// Streaming delimited-text parser. Pure functions, no DOM, so it runs in the
// worker and in Node for testing.

export const NA_TOKENS = new Set([
  '', 'NA', 'na', 'N/A', 'n/a', '#N/A', 'NaN', 'nan', 'NAN', '.', '-',
  'null', 'NULL', 'None', 'none',
]);

const INITIAL_CAP = 4096;

// One column of parsed data. Every value is stored as a float (NaN when
// missing or non-numeric). String values are kept only once a column turns
// out to contain non-numeric tokens.
export class ColumnBuilder {
  constructor(name) {
    this.name = name;
    this.num = new Float64Array(INITIAL_CAP);
    this.strs = null;
    this.n = 0;
    this.na = 0;
    this.nonnum = 0;
    this.underflow = null; // Map(row -> -log10(value)) for values below double range
  }

  grow() {
    const next = new Float64Array(this.num.length * 2);
    next.set(this.num);
    this.num = next;
  }

  push(tok) {
    if (this.n === this.num.length) this.grow();
    const i = this.n++;
    if (NA_TOKENS.has(tok)) {
      this.num[i] = NaN;
      this.na++;
      if (this.strs) this.strs.push('');
      return;
    }
    let v = +tok;
    if (v !== v) {
      const lc = tok.toLowerCase();
      if (lc === 'inf' || lc === '+inf') v = Infinity;
      else if (lc === '-inf') v = -Infinity;
    }
    if (v !== v) {
      this.nonnum++;
      if (!this.strs) this.backfill(i);
      this.strs.push(tok);
      this.num[i] = NaN;
      return;
    }
    if (v === 0 && (tok.indexOf('e') > 0 || tok.indexOf('E') > 0)) {
      // e.g. "3.2e-450": parses to 0 but the real value is tiny and positive.
      const k = tok.search(/[eE]/);
      const mant = parseFloat(tok.slice(0, k));
      const ex = parseInt(tok.slice(k + 1), 10);
      if (mant > 0 && Number.isFinite(ex)) {
        if (!this.underflow) this.underflow = new Map();
        this.underflow.set(i, -(Math.log10(mant) + ex));
      }
    }
    this.num[i] = v;
    if (this.strs) this.strs.push(tok);
  }

  backfill(i) {
    this.strs = new Array(i);
    for (let j = 0; j < i; j++) {
      const v = this.num[j];
      this.strs[j] = v === v ? String(v) : '';
    }
  }

  finish() {
    const num = this.num.slice(0, this.n);
    return {
      name: this.name,
      num,
      strs: this.strs,
      na: this.na,
      nonnum: this.nonnum,
      underflow: this.underflow ? Array.from(this.underflow) : null,
    };
  }
}

function isNumericToken(t) {
  if (NA_TOKENS.has(t)) return true; // missing values don't argue against a numeric column
  const v = +t;
  return v === v;
}

function splitQuotedCSV(line, delim) {
  const out = [];
  let cur = '';
  let inQ = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (inQ) {
      if (c === '"') {
        if (line[i + 1] === '"') { cur += '"'; i++; } else inQ = false;
      } else cur += c;
    } else if (c === '"') inQ = true;
    else if (c === delim) { out.push(cur.trim()); cur = ''; }
    else cur += c;
  }
  out.push(cur.trim());
  return out;
}

export function detectDelimiter(line) {
  if (line.indexOf('\t') >= 0) return '\t';
  if (line.indexOf(',') >= 0) return ',';
  if (line.indexOf(';') >= 0 && !/\s/.test(line.trim())) return ';';
  return 'ws';
}

export function makeSplitter(delim) {
  if (delim === 'ws') {
    const re = /\s+/;
    return (line) => line.trim().split(re);
  }
  return (line) => {
    if (line.indexOf('"') >= 0) return splitQuotedCSV(line, delim);
    const parts = line.split(delim);
    for (let i = 0; i < parts.length; i++) parts[i] = parts[i].trim();
    return parts;
  };
}

function uniqueNames(fields) {
  const seen = new Map();
  return fields.map((f, j) => {
    let name = f.replace(/^["']|["']$/g, '').trim() || `V${j + 1}`;
    const base = name;
    let k = seen.get(base) || 0;
    if (k > 0) name = `${base}.${k}`;
    seen.set(base, k + 1);
    return name;
  });
}

// Feed text chunks with push(); call finish() for the result.
export class DelimitedParser {
  constructor() {
    this.carry = '';
    this.pending = []; // lines buffered until the header is decided
    this.cols = null;
    this.split = null;
    this.delim = null;
    this.hasHeader = true;
    this.rows = 0;
    this.ragged = 0;
    this.raggedExample = null;
    this.lineNo = 0;
    this.comments = 0;
  }

  push(chunk) {
    const text = this.carry + chunk;
    const lines = text.split('\n');
    this.carry = lines.pop();
    for (let i = 0; i < lines.length; i++) this.line(lines[i]);
  }

  line(raw) {
    this.lineNo++;
    let line = raw;
    if (line.charCodeAt(line.length - 1) === 13) line = line.slice(0, -1);
    if (this.lineNo === 1 && line.charCodeAt(0) === 0xfeff) line = line.slice(1);
    if (line.trim() === '') return;
    if (this.cols) {
      if (line.charCodeAt(0) === 35) { this.comments++; return; } // '#'
      this.row(this.split(line));
      return;
    }
    if (line.startsWith('##')) { this.comments++; return; }
    if (this.pending.length > 0 && line.charCodeAt(0) === 35) { this.comments++; return; }
    this.pending.push(line);
    if (this.pending.length >= 2) this.init();
  }

  init() {
    let first = this.pending[0];
    const hashHeader = first.charCodeAt(0) === 35;
    if (hashHeader) first = first.replace(/^#+\s*/, '');
    this.delim = detectDelimiter(first);
    this.split = makeSplitter(this.delim);
    const f1 = this.split(first);
    const f2 = this.pending.length > 1 ? this.split(this.pending[1]) : null;
    let header;
    if (hashHeader) header = true;
    else if (f2) header = f1.some((t, j) => !isNumericToken(t) && j < f2.length && isNumericToken(f2[j]));
    else header = f1.some((t) => !isNumericToken(t));
    this.hasHeader = header;
    const names = header ? uniqueNames(f1) : f1.map((_, j) => `V${j + 1}`);
    this.cols = names.map((n) => new ColumnBuilder(n));
    const rest = this.pending.slice(header ? 1 : 0);
    this.pending = [];
    if (!header) this.row(f1);
    for (let i = header ? 0 : 1; i < rest.length; i++) this.row(this.split(rest[i]));
  }

  row(fields) {
    const cols = this.cols;
    const nc = cols.length;
    if (fields.length !== nc) {
      // PLINK and friends sometimes pad lines; only count real mismatches.
      this.ragged++;
      if (this.raggedExample === null) this.raggedExample = { line: this.lineNo, fields: fields.length };
    }
    for (let j = 0; j < nc; j++) {
      const t = j < fields.length ? fields[j] : '';
      cols[j].push(t);
    }
    this.rows++;
  }

  finish() {
    if (this.carry) { this.line(this.carry); this.carry = ''; }
    if (!this.cols && this.pending.length) this.init();
    if (!this.cols) throw new Error('No data found. Is the file empty?');
    return {
      columns: this.cols.map((c) => c.finish()),
      nrows: this.rows,
      delimiter: this.delim,
      hasHeader: this.hasHeader,
      ragged: this.ragged,
      raggedExample: this.raggedExample,
      expectedFields: this.cols.length,
      comments: this.comments,
    };
  }
}

export function parseText(text) {
  const p = new DelimitedParser();
  p.push(text);
  return p.finish();
}
