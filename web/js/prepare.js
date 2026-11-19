// Turn parsed columns + a column mapping into plot-ready, sorted arrays, and
// report what we found along the way.
import { fmtInt, fmtPFromLog, lambdaGC } from './stats.js';

const NAMES = {
  chr: ['chr', 'chrom', 'chromosome', 'chr_name', 'chrom_name', 'hg19chrc', 'hg38chrc', 'chr_id', 'contig', 'seqnames', 'seqname', 'chromosome_name', 'chrname'],
  bp: ['bp', 'pos', 'position', 'base_pair_location', 'basepair', 'base_pair', 'genpos', 'bp_hg19', 'bp_hg38', 'pos_hg19', 'pos_hg38', 'chr_position', 'start', 'location', 'physical_position', 'ps'],
  p: ['p', 'pval', 'p_value', 'p.value', 'p-value', 'pvalue', 'p_val', 'pvals', 'p_bolt_lmm', 'p_bolt_lmm_inf', 'p_lrt', 'p_wald', 'p_score', 'frequentist_add_pvalue', 'gc_pvalue', 'p_gc', 'all_inv_var_meta_p', 'p-val', 'pr(>|t|)', 'pr(>|z|)', 'p_dgc', 'p_linreg', 'pval_nominal', 'p.val', 'neg_log_10_p_value', 'log10p', 'mlog10p', 'neglog10p', 'neglog10_p', 'log10_p', 'minus_log10_p', '-log10p', 'mlogp', 'logp', 'lp'],
  snp: ['snp', 'rsid', 'rs', 'rs_id', 'id', 'snpid', 'snp_id', 'marker', 'markername', 'marker_name', 'variant', 'variant_id', 'variantid', 'name', 'rsids', 'snp_name', 'snpname', 'varid', 'rsnumber', 'rs_number'],
};
const LOG_NAMES = /^(neg_?)?(minus_?)?-?(log10_?p|mlog10p?|log10p|neglog10_?p|mlogp|logp|lp|neg_log_10_p_value)$/i;

function norm(s) { return s.trim().toLowerCase().replace(/^#/, ''); }

export function guessMapping(raw) {
  const names = raw.columns.map((c) => norm(c.name));
  const pick = (key, avoid) => {
    for (const cand of NAMES[key]) {
      const j = names.indexOf(cand);
      if (j >= 0 && !avoid.includes(j)) return j;
    }
    return -1;
  };
  const used = [];
  const chr = pick('chr', used); used.push(chr);
  const bp = pick('bp', used); used.push(bp);
  let p = -1;
  // Prefer a plain p-value column over a -log10 one if both exist.
  for (const cand of NAMES.p) {
    const j = names.indexOf(cand);
    if (j >= 0 && !used.includes(j)) { p = j; break; }
  }
  used.push(p);
  const snp = pick('snp', used);
  const m = { chr, bp, p, snp, scale: 'p' };
  if (p >= 0 && LOG_NAMES.test(names[p])) m.scale = 'mlog10';
  // Headerless or unusual names: fall back to column contents.
  if (m.chr < 0 || m.bp < 0 || m.p < 0) fallbackByContent(raw, m);
  return m;
}

function fallbackByContent(raw, m) {
  const cols = raw.columns;
  const taken = () => [m.chr, m.bp, m.p, m.snp];
  const sample = (c) => {
    const n = Math.min(c.num.length, 2000);
    let numeric = 0, unit = 0, ints = 0, big = 0, small = 0, max = -Infinity;
    for (let i = 0; i < n; i++) {
      const v = c.num[i];
      if (v === v) {
        numeric++;
        if (v >= 0 && v <= 1) unit++;
        if (Number.isInteger(v)) ints++;
        if (v > 1e4) big++;
        if (v >= 1 && v <= 30) small++;
        if (v > max) max = v;
      }
    }
    return { n, numeric, unit, ints, big, small, max };
  };
  const stats = cols.map(sample);
  if (m.chr < 0) {
    const j = stats.findIndex((s, k) => !taken().includes(k) && s.ints === s.numeric && s.numeric > 0 && s.max <= 30);
    if (j >= 0) m.chr = j;
  }
  if (m.bp < 0) {
    // the all-integer column with the largest values
    let best = -1;
    stats.forEach((s, k) => {
      if (taken().includes(k) || s.numeric === 0 || s.ints !== s.numeric) return;
      if (best < 0 || s.max > stats[best].max) best = k;
    });
    if (best >= 0 && stats[best].max > 30) m.bp = best;
  }
  if (m.p < 0) {
    // Last all-in-[0,1] column is usually the p-value (after allele frequencies etc.).
    for (let j = cols.length - 1; j >= 0; j--) {
      const s = stats[j];
      if (!taken().includes(j) && s.numeric > 0 && s.unit === s.numeric && s.ints < s.numeric) { m.p = j; break; }
    }
  }
  if (m.snp < 0) {
    const j = cols.findIndex((c, k) => !taken().includes(k) && c.strs && c.nonnum > c.num.length * 0.5 && /^(rs|ss|\d+[:_])/i.test(c.strs[0] || ''));
    if (j >= 0) m.snp = j;
  }
}

// ---------------------------------------------------------------------------
// Chromosomes

const PLINK = { 23: 'X', 24: 'Y', 25: 'XY', 26: 'MT' };

function chromFromToken(tok) {
  const s = String(tok).trim();
  const t = s.replace(/^(chr(om(osome)?)?)[_\s-]?/i, '');
  if (/^\d+$/.test(t)) {
    const n = parseInt(t, 10);
    return { key: `n${n}`, rank: n, label: String(n), known: true, letter: false };
  }
  const u = t.toUpperCase();
  if (u === 'X') return { key: 'n23', rank: 23, label: 'X', known: true, letter: true };
  if (u === 'Y') return { key: 'n24', rank: 24, label: 'Y', known: true, letter: true };
  if (u === 'XY' || u === 'PAR') return { key: 'n25', rank: 25, label: 'XY', known: true, letter: true };
  if (u === 'M' || u === 'MT') return { key: 'n26', rank: 26, label: 'MT', known: true, letter: true };
  return { key: `s${t}`, rank: Infinity, label: t || s, known: false, letter: true };
}

export function chromLabel(ch, opts) {
  if (opts && opts.plinkNames && !ch.letter && PLINK[ch.rank]) return PLINK[ch.rank];
  return ch.label;
}

// ---------------------------------------------------------------------------

const C = (level, title, detail) => ({ level, title, detail });
// "1 row with ... was dropped" / "3 rows with ... were dropped"
const droppedMsg = (n, what) => `${fmtInt(n)} row${n === 1 ? '' : 's'} with ${what} ${n === 1 ? 'was' : 'were'} dropped`;

export function prepare(raw, map) {
  const checks = [];
  const cols = raw.columns;
  const n = raw.nrows;

  if (raw.ragged > 0) {
    const ex = raw.raggedExample;
    checks.push(C('warn', `${fmtInt(raw.ragged)} line${raw.ragged > 1 ? 's' : ''} had the wrong number of fields`,
      `Expected ${raw.expectedFields} fields per line; line ${ex.line} had ${ex.fields}. Missing fields were treated as NA.`));
  }
  if (!raw.hasHeader) {
    checks.push(C('info', 'No header row detected', 'Columns were named V1, V2, … Pick the right ones under Columns.'));
  }

  const need = [['chr', 'chromosome'], ['bp', 'position'], ['p', 'p-value']];
  for (const [k, label] of need) {
    if (map[k] < 0 || map[k] === undefined || map[k] === null) {
      checks.push(C('error', `Choose a ${label} column`, `Couldn't find a ${label} column automatically. Select one under Columns.`));
    }
  }
  if (checks.some((c) => c.level === 'error')) return { data: null, checks };
  const dupes = [map.chr, map.bp, map.p].filter((v, i, a) => a.indexOf(v) !== i);
  if (dupes.length) {
    checks.push(C('error', 'The same column is used twice', 'CHR, BP and P each need their own column.'));
    return { data: null, checks };
  }

  const chrCol = cols[map.chr], bpCol = cols[map.bp], pCol = cols[map.p];
  const snpCol = map.snp >= 0 ? cols[map.snp] : null;
  const scale = map.scale;

  // --- chromosomes -------------------------------------------------------
  const chrId = new Int32Array(n).fill(-1);
  const chromMap = new Map(); // key -> chrom
  const chromList = [];
  const tokenCache = new Map();
  let nonIntChr = 0;
  let chrNA = 0;
  const getChrom = (tok) => {
    let c = tokenCache.get(tok);
    if (c !== undefined) return c;
    const info = chromFromToken(tok);
    c = chromMap.get(info.key);
    if (!c) {
      c = { ...info, id: chromList.length, count: 0 };
      chromMap.set(info.key, c);
      chromList.push(c);
    } else if (info.letter && !c.letter) {
      c.letter = true; c.label = info.label;
    }
    tokenCache.set(tok, c);
    return c;
  };
  const cn = chrCol.num, cs = chrCol.strs;
  for (let i = 0; i < n; i++) {
    const v = cn[i];
    let c;
    if (v === v) {
      if (!Number.isInteger(v) || v < 0) { nonIntChr++; c = getChrom(String(v)); }
      else c = getChrom(v);
      if (typeof c.rank === 'number' && !Number.isInteger(v)) c.known = false;
    } else if (cs && cs[i] !== '') c = getChrom(cs[i]);
    else { chrNA++; continue; }
    chrId[i] = c.id;
  }
  // numeric tokens were cached by number; normalize their labels
  for (const c of chromList) if (!c.letter && c.known) c.label = String(c.rank);

  if (chromList.length === 0) {
    checks.push(C('error', `No usable values in ${chrCol.name}`, 'Every value in the chromosome column is missing.'));
    return { data: null, checks };
  }
  if (chromList.length > 2000) {
    checks.push(C('error', `${fmtInt(chromList.length)} distinct values in ${chrCol.name}`,
      `That's too many to be chromosomes. Did you pick the right column? (Examples: ${chromList.slice(0, 4).map((c) => c.label).join(', ')})`));
    return { data: null, checks };
  }
  if (nonIntChr > 0) {
    checks.push(C('warn', `${chrCol.name} has non-integer values`, `${fmtInt(nonIntChr)} values like "${chromList.find((c) => !c.known)?.label}" don't look like chromosome numbers. Check that you picked the right column.`));
  }
  if (chromList.length > 60) {
    checks.push(C('warn', `${fmtInt(chromList.length)} chromosomes`, 'That is a lot of chromosomes or contigs. If that is unexpected, check the CHR column.'));
  }
  const unknown = chromList.filter((c) => !c.known);
  if (unknown.length && nonIntChr === 0) {
    checks.push(C('info', `${unknown.length} non-standard chromosome name${unknown.length > 1 ? 's' : ''}`,
      `Plotted after the numbered chromosomes: ${unknown.slice(0, 5).map((c) => c.label).join(', ')}${unknown.length > 5 ? ', …' : ''}`));
  }

  chromList.sort((a, b) => (a.rank !== b.rank ? a.rank - b.rank : a.label.localeCompare(b.label, undefined, { numeric: true })));
  const order = new Int32Array(chromList.length);
  chromList.forEach((c, k) => { order[c.id] = k; });

  // --- positions and p-values -------------------------------------------
  const bn = bpCol.num, pn = pCol.num;
  const pUnder = pCol.underflow ? new Map(pCol.underflow) : null;
  const y = new Float64Array(n);
  let bpBad = 0, bpNeg = 0, bpFrac = 0;
  let pNA = 0, pNonNum = 0, pHigh = 0, pNeg = 0, pZero = 0, pUnderflow = 0, pInf = 0;
  const keep = new Uint8Array(n);
  let minPos = Infinity;
  for (let i = 0; i < n; i++) {
    if (chrId[i] < 0) continue;
    const b = bn[i];
    if (!(b === b) || !Number.isFinite(b)) { bpBad++; continue; }
    if (b < 0) { bpNeg++; continue; }
    if (!Number.isInteger(b)) bpFrac++;
    const v = pn[i];
    if (!(v === v)) {
      if (pCol.strs && pCol.strs[i] !== '') pNonNum++; else pNA++;
      continue;
    }
    let val;
    if (scale === 'p') {
      if (v > 1) { pHigh++; continue; }
      if (v < 0) { pNeg++; continue; }
      if (v === 0) {
        const u = pUnder && pUnder.get(i);
        if (u !== undefined) { val = u; pUnderflow++; } else { pZero++; val = Infinity; }
      } else val = -Math.log10(v) + 0; // + 0 turns -0 (p = 1) into 0
    } else if (scale === 'mlog10') {
      if (v < 0) { pNeg++; continue; }
      if (!Number.isFinite(v)) { pInf++; continue; }
      val = v;
    } else {
      if (!Number.isFinite(v)) { pInf++; continue; }
      val = v;
    }
    y[i] = val;
    keep[i] = 1;
    if (b < minPos) minPos = b;
  }

  if (chrNA) checks.push(C('warn', droppedMsg(chrNA, 'a missing chromosome'), ''));
  if (bpBad) checks.push(C('warn', droppedMsg(bpBad, 'a missing or non-numeric position'), `Column ${bpCol.name}.`));
  if (bpNeg) checks.push(C('warn', droppedMsg(bpNeg, 'a negative position'), `Column ${bpCol.name}.`));
  if (bpFrac && bpFrac > n * 0.5) checks.push(C('warn', `${bpCol.name} values aren't whole numbers`, 'Base-pair positions are usually integers. Check that you picked the right column.'));
  if (pNA) checks.push(C('warn', droppedMsg(pNA, 'a missing p-value'), `Values like NA, ".", or blank in ${pCol.name}.`));
  if (pNonNum) {
    const ex = pCol.strs.find((s, i) => s && !(pn[i] === pn[i]));
    checks.push(C('warn', droppedMsg(pNonNum, 'a non-numeric p-value'), `For example "${ex}".`));
  }
  if (pNeg) checks.push(C('warn', droppedMsg(pNeg, `a negative value in ${pCol.name}`), scale === 'p' ? 'p-values must be between 0 and 1.' : '-log10(p) values must be ≥ 0.'));
  if (pInf) checks.push(C('warn', droppedMsg(pInf, `an infinite value in ${pCol.name}`), ''));

  let kept = 0;
  for (let i = 0; i < n; i++) kept += keep[i];

  if (pHigh) {
    const frac = pHigh / Math.max(1, pHigh + kept);
    if (frac > 0.2) {
      checks.push(C('error', `${Math.round(frac * 100)}% of values in ${pCol.name} are greater than 1`,
        `p-values must be between 0 and 1. If ${pCol.name} already holds −log10(p), set "P scale" to −log10(p). If it holds scores or test statistics, use "Raw value".`));
      return { data: null, checks };
    }
    checks.push(C('warn', droppedMsg(pHigh, 'p > 1'), 'p-values must be between 0 and 1.'));
  }
  if (kept === 0) {
    checks.push(C('error', 'No rows left to plot', 'Every row was dropped by the checks above.'));
    return { data: null, checks };
  }

  // --- sort by chromosome then position (counting sort on chromosome) ----
  const nch = chromList.length;
  const counts = new Int32Array(nch);
  for (let i = 0; i < n; i++) if (keep[i]) counts[order[chrId[i]]]++;
  const starts = new Int32Array(nch + 1);
  for (let k = 0; k < nch; k++) starts[k + 1] = starts[k] + counts[k];
  const fill = starts.slice(0, nch);
  const orig = new Uint32Array(kept);
  for (let i = 0; i < n; i++) if (keep[i]) orig[fill[order[chrId[i]]]++] = i;
  let resorted = 0;
  for (let k = 0; k < nch; k++) {
    const a = starts[k], b = starts[k + 1];
    let sorted = true;
    for (let i = a + 1; i < b; i++) if (bn[orig[i]] < bn[orig[i - 1]]) { sorted = false; break; }
    if (!sorted) {
      resorted++;
      const seg = Array.from(orig.subarray(a, b));
      seg.sort((u, v) => bn[u] - bn[v]);
      orig.set(seg, a);
    }
  }

  const D = {
    n: kept,
    scale,
    orig,
    bp: new Float64Array(kept),
    y: new Float64Array(kept),
    pos: new Float64Array(kept),
    chrIdx: new Uint16Array(kept),
    chroms: [],
    snpCol: map.snp >= 0 ? map.snp : -1,
    columns: cols,
  };
  for (let k = 0; k < nch; k++) {
    const c = chromList[k];
    const a = starts[k], b = starts[k + 1];
    D.chroms.push({ label: c.label, rank: c.rank, letter: c.letter, known: c.known, from: a, to: b, count: b - a, minbp: bn[orig[a]], maxbp: bn[orig[b - 1]] });
  }
  let maxFinite = -Infinity, minY = Infinity;
  for (let i = 0; i < kept; i++) {
    const o = orig[i];
    D.bp[i] = bn[o];
    D.y[i] = y[o];
    if (y[o] !== Infinity && y[o] > maxFinite) maxFinite = y[o];
    if (y[o] < minY) minY = y[o];
  }
  for (let k = 0; k < nch; k++) D.chrIdx.fill(k, D.chroms[k].from, D.chroms[k].to);

  if (pZero) {
    const cap = Number.isFinite(maxFinite) ? Math.max(maxFinite, 1) : 300;
    for (let i = 0; i < kept; i++) if (D.y[i] === Infinity) D.y[i] = cap;
    checks.push(C('warn', `${fmtInt(pZero)} p-value${pZero > 1 ? 's are' : ' is'} exactly 0`,
      `−log10(0) is infinite, so ${pZero > 1 ? 'these were' : 'it was'} plotted at the smallest non-zero p-value in the data (${fmtPFromLog(cap)}). Your software probably rounded a very small p-value; check the original output.`));
  }
  if (pUnderflow) {
    checks.push(C('info', `${fmtInt(pUnderflow)} p-value${pUnderflow > 1 ? 's' : ''} below 1e-308`,
      'These are smaller than a double can hold. They were read from the text directly, so their −log10(p) is exact.'));
  }

  if (map.snp < 0) {
    checks.push(C('info', 'No SNP column selected', 'Hover and labels will show chromosome:position instead. Highlighting by SNP name needs a SNP column.'));
  } else {
    const sc = cols[map.snp];
    if (sc.strs || sc.num.length) {
      const seen = new Set();
      let dup = 0, ex = null;
      for (let i = 0; i < kept; i++) {
        const s = snpName(D, i);
        if (s === '' || s === '.') continue;
        if (seen.has(s)) { dup++; if (ex === null) ex = s; } else seen.add(s);
      }
      if (dup) checks.push(C('info', `${fmtInt(dup)} duplicated SNP name${dup > 1 ? 's' : ''}`, `For example ${ex}. Multi-allelic sites often do this. It only matters for highlighting.`));
    }
  }
  if (resorted) checks.push(C('info', 'Rows were sorted by position', `Positions within ${resorted} chromosome${resorted > 1 ? 's' : ''} weren't in order, so they were sorted.`));

  const dropped = n - kept;
  const summary = `${fmtInt(kept)} of ${fmtInt(n)} rows plotted${dropped ? `, ${fmtInt(dropped)} dropped` : ''}.`;
  const range = chromList.length === 1 ? `Chromosome ${chromList[0].label}` : `${chromList.length} chromosomes: ${chromSummary(D.chroms)}`;
  checks.unshift(C('ok', summary, range));
  if (scale === 'p') {
    checks.push(C('ok', 'p-values are in range', `Smallest p = ${fmtPFromLog(maxFinite)}${pZero ? ' (excluding zeros)' : ''}.`));
  } else if (scale === 'mlog10') {
    if (maxFinite <= 1 && minY >= 0) checks.push(C('warn', `${pCol.name} values are all between 0 and 1`, 'These look like p-values rather than −log10(p). If so, set "P scale" to p-value.'));
    else checks.push(C('ok', 'Using values as −log10(p)', `Range ${minY.toFixed(2)} to ${maxFinite.toFixed(2)}, smallest p = ${fmtPFromLog(maxFinite)}.`));
  } else {
    checks.push(C('info', 'Plotting raw values', `No −log10 transform (like qqman's logp = FALSE). Range ${+minY.toPrecision(4)} to ${+maxFinite.toPrecision(4)}. The Q-Q plot needs p-values, so it's turned off.`));
  }

  D.maxY = maxFinite;
  D.minY = minY;
  computeSummary(D);
  if (kept > 100) {
    const s = D.sortedY;
    const tenth = s[kept - 10];
    if (s[kept - 1] > 40 && s[kept - 1] > 3 * Math.max(tenth, 8)) {
      checks.push(C('info', 'A few extreme values squash the plot',
        `The top value is ${s[kept - 1].toFixed(1)} on the y axis but the 10th highest is ${tenth.toFixed(1)}. Set "Y axis max" under Axes and text to cap the axis; points above the cap are drawn as triangles at the top.`));
    }
  }
  return { data: D, checks };
}

function chromSummary(chroms) {
  const labels = chroms.map((c) => c.label);
  if (labels.length <= 6) return labels.join(', ');
  return `${labels.slice(0, 3).join(', ')}, …, ${labels.slice(-3).join(', ')}`;
}

export function snpName(D, i) {
  if (D.snpCol < 0) return '';
  const c = D.columns[D.snpCol];
  const o = D.orig[i];
  if (c.strs) return c.strs[o];
  const v = c.num[o];
  return v === v ? String(v) : '';
}

export function cellText(col, o) {
  if (col.strs && col.strs[o] !== undefined && col.num[o] !== col.num[o]) return col.strs[o];
  const v = col.num[o];
  if (v !== v) return col.strs ? col.strs[o] || 'NA' : 'NA';
  return null; // numeric: caller formats
}

// Genome coordinates for each SNP, with an optional gap between chromosomes.
export function computePositions(D, gapFrac) {
  let span = 0;
  for (const c of D.chroms) span += c.maxbp - c.minbp;
  const minW = Math.max(1, span * 0.002);
  const gap = span * (gapFrac || 0);
  let off = 0;
  for (const c of D.chroms) {
    const w = Math.max(c.maxbp - c.minbp, minW);
    c.offset = off;
    c.start = off;
    c.end = off + w;
    const base = off - c.minbp + (w - (c.maxbp - c.minbp)) / 2;
    c.base = base;
    for (let i = c.from; i < c.to; i++) D.pos[i] = D.bp[i] + base;
    off += w + gap;
  }
  D.genomeStart = 0;
  D.genomeEnd = off - gap;
}

function computeSummary(D) {
  const sorted = D.y.slice().sort(); // ascending
  D.sortedY = sorted;
  const n = D.n;
  if (D.scale !== 'raw') {
    const mid = n % 2 ? sorted[(n - 1) / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
    D.lambda = lambdaGC(10 ** -mid);
  } else D.lambda = null;
  // Top hits: indices of the K largest values, sorted descending.
  const K = Math.min(n, 2000);
  const thr = sorted[n - K];
  const idx = [];
  for (let i = 0; i < n; i++) if (D.y[i] >= thr) idx.push(i);
  idx.sort((a, b) => D.y[b] - D.y[a] || a - b);
  D.top = idx.slice(0, K);
}
