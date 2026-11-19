// Parsing and simulation run here so the page stays responsive.
import { DelimitedParser, parseText } from './parse.js';
import { mulberry32, normalGen, zToMlog10p } from './stats.js';

let currentId = 0;

self.onmessage = async (ev) => {
  const msg = ev.data;
  currentId = msg.id;
  const t0 = performance.now();
  try {
    let result;
    if (msg.kind === 'file') result = await parseFile(msg.file, msg.id);
    else if (msg.kind === 'url') result = await parseUrl(msg.url, msg.id);
    else if (msg.kind === 'text') result = parseText(msg.text);
    else if (msg.kind === 'sim') result = msg.which === 'messy' ? parseText(simulateMessy()) : simulateGWAS(msg.n || 1_000_000, msg.id);
    else throw new Error(`Unknown request: ${msg.kind}`);
    if (msg.id !== currentId) return; // superseded by a newer request
    result.ms = performance.now() - t0;
    post(msg.id, 'done', { result }, result.columns.map((c) => c.num.buffer));
  } catch (err) {
    post(msg.id, 'error', { message: err && err.message ? err.message : String(err) });
  }
};

function post(id, type, payload, transfer) {
  self.postMessage({ id, type, ...payload }, transfer || []);
}

function progress(id, frac, rows) {
  post(id, 'progress', { frac, rows });
}

async function isGzip(blob) {
  const head = new Uint8Array(await blob.slice(0, 2).arrayBuffer());
  return head[0] === 0x1f && head[1] === 0x8b;
}

async function parseFile(file, id) {
  let stream = file.stream();
  let loaded = 0;
  const total = file.size || 1;
  const counter = new TransformStream({
    transform(chunk, ctrl) { loaded += chunk.byteLength; ctrl.enqueue(chunk); },
  });
  stream = stream.pipeThrough(counter);
  if (await isGzip(file)) {
    if (typeof DecompressionStream === 'undefined') throw new Error('This browser cannot read .gz files. Decompress the file first.');
    stream = stream.pipeThrough(new DecompressionStream('gzip'));
  }
  if (/\.(zip|xlsx?|parquet|bgen|bed|bim)$/i.test(file.name)) {
    throw new Error(`Can't read ${file.name.split('.').pop()} files. Use a plain text, CSV/TSV, or gzipped text file.`);
  }
  return parseStream(stream, () => loaded / total, id);
}

async function parseUrl(url, id) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Could not fetch ${url} (${res.status})`);
  const total = +res.headers.get('content-length') || 0;
  let loaded = 0;
  let stream = res.body.pipeThrough(new TransformStream({
    transform(chunk, ctrl) { loaded += chunk.byteLength; ctrl.enqueue(chunk); },
  }));
  if (/\.gz$/i.test(url)) stream = stream.pipeThrough(new DecompressionStream('gzip'));
  return parseStream(stream, () => (total ? loaded / total : 0), id);
}

async function parseStream(byteStream, frac, id) {
  const parser = new DelimitedParser();
  const reader = byteStream.pipeThrough(new TextDecoderStream()).getReader();
  let last = 0;
  for (;;) {
    const { value, done } = await reader.read();
    if (done) break;
    parser.push(value);
    const now = performance.now();
    if (now - last > 80) {
      last = now;
      if (id !== currentId) { reader.cancel(); throw new Error('cancelled'); }
      progress(id, frac(), parser.rows);
    }
  }
  return parser.finish();
}

// ---------------------------------------------------------------------------
// Simulated data

const HG38 = [
  ['1', 248956422], ['2', 242193529], ['3', 198295559], ['4', 190214555], ['5', 181538259],
  ['6', 170805979], ['7', 159345973], ['8', 145138636], ['9', 138394717], ['10', 133797422],
  ['11', 135086622], ['12', 133275309], ['13', 114364328], ['14', 107043718], ['15', 101991189],
  ['16', 90338345], ['17', 83257441], ['18', 80373285], ['19', 58617616], ['20', 64444167],
  ['21', 46709983], ['22', 50818468], ['X', 156040895],
];

function numCol(name, num) {
  return { name, num, strs: null, na: 0, nonnum: 0, underflow: null };
}

// A realistic-looking GWAS: ~1M SNPs, mild inflation, a dozen association
// peaks with LD-like decay around each lead SNP.
function simulateGWAS(n, id) {
  const rng = mulberry32(20140514);
  const rnorm = normalGen(rng);
  const total = HG38.reduce((s, c) => s + c[1], 0);
  const loci = [
    ['2', 60.7e6, 9.5], ['3', 141.1e6, 6.2], ['6', 32.6e6, 14.5], ['6', 26.1e6, 7.4],
    ['8', 127.4e6, 8.2], ['10', 112.9e6, 10.8], ['11', 2.1e6, 6.9], ['12', 111.8e6, 5.6],
    ['16', 53.8e6, 8.9], ['17', 45.9e6, 6.4], ['19', 44.9e6, 12.4], ['20', 44.4e6, 5.3],
    ['22', 30.4e6, 6.0], ['X', 67.5e6, 6.6], ['14', 75.0e6, 5.0],
  ].map(([chr, pos, z]) => ({ chr, pos, z, w: 60e3 + rng() * 240e3 }));

  const SNP = new Array(n);
  const CHR = new Array(n);
  const chrNum = new Float64Array(n);
  const BP = new Float64Array(n);
  const A1 = new Array(n);
  const A2 = new Array(n);
  const BETA = new Float64Array(n);
  const SE = new Float64Array(n);
  const P = new Float64Array(n);
  const bases = ['A', 'C', 'G', 'T'];
  const sdNull = Math.sqrt(1.035);
  let k = 0;
  let nonnum = 0;
  for (let c = 0; c < HG38.length; c++) {
    const [chr, len] = HG38[c];
    const m = c === HG38.length - 1 ? n - k : Math.round((n * len) / total);
    const here = loci.filter((l) => l.chr === chr);
    const step = (len - 2e5) / m;
    for (let i = 0; i < m && k < n; i++, k++) {
      const pos = Math.floor(1e5 + (i + rng()) * step);
      let z = rnorm() * sdNull;
      for (const l of here) {
        const d = Math.abs(pos - l.pos);
        if (d < l.w * 8) {
          // r^2 with the lead SNP decays with distance; the expected z scales with r
          const r2 = Math.exp(-d / l.w) * (0.15 + 0.85 * rng());
          z += l.z * Math.sqrt(r2);
        }
      }
      const maf = 0.01 + rng() * 0.49;
      const se = 1 / Math.sqrt(2 * 60000 * maf * (1 - maf));
      const a = (rng() * 4) | 0;
      SNP[k] = `rs${1000 + k * 53 + ((rng() * 53) | 0)}`;
      CHR[k] = chr;
      chrNum[k] = chr === 'X' ? NaN : +chr;
      if (chr === 'X') nonnum++;
      BP[k] = pos;
      A1[k] = bases[a];
      A2[k] = bases[(a + 1 + ((rng() * 3) | 0)) % 4];
      SE[k] = +se.toPrecision(4);
      BETA[k] = +(z * se).toPrecision(4);
      const L = zToMlog10p(z);
      P[k] = 10 ** -L;
      if (k % 100000 === 0) progress(id, k / n, k);
    }
  }
  const len = k;
  const sl = (a) => (a instanceof Float64Array ? a.slice(0, len) : a.slice(0, len));
  const strCol = (name, arr) => ({ name, num: new Float64Array(len).fill(NaN), strs: sl(arr), na: 0, nonnum: len, underflow: null });
  const chrCol = { name: 'CHR', num: sl(chrNum), strs: sl(CHR), na: 0, nonnum, underflow: null };
  return {
    columns: [strCol('SNP', SNP), chrCol, numCol('BP', sl(BP)), strCol('A1', A1), strCol('A2', A2), numCol('BETA', sl(BETA)), numCol('SE', sl(SE)), numCol('P', sl(P))],
    nrows: len, delimiter: 'sim', hasHeader: true, ragged: 0, raggedExample: null, expectedFields: 8, comments: 0,
  };
}

// A small file with the usual problems: chr prefixes, contigs, missing and
// out-of-range p-values, p = 0, an underflowing p-value, duplicate IDs, and a
// couple of truncated lines.
function simulateMessy() {
  const rng = mulberry32(7);
  const rnorm = normalGen(rng);
  const chroms = [...HG38.map((c) => c[0]), 'Y', 'M'];
  const lines = ['#CHROM\tPOS\tID\tREF\tALT\tPVAL'];
  let id = 1;
  for (const c of chroms) {
    const m = c === 'Y' || c === 'M' ? 25 : 220;
    const len = (HG38.find((h) => h[0] === c) || [0, 16e6])[1];
    for (let i = 0; i < m; i++) {
      const pos = Math.floor(((i + rng()) / m) * len);
      let z = rnorm();
      if (c === '6' && Math.abs(i - 60) < 8) z += 7 - Math.abs(i - 60) * 0.7;
      const p = 10 ** -zToMlog10p(z);
      lines.push(`chr${c}\t${pos}\trs${(id++) * 7}\tA\tG\t${p.toPrecision(4)}`);
    }
  }
  const at = (i) => lines.length > i ? i : lines.length - 1;
  const set = (i, f, v) => { const parts = lines[at(i)].split('\t'); parts[f] = v; lines[at(i)] = parts.join('\t'); };
  for (let i = 0; i < 24; i++) set(30 + i * 97, 5, i % 3 === 0 ? 'NA' : i % 3 === 1 ? '.' : '');
  set(210, 5, '1.7'); set(1210, 5, '2.31'); set(2410, 5, '-0.2');
  set(1320, 5, '0'); set(3020, 5, 'error');
  set(1324, 5, '3.1e-412');
  set(801, 1, 'NA'); set(1702, 1, '-500');
  for (let i = 0; i < 6; i++) set(400 + i * 300, 2, 'rs7');
  lines.splice(2500, 0, 'chr12\t1000');
  lines.splice(3500, 0, 'chr13\t52000\trs999999');
  for (let i = 0; i < 18; i++) lines.push(`chrUn_KI2707${30 + i}v1\t${1000 + i * 77}\t.\tC\tT\t${rng().toPrecision(3)}`);
  return lines.join('\n') + '\n';
}
