import { prepare, guessMapping, computePositions, snpName, chromLabel, cellText } from './prepare.js';
import { renderManhattan, renderQQ, CanvasPainter, SvgPainter, thresholdY, lineSpec, plainText, bpUnit } from './plot.js';
import { fmtPFromLog, fmtNum, fmtInt } from './stats.js';
import { rCode, highlightR } from './rcode.js';

const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => Array.from(el.querySelectorAll(s));

// ---------------------------------------------------------------------------
// Options

const PALETTES = {
  qqman: ['#1a1a1a', '#999999'],
  slate: ['#1e3a5f', '#7fb3e0'],
  navy: ['#1f3b73', '#e0812b'],
  okabe: ['#0072b2', '#e69f00', '#009e73', '#cc79a7', '#56b4e9', '#d55e00', '#000000'],
  viridis: ['#440154', '#3b528b', '#21918c', '#5ec962', '#c9b200'],
  rainbow: Array.from({ length: 12 }, (_, i) => hsl((i * 360) / 12, 70, 45)),
};

function hsl(h, s, l) {
  s /= 100; l /= 100;
  const k = (n) => (n + h / 30) % 12;
  const a = s * Math.min(l, 1 - l);
  const f = (n) => l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
  return `#${[f(0), f(8), f(4)].map((x) => Math.round(x * 255).toString(16).padStart(2, '0')).join('')}`;
}

const DEFAULTS = {
  palette: 'qqman',
  colors: PALETTES.qqman.slice(),
  pointSize: 4,
  opacity: 1,
  chrGap: 0,
  shade: false,
  shadeColor: '#f0f0f0',
  suggestive: { on: true, p: '1e-5', raw: '', color: '#0000ff' },
  genomewide: { on: true, p: '5e-8', raw: '', color: '#ff0000' },
  lineStyle: 'solid',
  lineWidth: 1,
  highlightText: '',
  highlightColor: '#00cd00',
  highlightScale: 1,
  annotate: false,
  annotateP: '1e-5',
  annotateRaw: '',
  annotateTop: true,
  labelCol: '',
  labelSize: 10,
  title: '',
  xlab: '',
  ylab: '',
  yMax: '',
  fontSize: 12,
  fontFamily: 'sans',
  chrLabels: '',
  plinkNames: false,
  grid: false,
  width: 1000,
  height: 480,
  margin: { t: 20, r: 20, b: 52, l: 62 },
  bg: '#ffffff',
  transparent: false,
  exportScale: 2,
  qq: { color: '#1a1a1a', lineColor: '#ff0000', ci: true, ciColor: '#dedede', showLambda: true, highlight: false, title: '', xlab: '', ylab: '', width: 520, height: 520, yMax: '' },
};

const clone = (o) => JSON.parse(JSON.stringify(o));
function merge(base, over) {
  if (!over || typeof over !== 'object') return base;
  for (const k of Object.keys(base)) {
    if (!(k in over)) continue;
    const b = base[k], v = over[k];
    if (b && typeof b === 'object' && !Array.isArray(b)) merge(b, v);
    else if (Array.isArray(b) ? Array.isArray(v) : typeof v === typeof b) base[k] = v;
  }
  return base;
}
const getPath = (o, path) => path.split('.').reduce((a, k) => (a ? a[k] : undefined), o);
function setPath(o, path, v) {
  const ks = path.split('.');
  const last = ks.pop();
  ks.reduce((a, k) => a[k], o)[last] = v;
}

// Settings that only make sense for one dataset: reset when new data loads and
// never remembered between visits. Everything else is remembered.
const DATASET_OPTS = ['highlightText', 'yMax', 'chrLabels', 'qq.yMax'];
function resetDatasetOpts(o) {
  for (const k of DATASET_OPTS) setPath(o, k, getPath(DEFAULTS, k));
  return o;
}

const STORE = 'qqman-web-opts-v1';
function loadOpts() {
  const o = clone(DEFAULTS);
  try { merge(o, JSON.parse(localStorage.getItem(STORE) || 'null')); } catch (e) { /* ignore */ }
  return resetDatasetOpts(o);
}
let saveTimer = 0;
function saveOpts() {
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => {
    try { localStorage.setItem(STORE, JSON.stringify(resetDatasetOpts(clone(opts)))); } catch (e) { /* ignore */ }
  }, 300);
}

let opts = loadOpts();

// ---------------------------------------------------------------------------
// State

const state = {
  raw: null,
  map: null,
  data: null,
  checks: [],
  view: null,
  pinned: new Set(),
  highlight: [],
  highlightNames: [],
  tab: 'manhattan',
  example: null,
  dirty: { manhattan: true, qq: true },
  mhL: null,
  qqL: null,
  hover: null,
  soi: null,
};

// ---------------------------------------------------------------------------
// Loading

const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
let reqId = 0;
let reqMeta = null;

worker.onmessage = (ev) => {
  const m = ev.data;
  if (m.id !== reqId) return;
  if (m.type === 'progress') setProgress(m.frac, m.rows);
  else if (m.type === 'error') {
    hideProgress();
    if (m.message !== 'cancelled') toast(m.message, true);
  } else if (m.type === 'done') onParsed(m.result, reqMeta);
};
worker.onerror = (e) => { hideProgress(); toast(`Something went wrong reading the data: ${e.message || 'worker error'}`, true); };

function request(kind, payload, meta) {
  reqId++;
  reqMeta = meta;
  showProgress(meta.label || meta.name);
  worker.postMessage({ id: reqId, kind, ...payload });
}

function loadFile(file) {
  if (!file) return;
  if (file.size === 0) { toast('That file is empty.', true); return; }
  request('file', { file }, { name: file.name, size: file.size, label: `Reading ${file.name}` });
}

function loadText(text, name = 'Pasted data') {
  if (!text || !text.trim()) { toast('Nothing to plot yet. Paste some data first.', true); return; }
  request('text', { text }, { name, size: text.length, label: 'Reading pasted data' });
}

function loadExample(which) {
  if (which === 'gwasResults') request('url', { url: new URL('examples/gwasResults.tsv', location.href).href }, { name: 'gwasResults', example: 'gwasResults', label: 'Loading gwasResults' });
  else if (which === 'big') request('sim', { which: 'big', n: 1_000_000 }, { name: 'Simulated GWAS (1M SNPs)', example: 'big', label: 'Simulating 1,000,000 SNPs' });
  else if (which === 'messy') request('sim', { which: 'messy' }, { name: 'messy_example.tsv', example: 'messy', label: 'Loading messy example' });
}

function onParsed(result, meta) {
  setProgress(1, result.nrows, 'Preparing plot…');
  // Let the progress text paint before the (synchronous) prepare step.
  requestAnimationFrame(() => setTimeout(() => {
    try { finishLoad(result, meta); } catch (err) {
      console.error(err);
      hideProgress();
      toast(`Couldn't prepare the plot: ${err.message}`, true);
    }
  }, 0));
}

function finishLoad(result, meta) {
  state.raw = { ...result, name: meta.name, size: meta.size };
  state.example = meta.example || null;
  state.map = guessMapping(state.raw);
  if (state.example === 'gwasResults' && !state.soi) {
    fetch('examples/snpsOfInterest.txt').then((r) => r.text()).then((t) => { state.soi = t.split(/\s+/).filter(Boolean); updateSoiButton(); }).catch(() => {});
  }
  resetDatasetOpts(opts);
  syncOptionInputs();
  fillColumnSelects();
  runPrepare(true);
  hideProgress();
  showFileInfo();
  $$('.example').forEach((b) => b.classList.toggle('active', b.dataset.example === state.example));
}

function runPrepare(fresh = false) {
  const t0 = performance.now();
  const { data, checks } = prepare(state.raw, state.map);
  state.data = data;
  state.checks = checks;
  state.view = null;
  state.pinned = new Set();
  state.prepMs = performance.now() - t0;
  if (data) computePositions(data, opts.chrGap);
  computeHighlight();
  $('#empty').hidden = true;
  $('#workspace').hidden = false;
  document.body.classList.add('has-data');
  $('#sec-cols').setAttribute('aria-disabled', 'false');
  updateRegionSelect();
  updateStats();
  updateChecks();
  updateAlert();
  updateLabelColSelect();
  updateHits();
  updateRCode();
  updateSoiButton();
  updatePinsButton();
  if (fresh && !data) setTab('checks');
  else if (fresh && state.tab === 'checks' && data) setTab('manhattan');
  markDirty();
}

// ---------------------------------------------------------------------------
// Progress, toasts, file info

function showProgress(label) {
  const p = $('#progress');
  p.hidden = false;
  $('.progress-text', p).textContent = label;
  const bar = $('.progress-bar', p);
  bar.classList.add('indeterminate');
  $('i', bar).style.width = '0';
}
function setProgress(frac, rows, text) {
  const p = $('#progress');
  const bar = $('.progress-bar', p);
  if (frac > 0) { bar.classList.remove('indeterminate'); $('i', bar).style.width = `${Math.min(100, frac * 100).toFixed(1)}%`; }
  $('.progress-text', p).textContent = text || `${fmtInt(rows)} rows read`;
}
function hideProgress() { $('#progress').hidden = true; }

let toastTimer = 0;
function toast(msg, error = false) {
  const t = $('#toast');
  t.textContent = msg;
  t.classList.toggle('error', error);
  t.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { t.hidden = true; }, error ? 6000 : 2600);
}

function fmtBytes(b) {
  if (!b) return '';
  if (b < 1024) return `${b} B`;
  if (b < 1024 ** 2) return `${(b / 1024).toFixed(0)} KB`;
  if (b < 1024 ** 3) return `${(b / 1024 ** 2).toFixed(1)} MB`;
  return `${(b / 1024 ** 3).toFixed(2)} GB`;
}

function showFileInfo() {
  const r = state.raw;
  const el = $('#file-info');
  el.hidden = false;
  const size = r.size && !state.example ? ` · ${fmtBytes(r.size)}` : '';
  el.innerHTML = `<div class="fi-main"><div class="fi-name"></div><div class="fi-meta">${fmtInt(r.nrows)} rows · ${r.columns.length} columns${size} · ${((r.ms || 0) + (state.prepMs || 0)).toFixed(0)} ms</div></div><button class="btn btn-sm btn-ghost" type="button" id="clear-data">Clear</button>`;
  $('.fi-name', el).textContent = r.name;
  $('#clear-data').onclick = clearData;
}

function clearData() {
  state.raw = state.data = null;
  state.example = null;
  $('#workspace').hidden = true;
  $('#empty').hidden = false;
  document.body.classList.remove('has-data');
  $('#file-info').hidden = true;
  $('#sec-cols').setAttribute('aria-disabled', 'true');
  $$('#sec-cols select[data-map]').forEach((s) => { if (s.dataset.map !== 'scale') s.innerHTML = ''; });
  $$('.example').forEach((b) => b.classList.remove('active'));
  $('#file-input').value = '';
  updateSoiButton();
}

// ---------------------------------------------------------------------------
// Columns panel

function fillColumnSelects() {
  const cols = state.raw.columns;
  for (const sel of $$('#sec-cols select[data-map]')) {
    const key = sel.dataset.map;
    if (key === 'scale') { sel.value = state.map.scale; continue; }
    sel.innerHTML = '';
    const none = document.createElement('option');
    none.value = '-1';
    none.textContent = key === 'snp' ? 'None' : 'Choose a column…';
    sel.appendChild(none);
    cols.forEach((c, j) => {
      const o = document.createElement('option');
      o.value = String(j);
      o.textContent = c.name;
      sel.appendChild(o);
    });
    sel.value = String(state.map[key]);
  }
}

function updateLabelColSelect() {
  const sel = $('#opt-labelCol');
  sel.innerHTML = '';
  const auto = document.createElement('option');
  auto.value = '';
  auto.textContent = state.map && state.map.snp >= 0 ? `SNP column (${state.raw.columns[state.map.snp].name})` : 'chr:position';
  sel.appendChild(auto);
  if (state.raw) {
    state.raw.columns.forEach((c) => {
      const o = document.createElement('option');
      o.value = c.name;
      o.textContent = c.name;
      sel.appendChild(o);
    });
  }
  sel.value = state.raw && state.raw.columns.some((c) => c.name === opts.labelCol) ? opts.labelCol : '';
}

$$('#sec-cols select[data-map]').forEach((sel) => {
  sel.addEventListener('change', () => {
    if (!state.raw) return;
    const key = sel.dataset.map;
    state.map[key] = key === 'scale' ? sel.value : parseInt(sel.value, 10);
    runPrepare();
    updatePlaceholders();
    updateStats();
  });
});

// ---------------------------------------------------------------------------
// Options panel

function syncOptionInputs() {
  for (const el of $$('[data-opt]')) {
    const v = getPath(opts, el.dataset.opt);
    if (el.type === 'checkbox') el.checked = !!v;
    else el.value = v ?? '';
  }
  $('#palette').value = opts.palette;
  $('#export-scale').value = String(opts.exportScale);
  $('#export-scale-qq').value = String(opts.exportScale);
  renderSwatches();
  updateOutputs();
  updatePlaceholders();
}

function updateOutputs() {
  for (const out of $$('output[data-for]')) {
    const k = out.dataset.for;
    const v = getPath(opts, k);
    out.textContent = k === 'chrGap' ? `${(v * 100).toFixed(1)}%` : k === 'opacity' ? `${Math.round(v * 100)}%` : `${v}px`;
  }
}

function updatePlaceholders() {
  const L = state.mhL;
  $('#opt-xlab').placeholder = L ? L.autoXlab : 'Chromosome';
  $('#opt-ylab').placeholder = plainText(autoYlab()).replace('log10', 'log₁₀');
  const raw = state.map && state.map.scale === 'raw';
  $('#line-hint').textContent = raw ? 'Raw-value mode: thresholds are y-axis values. Leave blank for no line.' : 'Thresholds are p-values. Lines are drawn at −log10(p).';
  // threshold inputs edit the p-value or the raw-value setting depending on the scale
  for (const [id, key] of [['#thr-sugg', 'suggestive'], ['#thr-gw', 'genomewide']]) {
    const el = $(id);
    const path = `${key}.${raw ? 'raw' : 'p'}`;
    if (el.dataset.opt !== path) { el.dataset.opt = path; el.value = getPath(opts, path); }
    el.placeholder = raw ? 'none' : '';
  }
  const an = $('#thr-annot');
  const apath = raw ? 'annotateRaw' : 'annotateP';
  if (an.dataset.opt !== apath) { an.dataset.opt = apath; an.value = opts[apath]; }
  an.placeholder = raw ? 'any' : '';
  $('#annot-label').textContent = raw ? 'Label values ≥' : 'Label SNPs with p <';
}

function autoYlab() {
  if (state.map && state.map.scale === 'raw' && state.raw) return state.raw.columns[state.map.p].name;
  return '−log_{10}(*p*)';
}

const NUMERIC = new Set(['pointSize', 'opacity', 'chrGap', 'lineWidth', 'highlightScale', 'labelSize', 'fontSize', 'width', 'height', 'margin.t', 'margin.r', 'margin.b', 'margin.l', 'qq.width', 'qq.height']);
const LIMITS = { width: [200, 6000], height: [150, 4000], 'qq.width': [150, 4000], 'qq.height': [150, 4000], fontSize: [4, 48], labelSize: [4, 40], pointSize: [0.5, 30], lineWidth: [0.1, 10], highlightScale: [0.25, 6] };

let hlTimer = 0;
function onOptInput(el, commit) {
  const path = el.dataset.opt;
  let v;
  if (el.type === 'checkbox') v = el.checked;
  else if (NUMERIC.has(path)) {
    v = parseFloat(el.value);
    if (!Number.isFinite(v)) return;
    const lim = LIMITS[path];
    if (lim) {
      if (!commit && (v < lim[0] || v > lim[1])) return; // wait for blur before clamping
      v = Math.min(lim[1], Math.max(lim[0], v));
      if (commit) el.value = v;
    }
    if (path.startsWith('margin.')) v = Math.max(0, Math.min(1000, v));
  } else v = el.value;
  if (getPath(opts, path) === v) return;
  setPath(opts, path, v);
  saveOpts();
  updateOutputs();
  if (path === 'highlightText') {
    clearTimeout(hlTimer);
    hlTimer = setTimeout(() => { computeHighlight(); updateRCode(); markDirty(); }, 200);
    return;
  }
  if (path === 'chrGap' && state.data) {
    computePositions(state.data, opts.chrGap);
    state.view = null;
    updateRegionSelect();
  }
  if (path === 'plinkNames' || path === 'chrLabels') updateRegionSelect();
  if (path.startsWith('genomewide') || path.startsWith('suggestive')) updateStats();
  if (path === 'transparent') updateTransparentBg();
  updateRCode();
  markDirty();
}

for (const el of $$('[data-opt]')) {
  const evt = el.type === 'checkbox' || el.tagName === 'SELECT' || el.type === 'color' ? 'change' : 'input';
  el.addEventListener(evt, () => onOptInput(el, false));
  if (el.type === 'color') el.addEventListener('input', () => onOptInput(el, false));
  if (el.type === 'number') el.addEventListener('change', () => onOptInput(el, true));
}

$('#palette').addEventListener('change', (e) => {
  const v = e.target.value;
  opts.palette = v;
  if (PALETTES[v]) opts.colors = PALETTES[v].slice();
  renderSwatches();
  saveOpts();
  updateRCode();
  markDirty();
});

function renderSwatches() {
  const box = $('#swatches');
  box.innerHTML = '';
  opts.colors.forEach((c, i) => {
    const w = document.createElement('span');
    w.className = 'swatch';
    const inp = document.createElement('input');
    inp.type = 'color';
    inp.value = c;
    inp.setAttribute('aria-label', `Color ${i + 1}`);
    inp.addEventListener('input', () => {
      opts.colors[i] = inp.value;
      opts.palette = 'custom';
      $('#palette').value = 'custom';
      saveOpts(); updateRCode(); markDirty();
    });
    w.appendChild(inp);
    if (opts.colors.length > 1) {
      const rm = document.createElement('button');
      rm.type = 'button';
      rm.className = 'rm';
      rm.textContent = '×';
      rm.setAttribute('aria-label', `Remove color ${i + 1}`);
      rm.onclick = () => {
        opts.colors.splice(i, 1);
        opts.palette = 'custom';
        $('#palette').value = 'custom';
        renderSwatches(); saveOpts(); updateRCode(); markDirty();
      };
      w.appendChild(rm);
    }
    box.appendChild(w);
  });
  if (opts.colors.length < 30) {
    const add = document.createElement('button');
    add.type = 'button';
    add.className = 'swatch-add';
    add.textContent = '+';
    add.title = 'Add a color';
    add.setAttribute('aria-label', 'Add a color');
    add.onclick = () => {
      opts.colors.push(hsl(Math.random() * 360, 60, 45));
      opts.palette = 'custom';
      $('#palette').value = 'custom';
      renderSwatches(); saveOpts(); updateRCode(); markDirty();
    };
    box.appendChild(add);
  }
}

$('#reset-opts').addEventListener('click', () => {
  const hl = opts.highlightText;
  opts = clone(DEFAULTS);
  opts.highlightText = hl;
  saveOpts();
  syncOptionInputs();
  if (state.data) { computePositions(state.data, opts.chrGap); state.view = null; }
  updateRegionSelect();
  updateLabelColSelect();
  updateStats();
  updateRCode();
  updateTransparentBg();
  markDirty();
  toast('Options reset to defaults');
});

for (const sel of [$('#export-scale'), $('#export-scale-qq')]) {
  sel.addEventListener('change', () => {
    opts.exportScale = parseInt(sel.value, 10);
    $('#export-scale').value = $('#export-scale-qq').value = sel.value;
    saveOpts();
  });
}

function updateTransparentBg() {
  for (const w of [$('#mh-wrap'), $('#qq-wrap')]) w.classList.toggle('transparent-bg', !!opts.transparent);
}

// ---------------------------------------------------------------------------
// Highlights and labels

function computeHighlight() {
  const status = $('#highlight-status');
  const names = opts.highlightText.split(/[\s,;]+/).filter(Boolean);
  state.highlightNames = names;
  state.highlight = [];
  if (!names.length) { status.textContent = ''; return; }
  const D = state.data;
  if (!D) { status.textContent = ''; return; }
  if (D.snpCol < 0) { status.textContent = 'Choose a SNP column to highlight by name.'; return; }
  const want = new Set(names);
  const found = new Set();
  for (let i = 0; i < D.n; i++) {
    const s = snpName(D, i);
    if (want.has(s)) { state.highlight.push(i); found.add(s); }
  }
  const miss = want.size - found.size;
  status.textContent = miss
    ? `${fmtInt(found.size)} of ${fmtInt(want.size)} found. Not in the data: ${[...want].filter((s) => !found.has(s)).slice(0, 3).join(', ')}${miss > 3 ? ', …' : ''}`
    : `All ${fmtInt(want.size)} found.`;
}

function labelColumnIndex() {
  if (!state.raw) return -1;
  if (opts.labelCol) {
    const j = state.raw.columns.findIndex((c) => c.name === opts.labelCol);
    if (j >= 0) return j;
  }
  return state.map.snp;
}

function labelText(i) {
  const D = state.data;
  const j = labelColumnIndex();
  if (j >= 0) {
    const col = state.raw.columns[j];
    const o = D.orig[i];
    const t = cellText(col, o);
    if (t !== null && t !== '' && t !== '.') return t;
    if (t === null) return fmtNum(col.num[o]);
  }
  return `${chromLabel(D.chroms[D.chrIdx[i]], opts)}:${fmtInt(D.bp[i])}`;
}

function computeLabels() {
  const D = state.data;
  const out = [];
  const seen = new Set();
  if (opts.annotate) {
    const spec = D.scale === 'raw' ? opts.annotateRaw : opts.annotateP;
    const thr = D.scale === 'raw' && spec === '' ? -Infinity : thresholdY(spec, D.scale);
    if (thr !== null) {
      let cands = [];
      if (D.top.length && D.y[D.top[D.top.length - 1]] < thr) {
        for (const i of D.top) { if (D.y[i] < thr) break; cands.push(i); }
      } else {
        for (let i = 0; i < D.n; i++) if (D.y[i] >= thr) cands.push(i);
        cands.sort((a, b) => D.y[b] - D.y[a]);
      }
      if (opts.annotateTop) {
        const per = new Set();
        cands = cands.filter((i) => { const c = D.chrIdx[i]; if (per.has(c)) return false; per.add(c); return true; });
      }
      state.labelsCapped = cands.length > 400;
      for (const i of cands.slice(0, 400)) { out.push({ i, text: labelText(i), pinned: false }); seen.add(i); }
    }
  }
  for (const i of state.pinned) if (!seen.has(i)) out.push({ i, text: labelText(i), pinned: true });
  return out;
}

function updateSoiButton() {
  $('#use-soi').hidden = !(state.example === 'gwasResults' && state.soi && !opts.highlightText.trim());
}
$('#use-soi').addEventListener('click', () => {
  if (!state.soi) return;
  opts.highlightText = state.soi.join(' ');
  $('#opt-highlight').value = opts.highlightText;
  computeHighlight();
  updateSoiButton();
  updateRCode();
  markDirty();
});
$('#opt-highlight').addEventListener('input', updateSoiButton);

function updatePinsButton() { $('#clear-pins').hidden = state.pinned.size === 0; }
$('#clear-pins').addEventListener('click', () => { state.pinned.clear(); updatePinsButton(); updateHits(); markDirty(); });

// ---------------------------------------------------------------------------
// Rendering

let rafPending = false;
function markDirty() {
  state.dirty.manhattan = state.dirty.qq = true;
  if (rafPending) return;
  rafPending = true;
  requestAnimationFrame(() => { rafPending = false; renderVisible(); });
}

function renderVisible() {
  if (!state.data) {
    for (const id of ['mh', 'qq']) { const c = $(`#${id}`); c.width = c.height = 0; }
    $('#mh-foot').textContent = '';
    return;
  }
  if (state.tab === 'manhattan' && state.dirty.manhattan) drawManhattan();
  if (state.tab === 'qq' && state.dirty.qq) drawQQ();
}

function canvasRes(W, wrap) {
  const dpr = window.devicePixelRatio || 1;
  const avail = (wrap.parentElement.clientWidth || W) - 0;
  const shown = Math.min(W, avail);
  return Math.max(0.5, Math.min(dpr, (dpr * shown) / W) * 1.0);
}

function sizeCanvas(canvas, overlay, wrap, W, H) {
  const res = canvasRes(W, wrap);
  const cw = Math.round(W * res), ch = Math.round(H * res);
  if (canvas.width !== cw || canvas.height !== ch) { canvas.width = cw; canvas.height = ch; }
  if (overlay.width !== cw || overlay.height !== ch) { overlay.width = cw; overlay.height = ch; }
  wrap.style.width = `${W}px`;
  canvas.style.width = '100%';
  return cw / W;
}

function drawManhattan() {
  const D = state.data;
  const W = opts.width, H = opts.height;
  const canvas = $('#mh'), over = $('#mh-over'), wrap = $('#mh-wrap');
  const res = sizeCanvas(canvas, over, wrap, W, H);
  const t0 = performance.now();
  const P = new CanvasPainter(canvas.getContext('2d'), res);
  P.fast = true;
  const labels = computeLabels();
  state.mhL = renderManhattan(P, D, opts, state.view, W, H, { highlight: state.highlight, labels, autoYlab: autoYlab() });
  state.mhL.res = res;
  const ms = performance.now() - t0;
  state.dirty.manhattan = false;
  drawOverlay();
  updatePlaceholders();
  updateZoomUI();
  const L = state.mhL;
  const shown = L.hi - L.lo;
  const parts = [];
  if (state.view) parts.push(`Showing ${regionText(L)} · ${fmtInt(shown)} variants`);
  else parts.push(`${fmtInt(D.n)} variants`);
  parts.push(`drawn in ${ms.toFixed(0)} ms`);
  if (L.capped) parts.push(`${fmtInt(L.capped)} above the y-axis max shown as triangles`);
  if (L.labelsSkipped) parts.push(`${L.labelsSkipped} label${L.labelsSkipped > 1 ? 's' : ''} hidden to avoid overlap`);
  if (state.labelsCapped && opts.annotate) parts.push('labels capped at 400; try a stricter threshold');
  const custom = opts.chrLabels.trim() ? opts.chrLabels.split(/[,\s]+/).filter(Boolean).length : 0;
  $('#chrlabels-status').textContent = custom && custom !== D.chroms.length ? `You gave ${custom} labels but the data has ${D.chroms.length} chromosomes, so the defaults are used.` : '';
  $('#mh-foot').textContent = parts.join(' · ');
}

function regionText(L) {
  const D = state.data;
  if (L.single) {
    const c = L.single;
    const a = Math.max(c.minbp, L.d0 - c.base), b = Math.min(c.maxbp, L.d1 - c.base);
    const [div, unit] = bpUnit(b - a);
    const dec = div === 1 ? 0 : 2;
    return `chr${chromLabel(c, opts)}:${(a / div).toFixed(dec)}–${(b / div).toFixed(dec)} ${unit}`;
  }
  const cs = D.chroms.filter((c) => c.end >= L.d0 && c.start <= L.d1);
  if (!cs.length) return 'an empty region';
  return `chromosomes ${chromLabel(cs[0], opts)}–${chromLabel(cs[cs.length - 1], opts)}`;
}

function drawQQ() {
  const D = state.data;
  const off = D.scale === 'raw';
  $('#qq-off').hidden = !off;
  $('#qq-wrap').hidden = off;
  $$('[data-export^="qq"]').forEach((b) => { b.disabled = off; });
  state.dirty.qq = false;
  if (off) { $('#qq-note').textContent = ''; return; }
  const W = opts.qq.width, H = opts.qq.height;
  const canvas = $('#qq'), over = $('#qq-over'), wrap = $('#qq-wrap');
  const res = sizeCanvas(canvas, over, wrap, W, H);
  const P = new CanvasPainter(canvas.getContext('2d'), res);
  P.fast = true;
  state.qqL = renderQQ(P, D, opts, W, H, { highlightRanks: qqHighlightRanks() });
  state.qqL.res = res;
  clearOverlay($('#qq-over'));
  $('#qq-note').textContent = D.lambda !== null ? `λGC = ${D.lambda.toFixed(3)} (median χ² / 0.455)${D.lambda > 1.1 ? '. That is high: check for population structure or relatedness.' : ''}` : '';
}

function qqHighlightRanks() {
  if (!opts.qq.highlight || !state.highlight.length) return null;
  const D = state.data, s = D.sortedY, n = D.n;
  return state.highlight.map((i) => {
    const v = D.y[i];
    let lo = 0, hi = n;
    while (lo < hi) { const m = (lo + hi) >>> 1; if (s[m] <= v) lo = m + 1; else hi = m; }
    return [n - lo, v];
  });
}

window.addEventListener('resize', () => { clearTimeout(window.__rz); window.__rz = setTimeout(markDirty, 120); });

// ---------------------------------------------------------------------------
// Interaction: hover, zoom, pin

function plotCoords(ev, canvas, W) {
  const r = canvas.getBoundingClientRect();
  const k = W / r.width;
  return { x: (ev.clientX - r.left) * k, y: (ev.clientY - r.top) * k, k, rect: r };
}

function lowerBound(arr, v) {
  let lo = 0, hi = arr.length;
  while (lo < hi) { const mid = (lo + hi) >>> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; }
  return lo;
}

function nearestManhattan(x, y) {
  const L = state.mhL, D = state.data;
  if (!L) return null;
  const R = Math.max(6, opts.pointSize / 2 + 4);
  const a = Math.max(L.lo, lowerBound(D.pos, L.ix(x - R)));
  const b = Math.min(L.hi, lowerBound(D.pos, L.ix(x + R)));
  let best = -1, bd = R * R;
  for (let i = a; i < b; i++) {
    const v = D.y[i];
    if (v < L.ymin) continue;
    const dx = L.sx(D.pos[i]) - x, dy = L.syc(v) - y;
    const d = dx * dx + dy * dy;
    if (d < bd || (d === bd && best >= 0 && v > D.y[best])) { bd = d; best = i; }
  }
  return best >= 0 ? best : null;
}

function nearestQQ(x, y) {
  const L = state.qqL, D = state.data;
  if (!L) return null;
  const R = Math.max(6, opts.pointSize / 2 + 4);
  let best = null, bd = R * R;
  const top = D.top;
  for (let j = 0; j < top.length; j++) {
    const i = top[j];
    const e = -Math.log10((j + 1 - L.a) / L.denom);
    const dx = L.sx(e) - x, dy = L.sy(D.y[i]) - y;
    const d = dx * dx + dy * dy;
    if (d < bd) { bd = d; best = { i, e }; }
  }
  return best;
}

function tooltipHTML(i, extraRows = []) {
  const D = state.data, raw = state.raw, m = state.map;
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const name = D.snpCol >= 0 ? snpName(D, i) : '';
  const c = D.chroms[D.chrIdx[i]];
  const rows = [];
  rows.push(['Position', `${esc(chromLabel(c, opts))}:${fmtInt(D.bp[i])}`]);
  if (D.scale === 'raw') rows.push([esc(raw.columns[m.p].name), fmtNum(D.y[i])]);
  else { rows.push(['p', fmtPFromLog(D.y[i])]); rows.push(['−log10 p', D.y[i].toFixed(2)]); }
  rows.push(...extraRows);
  const used = new Set([m.chr, m.bp, m.p, m.snp]);
  let extra = 0;
  raw.columns.forEach((col, j) => {
    if (used.has(j) || extra >= 5) return;
    const o = D.orig[i];
    const t = cellText(col, o);
    rows.push([esc(col.name), t === null ? fmtNum(col.num[o]) : esc(t || 'NA')]);
    extra++;
  });
  const pinned = state.pinned.has(i);
  return `<b>${esc(name || `${chromLabel(c, opts)}:${fmtInt(D.bp[i])}`)}</b><div class="tt-grid">${rows.map(([k, v]) => `<span>${k}</span><span>${v}</span>`).join('')}</div>`
    + (state.tab === 'manhattan' ? `<div class="tt-hint">Click to ${pinned ? 'remove the' : 'add a'} label</div>` : '');
}

function showTip(tip, canvas, W, px, py, html) {
  const r = canvas.getBoundingClientRect();
  const k = r.width / W;
  tip.innerHTML = html;
  tip.hidden = false;
  const x = px * k;
  const half = tip.offsetHeight / 2;
  const y = Math.min(Math.max(py * k, half + 4), r.height - half - 4); // keep it inside the plot
  tip.style.left = `${x}px`;
  tip.style.top = `${y}px`;
  tip.classList.toggle('left', x > r.width * 0.62);
}

function clearOverlay(over) {
  const ctx = over.getContext('2d');
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.clearRect(0, 0, over.width, over.height);
}

function ring(over, res, x, y, r) {
  const ctx = over.getContext('2d');
  ctx.setTransform(res, 0, 0, res, 0, 0);
  ctx.beginPath();
  ctx.arc(x, y, r + 3, 0, Math.PI * 2);
  ctx.lineWidth = 3.5; ctx.strokeStyle = 'rgba(255,255,255,0.95)'; ctx.stroke();
  ctx.lineWidth = 1.8; ctx.strokeStyle = '#0b6bcb'; ctx.stroke();
}

let brush = null;
function drawOverlay() {
  const over = $('#mh-over');
  clearOverlay(over);
  const L = state.mhL;
  if (!L) return;
  const ctx = over.getContext('2d');
  if (brush && brush.active) {
    ctx.setTransform(L.res, 0, 0, L.res, 0, 0);
    const a = Math.max(L.x0, Math.min(brush.x0, brush.x1)), b = Math.min(L.x1, Math.max(brush.x0, brush.x1));
    ctx.fillStyle = 'rgba(11,107,203,0.12)';
    ctx.fillRect(a, L.y0, b - a, L.y1 - L.y0);
    ctx.strokeStyle = 'rgba(11,107,203,0.7)';
    ctx.lineWidth = 1;
    ctx.strokeRect(a + 0.5, L.y0 + 0.5, b - a - 1, L.y1 - L.y0 - 1);
  }
  if (state.hover !== null && state.hover !== undefined) {
    const i = state.hover, D = state.data;
    ring(over, L.res, L.sx(D.pos[i]), L.syc(D.y[i]), opts.pointSize / 2);
  }
}

(function manhattanEvents() {
  const over = $('#mh-over'), tip = $('#mh-tip');
  over.addEventListener('mousemove', (ev) => {
    if (!state.data || !state.mhL) return;
    const { x, y } = plotCoords(ev, over, opts.width);
    if (brush) {
      brush.x1 = x;
      if (Math.abs(brush.x1 - brush.x0) > 4) brush.active = true;
      if (brush.active) { state.hover = null; tip.hidden = true; drawOverlay(); return; }
    }
    const i = nearestManhattan(x, y);
    if (i !== state.hover) {
      state.hover = i;
      drawOverlay();
    }
    if (i === null) { tip.hidden = true; over.style.cursor = 'crosshair'; return; }
    over.style.cursor = 'pointer';
    const L = state.mhL;
    showTip(tip, over, opts.width, L.sx(state.data.pos[i]), L.syc(state.data.y[i]), tooltipHTML(i));
  });
  over.addEventListener('mouseleave', () => {
    state.hover = null; tip.hidden = true;
    if (brush && !brush.active) brush = null;
    drawOverlay();
  });
  over.addEventListener('mousedown', (ev) => {
    if (ev.button !== 0 || !state.mhL) return;
    const { x } = plotCoords(ev, over, opts.width);
    brush = { x0: x, x1: x, active: false };
    ev.preventDefault();
  });
  window.addEventListener('mouseup', (ev) => {
    if (!brush) return;
    const b = brush;
    brush = null;
    const L = state.mhL;
    if (b.active) {
      const a = Math.max(L.x0, Math.min(b.x0, b.x1)), c = Math.min(L.x1, Math.max(b.x0, b.x1));
      const d0 = L.ix(a), d1 = L.ix(c);
      if (d1 - d0 > 50) { state.view = { x0: d0, x1: d1 }; state.hover = null; markDirty(); }
      drawOverlay();
    } else if (ev.target === over && state.hover !== null) {
      const i = state.hover;
      if (state.pinned.has(i)) state.pinned.delete(i); else state.pinned.add(i);
      updatePinsButton();
      updateHits();
      markDirty();
      const L2 = state.mhL;
      setTimeout(() => { if (state.hover === i) showTip(tip, over, opts.width, L2.sx(state.data.pos[i]), L2.syc(state.data.y[i]), tooltipHTML(i)); }, 30);
    }
  });
  over.addEventListener('dblclick', () => { if (state.view) { state.view = null; markDirty(); } });
})();

(function qqEvents() {
  const over = $('#qq-over'), tip = $('#qq-tip');
  let last = null;
  over.addEventListener('mousemove', (ev) => {
    if (!state.data || !state.qqL) return;
    const { x, y } = plotCoords(ev, over, opts.qq.width);
    const hit = nearestQQ(x, y);
    const L = state.qqL;
    if ((hit && hit.i) !== (last && last.i)) {
      clearOverlay(over);
      if (hit) ring(over, L.res, L.sx(hit.e), L.sy(state.data.y[hit.i]), opts.pointSize / 2);
      last = hit;
    }
    if (!hit) { tip.hidden = true; return; }
    showTip(tip, over, opts.qq.width, L.sx(hit.e), L.sy(state.data.y[hit.i]), tooltipHTML(hit.i, [['Expected', hit.e.toFixed(2)]]));
  });
  over.addEventListener('mouseleave', () => { tip.hidden = true; clearOverlay(over); last = null; });
})();

// Region / zoom
function updateRegionSelect() {
  const sel = $('#region');
  const D = state.data;
  sel.innerHTML = '';
  if (!D) return;
  const add = (v, t) => { const o = document.createElement('option'); o.value = v; o.textContent = t; sel.appendChild(o); };
  add('all', D.chroms.length > 1 ? 'Whole genome' : `Chromosome ${chromLabel(D.chroms[0], opts)}`);
  if (D.chroms.length > 1) D.chroms.forEach((c, k) => add(String(k), `Chromosome ${chromLabel(c, opts)}`));
  add('custom', 'Zoomed region');
  updateZoomUI();
}

function updateZoomUI() {
  const sel = $('#region');
  const D = state.data;
  if (!D) return;
  let v = 'all';
  if (state.view) {
    v = 'custom';
    if (state.view.chrom !== undefined) v = String(state.view.chrom);
  }
  $('option[value="custom"]', sel).hidden = v !== 'custom';
  sel.value = v;
  $('#zoom-reset').hidden = !state.view;
}

$('#region').addEventListener('change', (e) => {
  const D = state.data;
  const v = e.target.value;
  if (v === 'all') state.view = null;
  else if (v !== 'custom') {
    const c = D.chroms[+v];
    const pad = (c.end - c.start) * 0.02 || 1;
    state.view = { x0: c.start - pad, x1: c.end + pad, chrom: +v };
  }
  state.hover = null;
  markDirty();
});
$('#zoom-reset').addEventListener('click', () => { state.view = null; markDirty(); });

function zoomTo(i) {
  const D = state.data;
  const c = D.chroms[D.chrIdx[i]];
  const half = Math.max(1e6, (c.maxbp - c.minbp) * 0.004);
  const p = D.pos[i];
  state.view = { x0: Math.max(c.start, p - half), x1: Math.min(c.end, p + half) };
  if (state.view.x1 - state.view.x0 < 1000) state.view = { x0: p - 5e5, x1: p + 5e5 };
}

// ---------------------------------------------------------------------------
// Tabs

function setTab(tab) {
  state.tab = tab;
  for (const t of $$('.tab')) t.setAttribute('aria-selected', String(t.dataset.tab === tab));
  for (const p of $$('.panel')) p.hidden = p.dataset.panel !== tab;
  renderVisible();
}
$$('.tab').forEach((t) => t.addEventListener('click', () => setTab(t.dataset.tab)));

// ---------------------------------------------------------------------------
// Stats, alert, checks

function updateStats() {
  const el = $('#stats');
  const D = state.data;
  if (!D) { el.innerHTML = ''; return; }
  const tiles = [];
  tiles.push(['Variants', fmtInt(D.n)]);
  tiles.push(['Chromosomes', String(D.chroms.length)]);
  if (D.scale === 'raw') tiles.push(['Max value', fmtNum(D.maxY)]);
  else tiles.push(['Smallest p', fmtPFromLog(D.maxY)]);
  if (D.lambda !== null) tiles.push(['λGC', D.lambda.toFixed(3), D.lambda > 1.1 ? 'warn' : '']);
  const spec = lineSpec(opts, 'genomewide', D.scale);
  const gw = thresholdY(spec, D.scale);
  if (gw !== null) {
    let k = 0;
    for (let i = 0; i < D.n; i++) if (D.y[i] >= gw) k++;
    const label = D.scale === 'raw' ? `≥ ${spec}` : `p < ${spec}`;
    tiles.push([label, fmtInt(k), '', `${k === 1 ? 'variant' : 'variants'}`]);
  }
  el.innerHTML = tiles.map(([k, v, cls, small]) => `<div class="stat"><div class="k">${k}</div><div class="v ${cls || ''}">${v}${small ? `<small>${small}</small>` : ''}</div></div>`).join('');
}

function updateAlert() {
  const el = $('#alert');
  const errs = state.checks.filter((c) => c.level === 'error');
  const warns = state.checks.filter((c) => c.level === 'warn');
  const esc = (s) => s.replace(/[&<>]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
  if (errs.length) {
    el.className = 'alert error';
    el.innerHTML = `<div class="alert-body"><b>Can't plot this yet.</b><ul>${errs.map((e) => `<li><b>${esc(e.title)}.</b> ${esc(e.detail)}</li>`).join('')}</ul></div>`;
    el.hidden = false;
  } else if (warns.length) {
    el.className = 'alert warn';
    el.innerHTML = `<div class="alert-body"><b>${warns.length} warning${warns.length > 1 ? 's' : ''} about your data.</b> ${esc(warns[0].title)}${warns.length > 1 ? ', and more' : ''}. <a href="#" id="see-checks">See the data checks</a></div>`;
    el.hidden = false;
    $('#see-checks').onclick = (e) => { e.preventDefault(); setTab('checks'); };
  } else el.hidden = true;
  const badge = $('#checks-badge');
  const n = errs.length || warns.length;
  badge.hidden = !n;
  badge.textContent = String(n);
  badge.className = `badge ${errs.length ? 'error' : 'warn'}`;
  const hasData = !!state.data;
  for (const t of $$('.tab')) if (t.dataset.tab !== 'checks') t.disabled = !hasData;
  if (!hasData && state.tab !== 'checks') setTab('checks');
}

function updateChecks() {
  const icon = { ok: '✓', warn: '!', error: '×', info: 'i' };
  // errors, then the row-count summary, then warnings, other passes and notes
  const order = { error: 0, warn: 2, ok: 3, info: 4 };
  const summary = state.checks.find((c) => c.level === 'ok');
  const rank = (c) => (c === summary ? 1 : order[c.level]);
  const list = state.checks.slice().sort((a, b) => rank(a) - rank(b));
  const ul = $('#checks');
  ul.innerHTML = '';
  for (const c of list) {
    const li = document.createElement('li');
    li.className = `check-item ${c.level}`;
    li.innerHTML = `<span class="ic" aria-hidden="true">${icon[c.level]}</span><div><div class="t"></div><div class="d"></div></div>`;
    $('.t', li).textContent = c.title;
    $('.d', li).textContent = c.detail;
    ul.appendChild(li);
  }
  updatePreview();
}

function updatePreview() {
  const raw = state.raw;
  const tbl = $('#preview-table');
  if (!raw) { tbl.innerHTML = ''; return; }
  const m = state.map;
  const tags = new Map([[m.chr, 'CHR'], [m.bp, 'BP'], [m.p, 'P'], [m.snp, 'SNP']]);
  const n = Math.min(raw.nrows, 50);
  const delim = { '\t': 'tab', ',': 'comma', ';': 'semicolon', ws: 'whitespace', sim: 'simulated' }[raw.delimiter] || raw.delimiter;
  $('#preview-note').textContent = `First ${n} of ${fmtInt(raw.nrows)} rows as read (${delim}-delimited${raw.hasHeader ? ', with header' : ', no header'}${raw.comments ? `, ${fmtInt(raw.comments)} comment lines skipped` : ''}).`;
  const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
  let h = '<thead><tr>';
  raw.columns.forEach((c, j) => { h += `<th>${esc(c.name)}${tags.has(j) ? `<span class="map-tag">${tags.get(j)}</span>` : ''}</th>`; });
  h += '</tr></thead><tbody>';
  for (let i = 0; i < n; i++) {
    h += '<tr>';
    for (const c of raw.columns) {
      const t = cellText(c, i);
      if (t === null) h += `<td class="num">${esc(fmtNum(c.num[i], 6))}</td>`;
      else h += `<td${t === '' || t === 'NA' ? ' class="na"' : ''}>${esc(t === '' ? 'NA' : t)}</td>`;
    }
    h += '</tr>';
  }
  tbl.innerHTML = h + '</tbody>';
}

// ---------------------------------------------------------------------------
// Top hits

const HITS_SHOWN = 200;
function hitRows() {
  const D = state.data;
  return D.top.slice(0, HITS_SHOWN);
}

function updateHits() {
  const D = state.data;
  const tbl = $('#hits-table');
  if (!D) { tbl.innerHTML = ''; return; }
  const raw = state.raw, m = state.map;
  const used = new Set([m.chr, m.bp, m.p, m.snp]);
  const extras = raw.columns.map((c, j) => j).filter((j) => !used.has(j)).slice(0, 8);
  const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
  const rows = hitRows();
  const raws = D.scale === 'raw';
  let h = `<thead><tr><th class="num">#</th>${D.snpCol >= 0 ? '<th>SNP</th>' : ''}<th>CHR</th><th class="num">BP</th>${raws ? `<th class="num">${esc(raw.columns[m.p].name)}</th>` : '<th class="num">P</th><th class="num">−log10 P</th>'}`;
  for (const j of extras) h += `<th>${esc(raw.columns[j].name)}</th>`;
  h += '</tr></thead><tbody>';
  rows.forEach((i, k) => {
    const o = D.orig[i];
    h += `<tr class="clickable${state.pinned.has(i) ? ' pinned' : ''}" data-i="${i}"><td class="num">${k + 1}</td>`;
    if (D.snpCol >= 0) h += `<td>${esc(snpName(D, i))}</td>`;
    h += `<td>${esc(chromLabel(D.chroms[D.chrIdx[i]], opts))}</td><td class="num">${fmtInt(D.bp[i])}</td>`;
    h += raws ? `<td class="num">${fmtNum(D.y[i])}</td>` : `<td class="num">${fmtPFromLog(D.y[i])}</td><td class="num">${D.y[i].toFixed(2)}</td>`;
    for (const j of extras) {
      const c = raw.columns[j];
      const t = cellText(c, o);
      h += t === null ? `<td class="num">${esc(fmtNum(c.num[o]))}</td>` : `<td>${esc(t || 'NA')}</td>`;
    }
    h += '</tr>';
  });
  tbl.innerHTML = h + '</tbody>';
  $('#hits-note').textContent = `Top ${Math.min(HITS_SHOWN, D.n)} variants by ${raws ? 'value' : 'p-value'}. Click a row to zoom to it on the Manhattan plot and label it.`;
}

$('#hits-table').addEventListener('click', (e) => {
  const tr = e.target.closest('tr[data-i]');
  if (!tr) return;
  const i = +tr.dataset.i;
  state.pinned.add(i);
  updatePinsButton();
  zoomTo(i);
  setTab('manhattan');
  updateHits();
  markDirty();
});

$('#hits-csv').addEventListener('click', () => {
  const D = state.data;
  if (!D) return;
  const raw = state.raw;
  const quote = (s) => (/[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  const lines = [raw.columns.map((c) => quote(c.name)).join(',')];
  for (const i of hitRows()) {
    const o = D.orig[i];
    lines.push(raw.columns.map((c, j) => {
      if (j === state.map.p && D.scale === 'p') return fmtPFromLog(D.y[i]);
      const t = cellText(c, o);
      return quote(t === null ? String(c.num[o]) : t);
    }).join(','));
  }
  download(new Blob([lines.join('\n') + '\n'], { type: 'text/csv' }), `${baseName()}_top_hits.csv`);
});

// ---------------------------------------------------------------------------
// R code

function updateRCode() {
  const code = rCode(state, opts, state.soi);
  $('#rcode').innerHTML = highlightR(code);
  $('#rcode').dataset.raw = code;
}
$('#copy-r').addEventListener('click', async () => {
  try {
    await navigator.clipboard.writeText($('#rcode').dataset.raw || '');
    toast('R code copied');
  } catch (e) {
    toast('Copy failed. Select the code and copy it manually.', true);
  }
});

// ---------------------------------------------------------------------------
// Export

function baseName() {
  const n = (state.raw && state.raw.name) || 'qqman';
  return n.replace(/\.(gz)$/i, '').replace(/\.[a-z0-9]+$/i, '').replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '') || 'qqman';
}

function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}

function exportPlot(which, fmt) {
  const D = state.data;
  if (!D) return;
  const W = which === 'qq' ? opts.qq.width : opts.width;
  const H = which === 'qq' ? opts.qq.height : opts.height;
  const extra = which === 'qq'
    ? { highlightRanks: qqHighlightRanks() }
    : { highlight: state.highlight, labels: computeLabels(), autoYlab: autoYlab() };
  const draw = (P) => (which === 'qq' ? renderQQ(P, D, opts, W, H, extra) : renderManhattan(P, D, opts, state.view, W, H, extra));
  const name = `${baseName()}_${which === 'qq' ? 'qq' : 'manhattan'}`;
  if (fmt === 'svg') {
    const P = new SvgPainter(W, H);
    const L = draw(P);
    download(new Blob([P.toString()], { type: 'image/svg+xml' }), `${name}.svg`);
    if (L.drawn > 250000) toast(`SVG saved with ${fmtInt(L.drawn)} points. Some editors are slow with that many.`);
    return;
  }
  const scale = opts.exportScale || 2;
  const c = document.createElement('canvas');
  c.width = Math.round(W * scale);
  c.height = Math.round(H * scale);
  if (c.width * c.height > 268e6) { toast('That image is too large to export. Lower the size or resolution.', true); return; }
  const P = new CanvasPainter(c.getContext('2d'), scale);
  draw(P);
  c.toBlob((blob) => {
    if (!blob) { toast('PNG export failed. Try a lower resolution.', true); return; }
    download(blob, `${name}${scale > 1 ? `@${scale}x` : ''}.png`);
  }, 'image/png');
}
$$('[data-export]').forEach((b) => b.addEventListener('click', () => {
  const [which, fmt] = b.dataset.export.split(':');
  exportPlot(which, fmt);
}));

// ---------------------------------------------------------------------------
// Data input wiring

$('#file-input').addEventListener('change', (e) => loadFile(e.target.files[0]));
$('#hero-browse').addEventListener('click', () => $('#file-input').click());
$('#dropzone').addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); $('#file-input').click(); } });
$$('[data-example]').forEach((b) => b.addEventListener('click', () => loadExample(b.dataset.example)));

$('#paste-btn').addEventListener('click', () => {
  const box = $('#paste-box');
  box.hidden = !box.hidden;
  $('#paste-btn').setAttribute('aria-expanded', String(!box.hidden));
  $('#paste-btn').textContent = box.hidden ? 'Paste data instead' : 'Hide paste box';
  if (!box.hidden) $('#paste-text').focus();
});
$('#paste-go').addEventListener('click', () => loadText($('#paste-text').value));
$('#paste-text').addEventListener('keydown', (e) => { if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) loadText($('#paste-text').value); });

// Paste anywhere on the page (outside form fields) to load data.
document.addEventListener('paste', (e) => {
  const t = e.target;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.isContentEditable)) return;
  const text = e.clipboardData && e.clipboardData.getData('text');
  if (text && text.split('\n').length > 2) { e.preventDefault(); $('#paste-text').value = text.length < 2e6 ? text : ''; loadText(text); }
});

// Drag and drop anywhere.
let dragDepth = 0;
const hasFiles = (e) => e.dataTransfer && Array.from(e.dataTransfer.types || []).includes('Files');
window.addEventListener('dragenter', (e) => { if (!hasFiles(e)) return; dragDepth++; $('#drop-overlay').hidden = false; $('#dropzone').classList.add('over'); });
window.addEventListener('dragleave', () => { dragDepth = Math.max(0, dragDepth - 1); if (!dragDepth) { $('#drop-overlay').hidden = true; $('#dropzone').classList.remove('over'); } });
window.addEventListener('dragover', (e) => { if (hasFiles(e)) e.preventDefault(); });
window.addEventListener('drop', (e) => {
  if (!hasFiles(e)) return;
  e.preventDefault();
  dragDepth = 0;
  $('#drop-overlay').hidden = true;
  $('#dropzone').classList.remove('over');
  const f = e.dataTransfer.files[0];
  if (f) loadFile(f);
});

document.addEventListener('keydown', (e) => {
  if (e.key === 'Escape' && state.view && state.tab === 'manhattan' && !/INPUT|TEXTAREA|SELECT/.test(document.activeElement.tagName)) { state.view = null; markDirty(); }
});

// Theme
$('#theme-toggle').addEventListener('click', () => {
  const root = document.documentElement;
  const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
  root.dataset.theme = dark ? 'light' : 'dark';
  try { localStorage.setItem('qqman-theme', root.dataset.theme); } catch (err) { /* ignore */ }
});

// ---------------------------------------------------------------------------
// Decorative hero plot

function drawHero() {
  const svg = $('#hero-plot');
  if (!svg) return;
  let seed = 11;
  const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
  const parts = [];
  const chroms = 22, base = 204, x0 = 14, x1 = 586, line = base - 118;
  const peaks = { 2: 0.75, 5: 1.0, 9: 0.5, 15: 0.62, 18: 0.9 };
  const widths = Array.from({ length: chroms }, (_, k) => 1.3 - k * 0.035);
  const total = widths.reduce((a, b) => a + b, 0);
  let x = x0;
  for (let k = 0; k < chroms; k++) {
    const w = ((x1 - x0) * widths[k]) / total;
    const fill = k % 2 ? 'var(--text-3)' : 'var(--text)';
    const pts = [];
    for (let i = 0, n = Math.round(w * 2.2); i < n; i++) pts.push([x + rnd() * w, -Math.log10(rnd()) * 11]);
    if (peaks[k]) {
      const cx = x + w * (0.3 + rnd() * 0.4);
      for (let i = 0; i < 26; i++) {
        const d = (rnd() - 0.5) * w * 0.12;
        const h = 165 * peaks[k] * Math.exp(-Math.abs(d) / (w * 0.025)) * (0.35 + 0.65 * rnd());
        pts.push([cx + d, h]);
      }
    }
    for (const [px, h] of pts) {
      const hot = base - h < line;
      parts.push(`<circle cx="${px.toFixed(1)}" cy="${Math.max(8, base - h).toFixed(1)}" r="2.1" fill="${hot ? 'var(--accent)' : fill}" opacity="${hot ? 1 : 0.8}"/>`);
    }
    x += w;
  }
  parts.push(`<line x1="${x0}" x2="${x1}" y1="${line}" y2="${line}" stroke="#e5484d" stroke-width="1.2" stroke-dasharray="5 4"/>`);
  svg.innerHTML = parts.join('');
}

// ---------------------------------------------------------------------------
// Init

syncOptionInputs();
updateTransparentBg();
drawHero();
if (location.hash === '#example') loadExample('gwasResults');
