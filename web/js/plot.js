// Manhattan and Q-Q rendering against a small painter interface, so the same
// code draws to the screen (canvas), PNG export (canvas) and SVG export.
import { niceTicks, fmtTick, niceStep, qOrderStat } from './stats.js';
import { chromLabel } from './prepare.js';

export const FONTS = {
  sans: 'Helvetica, Arial, "Liberation Sans", sans-serif',
  serif: '"Times New Roman", Times, "Liberation Serif", serif',
  mono: 'Menlo, Consolas, "Liberation Mono", monospace',
};

const INK = '#1f1f1f';
const INK_SOFT = '#555555';

// ---------------------------------------------------------------------------
// Rich text: _{sub}, ^{sup}, *italic*

export function parseRich(s) {
  const segs = [];
  let cur = '', italic = false;
  const flush = () => { if (cur) segs.push({ t: cur, it: italic, sh: 0 }); cur = ''; };
  for (let i = 0; i < s.length;) {
    const ch = s[i];
    if ((ch === '_' || ch === '^') && s[i + 1] === '{') {
      const j = s.indexOf('}', i + 2);
      if (j > 0) { flush(); segs.push({ t: s.slice(i + 2, j), it: italic, sh: ch === '_' ? 1 : -1 }); i = j + 1; continue; }
    }
    if (ch === '*') { flush(); italic = !italic; i++; continue; }
    if (ch === '\\' && i + 1 < s.length) { cur += s[i + 1]; i += 2; continue; }
    cur += ch; i++;
  }
  flush();
  return segs;
}

export function plainText(s) {
  return parseRich(s).map((g) => g.t).join('');
}

let mctx = null;
function measureCtx() {
  if (!mctx) {
    const c = typeof OffscreenCanvas !== 'undefined' ? new OffscreenCanvas(4, 4) : document.createElement('canvas');
    mctx = c.getContext('2d');
  }
  return mctx;
}
const fontStr = (size, family, weight = 400, italic = false) => `${italic ? 'italic ' : ''}${weight} ${size}px ${family}`;

function segSize(g, size) { return g.sh ? size * 0.72 : size; }

export function measureText(text, size, family, weight = 400) {
  const segs = typeof text === 'string' ? [{ t: text, it: false, sh: 0 }] : text;
  const ctx = measureCtx();
  let w = 0;
  for (const g of segs) {
    ctx.font = fontStr(segSize(g, size), family, weight, g.it);
    w += ctx.measureText(g.t).width;
  }
  return w;
}

// ---------------------------------------------------------------------------
// Painters

const spriteCache = new Map();
function sprite(color, r, res) {
  const key = `${color}|${r}|${res}`;
  let s = spriteCache.get(key);
  if (s) return s;
  const R = r * res;
  const size = Math.ceil(R * 2 + 2);
  const c = typeof OffscreenCanvas !== 'undefined' ? new OffscreenCanvas(size, size) : Object.assign(document.createElement('canvas'), { width: size, height: size });
  const ctx = c.getContext('2d');
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(size / 2, size / 2, R, 0, Math.PI * 2);
  ctx.fill();
  s = { c, half: size / 2 };
  if (spriteCache.size > 200) spriteCache.clear();
  spriteCache.set(key, s);
  return s;
}

export class CanvasPainter {
  constructor(ctx, res) {
    this.ctx = ctx;
    this.res = res;
    this.reset();
  }
  reset() { this.ctx.setTransform(this.res, 0, 0, this.res, 0, 0); }
  clear(w, h, bg) {
    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    this.ctx.clearRect(0, 0, this.ctx.canvas.width, this.ctx.canvas.height);
    this.reset();
    if (bg) this.rect(0, 0, w, h, bg);
  }
  rect(x, y, w, h, fill) { this.ctx.fillStyle = fill; this.ctx.fillRect(x, y, w, h); }
  strokeRect(x, y, w, h, stroke, width = 1) {
    const c = this.ctx;
    c.strokeStyle = stroke; c.lineWidth = width; c.setLineDash([]);
    c.strokeRect(x, y, w, h);
  }
  line(x1, y1, x2, y2, stroke, width = 1, dash = null) {
    const c = this.ctx;
    c.strokeStyle = stroke; c.lineWidth = width;
    c.setLineDash(dash || []);
    c.beginPath(); c.moveTo(x1, y1); c.lineTo(x2, y2); c.stroke();
    c.setLineDash([]);
  }
  polygon(pts, fill) {
    const c = this.ctx;
    c.fillStyle = fill;
    c.beginPath();
    pts.forEach(([x, y], i) => (i ? c.lineTo(x, y) : c.moveTo(x, y)));
    c.closePath(); c.fill();
  }
  clip(x, y, w, h) { this.ctx.save(); this.ctx.beginPath(); this.ctx.rect(x, y, w, h); this.ctx.clip(); }
  unclip() { this.ctx.restore(); }
  text(text, x, y, o) {
    const c = this.ctx;
    const segs = typeof text === 'string' ? [{ t: text, it: false, sh: 0 }] : text;
    const w = measureText(segs, o.size, o.family, o.weight);
    c.save();
    c.translate(x, y);
    if (o.rotate) c.rotate((o.rotate * Math.PI) / 180);
    let cx = o.align === 'center' ? -w / 2 : o.align === 'right' ? -w : 0;
    const by = baselineShift(o);
    c.fillStyle = o.color || INK;
    c.textBaseline = 'alphabetic';
    for (const g of segs) {
      const sz = segSize(g, o.size);
      c.font = fontStr(sz, o.family, o.weight, g.it);
      c.fillText(g.t, cx, by + g.sh * o.size * 0.28);
      cx += c.measureText(g.t).width;
    }
    c.restore();
  }
  // xs/ys in CSS px; ci indexes into colors.
  points(xs, ys, ci, count, colors, r, alpha) {
    const c = this.ctx;
    const res = this.res;
    c.save();
    c.setTransform(1, 0, 0, 1, 0, 0);
    c.globalAlpha = alpha;
    const sprites = colors.map((col) => sprite(col, r, res));
    for (let k = 0; k < count; k++) {
      const s = sprites[ci[k]];
      c.drawImage(s.c, xs[k] * res - s.half, ys[k] * res - s.half);
    }
    c.restore();
  }
}

function baselineShift(o) {
  if (o.baseline === 'middle') return o.size * 0.36;
  if (o.baseline === 'top') return o.size * 0.78;
  return 0;
}

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const r1 = (v) => Math.round(v * 10) / 10;

export class SvgPainter {
  constructor(w, h) {
    this.w = w; this.h = h;
    this.res = 2; // dedupe resolution for vector output
    this.parts = [];
    this.clipId = 0;
    this.clipOpen = false;
  }
  clear(w, h, bg) { if (bg) this.rect(0, 0, w, h, bg); }
  rect(x, y, w, h, fill) { this.parts.push(`<rect x="${r1(x)}" y="${r1(y)}" width="${r1(w)}" height="${r1(h)}" fill="${fill}"/>`); }
  strokeRect(x, y, w, h, stroke, width = 1) { this.parts.push(`<rect x="${r1(x)}" y="${r1(y)}" width="${r1(w)}" height="${r1(h)}" fill="none" stroke="${stroke}" stroke-width="${width}"/>`); }
  line(x1, y1, x2, y2, stroke, width = 1, dash = null) {
    this.parts.push(`<line x1="${r1(x1)}" y1="${r1(y1)}" x2="${r1(x2)}" y2="${r1(y2)}" stroke="${stroke}" stroke-width="${width}"${dash ? ` stroke-dasharray="${dash.join(' ')}"` : ''}/>`);
  }
  polygon(pts, fill) { this.parts.push(`<polygon fill="${fill}" points="${pts.map(([x, y]) => `${r1(x)},${r1(y)}`).join(' ')}"/>`); }
  clip(x, y, w, h) {
    const id = `clip${this.clipId++}`;
    this.parts.push(`<clipPath id="${id}"><rect x="${r1(x)}" y="${r1(y)}" width="${r1(w)}" height="${r1(h)}"/></clipPath><g clip-path="url(#${id})">`);
  }
  unclip() { this.parts.push('</g>'); }
  text(text, x, y, o) {
    const segs = typeof text === 'string' ? [{ t: text, it: false, sh: 0 }] : text;
    const anchor = o.align === 'center' ? 'middle' : o.align === 'right' ? 'end' : 'start';
    const by = baselineShift(o);
    const tr = o.rotate ? ` transform="rotate(${o.rotate} ${r1(x)} ${r1(y)})"` : '';
    let inner = '';
    let shift = 0;
    for (const g of segs) {
      const target = g.sh * o.size * 0.28;
      const dy = target - shift;
      shift = target;
      const attrs = [];
      if (dy) attrs.push(`dy="${r1(dy)}"`);
      if (g.sh) attrs.push(`font-size="${r1(o.size * 0.72)}"`);
      if (g.it) attrs.push('font-style="italic"');
      inner += attrs.length ? `<tspan ${attrs.join(' ')}>${esc(g.t)}</tspan>` : esc(g.t);
    }
    this.parts.push(`<text x="${r1(x)}" y="${r1(y + by)}" font-family="${esc(o.family)}" font-size="${o.size}"${o.weight && o.weight !== 400 ? ` font-weight="${o.weight}"` : ''} fill="${o.color || INK}" text-anchor="${anchor}"${tr}>${inner}</text>`);
  }
  points(xs, ys, ci, count, colors, r, alpha) {
    let cur = -1;
    const out = [];
    const rr = r1(r);
    for (let k = 0; k < count; k++) {
      if (ci[k] !== cur) {
        if (cur >= 0) out.push('</g>');
        cur = ci[k];
        out.push(`<g fill="${colors[cur]}"${alpha < 1 ? ` fill-opacity="${alpha}"` : ''}>`);
      }
      out.push(`<circle cx="${r1(xs[k])}" cy="${r1(ys[k])}" r="${rr}"/>`);
    }
    if (cur >= 0) out.push('</g>');
    this.parts.push(out.join(''));
  }
  toString() {
    return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${this.w}" height="${this.h}" viewBox="0 0 ${this.w} ${this.h}">\n${this.parts.join('\n')}\n</svg>\n`;
  }
}

// ---------------------------------------------------------------------------
// Shared helpers

let pool = { n: 0, xs: null, ys: null, ci: null };
function buffers(n) {
  if (pool.n < n) pool = { n, xs: new Float32Array(n), ys: new Float32Array(n), ci: new Uint8Array(n) };
  return pool;
}

// Drop points that land on a cell already painted with the same color. Output
// is visually identical and much cheaper to draw (and to store as SVG).
function dedupe(b, count, W, H, res) {
  res = Math.max(0.05, res);
  const gw = Math.ceil(W * res) + 2, gh = Math.ceil(H * res) + 2;
  const grid = new Uint8Array(gw * gh);
  const { xs, ys, ci } = b;
  let m = 0;
  for (let k = 0; k < count; k++) {
    const gx = (xs[k] * res + 1) | 0, gy = (ys[k] * res + 1) | 0;
    const c = ci[k] + 1;
    if (gx >= 0 && gy >= 0 && gx < gw && gy < gh) {
      const key = gy * gw + gx;
      if (grid[key] === c) continue;
      grid[key] = c;
    }
    xs[m] = xs[k]; ys[m] = ys[k]; ci[m] = ci[k];
    m++;
  }
  return m;
}

function lowerBound(arr, v, lo = 0, hi = arr.length) {
  while (lo < hi) { const mid = (lo + hi) >>> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; }
  return lo;
}

const DASH = { solid: null, dashed: [6, 4], dotted: [1.5, 3] };

// Cells per CSS px for de-duplication. Exports use one cell per device pixel.
// The on-screen preview (painter.fast) uses cells of about half a point radius:
// opaque dots that close together look the same, and it is much cheaper.
function gridRes(P, r) {
  if (!P.fast) return P.res;
  return Math.min(P.res, 2 / Math.max(r, 0.5));
}

// Threshold settings are kept separately for p-values and for raw values.
export function lineSpec(o, key, scale) {
  return scale === 'raw' ? o[key].raw : o[key].p;
}

export function thresholdY(spec, scale) {
  if (spec === '' || spec === null || spec === undefined) return null;
  const v = typeof spec === 'number' ? spec : parseFloat(spec);
  if (!Number.isFinite(v)) return null;
  if (scale === 'raw') return v;
  if (v <= 0 || v >= 1) return null;
  return -Math.log10(v);
}

function axisFont(o) { return { size: o.fontSize, family: FONTS[o.fontFamily] || FONTS.sans, color: INK }; }

function drawYAxis(P, L, o, ticks, step) {
  const f = axisFont(o);
  const tl = 5;
  const tsz = f.size * 0.92;
  let maxW = 0;
  const x = L.x0;
  P.line(x, L.sy(ticks[0]), x, L.sy(ticks[ticks.length - 1]), INK, 1);
  for (const t of ticks) {
    const y = L.sy(t);
    P.line(x - tl, y, x, y, INK, 1);
    const s = fmtTick(t, step);
    maxW = Math.max(maxW, measureText(s, tsz, f.family));
    P.text(s, x - tl - 3, y, { ...f, size: tsz, align: 'right', baseline: 'middle' });
  }
  return maxW;
}

function drawYTitle(P, L, o, label, tickW) {
  const f = axisFont(o);
  const x = L.x0 - 5 - 3 - tickW - 8;
  P.text(parseRich(label), x, (L.y0 + L.y1) / 2, { ...f, align: 'center', rotate: -90 });
}

function drawTitle(P, L, o, W, title) {
  if (!title) return;
  const f = axisFont(o);
  P.text(parseRich(title), (L.x0 + L.x1) / 2, o.margin.t + f.size * 1.3, { ...f, size: f.size * 1.3, weight: 700, align: 'center' });
}

function titleSpace(o, title) { return title ? o.fontSize * 1.3 + 10 : 0; }

// ---------------------------------------------------------------------------
// Manhattan

export function manhattanLayout(D, o, view, W, H, hasLabels) {
  const m = o.margin;
  const L = { W, H, x0: m.l, x1: W - m.r, y0: m.t + titleSpace(o, o.title), y1: H - m.b };
  if (L.x1 - L.x0 < 20) L.x1 = L.x0 + 20;
  if (L.y1 - L.y0 < 20) L.y1 = L.y0 + 20;
  let d0, d1;
  if (view) { d0 = view.x0; d1 = view.x1; }
  else {
    const span = D.genomeEnd - D.genomeStart || 1;
    const pad = span * 0.015;
    d0 = D.genomeStart - pad; d1 = D.genomeEnd + pad;
  }
  L.d0 = d0; L.d1 = d1;
  L.lo = lowerBound(D.pos, d0);
  L.hi = lowerBound(D.pos, d1 + 1e-9);

  let maxV = -Infinity, minV = Infinity;
  if (!view) { maxV = D.maxY; minV = D.minY; }
  else for (let i = L.lo; i < L.hi; i++) { const v = D.y[i]; if (v > maxV) maxV = v; if (v < minV) minV = v; }
  if (!Number.isFinite(maxV)) { maxV = 1; minV = 0; }

  const lines = [];
  for (const k of ['suggestive', 'genomewide']) {
    if (!o[k].on) continue;
    const v = thresholdY(lineSpec(o, k, D.scale), D.scale);
    if (v !== null) lines.push({ v, color: o[k].color });
  }
  L.lines = lines;

  let top = maxV, bottom = minV;
  for (const l of lines) {
    if (l.v > top && (!view || l.v < maxV * 3 + 1)) top = l.v;
    if (l.v < bottom) bottom = l.v;
  }
  let ymin = 0, ymax;
  if (D.scale === 'raw') {
    // raw values can be negative and on any scale: round outward to a nice step
    const step = niceStep((top - Math.min(0, bottom)) || 1, 6);
    ymin = bottom < 0 ? Math.floor(bottom / step) * step : 0;
    ymax = Math.ceil((top + (hasLabels ? (top - ymin) * 0.06 : 0)) / step) * step;
    if (ymax <= ymin) ymax = ymin + step;
  } else {
    ymax = top > 0 ? Math.ceil(top * (hasLabels ? 1.06 : 1) + 1e-9) : 1;
  }
  const custom = parseFloat(o.yMax);
  if (Number.isFinite(custom) && custom > ymin) ymax = custom;
  L.ymin = ymin; L.ymax = ymax;

  const xs = (L.x1 - L.x0) / (d1 - d0);
  const ys = (L.y1 - L.y0) / (ymax - ymin);
  L.sx = (p) => L.x0 + (p - d0) * xs;
  L.sy = (v) => L.y1 - (v - ymin) * ys;
  L.syc = (v) => L.sy(Math.min(v, ymax)); // clamped: points above a custom max sit at the top
  L.ix = (px) => d0 + (px - L.x0) / xs;

  // x ticks
  const visibleChroms = D.chroms.filter((c) => c.end >= d0 && c.start <= d1);
  const mid = (d0 + d1) / 2;
  let single = D.chroms.length === 1 ? D.chroms[0] : null;
  if (!single && view && visibleChroms.length === 1) single = visibleChroms[0];
  if (!single && view) {
    const c = D.chroms.find((ch) => mid >= ch.start && mid <= ch.end);
    if (c && d0 >= c.start - (c.end - c.start) * 0.05 && d1 <= c.end + (c.end - c.start) * 0.05) single = c;
  }
  L.single = single;
  const labels = customChromLabels(o, D);
  if (single) {
    if (view) { // don't let neighboring chromosomes leak in at the edges
      L.lo = Math.max(L.lo, single.from);
      L.hi = Math.min(L.hi, single.to);
    }
    const b0 = d0 - single.base, b1 = d1 - single.base;
    const n = Math.max(2, Math.round((L.x1 - L.x0) / 90));
    const step = niceStep(b1 - b0, n);
    const ticks = niceTicks(Math.max(0, b0), b1, n);
    const [div, unit] = bpUnit(b1 - b0);
    L.xticks = ticks.map((b) => ({ x: L.sx(b + single.base), label: fmtTick(b / div, step / div) }));
    L.autoXlab = `Chromosome ${labels ? labels[D.chroms.indexOf(single)] : chromLabel(single, o)} position (${unit})`;
    L.xAxisFull = false;
  } else {
    L.xticks = [];
    D.chroms.forEach((c, k) => {
      const center = (c.start + c.end) / 2;
      if (center >= d0 && center <= d1) L.xticks.push({ x: L.sx(center), label: labels ? labels[k] : chromLabel(c, o) });
    });
    L.autoXlab = 'Chromosome';
  }
  return L;
}

// Pick a unit so position ticks stay readable.
export function bpUnit(span) {
  if (span >= 2e6) return [1e6, 'Mb'];
  if (span >= 2e3) return [1e3, 'kb'];
  return [1, 'bp'];
}

function customChromLabels(o, D) {
  if (!o.chrLabels || !o.chrLabels.trim()) return null;
  const labs = o.chrLabels.split(/[,\s]+/).filter(Boolean);
  return labs.length === D.chroms.length ? labs : null;
}

export function chromColors(o) {
  return o.colors && o.colors.length ? o.colors : ['#1a1a1a', '#999999'];
}

// extra: { highlight: number[] (sorted indices), labels: [{i, text}], autoYlab }
export function renderManhattan(P, D, o, view, W, H, extra = {}) {
  const hasLabels = extra.labels && extra.labels.length > 0;
  const L = manhattanLayout(D, o, view, W, H, hasLabels);
  const f = axisFont(o);
  P.clear(W, H, o.transparent ? null : o.bg);

  // alternating chromosome shading
  if (o.shade && !L.single) {
    P.clip(L.x0, L.y0, L.x1 - L.x0, L.y1 - L.y0);
    D.chroms.forEach((c, k) => {
      if (k % 2 === 1) {
        const prev = D.chroms[k - 1], next = D.chroms[k + 1];
        const a = prev ? (prev.end + c.start) / 2 : c.start;
        const b = next ? (c.end + next.start) / 2 : c.end;
        P.rect(L.sx(a), L.y0, L.sx(b) - L.sx(a), L.y1 - L.y0, o.shadeColor || '#f1f1f1');
      }
    });
    P.unclip();
  }

  const yTickCount = Math.max(2, Math.round((L.y1 - L.y0) / 48));
  const integer = D.scale !== 'raw' && L.ymax - L.ymin >= 3;
  const yticks = niceTicks(L.ymin, L.ymax, yTickCount, integer);
  const ystep = yticks.length > 1 ? yticks[1] - yticks[0] : 1;
  if (o.grid) for (const t of yticks) P.line(L.x0, L.sy(t), L.x1, L.sy(t), '#e4e4e4', 1);

  // points
  const r = o.pointSize / 2;
  const colors = chromColors(o);
  const nc = colors.length;
  const b = buffers(Math.max(1, L.hi - L.lo));
  let count = 0;
  const { xs, ys, ci } = b;
  const ymaxClip = L.ymax + (L.ymax - L.ymin) * 0.002;
  const capped = [];
  for (let i = L.lo; i < L.hi; i++) {
    const v = D.y[i];
    if (v > ymaxClip) { if (capped.length < 20000) capped.push(i); continue; }
    if (v < L.ymin) continue;
    xs[count] = L.sx(D.pos[i]);
    ys[count] = L.sy(v);
    ci[count] = D.chrIdx[i] % nc;
    count++;
  }
  count = dedupe(b, count, W, H, gridRes(P, r));
  P.points(xs, ys, ci, count, colors, r, o.opacity);
  L.drawn = count;
  L.capped = capped.length;
  if (capped.length) {
    const t = Math.max(2.5, r * 1.3);
    const yt = L.sy(L.ymax);
    let lastX = -Infinity, lastC = -1;
    for (const i of capped) {
      const x = L.sx(D.pos[i]);
      const c = D.chrIdx[i] % nc;
      if (Math.abs(x - lastX) < 0.5 && c === lastC) continue;
      lastX = x; lastC = c;
      P.polygon([[x, yt - t * 1.2], [x - t, yt + t * 0.6], [x + t, yt + t * 0.6]], colors[c]);
    }
  }

  // significance lines
  const dash = DASH[o.lineStyle] || null;
  for (const l of L.lines) {
    if (l.v >= L.ymin && l.v <= L.ymax) P.line(L.x0, L.sy(l.v), L.x1, L.sy(l.v), l.color, o.lineWidth, dash);
  }

  // highlighted SNPs
  if (extra.highlight && extra.highlight.length) {
    let hc = 0;
    const hb = { xs: new Float32Array(extra.highlight.length), ys: new Float32Array(extra.highlight.length), ci: new Uint8Array(extra.highlight.length) };
    for (const i of extra.highlight) {
      if (i < L.lo || i >= L.hi) continue;
      const v = D.y[i];
      if (v < L.ymin) continue;
      hb.xs[hc] = L.sx(D.pos[i]); hb.ys[hc] = L.syc(v); hc++;
    }
    P.points(hb.xs, hb.ys, hb.ci, hc, [o.highlightColor], r * (o.highlightScale || 1), Math.max(o.opacity, 0.85));
  }

  // labels
  if (hasLabels) {
    const items = [];
    for (const lab of extra.labels) {
      const i = lab.i;
      if (i < L.lo || i >= L.hi) continue;
      const v = D.y[i];
      if (v < L.ymin) continue;
      items.push({ px: L.sx(D.pos[i]), py: L.syc(v) - (v > ymaxClip ? r : 0), text: lab.text, v, pinned: lab.pinned });
    }
    L.labelsSkipped = placeLabels(P, items, L, o, r);
  }

  // axes
  const tickW = drawYAxis(P, L, o, yticks, ystep);
  const axisY = L.y1 + Math.max(3, r + 1);
  const tl = 5;
  if (L.xticks.length) {
    const first = L.xticks[0].x, last = L.xticks[L.xticks.length - 1].x;
    P.line(L.single ? L.x0 : first, axisY, L.single ? L.x1 : last, axisY, INK, 1);
    const tsz = f.size * 0.92;
    let lastRight = -Infinity;
    for (const t of L.xticks) {
      P.line(t.x, axisY, t.x, axisY + tl, INK, 1);
      const w = measureText(t.label, tsz, f.family);
      if (t.x - w / 2 > lastRight + 3) {
        P.text(t.label, t.x, axisY + tl + 3 + tsz * 0.8, { ...f, size: tsz, align: 'center' });
        lastRight = t.x + w / 2;
      }
    }
  }
  const xlab = o.xlab || L.autoXlab;
  P.text(parseRich(xlab), (L.x0 + L.x1) / 2, axisY + tl + 3 + f.size * 0.92 + 8 + f.size, { ...f, align: 'center' });
  const ylab = o.ylab || extra.autoYlab || '−log_{10}(*p*)';
  drawYTitle(P, L, o, ylab, tickW);
  drawTitle(P, L, o, W, o.title);
  return L;
}

function placeLabels(P, items, L, o, r) {
  const fam = FONTS[o.fontFamily] || FONTS.sans;
  const sz = o.labelSize;
  const h = sz * 1.05;
  const placed = [];
  const overlaps = (a) => placed.some((b) => a.x < b.x + b.w && a.x + a.w > b.x && a.y < b.y + b.h && a.y + a.h > b.y);
  items.sort((a, b) => (b.pinned - a.pinned) || (b.v - a.v));
  let skipped = 0;
  for (const it of items) {
    const w = measureText(it.text, sz, fam) + 2;
    // Candidate spots: just above the point, beside it, then higher up (with a
    // leader line), straight up first and then up and to either side.
    const side = w / 2 + r + 2;
    const cands = [{ dx: 0, lift: 0 }, { dx: side, beside: true }, { dx: -side, beside: true }];
    for (let k = 1; k < 14; k++) cands.push({ dx: 0, lift: k });
    for (let k = 1; k < 8; k++) cands.push({ dx: side, lift: k }, { dx: -side, lift: k });
    let ok = null;
    for (const c of cands) {
      const base = c.beside ? it.py + sz * 0.35 : it.py - r - 3 - c.lift * (h + 1);
      const x = Math.min(Math.max(it.px + c.dx - w / 2, 1), L.W - w - 1);
      const rect = { x, y: base - sz * 0.8, w, h };
      if (rect.y < 1 || rect.y + h > L.y1) continue;
      if (!overlaps(rect)) { ok = { rect, base, leader: !c.beside && c.lift > 0 }; break; }
    }
    if (!ok) { skipped++; continue; }
    placed.push(ok.rect);
    if (ok.leader) P.line(it.px, it.py - r - 1, ok.rect.x + w / 2, ok.rect.y + h - 1, '#9a9a9a', 0.75);
    P.text(it.text, ok.rect.x + w / 2, ok.base, { size: sz, family: fam, color: '#262626', align: 'center' });
  }
  return skipped;
}

// ---------------------------------------------------------------------------
// Q-Q

let bandCache = { n: 0, pts: null };
function ciBand(n, a) {
  if (bandCache.n === n) return bandCache.pts;
  const denom = n + 1 - 2 * a;
  const emax = -Math.log10((1 - a) / denom);
  const pts = [];
  const steps = 140;
  let lastI = -1;
  for (let s = 0; s <= steps; s++) {
    const e = (emax * s) / steps;
    let i = Math.round(10 ** -e * denom + a);
    i = Math.min(n, Math.max(1, i));
    if (i === lastI) continue;
    lastI = i;
    const ex = -Math.log10((i - a) / denom);
    const lo = -Math.log10(qOrderStat(0.975, i, n));
    const hi = -Math.log10(qOrderStat(0.025, i, n));
    pts.push([ex, lo, hi]);
  }
  bandCache = { n, pts };
  return pts;
}

export function qqLayout(D, q, o, W, H) {
  const m = o.margin;
  const L = { W, H, x0: m.l, x1: W - m.r, y0: m.t + titleSpace(o, q.title), y1: H - m.b };
  const n = D.n;
  const a = n <= 10 ? 3 / 8 : 1 / 2;
  L.a = a;
  L.denom = n + 1 - 2 * a;
  const emax = -Math.log10((1 - a) / L.denom);
  const omax = D.sortedY[n - 1];
  const custom = parseFloat(q.yMax);
  const ytop = Number.isFinite(custom) && custom > 0 ? custom : Math.max(omax, 1e-3);
  // R's default xaxs/yaxs = "r": extend the range by 4% on each side
  L.dx0 = -emax * 0.04; L.dx1 = emax * 1.04;
  L.dy0 = -ytop * 0.04; L.dy1 = ytop * 1.04;
  const xs = (L.x1 - L.x0) / (L.dx1 - L.dx0);
  const ys = (L.y1 - L.y0) / (L.dy1 - L.dy0);
  L.sx = (v) => L.x0 + (v - L.dx0) * xs;
  L.sy = (v) => L.y1 - (v - L.dy0) * ys;
  L.emax = emax; L.omax = ytop;
  return L;
}

export function expectedAt(L, k, n) {
  // k-th smallest p (0-based) -> expected -log10 p
  return -Math.log10((k + 1 - L.a) / L.denom);
}

export function renderQQ(P, D, o, W, H, extra = {}) {
  const q = o.qq;
  const L = qqLayout(D, q, o, W, H);
  const f = axisFont(o);
  P.clear(W, H, o.transparent ? null : o.bg);
  const n = D.n;

  P.clip(L.x0, L.y0, L.x1 - L.x0, L.y1 - L.y0);
  if (q.ci && n >= 3) {
    const band = ciBand(n, L.a);
    const pts = band.map(([e, lo]) => [L.sx(e), L.sy(lo)]).concat(band.slice().reverse().map(([e, , hi]) => [L.sx(e), L.sy(hi)]));
    P.polygon(pts, q.ciColor);
  }
  const big = Math.max(L.dx1, L.dy1) * 2;
  P.line(L.sx(-big), L.sy(-big), L.sx(big), L.sy(big), q.lineColor, o.lineWidth);

  const b = buffers(n);
  const { xs, ys, ci } = b;
  let count = 0;
  const sorted = D.sortedY;
  for (let k = 0; k < n; k++) {
    // sortedY ascending in -log10 p, so index k has the (n-k)-th smallest p
    const obs = sorted[k];
    if (obs > L.dy1) continue;
    const e = -Math.log10((n - k - L.a) / L.denom);
    xs[count] = L.sx(e); ys[count] = L.sy(obs); ci[count] = 0; count++;
  }
  count = dedupe(b, count, W, H, gridRes(P, o.pointSize / 2));
  P.points(xs, ys, ci, count, [q.color], o.pointSize / 2, o.opacity);
  P.unclip();

  // highlighted SNPs on the Q-Q plot
  if (extra.highlightRanks && extra.highlightRanks.length) {
    const hr = extra.highlightRanks;
    const hb = { xs: new Float32Array(hr.length), ys: new Float32Array(hr.length), ci: new Uint8Array(hr.length) };
    let hc = 0;
    for (const [rank, v] of hr) {
      if (v > L.dy1) continue;
      hb.xs[hc] = L.sx(-Math.log10((rank + 1 - L.a) / L.denom)); hb.ys[hc] = L.sy(v); hc++;
    }
    P.points(hb.xs, hb.ys, hb.ci, hc, [o.highlightColor], o.pointSize / 2, 1);
  }

  P.strokeRect(L.x0, L.y0, L.x1 - L.x0, L.y1 - L.y0, INK, 1);
  const tsz = f.size * 0.92;
  const xt = niceTicks(0, L.dx1, Math.max(2, Math.round((L.x1 - L.x0) / 70)), L.emax >= 3);
  const xstep = xt.length > 1 ? xt[1] - xt[0] : 1;
  for (const t of xt) {
    const x = L.sx(t);
    P.line(x, L.y1, x, L.y1 + 5, INK, 1);
    P.text(fmtTick(t, xstep), x, L.y1 + 8 + tsz * 0.8, { ...f, size: tsz, align: 'center' });
  }
  const yt = niceTicks(0, L.dy1, Math.max(2, Math.round((L.y1 - L.y0) / 50)), L.omax >= 3);
  const ystep = yt.length > 1 ? yt[1] - yt[0] : 1;
  let tickW = 0;
  for (const t of yt) {
    const y = L.sy(t);
    P.line(L.x0 - 5, y, L.x0, y, INK, 1);
    const s = fmtTick(t, ystep);
    tickW = Math.max(tickW, measureText(s, tsz, f.family));
    P.text(s, L.x0 - 8, y, { ...f, size: tsz, align: 'right', baseline: 'middle' });
  }
  P.text(parseRich(q.xlab || 'Expected −log_{10}(*p*)'), (L.x0 + L.x1) / 2, L.y1 + 8 + tsz + 8 + f.size, { ...f, align: 'center' });
  drawYTitle(P, L, o, q.ylab || 'Observed −log_{10}(*p*)', tickW);
  drawTitle(P, L, o, W, q.title);

  if (q.showLambda && D.lambda !== null && Number.isFinite(D.lambda)) {
    P.text(parseRich(`*λ*_{GC} = ${D.lambda.toFixed(3)}`), L.x0 + 10, L.y0 + 10 + f.size * 0.8, { ...f, align: 'left' });
  }
  return L;
}
