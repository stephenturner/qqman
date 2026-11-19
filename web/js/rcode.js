// Build the closest qqman R call for the current settings.
import { plainText } from './plot.js';

const DEFAULT_COLS = ['#1a1a1a', '#999999'];
const RNAMES = { '#1a1a1a': 'gray10', '#999999': 'gray60', '#000000': 'black', '#ff0000': 'red', '#0000ff': 'blue', '#00cd00': 'green3', '#bebebe': 'gray' };

const q = (s) => JSON.stringify(String(s));
const rcol = (c) => q(RNAMES[c.toLowerCase()] || c.toUpperCase());
const isName = (s) => /^[A-Za-z.][A-Za-z0-9._]*$/.test(s) && !/^\.\d/.test(s);
const ref = (obj, col) => (isName(col) ? `${obj}$${col}` : `${obj}[[${q(col)}]]`);
const num = (v) => String(+(+v).toPrecision(6));

function call(fn, args, indent = '') {
  if (args.length === 0) return `${fn}()`;
  const one = `${fn}(${args.join(', ')})`;
  if (one.length <= 78 && args.length <= 3) return one;
  const pad = ' '.repeat(fn.length + 1);
  return `${fn}(${args.join(`,\n${indent}${pad}`)})`;
}

export function rCode(state, o, soi) {
  const { raw, map, data: D } = state;
  if (!raw || !D) return '# Load some data to see the equivalent R code.';
  const cols = raw.columns;
  const nm = (j) => cols[j].name;
  const out = ['library(qqman)', ''];
  let obj = 'gwas';
  if (state.example === 'gwasResults') {
    obj = 'gwasResults';
    out.push('# gwasResults is the example data that comes with qqman');
  } else {
    // simulated and pasted data have no file on disk to point at
    const real = raw.name && !state.example && !raw.name.startsWith('Pasted');
    if (!real) out.push('# point this at your own results file');
    const file = real ? raw.name : 'results.txt';
    if (raw.delimiter === ',') out.push(`${obj} <- read.csv(${q(file)})`);
    else if (raw.delimiter === '\t') out.push(`${obj} <- read.delim(${q(file)})`);
    else out.push(`${obj} <- read.table(${q(file)}, header = TRUE)`);
  }
  const chr = nm(map.chr), bp = nm(map.bp), p = nm(map.p);
  const snp = map.snp >= 0 ? nm(map.snp) : null;

  if (cols[map.chr].strs) {
    out.push('', '# qqman needs a numeric chromosome column');
    const c = ref(obj, chr);
    out.push(`${c} <- sub("^chr", "", ${c}, ignore.case = TRUE)`);
    out.push(`${c}[${c} == "X"] <- 23`);
    out.push(`${c}[${c} == "Y"] <- 24`);
    out.push(`${c}[${c} %in% c("M", "MT")] <- 26`);
    out.push(`${c} <- suppressWarnings(as.numeric(${c}))`);
  }
  const dropped = raw.nrows - D.n;
  if (dropped > 0) {
    out.push('', `# ${dropped.toLocaleString('en-US')} rows were dropped by the data checks`);
    const cond = [`!is.na(${ref(obj, chr)})`, `!is.na(${ref(obj, bp)})`, `!is.na(${ref(obj, p)})`];
    if (map.scale === 'p') cond.push(`${ref(obj, p)} > 0`, `${ref(obj, p)} <= 1`);
    out.push(`${obj} <- ${obj}[${cond.join(' &\n' + ' '.repeat(obj.length * 2 + 5))}, ]`);
  }

  const args = [obj];
  if (chr !== 'CHR') args.push(`chr = ${q(chr)}`);
  if (bp !== 'BP') args.push(`bp = ${q(bp)}`);
  if (p !== 'P') args.push(`p = ${q(p)}`);
  if (snp && snp !== 'SNP') args.push(`snp = ${q(snp)}`);
  const colors = o.colors.map((c) => c.toLowerCase());
  if (colors.join() !== DEFAULT_COLS.join()) args.push(`col = c(${colors.map(rcol).join(', ')})`);
  // With logp = FALSE, qqman draws the lines at the raw values given.
  const lineArg = (key, name, def) => {
    const spec = map.scale === 'raw' ? o[key].raw : o[key].p;
    const v = parseFloat(spec);
    if (!o[key].on || !Number.isFinite(v)) { args.push(`${name} = FALSE`); return; }
    if (map.scale === 'raw') args.push(`${name} = ${num(v)}`);
    else if (map.scale === 'mlog10') args.push(`${name} = ${num(-Math.log10(v))}`);
    else if (v !== def) args.push(`${name} = -log10(${num(v)})`);
  };
  lineArg('suggestive', 'suggestiveline', 1e-5);
  lineArg('genomewide', 'genomewideline', 5e-8);

  let pre = [];
  const hl = state.highlightNames || [];
  if (hl.length) {
    if (state.example === 'gwasResults' && soi && hl.length === soi.length && hl.every((s, i) => s === soi[i])) {
      args.push('highlight = snpsOfInterest');
    } else if (hl.length <= 12) {
      args.push(`highlight = c(${hl.map(q).join(', ')})`);
    } else {
      pre.push(`snps <- c(${hl.slice(0, 300).map(q).join(', ')}${hl.length > 300 ? ' # …' : ''})`);
      args.push('highlight = snps');
    }
  }
  if (map.scale !== 'p') args.push('logp = FALSE');
  if (o.annotate) {
    const v = parseFloat(map.scale === 'raw' ? o.annotateRaw : o.annotateP);
    if (Number.isFinite(v)) {
      args.push(`annotatePval = ${map.scale === 'mlog10' ? num(-Math.log10(v)) : num(v)}`);
      if (!o.annotateTop) args.push('annotateTop = FALSE');
    }
  }
  if (o.chrLabels && o.chrLabels.trim()) {
    const labs = o.chrLabels.split(/[,\s]+/).filter(Boolean);
    if (labs.length === D.chroms.length) args.push(`chrlabs = c(${labs.map(q).join(', ')})`);
  } else if (o.plinkNames && D.chroms.some((c) => !c.letter && c.rank >= 23 && c.rank <= 26)) {
    const map26 = { 23: 'X', 24: 'Y', 25: 'XY', 26: 'MT' };
    args.push(`chrlabs = c(${D.chroms.map((c) => q(!c.letter && map26[c.rank] ? map26[c.rank] : c.label)).join(', ')})`);
  }
  if (o.title) args.push(`main = ${q(plainText(o.title))}`);
  if (o.xlab) args.push(`xlab = ${q(plainText(o.xlab))}`);
  if (o.ylab) args.push(`ylab = ${q(plainText(o.ylab))}`);
  else if (map.scale === 'mlog10') args.push('ylab = expression(-log[10](italic(p)))');
  const ym = parseFloat(o.yMax);
  if (Number.isFinite(ym)) args.push(`ylim = c(0, ${num(ym)})`);
  if (o.pointSize !== 4) args.push(`cex = ${num(o.pointSize / 4)}`);
  if (o.fontSize !== 12) args.push(`cex.axis = ${num(o.fontSize / 12)}`, `cex.lab = ${num(o.fontSize / 12)}`);

  out.push('');
  if (pre.length) out.push(...pre, '');
  out.push(`png("manhattan.png", width = ${o.width}, height = ${o.height})`);
  out.push(call('manhattan', args));
  out.push('dev.off()');
  if (o.highlightColor.toLowerCase() !== '#00cd00' && hl.length) out.push('# note: qqman always draws highlighted SNPs in green3');

  if (map.scale !== 'raw') {
    const pv = map.scale === 'mlog10' ? `10^-${ref(obj, p)}` : ref(obj, p);
    const qargs = [pv];
    if (o.qq.title) qargs.push(`main = ${q(plainText(o.qq.title))}`);
    out.push('', `png("qq.png", width = ${o.qq.width}, height = ${o.qq.height})`);
    out.push(call('qq', qargs));
    out.push('dev.off()');
  }
  return out.join('\n');
}

export function highlightR(code) {
  const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  return code.split('\n').map((line) => {
    const ci = line.search(/#(?=(?:[^"]*"[^"]*")*[^"]*$)/);
    const body = ci >= 0 ? line.slice(0, ci) : line;
    const comment = ci >= 0 ? line.slice(ci) : '';
    let h = esc(body)
      .replace(/("(?:[^"\\]|\\.)*")/g, '<span class="s">$1</span>')
      .replace(/\b([A-Za-z_.][A-Za-z0-9_.]*)(?=\()/g, '<span class="f">$1</span>')
      .replace(/\b(TRUE|FALSE|\d+(?:\.\d+)?(?:e-?\d+)?)\b(?![^<]*<\/span>)/g, '<span class="n">$1</span>');
    if (comment) h += `<span class="c">${esc(comment)}</span>`;
    return h;
  }).join('\n');
}
