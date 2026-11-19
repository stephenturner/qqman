// Numerical helpers: normal and beta distributions, ticks, formatting.

// Inverse standard normal CDF (Acklam), refined with one Halley step.
export function qnorm(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  const a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.38357751867269e2, -3.066479806614716e1, 2.506628277459239];
  const b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1];
  const c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416];
  const pl = 0.02425;
  let q, r, x;
  if (p < pl) {
    q = Math.sqrt(-2 * Math.log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p <= 1 - pl) {
    q = p - 0.5;
    r = q * q;
    x = ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q) / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  const e = 0.5 * erfc(-x / Math.SQRT2) - p;
  const u = e * Math.sqrt(2 * Math.PI) * Math.exp((x * x) / 2);
  return x - u / (1 + (x * u) / 2);
}

// Complementary error function, relative error < 1.2e-7 (Numerical Recipes erfcc).
export function erfc(x) {
  const z = Math.abs(x);
  const t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}

// -log10 of the two-sided normal p-value for a z statistic, stable for large |z|.
export function zToMlog10p(z) {
  const a = Math.abs(z);
  if (a < 8) return -Math.log10(erfc(a / Math.SQRT2));
  // Asymptotic tail: P(|Z|>a) ~ 2 phi(a)/a * (1 - 1/a^2 + 3/a^4)
  const lnp = Math.log(2) - (a * a) / 2 - Math.log(a) - 0.5 * Math.log(2 * Math.PI) + Math.log(1 - 1 / (a * a) + 3 / (a ** 4));
  return -lnp / Math.LN10;
}

export function lambdaGC(medianP) {
  const z = qnorm(medianP / 2);
  return (z * z) / 0.4549364231195724;
}

export function lgamma(x) {
  const g = 7;
  const cof = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
  x -= 1;
  let a = cof[0];
  const t = x + g + 0.5;
  for (let i = 1; i < g + 2; i++) a += cof[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

function betacf(a, b, x) {
  const MAXIT = 5000, EPS = 1e-12, FPMIN = 1e-300;
  const qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = 1 - (qab * x) / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}

// log B(a, b) for positive integers; exact-ish even when one argument is huge.
function lbetaInt(a, b) {
  const m = Math.min(a, b), M = Math.max(a, b);
  if (m > 500) return lgamma(a) + lgamma(b) - lgamma(a + b);
  let s = lgamma(m);
  for (let k = 0; k < m; k++) s -= Math.log(M + k);
  return s;
}

// Regularized incomplete beta I_x(a, b).
export function ibeta(x, a, b, lbeta) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const lb = lbeta === undefined ? lgamma(a) + lgamma(b) - lgamma(a + b) : lbeta;
  const lbt = -lb + a * Math.log(x) + b * Math.log1p(-x);
  if (x < (a + 1) / (a + b + 2)) return (Math.exp(lbt) * betacf(a, b, x)) / a;
  return 1 - (Math.exp(lbt) * betacf(b, a, 1 - x)) / b;
}

// Quantile of the i-th smallest of n uniform order statistics, Beta(i, n-i+1).
export function qOrderStat(p, i, n) {
  const a = i, b = n - i + 1;
  if (Math.min(a, b) > 300) {
    const s = a + b;
    const mu = a / s;
    const sd = Math.sqrt((a * b) / (s * s * (s + 1)));
    return Math.min(1, Math.max(1e-300, mu + qnorm(p) * sd));
  }
  const lb = lbetaInt(a, b);
  let lo = -745, hi = 0; // bisection on log(x)
  let mid = Math.log(a / (a + b));
  for (let it = 0; it < 90; it++) {
    const v = ibeta(Math.exp(mid), a, b, lb);
    if (v < p) lo = mid; else hi = mid;
    mid = (lo + hi) / 2;
    if (hi - lo < 1e-6) break;
  }
  return Math.exp(mid);
}

// "Nice" axis ticks covering [min, max].
export function niceStep(range, count) {
  const raw = range / Math.max(1, count);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10;
  return step * mag;
}

export function niceTicks(min, max, count, integer = false) {
  if (!(max > min)) return [min];
  let step = niceStep(max - min, count);
  if (integer) step = Math.max(1, Math.round(step));
  const start = Math.ceil(min / step - 1e-9) * step;
  const ticks = [];
  for (let v = start; v <= max + step * 1e-9; v += step) ticks.push(+v.toFixed(10));
  return ticks;
}

export function fmtTick(v, step) {
  const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return v.toFixed(dec);
}

// Format a p-value from its -log10, which survives values below 1e-308.
export function fmtPFromLog(L) {
  if (!Number.isFinite(L)) return L > 0 ? '0' : String(L);
  if (L < 3) return String(+(10 ** -L).toPrecision(3));
  const e = Math.floor(-L);
  let m = 10 ** (-L - e);
  let ex = e;
  if (m >= 9.995) { m /= 10; ex += 1; }
  return `${m.toFixed(2)}e${ex}`;
}

export function fmtNum(v, digits = 4) {
  if (v === null || v === undefined || v !== v) return 'NA';
  if (Number.isInteger(v) && Math.abs(v) < 1e15) return String(v);
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return v.toExponential(digits - 1);
  return String(+v.toPrecision(digits));
}

export const fmtInt = (v) => Math.round(v).toLocaleString('en-US');

// Seeded PRNG.
export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function normalGen(rng) {
  let spare = null;
  return function () {
    if (spare !== null) { const s = spare; spare = null; return s; }
    let u, v, s;
    do { u = rng() * 2 - 1; v = rng() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
    const m = Math.sqrt((-2 * Math.log(s)) / s);
    spare = v * m;
    return u * m;
  };
}
