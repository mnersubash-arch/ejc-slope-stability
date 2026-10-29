/* EJC slope-stability engine: a line-for-line JavaScript port of ejc_app/engine.py (mnersubash-arch/ejc-slope-stability).
 *
 * Hungr's simplified 3D method of columns and Bishop's 2D method on a prescribed spherical slip surface, with the
 * Extended Jiang Criterion (Route 2), GHB, Mohr-Coulomb, the thesis Section 4.6 projection, the native v15 path,
 * C0/T0 sensitivities, F_M, r_u pore pressure and base circles. Units: strength in MPa, column forces in kN,
 * stresses in kPa, areas in m^2. Equation numbers refer to the submitted thesis.
 *
 * Verified against the Python engine (see verify_engine.js); no DOM access, so it also runs under Node.
 */
(function (root) {
  "use strict";

  const GAMMA_W = 9.81;
  const TOL_EJC = 1e-8, TOL_CLOSURE = 1e-8, TOL_DF = 1e-12, TOL_DSN = 1e-9, TOL_LOCAL = 1e-9, TOL_GLOBAL = 1e-8;
  const MAX_ITER = 2000, N_SCAN = 4000, SCAN_DECADES = 14, N_BISECT = 200;
  const DEG = Math.PI / 180, RAD = 180 / Math.PI;      // as math.radians / math.degrees

  const MODELS = {
    ejc_moscarda: "EJC, Moscarda closure σ2 = ν(σ1 + σ3) (adopted)",
    ejc_original: "EJC, original closure σ2 = σ3/ν − σ1 (admissibility check only)",
    ghb: "Conventional generalised Hoek-Brown (σ2 = σ3)",
    ejc_sigma2_eq_sigma3: "EJC with the closure σ2 = σ3 (must reproduce GHB)",
    mc: "Mohr-Coulomb, constant c′ and φ′",
  };
  const OK = 0, DOMAIN = 1, EJC_RES = 2, CLOSURE_RES = 3, ORDERING = 4, GRADIENT = 5, SHEAR = 6, NO_STATE = 7,
    BISECTION = 8;
  const REASONS = {
    [OK]: "ACCEPTED",
    [DOMAIN]: "REJECTED: DOMAIN (Q <= 0, complex root or non-finite state)",
    [EJC_RES]: "REJECTED: EJC residual |f| >= 1e-8 MPa",
    [CLOSURE_RES]: "REJECTED: closure residual >= 1e-8 MPa",
    [ORDERING]: "REJECTED: ORDERING (sigma1 >= sigma2 >= sigma3 fails)",
    [GRADIENT]: "REJECTED: GRADIENT (M < 1 or not finite)",
    [SHEAR]: "REJECTED: SHEAR (tau_s < 0)",
    [NO_STATE]: "REJECTED: no admissible state at this sigma_n (outside the admissible range of the closure)",
    [BISECTION]: "REJECTED: bisection left the admissible interval",
  };
  const ITERATIONS = {
    default: "Default (brief §1.4): F from F_H/F_M at the current strength, then σ_n",
    native_v15: "Native v15 iteration path (start F = 1, Hungr explicit form, projection)",
  };

  // ---------------------------------------------------------------------------------------------------------------
  // numerics helpers (numpy's pairwise summation, linspace, searchsorted)
  // ---------------------------------------------------------------------------------------------------------------
  function pairwise(a, i0, n) {
    if (n < 8) {
      let res = 0;
      for (let i = 0; i < n; i++) res += a[i0 + i];
      return res;
    }
    if (n <= 128) {
      const r = [a[i0], a[i0 + 1], a[i0 + 2], a[i0 + 3], a[i0 + 4], a[i0 + 5], a[i0 + 6], a[i0 + 7]];
      let i;
      for (i = 8; i < n - (n % 8); i += 8) {
        for (let j = 0; j < 8; j++) r[j] += a[i0 + i + j];
      }
      let res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
      for (; i < n; i++) res += a[i0 + i];
      return res;
    }
    let n2 = Math.floor(n / 2);
    n2 -= n2 % 8;
    return pairwise(a, i0, n2) + pairwise(a, i0 + n2, n - n2);
  }
  const sum = (a) => pairwise(a, 0, a.length);
  const map = (n, f) => { const o = new Float64Array(n); for (let i = 0; i < n; i++) o[i] = f(i); return o; };
  function linspace(start, stop, num) {
    const step = (stop - start) / (num - 1);
    const o = new Float64Array(num);
    for (let k = 0; k < num; k++) o[k] = k * step + start;
    o[num - 1] = stop;
    return o;
  }
  function searchsortedLeft(w, t) {
    let lo = 0, hi = w.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (w[mid] < t) lo = mid + 1; else hi = mid; }
    return lo;
  }
  function maxAbsDiff(a, b) {
    let m = -Infinity;
    for (let i = 0; i < a.length; i++) { const d = Math.abs(a[i] - b[i]); if (Number.isNaN(d)) return NaN; if (d > m) m = d; }
    return m;
  }
  function fmtG(x, p) {                                     // Python's format(x, ".pg")
    if (Number.isNaN(x)) return "nan";
    if (!Number.isFinite(x)) return x > 0 ? "inf" : "-inf";
    if (x === 0) return "0";
    const e = parseInt(x.toExponential(p - 1).split("e")[1], 10);      // exponent after rounding to p digits
    const strip = (t) => (t.includes(".") ? t.replace(/0+$/, "").replace(/\.$/, "") : t);
    if (e < -4 || e >= p) {
      const [mant, ex] = x.toExponential(p - 1).split("e");
      const n = parseInt(ex, 10);
      return `${strip(mant)}e${n < 0 ? "-" : "+"}${String(Math.abs(n)).padStart(2, "0")}`;
    }
    return strip(x.toFixed(Math.max(0, p - 1 - e)));
  }
  function arrMin(a) { let v = Infinity; for (let k = 0; k < a.length; k++) { if (Number.isNaN(a[k])) return NaN; if (a[k] < v) v = a[k]; } return v; }

  // ---------------------------------------------------------------------------------------------------------------
  // eq. 2.5 rock-mass parameters, strength model
  // ---------------------------------------------------------------------------------------------------------------
  function hoekBrownParameters(gsi, mi, d) {
    return [mi * Math.exp((gsi - 100) / (28 - 14 * d)), Math.exp((gsi - 100) / (9 - 3 * d)),
      0.5 + (Math.exp(-gsi / 15) - Math.exp(-20 / 3)) / 6];
  }

  function makeStrength(sigma_ci, mb, s, a, model = "ejc_moscarda", nu = 0.271, c_mc = 0, phi_mc = 0) {
    if (!(model in MODELS)) throw new Error(`unknown strength model ${model}`);
    if (model === "mc") {
      if (!(c_mc >= 0 && phi_mc >= 0 && phi_mc < 90)) throw new Error("require c' >= 0 and 0 <= phi' < 90 deg");
    } else {
      if (!(sigma_ci > 0 && mb > 0 && s > 0 && a > 0 && a <= 1))
        throw new Error("require sigma_ci > 0, m_b > 0, s > 0 and 0 < a <= 1");
      if ((model === "ejc_moscarda" || model === "ejc_original") && !(nu > 0 && nu < 0.5))
        throw new Error("require 0 < nu < 0.5");
    }
    let p = 0, q = 1;                                       // ghb and ejc_sigma2_eq_sigma3
    if (model === "ejc_moscarda") { p = nu; q = nu; } else if (model === "ejc_original") { p = -1; q = 1 / nu; }
    return { sigma_ci, mb, s, a, model, nu, c_mc, phi_mc, p, q, s3t: -s * sigma_ci / mb,
      key: [sigma_ci, mb, s, a, model, nu, c_mc, phi_mc].join("|") };
  }

  /** Squared EJC on sigma2 = p sigma1 + q sigma3 (eqs 4.2, 4.4): P = p^2 - p + 1, R = 2pq - p - q - 1, S = q^2 - q + 1. */
  function closureConstants(p, q) { return [p * p - p + 1, 2 * p * q - p - q - 1, q * q - q + 1]; }

  /** Eq. 2.8 (unsquared EJC) residual; NaN if not real. */
  function ejcResidual(s1, s2, s3, Q, sci, a) {
    const d = s2 - s3;
    const rad = 4 * (sci * sci) * Math.pow(Q, 2 * a) - 3 * (d * d);
    return rad >= 0 ? 2 * s1 - s2 - s3 - Math.sqrt(rad) : NaN;
  }

  /** Principal state, gradient M and tangent parameters at sigma3 [MPa] (eqs 2.6-2.8, 4.1-4.7). */
  function principalState(s3, st) {
    const sci = st.sigma_ci, mb = st.mb, s = st.s, a = st.a, p = st.p, q = st.q;
    const Q = mb * s3 / sci + s;
    const qpos = Q > 0;
    const Qp = qpos ? Q : NaN;
    let s1, s2, M, real;
    if (st.model === "ghb") {
      s1 = s3 + sci * Math.pow(Qp, a);
      s2 = s3;
      M = 1 + a * mb * Math.pow(Qp, a - 1);
      real = qpos;
    } else {
      const [cP, cR, cS] = closureConstants(p, q);
      const K = (sci * sci) * Math.pow(Qp, 2 * a);
      const disc = cR * cR * s3 * s3 - 4 * cP * (cS * s3 * s3 - K);
      const sq = disc >= 0 ? Math.sqrt(disc) : NaN;
      const branch = (sign) => {
        const b1 = (-cR * s3 + sign * sq) / (2 * cP);
        const b2 = p * b1 + q * s3;
        const bf = ejcResidual(b1, b2, s3, Qp, sci, a);
        return [b1, qpos && disc >= 0 && Number.isFinite(bf) && Math.abs(bf) < TOL_EJC && b1 >= b2 && b2 >= s3];
      };
      const [s1p, okp] = branch(1), [s1m, okm] = branch(-1);
      s1 = (okp || !okm) ? s1p : s1m;
      s2 = p * s1 + q * s3;
      real = qpos && disc >= 0;
      M = (2 * a * mb * sci * Math.pow(Qp, 2 * a - 1) - cR * s1 - 2 * cS * s3) / (2 * cP * s1 + cR * s3);
    }
    const f = ejcResidual(s1, s2, s3, Qp, sci, a);
    const clo = s2 - (p * s1 + q * s3);
    const phi = Math.asin((M - 1) / (M + 1));
    const sn = s3 + (s1 - s3) / (M + 1);
    const tau = (s1 - s3) * Math.sqrt(M) / (M + 1);
    const c = tau - sn * Math.tan(phi);
    let code = OK;
    if (!(real && Number.isFinite(f) && Number.isFinite(s1))) code = DOMAIN;
    else if (!(Math.abs(f) < TOL_EJC)) code = EJC_RES;
    else if (!(Math.abs(clo) < TOL_CLOSURE)) code = CLOSURE_RES;
    else if (!(s1 >= s2 && s2 >= s3)) code = ORDERING;
    else if (!(Number.isFinite(M) && M >= 1)) code = GRADIENT;
    else if (!(tau >= 0)) code = SHEAR;
    return { s1, s2, s3, Q, M, f, clo, sn, tau, phi, c, code };
  }

  // ---------------------------------------------------------------------------------------------------------------
  // Route 2: scan sigma3 for the admissible interval, bracket sigma_n(sigma3) - target and bisect
  // ---------------------------------------------------------------------------------------------------------------
  function monotonePieces(ok, v) {
    const pieces = [], n = v.length;
    let i = 0;
    while (i < n - 1) {
      if (!(ok[i] && ok[i + 1])) { i++; continue; }
      const i0 = i, dir = Math.sign(v[i + 1] - v[i]);
      while (i < n - 1 && ok[i + 1]) {
        const sg = Math.sign(v[i + 1] - v[i]);
        if (!(sg === dir || sg === 0)) break;
        i++;
      }
      pieces.push([i0, i]);
    }
    return pieces;
  }

  function scanSigma3(st, top, nScan = N_SCAN) {
    const s3t = st.s3t, expo = linspace(-SCAN_DECADES, 0, nScan + 1);
    let grid = map(nScan + 1, (k) => s3t + (top - s3t) * Math.pow(10, expo[k]));
    let states = Array.from(grid, (x) => principalState(x, st));
    const edges = [];
    for (let k = 0; k < nScan; k++) {
      const okA = states[k].code === OK, okB = states[k + 1].code === OK;
      if (okA === okB) continue;
      let a = grid[k], b = grid[k + 1];
      for (let it = 0; it < 80; it++) {
        const mid = 0.5 * (a + b);
        if (!(mid > Math.min(a, b) && mid < Math.max(a, b))) break;
        if ((principalState(mid, st).code === OK) === okA) a = mid; else b = mid;
      }
      edges.push(okA ? a : b);
    }
    if (edges.length) {
      const all = Array.from(grid).concat(edges).sort((x, y) => x - y);
      grid = Float64Array.from(all.filter((x, i) => i === 0 || x !== all[i - 1]));
      states = Array.from(grid, (x) => principalState(x, st));
    }
    return { grid, states, ok: states.map((s) => s.code === OK), sn: Float64Array.from(states, (s) => s.sn) };
  }

  const scanCache = new Map();
  function cachedScan(st, spanExp, nScan) {
    const key = `${st.key}|${spanExp}|${nScan}`;
    if (scanCache.has(key)) { const v = scanCache.get(key); scanCache.delete(key); scanCache.set(key, v); return v; }
    const v = scanSigma3(st, st.s3t + Math.pow(2, spanExp), nScan);
    scanCache.set(key, v);
    if (scanCache.size > 64) scanCache.delete(scanCache.keys().next().value);
    return v;
  }

  function admissibleSnMin(st, top) {
    const sc = scanSigma3(st, top === undefined ? 10 * st.sigma_ci : top);
    let m = Infinity;
    for (let i = 0; i < sc.sn.length; i++) if (sc.ok[i] && sc.sn[i] < m) m = sc.sn[i];
    return Number.isFinite(m) ? m * 1000 : NaN;
  }

  function mohrCoulomb(snKpa, st) {
    const n = snKpa.length, phi0 = st.phi_mc * DEG, tanP = Math.tan(phi0);
    const out = emptyStrength(n);
    for (let i = 0; i < n; i++) {
      const t = snKpa[i], tau = st.c_mc + t * tanP, centre = t + tau * tanP, radius = tau / Math.cos(phi0);
      out.sn_target[i] = t; out.sn[i] = t; out.phi[i] = phi0; out.c[i] = st.c_mc; out.tau[i] = tau;
      out.M[i] = (1 + Math.sin(phi0)) / (1 - Math.sin(phi0));
      out.s1[i] = (centre + radius) / 1000; out.s3[i] = (centre - radius) / 1000; out.s2[i] = out.s3[i];
      out.code[i] = !Number.isFinite(tau) ? DOMAIN : (tau >= 0 ? OK : SHEAR);
      out.n_brackets[i] = 1;
    }
    return out;
  }

  function emptyStrength(n) {
    const nan = () => new Float64Array(n).fill(NaN);
    return { sn_target: new Float64Array(n), sn: nan(), phi: nan(), c: nan(), tau: nan(), M: nan(), s1: nan(), s2: nan(),
      s3: nan(), f: nan(), closure: nan(), code: new Int32Array(n).fill(NO_STATE), n_brackets: new Int32Array(n) };
  }

  /** Route 2 (brief §1.3) for target normal stresses [kPa]. */
  function solveStrength(snKpa, st) {
    if (st.model === "mc") return mohrCoulomb(snKpa, st);
    const n = snKpa.length, out = emptyStrength(n), s3t = st.s3t;
    const t = map(n, (i) => snKpa[i] / 1000);
    const lo = new Float64Array(n).fill(NaN), hi = new Float64Array(n).fill(NaN), glo = new Float64Array(n).fill(NaN);
    const code = new Int32Array(n).fill(NO_STATE);
    const cand = Array.from(t, (v) => Number.isFinite(v) && v > s3t);
    let tmax = -Infinity;
    for (let i = 0; i < n; i++) { out.sn_target[i] = snKpa[i]; if (cand[i] && t[i] > tmax) tmax = t[i]; }
    if (Number.isFinite(tmax)) {
      const sc = cachedScan(st, Math.ceil(Math.log2(tmax - s3t)), N_SCAN);
      for (const [i0, i1] of monotonePieces(sc.ok, sc.sn)) {
        const v = sc.sn.slice(i0, i1 + 1), inc = v[v.length - 1] >= v[0];
        const w = inc ? v : v.slice().reverse();
        for (let idx = 0; idx < n; idx++) {
          if (!(cand[idx] && t[idx] >= w[0] && t[idx] <= w[w.length - 1])) continue;
          const k = Math.min(Math.max(searchsortedLeft(w, t[idx]), 1), w.length - 1);
          const ja = inc ? i0 + k - 1 : i1 - k, jb = inc ? i0 + k : i1 - k + 1;
          out.n_brackets[idx]++;
          if (Number.isNaN(lo[idx])) { lo[idx] = sc.grid[ja]; hi[idx] = sc.grid[jb]; glo[idx] = sc.sn[ja] - t[idx]; }
        }
      }
    }
    for (let i = 0; i < n; i++) {
      if (!Number.isFinite(lo[i])) continue;
      code[i] = OK;
      for (let it = 0; it < N_BISECT; it++) {
        const mid = 0.5 * (lo[i] + hi[i]);
        if (!(mid > lo[i] && mid < hi[i])) break;
        const ms = principalState(mid, st);
        if (ms.code !== OK) { code[i] = BISECTION; break; }
        const gm = ms.sn - t[i];
        if (glo[i] * gm <= 0) hi[i] = mid; else { lo[i] = mid; glo[i] = gm; }
      }
      let s3 = lo[i];
      if (code[i] === OK) {
        const eLo = principalState(lo[i], st), eHi = principalState(hi[i], st);
        if (Math.abs(eHi.sn - t[i]) < Math.abs(eLo.sn - t[i])) s3 = hi[i];
      }
      const fs = principalState(s3, st);
      out.code[i] = code[i] === OK ? fs.code : code[i];
      out.sn[i] = fs.sn * 1000; out.phi[i] = fs.phi; out.c[i] = fs.c * 1000; out.tau[i] = fs.tau * 1000; out.M[i] = fs.M;
      out.s1[i] = fs.s1; out.s2[i] = fs.s2; out.s3[i] = fs.s3; out.f[i] = fs.f; out.closure[i] = fs.clo;
    }
    return out;
  }

  // ---------------------------------------------------------------------------------------------------------------
  // original closure (eqs 4.1, 4.2): minimum admissible sigma_n (Table 5.2); never a FoS
  // ---------------------------------------------------------------------------------------------------------------
  function originalClosureMinimum(sigma_ci, mb, s, a, nu, lo = 1e-6, hi = 1e6, n = 200000) {
    const st = makeStrength(sigma_ci, mb, s, a, "ejc_original", nu);
    const evaluate = (x) => {
      const ps = principalState(x, st);
      const ok = ps.code === OK && ps.s1 - ps.s3 >= sigma_ci * Math.pow(ps.Q > 0 ? ps.Q : NaN, a) &&
        ps.s1 >= ps.s3 / (2 * nu) && ps.s1 <= ps.s3 * (1 - nu) / nu;
      return [ps, ok];
    };
    const expo = linspace(Math.log10(lo), Math.log10(hi), n + 1);
    let best = null, k = -1;
    const grid = map(n + 1, (i) => Math.pow(10, expo[i]));
    for (let i = 0; i <= n; i++) {
      const [ps, ok] = evaluate(grid[i]);
      if (ok && (best === null || ps.sn < best[0])) { best = [ps.sn, grid[i], ps.phi]; k = i; }
    }
    if (best === null) return { nu, a, sigma3_min: NaN, sn_min: NaN, phi_deg: NaN, found: false };
    let a_ = grid[Math.max(k - 1, 0)], b_ = grid[Math.min(k + 1, n)];
    const gr = (Math.sqrt(5) - 1) / 2;
    const fval = (x) => { const [ps, ok] = evaluate(x); return ok ? ps.sn : Infinity; };
    for (let it = 0; it < 200; it++) {
      const c1 = b_ - gr * (b_ - a_), c2 = a_ + gr * (b_ - a_);
      if (fval(c1) < fval(c2)) b_ = c2; else a_ = c1;
    }
    // The minimum can lie on the edge of the admissible set (sigma2 = sigma3), where the bracket closes on the edge and
    // its midpoint can fall just outside it: keep the best admissible point among the midpoint and both bracket ends.
    for (const x of [0.5 * (a_ + b_), a_, b_]) {
      const [ps, ok] = evaluate(x);
      if (ok && ps.sn <= best[0]) best = [ps.sn, x, ps.phi];
    }
    return { nu, a, sigma3_min: best[1], sn_min: best[0], phi_deg: best[2] * RAD, found: true };
  }

  // ---------------------------------------------------------------------------------------------------------------
  // inputs
  // ---------------------------------------------------------------------------------------------------------------
  const DEFAULT_INPUTS = Object.freeze({
    H: 360, beta: 50, xc: 355.7423115, zc: 644.4192949, R: 725.2669975, gamma: 27, sigma_ci: 77.7, mi: 10, gsi: 42,
    D: 0, a_override: 0.5, nu: 0.271, model: "ejc_moscarda", wet: false, Hw: 360, beta_w: 50, n3: 20, n2: 40,
    form: "FH", tension: "none", iteration: "default", safeguard: "project", c_mc: 0, phi_mc: 0, ru: null,
    base_circle: false, n3y: null,
  });
  const inputs = (kw = {}) => Object.freeze(Object.assign({}, DEFAULT_INPUTS, kw));
  function hb(inp) {
    const [mb, s, a] = hoekBrownParameters(inp.gsi, inp.mi, inp.D);
    return [mb, s, inp.a_override === null || inp.a_override === undefined ? a : Number(inp.a_override)];
  }
  function strengthOf(inp) {
    const [mb, s, a] = hb(inp);
    return makeStrength(inp.sigma_ci, mb, s, a, inp.model, inp.nu, inp.c_mc, inp.phi_mc);
  }

  // ---------------------------------------------------------------------------------------------------------------
  // mesh and column quantities (eqs 4.10, 4.11)
  // ---------------------------------------------------------------------------------------------------------------
  const groundZ = (x, H, beta) => (x <= 0 ? Math.min(H, -x * Math.tan(beta * DEG)) : 0);
  const waterZ = (x, Hw, bw) => (x <= 0 ? Math.min(Hw, -x * Math.tan(bw * DEG)) : 0);
  function sphereZ(x, y, xc, zc, R) {
    const dx = x - xc, q = R * R - dx * dx - y * y;
    return q >= 0 ? zc - Math.sqrt(q) : NaN;
  }

  /** Entry and exit of the slip circle on y = 0 (ends of the interval where the ground is above the circle). */
  function slipExtent(inp, n = 20001) {
    const { xc, zc, R } = inp;
    const above = (x) => groundZ(x, inp.H, inp.beta) - sphereZ(x, 0, xc, zc, R) > 0;
    const xs = linspace(xc - R, xc + R, n);
    let first = -1, last = -1;
    for (let i = 0; i < n; i++) if (above(xs[i])) { if (first < 0) first = i; last = i; }
    if (first < 0) throw new Error("the slip circle does not cut the ground: there is no sliding mass");
    const refine = (a, b) => {
      for (let it = 0; it < 200; it++) {
        const mid = 0.5 * (a + b);
        if (mid === a || mid === b) break;
        if (above(mid)) b = mid; else a = mid;
      }
      return 0.5 * (a + b);
    };
    return [first === 0 ? xs[0] : refine(xs[first - 1], xs[first]),
      last === n - 1 ? xs[n - 1] : refine(xs[last + 1], xs[last])];
  }

  /** Heights of the vertical faces that bound the sliding mass where the slip circle on y = 0 ends below the ground:
   *  at its back extreme x_c - R (z = z_c), and with the base-circle grid also at its front extreme x_c + R (the
   *  default grid stops at the toe). No force acts on such a face: at the back it works as an open tension crack.
   *  Heights below 1e-9 H (rounding) count as zero. Same as vertical_faces in ejc_app/engine.py. */
  function verticalFaces(inp) {
    const { xc, zc, R, H, beta } = inp, tol = 1e-9 * H;
    const back = inp.base_circle || xc - R <= 0 ? groundZ(xc - R, H, beta) - zc : 0;
    const front = inp.base_circle ? groundZ(xc + R, H, beta) - zc : 0;
    return { back: back > tol ? back : 0, front: front > tol ? front : 0 };
  }

  /** The grid's x range on the central section, [x0, x1] with x0 > x1: brief §1.4, x0 = x_c - sqrt(R^2 - z_c^2) and
   *  x1 = x_c - R; with the base-circle option, the circle's exit and entry (slipExtent). Throws for a circle the
   *  default grid refuses. Same as grid_x_range in ejc_app/engine.py. */
  function gridXRange(inp) {
    const { xc, zc, R } = inp;
    if (inp.base_circle) {
      const [x1, x0] = slipExtent(inp);
      return [x0, x1];
    }
    if (!(R > zc)) throw new Error("require R > z_c (the sphere must cut the toe level); for a face circle, with its lowest point " +
      "above the toe level, use the base-circle option");
    const x0 = xc - Math.sqrt(R * R - zc * zc);
    if (x0 < -1e-5 * R)
      throw new Error(`the circle crosses the toe level ${fmtG(-x0, 4)} m behind the toe, so the default grid would ` +
        "cut off the toe region; use the base-circle option");
    return [x0, xc - R];
  }

  /** 3D columns across y for a mesh with n columns along x: n for the brief's N x N grid (n3y null); otherwise n3y
   *  scaled by n/n3 and rounded half up (at least 1). At n = n3 it is n3y. Same as columns_across in ejc_app/engine.py. */
  function columnsAcross(inp, n) {
    if (inp.n3y === null || inp.n3y === undefined) return n;
    if (!(Number.isInteger(inp.n3y) && inp.n3y >= 1 && Number.isInteger(inp.n3) && inp.n3 >= 1))
      throw new Error("the number of columns across y must be a positive integer");
    return Math.max(1, Math.floor((2 * inp.n3y * n + inp.n3) / (2 * inp.n3)));
  }
  /** The inputs with n columns along x and the columns across y scaled to keep their shape (mesh studies). */
  const withN3 = (inp, n) => inputs(Object.assign({}, inp,
    { n3: n, n3y: inp.n3y === null || inp.n3y === undefined ? null : columnsAcross(inp, n) }));
  /** Plan size of the 3D columns, [dx, dy] in m: dx = (x0 - x1)/N along x, dy = R/N_y across y. */
  function columnSize(inp) {
    const [x0, x1] = gridXRange(inp);
    return [(x0 - x1) / inp.n3, inp.R / columnsAcross(inp, inp.n3)];
  }
  /** Columns along x and across y whose plan size is at most dx by dy (m): N = ceil((x0 - x1)/dx), N_y = ceil(R/dy),
   *  each at least 1 (a relative 1e-9 below a whole number counts as it). Same as counts_for_size. */
  function countsForSize(inp, dx, dy) {
    if (!(dx > 0 && dy > 0)) throw new Error("the column length and width must be positive");
    const [x0, x1] = gridXRange(inp);
    return [Math.max(1, Math.ceil((x0 - x1) / dx * (1 - 1e-9))), Math.max(1, Math.ceil(inp.R / dy * (1 - 1e-9)))];
  }

  function buildMesh(inp, dim, n) {
    if (dim !== 2 && dim !== 3) throw new Error("dim must be 2 or 3");
    if (!(Number.isInteger(n) && n >= 1)) throw new Error("mesh size must be a positive integer");
    if (inp.ru !== null && inp.ru !== undefined) {
      if (inp.wet) throw new Error("choose a water table (wet) or r_u, not both");
      if (!(inp.ru >= 0 && inp.ru < 1)) throw new Error("require 0 <= r_u < 1");
    }
    const { xc, zc, R } = inp;
    const [x0, x1] = gridXRange(inp);
    const dx = (x0 - x1) / n;
    const ny = dim === 3 ? columnsAcross(inp, n) : 1, dy = dim === 3 ? R / ny : 1;
    const cols = { i: [], j: [], x: [], y: [], zb: [], zt: [], h: [], q: [] };
    for (let i = 1; i <= n; i++) {
      for (let j = 1; j <= ny; j++) {
        const x = x0 - (i - 0.5) * dx, y = dim === 3 ? (j - 0.5) * dy : 0;
        const ddx = x - xc, q = R * R - ddx * ddx - y * y;
        if (!(q > 0) || (!inp.base_circle && !(x <= 0))) continue;
        const zb = zc - Math.sqrt(q), zt = groundZ(x, inp.H, inp.beta), h = zt - zb;
        if (!(h > 0)) continue;
        cols.i.push(i); cols.j.push(j); cols.x.push(x); cols.y.push(y); cols.zb.push(zb); cols.zt.push(zt);
        cols.h.push(h); cols.q.push(q);
      }
    }
    const cnt = cols.x.length, f64 = (a) => Float64Array.from(a);
    const x = f64(cols.x), y = f64(cols.y), zb = f64(cols.zb), h = f64(cols.h);
    const ah = new Float64Array(cnt).fill(dx * dy);
    const g = map(cnt, (k) => Math.sqrt(cols.q[k]) / R);
    const rho = map(cnt, (k) => (dim === 3 ? Math.sqrt(R * R - y[k] * y[k]) : R));
    const kap = map(cnt, (k) => (xc - x[k]) / rho[k]);
    let u, zw;
    if (inp.ru !== null && inp.ru !== undefined) {
      zw = new Float64Array(cnt).fill(NaN); u = map(cnt, (k) => inp.ru * inp.gamma * h[k]);
    } else if (inp.wet) {
      zw = map(cnt, (k) => waterZ(x[k], inp.Hw, inp.beta_w)); u = map(cnt, (k) => GAMMA_W * Math.max(0, zw[k] - zb[k]));
    } else {
      zw = new Float64Array(cnt).fill(NaN); u = new Float64Array(cnt);
    }
    const mesh = { dim, n, ny, x0, x1, dx, dy, i: Int32Array.from(cols.i), j: Int32Array.from(cols.j), x, y, zb,
      zt: f64(cols.zt), h, ah, g, area: map(cnt, (k) => ah[k] / g[k]), rho, kap,
      w: map(cnt, (k) => inp.gamma * h[k] * ah[k]), u, zw, count: cnt };
    mesh.volume = (dim === 3 ? 2 : 1) * sum(map(cnt, (k) => h[k] * ah[k]));
    return mesh;
  }

  // ---------------------------------------------------------------------------------------------------------------
  // iteration (eqs 4.12-4.16) and final checks
  // ---------------------------------------------------------------------------------------------------------------
  function sensitivity(c, t, tau, sn, tension, tens) {
    const ce = Float64Array.from(c), te = Float64Array.from(t), taue = Float64Array.from(tau);
    for (let k = 0; k < c.length; k++) {
      if (!tens[k]) continue;
      if (tension === "C0") { ce[k] = 0; taue[k] = sn[k] * t[k]; } else if (tension === "T0") { ce[k] = 0; te[k] = 0; taue[k] = 0; }
    }
    return [ce, te, taue];
  }
  function globalF(m, taue, xc, form) {
    if (form === "FM") return sum(map(m.count, (k) => taue[k] * m.area[k] * m.rho[k])) / sum(map(m.count, (k) => m.w[k] * (xc - m.x[k])));
    return sum(map(m.count, (k) => taue[k] * m.area[k])) / sum(map(m.count, (k) => m.w[k] * m.kap[k]));
  }
  function localSn(m, ce, te, F) {
    const ma = map(m.count, (k) => m.g[k] + m.kap[k] * te[k] / F);
    const sn = map(m.count, (k) => (m.w[k] - m.u[k] * m.ah[k] - ce[k] * m.area[k] * m.kap[k] / F) / (m.area[k] * ma[k]));
    return [sn, ma];
  }

  function projectionLimit(st) {
    if (st.model === "mc") return null;
    const tensile = st.s3t * 1000;
    if (st.model === "ghb" || st.model === "ejc_sigma2_eq_sigma3") return tensile;
    const edge = admissibleSnMin(st);
    return Number.isFinite(edge) ? Math.max(tensile, edge) : tensile;
  }

  function solve(m, st, xc, gamma, form = "FH", tension = "none", iteration = "default", safeguard = "project",
    maxIter = MAX_ITER) {
    if (form !== "FH" && form !== "FM") throw new Error("form must be FH or FM");
    if (!["none", "C0", "T0"].includes(tension)) throw new Error("tension must be none, C0 or T0");
    if (!(iteration in ITERATIONS)) throw new Error("iteration must be default or native_v15");
    if (safeguard !== "none" && safeguard !== "project") throw new Error("safeguard must be none or project");
    if (iteration === "native_v15") {
      if (form !== "FH" || tension !== "none") throw new Error("the native v15 iteration path supports F_H only, without C0/T0");
      safeguard = "project";
    }
    const res = { dim: m.dim, status: "REJECTED", reason: "", fos: null, f_trial: null, iterations: 0, form, tension,
      iteration, mesh: m, strength: null, sn: null, c_used: null, tan_used: null, tau_used: null, tensile: null,
      m_alpha: null, local_res: null, global_res: NaN, bishop_res: NaN, checks: [], history: [], projections: 0,
      rejected_index: null, safeguard, projection_limit: null, projection_value: null };
    if (st.model === "ejc_original") {
      res.status = "NOT APPLICABLE";
      res.reason = "The original closure is an admissibility check only and never produces a FoS.";
      return res;
    }
    if (m.count === 0) { res.reason = "No active columns."; return res; }
    const drive = form === "FM" ? sum(map(m.count, (k) => m.w[k] * (xc - m.x[k]))) : sum(map(m.count, (k) => m.w[k] * m.kap[k]));
    if (!(drive > 0)) { res.reason = "Non-positive driving sum."; return res; }
    let sn = map(m.count, (k) => (gamma * m.h[k] - m.u[k]) * (m.g[k] * m.g[k]));
    const limit = safeguard === "project" ? projectionLimit(st) : null;
    let project = (v) => [v, 0];
    if (limit !== null) {
      const value = limit + Math.max(1e-9, 0.001 * Math.abs(limit));
      res.projection_limit = limit; res.projection_value = value;
      project = (v) => { let c = 0; const o = map(v.length, (k) => (v[k] <= limit ? (c++, value) : v[k])); return [o, c]; };
    }
    [sn, res.projections] = project(sn);
    let fOld = iteration === "native_v15" ? 1.0 : null;
    let tens = new Uint8Array(m.count);
    for (let it = 1; it <= maxIter; it++) {
      const stg = solveStrength(sn, st);
      res.iterations = it;
      const bad = [];
      for (let k = 0; k < m.count; k++) {
        if (stg.code[k] !== OK && !(st.model === "mc" && stg.code[k] === SHEAR)) bad.push(k);
      }
      if (bad.length) {
        res.strength = stg; res.sn = sn; res.rejected_index = bad;
        res.reason = `${bad.length} of ${m.count} bases have no accepted local state at iteration ${it}; first: column ` +
          `${bad[0]} at sigma_n = ${fmtG(sn[bad[0]], 6)} kPa, ${REASONS[stg.code[bad[0]]]}`;
        return res;
      }
      const t = map(m.count, (k) => Math.tan(stg.phi[k]));
      if (tension !== "none") { const nt = new Uint8Array(m.count); for (let k = 0; k < m.count; k++) nt[k] = tens[k] || sn[k] < 0 ? 1 : 0; tens = nt; }
      const [ce, te, taue] = sensitivity(stg.c, t, stg.tau, sn, tension, tens);
      let F, nw, ma, nProj;
      if (iteration === "default") {
        F = globalF(m, taue, xc, form);
        if (!(Number.isFinite(F) && F > 0)) { res.reason = `Non-positive or non-finite trial F at iteration ${it}.`; return res; }
        [nw, ma] = localSn(m, ce, te, F);
        [nw, nProj] = project(nw);
      } else {
        [nw, ma] = localSn(m, ce, te, fOld);
        [nw, nProj] = project(nw);
        F = globalF(m, map(m.count, (k) => ce[k] + nw[k] * te[k]), xc, form);
        if (!(Number.isFinite(F) && F > 0)) { res.reason = `Non-positive or non-finite trial F at iteration ${it}.`; return res; }
      }
      res.projections += nProj;
      if (ma.some((v) => v <= 0)) { res.reason = `Non-positive m_alpha at iteration ${it}.`; return res; }
      const prev = iteration === "native_v15" ? fOld : (res.history.length ? res.history[res.history.length - 1][1] : null);
      const dF = prev !== null ? Math.abs(F - prev) / F : Infinity;
      const move = maxAbsDiff(nw, sn);
      let ntens = 0; for (let k = 0; k < m.count; k++) ntens += tens[k];
      res.history.push([it, F, dF, move, ntens, nProj]);
      res.f_trial = F;
      sn = nw;
      fOld = F;
      if (dF < TOL_DF && move < TOL_DSN) {
        if (iteration === "native_v15") [sn] = localSn(m, ce, te, F);
        return finalChecks(res, m, st, xc, F, sn, tens, dF, move);
      }
    }
    res.status = "NOT CONVERGED";
    res.reason = `No convergence in ${maxIter} iterations; no accepted FoS.`;
    res.sn = sn;
    return res;
  }

  function finalChecks(res, m, st, xc, F, sn, tensUsed, dF, move) {
    const stg = solveStrength(sn, st), n = m.count;
    const t = map(n, (k) => Math.tan(stg.phi[k]));
    const tens = map(n, (k) => (res.tension !== "none" ? (tensUsed[k] || sn[k] < 0 ? 1 : 0) : tensUsed[k]));
    const [ce, te, taue] = sensitivity(stg.c, t, stg.tau, sn, res.tension, tensUsed);
    const ma = map(n, (k) => m.g[k] + m.kap[k] * te[k] / F);
    const local = map(n, (k) => Math.abs(m.w[k] - (sn[k] + m.u[k]) * m.ah[k] - taue[k] * m.area[k] * m.kap[k] / F) / m.w[k]);
    const glob = Math.abs(globalF(m, taue, xc, res.form) - F) / F;
    const bishop = sum(map(n, (k) => (ce[k] * m.ah[k] + (m.w[k] - m.u[k] * m.ah[k]) * te[k]) / ma[k])) /
      sum(map(n, (k) => m.w[k] * m.kap[k]));
    const bishopRes = Math.abs(bishop - F) / F;
    const acc = Array.from(stg.code, (c) => c === OK);
    const nanMax = (a) => { let v = -Infinity, any = false; for (const x of a) if (!Number.isNaN(x)) { any = true; if (x > v) v = x; } return any ? v : NaN; };
    const nanMin = (a) => { let v = Infinity, any = false; for (const x of a) if (!Number.isNaN(x)) { any = true; if (x < v) v = x; } return any ? v : NaN; };
    const all = (a, f) => { for (let k = 0; k < a.length; k++) if (!f(a[k], k)) return false; return true; };
    const chk = (name, value, limit, passed) => ({ name, value, limit, passed: !!passed });
    const checks = [
      chk("Convergence |dF|/F", dF, "< 1e-12", dF < TOL_DF),
      chk("Convergence max |d sigma_n| (kPa)", move, "< 1e-9", move < TOL_DSN),
      chk("Local states accepted (count rejected)", acc.filter((v) => !v).length, "0", acc.every((v) => v)),
    ];
    if (st.model !== "mc") {
      const absF = Array.from(stg.f, Math.abs), absC = Array.from(stg.closure, Math.abs);
      checks.push(chk("Max |f| unsquared EJC (MPa)", absF.some(Number.isFinite) ? nanMax(absF) : NaN, "< 1e-8",
        all(stg.f, Number.isFinite) && all(absF, (v) => v < TOL_EJC)));
      checks.push(chk("Max closure residual (MPa)", nanMax(absC), "< 1e-8", all(absC, (v) => v < TOL_CLOSURE)));
    }
    const ordered = (k) => stg.s1[k] >= stg.s2[k] && stg.s2[k] >= stg.s3[k];
    let nOrd = 0; for (let k = 0; k < n; k++) if (!ordered(k)) nOrd++;
    checks.push(
      chk("Ordering sigma1 >= sigma2 >= sigma3 (failures)", nOrd, "0", nOrd === 0),
      chk("Min M (gradient)", nanMin(stg.M), ">= 1, finite", all(stg.M, (v) => Number.isFinite(v) && v >= 1)),
      chk("Min tau_s (kPa)", nanMin(stg.tau), ">= 0", all(stg.tau, (v) => v >= 0)),
      chk("Min m_alpha", arrMin(ma), "> 0", all(ma, (v) => v > 0)),
      chk("Max relative local residual", nanMax(local), "<= 1e-9", all(local, (v) => v <= TOL_LOCAL)),
      chk(`Relative global residual (${res.form})`, glob, "<= 1e-8", glob <= TOL_GLOBAL));
    if (m.dim === 2) checks.push(chk("Relative Bishop residual", bishopRes, "<= 1e-8", bishopRes <= TOL_GLOBAL));
    if (res.tension === "C0" || res.tension === "T0") {
      let changes = 0; for (let k = 0; k < n; k++) if (tens[k] !== tensUsed[k]) changes++;
      checks.push(chk("New tensile bases at the final state", changes, "0", changes === 0));
    }
    if (res.projection_limit !== null) {
      const last = res.history.length ? res.history[res.history.length - 1][5] : 0;
      checks.push(chk("Projected bases at convergence", last, "0", last === 0));
      checks.push(chk("Min sigma_n above the projection limit (kPa)", arrMin(sn) - res.projection_limit, "> 0",
        all(sn, (v) => v > res.projection_limit)));
    }
    Object.assign(res, { strength: stg, sn, c_used: ce, tan_used: te, tau_used: taue, tensile: Uint8Array.from(tensUsed),
      m_alpha: ma, local_res: local, global_res: glob, bishop_res: bishopRes, checks, f_trial: F });
    if (checks.every((c) => c.passed)) {
      res.status = "ACCEPTED"; res.fos = F;
      res.reason = "Every local state, local and global equilibrium check passed.";
      if (res.tension !== "none") res.reason += ` ${res.tension}: sensitivity, not a calibrated contact law.`;
    } else {
      res.status = "REJECTED";
      res.reason = "Final-state check failed: " + checks.filter((c) => !c.passed).map((c) => c.name).join("; ");
    }
    return res;
  }

  // ---------------------------------------------------------------------------------------------------------------
  // one case = 3D + 2D; helpers for the pages
  // ---------------------------------------------------------------------------------------------------------------
  const run3d = (inp) => solve(buildMesh(inp, 3, inp.n3), strengthOf(inp), inp.xc, inp.gamma, inp.form, inp.tension,
    inp.iteration, inp.safeguard);
  const run2d = (inp) => solve(buildMesh(inp, 2, inp.n2), strengthOf(inp), inp.xc, inp.gamma, "FH", inp.tension,
    inp.iteration, inp.safeguard);
  function run(inp) {
    const r3 = run3d(inp), r2 = run2d(inp);
    return { inputs: inp, r3, r2, ratio: r3.fos !== null && r2.fos !== null ? r3.fos / r2.fos : null };
  }
  const gravityCap = (m, gamma) => map(m.count, (k) => gamma * m.h[k] - m.u[k]);

  /** Forces on element k at the reported F (thesis eq. 4.12): N' = sigma_n A and U = u A along the base normal n,
   *  S = tau_s A / F along s = (-(z_c - z_b), 0, x_c - x)/rho; vertical balance W = (N' + U) g_z + S kappa. */
  function columnForces(inp, r, k) {
    const m = r.mesh, F = r.fos !== null ? r.fos : r.f_trial;
    if (!r.sn || !r.tau_used || F === null) throw new Error("no converged state for this element");
    const x = m.x[k], y = m.y[k], zb = m.zb[k];
    const nvec = [(inp.xc - x) / inp.R, -y / inp.R, (inp.zc - zb) / inp.R];
    const svec = [-(inp.zc - zb) / m.rho[k], 0, (inp.xc - x) / m.rho[k]];
    const N = r.sn[k] * m.area[k], U = m.u[k] * m.area[k], S = r.tau_used[k] * m.area[k] / F, W = m.w[k];
    const resid = W - (N + U) * m.g[k] - S * m.kap[k];
    const foot = [inp.xc, y, inp.zc];
    return { W, N, U, S, F, n: nvec, s: svec, g: m.g[k], kap: m.kap[k], rho: m.rho[k], area: m.area[k],
      residual: resid, residual_rel: Math.abs(resid) / W, base: [x, y, zb], top: [x, y, m.zt[k]], lever_foot: foot,
      lever: Math.hypot(x - foot[0], y - foot[1], zb - foot[2]) };
  }

  /** All forces on the sliding mass at the reported F (free-body diagram of the whole mass): for every column or slice
   *  W (vertical), N' = sigma_n A and U = u A along the base normal n, S = tau_s A / F along s, as columnForces; their
   *  totals for the full mass (both halves in 3D, per metre in 2D) and their moments about the rotation axis (parallel
   *  to y through (x_c, z_c)), positive in the sense in which the weight drives. N' and U act along n, through the
   *  sphere centre, so their moments are zero; S has lever arm rho. Vertical balance holds element by element
   *  (eq. 4.12); horizontal force balance is not part of the method. The method's own global equation is
   *  sum S = sum W kappa for F_H (eq. 4.14) and sum S rho = sum W (x_c - x) for F_M (eq. 4.15). Same as
   *  ejc_app/engine.py mass_forces. */
  function massForces(inp, r) {
    const m = r.mesh, F = r.fos !== null ? r.fos : r.f_trial;
    if (!r.sn || !r.tau_used || F === null) throw new Error("no converged state for this run");
    const n = m.count, mult = m.dim === 3 ? 2 : 1, xc = inp.xc, zc = inp.zc, R = inp.R;
    const nx = map(n, (k) => (xc - m.x[k]) / R), ny = map(n, (k) => -m.y[k] / R), nz = Float64Array.from(m.g);
    const sx = map(n, (k) => -(zc - m.zb[k]) / m.rho[k]), sz = map(n, (k) => (xc - m.x[k]) / m.rho[k]);
    const W = Float64Array.from(m.w), N = map(n, (k) => r.sn[k] * m.area[k]), U = map(n, (k) => m.u[k] * m.area[k]);
    const S = map(n, (k) => r.tau_used[k] * m.area[k] / F);
    const tot = (f, d) => mult * sum(map(n, (k) => f[k] * d[k]));
    const vert = map(n, (k) => W[k] - (N[k] + U[k]) * nz[k] - S[k] * sz[k]);
    const totals = { W: { x: 0, z: -mult * sum(W) }, N: { x: tot(N, nx), z: tot(N, nz) }, U: { x: tot(U, nx), z: tot(U, nz) },
      S: { x: tot(S, sx), z: tot(S, sz) },
      net: { x: mult * sum(map(n, (k) => (N[k] + U[k]) * nx[k] + S[k] * sx[k])), z: -mult * sum(vert) } };
    const mW = mult * sum(map(n, (k) => W[k] * (xc - m.x[k]))), mS = -mult * sum(map(n, (k) => S[k] * m.rho[k]));
    const fm = r.form === "FM";
    const lhs = fm ? -mS : mult * sum(S), rhs = fm ? mW : mult * sum(map(n, (k) => W[k] * m.kap[k]));
    let vmax = 0;
    for (let k = 0; k < n; k++) vmax = Math.max(vmax, Math.abs(vert[k]) / W[k]);
    return { F, dim: m.dim, mult, count: n, form: r.form, weight: -totals.W.z, totals,
      moments: { W: mW, N: 0, U: 0, S: mS, net: mW + mS },
      method: { form: r.form, lhs, rhs, residual: (lhs - rhs) / rhs }, vertical_max_rel: vmax,
      el: { x: m.x, y: m.y, zb: m.zb, zm: map(n, (k) => (m.zb[k] + m.zt[k]) / 2), W, N, U, S, nx, ny, nz, sx, sz } };
  }

  // ---------------------------------------------------------------------------------------------------------------
  // critical-surface search (QA/QC finding R1)
  // ---------------------------------------------------------------------------------------------------------------
  /** Settings of criticalSearch: limits of the exits (multiples of H behind the crest edge, and ahead of the toe for
   *  base circles), grid of starting surfaces (front exit x back exit x curvature; grid_base with the base-circle
   *  option), Nelder-Mead starts, tolerance on the normalised parameters, evaluations per run, restarts, and the
   *  least share of the default grid's span that the mass must cover. Same as SEARCH in ejc_app/engine.py. */
  const SEARCH = Object.freeze({ back: 2, front: 2, grid: Object.freeze([4, 9, 4]), grid_base: Object.freeze([6, 9, 4]),
    starts: 2, tol: 1e-3, max_eval: 100, restarts: 2, coverage: 0.5 });

  /** Point of the ground surface at s (m): s <= 0 on the pit floor ahead of the toe (x = -s), 0 < s <= H up the face
   *  by height (z = s), s > H on the crest at s - H behind the crest edge. */
  function groundAt(s, H, tanB) {
    if (s <= 0) return [-s, 0];
    if (s <= H) return [-s / tanB, s];
    return [-H / tanB - (s - H), H];
  }

  /** Circle through the front exit pf and the back exit pb (x_b < x_f, z_b >= z_f) whose lower arc joins them. The
   *  centre lies on the chord's perpendicular bisector, d = (c/2)(dz/dx + 1/t - 1) above the chord's midpoint, so that
   *  z_c - z_b = (1/t - 1) dx/2: t = 1 makes the arc vertical at pb, t -> 0 flattens it towards the chord. For
   *  t <= 1 the arc meets the ground at pb before its back extreme x_c - R, so the mass has no vertical back face. */
  function circleFromExits(pf, pb, t) {
    const dx = pf[0] - pb[0], dz = pb[1] - pf[1];
    const c = Math.sqrt(dx * dx + dz * dz), h = c / 2;
    const d = h * (dz / dx + 1 / t - 1);
    return { xc: (pf[0] + pb[0]) / 2 + d * dz / c, zc: (pf[1] + pb[1]) / 2 + d * dx / c, R: Math.sqrt(h * h + d * d) };
  }

  /** The surface of search parameters p = [a, b, t]: front exit at s = aH, back exit at s = bH (groundAt), curvature
   *  t (circleFromExits). */
  function searchSurface(inp, p) {
    const H = inp.H, tanB = Math.tan(inp.beta * DEG);
    const front = groundAt(p[0] * H, H, tanB), back = groundAt(p[1] * H, H, tanB);
    return Object.assign(circleFromExits(front, back, p[2]), { front, back });
  }

  /** True if the ground is above the circle on y = 0 over one interval that ends inside (x_c - R, x_c + R), sampled as
   *  slipExtent samples it: one sliding mass, with no vertical face at the circle's extremes (base-circle grid). */
  function singleMass(inp) {
    const { xc, zc, R } = inp, n = 20001, xs = linspace(xc - R, xc + R, n);
    let runs = 0, prev = false, first = -1, last = -1;
    for (let i = 0; i < n; i++) {
      const on = groundZ(xs[i], inp.H, inp.beta) - sphereZ(xs[i], 0, xc, zc, R) > 0;
      if (on) { if (!prev) runs++; if (first < 0) first = i; last = i; }
      prev = on;
    }
    return runs === 1 && first > 0 && last < n - 1;
  }

  /** Nelder-Mead (reflection 1, expansion 2, contractions 1/2, shrink 1/2) from x0 and x0 + step_i e_i, stopped when
   *  every vertex lies within tol of the best in every parameter or after maxEval calls; ties go to the lower index.
   *  Same as _nelder_mead in ejc_app/engine.py. */
  function nelderMead(f, x0, step, tol, maxEval) {
    const n = x0.length;
    let pts = [x0.slice()];
    for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] = x0[i] + step[i]; pts.push(p); }
    let vals = pts.map(f), evals = n + 1;
    while (evals < maxEval) {
      const idx = vals.map((v, i) => i).sort((a, b) => (vals[a] < vals[b] ? -1 : vals[a] > vals[b] ? 1 : a - b));
      pts = idx.map((i) => pts[i]); vals = idx.map((i) => vals[i]);
      let spread = 0;
      for (let i = 1; i <= n; i++) for (let j = 0; j < n; j++) spread = Math.max(spread, Math.abs(pts[i][j] - pts[0][j]));
      if (spread < tol) break;
      const c = [];
      for (let j = 0; j < n; j++) { let s = 0; for (let i = 0; i < n; i++) s += pts[i][j]; c.push(s / n); }
      const at = (t) => c.map((cj, j) => cj + t * (pts[n][j] - cj));
      const xr = at(-1), fr = f(xr); evals++;
      if (fr < vals[0]) {
        const xe = at(-2), fe = f(xe); evals++;
        if (fe < fr) { pts[n] = xe; vals[n] = fe; } else { pts[n] = xr; vals[n] = fr; }
      } else if (fr < vals[n - 1]) { pts[n] = xr; vals[n] = fr; }
      else {
        const outside = fr < vals[n], xk = at(outside ? -0.5 : 0.5), fk = f(xk); evals++;
        if (outside ? fk <= fr : fk < vals[n]) { pts[n] = xk; vals[n] = fk; }
        else {
          for (let i = 1; i <= n; i++) { pts[i] = pts[i].map((v, j) => pts[0][j] + 0.5 * (v - pts[0][j])); vals[i] = f(pts[i]); evals++; }
        }
      }
    }
    let b = 0;
    for (let i = 1; i <= n; i++) if (vals[i] < vals[b]) b = i;
    return { x: pts[b], f: vals[b] };
  }

  /** Critical-surface search: the lowest accepted FoS over circles (2D, Bishop) and spheres centred on y = 0 (3D,
   *  Hungr) that leave the ground at a front exit (the toe or the face; with the base-circle option also the pit
   *  floor, up to front H ahead of the toe) and a back exit (face or crest, up to back H behind the crest edge), with
   *  the lower arc between them (circleFromExits, so no vertical back face). Each surface is solved with the current
   *  inputs, mesh and grid option; it is skipped when the grid option cannot analyse it (the default grid refuses
   *  it, or its mass covers less than `coverage` of the default grid's span from x_c - R to the toe-level crossing,
   *  so that few slices or columns resolve it; with the base-circle option, the ground is not above it over one
   *  interval), and a result that is not ACCEPTED is withheld. 2D: a grid of starting surfaces (front exit x back
   *  exit x curvature), then Nelder-Mead from the `starts` best grid surfaces that are not neighbours, each restarted
   *  from its minimum with half the initial step until it stops improving. 3D: Nelder-Mead from the 2D minima.
   *  Then each critical surface is probed one step 2 tol along every parameter (to name the limits it lies on) and
   *  re-solved on a finer mesh. A local search: the lowest FoS among the surfaces tried, not a proof of the global
   *  minimum. Same as critical_search in ejc_app/engine.py. */
  function criticalSearch(inp, opts = {}, progress = null) {
    const o = Object.assign({}, SEARCH), dims = opts.dims || [2];
    for (const [k, v] of Object.entries(opts)) if (v !== undefined && v !== null) o[k] = v;
    const t0 = Date.now(), H = inp.H, base = !!inp.base_circle;
    const [na, nb, nt] = base ? o.grid_base : o.grid;
    const lo0 = base ? -o.front : 0, hi1 = 1 + o.back;
    const sp = [(1 - lo0) / na, hi1 / nb, 1 / nt], step = sp.map((v) => v / 2);
    const local = o.max_eval * (1 + o.restarts);
    const budget = { 2: na * nb * nt + o.starts * local + 8, 3: o.starts * local + 8 };
    const total = dims.reduce((s, d) => s + budget[d], 0);
    let done = 0;
    const say = (label) => { if (progress) progress({ done: Math.min(done, total), total, label }); };
    const faces = verticalFaces(inp);
    const out = { settings: { back: o.back, front: o.front, grid: [na, nb, nt], starts: o.starts, tol: o.tol,
      max_eval: o.max_eval, restarts: o.restarts, coverage: o.coverage, dims: dims.slice(), base_circle: base },
    domain: { front: [lo0, 1], back: [0, hi1], curvature: [0, 1] },
    prescribed: { xc: inp.xc, zc: inp.zc, R: inp.R, vertical_back: faces.back, vertical_front: faces.front }, r2: null,
      r3: null, ms: 0 };
    const solveAt = (dim, kw) => {
      try {
        const r = (dim === 3 ? run3d : run2d)(inputs(Object.assign({}, inp, kw)));
        return { fos: r.fos, status: r.status, reason: r.fos !== null ? "" : r.reason };
      } catch (e) { return { fos: null, status: "INPUT ERROR", reason: String((e && e.message) || e) }; }
    };
    for (const dim of dims) {
      const pr = solveAt(dim, {});
      out.prescribed["fos" + dim] = pr.fos; out.prescribed["status" + dim] = pr.status;
      out.prescribed["reason" + dim] = pr.reason;
    }
    const stageEnd = { 2: budget[2], 3: (dims.includes(2) ? budget[2] : 0) + budget[3] };
    let seeds = null;
    for (const dim of dims) {
      const td = Date.now(), cache = new Map();
      const counts = { solved: 0, withheld: 0, skipped: 0 }, why = { coverage: 0, masses: 0, grid: 0 };
      const cloud = { xc: [], zc: [], R: [], fos: [] };
      let stage = "";
      const evaluate = (p) => {
        const [a, b, t] = p;
        if (!(a >= lo0 && a < 1 && b > Math.max(a, 0) && b <= hi1 && t > 0 && t <= 1)) return { v: Infinity, kind: "outside" };
        const s = searchSurface(inp, p), kw = { xc: s.xc, zc: s.zc, R: s.R };
        let skip = null;
        if (base) { if (!singleMass(inputs(Object.assign({}, inp, kw)))) skip = "masses"; }
        else if (s.R > s.zc) {
          const span = (s.xc - Math.sqrt(s.R * s.R - s.zc * s.zc)) - (s.xc - s.R);
          if (s.front[0] - s.back[0] < o.coverage * span) skip = "coverage";
        }
        let e;
        if (skip) e = { v: Infinity, kind: "skipped", why: skip };
        else {
          const r = solveAt(dim, kw);
          if (r.status === "INPUT ERROR") e = { v: Infinity, kind: "skipped", why: "grid" };
          else if (r.fos === null) e = { v: Infinity, kind: "withheld" };
          else e = { v: r.fos, kind: "solved" };
        }
        counts[e.kind]++;
        if (e.kind === "skipped") why[e.why]++;
        cloud.xc.push(s.xc); cloud.zc.push(s.zc); cloud.R.push(s.R); cloud.fos.push(e.kind === "solved" ? e.v : null);
        done++;
        if (progress && (done % 5 === 0)) say(stage);
        return e;
      };
      const look = (p) => {
        const key = p.join(",");
        let e = cache.get(key);
        if (e === undefined) { e = evaluate(p); cache.set(key, e); }
        return e;
      };
      const f = (p) => look(p).v;
      const runs = [];
      const res = { dim, fos: null, message: "", runs, counts, skipped_why: why, evaluations: 0, edges: [], check: null,
        cloud };
      let starts = [];
      if (dim === 2) {
        stage = "2D: grid of starting surfaces";
        const cells = [];
        for (let i = 0; i < na; i++) for (let j = 1; j <= nb; j++) for (let k = 0; k < nt; k++) {
          const p = [lo0 + sp[0] * i, sp[1] * j, sp[2] * (k + 0.5)];
          cells.push({ p, ijk: [i, j, k], v: f(p) });
        }
        const order = cells.map((c, i) => i).sort((a, b) => (cells[a].v < cells[b].v ? -1 : cells[a].v > cells[b].v ? 1 : a - b));
        for (const i of order) {
          if (!(cells[i].v < Infinity) || starts.length >= o.starts) break;
          const ijk = cells[i].ijk;
          if (starts.some((c) => Math.max(Math.abs(c.ijk[0] - ijk[0]), Math.abs(c.ijk[1] - ijk[1]), Math.abs(c.ijk[2] - ijk[2])) < 2)) continue;
          starts.push({ p: cells[i].p, ijk, v: cells[i].v });
        }
        done = Math.max(done, na * nb * nt);
        if (!starts.length) res.message = "No surface of the starting grid gave an accepted FoS.";
      } else if (seeds && seeds.length) starts = seeds.map((p) => ({ p, v: f(p) }));
      else res.message = "The 2D search found no critical circle to start the 3D search from.";
      let best = null;
      starts.forEach((c, q) => {
        stage = `${dim}D: Nelder–Mead from start ${q + 1} of ${starts.length}`;
        say(stage);
        const from = done;
        let r = nelderMead(f, c.p, step, o.tol, o.max_eval);
        for (let k = 0; k < o.restarts; k++) {
          const r2 = nelderMead(f, r.x, step.map((v) => v / 2), o.tol, o.max_eval);
          if (!(r2.f < r.f)) break;
          r = r2;
        }
        done = Math.max(done, from + local);
        runs.push({ p0: c.p, f0: c.v < Infinity ? c.v : null, p: r.x, f: r.f < Infinity ? r.f : null });
        if (r.f < Infinity && (!best || r.f < best.f)) best = r;
      });
      if (dim === 2) seeds = runs.filter((r) => r.f !== null).map((r) => r.p);
      if (best) {
        stage = `${dim}D: limits and mesh check`;
        say(stage);
        const names = ["front", "back", "curvature"];
        let top = best;
        for (let i = 0; i < 3; i++) for (const sg of [-1, 1]) {
          const q = best.x.slice();
          q[i] = best.x[i] + sg * 2 * o.tol;
          const e = look(q);
          if (e.kind === "solved") { if (e.v < top.f) top = { x: q, f: e.v }; continue; }
          let w = e.kind === "skipped" ? e.why : e.kind;
          if (e.kind === "outside") {
            w = q[0] < lo0 ? (base ? "front limit" : "toe") : q[0] >= 1 ? "crest" : q[1] <= Math.max(q[0], 0) ? "exits meet"
              : q[1] > hi1 ? "back limit" : q[2] > 1 ? "vertical back" : "flat";
          }
          res.edges.push({ param: names[i], side: sg, why: w });
        }
        best = top;
        const s = searchSurface(inp, best.x);
        Object.assign(res, { fos: best.f, p: best.x, xc: s.xc, zc: s.zc, R: s.R, front: s.front, back: s.back });
        if (!base && s.R > s.zc) res.coverage = (s.front[0] - s.back[0]) / ((s.xc - Math.sqrt(s.R * s.R - s.zc * s.zc)) - (s.xc - s.R));
        const n = dim === 3 ? Math.min(2 * inp.n3, 320) : Math.min(4 * inp.n2, 10000);
        const kw = dim === 3 ? { xc: s.xc, zc: s.zc, R: s.R, n3: n } : { xc: s.xc, zc: s.zc, R: s.R, n2: n };
        let ny = null;
        if (dim === 3 && inp.n3y !== null && inp.n3y !== undefined) {   // columns across y in the same ratio, at most 320
          ny = Math.min(columnsAcross(inp, n), 320);
          kw.n3y = ny;
        }
        const ck = solveAt(dim, kw);
        res.check = { n, ny: ny !== null ? ny : (dim === 3 ? n : 1), fos: ck.fos, status: ck.status, reason: ck.reason,
          change: ck.fos !== null ? (ck.fos - best.f) / ck.fos : null };
      }
      res.evaluations = counts.solved + counts.withheld + counts.skipped;
      res.ms = Date.now() - td;
      out["r" + dim] = res;
      done = stageEnd[dim];
      say(`${dim}D: done`);
    }
    out.ms = Date.now() - t0;
    return out;
  }

  const api = { GAMMA_W, MODELS, REASONS, ITERATIONS, OK, SHEAR, NO_STATE, DOMAIN, hoekBrownParameters, makeStrength,
    closureConstants, principalState, scanSigma3, admissibleSnMin, solveStrength, originalClosureMinimum,
    DEFAULT_INPUTS, inputs, hb, strengthOf, groundZ, waterZ, sphereZ, slipExtent, verticalFaces, gridXRange,
    columnsAcross, withN3, columnSize, countsForSize, buildMesh,
    projectionLimit, solve,
    run3d, run2d, run, gravityCap, columnForces, massForces, SEARCH, groundAt, circleFromExits, searchSurface,
    singleMass, nelderMead, criticalSearch, sum, fmtG, DEG, RAD };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.EJC = api;
})(typeof window !== "undefined" ? window : globalThis);
