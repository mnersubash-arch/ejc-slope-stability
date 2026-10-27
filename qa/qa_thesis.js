// Recompute every number reported in the thesis with the shipped web engine (webapp/engine.js).
// Writes qa_thesis.json: { key: value } for compare_thesis.py.
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "webapp", "engine.js"));
const OUT = path.join(__dirname, "out");                 // results (git-ignored)
fs.mkdirSync(OUT, { recursive: true });

const out = {};
const put = (k, v) => { out[k] = v; };
const t0 = Date.now();
const log = (s) => console.error(`[${((Date.now() - t0) / 1000).toFixed(1)} s] ${s}`);
const MAIN = E.inputs();
const WET = E.inputs({ wet: true });
const CT = E.inputs({ xc: 207.28, zc: 586.53, R: 622.08 });
const cache = new Map();
function r3(kw) { const k = "3" + JSON.stringify(kw); if (!cache.has(k)) cache.set(k, E.run3d(E.inputs(kw))); return cache.get(k); }
function r2(kw) { const k = "2" + JSON.stringify(kw); if (!cache.has(k)) cache.set(k, E.run2d(E.inputs(kw))); return cache.get(k); }
const count = (arr, f) => { let n = 0; for (let k = 0; k < arr.length; k++) if (f(arr[k], k)) n++; return n; };

// ---- Table 4.1 and eq. 2.5
const [mb, s, aG] = E.hoekBrownParameters(42, 10, 0);
put("mb", mb); put("s", s); put("a_gsi", aG);

// ---- Figure 4.2: entry and exit of the central-section circle, active counts
{
  const { xc, zc, R, H, beta } = MAIN;
  const g = (x) => E.groundZ(x, H, beta) - E.sphereZ(x, 0, xc, zc, R);
  const bis = (a, b) => { for (let i = 0; i < 200; i++) { const m = 0.5 * (a + b); if ((g(a) > 0) === (g(m) > 0)) a = m; else b = m; } return 0.5 * (a + b); };
  const xEntry = bis(-400, -300), xExit = bis(-100, -1);
  put("entry_x", xEntry); put("entry_z", E.groundZ(xEntry, H, beta));
  put("exit_x", xExit); put("exit_z", E.groundZ(xExit, H, beta));
  put("active_3d_20", E.buildMesh(MAIN, 3, 20).count); put("active_2d_40", E.buildMesh(MAIN, 2, 40).count);
  put("active_3d_10", E.buildMesh(MAIN, 3, 10).count);
}
log("geometry");

// ---- Section 5.1 / Table C.2: baseline (default path) and native v15 path
for (const [tag, kw] of [["dry", {}], ["wet", { wet: true }]]) {
  const d3 = r3(kw), d2 = r2(kw), n3 = r3(Object.assign({ iteration: "native_v15" }, kw)), n2 = r2(Object.assign({ iteration: "native_v15" }, kw));
  put(`fos3_${tag}`, d3.fos); put(`fos2_${tag}`, d2.fos);
  put(`native_fos3_${tag}`, n3.fos); put(`native_iter3_${tag}`, n3.iterations); put(`native_fos2_${tag}`, n2.fos); put(`native_iter2_${tag}`, n2.iterations);
  put(`default_iter3_${tag}`, d3.iterations);
  const chk = (r, name) => r.checks.find((c) => c.name.startsWith(name)).value;
  put(`local_res3_${tag}`, chk(d3, "Max relative local")); put(`native_local_res3_${tag}`, chk(n3, "Max relative local"));
  put(`ejc_res3_${tag}`, chk(d3, "Max |f|")); put(`native_ejc_res3_${tag}`, chk(n3, "Max |f|"));
  // Section 5.8: negative sigma_n and sigma3 (baseline 20x20)
  put(`neg_sn3_${tag}`, count(d3.sn, (v) => v < 0)); put(`neg_s3_${tag}`, count(d3.strength.s3, (v) => v < 0));
  put(`native_neg_sn3_${tag}`, count(n3.sn, (v) => v < 0)); put(`native_neg_s3_${tag}`, count(n3.strength.s3, (v) => v < 0));
  put(`neg_sn2_${tag}`, count(d2.sn, (v) => v < 0));
  // Section 5.4: cohesion share of the 3D resistance sum(c A)/sum(tau A)
  let cA = 0, tA = 0;
  for (let k = 0; k < d3.mesh.count; k++) { cA += d3.c_used[k] * d3.mesh.area[k]; tA += d3.tau_used[k] * d3.mesh.area[k]; }
  put(`coh_share3_${tag}`, cA / tA);
  let cA2 = 0, tA2 = 0;
  for (let k = 0; k < d2.mesh.count; k++) { cA2 += d2.c_used[k] * d2.mesh.area[k]; tA2 += d2.tau_used[k] * d2.mesh.area[k]; }
  put(`coh_share2_${tag}`, cA2 / tA2);
}
// Figure 5.1: instantaneous parameters along the dry central section
{
  const d2 = r2({});
  const phi = Array.from(d2.strength.phi, (v) => v * E.RAD), c = Array.from(d2.c_used);
  put("fig51_phi_min", Math.min(...phi)); put("fig51_phi_max", Math.max(...phi));
  put("fig51_c_min", Math.min(...c)); put("fig51_c_max", Math.max(...c));
  put("fig51_s1_max", Math.max(...d2.strength.s1)); put("fig51_s2_max", Math.max(...d2.strength.s2)); put("fig51_sn_max", Math.max(...d2.sn) / 1000);
  put("projection_limit", r3({}).projection_limit);
}
// sigma2 = sigma3 closure against GHB: eight matched cases
{
  let worst = 0;
  for (const wet of [false, true]) for (const a of [0.5, 0.509923]) for (const dim of [3, 2]) {
    const run = dim === 3 ? r3 : r2;
    const e = run({ wet, a_override: a, model: "ejc_sigma2_eq_sigma3" }).fos, h = run({ wet, a_override: a, model: "ghb" }).fos;
    worst = Math.max(worst, Math.abs(e - h) / h);
  }
  put("s2s3_vs_ghb_worst_rel", worst);
}
log("baseline");

// ---- Section 5.1: Carranza-Torres circle
put("ct_ghb_2d", r2({ xc: 207.28, zc: 586.53, R: 622.08, model: "ghb" }).fos);
put("ct_ghb_3d", r3({ xc: 207.28, zc: 586.53, R: 622.08, model: "ghb" }).fos);

// ---- local states: Table 5.1 (Route 2), Appendix B.5 and A.5 fixtures
{
  const st = E.strengthOf(MAIN);
  const r = E.solveStrength(Float64Array.from([50, 400, 1600]), st);
  [50, 400, 1600].forEach((sn, k) => { put(`t51_phi_${sn}`, r.phi[k] * E.RAD); put(`t51_c_${sn}`, r.c[k]); put(`t51_tau_${sn}`, r.tau[k]); });
  // B.5: sigma_n = 0.400 MPa
  put("b5_s1", r.s1[1]); put("b5_s2", r.s2[1]); put("b5_s3", r.s3[1]); put("b5_M", r.M[1]); put("b5_phi", r.phi[1] * E.RAD);
  put("b5_tau", r.tau[1] / 1000); put("b5_c", r.c[1] / 1000); put("b5_ejc_res", Math.abs(r.f[1]));
  // A.5: original closure, sigma_ci 120, m_b 3.5, s 0.002, a 0.505, nu 0.05, sigma_n 5.7403 MPa
  const so = E.makeStrength(120, 3.5, 0.002, 0.505, "ejc_original", 0.05);
  const o = E.solveStrength(Float64Array.from([5740.3]), so);
  put("a5_s1", o.s1[0]); put("a5_s2", o.s2[0]); put("a5_s3", o.s3[0]); put("a5_M", o.M[0]); put("a5_phi", o.phi[0] * E.RAD);
  put("a5_tau", o.tau[0] / 1000); put("a5_c", o.c[0] / 1000); put("a5_ejc_res", Math.abs(o.f[0]));
  // the same strength with the Moscarda closure has no admissible state there
  put("a5_moscarda_code", E.REASONS[E.solveStrength(Float64Array.from([5740.3]), E.makeStrength(120, 3.5, 0.002, 0.505, "ejc_moscarda", 0.05)).code[0]]);
}
// Figure 2.2: Mohr-Coulomb circle, c' = 20 MPa, phi' = 30 deg, tangency at sigma_n = 66.821 MPa
{
  const st = E.makeStrength(0, 0, 0, 0, "mc", 0.271, 20000, 30);
  const r = E.solveStrength(Float64Array.from([66820.5]), st);
  put("fig22_s1", r.s1[0]); put("fig22_s3", r.s3[0]); put("fig22_tau", r.tau[0] / 1000); put("fig22_centre", (r.s1[0] + r.s3[0]) / 2);
}
// Figure 2.5: Lode parameter b and the EJC/GHB principal-difference gain along the Moscarda closure
{
  const st = E.strengthOf(MAIN);
  // closure state at a given mean stress p
  const stateAtP = (p) => {
    let lo = st.s3t * 0.999999, hi = 50;
    for (let i = 0; i < 200; i++) {
      const m = 0.5 * (lo + hi), ps = E.principalState(m, st);
      const pm = (ps.s1 + ps.s2 + ps.s3) / 3;
      if (ps.code !== E.OK || pm < p) lo = m; else hi = m;
    }
    return E.principalState(0.5 * (lo + hi), st);
  };
  for (const p of [1.810, 2.906]) {
    const ps = stateAtP(p);
    const b = (ps.s2 - ps.s3) / (ps.s1 - ps.s3);
    // deviatoric radius of the EJC state and of the GHB surface at the same p and Lode angle
    const dev = [ps.s1 - p, ps.s2 - p, ps.s3 - p];
    const r = Math.hypot(...dev), u = dev.map((v) => v / r);
    const ghbF = (rr) => { const s1 = p + rr * u[0], s3 = p + rr * u[2]; const Q = st.mb * s3 / st.sigma_ci + st.s; return Q <= 0 ? 1 : (s1 - s3) - st.sigma_ci * Math.pow(Q, st.a); };
    let lo = 0, hi = r * 2;
    for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); if (ghbF(m) < 0) lo = m; else hi = m; }
    put(`fig25_b_${p}`, b); put(`fig25_radial_gain_${p}`, r / lo - 1);
    put(`fig25_gain_fixed_s3_${p}`, 1 / Math.sqrt(1 - b + b * b) - 1);
  }
  put("closure_gain_nu", Math.pow(1 - 0.271 + 0.271 * 0.271, -0.5) - 1);
}
log("local states");

// ---- Section 5.3 / Table 5.2: original closure
{
  const res = {};
  for (const nu of [0.05, 0.10, 0.15, 0.20, 0.25, 0.271, 0.30, 0.35, 0.40, 0.45, 0.49]) {
    const m = E.originalClosureMinimum(77.7, mb, s, 0.5, nu);
    res[nu] = m;
    put(`t52_snmin_${nu}`, m.sn_min); put(`t52_s3min_${nu}`, m.sigma3_min);
  }
  const md = E.buildMesh(MAIN, 3, 20), mw = E.buildMesh(WET, 3, 20);
  const cd = E.gravityCap(md, 27), cw = E.gravityCap(mw, 27);
  put("cap_dry", Math.max(...cd) / 1000); put("cap_wet", Math.max(...cw) / 1000);
  for (const nu of [0.05, 0.10, 0.20, 0.271, 0.35, 0.45, 0.49]) {
    put(`t52_dry_below_${nu}`, count(cd, (v) => v < res[nu].sn_min * 1000));
    put(`t52_wet_below_${nu}`, count(cw, (v) => v < res[nu].sn_min * 1000));
  }
  put("orig_ratio_to_dry_cap", res[0.271].sn_min / (Math.max(...cd) / 1000));
  // Figure 5.2: end of the admissible original-closure branch, where M = 1
  const so = E.makeStrength(77.7, mb, s, 0.5, "ejc_original", 0.271);
  const ok = (x) => { const ps = E.principalState(x, so); return ps.code === E.OK; };
  let lo = 34.376, hi = 1000;           // admissible at the minimum, not at 1000 MPa
  for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); if (ok(m)) lo = m; else hi = m; }
  put("orig_branch_end_s3", lo); put("orig_branch_end_M", E.principalState(lo, so).M);
}
log("original closure");

// ---- Tables 5.3 and 5.4
for (const wet of [false, true]) for (const a of [0.5, 0.509923]) {
  const tag = `${wet ? "wet" : "dry"}_${a}`;
  put(`t53_hb3_${tag}`, r3({ wet, a_override: a, model: "ghb" }).fos); put(`t53_hb2_${tag}`, r2({ wet, a_override: a, model: "ghb" }).fos);
  put(`t53_ejc3_${tag}`, r3({ wet, a_override: a }).fos); put(`t53_ejc2_${tag}`, r2({ wet, a_override: a }).fos);
}
// tangent shear-strength gain EJC over GHB at the same sigma_n, over the converged column states
{
  const ste = E.strengthOf(MAIN), sth = E.strengthOf(E.inputs({ model: "ghb" }));
  let lo = Infinity, hi = -Infinity;
  for (const kw of [{}, { wet: true }]) {
    const r = r3(kw);
    const sn = Float64Array.from(r.sn).filter((v) => v > 0);
    const te = E.solveStrength(sn, ste), th = E.solveStrength(sn, sth);
    for (let k = 0; k < sn.length; k++) { const gq = te.tau[k] / th.tau[k] - 1; lo = Math.min(lo, gq); hi = Math.max(hi, gq); }
  }
  put("tau_gain_min", lo); put("tau_gain_max", hi);
}
log("tables 5.3-5.4");

// ---- Table 5.5: F_H and F_M at 20x20 and 240x240
for (const n of [20, 240]) for (const wet of [false, true]) {
  const tag = `${n}_${wet ? "wet" : "dry"}`;
  put(`t55_fh_${tag}`, r3({ n3: n, wet }).fos); put(`t55_fm_${tag}`, r3({ n3: n, wet, form: "FM" }).fos);
  log(`F_M ${tag}`);
}

// ---- Table 5.6 and Figure 5.5: mesh series
const N3 = [10, 15, 20, 30, 40, 60, 80, 120, 160, 240], N2 = [20, 40, 80, 160, 320, 640, 1280, 2560, 5120];
for (const a of [0.5, 0.509923]) for (const wet of [false, true]) {
  const tag = `${a}_${wet ? "wet" : "dry"}`;
  put(`mesh3_${tag}`, N3.map((n) => r3({ n3: n, wet, a_override: a }).fos));
  put(`mesh2_${tag}`, N2.map((n) => r2({ n2: n, wet, a_override: a }).fos));
  log(`mesh ${tag}`);
}
// Section 5.8 at 240x240 and 80x80
for (const n of [80, 240]) for (const wet of [false, true]) {
  const r = r3({ n3: n, wet });
  put(`neg_sn_frac_${n}_${wet ? "wet" : "dry"}`, count(r.sn, (v) => v < 0) / r.mesh.count);
  put(`neg_sn_count_${n}_${wet ? "wet" : "dry"}`, count(r.sn, (v) => v < 0));
  put(`active_${n}`, r.mesh.count);
}
log("tension counts");

// ---- Table 5.7: tension sensitivity (20x20, 80x80, 2D 40)
for (const [tag, run, kw] of [["3d20", r3, {}], ["3d80", r3, { n3: 80 }], ["2d40", r2, {}]]) {
  for (const wet of [false, true]) for (const tension of ["none", "C0", "T0"]) {
    const r = run(Object.assign({ wet, tension }, kw));
    put(`t57_${tag}_${wet ? "wet" : "dry"}_${tension}`, r.fos === null ? `${r.status}: ${r.reason.slice(0, 90)}` : r.fos);
    if (tension === "C0") put(`t57_${tag}_${wet ? "wet" : "dry"}_C0_set`, count(r.tensile || [], (v) => v));
  }
  log(`table 5.7 ${tag}`);
}

// ---- Figure 5.6: nu sweep (a = 0.5 and 0.509923), ends and the reference value
for (const a of [0.5, 0.509923]) for (const wet of [false, true]) {
  const nus = [0.05, 0.1, 0.15, 0.2, 0.25, 0.271, 0.3, 0.35, 0.4, 0.45, 0.49];
  put(`nu_${a}_${wet ? "wet" : "dry"}_3d`, nus.map((nu) => r3({ nu, wet, a_override: a }).fos));
  put(`nu_${a}_${wet ? "wet" : "dry"}_2d`, nus.map((nu) => r2({ nu, wet, a_override: a }).fos));
}
log("nu sweep");

fs.writeFileSync(path.join(OUT, "qa_thesis.json"), JSON.stringify(out, null, 1));
log("done");
