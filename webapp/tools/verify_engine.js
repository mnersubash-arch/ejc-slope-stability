// Compare engine.js with the Python engine's reference results (reference.json from dump_reference.py).
// Usage (from webapp/): python3 tools/dump_reference.py && node tools/verify_engine.js [tools/ref/reference.json]
// Passes when every analysis matches: FoS and trial F to 1e-12 (relative), the same status, iteration count,
// projection count, mesh and check results, sigma_n and tau_s arrays to 1e-9; strength fixtures to 1e-11.
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "engine.js"));

const ref = JSON.parse(fs.readFileSync(process.argv[2] || path.join(__dirname, "ref", "reference.json"), "utf8"));
const val = (v) => (v === "nan" ? NaN : v === "inf" ? Infinity : v === "-inf" ? -Infinity : v);
let fails = 0, compared = 0, worstFos = 0, worstArr = 0;
const fail = (msg) => { fails++; console.log("FAIL", msg); };

function rel(a, b) {
  a = val(a); b = val(b);
  if (a === null || b === null) return a === b ? 0 : Infinity;
  if (Number.isNaN(a) && Number.isNaN(b)) return 0;
  if (a === b) return 0;
  return Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
}
function arrDiff(js, py, absFloor) {
  if (py === null || py === undefined) return js === null || js === undefined ? 0 : Infinity;
  if (!js || js.length !== py.length) return Infinity;
  let w = 0;
  for (let k = 0; k < py.length; k++) {
    const a = js[k], b = val(py[k]);
    if (Number.isNaN(a) && Number.isNaN(b)) continue;
    const d = Math.abs(a - b) / Math.max(Math.abs(b), absFloor);
    if (!(d <= w)) w = d;                                   // NaN-safe
  }
  return w;
}
const inputsOf = (d) => E.inputs(Object.fromEntries(Object.entries(d).map(([k, v]) => [k, val(v)])));

for (const c of ref.cases) {
  const inp = inputsOf(c.inputs);
  for (const dim of [3, 2]) {
    const py = c[`r${dim}`];
    if (!py) continue;
    compared++;
    const tag = `${c.name} ${dim}D`;
    let js;
    try { js = dim === 3 ? E.run3d(inp) : E.run2d(inp); } catch (e) {
      if (py.error) { if (e.message !== py.error) fail(`${tag}: error text\n  js ${e.message}\n  py ${py.error}`); continue; }
      fail(`${tag}: js threw ${e.message}`); continue;
    }
    if (py.error) { fail(`${tag}: python raised "${py.error}" but js gave ${js.status}`); continue; }
    if (js.status !== py.status) fail(`${tag}: status js ${js.status} py ${py.status} (${js.reason} | ${py.reason})`);
    if (js.reason !== py.reason) console.log(`note ${tag}: reason text differs\n  js ${js.reason}\n  py ${py.reason}`);
    const dF = rel(js.fos, py.fos);
    worstFos = Math.max(worstFos, dF === Infinity ? 0 : dF);
    if (!(dF <= 1e-12)) fail(`${tag}: fos js ${js.fos} py ${py.fos} rel ${dF}`);
    if (rel(js.f_trial, py.f_trial) > 1e-12) fail(`${tag}: f_trial js ${js.f_trial} py ${py.f_trial}`);
    for (const k of ["iterations", "projections"]) if (js[k] !== py[k]) fail(`${tag}: ${k} js ${js[k]} py ${py[k]}`);
    if (js.mesh.count !== py.count) fail(`${tag}: count js ${js.mesh.count} py ${py.count}`);
    if (rel(js.mesh.volume, py.volume) > 1e-13) fail(`${tag}: volume js ${js.mesh.volume} py ${py.volume}`);
    if (rel(js.mesh.x0, py.x0) > 1e-14 || rel(js.mesh.x1, py.x1) > 1e-14) fail(`${tag}: grid ends differ`);
    if (rel(js.projection_limit, py.projection_limit) > 1e-12) fail(`${tag}: projection limit js ${js.projection_limit} py ${py.projection_limit}`);
    if (js.history.length !== py.history.length) fail(`${tag}: history length js ${js.history.length} py ${py.history.length}`);
    if (js.checks.length !== py.checks.length) fail(`${tag}: ${js.checks.length} checks in js, ${py.checks.length} in py`);
    py.checks.forEach((pc, k) => {
      const jc = js.checks[k];
      if (!jc || jc.name !== pc.name || jc.passed !== pc.passed) fail(`${tag}: check ${pc.name} js ${jc && jc.name} ${jc && jc.passed} py ${pc.passed}`);
    });
    const dsn = arrDiff(js.sn, py.sn, 1);
    worstArr = Math.max(worstArr, dsn === Infinity ? 0 : dsn);
    if (!(dsn <= 1e-9)) fail(`${tag}: sigma_n arrays differ by ${dsn}`);
    if (!(arrDiff(js.tau_used, py.tau_used, 1) <= 1e-9)) fail(`${tag}: tau_used arrays differ`);
  }
}

let fx = 0;
for (const f of ref.fixtures) {
  const s = f.strength;
  const st = E.makeStrength(s.sigma_ci, s.mb, s.s, s.a, s.model, s.nu, s.c_mc, s.phi_mc);
  const r = E.solveStrength(Float64Array.from(f.targets), st);
  for (const k of ["sn", "phi", "c", "tau", "M", "s1", "s3"]) {
    const d = arrDiff(r[k], f[k], 1e-9);
    fx++;
    if (!(d <= 1e-11)) fail(`fixture ${f.name}: ${k} differs by ${d}\n  js ${Array.from(r[k])}\n  py ${f[k]}`);
  }
  if (Array.from(r.code).join() !== f.code.join()) fail(`fixture ${f.name}: codes js ${Array.from(r.code)} py ${f.code}`);
  if (Array.from(r.n_brackets).join() !== f.n_brackets.join()) fail(`fixture ${f.name}: brackets js ${Array.from(r.n_brackets)} py ${f.n_brackets}`);
  if (rel(E.projectionLimit(st), f.projection_limit) > 1e-12) fail(`fixture ${f.name}: projection limit`);
  if (f.ps) {
    const ps = f.ps.s3.map((x) => E.principalState(x, st));
    for (const k of ["s1", "M", "sn", "tau"]) {
      const d = arrDiff(ps.map((p) => p[k]), f.ps[k], 1e-9);
      if (!(d <= 1e-12)) fail(`fixture ${f.name}: principal ${k} differs by ${d}`);
    }
    if (ps.map((p) => p.code).join() !== f.ps.code.join()) fail(`fixture ${f.name}: principal codes`);
  }
}

for (const o of ref.original_closure) {
  const [mb, s] = E.hoekBrownParameters(42, 10, 0);
  const r = E.originalClosureMinimum(77.7, mb, s, 0.5, o.nu);
  for (const k of ["sigma3_min", "sn_min", "phi_deg"]) if (rel(r[k], o[k]) > 1e-9) fail(`original closure nu ${o.nu}: ${k} js ${r[k]} py ${o[k]}`);
}
for (const x of ref.extents) {
  const e = E.slipExtent(inputsOf(x.inputs));
  if (rel(e[0], x.extent[0]) > 1e-14 || rel(e[1], x.extent[1]) > 1e-14) fail(`extent ${x.key}: js ${e} py ${x.extent}`);
}
const hb = E.hoekBrownParameters(42, 10, 0);
hb.forEach((v, k) => { if (rel(v, ref.hb[k]) > 1e-15) fail(`hb ${k}`); });

console.log(`${compared} runs, ${fx} fixture arrays, ${ref.original_closure.length} original-closure minima, ` +
  `${ref.extents.length} extents compared; worst FoS rel diff ${worstFos.toExponential(2)}, worst sigma_n rel diff ` +
  `${worstArr.toExponential(2)}; ${fails} failures`);
process.exit(fails ? 1 : 0);
