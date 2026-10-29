// Compare engine.js with the Python engine's reference results (reference.json from dump_reference.py).
// Usage (from webapp/): python3 tools/dump_reference.py && node tools/verify_engine.js [tools/ref/reference.json]
// Passes when every analysis matches: FoS and trial F to 1e-12 (relative), the same status, iteration count,
// projection count, mesh and check results, sigma_n and tau_s arrays to 1e-9; strength fixtures to 1e-11; the
// critical-surface searches field by field (every surface tried, numbers to 1e-12).
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "engine.js"));

const ref = JSON.parse(fs.readFileSync(process.argv[2] || path.join(__dirname, "ref", "reference.json"), "utf8"));
const val = (v) => (v === "nan" ? NaN : v === "inf" ? Infinity : v === "-inf" ? -Infinity : v);
let fails = 0, compared = 0, worstFos = 0, worstArr = 0, masses = 0, worstMass = 0;
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
    if (js.mesh.n !== py.n || js.mesh.ny !== py.ny || rel(js.mesh.dx, py.dx) > 1e-14 || rel(js.mesh.dy, py.dy) > 1e-14)
      fail(`${tag}: grid js ${js.mesh.n} x ${js.mesh.ny} (${js.mesh.dx} by ${js.mesh.dy} m) py ${py.n} x ${py.ny} (${py.dx} by ${py.dy} m)`);
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
    const pm = c[`m${dim}`];                     // whole-mass forces: totals to 1e-12 of the weight, moments to 1e-12
    if (pm) {
      let jm = null;
      try { jm = E.massForces(inp, js); } catch (e) { fail(`${tag}: massForces threw ${e.message}`); }
      if (jm) {
        masses++;
        if (rel(jm.weight, pm.weight) > 1e-13) fail(`${tag}: mass weight js ${jm.weight} py ${pm.weight}`);
        for (const k of ["W", "N", "U", "S", "net"]) for (const ax of ["x", "z"]) {
          const d = Math.abs(jm.totals[k][ax] - val(pm.totals[k][ax])) / pm.weight;
          worstMass = Math.max(worstMass, d);
          if (!(d <= 1e-12)) fail(`${tag}: mass total ${k}.${ax} js ${jm.totals[k][ax]} py ${pm.totals[k][ax]}`);
        }
        for (const k of ["W", "S", "net"]) if (!(Math.abs(jm.moments[k] - val(pm.moments[k])) <= 1e-12 * Math.abs(val(pm.moments.W)))) fail(`${tag}: mass moment ${k} js ${jm.moments[k]} py ${pm.moments[k]}`);
        if (jm.method.form !== pm.method.form || !(Math.abs(jm.method.residual - val(pm.method.residual)) <= 1e-12)) fail(`${tag}: method equation js ${jm.method.residual} py ${pm.method.residual}`);
        if (!(Math.abs(jm.vertical_max_rel - val(pm.vertical_max_rel)) <= 1e-12)) fail(`${tag}: vertical balance js ${jm.vertical_max_rel} py ${pm.vertical_max_rel}`);
      }
    }
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
// The same minima against 50-digit values on the edge sigma2 = sigma3 of the admissible set (QA/QC finding S1; see
// ejc_app/tests/test_strength.py): [sigma3_min, sigma_n,min] in MPa.
const EDGE_MINIMA = { 0.05: [0.380094742856100265, 0.915855437380532255], 0.1: [1.62220790622158018, 3.42932759983678206],
  0.2: [10.975624118170999, 19.2255957374544289], 0.271: [34.3759758649727541, 52.9663755316253943],
  0.35: [133.35940221633931, 177.816163701316242], 0.45: [1982.70231915863376, 2187.80980520517984],
  0.49: [58768.4063005811189, 59951.6628165205653] };
for (const [nu, [s3, sn]] of Object.entries(EDGE_MINIMA)) {
  const [mb, s] = E.hoekBrownParameters(42, 10, 0);
  const r = E.originalClosureMinimum(77.7, mb, s, 0.5, Number(nu));
  if (!(rel(r.sigma3_min, s3) < 1e-9 && rel(r.sn_min, sn) < 1e-9)) fail(`original closure nu ${nu}: js ${r.sigma3_min}, ${r.sn_min}; 50-digit ${s3}, ${sn}`);
}
for (const x of ref.extents) {
  const e = E.slipExtent(inputsOf(x.inputs));
  if (rel(e[0], x.extent[0]) > 1e-14 || rel(e[1], x.extent[1]) > 1e-14) fail(`extent ${x.key}: js ${e} py ${x.extent}`);
}
const hb = E.hoekBrownParameters(42, 10, 0);
hb.forEach((v, k) => { if (rel(v, ref.hb[k]) > 1e-15) fail(`hb ${k}`); });

for (const x of ref.faces || []) {                    // vertical faces where the circle ends below the ground
  const f = E.verticalFaces(inputsOf(x.inputs));
  for (const k of ["back", "front"]) if (!(Math.abs(f[k] - x.faces[k]) <= 1e-12 * x.inputs.H)) fail(`vertical ${k} face: js ${f[k]} py ${x.faces[k]} (${JSON.stringify(x.inputs)})`);
}

for (const x of ref.columns || []) {                  // 3D columns across y, plan size, counts for a size
  const inp = inputsOf(x.inputs), tag = `columns ${JSON.stringify({ n3: x.inputs.n3, n3y: x.inputs.n3y })}`;
  if (x.counts) {
    const c = E.countsForSize(inp, val(x.dx), val(x.dy));
    if (c.join() !== x.counts.join()) fail(`${tag}: counts for ${x.dx} by ${x.dy} m js ${c} py ${x.counts}`);
    continue;
  }
  const a = E.columnsAcross(inp, x.n), sz = E.columnSize(inp);
  if (a !== x.across) fail(`${tag}: columns across for n = ${x.n} js ${a} py ${x.across}`);
  if (rel(sz[0], x.size[0]) > 1e-14 || rel(sz[1], x.size[1]) > 1e-14) fail(`${tag}: column size js ${sz} py ${x.size}`);
  const w = E.withN3(inp, x.n);
  if (w.n3 !== x.n || (inp.n3y === null ? w.n3y !== null : w.n3y !== x.across)) fail(`${tag}: withN3 js ${w.n3} x ${w.n3y}`);
}

// critical-surface searches: every field, the clouds of surfaces tried included; FoS and other numbers to 1e-12
// (relative), lengths and search parameters to 1e-12 of H (V8's tan and the C library's can differ in the last bit, so
// a coordinate near zero can differ by more than 1e-12 of itself)
let searches = 0, surfaces = 0, worstLen = 0;
const sameNum = (a, b) => rel(a === undefined ? null : a, b === undefined ? null : b) <= 1e-12;
const sameArr = (a, b) => (a === undefined || a === null ? b === null || b === undefined
  : !!b && a.length === b.length && a.every((x, k) => sameNum(x, b[k])));
for (const s of ref.searches || []) {
  const res = E.criticalSearch(inputsOf(s.inputs), Object.assign({ dims: s.dims }, s.opts));
  const py = s.result, tag = `search ${s.name}`, H = s.inputs.H;
  const sameLen = (a, b) => {
    if (a === undefined || a === null || b === undefined || b === null) return (a === undefined || a === null) === (b === undefined || b === null);
    const d = Math.abs(a - b) / H;
    worstLen = Math.max(worstLen, d);
    return d <= 1e-12;
  };
  const sameLens = (a, b) => (a === undefined || a === null ? b === null || b === undefined
    : !!b && a.length === b.length && a.every((x, k) => sameLen(x, b[k])));
  searches++;
  for (const [k, pv] of Object.entries(py.prescribed)) {
    const jv = res.prescribed[k];
    if (typeof pv === "string" ? jv !== pv : !sameNum(jv, pv)) fail(`${tag}: prescribed ${k} js ${jv} py ${pv}`);
  }
  for (const k of ["back", "front", "starts", "tol", "max_eval", "restarts", "coverage", "base_circle"]) if (res.settings[k] !== py.settings[k]) fail(`${tag}: setting ${k}`);
  if (res.settings.grid.join() !== py.settings.grid.join()) fail(`${tag}: grid`);
  for (const d of s.dims) {
    const j = res["r" + d], p = py["r" + d], t = `${tag} ${d}D`;
    surfaces += p.evaluations;
    for (const k of ["fos", "coverage", "evaluations"]) if (!sameNum(j[k], p[k])) fail(`${t}: ${k} js ${j[k]} py ${p[k]}`);
    for (const k of ["xc", "zc", "R"]) if (!sameLen(j[k], p[k])) fail(`${t}: ${k} js ${j[k]} py ${p[k]}`);
    for (const k of ["p", "front", "back"]) if (!sameLens(j[k], p[k])) fail(`${t}: ${k} js ${j[k]} py ${p[k]}`);
    for (const k of ["counts", "skipped_why", "edges", "message"]) if (JSON.stringify(j[k]) !== JSON.stringify(p[k])) fail(`${t}: ${k} js ${JSON.stringify(j[k])} py ${JSON.stringify(p[k])}`);
    if (!!j.check !== !!p.check) fail(`${t}: mesh check present js ${!!j.check} py ${!!p.check}`);
    else if (p.check) for (const [k, pv] of Object.entries(p.check)) {
      const jv = j.check[k];
      if (typeof pv === "string" ? jv !== pv : !sameNum(jv, pv)) fail(`${t}: mesh check ${k} js ${jv} py ${pv}`);
    }
    if (j.runs.length !== p.runs.length) fail(`${t}: ${j.runs.length} runs in js, ${p.runs.length} in py`);
    else p.runs.forEach((q, k) => {
      const r = j.runs[k];
      if (!sameLens(r.p0, q.p0) || !sameLens(r.p, q.p) || !sameNum(r.f0, q.f0) || !sameNum(r.f, q.f)) fail(`${t}: run ${k} js ${JSON.stringify(r)} py ${JSON.stringify(q)}`);
    });
    for (const k of ["xc", "zc", "R"]) if (!sameLens(j.cloud[k], p.cloud[k])) fail(`${t}: cloud ${k} differs`);
    if (!sameArr(j.cloud.fos, p.cloud.fos)) fail(`${t}: cloud fos differs`);
  }
}

console.log(`${compared} runs, ${fx} fixture arrays, ${ref.original_closure.length} original-closure minima ` +
  `(and ${Object.keys(EDGE_MINIMA).length} against 50-digit values), ` +
  `${ref.extents.length} extents, ${(ref.faces || []).length} vertical-face cases, ${(ref.columns || []).length} column-shape cases, ${masses} whole-mass force sets compared (worst ${worstMass.toExponential(1)} of the weight), ` +
  `${searches} critical-surface searches (${surfaces} surfaces, every field; worst length diff ${worstLen.toExponential(1)} of H); ` +
  `worst FoS rel diff ${worstFos.toExponential(2)}, worst sigma_n rel diff ` +
  `${worstArr.toExponential(2)}; ${fails} failures`);
process.exit(fails ? 1 : 0);
