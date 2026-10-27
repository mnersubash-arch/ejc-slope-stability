// Run a batch of analyses with the shipped web engine. Usage: node run_engine.js jobs.json out.json
// jobs: [{id, dim: 2|3, inputs: {...}, arrays: bool}] -> results: {id: {fos, status, reason, iterations, count, ...}}
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "webapp", "engine.js"));
const jobs = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
const out = {};
for (const j of jobs) {
  let r;
  try {
    const inp = E.inputs(j.inputs);
    r = j.dim === 3 ? E.run3d(inp) : E.run2d(inp);
  } catch (e) {
    out[j.id] = { error: String(e.message || e) };
    continue;
  }
  const o = { fos: r.fos, f_trial: r.f_trial, status: r.status, reason: r.reason, iterations: r.iterations, count: r.mesh.count,
    projections: r.projections, checks: r.checks.map((c) => [c.name, c.value, c.passed]) };
  if (j.arrays && r.sn) {
    const m = r.mesh;
    o.mesh = { x: Array.from(m.x), y: Array.from(m.y), zb: Array.from(m.zb), zt: Array.from(m.zt), h: Array.from(m.h), ah: Array.from(m.ah),
      area: Array.from(m.area), g: Array.from(m.g), rho: Array.from(m.rho), kap: Array.from(m.kap), w: Array.from(m.w), u: Array.from(m.u),
      dx: m.dx, dy: m.dy, x0: m.x0, x1: m.x1 };
    o.sn = Array.from(r.sn);
    if (r.tau_used) {
      o.tau_used = Array.from(r.tau_used); o.c_used = Array.from(r.c_used); o.tan_used = Array.from(r.tan_used);
      const s = r.strength;
      o.s1 = Array.from(s.s1); o.s2 = Array.from(s.s2); o.s3 = Array.from(s.s3); o.phi = Array.from(s.phi); o.c = Array.from(s.c);
      o.tau = Array.from(s.tau); o.M = Array.from(s.M); o.code = Array.from(s.code);
    }
  }
  out[j.id] = o;
}
fs.writeFileSync(process.argv[3], JSON.stringify(out));
