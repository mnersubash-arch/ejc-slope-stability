// D'Esterre (2024): is each published circle near the app's own critical circle? Nelder-Mead over (x_c, z_c, R) from
// the published circle, base-circle grid, Bishop 2D with 640 slices. If the published circles are Slide2's critical
// circles and the app agrees with Slide2, the app's local minimum should sit close to the published circle and FoS.
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "webapp", "engine.js"));
const P = require(path.join(__dirname, "..", "webapp", "presets.js"));
const OUT = path.join(__dirname, "out");                 // results (git-ignored)
fs.mkdirSync(OUT, { recursive: true });

function fos(inp, p) {
  const [xc, zc, R] = p;
  if (!(R > 0.1)) return 1e3;
  try { const r = E.run2d(E.inputs(Object.assign({}, inp, { xc, zc, R, n2: 640 }))); return r.fos !== null ? r.fos : 1e3 - 1; } catch (e) { return 1e3; }
}
function nm(f, x0, step, iters) {
  const n = x0.length; let pts = [x0.slice()];
  for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += step[i]; pts.push(p); }
  let vals = pts.map(f);
  for (let it = 0; it < iters; it++) {
    const idx = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]); pts = idx.map((i) => pts[i]); vals = idx.map((i) => vals[i]);
    const c = Array(n).fill(0); for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += pts[i][j] / n;
    const at = (t) => c.map((cj, j) => cj + t * (pts[n][j] - cj));
    const xr = at(-1), fr = f(xr);
    if (fr < vals[0]) { const xe = at(-2), fe = f(xe); if (fe < fr) { pts[n] = xe; vals[n] = fe; } else { pts[n] = xr; vals[n] = fr; } }
    else if (fr < vals[n - 1]) { pts[n] = xr; vals[n] = fr; }
    else { const xk = fr < vals[n] ? at(-0.5) : at(0.5), fk = f(xk); if (fk < Math.min(fr, vals[n])) { pts[n] = xk; vals[n] = fk; } else { for (let i = 1; i <= n; i++) { pts[i] = pts[i].map((v, j) => pts[0][j] + 0.5 * (v - pts[0][j])); vals[i] = f(pts[i]); } } }
    if (Math.max(...pts.map((p) => Math.max(...p.map((v, j) => Math.abs(v - pts[0][j]))))) < 1e-4) break;
  }
  const i = vals.indexOf(Math.min(...vals)); return { x: pts[i], f: vals[i] };
}
const out = [];
for (const p of P.filter((q) => q.group === "D'Esterre (2024)")) {
  const d = p.data, inp = p.inputs, x0 = [inp.xc, inp.zc, inp.R];
  const onPub = fos(inp, x0);
  let best = null;
  for (const st of [0.3, 1.0]) { const r = nm((q) => fos(inp, q), x0, [st, st, st], 600); if (!best || r.f < best.f) best = r; }
  const shift = Math.hypot(best.x[0] - x0[0], best.x[1] - x0[1]);
  out.push({ id: d.id, slide2: d.reported.Slide2_2D, app_published: onPub, app_min: best.f, min_circle: best.x, centre_shift_m: shift, dR: best.x[2] - x0[2] });
  console.log(`${d.id.padEnd(16)} Slide2 ${d.reported.Slide2_2D.toFixed(2)} | app on published ${onPub.toFixed(3)} (${((onPub / d.reported.Slide2_2D - 1) * 100).toFixed(1)}%) | app minimum ${best.f.toFixed(3)} (${((best.f / d.reported.Slide2_2D - 1) * 100).toFixed(1)}%) at (${best.x.map((v) => v.toFixed(2)).join(", ")}), centre moved ${shift.toFixed(2)} m, dR ${(best.x[2] - x0[2]).toFixed(2)} m`);
}
fs.writeFileSync(path.join(OUT, "qa_desterre_search.json"), JSON.stringify(out, null, 1));
