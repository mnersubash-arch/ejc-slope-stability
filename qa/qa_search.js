// Critical-surface searches with the shipped web engine (the app itself evaluates one prescribed surface).
// 2D: the circle with the lowest Bishop FoS; 3D: the sphere with the lowest Hungr FoS. Nelder-Mead over (x_c, z_c, R)
// from several starts, default grid (circles that cross the toe level at or ahead of the toe). Then the minima are
// re-solved on finer meshes. Compares the GHB 2D minimum with Carranza-Torres & Hormazabal (2020): circle
// (207.28, 586.53, 622.08), FoS 2.01.
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "webapp", "engine.js"));
const OUT = path.join(__dirname, "out");                 // results (git-ignored)
fs.mkdirSync(OUT, { recursive: true });
const t0 = Date.now();
const log = (s) => console.error(`[${((Date.now() - t0) / 1000).toFixed(0)} s] ${s}`);

function fos(dim, p, extra) {
  const [xc, zc, R] = p;
  if (!(R > zc + 1 && zc > 0 && R < 3000)) return 1e3;
  try {
    const inp = E.inputs(Object.assign({ xc, zc, R }, extra));
    const r = dim === 3 ? E.run3d(inp) : E.run2d(inp);
    return r.fos !== null ? r.fos : 1e3 - 1;
  } catch (e) { return 1e3; }
}
function nelderMead(f, x0, step, iters) {
  const n = x0.length;
  let pts = [x0.slice()];
  for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += step[i]; pts.push(p); }
  let vals = pts.map(f);
  let evals = n + 1;
  for (let it = 0; it < iters; it++) {
    const idx = vals.map((v, i) => i).sort((a, b) => vals[a] - vals[b]);
    pts = idx.map((i) => pts[i]); vals = idx.map((i) => vals[i]);
    const c = Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) c[j] += pts[i][j] / n;
    const at = (t) => c.map((cj, j) => cj + t * (pts[n][j] - cj));
    const xr = at(-1), fr = f(xr); evals++;
    if (fr < vals[0]) { const xe = at(-2), fe = f(xe); evals++; if (fe < fr) { pts[n] = xe; vals[n] = fe; } else { pts[n] = xr; vals[n] = fr; } }
    else if (fr < vals[n - 1]) { pts[n] = xr; vals[n] = fr; }
    else {
      const xk = fr < vals[n] ? at(-0.5) : at(0.5), fk = f(xk); evals++;
      if (fk < Math.min(fr, vals[n])) { pts[n] = xk; vals[n] = fk; }
      else { for (let i = 1; i <= n; i++) { pts[i] = pts[i].map((v, j) => pts[0][j] + 0.5 * (v - pts[0][j])); vals[i] = f(pts[i]); evals++; } }
    }
    const spread = Math.max(...pts.map((p) => Math.max(...p.map((v, j) => Math.abs(v - pts[0][j])))));
    if (spread < 0.05) break;
  }
  const i = vals.indexOf(Math.min(...vals));
  return { x: pts[i], f: vals[i], evals };
}

const out = {};
const best2 = {};
const starts = [[355.7423115, 644.4192949, 725.2669975], [207.28, 586.53, 622.08], [150, 450, 480], [300, 800, 860]];
for (const model of ["ghb", "ejc_moscarda"]) for (const wet of [false, true]) {
  const extra = { model, wet };
  const tag = `${model}_${wet ? "wet" : "dry"}`;
  for (const dim of [2, 3]) {
    let best = null;
    for (const s of (dim === 2 ? starts : [best2[tag].x, starts[0]])) {
      const r = nelderMead((p) => fos(dim, p, extra), s, [25, 25, 25], dim === 2 ? 400 : 150);
      if (!best || r.f < best.f) best = r;
    }
    if (dim === 2) best2[tag] = best;
    const fine = dim === 2 ? fos(2, best.x, Object.assign({ n2: 640 }, extra)) : fos(3, best.x, Object.assign({ n3: 80 }, extra));
    const onMain = fos(dim, starts[0], extra);
    out[`${tag}_${dim}d`] = { centre_radius: best.x, fos: best.f, fos_fine: fine, fos_main_surface: onMain, evals: best.evals };
    log(`${tag} ${dim}D: min ${best.f.toFixed(4)} at (${best.x.map((v) => v.toFixed(1)).join(", ")}); fine mesh ${fine.toFixed(4)}; main surface ${onMain.toFixed(4)}`);
  }
  // the 2D critical circle as a sphere, and the 3D critical sphere's central section, for the 3D/2D ratio
  const c2 = out[`${tag}_2d`].centre_radius, c3 = out[`${tag}_3d`].centre_radius;
  out[`${tag}_ratio`] = { sphere_on_2d_critical: fos(3, c2, extra), section_of_3d_critical: fos(2, c3, extra),
    ratio_minima: out[`${tag}_3d`].fos / out[`${tag}_2d`].fos };
}
fs.writeFileSync(path.join(OUT, "qa_search.json"), JSON.stringify(out, null, 1));
log("done");
