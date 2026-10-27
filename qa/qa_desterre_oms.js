// Ordinary method of slices on the app's own slices for all 24 D'Esterre circles, two pore-pressure forms:
// N' = W cos(a) - u l  and  N' = (W - u b) cos(a). Same geometry, weights and pore pressures as the app's Bishop run.
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "webapp", "engine.js"));
const P = require(path.join(__dirname, "..", "webapp", "presets.js"));
const OUT = path.join(__dirname, "out");                 // results (git-ignored)
fs.mkdirSync(OUT, { recursive: true });
const out = [];
for (const p of P.filter((q) => q.group === "D'Esterre (2024)")) {
  const d = p.data, inp = Object.assign({}, p.inputs, { n2: 640 });
  const r = E.run2d(E.inputs(inp)), m = r.mesh, c = inp.c_mc, t = Math.tan(inp.phi_mc * Math.PI / 180);
  let na = 0, nb = 0, den = 0;
  for (let k = 0; k < m.count; k++) {
    na += c * m.area[k] + (m.w[k] * m.g[k] - m.u[k] * m.area[k]) * t;
    nb += c * m.area[k] + (m.w[k] - m.u[k] * m.ah[k]) * m.g[k] * t;
    den += m.w[k] * m.kap[k];
  }
  out.push({ id: d.id, slide2: d.reported.Slide2_2D, bishop: r.fos, oms_wcos_ul: na / den, oms_wub_cos: nb / den });
}
fs.writeFileSync(path.join(OUT, "qa_desterre_oms.json"), JSON.stringify(out, null, 1));
console.log(out.map((o) => `${o.id} ${((o.bishop / o.slide2 - 1) * 100).toFixed(1)} ${((o.oms_wub_cos / o.slide2 - 1) * 100).toFixed(1)} ${((o.oms_wcos_ul / o.slide2 - 1) * 100).toFixed(1)}`).join("\n"));
