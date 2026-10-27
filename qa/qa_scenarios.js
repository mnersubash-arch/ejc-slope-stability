// One-at-a-time sweeps over realistic open-pit inputs on the thesis slope and main sphere (EJC, Moscarda closure,
// 20x20 / 40 slices): does the FoS move in the physically expected direction, and when is it withheld?
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "webapp", "engine.js"));
const OUT = path.join(__dirname, "out");                 // results (git-ignored)
fs.mkdirSync(OUT, { recursive: true });

const BASE = { a_override: null };                 // a from GSI for realism (the thesis baseline fixes a = 0.5)
const sweeps = {
  gsi: { values: [20, 30, 42, 55, 70, 85], expect: "up" },
  sigma_ci: { values: [10, 25, 50, 77.7, 120, 200], expect: "up" },
  mi: { values: [4, 7, 10, 17, 25, 32], expect: "up" },
  D: { values: [0, 0.3, 0.5, 0.7, 0.85, 1.0], expect: "down" },
  gamma: { values: [22, 24, 26, 27, 28, 30], expect: "down" },
  Hw: { values: [0, 72, 144, 216, 288, 360], expect: "down", extra: { wet: true } },
  ru: { values: [0, 0.1, 0.2, 0.3, 0.4, 0.5], expect: "down" },
  scale: { values: [100, 200, 300, 360, 450, 600], expect: "down" },  // slope height, sphere scaled with it
};
const out = {};
for (const [name, sw] of Object.entries(sweeps)) {
  const rows = [];
  for (const v of sw.values) {
    let kw = Object.assign({}, BASE, sw.extra || {});
    if (name === "scale") {
      const k = v / 360;
      kw = Object.assign(kw, { H: v, xc: 355.7423115 * k, zc: 644.4192949 * k, R: 725.2669975 * k });
    } else kw[name] = v;
    const row = { value: v };
    for (const dim of [3, 2]) {
      const t = Date.now();
      let r;
      try { r = dim === 3 ? E.run3d(E.inputs(kw)) : E.run2d(E.inputs(kw)); } catch (e) { row[`err${dim}`] = e.message; continue; }
      row[`fos${dim}`] = r.fos; row[`status${dim}`] = r.status; row[`reason${dim}`] = r.fos === null ? r.reason.slice(0, 110) : "";
      row[`ms${dim}`] = Date.now() - t; row[`proj${dim}`] = r.projections; row[`iter${dim}`] = r.iterations;
      if (r.sn) row[`tensile${dim}`] = Array.from(r.sn).filter((x) => x < 0).length / r.mesh.count;
    }
    row.ratio = row.fos3 && row.fos2 ? row.fos3 / row.fos2 : null;
    rows.push(row);
  }
  // monotonic in the expected direction over the accepted values?
  const f3 = rows.map((r) => r.fos3).filter((x) => x != null), f2 = rows.map((r) => r.fos2).filter((x) => x != null);
  const mono = (a) => a.every((x, i) => i === 0 || (sw.expect === "up" ? x >= a[i - 1] : x <= a[i - 1]));
  out[name] = { rows, expect: sw.expect, monotonic3: mono(f3), monotonic2: mono(f2), accepted3: f3.length, accepted2: f2.length };
  console.log(`${name.padEnd(9)} (${sw.expect}) 3D: ${rows.map((r) => (r.fos3 != null ? r.fos3.toFixed(3) : "withheld")).join(" ")} | 2D: ${rows.map((r) => (r.fos2 != null ? r.fos2.toFixed(3) : "withheld")).join(" ")} | monotonic ${out[name].monotonic3}/${out[name].monotonic2} | ratio ${rows.map((r) => (r.ratio ? r.ratio.toFixed(3) : "-")).join(" ")}`);
  for (const r of rows) if (r.fos3 == null || r.fos2 == null) console.log(`   ${name} = ${r.value}: ${r.reason3 || r.err3 || ""} | ${r.reason2 || r.err2 || ""}`);
}
fs.writeFileSync(path.join(OUT, "qa_scenarios.json"), JSON.stringify(out, null, 1));
