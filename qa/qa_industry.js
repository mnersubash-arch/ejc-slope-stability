// The web engine against published/commercial results on the same surfaces: D'Esterre (2024) Slide2/Slide3 and his
// Hungr spreadsheet (24 Mohr-Coulomb cases), Law (2024) Slide2 and Excel (GHB), Carranza-Torres & Hormazabal (2020).
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "webapp", "engine.js"));
const P = require(path.join(__dirname, "..", "webapp", "presets.js"));
const OUT = path.join(__dirname, "out");                 // results (git-ignored)
fs.mkdirSync(OUT, { recursive: true });

const out = { desterre: [], law: [] };
const fos = (dim, inp, kw) => { try { const r = dim === 3 ? E.run3d(E.inputs(Object.assign({}, inp, kw))) : E.run2d(E.inputs(Object.assign({}, inp, kw))); return r.fos; } catch (e) { return null; } };
for (const p of P.filter((q) => q.group === "D'Esterre (2024)")) {
  const d = p.data, rep = d.reported;
  const row = { id: d.id, slope: d.slope_H_to_L, ru: p.inputs.ru, floor: d.needs_floor_columns, lowest: d.lowest_point_z_m,
    slide2: rep.Slide2_2D, slide3: rep.Slide3_3D, sheet3d: rep.DEsterre_spreadsheet_3D_Hungr,
    app2_40: fos(2, p.inputs, {}), app2_640: fos(2, p.inputs, { n2: 640 }), app3_20: fos(3, p.inputs, {}), app3_80: fos(3, p.inputs, { n3: 80 }) };
  out.desterre.push(row);
}
for (const p of P.filter((q) => q.group === "Law (2024)")) {
  const d = p.data;
  out.law.push({ id: d.id, water: d.water.split(":")[0], reported: d.reported, app2_40: fos(2, p.inputs, {}), app2_640: fos(2, p.inputs, { n2: 640 }),
    app3_20: fos(3, p.inputs, {}), app3_80: fos(3, p.inputs, { n3: 80 }) });
}
const ct = { xc: 207.28, zc: 586.53, R: 622.08, model: "ghb" };
out.ct = { published: 2.01, app2_40: fos(2, ct, {}), app2_640: fos(2, ct, { n2: 640 }), app2_5120: fos(2, ct, { n2: 5120 }) };
fs.writeFileSync(path.join(OUT, "qa_industry.json"), JSON.stringify(out, null, 1));
const pct = (a, b) => (a == null || b == null ? "  n/a " : `${((a / b - 1) * 100).toFixed(1).padStart(5)}%`);
console.log("D'Esterre: id | Slide2 | app 2D 40/640 | diff vs Slide2 | Slide3 | app 3D 20/80 | diff vs Slide3 | his sheet 3D | diff");
for (const r of out.desterre) console.log(`${r.id.padEnd(17)} ${r.slide2.toFixed(2)} | ${r.app2_40.toFixed(3)}/${r.app2_640.toFixed(3)} | ${pct(r.app2_640, r.slide2)} | ${r.slide3.toFixed(2)} | ${r.app3_20.toFixed(3)}/${r.app3_80.toFixed(3)} | ${pct(r.app3_80, r.slide3)} | ${r.sheet3d.toFixed(2)} | ${pct(r.app3_80, r.sheet3d)}`);
for (const r of out.law) console.log(`${r.id} (${r.water}): reported ${JSON.stringify(r.reported)} | app 2D ${r.app2_40.toFixed(4)}/${r.app2_640.toFixed(4)} | app 3D ${r.app3_20.toFixed(4)}/${r.app3_80.toFixed(4)}`);
console.log("Carranza-Torres:", JSON.stringify(out.ct));
