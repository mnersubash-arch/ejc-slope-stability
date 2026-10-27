// Compare the verification table of checks.js with the Python table (tools/ref/verification_rows.json from
// dump_reference.py), field by field. Usage (from webapp/): node tools/verify_checks.js [--slow]
// With --slow the slow rows (320x320, 240x240, 5120 slices; about a minute) are computed too and printed, for
// comparison with ejc_app/VERIFICATION.md; the other rows are still compared with the Python table.
"use strict";
const fs = require("fs");
const path = require("path");
const E = require(path.join(__dirname, "..", "engine.js"));
const PRESETS = require(path.join(__dirname, "..", "presets.js"));
const C = require(path.join(__dirname, "..", "checks.js")).make(E, PRESETS);

const slow = process.argv.includes("--slow");
let py = JSON.parse(fs.readFileSync(path.join(__dirname, "ref", "verification_rows.json"), "utf8"));
const t = Date.now();
let js = C.verificationRows({ slow, progress: (s) => process.stderr.write(`... ${s}\n`) });
console.log(`js table in ${((Date.now() - t) / 1000).toFixed(1)} s:`, C.summary(js));
if (slow) {
  const skipped = new Set(py.filter((r) => r.status === "SKIPPED").map((r) => r.group + "|" + r.item));
  for (const r of js) if (skipped.has(r.group + "|" + r.item)) console.log(`slow: ${r.item} | target ${r.target} | ${r.computed} | ${r.difference} | ${r.status} | ${r.note}`);
  js = js.filter((r) => !skipped.has(r.group + "|" + r.item));
  py = py.filter((r) => !skipped.has(r.group + "|" + r.item));
}
let diffs = 0;
if (js.length !== py.length) { console.log(`row count js ${js.length} py ${py.length}`); diffs++; }
for (let k = 0; k < Math.max(js.length, py.length); k++) {
  const a = js[k] || {}, b = py[k] || {};
  for (const f of ["group", "item", "target", "computed", "difference", "tolerance", "status", "note"]) {
    if (a[f] !== b[f]) { diffs++; console.log(`row ${k} ${f}:\n  js ${JSON.stringify(a[f])}\n  py ${JSON.stringify(b[f])}`); }
  }
}
console.log(diffs ? `${diffs} differences` : "identical to the Python table");
process.exit(diffs ? 1 : 0);
