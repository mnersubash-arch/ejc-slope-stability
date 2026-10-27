/* Verification targets of ejc_app/verification.py (brief section 3), recomputed with engine.js.
 *
 * Same rows, targets, tolerances and wording as the Python table, so the two can be compared line by line
 * (verify_checks.js). Tolerances: FoS given to 9 decimals must match to 1e-6; values given to 3-4 decimals must
 * match to rounding (|computed - target| <= half a unit in the last decimal shown). INFO rows are reference values,
 * not targets. Without `slow`, the fine-mesh rows are listed as SKIPPED.
 */
(function (root) {
  "use strict";

  function make(E, PRESETS) {
    const MAIN = E.inputs();
    const CARRANZA = E.inputs({ xc: 207.28, zc: 586.53, R: 622.08 });
    const fmtG = E.fmtG;
    const fmtE = (x, p, plus) => {                          // Python format(x, "+.pe") / ".pe"
      const [m, e] = x.toExponential(p).split("e");
      const n = parseInt(e, 10);
      const sign = plus && x >= 0 ? "+" : "";
      return `${sign}${m}e${n < 0 ? "-" : "+"}${String(Math.abs(n)).padStart(2, "0")}`;
    };
    const fixed = (v, d) => v.toFixed(d);
    const fmt = (v, d) => {                                 // verification._fmt
      if (v === null || v === undefined || (typeof v === "number" && !Number.isFinite(v))) return "none";
      if (typeof v === "string") return v;
      return d !== undefined ? fixed(v, d) : fmtG(v, 10);
    };
    const pyStr = (v) => (typeof v === "boolean" ? (v ? "True" : "False") : String(v));
    const row = (group, item, target, computed, difference, tolerance, status, note = "") =>
      ({ group, item, target, computed, difference, tolerance, status, note });

    function byRounding(group, item, computed, target, decimals, note = "") {
      const tol = 0.5 * Math.pow(10, -decimals);
      if (computed === null || computed === undefined || !Number.isFinite(computed))
        return row(group, item, fixed(target, decimals), "none", "n/a", `rounding to ${decimals} dp`, "FAIL", note);
      const diff = computed - target;
      const ok = Math.abs(diff) <= tol * (1 + 1e-12);
      return row(group, item, fixed(target, decimals), fmt(computed, decimals + 3), fmtE(diff, 2, true),
        `rounding to ${decimals} dp`, ok ? "PASS" : "FAIL", note);
    }
    function byAbs(group, item, computed, target, tol, decimals = 9, note = "") {
      if (computed === null || computed === undefined || !Number.isFinite(computed))
        return row(group, item, fixed(target, decimals), "none", "n/a", `abs <= ${fmtG(tol, 6)}`, "FAIL", note);
      const diff = computed - target;
      return row(group, item, fixed(target, decimals), fmt(computed, 12), fmtE(diff, 2, true), `abs <= ${fmtG(tol, 6)}`,
        Math.abs(diff) <= tol ? "PASS" : "FAIL", note);
    }
    function byRel(group, item, computed, reference, tol, note = "") {
      if (computed === null || reference === null)
        return row(group, item, "relative difference", "none", "n/a", `rel <= ${fmtG(tol, 6)}`, "FAIL", note);
      const rel = Math.abs(computed - reference) / Math.abs(reference);
      return row(group, item, fmt(reference, 12), fmt(computed, 12), `${fmtE(rel, 2, false)} (rel)`,
        `rel <= ${fmtG(tol, 6)}`, rel <= tol ? "PASS" : "FAIL", note);
    }
    function byExact(group, item, computed, target, note = "") {           // integers or strings
      const ok = computed === target;
      const diff = ok ? "0" : (typeof computed === "number" && typeof target === "number"
        ? `${computed - target >= 0 ? "+" : ""}${computed - target}` : "differs");
      return row(group, item, String(target), String(computed), diff, "exact", ok ? "PASS" : "FAIL", note);
    }
    function noFos(group, item, r, note = "") {
      const ok = r.fos === null;
      return row(group, item, "no accepted FoS", ok ? r.status : fmt(r.fos, 6), "n/a", "status", ok ? "PASS" : "FAIL",
        note || r.reason.slice(0, 160));
    }
    const info = (group, item, computed, target = "", note = "") =>
      row(group, item, target, computed, "", "reference only", "INFO", note);
    const SKIP_TEXT = "skipped (slow), see VERIFICATION.md";
    const skipped = (group, item, target, tolerance) => row(group, item, target, SKIP_TEXT, "", tolerance, "SKIPPED");

    /** Runs and caches the analyses the targets need. */
    function Runner() {
      const cache = new Map();
      const with_ = (base, kw) => E.inputs(Object.assign({}, base, kw || {}));
      return {
        r3(base = MAIN, kw = {}) {
          const inp = with_(base, kw), key = "3" + JSON.stringify(inp);
          if (!cache.has(key)) cache.set(key, E.run3d(inp));
          return cache.get(key);
        },
        r2(base = MAIN, kw = {}) {
          const inp = with_(base, kw), key = "2" + JSON.stringify(Object.assign({}, inp, { form: "FH", n3: 0 }));
          if (!cache.has(key)) cache.set(key, E.run2d(inp));
          return cache.get(key);
        },
      };
    }

    const presetInputs = (p) => E.inputs(p.inputs);
    const lawPresets = () => PRESETS.filter((p) => p.group === "Law (2024)");
    const desterrePresets = () => PRESETS.filter((p) => p.group === "D'Esterre (2024)");
    const valuesText = (d) => Object.entries(d).map(([k, v]) => `${k.replace(/_/g, " ")} ${pyStr(v)}`).join("; ");

    const CT_SLOW = [["2D EJC, 5120 slices (slow)", 2.121420, 6, (R) => R.r2(CARRANZA, { n2: 5120 })],
      ["3D GHB, 320x320 (slow)", 2.140776, 6, (R) => R.r3(CARRANZA, { model: "ghb", n3: 320 })],
      ["3D EJC, 320x320 (slow)", 2.267835, 6, (R) => R.r3(CARRANZA, { n3: 320 })]];
    const FINEST = [["Dry: F 3D 240x240", 2.7112, 4, (R) => R.r3(MAIN, { n3: 240 })],
      ["Dry: F 2D 5120", 2.4422, 4, (R) => R.r2(MAIN, { n2: 5120 })],
      ["Wet: F 3D 240x240", 1.2776, 4, (R) => R.r3(MAIN, { n3: 240, wet: true })],
      ["Wet: F 2D 5120", 1.0916, 4, (R) => R.r2(MAIN, { n2: 5120, wet: true })]];
    const COARSE_INFO = [["dry", false, "about 2.6"], ["wet", true, "about 4.4"]];

    /** The rows of verification_rows(); `progress(label)` is called before each group. */
    function verificationRows(opts = {}) {
      const slow = !!opts.slow, say = opts.progress || (() => {});
      const R = Runner();
      const rows = [];
      const push = (...r) => rows.push(...r);
      let g;

      g = "1.1 Rock mass (GSI 42, m_i 10, D 0)";
      say("rock mass");
      let [mb, s, aGsi] = E.hoekBrownParameters(42.0, 10.0, 0.0);
      push(byRounding(g, "m_b", mb, 1.260056, 6), byRounding(g, "s", s, 0.0015893, 7),
        byRounding(g, "a from GSI", aGsi, 0.509923, 6));

      say("strength fixtures");
      g = "3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2)";
      const st = E.strengthOf(MAIN);
      const fx = E.solveStrength(Float64Array.from([50.0, 400.0, 1600.0]), st);
      [[50, 67.45, 362.5, 482.9], [400, 59.90, 487.2, 1177.2], [1600, 50.21, 935.3, 2856.0]].forEach(([sn, phi, c, tau], k) => {
        push(byRounding(g, `sigma_n = ${sn} kPa: phi_i (deg)`, fx.phi[k] * E.RAD, phi, 2),
          byRounding(g, `sigma_n = ${sn} kPa: c_i (kPa)`, fx.c[k], c, 1),
          byRounding(g, `sigma_n = ${sn} kPa: tau_s (kPa)`, fx.tau[k], tau, 1));
      });
      const gh = E.solveStrength(Float64Array.from([400.0]), E.strengthOf(E.inputs({ model: "ghb" })));
      push(byRounding(g, "GHB, sigma_n = 400 kPa: phi_i (deg)", gh.phi[0] * E.RAD, 58.45, 2),
        byRounding(g, "GHB, sigma_n = 400 kPa: c_i (kPa)", gh.c[0], 460.0, 0),
        byRounding(g, "GHB, sigma_n = 400 kPa: tau_s (kPa)", gh.tau[0], 1111.0, 0));
      const xo = E.solveStrength(Float64Array.from([5740.3]), E.makeStrength(120.0, 3.5, 0.002, 0.505, "ejc_original", 0.05));
      const xm = E.solveStrength(Float64Array.from([5740.3]), E.makeStrength(120.0, 3.5, 0.002, 0.505, "ejc_moscarda", 0.05));
      let note = "sigma_ci 120 MPa, m_b 3.5, s 0.002, a 0.505, sigma_n 5.7403 MPa";
      push(byRounding(g, "Extra fixture, original closure (nu = 0.05): phi_i (deg)", xo.phi[0] * E.RAD, 55.0720, 4, note),
        byRounding(g, "Extra fixture, original closure (nu = 0.05): c_i (kPa)", xo.c[0], 2079.5, 1, note),
        info(g, "Extra fixture with the Moscarda closure", E.REASONS[xm.code[0]], "",
          "the Moscarda closure has no admissible state here"));

      say("sigma2 = sigma3 equivalence");
      g = "3 EJC route at sigma2 = sigma3 reproduces GHB";
      for (const wet of [false, true]) {
        for (const [dim, fn] of [["3D", "r3"], ["2D", "r2"]]) {
          const e = R[fn](MAIN, { model: "ejc_sigma2_eq_sigma3", wet }), h = R[fn](MAIN, { model: "ghb", wet });
          push(byRel(g, `${dim} ${wet ? "wet" : "dry"}, a = 0.5: F(EJC, sigma2 = sigma3) vs F(GHB)`, e.fos, h.fos, 1e-11));
        }
      }

      say("main cases");
      g = "3 Main sphere, Moscarda closure, EJC (20x20, 40 slices)";
      const d3 = R.r3(), d2 = R.r2();
      push(byExact(g, "Active half-columns (3D 20x20)", d3.mesh.count, 76), byExact(g, "Active slices (2D 40)", d2.mesh.count, 30));
      const w3 = R.r3(MAIN, { wet: true }), w2 = R.r2(MAIN, { wet: true });
      push(byAbs(g, "Dry, a = 0.5: F 3D", d3.fos, 2.706024342, 1e-6), byAbs(g, "Dry, a = 0.5: F 2D", d2.fos, 2.447910895, 1e-6),
        byAbs(g, "Wet, a = 0.5: F 3D", w3.fos, 1.272331405, 1e-6), byAbs(g, "Wet, a = 0.5: F 2D", w2.fos, 1.098069496, 1e-6),
        info(g, "Dry / wet, a = 0.5: 3D iterations and projections (default path)",
          `${d3.iterations} / ${w3.iterations} iterations; ${d3.projections} / ${w3.projections} projections`, "shown, not tested"));
      const gn = "3 Native v15 iteration path (start F = 1, Hungr explicit form, projection)";
      const n3d = R.r3(MAIN, { iteration: "native_v15" }), n3w = R.r3(MAIN, { iteration: "native_v15", wet: true });
      push(byAbs(gn, "Dry, a = 0.5: F 3D", n3d.fos, 2.706024342, 1e-6), byExact(gn, "Dry, a = 0.5: 3D iterations", n3d.iterations, 34),
        byAbs(gn, "Wet, a = 0.5: F 3D", n3w.fos, 1.272331405, 1e-6), byExact(gn, "Wet, a = 0.5: 3D iterations", n3w.iterations, 63),
        info(gn, "Projections dry / wet", `${n3d.projections} / ${n3w.projections}`, ""));
      for (const [wet, a, t3, t2, g3, g2] of [[false, 0.5, null, null, 2.554, 2.310], [true, 0.5, null, null, 1.200, 1.036],
        [false, 0.509923, 2.624, 2.376, 2.476, 2.242], [true, 0.509923, 1.228, 1.061, 1.158, 1.000]]) {
        const lab = `${wet ? "Wet" : "Dry"}, a = ${fmtG(a, 6)}`;
        if (t3 !== null) {
          push(byRounding(g, `${lab}: F 3D`, R.r3(MAIN, { wet, a_override: a }).fos, t3, 3),
            byRounding(g, `${lab}: F 2D`, R.r2(MAIN, { wet, a_override: a }).fos, t2, 3));
        }
        push(byRounding(g, `${lab}: GHB F 3D`, R.r3(MAIN, { wet, a_override: a, model: "ghb" }).fos, g3, 3),
          byRounding(g, `${lab}: GHB F 2D`, R.r2(MAIN, { wet, a_override: a, model: "ghb" }).fos, g2, 3));
      }

      say("moment form");
      g = "3 Moment form F_M (20x20, own iteration)";
      const fmD = R.r3(MAIN, { form: "FM" }), fmW = R.r3(MAIN, { form: "FM", wet: true });
      push(byRounding(g, "Dry: F_M 3D", fmD.fos, 2.7022, 4), byRounding(g, "Wet: F_M 3D", fmW.fos, 1.2695, 4));

      say("tension sensitivity");
      g = "3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law)";
      push(byRounding(g, "Dry: C0 3D", R.r3(MAIN, { tension: "C0" }).fos, 2.637, 3),
        byRounding(g, "Dry: T0 3D", R.r3(MAIN, { tension: "T0" }).fos, 2.597, 3),
        byRounding(g, "Dry: C0 2D", R.r2(MAIN, { tension: "C0" }).fos, 2.430, 3),
        byRounding(g, "Dry: T0 2D", R.r2(MAIN, { tension: "T0" }).fos, 2.429, 3),
        byRounding(g, "Wet: C0 3D", R.r3(MAIN, { tension: "C0", wet: true }).fos, 0.964, 3),
        noFos(g, "Wet: T0 3D", R.r3(MAIN, { tension: "T0", wet: true })),
        byRounding(g, "Wet: C0 2D", R.r2(MAIN, { tension: "C0", wet: true }).fos, 0.978, 3),
        noFos(g, "Wet: T0 2D", R.r2(MAIN, { tension: "T0", wet: true })));

      say("original closure");
      g = "1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5)";
      [mb, s] = E.hb(MAIN);
      const m0 = E.originalClosureMinimum(77.7, mb, s, 0.5, 0.271);
      push(byRounding(g, "nu = 0.271: sigma3,min (MPa)", m0.sigma3_min, 34.376, 3),
        byRounding(g, "nu = 0.271: sigma_n,min (MPa)", m0.sn_min, 52.97, 2));
      const md = E.buildMesh(MAIN, 3, 20), mw = E.buildMesh(E.inputs({ wet: true }), 3, 20);
      const capD = E.gravityCap(md, MAIN.gamma), capW = E.gravityCap(mw, MAIN.gamma);
      const amax = (a) => a.reduce((m, v) => (v > m ? v : m), -Infinity);
      push(byRounding(g, "Gravity cap max(gamma h - u), dry (MPa, 3D 20x20)", amax(capD) / 1000, 1.61, 2),
        byRounding(g, "Gravity cap max(gamma h - u), wet (MPa, 3D 20x20)", amax(capW) / 1000, 1.02, 2));
      for (const [nu, tgt, d] of [[0.05, 0.92, 2], [0.10, 3.43, 2], [0.20, 19.23, 2], [0.271, 52.97, 2], [0.35, 177.8, 1], [0.45, 2188.0, 0]]) {
        const mn = nu === 0.271 ? m0 : E.originalClosureMinimum(77.7, mb, s, 0.5, nu);
        push(byRounding(g, `Table 5.2: min sigma_n at nu = ${fmtG(nu, 6)} (MPa)`, mn.sn_min, tgt, d));
      }
      const thr = E.originalClosureMinimum(77.7, mb, s, 0.5, 0.05).sn_min * 1000.0;
      push(byExact(g, "nu = 0.05: dry columns with gamma h - u below the threshold (of 76)", capD.filter((v) => v < thr).length, 43));
      const orig = R.r3(MAIN, { model: "ejc_original" });
      push(row(g, "Original closure never produces a FoS", "NOT APPLICABLE", orig.status, "n/a", "status",
        orig.fos === null && orig.status === "NOT APPLICABLE" ? "PASS" : "FAIL"));

      say("Carranza-Torres & Hormazabal (2020)");
      g = "3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5";
      const ct2 = R.r2(CARRANZA, { model: "ghb" }), ct3 = R.r3(CARRANZA, { model: "ghb" });
      const ce2 = R.r2(CARRANZA), ce3 = R.r3(CARRANZA);
      push(byRounding(g, "2D GHB, 40 slices (brief)", ct2.fos, 2.0039, 4),
        byRounding(g, "2D GHB, 40 slices (thesis record)", ct2.fos, 2.003873, 6),
        byRounding(g, "3D GHB, 20x20 (thesis value)", ct3.fos, 2.142, 3),
        byRounding(g, "3D GHB, 20x20 (thesis record)", ct3.fos, 2.141937, 6),
        byRounding(g, "2D EJC, 40 slices", ce2.fos, 2.122598, 6), byRounding(g, "3D EJC, 20x20", ce3.fos, 2.269068, 6),
        info(g, "Published FoS (their method)", "", "2.01", "reference only"));
      if (slow) say("Carranza-Torres fine meshes (slow)");
      for (const [item, tgt, d, fn] of CT_SLOW)
        push(slow ? byRounding(g, item, fn(R).fos, tgt, d) : skipped(g, item, fixed(tgt, d), `rounding to ${d} dp`));
      g = "3 Moscarda (2025), main sphere, nu = 0.271 (reference values, not targets)";
      for (const [cs, wet, a, rep, sheet] of [[1, false, 0.5, 2.24, 2.245], [2, true, 0.5, 1.57, 1.567],
        [3, false, 0.509923, 2.04, 2.041], [4, true, 0.509923, 1.13, 1.132]]) {
        const ours = R.r2(MAIN, { wet, a_override: a }).fos;
        push(info(g, `Case ${cs}: ${wet ? "wet" : "dry"}, a = ${fmtG(a, 6)}, 2D FoS`, `this engine (Route 2, 40 slices): ${fmt(ours, 4)}`,
          `reported ${rep.toFixed(2)} (sheet at nu = 0.271: ${sheet.toFixed(3)})`, "his own slice model and the inherited route"));
      }

      say("Law (2024)");
      g = "3 Cross-check: Law (2024) presets, GHB 2D against the reference code lem_ref.py";
      for (const p of lawPresets()) {
        for (const [key, value] of Object.entries(p.data.engine_check)) {
          const n2 = parseInt(key.slice(key.lastIndexOf("_") + 1), 10);
          push(byRounding(g, `${p.data.id}: GHB 2D, ${n2} slices`, R.r2(presetInputs(p), { n2 }).fos, value, 4,
            value === 1.953 ? "read to 4 decimals (1.953 = 1.9530)" : ""));
        }
      }
      g = "3 Law (2024), reported values (reference only, not targets)";
      for (const p of lawPresets()) {
        push(info(g, `${p.data.id}: 2D / 3D FoS`, `this engine (GHB, 40 slices / 20x20): ${fmt(R.r2(presetInputs(p)).fos, 4)} / ` +
          `${fmt(R.r3(presetInputs(p)).fos, 4)}`, valuesText(p.data.reported), "her own models"));
      }

      say("Mohr-Coulomb");
      g = "3 Cross-check: Mohr-Coulomb option (c' = 100 kPa, phi' = 0, main sphere, dry) against the slope3d record";
      const mc = { model: "mc", c_mc: 100.0, phi_mc: 0.0 };
      push(byAbs(g, "F 3D, 20x20", R.r3(MAIN, mc).fos, 0.255354533635107, 1e-6, 15),
        byAbs(g, "F 2D, 40 slices", R.r2(MAIN, mc).fos, 0.190263799071471, 1e-6, 15));

      say("D'Esterre (2024)");
      g = "3 D'Esterre (2024): this engine against his reported values (reference only, not targets)";
      for (const p of desterrePresets()) {
        const rep = p.data.reported;
        push(info(g, p.data.id, `2D ${fmt(R.r2(presetInputs(p)).fos, 3)} (40 slices); 3D ${fmt(R.r3(presetInputs(p)).fos, 3)} (20x20)`,
          `Slide2 2D ${pyStr(rep.Slide2_2D)}; spreadsheet 3D ${pyStr(rep.DEsterre_spreadsheet_3D_Hungr)}; Slide3 3D ${pyStr(rep.Slide3_3D)}`,
          `needs_floor_columns ${String(p.data.needs_floor_columns).toLowerCase()}`));
      }

      g = "6 Limitations statements (about)";
      push(info(g, "F_H vs F_M difference, dry / wet (%)",
        [d3.fos, fmD.fos, w3.fos, fmW.fos].includes(null) ? "none"
          : `${((d3.fos - fmD.fos) / d3.fos * 100).toFixed(3)} / ${((w3.fos - fmW.fos) / w3.fos * 100).toFixed(3)}`, "0.14-0.27"));
      const c0w = R.r3(MAIN, { tension: "C0", wet: true }).fos;
      push(info(g, "Wet 3D FoS drop with C0 (%)", [w3.fos, c0w].includes(null) ? "none" : ((w3.fos - c0w) / w3.fos * 100).toFixed(1), "about 25"));

      g = "3 Finest mesh (240x240, 5120 slices; slow)";
      const g6 = "6 Limitations statements (about)";
      if (slow) say("finest mesh (slow)");
      for (const [item, tgt, d, fn] of FINEST) {
        if (slow) {
          const r = fn(R);
          note = r.fos ? `${r.projections} projections on the path` : r.reason.slice(0, 160);
          push(byRounding(g, item, r.fos, tgt, d, note));
        } else push(skipped(g, item, fixed(tgt, d), `rounding to ${d} dp`));
      }
      for (const [lab, wet, tgt] of COARSE_INFO) {
        const item = `10x10 below 240x240, ${lab} (%)`;
        if (slow) {
          const fin = R.r3(MAIN, { n3: 240, wet }), c10 = R.r3(MAIN, { n3: 10, wet });
          push(info(g6, item, [fin.fos, c10.fos].includes(null) ? "none" : ((fin.fos - c10.fos) / fin.fos * 100).toFixed(2), tgt));
        } else push(skipped(g6, item, tgt, "reference only"));
      }
      return rows;
    }

    const summary = (rows) => Object.fromEntries(["PASS", "FAIL", "SKIPPED", "INFO"].map((k) => [k, rows.filter((r) => r.status === k).length]));
    return { verificationRows, summary, SKIP_TEXT };
  }

  const api = { make };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.EJC_CHECKS = api;
})(typeof window !== "undefined" ? window : globalThis);
