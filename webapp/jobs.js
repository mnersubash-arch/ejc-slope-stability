/* Jobs for the page: one analysis (3D + 2D, or the original-closure check), the mesh study, the nu sweep, the
 * critical-surface search and the verification table. Runs in the Web Worker (worker.js) or, if the worker cannot
 * start, on the page itself. */
(function (root) {
  "use strict";

  function make(E, PRESETS, CHECKS) {
    const attempt = (fn) => {
      try { return { r: fn() }; } catch (e) { return { error: String((e && e.message) || e) }; }
    };

    /** Brief §1.2(c): the original closure never gives a FoS; compare its minimum admissible sigma_n with the gravity
     *  cap sigma_n <= gamma h - u of the dry and wet meshes (as the Calculator page does). */
    function originalCheck(inp) {
      const [mb, s, a] = E.hb(inp);
      const r = E.originalClosureMinimum(inp.sigma_ci, mb, s, a, inp.nu);
      const out = { found: r.found, sigma3_min: r.sigma3_min, sn_min: r.sn_min, phi_deg: r.phi_deg, nu: inp.nu, a, n3: inp.n3 };
      try {
        const md = E.buildMesh(E.inputs(Object.assign({}, inp, { wet: false })), 3, inp.n3);
        const mw = E.buildMesh(E.inputs(Object.assign({}, inp, { wet: true })), 3, inp.n3);
        const cd = E.gravityCap(md, inp.gamma), cw = E.gravityCap(mw, inp.gamma);
        out.ny = md.ny;
        const max = (v) => v.reduce((m, x) => (x > m ? x : m), -Infinity);
        Object.assign(out, { cap_dry: max(cd), cap_wet: max(cw), count_dry: md.count, count_wet: mw.count,
          below_dry: cd.filter((v) => v < r.sn_min * 1000).length, below_wet: cw.filter((v) => v < r.sn_min * 1000).length });
      } catch (e) { out.mesh_error = String(e.message || e); }
      return out;
    }

    function handle(job, progress) {
      const t0 = Date.now();
      const say = progress || (() => {});
      if (job.type === "run") {
        const inp = E.inputs(job.inputs);
        const a3 = attempt(() => E.run3d(inp)), a2 = attempt(() => E.run2d(inp));
        const out = { type: "run", id: job.id, key: job.key, r3: a3.r || null, e3: a3.error || null, r2: a2.r || null,
          e2: a2.error || null };
        if (inp.model === "ejc_original") {
          const oc = attempt(() => originalCheck(inp));
          out.orig = oc.r || null; out.orig_error = oc.error || null;
        }
        out.ms = Date.now() - t0;
        return out;
      }
      if (job.type === "nusweep") {
        const rows = [];
        job.nus.forEach((nu, q) => {
          const b = E.inputs(Object.assign({}, job.inputs, { model: "ejc_moscarda", nu }));
          const a3 = attempt(() => E.run3d(b)), a2 = attempt(() => E.run2d(b));
          const r3 = a3.r, r2 = a2.r;
          rows.push({ nu, fos3: r3 ? r3.fos : null, fos2: r2 ? r2.fos : null,
            status3: r3 ? r3.status : "INPUT ERROR", status2: r2 ? r2.status : "INPUT ERROR",
            reason3: r3 ? (r3.fos !== null ? "" : r3.reason.slice(0, 120)) : a3.error });
          say({ done: q + 1, total: job.nus.length, label: `ν = ${nu}` });
        });
        return { type: "nusweep", id: job.id, rows, ms: Date.now() - t0 };
      }
      if (job.type === "mesh") {
        const jobs = job.n3s.map((n) => [3, n]).concat(job.n2s.map((n) => [2, n]));
        const rows = [];
        jobs.forEach(([dim, n], q) => {
          // 3D: columns across y scaled with n, so the columns keep their shape (withN3)
          const b = attempt(() => (dim === 3 ? E.withN3(E.inputs(job.inputs), n) : E.inputs(Object.assign({}, job.inputs, { n2: n }))));
          const a = b.r ? attempt(() => (dim === 3 ? E.run3d(b.r) : E.run2d(b.r))) : b;
          const r = a.r, ny = dim === 3 && b.r ? E.columnsAcross(b.r, n) : null;
          rows.push({ dim: `${dim}D`, mesh: dim === 3 ? `${n}x${ny !== null ? ny : n}` : `${n} slices`, n, ny,
            elements: r ? r.mesh.count : null,
            status: r ? r.status : "INPUT ERROR", fos: r ? r.fos : null, iterations: r ? r.iterations : null,
            projections: r ? r.projections : null, reason: r ? (r.fos !== null ? "" : r.reason.slice(0, 120)) : a.error });
          say({ done: q + 1, total: jobs.length, label: `${dim}D ${n}` });
        });
        return { type: "mesh", id: job.id, rows, ms: Date.now() - t0 };
      }
      if (job.type === "search") {
        const res = E.criticalSearch(E.inputs(job.inputs), { dims: job.dims, back: job.back, front: job.front }, say);
        return Object.assign(res, { type: "search", id: job.id, ms: Date.now() - t0 });
      }
      if (job.type === "verify") {
        let step = 0;
        const rows = CHECKS.verificationRows({ slow: !!job.slow, progress: (label) => say({ done: ++step, total: job.slow ? 13 : 11, label }) });
        return { type: "verify", id: job.id, rows, summary: CHECKS.summary(rows), slow: !!job.slow, ms: Date.now() - t0 };
      }
      return { type: "error", id: job.id, error: `unknown job ${job.type}` };
    }
    return { handle, originalCheck };
  }

  const api = { make };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.EJC_JOBS = api;
})(typeof window !== "undefined" ? window : globalThis);
