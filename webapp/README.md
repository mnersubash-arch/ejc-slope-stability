# webapp: EJC Slope Stability in the browser

The `ejc_app` analysis as a single web page: Hungr 3D columns and Bishop 2D slices on a prescribed spherical slip
surface, with the Extended Jiang Criterion (Route 2), GHB and Mohr-Coulomb, and the nine additions described under
"Additions" in `ejc_app/README.md`. The engine is a JavaScript port of
`ejc_app/engine.py`, checked against it (below). Published privately on claude.ai at
<https://claude.ai/artifact/HaR96PNS4QgJLWSewwrfGg>; only the owner can open it unless it is shared.

Research tool for the prescribed-surface comparison in Mahat (2026) thesis; not for design.

## What the page does

- **Inputs:** a number box and a slider for every input, the 56 presets of `ejc_app/presets.py`, and a "Reset inputs
  to preset" button. The number box takes any value in the valid range; the slider spans a window around the preset's
  value (for example ±H for the sphere centre and radius), re-centred when a typed value falls outside it.
- **3D columns:** N × N (the brief's grid, the default), N × N_y, or a plan size Δx × Δy in metres that becomes
  counts on the current circle (N = ⌈(x0 − x1)/Δx⌉, N_y = ⌈R/Δy⌉, each 2 to 320). The rail shows each column's plan
  size and the 2D slice width; the summary gives the grid, Δx × Δy and the range of column heights h, and the 3D view
  and the free-body diagram give the selected column's Δx × Δy × h. The column height is set by the ground and the
  slip surface. See "Column shape" in `ejc_app/README.md`.
- **Additions** (all off by default; formulations and sources in `ejc_app/README.md`):
  - Slope: a planar face, benches or typed points; ground rising behind the crest; a wall concave or convex in plan.
  - Slip surface: ellipsoidal end caps (λ) around a cylindrical middle (L).
  - Strength: D only within T of the face.
  - Water: a piezometric line through points.
  - Loads: a tension crack with water, and a seismic coefficient k_h.

  The plan radius, T and the crack depth have no default: the run waits until they are entered. The views draw
  them: the ground and the slip surface in 3D, the crack and its water on the cross-section, and k_h W and P in the
  free-body diagrams.
- **Updating:** results update live up to 40 × 40 columns and N_2D ≤ 640; above that, press Run. The engine runs in
  a Web Worker (`worker.js`), or on the page itself if the worker cannot start.
- **Summary:** both factors of safety and their ratio, iterations, projections and the bases with tensile effective
  normal stress (σ′n < 0), with a note when there are any and a warning above a quarter of the bases (thesis §5.8).
  A warning appears when the slip circle ends below the ground at its back extreme (or, with the base-circle option,
  its front extreme): the mass then has a vertical face on which no force acts (`verticalFaces`).
- **Views:** 3D view (columns coloured by σn, φi, ci, τs, u, τs A or the tension flag, with a mirror toggle and a plan
  view to pick a column); cross-section (slices with hover values; a click opens the slice's free-body diagram);
  free-body diagram (every force on the sliding mass, on the 2D central section or the 3D full mass: the weight W of
  each slice or column and the N′, U and S on its base, the force polygon, the totals with their moments about the
  rotation axis, and the vertical, horizontal and moment balance; then one 3D column with W, N′, U and S, the lever
  arm ρ to the rotation axis, the force table and the vertical balance of eq. 4.12, and the matching 2D slice); final-state checks with the iteration history and a
  copy-as-CSV results table; studies; the verification table.
- **Studies:**
  - the critical-surface search (below), the mesh study (eq. 4.17) and the ν sweep;
  - EJC against GHB, in 3D and 2D, with the ratios;
  - a sweep of any input in use, for EJC and GHB, in 3D and 2D;
  - a back-analysis: the value of one input that gives a target FoS, by bisection;
  - Monte Carlo: truncated normal inputs with your standard deviations and a seed; the result gives P_f with its
    standard error (or the rule-of-three bound when no sample fails), β and the percentiles.
- **Cases and report** (the box under the preset):
  - save the case as JSON, open a case file, or copy and paste a case code; the Streamlit app reads and writes the
    same format;
  - where the page runs as an artifact with storage, its owner can share named cases, which everyone who can open
    the page can open;
  - **Report** gives one self-contained HTML page: the FoS, the cross-section, the notes, the inputs in use and the
    studies run on them. It needs results for the inputs on screen. Where files cannot be saved, it opens in a
    window.
- **Critical-surface search** (QA/QC finding R1; `criticalSearch`, the same algorithm as `critical_search` in
  `ejc_app/engine.py`, described in `ejc_app/README.md`): the circle (2D) and the sphere centred on y = 0 (3D) with
  the lowest accepted FoS for the current inputs, mesh and grid option, among surfaces that meet the ground at both
  exits. It shows the prescribed and critical FoS with the change, the critical circle and its exits, a mesh check on
  a finer mesh, the surfaces tried (solved, withheld, skipped), the limits the critical surface lies on, and a chart
  of the surfaces with the centres tried coloured by FoS. **Use this circle** or **Use this sphere** makes the
  critical surface the prescribed one and reruns the analysis. On the thesis slope it takes about 7 s (dry) to 17 s
  (wet) for 2D and 17 to 37 s more for 3D.
- Inputs are remembered in the browser (localStorage) only.

## Files

| File | What it is |
|---|---|
| `index.html` | The page. It is claude.ai artifact source, so it has no `<html>`, `<head>` or `<body>`: the publish step wraps it in a document skeleton (`tools/make_preview.py` does the same locally). |
| `engine.js` | Line-for-line port of `ejc_app/engine.py` (including `massForces` = `mass_forces` and `criticalSearch` = `critical_search`), plus `columnForces` from `ejc_app/plots.py`. Same equations, tolerances, σ3 scan, bisection and checks. |
| `presets.js` | Generated by `tools/make_presets.py` from `ejc_app/presets.py`; do not edit by hand. |
| `checks.js` | The rows of `ejc_app/verification.py` (same targets, tolerances and wording). |
| `jobs.js` | One analysis (3D + 2D, or the original-closure check), the mesh study, the ν sweep, the critical-surface search, EJC against GHB, the parameter sweep, the back-analysis, Monte Carlo and the verification table. |
| `worker.js` | The Web Worker that runs `jobs.js`. |
| `tools/` | The checks below. |

Plotly.js 2.35.2 loads from cdnjs (jsDelivr as a fallback); fonts from Google Fonts.

## How the port was checked

Run from this folder:

1. `python3 tools/dump_reference.py` writes the Python engine's results to `tools/ref/` (about 4 minutes).
2. `node tools/verify_engine.js` compares `engine.js` with them. Result: 0 failures. It compares:
   - 248 analyses: every main case, GHB, the σ2 = σ3 closure, F_M, C0/T0, the native v15 path, the safeguard off,
     other meshes, N × N_y columns, ν, GSI, D, r_u, base and face circles, Mohr-Coulomb, the Carranza-Torres, Law,
     D'Esterre and Baskoyak presets, the input errors, and the additions (the Baskoyak real section as a profile,
     concave and convex walls, ellipsoids and cylindrical middles, the D zone, the piezometric line, tension cracks
     with and without water, k_h with F_M, the native path and C0, all additions together, and their input errors);
   - the whole-mass forces of the 201 analyses that reach a final state;
   - 63 strength-fixture arrays, 7 original-closure minima (also checked against 50-digit values) and 24 slip
     extents;
   - 10 vertical-face cases and 17 column-shape cases (columns across y for other N, plan sizes, counts for a size);
   - 15 study checks: EJC against GHB, sweeps, back-analyses, Monte Carlo and the random stream;
   - 11 critical-surface searches, 7,732 surfaces with every field compared: Taylor's φ = 0 slopes; D'Esterre slopes
     with the base-circle option, among them two with N × N_y columns, one checked on 16 × 320, and one with k_h
     and a tension crack; the thesis slope with GHB in 2D and 3D; the Baskoyak failure with GHB in 2D and 3D, and on
     the real section; two searches that stop on a domain limit.

   The largest FoS difference is 4.4e-16 (relative), and the whole-mass totals agree within 4e-16 of the weight. The
   statuses, iteration counts, projection counts, grids, error texts and check results are the same. The searches
   try the same surfaces in the same order and find the same critical surfaces, with FoS within 1e-12 and lengths
   within 7e-15 of H. The results agree to rounding rather than bit for bit because the
   JavaScript `pow`, `tan` and `asin` can differ from NumPy's in the last bit.
3. `node tools/verify_checks.js` compares the page's verification table with `verification.py`: identical text, row
   for row (82 PASS, 0 FAIL, 38 INFO, 9 SKIPPED). `--slow` also computes the slow rows: 89 PASS, 0 FAIL, 40 INFO,
   with the values in `ejc_app/VERIFICATION.md`.
4. `node tools/ui_test.js` drives the page in Chromium (Playwright): the default case gives 2.706024 / 2.447911, the
   summary counts the tensile bases, every tab draws without console errors, the whole-mass free-body diagram
   balances (2D and 3D), the critical-surface search finds 2.118587 in 2D and Use this circle reruns the analysis on
   it, presets (Baskoyak among them) and inputs update the results, the column option gives 2.720204 with 20 × 40
   columns and 2.713457 with columns of at most 20 m by 20 m (20 × 37) and flags a size that needs more than 320
   columns, and the layout holds in dark mode and at 390 px. It also checks the additions:
   - one bench gives the planar face, four benches give 2.192610 / 2.022105, and points that run back are refused;
   - a concave wall of R_p = 800 m gives 2.755516;
   - the default piezometric line gives the water table's result;
   - a 20 m crack half full with k_h = 0.1 gives 2.339066 / 2.101695;
   - the EJC/GHB comparison and a back-analysis (D = 0.417641 for F = 2 in 2D);
   - a saved case reopens from its code, and the report holds the cross-section and the studies.

   46 checks in all.
   See the file's header for the offline option.

The QA/QC of 27 October 2026 (`qa/` at the repository root) checked the engines against the thesis, exact solutions,
randomised tests and published case studies; its findings S1 to S3 are fixed in this version and R1 is addressed by
the critical-surface search. `qa/qa_baskoyak.py` (29 October 2026) checks both apps on a real slope failure, the
Baskoyak presets (see `qa/README.md`).

## Run it locally

```bash
cd webapp
python3 tools/make_preview.py                    # writes preview.html (index.html in the publish skeleton)
python3 -m http.server 8765 --bind 127.0.0.1     # then open http://127.0.0.1:8765/preview.html
```

Opening `index.html` straight from disk does not work: the worker needs `http://`, and the page relies on the
skeleton's `[hidden]` rule. For a program that does all this on a desktop (and offline), see `desktop/README.md`.

## Changing it

- After a change to `ejc_app/engine.py`, make the same change in `engine.js` and run steps 1 to 3.
- After a change to the presets, run `python3 tools/make_presets.py`.
- To update the published page, republish `index.html` with the five `.js` files to the same artifact URL from a
  Claude session.
