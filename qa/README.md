# qa: QA/QC of the EJC engines

Independent checks of `webapp/engine.js` and `ejc_app/engine.py` against the thesis (Mahat 2026, GENG5512), exact
solutions, randomised tests and published case studies. They produced the QA/QC record of 27 October 2026,
`qa-report-2026-10-27.html` in this folder (also published privately on claude.ai at
<https://claude.ai/artifact/Q36P1e5CrdfqKN23MLc36T>). Open the HTML file in a browser; its numbers are those of the
run before the fixes below.

Results go to `qa/out/`, which git ignores. Paths are relative to this folder, so the scripts run from anywhere.
They need Node 18 or newer (run with 22.22) and Python 3.11 or newer with numpy, scipy and mpmath
(`pip install numpy scipy mpmath`), plus `ejc_app/requirements.txt` for the Python engine.

| Run | What it checks | Time |
|---|---|---|
| `node qa_thesis.js` | Recomputes every thesis value the app can compute, with the web engine. | 3 min |
| `python3 qa_extra.py` | Table 5.2 at ν = 0.49 at 50 digits, the tangent shear-strength gain by σ′n, the root order of the squared EJC, and an independent rebuild of Route 1 (eqs 4.8, 4.9, 4.5) for Table 5.1. | 1 min |
| `python3 compare_thesis.py` | Compares 193 thesis values with the two runs above: a value matches when it lies within half a unit of its last printed digit; counts and statuses must match exactly. Needs both runs above. | seconds |
| `python3 qa_analytic.py` | Exact φ = 0 solutions by adaptive quadrature (Bishop, F_M and the F_H limit) for three surfaces, Taylor's critical toe circles, and metamorphic checks (scaling, pore pressure at φ = 0, linearity in c). | 3 min |
| `python3 qa_random.py` | 1,600 random runs (realistic and extreme inputs, 3D and 2D): JavaScript against Python, an independent re-check of every accepted result with numpy, and local strength against 50-digit solutions. | 10 min |
| `python3 find_engine_only.py` | Replays `qa_random.py`'s random stream and prints the one local state that only the engine solved, re-evaluated at 50 digits. | 3 min |
| `node qa_industry.js` | D'Esterre (2024) against Slide2 and Slide3, Law (2024) against Slide2 and Excel, Carranza-Torres & Hormazábal (2020). | 1 min |
| `node qa_desterre_search.js` | The app's own critical circle near each D'Esterre circle (finding B1). | 1 min |
| `node qa_desterre_oms.js` | The Ordinary method of slices on the app's slices for the D'Esterre circles (finding B1). | seconds |
| `node qa_search.js` | Critical-surface searches on the thesis slope, GHB and EJC, dry and wet, 2D and 3D (finding R1). | 3 min |
| `node qa_scenarios.js` | One input at a time on the thesis slope: GSI, σci, mi, D, γ, Hw, r_u and slope height. | 1 min |

`run_engine.js` is the batch runner the Python scripts use (`node run_engine.js jobs.json out.json`).

## Findings fixed since the record

- **S1**: the original-closure search (`originalClosureMinimum` / `original_closure_minimum`) now keeps the best
  admissible point of its final bracket. Before, at ν = 0.49 it returned 59,953.87 MPa, one grid step inside the
  edge of the admissible set; it now returns 59,951.6628 MPa, the 50-digit value, and agrees with the 50-digit edge
  minima at every ν of Table 5.2 to about 1e-13 (tests in `ejc_app/tests/test_strength.py` and
  `webapp/tools/verify_engine.js`). With this fix `compare_thesis.py` finds 184 of the 193 values matching instead of
  183; the other 9 are the thesis items of the record.
- **S2**: the default grid's "require R > z_c" error now names the base-circle option, which solves face circles.
- **S3**: the results summary counts the bases with tensile effective normal stress (σ′n < 0), with a note when
  there are any and a warning above a quarter of the bases.

The record's other findings concern the thesis text (T1 to T4), published benchmarks (B1 to B3) and real-life use
(R1 to R3); they need no change to the code.
