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
| `node qa_search.js` | Critical-surface searches on the thesis slope, GHB and EJC, dry and wet, 2D and 3D (finding R1; the apps now have their own search, below). | 3 min |
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

## Finding addressed since the record

- **R1**: both apps now search for the critical surface (Studies: `criticalSearch` in `webapp/engine.js`,
  `critical_search` in `ejc_app/engine.py`, the same algorithm with the same results; method and limits in
  `ejc_app/README.md`). On the thesis slope at the default meshes its minima are at or below those of
  `qa_search.js`, which the record used:

  | Case | 2D: app search (record) | 3D: app search (record) |
  |---|---|---|
  | EJC dry | 2.1186 (2.1212) | 2.2600 (2.2635) |
  | EJC wet | 0.9929 (0.9943) | 1.0770 (1.0821) |
  | GHB dry | 2.0002 (2.0025) | 2.1359 (2.1368) |
  | GHB wet | 0.9369 (0.9373) | 1.0214 (1.0224) |

  Building it showed three things that the record did not cover:
  - The engines accept a circle whose arc ends below the ground at its back extreme x_c − R. The mass then has a
    vertical back face that carries no force, as a tension crack would. On Taylor's 75° slope (H = 10 m, c = 100 kPa,
    φ = 0) the circle (0.40, 7.20, 7.21) has a 2.8 m face and an exact FoS of 1.6245 (by `qa_analytic.py`'s
    quadrature), below the 1.6906 of the critical circle. The search tries only circles that meet the ground at both
    exits, and the results summary of both apps warns about such a face (its FoS is unchanged).
  - A search on a coarse mesh favours surfaces whose FoS the mesh underestimates: on that slope with 40 slices it
    reports 1.6925 on a circle that gives 1.7503 on 160 slices. The search re-solves its critical surfaces on a finer
    mesh and warns above the 0.5 % threshold of eq. 4.17. The thesis's wet cases exceed it at the default meshes
    (+0.8 % to +1.8 %), so their critical FoS needs a finer search mesh.
  - On the default grid, slices are spread from x_c − R to the toe-level crossing, so a nearly flat circle can leave
    its mass on one slice: on Taylor's 60° slope, R = 994 m gives 1.34 with 40 slices (2.55 with 640). The search
    skips circles whose mass covers less than half of that span.

The record's other findings concern the thesis text (T1 to T4), published benchmarks (B1 to B3) and real-life use
(R2, R3); they need no change to the code.
