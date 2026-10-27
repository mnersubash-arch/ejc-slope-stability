# Verification of slope3d against the thesis workbook

**Reference:** `GENG5512_Thesis_v14_Final.xlsm` (Excel 16.0 records of 2026-10-07: 22/22 strength fixtures, 10/10 worksheet checks,
11/11 engine checks). **Software:** slope3d 0.1.0, Python 3.13, numpy 2.5. Values are the workbook's *cached* cell values,
read with openpyxl; no Excel is needed. Reproduce this page with
`python tools/compare_with_workbook.py --summary` (about 20 s), `pytest`, `pytest -m slow` and `slope3d validate --slow`.
Relative difference = |workbook - slope3d| / |workbook|. Acceptance tolerance: 1e-9.

## 1. Every recorded value in the workbook QA sheets

| Workbook sheet | Comparisons | Reproduced | Failed | Open | Worst relative difference | Where |
|---|---|---|---|---|---|---|
| VBA Validation (22 strength fixtures) | 60 | 60 | 0 | 0 | 1.1e-10 | T4, phi = 89.95 deg |
| Native Workbook QA (9 worksheet cases) | 48 | 48 | 0 | 0 | 1.5e-10 | Moscarda diagnostic F |
| Native Engine QA (11 engine cases) | 58 | 58 | 0 | 0 | 3.1e-12 | recorded "expected" F |
| VBA Columns (76 elements x 23 fields) | 4 | 4 | 0 | 0 | 7.8e-13 | element-table cells |
| VBA Slices (30 elements x 23 fields) | 4 | 4 | 0 | 0 | 7.7e-13 | element-table cells |
| VBA Run | 5 | 4 | 0 | **1** | n.a. | see section 4 |
| Interface (slip-surface exit, F39:G39) | 2 | 2 | 0 | 0 | 1.5e-14 | exit x |
| Refined Cases (finest-mesh F) | 10 | 10 | 0 | 0 | 1.5e-12 | wet a = 0.509923, 2D |
| **All** | **191** | **190** | **0** | **1** | **1.5e-10** | |

A comparison is a recorded number (F, angle, residual, coordinate, table cell), an acceptance label, a status, a count or a
reason text. Strength-fixture residuals are recorded as 4-significant-digit text and are compared to that precision.

## 2. Where the numbers are tightest and loosest

* **Factor of safety, six benchmark cases, 3D 20x20 and 2D 40 slices (24 values):** worst 3.1e-12 against the recorded
  "expected" F (these were written by the worksheet engine) and **2.8e-15 against the VBA engine's own recorded output**
  (`Native Engine QA`, columns C and E). The Chile thesis closure is NOT ACCEPTED with 76 columns / 30 slices, as recorded.
* **Strength fixtures T1-T22:** all 22 pass; returned angles agree to 1.1e-10 or better (worst at phi = 89.95 deg, where the
  eq. 4.21 root is ill-conditioned; both codes stop at the same 1e-10 rad tolerance). All 22 acceptance labels agree.
* **Element tables of the workbook's last run** (Chile, thesis closure, NOT ACCEPTED): every cell of the 76 x 23 and 30 x 23
  tables agrees to 7.8e-13, including the "NOT EVALUATED" fields and status text.
* **Diagnostic mode** (Moscarda closure, wet, a = 0.509923, checks = 0): internal F 1.0432774778 (3D) and 0.9017605438 (2D)
  against the recorded 1.04327747791523 and 0.90176054371878: 1.5e-10 and 1.4e-10. The difference is the convergence
  criterion (|dF|/F <= 1e-10); the worksheet engine ran a fixed 34 iterations.
* **Slip-surface exit** on the face (x = -20.4057 m, 24.3186 m above the toe) agrees with the workbook's chart helper to 1.5e-14.

## 3. Independent records and thresholds

| Check | Result |
|---|---|
| 95 mesh-study records (5 cases; 3D 10x10 to 320x320; 2D 20 to 2,560 slices) | all reproduced: F worst 1.5e-12, volume worst 8.7e-15, element counts 95/95, acceptance 95/95 |
| Exact-EJC admissibility threshold at nu = 0.271 | 52.966376 MPa (a = 0.5), 52.283035 MPa (a = 0.509923), as in `nu_exact.txt`; the whole nu table (10 rows x 2) also agrees to 5e-7 MPa |
| Refinement-study finest F (`Refined Cases`) | 10 of 10, worst 1.5e-12 |

## 4. The one open item

**`VBA Run!B18`** (the "2D FoS" display cell) holds the cached number **3**, while its formula
`=IFERROR(IF(C8="CURRENT",D18,"WITHHELD"),"WITHHELD")` refers to `D18` = "WITHHELD" and the 2D run is NOT ACCEPTED. This looks like a
stale cached value in the workbook, not a calculation: slope3d returns WITHHELD. Recalculating the workbook in Excel should clear
it; until then it is recorded as OPEN and is not counted as reproduced.

## 5. Software tests

| Suite | Content | Result |
|---|---|---|
| A | 22 strength fixtures T1-T22 (assertions follow the VBA) | all pass |
| B | engine regressions, rejections, perturbation (+1 kPa), geometry refusals, mesh limits, diagnostic mode | all pass |
| equivalence | numpy path against the scalar port (random states, 3 closures) | all pass |
| report / GUI | slip-surface exit, six plots, PDF/HTML report, Streamlit app (headless), labelled diagnostic F | all pass |
| C (slow) | 95 mesh-study records, nu thresholds | all pass |
| `slope3d validate --slow` | the same regression checks without pytest | ALL PASS |

Counts at the time of writing: `pytest` 84 passed; `pytest -m slow` 5 passed; `slope3d validate --slow` 31 of 31; `tools/compare_with_workbook.py` 190 of 191 reproduced, 0 failed, 1 open.

## 6. What this does and does not show

It shows that slope3d reproduces the workbook's recorded numbers, statuses and withheld results; the VBA engine to machine
precision, the worksheet engine and the independent mesh study to 1e-9 or better. It does not validate the method against
field data or an independent 3D code; the model is a prescribed spherical surface, a symmetric half-model, Hungr's
vertical-equilibrium simplification and the positive-centre, toe-cutoff geometry. Deliberate differences from the
VBA, each tested, are listed in `docs/SPEC.md` section 0 (the diagnostic mode is the only added behaviour).
