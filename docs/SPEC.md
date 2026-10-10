# slope3d specification (derived line by line from the VBA)

Primary specification: `reference/EJC_Eq421.bas` (cited **EJC:line**) and `reference/Thesis_Columns.bas` (cited **COL:line**).
`reference/Thesis_Final_QA.bas` is cited **QA:...** (test ideas only). `reference/code/*.py` are the Python mirrors (**PY**).
Where the VBA and the brief or the mirror disagree, the VBA is followed and the difference is listed in section 0.

Code locations: `slope3d/strength.py` (EJC), `slope3d/geometry.py`, `slope3d/engine.py`, `slope3d/models.py` (COL).
`slope3d/_vector.py` is a numpy acceleration of the same algorithm, tested elementwise against the scalar code.

---

## 0. Differences and open questions (read this first)

| # | Finding | What slope3d does |
|---|---|---|
| **D1** | **The Moscarda "checks = 0" result in the brief (F = 1.04327747791523 3D / 0.90176054371878 2D) is not produced by the VBA engine.** It is recorded by the *worksheet formula engine* (`Native Workbook QA` row 11, QA test 8: `P5=2, P19=0, a=.509923, wet`). `Thesis_Columns` always requires an ADMISSIBLE strength for every element (COL:199-215, COL:279-285): `checks` only changes the final label (COL:309). Verified: the VBA engine and its mirror give NOT ACCEPTED, 76 of 76 / 30 of 30 failed strengths, for `checks` = 0 and 1. The workbook documents `P19 = 0` as "use converged eq. 4.21 roots even if they fail the checks" (`Interface!R19`). Experiment: accepting the first converged eq. 4.21 root of every element reproduces 1.0432774779 / 0.9017605437 (relative difference 1.5e-11 after 200 iterations). | **Not implemented, because it changes the engine.** The VBA behaviour is kept. One test is marked `xfail` and `tools/compare_with_workbook.py` lists the two values as OPEN. **Decision needed from you** (see the end of the chat summary). |
| D2 | The brief calls the Moscarda rejection test "wet, a = 0.509923". The VBA engine's Moscarda test (`Native Engine QA` row 11) is **dry, a = 0.5**; the wet a = 0.509923 case is the +1 kPa perturbation test. | Both Moscarda variants are tested (NOT ACCEPTED, F withheld). |
| D3 | Failure text. The brief gives "k failed strengths" (mirror) and "76 of 76 elements have no admissible eq. 4.21 strength". COL:283 is `"<k> of <n> elements have no admissible eq.4.21 strength; no FoS."` (no space in `eq.4.21`). | VBA text. |
| D4 | Mirror `fastphi` measures the angular residual with `ejc_check`, which also demands radicand >= 0, 0 <= phi < 90 deg and recomputes M from sigma1, sigma3 after a degree round trip. VBA uses `AngRes` (EJC:55-66): GState, M+1 != 0, \|z\| <= 1 only. | VBA (`strength.ang_res`). Effect on results < 1e-12. |
| D5 | Mirror `ejc_accept(non-array)` returns `REJECTED: DOMAIN`; VBA returns `REJECTED: INPUT/DOMAIN` (EJC:285). | VBA text. |
| D6 | pi. EJC uses `PI_ = 3.14159265358979` (15 digits, EJC:21); `Thesis_Columns` uses `3.141592653589793` (COL:10). The mirrors use `math.pi` everywhere. | Both VBA constants are kept (difference 7.6e-16 relative; no visible effect). |
| D7 | `ReadConfig` (COL:115-142) refuses tolerances looser than the defaults (tolA <= 1e-10, tolE <= 1e-8, tolF <= 1e-8), nu outside (0, 0.5), a outside (0, 1], etc. The mirror has none of this. | Implemented in `Config.input_issue()`; a refused input gives NOT ACCEPTED, "INVALID INPUT: ...", no number. You can tighten tolerances but not loosen them. |
| D8 | The geometry guard in the VBA (COL:98-113) has no finiteness test; the brief and the mirror do. | Finiteness test included. |
| D9 | `Native Workbook QA` expected F were recorded by the worksheet engine; the VBA engine's own recorded F (`Native Engine QA`, column C/E) differ from them by about 1e-12. | slope3d matches the VBA engine's recorded values to better than 1e-12 and the "expected" values to 1e-9. |
| D10 | **Workbook display anomaly.** `VBA Run!B18` (2D FoS display cell) holds the cached number **3** although its formula returns `D18` = "WITHHELD". It is a stale cached value (the 2D run is NOT ACCEPTED). A number is displayed for a withheld result. | Reported by the compare tool as OPEN. Recalculating in Excel (Ctrl+Alt+F9) should clear it. |
| D11 | Observation, not an error. For the Appendix E circle the slip surface passes **12.39 m above toe level at x = 0** (the circle meets z = 0 at x = +22.97 m, which is cut off by the `x <= 0` rule, COL:170). The mass therefore ends at a vertical face at x = 0 with no resistance. The published circle meets the toe to within 2.6 mm. | Documented as a model limit. |
| D12 | Observation. Accepted conventional-HB runs contain tension elements (sigma_n' < 0): 6 of 76 columns at 3D 20x20 dry a = 0.5 (min -50.26 kPa), all inside the HB tensile limit (-98.0 kPa). They are reported in the Summary, not hidden. | Reported (`tension`, `min_normal`). |

---

## 1. Units and constants

| # | Rule | Source |
|---|---|---|
| 1.1 | Lengths m; unit weight kN/m^3; stresses, c', u in kPa inside the engine; sigma_ci in MPa; principal stresses in MPa in the strength routines (sigma_n' converted kPa -> MPa by /1000). Angles entered/reported in degrees, radians internally. | EJC:8, COL:342 |
| 1.2 | `PI_ = 3.14159265358979` in the strength library, `3.141592653589793` in the engine. | EJC:21, COL:10 |
| 1.3 | `HI_DEG = 89.9999`, `NSCAN = 1800`. | EJC:22-23 |

## 2. Inputs and defaults

| # | Input | Default | Source |
|---|---|---|---|
| 2.1 | Circle centre xc, zc; radius R (m) | 355.7423115464064, 644.4192948645371; 725.2669975187362 | COL:467 |
| 2.2 | Sphere centre y = 0 (symmetry model; not an input) | 0 | COL:118 |
| 2.3 | Slope height H (m); slope = tan(face angle) = D16/E16 | 360; tan 50 deg | COL:120, 468 |
| 2.4 | Unit weight gamma (kN/m^3) | 27 | COL:121 |
| 2.5 | Water: waterAngle (deg), waterHeight (m), uExtra (kPa) | 0, 0, 0 (wet: 50, 360) | COL:121-122 |
| 2.6 | strength: 1 = HB/EJC, 0 = MC; MC c (kPa), phi (deg) | 1; workbook default 7.2, 30; numerical check 100, 0 | COL:123, 133 |
| 2.7 | sigma_ci (MPa), mi, GSI, D | 77.7, 10, 42, 0 | COL:124-126 |
| 2.8 | a override (blank = from GSI) | 0.5 (workbook `P14`) | COL:127-132 |
| 2.9 | closure 1/2/3; nu | 1; 0.271 | COL:123 |
| 2.10 | tolA = 1e-10 rad, tolE = 1e-8 MPa, tolF = 1e-8; checks = 1; maxIter = 2000 (1-50000) | | COL:134-135, 433 |
| 2.11 | Mesh: 3D n x n half-model 2-400 (default 20); 2D slices 2-20000 (default 40) | | COL:157-158 |

## 3. Derived rock-mass parameters (COL:125-132)

* 3.1 `mb = mi * exp((GSI - 100)/(28 - 14 D))`
* 3.2 `s = exp((GSI - 100)/(9 - 3 D))`
* 3.3 `a = 1/2 + (exp(-GSI/15) - exp(-20/3))/6` unless overridden (0.5 or 0.509923).

## 4. Closure (COL:144-150) and constants (EJC:37)

* 4.1 `sigma2' = p sigma1' + q sigma3'`. Closure 1 (thesis): p = -1, q = 1/nu. Closure 2 (Moscarda): p = q = nu. Closure 3 (conventional 2D HB comparator): p = 0, q = 1.
* 4.2 `P = p^2 - p + 1`, `R = 2pq - p - q - 1`, `S = q^2 - q + 1`.
* 4.3 Closure 3 is labelled "Conventional HB comparator (3D geometry, 2D strength)" and is never a 3D-EJC result. A **thesis-closure result** needs HB on, closure 1, checks 1, and an ACCEPTED dimension (COL:442).

## 5. Input validation (COL:136-141) and geometry guard (COL:98-113)

* 5.1 checks must be 0 or 1; strength 0/1; closure 1-3.
* 5.2 `sci > 0, mb > 0, s >= 0, 0 < a <= 1, 0 < nu < 0.5`.
* 5.3 `cMC >= 0, 0 <= phiMC < 90, uExtra >= 0, waterHeight >= 0, 0 <= waterAngle < 90`.
* 5.4 `0 < tolA <= 1e-10, 0 < tolE <= 1e-8, 0 < tolF <= 1e-8`.
* 5.5 **Geometry guard** (refuse with `UNSUPPORTED GEOMETRY: ...`, F absent, zero active elements): every input finite; `R > zc > 0`, `H > 0`, `slope > 0`, `gamma > 0`; `0 < xc < R`; `x_toe = xc - sqrt(R^2 - zc^2) >= -1e-5 R`.

## 6. Elements (COL:152-197)

* 6.1 `x0 = xc - sqrt(R^2 - zc^2)`, `x1 = xc - R`, `dx = (x0 - x1)/n`; 3D: `dy = R/n`, 2D: `dy = 1`.
* 6.2 Centres `x = x0 - (i + 1/2) dx` (i = 0..n-1); 3D `y = (j + 1/2) dy`, 2D `y = 0`. Loop order: x outer, y inner.
* 6.3 `q = R^2 - (x - xc)^2 - y^2`. Skip if `q <= 0` or `x > 0`.
* 6.4 `z_b = zc - sqrt(q)`; `h = min(H, -x slope) - z_b`; skip if `h <= 0`.
* 6.5 `A_h = dx dy`, `g = sqrt(q)/R`, `rho = sqrt(R^2 - y^2)`, `A = A_h/g`, `kappa = (xc - x)/rho`, `W = gamma h A_h`.
* 6.6 `u = uExtra + 9.81 max(0, min(waterHeight, -x tan(waterAngle)) - z_b)` when `waterHeight > 0`, else `uExtra`.
* 6.7 Initial `sigma_n' = (gamma h - u) g^2`; initial `phi = -1`.
* 6.8 3D volume = 2 * sum(h A_h) (half-model). No active element: "No active elements."

## 7. Strength: eq. 4.21 (EJC:25-215)

* 7.1 `ValidInputs`: `sci > 0, mb > 0, s >= 0, 0 < a <= 1, 0 < tolA <= 1e-6` (EJC:26).
* 7.2 **g(phi)** (EJC:29-53), sigma_n in MPa. `B = mb sn/sci + s` must be > 0.
  `H = (B/(1 + sin phi/a))^a`, `dH = -H cos phi/(1 + sin phi/a)`;
  `sig3 = sn - (sci/2) H (1 - sin phi)`, `sig1 = sn + (sci/2) H (1 + sin phi)` and their exact derivatives;
  `Q = mb sig3/sci + s` must be > 0;
  `N = 2 a mb sci Q^(2a-1) - R sig1 - 2 S sig3`; `D = 2 P sig1 + R sig3`; both with exact derivatives;
  `g = sin phi - 1 + 2D/(N + D)` (= `sin phi - (M-1)/(M+1)`, `M = N/D`); `dg = cos phi - 2 (dN D - N dD)/(N + D)^2`.
  Undefined (None) if `N + D = 0` or `D = 0`.
* 7.3 **Unclipped angular residual** `phi - asin((M-1)/(M+1))`; undefined if `M + 1 = 0` or `|(M-1)/(M+1)| > 1` (EJC:55-66).
* 7.4 **phi_b** (EJC:68-87): `T = s sci/mb`, `K = sn + T` (fail if <= 0); `Z = 2K/(sci (mb K/sci)^a)`; if `Z >= 1` -> 0; start `p = (1+Z^2) - sqrt((1+Z^2)^2 - (1-Z^2))`; 8 Newton steps on `(1-p)(1 + p/a)^(-a) - Z`, clamped to `[0, 1 - 1e-15]`; `phi_b = asin p`.
* 7.5 **Scan** (EJC:89-135): from `max(phi_b + 0.001 deg, 0.01 deg)` to `89.9999 deg`, 1800 intervals; sign change (`g_prev g <= 0`) only between points where g is defined; 30 bisections of each bracket; discard if g becomes undefined or if both end values still have `|g| > 1e-3` (pole).
* 7.6 **Root** (EJC:137-165): `tol <= 0`, `rootIndex < 1` or invalid inputs -> `#VALUE!`; fewer brackets than the index -> `#N/A`; undefined g -> `#NUM!`. Newton from the bracket mid-point, bisection when the step leaves (lo, hi); 200 iterations. **Converged** only when the unclipped residual `<= tol` and the last step `<= tol`; otherwise `#NUM!`. No fallback.
* 7.7 **Fast path** `EJC_FASTPHI` (EJC:167-208), closure 3 only: bracket `[phi_b + 1e-9, 89.9999 deg]`; warm start from the previous phi if inside; converged when residual `<= tol` and `|g/dg| <= tol` (before a step) or residual and step `<= tol` (after); any failure -> the full scan. Other closures call `EJC_PHI` directly.
* 7.8 In the engine, `solverTol = tolA`, and for closure 3 `min(1e-12, tolA/10)` (COL:204-205).

## 8. State checks and acceptance (EJC:217-305)

* 8.1 **Force-used state** (EJC:270-281): `tau = c + sn tan phi`, `centre = sn + tau tan phi`, `radius = tau/cos phi`, `sig1/3 = (centre +/- radius)/1000` (MPa), `sig2 = p sig1 + q sig3`. phi outside `[0, 90)` deg -> controlled domain failure.
* 8.2 **StateChecks** (EJC:224-254): closure residual `sig2 - (p sig1 + q sig3)`; ordering OK if `sig1 >= sig2 - 1e-9` and `sig2 >= sig3 - 1e-9`; domain = 1 if `hbQ = mb sig3/sci + s > 0`, `denM = 2P sig1 + R sig3 != 0`, `M + 1 != 0`, `|z| <= 1`, `rad = 4 sci^2 hbQ^(2a) - 3 (sig2 - sig3)^2 >= 0`; then `angular = phi - asin z`, `EJC residual = sig1 - (sig2 + sig3 + sqrt(rad))/2` (original, unsquared); branch = 1 if `0 <= phi < pi/2` and `asin z >= -tolA`.
* 8.3 **AcceptText** (EJC:283-293), first failing test wins, in this order: invalid tolerances -> `REJECTED: INPUT` (`0 < tolA <= 1e-6`, `0 < tolE <= 1e-4`); not an array -> `REJECTED: INPUT/DOMAIN`; domain != 1 -> `DOMAIN`; branch != 1 -> `BRANCH`; `|angular| > tolA` -> `TANGENCY`; `|EJC| > tolE` -> `EJC`; `|closure| > tolE` -> `CLOSURE`; ordering != 1 -> `ORDERING`; else `ADMISSIBLE`.
* 8.4 **Instantaneous cohesion** (COL:212): `c' = 1000 sci H cos phi/2 - sn tan phi` kPa, with `H = (b/(1 + sin phi/a))^a`, `b = mb sn/(1000 sci) + s`; `b <= 0` -> `OUTSIDE HB DOMAIN`; error from the solver -> `NO CONVERGED EQ.4.21 ROOT`.
* 8.5 `SolveStrength` is true only if the label is `ADMISSIBLE`. Mohr-Coulomb: `phi = phiMC`, `c = cMC`, always true (COL:201-203).

## 9. Iteration (COL:260-320)

* 9.1 Start `F = 1`.
* 9.2 **Tensile projection** (HB only): `limit = -1000 s sci/mb` kPa. Any `sn <= limit` becomes `limit + max(1e-9, 0.001 |limit|)`; each projection is counted (initial elements and every update).
* 9.3 Each iteration: (a) strength of every element; `k` failures -> NOT ACCEPTED, reason of D3, F withheld; (b) `drive = sum(W kappa)` must be > 0; (c) `t = tan(phi)`, `m_alpha = g + kappa t/F` must be > 0; (d) `SF = (c A_h + (W - u A_h) t)/m_alpha`, `F_new = sum(SF)/sum(W kappa)` must be > 0; (e) `sn = (W - u A_h - c A kappa/F)/(A m_alpha)` with the projection rule; track the largest `|sn change|`.
* 9.4 **Converged** when `|F_new - F|/F_new <= min(1e-12, tolF/10)` and the largest sn change `<= 1e-9` kPa. Otherwise `F = F_new`. Iteration limit -> NOT ACCEPTED.
* 9.5 Sums are accumulated left to right in element order (`engine.seqsum`).

## 10. Final-state validation (COL:217-258)

On the **actual final state** at the converged F (the stress is rebuilt at F_new, not reused from the update):

* 10.1 `final_n = (W - u A_h - c A kappa/F)/(A m_alpha)`; `tau = c + final_n tan phi`.
* 10.2 Per-element equilibrium `|W - ((n + u) A_h + tau A kappa/F)|/W <= 1e-9`, and `tau >= 0`, else `REJECTED: EQUILIBRIUM/SHEAR`. `m_alpha <= 0`, `g <= 0` or `W <= 0` -> `INVALID EQUILIBRIUM`, not evaluated.
* 10.3 HB: `EJC_ACCEPT_STATE` on `(final_n, c, phi)` for every element; any label other than ADMISSIBLE replaces the equilibrium label.
* 10.4 Global: `|sum(tau A)/sum(W kappa) - F|/F <= tolF`.
* 10.5 Diagnostic only: spherical moment form `sum(tau A rho)/sum(W kappa rho)`.
* 10.6 Statistics: max |angular|, max |EJC|, max equilibrium residual, minimum final normal stress, number of tension elements (`n < 0`).

## 11. Result rules (non-negotiable)

* 11.1 **ACCEPTED** only if every check of section 10 passes and `checks = 1`. With `checks = 0` the result is **DIAGNOSTIC ONLY** (COL:309) and F is withheld.
* 11.2 A factor of safety exists only for an ACCEPTED dimension. Otherwise "WITHHELD" with the reason. `ModelResult.fos` is `None` unless ACCEPTED; `internal_f` exists for diagnostics and tests and is never displayed or exported for a non-accepted run (iteration logs show WITHHELD).
* 11.3 3D/2D ratio only when both are ACCEPTED.
* 11.4 No fallback strength, no clipped asin, no substitute C1.
* 11.5 Unsupported geometry: NOT ACCEPTED, F absent, 0 active elements, not evaluated (COL:511-514).

## 12. Tests that realise this specification

| Suite | File | Content |
|---|---|---|
| A | `tests/test_strength.py` | T1-T22 (EJC:342-447) |
| B | `tests/test_engine.py` | six benchmark pairs (1e-9, and the VBA engine's own recorded values to 1e-12), Chile thesis closure (76/30), Moscarda, checks = 0, +1 kPa perturbation, geometry refusals, guard rules, mesh limits |
| - | `tests/test_vector_equivalence.py` | numpy path = scalar path |
| C (slow) | `tests/test_slow.py` | 95 mesh-study records (1e-9), nu thresholds 52.966376 / 52.283035 MPa and the whole `nu_exact.txt` table |
| D | `tools/compare_with_workbook.py` | `VBA Validation`, `Native Workbook QA`, `Native Engine QA`, `Refined Cases`, `VBA Columns`/`VBA Slices` |
| - | `slope3d validate [--slow]` | the same checks without pytest |
