# Verification of the engine

Generated with `python verification.py --slow` (the same table is on the app's Verification page).
Default path with the Section 4.6 projection; the native v15 iteration path only where a row says so.
Main sphere x_c = 355.7423115, z_c = 644.4192949, R = 725.2669975 m; 3D 20x20 and 2D 40 slices unless stated.
Tolerances: FoS given to 9 decimals to 1e-6; values given to 3-6 decimals to rounding; fixtures to the decimals
shown. INFO rows are reference values, not targets.

**Result: {'PASS': 78, 'FAIL': 0, 'INFO': 12}**

| Group | Item | Target | Computed | Difference | Tolerance | Status | Note |
|---|---|---|---|---|---|---|---|
| 1.1 Rock mass (GSI 42, m_i 10, D 0) | m_b | 1.260056 | 1.260056450 | +4.50e-07 | rounding to 6 dp | **PASS** |  |
| 1.1 Rock mass (GSI 42, m_i 10, D 0) | s | 0.0015893 | 0.0015893273 | +2.73e-08 | rounding to 7 dp | **PASS** |  |
| 1.1 Rock mass (GSI 42, m_i 10, D 0) | a from GSI | 0.509923 | 0.509922905 | -9.52e-08 | rounding to 6 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 50 kPa: phi_i (deg) | 67.45 | 67.44917 | -8.29e-04 | rounding to 2 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 50 kPa: c_i (kPa) | 362.5 | 362.4708 | -2.92e-02 | rounding to 1 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 50 kPa: tau_s (kPa) | 482.9 | 482.8792 | -2.08e-02 | rounding to 1 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 400 kPa: phi_i (deg) | 59.90 | 59.89885 | -1.15e-03 | rounding to 2 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 400 kPa: c_i (kPa) | 487.2 | 487.1667 | -3.33e-02 | rounding to 1 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 400 kPa: tau_s (kPa) | 1177.2 | 1177.1712 | -2.88e-02 | rounding to 1 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 1600 kPa: phi_i (deg) | 50.21 | 50.20528 | -4.72e-03 | rounding to 2 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 1600 kPa: c_i (kPa) | 935.3 | 935.2630 | -3.70e-02 | rounding to 1 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | sigma_n = 1600 kPa: tau_s (kPa) | 2856.0 | 2856.0026 | +2.61e-03 | rounding to 1 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | GHB, sigma_n = 400 kPa: phi_i (deg) | 58.45 | 58.44835 | -1.65e-03 | rounding to 2 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | GHB, sigma_n = 400 kPa: c_i (kPa) | 460 | 459.504 | -4.96e-01 | rounding to 0 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | GHB, sigma_n = 400 kPa: tau_s (kPa) | 1111 | 1110.925 | -7.47e-02 | rounding to 0 dp | **PASS** |  |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | Extra fixture, original closure (nu = 0.05): phi_i (deg) | 55.0720 | 55.0720172 | +1.72e-05 | rounding to 4 dp | **PASS** | sigma_ci 120 MPa, m_b 3.5, s 0.002, a 0.505, sigma_n 5.7403 MPa |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | Extra fixture, original closure (nu = 0.05): c_i (kPa) | 2079.5 | 2079.5128 | +1.28e-02 | rounding to 1 dp | **PASS** | sigma_ci 120 MPa, m_b 3.5, s 0.002, a 0.505, sigma_n 5.7403 MPa |
| 3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2) | Extra fixture with the Moscarda closure |  | REJECTED: no admissible state at this sigma_n (outside the admissible range of the closure) |  | reference only | **INFO** | the Moscarda closure has no admissible state here |
| 3 EJC route at sigma2 = sigma3 reproduces GHB | 3D dry, a = 0.5: F(EJC, sigma2 = sigma3) vs F(GHB) | 2.553668420403 | 2.553668420403 | 0.00e+00 (rel) | rel <= 1e-11 | **PASS** |  |
| 3 EJC route at sigma2 = sigma3 reproduces GHB | 2D dry, a = 0.5: F(EJC, sigma2 = sigma3) vs F(GHB) | 2.310255027833 | 2.310255027833 | 1.92e-16 (rel) | rel <= 1e-11 | **PASS** |  |
| 3 EJC route at sigma2 = sigma3 reproduces GHB | 3D wet, a = 0.5: F(EJC, sigma2 = sigma3) vs F(GHB) | 1.200121845050 | 1.200121845050 | 1.85e-16 (rel) | rel <= 1e-11 | **PASS** |  |
| 3 EJC route at sigma2 = sigma3 reproduces GHB | 2D wet, a = 0.5: F(EJC, sigma2 = sigma3) vs F(GHB) | 1.035870832725 | 1.035870832725 | 0.00e+00 (rel) | rel <= 1e-11 | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Active half-columns (3D 20x20) | 76 | 76 | 0 | exact | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Active slices (2D 40) | 30 | 30 | 0 | exact | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.5: F 3D | 2.706024342 | 2.706024342432 | +4.32e-10 | abs <= 1e-06 | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.5: F 2D | 2.447910895 | 2.447910895515 | +5.15e-10 | abs <= 1e-06 | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.5: F 3D | 1.272331405 | 1.272331405186 | +1.86e-10 | abs <= 1e-06 | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.5: F 2D | 1.098069496 | 1.098069496535 | +5.35e-10 | abs <= 1e-06 | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry / wet, a = 0.5: 3D iterations and projections (default path) | shown, not tested | 30 / 65 iterations; 0 / 0 projections |  | reference only | **INFO** |  |
| 3 Native v15 iteration path (start F = 1, Hungr explicit form, projection) | Dry, a = 0.5: F 3D | 2.706024342 | 2.706024342432 | +4.32e-10 | abs <= 1e-06 | **PASS** |  |
| 3 Native v15 iteration path (start F = 1, Hungr explicit form, projection) | Dry, a = 0.5: 3D iterations | 34 | 34 | 0 | exact | **PASS** |  |
| 3 Native v15 iteration path (start F = 1, Hungr explicit form, projection) | Wet, a = 0.5: F 3D | 1.272331405 | 1.272331405181 | +1.81e-10 | abs <= 1e-06 | **PASS** |  |
| 3 Native v15 iteration path (start F = 1, Hungr explicit form, projection) | Wet, a = 0.5: 3D iterations | 63 | 63 | 0 | exact | **PASS** |  |
| 3 Native v15 iteration path (start F = 1, Hungr explicit form, projection) | Projections dry / wet |  | 1 / 2 |  | reference only | **INFO** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.5: GHB F 3D | 2.554 | 2.553668 | -3.32e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.5: GHB F 2D | 2.310 | 2.310255 | +2.55e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.5: GHB F 3D | 1.200 | 1.200122 | +1.22e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.5: GHB F 2D | 1.036 | 1.035871 | -1.29e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.509923: F 3D | 2.624 | 2.623905 | -9.51e-05 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.509923: F 2D | 2.376 | 2.375733 | -2.67e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.509923: GHB F 3D | 2.476 | 2.476038 | +3.81e-05 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Dry, a = 0.509923: GHB F 2D | 2.242 | 2.242063 | +6.33e-05 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.509923: F 3D | 1.228 | 1.228328 | +3.28e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.509923: F 2D | 1.061 | 1.060695 | -3.05e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.509923: GHB F 3D | 1.158 | 1.158427 | +4.27e-04 | rounding to 3 dp | **PASS** |  |
| 3 Main sphere, Moscarda closure, EJC (20x20, 40 slices) | Wet, a = 0.509923: GHB F 2D | 1.000 | 1.000466 | +4.66e-04 | rounding to 3 dp | **PASS** |  |
| 3 Moment form F_M (20x20, own iteration) | Dry: F_M 3D | 2.7022 | 2.7021893 | -1.07e-05 | rounding to 4 dp | **PASS** |  |
| 3 Moment form F_M (20x20, own iteration) | Wet: F_M 3D | 1.2695 | 1.2694781 | -2.19e-05 | rounding to 4 dp | **PASS** |  |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Dry: C0 3D | 2.637 | 2.637135 | +1.35e-04 | rounding to 3 dp | **PASS** |  |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Dry: T0 3D | 2.597 | 2.597232 | +2.32e-04 | rounding to 3 dp | **PASS** |  |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Dry: C0 2D | 2.430 | 2.430370 | +3.70e-04 | rounding to 3 dp | **PASS** |  |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Dry: T0 2D | 2.429 | 2.428723 | -2.77e-04 | rounding to 3 dp | **PASS** |  |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Wet: C0 3D | 0.964 | 0.963803 | -1.97e-04 | rounding to 3 dp | **PASS** |  |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Wet: T0 3D | no accepted FoS | REJECTED | n/a | status | **PASS** | Non-positive or non-finite trial F at iteration 9. |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Wet: C0 2D | 0.978 | 0.978391 | +3.91e-04 | rounding to 3 dp | **PASS** |  |
| 3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law) | Wet: T0 2D | no accepted FoS | REJECTED | n/a | status | **PASS** | Non-positive or non-finite trial F at iteration 16. |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | nu = 0.271: sigma3,min (MPa) | 34.376 | 34.375976 | -2.41e-05 | rounding to 3 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | nu = 0.271: sigma_n,min (MPa) | 52.97 | 52.96638 | -3.62e-03 | rounding to 2 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Gravity cap max(gamma h - u), dry (MPa, 3D 20x20) | 1.61 | 1.60817 | -1.83e-03 | rounding to 2 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Gravity cap max(gamma h - u), wet (MPa, 3D 20x20) | 1.02 | 1.02387 | +3.87e-03 | rounding to 2 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Table 5.2: min sigma_n at nu = 0.05 (MPa) | 0.92 | 0.91586 | -4.14e-03 | rounding to 2 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Table 5.2: min sigma_n at nu = 0.1 (MPa) | 3.43 | 3.42933 | -6.72e-04 | rounding to 2 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Table 5.2: min sigma_n at nu = 0.2 (MPa) | 19.23 | 19.22560 | -4.40e-03 | rounding to 2 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Table 5.2: min sigma_n at nu = 0.271 (MPa) | 52.97 | 52.96638 | -3.62e-03 | rounding to 2 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Table 5.2: min sigma_n at nu = 0.35 (MPa) | 177.8 | 177.8162 | +1.62e-02 | rounding to 1 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Table 5.2: min sigma_n at nu = 0.45 (MPa) | 2188 | 2187.810 | -1.90e-01 | rounding to 0 dp | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | nu = 0.05: dry columns with gamma h - u below the threshold (of 76) | 43 | 43 | 0 | exact | **PASS** |  |
| 1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5) | Original closure never produces a FoS | NOT APPLICABLE | NOT APPLICABLE | n/a | status | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 2D GHB, 40 slices (brief) | 2.0039 | 2.0038734 | -2.66e-05 | rounding to 4 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 2D GHB, 40 slices (thesis record) | 2.003873 | 2.003873411 | +4.11e-07 | rounding to 6 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 3D GHB, 20x20 (thesis value) | 2.142 | 2.141937 | -6.26e-05 | rounding to 3 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 3D GHB, 20x20 (thesis record) | 2.141937 | 2.141937395 | +3.95e-07 | rounding to 6 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 2D EJC, 40 slices | 2.122598 | 2.122598022 | +2.22e-08 | rounding to 6 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 3D EJC, 20x20 | 2.269068 | 2.269068429 | +4.29e-07 | rounding to 6 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | Published FoS (their method) | 2.01 |  |  | reference only | **INFO** | reference only |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 2D EJC, 5120 slices (slow) | 2.121420 | 2.121420104 | +1.04e-07 | rounding to 6 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 3D GHB, 320x320 (slow) | 2.140776 | 2.140775815 | -1.85e-07 | rounding to 6 dp | **PASS** |  |
| 3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5 | 3D EJC, 320x320 (slow) | 2.267835 | 2.267835497 | +4.97e-07 | rounding to 6 dp | **PASS** |  |
| 3 Moscarda (2025), main sphere, nu = 0.271 (reference values, not targets) | Case 1: dry, a = 0.5, 2D FoS | reported 2.24 (sheet at nu = 0.271: 2.245) | this engine (Route 2, 40 slices): 2.4479 |  | reference only | **INFO** | his own slice model and the inherited route |
| 3 Moscarda (2025), main sphere, nu = 0.271 (reference values, not targets) | Case 2: wet, a = 0.5, 2D FoS | reported 1.57 (sheet at nu = 0.271: 1.567) | this engine (Route 2, 40 slices): 1.0981 |  | reference only | **INFO** | his own slice model and the inherited route |
| 3 Moscarda (2025), main sphere, nu = 0.271 (reference values, not targets) | Case 3: dry, a = 0.509923, 2D FoS | reported 2.04 (sheet at nu = 0.271: 2.041) | this engine (Route 2, 40 slices): 2.3757 |  | reference only | **INFO** | his own slice model and the inherited route |
| 3 Moscarda (2025), main sphere, nu = 0.271 (reference values, not targets) | Case 4: wet, a = 0.509923, 2D FoS | reported 1.13 (sheet at nu = 0.271: 1.132) | this engine (Route 2, 40 slices): 1.0607 |  | reference only | **INFO** | his own slice model and the inherited route |
| 6 Limitations statements (about) | F_H vs F_M difference, dry / wet (%) | 0.14-0.27 | 0.142 / 0.224 |  | reference only | **INFO** |  |
| 6 Limitations statements (about) | Wet 3D FoS drop with C0 (%) | about 25 | 24.2 |  | reference only | **INFO** |  |
| 3 Finest mesh (240x240, 5120 slices; slow) | Dry: F 3D 240x240 | 2.7112 | 2.7112235 | +2.35e-05 | rounding to 4 dp | **PASS** |  |
| 3 Finest mesh (240x240, 5120 slices; slow) | Dry: F 2D 5120 | 2.4422 | 2.4422059 | +5.87e-06 | rounding to 4 dp | **PASS** |  |
| 3 Finest mesh (240x240, 5120 slices; slow) | Wet: F 3D 240x240 | 1.2776 | 1.2775700 | -3.00e-05 | rounding to 4 dp | **PASS** | 3 projections on the path |
| 3 Finest mesh (240x240, 5120 slices; slow) | Wet: F 2D 5120 | 1.0916 | 1.0916067 | +6.70e-06 | rounding to 4 dp | **PASS** | 7 projections on the path |
| 6 Limitations statements (about) | 10x10 below 240x240, dry (%) | about 2.6 | 2.55 |  | reference only | **INFO** |  |
| 6 Limitations statements (about) | 10x10 below 240x240, wet (%) | about 4.4 | 4.43 |  | reference only | **INFO** |  |
