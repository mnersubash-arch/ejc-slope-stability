# EJC 3D/2D rock-slope stability app

Interactive app for the UWA GENG5512 thesis *Developing a fully 3D limit equilibrium method, using a 3D Hoek–Brown
failure criterion, to analyse rock slope stability, and comparison to 2D methods*. It reproduces the thesis
calculation: Hungr's simplified 3D method of columns and Bishop's 2D method on a prescribed spherical slip surface,
with the Extended Jiang Criterion (EJC) and the Moscarda closure, solved by Route 2.

This folder stands alone. The `slope3d` package at the repository root is the earlier port of the thesis workbook
(VBA). The two share no code; this engine's GHB results reproduce that package's closure-3 benchmarks.

Research tool for the prescribed-surface comparison in Mahat (2026) thesis; not for design. Read the
Limitations page.

## Run

```bash
cd ejc_app
pip install -r requirements.txt
streamlit run app.py              # opens http://localhost:8501
pytest                            # about 70 s, includes the fine-mesh runs
pytest -m "not slow"              # skips the fine-mesh runs
python verification.py --slow     # the verification table as Markdown
```

The app is for local use. `.streamlit/config.toml` hides Streamlit's deploy button and turns off usage statistics.
Nothing in this folder deploys it to a public host.

## Pages

1. **Calculator.** Sidebar inputs: each has a number box and a slider that stay in sync, with units, valid range
   and reset-to-preset. Outputs:
   - FoS 3D and 2D, the 3D/2D ratio, iterations, active columns and slices, and the projection count;
   - PASS/FAIL for every final check, with the largest residuals;
   - CSV downloads of the column and slice results;
   - studies: a mesh study (eq. 4.17) and a ν sweep.

   Updates live for N ≤ 40; for larger N press **Run**. With the original closure selected, the page shows the
   admissibility check instead of a FoS.
2. **Cross-section.** Central section y = 0: ground, slip surface, water table, Bishop slices (hover for σ_n,
   φ_i, c_i, τ_s and u) and the circle centre.
3. **3D view.** Ground, spherical slip surface and column prisms, coloured by σ_n, φ_i, c_i, τ_s, u, contribution
   to resistance τ_s A, or the tension flag. A toggle mirrors the half-model into the full mass, and a plan view
   highlights a column.
4. **Free-body diagram.** Choose a column by clicking the plan view or by its index. The page shows:
   - W, N′ = σ_n A, U = u A, S = τ_s A / F, the base normal, the lever arm ρ and the rotation axis;
   - the force values and the vertical-balance residual;
   - the matching 2D slice diagram.
5. **Strength explorer.** σ1–σ3 curves for EJC + Moscarda, EJC + original closure (with its ordering band shaded)
   and GHB, plus the Mohr circle, tangent and envelope at a chosen σ_n.
6. **Limitations.**
7. **Verification.** Runs every target in one go and shows target, computed value, difference and PASS/FAIL.

## Presets

Every preset fills all the inputs and shows its source:

- **This thesis.** The main sphere (Moscarda 2025, Appendix E, mirrored in x) in four cases: dry or wet,
  a = 0.5 or 0.509923. Also: finest mesh (slow), F_M, C0/T0 (Table 5.7), the native v15 iteration path, the
  ν sweep, and the original-closure checks at ν = 0.271 and 0.05 (Table 5.2).
- **Carranza-Torres & Hormazábal (2020).** The published circle with GHB and EJC, at 20×20 / 40 slices and at
  320×320 / 5120 slices (slow).
- **Moscarda (2025), Cases 1–4.** His reported 2D FoS (2.24, 1.57, 2.04, 1.13) are shown as reference values,
  not targets. They come from his own slice model and the inherited route. ν = 0.271 is used, because his reported
  values are what his spreadsheet gives at ν = 0.271.
- **Law (2024) and D'Esterre (2024).** Empty slots; the inputs are still to be entered from the sources.

## Equations (thesis equation numbers)

Conventions: compression is positive and all stresses are effective. Strength is evaluated in MPa; column forces
are in kN, stresses in kPa and areas in m². The toe is at the origin and the slope rises towards −x:
$z_t = \min(H, -x\tan\beta)$ for $x \le 0$, and $z_t = 0$ on the pit floor.

| Eq. | Content | In `engine.py` |
|---|---|---|
| 2.5 | $m_b = m_i e^{(GSI-100)/(28-14D)}$, $s = e^{(GSI-100)/(9-3D)}$, $a = \tfrac12 + \tfrac16(e^{-GSI/15} - e^{-20/3})$; $a$ can be overridden | `hoek_brown_parameters` |
| 2.6 | Tangent geometry: $\sin\phi_i = \frac{M-1}{M+1}$, $\sigma_n = \sigma_3 + \frac{\sigma_1-\sigma_3}{M+1}$, $\tau_s = \frac{(\sigma_1-\sigma_3)\sqrt M}{M+1}$ with $M = d\sigma_1/d\sigma_3$ | `principal_states` |
| 2.7 | GHB: $\sigma_1 = \sigma_3 + \sigma_{ci}Q^a$, $Q = m_b\sigma_3/\sigma_{ci} + s$, $M = 1 + a m_b Q^{a-1}$ | `principal_states` (`ghb`) |
| 2.8 | Unsquared EJC: $f = 2\sigma_1 - \sigma_2 - \sigma_3 - \sqrt{4\sigma_{ci}^2 Q^{2a} - 3(\sigma_2-\sigma_3)^2} = 0$ (governs acceptance) | `ejc_residual` |
| 4.1 | Original closure $\sigma_2 = \sigma_3/\nu - \sigma_1$ (admissibility check only) | `Strength.pq` |
| 4.2 | Its squared form $3\sigma_1^2 - \frac3\nu\sigma_1\sigma_3 + (1 - \frac1\nu + \frac1{\nu^2})\sigma_3^2 = \sigma_{ci}^2Q^{2a}$ | `closure_constants` |
| 4.3 | Moscarda closure $\sigma_2 = \nu(\sigma_1 + \sigma_3)$ (adopted) | `Strength.pq` |
| 4.4 | Its squared form $C(\sigma_1^2 + \sigma_3^2) - B\sigma_1\sigma_3 = \sigma_{ci}^2Q^{2a}$, $C = 1 - \nu + \nu^2$, $B = 1 + 2\nu - 2\nu^2$ | `closure_constants` |
| 4.5 | Gradient $M = \frac{2a m_b\sigma_{ci}Q^{2a-1} + B\sigma_1 - 2C\sigma_3}{2C\sigma_1 - B\sigma_3}$ (the general $\sigma_2 = p\sigma_1 + q\sigma_3$ form is used, with $P = p^2-p+1$, $R = 2pq-p-q-1$, $S = q^2-q+1$) | `principal_states` |
| 4.6 | $\phi_i$, $\sigma_n$, $\tau_s$ from eq. 2.6 with $M$ of eq. 4.5 | `principal_states` |
| 4.7 | $c_i = \tau_s - \sigma_n\tan\phi_i$ | `principal_states` |
| 4.8 | Inherited GHB $\tau_s$ (Moscarda's route) | not used (Moscarda presets: reference only) |
| 4.9 | Route 1 angular residual | not used (the app uses Route 2) |
| 4.10 | Sphere $z_b = z_c - \sqrt{R^2 - (x-x_c)^2 - y^2}$ | `sphere_z`, `build_mesh` |
| 4.11 | $A_h = \Delta x\Delta y$, $g_z = \sqrt{R^2-(x-x_c)^2-y^2}/R$, $A = A_h/g_z$, $W = \gamma h A_h$, $\rho = \sqrt{R^2-y^2}$, $\kappa = (x_c - x)/\rho$ | `build_mesh` |
| 4.12 | Local vertical equilibrium $W = (\sigma_n + u)A_h + \tau_s A\kappa/F$, $m_\alpha = g_z + \kappa\tan\phi_i/F$ | `_local_sn`, final checks |
| 4.13 | $\sigma_n = (W - uA_h - c_iA\kappa/F)/(A\,m_\alpha)$ | `_local_sn` |
| 4.14 | $F_H = \sum\tau_sA / \sum W\kappa$ (governing) | `_global_f` |
| 4.15 | $F_M = \sum\tau_sA\rho / \sum W(x_c - x)$, re-solved with its own iteration | `_global_f` |
| 4.16 | Bishop 2D: $F = \sum\frac{c_ib + (W - ub)\tan\phi_i}{\cos\alpha + \sin\alpha\tan\phi_i/F} / \sum W\sin\alpha$ (unit-thickness slices, $\rho = R$, $g_z = \cos\alpha$, $\kappa = \sin\alpha$) | `run_2d`, Bishop residual check |
| 4.17 | $\varepsilon_{mesh} = (F_{max} - F_{min})/F_{finest} \times 100\%$; 0.5 % screening threshold | `mesh_error` |

Water: $u = \gamma_w\max(0, z_w - z_b)$, $z_w = \min(H_w, -x\tan\beta_w)$, $\gamma_w = 9.81$ kN/m³.

### Mesh

The half-model (y ≥ 0) uses N × N nominal cells:

- x runs from $x_0 = x_c - \sqrt{R^2 - z_c^2}$ down to $x_1 = x_c - R$, with $x_i = x_0 - (i - \tfrac12)\Delta x$;
- y runs from 0 to R in steps of $\Delta y = R/N$.

A column is active when $z_b$ is real, $x \le 0$ and $h = z_t - z_b > 0$. The 2D section uses the same x grid
with $N_{2D}$ slices.

### Route 2 local strength

For a target σ_n, σ3 is scanned over (σ3,t, σ_n] with σ3,t = −sσ_ci/m_b. The scan is geometric towards
σ3,t, and the ends of the admissible interval are located by bisection. At each σ3:

1. σ1 is the root of the squared form (eq. 4.2 or 4.4) that satisfies the unsquared EJC (eq. 2.8), Q > 0, a real
   radicand and σ1 ≥ σ2 ≥ σ3.
2. M (eq. 4.5) and σ_n(σ3) (eq. 4.6) are evaluated.
3. σ_n(σ3) − σ_n,target is bracketed and bisected to machine precision.

A local state is accepted only if:

- |f| < 1e-8 MPa;
- the closure residual is < 1e-8 MPa;
- the stresses are ordered;
- the state is real;
- M ≥ 1 is finite;
- τ_s ≥ 0.

Otherwise it is REJECTED with its reason. GHB is never substituted. The scan does not prove that the root is
unique; the number of brackets seen is reported for each base (`n_brackets` in the CSV).

### Iteration

**Default path (brief §1.4).**

1. Start from σ_n0 = (γh − u) g_z².
2. Each iteration takes the strength at σ_n, then F from F_H (or F_M), then σ_n from eq. 4.13 at that F.
3. Stop when |ΔF|/F < 1e-12 and max |Δσ_n| < 1e-9 kPa, or after 2000 iterations.
4. Check the final state independently: fresh strength at the final σ_n, relative local residual ≤ 1e-9 and
   relative global residual ≤ 1e-8. In 2D, also check Bishop's expression (eq. 4.16).

A FoS is reported only if every check passes. The default path converges in 30 iterations (dry) and 65 (wet) for
the main cases; these counts are shown, not tested.

**Projection (thesis Section 4.6), always on.**

- limit = max(−1000 s σ_ci/m_b kPa, σ_n at the admissible edge of the closure). For the Moscarda closure at
  a = 0.5 the limit is −97.916 kPa.
- A trial σ_n ≤ limit is set to limit + max(1e-9, 0.001 |limit|) and counted.
- A solution is accepted only if no base is projected at convergence, so the projection changes the path, not
  the fixed point. The projection count is reported.

**Native v15 iteration path (option).** The workbook path: start at F = 1, update σ_n with the previous F, and
take the new F from the Hungr explicit form Σ[(c A_h + (W − uA_h) tanφ)/m_α]/ΣWκ, with projection.

- It reaches the same F as the default path (to about 1e-11).
- The main cases converge in 34 (dry) and 63 (wet) iterations; these are the tested counts.
- It supports F_H only and does not offer C0/T0: the Table 5.7 values come from the default path.

### Tension sensitivities C0 / T0

- **C0** removes c_i on bases in the tensile set.
- **T0** removes all basal shear there.
- Both are sensitivities, not a calibrated contact law.

The tensile set is cumulative. A base joins it the first time its σ_n is negative during the iteration, and stays
in it. With c_i removed, its σ_n turns positive, so a set rebuilt from the current σ_n alone flips between 7 and 0
bases and never converges. The final check requires that no new base is tensile at the converged state.

The thesis reference code grew the set across full re-solves until no new base went tensile. This version gives
the same values (Table 5.7): dry C0 2.637 / 2.430 and T0 2.597 / 2.429; wet C0 0.964 / 0.978; wet T0 has no
accepted FoS.

## Verification

`python verification.py --slow` (or the Verification page with the slow box ticked) gives **78 PASS, 0 FAIL,
12 INFO**. See `VERIFICATION.md`. Without the slow targets, the page lists the 9 fine-mesh rows (7 targets, 2
reference rows) as "skipped (slow), see VERIFICATION.md", so the totals still add up to 78 and 12. The targets are:

- the strength fixtures;
- the four main cases (3D and 2D, EJC and GHB);
- EJC at σ2 = σ3 against GHB (≤ 1e-11 relative; observed ≤ 2e-16);
- F_M, C0/T0 and the native-path iteration counts;
- the original-closure table (Table 5.2) and the 43-of-76 count;
- the Carranza-Torres records at 20×20 / 40 and at 320×320 / 5120;
- the finest meshes.

INFO rows are reference values: the published 2.01, Moscarda's reported values, and the "about" figures on the
Limitations page.

## Known limitations

- The surface is prescribed (spherical, symmetric); there is no critical-surface search.
- σ2 is imposed by the closure rather than solved from a 3D stress field, and ν is inherited, not measured.
- Hungr's simplified method neglects vertical intercolumn shear and does not satisfy all force and moment
  conditions. F_H and F_M differ by about 0.14–0.27 %.
- The rock mass is an equivalent isotropic medium, with no discontinuities.
- The water surface is hydrostatic, with no seepage solution.
- Tensile basal states are allowed. C0/T0 are only sensitivities, and the wet margin depends strongly on them
  (about 25 %).
- Coarse meshes underestimate FoS (10×10 is about 2.6 % low dry and 4.4 % low wet), and convergence is not an
  error bound.
- The bisection scan does not prove that the root is globally unique.
- Results are a research tool, not a design tool.

## Files

| File | Content |
|---|---|
| `app.py`, `views/` | Streamlit entry point and the seven pages |
| `ui.py` | Sidebar, input state, presets, running and caching |
| `engine.py` | Calculation engine (no Streamlit) |
| `plots.py` | Plotly figures (no Streamlit) |
| `presets.py` | Preset cases and their sources |
| `verification.py` | Targets and the verification runner |
| `tests/` | pytest suite: strength, LEM, targets, figures, headless pages |
| `VERIFICATION.md` | The verification table |
