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
pytest                            # about 3 minutes, includes the fine-mesh runs and a search
pytest -m "not slow"              # skips the fine-mesh runs and that search
python verification.py --slow     # the verification table as Markdown
```

The app is for local use. `.streamlit/config.toml` hides Streamlit's deploy button and turns off usage statistics.
Nothing in this folder deploys it to a public host.

## Pages

1. **Calculator.** Sidebar inputs: each has a number box and a slider that stay in sync, with units, valid range
   and reset-to-preset. Outputs:
   - FoS 3D and 2D, the 3D/2D ratio, iterations, active columns and slices, and the projection count;
   - the bases with tensile effective normal stress (σ_n < 0), with a note when there are any and a warning above
     a quarter of the bases (thesis §5.8);
   - a warning when the slip circle ends below the ground at its back extreme x_c − R (or, with the base-circle
     option, its front extreme x_c + R): the mass then has a vertical face on which no force acts, so the FoS is
     for that mass, not for a circular slip surface (`engine.vertical_faces`; see below);
   - PASS/FAIL for every final check, with the largest residuals;
   - CSV downloads of the column and slice results;
   - a one-page HTML report (`report.py`): the FoS, the cross-section, the notes, the inputs in use and the studies
     run on them;
   - studies: a critical-surface search (below), a mesh study (eq. 4.17), a ν sweep, EJC against GHB, a parameter
     sweep, a back-analysis and Monte Carlo (see "Additions").

   The sidebar also has **Save, open or share a case** (`casefile.py`): save the inputs as a JSON file, open a file
   or a pasted case code, or put the case in the page's link to share it. The web app reads and writes the same
   files and codes, so a case opens in either app with the same inputs. Its inputs for the additions (slope face,
   wall in plan, slip-surface shape, D zone, piezometric line, tension crack, k_h) are described under "Additions".

   The Mesh panel sets the 3D columns as N × N (the default), N × N_y, or a plan size Δx × Δy in metres (see
   "Column shape" below), and shows each column's plan size; the results give the grid and the range of column
   heights. Updates live up to 40 × 40 columns; for larger meshes press **Run**. With the original closure selected,
   the page shows the admissibility check instead of a FoS.
2. **Cross-section.** Central section y = 0: ground, slip surface, water table, Bishop slices (hover for σ_n,
   φ_i, c_i, τ_s and u) and the circle centre.
3. **3D view.** Ground, spherical slip surface and column prisms, coloured by σ_n, φ_i, c_i, τ_s, u, contribution
   to resistance τ_s A, or the tension flag. A toggle mirrors the half-model into the full mass, and a plan view
   highlights a column.
4. **Free-body diagram.** First every force on the sliding mass, on the 2D central section or the 3D full mass:
   - the weight W of each slice or column and the N′, U and S on its base (`engine.mass_forces`);
   - the force polygon, the totals with their moments about the rotation axis, and the vertical (eq. 4.12),
     horizontal and moment balance (F_H balances ΣS = ΣWκ; F_M and Bishop balance the moment).

   Then one column, chosen by clicking the plan view or by its index:
   - W, N′ = σ_n A, U = u A, S = τ_s A / F, the base normal, the lever arm ρ and the rotation axis;
   - the force values and the vertical-balance residual;
   - the matching 2D slice diagram.
5. **Strength explorer.** σ1–σ3 curves for EJC + Moscarda, EJC + original closure (with its ordering band shaded)
   and GHB, plus the Mohr circle, tangent and envelope at a chosen σ_n.
6. **Limitations.**
7. **Verification.** Runs every target in one go and shows target, computed value, difference and PASS/FAIL.

## Critical-surface search

QA/QC finding R1: the thesis's main sphere is not the critical surface. The Calculator's **Critical surface (search)**
study (`engine.critical_search`, also in the web app) finds the circle (2D, Bishop) and the sphere centred on y = 0
(3D, Hungr) with the lowest accepted FoS for the current inputs, mesh and grid option.

- **Surfaces tried.** Each surface is set by its front exit on the ground (the toe or the face; with the base-circle
  option also the pit floor, up to 2 H ahead of the toe), its back exit (the face or the crest, up to 2 H behind the
  crest edge) and the curvature of the arc between them. The arc meets the ground at both exits, so the mass never
  has a vertical back face. A circle whose centre lies below the ground at its back extreme x_c − R would have one:
  the engine accepts such a circle as a prescribed surface and gives the mass a back face that carries no force, as
  a tension crack would (for Taylor's 75° slope, a circle with a 2.8 m face gives 1.6245 against 1.6906 for the
  critical circle). The results summary warns about such a face, and the search does not try such circles.
- **Skipped surfaces.** A surface is skipped when the grid option cannot analyse it: the default grid refuses it, or
  its mass covers less than half of the default grid's span from x_c − R to the toe-level crossing, or (base-circle
  option) the ground is not above it over one interval. With the coverage rule, nearly flat circles resolved by a
  few slices are not tried: on Taylor's 60° slope a circle with R = 994 m, cut by one of 40 slices, gives 1.34 (2.55
  on 640 slices), against 1.94 for the critical circle. A result that is not ACCEPTED is withheld and does not count.
- **Method.** 2D: a grid of 4 × 9 × 4 = 144 starting circles (6 × 9 × 4 with the base-circle option), then
  Nelder–Mead from the two best grid circles that are not neighbours, each restarted from its minimum until it stops
  improving. 3D: Nelder–Mead from the 2D minima. The critical surface is the lowest accepted FoS among the surfaces
  tried. The result names any limit the surface lies on (a domain limit to widen, the steepest arc, the coverage rule,
  the default grid's limit) and re-solves it on a finer mesh (4 × N_2D slices, 2N × 2N columns; with N × N_y
  columns, 2N × 2N_y, each at most 320), with a warning when the FoS changes by more than the 0.5 % screening
  threshold of eq. 4.17. Every surface is solved with the same counts N and N_y, so the column size changes with the
  surface.
- **Mesh dependence.** The search minimises the FoS on the current mesh, so it favours surfaces whose FoS the mesh
  underestimates. On Taylor's 75° slope (φ = 0) with 40 slices it reports 1.6925, on a circle that gives 1.7503 on
  160 slices, and the mesh check warns; with 640 slices it reports 1.6889, against the exact 1.6906. On a 20 × 20
  mesh the 3D FoS changes in steps of a few tenths of a percent as columns enter or leave the mass, so searches from
  different starts can stop at FoS that differ by that much. Check the critical surface with the mesh check or a
  mesh study.
- **Results on the thesis slope** (40 slices and 20 × 20 columns; the main sphere in brackets; mesh check change):

  | Case | 2D critical FoS | 3D critical FoS |
  |---|---|---|
  | EJC dry | 2.1186 (2.4479), +0.3 % | 2.2600 (2.7060), +0.4 % |
  | EJC wet | 0.9929 (1.0981), +0.9 % | 1.0770 (1.2723), +1.8 % |
  | GHB dry | 2.0002 (2.3103), +0.3 % | 2.1359 (2.5537), +0.2 % |
  | GHB wet | 0.9369 (1.0359), +0.8 % | 1.0214 (1.2001), +1.1 % |

  The GHB dry circle's exits lie within 2 m of those of the Carranza-Torres & Hormazábal critical circle (FoS 2.01);
  the centres differ by 33 m, as circles through nearly the same exits differ only in curvature. The wet cases exceed
  the 0.5 % threshold on the finer mesh, so their critical FoS need a finer search mesh.
- **Time.** In the web app about 7 s (dry) to 17 s (wet) for 2D, and 17 to 37 s more for 3D. This engine is about
  20 times slower: a 2D search on the thesis slope takes about 1.5 minutes with GHB and 2.5 minutes with the EJC
  (dry), longer when wet.
- **Use this surface** makes the critical circle or sphere the prescribed surface and reruns the analysis.

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
- **Law (2024), Cases 1–4.** GHB on the overall 50° profile, with R exactly as given (through the toe), not the
  rounded reported R. They are cross-checked against the thesis reference code. Law's own and the Slide2 values
  are reference only.
- **D'Esterre (2024), 24 cases.** Three slopes (1:1.5, 1:2.5, 1:3.5; H = 2.1 m) × four c′/φ′ combinations ×
  r_u = 0 or 0.5. They use Mohr-Coulomb, r_u and the base-circle grid. The preset information shows
  `needs_floor_columns` (whether the circle exits ahead of the toe) and his Slide2, spreadsheet-3D and Slide3
  values, as reference only.
- **Baskoyak (real failure), 3 cases.** A circular failure through one bench in heavily jointed, dry schist at the
  Baskoyak barite pit, Turkey (Sonmez & Ulusay 1999, Case 2): γ = 22.2 kN/m³, σci = 5.2 MPa, m_i = 7. The section
  (H = 18.99 m, overall face 31.92°, planar here) and the failure circle (x_c = −2.19, z_c = 31.52, R = 32.40 m;
  it exits on the pit floor, so the base-circle grid is used) are traced from their Fig. 7. Three input sets: GHB
  with the SLIDE back-analysis of Li, Lo & Cassidy (2016) (GSI 16, D 0.68; SLIDE F = 1.007, reference only); GHB
  with site and guideline inputs (GSI 16.5, D 0.7; nothing fitted to the failure); and the EJC with ν = 0.45, a
  rock-mass value derived from published relations, not measured. `qa/qa_baskoyak.py` checks them.

The Law and D'Esterre cases are read from `data/presets_law_desterre.json` (already mirrored to engine
coordinates), the Baskoyak cases from `data/presets_baskoyak.json` (sources, traced geometry, derivation of ν).

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

The half-model (y ≥ 0) uses N × N nominal cells (N × N_y with the column-shape option below):

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

### Options for the D'Esterre (2024) presets (all off by default)

These were added for the D'Esterre soil-slope cases. With them off, every thesis result is unchanged.

- **Mohr-Coulomb strength** (`model="mc"`, `c_mc` kPa, `phi_mc` deg): τ_s = c′ + σ_n tan φ′.
  - There is no projection.
  - τ_s < 0 is allowed on the iteration path (as in the thesis workbook) but must be ≥ 0 in the final state.
  - The EJC residual checks do not apply.
- **Pore pressure by r_u** (`ru`): u = r_u γ h for every column or slice. It cannot be combined with the water
  table.
- **Base circles** (`base_circle=True`): the x grid spans the circle's exit-to-entry extent on the central section.
  - The exit and entry are found where the ground (pit floor z_t = 0 included) is above the circle, by sampling
    and bisection.
  - Floor columns are allowed: x > 0, z_t = 0, h = −z_b > 0. They appear only where the circle exits ahead of
    the toe.
  - The cell numbering, the y grid (0 to R, Δy = R/N) and all column equations are unchanged.

The default grid (x0 = x_c − √(R² − z_c²)) needs the circle to cross the toe level at or ahead of the toe. If it
crosses more than 1e-5 R behind the toe, the grid would cut off the toe region, so the default grid now refuses
such a circle and points to the base-circle option. The tolerance is the thesis workbook's; the Law and
Carranza-Torres circles are 0.3–2.6 mm behind the toe, within that tolerance. No existing preset is affected.

### Column shape (off by default)

The brief's 3D mesh has N × N cells (§1.4), so the column plan size is set by N alone: Δx = (x0 − x1)/N along x and
Δy = R/N across y. This option sets the two independently. With it off, every thesis result is unchanged.

- **Columns across y** (`n3y`, N_y): the y grid has N_y cells, Δy = R/N_y, y_j = (j − ½)Δy for j = 1 to N_y. The
  x grid, the active-column rule and every column equation (eq. 4.11, with A_h = Δx Δy) are unchanged. `n3y=None`
  (the default) is the N × N grid, and `n3y=N` gives it bit for bit (tested).
- **Plan size in metres** (the apps): a column length Δx and width Δy, both upper limits, become counts on the
  current circle, N = ⌈(x0 − x1)/Δx⌉ and N_y = ⌈R/Δy⌉ (`counts_for_size`); they must lie in 2 to 320, as N. When the
  circle changes, the counts follow. The third dimension, the column height h = z_t − z_b, is not an input: the
  ground and the slip surface set it. The apps show each column's Δx × Δy × h.
- **Studies.** The mesh study and the search's finer mesh scale N_y in the same ratio as N (`with_n3`,
  `columns_across`: N_y n/N rounded half up, at least 1), so the columns keep their shape.
- **Convergence.** Columns of very different shape converge to one result: on the thesis slope 160 × 40 columns
  (2.45 m by 18.1 m), 40 × 160 (9.81 m by 4.53 m) and 160 × 160 agree within 0.12 % dry and 0.26 % wet, inside the
  0.5 % threshold of eq. 4.17 (`test_rectangular_columns_converge`). At 20 × 20 columns (19.6 m by 36.3 m) the dry
  FoS is 2.7060; with 20 × 40 (19.6 m by 18.1 m) it is 2.7202, and with columns of at most 20 m by 20 m (20 × 37)
  2.7135.

### Additions: geometry, water, loading and studies (all off by default)

Nine additions extend the method beyond the thesis, added on request on 29 October 2026. Each is off by default, and
with all of them off the engine takes the thesis code path: the verification table is unchanged (89 PASS, 0 FAIL),
and 144 runs of the thesis, GHB, Law, D'Esterre and Baskoyak cases gave the same results, bit for bit, before and
after the change. Each formulation is standard limit-equilibrium practice or plain geometry, written out here so you
can check it; nothing was fitted to a result. `engine.py` and `webapp/engine.js` implement them line for line
(`webapp/README.md`), and `tests/test_additions.py` checks them against independent calculations where one exists.

**Studies (no new equations)**

- **EJC against GHB** (`compare_models`): the EJC with the Moscarda closure at the current ν, and GHB, on the same
  surface and mesh, in 3D and 2D, with the 3D/2D and EJC/GHB ratios.
- **Parameter sweep** (`parameter_sweep`): the FoS at evenly spaced values of one input, for the chosen models and
  dimensions; the other inputs stay as they are.
- **Back-analysis** (`back_analysis`): the value of one input at which the FoS equals a target, by bisection. It needs
  an accepted FoS at both ends of [lo, hi] and a sign change of F − target; it halves the bracket until it is
  narrower than 10⁻⁶ (hi − lo), then solves once more at the linear interpolation across the last bracket. If F steps
  (a column entering the mass), the result is the step. On the Baskoyak failure circle (GHB, 2D), F = 1 at
  D = 0.72087, where `qa/qa_baskoyak.py` found 0.721 from a cubic fit.
- **Monte Carlo** (`monte_carlo`): each chosen input is drawn from a normal distribution truncated to its valid
  range. The mean is the current value and the standard deviation is yours to enter: there is no default spread.
  Normal deviates come from the Box–Muller transform (Box and Muller 1958) on a seeded Mulberry32 stream, the same in
  both apps; draws outside the range are rejected. For each model and dimension the result gives:
  - the accepted count, mean and standard deviation (n − 1);
  - the reliability index β = (mean − 1)/SD;
  - the probability of failure P_f = P(F < 1) among the accepted results, with its standard error
    √(P_f(1 − P_f)/n); when no sample fails, the apps give the rule-of-three bound P_f < 3/n at 95 % (Hanley and
    Lippman-Hand 1983);
  - the 5, 50 and 95 % values.

  Withheld runs are counted and left out. The inputs are sampled independently: correlation (between GSI and σci,
  say) is not modelled.

**Slope geometry**

- **Slope profile** (`profile`): the face as points (x, z) from the toe (0, 0) to the crest edge (x_e, h_e), with x
  decreasing and z not decreasing. The ground is z = 0 on the pit floor ahead of the toe, the face by linear
  interpolation, and the crest behind the edge. The apps also build the points from benches: n benches of height h_b
  at the bench face angle β_b, separated by berms w_b wide (the bench terms of Read and Stacey 2009). H becomes h_e
  and β the overall angle from the toe. The critical-surface search follows the profile for its exits
  (`ground_point`). On the traced Baskoyak section the engine matches the independent Bishop code of
  `qa/qa_baskoyak.py` to 1e-12 (2000 slices).
- **Ground behind the crest** (`crest_slope`, α): z = h_e + (x_e − x) tan α behind the crest edge; α = 0 is flat.
- **Wall curved in plan** (`plan_radius`, R_p; 3D only): the toe line is a circle of radius |R_p| about a vertical
  axis at (R_p, 0). A point (x, y) lies on the section at x′, its signed distance from that line:
  - concave wall (R_p > 0, the pit inside the circle): x′ = R_p − r;
  - convex nose (R_p < 0): x′ = r − |R_p|;
  - with r = √((x − R_p)² + y²).

  Ground and water are taken at x′. On y = 0, x′ = x, so the 2D section is unchanged. A convex axis must lie behind
  the grid (|R_p| more than the distance from the toe to the back of the grid).

  On the thesis slope (80 × 80 columns, dry), a concave wall raises the FoS steadily: 2.7116 straight; 2.7146,
  2.7160, 2.7235 and 2.7568 at R_p = 6400, 3200, 1600 and 800 m. A convex nose does not change it steadily on a fixed
  sphere: 2.7100, 2.7086, 2.7105 and 2.7183 at R_p = −6400, −3000, −1500 and −800 m. Before reading a trend from
  convex walls, search for the critical sphere at each radius.

**Slip-surface shape** (`elong`, λ; `mid_length`, L; 3D only)

z_b = z_c − √q with q = R² − (x − x_c)² − (y_c/λ)² and y_c = max(|y| − L/2, 0): ellipsoidal end caps of half-width
λR around a cylindrical middle part of length L. The y grid spans the half-width Y = λR + L/2 (Δy = Y/N_y); λ = 1,
L = 0 is the sphere. At each y the surface is a circle about the rotation axis through (x_c, z_c), of radius
ρ = √(R² − (y_c/λ)²). The base normal, n ∝ (x_c − x, −y_c/λ², z_c − z_b), has its x–z part pointing at the axis and
its y part parallel to it, so N′ and U have no moment about the axis. Hungr's equations therefore keep their form,
with ρ per column. The central section stays the circle (x_c, z_c, R), so the 2D analysis is unchanged.

On the thesis slope (40 × 40 columns; 2D with 40 slices), 3D/2D is 1.110 for the sphere and 1.104, 1.102 and 1.100
at λ = 1.5, 2 and 3: wider caps alone keep the end effect. A cylindrical middle removes it: 3D/2D is 1.038, 1.013 and
1.003 at L = 500, 2000 and 10,000 m. The 3D FoS tends to the 2D FoS, dry and wet, within the 0.5 % of eq. 4.17 at
L = 100R (`test_a_long_slide_tends_to_the_2d_result`).

**Rock mass**

- **Disturbance zone** (`d_zone`, T): bases closer than T to the face take m_b, s and a with D; the others take
  D = 0. The distance is the shortest one from the base point (x, z_b) on its section to the face line
  (`face_distance`). Each zone has its own projection limit (Section 4.6). The option is not used with Mohr-Coulomb.
  - D is the blast-damage factor of Hoek et al. (2002). Hoek (2012, a RocNews technical note) recommends applying D
    only within a damaged zone of thickness T behind the face.
  - That note could not be opened from the environment these additions were written in. It is cited from memory,
    and its suggested thicknesses are not used: T has no default, so enter it from your site.
  - T = 0 gives the D = 0 result, and T beyond every base the D-everywhere result, both exactly.

**Water**

- **Piezometric line** (`piezo`): u = γ_w max(0, z_w − z_b), with z_w interpolated along the points (x decreasing)
  and constant beyond the end points. On a curved wall it is swept like the ground. It replaces H_w and β_w. The
  line (0, 0), (−H_w/tan β_w, H_w) reproduces the water table exactly.
- **Tension crack** (`crack_depth`, z_t; `crack_water`, f_w):
  - The crack is vertical, at x_tc, the first point from the back of the central section where the slip circle lies
    z_t below the ground. The mass ends there: columns and slices with x < x_tc are removed.
  - Water stands f_w z_t above the crack bottom on the central section, at the same level across y.
  - Each row of columns (in 2D the section, per metre) carries the hydrostatic force P = ½ γ_w h_w² (times Δy),
    horizontal towards the pit, at h_w/3 above the crack bottom, with lever arm e = z_c − z_P about the axis.
  - This is the usual treatment of a water-filled tension crack in limit equilibrium (Hoek and Bray 1981; Duncan and
    Wright 2005).

**Loading**

- **Seismic load** (`kh`): the pseudo-static method (Terzaghi 1950; Kramer 1996). A horizontal force k_h W acts at
  each column's or slice's centroid, at height (z_b + z_t)/2, towards the pit.

**The global equation with these forces.** A horizontal force Q at height z_Q has the moment Q e about the rotation
axis, e = z_c − z_Q, in the same sense as the weight's W (x_c − x):

- F_M (eq. 4.15) adds the moments: F = Σ τ A ρ / (Σ W (x_c − x) + Σ Q e);
- F_H (eq. 4.14) adds each moment divided by the lever arm ρ of its column or row, as W κ = W (x_c − x)/ρ:
  F = Σ τ A / (Σ W κ + Σ Q e/ρ);
- Bishop's 2D equation does the same per metre.

Local vertical balance (eq. 4.12) is unchanged, since the forces are horizontal. In 2D the engine matches an
independent Bishop calculation with these moment terms to 1e-10. The whole-mass free-body diagram shows E = k_h W and
P, and its moment balance includes them.

References for the additions:

- Box, G. E. P. and Muller, M. E. (1958). A note on the generation of random normal deviates. *Annals of Mathematical
  Statistics* 29(2), 610–611.
- Duncan, J. M. and Wright, S. G. (2005). *Soil Strength and Slope Stability*. Wiley.
- Hanley, J. A. and Lippman-Hand, A. (1983). If nothing goes wrong, is everything all right? Interpreting zero
  numerators. *JAMA* 249(13), 1743–1745.
- Hoek, E. (2012). Blast damage factor D. Technical note for RocNews, February 2012.
- Hoek, E., Carranza-Torres, C. and Corkum, B. (2002). Hoek-Brown failure criterion – 2002 edition. *Proc. NARMS-TAC
  2002*, Toronto, 267–273.
- Hoek, E. and Bray, J. W. (1981). *Rock Slope Engineering*, 3rd edn. Institution of Mining and Metallurgy, London.
- Kramer, S. L. (1996). *Geotechnical Earthquake Engineering*. Prentice Hall.
- Read, J. and Stacey, P. (eds) (2009). *Guidelines for Open Pit Slope Design*. CSIRO Publishing.
- Terzaghi, K. (1950). Mechanism of landslides. In *Application of Geology to Engineering Practice* (Berkey volume),
  Geological Society of America, 83–123.

## Verification

`python verification.py --slow` (or the Verification page with the slow box ticked) gives **89 PASS, 0 FAIL,
40 INFO**. See `VERIFICATION.md`. Without the slow targets, the page lists the 9 fine-mesh rows (7 targets, 2
reference rows) as "skipped (slow), see VERIFICATION.md", so the totals still add up. The PASS/FAIL rows are:

- the strength fixtures;
- the four main cases (3D and 2D, EJC and GHB);
- EJC at σ2 = σ3 against GHB (≤ 1e-11 relative; observed ≤ 2e-16);
- F_M, C0/T0 and the native-path iteration counts;
- the original-closure table (Table 5.2, ν = 0.05 to 0.49) and the 43-of-76 count;
- the Carranza-Torres records at 20×20 / 40 and at 320×320 / 5120;
- the finest meshes;
- the Law (2024) presets against the thesis reference code `lem_ref.py` (GHB 2D, 40 and 640 slices, 4 decimals);
- the Mohr-Coulomb option against the `slope3d` package's record (c′ = 100 kPa, φ′ = 0: 0.2553545 / 0.1902638).

The first seven groups are the 78 targets of the brief and the thesis records.

INFO rows are reference values, never targets:

- the published 2.01;
- the reported values of Moscarda (2025), Law (2024) and D'Esterre (2024);
- the "about" figures on the Limitations page.

## Known limitations

- Each analysis uses a prescribed sphere, symmetric about y = 0. The critical-surface search tries only circles
  and spheres centred on y = 0; it is a local search on the current mesh, not a proof of the global minimum.
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
- The additions are standard formulations implemented here and checked against independent calculations where one
  exists, but they are not part of the thesis, and no published 3D case checks them.
  - Curved walls and slip-surface shapes are 3D only; the 2D section is y = 0.
  - The tension crack is placed on the central section and runs straight across y.
  - The seismic load is pseudo-static, with one constant k_h.
  - Monte Carlo samples its inputs independently, without correlation.
- Results are a research tool, not a design tool.

## Files

| File | Content |
|---|---|
| `app.py`, `views/` | Streamlit entry point and the seven pages |
| `ui.py` | Sidebar, input state, presets, cases, running and caching |
| `study_panels.py` | Calculator studies: EJC against GHB, parameter sweep, back-analysis, Monte Carlo |
| `casefile.py` | Case files and the page link, in the web app's format (no Streamlit) |
| `report.py` | The one-page HTML report, with the cross-section as SVG (no Streamlit) |
| `engine.py` | Calculation engine (no Streamlit) |
| `plots.py` | Plotly figures (no Streamlit) |
| `presets.py` | Preset cases and their sources |
| `data/presets_law_desterre.json` | Law (2024) and D'Esterre (2024) cases in engine coordinates (source of those presets) |
| `data/presets_baskoyak.json` | Baskoyak case: sources, site data, traced section and failure circle, derivation of ν (source of those presets) |
| `verification.py` | Targets and the verification runner |
| `tests/` | pytest suite: strength, LEM, targets, figures, headless pages, the additions (`test_additions.py`, `test_app_additions.py`) |
| `VERIFICATION.md` | The verification table |
