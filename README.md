# ejc-slope-stability

GENG5512 thesis code, brought over from `mnersubash-arch/slope3D` (branch `claude/awesome-ride-ukq1wd`) with its full
history:

- `ejc_app/`: the interactive Streamlit app for the submitted thesis method (EJC with the Moscarda closure, Route 2,
  Hungr 3D / Bishop 2D, with the Law and D'Esterre presets and a critical-surface search). Run it with
  `cd ejc_app && streamlit run app.py`; see `ejc_app/README.md` and `ejc_app/VERIFICATION.md`.
- `webapp/`: the same app as a web page, with the engine ported to JavaScript and checked against
  `ejc_app/engine.py`; published privately on claude.ai. See `webapp/README.md`.
- `qa/`: the QA/QC of both engines against the thesis, exact solutions, randomised tests and published case
  studies (27 October 2026): the scripts and the record. See `qa/README.md`.
- `slope3d/` and the other files described below: the earlier port of the thesis workbook (VBA).

Research tool for the prescribed-surface comparison in Mahat (2026) thesis; not for design.

---

# slope3d

> **See also `ejc_app/`:** the interactive Streamlit app for the submitted thesis method (EJC with the Moscarda
> closure, Route 2, Hungr 3D / Bishop 2D). It is a separate engine; run it with `cd ejc_app && streamlit run app.py`.

Standalone 3D limit-equilibrium rock-slope stability with a 3D Hoek-Brown (EJC eq. 4.21) strength model: a port of the
verified GENG5512 thesis workbook (`reference/`) with a command line, a Streamlit interface, plots and a run report.
**The science is unchanged**: the VBA is the specification (`docs/SPEC.md` lists every rule with its VBA line, and the
few places where the VBA, the Python mirror and the brief differed).

* 3D Hungr (1987) simplified columns (n x n half-model) and 2D Bishop slices on a prescribed spherical slip surface.
* A factor of safety is shown **only** for an ACCEPTED dimension; a NOT ACCEPTED run shows `FoS withheld` and the reason.
* Closure 3 is always labelled *Conventional HB comparator (3D geometry, 2D strength)*, never a 3D-EJC result.

## Install on Windows

1. Install Python 3.11 or newer from python.org (tick "Add python.exe to PATH").
2. In PowerShell, in this folder:

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[gui]"          # engine + CLI + GUI + plots.  Add ",dev" instead of ",gui" to run the tests
```

## Run

```powershell
slope3d gui                                       # opens the interface in your browser (http://localhost:8501)
slope3d presets                                   # list the built-in benchmark cases
slope3d run cases\chile_dry_conventional_hb.yaml --mesh3d 20 --mesh2d 40 --out out --format csv,xlsx,json,report
slope3d run cases\chile_thesis_default.yaml       # thesis closure: both NOT ACCEPTED, FoS withheld (exit status 1)
slope3d validate                                  # every regression check, PASS/FAIL
slope3d validate --slow                           # + 95-record mesh study and nu thresholds (about 30 s)
slope3d mesh-study --out mesh.csv --compare reference\mesh_study\mesh_study_records.csv
```

`run` exit status: 0 both dimensions ACCEPTED, 1 at least one not accepted (a valid outcome), 2 input error.
`--format` accepts `csv`, `xlsx`, `json`, `pdf`, `html` or `report` (= pdf + html).

### The interface

* **Sidebar:** preset buttons (the six benchmark cases and the thesis closure), then the inputs grouped as on the workbook's
  Interface sheet: geometry, slope, water, rock mass, strength model, closure, tolerances, mesh.
* **Status panel:** a large ACCEPTED / NOT ACCEPTED / DIAGNOSTIC ONLY label for 3D and for 2D, with the reason.
* **Tabs:** Summary (F 3D, F 2D, ratio, iterations, slip-surface entry/exit, number and minimum of elements with
  sigma_n' < 0), Element table (same 23 fields as the VBA `VBA Columns` sheet), Plots (six, interactive) and Exports
  (CSV, XLSX, JSON, case file, PDF/HTML report).

### Case files

YAML or JSON, flat or grouped like the Interface sheet; see `cases\`. Unknown keys are an error.
The slope can be given as `face_angle_deg` or as `slope` (= tan of the angle).

### Acceptance checks off = diagnostic mode

With `checks: 0` an HB run uses the first converged eq. 4.21 root of every element with no admissibility test and checks
only local equilibrium and convergence of F. The result is always **DIAGNOSTIC ONLY**; its value is shown only as "Diagnostic F (unchecked roots, not an accepted result)" (grey card, never the ACCEPTED style), and it is not a factor of safety. NOT ACCEPTED runs still show "FoS withheld". This reproduces
the workbook's worksheet engine (`Interface!P19 = 0`); it is not part of the VBA engine. The default is `checks: 1`.

## Tests

```powershell
pip install -e ".[dev]"
pytest                    # suites A (22 strength fixtures), B (engine), vector equivalence, report/GUI  (about 90 s)
pytest -m slow            # suite C: 95 mesh-study records and the exact-EJC thresholds (about 30 s)
python tools\compare_with_workbook.py     # suite D: every recorded value in the .xlsm QA sheets (about 20 s)
```

## Single .exe (optional, untested on Windows by the author of this port)

```powershell
pip install -e ".[exe]"
pyinstaller packaging\slope3d.spec --noconfirm
dist\slope3d.exe                  # double-click: GUI;   dist\slope3d.exe run case.yaml  /  validate  also work
```

## Limits

* Prescribed spherical slip surface; no search for the critical surface.
* Symmetric half-model (y_c = 0).
* Hungr's vertical-equilibrium simplification.
* Only the positive-centre, toe-cutoff geometry is supported; other geometry is refused (`UNSUPPORTED GEOMETRY`).
* Closures 1 and 2 scan for roots element by element and are slow on very fine 3D meshes (the default 20x20 takes under a second).
* Tolerances can be tightened but not loosened (as in the VBA), and nu must lie in (0, 0.5).
