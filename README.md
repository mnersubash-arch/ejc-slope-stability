# slope3d

Standalone 3D limit-equilibrium rock-slope stability with a 3D Hoek-Brown (EJC eq. 4.21) strength model:
a port of the verified GENG5512 thesis workbook (`reference/`). The science is unchanged; the VBA is the specification
(`docs/SPEC.md`).

**Status: engine, strength library, CLI and tests are done. The Streamlit GUI, plots/report and the PyInstaller packaging
are waiting for your review of `docs/SPEC.md`.**

## Install (Windows, PowerShell)

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
```

## Use

```powershell
slope3d presets                                   # list the built-in benchmark cases
slope3d run cases\chile_dry_conventional_hb.yaml --mesh3d 20 --mesh2d 40 --out out
slope3d run cases\chile_thesis_default.yaml       # thesis closure: both NOT ACCEPTED, FoS withheld
slope3d validate                                  # every regression check, PASS/FAIL
slope3d validate --slow                           # + 95-record mesh study and nu thresholds (~30 s)
slope3d mesh-study --out mesh.csv --compare reference\mesh_study\mesh_study_records.csv
pytest                                            # fast suites;  pytest -m slow  for suite C
python tools\compare_with_workbook.py             # cross-check against the .xlsm's recorded values
```

A factor of safety is shown only for an ACCEPTED dimension; otherwise the output says `FoS withheld`.
