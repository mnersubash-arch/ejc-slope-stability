# EJC 3D/2D rock-slope app (work in progress)

Interactive limit-equilibrium app for the GENG5512 thesis: Hungr's simplified 3D method of columns and Bishop's 2D
method on a prescribed spherical slip surface, with the Extended Jiang Criterion (EJC) and the Moscarda closure.

**Status:** step 1-2 of the brief (engine, tests, verification table). The Streamlit pages are not built yet.

This folder is separate from the `slope3d` package at the repository root, which is the port of the thesis
workbook (VBA). The two do not share code. The GHB results of this engine reproduce that package's closure-3
benchmarks (for example 2.5536684 / 2.3102550 dry, a = 0.5).

## Files

| File | Content |
|---|---|
| `engine.py` | Calculation engine (no Streamlit): rock-mass parameters, EJC states on the closures, Route 2 strength solver, mesh, Hungr/Bishop iteration, final checks, original-closure admissibility minimum |
| `verification.py` | The brief's §3 targets and a runner that prints target, computed, difference and PASS/FAIL/INFO |
| `tests/` | pytest suite (strength, LEM, targets) |
| `VERIFICATION.md` | The verification tables |

## Run

```bash
pip install -r requirements.txt
pytest                                   # about 40 s, includes the finest-mesh runs
python verification.py --slow            # verification table (Markdown)
python verification.py --slow --safeguard            # with the trial-stress projection
python verification.py --slow --safeguard --workbook # workbook iteration path
```

## Engine options that are not in the brief

The defaults follow brief §1.4 as written. Two options are needed to reproduce some recorded values. Both are
off by default and are listed in `VERIFICATION.md`:

* `safeguard="project"`: a trial sigma_n below the lowest admissible sigma_n of the closure is moved to just above
  it during the iteration. The run is accepted only if no base is projected in the final iteration, so F is the
  same fixed point.
* `iteration="workbook"`: the workbook/VBA path (start at F = 1, Hungr explicit form). It gives the same F to
  about 1e-11 and reproduces the recorded iteration counts.
