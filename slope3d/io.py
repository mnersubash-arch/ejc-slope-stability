"""Load and save cases (JSON/YAML) and export results (CSV/XLSX).

Case files may be flat (every :class:`~slope3d.models.Config` field at the top level) or grouped like the Interface
sheet (``geometry``, ``slope``, ``water``, ``rock_mass``, ``strength``, ``closure``, ``tolerances``, ``mesh``);
groups are flattened on load. ``mesh3d``/``mesh2d`` may sit at the top level or under ``mesh``.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import yaml

from .geometry import slip_surface_ends
from .models import CaseResult, Config, Elements, ModelRun

GROUPS = ("geometry", "slope", "water", "rock_mass", "strength_model", "closure_inputs", "tolerances", "mesh",
          "strength", "closure", "run")
# table columns exactly as the VBA 'VBA Columns' / 'VBA Slices' sheets (WriteElements, line 344)
ELEMENT_COLUMNS = ["x", "y", "z_base", "height", "A_base", "W", "u", "g", "kappa", "rho", "trial normal",
                   "actual final normal", "phi deg", "c kPa", "tau kPa", "sigma1 MPa", "sigma2 MPa", "sigma3 MPa",
                   "angle rad", "EJC MPa", "closure MPa", "eq/W", "status"]
ELEMENT_UNITS = ("x/y/z m, area m2, forces kN, normal/c/shear kPa, principal stresses MPa "
                 "(columns: 'trial normal' = last iterate, 'actual final normal' = rebuilt at the reported F)")
NOT_EVALUATED = "NOT EVALUATED"


# --------------------------------------------------------------------------------------------------
# case files
# --------------------------------------------------------------------------------------------------
def _flatten(data: dict[str, Any]) -> dict[str, Any]:
    flat: dict[str, Any] = {}
    for key, val in data.items():
        if isinstance(val, dict) and key in GROUPS:
            for k2, v2 in val.items():
                if k2 in flat:
                    raise ValueError(f"input '{k2}' given twice")
                flat[k2] = v2
        else:
            if key in flat:
                raise ValueError(f"input '{key}' given twice")
            flat[key] = val
    return flat


def parse_case(data: dict[str, Any]) -> tuple[Config, int, int]:
    """Build (Config, mesh3d, mesh2d) from a parsed case document. Unknown keys are an error, not ignored."""
    flat = _flatten(dict(data))
    flat.pop("name", None)
    flat.pop("description", None)
    mesh3d = int(flat.pop("mesh3d", 20))
    mesh2d = int(flat.pop("mesh2d", 40))
    return Config.from_dict(flat), mesh3d, mesh2d


def load_case(path: str | Path) -> tuple[Config, int, int]:
    """Load a ``.yaml``/``.yml``/``.json`` case. Returns (Config, mesh3d, mesh2d)."""
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    data = json.loads(text) if p.suffix.lower() == ".json" else yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError(f"{p}: a case file must contain a mapping of inputs")
    return parse_case(data)


def save_case(path: str | Path, cfg: Config, mesh3d: int = 20, mesh2d: int = 40, name: str = "") -> None:
    """Save a case as grouped JSON or YAML (by extension); round-trips through :func:`load_case`."""
    d = cfg.to_dict()
    pick = lambda *keys: {k: d[k] for k in keys}   # noqa: E731
    doc = {
        "name": name or "slope3d case",
        "geometry": pick("xc", "zc", "radius", "height", "gamma"),
        "slope": {},
        "water": pick("water_angle", "water_height", "u_extra"),
        "rock_mass": pick("sci", "mi", "gsi", "d", "a_override"),
        "strength_model": pick("strength", "c_mc", "phi_mc"),
        "closure_inputs": pick("closure", "nu"),
        "tolerances": pick("tol_a", "tol_e", "tol_f", "checks", "max_iter"),
        "mesh": {"mesh3d": mesh3d, "mesh2d": mesh2d},
    }
    doc["slope"] = {"slope": cfg.slope}      # exact tan value (face_angle_deg would be a second source of the same input)
    p = Path(path)
    if p.suffix.lower() == ".json":
        p.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    else:
        p.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")


# --------------------------------------------------------------------------------------------------
# result tables
# --------------------------------------------------------------------------------------------------
def elements_table(run: ModelRun, cfg: Config) -> pd.DataFrame:
    """The element table with the same 23 fields as the VBA ``VBA Columns``/``VBA Slices`` sheets.

    Mirrors ``WriteElements``: final normal and shear are 'NOT EVALUATED' unless the final state was evaluated,
    and the principal-stress/residual columns exist only for Hoek-Brown runs.
    """
    res, els = run.result, run.elements
    if els is None or res.active < 1:
        return pd.DataFrame(columns=ELEMENT_COLUMNS)
    n = els.n
    na = np.full(n, NOT_EVALUATED, dtype=object)
    cols: dict[str, Any] = {
        "x": els.x, "y": els.y, "z_base": els.zb, "height": els.h, "A_base": els.area, "W": els.weight,
        "u": els.pressure, "g": els.g, "kappa": els.kap, "rho": els.rho, "trial normal": els.sn,
        "actual final normal": els.final_n if res.evaluated else na,
        "phi deg": els.phi, "c kPa": els.cohesion,
        "tau kPa": els.tau if res.evaluated else na,
    }
    blank = np.full(n, None, dtype=object)
    hb = res.evaluated and cfg.strength == 1
    cols["sigma1 MPa"] = els.s1 if hb else blank
    cols["sigma2 MPa"] = els.s2 if hb else blank
    cols["sigma3 MPa"] = els.s3 if hb else blank
    cols["angle rad"] = els.angular if hb else blank
    cols["EJC MPa"] = els.ejc if hb else blank
    cols["closure MPa"] = els.closure_res if hb else blank
    cols["eq/W"] = els.equilibrium if res.evaluated else blank
    cols["status"] = els.label
    return pd.DataFrame(cols, columns=ELEMENT_COLUMNS)


def summary_table(case: CaseResult) -> pd.DataFrame:
    """Key/value summary. A withheld FoS is written as the word WITHHELD, never as a number."""
    rows: list[tuple[str, Any]] = [("result identity", case.label)]
    for tag, run, mesh in (("3D", case.run3d, case.mesh3d), ("2D", case.run2d, case.mesh2d)):
        r = run.result
        rows += [(f"{tag} mesh", mesh), (f"{tag} status", r.status.value), (f"{tag} FoS", r.fos if r.fos is not None else "WITHHELD"),
                 (f"{tag} active elements", r.active), (f"{tag} iterations", r.iterations), (f"{tag} reason", r.reason),
                 (f"{tag} min final normal kPa", r.min_normal if r.evaluated else NOT_EVALUATED),
                 (f"{tag} tension elements", r.tension if r.evaluated else NOT_EVALUATED)]
    ends = slip_surface_ends(case.config)
    if ends.get("supported"):
        rows += [("slip surface entry (crest side) x, z [m]", f"{ends['entry_x']:.2f}, {ends['entry_z']:.2f}"),
                 ("slip surface exit x, z above toe [m]", f"{ends['exit_x']:.2f}, {ends['exit_z']:.2f}"
                  + ("" if ends["exit_on_face"] else " (reaches the toe plane)")),
                 ("exit behind toe / below crest edge [m]", f"{ends['exit_behind_toe']:.2f} / {ends['exit_distance_from_crest']:.2f}")]
    rows.append(("3D/2D ratio", case.ratio if case.ratio is not None else "n.a."))
    m = case.run3d.result.moment_diagnostic
    rows.append(("spherical moment ratio (diagnostic only)", m if m is not None else "NOT EVALUATED"))
    return pd.DataFrame(rows, columns=["item", "value"])


def inputs_table(cfg: Config, mesh3d: int, mesh2d: int) -> pd.DataFrame:
    d = cfg.to_dict()
    d.update(mesh3d=mesh3d, mesh2d=mesh2d, mb=cfg.mb, s=cfg.s, a_used=cfg.a)
    return pd.DataFrame(list(d.items()), columns=["input", "value"])


def iterations_table(case: CaseResult) -> pd.DataFrame:
    """Iteration log (VBA 'VBA Iterations'). The trial F of a run that is not ACCEPTED is written as WITHHELD."""
    rows = []
    for dim, run in ((3, case.run3d), (2, case.run2d)):
        for rec in run.result.history:
            rows.append((dim, rec.iteration, rec.f if run.result.accepted else "WITHHELD", rec.df_rel,
                         rec.max_move_kpa, rec.projections))
    return pd.DataFrame(rows, columns=["dimension", "iteration", "F", "dF/F", "max normal move kPa", "path projections"])


def export_csv(case: CaseResult, outdir: str | Path, stem: str = "slope3d") -> list[Path]:
    """Write ``<stem>_summary.csv``, ``_columns3d.csv``, ``_slices2d.csv`` and ``_inputs.csv``. Returns the paths."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = {
        f"{stem}_summary.csv": summary_table(case),
        f"{stem}_inputs.csv": inputs_table(case.config, case.mesh3d, case.mesh2d),
        f"{stem}_columns3d.csv": elements_table(case.run3d, case.config),
        f"{stem}_slices2d.csv": elements_table(case.run2d, case.config),
    }
    paths = []
    for name, df in files.items():
        p = out / name
        df.to_csv(p, index=False)
        paths.append(p)
    return paths


def export_xlsx(case: CaseResult, path: str | Path) -> Path:
    """One workbook: Summary, Inputs, Columns 3D, Slices 2D (same fields as the VBA sheets), Iterations."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(p, engine="openpyxl") as xw:
        summary_table(case).to_excel(xw, sheet_name="Summary", index=False)
        inputs_table(case.config, case.mesh3d, case.mesh2d).to_excel(xw, sheet_name="Inputs", index=False)
        elements_table(case.run3d, case.config).to_excel(xw, sheet_name="Columns 3D", index=False, startrow=3)
        elements_table(case.run2d, case.config).to_excel(xw, sheet_name="Slices 2D", index=False, startrow=3)
        for sheet, run in (("Columns 3D", case.run3d), ("Slices 2D", case.run2d)):
            ws = xw.sheets[sheet]
            ws["A1"] = ELEMENT_UNITS
            ws["A2"] = f"{run.result.status.value}: {run.result.reason}"
        it = iterations_table(case)
        it.to_excel(xw, sheet_name="Iterations", index=False)
    return p


def result_to_dict(case: CaseResult) -> dict[str, Any]:
    """JSON-safe summary of a run (F is ``null`` unless ACCEPTED)."""
    def one(run: ModelRun, mesh: int) -> dict[str, Any]:
        r = run.result
        return {"mesh": mesh, "status": r.status.value, "fos": r.fos, "reason": r.reason, "active": r.active,
                "iterations": r.iterations, "admissible": r.admissible,
                "min_normal_kpa": r.min_normal if r.evaluated else None, "tension_elements": r.tension if r.evaluated else None}
    return {"identity": case.label, "inputs": case.config.to_dict(), "3d": one(case.run3d, case.mesh3d),
            "2d": one(case.run2d, case.mesh2d), "ratio_3d_2d": case.ratio}
