"""Mesh-convergence study and the exact-EJC admissibility threshold.

``run_mesh_study`` repeats the study of ``reference/mesh_study/hungr_ejc_mesh_study.py`` (conventional-HB comparator,
closure 3) with the slope3d engine. ``exact_ejc_threshold`` ports ``reference/mesh_study/nu_exact.py``.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from .engine import solve_model
from .models import Config

SURFACES = {"AppendixE": (355.742311546406, 644.419294864537, 725.2669975187362),
            "Published": (207.28, 586.53, 622.08)}
CASES = [("Dry, a = 0.500000", "AppendixE", 0.5, False), ("Wet, a = 0.500000", "AppendixE", 0.5, True),
         ("Dry, a = 0.509923", "AppendixE", 0.509923, False), ("Wet, a = 0.509923", "AppendixE", 0.509923, True),
         ("Dry, a = 0.500000, published circle", "Published", 0.5, False)]
MESH = {3: (10, 15, 20, 30, 40, 60, 80, 120, 160, 240, 320), 2: (20, 40, 80, 160, 320, 640, 1280, 2560)}


def case_config(surface: str, a: float, wet: bool) -> Config:
    xc, zc, r = SURFACES[surface]
    return Config(xc=xc, zc=zc, radius=r, closure=3, a_override=a,
                  water_angle=50.0 if wet else 0.0, water_height=360.0 if wet else 0.0)


def run_mesh_study(meshes3=MESH[3], meshes2=MESH[2], cases=CASES,
                   progress: Optional[Callable[[dict], None]] = None) -> list[dict]:
    """One row per case/dimension/mesh. ``F`` is ``None`` unless the run is ACCEPTED."""
    rows: list[dict] = []
    for label, surface, a, wet in cases:
        cfg = case_config(surface, a, wet)
        for dim, meshes in ((3, meshes3), (2, meshes2)):
            prev = None
            for n in meshes:
                res = solve_model(cfg, dim, n).result
                f = res.fos
                change = None if (prev is None or f is None) else abs(f - prev) / f * 100
                prev = f if f is not None else prev
                row = dict(case=label, surface=surface, a=a, wet=wet, dim=f"{dim}D",
                           mesh=f"{n}x{n}" if dim == 3 else f"{n} slices", elements=res.active,
                           volume_m3=res.volume, change_from_previous_pct=change, F=f, iters=res.iterations,
                           proj_path=res.projections, min_sn=res.min_normal if res.evaluated else None,
                           n_tension=res.tension if res.evaluated else None, accepted=res.accepted,
                           status=res.status.value)
                rows.append(row)
                if progress:
                    progress(row)
    return rows


def compare_with_records(rows: list[dict], records_csv: str | Path, rel_tol: float = 1e-9) -> list[str]:
    """Return a list of discrepancies against ``mesh_study_records.csv`` (empty list = all reproduced)."""
    ref = {(r["case"], r["dim"], r["mesh"]): r for r in csv.DictReader(open(records_csv, newline=""))}
    problems = []
    for row in rows:
        key = (row["case"], row["dim"], row["mesh"])
        r = ref.get(key)
        if r is None:
            continue
        f_ref = float(r["F"])
        if row["F"] is None or abs(row["F"] - f_ref) > rel_tol * abs(f_ref):
            problems.append(f"{key}: F {row['F']} vs record {f_ref}")
        if bool(row["accepted"]) != (r["accepted"] == "True"):
            problems.append(f"{key}: accepted {row['accepted']} vs record {r['accepted']}")
        if row["elements"] != int(r["elements"]):
            problems.append(f"{key}: elements {row['elements']} vs record {r['elements']}")
        if abs(row["volume_m3"] - float(r["volume_m3"])) > rel_tol * float(r["volume_m3"]):
            problems.append(f"{key}: volume {row['volume_m3']} vs record {r['volume_m3']}")
    return problems


def write_rows(rows: list[dict], path: str | Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


# --------------------------------------------------------------------------------------------------
# exact-EJC admissibility threshold (nu_exact.py)
# --------------------------------------------------------------------------------------------------
def _state(s3, nu, a, sci, mb, s):
    """Exact thesis-closure EJC tangent state at sigma3' (MPa): (sigma_n', phi deg) or ``None`` if not admissible."""
    q = mb * s3 / sci + s
    k = 1 / (2 * nu) - 1
    rad = sci ** 2 * q ** (2 * a) - k * k * s3 * s3
    if rad < 0:
        return None
    s1 = s3 / (2 * nu) + math.sqrt(rad) / math.sqrt(3)
    s2 = s3 / nu - s1
    if not (s1 >= s2 >= s3):
        return None
    m = (3 * s1 / nu - 2 * (1 - 1 / nu + 1 / nu ** 2) * s3 + 2 * a * mb * sci * q ** (2 * a - 1)) / (6 * s1 - 3 * s3 / nu)
    if m < 1:
        return None
    p = (m - 1) / (m + 1)
    return (s1 + s3) / 2 - (s1 - s3) / 2 * p, math.degrees(math.asin(p))


def exact_ejc_threshold(nu: float, a: float, sci: float = 77.7, mi: float = 10.0, gsi: float = 42.0, d: float = 0.0,
                        lo: float = 1e-6, hi: float = 1e6, npts: int = 200000):
    """Smallest sigma_n' (MPa) at which the exact thesis-closure EJC has an admissible state, searched over sigma3' in
    [lo, hi] MPa (log-spaced) and refined by golden section. Returns (sigma_n', phi deg, sigma3') or ``None``.
    """
    mb = mi * math.exp((gsi - 100) / (28 - 14 * d))
    s = math.exp((gsi - 100) / (9 - 3 * d))
    best = None
    l0, l1 = math.log10(lo), math.log10(hi)
    for i in range(npts + 1):
        s3 = 10 ** (l0 + (l1 - l0) * i / npts)
        st = _state(s3, nu, a, sci, mb, s)
        if st and (best is None or st[0] < best[0]):
            best = (st[0], st[1], s3)
    if best is None:
        return None
    a_, b_ = best[2] / 1.01, best[2] * 1.01
    g = (math.sqrt(5) - 1) / 2
    f = lambda x: (_state(x, nu, a, sci, mb, s) or (1e300,))[0]   # noqa: E731
    for _ in range(200):
        c1 = b_ - g * (b_ - a_)
        c2 = a_ + g * (b_ - a_)
        if f(c1) < f(c2):
            b_ = c2
        else:
            a_ = c1
    x = (a_ + b_) / 2
    st = _state(x, nu, a, sci, mb, s)
    return (st[0], st[1], x) if st and st[0] <= best[0] else best
