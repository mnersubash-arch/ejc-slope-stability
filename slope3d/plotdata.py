"""Plain-numpy data behind every plot (shared by the matplotlib report and the plotly GUI).

Nothing here recomputes the engine: it only arranges the final element states. A quantity that is withheld
(F of a run that is not ACCEPTED) or not evaluated (final stresses) is returned as ``None``.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from .engine import solve_model
from .geometry import base_height, ground_height, slip_surface_ends
from .models import CaseResult, Config, ModelRun


def section_data(case: CaseResult) -> dict:
    """Central cross-section: ground line, slip circle arc and the 2D slice boundaries."""
    cfg = case.config
    ends = slip_surface_ends(cfg)
    if not ends.get("supported"):
        return {"supported": False}
    x_left = min(ends["entry_x"] - 0.04 * cfg.radius, -cfg.height / cfg.slope - 0.1 * cfg.height)
    gx = np.array([x_left, -cfg.height / cfg.slope, 0.0, max(0.12 * cfg.height, 10.0)])
    gz = np.array([cfg.height, cfg.height, 0.0, 0.0])
    xs = np.linspace(ends["entry_x"], ends["exit_x"], 400)
    arc = (xs, base_height(cfg, xs))
    # the full circle outside the ground, for context
    xc_full = np.linspace(cfg.xc - cfg.radius, min(cfg.xc + cfg.radius, ends["exit_x"] + 0.5 * cfg.radius), 600)
    out = {"supported": True, "ground": (gx, gz), "arc": arc, "ends": ends,
           "circle_full": (xc_full, base_height(cfg, xc_full))}
    els = case.run2d.elements
    if els is not None:
        dx = (cfg.xc - np.sqrt(cfg.radius ** 2 - cfg.zc ** 2) - (cfg.xc - cfg.radius)) / case.mesh2d
        out["slices"] = dict(x=els.x, zb=els.zb, top=els.zb + els.h, dx=dx)
    return out


def _mirror(a: np.ndarray) -> np.ndarray:
    return np.concatenate([a, a])


def plan_data(case: CaseResult) -> Optional[dict]:
    """Final base normal stress per column, half-model mirrored about y = 0. ``None`` if not evaluated."""
    run = case.run3d
    if run.elements is None or not run.result.evaluated:
        return None
    e, cfg = run.elements, case.config
    dx = (cfg.xc - np.sqrt(cfg.radius ** 2 - cfg.zc ** 2) - (cfg.xc - cfg.radius)) / case.mesh3d
    dy = cfg.radius / case.mesh3d
    return dict(x=_mirror(e.x), y=np.concatenate([e.y, -e.y]), value=_mirror(e.final_n), dx=dx, dy=dy)


def surface_grid(case: CaseResult) -> Optional[dict]:
    """Slip-surface elevation z_b on the n x n grid, mirrored (NaN where there is no mass). Needs only the geometry."""
    run = case.run3d
    if run.elements is None:
        return None
    e, n = run.elements, case.mesh3d
    xs = np.full((n, n), np.nan)
    ys = np.full((n, n), np.nan)
    zs = np.full((n, n), np.nan)
    xs[e.ix, e.iy], ys[e.ix, e.iy], zs[e.ix, e.iy] = e.x, e.y, e.zb
    return dict(x=np.concatenate([xs[:, ::-1], xs], axis=1), y=np.concatenate([-ys[:, ::-1], ys], axis=1),
                z=np.concatenate([zs[:, ::-1], zs], axis=1))


def convergence_data(case: CaseResult) -> dict:
    """F against iteration, only for ACCEPTED dimensions (otherwise withheld: ``None``)."""
    out = {}
    for tag, run in (("3D", case.run3d), ("2D", case.run2d)):
        h = run.result.history
        out[tag] = (np.array([r.iteration for r in h]), np.array([r.f for r in h])) if (run.result.accepted and h) else None
    return out


def strength_section(case: CaseResult) -> dict:
    """phi'_i and c'_i along the central section: 2D slices and the 3D column row nearest to y = 0."""
    out = {}
    r2 = case.run2d
    if r2.elements is not None and r2.result.evaluated:
        e = r2.elements
        out["2D"] = (e.x, e.phi, e.cohesion)
    r3 = case.run3d
    if r3.elements is not None and r3.result.evaluated:
        e = r3.elements
        m = e.iy == 0
        out["3D"] = (e.x[m], e.phi[m], e.cohesion[m])
    return out


def mesh_convergence(cfg: Config, meshes3=(10, 15, 20, 30, 40, 60, 80), meshes2=(20, 40, 80, 160, 320, 640),
                     progress=None) -> list[dict]:
    """F against mesh for the current configuration. F is ``None`` unless that run is ACCEPTED."""
    rows = []
    for dim, meshes in ((3, meshes3), (2, meshes2)):
        for n in meshes:
            r = solve_model(cfg, dim, n).result
            row = dict(dim=f"{dim}D", n=n, F=r.fos, status=r.status.value, elements=r.active)
            rows.append(row)
            if progress:
                progress(row)
    return rows
