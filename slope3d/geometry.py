"""Geometry guard and element builder (``Thesis_Columns.GeometryIssue`` and ``BuildElements``)."""
from __future__ import annotations

import math

import numpy as np

from .models import MESH2D_RANGE, MESH3D_RANGE, Config, Elements

PI_C = 3.141592653589793        # Thesis_Columns.bas line 10


class GeometryError(ValueError):
    """Raised for refused geometry, an invalid mesh, or a mesh with no active element."""


def geometry_issue(cfg: Config) -> str:
    """Return ``"UNSUPPORTED GEOMETRY: ..."`` if the case is outside the supported model, else ``""``.

    The model is the documented positive-centre, driving-side, toe-cutoff sphere. Refuse rather than compute wrong numbers.
    """
    vals = (cfg.xc, cfg.zc, cfg.radius, cfg.height, cfg.slope, cfg.gamma)
    if not all(isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool) and math.isfinite(v)
               for v in vals):
        return "UNSUPPORTED GEOMETRY: finite numeric geometry inputs required."
    if cfg.radius <= cfg.zc or cfg.zc <= 0 or cfg.height <= 0 or cfg.slope <= 0 or cfg.gamma <= 0:
        return "UNSUPPORTED GEOMETRY: require R > zc > 0 and positive height, slope and unit weight."
    if cfg.xc <= 0 or cfg.xc >= cfg.radius:
        return "UNSUPPORTED GEOMETRY: this driving-side toe-cutoff model requires 0 < xc < R."
    x_toe = cfg.xc - math.sqrt(cfg.radius ** 2 - cfg.zc ** 2)
    # Permit the published reference circle's 2.62 mm rounding mismatch; do not move its grid.
    if x_toe < -0.00001 * cfg.radius:
        return "UNSUPPORTED GEOMETRY: the adopted grid would omit slope-side mass between its right bound and the toe."
    return ""


def mesh_issue(dim: int, n: int) -> str:
    """Mesh limits: 3D n x n half-model with 2 <= n <= 400; 2D 2 <= n <= 20,000 slices."""
    if dim not in (2, 3) or isinstance(n, bool) or not isinstance(n, (int, np.integer)):
        return "invalid mesh (dimension must be 2 or 3 and the mesh an integer)"
    lo, hi = MESH3D_RANGE if dim == 3 else MESH2D_RANGE
    if n < 2 or n > 20000:
        return f"invalid mesh: {dim}D mesh must be {lo}-{hi}"
    if dim == 3 and n > 400:
        return "3D limit is 400 x 400; reduce mesh."
    return ""


def build_elements(cfg: Config, dim: int, n: int) -> tuple[Elements, float]:
    """Build the active elements in VBA order (x outer, y inner). Returns (elements, volume m^3).

    x0 = xc - sqrt(R^2 - zc^2), x1 = xc - R; dx = (x0 - x1)/n; dy = R/n (3D) or 1 (2D).
    Element centres x = x0 - (i + 1/2) dx, y = (j + 1/2) dy (3D) or 0 (2D). Skipped where
    q = R^2 - (x - xc)^2 - y^2 <= 0, or x > 0, or h <= 0. 3D volume is doubled (half-model).
    """
    issue = geometry_issue(cfg)
    if issue:
        raise GeometryError(issue)
    issue = mesh_issue(dim, n)
    if issue:
        raise GeometryError(issue)
    r = cfg.radius
    x0 = cfg.xc - math.sqrt(r ** 2 - cfg.zc ** 2)
    x1 = cfg.xc - r
    dx = (x0 - x1) / n
    dy = r / n if dim == 3 else 1.0
    ny = n if dim == 3 else 1
    ii, jj = np.meshgrid(np.arange(n), np.arange(ny), indexing="ij")     # i outer, j inner (row-major ravel)
    ii, jj = ii.ravel(), jj.ravel()
    xx = x0 - (ii + 0.5) * dx
    yy = (jj + 0.5) * dy if dim == 3 else np.zeros(len(ii))
    with np.errstate(invalid="ignore"):
        q = r ** 2 - (xx - cfg.xc) ** 2 - yy ** 2
        keep = (q > 0) & (xx <= 0)
        ii, jj, xx, yy, q = ii[keep], jj[keep], xx[keep], yy[keep], q[keep]
        sq = np.sqrt(q)
        base = cfg.zc - sq
        top = np.minimum(cfg.height, -xx * cfg.slope)
        thick = top - base
        keep = thick > 0
    ii, jj, xx, yy, q, sq, base, thick = (v[keep] for v in (ii, jj, xx, yy, q, sq, base, thick))
    if len(xx) == 0:
        raise GeometryError("No active elements.")
    ah = dx * dy
    g = sq / r
    rho = np.sqrt(r ** 2 - yy ** 2)
    pressure = np.full(len(xx), cfg.u_extra, dtype=float)
    if cfg.water_height > 0:
        water = np.minimum(-xx * math.tan(cfg.water_angle * PI_C / 180), cfg.water_height)
        pressure = pressure + np.where(water > base, 9.81 * (water - base), 0.0)
    els = Elements(
        x=xx, y=yy, zb=base, h=thick, ah=np.full(len(xx), ah), g=g, area=ah / g, rho=rho,
        kap=(cfg.xc - xx) / rho, weight=cfg.gamma * thick * ah, pressure=pressure,
        sn=(cfg.gamma * thick - pressure) * g ** 2, ix=ii, iy=jj)
    volume = float(np.cumsum(thick * ah)[-1])        # sequential sum, as the VBA accumulates
    if dim == 3:
        volume *= 2
    return els, volume
