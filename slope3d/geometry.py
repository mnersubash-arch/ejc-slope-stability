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


def ground_height(cfg: Config, x):
    """Ground surface on the central section: z = min(H, -x * slope) for x <= 0 (the slope rises toward -x)."""
    return np.minimum(cfg.height, -np.asarray(x, dtype=float) * cfg.slope)


def base_height(cfg: Config, x, y=0.0):
    """Slip-sphere base z_b = zc - sqrt(R^2 - (x - xc)^2 - y^2); NaN outside the sphere."""
    x = np.asarray(x, dtype=float)
    with np.errstate(invalid="ignore"):
        return cfg.zc - np.sqrt(cfg.radius ** 2 - (x - cfg.xc) ** 2 - y ** 2)


def slip_surface_ends(cfg: Config) -> dict:
    """Where the slip surface enters and exits the ground on the central section (y = 0), x <= 0.

    The mass exists where the ground is above the sphere base. ``entry`` is the first crossing from the back
    (z at the crest level or on the face), ``exit`` the last one before the toe-side cutoff at x = 0.
    Heights are above the toe level (z = 0). ``base_at_toe_plane`` is z_b at x = 0; a value above 0 means the
    circle passes above the toe plane, so the mass pinches out on the face (``exit_on_face``) before x = 0.
    """
    out: dict = {"supported": False}
    if geometry_issue(cfg):
        return out
    x1 = cfg.xc - cfg.radius
    xs = np.linspace(x1, 0.0, 200001)
    diff = ground_height(cfg, xs) - base_height(cfg, xs)       # > 0 where there is mass
    ok = np.isfinite(diff)
    pos = ok & (diff > 0)
    out["base_at_toe_plane"] = float(base_height(cfg, 0.0))
    if not pos.any():
        return out

    def refine(i0: int, i1: int) -> float:
        a, b = xs[i0], xs[i1]
        fa = float(ground_height(cfg, a) - base_height(cfg, a))
        for _ in range(80):
            m = (a + b) / 2
            fm = float(ground_height(cfg, m) - base_height(cfg, m))
            if (fm > 0) == (fa > 0):
                a, fa = m, fm
            else:
                b = m
        return (a + b) / 2

    first, last = int(np.argmax(pos)), len(pos) - 1 - int(np.argmax(pos[::-1]))
    x_in = refine(max(first - 1, 0), first) if first > 0 else float(xs[0])
    if last < len(xs) - 1:
        x_out = refine(last, last + 1)
        reaches_face = True
    else:
        x_out, reaches_face = 0.0, False                        # mass still present at the toe plane
    out.update(
        supported=True,
        entry_x=float(x_in), entry_z=float(ground_height(cfg, x_in)),
        exit_x=float(x_out), exit_z=float(base_height(cfg, x_out)),
        exit_on_face=reaches_face,
        crest_x=-cfg.height / cfg.slope,
    )
    out["exit_height_above_toe"] = out["exit_z"]
    out["exit_behind_toe"] = -float(x_out)                      # horizontal distance of the exit behind the toe (x = 0)
    out["exit_distance_from_crest"] = float(x_out - out["crest_x"])  # horizontal distance below the crest edge
    return out
