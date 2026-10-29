"""Plotly figures for the app pages (no Streamlit imports, so they can be tested on their own).

Lengths in m, stresses in kPa (strength explorer: MPa), forces in kN (kN/m for 2D slices).
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

import engine as E

COLORSCALE = "Viridis"
GROUND = "#8c7b6b"
ROCK = "rgba(166, 146, 120, 0.25)"
WATER = "#2b7bba"
SLIP = "#c0392b"
FORCE = {"W": "#444444", "N": "#1f77b4", "U": WATER, "S": "#d62728", "n": "#7f7f7f", "E": "#7f7f7f", "P": WATER}

FIELDS = {
    "sigma_n": ("sigma_n (kPa)", "kPa"),
    "phi": ("phi_i (deg)", "deg"),
    "c": ("c_i (kPa)", "kPa"),
    "tau": ("tau_s (kPa)", "kPa"),
    "u": ("u (kPa)", "kPa"),
    "resistance": ("Contribution to resistance tau_s A (kN)", "kN"),
    "tension": ("Tension flag (sigma_n < 0)", ""),
}


def field_values(r: E.LemResult, name: str) -> np.ndarray:
    """Per-element value of a display field (NaN where the strength is not available)."""
    m = r.mesh
    nan = np.full(m.count, np.nan)
    s = r.strength
    if name == "u":
        return m.u.astype(float)
    if name == "sigma_n":
        return r.sn if r.sn is not None else nan
    if s is None:
        return nan
    if name == "phi":
        return np.degrees(s.phi)
    if name == "c":
        return r.c_used if r.c_used is not None else s.c
    if name == "tau":
        return r.tau_used if r.tau_used is not None else s.tau
    if name == "resistance":
        tau = r.tau_used if r.tau_used is not None else s.tau
        return tau * m.area
    if name == "tension":
        sn = r.sn if r.sn is not None else nan
        flag = sn < 0
        if r.tensile is not None:
            flag = flag | r.tensile
        return flag.astype(float)
    raise ValueError(name)


def _colors(values: np.ndarray, vmin: float, vmax: float) -> list[str]:
    span = vmax - vmin if vmax > vmin else 1.0
    t = np.clip((np.nan_to_num(values, nan=vmin) - vmin) / span, 0, 1)
    return sample_colorscale(COLORSCALE, list(t))


def _hover_rows(r: E.LemResult, k: int, two_d: bool) -> str:
    m = r.mesh
    s = r.strength
    lines = [f"{'slice' if two_d else 'column'} {k}  (i = {int(m.i[k])}" + ("" if two_d else f", j = {int(m.j[k])}") + ")",
             f"x = {m.x[k]:.2f} m" + ("" if two_d else f", y = {m.y[k]:.2f} m"),
             f"h = {m.h[k]:.2f} m, W = {m.w[k]:.1f} kN" + ("/m" if two_d else ""),
             f"u = {m.u[k]:.2f} kPa"]
    if r.sn is not None:
        lines.append(f"sigma_n = {r.sn[k]:.2f} kPa")
    if s is not None and np.isfinite(s.phi[k]):
        tau = r.tau_used[k] if r.tau_used is not None else s.tau[k]
        c = r.c_used[k] if r.c_used is not None else s.c[k]
        lines += [f"phi_i = {math.degrees(s.phi[k]):.3f} deg", f"c_i = {c:.2f} kPa", f"tau_s = {tau:.2f} kPa"]
    elif s is not None:
        lines.append(s.reason(k))
    return "<br>".join(lines)


def _frame(inp: E.Inputs, m: E.Mesh) -> tuple[float, float, float]:
    """Plot extent (x from, x to, padding), scaled to the problem: metres for a small slope, hundreds for a pit."""
    crest = -inp.H / math.tan(math.radians(inp.beta))
    left, right = min(m.x1, crest), max(m.x0, 0.0)
    pad = 0.1 * max(right - left, inp.H)
    return left - pad, right + 2 * pad, pad


def _slip_x(inp: E.Inputs, m: E.Mesh, n: int) -> np.ndarray:
    """x samples of the slip surface: up to the toe, or up to the exit on the floor for a base circle."""
    return np.linspace(m.x1, m.x0 if inp.base_circle else min(m.x0, 0.0), n)


# ---------------------------------------------------------------------------------------------------------------
# page 2: cross-section y = 0
# ---------------------------------------------------------------------------------------------------------------
def cross_section(inp: E.Inputs, r2: E.LemResult, field: str = "sigma_n", max_polygons: int = 400) -> go.Figure:
    m = r2.mesh
    xa, xb, pad = _frame(inp, m)
    xs = np.linspace(xa, xb, 1200)
    zt = E.ground(inp, xs)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.r_[xs, xs[::-1]], y=np.r_[zt, np.full(len(xs), min(-pad / 2, inp.zc - inp.R - pad / 2))],
                             fill="toself", fillcolor=ROCK, line=dict(width=0), hoverinfo="skip", name="rock mass"))
    fig.add_trace(go.Scatter(x=xs, y=zt, mode="lines", line=dict(color=GROUND, width=2), name="ground surface",
                             hoverinfo="skip"))
    # slip circle (the part inside the mass solid, the rest dotted)
    th = np.linspace(0, 2 * np.pi, 1441)
    cx, cz = inp.xc + inp.R * np.cos(th), inp.zc + inp.R * np.sin(th)
    fig.add_trace(go.Scatter(x=cx, y=cz, mode="lines", line=dict(color=SLIP, width=1, dash="dot"),
                             name="slip circle (full)", hoverinfo="skip", visible="legendonly"))
    xb_ = _slip_x(inp, m, 1200)
    zb_ = E.sphere_z(xb_, 0.0, inp.xc, inp.zc, inp.R)
    with np.errstate(invalid="ignore"):
        inside = zb_ < E.ground(inp, xb_)
    fig.add_trace(go.Scatter(x=xb_, y=np.where(inside, zb_, np.nan), mode="lines", line=dict(color=SLIP, width=3),
                             name="slip surface", hoverinfo="skip"))
    if inp.wet:
        zw = np.minimum(E.water(inp, xs), zt)
        fig.add_trace(go.Scatter(x=xs, y=zw, mode="lines", line=dict(color=WATER, width=2, dash="dash"),
                                 name="piezometric line" if inp.piezo is not None else "water table", hoverinfo="skip"))
    notes = []
    if inp.ru is not None:
        notes.append(f"pore pressure u = r_u γ h, r_u = {inp.ru:g}")
    if inp.kh > 0:
        notes.append(f"seismic k_h = {inp.kh:g}: k_h W at each centroid, towards the pit")
    if notes:
        fig.add_annotation(xref="paper", yref="paper", x=0.01, y=0.02, showarrow=False, align="left",
                           text="<br>".join(notes), font=dict(color="#333333"))
    ck = m.crack
    if ck and ck["rows"]:
        rw = ck["rows"][0]
        fig.add_trace(go.Scatter(x=[ck["x"], ck["x"]], y=[rw["zbot"], rw["ztop"]], mode="lines",
                                 line=dict(color="#222222", width=2.5), hoverinfo="name",
                                 name=f"tension crack, {rw['ztop'] - rw['zbot']:.2f} m deep"))
        if rw["hw"] > 0:
            fig.add_trace(go.Scatter(x=[ck["x"], ck["x"]], y=[rw["zbot"], rw["zbot"] + rw["hw"]], mode="lines",
                                     line=dict(color=WATER, width=7), hoverinfo="name",
                                     name=f"water in the crack: {rw['hw']:.2f} m, P = {_big(ck['P'])} kN/m"))
    # slices
    vals = field_values(r2, field)
    finite = np.isfinite(vals)
    vmin = float(np.nanmin(vals)) if finite.any() else 0.0
    vmax = float(np.nanmax(vals)) if finite.any() else 1.0
    cols = _colors(vals, vmin, vmax)
    hover = [_hover_rows(r2, k, True) for k in range(m.count)]
    if m.count <= max_polygons:
        for k in range(m.count):
            xl, xr = m.x[k] - m.dx / 2, m.x[k] + m.dx / 2
            zl = E.sphere_z(xl, 0.0, inp.xc, inp.zc, inp.R)
            zr = E.sphere_z(xr, 0.0, inp.xc, inp.zc, inp.R)
            zl, zr = (m.zb[k] if not np.isfinite(v) else v for v in (zl, zr))
            gl, gr = E.ground(inp, xl), E.ground(inp, xr)
            fig.add_trace(go.Scatter(x=[xl, xr, xr, xl, xl], y=[zl, zr, gr, gl, zl], fill="toself",
                                     fillcolor=cols[k] if finite[k] else "rgba(200,200,200,0.5)",
                                     line=dict(color="rgba(60,60,60,0.6)", width=0.6), mode="lines",
                                     hoveron="fills", text=hover[k], hoverinfo="text", showlegend=False))
    fig.add_trace(go.Scatter(x=m.x, y=m.zb, mode="markers", name="slice bases",
                             marker=dict(size=6, color=np.where(finite, vals, np.nan), colorscale=COLORSCALE,
                                         cmin=vmin, cmax=vmax, line=dict(width=0.5, color="white"),
                                         colorbar=dict(title=FIELDS[field][0], len=0.75)),
                             text=hover, hoverinfo="text"))
    # centre and radii to the ends of the slip surface
    fig.add_trace(go.Scatter(x=[inp.xc], y=[inp.zc], mode="markers+text", marker=dict(symbol="x", size=11, color=SLIP),
                             text=["centre"], textposition="middle right", name="circle centre",
                             hovertext=[f"centre ({inp.xc:.3f}, {inp.zc:.3f}), R = {inp.R:.3f} m"], hoverinfo="text"))
    xi = xb_[inside]
    if len(xi):
        for xe in (xi[0], xi[-1]):
            ze = float(E.sphere_z(xe, 0.0, inp.xc, inp.zc, inp.R))
            fig.add_trace(go.Scatter(x=[inp.xc, xe], y=[inp.zc, ze], mode="lines", showlegend=False, hoverinfo="skip",
                                     line=dict(color=SLIP, width=0.8, dash="dash")))
    fig.update_layout(xaxis_title="x (m), toe at 0, slope rises to -x", yaxis_title="z (m)",
                      yaxis=dict(scaleanchor="x", scaleratio=1), height=620, margin=dict(l=10, r=10, t=30, b=10),
                      legend=dict(orientation="h", y=1.02, x=0), hovermode="closest")
    return fig


# ---------------------------------------------------------------------------------------------------------------
# page 3: 3D view and plan
# ---------------------------------------------------------------------------------------------------------------
def _prisms(x, y, dx, dy, zb, zt, values):
    """Mesh3d vertex/face arrays for vertical prisms; one value per face (12 faces per prism)."""
    n = len(x)
    ox = np.array([-1, 1, 1, -1, -1, 1, 1, -1]) * dx / 2
    oy = np.array([-1, -1, 1, 1, -1, -1, 1, 1]) * dy / 2
    vx = (x[:, None] + ox[None, :]).ravel()
    vy = (y[:, None] + oy[None, :]).ravel()
    vz = np.c_[np.repeat(zb[:, None], 4, 1), np.repeat(zt[:, None], 4, 1)].ravel()
    tri = np.array([[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
                    [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    base = (np.arange(n) * 8)[:, None, None]
    faces = (tri[None, :, :] + base).reshape(-1, 3)
    return vx, vy, vz, faces[:, 0], faces[:, 1], faces[:, 2], np.repeat(values, 12)


def view3d(inp: E.Inputs, r3: E.LemResult, field: str = "sigma_n", mirror: bool = True,
           highlight: Optional[int] = None, max_prisms: int = 4000) -> go.Figure:
    m = r3.mesh
    vals = field_values(r3, field)
    finite = np.isfinite(vals)
    vmin = float(np.nanmin(vals)) if finite.any() else 0.0
    vmax = float(np.nanmax(vals)) if finite.any() else 1.0
    fig = go.Figure()
    xa, xb, _ = _frame(inp, m)
    ymass = float(np.max(m.y + m.dy / 2)) if m.count else 0.0
    yext = min(E.y_extent(inp), max(1.3 * ymass, 0.25 * (xb - xa)))
    ylo = -yext if mirror else 0.0
    gx = np.linspace(xa, xb, 90)
    gy = np.linspace(ylo, yext, 60)
    GX, GY = np.meshgrid(gx, gy)
    GZ = E.ground(inp, GX, GY)
    fig.add_trace(go.Surface(x=GX, y=GY, z=GZ, colorscale=[[0, "#bdb2a5"], [1, "#bdb2a5"]], showscale=False,
                             opacity=0.35, name="ground surface", hoverinfo="skip", showlegend=True))
    sx = _slip_x(inp, m, 90)
    sy = np.linspace(ylo, yext, 90)
    SX, SY = np.meshgrid(sx, sy)
    SZ = E.slip_z(inp, SX, SY)
    with np.errstate(invalid="ignore"):
        SZ = np.where(SZ < E.ground(inp, SX, SY), SZ, np.nan)
    fig.add_trace(go.Surface(x=SX, y=SY, z=SZ, colorscale=[[0, "#e8a598"], [1, "#e8a598"]], showscale=False,
                             opacity=0.55, name="slip surface", hoverinfo="skip", showlegend=True))
    xs, ys, zb, zt, vv = m.x, m.y, m.zb, m.zt, np.nan_to_num(vals, nan=vmin)
    if mirror:
        xs, ys, zb, zt, vv = (np.r_[xs, xs], np.r_[ys, -ys], np.r_[zb, zb], np.r_[zt, zt], np.r_[vv, vv])
    if len(xs) <= max_prisms:
        vx, vy, vz, i, j, k, fv = _prisms(xs, ys, m.dx, m.dy, zb, zt, vv)
        fig.add_trace(go.Mesh3d(x=vx, y=vy, z=vz, i=i, j=j, k=k, intensity=fv, intensitymode="cell",
                                colorscale=COLORSCALE, cmin=vmin, cmax=vmax, flatshading=True, opacity=1.0,
                                colorbar=dict(title=FIELDS[field][0], len=0.7), name="columns", hoverinfo="skip"))
    else:                                              # too many prisms to draw: coloured base cells only
        vx, vy, vz, i, j, k, fv = _prisms(xs, ys, m.dx, m.dy, zb, zb + 0.01, vv)
        fig.add_trace(go.Mesh3d(x=vx, y=vy, z=vz, i=i, j=j, k=k, intensity=fv, intensitymode="cell",
                                colorscale=COLORSCALE, cmin=vmin, cmax=vmax, flatshading=True,
                                colorbar=dict(title=FIELDS[field][0], len=0.7), name="column bases", hoverinfo="skip"))
    hover = [_hover_rows(r3, q, False) for q in range(m.count)]
    fig.add_trace(go.Scatter3d(x=m.x, y=m.y, z=m.zb, mode="markers", marker=dict(size=2, color="black"),
                               text=hover, hoverinfo="text", name="base centres (half-model)"))
    if highlight is not None and 0 <= highlight < m.count:
        q = highlight
        fig.add_trace(go.Scatter3d(x=[m.x[q], m.x[q]], y=[m.y[q], m.y[q]], z=[m.zb[q], m.zt[q]], mode="lines+markers",
                                   line=dict(color="red", width=8), marker=dict(size=4, color="red"),
                                   name=f"column {q}"))
    fig.update_layout(scene=dict(aspectmode="data", xaxis_title="x (m)", yaxis_title="y (m)", zaxis_title="z (m)",
                                 camera=dict(eye=dict(x=1.9, y=-1.1, z=0.8))),
                      height=680, margin=dict(l=0, r=0, t=30, b=0), legend=dict(orientation="h", y=1.02))
    return fig


def plan_view(r3: E.LemResult, field: str = "sigma_n", highlight: Optional[int] = None,
              mirror: bool = False) -> go.Figure:
    """Plan view of the columns (cells as a heatmap) with clickable centres and the highlighted column."""
    m = r3.mesh
    vals = field_values(r3, field)
    n, ny = m.n, m.ny
    grid = np.full((ny, n), np.nan)
    grid[m.j - 1, m.i - 1] = vals
    xs = m.x0 - (np.arange(1, n + 1) - 0.5) * m.dx
    ys = (np.arange(1, ny + 1) - 0.5) * m.dy
    if mirror:
        grid = np.r_[grid[::-1], grid]
        ys = np.r_[-ys[::-1], ys]
    fig = go.Figure(go.Heatmap(x=xs, y=ys, z=grid, colorscale=COLORSCALE, hoverinfo="skip",
                               colorbar=dict(title=FIELDS[field][0], len=0.8)))
    hover = [_hover_rows(r3, q, False) for q in range(m.count)]
    fig.add_trace(go.Scatter(x=m.x, y=m.y, mode="markers", marker=dict(size=7, color="rgba(255,255,255,0.75)",
                                                                       line=dict(width=0.5, color="black")),
                             text=hover, hoverinfo="text", customdata=np.arange(m.count), name="columns"))
    if highlight is not None and 0 <= highlight < m.count:
        q = highlight
        fig.add_shape(type="rect", x0=m.x[q] - m.dx / 2, x1=m.x[q] + m.dx / 2, y0=m.y[q] - m.dy / 2,
                      y1=m.y[q] + m.dy / 2, line=dict(color="red", width=3))
    ymax = float(m.y.max() + m.dy)
    fig.update_layout(xaxis_title="x (m)", yaxis_title="y (m)", height=480, margin=dict(l=10, r=10, t=30, b=10),
                      xaxis=dict(range=[float(m.x.min() - m.dx), float(m.x.max() + m.dx)]),
                      yaxis=dict(scaleanchor="x", scaleratio=1, range=[-ymax if mirror else 0.0, ymax]),
                      showlegend=False, clickmode="event+select")
    return fig


# ---------------------------------------------------------------------------------------------------------------
# page 4: free-body diagrams
# ---------------------------------------------------------------------------------------------------------------
def column_forces(inp: E.Inputs, r: E.LemResult, k: int) -> dict:
    """Forces on element k at the reported F (3D column, or 2D slice per metre).

    N' = sigma_n A, U = u A along the base normal n (into the mass); S = tau_s A / F along the resisting shear
    direction s = (-(z_c - z_b), 0, x_c - x)/rho, which lies in the base plane and in the x-z plane, so its vertical
    component is kappa. Vertical balance: W = (N' + U) g_z + S kappa.
    """
    m = r.mesh
    F = r.fos if r.fos is not None else r.f_trial
    if r.sn is None or r.tau_used is None or F is None:
        raise ValueError("no converged state for this element")
    x, y, zb = m.x[k], m.y[k], m.zb[k]
    if inp.elong == 1 and not inp.mid_length > 0:
        n = np.array([inp.xc - x, -y, inp.zc - zb]) / inp.R
    else:                                    # ellipsoidal caps or a cylindrical middle: the mesh's base normal
        n = np.array([m.n_x[k], m.n_y[k], m.g[k]])
    s_dir = np.array([-(inp.zc - zb), 0.0, inp.xc - x]) / m.rho[k]
    N = r.sn[k] * m.area[k]
    U = m.u[k] * m.area[k]
    S = r.tau_used[k] * m.area[k] / F
    W = m.w[k]
    resid = W - (N + U) * m.g[k] - S * m.kap[k]
    lever_foot = np.array([inp.xc, y, inp.zc])
    return dict(W=W, N=N, U=U, S=S, F=F, n=n, s=s_dir, g=m.g[k], kap=m.kap[k], rho=m.rho[k], area=m.area[k],
                residual=resid, residual_rel=abs(resid) / W, base=np.array([x, y, zb]),
                top=np.array([x, y, m.zt[k]]), lever_foot=lever_foot,
                lever=float(np.linalg.norm(np.array([x, y, zb]) - lever_foot)))


def _arrow3d(fig, start, vec, color, name, width=6):
    end = start + vec
    fig.add_trace(go.Scatter3d(x=[start[0], end[0]], y=[start[1], end[1]], z=[start[2], end[2]], mode="lines",
                               line=dict(color=color, width=width), name=name, hoverinfo="name"))
    L = float(np.linalg.norm(vec))
    if L > 0:
        fig.add_trace(go.Cone(x=[end[0]], y=[end[1]], z=[end[2]], u=[vec[0] / L], v=[vec[1] / L], w=[vec[2] / L],
                              anchor="tip", sizemode="absolute", sizeref=L * 0.18, showscale=False,
                              colorscale=[[0, color], [1, color]], hoverinfo="skip", showlegend=False))


def fbd3d(inp: E.Inputs, r3: E.LemResult, k: int) -> go.Figure:
    """Local free-body diagram of column k: prism, W, N', U, S and the base normal (forces scaled to W)."""
    m = r3.mesh
    f = column_forces(inp, r3, k)
    fig = go.Figure()
    vx, vy, vz, i, j, kk, _ = _prisms(np.array([m.x[k]]), np.array([m.y[k]]), m.dx, m.dy, np.array([m.zb[k]]),
                                      np.array([m.zt[k]]), np.array([0.0]))
    fig.add_trace(go.Mesh3d(x=vx, y=vy, z=vz, i=i, j=j, k=kk, color="#d9cbb8", opacity=0.45, name="column",
                            hoverinfo="skip", flatshading=True))
    L = max(m.h[k], m.dx, m.dy) * 1.2                  # length of W
    sc = L / f["W"]
    cen = (f["base"] + f["top"]) / 2
    _arrow3d(fig, cen, np.array([0, 0, -f["W"] * sc]), FORCE["W"], f"W = {f['W']:.1f} kN")
    _arrow3d(fig, f["base"] - f["n"] * f["N"] * sc, f["n"] * f["N"] * sc, FORCE["N"], f"N' = {f['N']:.1f} kN")
    if f["U"] > 0:
        off = f["s"] * 0.15 * L
        _arrow3d(fig, f["base"] + off - f["n"] * f["U"] * sc, f["n"] * f["U"] * sc, FORCE["U"], f"U = {f['U']:.1f} kN")
    _arrow3d(fig, f["base"], f["s"] * f["S"] * sc, FORCE["S"], f"S = {f['S']:.1f} kN")
    _arrow3d(fig, f["base"], f["n"] * 0.5 * L, FORCE["n"], "base normal n", width=3)
    fig.update_layout(scene=dict(aspectmode="data", xaxis_title="x (m)", yaxis_title="y (m)", zaxis_title="z (m)"),
                      height=560, margin=dict(l=0, r=0, t=30, b=0), legend=dict(orientation="h", y=1.05))
    return fig


def lever_view(inp: E.Inputs, r3: E.LemResult, k: int) -> go.Figure:
    """Elevation (x-z) at the column's y: the sphere section (radius rho), the rotation axis (parallel to y through
    (x_c, z_c)), the lever arm rho to the base point and the column."""
    m = r3.mesh
    f = column_forces(inp, r3, k)
    rho = f["rho"]
    th = np.linspace(np.pi, 2 * np.pi, 400)
    fig = go.Figure()
    xa, xb, _ = _frame(inp, m)
    xs = np.linspace(xa, xb, 500)
    fig.add_trace(go.Scatter(x=xs, y=E.ground(inp, xs, np.full_like(xs, m.y[k])), mode="lines",
                             line=dict(color=GROUND, width=2), name=f"ground at y = {m.y[k]:.1f} m", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=inp.xc + rho * np.cos(th), y=inp.zc + rho * np.sin(th), mode="lines",
                             line=dict(color=SLIP, width=2), name=f"sphere section at y = {m.y[k]:.1f} m",
                             hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=[inp.xc], y=[inp.zc], mode="markers+text", marker=dict(symbol="circle-x", size=13,
                                                                                         color=SLIP),
                             text=["rotation axis (parallel to y)"], textposition="top center", name="rotation axis",
                             hoverinfo="name"))
    fig.add_trace(go.Scatter(x=[inp.xc, m.x[k]], y=[inp.zc, m.zb[k]], mode="lines+text", line=dict(color="black",
                                                                                                    dash="dash"),
                             text=["", f"rho = {rho:.1f} m"], textposition="bottom left", name="lever arm rho"))
    fig.add_trace(go.Scatter(x=[m.x[k] - m.dx / 2, m.x[k] + m.dx / 2, m.x[k] + m.dx / 2, m.x[k] - m.dx / 2,
                                m.x[k] - m.dx / 2], y=[m.zb[k], m.zb[k], m.zt[k], m.zt[k], m.zb[k]], fill="toself",
                             fillcolor="rgba(214,39,40,0.3)", line=dict(color="red"), name=f"column {k}"))
    fig.update_layout(xaxis_title="x (m)", yaxis_title="z (m)", yaxis=dict(scaleanchor="x", scaleratio=1),
                      height=480, margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.05))
    return fig


def fbd2d(inp: E.Inputs, r2: E.LemResult, k: int) -> go.Figure:
    """Free-body diagram of slice k (per metre): W, N', U, S, base normal."""
    m = r2.mesh
    f = column_forces(inp, r2, k)
    xl, xr = m.x[k] - m.dx / 2, m.x[k] + m.dx / 2
    zl, zr = (float(E.sphere_z(v, 0.0, inp.xc, inp.zc, inp.R)) for v in (xl, xr))
    gl, gr = (float(E.ground(inp, v)) for v in (xl, xr))
    fig = go.Figure(go.Scatter(x=[xl, xr, xr, xl, xl], y=[zl, zr, gr, gl, zl], fill="toself",
                               fillcolor="rgba(217,203,184,0.6)", line=dict(color="black"), name=f"slice {k}",
                               hoverinfo="skip"))
    L = max(m.h[k], m.dx) * 1.0
    sc = L / f["W"]
    base = np.array([m.x[k], m.zb[k]])
    cen = np.array([m.x[k], (m.zb[k] + m.zt[k]) / 2])
    n2 = np.array([f["n"][0], f["n"][2]])
    s2 = np.array([f["s"][0], f["s"][2]])
    arrows = [("W", cen, np.array([0, -f["W"] * sc]), FORCE["W"]),
              ("N'", base - n2 * f["N"] * sc, n2 * f["N"] * sc, FORCE["N"]),
              ("S", base, s2 * f["S"] * sc, FORCE["S"])]
    if f["U"] > 0:
        arrows.append(("U", base + s2 * 0.12 * L - n2 * f["U"] * sc, n2 * f["U"] * sc, FORCE["U"]))
    for name, st_, vec, col in arrows:
        end = st_ + vec
        fig.add_annotation(x=end[0], y=end[1], ax=st_[0], ay=st_[1], xref="x", yref="y", axref="x", ayref="y",
                           showarrow=True, arrowhead=3, arrowsize=1.2, arrowwidth=2.5, arrowcolor=col, text="")
        fig.add_annotation(x=end[0], y=end[1], text=f"<b>{name}</b>", showarrow=False, font=dict(color=col, size=14),
                           xshift=12, yshift=8)
    fig.add_trace(go.Scatter(x=[base[0], base[0] + n2[0] * 0.5 * L], y=[base[1], base[1] + n2[1] * 0.5 * L],
                             mode="lines", line=dict(color=FORCE["n"], dash="dot"), name="base normal"))
    fig.update_layout(xaxis_title="x (m)", yaxis_title="z (m)", yaxis=dict(scaleanchor="x", scaleratio=1),
                      height=480, margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.05))
    return fig


# ---------------------------------------------------------------------------------------------------------------
# page 4: free-body diagram of the whole mass (E.mass_forces); as renderMass in webapp/index.html
# ---------------------------------------------------------------------------------------------------------------
MASS_NAMES = {"W": "W, weight", "N": "N′, effective normal force", "U": "U, pore-water force", "S": "S, mobilised shear",
              "E": "E = k_h W, seismic force", "P": "P, water in the tension crack"}


def _big(v: float) -> str:
    """Whole kN or kN·m with thousands separators; smaller values to 6 significant digits."""
    return f"{round(v):,}" if abs(v) >= 1e5 else f"{v:,.6g}"


def mass_groups(mf: dict, m: E.Mesh, limit: int) -> tuple[list[dict], int]:
    """Element forces summed over groups of neighbouring slices (2D) or blocks of b x b columns (3D) so that at most
    `limit` arrows of each kind are drawn. Display only: the table and the balance use every element."""
    n, el = mf["count"], mf["el"]
    if m.dim == 2:
        size = max(1, math.ceil(n / limit))
        keys = [str(k // size) for k in range(n)]
    else:
        size = 1

        def key(k, b):
            return f"{(int(m.i[k]) - 1) // b},{(int(m.j[k]) - 1) // b}"
        while len({key(k, size) for k in range(n)}) > limit:
            size += 1
        keys = [key(k, size) for k in range(n)]
    by: dict[str, list[int]] = {}
    for k, q in enumerate(keys):
        by.setdefault(q, []).append(k)
    groups = []
    for ks in by.values():
        ks = np.array(ks)
        w = el["W"][ks]
        W = float(w.sum())
        cen = np.array([np.sum(w * el["x"][ks]), np.sum(w * el["y"][ks]), np.sum(w * el["zm"][ks])]) / W
        N = np.array([np.sum(el["N"][ks] * el[c][ks]) for c in ("nx", "ny", "nz")])
        U = np.array([np.sum(el["U"][ks] * el[c][ks]) for c in ("nx", "ny", "nz")])
        S = np.array([np.sum(el["S"][ks] * el["sx"][ks]), 0.0, np.sum(el["S"][ks] * el["sz"][ks])])
        Ev = np.array([float(np.sum(el["E"][ks])) if "E" in el else 0.0, 0.0, 0.0])
        b = int(ks[np.argmin((el["x"][ks] - cen[0]) ** 2 + (el["y"][ks] - cen[1]) ** 2)])
        if m.dim == 2:
            name = f"slice {ks[0]}" if len(ks) == 1 else f"slices {ks[0]}–{ks[-1]}"
        else:
            name = (f"column {ks[0]} (i = {m.i[ks[0]]}, j = {m.j[ks[0]]})" if len(ks) == 1
                    else f"{len(ks)} columns around i = {m.i[b]}, j = {m.j[b]}")
        groups.append(dict(name=name, W=W, cen=cen, base=np.array([el["x"][b], el["y"][b], el["zb"][b]]),
                           tan=np.array([el["sx"][b], 0.0, el["sz"][b]]),
                           nrm=np.array([el["nx"][b], el["ny"][b], el["nz"][b]]), N=N, U=U, S=S, E=Ev))
    return groups, size


def _mass_arrows(groups: list[dict], kind: str, sc: float, off: float, unit: str) -> list[tuple]:
    """(start, end, hover text) of every arrow of one kind: W and E from the centre of gravity, S from the base point,
    N' and U ending at the base point (pointing into the mass when compressive)."""
    out = []
    for g in groups:
        v = np.array([0.0, 0.0, -g["W"]]) if kind == "W" else g[kind]
        mag = float(np.linalg.norm(v))
        if not mag > 0:
            continue
        vec = v * sc
        if kind in ("W", "E"):
            a, b = g["cen"], g["cen"] + vec
        elif kind == "S":
            a, b = g["base"], g["base"] + vec
        else:
            p = g["base"] + g["tan"] * off if kind == "U" else g["base"]
            a, b = p - vec, p
        tensile = kind == "N" and float(v @ g["nrm"]) < 0
        out.append((a, b, f"{g['name']}<br>{'N′' if kind == 'N' else kind} = {_big(mag)} {unit}"
                          + (" (tensile: pulls on the mass)" if tensile else "")))
    return out


def _mass_setup(mf: dict, m: E.Mesh, factor: float, limit: int):
    groups, size = mass_groups(mf, m, limit)
    span = max(float(m.x.max() - m.x.min() + m.dx), float(m.zt.max() - m.zb.min()), 2 * float(m.y.max()))
    peak = max(max(g["W"], *(float(np.linalg.norm(g[k])) for k in ("N", "U", "S", "E"))) for g in groups)
    t = mf["totals"]
    kinds = [k for k in ("W", "N", "U", "S", "E") if (k != "U" or t["U"]["x"] != 0 or t["U"]["z"] != 0)
             and (k != "E" or ("E" in t and t["E"]["x"] != 0))]
    return groups, size, span, factor * span / peak, 0.025 * span, kinds


def _mass_legend(mf: dict, k: str, unit: str) -> str:
    t = mf["weight"] if k == "W" else math.hypot(mf["totals"][k]["x"], mf["totals"][k]["z"])
    return f"{MASS_NAMES[k]} (Σ = {_big(t)} {unit})"


def mass_fbd2d(inp: E.Inputs, r2: E.LemResult, mf: dict) -> tuple[go.Figure, int]:
    """Every slice's W, N', U and S on the central section (per metre), the mass, the water table and the centre O.
    Returns the figure and the display group size (1 = every slice drawn)."""
    m = r2.mesh
    groups, size, span, sc, off, kinds = _mass_setup(mf, m, 0.14, 50)
    xa, xb, _ = _frame(inp, m)
    xs = np.linspace(xa, xb, 400)
    xe = _slip_x(inp, m, 600)
    arc = E.sphere_z(xe, 0.0, inp.xc, inp.zc, inp.R)
    with np.errstate(invalid="ignore"):
        keep = E.ground(inp, xe) > arc                                        # ground above the arc
    xe, arc = xe[keep], arc[keep]
    gz = E.ground(inp, xe)
    fig = go.Figure(go.Scatter(x=np.r_[xe, xe[::-1]], y=np.r_[gz, arc[::-1]], fill="toself", fillcolor=ROCK,
                               line=dict(color=GROUND, width=1), mode="lines", name="sliding mass", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=xs, y=E.ground(inp, xs), mode="lines", line=dict(color=GROUND, width=2),
                             name="ground", hoverinfo="skip", showlegend=False))
    if inp.wet:
        fig.add_trace(go.Scatter(x=xs, y=np.minimum(E.water(inp, xs), E.ground(inp, xs)),
                                 mode="lines", line=dict(color=WATER, width=1.5, dash="dash"),
                                 name="piezometric line" if inp.piezo is not None else "water table", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=xe, y=arc, mode="lines", line=dict(color=SLIP, width=2), name="slip surface",
                             hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=[inp.xc], y=[inp.zc], mode="markers+text", marker=dict(symbol="circle-x", size=12, color=SLIP),
                             text=["O, rotation axis"], textposition="middle right", name="centre", showlegend=False,
                             hovertext=[f"centre ({inp.xc:.2f}, {inp.zc:.2f}), R = {inp.R:.2f} m"], hoverinfo="text"))
    pts = [np.c_[xe, gz], np.c_[xe, arc]]
    for k in kinds:
        x, y, size_, text = [], [], [], []
        for a, b, t in _mass_arrows(groups, k, sc, off, "kN/m"):
            x += [a[0], b[0], None]
            y += [a[2], b[2], None]
            size_ += [0, 10, 0]
            text += [t, t, ""]
            pts.append(np.array([[a[0], a[2]], [b[0], b[2]]]))
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", text=text, hoverinfo="text", name=_mass_legend(mf, k, "kN/m"),
                                 line=dict(color=FORCE[k], width=1.6),
                                 marker=dict(symbol="arrow", size=size_, angleref="previous", color=FORCE[k])))
    ck = mf.get("crack")
    if ck and ck["P"] > 0:                     # water in the tension crack: one force on its face
        rw, vec = ck["rows"][0], ck["P"] * sc
        fig.add_trace(go.Scatter(x=[ck["x"] - vec, ck["x"]], y=[rw["zP"], rw["zP"]], mode="lines+markers",
                                 text=[f"P = {_big(ck['P'])} kN/m"] * 2, hoverinfo="text",
                                 name=f"{MASS_NAMES['P']} (Σ = {_big(ck['P'] * mf['mult'])} kN/m)",
                                 line=dict(color=FORCE["P"], width=2.4),
                                 marker=dict(symbol="arrow", size=[0, 12], angleref="previous", color=FORCE["P"])))
        pts.append(np.array([[ck["x"] - vec, rw["zP"]], [ck["x"], rw["zP"]]]))
    p = np.vstack(pts)
    lo, hi = p.min(axis=0), p.max(axis=0)
    pad = 0.06 * (hi - lo)
    fig.update_layout(xaxis=dict(title="x (m), toe at 0", range=[lo[0] - pad[0], hi[0] + pad[0]]),
                      yaxis=dict(title="z (m)", scaleanchor="x", scaleratio=1, range=[lo[1] - pad[1], hi[1] + pad[1]]),
                      height=520, margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.02, yanchor="bottom"))
    return fig, size


def mass_fbd3d(inp: E.Inputs, r3: E.LemResult, mf: dict) -> tuple[go.Figure, int]:
    """The same forces on every column of the full mass (both halves), with the slip and ground surfaces."""
    m = r3.mesh
    groups, size, span, sc, off, kinds = _mass_setup(mf, m, 0.25, 120)
    bx, by = np.r_[m.x, m.x], np.r_[m.y, -m.y]
    fig = go.Figure(go.Mesh3d(x=bx, y=by, z=np.r_[m.zb, m.zb], color=SLIP, opacity=0.35, name="slip surface",
                              hoverinfo="skip", showlegend=True, delaunayaxis="z"))
    fig.add_trace(go.Mesh3d(x=bx, y=by, z=np.r_[m.zt, m.zt], color=GROUND, opacity=0.25, name="ground over the mass",
                            hoverinfo="skip", showlegend=True, delaunayaxis="z"))
    for k in kinds:
        x, y, z, text, tips, dirs = [], [], [], [], [], []
        for a, b, t in _mass_arrows(groups, k, sc, off, "kN"):
            for sg in (1.0, -1.0):
                A, B = np.array([a[0], sg * a[1], a[2]]), np.array([b[0], sg * b[1], b[2]])
                d = B - A
                x += [A[0], B[0], None]
                y += [A[1], B[1], None]
                z += [A[2], B[2], None]
                text += [t, t, ""]
                tips.append(B)
                dirs.append(d / np.linalg.norm(d))
        fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode="lines", text=text, hoverinfo="text", name=_mass_legend(mf, k, "kN"),
                                   legendgroup=k, line=dict(color=FORCE[k], width=3)))
        tips, dirs = np.array(tips), np.array(dirs)
        fig.add_trace(go.Cone(x=tips[:, 0], y=tips[:, 1], z=tips[:, 2], u=dirs[:, 0], v=dirs[:, 1], w=dirs[:, 2],
                              anchor="tip", sizemode="raw", sizeref=0.035 * span, showscale=False,
                              colorscale=[[0, FORCE[k]], [1, FORCE[k]]], hoverinfo="skip", showlegend=False, legendgroup=k))
    fig.update_layout(scene=dict(aspectmode="data", xaxis_title="x (m)", yaxis_title="y (m)", zaxis_title="z (m)"),
                      height=560, margin=dict(l=0, r=0, t=30, b=0), legend=dict(orientation="h", y=1.02, yanchor="bottom"))
    return fig, size


def force_polygon(mf: dict) -> go.Figure:
    """The totals head to tail (W, U, N', S); the dashed arrow is their resultant (zero for full equilibrium)."""
    unit = "kN" if mf["dim"] == 3 else "kN/m"
    fig = go.Figure()
    p = np.zeros(2)
    for k in ("W", "U", "N", "S", "E", "P"):
        t = mf["totals"].get(k)
        if t is None or not math.hypot(t["x"], t["z"]) > 1e-12 * mf["weight"]:
            continue
        q = p + np.array([t["x"], t["z"]])
        lab = f"Σ{'N′' if k == 'N' else k}"
        fig.add_trace(go.Scatter(x=[p[0], (p[0] + q[0]) / 2, q[0]], y=[p[1], (p[1] + q[1]) / 2, q[1]], mode="lines+markers+text",
                                 text=["", lab, ""], textposition="middle right", line=dict(color=FORCE[k], width=3),
                                 marker=dict(symbol="arrow", size=[0, 0, 14], angleref="previous", color=FORCE[k]),
                                 hovertext=[f"{lab}: x {_big(t['x'])}, z {_big(t['z'])} {unit}"] * 3, hoverinfo="text",
                                 showlegend=False))
        p = q
    net = mf["totals"]["net"]
    fig.add_trace(go.Scatter(x=[0.0, net["x"]], y=[0.0, net["z"]], mode="lines+markers",
                             line=dict(color="#555555", width=2, dash="dash"),
                             marker=dict(symbol="arrow", size=[0, 12], angleref="previous", color="#555555"),
                             name=f"resultant: x {_big(net['x'])}, z {net['z']:.2e} {unit}", hoverinfo="name"))
    fig.update_layout(title=dict(text="Force polygon (totals, head to tail)", font=dict(size=13)),
                      xaxis_title=f"horizontal x ({unit})", yaxis=dict(title=f"vertical z ({unit})", scaleanchor="x", scaleratio=1),
                      height=520, margin=dict(l=10, r=10, t=40, b=10), legend=dict(orientation="h", y=-0.15, yanchor="top"))
    return fig


def mass_rows(mf: dict) -> list[dict]:
    """The totals and their moments about the rotation axis, for a table."""
    t, M = mf["totals"], mf["moments"]
    uF, uM = ("kN", "kN·m") if mf["dim"] == 3 else ("kN/m", "kN·m/m")

    def f(v):
        return "0" if v == 0 else f"{v:.3e}" if abs(v) < 1e-6 * mf["weight"] else _big(v)
    rows = [("W, weight", t["W"], _big(M["W"]), "vertical, down, at each element's centre of gravity; drives the rotation"),
            ("N′ = σ_n A, effective normal force", t["N"], "0", "normal to the slip surface, through the sphere centre"),
            ("U = u A, pore-water force", t["U"], "0", "normal to the slip surface, through the sphere centre"),
            ("S = τ_s A / F, mobilised shear", t["S"], _big(M["S"]), f"along the slip surface, against the rotation (F = {mf['F']:.6f})")]
    if "E" in t:
        rows += [("E = k_h W, seismic force", t["E"], _big(M["E"]), "horizontal, towards the pit, at each element's centre of gravity"),
                 ("P, water in the tension crack", t["P"], _big(M["P"]), "horizontal, towards the pit, on the crack face at h_w/3")]
    rows.append(("Sum of all forces", t["net"], f"{M['net']:.3e}", "resultant on the mass"))
    return [{"Force": n, f"Horizontal x ({uF})": f(v["x"]), f"Vertical z ({uF})": f(v["z"]),
             f"Magnitude ({uF})": f(math.hypot(v["x"], v["z"])), f"Moment about the axis ({uM})": mom, "Direction": d}
            for n, v, mom, d in rows]


def mass_balance(mf: dict) -> list[str]:
    """The vertical, horizontal and moment balance of the whole mass, as markdown lines (as massBalance in the web
    app)."""
    t, M, W = mf["totals"], mf["moments"], mf["weight"]
    uF, uM = ("kN", "kN·m") if mf["dim"] == 3 else ("kN/m", "kN·m/m")
    el = "column" if mf["dim"] == 3 else "slice"
    up = t["N"]["z"] + t["U"]["z"] + t["S"]["z"]
    v_ok = abs(t["net"]["z"]) <= 1e-9 * W and mf["vertical_max_rel"] <= 1e-9
    lines = [f"**Vertical:** ΣW = {_big(W)} {uF} down; the vertical parts of N′, U and S add to {_big(up)} {uF} up. Net "
             f"{t['net']['z']:.3e} {uF} ({abs(t['net']['z']) / W:.2e} of the weight): every {el} balances vertically "
             f"(eq. 4.12; largest relative residual {mf['vertical_max_rel']:.2e}) → **{'PASS' if v_ok else 'FAIL'}**",
             f"**Horizontal:** net {_big(t['net']['x'])} {uF}, {100 * abs(t['net']['x']) / W:.2f}% of the weight. "
             f"{'Hungr' if mf['dim'] == 3 else 'Bishop'}'s {'method' if mf['dim'] == 3 else 'simplified method'} does not "
             "enforce horizontal force equilibrium (thesis §5.5, §5.9), so this part stays unbalanced. The forces between "
             f"{el}s act inside the mass and cancel in these sums; the method takes their vertical (shear) parts as zero."]
    moment = (f"**Moments about the rotation axis** (through x_c, z_c, parallel to y): the weight drives {_big(M['W'])} {uM}"
              + (f", the seismic forces {_big(M['E'])} {uM} and the crack water {_big(M['P'])} {uM}" if "E" in M else "")
              + f"; the mobilised shear resists {_big(-M['S'])} {uM}; N′ and U act through the axis. Net {M['net']:.3e} {uM} "
              f"({abs(M['net']) / M['W']:.2e} of the weight's moment). ")
    if mf["dim"] == 2 or mf["form"] == "FM":
        ok = abs(M["net"]) <= 1e-8 * M["W"]
        moment += (f"{'Bishop' + chr(39) + 's method (eq. 4.16)' if mf['dim'] == 2 else 'The moment form F_M (eq. 4.15)'} "
                   f"balances this moment → **{'PASS' if ok else 'FAIL'}**")
    else:
        moment += (f"F_H balances ΣS = ΣWκ instead (eq. 4.14; relative residual {abs(mf['method']['residual']):.2e}); the "
                   "moment form F_M balances this moment (eq. 4.15, thesis §5.5).")
    return lines + [moment, "Not in the model: surcharge and water pressure against the face. Seismic load and water in a "
                            "tension crack are included when set (inputs)."]


# ---------------------------------------------------------------------------------------------------------------
# page 1: critical-surface search (engine.critical_search; QA/QC finding R1)
# ---------------------------------------------------------------------------------------------------------------
SEARCH_RAMP = ["#0d366b", "#104281", "#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5", "#5598e7", "#6da7ec",
               "#86b6ef", "#9ec5f4", "#b7d3f6", "#cde2fb"]          # the web app's blue ramp, darkest = lowest FoS
CRITICAL_3D = "#1f77b4"


def _f2(v: float) -> str:
    return "0.00" if abs(v) < 0.005 else f"{v:.2f}"


def _exit_text(inp: E.Inputs, p) -> str:
    crest = -inp.H / math.tan(math.radians(inp.beta))
    where = ("toe" if abs(p[0]) < 0.005 and abs(p[1]) < 0.005 else f"crest, {_f2(crest - p[0])} m behind the edge"
             if p[1] >= inp.H else "face" if p[1] > 0 else f"floor, {_f2(p[0])} m ahead of the toe")
    return f"({_f2(p[0])}, {_f2(p[1])}) {where}"


def _pct(a, b) -> str:
    if a is None or b is None:
        return "–"
    return f"{'+' if a >= b else '−'}{abs((a / b - 1) * 100):.1f} %"


def search_rows(inp: E.Inputs, res: dict) -> list[dict]:
    """The search results as table rows (as the web app's search table): FoS on the prescribed and the critical
    surface, the change, the critical circle and its exits, the mesh check and the surfaces tried."""
    P, dims = res["prescribed"], res["settings"]["dims"]
    cols = {2: "2D, Bishop", 3: "3D, Hungr"}

    def mesh(d, c):
        return f"{c['n']} x {c['ny']} columns" if d == 3 else f"{c['n']} slices"

    def fos_or_status(v, s):
        return f"{v:.6f}" if v is not None else (s or "withheld")

    def cell(d, what):
        r = res[f"r{d}"]
        found = r["fos"] is not None
        if what == "prescribed":
            return fos_or_status(P.get(f"fos{d}"), P.get(f"status{d}"))
        if what == "critical":
            return f"{r['fos']:.6f}" if found else "none found"
        if what == "change":
            return _pct(r["fos"], P.get(f"fos{d}")) if found else "–"
        if what == "circle":
            return f"{r['xc']:.3f}, {r['zc']:.3f}, {r['R']:.3f}" if found else "–"
        if what in ("front", "back"):
            return _exit_text(inp, r[what]) if found else "–"
        if what == "check":
            c = r["check"]
            if not c:
                return "–"
            return (f"{mesh(d, c)}: {fos_or_status(c['fos'], c['status'])}"
                    + (f" ({_pct(c['fos'], r['fos'])})" if c["change"] is not None else ""))
        if what == "tried":
            n = r["counts"]
            return f"{r['evaluations']}: {n['solved']} solved, {n['withheld']} withheld, {n['skipped']} skipped"
        return f"{r['ms'] / 1000:.1f} s"

    labels = [("prescribed", "FoS on the prescribed surface"), ("critical", "FoS on the critical surface"),
              ("change", "Change"), ("circle", "x_c, z_c, R (m)"), ("front", "Front exit (x, z) (m)"),
              ("back", "Back exit (x, z) (m)"), ("check", "Mesh check (finer mesh)"), ("tried", "Surfaces tried"),
              ("time", "Time")]
    return [{"": lab, **{cols[d]: cell(d, what) for d in dims}} for what, lab in labels]


def search_figure(inp: E.Inputs, res: dict) -> go.Figure:
    """The prescribed surface, the 2D critical circle and the 3D critical sphere's section y = 0 on the slope, with the
    centres of the circles the 2D search tried, coloured by their FoS (as the web app's search chart)."""
    P = res["prescribed"]
    crest = -inp.H / math.tan(math.radians(inp.beta))
    surfs = [("prescribed surface", P, "#7f7f7f", "dash")]
    for d, name, col in ((2, "2D critical circle", SLIP), (3, "3D critical sphere, section y = 0", CRITICAL_3D)):
        r = res.get(f"r{d}")
        if r and r["fos"] is not None:
            surfs.append((name, r, col, "solid"))
    xs, zs, arcs = [crest, 0.0], [0.0, inp.H], []
    for name, c, col, dash in surfs:
        x = np.linspace(c["xc"] - c["R"], c["xc"] + c["R"], 1601)
        z = E.sphere_z(x, 0.0, c["xc"], c["zc"], c["R"])
        with np.errstate(invalid="ignore"):
            ok = np.isfinite(z) & (z < E.ground(inp, x)) & ((x <= 0) | inp.base_circle)
        xs += list(x[ok]) + [c["xc"]]
        zs += list(z[ok]) + [c["zc"]]
        arcs.append(go.Scatter(x=x, y=np.where(ok, z, np.nan), mode="lines", name=name, hoverinfo="skip",
                               line=dict(color=col, width=2 if c is P else 3, dash=dash)))
    xa, xb, za, zb = min(xs), max(xs), min(zs), max(zs)
    pad = 0.08 * max(xb - xa, zb - za, inp.H)
    gx = np.r_[xa - 3 * (xb - xa) - pad, np.linspace(min(xa - pad, crest - 1e-9), max(xb + pad, 1e-9), 600),
               xb + 3 * (xb - xa) + pad]
    gz = E.ground(inp, gx)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.r_[gx, gx[::-1]], y=np.r_[gz, np.full(len(gx), za - 3 * (zb - za) - pad)],
                             fill="toself", fillcolor=ROCK, line=dict(width=0), hoverinfo="skip", name="rock mass",
                             showlegend=False))
    fig.add_trace(go.Scatter(x=gx, y=gz, mode="lines", line=dict(color=GROUND, width=2), name="ground surface",
                             hoverinfo="skip", showlegend=False))
    for d, name in ((2, "centres tried (2D)"), (3, "centres tried (3D)")):
        r = res.get(f"r{d}")
        if not r:
            continue
        c = r["cloud"]
        keep = [q for q, v in enumerate(c["fos"]) if v is not None]
        text = [f"centre ({c['xc'][q]:.2f}, {c['zc'][q]:.2f}), R = {c['R'][q]:.2f} m<br>FoS {d}D {c['fos'][q]:.4f}"
                for q in keep]
        if d == 2:
            best = r["fos"]
            marker = dict(size=6, color=[c["fos"][q] for q in keep], colorscale=[[i / (len(SEARCH_RAMP) - 1), col]
                                                                               for i, col in enumerate(SEARCH_RAMP)],
                          cmin=best, cmax=1.5 * best if best is not None else None, line=dict(width=0.5, color="white"),
                          colorbar=dict(title="FoS 2D", len=0.6))
            fig.add_trace(go.Scatter(x=[c["xc"][q] for q in keep], y=[c["zc"][q] for q in keep], mode="markers",
                                     name=name, marker=marker, text=text, hoverinfo="text"))
        else:
            fig.add_trace(go.Scatter(x=[c["xc"][q] for q in keep], y=[c["zc"][q] for q in keep], mode="markers",
                                     name=name, visible="legendonly", text=text, hoverinfo="text",
                                     marker=dict(size=6, color="rgba(0,0,0,0)", line=dict(width=1.2, color=CRITICAL_3D))))
    for a in arcs:
        fig.add_trace(a)
    for name, c, col, _ in surfs:
        fos = c.get("fos")
        fig.add_trace(go.Scatter(x=[c["xc"]], y=[c["zc"]], mode="markers", showlegend=False, name=name,
                                 marker=dict(symbol="x", size=9 if c is P else 11, color=col), hoverinfo="text",
                                 hovertext=[f"{name}: centre ({c['xc']:.3f}, {c['zc']:.3f}), R = {c['R']:.3f} m"
                                            + (f", FoS {fos:.4f}" if fos is not None else "")]))
    fig.update_layout(xaxis_title="x (m), toe at 0, slope rises to -x", yaxis_title="z (m)",
                      xaxis=dict(range=[xa - pad, xb + pad]),
                      yaxis=dict(scaleanchor="x", scaleratio=1, range=[za - pad, zb + pad]), height=560,
                      margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.02, x=0),
                      hovermode="closest")
    return fig


# ---------------------------------------------------------------------------------------------------------------
# page 5: strength explorer
# ---------------------------------------------------------------------------------------------------------------
def strength_curves(sigma_ci, mb, s, a, nu, s3_max, n=1500) -> go.Figure:
    """sigma1-sigma3 curves: EJC + Moscarda closure, EJC + original closure (admissible band shaded), GHB."""
    s3t = -s * sigma_ci / mb
    s3 = np.r_[s3t + (0 - s3t) * np.logspace(-10, 0, 300), np.linspace(0, s3_max, n)[1:]]
    fig = go.Figure()
    band_x = np.linspace(0, s3_max, 200)
    lo, hi = band_x / (2 * nu), band_x * (1 - nu) / nu
    fig.add_trace(go.Scatter(x=np.r_[band_x, band_x[::-1]], y=np.r_[lo, hi[::-1]], fill="toself",
                             fillcolor="rgba(255,165,0,0.15)", line=dict(width=0), hoverinfo="skip",
                             name="original closure: ordering band sigma3/(2 nu) <= sigma1 <= sigma3 (1-nu)/nu"))
    for model, name, color, dash in (("ejc_moscarda", "EJC, Moscarda closure", "#1f77b4", "solid"),
                                     ("ejc_original", "EJC, original closure (admissible part)", "#ff7f0e", "solid"),
                                     ("ghb", "GHB (sigma2 = sigma3)", "#2ca02c", "dash")):
        ps = E.principal_states(s3, E.Strength(sigma_ci, mb, s, a, model, nu))
        ok = ps.code == E.OK
        fig.add_trace(go.Scatter(x=s3, y=np.where(ok, ps.s1, np.nan), mode="lines", name=name,
                                 line=dict(color=color, dash=dash, width=2.5),
                                 hovertemplate="sigma3 = %{x:.4f} MPa<br>sigma1 = %{y:.4f} MPa"))
    fig.add_trace(go.Scatter(x=[s3t, s3_max], y=[s3t, s3_max], mode="lines", line=dict(color="grey", width=1),
                             name="sigma1 = sigma3"))
    fig.update_layout(xaxis_title="sigma3 (MPa)", yaxis_title="sigma1 (MPa)", height=520,
                      margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.08))
    return fig


def mohr(sigma_ci, mb, s, a, nu, sn_kpa: float, models=("ejc_moscarda", "ghb", "ejc_original")) -> tuple[go.Figure, list]:
    """Mohr circles and tangent lines at sigma_n (kPa) for each model, with each model's envelope. Stresses in MPa."""
    colors = {"ejc_moscarda": "#1f77b4", "ghb": "#2ca02c", "ejc_original": "#ff7f0e"}
    names = {"ejc_moscarda": "EJC Moscarda", "ghb": "GHB", "ejc_original": "EJC original"}
    fig = go.Figure()
    rows = []
    xmax = 0.0
    for model in models:
        st = E.Strength(sigma_ci, mb, s, a, model, nu)
        r = E.solve_strength([sn_kpa], st)
        row = dict(model=names[model], state=r.reason(0))
        if r.accepted[0]:
            s1, s3 = float(r.s1[0]), float(r.s3[0])
            cen, rad = (s1 + s3) / 2, (s1 - s3) / 2
            th = np.linspace(0, np.pi, 200)
            fig.add_trace(go.Scatter(x=cen + rad * np.cos(th), y=rad * np.sin(th), mode="lines",
                                     line=dict(color=colors[model], width=1.5), name=f"{names[model]}: Mohr circle"))
            phi, c, tau, sn = float(r.phi[0]), float(r.c[0]) / 1000, float(r.tau[0]) / 1000, float(r.sn[0]) / 1000
            xl = np.array([min(s3, sn) - 0.1 * rad, s1 + 0.1 * rad])
            fig.add_trace(go.Scatter(x=xl, y=c + xl * math.tan(phi), mode="lines",
                                     line=dict(color=colors[model], dash="dash", width=1), name=f"{names[model]}: tangent"))
            fig.add_trace(go.Scatter(x=[sn], y=[tau], mode="markers", marker=dict(color=colors[model], size=9),
                                     name=f"{names[model]}: tangent point"))
            xmax = max(xmax, s1)
            row.update(sigma1_MPa=s1, sigma2_MPa=float(r.s2[0]), sigma3_MPa=s3, M=float(r.M[0]),
                       phi_deg=math.degrees(phi), c_kPa=c * 1000, tau_kPa=tau * 1000)
        rows.append(row)
    # envelopes tau_s(sigma_n) over the visible range
    if xmax > 0:
        for model in models:
            st = E.Strength(sigma_ci, mb, s, a, model, nu)
            s3 = st.sigma3_tensile + (xmax - st.sigma3_tensile) * np.logspace(-12, 0, 800)
            ps = E.principal_states(s3, st)
            ok = (ps.code == E.OK) & (ps.sn <= xmax)
            if ok.any():
                fig.add_trace(go.Scatter(x=ps.sn[ok], y=ps.tau[ok], mode="lines", line=dict(color=colors[model],
                                                                                              width=3),
                                         opacity=0.35, name=f"{names[model]}: envelope"))
    fig.add_vline(x=sn_kpa / 1000, line=dict(color="black", dash="dot", width=1))
    fig.update_layout(xaxis_title="sigma (MPa)", yaxis_title="tau (MPa)", yaxis=dict(scaleanchor="x", scaleratio=1),
                      height=560, margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.1))
    return fig, rows


def history(r: E.LemResult) -> go.Figure:
    h = np.array([row[:4] for row in r.history]) if r.history else np.zeros((0, 4))
    fig = go.Figure()
    if len(h):
        fig.add_trace(go.Scatter(x=h[:, 0], y=h[:, 1], mode="lines+markers", name="F"))
    fig.update_layout(xaxis_title="iteration", yaxis_title="trial F", height=260, margin=dict(l=10, r=10, t=10, b=10))
    return fig
