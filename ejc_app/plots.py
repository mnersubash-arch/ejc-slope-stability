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
FORCE = {"W": "#444444", "N": "#1f77b4", "U": WATER, "S": "#d62728", "n": "#7f7f7f"}

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


# ---------------------------------------------------------------------------------------------------------------
# page 2: cross-section y = 0
# ---------------------------------------------------------------------------------------------------------------
def cross_section(inp: E.Inputs, r2: E.LemResult, field: str = "sigma_n", max_polygons: int = 400) -> go.Figure:
    m = r2.mesh
    crest = -inp.H / math.tan(math.radians(inp.beta))
    xa = min(m.x1, crest) - 40.0
    xb = max(m.x0, 0.0) + 80.0
    xs = np.linspace(xa, xb, 1200)
    zt = E.ground_z(xs, inp.H, inp.beta)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.r_[xs, xs[::-1]], y=np.r_[zt, np.full(len(xs), min(-20.0, inp.zc - inp.R - 20))],
                             fill="toself", fillcolor=ROCK, line=dict(width=0), hoverinfo="skip", name="rock mass"))
    fig.add_trace(go.Scatter(x=xs, y=zt, mode="lines", line=dict(color=GROUND, width=2), name="ground surface",
                             hoverinfo="skip"))
    # slip circle (the part inside the mass solid, the rest dotted)
    th = np.linspace(0, 2 * np.pi, 1441)
    cx, cz = inp.xc + inp.R * np.cos(th), inp.zc + inp.R * np.sin(th)
    fig.add_trace(go.Scatter(x=cx, y=cz, mode="lines", line=dict(color=SLIP, width=1, dash="dot"),
                             name="slip circle (full)", hoverinfo="skip", visible="legendonly"))
    xb_ = np.linspace(m.x1, min(m.x0, 0.0), 1200)
    zb_ = E.sphere_z(xb_, 0.0, inp.xc, inp.zc, inp.R)
    inside = zb_ < E.ground_z(xb_, inp.H, inp.beta)
    fig.add_trace(go.Scatter(x=xb_, y=np.where(inside, zb_, np.nan), mode="lines", line=dict(color=SLIP, width=3),
                             name="slip surface", hoverinfo="skip"))
    if inp.wet:
        zw = np.minimum(E.water_z(xs, inp.Hw, inp.beta_w), zt)
        fig.add_trace(go.Scatter(x=xs, y=zw, mode="lines", line=dict(color=WATER, width=2, dash="dash"),
                                 name="water table", hoverinfo="skip"))
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
            gl, gr = E.ground_z(xl, inp.H, inp.beta), E.ground_z(xr, inp.H, inp.beta)
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
    ylo = -inp.R if mirror else 0.0
    gx = np.linspace(m.x1 - 30, max(m.x0, 0) + 60, 90)
    gy = np.linspace(ylo, inp.R, 60)
    GX, GY = np.meshgrid(gx, gy)
    GZ = E.ground_z(GX, inp.H, inp.beta)
    fig.add_trace(go.Surface(x=GX, y=GY, z=GZ, colorscale=[[0, "#bdb2a5"], [1, "#bdb2a5"]], showscale=False,
                             opacity=0.35, name="ground surface", hoverinfo="skip", showlegend=True))
    sx = np.linspace(m.x1, min(m.x0, 0.0), 90)
    sy = np.linspace(ylo, inp.R, 90)
    SX, SY = np.meshgrid(sx, sy)
    SZ = E.sphere_z(SX, SY, inp.xc, inp.zc, inp.R)
    SZ = np.where(SZ < E.ground_z(SX, inp.H, inp.beta), SZ, np.nan)
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
    n = m.n
    grid = np.full((n, n), np.nan)
    grid[m.j - 1, m.i - 1] = vals
    xs = m.x0 - (np.arange(1, n + 1) - 0.5) * m.dx
    ys = (np.arange(1, n + 1) - 0.5) * m.dy
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
    n = np.array([inp.xc - x, -y, inp.zc - zb]) / inp.R
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
    xs = np.linspace(min(m.x1, -inp.H / math.tan(math.radians(inp.beta))) - 30, max(m.x0, 0) + 60, 500)
    fig.add_trace(go.Scatter(x=xs, y=E.ground_z(xs, inp.H, inp.beta), mode="lines", line=dict(color=GROUND, width=2),
                             name="ground (any y)", hoverinfo="skip"))
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
    gl, gr = (float(E.ground_z(v, inp.H, inp.beta)) for v in (xl, xr))
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
