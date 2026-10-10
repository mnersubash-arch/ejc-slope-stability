"""Interactive plotly versions of the six plots (used by the GUI). Same data and withholding rules as ``report.py``."""
from __future__ import annotations

from typing import Optional

import numpy as np
import plotly.graph_objects as go

from .models import CaseResult
from .plotdata import convergence_data, plan_data, section_data, strength_section, surface_grid

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#8a8984"


def _empty(title: str, message: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=message, x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False, font=dict(size=16))
    fig.update_layout(title=title, xaxis=dict(visible=False), yaxis=dict(visible=False), height=320)
    return fig


def _layout(fig: go.Figure, title: str, xt: str, yt: str, height: int = 460) -> go.Figure:
    fig.update_layout(title=dict(text=title, x=0), xaxis_title=xt, yaxis_title=yt, height=height,
                      margin=dict(l=60, r=20, t=100, b=50), legend=dict(orientation="h", y=1.0, yanchor="bottom", x=0))
    return fig


def section(case: CaseResult) -> go.Figure:
    d = section_data(case)
    if not d["supported"]:
        return _empty("Central cross-section", "Unsupported geometry - nothing to draw")
    fig = go.Figure()
    gx, gz = d["ground"]
    fig.add_scatter(x=gx, y=gz, mode="lines", name="ground surface", line=dict(color="#0b0b0b", width=2))
    fig.add_scatter(x=d["circle_full"][0], y=d["circle_full"][1], mode="lines", name="slip circle (full)",
                    line=dict(color=GREY, width=1, dash="dash"))
    s = d.get("slices")
    if s is not None:
        for i, (x, zb, top) in enumerate(zip(s["x"], s["zb"], s["top"])):
            fig.add_shape(type="rect", x0=x - s["dx"] / 2, x1=x + s["dx"] / 2, y0=zb, y1=top, line=dict(color="white", width=1),
                          fillcolor="rgba(42,120,214,0.35)", layer="below")
    fig.add_scatter(x=d["arc"][0], y=d["arc"][1], mode="lines", name="slip surface", line=dict(color=ORANGE, width=3))
    e = d["ends"]
    fig.add_scatter(x=[e["entry_x"], e["exit_x"]], y=[e["entry_z"], e["exit_z"]], mode="markers+text", name="entry / exit",
                    marker=dict(color=ORANGE, size=10, line=dict(color="white", width=2)),
                    text=[f"entry x={e['entry_x']:.1f}", f"exit x={e['exit_x']:.1f}, {e['exit_z']:.1f} m above toe"],
                    textposition=["top right", "top right"])
    _layout(fig, f"Central cross-section, {case.mesh2d} slices", "x [m] (slope rises toward -x)", "z [m]", 520)
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


def plan(case: CaseResult) -> go.Figure:
    d = plan_data(case)
    if d is None:
        return _empty("Final base normal stress (plan)", f"sigma_n' not evaluated - 3D run is {case.run3d.result.status.value}")
    v = d["value"]
    fig = go.Figure()
    fig.add_scatter(x=d["x"], y=d["y"], mode="markers", name="columns",
                    marker=dict(symbol="square", size=16, color=v, colorscale="RdBu", cmid=0 if v.min() < 0 else None,
                                colorbar=dict(title="sigma_n' [kPa]<br>(negative = tension)"), line=dict(width=0)),
                    hovertemplate="x=%{x:.1f} m<br>y=%{y:.1f} m<br>sigma_n'=%{marker.color:.2f} kPa<extra></extra>")
    r3 = case.run3d.result
    _layout(fig, f"Final base normal stress (half-model mirrored): {r3.tension} columns with sigma_n' < 0, min {r3.min_normal:.1f} kPa",
            "x [m]", "y [m]", 560)
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


def surface3d(case: CaseResult) -> go.Figure:
    g = surface_grid(case)
    if g is None:
        return _empty("Slip surface (3D)", "No elements - nothing to draw")
    fig = go.Figure(go.Surface(x=g["x"], y=g["y"], z=g["z"], colorscale="Blues", reversescale=True, showscale=False))
    fig.update_layout(title=dict(text="Slip surface (spherical base of the sliding mass)", x=0), height=560,
                      scene=dict(xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]", aspectmode="data"),
                      margin=dict(l=0, r=0, t=60, b=0))
    return fig


def convergence(case: CaseResult) -> go.Figure:
    data = convergence_data(case)
    if all(v is None for v in data.values()):
        return _empty("F against iteration", "FoS withheld - no ACCEPTED dimension")
    fig = go.Figure()
    for tag, col in (("3D", BLUE), ("2D", ORANGE)):
        if data[tag] is not None:
            it, f = data[tag]
            fig.add_scatter(x=it, y=f, mode="lines+markers", name=f"{tag}  F = {f[-1]:.4f}", line=dict(color=col, width=2))
        else:
            fig.add_scatter(x=[None], y=[None], mode="lines", name=f"{tag}  FoS withheld", line=dict(color=col))
    return _layout(fig, "F against iteration (start F = 1)", "iteration", "factor of safety F", 420)


def strength(case: CaseResult) -> go.Figure:
    d = strength_section(case)
    if not d:
        return _empty("phi'_i and c'_i along the central section", "Final strengths not evaluated")
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08)
    for tag, col, sym, lab in (("2D", ORANGE, "circle", "2D slices (y = 0)"), ("3D", BLUE, "square", "3D columns nearest y = 0")):
        if tag in d:
            x, phi, c = d[tag]
            o = np.argsort(x)
            fig.add_scatter(x=x[o], y=phi[o], name=lab, mode="lines+markers", marker=dict(symbol=sym), line=dict(color=col), row=1, col=1)
            fig.add_scatter(x=x[o], y=c[o], name=lab, mode="lines+markers", marker=dict(symbol=sym), line=dict(color=col),
                            showlegend=False, row=2, col=1)
    fig.update_yaxes(title_text="phi'_i [deg]", row=1, col=1)
    fig.update_yaxes(title_text="c'_i [kPa]", row=2, col=1)
    fig.update_xaxes(title_text="x [m]", row=2, col=1)
    fig.update_layout(title=dict(text="Instantaneous friction angle and cohesion along the central section", x=0), height=560,
                      margin=dict(l=60, r=20, t=100, b=50), legend=dict(orientation="h", y=1.0, yanchor="bottom", x=0))
    return fig


def mesh(rows: Optional[list[dict]]) -> go.Figure:
    if not rows:
        return _empty("Mesh convergence", "Not run yet - press 'Run mesh convergence'")
    fig = go.Figure()
    for tag, col in (("3D", BLUE), ("2D", ORANGE)):
        pts = [(r["n"], r["F"]) for r in rows if r["dim"] == tag and r["F"] is not None]
        if pts:
            n, f = zip(*pts)
            fig.add_scatter(x=n, y=f, mode="lines+markers", name=f"{tag} ({'n x n columns' if tag == '3D' else 'slices'})", line=dict(color=col))
        else:
            fig.add_scatter(x=[None], y=[None], name=f"{tag}: FoS withheld", line=dict(color=col))
    fig.update_xaxes(type="log")
    return _layout(fig, "Mesh convergence of F (ACCEPTED runs only)", "mesh (n, log scale)", "F", 420)
