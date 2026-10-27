"""Matplotlib figures and the PDF / HTML run report.

Rule kept everywhere: a factor of safety is drawn or printed only for an ACCEPTED dimension. Withheld results
show the words "FoS withheld"; stresses that were not evaluated show "not evaluated".
"""
from __future__ import annotations

import base64
import html
import io as _io
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages                    # noqa: E402
from matplotlib.colors import TwoSlopeNorm                              # noqa: E402
from matplotlib.patches import Rectangle                                # noqa: E402
from mpl_toolkits.mplot3d import Axes3D                                 # noqa: E402,F401
import numpy as np                                                      # noqa: E402

from . import __version__                                               # noqa: E402
from .io import inputs_table, summary_table                             # noqa: E402
from .models import CaseResult, Status                                  # noqa: E402
from .plotdata import (convergence_data, plan_data, section_data, strength_section,   # noqa: E402
                       surface_grid)

BLUE, ORANGE, GREY, INK = "#2a78d6", "#eb6834", "#8a8984", "#0b0b0b"
STATUS_COLOUR = {Status.ACCEPTED: "#008300", Status.NOT_ACCEPTED: "#e34948", Status.DIAGNOSTIC_ONLY: "#6b6a66"}


def _style(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, color="#e6e5e1", lw=0.6)
    ax.set_axisbelow(True)


def _placeholder(title: str, message: str):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=11)
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=13, color="#52514e", transform=ax.transAxes)
    return fig


def fig_section(case: CaseResult):
    """Central cross-section with the slip circle, the 2D slices and the entry / exit points."""
    d = section_data(case)
    if not d["supported"]:
        return _placeholder("Central cross-section", "Unsupported geometry - nothing to draw")
    fig, ax = plt.subplots(figsize=(9, 5))
    gx, gz = d["ground"]
    ax.fill_between(gx, gz, min(0.0, d["arc"][1].min()) - 0.05 * case.config.height, color="#efeee9", zorder=0)
    ax.plot(gx, gz, color=INK, lw=1.6, label="ground surface")
    ax.plot(*d["circle_full"], color=GREY, lw=0.8, ls="--", label="slip circle (full)")
    ax.plot(*d["arc"], color=ORANGE, lw=2.2, label="slip surface")
    s = d.get("slices")
    if s is not None:
        for x, zb, top in zip(s["x"], s["zb"], s["top"]):
            ax.add_patch(Rectangle((x - s["dx"] / 2, zb), s["dx"], top - zb, fc=BLUE, ec="white", lw=0.6, alpha=0.35))
    e = d["ends"]
    ax.plot([e["entry_x"]], [e["entry_z"]], "o", color=ORANGE, mec="white", ms=7)
    ax.plot([e["exit_x"]], [e["exit_z"]], "o", color=ORANGE, mec="white", ms=7)
    ax.annotate(f"exit  x = {e['exit_x']:.1f} m\n{e['exit_z']:.1f} m above toe", (e["exit_x"], e["exit_z"]),
                (e["exit_x"] + 0.08 * case.config.height, e["exit_z"] + 0.18 * case.config.height),
                arrowprops=dict(arrowstyle="-", color=GREY), fontsize=9)
    ax.annotate(f"entry  x = {e['entry_x']:.1f} m", (e["entry_x"], e["entry_z"]),
                (e["entry_x"] + 0.04 * case.config.height, e["entry_z"] + 0.2 * case.config.height),
                arrowprops=dict(arrowstyle="-", color=GREY), fontsize=9, ha="left")
    ax.set_aspect("equal")
    ax.set_xlabel("x [m] (slope rises toward -x)")
    ax.set_ylabel("z [m]")
    ax.set_title(f"Central cross-section, {case.mesh2d} slices", loc="left", fontsize=11)
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    _style(ax)
    fig.tight_layout()
    return fig


def fig_plan(case: CaseResult):
    """Plan heat-map of the final base normal stress (half-model mirrored about y = 0)."""
    d = plan_data(case)
    if d is None:
        return _placeholder("Final base normal stress (plan)",
                            f"sigma_n' not evaluated - 3D run is {case.run3d.result.status.value}")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    v = d["value"]
    norm = TwoSlopeNorm(vcenter=0.0, vmin=min(v.min(), -1e-9), vmax=max(v.max(), 1e-9)) if v.min() < 0 else None
    sc = ax.scatter(d["x"], d["y"], c=v, s=_marker_area(ax, fig, d), marker="s", cmap="RdBu" if norm else "Blues",
                    norm=norm, linewidths=0)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("final base normal stress sigma_n' [kPa]  (negative = tension)")
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]  (mirrored about y = 0)")
    r3 = case.run3d.result
    ax.set_title(f"Final base normal stress, {case.mesh3d}x{case.mesh3d} columns\n{r3.tension} columns with sigma_n' < 0, "
                 f"minimum {r3.min_normal:.1f} kPa", loc="left", fontsize=10)
    _style(ax)
    fig.tight_layout()
    return fig


def _marker_area(ax, fig, d) -> float:
    """Square markers sized to the cell: cell width in points, squared."""
    ax.set_xlim(d["x"].min() - d["dx"], d["x"].max() + d["dx"])
    ax.set_ylim(d["y"].min() - d["dy"], d["y"].max() + d["dy"])
    ax.set_aspect("equal")
    fig.canvas.draw()
    bbox = ax.get_window_extent()
    px_per_m = min(bbox.width / (ax.get_xlim()[1] - ax.get_xlim()[0]), bbox.height / (ax.get_ylim()[1] - ax.get_ylim()[0]))
    pts = max(px_per_m * min(d["dx"], d["dy"]) * 72 / fig.dpi, 1.0)
    return pts ** 2


def fig_surface3d(case: CaseResult):
    """3D view of the slip surface (mirrored half-model)."""
    g = surface_grid(case)
    if g is None:
        return _placeholder("Slip surface (3D)", "No elements - nothing to draw")
    fig = plt.figure(figsize=(8, 5.5))
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(g["x"], g["y"], g["z"], cmap="Blues_r", edgecolor="#cfcdc7", linewidth=0.2, antialiased=True, alpha=0.95)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.view_init(elev=28, azim=-60)
    ax.set_title("Slip surface (spherical base of the sliding mass)", loc="left", fontsize=11)
    fig.tight_layout()
    return fig


def fig_convergence(case: CaseResult):
    """F against iteration. Drawn only for ACCEPTED dimensions; a withheld F is never plotted."""
    data = convergence_data(case)
    if all(v is None for v in data.values()):
        return _placeholder("F against iteration", "FoS withheld - no ACCEPTED dimension")
    fig, ax = plt.subplots(figsize=(8, 4))
    for tag, col in (("3D", BLUE), ("2D", ORANGE)):
        if data[tag] is not None:
            it, f = data[tag]
            ax.plot(it, f, color=col, lw=2, marker="o", ms=3.5)
            ax.annotate(f"{tag}: F = {f[-1]:.4f}", (it[-1], f[-1]), (6, 6 if tag == "3D" else -14),
                        textcoords="offset points", color=INK, fontsize=9)
    for tag in ("3D", "2D"):
        if data[tag] is None:
            ax.text(0.98, 0.06 if tag == "2D" else 0.12, f"{tag}: FoS withheld", transform=ax.transAxes, ha="right", color="#52514e")
    ax.set_xlabel("iteration")
    ax.set_ylabel("factor of safety F")
    ax.set_title("F against iteration (start F = 1)", loc="left", fontsize=11)
    _style(ax)
    fig.tight_layout()
    return fig


def fig_strength(case: CaseResult):
    """phi'_i and c'_i along the central section."""
    d = strength_section(case)
    if not d:
        return _placeholder("phi'_i and c'_i along the central section",
                            "Final strengths not evaluated (no ACCEPTED or DIAGNOSTIC final state)")
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for tag, col, mk in (("2D", ORANGE, "o"), ("3D", BLUE, "s")):
        if tag in d:
            x, phi, c = d[tag]
            o = np.argsort(x)
            lab = "2D slices (y = 0)" if tag == "2D" else "3D columns nearest y = 0"
            a1.plot(x[o], phi[o], color=col, marker=mk, ms=4, lw=1.6, label=lab)
            a2.plot(x[o], c[o], color=col, marker=mk, ms=4, lw=1.6, label=lab)
    a1.set_ylabel("phi'_i [deg]")
    a2.set_ylabel("c'_i [kPa]")
    a2.set_xlabel("x [m]")
    a1.set_title("Instantaneous friction angle and cohesion along the central section", loc="left", fontsize=11)
    a1.legend(frameon=False, fontsize=8)
    for a in (a1, a2):
        _style(a)
    fig.tight_layout()
    return fig


def fig_mesh(rows: Optional[list[dict]]):
    """Mesh convergence of F (ACCEPTED points only). ``rows`` from ``plotdata.mesh_convergence``."""
    if not rows:
        return _placeholder("Mesh convergence", "Not run (use --mesh-study or the 'Mesh convergence' button)")
    fig, (a3, a2) = plt.subplots(1, 2, figsize=(10, 4))
    for ax, tag, col, label in ((a3, "3D", BLUE, "columns per direction n (n x n half-model)"), (a2, "2D", ORANGE, "slices")):
        pts = [(r["n"], r["F"]) for r in rows if r["dim"] == tag and r["F"] is not None]
        if pts:
            n, f = zip(*pts)
            ax.plot(n, f, color=col, marker="o", lw=2, ms=4)
            ax.annotate(f"{f[-1]:.4f}", (n[-1], f[-1]), (4, 6), textcoords="offset points", fontsize=9)
        else:
            ax.text(0.5, 0.5, "FoS withheld", transform=ax.transAxes, ha="center", color="#52514e")
        ax.set_xscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("F")
        ax.set_title(f"{tag} mesh convergence", loc="left", fontsize=11)
        _style(ax)
    fig.tight_layout()
    return fig


FIGURES = (("section", "Central cross-section", fig_section), ("plan", "Final base normal stress (plan)", fig_plan),
           ("surface3d", "Slip surface (3D)", fig_surface3d), ("convergence", "F against iteration", fig_convergence),
           ("strength", "phi' and c' along the central section", fig_strength))


def all_figures(case: CaseResult, mesh_rows: Optional[list[dict]] = None) -> list[tuple[str, str, "plt.Figure"]]:
    figs = [(key, title, fn(case)) for key, title, fn in FIGURES]
    figs.append(("mesh", "Mesh convergence", fig_mesh(mesh_rows)))
    return figs


# --------------------------------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------------------------------
def _summary_lines(case: CaseResult) -> list[str]:
    df = summary_table(case)
    return [f"{k}: {v}" for k, v in zip(df["item"], df["value"])]


def _png(fig) -> bytes:
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    return buf.getvalue()


def write_report(case: CaseResult, outdir: str | Path, stem: str = "slope3d", pdf: bool = True, html: bool = True,
                 mesh_rows: Optional[list[dict]] = None) -> list[Path]:
    """Write ``<stem>_report.pdf`` and/or ``<stem>_report.html`` into ``outdir``. Returns the paths."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    figs = all_figures(case, mesh_rows)
    written: list[Path] = []
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    if pdf:
        p = out / f"{stem}_report.pdf"
        with PdfPages(p) as pp:
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.07, 0.95, "slope3d run report", fontsize=17, weight="bold")
            fig.text(0.07, 0.925, f"{case.label}   |   slope3d {__version__}   |   {stamp}", fontsize=9, color="#52514e")
            y = 0.89
            for tag, run in (("3D", case.run3d), ("2D", case.run2d)):
                r = run.result
                fig.text(0.07, y, f"{tag}: {r.status.value}", fontsize=15, weight="bold", color=STATUS_COLOUR[r.status])
                fig.text(0.30, y, r.f_display, fontsize=15 if r.diagnostic_f is None else 9.5, color=INK if r.diagnostic_f is None else "#52514e")
                fig.text(0.07, y - 0.022, r.reason, fontsize=8.5, color="#52514e", wrap=True)
                y -= 0.06
            y -= 0.01
            for line in _summary_lines(case):
                fig.text(0.07, y, line, fontsize=8.5)
                y -= 0.0175
            y -= 0.015
            fig.text(0.07, y, "Inputs", fontsize=11, weight="bold")
            y -= 0.02
            tab = inputs_table(case.config, case.mesh3d, case.mesh2d)
            half = (len(tab) + 1) // 2
            for i, (k, v) in enumerate(zip(tab["input"], tab["value"])):
                col, row = divmod(i, half)
                fig.text(0.07 + 0.45 * col, y - row * 0.0165, f"{k} = {v:.10g}" if isinstance(v, float) else f"{k} = {v}", fontsize=7.5)
            pp.savefig(fig)
            plt.close(fig)
            for _, _, f in figs:
                pp.savefig(f, bbox_inches="tight")
        written.append(p)
    if html:
        p = out / f"{stem}_report.html"
        parts = [f"<!doctype html><meta charset='utf-8'><title>slope3d report</title>"
                 "<style>body{font:14px system-ui,sans-serif;max-width:960px;margin:2rem auto;padding:0 1rem;color:#0b0b0b}"
                 "h1{margin-bottom:0}.s{display:inline-block;padding:.2rem .7rem;border-radius:4px;color:#fff;font-weight:700}"
                 "table{border-collapse:collapse}td{border-bottom:1px solid #e6e5e1;padding:.2rem .8rem}img{max-width:100%}</style>",
                 f"<h1>slope3d run report</h1><p>{html_escape(case.label)} &middot; slope3d {__version__} &middot; {stamp}</p>"]
        for tag, run in (("3D", case.run3d), ("2D", case.run2d)):
            r = run.result
            val = html_escape(r.f_display)
            parts.append(f"<h2>{tag}: <span class='s' style='background:{STATUS_COLOUR[r.status]}'>{r.status.value}</span> {val}</h2>"
                         f"<p>{html_escape(r.reason)}</p>")
        parts.append("<h2>Summary</h2><table>" + "".join(f"<tr><td>{html_escape(str(k))}</td><td>{html_escape(str(v))}</td></tr>"
                     for k, v in zip(summary_table(case)["item"], summary_table(case)["value"])) + "</table>")
        parts.append("<h2>Inputs</h2><table>" + "".join(f"<tr><td>{html_escape(str(k))}</td><td>{html_escape(str(v))}</td></tr>"
                     for k, v in zip(*[inputs_table(case.config, case.mesh3d, case.mesh2d)[c] for c in ("input", "value")])) + "</table>")
        for _, title, f in figs:
            b64 = base64.b64encode(_png(f)).decode()
            parts.append(f"<h2>{html_escape(title)}</h2><img alt='{html_escape(title)}' src='data:image/png;base64,{b64}'>")
        p.write_text("\n".join(parts), encoding="utf-8")
        written.append(p)
    for _, _, f in figs:
        plt.close(f)
    return written


def html_escape(s: str) -> str:
    return html.escape(s, quote=True)
