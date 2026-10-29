"""The one-page report of an analysis: the FoS, the cross-section, the analysis notes, the inputs in use and the
studies run on the same inputs, as one self-contained HTML file (the web app's Report, with the same sections). The
cross-section is drawn here as inline SVG, so the file needs nothing else to open or print. No Streamlit here.
"""
from __future__ import annotations

import html
import math
from datetime import date
from typing import Optional

import numpy as np

import engine as E

REPORT_CSS = """
.rpt { font: 14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; color: #111; background: #fff; max-width: 920px; margin: 0 auto; padding: 24px 16px 40px; }
.rpt .eyebrow { margin: 0; font-size: 12px; letter-spacing: 0.06em; text-transform: uppercase; color: #555; }
.rpt h1 { font-size: 22px; margin: 2px 0; }
.rpt h2 { font-size: 15px; margin: 20px 0 6px; padding-bottom: 3px; border-bottom: 1px solid #c8c8c8; }
.rpt h3 { font-size: 13.5px; margin: 12px 0 4px; }
.rpt p { margin: 6px 0; }
.rpt .sub { margin: 0 0 12px; color: #555; }
.rpt .fos { display: flex; flex-wrap: wrap; gap: 10px; margin: 12px 0; }
.rpt .fos > div { border: 1px solid #c8c8c8; border-radius: 6px; padding: 8px 12px; min-width: 160px; }
.rpt .fos span.lab { display: block; font-size: 12.5px; color: #555; }
.rpt .fos b { display: block; font-size: 24px; font-variant-numeric: tabular-nums; }
.rpt figure { margin: 12px 0; }
.rpt svg.fig { width: 100%; height: auto; border: 1px solid #e3e3e3; }
.rpt figcaption { font-size: 12px; color: #555; }
.rpt table { border-collapse: collapse; width: 100%; font-size: 12.5px; margin: 4px 0; }
.rpt th, .rpt td { text-align: left; padding: 3px 6px; border-bottom: 1px solid #e3e3e3; vertical-align: top; }
.rpt thead th { border-bottom: 1px solid #999; }
.rpt td.n { text-align: right; font-variant-numeric: tabular-nums; }
.rpt dl.facts { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 0 18px; font-size: 12.5px; margin: 8px 0; }
.rpt dl.facts > div { display: flex; justify-content: space-between; gap: 10px; border-bottom: 1px solid #eee; padding: 2px 0; }
.rpt dt { color: #555; }
.rpt dd { margin: 0; text-align: right; font-variant-numeric: tabular-nums; }
.rpt .note { display: flex; gap: 8px; border-left: 3px solid #999; padding: 4px 8px; margin: 6px 0; font-size: 12.5px; background: #f6f6f6; }
.rpt .note.critical { border-color: #d03b3b; }
.rpt .note.warning { border-color: #c98a00; }
.rpt .note.info { border-color: #1c5cab; }
.rpt .note > i, .rpt .chip i { font-style: normal; font-weight: 700; }
.rpt .chip { display: inline-flex; gap: 4px; font-size: 12px; }
.rpt .chip.good i { color: #0a7d0a; }
.rpt .chip.critical i { color: #d03b3b; }
.rpt .chip.warning i { color: #a87300; }
.rpt .hint { color: #555; font-size: 12px; }
.rpt footer { margin-top: 22px; padding-top: 8px; border-top: 1px solid #c8c8c8; font-size: 11.5px; color: #555; }
@media (max-width: 640px) { .rpt dl.facts { grid-template-columns: minmax(0, 1fr); } }
@media print { .rpt { padding: 0; max-width: none; } .rpt h2, .rpt h3 { break-after: avoid; } .rpt figure, .rpt tr, .rpt .fos { break-inside: avoid; } }
"""
STATUS = {"ACCEPTED": ("good", "✓", "Accepted"), "REJECTED": ("critical", "✕", "Rejected"),
          "NOT CONVERGED": ("warning", "!", "Not converged"), "NOT APPLICABLE": ("neutral", "–", "Not applicable"),
          "INPUT ERROR": ("critical", "✕", "Input error")}
esc = html.escape


def chip(status: str) -> str:
    c, i, t = STATUS.get(status, ("neutral", "–", status))
    return f'<span class="chip {c}"><i aria-hidden="true">{i}</i>{esc(t)}</span>'


def note(kind: str, text_html: str) -> str:
    icon = {"critical": "✕", "warning": "!", "info": "i", "good": "✓"}[kind]
    return f'<div class="note {kind}"><i aria-hidden="true">{icon}</i><div>{text_html}</div></div>'


def g(v: float, p: int = 6) -> str:
    return f"{v:.{p}g}"


def gc(v: float, p: int = 4) -> str:
    return f"{v:,.{p}g}"


def big(v: float) -> str:
    """Whole kN or kN·m with thousands separators; smaller values to 6 significant digits."""
    return f"{round(v):,}" if abs(v) >= 1e5 else f"{v:,.6g}"


def _tensile(r: E.LemResult) -> str:
    if r is None or r.fos is None or r.sn is None:
        return "–"
    n, total = int((r.sn < 0).sum()), len(r.sn)
    share = "0%" if n == 0 else "<1%" if 100 * n / total < 1 else f"{round(100 * n / total)}%"
    return f"{n} of {total} ({share})"


# ---------------------------------------------------------------------------------------------------------------
# what the results summary says (the web app's caption and facts)
# ---------------------------------------------------------------------------------------------------------------
def caption_text(A: E.Analysis) -> str:
    """Strength, water, iteration path, grid, projection limit and additions of the run, as one paragraph."""
    inp, r3, r2 = A.inputs, A.r3, A.r2
    mb, s, a = inp.hb
    if inp.model == "mc":
        strength = f"c′ = {g(inp.c_mc)} kPa, φ′ = {g(inp.phi_mc)}°"
    else:
        strength = (f"m_b = {mb:.6f}, s = {s:.7f}, a = {a:.6f}"
                    + (f", ν = {g(inp.nu)}" if inp.model in ("ejc_moscarda", "ejc_original") else ""))
    water = (f"u = r_u γ h with r_u = {g(inp.ru)}" if inp.ru is not None else
             ("piezometric line" if inp.piezo is not None else f"water table (H_w = {g(inp.Hw)} m, β_w = {g(inp.beta_w)}°)")
             if inp.wet else "dry")
    cap = (f"Strength: {E.MODELS[inp.model]}. {strength}. Water: {water}. "
           f"Iteration path: {E.ITERATIONS[inp.iteration]}.")
    if inp.base_circle:
        cap += (f" Base-circle grid (exit to entry, floor columns allowed): {int((r3.mesh.x > 0).sum())} floor columns, "
                f"{int((r2.mesh.x > 0).sum())} floor slices.")
    if r3.projection_limits is not None:
        cap += (f" Projection limits {r3.projection_limits[0]:.4f} kPa within T of the face (D = {g(inp.D, 4)}), "
                f"{r3.projection_limits[1]:.4f} kPa beyond (D = 0).")
    elif r3.projection_limit is not None:
        cap += f" Projection limit {r3.projection_limit:.4f} kPa."
    adds = []
    if inp.profile is not None:
        adds.append(f"slope face through {len(inp.profile)} points, crest edge at ({inp.profile[-1][0]:.2f}, "
                    f"{inp.H:.2f}) m")
    if inp.crest_slope:
        adds.append(f"ground behind the crest rising at {g(inp.crest_slope, 4)}°")
    if inp.plan_radius is not None:
        adds.append(f"{'concave' if inp.plan_radius > 0 else 'convex'} wall in plan, R_p = {g(abs(inp.plan_radius))} m")
    if inp.elong != 1 or inp.mid_length > 0:
        adds.append(f"slip surface: ellipsoidal caps λ = {g(inp.elong, 4)}"
                    + (f", cylindrical middle L = {g(inp.mid_length)} m" if inp.mid_length > 0 else ""))
    if inp.piezo is not None:
        adds.append(f"piezometric line through {len(inp.piezo)} points")
    if inp.crack_depth is not None:
        adds.append(f"tension crack {g(inp.crack_depth, 4)} m deep, {g(100 * inp.crack_water, 3)}% full of water")
    if inp.kh > 0:
        adds.append(f"seismic k_h = {g(inp.kh, 4)}")
    if adds:
        cap += f" Additions: {'; '.join(adds)}."
    return cap


def addition_facts(A: E.Analysis) -> list[tuple[str, str]]:
    """The disturbed-zone, tension-crack and seismic rows of the summary (plain text), when they apply."""
    inp, r3, r2 = A.inputs, A.r3, A.r2
    out = []
    if inp.d_zone is not None and inp.model != "mc":
        def zone(r):
            return f"{int(np.sum(r.mesh.zone))} of {r.mesh.count}" if r.mesh.zone is not None else "–"
        out.append((f"Bases within T = {g(inp.d_zone, 4)} m of the face (D applies) 3D / 2D", f"{zone(r3)} / {zone(r2)}"))
    ck = r2.mesh.crack or r3.mesh.crack
    if ck and ck["rows"]:
        rw = ck["rows"][0]
        out.append(("Tension crack", f"x = {ck['x']:.2f} m, {rw['ztop'] - rw['zbot']:.2f} m deep, "
                    + (f"water {rw['hw']:.2f} m, P = {big(r2.mesh.crack['P'])} kN/m (2D)"
                       if r2.mesh.crack and r2.mesh.crack["P"] > 0 else "dry")))
    if inp.kh > 0:
        ex = r3.mesh.ext["Ex"] if r3.mesh.ext is not None else 0.0
        out.append(("Seismic k_h", f"{g(inp.kh)}: k_h ΣW = {big(2 * ex)} kN (3D)"))
    return out


def facts(A: E.Analysis) -> list[tuple[str, str]]:
    r3, r2 = A.r3, A.r2
    out = [("Iterations 3D / 2D", f"{r3.iterations} / {r2.iterations}"),
           ("Active columns / slices", f"{r3.mesh.count} / {r2.mesh.count}"),
           ("Projections 3D / 2D", f"{r3.projections} / {r2.projections}"),
           ("Tensile bases (σ′_n < 0) 3D / 2D", f"{_tensile(r3)} / {_tensile(r2)}")]
    if r3.mesh.count:
        m = r3.mesh
        out += [("3D columns, Δx × Δy", f"{m.n} × {m.ny}, {m.dx:.2f} × {m.dy:.2f} m"),
                ("Column height h", f"{m.h.min():.2f} to {m.h.max():.2f} m")]
    if r2.mesh.count:
        out.append(("2D slice width", f"{r2.mesh.dx:.2f} m"))
    out += [("Sliding volume (full mass)", f"{gc(r3.mesh.volume)} m³"), ("2D area", f"{gc(r2.mesh.volume)} m²")]
    return out + addition_facts(A)


def input_rows(inp: E.Inputs) -> list[tuple[str, str, str]]:
    """The inputs of a run as (label, value, unit): those in use only (as inputRows in the web app)."""
    def pts(q):
        return " ".join(f"({g(u[0])}, {g(u[1])})" for u in q)
    rows = [("Slope height H", g(inp.H), "m")]
    rows.append(("Face angle β", g(inp.beta), "°") if inp.profile is None else
                ("Slope face points (x, z), from the toe", pts(inp.profile), "m"))
    if inp.crest_slope:
        rows.append(("Ground behind the crest rises at α", g(inp.crest_slope), "°"))
    if inp.plan_radius is not None:
        rows.append((f"Wall in plan, {'concave' if inp.plan_radius > 0 else 'convex'}: radius R_p",
                     g(abs(inp.plan_radius)), "m"))
    rows += [("Slip circle centre (x_c, z_c)", f"({g(inp.xc)}, {g(inp.zc)})", "m"), ("Radius R", g(inp.R), "m")]
    if inp.elong != 1 or inp.mid_length:
        rows.append(("3D end caps λ; cylindrical middle L", f"{g(inp.elong)}; {g(inp.mid_length)} m", ""))
    rows += [("Unit weight γ", g(inp.gamma), "kN/m³"), ("Strength", E.MODELS.get(inp.model, inp.model), "")]
    if inp.model != "mc":
        rows += [("σci", g(inp.sigma_ci), "MPa"), ("m_i", g(inp.mi), ""), ("GSI", g(inp.gsi), ""), ("D", g(inp.D), "")]
        if inp.d_zone is not None:
            rows.append(("D applies within T of the face only: T", g(inp.d_zone), "m"))
        if inp.a_override is not None:
            rows.append(("a (set)", g(inp.a_override), ""))
        if inp.model.startswith("ejc"):
            rows.append(("ν", g(inp.nu), ""))
    else:
        rows += [("c′", g(inp.c_mc), "kPa"), ("φ′", g(inp.phi_mc), "°")]
    rows.append(("Water", f"r_u = {g(inp.ru)}" if inp.ru is not None else "dry" if not inp.wet else
                 f"piezometric line {pts(inp.piezo)} m" if inp.piezo is not None else
                 f"water table H_w = {g(inp.Hw)} m, β_w = {g(inp.beta_w)}°", ""))
    if inp.crack_depth is not None:
        rows.append(("Tension crack depth; water fill", f"{g(inp.crack_depth)} m; {g(100 * inp.crack_water)} %", ""))
    if inp.kh:
        rows.append(("Seismic coefficient k_h", g(inp.kh), ""))
    rows += [("3D columns N; N_y; 2D slices", f"{inp.n3}; {'= N' if inp.n3y is None else inp.n3y}; {inp.n2}", ""),
             ("3D equation; tension option",
              f"{'F_M (moments)' if inp.form == 'FM' else 'F_H (horizontal forces)'}; {inp.tension}", ""),
             ("Iteration", E.ITERATIONS.get(inp.iteration, inp.iteration), "")]
    if inp.base_circle:
        rows.append(("Grid", "base circle: exit to entry, floor columns allowed", ""))
    return rows


# ---------------------------------------------------------------------------------------------------------------
# the cross-section y = 0 as SVG: ground, water, slices, slip surface, crack, centre
# ---------------------------------------------------------------------------------------------------------------
def _step(span: float) -> float:
    raw = span / 6.0
    p = 10.0 ** math.floor(math.log10(raw))
    return next(k * p for k in (1, 2, 5, 10) if raw <= k * p)


def section_svg(inp: E.Inputs, r2: E.LemResult, width: int = 1000, max_height: int = 640) -> str:
    m = r2.mesh
    crest = -inp.H / math.tan(math.radians(inp.beta))
    if m.count:
        left, right = min(m.x1, crest), max(m.x0, 0.0)
    else:
        x0 = inp.xc - math.sqrt(inp.R * inp.R - inp.zc * inp.zc) if inp.R > inp.zc else inp.xc
        left, right = min(inp.xc - inp.R, crest), max(x0, 0.0)
    pad = 0.1 * max(right - left, inp.H)
    xa, xb = left - pad, right + 2 * pad
    span = xb - xa
    x0, x1 = min(xa, inp.xc - 0.05 * span), max(xb, inp.xc + 0.12 * span)
    z0, z1 = min(-pad / 2, inp.zc - inp.R - pad / 2), max(inp.H, inp.zc) + pad
    L, R_, T, B = 58, 16, 14, 46
    s = (width - L - R_) / (x1 - x0)
    if T + B + s * (z1 - z0) > max_height:
        s = (max_height - T - B) / (z1 - z0)
    pw, ph = s * (x1 - x0), s * (z1 - z0)
    height = int(math.ceil(T + B + ph))
    ox = L + ((width - L - R_) - pw) / 2

    def X(x):
        return ox + (x - x0) * s

    def Z(z):
        return T + (z1 - z) * s

    def poly(xs, zs):
        return " ".join(f"{X(x):.1f},{Z(z):.1f}" for x, z in zip(xs, zs))

    def lines(xs, zs):                            # polylines split where z is NaN
        out, cur = [], []
        for x, z in zip(xs, zs):
            if np.isfinite(z):
                cur.append((x, z))
            elif cur:
                out.append(cur)
                cur = []
        return out + ([cur] if cur else [])

    parts = [f'<svg class="fig" viewBox="0 0 {width} {height}" role="img" xmlns="http://www.w3.org/2000/svg" '
             f'aria-label="Cross-section on the central plane y = 0: the slope, the slip surface and the 2D slices.">',
             f'<defs><clipPath id="rpt-plot"><rect x="{ox:.1f}" y="{T}" width="{pw:.1f}" height="{ph:.1f}"/></clipPath></defs>',
             f'<rect x="{ox:.1f}" y="{T}" width="{pw:.1f}" height="{ph:.1f}" fill="#fcfcfb" stroke="#c8c8c8"/>',
             '<g clip-path="url(#rpt-plot)">']
    xs = np.linspace(x0, x1, 600)
    zt = E.ground(inp, xs)
    parts.append(f'<polygon points="{poly(np.r_[xs, x1, x0], np.r_[zt, z0, z0])}" fill="#e7e2d8"/>')
    # grid
    for v in np.arange(math.ceil(x0 / _step(x1 - x0)) * _step(x1 - x0), x1, _step(x1 - x0)):
        parts.append(f'<line x1="{X(v):.1f}" y1="{T}" x2="{X(v):.1f}" y2="{T + ph:.1f}" stroke="#000" stroke-opacity="0.07"/>')
    for v in np.arange(math.ceil(z0 / _step(z1 - z0)) * _step(z1 - z0), z1, _step(z1 - z0)):
        parts.append(f'<line x1="{ox:.1f}" y1="{Z(v):.1f}" x2="{ox + pw:.1f}" y2="{Z(v):.1f}" stroke="#000" stroke-opacity="0.07"/>')
    if m.count and m.count <= 400:
        for k in range(m.count):
            xl, xr = m.x[k] - m.dx / 2, m.x[k] + m.dx / 2
            zl, zr = (float(E.sphere_z(v, 0.0, inp.xc, inp.zc, inp.R)) for v in (xl, xr))
            zl, zr = (v if math.isfinite(v) else float(m.zb[k]) for v in (zl, zr))
            gl, gr = float(E.ground(inp, xl)), float(E.ground(inp, xr))
            parts.append(f'<polygon points="{poly([xl, xr, xr, xl], [zl, zr, gr, gl])}" fill="#2a78d6" fill-opacity="0.13" '
                         'stroke="#3c3c3c" stroke-opacity="0.45" stroke-width="0.6"/>')
    parts.append(f'<polyline points="{poly(xs, zt)}" fill="none" stroke="#8a7f70" stroke-width="2"/>')
    if inp.wet:
        zw = np.minimum(E.water(inp, xs), zt)
        parts.append(f'<polyline points="{poly(xs, zw)}" fill="none" stroke="#1baf7a" stroke-width="2" stroke-dasharray="7 5"/>')
    ends = []
    if m.count:
        xb_ = np.linspace(m.x1, m.x0 if inp.base_circle else min(m.x0, 0.0), 800)
        zb_ = E.sphere_z(xb_, 0.0, inp.xc, inp.zc, inp.R)
        with np.errstate(invalid="ignore"):
            zin = np.where(zb_ < E.ground(inp, xb_), zb_, np.nan)
        for seg in lines(xb_, zin):
            parts.append(f'<polyline points="{poly(*zip(*seg))}" fill="none" stroke="#eb6834" stroke-width="3"/>')
            ends += [seg[0], seg[-1]]
    ck = m.crack
    if ck and ck["rows"]:
        rw = ck["rows"][0]
        parts.append(f'<line x1="{X(ck["x"]):.1f}" y1="{Z(rw["zbot"]):.1f}" x2="{X(ck["x"]):.1f}" y2="{Z(rw["ztop"]):.1f}" '
                     'stroke="#222" stroke-width="2.5"/>')
        if rw["hw"] > 0:
            parts.append(f'<line x1="{X(ck["x"]):.1f}" y1="{Z(rw["zbot"]):.1f}" x2="{X(ck["x"]):.1f}" '
                         f'y2="{Z(rw["zbot"] + rw["hw"]):.1f}" stroke="#1baf7a" stroke-width="6"/>')
    if ends:
        for x, z in (ends[0], ends[-1]):
            parts.append(f'<line x1="{X(inp.xc):.1f}" y1="{Z(inp.zc):.1f}" x2="{X(x):.1f}" y2="{Z(z):.1f}" stroke="#eb6834" '
                         'stroke-width="0.8" stroke-dasharray="5 4"/>')
    cx, cz = X(inp.xc), Z(inp.zc)
    parts.append(f'<path d="M{cx - 6:.1f},{cz - 6:.1f}L{cx + 6:.1f},{cz + 6:.1f}M{cx - 6:.1f},{cz + 6:.1f}L{cx + 6:.1f},{cz - 6:.1f}" '
                 f'stroke="#eb6834" stroke-width="2.5"/><text x="{cx + 10:.1f}" y="{cz + 4:.1f}" font-size="12" fill="#111">centre</text>')
    parts.append("</g>")
    # axes
    st = _step(x1 - x0)
    for v in np.arange(math.ceil(x0 / st) * st, x1 + 1e-9, st):
        parts.append(f'<text x="{X(v):.1f}" y="{T + ph + 16:.1f}" font-size="11" text-anchor="middle" fill="#444">{v:g}</text>')
    st = _step(z1 - z0)
    for v in np.arange(math.ceil(z0 / st) * st, z1 + 1e-9, st):
        parts.append(f'<text x="{ox - 6:.1f}" y="{Z(v) + 4:.1f}" font-size="11" text-anchor="end" fill="#444">{v:g}</text>')
    parts.append(f'<text x="{ox + pw / 2:.1f}" y="{height - 8}" font-size="12" text-anchor="middle" fill="#444">'
                 'x (m), toe at 0, slope rises to −x</text>')
    parts.append(f'<text x="14" y="{T + ph / 2:.1f}" font-size="12" text-anchor="middle" fill="#444" '
                 f'transform="rotate(-90 14 {T + ph / 2:.1f})">z (m)</text>')
    legend = [("#8a7f70", "ground", ""), ("#eb6834", "slip surface", "")]
    if inp.wet:
        legend.append(("#1baf7a", "piezometric line" if inp.piezo is not None else "water table", "7 5"))
    if ck and ck["rows"]:
        legend.append(("#222", "tension crack", ""))
    lx = ox + 10
    for col, lab, dash in legend:
        parts.append(f'<line x1="{lx:.1f}" y1="{T + 14}" x2="{lx + 22:.1f}" y2="{T + 14}" stroke="{col}" stroke-width="3"'
                     + (f' stroke-dasharray="{dash}"' if dash else "") + "/>"
                     f'<text x="{lx + 28:.1f}" y="{T + 18}" font-size="12" fill="#111">{esc(lab)}</text>')
        lx += 46 + 7 * len(lab)
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------------------------------------------
# studies as tables
# ---------------------------------------------------------------------------------------------------------------
MODEL_LAB = {"ejc_moscarda": "EJC", "ghb": "GHB", "mc": "MC"}


def _f(v, d=6):
    return "withheld" if v is None else f"{v:.{d}f}"


def compare_table(res: dict, nu: float) -> str:
    e, gm = res["models"]["ejc"], res["models"]["ghb"]
    r4 = (lambda v: "–" if v is None else f"{v:.4f}")
    return ('<table><thead><tr><th></th><th>FoS 3D</th><th>FoS 2D</th><th>3D / 2D</th></tr></thead><tbody>'
            f'<tr><td>EJC, Moscarda (ν = {g(nu, 4)})</td><td class="n">{_f(e["fos3"])}</td><td class="n">{_f(e["fos2"])}</td>'
            f'<td class="n">{r4(e["ratio"])}</td></tr>'
            f'<tr><td>GHB</td><td class="n">{_f(gm["fos3"])}</td><td class="n">{_f(gm["fos2"])}</td><td class="n">{r4(gm["ratio"])}</td></tr>'
            f'<tr><td><b>EJC / GHB</b></td><td class="n"><b>{r4(res["ratio"]["ejc_ghb3"])}</b></td>'
            f'<td class="n"><b>{r4(res["ratio"]["ejc_ghb2"])}</b></td><td></td></tr></tbody></table>')


def back_text(res: dict, label: str) -> str:
    if res["value"] is None:
        return note("warning", esc(res["message"][:1].upper() + res["message"][1:]) + ".")
    return (f'<p><b>{esc(label)} = {g(res["value"])}</b> gives F = {_f(res["fos"])} ({res["dim"]}D), for the target '
            f'{g(res["target"])}; bracket {g(res["bracket"][0], 8)} to {g(res["bracket"][1], 8)}, {len(res["steps"])} analyses.</p>'
            + (note("warning", esc(res["message"])) if res["message"] else ""))


def sweep_table(rows: list, label: str, models, dims) -> str:
    series = [(m, d) for m in models for d in dims]
    head = "".join(f"<th>{MODEL_LAB.get(m, m)} {d}D</th>" for m, d in series)
    body = "".join(f'<tr><td class="n">{g(r["value"])}</td>'
                   + "".join(f'<td class="n">{_f(r[m][f"fos{d}"])}</td>' for m, d in series) + "</tr>" for r in rows)
    return f"<table><thead><tr><th>{esc(label)}</th>{head}</tr></thead><tbody>{body}</tbody></table>"


def mc_table(res: dict, labels: Optional[list] = None) -> str:
    f4 = (lambda v: "–" if v is None else f"{v:.4f}")
    rows = []
    for m in res["models"]:
        for d in res["dims"]:
            st = res["stats"][f"{m}{d}"]
            if st["pof"] is None:
                pf = "–"
            elif st["pof"] > 0:
                pf = (f'{100 * st["pof"]:.1f} % ± {100 * st["pof_se"]:.1f} % '
                      f'({round(st["pof"] * st["accepted"])} failed)')
            else:
                pf = f'0 of {st["accepted"]} failed: below {300 / st["accepted"]:.2g} % (95 %)'
            rows.append(f'<tr><td>{MODEL_LAB.get(m, m)} {d}D</td><td class="n">{st["accepted"]} of {st["n"]}</td>'
                        f'<td class="n">{f4(st["mean"])}</td><td class="n">{f4(st["sd"])}</td><td class="n">{f4(st["beta"])}</td>'
                        f'<td class="n">{pf}</td><td class="n">{f4(st["p05"])} / {f4(st["p50"])} / {f4(st["p95"])}</td></tr>')
    varied = ", ".join(f"{esc((labels or [])[q] if labels and q < len(labels) else v[0])} (mean {g(v[1])}, SD {g(v[2])})"
                       for q, v in enumerate(res["variables"]))
    return ('<table><thead><tr><th></th><th>Accepted</th><th>Mean F</th><th>SD</th><th>β = (mean − 1)/SD</th>'
            '<th>P<sub>f</sub> = P(F &lt; 1)</th><th>5 % / 50 % / 95 %</th></tr></thead><tbody>' + "".join(rows)
            + f'</tbody></table><p class="hint">{len(res["samples"])} samples, seed {res["seed"]}; varied: {varied}. '
            "The ± is one standard error of P<sub>f</sub>, √(P<sub>f</sub>(1 − P<sub>f</sub>)/n). When no sample fails, the "
            "95 % upper bound on P<sub>f</sub> is 3/n (the rule of three).</p>")


# ---------------------------------------------------------------------------------------------------------------
# the report
# ---------------------------------------------------------------------------------------------------------------
def report_html(A: E.Analysis, *, title: str, subtitle: str = "", notes: Optional[list] = None,
                studies: Optional[list] = None, day: Optional[str] = None) -> str:
    """The report as a complete HTML document. ``notes`` = [(kind, text)], ``studies`` = [(heading, html)] run on
    A.inputs."""
    inp, r3, r2 = A.inputs, A.r3, A.r2
    day = day or date.today().isoformat()

    def box(lab, r):
        return (f'<div><span class="lab">{lab}</span><b>{_f(r.fos) if r is not None else "withheld"}</b>'
                f'{chip(r.status if r is not None else "INPUT ERROR")}</div>')
    ratio = f"{A.ratio:.4f}" if A.ratio is not None else "withheld"
    body = [f'<header><p class="eyebrow">EJC Slope Stability · case report · {esc(day)}</p><h1>{esc(title)}</h1>'
            + (f'<p class="sub">{esc(subtitle)}</p>' if subtitle else "") + "</header>"]
    if inp.model == "ejc_original":
        body.append("<h2>Result</h2><p>Original closure: an admissibility check only; it gives no FoS.</p>")
    else:
        body.append(f'<div class="fos">{box("FoS 3D, Hungr " + ("F<sub>M</sub>" if inp.form == "FM" else "F<sub>H</sub>"), r3)}'
                    f'{box("FoS 2D, Bishop", r2)}<div><span class="lab">3D / 2D</span><b>{ratio}</b></div></div>')
        body.append(f"<figure>{section_svg(inp, r2)}<figcaption>Cross-section on the central plane (y = 0) with the 2D "
                    "slices.</figcaption></figure>")
    body.append(f"<h2>Analysis</h2><p>{esc(caption_text(A))}</p><dl class=\"facts\">"
                + "".join(f"<div><dt>{esc(k)}</dt><dd>{esc(v)}</dd></div>" for k, v in facts(A)) + "</dl>")
    for lab, r in (("3D", r3), ("2D", r2)):
        if r.status != "ACCEPTED" and inp.model != "ejc_original":
            body.append(note("warning" if r.status == "NOT CONVERGED" else "critical",
                             f"<b>{lab}: {esc(r.status)}, no accepted FoS.</b> {esc(r.reason)}"))
    for kind, text in notes or []:
        body.append(note(kind, esc(text)))
    body.append("<h2>Inputs</h2><table><tbody>" + "".join(
        f'<tr><th scope="row">{esc(lab)}</th><td>{esc(v)}</td><td>{esc(u)}</td></tr>' for lab, v, u in input_rows(inp))
        + "</tbody></table>")
    if studies:
        body.append("<h2>Studies on these inputs</h2>" + "".join(f"<h3>{esc(h)}</h3>{t}" for h, t in studies))
    body.append("<footer><p>EJC Slope Stability: fully 3D limit equilibrium with Hungr's columns and Bishop's 2D slices on "
                "the central section, for a prescribed slip surface; strength from the Extended Jiang Criterion solved by "
                "Route 2, with GHB and Mohr-Coulomb for comparison. Research tool for the prescribed-surface comparison in "
                "Mahat (2026), UWA GENG5512 thesis; not for design.</p><p>Solved by ejc_app/engine.py. Compression "
                "positive, effective stresses.</p></footer>")
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f"<title>{esc(title)}: EJC Slope Stability report</title><style>body{{margin:0;background:#fff}}{REPORT_CSS}"
            f'</style></head><body><main class="rpt">{"".join(body)}</main></body></html>')
