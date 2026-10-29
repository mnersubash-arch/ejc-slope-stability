"""Shared Streamlit helpers: input state, the sidebar (number box + slider kept in sync), presets and running.

Input values live in ``st.session_state["inp"]`` (a plain dict, kept across pages). Each widget has its own key,
initialised from that dict when it is first drawn on a page and copied back to it in an ``on_change`` callback.
"""
from __future__ import annotations

import math
import re
from dataclasses import asdict
from typing import Optional

import pandas as pd
import streamlit as st

import casefile as CF
import engine as E
import presets as PR

LIVE_MAX_N = 40          # live update only for N <= 40 (brief, page 1); N x N_y <= 40 x 40 with N_y
DISCLAIMER = "Research tool for the prescribed-surface comparison in Mahat (2026) thesis; not for design."

# name: (label, unit, min, max, step, format)
NUM = {
    "H": ("H, slope height", "m", 0.5, 1000.0, 0.1, "%.2f"),
    "beta": ("β, slope face angle", "deg", 10.0, 89.0, 0.5, "%.3f"),
    "xc": ("x_c, sphere centre", "m", -1000.0, 2000.0, 0.01, "%.7f"),
    "zc": ("z_c, sphere centre", "m", 1.0, 3000.0, 0.01, "%.7f"),
    "R": ("R, sphere radius", "m", 1.0, 4000.0, 0.01, "%.7f"),
    "gamma": ("γ, unit weight", "kN/m³", 10.0, 40.0, 0.1, "%.2f"),
    "sigma_ci": ("σ_ci, intact UCS", "MPa", 1.0, 300.0, 0.1, "%.2f"),
    "mi": ("m_i", "-", 1.0, 40.0, 0.5, "%.2f"),
    "gsi": ("GSI", "-", 10.0, 100.0, 1.0, "%.1f"),
    "D": ("D, disturbance factor", "-", 0.0, 1.0, 0.05, "%.2f"),
    "a_value": ("a (override value)", "-", 0.30, 1.00, 0.001, "%.6f"),
    "nu": ("ν, Poisson's ratio", "-", 0.01, 0.49, 0.001, "%.3f"),
    "c_mc": ("c′, cohesion (Mohr-Coulomb)", "kPa", 0.0, 1000.0, 0.1, "%.2f"),
    "phi_mc": ("φ′, friction angle (Mohr-Coulomb)", "deg", 0.0, 60.0, 0.1, "%.2f"),
    "ru_value": ("r_u, pore-pressure ratio", "-", 0.0, 0.99, 0.01, "%.3f"),
    "Hw": ("H_w, water-table height", "m", 0.0, 1000.0, 1.0, "%.1f"),
    "beta_w": ("β_w, water-table angle", "deg", 1.0, 89.0, 0.5, "%.1f"),
    "n3": ("N, 3D mesh (N x N half-model)", "-", 2, 320, 1, None),
    "n3y_value": ("N_y, 3D columns across y (0 to R)", "-", 2, 320, 1, None),
    "dx_value": ("Δx, column length along x, at most", "m", 0.01, 4000.0, 0.01, "%.4f"),
    "dy_value": ("Δy, column width across y, at most", "m", 0.01, 4000.0, 0.01, "%.4f"),
    "n2": ("N_2D, 2D slices", "-", 2, 10000, 1, None),
    # additions (all off by default)
    "bench_n": ("n, benches", "-", 1, 50, 1, None),
    "bench_h": ("h_b, bench height", "m", 0.5, 500.0, 0.1, "%.2f"),
    "bench_beta": ("β_b, bench face angle", "deg", 10.0, 89.0, 0.1, "%.2f"),
    "berm_w": ("w_b, berm width", "m", 0.0, 500.0, 0.1, "%.2f"),
    "crest_slope": ("α, ground behind the crest edge", "deg", 0.0, 60.0, 0.1, "%.2f"),
    "plan_r": ("R_p, radius of the toe line in plan", "m", 1.0, 100000.0, 1.0, "%.2f"),
    "elong": ("λ, end caps: half-width / R (1 = sphere)", "-", 0.1, 20.0, 0.01, "%.3f"),
    "mid_length": ("L, cylindrical middle part, length", "m", 0.0, 100000.0, 1.0, "%.2f"),
    "d_zone": ("T, disturbed zone, distance from the face", "m", 0.0, 10000.0, 0.5, "%.2f"),
    "crack_depth": ("z_t, tension crack depth below the ground", "m", 0.01, 5000.0, 0.5, "%.2f"),
    "crack_water": ("f_w, share of the crack filled with water", "-", 0.0, 1.0, 0.01, "%.2f"),
    "kh": ("k_h, horizontal seismic coefficient", "-", 0.0, 0.99, 0.005, "%.3f"),
}
# inputs with no default: the box starts empty and the run waits until a value is entered (never invented)
OPTIONAL = ("plan_r", "d_zone", "crack_depth")
STRENGTH = {"EJC": "EJC (Extended Jiang Criterion)", "GHB": "GHB (conventional generalised Hoek-Brown)",
            "MC": "Mohr-Coulomb (constant c′, φ′)"}
WATER = {"Dry": "Dry (u = 0)", "Wet": "Water table", "Piezo": "Piezometric line", "ru": "r_u (u = r_u γ h)"}
FACE = {"planar": "Planar", "benches": "Benches", "profile": "Points"}
WALL = {"straight": "Straight", "concave": "Concave", "convex": "Convex"}
CLOSURE = {"moscarda": "Moscarda: σ2 = ν(σ1 + σ3) (adopted)",
           "original": "Original: σ2 = σ3/ν − σ1 (admissibility check only)",
           "s2s3": "σ2 = σ3 (check: must reproduce GHB)"}
FORM = {"FH": "F_H, force form (governing, eq. 4.14)", "FM": "F_M, moment form (eq. 4.15, own iteration, 3D)"}
TENSION = {"none": "None", "C0": "C0: remove c_i on bases with σ_n < 0", "T0": "T0: remove all basal shear there"}
ITER = {"default": "Default (§1.4): F from F_H/F_M, then σ_n", "native_v15": "Native v15 iteration path"}
COLS = {"square": "N × N", "counts": "N × N_y", "size": "Δx, Δy (m)"}
CHOICES = {"strength": STRENGTH, "closure": CLOSURE, "form": FORM, "tension": TENSION, "iteration": ITER}


# ---------------------------------------------------------------------------------------------------------------
# Inputs <-> state dict
# ---------------------------------------------------------------------------------------------------------------
def to_state(inp: E.Inputs) -> dict:
    d = asdict(inp)
    d["a_on"] = inp.a_override is not None
    d["a_value"] = 0.5 if inp.a_override is None else float(inp.a_override)
    d["strength"] = {"ghb": "GHB", "mc": "MC"}.get(inp.model, "EJC")
    d["closure"] = {"ejc_original": "original", "ejc_sigma2_eq_sigma3": "s2s3"}.get(inp.model, "moscarda")
    d["water"] = "ru" if inp.ru is not None else ("Wet" if inp.wet else "Dry")
    d["ru_value"] = 0.0 if inp.ru is None else float(inp.ru)
    d["cols"] = "square" if inp.n3y is None else "counts"
    d["n3y_value"] = inp.n3 if inp.n3y is None else int(inp.n3y)
    # additions: the slope face, wall in plan, D zone, piezometric line, tension crack (off: empty or the planar values)
    tb, tw = math.tan(math.radians(inp.beta)), math.tan(math.radians(inp.beta_w))
    d["face"] = "planar" if inp.profile is None else "profile"
    d["profile_text"] = points_text(((0.0, 0.0), (-inp.H / tb, inp.H)) if inp.profile is None else inp.profile)
    d["bench_n"], d["bench_h"], d["bench_beta"], d["berm_w"] = 1, inp.H, inp.beta, 0.0
    d["wall"] = "straight" if inp.plan_radius is None else ("concave" if inp.plan_radius > 0 else "convex")
    d["plan_r"] = None if inp.plan_radius is None else abs(float(inp.plan_radius))
    d["dz_on"] = inp.d_zone is not None
    if inp.piezo is not None and inp.wet:
        d["water"] = "Piezo"
    d["piezo_text"] = points_text(((0.0, 0.0), (-inp.Hw / tw, inp.Hw)) if inp.piezo is None else inp.piezo)
    d["crack_on"] = inp.crack_depth is not None
    d["profile"] = d["piezo"] = None       # the points live in the text boxes
    sz = sizes_for(d, inp.n3, d["n3y_value"])
    d["dx_value"], d["dy_value"] = sz if sz else (inp.R / inp.n3, inp.R / inp.n3)
    return d


def parse_points(text: str) -> list[tuple[float, float]]:
    """Points [(x, z), ...] from text: one "x, z" pair per line (or separated by ";"). Raises ValueError naming the
    line. Same as parsePoints in webapp/index.html."""
    pts = []
    lines = [ln.strip() for ln in re.split(r"[\n;]+", str(text or "")) if ln.strip()]
    for q, ln in enumerate(lines):
        v = [t for t in re.split(r"[\s,]+", ln) if t]
        try:
            x, z = (float(t) for t in v) if len(v) == 2 else (math.nan, math.nan)
        except ValueError:
            x = z = math.nan
        if not (math.isfinite(x) and math.isfinite(z)):
            raise ValueError(f'line {q + 1} ("{ln}") is not two numbers x, z')
        pts.append((x, z))
    return pts


def points_text(pts) -> str:
    return "\n".join(f"{float(x):.10g}, {float(z):.10g}" for x, z in pts)


def bench_points(d: dict) -> list[tuple[float, float]]:
    """Bench face from the toe: n benches of height h_b at beta_b, separated by berms w_b wide."""
    n, h, w = int(d["bench_n"]), float(d["bench_h"]), float(d["berm_w"])
    run = h / math.tan(math.radians(float(d["bench_beta"])))
    pts, x, z = [(0.0, 0.0)], 0.0, 0.0
    for k in range(n):
        x -= run
        z += h
        pts.append((x, z))
        if k < n - 1 and w > 0:
            x -= w
            pts.append((x, z))
    return pts


def face_fields(d: dict) -> tuple[float, float, Optional[tuple]]:
    """(H, beta, profile) of the slope-face option: the planar H and beta, or the points (benches or typed) with
    H = the crest edge's z and beta the overall angle from the toe. Same as faceFields in webapp/index.html."""
    try:
        pts = (bench_points(d) if d.get("face") == "benches" else
               parse_points(d.get("profile_text", "")) if d.get("face") == "profile" else None)
    except ValueError:
        pts = None
    if not pts or len(pts) < 2:
        return float(d["H"]), float(d["beta"]), None
    xe, he = pts[-1]
    beta = math.degrees(math.atan(he / -xe)) if xe < 0 and he > 0 else float(d["beta"])
    return float(he), beta, tuple((float(x), float(z)) for x, z in pts)


def _opt(d: dict, name: str):
    v = d.get(name)
    return None if v is None or not math.isfinite(float(v)) else float(v)


# 3D columns: N x N (the brief's grid), N x N_y, or plan sizes in metres converted to counts on the current circle
def _geom(d: dict) -> E.Inputs:
    H, beta, profile = face_fields(d)
    return E.Inputs(H=H, beta=beta, xc=float(d["xc"]), zc=float(d["zc"]), R=float(d["R"]),
                    base_circle=bool(d["base_circle"]), crest_slope=float(d.get("crest_slope", 0.0)),
                    elong=float(d.get("elong", 1.0)), mid_length=float(d.get("mid_length", 0.0)), profile=profile)


def _ceil4(v: float) -> float:
    q = 10.0 ** (math.floor(math.log10(v)) - 3)
    return float(f"{math.ceil(v / q - 1e-9) * q:.4g}")


def _floor4(v: float) -> float:
    q = 10.0 ** (math.floor(math.log10(v)) - 3)
    return float(f"{math.floor(v / q + 1e-9) * q:.4g}")


def span_of(d: dict):
    """Grid span x0 - x1 (m) of the circle of d, or None when the default grid refuses it (the run then says why)."""
    try:
        x0, x1 = E.grid_x_range(_geom(d))
    except ValueError:
        return None
    return x0 - x1


def sizes_for(d: dict, n: int, ny: int):
    """Plan size (dx, dy) in m of n x ny columns on the circle of d, rounded up to 4 significant digits (which gives
    the same counts back), or None when the circle has no grid."""
    span = span_of(d)
    if span is None:
        return None
    return _clamp("dx_value", _ceil4(span / n)), _clamp("dy_value", _ceil4(E.y_extent(_geom(d)) / ny))


def col_counts(d: dict):
    """Counts (N, N_y) of the column option of d; N_y is None for the N x N grid. Sizes in metres become counts with
    engine.counts_for_size on the current circle; if the circle has no grid, the run reports why."""
    cols = d.get("cols", "square")
    if cols == "counts":
        return int(d["n3"]), int(d["n3y_value"])
    if cols == "size":
        try:
            return E.counts_for_size(_geom(d), float(d["dx_value"]), float(d["dy_value"]))
        except ValueError:
            return int(d["n3"]), None
    return int(d["n3"]), None


def col_grid(d: dict) -> tuple[int, int]:
    n, ny = col_counts(d)
    return n, n if ny is None else ny


def mesh_problem(d: dict):
    """With sizes in metres, the counts they give must lie in 2 to 320, as N: the message, or None."""
    span = span_of(d) if d.get("cols") == "size" else None
    if span is None:
        return None
    n, ny = col_counts(d)
    out = []
    for sym, v, cnt, csym, axis, L in (("Δx", d["dx_value"], n, "N", "along x", span),
                                        ("Δy", d["dy_value"], ny, "N_y", "across y", E.y_extent(_geom(d)))):
        if not 2 <= cnt <= 320:
            out.append(f"{sym} = {float(v):g} m gives {csym} = {cnt} column{'' if cnt == 1 else 's'} {axis} on this "
                       f"circle; {csym} must be from 2 to 320, so use {sym} from {_ceil4(L / 320):g} to "
                       f"{_floor4(L / 2):g} m.")
    return " ".join(out) or None


def cols_text(d: dict) -> str:
    span = span_of(d)
    if span is None:
        return ""
    n, ny = col_grid(d)
    ysym = "N" if d.get("cols", "square") == "square" else "N_y"
    Y, R = E.y_extent(_geom(d)), float(d["R"])
    dy = (f"Δy = R/{ysym} = {R / ny:.3f} m" if Y == R else
          f"Δy = Y/{ysym} = {Y / ny:.3f} m (Y = {Y:.2f} m, the half-width across y)")
    return ((f"On this circle: N = {n} along x and N_y = {ny} across y. " if d.get("cols") == "size" else "")
            + f"Columns {n} × {ny}: Δx = (x0 − x1)/N = {span / n:.3f} m by {dy}. "
            "The height h = z_t − z_b of each column is set by the ground and the slip surface.")


def input_problems(d: dict) -> list[str]:
    """What stops a run: unparsable points, a required value not entered, the engine's checks of the additions, and
    column sizes whose counts fall outside 2 to 320. Empty when the inputs can be solved."""
    out = []
    if d.get("face") == "profile":
        try:
            parse_points(d.get("profile_text", ""))
        except ValueError as exc:
            out.append(f"Face points: {exc}.")
    if d.get("water") == "Piezo":
        try:
            parse_points(d.get("piezo_text", ""))
        except ValueError as exc:
            out.append(f"Piezometric line: {exc}.")
    need = [("plan_r", d.get("wall", "straight") != "straight"), ("d_zone", d.get("dz_on") and d["strength"] != "MC"),
            ("crack_depth", d.get("crack_on"))]
    out += [f"Enter {NUM[k][0]} ({NUM[k][1]}): it has no default." for k, on in need if on and _opt(d, k) is None]
    if not out:
        try:
            E.check_additions(to_inputs(d))
        except ValueError as exc:
            out.append(f"{str(exc)[:1].upper()}{str(exc)[1:]}.")
    problem = mesh_problem(d)
    if problem:
        out.append(problem)
    return out


def to_inputs(d: dict) -> E.Inputs:
    if d["strength"] == "GHB":
        model = "ghb"
    elif d["strength"] == "MC":
        model = "mc"
    else:
        model = {"moscarda": "ejc_moscarda", "original": "ejc_original", "s2s3": "ejc_sigma2_eq_sigma3"}[d["closure"]]
    native = d["iteration"] == "native_v15"
    n3, ny = col_counts(d)                 # N_y = N is the N x N grid
    H, beta, profile = face_fields(d)
    wall = {"concave": 1.0, "convex": -1.0}.get(d.get("wall", "straight"), 0.0)
    piezo = None
    if d["water"] == "Piezo":
        try:
            piezo = tuple(parse_points(d.get("piezo_text", "")))
        except ValueError:
            piezo = None
    plan_r, d_zone, crack = _opt(d, "plan_r"), _opt(d, "d_zone"), _opt(d, "crack_depth")
    return E.Inputs(H=H, beta=beta, xc=float(d["xc"]), zc=float(d["zc"]), R=float(d["R"]),
                    gamma=float(d["gamma"]), sigma_ci=float(d["sigma_ci"]), mi=float(d["mi"]), gsi=float(d["gsi"]),
                    D=float(d["D"]), a_override=float(d["a_value"]) if d["a_on"] else None, nu=float(d["nu"]),
                    model=model, wet=d["water"] in ("Wet", "Piezo"), Hw=float(d["Hw"]), beta_w=float(d["beta_w"]),
                    n3=n3, n3y=None if ny is None or ny == n3 else ny, n2=int(d["n2"]),
                    form="FH" if native else d["form"],
                    tension="none" if native else d["tension"], iteration=d["iteration"], c_mc=float(d["c_mc"]),
                    phi_mc=float(d["phi_mc"]), ru=float(d["ru_value"]) if d["water"] == "ru" else None,
                    base_circle=bool(d["base_circle"]), profile=profile, crest_slope=float(d.get("crest_slope", 0.0)),
                    plan_radius=wall * plan_r if wall and plan_r is not None else None,
                    elong=float(d.get("elong", 1.0)), mid_length=float(d.get("mid_length", 0.0)),
                    d_zone=d_zone if d.get("dz_on") and d["strength"] != "MC" else None, piezo=piezo,
                    crack_depth=crack if d.get("crack_on") else None,
                    crack_water=float(d.get("crack_water", 0.0)) if d.get("crack_on") else 0.0, kh=float(d.get("kh", 0.0)))


WIDGET_NAMES = list(NUM) + ["a_on", "strength", "closure", "water", "form", "tension", "iteration", "base_circle",
                             "cols", "face", "wall", "dz_on", "crack_on", "profile_text", "piezo_text"]


def _wkeys(name: str) -> list[str]:
    if name in OPTIONAL:
        return [f"nb_{name}"]
    return [f"nb_{name}", f"sl_{name}"] if name in NUM else [f"w_{name}"]


def _clamp(name, v):
    lo, hi = NUM[name][2], NUM[name][3]
    return min(max(v, lo), hi)


def _write_widgets():
    """Copy the state dict into every widget key (called from callbacks, before widgets are drawn)."""
    d = st.session_state.inp
    for name in WIDGET_NAMES:
        for k in _wkeys(name):
            if name in NUM:
                st.session_state[k] = _clamp(name, d[name]) if k.startswith("sl_") else d[name]
            else:
                st.session_state[k] = d[name]


def init_state():
    if "inp" not in st.session_state:
        st.session_state.inp = to_state(PR.BY_KEY[PR.DEFAULT_KEY].inputs)
        st.session_state.preset = PR.DEFAULT_KEY
        st.session_state.preset_note = ""
        code = st.query_params.get("case")
        if code:                           # a case in the page link (Put this case in the page link)
            try:
                open_case(CF.decode_link(code))
                st.session_state.case_msg = ("success", "Opened the case in this page's link.")
            except ValueError as exc:
                st.session_state.case_msg = ("error", f"The case in this page's link was not opened: {exc}.")
    d = st.session_state.inp
    if "face" not in d:                    # a session that started before the additions
        base = to_state(E.Inputs())
        for k, v in base.items():
            d.setdefault(k, v)
    if "cols" not in d:                    # a session that started before the column option was added
        n3 = int(d["n3"])
        sz = sizes_for(d, n3, n3)
        d.update(cols="square", n3y=None, n3y_value=n3, dx_value=sz[0] if sz else float(d["R"]) / n3,
                 dy_value=sz[1] if sz else float(d["R"]) / n3)


def apply_preset(key: str):
    p = PR.BY_KEY[key]
    st.session_state.preset = key
    st.session_state.case_name = ""
    if p.inputs is None:
        st.session_state.preset_note = "empty"
        return
    st.session_state.preset_note = ""
    st.session_state.inp = to_state(p.inputs)
    _write_widgets()
    if p.action == "nu_sweep":
        st.session_state.open_nu_sweep = True


def _on_preset():
    apply_preset(st.session_state.w_preset)


def _on_reset():
    apply_preset(st.session_state.preset)


# ---------------------------------------------------------------------------------------------------------------
# cases: save, open (a file or a case code, from this app or the web app) and the page link (casefile)
# ---------------------------------------------------------------------------------------------------------------
CHOICE_KEYS = {"strength": STRENGTH, "closure": CLOSURE, "form": FORM, "tension": TENSION, "iteration": ITER,
               "cols": COLS, "face": FACE, "wall": WALL, "water": WATER}


def _state_value(k: str, base, v):
    """v for the state key k if it has the right type (else the base value): choices from their lists, free text
    for the points, numbers for numbers."""
    num = isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
    if k in CHOICE_KEYS:
        return v if isinstance(v, str) and v in CHOICE_KEYS[k] else base
    if k in ("profile_text", "piezo_text"):
        return v if isinstance(v, str) and len(v) <= 20000 else base
    if k in OPTIONAL:
        return None if v is None else float(v) if num else base
    if isinstance(base, bool):
        return v if isinstance(v, bool) else base
    if isinstance(base, int):
        return int(v) if num and float(v).is_integer() else base
    if isinstance(base, float):
        return float(v) if num else base
    return base                            # engine fields such as model: set from the inputs


def open_case(c: dict):
    """Open a case: its engine inputs, with this app's page state (benches, sizes in metres...) where that gives the
    same inputs. Raises ValueError if the case cannot be read."""
    inp = CF.inputs_from_case(c)
    d = to_state(inp)
    ss = c.get("streamlit_state")
    if isinstance(ss, dict):
        cand = dict(d)
        for k, v in ss.items():
            if k in cand:
                cand[k] = _state_value(k, cand[k], v)
        try:
            if to_inputs(cand) == inp:
                d = cand
        except (ValueError, KeyError, TypeError, ZeroDivisionError):
            pass
    st.session_state.inp = d
    if c.get("preset") in PR.BY_KEY:
        st.session_state.preset = c["preset"]
    st.session_state.preset_note = ""
    st.session_state.case_name = str(c.get("name") or "")[:80]
    _write_widgets()


def _open_file():
    key = f"case_file_{st.session_state.get('case_n', 0)}"
    f = st.session_state.get(key)
    if f is None:
        return
    try:
        data = f.getvalue()
        if len(data) > 2_000_000:
            raise ValueError("the file is over 2 MB, too large for a case")
        open_case(CF.read_case(data))
        st.session_state.case_msg = ("success", f"Opened {f.name}.")
    except ValueError as exc:
        st.session_state.case_msg = ("error", f"Not opened: {exc}.")
    st.session_state.case_n = st.session_state.get("case_n", 0) + 1     # an empty uploader next time


def _open_code():
    text = str(st.session_state.get("case_code", "")).strip()
    if not text:
        st.session_state.case_msg = ("warning", "Paste a case code into the box first.")
        return
    try:
        c = CF.read_case(text)
        open_case(c)
        st.session_state.case_msg = ("success", f"Opened {'“' + c['name'] + '”' if c.get('name') else 'the case'}.")
    except ValueError as exc:
        st.session_state.case_msg = ("error", f"Not opened: {exc}.")


def _show_code(code: str):
    st.session_state.case_code = code
    st.session_state.case_msg = ("info", "The case code is in the box: copy it.")


def _link(code: str):
    st.query_params["case"] = code
    st.session_state.case_msg = ("info", "This page's address now holds the case: copy it from the address bar to share "
                                         "it with people who can open this app.")


def case_box(inp: E.Inputs):
    """Sidebar: save the case, open a case file or code, put the case in the page link."""
    A = st.session_state.get("analysis")
    if (A is None or A.inputs != inp) and not input_problems(st.session_state.inp) \
            and inp.n3 * E.columns_across(inp, inp.n3) <= LIVE_MAX_N * LIVE_MAX_N:
        A = _safe_run(inp)[0]              # live mesh: the run the page makes next (cached)
    res = CF.results_of(A) if A is not None and A.inputs == inp else None
    c = CF.make_case(inp, name=st.session_state.get("case_name", ""), preset=st.session_state.preset,
                     ui_state=st.session_state.inp, results=res)
    with st.sidebar.expander("Save, open or share a case", expanded=False):
        st.download_button("Save case (.json)", CF.case_json(c), CF.file_name(c), "application/json", key="case_save",
                           width="stretch")
        st.file_uploader("Open a case file (saved by this app or the web app)", type=["json"],
                         key=f"case_file_{st.session_state.get('case_n', 0)}", on_change=_open_file)
        st.text_area("Case code: copy it to keep or send this case; paste one and press Open code.", key="case_code",
                     height=90)
        c1, c2 = st.columns(2)
        c1.button("Show this case's code", key="case_code_make", on_click=_show_code, args=(CF.case_json(c),))
        c2.button("Open code", key="case_code_open", on_click=_open_code)
        st.button("Put this case in the page link", key="case_link", on_click=_link, args=(CF.encode_link(c),),
                  width="stretch")
        st.caption("A case holds the inputs (and the FoS when the results are for them). The web app opens the same "
                   "files and codes. The report is on the Calculator page.")
    msg = st.session_state.pop("case_msg", None)
    if msg:
        getattr(st.sidebar, msg[0])(msg[1])


# ---------------------------------------------------------------------------------------------------------------
# widgets
# ---------------------------------------------------------------------------------------------------------------
def _num(container, name: str, disabled: bool = False, label: str = ""):
    label, unit, lo, hi, step, fmt = (label or NUM[name][0],) + NUM[name][1:]
    d = st.session_state.inp
    nk, sk = f"nb_{name}", f"sl_{name}"
    if nk not in st.session_state:
        st.session_state[nk] = d[name]
    if sk not in st.session_state:
        st.session_state[sk] = _clamp(name, d[name])

    def from_nb():
        v = st.session_state[nk]
        st.session_state.inp[name] = v
        st.session_state[sk] = _clamp(name, v)

    def from_sl():
        v = st.session_state[sk]
        st.session_state.inp[name] = v
        st.session_state[nk] = v

    rng = f"{lo:g} to {hi:g}"
    unit_txt = "" if unit == "-" else f" ({unit})"
    kw = dict(min_value=lo, max_value=hi, step=step, key=nk, on_change=from_nb, disabled=disabled,
              help=f"Valid range {rng}{unit_txt}.")
    if fmt:
        kw["format"] = fmt
    container.number_input(f"{label}{unit_txt} · {rng}", **kw)
    container.slider(f"{label} slider", min_value=lo, max_value=hi, step=step, key=sk, on_change=from_sl,
                     disabled=disabled, label_visibility="collapsed")


def _num_opt(container, name: str, disabled: bool = False):
    """A number with no default: an empty box until a value is entered (no slider)."""
    label, unit, lo, hi, step, fmt = NUM[name]
    k = f"nb_{name}"
    if k not in st.session_state:
        st.session_state[k] = st.session_state.inp[name]

    def cb():
        st.session_state.inp[name] = st.session_state[k]

    unit_txt = "" if unit == "-" else f" ({unit})"
    container.number_input(f"{label}{unit_txt} · {lo:g} to {hi:g}", min_value=lo, max_value=hi, step=step, format=fmt,
                           key=k, on_change=cb, disabled=disabled, placeholder="enter a value (no default)",
                           help=f"Valid range {lo:g} to {hi:g}{unit_txt}. There is no default: enter your value.")


def _text(container, name: str, label: str, height: int = 120):
    k = f"w_{name}"
    if k not in st.session_state:
        st.session_state[k] = st.session_state.inp[name]

    def cb():
        st.session_state.inp[name] = st.session_state[k]

    container.text_area(label, key=k, on_change=cb, height=height)


def _choice(container, name: str, label: str, options: dict, disabled=False, horizontal=False, on_change=None):
    k = f"w_{name}"
    if k not in st.session_state:
        st.session_state[k] = st.session_state.inp[name]

    def cb():
        st.session_state.inp[name] = st.session_state[k]
        if on_change:
            on_change()

    container.radio(label, list(options), format_func=options.get, key=k, on_change=cb, disabled=disabled,
                    horizontal=horizontal)


def _check(container, name: str, label: str, help: str = "", disabled: bool = False):
    k = f"w_{name}"
    if k not in st.session_state:
        st.session_state[k] = st.session_state.inp[name]

    def cb():
        st.session_state.inp[name] = st.session_state[k]

    container.checkbox(label, key=k, on_change=cb, disabled=disabled, help=help or None)


def _on_cols():
    """Switch the 3D column option and keep the current columns: the counts carry over, and sizes in metres start
    from the current plan size (rounded up to 4 digits, which gives the same counts)."""
    d = st.session_state.inp
    new = st.session_state["w_cols"]
    if new == d["cols"]:
        return
    n, ny = col_grid(d)
    if new == "size":
        sz = sizes_for(d, n, ny)
        if sz:
            d["dx_value"], d["dy_value"] = sz
    else:                                  # counts outside 2 to 320 (sizes the sidebar flags) are clamped
        d["n3"] = _clamp("n3", n)
        if new == "counts":
            d["n3y_value"] = _clamp("n3y_value", ny)
    d["cols"] = new
    _write_widgets()


def _on_face():
    """Switch the slope-face option. To points: the current face as points. To planar: keep H and beta unless the
    face was changed to a valid one (one bench, or points other than the planar face's), whose height and overall
    angle are then used. Same as setFace in webapp/index.html."""
    d = st.session_state.inp
    new = st.session_state["w_face"]
    if new == d["face"]:
        return
    H, beta, prof = face_fields(d)
    if new == "profile":
        tb = math.tan(math.radians(beta))
        d["profile_text"] = points_text(prof or ((0.0, 0.0), (-H / tb, H)))
    elif new == "planar":
        keep = None
        if d["face"] == "benches" and int(d["bench_n"]) == 1:
            keep = (float(d["bench_h"]), float(d["bench_beta"]))
        elif prof:
            tb = math.tan(math.radians(float(d["beta"])))
            ok = not (d["face"] == "profile" and str(d["profile_text"]).strip()
                      == points_text(((0.0, 0.0), (-float(d["H"]) / tb, float(d["H"])))))
            if ok:
                try:
                    E.check_additions(E.Inputs(profile=prof))
                except ValueError:
                    ok = False
            if ok:
                keep = (H, float(f"{beta:.12g}"))
        if keep:
            d["H"], d["beta"] = keep
    d["face"] = new
    _write_widgets()


def _on_iteration():
    if st.session_state.inp["iteration"] == "native_v15":
        st.session_state.inp["form"] = "FH"
        st.session_state.inp["tension"] = "none"
        st.session_state["w_form"] = "FH"
        st.session_state["w_tension"] = "none"


def sidebar() -> E.Inputs:
    """Draw the input sidebar and return the current inputs."""
    init_state()
    # Each page draws its own sidebar widgets, and Streamlit resets widgets it sees for the first time on a page to
    # their defaults. The inputs live in st.session_state.inp, so every widget is set from them before it is drawn
    # (on_change callbacks keep st.session_state.inp up to date).
    _write_widgets()
    st.session_state.w_preset = st.session_state.preset
    sb = st.sidebar
    labels = {p.key: f"{p.group}: {p.label}" for p in PR.PRESETS}
    sb.selectbox("Preset", list(labels), format_func=labels.get, key="w_preset", on_change=_on_preset)
    sb.button("Reset inputs to preset", on_click=_on_reset, width="stretch")
    d = st.session_state.inp
    case_box(to_inputs(d))
    with sb.expander("Slope", expanded=False):
        if "w_face" not in st.session_state:
            st.session_state["w_face"] = d["face"]
        st.radio("Slope face", list(FACE), format_func=FACE.get, key="w_face", on_change=_on_face, horizontal=True)
        if d["face"] == "planar":
            _num(st, "H")
            _num(st, "beta")
        elif d["face"] == "benches":
            for n in ("bench_n", "bench_h", "bench_beta", "berm_w"):
                _num(st, n)
        else:
            _text(st, "profile_text", "Face points x, z (m), one per line, from the toe (0, 0) up to the crest edge")
        if d["face"] != "planar":
            H, beta, prof = face_fields(d)
            if prof:
                st.caption(f"Crest edge at x = {prof[-1][0]:.2f} m, H = {H:.2f} m; overall angle {beta:.2f}° from the "
                           f"toe ({len(prof)} points).")
        _num(st, "crest_slope")
        _choice(st, "wall", "Wall in plan", WALL, horizontal=True)
        if d["wall"] != "straight":
            _num_opt(st, "plan_r")
        st.caption(("Concave wall: the toe line is a circle of radius R_p around the pit, which lies inside it; the "
                    "section is swept around its axis. " if d["wall"] == "concave" else
                    "Convex wall (a nose): the toe line is a circle of radius R_p around an axis behind the slope. "
                    if d["wall"] == "convex" else "")
                   + "Ground behind the crest edge: flat (α = 0) or rising at α. The 2D section is the central one, "
                   "y = 0.")
    with sb.expander("Slip surface", expanded=False):
        for n in ("xc", "zc", "R", "elong", "mid_length"):
            _num(st, n)
        try:
            half = E.y_extent(_geom(d))
            st.caption(f"Half-width across y: {'L/2 + ' if float(d['mid_length']) > 0 else ''}λR = {half:.2f} m"
                       + (" (the sphere)." if float(d["elong"]) == 1 and float(d["mid_length"]) == 0 else
                          ". The central section y = 0 stays the circle (x_c, z_c, R), as in 2D."))
        except ValueError:
            pass
        _check(st, "base_circle", "Base circle: grid from exit to entry, floor columns allowed",
               "Off by default (thesis grid). Needed for circles that cross the toe level behind the toe (x_c ≤ 0), "
               "exit on the floor ahead of the toe, or stay above the toe level (face circles, R ≤ z_c).")
    hb_off = d["strength"] == "MC"
    with sb.expander("Rock mass", expanded=False):
        _num(st, "gamma")
        for n in ("sigma_ci", "mi", "gsi", "D"):
            _num(st, n, disabled=hb_off)
        _check(st, "a_on", "Override a (otherwise a from GSI, eq. 2.5)", disabled=hb_off)
        _num(st, "a_value", disabled=hb_off or not d["a_on"])
        _check(st, "dz_on", "Apply D only within T of the face (D = 0 beyond)", disabled=hb_off,
               help="Hoek et al. (2002): D is for the blast-damaged zone near the face. Bases closer than T to the face "
                    "use D; the rest use D = 0. T has no default: enter it from your site.")
        if d["dz_on"] and not hb_off:
            _num_opt(st, "d_zone")
        mb, s, a = to_inputs(d).hb
        st.caption(f"m_b = {mb:.6f}, s = {s:.7f}, a used = {a:.6f}" + (" (not used by Mohr-Coulomb)" if hb_off else ""))
    with sb.expander("Strength model", expanded=True):
        _choice(st, "strength", "Strength model", STRENGTH)
        _choice(st, "closure", "Closure (EJC)", CLOSURE, disabled=d["strength"] != "EJC")
        _num(st, "nu", disabled=d["strength"] != "EJC")
        _num(st, "c_mc", disabled=not hb_off)
        _num(st, "phi_mc", disabled=not hb_off)
    with sb.expander("Water, tension crack and seismic load", expanded=False):
        _choice(st, "water", "Pore pressure", WATER, horizontal=True)
        _num(st, "Hw", disabled=d["water"] != "Wet")
        _num(st, "beta_w", disabled=d["water"] != "Wet")
        if d["water"] == "Piezo":
            _text(st, "piezo_text", "Piezometric line x, z (m), one point per line, x decreasing (into the slope)", 100)
        _num(st, "ru_value", disabled=d["water"] != "ru")
        st.caption({"Wet": "u = γ_w max(0, z_w − z_b), γ_w = 9.81 kN/m³; z_w = min(H_w, −x tan β_w).",
                    "Piezo": "u = γ_w max(0, z_w − z_b) with z_w on the line through the points (constant beyond its "
                             "ends).",
                    "ru": "u = r_u γ h for every column and slice (h = column height)."}.get(d["water"], "u = 0."))
        _check(st, "crack_on", "Tension crack at the back of the mass")
        if d["crack_on"]:
            _num_opt(st, "crack_depth")
            _num(st, "crack_water")
        _num(st, "kh")
        st.caption(("The crack is vertical, where the slip circle on y = 0 first lies z_t below the ground from the "
                    "back; the mass ends there, and water in it pushes on its face with γ_w h_w²/2 at h_w/3. "
                    if d["crack_on"] else "")
                   + "Seismic: a pseudo-static horizontal force k_h W at each column's centroid, towards the pit "
                   "(0 = none).")
    with sb.expander("Mesh", expanded=True):
        if "w_cols" not in st.session_state:
            st.session_state["w_cols"] = d["cols"]
        st.radio("3D columns: plan size Δx × Δy", list(COLS), format_func=COLS.get, key="w_cols", on_change=_on_cols,
                 horizontal=True)
        if d["cols"] != "size":
            _num(st, "n3", label="N, 3D columns along x" if d["cols"] == "counts" else "")
        if d["cols"] == "counts":
            _num(st, "n3y_value")
        if d["cols"] == "size":
            _num(st, "dx_value")
            _num(st, "dy_value")
        problem = mesh_problem(d)
        if problem:
            st.error(problem)
        st.caption(cols_text(d))
        _num(st, "n2")
        span = span_of(d)
        if span is not None:
            st.caption(f"2D slices: Δx = (x0 − x1)/N_2D = {span / int(d['n2']):.3f} m, unit thickness.")
        st.caption(f"Live update up to {LIVE_MAX_N} × {LIVE_MAX_N} columns; otherwise press Run.")
    with sb.expander("Options", expanded=False):
        native = d["iteration"] == "native_v15"
        _choice(st, "iteration", "Iteration path", ITER, on_change=_on_iteration)
        _choice(st, "form", "Global equilibrium", FORM, disabled=native)
        _choice(st, "tension", "Tension sensitivity", TENSION, disabled=native)
        if native:
            st.caption("The native v15 path supports F_H only, without C0/T0 (Table 5.7 values come from the default "
                       "path).")
        st.caption("Trial-stress projection (thesis Section 4.6) is always on; a run is accepted only if no base is "
                   "projected at convergence.")
    for msg in input_problems(d):
        if msg != mesh_problem(d):
            sb.error(msg)
    return to_inputs(d)


# ---------------------------------------------------------------------------------------------------------------
# running
# ---------------------------------------------------------------------------------------------------------------
@st.cache_data(show_spinner=False, max_entries=64)
def cached_run(inp: E.Inputs) -> E.Analysis:
    return E.run(inp)


@st.cache_data(show_spinner=False, max_entries=64)
def cached_run_3d(inp: E.Inputs) -> E.LemResult:
    return E.run_3d(inp)


@st.cache_data(show_spinner=False, max_entries=64)
def cached_run_2d(inp: E.Inputs) -> E.LemResult:
    return E.run_2d(inp)


def _safe_run(inp: E.Inputs):
    try:
        return cached_run(inp), None
    except ValueError as exc:
        return None, str(exc)


def analysis(inp: E.Inputs):
    """Run (live up to 40 x 40 columns, otherwise on the Run button) and return (analysis or None, stale flag)."""
    problems = input_problems(st.session_state.inp)
    if problems:
        st.error("Input error: " + " ".join(problems))
        return None, False
    n, ny = inp.n3, E.columns_across(inp, inp.n3)
    live = n * ny <= LIVE_MAX_N * LIVE_MAX_N
    p = PR.BY_KEY[st.session_state.preset]
    if live:
        st.sidebar.caption(f"Live update is on (up to {LIVE_MAX_N} × {LIVE_MAX_N} columns).")
        A, err = _safe_run(inp)
        st.session_state.analysis, st.session_state.analysis_error = A, err
    else:
        msg = "Run (slow preset, may take a minute)" if p.slow else "Run"
        if st.sidebar.button(msg, type="primary", width="stretch"):
            with st.spinner(f"Solving 3D {n}x{ny} columns and 2D {inp.n2} slices..."):
                A, err = _safe_run(inp)
            st.session_state.analysis, st.session_state.analysis_error = A, err
    A = st.session_state.get("analysis")
    err = st.session_state.get("analysis_error")
    if err:
        st.error(f"Input error: {err}")
        return None, False
    stale = A is not None and A.inputs != inp
    if A is None:
        st.info(f"The mesh ({n} × {ny} columns) is above the live limit of {LIVE_MAX_N} × {LIVE_MAX_N}: press **Run** "
                "in the sidebar to solve.")
    elif stale:
        st.warning("The inputs have changed since the last run; the results below are for the previous inputs. "
                   "Press **Run** to update.")
    return A, stale


def preset_box():
    p = PR.BY_KEY[st.session_state.preset]
    if p.inputs is None:
        st.warning(f"**{p.label}**: {p.source} The inputs below are unchanged.")
        return
    lines = [f"**Preset: {p.group}, {p.label}.** {p.source}"]
    if p.info:
        lines.append(p.info)
    if p.targets:
        lines.append(f"Targets: {p.targets}")
    if p.reference:
        lines.append(f"Reference only (not a target): {p.reference}")
    if p.slow:
        lines.append("This preset uses fine meshes and is slow.")
    st.info("  \n".join(lines))


def fmt(v, d=6):
    return "withheld" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f"{v:.{d}f}"


def status_badge(r: E.LemResult, label: str):
    if r.status == "ACCEPTED":
        st.success(f"**{label}: ACCEPTED.** {r.reason}")
    elif r.status == "NOT APPLICABLE":
        st.info(f"**{label}: NOT APPLICABLE.** {r.reason}")
    else:
        st.error(f"**{label}: {r.status}, no accepted FoS.** {r.reason}")


# Bases with tensile effective normal stress at an accepted state (QA/QC finding S3). The summary warns above a quarter
# of the bases: the thesis's wet baseline has 30% and C0 lowers its FoS by 24%, the dry baseline has 8% and C0 lowers
# it by 2.5% (thesis §5.8). The same rule is in webapp/index.html.
TENSILE_WARN = 0.25


def tensile_count(r: E.LemResult):
    """(number of bases with sigma_n < 0, number of bases) at an accepted state, else None."""
    if r is None or r.fos is None or r.sn is None:
        return None
    return int((r.sn < 0).sum()), len(r.sn)


def tensile_text(t) -> str:
    if t is None:
        return "–"
    n, total = t
    share = "0%" if n == 0 else "<1%" if 100 * n / total < 1 else f"{round(100 * n / total)}%"
    return f"{n} of {total} ({share})"


def tensile_note(inp: E.Inputs, r3: E.LemResult, r2: E.LemResult):
    """(kind, markdown) for the summary, or None: as tensileNote in webapp/index.html."""
    if inp.tension != "none":
        return None
    found = [(lab, t) for lab, t in (("3D", tensile_count(r3)), ("2D", tensile_count(r2))) if t and t[0] > 0]
    if not found:
        return None
    lst = ", ".join(f"{lab} {tensile_text(t)}" for lab, t in found)
    if any(n / total > TENSILE_WARN for _, (n, total) in found):
        return "warning", (f"**Tensile bases: {lst}.** More than a quarter of the bases carry tensile effective normal "
                           "stress (σ′_n < 0) and the FoS counts their resistance. Compare the C0 and T0 options before "
                           "relying on it (thesis §5.8, Table 5.7).")
    return "info", (f"**Tensile bases: {lst}.** These bases carry tensile effective normal stress (σ′_n < 0) and the FoS "
                    "counts their resistance; the C0 and T0 options show how much it depends on them (thesis §5.8).")


def _len2(v: float) -> str:
    return f"{v:.2f}" if v == 0 or abs(v) >= 0.01 else f"{v:.2g}"


def vertical_face_note(inp: E.Inputs):
    """Markdown warning when the slip circle ends below the ground (engine.vertical_faces), else None: no force acts on
    such a face, so the FoS is for a mass with an open crack at the back or a free face at the front, not for a
    circular slip surface. As faceNote in webapp/index.html."""
    f = E.vertical_faces(inp)
    back, front = f["back"] > 0, f["front"] > 0
    if not back and not front:
        return None
    bk, fr = f"a vertical back face {_len2(f['back'])} m high", f"a vertical front face {_len2(f['front'])} m high"
    title = "Vertical back and front faces" if back and front else "Vertical back face" if back else "Vertical front face"
    if back and front:
        text = (f"The slip circle ends below the ground at both extremes (on y = 0), so the analysed mass has {bk} at "
                f"x_c − R = {inp.xc - inp.R:.2f} m and {fr} at x_c + R = {inp.xc + inp.R:.2f} m. No force acts on them: the "
                "back face works as an open tension crack and the front face as if the ground ahead of it were removed, so "
                "the FoS is for that mass, not for a circular slip surface. The arc meets the ground before its extremes "
                f"only when z_c is at least the ground height there ({_len2(inp.zc + f['back'])} m and "
                f"{_len2(inp.zc + f['front'])} m).")
    else:
        end = "back" if back else "front"
        x_end = inp.xc - inp.R if back else inp.xc + inp.R
        text = (f"The slip circle ends below the ground at its {end} extreme (x_c {'−' if back else '+'} R = {x_end:.2f} m "
                f"on y = 0), so the analysed mass has {bk if back else fr}. No force acts on it: it works "
                f"{'as an open tension crack' if back else 'as if the ground ahead of it were removed'}, and the FoS is for "
                f"that mass, not for a circular slip surface. The arc meets the ground before its {end} extreme only when "
                f"z_c is at least the ground height there ({_len2(inp.zc + (f['back'] if back else f['front']))} m).")
    return (f"**{title}.** {text} The critical-surface search (Studies) tries only circles that meet the ground at both "
            "exits.")


def checks_frame(r: E.LemResult) -> pd.DataFrame:
    return pd.DataFrame([{"Check": c.name, "Value": f"{c.value:.3e}" if abs(c.value) < 1e-3 or abs(c.value) >= 1e4
                          else f"{c.value:.6g}", "Limit": c.limit, "Status": "PASS" if c.passed else "FAIL"}
                         for c in r.checks])


def color_status(v):
    return {"PASS": "background-color: #d4edda", "FAIL": "background-color: #f8d7da",
            "INFO": "background-color: #e9ecef", "SKIPPED": "background-color: #fff3cd"}.get(v, "")


# ---------------------------------------------------------------------------------------------------------------
# critical-surface search (engine.critical_search; QA/QC finding R1). Messages as in webapp/index.html.
# ---------------------------------------------------------------------------------------------------------------
SEARCH_EDGE = {
    "front limit": lambda s: f"lies on the front limit ({s['front']:g} H ahead of the toe): widen it and search again",
    "back limit": lambda s: f"lies on the back limit ({s['back']:g} H behind the crest): widen it and search again",
    "crest": lambda s: "has its front exit at the crest edge",
    "exits meet": lambda s: "has its two exits next to each other",
    "vertical back": lambda s: "has an arc that is vertical at the back exit, the steepest arc the search tries",
    "flat": lambda s: "is nearly flat, the flattest arc the search tries",
    "coverage": lambda s: ("lies on the default grid's coverage limit: a flatter or longer circle would cover less than "
                           f"{round(s['coverage'] * 100)} % of the grid's span from x_c − R to the toe-level crossing, "
                           "so few slices or columns would resolve it, and the search skips it (the base-circle option "
                           "meshes the mass itself)"),
    "masses": lambda s: "lies on the base-circle grid's limit of one sliding mass",
    "grid": lambda s: ("lies on the default grid's limit (it analyses circles that cross the toe level at or ahead of the "
                       "toe); the base-circle option also analyses circles that pass below the toe, or stay above the "
                       "toe level"),
    "withheld": lambda s: "lies next to surfaces whose results were withheld (not accepted)",
}
SURFACE_KEYS = ("xc", "zc", "R")


def _without_surface(inp: E.Inputs) -> E.Inputs:
    return inp.with_(xc=0.0, zc=0.0, R=0.0)


def surface_in_limits(r: dict) -> bool:
    """True if the surface's x_c, z_c and R lie within the sidebar's input limits."""
    return all(NUM[k][2] <= r[k] <= NUM[k][3] for k in SURFACE_KEYS)


def search_notes(current: E.Inputs, searched: E.Inputs, res: dict) -> list[tuple[str, str]]:
    """(kind, markdown) notes on a search result: stale inputs, the prescribed circle's vertical back face, the limits
    the critical surfaces lie on and the mesh check (as the web app's search panel)."""
    notes = []
    if _without_surface(current) != _without_surface(searched):
        notes.append(("warning", "The inputs have changed since this search (other than the surface). Search again for "
                                 "the current inputs."))
    P, s = res["prescribed"], res["settings"]
    if P["vertical_back"] > 0 or P["vertical_front"] > 0:
        faces = " and ".join(t for t in (f"a vertical back face {_len2(P['vertical_back'])} m high" if P["vertical_back"] > 0
                                         else "", f"a vertical front face {_len2(P['vertical_front'])} m high"
                                         if P["vertical_front"] > 0 else "") if t)
        both = P["vertical_back"] > 0 and P["vertical_front"] > 0
        notes.append(("info", f"The prescribed circle ends below the ground on y = 0, so the mass analysed for it has {faces} "
                              f"that {'carry' if both else 'carries'} no force (see the results summary). The search does "
                              "not try such circles."))
    for d in s["dims"]:
        r = res[f"r{d}"]
        if r["message"]:
            notes.append(("warning", f"{d}D: {r['message']}"))
        if r["fos"] is None:
            continue
        if not surface_in_limits(r):
            notes.append(("warning", f"{d}D: the critical surface lies outside the sidebar's input limits (x_c "
                                     f"{NUM['xc'][2]:g} to {NUM['xc'][3]:g} m, z_c {NUM['zc'][2]:g} to {NUM['zc'][3]:g} m, "
                                     f"R {NUM['R'][2]:g} to {NUM['R'][3]:g} m), so it cannot be used as the prescribed "
                                     "surface."))
        edges = []
        for e in r["edges"]:
            text = SEARCH_EDGE.get(e["why"])
            if e["why"] != "toe" and text and text(s) not in edges:
                edges.append(text(s))
        if edges:
            notes.append(("warning", f"{d}D: the critical {'circle' if d == 2 else 'sphere'} " + "; it ".join(edges) + "."))
        c = r["check"]
        mesh = f"{c['n']} x {c['ny']} columns" if d == 3 else f"{c['n']} slices"
        if c["change"] is not None and abs(c["change"]) > 0.005:
            change = f"{'+' if c['fos'] >= r['fos'] else '−'}{abs((c['fos'] / r['fos'] - 1) * 100):.1f} %"
            notes.append(("warning", f"{d}D: on {mesh} the critical surface's FoS changes by {change}, more than the "
                                     "0.5 % screening threshold of eq. 4.17. The current mesh under-resolves this surface, "
                                     "and the search may favour surfaces whose FoS the mesh underestimates. Refine the "
                                     f"mesh ({'N_2D' if d == 2 else 'N' if searched.n3y is None else 'N and N_y'}) "
                                     "and search again."))
        elif c["fos"] is None:
            notes.append(("warning", f"{d}D: on {mesh} the critical surface gave no accepted FoS ({c['status']})."))
    return notes


def use_surface(xc: float, zc: float, R: float):
    """Button callback: make (x_c, z_c, R) the prescribed surface."""
    d = st.session_state.inp
    d["xc"], d["zc"], d["R"] = float(xc), float(zc), float(R)
    _write_widgets()
