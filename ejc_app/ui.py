"""Shared Streamlit helpers: input state, the sidebar (number box + slider kept in sync), presets and running.

Input values live in ``st.session_state["inp"]`` (a plain dict, kept across pages). Each widget has its own key,
initialised from that dict when it is first drawn on a page and copied back to it in an ``on_change`` callback.
"""
from __future__ import annotations

import math
from dataclasses import asdict

import pandas as pd
import streamlit as st

import engine as E
import presets as PR

LIVE_MAX_N = 40          # live update only for N <= 40 (brief, page 1)
DISCLAIMER = "Research tool for the prescribed-surface comparison in Mahat (2026) thesis; not for design."

# name: (label, unit, min, max, step, format)
NUM = {
    "H": ("H, slope height", "m", 10.0, 1000.0, 1.0, "%.2f"),
    "beta": ("β, slope face angle", "deg", 10.0, 89.0, 0.5, "%.2f"),
    "xc": ("x_c, sphere centre", "m", -1000.0, 2000.0, 0.5, "%.7f"),
    "zc": ("z_c, sphere centre", "m", 1.0, 3000.0, 0.5, "%.7f"),
    "R": ("R, sphere radius", "m", 1.0, 4000.0, 0.5, "%.7f"),
    "gamma": ("γ, unit weight", "kN/m³", 10.0, 40.0, 0.1, "%.2f"),
    "sigma_ci": ("σ_ci, intact UCS", "MPa", 1.0, 300.0, 0.1, "%.2f"),
    "mi": ("m_i", "-", 1.0, 40.0, 0.5, "%.2f"),
    "gsi": ("GSI", "-", 10.0, 100.0, 1.0, "%.1f"),
    "D": ("D, disturbance factor", "-", 0.0, 1.0, 0.05, "%.2f"),
    "a_value": ("a (override value)", "-", 0.30, 1.00, 0.001, "%.6f"),
    "nu": ("ν, Poisson's ratio", "-", 0.01, 0.49, 0.001, "%.3f"),
    "Hw": ("H_w, water-table height", "m", 0.0, 1000.0, 1.0, "%.1f"),
    "beta_w": ("β_w, water-table angle", "deg", 1.0, 89.0, 0.5, "%.1f"),
    "n3": ("N, 3D mesh (N x N half-model)", "-", 2, 320, 1, None),
    "n2": ("N_2D, 2D slices", "-", 2, 10000, 1, None),
}
STRENGTH = {"EJC": "EJC (Extended Jiang Criterion)", "GHB": "GHB (conventional generalised Hoek-Brown)"}
CLOSURE = {"moscarda": "Moscarda: σ2 = ν(σ1 + σ3) (adopted)",
           "original": "Original: σ2 = σ3/ν − σ1 (admissibility check only)",
           "s2s3": "σ2 = σ3 (check: must reproduce GHB)"}
FORM = {"FH": "F_H, force form (governing, eq. 4.14)", "FM": "F_M, moment form (eq. 4.15, own iteration, 3D)"}
TENSION = {"none": "None", "C0": "C0: remove c_i on bases with σ_n < 0", "T0": "T0: remove all basal shear there"}
ITER = {"default": "Default (§1.4): F from F_H/F_M, then σ_n", "native_v15": "Native v15 iteration path"}
CHOICES = {"strength": STRENGTH, "closure": CLOSURE, "form": FORM, "tension": TENSION, "iteration": ITER}


# ---------------------------------------------------------------------------------------------------------------
# Inputs <-> state dict
# ---------------------------------------------------------------------------------------------------------------
def to_state(inp: E.Inputs) -> dict:
    d = asdict(inp)
    d["a_on"] = inp.a_override is not None
    d["a_value"] = 0.5 if inp.a_override is None else float(inp.a_override)
    d["strength"] = "GHB" if inp.model == "ghb" else "EJC"
    d["closure"] = {"ejc_original": "original", "ejc_sigma2_eq_sigma3": "s2s3"}.get(inp.model, "moscarda")
    d["water"] = "Wet" if inp.wet else "Dry"
    return d


def to_inputs(d: dict) -> E.Inputs:
    if d["strength"] == "GHB":
        model = "ghb"
    else:
        model = {"moscarda": "ejc_moscarda", "original": "ejc_original", "s2s3": "ejc_sigma2_eq_sigma3"}[d["closure"]]
    native = d["iteration"] == "native_v15"
    return E.Inputs(H=float(d["H"]), beta=float(d["beta"]), xc=float(d["xc"]), zc=float(d["zc"]), R=float(d["R"]),
                    gamma=float(d["gamma"]), sigma_ci=float(d["sigma_ci"]), mi=float(d["mi"]), gsi=float(d["gsi"]),
                    D=float(d["D"]), a_override=float(d["a_value"]) if d["a_on"] else None, nu=float(d["nu"]),
                    model=model, wet=d["water"] == "Wet", Hw=float(d["Hw"]), beta_w=float(d["beta_w"]),
                    n3=int(d["n3"]), n2=int(d["n2"]), form="FH" if native else d["form"],
                    tension="none" if native else d["tension"], iteration=d["iteration"])


WIDGET_NAMES = list(NUM) + ["a_on", "strength", "closure", "water", "form", "tension", "iteration"]


def _wkeys(name: str) -> list[str]:
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


def apply_preset(key: str):
    p = PR.BY_KEY[key]
    st.session_state.preset = key
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
# widgets
# ---------------------------------------------------------------------------------------------------------------
def _num(container, name: str, disabled: bool = False):
    label, unit, lo, hi, step, fmt = NUM[name]
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


def _on_iteration():
    if st.session_state.inp["iteration"] == "native_v15":
        st.session_state.inp["form"] = "FH"
        st.session_state.inp["tension"] = "none"
        st.session_state["w_form"] = "FH"
        st.session_state["w_tension"] = "none"


def sidebar() -> E.Inputs:
    """Draw the input sidebar and return the current inputs."""
    init_state()
    sb = st.sidebar
    labels = {p.key: f"{p.group}: {p.label}" for p in PR.PRESETS}
    if "w_preset" not in st.session_state:
        st.session_state.w_preset = st.session_state.preset
    sb.selectbox("Preset", list(labels), format_func=labels.get, key="w_preset", on_change=_on_preset)
    sb.button("Reset inputs to preset", on_click=_on_reset, width="stretch")
    d = st.session_state.inp
    with sb.expander("Geometry", expanded=False):
        for n in ("H", "beta", "xc", "zc", "R"):
            _num(st, n)
    with sb.expander("Rock mass", expanded=False):
        for n in ("gamma", "sigma_ci", "mi", "gsi", "D"):
            _num(st, n)
        if "w_a_on" not in st.session_state:
            st.session_state.w_a_on = d["a_on"]

        def a_cb():
            st.session_state.inp["a_on"] = st.session_state.w_a_on
        st.checkbox("Override a (otherwise a from GSI, eq. 2.5)", key="w_a_on", on_change=a_cb)
        _num(st, "a_value", disabled=not d["a_on"])
        mb, s, a = to_inputs(d).hb
        st.caption(f"m_b = {mb:.6f}, s = {s:.7f}, a used = {a:.6f}")
    with sb.expander("Strength model", expanded=True):
        _choice(st, "strength", "Strength model", STRENGTH)
        _choice(st, "closure", "Closure (EJC)", CLOSURE, disabled=d["strength"] == "GHB")
        _num(st, "nu", disabled=d["strength"] == "GHB")
    with sb.expander("Water", expanded=False):
        _choice(st, "water", "Condition", {"Dry": "Dry (u = 0)", "Wet": "Wet (hydrostatic)"}, horizontal=True)
        _num(st, "Hw", disabled=d["water"] == "Dry")
        _num(st, "beta_w", disabled=d["water"] == "Dry")
    with sb.expander("Mesh", expanded=True):
        _num(st, "n3")
        _num(st, "n2")
        st.caption(f"Live update when N ≤ {LIVE_MAX_N}; otherwise press Run.")
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
    """Run (live for N <= 40, otherwise on the Run button) and return (analysis or None, stale flag)."""
    live = inp.n3 <= LIVE_MAX_N
    p = PR.BY_KEY[st.session_state.preset]
    if live:
        st.sidebar.caption("Live update is on (N ≤ 40).")
        A, err = _safe_run(inp)
        st.session_state.analysis, st.session_state.analysis_error = A, err
    else:
        msg = "Run (slow preset, may take a minute)" if p.slow else "Run"
        if st.sidebar.button(msg, type="primary", width="stretch"):
            with st.spinner(f"Solving 3D {inp.n3}x{inp.n3} and 2D {inp.n2} slices..."):
                A, err = _safe_run(inp)
            st.session_state.analysis, st.session_state.analysis_error = A, err
    A = st.session_state.get("analysis")
    err = st.session_state.get("analysis_error")
    if err:
        st.error(f"Input error: {err}")
        return None, False
    stale = A is not None and A.inputs != inp
    if A is None:
        st.info(f"N = {inp.n3} > {LIVE_MAX_N}: press **Run** in the sidebar to solve.")
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


def checks_frame(r: E.LemResult) -> pd.DataFrame:
    return pd.DataFrame([{"Check": c.name, "Value": f"{c.value:.3e}" if abs(c.value) < 1e-3 or abs(c.value) >= 1e4
                          else f"{c.value:.6g}", "Limit": c.limit, "Status": "PASS" if c.passed else "FAIL"}
                         for c in r.checks])


def color_status(v):
    return {"PASS": "background-color: #d4edda", "FAIL": "background-color: #f8d7da",
            "INFO": "background-color: #e9ecef", "SKIPPED": "background-color: #fff3cd"}.get(v, "")
