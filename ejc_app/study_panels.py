"""Studies on the Calculator page (additions): EJC against GHB, a parameter sweep, a back-analysis and Monte Carlo,
as the web app's Studies panels. Each result is kept with the inputs it ran on; the report includes only the studies
run on the inputs it reports."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import engine as E
import report as RP
import ui

# inputs a study can vary: (state key, engine input, label)
KEYS = [("gsi", "gsi", "GSI"), ("D", "D", "D, disturbance factor"), ("sigma_ci", "sigma_ci", "σci (MPa)"),
        ("mi", "mi", "m_i"), ("nu", "nu", "ν (EJC)"), ("gamma", "gamma", "γ (kN/m³)"), ("c_mc", "c_mc", "c′ (kPa)"),
        ("phi_mc", "phi_mc", "φ′ (°)"), ("Hw", "Hw", "H_w (m)"), ("beta_w", "beta_w", "β_w (°)"), ("ru_value", "ru", "r_u"),
        ("kh", "kh", "k_h, seismic"), ("crack_depth", "crack_depth", "tension crack depth (m)"),
        ("crack_water", "crack_water", "crack water share"), ("d_zone", "d_zone", "T, disturbed zone (m)"),
        ("H", "H", "H (m)"), ("beta", "beta", "β (°)"), ("crest_slope", "crest_slope", "α, ground behind the crest (°)"),
        ("elong", "elong", "λ, end caps (3D)"), ("mid_length", "mid_length", "L, cylindrical middle (m, 3D)")]
GEOMETRY = {"H", "beta", "crest_slope", "elong", "mid_length"}
COLOR = {3: "#2a78d6", 2: "#eb6834"}                 # 3D blue, 2D orange, as the ν sweep and the web app
DASH = {"ejc_moscarda": "solid", "ghb": "dash", "mc": "solid"}
SPEED = ("This engine takes about 0.1 to 1 s per analysis, so a study of many analyses takes minutes; keep this page "
         "open while it runs. The web app (webapp/) runs the same studies, with the same results, faster.")


def _relevant(k: str, d: dict) -> bool:
    if k in ("sigma_ci", "mi", "gsi", "D"):
        return d["strength"] != "MC"
    if k == "nu":
        return d["strength"] == "EJC"
    if k in ("c_mc", "phi_mc"):
        return d["strength"] == "MC"
    if k in ("Hw", "beta_w"):
        return d["water"] == "Wet"
    if k == "ru_value":
        return d["water"] == "ru"
    if k in ("H", "beta"):
        return d.get("face", "planar") == "planar"
    if k == "d_zone":
        return bool(d.get("dz_on")) and d["strength"] != "MC"
    if k in ("crack_depth", "crack_water"):
        return bool(d.get("crack_on"))
    return True


def study_keys(d: dict) -> list[tuple[str, str, str]]:
    return [kk for kk in KEYS if _relevant(kk[0], d)]


def _range(k: str, v: float) -> tuple[float, float]:
    """Default range of a sweep or back-analysis: the current value ± 50 %, inside the valid range."""
    lo, hi = ui.NUM[k][2], ui.NUM[k][3]
    a = lo if v == 0 else max(lo, v - 0.5 * abs(v))
    b = lo + (hi - lo) / 10 if v == 0 else min(hi, v + 0.5 * abs(v))
    return float(f"{a:.6g}"), float(f"{b:.6g}")


def _stale(found: dict, inp: E.Inputs, button: str):
    if found["inputs"] != inp:
        st.caption(f"For the previous inputs; press {button} to update.")


def _bar(text: str):
    bar = st.progress(0.0, text=text)

    def progress(done, total, label):
        bar.progress(min(done / total, 1.0), text=f"{label} ({done} of {total})")
    return bar, progress


def compare_panel(inp: E.Inputs):
    st.write("The EJC (Moscarda closure, at the current ν) and GHB on the current surface, mesh and inputs, in 3D and "
             "2D, with their ratios.")
    if st.button("Compare", key="cmp_run"):
        with st.spinner("Solving EJC and GHB, 3D and 2D…"):
            st.session_state.cmp = {"inputs": inp, "res": E.compare_models(inp)}
    found = st.session_state.get("cmp")
    if not found:
        return
    _stale(found, inp, "Compare")
    res, nu = found["res"], found["inputs"].nu
    e, g = res["models"]["ejc"], res["models"]["ghb"]
    st.dataframe(pd.DataFrame([
        {"": f"EJC, Moscarda (ν = {nu:.4g})", "FoS 3D": e["fos3"], "FoS 2D": e["fos2"], "3D / 2D": e["ratio"]},
        {"": "GHB", "FoS 3D": g["fos3"], "FoS 2D": g["fos2"], "3D / 2D": g["ratio"]},
        {"": "EJC / GHB", "FoS 3D": res["ratio"]["ejc_ghb3"], "FoS 2D": res["ratio"]["ejc_ghb2"], "3D / 2D": None}]),
        hide_index=True, width="stretch",
        column_config={c: st.column_config.NumberColumn(format="%.6f") for c in ("FoS 3D", "FoS 2D")}
        | {"3D / 2D": st.column_config.NumberColumn(format="%.4f")})
    for lab, r in (("EJC", e), ("GHB", g)):
        for d in (3, 2):
            if r[f"fos{d}"] is None:
                st.warning(f"{lab} {d}D: {r[f'status{d}']}. {r[f'reason{d}']}")


def sweep_panel(inp: E.Inputs, d: dict):
    keys = study_keys(d)
    labels = {k: lab for k, _, lab in keys}
    eng = {k: n for k, n, _ in keys}
    st.write("The FoS as one input varies over evenly spaced values, for the strength models and dimensions ticked. "
             "The other inputs stay as they are.")
    c = st.columns([2, 1, 1, 1])
    k = c[0].selectbox("Input", list(labels), format_func=labels.get, key="sweep_name")
    lo0, hi0 = _range(k, float(getattr(inp, eng[k])))
    lo = c[1].number_input("from", value=lo0, key=f"sweep_lo_{k}", format="%.6g")
    hi = c[2].number_input("to", value=hi0, key=f"sweep_hi_{k}", format="%.6g")
    n = c[3].number_input("values", 2, 101, 9, key="sweep_n")
    c = st.columns(4)
    mc = inp.model == "mc"
    use = {"ejc_moscarda": c[0].checkbox("EJC", True, key="sweep_ejc", disabled=mc),
           "ghb": c[1].checkbox("GHB", True, key="sweep_ghb", disabled=mc)}
    dims = [dd for dd, col in ((3, c[2]), (2, c[3])) if col.checkbox(f"{dd}D", True, key=f"sweep_{dd}")]
    models = ["mc"] if mc else [m for m, on in use.items() if on]
    st.caption(SPEED)
    if st.button("Run the sweep", key="sweep_run"):
        sp = ui.NUM[k]
        if not (sp[2] <= lo < hi <= sp[3]):
            st.error(f"Use {sp[2]:g} ≤ from < to ≤ {sp[3]:g}.")
        elif not models or not dims:
            st.warning("Tick at least one model and one dimension.")
        else:
            values = [float(f"{v:.12g}") for v in np.linspace(lo, hi, int(n))]
            bar, progress = _bar("Solving…")
            rows = E.parameter_sweep(inp, eng[k], values, models, dims, progress)
            bar.empty()
            st.session_state.sweep = {"inputs": inp, "rows": rows, "label": labels[k], "models": models, "dims": dims}
    found = st.session_state.get("sweep")
    if not found:
        return
    _stale(found, inp, "Run the sweep")
    rows, series = found["rows"], [(m, dd) for m in found["models"] for dd in found["dims"]]
    fig = go.Figure()
    for m, dd in series:
        fig.add_trace(go.Scatter(x=[r["value"] for r in rows], y=[r[m][f"fos{dd}"] for r in rows], mode="lines+markers",
                                 name=f"{RP.MODEL_LAB[m]} {dd}D", line=dict(color=COLOR[dd], width=2, dash=DASH[m]),
                                 marker=dict(size=7, symbol="square" if m == "ghb" else "circle")))
    fig.add_hline(y=1.0, line=dict(color="#76746e", width=1, dash="dot"), annotation_text="F = 1",
                  annotation_position="bottom left")
    fig.update_layout(xaxis_title=found["label"], yaxis_title="FoS", height=380, margin=dict(l=10, r=10, t=30, b=10),
                      legend=dict(orientation="h", y=1.02, yanchor="bottom"))
    st.plotly_chart(fig, width="stretch", key="sweep_fig")
    st.dataframe(pd.DataFrame([{found["label"]: r["value"], **{f"{RP.MODEL_LAB[m]} {dd}D": r[m][f"fos{dd}"]
                                                                for m, dd in series}} for r in rows]),
                 hide_index=True, width="stretch")
    st.caption("Empty cells are runs with no accepted FoS.")


def back_panel(inp: E.Inputs, d: dict):
    keys = study_keys(d)
    labels = {k: lab for k, _, lab in keys}
    eng = {k: n for k, n, _ in keys}
    st.write("The value of one input at which the FoS equals a target (F = 1 for a failure), found by bisection between "
             "two values where F is above and below the target, on the prescribed surface.")
    c = st.columns([2, 1, 1, 1, 1])
    k = c[0].selectbox("Input", list(labels), format_func=labels.get, key="back_name")
    lo0, hi0 = _range(k, float(getattr(inp, eng[k])))
    lo = c[1].number_input("from", value=lo0, key=f"back_lo_{k}", format="%.6g")
    hi = c[2].number_input("to", value=hi0, key=f"back_hi_{k}", format="%.6g")
    target = c[3].number_input("target F", value=1.0, min_value=0.01, step=0.05, key="back_target")
    dim = c[4].selectbox("analysis", [2, 3], format_func=lambda v: f"{v}D", key="back_dim")
    st.caption(SPEED)
    if st.button("Back-analyse", key="back_run"):
        sp = ui.NUM[k]
        if not (sp[2] <= lo < hi <= sp[3]):
            st.error(f"Use {sp[2]:g} ≤ from < to ≤ {sp[3]:g}.")
        else:
            bar, progress = _bar("Solving…")
            res = E.back_analysis(inp, eng[k], lo, hi, target, dim, 1e-6, 60, progress)
            bar.empty()
            st.session_state.back = {"inputs": inp, "res": res, "label": labels[k]}
    found = st.session_state.get("back")
    if not found:
        return
    _stale(found, inp, "Back-analyse")
    res = found["res"]
    if res["value"] is not None:
        st.success(f"**{found['label']} = {res['value']:.6g}** gives F = {RP._f(res['fos'])} ({res['dim']}D), for the "
                   f"target {res['target']:g}; bracket {res['bracket'][0]:.8g} to {res['bracket'][1]:.8g}, "
                   f"{len(res['steps'])} analyses.")
    if res["message"]:
        st.warning(res["message"][:1].upper() + res["message"][1:] + ".")
    with st.expander("Every analysis"):
        st.dataframe(pd.DataFrame(res["steps"]).rename(columns={"value": found["label"], "fos": f"FoS {res['dim']}D"}),
                     hide_index=True, width="stretch")


def mc_panel(inp: E.Inputs, d: dict):
    keys = [kk for kk in study_keys(d) if kk[0] not in GEOMETRY]
    st.write("Each input below is drawn from a normal distribution with the current value as its mean and the standard "
             "deviation you enter (0 keeps it fixed), cut off at its valid range. No standard deviation is assumed: "
             "enter them from your site data. The FoS of every sample is solved with the current model; "
             "P_f = P(F < 1) among the accepted results.")
    variables, labels = [], []
    cols = st.columns(3)
    for q, (k, name, lab) in enumerate(keys):
        mean = float(getattr(inp, name))
        lo, hi = ui.NUM[k][2], ui.NUM[k][3]
        sd = cols[q % 3].number_input(f"SD of {lab} (mean {mean:g}; {lo:g} to {hi:g})", min_value=0.0, value=0.0,
                                      key=f"mc_sd_{k}", format="%.6g")
        if sd > 0:
            variables.append((name, mean, float(sd), float(lo), float(hi)))
            labels.append(lab)
    c = st.columns(3)
    n = c[0].number_input("samples", 10, 5000, 100, key="mc_n")
    seed = c[1].number_input("seed", 0, 2**31 - 1, 1, key="mc_seed")
    dsel = c[2].selectbox("analysis", ["2D", "3D", "3D and 2D"], key="mc_dim")
    dims = {"2D": [2], "3D": [3], "3D and 2D": [3, 2]}[dsel]
    st.caption(SPEED)
    if st.button("Run Monte Carlo", key="mc_run"):
        if not variables:
            st.warning("Enter a standard deviation above 0 for at least one input.")
        else:
            bar, progress = _bar("Sampling…")
            res = E.monte_carlo(inp, variables, int(n), int(seed), None, dims, progress)
            bar.empty()
            st.session_state.mc = {"inputs": inp, "res": res, "labels": labels}
    found = st.session_state.get("mc")
    if not found:
        return
    _stale(found, inp, "Run Monte Carlo")
    res = found["res"]
    rows = []
    for m in res["models"]:
        for dd in res["dims"]:
            s = res["stats"][f"{m}{dd}"]
            if s["pof"] is None:
                pf = "–"
            elif s["pof"] > 0:
                pf = f"{100 * s['pof']:.1f} % ± {100 * s['pof_se']:.1f} % ({round(s['pof'] * s['accepted'])} failed)"
            else:
                pf = f"0 of {s['accepted']} failed: below {300 / s['accepted']:.2g} % (95 %)"
            rows.append({"": f"{RP.MODEL_LAB.get(m, m)} {dd}D", "Accepted": f"{s['accepted']} of {s['n']}",
                         "Mean F": s["mean"], "SD": s["sd"], "β = (mean − 1)/SD": s["beta"], "P_f = P(F < 1)": pf,
                         "5 % / 50 % / 95 %": " / ".join("–" if s[q] is None else f"{s[q]:.4f}" for q in ("p05", "p50", "p95"))})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch",
                 column_config={c: st.column_config.NumberColumn(format="%.4f") for c in ("Mean F", "SD", "β = (mean − 1)/SD")})
    st.caption(f"{len(res['samples'])} samples, seed {res['seed']}; varied: "
               + ", ".join(f"{found['labels'][q]} (mean {v[1]:g}, SD {v[2]:g})" for q, v in enumerate(res["variables"]))
               + ". The ± is one standard error of P_f, √(P_f(1 − P_f)/n). When no sample fails, the 95 % upper bound on "
               "P_f is 3/n (the rule of three). Withheld results (no accepted FoS) are left out and counted above.")
    fig = go.Figure()
    for m in res["models"]:
        for dd in res["dims"]:
            fig.add_trace(go.Histogram(x=[s[f"{m}{dd}"] for s in res["samples"] if s[f"{m}{dd}"] is not None],
                                       name=f"{RP.MODEL_LAB.get(m, m)} {dd}D", opacity=0.6,
                                       marker=dict(color=COLOR[dd], line=dict(color="white", width=1))))
    fig.add_vline(x=1.0, line=dict(color="#76746e", width=1, dash="dot"), annotation_text="F = 1")
    fig.update_layout(barmode="overlay", xaxis_title="FoS", yaxis_title="samples", height=340,
                      margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.02, yanchor="bottom"))
    st.plotly_chart(fig, width="stretch", key="mc_fig")


def report_studies(A: E.Analysis) -> list[tuple[str, str]]:
    """The studies run on A's inputs, as (heading, html) for the report."""
    out, ss = [], st.session_state
    if ss.get("cmp") and ss.cmp["inputs"] == A.inputs:
        out.append(("EJC against GHB", RP.compare_table(ss.cmp["res"], A.inputs.nu)))
    if ss.get("back") and ss.back["inputs"] == A.inputs:
        out.append(("Back-analysis", RP.back_text(ss.back["res"], ss.back["label"])))
    if ss.get("sweep") and ss.sweep["inputs"] == A.inputs:
        f = ss.sweep
        out.append((f"Parameter sweep: {f['label']}", RP.sweep_table(f["rows"], f["label"], f["models"], f["dims"])))
    if ss.get("mc") and ss.mc["inputs"] == A.inputs:
        out.append(("Monte Carlo", RP.mc_table(ss.mc["res"], ss.mc["labels"])))
    return out
