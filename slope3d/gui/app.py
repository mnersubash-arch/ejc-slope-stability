"""slope3d Streamlit app: input panel, presets, status panel, results tabs, plots and exports.

Rule enforced everywhere on screen: a factor of safety is shown only for an ACCEPTED dimension. A NOT ACCEPTED run says
"FoS withheld"; a DIAGNOSTIC ONLY run shows its F on a grey card labelled "Diagnostic F (unchecked roots, not an accepted result)". Closure 3 is always labelled as the conventional-HB comparator.
"""
from __future__ import annotations

import io
import math
import tempfile
from pathlib import Path

import streamlit as st

from slope3d import Config, __version__, run_case
from slope3d import io as s3io
from slope3d import plotly_figs as pf
from slope3d.geometry import slip_surface_ends
from slope3d.models import CLOSURE_NAMES, DIAGNOSTIC_F_LABEL, MESH2D_RANGE, MESH3D_RANGE, Status
from slope3d.plotdata import mesh_convergence
from slope3d.presets import PRESETS

STRENGTH_LABELS = {1: "Hoek-Brown / EJC eq. 4.21", 0: "Mohr-Coulomb (verification model)"}
CLOSURE_LABELS = {1: "1 - Thesis closure (sigma2' = sigma3'/nu - sigma1')",
                  2: "2 - Moscarda closure (sigma2' = nu (sigma1' + sigma3'))",
                  3: "3 - Conventional HB comparator (3D geometry, 2D strength)"}
STATUS_STYLE = {Status.ACCEPTED: ("#008300", "ACCEPTED"), Status.NOT_ACCEPTED: ("#c0392b", "NOT ACCEPTED"),
                Status.DIAGNOSTIC_ONLY: ("#6b6a66", "DIAGNOSTIC ONLY")}
NUM_FIELDS = ["xc", "zc", "radius", "height", "gamma", "water_angle", "water_height", "u_extra", "c_mc", "phi_mc",
              "sci", "mi", "gsi", "d", "nu", "tol_a", "tol_e", "tol_f"]


# --------------------------------------------------------------------------------------------------
# state <-> Config
# --------------------------------------------------------------------------------------------------
def load_into_state(cfg: Config, mesh3d: int, mesh2d: int) -> None:
    ss = st.session_state
    for f in NUM_FIELDS:
        ss[f"in_{f}"] = float(getattr(cfg, f))
    ss["in_face"] = float(cfg.face_angle_deg)
    ss["_face_src"], ss["_slope_src"] = ss["in_face"], cfg.slope
    ss["in_a_override"] = cfg.a_override is not None
    ss["in_a"] = float(cfg.a_override if cfg.a_override is not None else cfg.a_gsi)
    ss["in_strength"] = cfg.strength
    ss["in_closure"] = cfg.closure
    ss["in_checks"] = cfg.checks
    ss["in_max_iter"] = int(cfg.max_iter)
    ss["in_mesh3d"], ss["in_mesh2d"] = int(mesh3d), int(mesh2d)


def config_from_state() -> tuple[Config, int, int]:
    ss = st.session_state
    kw = {f: ss[f"in_{f}"] for f in NUM_FIELDS}
    face = ss["in_face"]
    slope = ss["_slope_src"] if face == ss["_face_src"] else math.tan(math.radians(face))
    cfg = Config(**kw, slope=slope, a_override=ss["in_a"] if ss["in_a_override"] else None, strength=ss["in_strength"],
                 closure=ss["in_closure"], checks=ss["in_checks"], max_iter=int(ss["in_max_iter"]))
    return cfg, int(ss["in_mesh3d"]), int(ss["in_mesh2d"])


def _apply_preset(name: str) -> None:
    p = PRESETS[name]
    load_into_state(p.config, p.mesh3d, p.mesh2d)
    st.session_state["_run_now"] = True
    st.session_state["_preset"] = name


# --------------------------------------------------------------------------------------------------
# widgets
# --------------------------------------------------------------------------------------------------
def _num(label: str, key: str, help_: str = "", fmt: str = "%.10g", **kw) -> None:
    st.number_input(label, key=key, help=help_, format=fmt, **kw)


def input_panel() -> None:
    st.sidebar.markdown("### Presets")
    names = list(PRESETS)
    for i, name in enumerate(names):
        st.sidebar.button(name, key=f"preset_{i}", on_click=_apply_preset, args=(name,), width="stretch")
    st.sidebar.caption("Presets load the inputs and run the case (3D 20x20, 2D 40 slices).")
    st.sidebar.markdown("### Inputs")
    with st.sidebar.expander("Geometry", expanded=False):
        _num("Centre x_c [m]", "in_xc", "Slip-sphere centre; 0 < x_c < R. The sphere centre y is fixed at 0 (symmetry).")
        _num("Centre z_c [m]", "in_zc", "R > z_c > 0")
        _num("Radius R [m]", "in_radius")
        _num("Unit weight gamma [kN/m3]", "in_gamma")
    with st.sidebar.expander("Slope"):
        _num("Slope height H [m]", "in_height")
        _num("Face angle [deg]", "in_face", "Entered as an angle; the engine uses slope = tan(angle).")
        st.caption(f"slope = tan(face angle) = {config_from_state()[0].slope:.12g}")
    with st.sidebar.expander("Water"):
        _num("Water table angle [deg]", "in_water_angle")
        _num("Water table height [m]", "in_water_height", "0 = dry. Wet benchmark cases: 50 deg and 360 m.")
        _num("Extra uniform pore pressure [kPa]", "in_u_extra")
    with st.sidebar.expander("Rock mass"):
        _num("sigma_ci [MPa]", "in_sci")
        _num("m_i", "in_mi")
        _num("GSI", "in_gsi")
        _num("Disturbance D", "in_d")
        st.checkbox("Override a (otherwise a from GSI)", key="in_a_override")
        _num("a", "in_a", "The thesis uses 0.5 or 0.509923.", disabled=not st.session_state["in_a_override"])
        c = config_from_state()[0]
        st.caption(f"m_b = {c.mb:.6f}, s = {c.s:.6g}, a from GSI = {c.a_gsi:.6f}, a used = {c.a:.6f}")
    with st.sidebar.expander("Strength model"):
        st.selectbox("Strength model", options=list(STRENGTH_LABELS), format_func=STRENGTH_LABELS.get, key="in_strength")
        _num("Mohr-Coulomb c [kPa]", "in_c_mc", disabled=st.session_state["in_strength"] == 1)
        _num("Mohr-Coulomb phi [deg]", "in_phi_mc", disabled=st.session_state["in_strength"] == 1)
    with st.sidebar.expander("Closure"):
        st.selectbox("sigma2' closure", options=list(CLOSURE_LABELS), format_func=CLOSURE_LABELS.get, key="in_closure",
                     disabled=st.session_state["in_strength"] == 0)
        _num("Poisson's ratio nu", "in_nu", "0 < nu < 0.5 (used by closures 1 and 2).")
    with st.sidebar.expander("Tolerances and run controls"):
        _num("Tangency tolerance [rad]", "in_tol_a", "<= 1e-10 (can be tightened, not loosened)", fmt="%.3e")
        _num("EJC / closure tolerance [MPa]", "in_tol_e", "<= 1e-8", fmt="%.3e")
        _num("F tolerance (relative)", "in_tol_f", "<= 1e-8", fmt="%.3e")
        st.radio("Acceptance checks", [1, 0], key="in_checks", horizontal=True,
                 format_func=lambda v: "1 - on (default)" if v == 1 else "0 - off: DIAGNOSTIC ONLY, F withheld")
        if st.session_state["in_checks"] == 0:
            st.warning("Diagnostic mode: the first converged eq. 4.21 root of every element is used without an admissibility "
                       "test. The result can never be ACCEPTED; its F is shown only as a labelled diagnostic value.")
        st.number_input("Maximum iterations", 1, 50000, key="in_max_iter", step=100)
    with st.sidebar.expander("Mesh"):
        st.number_input(f"3D mesh n (n x n half-model, {MESH3D_RANGE[0]}-{MESH3D_RANGE[1]})", MESH3D_RANGE[0], MESH3D_RANGE[1],
                        key="in_mesh3d", step=5)
        st.number_input(f"2D slices ({MESH2D_RANGE[0]}-{MESH2D_RANGE[1]})", MESH2D_RANGE[0], MESH2D_RANGE[1], key="in_mesh2d", step=10)
    st.sidebar.button("Run", type="primary", key="run_btn", width="stretch")


# --------------------------------------------------------------------------------------------------
# status panel and tabs
# --------------------------------------------------------------------------------------------------
def status_card(tag: str, run, mesh: int) -> None:
    r = run.result
    colour, text = STATUS_STYLE[r.status]
    if r.fos is not None:
        fos = f"<div style='font-size:1.4rem;font-weight:700'>F = {r.fos:.4f}</div>"
    elif r.diagnostic_f is not None:      # grey card, never the ACCEPTED style, never a bare "F ="
        fos = (f"<div style='font-size:0.85rem;font-weight:600;opacity:.8'>{DIAGNOSTIC_F_LABEL}</div>"
               f"<div style='font-size:1.3rem;font-weight:600;color:{colour}'>{r.diagnostic_f:.4f}</div>")
    else:
        fos = "<div style='font-size:1.4rem;font-weight:700'>FoS withheld</div>"
    bg = "background:rgba(128,128,128,0.14);" if r.status == Status.DIAGNOSTIC_ONLY else ""
    st.markdown(
        f"<div style='border:2px solid {colour};{bg}border-radius:8px;padding:12px 16px'>"
        f"<div style='font-size:0.9rem;opacity:.75'>{tag} &middot; {mesh}{'x' + str(mesh) if tag == '3D' else ''} mesh &middot; "
        f"{r.active} active {'columns' if tag == '3D' else 'slices'}</div>"
        f"<div style='font-size:2rem;font-weight:800;color:{colour};line-height:1.2'>{text}</div>"
        f"{fos}"
        f"<div style='font-size:0.85rem;opacity:.85;margin-top:4px'>{r.reason}</div></div>", unsafe_allow_html=True)


def identity_banner(case) -> None:
    c = case.config
    if c.strength == 0:
        st.info("Mohr-Coulomb verification model (numerical check, not a calibrated material).")
    elif c.is_thesis_closure:
        st.success("Thesis closure (3D EJC): a result here is a thesis-closure result only if that dimension is ACCEPTED.")
    elif c.closure == 3:
        st.warning("Conventional HB comparator (3D geometry, 2D strength) - not a 3D-EJC result.")
    elif c.checks == 0:
        st.warning(f"{CLOSURE_NAMES[c.closure]} - acceptance checks OFF: DIAGNOSTIC ONLY. The diagnostic F is not a factor of safety.")
    else:
        st.info(CLOSURE_NAMES[c.closure])


def tab_summary(case) -> None:
    r3, r2 = case.run3d.result, case.run2d.result
    a, b, c, d = st.columns(4)
    for col, tag, r in ((a, "3D", r3), (b, "2D", r2)):
        if r.diagnostic_f is not None:
            col.metric(f"Diagnostic F {tag} (unchecked, not accepted)", f"{r.diagnostic_f:.4f}")
        else:
            col.metric(f"F {tag}", f"{r.fos:.4f}" if r.fos is not None else "WITHHELD")
    c.metric("3D / 2D ratio", f"{case.ratio:.4f}" if case.ratio is not None else "n.a.")
    d.metric("Iterations 3D / 2D", f"{r3.iterations} / {r2.iterations}")
    ends = slip_surface_ends(case.config)
    if ends.get("supported"):
        st.markdown(f"**Slip surface:** enters at x = {ends['entry_x']:.2f} m (z = {ends['entry_z']:.2f} m); " +
                    (f"exits the slope face at x = {ends['exit_x']:.2f} m, **{ends['exit_z']:.2f} m above the toe** "
                     f"({ends['exit_behind_toe']:.2f} m behind the toe)." if ends["exit_on_face"] else
                     f"reaches the toe plane (x = 0) at z = {ends['exit_z']:.3f} m."))
    for tag, r in (("3D columns", r3), ("2D slices", r2)):
        if r.evaluated:
            st.markdown(f"**{tag}:** {r.tension} element(s) with sigma_n' < 0; minimum sigma_n' = {r.min_normal:.2f} kPa.")
        else:
            st.markdown(f"**{tag}:** final stresses not evaluated.")
    if r3.moment_diagnostic is not None:
        st.caption(f"Spherical moment-form ratio (diagnostic only, not a separate equilibrium solution): {r3.moment_diagnostic:.6f}")
    st.dataframe(s3io.summary_table(case).astype(str), hide_index=True, width="stretch")


def tab_elements(case) -> None:
    which = st.radio("Table", ["3D columns", "2D slices"], horizontal=True)
    run = case.run3d if which.startswith("3D") else case.run2d
    df = s3io.elements_table(run, case.config)
    st.caption(s3io.ELEMENT_UNITS)
    st.caption(f"{run.result.status.value}: {run.result.reason}")
    if df.empty:
        st.info("No elements for this run.")
    else:
        st.dataframe(df.astype(str) if (df.dtypes == object).any() else df, hide_index=True, width="stretch", height=520)


def tab_plots(case) -> None:
    st.plotly_chart(pf.section(case), width="stretch")
    c1, c2 = st.columns(2)
    c1.plotly_chart(pf.plan(case), width="stretch")
    c2.plotly_chart(pf.surface3d(case), width="stretch")
    c1, c2 = st.columns(2)
    c1.plotly_chart(pf.convergence(case), width="stretch")
    c2.plotly_chart(pf.strength(case), width="stretch")
    heavy = case.config.strength == 1 and case.config.closure != 3
    if st.button("Run mesh convergence", help="3D 10-80, 2D 20-640. " + ("Slow for closures 1 and 2 (full root scan)." if heavy else "")):
        with st.spinner("Running meshes..."):
            st.session_state["mesh_rows"] = mesh_convergence(case.config)
            st.session_state["mesh_for"] = repr(case.config)
    rows = st.session_state.get("mesh_rows") if st.session_state.get("mesh_for") == repr(case.config) else None
    st.plotly_chart(pf.mesh(rows), width="stretch")


def tab_exports(case) -> None:
    import json
    stem = "slope3d"
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Tables**")
        st.download_button("Summary (CSV)", s3io.summary_table(case).to_csv(index=False), f"{stem}_summary.csv", "text/csv")
        st.download_button("3D columns (CSV)", s3io.elements_table(case.run3d, case.config).to_csv(index=False), f"{stem}_columns3d.csv", "text/csv")
        st.download_button("2D slices (CSV)", s3io.elements_table(case.run2d, case.config).to_csv(index=False), f"{stem}_slices2d.csv", "text/csv")
        with tempfile.TemporaryDirectory() as td:
            p = s3io.export_xlsx(case, Path(td) / f"{stem}.xlsx")
            st.download_button("Everything (XLSX)", p.read_bytes(), f"{stem}.xlsx",
                               "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        st.download_button("Result summary (JSON)", json.dumps(s3io.result_to_dict(case), indent=2), f"{stem}.json", "application/json")
        with tempfile.TemporaryDirectory() as td:
            s3io.save_case(Path(td) / "case.yaml", case.config, case.mesh3d, case.mesh2d)
            st.download_button("Case file (YAML)", (Path(td) / "case.yaml").read_text(), "case.yaml", "text/yaml")
    with c2:
        st.markdown("**Run report** (all six plots, summary, inputs, slip-surface exit and tension statistics)")
        if st.button("Build report (PDF + HTML)"):
            from slope3d.report import write_report
            with st.spinner("Rendering plots..."), tempfile.TemporaryDirectory() as td:
                rows = st.session_state.get("mesh_rows") if st.session_state.get("mesh_for") == repr(case.config) else None
                paths = write_report(case, td, stem, mesh_rows=rows)
                st.session_state["report_files"] = {p.name: p.read_bytes() for p in paths}
        for name, data in st.session_state.get("report_files", {}).items():
            st.download_button(f"Download {name}", data, name, "application/pdf" if name.endswith("pdf") else "text/html")


# --------------------------------------------------------------------------------------------------
def main() -> None:
    st.set_page_config(page_title="slope3d", layout="wide")
    if "in_xc" not in st.session_state:
        load_into_state(Config(), 20, 40)
    st.title("slope3d")
    st.caption(f"3D limit-equilibrium rock-slope stability with a 3D Hoek-Brown (EJC eq. 4.21) strength model - v{__version__}. "
               "Prescribed spherical slip surface; no search for the critical surface.")
    input_panel()
    if st.session_state.pop("_run_now", False) or st.session_state.get("run_btn") or "case" not in st.session_state:
        cfg, m3, m2 = config_from_state()
        with st.spinner("Solving 3D columns and 2D slices..."):
            st.session_state["case"] = run_case(cfg, m3, m2)
            st.session_state.pop("report_files", None)
    case = st.session_state["case"]
    if st.session_state.get("_preset"):
        st.caption(f"Preset: {st.session_state['_preset']} - {PRESETS[st.session_state['_preset']].note}")
    identity_banner(case)
    c1, c2 = st.columns(2)
    with c1:
        status_card("3D", case.run3d, case.mesh3d)
    with c2:
        status_card("2D", case.run2d, case.mesh2d)
    t1, t2, t3, t4 = st.tabs(["Summary", "Element table", "Plots", "Exports"])
    with t1:
        tab_summary(case)
    with t2:
        tab_elements(case)
    with t3:
        tab_plots(case)
    with t4:
        tab_exports(case)


main()
