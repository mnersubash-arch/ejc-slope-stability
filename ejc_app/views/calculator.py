"""Page 1: Calculator."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import engine as E
import plots as P
import ui

st.title("Calculator")
inp = ui.sidebar()
ui.preset_box()


@st.cache_data(show_spinner=False)
def original_check(sigma_ci, mb, s, a, nu):
    return E.original_closure_minimum(sigma_ci, mb, s, a, nu)


def original_closure_panel(inp: E.Inputs):
    """Brief §1.2(c): the original closure never produces a FoS; compare its minimum admissible sigma_n with the
    gravity cap sigma_n <= gamma h - u."""
    st.subheader("Original closure: admissibility check (no FoS)")
    mb, s, a = inp.hb
    with st.spinner("Searching the admissible states..."):
        r = original_check(inp.sigma_ci, mb, s, a, inp.nu)
    try:
        md, mw = E.build_mesh(inp.with_(wet=False), 3, inp.n3), E.build_mesh(inp.with_(wet=True), 3, inp.n3)
    except ValueError as exc:
        st.error(str(exc))
        return
    cap_d, cap_w = E.gravity_cap(md, inp.gamma), E.gravity_cap(mw, inp.gamma)
    if not r.found:
        st.error("No admissible state of the original closure in sigma3 = 1e-6 to 1e6 MPa.")
        return
    c = st.columns(4)
    c[0].metric("σ3 at the minimum (MPa)", f"{r.sigma3_min:.3f}")
    c[1].metric("Minimum admissible σ_n (MPa)", f"{r.sn_min:.3f}")
    c[2].metric("Gravity cap γh − u, dry (MPa)", f"{cap_d.max() / 1000:.3f}")
    c[3].metric("Gravity cap γh − u, wet (MPa)", f"{cap_w.max() / 1000:.3f}")
    nd, nw = int(np.sum(cap_d < r.sn_min * 1000)), int(np.sum(cap_w < r.sn_min * 1000))
    st.write(f"At ν = {inp.nu:g}, a = {a:g}: **{nd} of {md.count}** dry columns and **{nw} of {mw.count}** wet "
             f"columns ({inp.n3}x{inp.n3} half-model) have γh − u below the minimum admissible σ_n "
             f"(φ_i there = {r.phi_deg:.3f}°). Those bases have no admissible state, so this closure gives no FoS.")
    verdict = ("The minimum admissible σ_n is above the gravity cap of every column."
               if r.sn_min * 1000 > max(cap_d.max(), cap_w.max()) else
               "Some columns lie above the minimum, but every base must have an admissible state for a FoS.")
    st.info(verdict + " Ordering requires σ3/(2ν) ≤ σ1 ≤ σ3(1 − ν)/ν, and the EJC requires σ1 − σ3 ≥ σ_ci Q^a.")


if inp.model == "ejc_original":
    original_closure_panel(inp)
    st.stop()

A, stale = ui.analysis(inp)
if A is None:
    st.stop()
r3, r2 = A.r3, A.r2

# ---- headline numbers
c = st.columns(6)
c[0].metric("FoS 3D" + (" (F_M)" if A.inputs.form == "FM" else ""), ui.fmt(r3.fos))
c[1].metric("FoS 2D (Bishop)", ui.fmt(r2.fos))
c[2].metric("3D / 2D ratio", ui.fmt(A.ratio, 4))
c[3].metric("Iterations 3D / 2D", f"{r3.iterations} / {r2.iterations}")
c[4].metric("Active columns / slices", f"{r3.mesh.count} / {r2.mesh.count}")
c[5].metric("Projections 3D / 2D", f"{r3.projections} / {r2.projections}",
            help="Trial σ_n moved above the projection limit during the iteration (thesis Section 4.6). A run is "
                 "accepted only if no base is projected at convergence.")
ui.status_badge(r3, "3D")
ui.status_badge(r2, "2D")
if A.inputs.tension != "none":
    st.caption(f"{A.inputs.tension}: sensitivity, not a calibrated contact law.")
mb, s, a = A.inputs.hb
st.caption(f"Strength: {E.MODELS[A.inputs.model]}. m_b = {mb:.6f}, s = {s:.7f}, a = {a:.6f}, ν = {A.inputs.nu:g}. "
           f"Iteration path: {E.ITERATIONS[A.inputs.iteration]}. Sliding volume (full mass) "
           f"{r3.mesh.volume:,.0f} m³; 2D area {r2.mesh.volume:,.1f} m²."
           + (f" Projection limit {r3.projection_limit:.4f} kPa." if r3.projection_limit is not None else ""))

# ---- checks
st.subheader("Checks")
c1, c2 = st.columns(2)
for col, r, lab in ((c1, r3, "3D Hungr columns"), (c2, r2, "2D Bishop slices")):
    with col:
        st.markdown(f"**{lab}**")
        if r.checks:
            st.dataframe(ui.checks_frame(r).style.map(ui.color_status, subset=["Status"]), hide_index=True,
                         width="stretch")
        else:
            st.write("Final checks not reached: " + r.reason)

with st.expander("Iteration history"):
    c1, c2 = st.columns(2)
    c1.plotly_chart(P.history(r3), width="stretch", key="h3")
    c2.plotly_chart(P.history(r2), width="stretch", key="h2")

# ---- downloads
st.subheader("Column and slice results")
t3, t2 = pd.DataFrame(E.columns_table(r3)), pd.DataFrame(E.columns_table(r2))
c1, c2 = st.columns(2)
c1.download_button("Download 3D column results (CSV)", t3.to_csv(index=False).encode(), "columns_3d.csv",
                   "text/csv", width="stretch")
c2.download_button("Download 2D slice results (CSV)", t2.to_csv(index=False).encode(), "slices_2d.csv",
                   "text/csv", width="stretch")
with st.expander("Show the 3D column table"):
    st.dataframe(t3, width="stretch", height=320)

# ---- studies
st.subheader("Studies")
tab_mesh, tab_nu = st.tabs(["Mesh study (eq. 4.17)", "ν sweep (Moscarda closure)"])
with tab_mesh:
    base = A.inputs
    c1, c2 = st.columns(2)
    m3 = c1.multiselect("3D meshes N", [10, 15, 20, 30, 40, 60, 80, 120, 160, 240, 320], default=[10, 20, 40, 80])
    m2 = c2.multiselect("2D slices N_2D", [20, 40, 80, 160, 320, 640, 1280, 2560, 5120], default=[20, 40, 80, 160, 320])
    if st.button("Run mesh study"):
        rows = []
        bar = st.progress(0.0)
        jobs = [(3, n) for n in sorted(m3)] + [(2, n) for n in sorted(m2)]
        for q, (dim, n) in enumerate(jobs):
            r = ui.cached_run_3d(base.with_(n3=n)) if dim == 3 else ui.cached_run_2d(base.with_(n2=n))
            rows.append(dict(dim=f"{dim}D", mesh=f"{n}x{n}" if dim == 3 else f"{n} slices", n=n, elements=r.mesh.count,
                             status=r.status, FoS=r.fos, iterations=r.iterations, projections=r.projections))
            bar.progress((q + 1) / len(jobs))
        st.session_state.mesh_study = rows
    if st.session_state.get("mesh_study"):
        df = pd.DataFrame(st.session_state.mesh_study)
        st.dataframe(df, hide_index=True, width="stretch")
        for dim in ("3D", "2D"):
            d = df[(df.dim == dim) & df.FoS.notna()]
            if len(d) >= 2:
                eps = E.mesh_error(list(d.FoS), float(d.sort_values("n").FoS.iloc[-1]))
                ok = "meets" if eps <= 0.5 else "does not meet"
                st.write(f"{dim}: ε_mesh = (F_max − F_min)/F_finest × 100 = **{eps:.3f} %** over "
                         f"{', '.join(d.mesh)}; {ok} the 0.5 % screening threshold. Convergence is not an error bound.")
with tab_nu:
    if st.session_state.pop("open_nu_sweep", False):
        st.info("ν sweep preset loaded: set the range and press Run ν sweep.")
    c1, c2, c3 = st.columns(3)
    nu0 = c1.number_input("ν from", 0.01, 0.49, 0.05, 0.01, format="%.3f")
    nu1 = c2.number_input("ν to", 0.01, 0.49, 0.49, 0.01, format="%.3f")
    dnu = c3.number_input("step", 0.005, 0.2, 0.02, 0.005, format="%.3f")
    if st.button("Run ν sweep"):
        nus = np.round(np.arange(nu0, nu1 + 1e-9, dnu), 6)
        rows = []
        bar = st.progress(0.0)
        for q, nu in enumerate(nus):
            b = A.inputs.with_(model="ejc_moscarda", nu=float(nu))
            r3n, r2n = ui.cached_run_3d(b), ui.cached_run_2d(b)
            rows.append(dict(nu=float(nu), FoS_3D=r3n.fos, FoS_2D=r2n.fos, status_3D=r3n.status, status_2D=r2n.status,
                             reason_3D="" if r3n.fos else r3n.reason[:120]))
            bar.progress((q + 1) / len(nus))
        st.session_state.nu_sweep = rows
    if st.session_state.get("nu_sweep"):
        df = pd.DataFrame(st.session_state.nu_sweep)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=df.nu, y=df.FoS_3D, mode="lines+markers", name="FoS 3D"))
        fig.add_trace(go.Scatter(x=df.nu, y=df.FoS_2D, mode="lines+markers", name="FoS 2D"))
        fig.update_layout(xaxis_title="ν", yaxis_title="FoS", height=360, margin=dict(l=10, r=10, t=10, b=10))
        st.plotly_chart(fig, width="stretch", key="nu_fig")
        st.caption("Gaps are runs with no accepted FoS (reason in the table).")
        st.dataframe(df, hide_index=True, width="stretch")
