"""Page 1: Calculator."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import engine as E
import plots as P
import ui

st.markdown(f"**{ui.DISCLAIMER}**")
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


def search_panel(inp: E.Inputs):
    """Critical-surface search (engine.critical_search; QA/QC finding R1), as the web app's Studies panel."""
    st.write("Finds the circle (2D) and the sphere centred on y = 0 (3D) with the lowest accepted FoS for the current "
             "inputs, mesh and grid option. The search tries surfaces that leave the ground at a front exit (the toe or "
             "the face; with the base-circle option also the pit floor) and at a back exit (the face or the crest), and "
             "that meet the ground at both exits, so the mass has no vertical back face.")
    c1, c2, c3 = st.columns(3)
    back = c1.number_input("Back exit up to (× H behind the crest)", 0.25, 10.0, E.SEARCH["back"], 0.25, key="srch_back")
    front = c2.number_input("Front exit up to (× H ahead of the toe)", 0.25, 10.0, E.SEARCH["front"], 0.25,
                            key="srch_front", disabled=not inp.base_circle,
                            help="Base-circle option only: with the default grid the front exit is the toe or the face.")
    with3d = c3.checkbox("Include 3D (slower: it starts from the 2D minima)", key="srch_3d")
    st.caption("This engine takes about 0.1 to 1 s per surface (Hoek-Brown and EJC; Mohr-Coulomb is faster), and a 2D "
               "search tries about 500 to 700 surfaces, so it takes minutes; 3D adds about 300 more. Keep this page "
               "open while it runs. The web app (webapp/) runs the same search, with the same results, about 15 times "
               "faster.")
    if st.button("Search", type="primary", key="srch_run"):
        bar = st.progress(0.0, text="Searching…")

        def progress(done, total, label):
            bar.progress(min(done / total, 1.0), text=f"{label}… ({done} of about {total} surfaces)")

        res = E.critical_search(inp, dims=(2, 3) if with3d else (2,), progress=progress, back=back, front=front)
        bar.empty()
        st.session_state.search = {"inputs": inp, "res": res}
    found = st.session_state.get("search")
    if not found:
        return
    res, searched = found["res"], found["inputs"]
    dims = res["settings"]["dims"]
    st.write(f"Finished in {res['ms'] / 1000:.1f} s: "
             + ", ".join(f"{res[f'r{d}']['evaluations']} surfaces in {d}D" for d in dims) + ".")
    st.dataframe(pd.DataFrame(P.search_rows(searched, res)), hide_index=True, width="stretch")
    cols = st.columns(len(dims))
    for col, d in zip(cols, dims):
        r = res[f"r{d}"]
        if r["fos"] is not None:
            col.button(f"Use this {'circle' if d == 2 else 'sphere'} ({d}D)", key=f"srch_use{d}", on_click=ui.use_surface,
                       args=(r["xc"], r["zc"], r["R"]), disabled=not ui.surface_in_limits(r), width="stretch")
    pr = res["prescribed"]
    if len(dims) == 2 and all(res[f"r{d}"]["fos"] is not None for d in dims) and pr.get("fos2") and pr.get("fos3"):
        st.write(f"3D/2D ratio: **{pr['fos3'] / pr['fos2']:.3f}** on the prescribed surface, "
                 f"**{res['r3']['fos'] / res['r2']['fos']:.3f}** on the critical surfaces (the 2D circle and the 3D "
                 "sphere each with its own lowest FoS).")
    for kind, text in ui.search_notes(inp, searched, res):
        (st.warning if kind == "warning" else st.info)(text)
    st.plotly_chart(P.search_figure(searched, res), width="stretch", key="search_fig")
    s = res["settings"]
    n_grid = s["grid"][0] * s["grid"][1] * s["grid"][2]
    st.caption(f"2D: a grid of {' × '.join(str(v) for v in s['grid'])} = {n_grid} starting circles (front exit × back "
               f"exit × curvature), then Nelder–Mead from the {s['starts']} best that are not neighbours, each restarted "
               "until it stops improving." + (" 3D: Nelder–Mead from the 2D minima." if 3 in dims else "")
               + " The critical surface is the lowest accepted FoS among the surfaces tried: a local search, not a proof "
               "of the global minimum. Withheld results and surfaces the grid option cannot analyse are skipped. Every "
               f"surface is solved on the current mesh ({searched.n2} slices, {searched.n3} × {searched.n3} columns); "
               "the mesh check re-solves the critical surface on a finer one. On a coarse 3D mesh the FoS changes in "
               "small steps as columns enter or leave the mass, so searches from different starts can stop at FoS that "
               "differ by a few tenths of a percent.")


if inp.model == "ejc_original":
    original_closure_panel(inp)
    st.stop()

A, stale = ui.analysis(inp)
if A is None:
    st.subheader("Critical surface (search)")
    search_panel(inp)
    st.stop()
r3, r2 = A.r3, A.r2

# ---- headline numbers
c = st.columns(7)
c[0].metric("FoS 3D" + (" (F_M)" if A.inputs.form == "FM" else ""), ui.fmt(r3.fos))
c[1].metric("FoS 2D (Bishop)", ui.fmt(r2.fos))
c[2].metric("3D / 2D ratio", ui.fmt(A.ratio, 4))
c[3].metric("Iterations 3D / 2D", f"{r3.iterations} / {r2.iterations}")
c[4].metric("Active columns / slices", f"{r3.mesh.count} / {r2.mesh.count}")
c[5].metric("Projections 3D / 2D", f"{r3.projections} / {r2.projections}",
            help="Trial σ_n moved above the projection limit during the iteration (thesis Section 4.6). A run is "
                 "accepted only if no base is projected at convergence.")
c[6].metric("Tensile bases 3D / 2D", f"{ui.tensile_text(ui.tensile_count(r3))} / {ui.tensile_text(ui.tensile_count(r2))}",
            help="Bases with tensile effective normal stress σ′_n < 0 at the accepted state (thesis Section 5.8).")
ui.status_badge(r3, "3D")
ui.status_badge(r2, "2D")
if r3.mesh.count or r2.mesh.count:
    vf = ui.vertical_face_note(A.inputs)
    if vf:
        st.warning(vf)
tn = ui.tensile_note(A.inputs, r3, r2)
if tn:
    (st.warning if tn[0] == "warning" else st.info)(tn[1])
if A.inputs.tension != "none":
    st.caption(f"{A.inputs.tension}: sensitivity, not a calibrated contact law.")
mb, s, a = A.inputs.hb
if A.inputs.model == "mc":
    strength_txt = f"c′ = {A.inputs.c_mc:g} kPa, φ′ = {A.inputs.phi_mc:g}°"
else:
    strength_txt = f"m_b = {mb:.6f}, s = {s:.7f}, a = {a:.6f}, ν = {A.inputs.nu:g}"
water_txt = (f"u = r_u γ h with r_u = {A.inputs.ru:g}" if A.inputs.ru is not None else
             "water table" if A.inputs.wet else "dry")
st.caption(f"Strength: {E.MODELS[A.inputs.model]}. {strength_txt}. Water: {water_txt}. "
           f"Iteration path: {E.ITERATIONS[A.inputs.iteration]}. "
           + ("Base-circle grid (exit to entry, floor columns allowed): "
              f"{int((r3.mesh.x > 0).sum())} floor columns, {int((r2.mesh.x > 0).sum())} floor slices. "
              if A.inputs.base_circle else "")
           + f"Sliding volume (full mass) {r3.mesh.volume:,.4g} m³; 2D area {r2.mesh.volume:,.4g} m²."
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
tab_search, tab_mesh, tab_nu = st.tabs(["Critical surface (search)", "Mesh study (eq. 4.17)",
                                        "ν sweep (Moscarda closure)"])
with tab_search:
    search_panel(inp)
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
