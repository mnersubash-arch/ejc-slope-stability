"""Page 4: free-body diagram of one column (3D) and of the matching slice (2D)."""
import numpy as np
import pandas as pd
import streamlit as st

import plots as P
import ui

st.title("Free-body diagram")
inp = ui.sidebar()
if inp.model == "ejc_original":
    st.info("The original closure is an admissibility check only; there are no column forces. See the Calculator.")
    st.stop()
A, stale = ui.analysis(inp)
if A is None:
    st.stop()
r3, r2 = A.r3, A.r2
m3, m2 = r3.mesh, r2.mesh
if r3.tau_used is None:
    ui.status_badge(r3, "3D")
    st.warning("The 3D run did not reach a final state, so there are no column forces to draw.")
    st.stop()
if r3.fos is None:
    st.warning(f"3D: {r3.status}. The diagram uses the last trial F = {r3.f_trial:.6f}, which is not a FoS.")

# ---- choose the column: click in the plan view or pick its index
if "fbd_col" not in st.session_state or st.session_state.fbd_col >= m3.count:
    st.session_state.fbd_col = int(np.argmax(m3.h))
st.markdown("Click a column in the plan view, or pick its index.")
ev = st.plotly_chart(P.plan_view(r3, "sigma_n", highlight=st.session_state.fbd_col), width="stretch",
                     key="fbd_plan", on_select="rerun", selection_mode="points")
pts = (ev.get("selection", {}) or {}).get("points", []) if isinstance(ev, dict) else getattr(ev.selection, "points", [])
if pts:
    sel = [p for p in pts if p.get("curve_number") == 1]
    if sel:
        k_new = int(sel[0].get("point_index", sel[0].get("point_number", 0)))
        if st.session_state.get("fbd_last_click") != k_new:
            st.session_state.fbd_last_click = k_new
            st.session_state.fbd_col = k_new
            st.rerun()
st.number_input(f"Column index (0 to {m3.count - 1})", 0, m3.count - 1, key="fbd_col")
k = int(st.session_state.fbd_col)

f = P.column_forces(A.inputs, r3, k)
st.subheader(f"Column {k} (i = {int(m3.i[k])}, j = {int(m3.j[k])}), x = {m3.x[k]:.2f} m, y = {m3.y[k]:.2f} m")
c1, c2 = st.columns([3, 2])
c1.plotly_chart(P.fbd3d(A.inputs, r3, k), width="stretch", key="fbd3d")
c2.plotly_chart(P.lever_view(A.inputs, r3, k), width="stretch", key="fbd_lever")
rows = [("W (weight)", f["W"], "kN", "vertical, down"),
        ("N' = σ_n A (effective normal)", f["N"], "kN", "along the base normal n, into the column"),
        ("U = u A (pore force)", f["U"], "kN", "along n"),
        ("S = τ_s A / F (mobilised shear)", f["S"], "kN", "in the base plane, in x-z, resisting sliding"),
        ("F used", f["F"], "-", "accepted FoS" if r3.fos is not None else "last trial F (not a FoS)"),
        ("σ_n", r3.sn[k], "kPa", ""), ("u", m3.u[k], "kPa", ""), ("τ_s used", r3.tau_used[k], "kPa", ""),
        ("A (base area)", f["area"], "m²", "A = A_h / g_z"), ("g_z (vertical component of n)", f["g"], "-", ""),
        ("κ (vertical component of the shear direction)", f["kap"], "-", "κ = (x_c − x)/ρ"),
        ("ρ (lever arm to the rotation axis)", f["rho"], "m", "ρ = √(R² − y²)")]
st.dataframe(pd.DataFrame(rows, columns=["Quantity", "Value", "Unit", "Note"]), hide_index=True,
             width="stretch")
ok = f["residual_rel"] <= 1e-9
st.write(f"Vertical balance W − (N' + U) g_z − S κ = **{f['residual']:.3e} kN** "
         f"(relative {f['residual_rel']:.2e}; limit 1e-9) → **{'PASS' if ok else 'FAIL'}**")

# ---- matching 2D slice
st.subheader("Matching 2D slice (central section, per metre)")
if r2.tau_used is None:
    ui.status_badge(r2, "2D")
    st.stop()
k2_default = int(np.argmin(np.abs(m2.x - m3.x[k])))
if st.session_state.get("fbd_slice_for") != k:
    st.session_state.fbd_slice_for = k
    st.session_state.fbd_slice = k2_default
st.number_input(f"Slice index (0 to {m2.count - 1}); default: the slice nearest to the column's x", 0,
                m2.count - 1, key="fbd_slice")
j = int(st.session_state.fbd_slice)
g = P.column_forces(A.inputs, r2, j)
c1, c2 = st.columns([3, 2])
c1.plotly_chart(P.fbd2d(A.inputs, r2, j), width="stretch", key="fbd2d")
alpha = float(np.degrees(np.arcsin(m2.kap[j])))
rows2 = [("W", g["W"], "kN/m"), ("N' = σ_n A", g["N"], "kN/m"), ("U = u A", g["U"], "kN/m"),
         ("S = τ_s A / F", g["S"], "kN/m"), ("F used", g["F"], "-"), ("α (base inclination)", alpha, "deg"),
         ("σ_n", r2.sn[j], "kPa"), ("u", m2.u[j], "kPa")]
c2.dataframe(pd.DataFrame(rows2, columns=["Quantity", "Value", "Unit"]), hide_index=True, width="stretch")
ok2 = g["residual_rel"] <= 1e-9
c2.write(f"W − (N' + U) cos α − S sin α = **{g['residual']:.3e} kN/m** (relative {g['residual_rel']:.2e}) → "
         f"**{'PASS' if ok2 else 'FAIL'}**")
