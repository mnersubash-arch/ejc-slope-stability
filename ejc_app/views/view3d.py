"""Page 3: 3D view of the columns, slip surface and ground, with a plan view."""
import streamlit as st

import plots as P
import ui

st.title("3D view")
inp = ui.sidebar()
if inp.model == "ejc_original":
    st.info("The original closure is an admissibility check only; there are no column results. See the Calculator.")
    st.stop()
A, stale = ui.analysis(inp)
if A is None:
    st.stop()
r3 = A.r3
m = r3.mesh
c1, c2, c3 = st.columns([3, 1, 1])
field = c1.selectbox("Colour columns by", list(P.FIELDS), format_func=lambda k: P.FIELDS[k][0], key="v3_field")
mirror = c2.toggle("Mirror (full mass)", value=True, key="v3_mirror")
if "v3_col" not in st.session_state or st.session_state.v3_col >= m.count:
    st.session_state.v3_col = 0
col = c3.number_input(f"Highlight column (0 to {m.count - 1})", 0, m.count - 1, key="v3_col")
ui.status_badge(r3, f"3D, {m.count} active half-columns")
st.plotly_chart(P.view3d(A.inputs, r3, field, mirror, highlight=int(col)), width="stretch", key="v3_fig")
st.caption("Column prisms run from the slip surface to the ground, coloured by the selected field (contribution to "
           "resistance = τ_s A; tension flag = σ_n < 0 or in the C0/T0 tensile set). Hover over the black base points "
           "for values. Above 4,000 prisms only the coloured bases are drawn.")
st.subheader("Plan view")
st.plotly_chart(P.plan_view(r3, field, highlight=int(col), mirror=mirror), width="stretch",
                key="v3_plan")
k = st.session_state.v3_col
st.write(f"Column {k}: i = {int(m.i[k])}, j = {int(m.j[k])}, x = {m.x[k]:.2f} m, y = {m.y[k]:.2f} m, "
         f"h = {m.h[k]:.2f} m.")
