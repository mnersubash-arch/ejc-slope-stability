"""Page 2: central cross-section y = 0 with the Bishop slices."""
import pandas as pd
import streamlit as st

import engine as E
import plots as P
import ui

st.title("Cross-section (y = 0)")
inp = ui.sidebar()
if inp.model == "ejc_original":
    st.info("The original closure is an admissibility check only; there are no slice results. See the Calculator.")
    st.stop()
A, stale = ui.analysis(inp)
if A is None:
    st.stop()
r2 = A.r2
names = {k: P.FIELDS[k][0] for k in ("sigma_n", "phi", "c", "tau", "u")}
field = st.radio("Colour slices by", list(names), format_func=names.get, horizontal=True, key="xs_field")
ui.status_badge(r2, f"2D, {r2.mesh.count} active slices")
st.plotly_chart(P.cross_section(A.inputs, r2, field), width="stretch", key="xs_fig")
st.caption("Hover over a slice for σ_n, φ_i, c_i, τ_s and u. The slip surface is drawn where it lies inside the mass; "
           "show the full circle from the legend. Dashed lines join the centre to the ends of the slip surface.")
with st.expander("Slice table"):
    st.dataframe(pd.DataFrame(E.columns_table(r2)), width="stretch", height=320)
