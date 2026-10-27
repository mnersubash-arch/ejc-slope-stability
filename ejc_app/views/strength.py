"""Page 5: strength explorer: sigma1-sigma3 curves and the Mohr circle with its tangent at a chosen sigma_n."""
import pandas as pd
import streamlit as st

import plots as P
import ui

st.title("Strength explorer")
inp = ui.sidebar()
mb, s, a = inp.hb
if inp.model == "mc":
    st.info(f"The strength explorer covers the Hoek-Brown models. The selected Mohr-Coulomb strength is the straight line "
            f"τ = c′ + σ_n tan φ′ with c′ = {inp.c_mc:g} kPa and φ′ = {inp.phi_mc:g}°; the curves below use the rock-mass "
            f"inputs of the sidebar.")
st.caption(f"Rock mass from the sidebar: σ_ci = {inp.sigma_ci:g} MPa, m_b = {mb:.6f}, s = {s:.7f}, a = {a:.6f}, "
           f"ν = {inp.nu:g}. Tensile limit σ3 = {-s * inp.sigma_ci / mb * 1000:.3f} kPa.")
c1, c2 = st.columns(2)
s3max = c1.number_input("σ3 range up to (MPa)", 0.1, 5000.0, 50.0, 1.0,
                        help="The original closure is admissible only from σ3 ≈ 34 MPa at ν = 0.271.")
sn = c2.number_input("σ_n for the Mohr circle (kPa)", -97.0, 1.0e6, 400.0, 50.0)
st.plotly_chart(P.strength_curves(inp.sigma_ci, mb, s, a, inp.nu, s3max), width="stretch", key="sx_curves")
st.caption("Each EJC curve is drawn only where the state is admissible: unsquared EJC, real state, ordering, M ≥ 1. "
           "The shaded band is the ordering range of the original closure.")
fig, rows = P.mohr(inp.sigma_ci, mb, s, a, inp.nu, sn)
st.plotly_chart(fig, width="stretch", key="sx_mohr")
st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
