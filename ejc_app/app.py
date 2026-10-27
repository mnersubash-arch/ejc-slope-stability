"""EJC 3D/2D rock-slope stability app (GENG5512). Run with:  streamlit run app.py"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import streamlit as st  # noqa: E402

st.set_page_config(page_title="EJC slope stability", page_icon=":material/landscape:", layout="wide")

pages = [
    st.Page("views/calculator.py", title="Calculator", icon=":material/calculate:", default=True),
    st.Page("views/cross_section.py", title="Cross-section", icon=":material/timeline:"),
    st.Page("views/view3d.py", title="3D view", icon=":material/view_in_ar:"),
    st.Page("views/fbd.py", title="Free-body diagram", icon=":material/call_split:"),
    st.Page("views/strength.py", title="Strength explorer", icon=":material/show_chart:"),
    st.Page("views/limitations.py", title="Limitations", icon=":material/warning:"),
    st.Page("views/verification.py", title="Verification", icon=":material/fact_check:"),
]
st.navigation(pages).run()
