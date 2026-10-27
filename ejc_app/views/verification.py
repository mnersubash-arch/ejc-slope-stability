"""Page 7: run every verification target."""
import pandas as pd
import streamlit as st

import verification as V
import ui

st.title("Verification")
st.write("Runs every target of the brief (section 3) on the engine and shows target, computed value, difference and "
         "PASS / FAIL. INFO rows are reference values, not targets. Tolerances: FoS given to 9 decimals to 1e-6; "
         "values given to 3-4 decimals to rounding; fixtures to the decimals shown.")
slow = st.checkbox("Include the slow targets (240x240, 320x320 and 5120-slice meshes; about 1-2 minutes)")
if st.button("Run verification", type="primary"):
    status = st.status("Running the verification targets...", expanded=True)
    rows = V.verification_rows(V.Runner(), slow=slow, progress=lambda s: status.write(s))
    status.update(label="Done", state="complete", expanded=False)
    st.session_state.verification = (slow, rows)
if "verification" in st.session_state:
    slow_used, rows = st.session_state.verification
    sm = V.summary(rows)
    c = st.columns(3)
    c[0].metric("PASS", sm["PASS"])
    c[1].metric("FAIL", sm["FAIL"])
    c[2].metric("INFO", sm["INFO"])
    if sm["FAIL"]:
        st.error(f"{sm['FAIL']} target(s) do not match. They are listed in red below.")
    else:
        st.success("Every target matches." + ("" if slow_used else " The slow targets were not run."))
    df = pd.DataFrame([r.__dict__ for r in rows]).rename(columns=str.capitalize)
    st.dataframe(df.style.map(ui.color_status, subset=["Status"]), hide_index=True, width="stretch",
                 height=700)
    st.download_button("Download the table (CSV)", df.to_csv(index=False).encode(), "verification.csv", "text/csv")
