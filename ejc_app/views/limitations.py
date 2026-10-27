"""Page 6: limitations (static)."""
import streamlit as st

import ui

st.markdown(f"**{ui.DISCLAIMER}**")
st.title("Limitations")
st.markdown("""
- **Spherical surfaces only.** Each analysis uses a prescribed sphere, symmetric about y = 0. The Calculator's
  critical-surface search looks for the lowest FoS among circles (2D) and spheres centred on y = 0 (3D) that meet
  the ground at both exits. It is a local search on the current mesh, not a proof of the global minimum, and it
  does not try non-spherical or asymmetric surfaces.
- **σ2 is imposed, not solved.** σ2 comes from the closure, not from a 3D stress field, and ν is inherited, not
  measured.
- **Hungr's simplified method** neglects vertical intercolumn shear and does not satisfy every force and moment
  condition. F_H and F_M differ by about 0.14–0.27 %.
- **Equivalent isotropic rock mass**, with no discontinuities.
- **Hydrostatic water surface**, with no seepage solution.
- **Tensile basal states are allowed.** C0/T0 are sensitivities only, not a calibrated contact law, and the wet
  margin depends strongly on them (about 25 %).
- **Coarse meshes underestimate FoS.** 10x10 is about 2.6 % low dry and 4.4 % low wet. Convergence is not an
  error bound.
- **Uniqueness is not proven.** The bisection scan does not prove that the root is globally unique.
- **Research tool, not a design tool.**

Notes on this implementation:

- The trial-stress projection (thesis Section 4.6) only changes the iteration path. A run is accepted only if
  no base is projected at convergence, and the projection count is reported.
- The native v15 iteration path (start at F = 1, Hungr explicit form) reaches the same F as the default path.
  Only its iteration count differs (34 / 63 against 30 / 65 for the main cases).
- The original closure (σ2 = σ3/ν − σ1) is an admissibility check only and never produces a FoS.
- Moscarda (2025) values are reference values from his own slice model and the inherited route, not targets.
""")
