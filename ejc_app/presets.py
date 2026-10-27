"""Preset cases: this thesis and previous studies. Each preset fills every input and names its source.

Targets are verification targets (see ``verification.py``); reference values are shown for comparison only.
No Streamlit imports.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import engine as E

MAIN = E.Inputs()
CARRANZA = MAIN.with_(xc=207.28, zc=586.53, R=622.08)

SRC_THESIS = ("This thesis (UWA GENG5512). Main sphere x_c = 355.7423115, z_c = 644.4192949, R = 725.2669975 m: "
              "Moscarda (2025) Appendix E, mirrored in x. H = 360 m, beta = 50 deg, gamma = 27 kN/m3, "
              "sigma_ci = 77.7 MPa, m_i = 10, GSI = 42, D = 0, nu = 0.271.")
SRC_CT = ("Carranza-Torres & Hormazabal (2020), published critical circle x_c = 207.28, z_c = 586.53, R = 622.08 m; "
          "dry, a = 0.5. Targets are the thesis records for this circle.")
SRC_MOSCARDA = ("Moscarda (2025). His reported 2D FoS comes from his own slice model and the inherited route. His "
                "spreadsheet cells hold nu = 0.275-0.276, but the reported values are what it gives at nu = 0.271, so "
                "nu = 0.271 is used here. Reference values only, not targets.")


@dataclass(frozen=True)
class Preset:
    key: str
    label: str
    group: str
    inputs: Optional[E.Inputs]          # None = empty slot (inputs to be entered from the source)
    source: str
    targets: str = ""
    reference: str = ""
    slow: bool = False
    action: str = ""                    # "nu_sweep" or "" (the original-closure check follows from the model)


def _p(*a, **k) -> Preset:
    return Preset(*a, **k)


PRESETS: list[Preset] = [
    # ---- this thesis, main sphere, Moscarda closure, EJC
    _p("thesis_dry_05", "Dry, a = 0.5", "This thesis", MAIN, SRC_THESIS,
       "3D 2.706024342, 2D 2.447910895 (20x20, 40 slices). GHB: 2.554 / 2.310."),
    _p("thesis_wet_05", "Wet, a = 0.5", "This thesis", MAIN.with_(wet=True), SRC_THESIS + " Wet: H_w = 360 m, beta_w = 50 deg.",
       "3D 1.272331405, 2D 1.098069496. GHB: 1.200 / 1.036."),
    _p("thesis_dry_gsi", "Dry, a = 0.509923 (from GSI)", "This thesis", MAIN.with_(a_override=0.509923), SRC_THESIS,
       "3D 2.624, 2D 2.376. GHB: 2.476 / 2.242."),
    _p("thesis_wet_gsi", "Wet, a = 0.509923 (from GSI)", "This thesis", MAIN.with_(wet=True, a_override=0.509923),
       SRC_THESIS, "3D 1.228, 2D 1.061. GHB: 1.158 / 1.000."),
    _p("thesis_native_dry", "Native v15 iteration path, dry", "This thesis", MAIN.with_(iteration="native_v15"),
       SRC_THESIS + " Native v15 iteration path of the thesis workbook (start F = 1, Hungr explicit form, projection).",
       "3D 2.706024342 in 34 iterations."),
    _p("thesis_native_wet", "Native v15 iteration path, wet", "This thesis",
       MAIN.with_(wet=True, iteration="native_v15"),
       SRC_THESIS + " Native v15 iteration path of the thesis workbook (start F = 1, Hungr explicit form, projection).",
       "3D 1.272331405 in 63 iterations."),
    _p("thesis_fm_dry", "Moment form F_M, dry", "This thesis", MAIN.with_(form="FM"), SRC_THESIS,
       "F_M 3D 2.7022 (re-solved with its own iteration)."),
    _p("thesis_fm_wet", "Moment form F_M, wet", "This thesis", MAIN.with_(form="FM", wet=True), SRC_THESIS,
       "F_M 3D 1.2695."),
    _p("thesis_c0_dry", "Tension sensitivity C0, dry", "This thesis", MAIN.with_(tension="C0"),
       SRC_THESIS + " Table 5.7 (Python, default path). Sensitivity, not a calibrated contact law.",
       "3D 2.637, 2D 2.430."),
    _p("thesis_t0_dry", "Tension sensitivity T0, dry", "This thesis", MAIN.with_(tension="T0"),
       SRC_THESIS + " Table 5.7 (Python, default path). Sensitivity, not a calibrated contact law.",
       "3D 2.597, 2D 2.429."),
    _p("thesis_c0_wet", "Tension sensitivity C0, wet", "This thesis", MAIN.with_(tension="C0", wet=True),
       SRC_THESIS + " Table 5.7 (Python, default path). Sensitivity, not a calibrated contact law.",
       "3D 0.964, 2D 0.978."),
    _p("thesis_t0_wet", "Tension sensitivity T0, wet", "This thesis", MAIN.with_(tension="T0", wet=True),
       SRC_THESIS + " Table 5.7 (Python, default path). Sensitivity, not a calibrated contact law.",
       "No accepted FoS in 3D or 2D."),
    _p("thesis_finest_dry", "Finest mesh, dry (240x240, 5120 slices; slow)", "This thesis",
       MAIN.with_(n3=240, n2=5120), SRC_THESIS, "3D 2.7112, 2D 2.4422.", slow=True),
    _p("thesis_finest_wet", "Finest mesh, wet (240x240, 5120 slices; slow)", "This thesis",
       MAIN.with_(n3=240, n2=5120, wet=True), SRC_THESIS, "3D 1.2776, 2D 1.0916.", slow=True),
    _p("thesis_nu_sweep", "nu sweep 0.05-0.49, a = 0.5", "This thesis", MAIN, SRC_THESIS,
       "Runs the Moscarda-closure FoS over nu (Calculator, Studies).", action="nu_sweep"),
    _p("thesis_orig_0271", "Original closure check, nu = 0.271 (Table 5.2)", "This thesis",
       MAIN.with_(model="ejc_original"), SRC_THESIS + " Table 5.2.",
       "sigma3,min = 34.376 MPa, sigma_n,min = 52.97 MPa against the gravity cap 1.61 MPa dry / 1.02 MPa wet. No FoS."),
    _p("thesis_orig_005", "Original closure check, nu = 0.05 (Table 5.2)", "This thesis",
       MAIN.with_(model="ejc_original", nu=0.05), SRC_THESIS + " Table 5.2.",
       "Minimum sigma_n 0.92 MPa; 43 of 76 dry columns have gamma h - u below it. No FoS."),
    # ---- previous studies
    _p("ct_ghb", "Carranza-Torres & Hormazabal (2020), GHB", "Previous studies", CARRANZA.with_(model="ghb"), SRC_CT,
       "2D GHB 40 slices 2.003873 (brief: 2.0039); 3D GHB 20x20 2.141937 (thesis value 2.142).",
       "Published FoS 2.01 (their method)."),
    _p("ct_ejc", "Carranza-Torres & Hormazabal (2020), EJC", "Previous studies", CARRANZA, SRC_CT,
       "2D EJC 40 slices 2.122598; 3D EJC 20x20 2.269068.", "Published FoS 2.01 (their method, GHB)."),
    _p("ct_ghb_fine", "Carranza-Torres & Hormazabal (2020), GHB, 320x320 (slow)", "Previous studies",
       CARRANZA.with_(model="ghb", n3=320), SRC_CT, "3D GHB 320x320 2.140776; 2D GHB 40 slices 2.003873.", slow=True),
    _p("ct_ejc_fine", "Carranza-Torres & Hormazabal (2020), EJC, 320x320 / 5120 (slow)", "Previous studies",
       CARRANZA.with_(n3=320, n2=5120), SRC_CT, "3D EJC 320x320 2.267835; 2D EJC 5120 slices 2.121420.", slow=True),
    _p("moscarda_1", "Moscarda (2025) Case 1: dry, a = 0.5", "Previous studies", MAIN, SRC_MOSCARDA,
       "", "Moscarda reported 2D FoS 2.24 (his spreadsheet at nu = 0.271: 2.245)."),
    _p("moscarda_2", "Moscarda (2025) Case 2: wet, a = 0.5", "Previous studies", MAIN.with_(wet=True), SRC_MOSCARDA,
       "", "Moscarda reported 2D FoS 1.57 (his spreadsheet at nu = 0.271: 1.567)."),
    _p("moscarda_3", "Moscarda (2025) Case 3: dry, a = 0.509923", "Previous studies",
       MAIN.with_(a_override=0.509923), SRC_MOSCARDA,
       "", "Moscarda reported 2D FoS 2.04 (his spreadsheet at nu = 0.271: 2.041)."),
    _p("moscarda_4", "Moscarda (2025) Case 4: wet, a = 0.509923", "Previous studies",
       MAIN.with_(wet=True, a_override=0.509923), SRC_MOSCARDA,
       "", "Moscarda reported 2D FoS 1.13 (his spreadsheet at nu = 0.271: 1.132)."),
    _p("law_2024", "Law (2024): enter inputs from source", "Previous studies", None,
       "Law (2024). Inputs not entered yet: enter the geometry and inputs from the source."),
    _p("desterre_2024", "D'Esterre (2024): enter inputs from source", "Previous studies", None,
       "D'Esterre (2024). Inputs not entered yet: enter the geometry and inputs from the source."),
]
BY_KEY = {p.key: p for p in PRESETS}
DEFAULT_KEY = "thesis_dry_05"
