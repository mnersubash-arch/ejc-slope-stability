"""Preset cases: this thesis and previous studies. Each preset fills every input and names its source.

Targets are verification targets (see ``verification.py``); reference values are shown for comparison only.
The Law (2024) and D'Esterre (2024) cases are read from ``data/presets_law_desterre.json`` (engine coordinates).
No Streamlit imports.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import engine as E

DATA_FILE = Path(__file__).resolve().parent / "data" / "presets_law_desterre.json"
DATA = json.loads(DATA_FILE.read_text())
BASKOYAK_FILE = Path(__file__).resolve().parent / "data" / "presets_baskoyak.json"
BASKOYAK = json.loads(BASKOYAK_FILE.read_text())

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
    info: str = ""                      # extra facts shown with the preset (e.g. needs_floor_columns)
    data: dict = field(default_factory=dict, compare=False)   # the source record (Law / D'Esterre JSON)


def _p(*a, **k) -> Preset:
    return Preset(*a, **k)


def values_text(d: dict) -> str:
    return "; ".join(f"{k.replace('_', ' ')} {v}" for k, v in d.items())


def law_presets() -> list[Preset]:
    """Law (2024): GHB, R exactly as given in the file (through the toe), default grid."""
    out = []
    for c in DATA["Law_2024"]:
        inp = E.Inputs(H=c["H_m"], beta=c["beta_deg"], xc=c["x_c_engine"], zc=c["z_c"], R=c["R"], gamma=c["gamma_kNm3"],
                       sigma_ci=c["sigma_ci_MPa"], mi=c["m_i"], gsi=c["GSI"], D=c["D"], a_override=c["a"], model="ghb",
                       wet=c["water"] != "dry")
        circle = f"x_c = {c['x_c_engine']}, z_c = {c['z_c']}, R = {c['R']} m (engine coordinates; R through the toe"
        circle += f", reported {c['R_reported']})" if "R_reported" in c else ")"
        notes = " ".join(s for s in (c.get("note", ""), c.get("geometry_note", "")) if s)
        out.append(_p(f"law_{c['id'][-1]}", f"{c['id']}: {c['label']}", "Law (2024)", inp,
                      f"Law (2024) GENG5512 report. {c['strength']}. {circle}. Water: {c['water']}. {notes}",
                      "Cross-check against the reference code (lem_ref.py): "
                      + "; ".join(f"{k.replace('_', ' ')} slices {v}" for k, v in c["engine_check"].items()) + ".",
                      "Reported (reference only): " + values_text(c["reported"]) + ".", data=c))
    return out


def desterre_presets() -> list[Preset]:
    """D'Esterre (2024): Mohr-Coulomb, u = r_u gamma h, base-circle grid (exit to entry, floor columns allowed)."""
    common = DATA["DEsterre_2024"]["common"]
    out = []
    for c in DATA["DEsterre_2024"]["cases"]:
        ru = float(re.search(r"ru = ([0-9.]+)", c["pore_pressure"]).group(1))
        inp = E.Inputs(H=c["H_m"], beta=c["beta_deg"], xc=c["x_c_engine"], zc=c["z_c"], R=c["R"], gamma=c["gamma_kNm3"],
                       model="mc", c_mc=c["c_kPa"], phi_mc=c["phi_deg"], ru=ru, base_circle=True)
        rep = c["reported"]
        out.append(_p(f"desterre_{c['id']}", c["label"], "D'Esterre (2024)", inp,
                      f"{c['source']}. Slope {c['slope_H_to_L']} (H = {c['H_m']} m), circle x_c = {c['x_c_engine']}, "
                      f"z_c = {c['z_c']}, R = {c['R']} m (engine coordinates; {common['rounding']}). "
                      f"γ = {common['note_gamma']}. His 3D model: {common['his_3D_model']}.",
                      "",
                      f"Reported (reference only): Slide2 2D {rep['Slide2_2D']}; D'Esterre spreadsheet 3D (Hungr) "
                      f"{rep['DEsterre_spreadsheet_3D_Hungr']}; Slide3 3D {rep['Slide3_3D']}; spreadsheet 3D / Slide2 "
                      f"{rep['ratio_spreadsheet3D_over_Slide2']}.",
                      info=f"needs_floor_columns: {str(c['needs_floor_columns']).lower()} (circle exits "
                           f"{c['exit_ahead_of_toe_m']} m ahead of the toe; lowest point z = {c['lowest_point_z_m']} m).",
                      data=c))
    return out


def baskoyak_presets() -> list[Preset]:
    """Baskoyak barite pit (a real failure): GHB with the SLIDE back-analysis inputs, GHB with site and guideline inputs
    (nothing fitted to the failure), and the EJC with a derived Poisson's ratio. All use the app's planar section and the
    traced failure circle with the base-circle grid, from data/presets_baskoyak.json; checked by qa/qa_baskoyak.py."""
    b = BASKOYAK
    site, sec, circ, nu = b["site"], b["geometry"]["app_section"], b["geometry"]["failure_circle"], b["poisson_ratio"]
    base = E.Inputs(H=sec["H_m"], beta=sec["beta_deg"], xc=circ["x_c_engine"], zc=circ["z_c"], R=circ["R"],
                    gamma=site["gamma_kNm3"], sigma_ci=site["sigma_ci_MPa"], mi=site["m_i"], a_override=None,
                    nu=nu["nu"], model="ghb", base_circle=True)
    source = (f"Baskoyak barite open pit, western Anatolia (Turkey): a circular failure through one bench in heavily "
              f"jointed, dry schist (Sonmez & Ulusay 1999, Case 2). gamma = {site['gamma_kNm3']} kN/m3, "
              f"sigma_ci = {site['sigma_ci_MPa']} MPa, m_i = {site['m_i']:g}. Section and failure circle traced from "
              f"Sonmez & Ulusay (1999) Fig. 7: H = {sec['H_m']} m, overall face {sec['beta_deg']} deg (planar here; the "
              f"real section has a 2.2 m bench at 8.5 m and ground rising 5.1 deg behind the crest); failure circle "
              f"x_c = {circ['x_c_engine']}, z_c = {circ['z_c']}, R = {circ['R']:.2f} m (engine coordinates), which exits "
              f"on the pit floor 5.3 m ahead of the toe (base-circle grid). Record: data/presets_baskoyak.json; checks: "
              f"qa/qa_baskoyak.py.")
    failed = "the slope failed on the traced circle (Sonmez & Ulusay 1999), so F = 1 at failure"
    slide = b["reported"]["LLC2016_SLIDE"]
    site_note = ("GSI 16.5 from the field survey (Sonmez & Ulusay 2002, their modified GSI chart); D = 0.7 for a pit "
                 "excavated mechanically without blasting (Hoek et al. 2002, as Li et al. 2016 note for this mine). "
                 "No input is fitted to the failure.")
    nu_note = (f"nu = {nu['nu']} is the rock-mass Poisson's ratio derived from Aydan, Tokashiki & Ulusay (2015) and "
               f"Aydan's intact-rock relation (range {nu['range'][0]:.2f}-{nu['range'][1]:.2f}; not measured for this "
               f"site). At "
               f"the thesis nu = 0.271 the EJC has no admissible state above sigma_n = 95 kPa for this rock mass.")
    return [
        _p("baskoyak_slide_ghb", "GHB, SLIDE back-analysis inputs (GSI 16, D 0.68)", "Baskoyak (real failure)",
           base.with_(gsi=slide["GSI"], D=slide["D"]),
           source + f" GSI {slide['GSI']:g} and D {slide['D']} are the SLIDE back-analysis of Li, Lo & Cassidy (2016).",
           "", f"Reported (reference only): SLIDE (Bishop) F = {slide['F']} on its own critical circle of the real "
               f"section (Li et al. 2016); {failed}.", data=b),
        _p("baskoyak_site_ghb", "GHB, site inputs, D 0.7 (forward prediction)", "Baskoyak (real failure)",
           base.with_(gsi=b["site"]["GSI"]["SU2002"], D=b["reported"]["LLC2016_D_adopted"]),
           source + " " + site_note, "", f"Reported (reference only): {failed}.", data=b),
        _p("baskoyak_site_ejc", f"EJC (Moscarda closure), site inputs, D 0.7, nu = {nu['nu']}", "Baskoyak (real failure)",
           base.with_(gsi=b["site"]["GSI"]["SU2002"], D=b["reported"]["LLC2016_D_adopted"], model="ejc_moscarda"),
           source + " " + site_note + " " + nu_note, "", f"Reported (reference only): {failed}.", data=b),
    ]


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
    *law_presets(),
    *desterre_presets(),
    *baskoyak_presets(),
]
BY_KEY = {p.key: p for p in PRESETS}
DEFAULT_KEY = "thesis_dry_05"
