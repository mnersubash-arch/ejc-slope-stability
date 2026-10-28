"""Presets, figures and a headless run of every Streamlit page."""
import math
from pathlib import Path

import numpy as np
import pytest

import plots as P
import presets as PR
from cases import MAIN, WET, r2, r3

APP = Path(__file__).resolve().parents[1] / "app.py"


def test_presets_fill_every_input_and_name_a_source():
    keys = [p.key for p in PR.PRESETS]
    assert len(keys) == len(set(keys)) and PR.DEFAULT_KEY in PR.BY_KEY
    for p in PR.PRESETS:
        assert p.source and p.inputs is not None
        p.inputs.strength()                                   # valid strength inputs
        if p.inputs.iteration == "native_v15":
            assert p.inputs.tension == "none" and p.inputs.form == "FH"
    assert all(PR.BY_KEY[k].inputs.nu == 0.271 for k in ("moscarda_1", "moscarda_2", "moscarda_3", "moscarda_4"))
    assert all(PR.BY_KEY[k].targets == "" and PR.BY_KEY[k].reference for k in ("moscarda_1", "moscarda_4"))
    law, des = PR.law_presets(), PR.desterre_presets()
    assert len(law) == 4 and len(des) == 24
    assert all(p.inputs.model == "ghb" and not p.inputs.base_circle and p.inputs.R == p.data["R"] for p in law)
    assert all(p.inputs.model == "mc" and p.inputs.base_circle and p.inputs.ru in (0.0, 0.5) for p in des)
    assert all("needs_floor_columns" in p.info for p in des)
    assert all(p.reference.startswith("Reported (reference only)") for p in law + des)


@pytest.mark.parametrize("page", ["calculator", "cross_section", "view3d", "fbd"])
def test_pages_with_a_desterre_preset(page):
    """Mohr-Coulomb, r_u and a base circle with floor columns, on every page that draws the result."""
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    at.sidebar.selectbox(key="w_preset").select("desterre_DE-1.5-c3-ru0.5").run()
    at.switch_page(f"views/{page}.py")
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    if page == "calculator":
        assert any("needs_floor_columns: true" in i.value for i in at.info)
        assert dict((m.label, m.value) for m in at.metric)["FoS 2D (Bishop)"] != "withheld"


@pytest.mark.parametrize("field", list(P.FIELDS))
def test_figures_build(field):
    for inp in (MAIN, WET):
        a3, a2 = r3(inp), r2(inp)
        P.cross_section(inp, a2, field)
        P.view3d(inp, a3, field, mirror=True, highlight=3)
        P.plan_view(a3, field, highlight=3, mirror=True)


def test_tensile_bases_summary():
    """QA/QC finding S3: the summary counts bases with sigma_n < 0 and warns above a quarter of them."""
    import ui
    dry, wet = (r3(MAIN), r2(MAIN)), (r3(WET), r2(WET))
    assert ui.tensile_text(ui.tensile_count(dry[0])) == "6 of 76 (8%)"           # thesis §5.8: 6 of 76 dry bases
    assert ui.tensile_text(ui.tensile_count(wet[0])) == "23 of 76 (30%)"         # and 23 of 76 wet bases
    assert ui.tensile_note(MAIN, *dry)[0] == "info" and ui.tensile_note(WET, *wet)[0] == "warning"
    c0 = MAIN.with_(tension="C0")
    assert ui.tensile_note(c0, r3(c0), r2(c0)) is None                          # C0/T0 already treat those bases
    assert ui.tensile_text(None) == "–"


def test_free_body_balances():
    for inp in (MAIN, WET):
        a3, a2 = r3(inp), r2(inp)
        for k in range(a3.mesh.count):
            assert P.column_forces(inp, a3, k)["residual_rel"] <= 1e-9
        for k in range(a2.mesh.count):
            assert P.column_forces(inp, a2, k)["residual_rel"] <= 1e-9
        P.fbd3d(inp, a3, 10)
        P.lever_view(inp, a3, 10)
        P.fbd2d(inp, a2, 10)


def test_lever_arm_is_rho():
    a3 = r3(MAIN)
    f = P.column_forces(MAIN, a3, 5)
    assert f["lever"] == pytest.approx(f["rho"], rel=1e-12)
    assert np.dot(f["n"], f["s"]) == pytest.approx(0.0, abs=1e-12)        # shear lies in the base plane
    assert f["s"][2] == pytest.approx(f["kap"]) and f["n"][2] == pytest.approx(f["g"])


def test_strength_figures():
    mb, s, _ = MAIN.hb
    P.strength_curves(77.7, mb, s, 0.5, 0.271, 50.0)
    fig, rows = P.mohr(77.7, mb, s, 0.5, 0.271, 400.0)
    assert rows[0]["state"] == "ACCEPTED" and abs(rows[0]["phi_deg"] - 59.90) < 0.005


@pytest.mark.parametrize("page", ["calculator", "cross_section", "view3d", "fbd", "strength", "limitations",
                                  "verification"])
def test_pages_run_headless(page):
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    at.switch_page(f"views/{page}.py")
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    if page == "calculator":
        assert dict((m.label, m.value) for m in at.metric)["FoS 3D"] == "2.706024"


@pytest.mark.parametrize("page", ["calculator", "view3d", "fbd"])
def test_pages_with_no_active_columns(page):
    """N = 2 leaves no active 3D column: the pages that use the 3D result report it instead of raising."""
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    at.number_input(key="nb_n3").set_value(2)
    at.run()
    at.switch_page(f"views/{page}.py")
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    assert any("No active columns" in e.value for e in at.error)


def test_sidebar_keeps_its_values_across_pages():
    """Streamlit resets widgets that are new to a page; the sidebar must show the stored inputs on every page."""
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    at.sidebar.selectbox(key="w_preset").select("desterre_DE-1.5-c3-ru0").run()
    at.number_input(key="nb_n3").set_value(24).run()
    for page in ("cross_section", "view3d", "fbd", "strength", "calculator"):
        at.switch_page(f"views/{page}.py")
        at.run()
        assert at.number_input(key="nb_n3").value == 24 and at.number_input(key="nb_c_mc").value == 21.6
        assert at.radio(key="w_strength").value == "MC" and at.checkbox(key="w_base_circle").value
        assert at.sidebar.selectbox(key="w_preset").value == "desterre_DE-1.5-c3-ru0"


def test_whole_mass_figures():
    """All forces on the sliding mass: every element drawn on the default meshes, groups on fine ones."""
    import engine as E
    for inp in (MAIN, WET, MAIN.with_(n3=40, n2=640)):
        a3, a2 = r3(inp), r2(inp)
        f2, size2 = P.mass_fbd2d(inp, a2, E.mass_forces(inp, a2))
        f3, size3 = P.mass_fbd3d(inp, a3, E.mass_forces(inp, a3))
        n_w2 = len(next(t for t in f2.data if t.name.startswith("W,")).x) // 3       # 3 points per arrow
        n_w3 = len(next(t for t in f3.data if t.name and t.name.startswith("W,")).x) // 6   # both halves
        assert size2 == max(1, math.ceil(a2.mesh.count / 50)) and n_w2 == math.ceil(a2.mesh.count / size2) <= 50
        assert n_w3 <= 120 and (size3 == 1) == (a3.mesh.count <= 120)
        if inp.n3 == 20:
            assert (size2, size3, n_w2, n_w3) == (1, 1, a2.mesh.count, a3.mesh.count)   # every element drawn
        assert any(t.name.startswith("N′") for t in f2.data) and any(t.name and t.name.startswith("S,") for t in f3.data)
        assert any(t.name.startswith("U,") for t in f2.data) == inp.wet
        for mf in (E.mass_forces(inp, a3), E.mass_forces(inp, a2)):
            P.force_polygon(mf)
            assert len(P.mass_rows(mf)) == 5 and "PASS" in P.mass_balance(mf)[0]


def test_fbd_page_whole_mass_views():
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    at.switch_page("views/fbd.py")
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    assert any("Vertical:" in m.value and "PASS" in m.value for m in at.markdown)
    at.radio(key="mass_view").set_value("3D full mass").run()
    assert not at.exception, [e.value for e in at.exception]
    assert any("F_H balances ΣS = ΣWκ instead" in m.value for m in at.markdown)


def test_calculator_search_tab():
    """Critical-surface search on a Mohr-Coulomb preset (fast), then Use this circle reruns the analysis on it."""
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=300)
    at.run()
    at.sidebar.selectbox(key="w_preset").select("desterre_DE-1.5-c1-ru0").run()
    at.button(key="srch_run").click().run()
    assert not at.exception, [e.value for e in at.exception]
    table = next(d.value for d in at.dataframe if "FoS on the critical surface" in list(d.value.iloc[:, 0]))
    rows = dict(zip(table.iloc[:, 0], table["2D, Bishop"]))
    assert rows["FoS on the prescribed surface"] == "2.610501" and rows["FoS on the critical surface"] == "2.605316"
    assert rows["Mesh check (finer mesh)"] == "160 slices: 2.606797 (+0.1 %)"
    at.button(key="srch_use2").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert dict((m.label, m.value) for m in at.metric)["FoS 2D (Bishop)"] == "2.605316"
    assert at.number_input(key="nb_R").value == pytest.approx(3.540311152206755, rel=1e-15)


def test_vertical_back_face_warning():
    """A circle whose arc ends below the ground at its back extreme: the Calculator warns about the vertical face."""
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    assert not any("Vertical" in w.value for w in at.warning)
    for k, v in (("xc", 25.0), ("zc", 250.0), ("R", 251.0)):
        at.number_input(key=f"nb_{k}").set_value(v).run()
    assert not at.exception, [e.value for e in at.exception]
    w = [x.value for x in at.warning if "Vertical back face" in x.value]
    assert len(w) == 1 and "19.34 m high" in w[0] and "(269.34 m)" in w[0]
    assert dict((m.label, m.value) for m in at.metric)["FoS 3D"] == "2.758886"           # the result itself is unchanged
