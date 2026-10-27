"""Presets, figures and a headless run of every Streamlit page."""
from pathlib import Path

import numpy as np
import pytest

import engine as E
import plots as P
import presets as PR
from cases import MAIN, WET, r2, r3

APP = Path(__file__).resolve().parents[1] / "app.py"


def test_presets_fill_every_input_and_name_a_source():
    keys = [p.key for p in PR.PRESETS]
    assert len(keys) == len(set(keys)) and PR.DEFAULT_KEY in PR.BY_KEY
    for p in PR.PRESETS:
        assert p.source
        if p.inputs is None:
            assert p.key in ("law_2024", "desterre_2024")
            continue
        p.inputs.strength()                                   # valid strength inputs
        if p.inputs.iteration == "native_v15":
            assert p.inputs.tension == "none" and p.inputs.form == "FH"
    assert all(PR.BY_KEY[k].inputs.nu == 0.271 for k in ("moscarda_1", "moscarda_2", "moscarda_3", "moscarda_4"))
    assert all(PR.BY_KEY[k].targets == "" and PR.BY_KEY[k].reference for k in ("moscarda_1", "moscarda_4"))


@pytest.mark.parametrize("field", list(P.FIELDS))
def test_figures_build(field):
    for inp in (MAIN, WET):
        a3, a2 = r3(inp), r2(inp)
        P.cross_section(inp, a2, field)
        P.view3d(inp, a3, field, mirror=True, highlight=3)
        P.plan_view(a3, field, highlight=3, mirror=True)


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
