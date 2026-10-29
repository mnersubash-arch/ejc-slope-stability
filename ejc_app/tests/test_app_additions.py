"""The additions in the Streamlit app (headless): sidebar inputs, studies, cases (file, code, page link) and the
report, with the values the web app gives for the same inputs (webapp/tools/ui_test.js)."""
import json
from pathlib import Path

import pytest

import casefile as CF
import engine as E
import report as RP

APP = Path(__file__).resolve().parents[1] / "app.py"


def _app():
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at = at_mod.AppTest.from_file(str(APP), default_timeout=600)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def _fos(at):
    m = dict((x.label, x.value) for x in at.metric)
    return m["FoS 3D"], m["FoS 2D (Bishop)"]


def _reset(at):
    [b for b in at.sidebar.button if b.label == "Reset inputs to preset"][0].click().run()


def test_sidebar_additions_match_the_web_app():
    at = _app()
    at.radio(key="w_face").set_value("benches").run()
    assert _fos(at) == ("2.706024", "2.447911")                       # one bench is the planar face
    for k, v in (("bench_n", 4), ("bench_h", 90.0), ("bench_beta", 60.0), ("berm_w", 15.0)):
        at.number_input(key=f"nb_{k}").set_value(v).run()
    assert _fos(at) == ("2.192610", "2.022105")
    at.radio(key="w_face").set_value("profile").run()
    at.text_area(key="w_profile_text").set_value("0, 0\n-10, 5\n-5, 8").run()
    assert any("x must decrease" in e.value for e in at.sidebar.error)
    _reset(at)
    at.radio(key="w_wall").set_value("concave").run()
    assert any("no default" in e.value for e in at.error) and at.number_input(key="nb_plan_r").value is None
    at.number_input(key="nb_plan_r").set_value(800.0).run()
    assert _fos(at) == ("2.755516", "2.447911")
    at.radio(key="w_wall").set_value("straight").run()
    at.radio(key="w_water").set_value("Wet").run()
    wet = _fos(at)
    at.radio(key="w_water").set_value("Piezo").run()
    assert _fos(at) == wet == ("1.272331", "1.098069")                # the default line is the water table
    at.radio(key="w_water").set_value("Dry").run()
    at.checkbox(key="w_crack_on").check().run()
    at.number_input(key="nb_crack_depth").set_value(20.0).run()
    at.number_input(key="nb_crack_water").set_value(0.5).run()
    at.number_input(key="nb_kh").set_value(0.1).run()
    assert _fos(at) == ("2.339066", "2.101695")
    assert any("Tension crack:" in c.value and "P = 490.5 kN/m" in c.value for c in at.caption)
    for page in ("cross_section", "view3d", "fbd"):
        at.switch_page(f"views/{page}.py")
        at.run()
        assert not at.exception, (page, [e.value for e in at.exception])


def test_studies_on_the_calculator():
    at = _app()
    at.button(key="cmp_run").click().run()
    df = [d.value for d in at.dataframe if "3D / 2D" in d.value.columns][0]
    assert [f"{v:.6f}" for v in df["FoS 3D"][:2]] == ["2.706024", "2.553668"]
    assert [f"{v:.6f}" for v in df["FoS 2D"][:2]] == ["2.447911", "2.310255"]
    at.selectbox(key="back_name").set_value("D").run()
    at.number_input(key="back_lo_D").set_value(0.0)
    at.number_input(key="back_hi_D").set_value(1.0)
    at.number_input(key="back_target").set_value(2.0)
    at.selectbox(key="back_dim").set_value(2)
    at.button(key="back_run").click().run()
    assert any("D, disturbance factor = 0.417641** gives F = 2.000000 (2D)" in s.value for s in at.success)
    at.selectbox(key="sweep_name").set_value("gsi").run()
    at.number_input(key="sweep_lo_gsi").set_value(40.0)
    at.number_input(key="sweep_hi_gsi").set_value(80.0)
    at.number_input(key="sweep_n").set_value(3)
    at.checkbox(key="sweep_3").uncheck()
    at.button(key="sweep_run").click().run()
    sw = [d.value for d in at.dataframe if "GSI" in d.value.columns][0]
    assert [f"{v:.6f}" for v in sw["EJC 2D"]] == ["2.357616", "3.925349", "9.542099"]
    assert [f"{v:.6f}" for v in sw["GHB 2D"]] == ["2.225298", "3.698480", "8.963219"]
    at.number_input(key="mc_sd_gsi").set_value(5.0)
    at.number_input(key="mc_n").set_value(10)
    at.button(key="mc_run").click().run()
    mc = [d.value for d in at.dataframe if "Accepted" in d.value.columns][0]
    ref = E.monte_carlo(E.Inputs(), [("gsi", 42.0, 5.0, 10.0, 100.0)], 10, 1, None, [2])["stats"]["ejc_moscarda2"]
    assert mc["Accepted"][0] == "10 of 10" and mc["Mean F"][0] == ref["mean"]
    assert not at.exception, [e.value for e in at.exception]
    assert [b.proto.label for b in at.get("download_button")].count("Download the report (HTML)") == 1


def test_cases_open_from_code_web_app_and_link():
    at = _app()
    at.button(key="case_code_make").click().run()
    code = at.text_area(key="case_code").value
    c = json.loads(code)
    assert c["format"] == CF.CASE_FORMAT and c["results"]["r3"]["status"] == "ACCEPTED" and "streamlit_state" in c
    at.number_input(key="nb_gsi").set_value(50.0).run()
    assert _fos(at)[0] == "3.315975"
    at.text_area(key="case_code").set_value(code).run()
    at.button(key="case_code_open").click().run()
    assert _fos(at) == ("2.706024", "2.447911") and any("Opened" in s.value for s in at.sidebar.success)
    # a case written by the web app: its page state is for the web page; the engine inputs are read here
    web = {"format": CF.CASE_FORMAT, "name": "from the web app", "preset": "thesis_dry_05",
           "state": {"face": "benches", "gsi": 1}, "inputs": dict(CF.inputs_dict(E.Inputs()), gsi=45)}
    at.text_area(key="case_code").set_value(json.dumps(web)).run()
    at.button(key="case_code_open").click().run()
    assert _fos(at) == ("2.902160", "2.601524") and at.radio(key="w_face").value == "planar"
    at.text_area(key="case_code").set_value("{").run()
    at.button(key="case_code_open").click().run()
    assert any("not valid JSON" in e.value for e in at.sidebar.error)
    # the page link
    at_mod = pytest.importorskip("streamlit.testing.v1")
    at2 = at_mod.AppTest.from_file(str(APP), default_timeout=600)
    at2.query_params["case"] = CF.encode_link(CF.make_case(E.Inputs(gsi=50.0, kh=0.05), name="Linked"))
    at2.run()
    assert not at2.exception and _fos(at2)[0] == f"{E.run_3d(E.Inputs(gsi=50.0, kh=0.05)).fos:.6f}"
    assert at2.number_input(key="nb_kh").value == 0.05


def test_case_file_round_trip_and_refusals():
    inp = E.Inputs(profile=((0.0, 0.0), (-150.0, 200.0), (-260.0, 360.0)), crest_slope=5.0, plan_radius=-3000.0,
                   elong=1.5, mid_length=100.0, d_zone=10.0, D=0.7, wet=True, piezo=((0.0, 0.0), (-500.0, 250.0)),
                   crack_depth=15.0, crack_water=0.3, kh=0.05, n3y=30)
    c = CF.read_case(CF.case_json(CF.make_case(inp, name="All additions", preset="thesis_dry_05")))
    assert CF.inputs_from_case(c) == inp
    assert CF.inputs_from_case(CF.decode_link(CF.encode_link(c))) == inp
    assert CF.file_name(c) == f"ejc-case-all-additions-{c['saved'][:10]}.json"
    for bad, why in (('{"format": "x", "inputs": {}}', "not a case"), ("[", "not valid JSON"),
                     (json.dumps({"format": CF.CASE_FORMAT, "inputs": {"H": "high"}}), "H must be a number"),
                     (json.dumps({"format": CF.CASE_FORMAT, "inputs": {"n3": 2.5}}), "whole number"),
                     (json.dumps({"format": CF.CASE_FORMAT, "inputs": {"model": "<b>"}}), "choices"),
                     (json.dumps({"format": CF.CASE_FORMAT, "inputs": {"profile": [[0, 0], [1]]}}), "points")):
        with pytest.raises(ValueError, match=why):
            CF.inputs_from_case(CF.read_case(bad))
    with pytest.raises(ValueError, match="does not hold a case"):
        CF.decode_link("AAAA")


def test_report_has_every_section():
    A = E.run(E.Inputs(crack_depth=20.0, crack_water=0.5, kh=0.1, d_zone=10.0, D=0.7))
    cmp = E.compare_models(A.inputs)
    doc = RP.report_html(A, title="Report <test>", subtitle="sub", notes=[("info", "a note")],
                         studies=[("EJC against GHB", RP.compare_table(cmp, A.inputs.nu))])
    assert doc.startswith("<!doctype html>") and "Report &lt;test&gt;" in doc and "<b>" not in doc.split("<main")[0]
    assert f"{A.r3.fos:.6f}" in doc and f"{A.r2.fos:.6f}" in doc and "<svg" in doc and "tension crack" in doc
    for text in ("Analysis", "Inputs", "Studies on these inputs", "Tension crack", "Seismic k_h", "within T = 10 m",
                 "D applies within T of the face only", "a note", f"{cmp['models']['ghb']['fos3']:.6f}"):
        assert text in doc, text
    assert RP.caption_text(A).endswith("tension crack 20 m deep, 50% full of water; seismic k_h = 0.1.")
