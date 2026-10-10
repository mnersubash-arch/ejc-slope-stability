"""Slip-surface exit, plots, run report and the Streamlit app (headless AppTest)."""
import re
from pathlib import Path

import pytest

from slope3d import Config, run_case
from slope3d import plotly_figs as pf
from slope3d import report
from slope3d.geometry import slip_surface_ends


def test_slip_surface_exit_matches_workbook_chart_helper():
    """Interface!F39:G39 of the workbook: (-20.40573098459213, 24.31860321039869)."""
    e = slip_surface_ends(Config())
    assert e["exit_on_face"]
    assert e["exit_x"] == pytest.approx(-20.40573098459213, abs=1e-9)
    assert e["exit_z"] == pytest.approx(24.31860321039869, abs=1e-9)
    assert e["exit_behind_toe"] == pytest.approx(20.40573098459213, abs=1e-9)


def test_slip_surface_ends_published_circle_reaches_toe_plane():
    e = slip_surface_ends(Config(xc=207.28, zc=586.53, radius=622.08))
    assert not e["exit_on_face"] and abs(e["exit_z"]) < 2e-3


def test_slip_surface_ends_unsupported():
    assert slip_surface_ends(Config(xc=-1.45, zc=3.24, radius=4.46))["supported"] is False


@pytest.mark.parametrize("cfg,accepted", [(Config(closure=3), True), (Config(), False)])
def test_report_and_figures(tmp_path, cfg, accepted):
    case = run_case(cfg, 20, 40)
    for fn in (report.fig_section, report.fig_plan, report.fig_surface3d, report.fig_convergence, report.fig_strength):
        report.plt.close(fn(case))
    for fn in (pf.section, pf.plan, pf.surface3d, pf.convergence, pf.strength):
        assert fn(case).to_plotly_json()["data"] is not None
    paths = report.write_report(case, tmp_path, "t")
    assert {p.suffix for p in paths} == {".pdf", ".html"} and all(p.stat().st_size > 5000 for p in paths)
    html = (tmp_path / "t_report.html").read_text()
    assert "exits the face" in html or "slip surface exit" in html
    assert "tension elements" in html
    if accepted:
        assert "F = 2.553668" in html
    else:
        assert "FoS withheld" in html and "F = " not in html.split("<h2>Summary")[0]


def test_withheld_runs_do_not_plot_f():
    case = run_case(Config(), 20, 40)
    from slope3d.plotdata import convergence_data, plan_data
    assert convergence_data(case) == {"3D": None, "2D": None}
    assert plan_data(case) is None


def test_report_via_cli(tmp_path):
    from slope3d.cli import main
    assert main(["run", "cases/chile_dry_conventional_hb.yaml", "--out", str(tmp_path), "--format", "report"]) == 0
    assert (tmp_path / "chile_dry_conventional_hb_report.pdf").exists()


APP = str(Path(__file__).resolve().parents[1] / "slope3d" / "gui" / "app.py")


def _app():
    from streamlit.testing.v1 import AppTest
    return AppTest.from_file(APP, default_timeout=180).run()


def _text(at):
    return " ".join([m.value for m in at.markdown] + [m.value for m in at.warning] + [m.value for m in at.success] + [m.value for m in at.info])


def test_gui_preset_shows_f_and_thesis_preset_withholds():
    at = _app()
    assert not at.exception
    at.sidebar.button(key="preset_0").click().run()                  # Dry, a = 0.5
    assert not at.exception
    assert [m.value for m in at.metric][:3] == ["2.5537", "2.3103", "1.1054"]
    assert "ACCEPTED" in _text(at) and "Conventional HB comparator" in _text(at)
    at.sidebar.button(key="preset_6").click().run()                  # thesis closure
    assert not at.exception
    assert [m.value for m in at.metric][:3] == ["WITHHELD", "WITHHELD", "n.a."]
    t = _text(at)
    assert "NOT ACCEPTED" in t and "FoS withheld" in t and "F = " not in t


def test_gui_has_all_six_presets_and_groups():
    at = _app()
    keys = [b.key for b in at.sidebar.button]
    assert sum(k.startswith("preset_") for k in keys) == 7
    labels = [e.label for e in at.sidebar.expander]
    assert labels == ["Geometry", "Slope", "Water", "Rock mass", "Strength model", "Closure", "Tolerances and run controls", "Mesh"]


def test_gui_checks_off_shows_labelled_diagnostic_f_on_a_grey_card():
    at = _app()
    at.sidebar.button(key="preset_0").click().run()
    at.sidebar.radio(key="in_checks").set_value(0)
    at.sidebar.button(key="run_btn").click().run()
    assert not at.exception
    cards = [m.value for m in at.markdown if "<div style='border:2px solid" in m.value]
    assert len(cards) == 2
    for c in cards:
        assert "DIAGNOSTIC ONLY" in c and "Diagnostic F (unchecked roots, not an accepted result)" in c
        assert "#6b6a66" in c and "#008300" not in c and "F = " not in c and "FoS withheld" not in c
        assert "rgba(128,128,128" in c
    assert "ACCEPTED" not in " ".join(re.sub(r"NOT ACCEPTED|DIAGNOSTIC ONLY", "", c) for c in cards)
    labels = [m.label for m in at.metric][:2]
    assert labels == ["Diagnostic F 3D (unchecked, not accepted)", "Diagnostic F 2D (unchecked, not accepted)"]
    assert [m.value for m in at.metric][:2] == ["2.5537", "2.3103"]


def test_gui_not_accepted_keeps_fos_withheld_card():
    at = _app()
    at.sidebar.button(key="preset_6").click().run()
    cards = [m.value for m in at.markdown if "<div style='border:2px solid" in m.value]
    assert all("FoS withheld" in c and "Diagnostic F" not in c for c in cards) and len(cards) == 2


def test_diagnostic_f_is_labelled_in_cli_tables_and_json(capsys, tmp_path):
    from slope3d import io
    from slope3d.cli import format_result
    case = run_case(Config(closure=3, checks=0), 20, 40)
    text = format_result(case)
    assert text.count("Diagnostic F (unchecked roots, not an accepted result): 2.") + text.count(": 2.3102") >= 2
    assert "F       : 2." not in text and "FoS withheld" not in text
    summ = io.summary_table(case)
    assert (summ["item"] == "3D FoS").any() and summ.loc[summ["item"] == "3D FoS", "value"].iloc[0] == "WITHHELD"
    assert summ["item"].str.contains("Diagnostic F").sum() == 2
    d = io.result_to_dict(case)
    assert d["3d"]["fos"] is None and d["3d"]["diagnostic_f"] == pytest.approx(2.5536684200376, rel=1e-9)
    # a NOT ACCEPTED run never carries a diagnostic value
    assert io.result_to_dict(run_case(Config(), 20, 40))["3d"]["diagnostic_f"] is None
    html = report.write_report(case, tmp_path, "d", pdf=False)[0].read_text()
    assert "Diagnostic F (unchecked roots, not an accepted result)" in html and "F = 2.55" not in html


def test_root_streamlit_app_entry_point_runs():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "streamlit_app.py"), default_timeout=180).run()
    assert not at.exception
    at.sidebar.button(key="preset_0").click().run()           # a rerun must work too (runpy re-executes the app)
    assert not at.exception and [m.value for m in at.metric][:2] == ["2.5537", "2.3103"]


def test_mesh_convergence_is_reduced_by_default_and_gated_for_slow_closures():
    at = _app()
    at.sidebar.button(key="preset_0").click().run()           # closure 3
    plots_buttons = [b.label for b in at.button]
    assert any("3D 10-40, 2D 20-320 (reduced)" in b for b in plots_buttons)
    assert at.checkbox(key="mesh_full").value is False
    at.sidebar.button(key="preset_6").click().run()           # thesis closure (scan closure)
    btn = [b for b in at.button if b.label.startswith("Run mesh convergence")][0]
    assert btn.proto.disabled and "3D 10-20, 2D 20-80" in btn.label
    assert any("slow" in w.value for w in at.warning)
    at.checkbox(key="mesh_ok").check().run()
    assert not [b for b in at.button if b.label.startswith("Run mesh convergence")][0].proto.disabled


def test_mesh_convergence_reduced_run_gives_accepted_points_only():
    from slope3d.plotdata import MESHES_REDUCED, mesh_convergence
    rows = mesh_convergence(Config(closure=3), (10, 20), (20, 40))
    assert len(rows) == 4 and all(r["F"] is not None for r in rows)
    rows = mesh_convergence(Config(), (10,), (20,))            # thesis closure: nothing plotted, nothing shown
    assert all(r["F"] is None for r in rows)
    assert MESHES_REDUCED[0][-1] <= 40 and MESHES_REDUCED[1][-1] <= 320
