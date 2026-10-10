import json

import pandas as pd
import pytest

from slope3d import Config, io, run_case
from slope3d.cli import main
from slope3d.selftest import run_all


def test_validate_command_all_pass(capsys):
    assert main(["validate"]) == 0
    out = capsys.readouterr().out
    assert "ALL PASS" in out and "FAIL " not in out


def test_selftest_programmatic():
    assert all(r.ok for r in run_all())


def test_run_exit_codes_and_withheld_output(capsys, tmp_path):
    assert main(["run", "cases/chile_thesis_default.yaml"]) == 1
    out = capsys.readouterr().out
    assert "FoS withheld" in out and out.count("NOT ACCEPTED") == 2
    assert main(["run", "cases/chile_dry_conventional_hb.yaml", "--mesh3d", "10", "--mesh2d", "20",
                 "--out", str(tmp_path), "--format", "csv,xlsx,json"]) == 0
    assert (tmp_path / "chile_dry_conventional_hb.xlsx").exists()
    data = json.loads((tmp_path / "chile_dry_conventional_hb.json").read_text())
    assert data["3d"]["status"] == "ACCEPTED" and data["3d"]["fos"] > 0


def test_withheld_runs_never_export_a_number(tmp_path):
    case = run_case(Config(), 20, 40)
    paths = io.export_csv(case, tmp_path)
    summary = pd.read_csv([p for p in paths if p.name.endswith("summary.csv")][0])
    assert summary.loc[summary["item"] == "3D FoS", "value"].iloc[0] == "WITHHELD"
    assert summary.loc[summary["item"] == "2D FoS", "value"].iloc[0] == "WITHHELD"
    assert io.result_to_dict(case)["3d"]["fos"] is None
    cols = io.elements_table(case.run3d, case.config)
    assert list(cols.columns) == io.ELEMENT_COLUMNS and len(cols) == 76
    assert (cols["actual final normal"] == "NOT EVALUATED").all()


def test_case_round_trip(tmp_path):
    cfg = Config(closure=3, a_override=None, water_angle=50.0, water_height=360.0).with_face_angle(55)
    for name in ("c.yaml", "c.json"):
        io.save_case(tmp_path / name, cfg, 30, 60)
        cfg2, m3, m2 = io.load_case(tmp_path / name)
        assert cfg2 == cfg and (m3, m2) == (30, 60)


def test_unknown_input_is_an_error(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("closure: 3\nradiuss: 5\n")
    with pytest.raises(ValueError):
        io.load_case(p)


def test_element_table_matches_workbook_fields():
    case = run_case(Config(closure=3), 20, 40)
    t = io.elements_table(case.run3d, case.config)
    assert len(t.columns) == 23 and t["status"].eq("ADMISSIBLE").all()
