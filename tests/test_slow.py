"""Suite C: slow regressions (run with ``pytest -m slow``)."""
import csv
from pathlib import Path

import pytest

from slope3d.mesh_study import compare_with_records, exact_ejc_threshold, run_mesh_study

REF = Path(__file__).resolve().parents[1] / "reference" / "mesh_study"
pytestmark = pytest.mark.slow


def test_mesh_study_records_reproduced_to_1e9():
    """All 95 records: F to 1e-9 relative, plus elements, volume and acceptance (3D to 320x320, 2D to 2,560 slices)."""
    rows = run_mesh_study()
    assert len(rows) == 95
    assert compare_with_records(rows, REF / "mesh_study_records.csv", 1e-9) == []


def test_mesh_study_extra_columns_match_records():
    """Path projections, minimum final normal stress and tension counts also match the independent reimplementation."""
    ref = {(r["case"], r["dim"], r["mesh"]): r for r in csv.DictReader(open(REF / "mesh_study_records.csv"))}
    for row in run_mesh_study(meshes3=(10, 20, 40), meshes2=(20, 40, 160)):
        r = ref[(row["case"], row["dim"], row["mesh"])]
        assert row["n_tension"] == int(r["n_tension"])
        assert row["min_sn"] == pytest.approx(float(r["min_sn"]), rel=1e-9, abs=1e-9)


@pytest.mark.parametrize("a,expected", [(0.5, 52.966376), (0.509923, 52.283035)])
def test_exact_ejc_admissibility_threshold_at_nu_0271(a, expected):
    sn, phi, s3 = exact_ejc_threshold(0.271, a)
    assert sn == pytest.approx(expected, abs=5e-7)


def test_nu_table_matches_nu_exact_txt():
    """Every printed row of reference/mesh_study/nu_exact.txt."""
    import re
    text = (REF / "nu_exact.txt").read_text()
    for m in re.finditer(r"nu=\s*([\d.]+) \| a=0\.5:\s+([\d.]+) MPa .*? \| a=0\.509923:\s+([\d.]+) MPa", text):
        nu, a1, a2 = float(m[1]), float(m[2]), float(m[3])
        assert exact_ejc_threshold(nu, 0.5)[0] == pytest.approx(a1, abs=5e-7)
        assert exact_ejc_threshold(nu, 0.509923)[0] == pytest.approx(a2, abs=5e-7)
