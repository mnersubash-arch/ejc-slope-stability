"""Suite B: engine regressions (3D 20x20, 2D 40 slices), rejection and robustness tests.

Expected values are those of Thesis_Columns.bas ``Thesis_RunEngineQA`` (also the workbook's recorded values).
"""
import math

import numpy as np
import pytest

from slope3d import Config, Status, run_case, solve_model, validate_final
from slope3d.geometry import build_elements, geometry_issue

REL = 1e-9

BENCH = [  # name, config overrides (closure 3 comparator unless stated), F3D, F2D
    ("dry a=0.5", {}, 2.5536684200376, 2.31025502757599),
    ("wet a=0.5", dict(water_angle=50.0, water_height=360.0), 1.20012184467095, 1.03587083245705),
    ("dry a=0.509923", dict(a_override=0.509923), 2.47603808311863, 2.24206326532995),
    ("wet a=0.509923", dict(a_override=0.509923, water_angle=50.0, water_height=360.0), 1.15842745018062, 1.00046615780691),
    ("published circle", dict(xc=207.28, zc=586.53, radius=622.08), 2.14193739499583, 2.00387341089056),
    ("MC c=100 phi=0", dict(strength=0, c_mc=100.0, phi_mc=0.0), 0.255354533635107, 0.190263799071471),
]

# actual values recorded by the VBA engine in Excel ('Native Engine QA'!C5:E10): the port should match these tighter
VBA_ACTUAL = {
    "dry a=0.5": (2.5536684200348696, 2.310255027573278),
    "wet a=0.5": (1.2001218446678463, 1.0358708324538346),
    "dry a=0.509923": (2.4760380831157867, 2.2420632653272645),
    "wet a=0.509923": (1.1584274501774219, 1.0004661578100102),
    "published circle": (2.1419373949945077, 2.003873410889333),
    "MC c=100 phi=0": (0.25535453363510774, 0.19026379907147156),
}


def cfg_for(over, **extra):
    return Config(closure=3, **{**over, **extra})


@pytest.mark.parametrize("name,over,f3,f2", BENCH, ids=[b[0] for b in BENCH])
def test_benchmark_regression(name, over, f3, f2):
    r = run_case(cfg_for(over), 20, 40)
    for run, expected in ((r.run3d, f3), (r.run2d, f2)):
        assert run.result.status is Status.ACCEPTED, run.result.reason
        assert run.result.fos == pytest.approx(expected, rel=REL)


@pytest.mark.parametrize("name,over,f3,f2", BENCH, ids=[b[0] for b in BENCH])
def test_matches_vba_engine_actuals(name, over, f3, f2):
    """The VBA engine's own recorded output (Excel 16.0) is reproduced to better than 1e-12."""
    r = run_case(cfg_for(over), 20, 40)
    v3, v2 = VBA_ACTUAL[name]
    assert r.run3d.result.fos == pytest.approx(v3, rel=1e-12)
    assert r.run2d.result.fos == pytest.approx(v2, rel=1e-12)


def test_thesis_closure_chile_is_rejected_and_withheld():
    r = run_case(Config(), 20, 40)                      # defaults: closure 1, Chile
    for run, n in ((r.run3d, 76), (r.run2d, 30)):
        res = run.result
        assert res.status is Status.NOT_ACCEPTED
        assert res.fos is None and res.fos_text == "WITHHELD"
        assert res.active == n
        assert res.reason == f"{n} of {n} elements have no admissible eq.4.21 strength; no FoS."
        assert not res.evaluated
    assert r.ratio is None


def test_moscarda_closure_is_not_accepted_dry_and_wet():
    for over in ({}, dict(a_override=0.509923, water_angle=50.0, water_height=360.0)):
        r = run_case(Config(closure=2, **over), 20, 40)
        for run in (r.run3d, r.run2d):
            assert run.result.status is Status.NOT_ACCEPTED and run.result.fos is None


def test_checks_off_gives_diagnostic_only_and_withholds_f():
    """VBA: the acceptance switch only changes the final label of an otherwise valid run."""
    r = run_case(Config(closure=3, checks=0), 20, 40)
    for run, f in ((r.run3d, 2.5536684200376), (r.run2d, 2.31025502757599)):
        assert run.result.status is Status.DIAGNOSTIC_ONLY
        assert run.result.fos is None
        assert run.result.internal_f == pytest.approx(f, rel=REL)


@pytest.mark.xfail(strict=True, reason="OPEN QUESTION for the author: these numbers come from the workbook's worksheet "
                   "engine (Native Workbook QA test 8, P19 = 0 uses unchecked eq.4.21 roots). The VBA engine "
                   "(Thesis_Columns) never accepts an inadmissible strength, so it returns NOT ACCEPTED. See docs/SPEC.md 0.")
def test_moscarda_wet_checks_off_diagnostic_values():
    cfg = Config(closure=2, a_override=0.509923, water_angle=50.0, water_height=360.0, checks=0)
    r = run_case(cfg, 20, 40)
    assert r.run3d.result.status is Status.DIAGNOSTIC_ONLY and r.run2d.result.status is Status.DIAGNOSTIC_ONLY
    assert r.run3d.result.internal_f == pytest.approx(1.04327747791523, rel=REL)
    assert r.run2d.result.internal_f == pytest.approx(0.90176054371878, rel=REL)


def test_perturbing_final_normal_stress_is_rejected():
    cfg = cfg_for(dict(a_override=0.509923, water_angle=50.0, water_height=360.0))
    run = solve_model(cfg, 3, 20)
    res, els = run.result, run.elements
    assert res.status is Status.ACCEPTED
    assert validate_final(cfg, els, res.internal_f, res, rebuild=False)        # untouched actual state still passes
    els.final_n[0] += 1.0                                                        # +1 kPa on one column
    assert not validate_final(cfg, els, res.internal_f, res, rebuild=False)
    assert res.admissible < res.active


@pytest.mark.parametrize("over", [
    dict(xc=-1.45, zc=3.24, radius=4.46, height=2.1, slope=1 / 3.5, gamma=19.0, strength=0, c_mc=7.2, phi_mc=30.0),
    dict(xc=207.28, zc=300.0, radius=622.08),
], ids=["original MC geometry", "cropped geometry"])
def test_unsupported_geometry_is_refused(over):
    r = run_case(Config(closure=3, **over), 20, 40)
    for run in (r.run3d, r.run2d):
        res = run.result
        assert res.status is Status.NOT_ACCEPTED
        assert res.fos is None and res.internal_f is None
        assert res.active == 0 and not res.evaluated
        assert res.reason.startswith("UNSUPPORTED GEOMETRY:")
        assert run.elements is None


@pytest.mark.parametrize("over,needle", [
    (dict(radius=float("nan")), "finite"),
    (dict(zc=-1.0), "R > zc > 0"),
    (dict(xc=800.0), "0 < xc < R"),
    (dict(height=0.0), "R > zc > 0"),
    (dict(slope=0.0), "R > zc > 0"),
    (dict(gamma=-1.0), "R > zc > 0"),
])
def test_geometry_guard_rules(over, needle):
    issue = geometry_issue(Config(**over))
    assert issue.startswith("UNSUPPORTED GEOMETRY:") and needle in issue


def test_geometry_guard_toe_tolerance():
    # the published circle's toe mismatch (2.62 mm) is allowed; anything beyond 1e-5 R is refused
    assert geometry_issue(Config(xc=207.28, zc=586.53, radius=622.08)) == ""


def test_element_counts_and_volume():
    els3, vol3 = build_elements(Config(), 3, 20)
    els2, vol2 = build_elements(Config(), 2, 40)
    assert (els3.n, els2.n) == (76, 30)
    assert vol3 > 0 and vol2 > 0
    # derived quantities
    assert np.allclose(els3.area * els3.g, els3.ah)
    assert np.allclose(els3.weight, 27.0 * els3.h * els3.ah)


@pytest.mark.parametrize("dim,n", [(3, 1), (3, 401), (2, 1), (2, 20001), (4, 20)])
def test_mesh_limits(dim, n):
    res = solve_model(Config(closure=3), dim, n).result
    assert res.status is Status.NOT_ACCEPTED and res.fos is None and res.reason.startswith("INVALID MESH")


def test_invalid_inputs_do_not_produce_numbers():
    res = solve_model(Config(nu=0.7), 3, 20).result
    assert res.status is Status.NOT_ACCEPTED and res.fos is None and res.reason.startswith("INVALID INPUT")
    res = solve_model(Config(closure=3, tol_f=1e-3), 2, 40).result
    assert res.fos is None


def test_ratio_and_moment_diagnostic():
    r = run_case(cfg_for({}), 20, 40)
    assert r.ratio == pytest.approx(2.5536684200376 / 2.31025502757599, rel=1e-9)
    assert r.run3d.result.moment_diagnostic is not None
    r = run_case(Config(), 20, 40)
    assert r.run3d.result.moment_diagnostic is None


def test_thesis_closure_label_rules():
    assert run_case(Config(), 4, 4).label == "Thesis closure (3D EJC)"
    assert "Conventional HB comparator (3D geometry, 2D strength)" in run_case(Config(closure=3), 4, 4).label
    assert not Config(closure=1, checks=0).is_thesis_closure
    assert not Config(closure=3).is_thesis_closure
