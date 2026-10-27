"""Critical-surface search (engine.critical_search; QA/QC finding R1)."""
import math

import numpy as np
import pytest

import engine as E
import presets as PR

# Exact critical toe circles of phi = 0 slopes (Taylor's problem), by adaptive quadrature in qa/qa_analytic.py:
# H = 10 m, gamma = 27 kN/m3, c = 100 kPa; stability numbers 0.2191 and 0.1906 (Taylor's chart: 0.219 and 0.191).
TAYLOR = {75.0: 1.6905691599108592, 60.0: 1.943466058585107}


def taylor(beta, **kw):
    return E.Inputs(H=10.0, beta=beta, gamma=27.0, model="mc", c_mc=100.0, phi_mc=0.0, n2=640, **kw)


def test_circle_from_exits():
    for pf, pb, t in (((-1.0, 1.2), (-12.0, 10.0), 0.4), ((0.0, 0.0), (-5.0, 3.0), 1.0), ((2.0, 0.0), (-8.0, 4.0), 0.05)):
        c = E.circle_from_exits(pf, pb, t)
        for x, z in (pf, pb):
            assert math.hypot(x - c["xc"], z - c["zc"]) == pytest.approx(c["R"], rel=1e-12)
            assert z <= c["zc"]                                        # both exits on the lower half
        assert c["zc"] - pb[1] == pytest.approx((1 / t - 1) * (pf[0] - pb[0]) / 2, rel=1e-9, abs=1e-12)
    c = E.circle_from_exits((0.0, 0.0), (-5.0, 3.0), 1.0)              # t = 1: the arc is vertical at the back exit
    assert c["xc"] - c["R"] == pytest.approx(-5.0, abs=1e-12) and c["zc"] == pytest.approx(3.0, abs=1e-12)


@pytest.mark.parametrize("base", [False, True])
def test_search_surfaces_meet_the_ground_at_both_exits(base):
    """Every surface of the search domain passes through its two exits on the ground surface and meets the ground at
    the back exit before its back extreme x_c - R, so the mass has no vertical back face."""
    inp = E.Inputs(base_circle=base)
    lo = -2.0 if base else 0.0
    for a in np.linspace(lo, 0.99, 7):
        for b in np.linspace(max(a, 0.0) + 0.01, 3.0, 9):
            for t in np.linspace(0.02, 1.0, 7):
                s = E.search_surface(inp, (a, b, t))
                for x, z in (s["front"], s["back"]):
                    assert float(E.ground_z(x, inp.H, inp.beta)) == pytest.approx(z, abs=1e-9 * inp.H)
                    assert math.hypot(x - s["xc"], z - s["zc"]) == pytest.approx(s["R"], abs=1e-9 * inp.H)
                    assert z <= s["zc"] + 1e-9 * inp.H                 # on the lower half of the circle
                assert float(E.ground_z(s["xc"] - s["R"], inp.H, inp.beta)) <= s["zc"] + 1e-9 * inp.H


def test_vertical_back_circle_is_reported_and_not_searched():
    """A circle whose centre lies below the crest ends at its back extreme below the ground: the engine analyses it
    with a vertical back face that carries no force (here 2.8 m, exact FoS 1.6245, below Taylor's 1.6906). The search
    reports it for the prescribed surface and does not try such circles."""
    res = E.critical_search(taylor(75.0, xc=0.40, zc=7.20, R=7.21))
    P = res["prescribed"]
    assert P["vertical_back"] == pytest.approx(10.0 - 7.20) and P["fos2"] < TAYLOR[75.0]
    assert res["r2"]["fos"] > P["fos2"]


@pytest.mark.parametrize("beta", [75.0, 60.0])
def test_search_finds_taylors_critical_toe_circle(beta):
    """On 640 slices the critical FoS lies within the 0.5 % screening threshold (eq. 4.17) of the exact value, on a
    toe circle (front exit within 1 % of H of the toe)."""
    r = E.critical_search(taylor(beta))["r2"]
    assert abs(r["fos"] / TAYLOR[beta] - 1) < 0.005
    assert r["front"][1] < 0.01 * 10.0
    assert r["check"]["n"] == 2560 and abs(r["check"]["fos"] / TAYLOR[beta] - 1) < 0.005


def test_back_limit_is_reported():
    """Taylor 75 deg with the back exit limited to 1 H behind the crest: the critical circle exits 0.85 H behind it, so
    the search stops at the limit and says so."""
    r = E.critical_search(taylor(75.0), back=1.0)["r2"]
    assert any(e["why"] == "back limit" for e in r["edges"])
    assert r["back"][0] == pytest.approx(-10.0 / math.tan(math.radians(75.0)) - 10.0, abs=0.01)


def test_search_result_is_the_lowest_fos_tried():
    """Base-circle option (D'Esterre 2024, DE-1.5-c1-ru0), 2D and 3D: the critical FoS is the lowest accepted FoS of all
    the surfaces tried, at most the published circle's, and every surface tried is counted once."""
    inp = PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs
    res = E.critical_search(inp, dims=(2, 3))
    for d in (2, 3):
        r = res[f"r{d}"]
        fos = [v for v in r["cloud"]["fos"] if v is not None]
        assert r["fos"] == min(fos) and r["fos"] <= res["prescribed"][f"fos{d}"]
        assert r["evaluations"] == sum(r["counts"].values()) == len(r["cloud"]["xc"])
        assert r["counts"]["solved"] == len(fos) and r["counts"]["skipped"] == sum(r["skipped_why"].values())
    assert res["r2"]["fos"] == pytest.approx(2.605316267142031, rel=1e-12)      # as webapp/engine.js (verify_engine.js)
    assert res["r3"]["fos"] == pytest.approx(2.9652854773051613, rel=1e-12)


def test_single_mass():
    """With the base-circle grid the thesis circle would leave two masses (it also dips below the pit floor ahead of
    the toe); a D'Esterre circle leaves one."""
    assert not E.single_mass(E.Inputs(base_circle=True))
    assert E.single_mass(PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs)


@pytest.mark.slow
def test_search_thesis_slope_ghb_dry():
    """The thesis slope with GHB, dry (about a minute): the 2D critical FoS is 13 % below the main sphere's (QA/QC
    finding R1); Carranza-Torres & Hormazabal (2020) report 2.01 for their critical circle on this slope."""
    res = E.critical_search(E.Inputs(model="ghb"))
    r = res["r2"]
    assert r["fos"] == pytest.approx(2.000179246888455, rel=1e-12)              # as webapp/engine.js
    assert r["check"]["fos"] == pytest.approx(2.0065732918732504, rel=1e-12)    # 160 slices
    assert r["fos"] / res["prescribed"]["fos2"] - 1 == pytest.approx(-0.134, abs=0.001)
