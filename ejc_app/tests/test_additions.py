"""Additions (all off by default): slope profile, crest slope, walls curved in plan, ellipsoidal slip surfaces with a
cylindrical middle, the disturbance zone, the piezometric line, the water-filled tension crack, pseudo-static seismic
load, and the studies (EJC against GHB, sweeps, back-analysis, Monte Carlo)."""
import math

import numpy as np
import pytest

import engine as E
import presets as PR
from cases import MAIN, WET, r2, r3

BAS = PR.BY_KEY["baskoyak_site_ghb"].inputs
SLIDE = PR.BY_KEY["baskoyak_slide_ghb"].inputs
REAL = ((0.0, 0.0), (-11.9, 8.5), (-14.1, 8.5), (-30.48, 18.99))      # Baskoyak real section (presets_baskoyak.json)
CREST = math.degrees(math.atan((20.68 - 18.99) / (49.3 - 30.48)))      # ground behind the crest, 5.13 deg
DE = PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs                         # Mohr-Coulomb, base circle


# ---- geometry ----------------------------------------------------------------------------------------------------
def test_defaults_are_off():
    d = E.Inputs()
    assert (d.profile, d.crest_slope, d.plan_radius, d.elong, d.mid_length, d.d_zone, d.piezo, d.crack_depth,
            d.crack_water, d.kh) == (None, 0.0, None, 1.0, 0.0, None, None, None, 0.0, 0.0)
    m = r3(MAIN).mesh
    assert m.zone is None and m.crack is None and m.ext is None


def test_profile_matches_the_independent_code_on_the_real_section():
    """qa/qa_baskoyak.py part 4: an independent Bishop/GHB code on the traced real section (bench, ground rising
    behind the crest), failure circle, 2000 slices: 1.0590378445227162 (SLIDE's inputs), 1.0506711619962341 (site)."""
    for inp, ref in ((SLIDE, 1.0590378445227162), (BAS, 1.0506711619962341)):
        r = E.run_2d(inp.with_(profile=REAL, crest_slope=CREST, n2=2000))
        assert r.accepted and abs(r.fos - ref) <= 1e-12 * ref


def test_planar_profile_and_flat_crest_give_the_planar_slope():
    tb = math.tan(math.radians(MAIN.beta))
    planar = MAIN.with_(profile=((0.0, 0.0), (-MAIN.H / tb, MAIN.H)))
    for run in (E.run_3d, E.run_2d):
        assert run(planar).fos == pytest.approx(run(MAIN).fos, rel=1e-12)
    xs = np.linspace(-800.0, 50.0, 401)
    assert np.allclose(E.ground(planar, xs), E.ground_z(xs, MAIN.H, MAIN.beta), rtol=0, atol=1e-9)
    assert np.array_equal(E.ground(MAIN.with_(crest_slope=0.0), xs), E.ground_z(xs, MAIN.H, MAIN.beta))


def test_ground_and_search_points_on_a_profile():
    inp = BAS.with_(profile=REAL, crest_slope=CREST)
    assert float(E.ground(inp, 1.0)) == 0.0 and float(E.ground(inp, -13.0)) == pytest.approx(8.5)
    assert float(E.ground(inp, -40.0)) == pytest.approx(18.99 + (40.0 - 30.48) * math.tan(math.radians(CREST)))
    assert E.face_height(inp) == 18.99 and E.ground_point(inp, 18.99) == pytest.approx((-30.48, 18.99))
    assert E.ground_point(inp, 0.0) == (0.0, 0.0) and E.ground_point(inp, -2.0) == (2.0, 0.0)
    x, z = E.ground_point(inp, 25.0)                                     # behind the crest edge
    assert x == pytest.approx(-30.48 - 6.01) and z == pytest.approx(float(E.ground(inp, x)))


def test_profile_errors():
    for prof, text in ((((1.0, 0.0), (-10.0, 5.0)), "starts at the toe"),
                       (((0.0, 0.0), (-10.0, 5.0), (-5.0, 8.0)), "x must decrease"),
                       (((0.0, 0.0), (-10.0, 5.0), (-20.0, 4.0)), "cannot descend"),
                       (((0.0, 0.0),), "at least two points")):
        with pytest.raises(ValueError, match=text):
            E.run_2d(MAIN.with_(profile=prof))


def test_walls_curved_in_plan():
    """The central section (2D) is unchanged; on the thesis sphere a concave wall, which surrounds the mass, raises
    the 3D FoS the more the tighter it is; a convex wall with its axis inside the mass is refused."""
    assert r2(MAIN.with_(plan_radius=800.0)).fos == r2(MAIN).fos
    f = [r3(MAIN.with_(n3=80, plan_radius=rp)).fos for rp in (None, 3000.0, 1500.0, 800.0)]
    assert f[0] < f[1] < f[2] < f[3]
    assert float(E.plan_x(MAIN.with_(plan_radius=800.0), -100.0, 0.0)) == pytest.approx(-100.0)
    assert float(E.plan_x(MAIN.with_(plan_radius=-800.0), 0.0, 300.0)) == pytest.approx(math.hypot(800, 300) - 800)
    with pytest.raises(ValueError, match="convex wall needs its plan axis behind"):
        E.run_3d(MAIN.with_(plan_radius=-200.0))


# ---- slip surface shape ---------------------------------------------------------------------------------------------
def test_ellipsoid_and_cylindrical_middle():
    """lambda = 1 and L = 0 are the sphere; the 2D section is the circle for any shape."""
    assert r3(MAIN.with_(elong=1.0, mid_length=0.0)).fos == r3(MAIN).fos
    assert r2(MAIN.with_(elong=3.0, mid_length=500.0)).fos == r2(MAIN).fos
    m = E.build_mesh(MAIN.with_(elong=2.0, n3y=40), 3, 20)
    assert m.dy == 2.0 * MAIN.R / 40
    n = np.sqrt(m.n_x ** 2 + m.n_y ** 2 + m.g ** 2)
    assert np.allclose(n, 1.0) and np.allclose(m.rho, np.sqrt(MAIN.R ** 2 - (m.y / 2.0) ** 2))


def test_a_long_slide_tends_to_the_2d_result():
    """A cylinder with spherical ends: as its middle part grows, the 3D FoS tends to the 2D FoS of the same section
    (the same x grid, n3 = n2), dry and wet; at L/2R = 50 within the 0.5 % threshold of eq. 4.17."""
    for inp in (MAIN, WET):
        f2 = r2(inp.with_(n2=40)).fos
        dev = [abs(r3(inp.with_(n3=40, n3y=80, mid_length=L)).fos / f2 - 1) for L in (725.0, 2900.0, 14500.0, 72500.0)]
        assert dev[0] > dev[1] > dev[2] > dev[3] and dev[3] < 0.005


# ---- disturbance zone ------------------------------------------------------------------------------------------------
def test_disturbance_zone_limits():
    """T = 0: no base disturbed, the D = 0 result; T beyond every base: the D-everywhere result, both exactly."""
    g = MAIN.with_(model="ghb", D=0.7, a_override=None)
    assert r3(g.with_(d_zone=0.0)).fos == r3(g.with_(D=0.0)).fos
    assert r3(g.with_(d_zone=1e6)).fos == r3(g).fos
    mid = r3(g.with_(d_zone=10.0))
    assert r3(g).fos < mid.fos < r3(g.with_(D=0.0)).fos and 0 < mid.mesh.zone.sum() < mid.mesh.count
    assert mid.projection_limits is not None and mid.projection_limits[0] > mid.projection_limits[1]
    m = mid.mesh
    assert np.array_equal(m.zone, E.face_distance(g, m.xsec, m.zb) <= 10.0)
    mc = MAIN.with_(model="mc", c_mc=50.0, phi_mc=35.0)
    assert r3(mc.with_(d_zone=10.0)).fos == r3(mc).fos                    # Mohr-Coulomb has no D


def test_face_distance():
    tb = math.tan(math.radians(MAIN.beta))
    assert float(E.face_distance(MAIN, 0.0, -10.0)) == pytest.approx(10.0)          # below the toe
    x, z = -100.0 / tb, 100.0                                                       # a point on the face, moved
    nx, nz = -math.sin(math.radians(MAIN.beta)), -math.cos(math.radians(MAIN.beta))  # 20 m into the slope
    assert float(E.face_distance(MAIN, x + 20 * nx, z + 20 * nz)) == pytest.approx(20.0)


# ---- water -----------------------------------------------------------------------------------------------------------
def test_piezometric_line_equal_to_the_water_table():
    pz = ((0.0, 0.0), (-360.0 / math.tan(math.radians(50.0)), 360.0), (-2000.0, 360.0))
    for run in (E.run_3d, E.run_2d):
        assert run(WET.with_(piezo=pz)).fos == pytest.approx(run(WET).fos, rel=1e-12)
    with pytest.raises(ValueError, match="piezometric line runs back"):
        E.run_2d(WET.with_(piezo=((0.0, 0.0), (10.0, 5.0))))


def _bishop_2d(inp, ext):
    """Bishop's simplified method written out (textbook form), Mohr-Coulomb, on the engine's slices: F = sum[(c b +
    (W - u b) tan phi)/m_a] / [sum W sin a + external moment / R], m_a = cos a + sin a tan phi / F."""
    m = E.build_mesh(inp, 2, inp.n2)
    sa, ca, b = m.kap, m.g, m.dx
    tp, c = math.tan(math.radians(inp.phi_mc)), inp.c_mc
    den = float(np.sum(m.w * sa)) + ext(m) / inp.R
    F = 1.0
    for _ in range(500):
        Fn = float(np.sum((c * b + (m.w - m.u * b) * tp) / (ca + sa * tp / F))) / den
        if abs(Fn - F) < 1e-15 * Fn:
            break
        F = Fn
    return Fn


def test_seismic_load_against_bishop_with_a_seismic_term():
    """2D: the engine against Bishop's method with the pseudo-static moment sum k_h W (z_c - z_g) about the centre."""
    for kh in (0.05, 0.15):
        inp = DE.with_(kh=kh)
        ref = _bishop_2d(inp, lambda m: float(np.sum(kh * m.w * (inp.zc - 0.5 * (m.zb + m.zt)))))
        assert E.run_2d(inp).fos == pytest.approx(ref, rel=1e-10)
    f = [r3(MAIN.with_(kh=kh)).fos for kh in (0.0, 0.05, 0.1, 0.2)]
    assert f[0] == r3(MAIN).fos and f[0] > f[1] > f[2] > f[3]


def test_tension_crack_against_bishop_with_the_crack_water_force():
    """The crack sits where the circle is crack_depth below the ground, the mass ends there, and water in it adds
    P (z_c - z_P) to the driving moment: P = gamma_w h_w^2/2 at h_w/3 above the crack bottom."""
    inp = DE.with_(crack_depth=0.5, crack_water=0.6)
    xt = E.crack_x(inp)
    depth = float(E.ground(inp, xt) - E.sphere_z(xt, 0.0, inp.xc, inp.zc, inp.R))
    assert depth == pytest.approx(0.5, abs=1e-9)
    hw = 0.6 * depth
    zbot = float(E.sphere_z(xt, 0.0, inp.xc, inp.zc, inp.R))
    P, e = 0.5 * E.GAMMA_W * hw * hw, inp.zc - (zbot + hw / 3)
    m = E.build_mesh(inp, 2, inp.n2)
    assert m.x.min() >= xt and m.crack["P"] == pytest.approx(P, rel=1e-12)
    assert E.run_2d(inp).fos == pytest.approx(_bishop_2d(inp, lambda m: P * e), rel=1e-10)
    dry = E.run_2d(inp.with_(crack_water=0.0)).fos
    assert E.run_2d(inp).fos < dry
    with pytest.raises(ValueError, match="deeper than the sliding mass"):
        E.run_2d(MAIN.with_(crack_depth=60.0))
    assert E.vertical_faces(MAIN.with_(xc=25.0, zc=250.0, R=251.0, crack_depth=5.0))["back"] == 0.0


def test_mass_forces_balance_with_seismic_load_and_crack_water():
    for inp in (MAIN.with_(kh=0.1), MAIN.with_(kh=0.1, form="FM"), MAIN.with_(crack_depth=20.0, crack_water=1.0),
                DE.with_(kh=0.1, crack_depth=0.5, crack_water=0.5)):
        for r in (r3(inp), r2(inp)):
            mf = E.mass_forces(inp, r)
            assert abs(mf["method"]["residual"]) <= 1e-9 and mf["vertical_max_rel"] <= 1e-9
            if r.mesh.ext is not None:
                assert mf["totals"]["E"]["x"] == pytest.approx(mf["mult"] * inp.kh * r.mesh.w.sum(), rel=1e-12)


def test_all_additions_together_in_both_dimensions():
    inp = BAS.with_(profile=REAL, crest_slope=5.0, plan_radius=300.0, elong=1.5, mid_length=20.0, d_zone=3.0, kh=0.05,
                    crack_depth=1.0, crack_water=0.5, n3y=25)
    a, b = r3(inp), r2(inp)
    assert a.accepted and b.accepted
    for r in (a, b):
        assert abs(E.mass_forces(inp, r)["method"]["residual"]) <= 1e-9


# ---- studies ---------------------------------------------------------------------------------------------------------
def test_compare_models():
    c = E.compare_models(MAIN)
    assert c["models"]["ejc"]["fos3"] == r3(MAIN).fos and c["models"]["ghb"]["fos2"] == r2(MAIN.with_(model="ghb")).fos
    assert c["ratio"]["ejc_ghb3"] == pytest.approx(c["models"]["ejc"]["fos3"] / c["models"]["ghb"]["fos3"])
    assert c["models"]["ejc"]["ratio"] == pytest.approx(r3(MAIN).fos / r2(MAIN).fos)


def test_parameter_sweep():
    rows = E.parameter_sweep(MAIN, "gsi", [30.0, 42.0], ("ghb",), (2,))
    assert [r["value"] for r in rows] == [30.0, 42.0]
    assert rows[1]["ghb"]["fos2"] == r2(MAIN.with_(model="ghb")).fos and rows[0]["ghb"]["fos2"] < rows[1]["ghb"]["fos2"]
    with pytest.raises(ValueError, match="unknown input"):
        E.parameter_sweep(MAIN, "nope", [1.0])


def test_back_analysis_of_d_on_the_baskoyak_failure_circle():
    """The disturbance factor that gives F = 1 on the failure circle, GHB 2D: 0.721 in qa/qa_baskoyak.py (a cubic fit
    of F against D); bisection here."""
    r = E.back_analysis(BAS, "D", 0.5, 0.9, 1.0, dim=2)
    assert r["message"] == "" and abs(r["fos"] - 1.0) < 1e-6 and 0.7205 < r["value"] < 0.7215
    bad = E.back_analysis(BAS, "D", 0.0, 0.3, 1.0, dim=2)
    assert bad["value"] is None and "does not cross" in bad["message"]


def test_monte_carlo_is_reproducible():
    var = [("gsi", 42.0, 5.0, 10.0, 100.0), ("sigma_ci", 77.7, 15.0, 1.0, 300.0)]
    a = E.monte_carlo(MAIN.with_(model="ghb"), var, 8, seed=3, dims=(2,))
    b = E.monte_carlo(MAIN.with_(model="ghb"), var, 8, seed=3, dims=(2,))
    assert a["samples"] == b["samples"]
    f = [s["ghb2"] for s in a["samples"]]
    st = a["stats"]["ghb2"]
    assert st["accepted"] == 8 and st["mean"] == pytest.approx(np.mean(f)) and st["sd"] == pytest.approx(np.std(f, ddof=1))
    assert all(10 <= s["values"]["gsi"] <= 100 for s in a["samples"]) and st["pof"] == 0.0
    r = E.mulberry32(1)
    u = [r() for _ in range(1000)]
    assert all(0 <= v < 1 for v in u) and abs(np.mean(u) - 0.5) < 0.05
