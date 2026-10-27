"""Rock-mass parameters, the EJC on its closures, the tangent parameters and the Route 2 solver (brief §1.1-1.3)."""
import math

import numpy as np
import pytest

import engine as E
from cases import MAIN

ST = MAIN.strength()                          # Moscarda, nu = 0.271, a = 0.5
GHB = MAIN.with_(model="ghb").strength()
NU = 0.271


def rounds_to(value, target, decimals):
    return abs(value - target) <= 0.5 * 10.0 ** (-decimals) * (1 + 1e-12)


# ---- §1.1 ------------------------------------------------------------------------------------------------------
def test_rock_mass_parameters():
    mb, s, a = E.hoek_brown_parameters(42.0, 10.0, 0.0)
    assert rounds_to(mb, 1.260056, 6)
    assert rounds_to(s, 0.0015893, 7)
    assert rounds_to(a, 0.509923, 6)


def test_a_override():
    assert MAIN.hb[2] == 0.5
    assert MAIN.with_(a_override=None).hb[2] == pytest.approx(0.5 + (math.exp(-42 / 15) - math.exp(-20 / 3)) / 6)


# ---- §1.2 closure algebra ----------------------------------------------------------------------------------------
def test_moscarda_constants_match_brief():
    P, R, S = E.closure_constants(NU, NU)
    C, B = 1 - NU + NU ** 2, 1 + 2 * NU - 2 * NU ** 2
    assert P == pytest.approx(C, abs=1e-15) and S == pytest.approx(C, abs=1e-15) and R == pytest.approx(-B, abs=1e-15)


def test_original_constants_match_brief():
    P, R, S = E.closure_constants(-1.0, 1.0 / NU)
    assert (P, R) == pytest.approx((3.0, -3.0 / NU)) and S == pytest.approx(1 - 1 / NU + 1 / NU ** 2)


@pytest.mark.parametrize("model", ["ejc_moscarda", "ejc_original", "ghb", "ejc_sigma2_eq_sigma3"])
def test_states_satisfy_unsquared_and_squared_forms(model):
    st = MAIN.with_(model=model, nu=0.05 if model == "ejc_original" else NU).strength()
    s3 = np.linspace(0.5, 30.0, 50) if model == "ejc_original" else np.linspace(-0.09, 30.0, 50)
    ps = E.principal_states(s3, st)
    ok = ps.code == E.OK
    assert ok.any()
    p, q = st.pq
    P, R, S = E.closure_constants(p, q)
    K = st.sigma_ci ** 2 * ps.Q[ok] ** (2 * st.a)
    s1, s2, s3o = ps.s1[ok], ps.s2[ok], s3[ok]
    assert np.allclose(P * s1 ** 2 + R * s1 * s3o + S * s3o ** 2, K, rtol=1e-12, atol=1e-10)
    lhs = 2 * s1 - s2 - s3o
    assert np.all(np.abs(lhs - np.sqrt(4 * K - 3 * (s2 - s3o) ** 2)) < 1e-8)
    assert np.all((s1 >= s2) & (s2 >= s3o))


def test_moscarda_gradient_matches_brief_and_finite_difference():
    s3 = np.array([-0.05, 0.08, 0.5, 2.0])
    ps = E.principal_states(s3, ST)
    C, B = 1 - NU + NU ** 2, 1 + 2 * NU - 2 * NU ** 2
    sci, mb, a = ST.sigma_ci, ST.mb, ST.a
    Q = mb * s3 / sci + ST.s
    M = (2 * a * mb * sci * Q ** (2 * a - 1) + B * ps.s1 - 2 * C * s3) / (2 * C * ps.s1 - B * s3)
    assert np.allclose(ps.M, M, rtol=1e-13)
    h = 1e-7
    fd = (E.principal_states(s3 + h, ST).s1 - E.principal_states(s3 - h, ST).s1) / (2 * h)
    assert np.allclose(ps.M, fd, rtol=1e-6)


def test_ghb_closed_form():
    s3 = np.array([-0.05, 0.1, 1.0])
    ps = E.principal_states(s3, GHB)
    Q = GHB.mb * s3 / GHB.sigma_ci + GHB.s
    assert np.allclose(ps.s1, s3 + GHB.sigma_ci * Q ** 0.5, rtol=1e-15)
    assert np.allclose(ps.M, 1 + 0.5 * GHB.mb * Q ** -0.5, rtol=1e-15)


# ---- §1.3 tangent parameters ---------------------------------------------------------------------------------------
def test_tangent_identities():
    ps = E.principal_states(np.array([-0.05, 0.08, 0.5]), ST)
    M = ps.M
    assert np.allclose(np.sin(ps.phi), (M - 1) / (M + 1))
    assert np.allclose(ps.sn, ps.s3 + (ps.s1 - ps.s3) / (M + 1))
    assert np.allclose(ps.tau, (ps.s1 - ps.s3) * np.sqrt(M) / (M + 1))
    assert np.allclose(ps.c, ps.tau - ps.sn * np.tan(ps.phi))
    # the tangent point lies on the Mohr circle of (sigma1, sigma3)
    centre, radius = (ps.s1 + ps.s3) / 2, (ps.s1 - ps.s3) / 2
    assert np.allclose((ps.sn - centre) ** 2 + ps.tau ** 2, radius ** 2, rtol=1e-12)


# ---- §3 strength fixtures ----------------------------------------------------------------------------------------
@pytest.mark.parametrize("sn,phi,c,tau", [(50, 67.45, 362.5, 482.9), (400, 59.90, 487.2, 1177.2),
                                          (1600, 50.21, 935.3, 2856.0)])
def test_moscarda_fixtures(sn, phi, c, tau):
    r = E.solve_strength([sn], ST)
    assert r.accepted[0], r.reason(0)
    assert rounds_to(math.degrees(r.phi[0]), phi, 2)
    assert rounds_to(r.c[0], c, 1)
    assert rounds_to(r.tau[0], tau, 1)
    assert abs(r.sn[0] - sn) < 1e-9
    assert r.n_brackets[0] == 1


def test_ghb_fixture():
    r = E.solve_strength([400.0], GHB)
    assert rounds_to(math.degrees(r.phi[0]), 58.45, 2)
    assert rounds_to(r.c[0], 460.0, 0)
    assert rounds_to(r.tau[0], 1111.0, 0)


def test_extra_fixture_original_closure():
    """Extra fixture, original closure (nu = 0.05): sigma_ci 120, m_b 3.5, s 0.002, a 0.505, sigma_n 5.7403 MPa.
    The Moscarda closure has no admissible state there."""
    r = E.solve_strength([5740.3], E.Strength(120.0, 3.5, 0.002, 0.505, "ejc_original", 0.05))
    assert r.accepted[0]
    assert rounds_to(math.degrees(r.phi[0]), 55.0720, 4)
    assert rounds_to(r.c[0], 2079.5, 1)
    m = E.solve_strength([5740.3], E.Strength(120.0, 3.5, 0.002, 0.505, "ejc_moscarda", 0.05))
    assert m.code[0] == E.NO_STATE


def test_ejc_route_at_sigma2_eq_sigma3_reproduces_ghb_strength():
    sn = np.array([-90.0, -20.0, 0.0, 50.0, 400.0, 1600.0, 5000.0])
    e = E.solve_strength(sn, MAIN.with_(model="ejc_sigma2_eq_sigma3").strength())
    g = E.solve_strength(sn, GHB)
    assert e.accepted.all() and g.accepted.all()
    assert np.allclose(e.phi, g.phi, rtol=1e-12, atol=0)
    assert np.allclose(e.c, g.c, rtol=1e-11, atol=0)


# ---- acceptance: rejections are reported, never replaced ----------------------------------------------------------
def test_below_tensile_limit_is_rejected():
    r = E.solve_strength([-98.1, -200.0], ST)
    assert not r.accepted.any()
    assert np.isnan(r.phi).all()


def test_moscarda_admissible_range_starts_above_tensile_limit():
    mn = E.admissible_sn_min(ST)
    assert ST.sigma3_tensile * 1000 < mn < -97.9
    assert E.solve_strength([mn], ST).accepted[0]
    assert E.solve_strength([mn - 1e-6], ST).code[0] == E.NO_STATE


def test_acceptance_codes():
    ps = E.principal_states(np.array([ST.sigma3_tensile - 0.01]), ST)
    assert ps.code[0] == E.DOMAIN
    # original closure at nu = 0.271 has no admissible state at low stress
    st = MAIN.with_(model="ejc_original").strength()
    assert E.solve_strength([1000.0], st).code[0] == E.NO_STATE


def test_original_closure_minimum():
    mb, s, _ = MAIN.hb
    r = E.original_closure_minimum(77.7, mb, s, 0.5, 0.271)
    assert rounds_to(r.sigma3_min, 34.376, 3)
    assert rounds_to(r.sn_min, 52.97, 2)


@pytest.mark.parametrize("nu,target,d", [(0.05, 0.92, 2), (0.10, 3.43, 2), (0.20, 19.23, 2), (0.271, 52.97, 2),
                                         (0.35, 177.8, 1), (0.45, 2188.0, 0)])
def test_table_5_2(nu, target, d):
    mb, s, _ = MAIN.hb
    assert rounds_to(E.original_closure_minimum(77.7, mb, s, 0.5, nu).sn_min, target, d)
