"""Suite A: the 22 strength fixtures of ``EJC_RunValidation`` (EJC_Eq421.bas lines 342-447).

Assertions follow the VBA, not the Python mirror, wherever the two differ (see docs/SPEC.md section 0).
"""
import math

import pytest

from slope3d.strength import (ErrorCode, accept_text, ejc_accept, ejc_accept_state, ejc_check,
                              ejc_check_state, ejc_phi)

MB = 10 * math.exp((42 - 100) / 28)
S = math.exp((42 - 100) / 9)
V = (120, 3.5, 0.002, 0.505)          # sigma_ci, mb, s, a of the T2/T7 counterexample
T7SN = 13265.01174607056
ADM = "ADMISSIBLE"


def chile(sn, a=0.5, p=0, q=1, **kw):
    return ejc_phi(sn, 77.7, MB, S, a, p, q, **kw)


def label(sn, ph, *args):
    return accept_text(ejc_check(sn, ph, *args), 1e-10, 1e-8)


def test_T1_closure3_chile():
    ph = chile(400)
    assert isinstance(ph, float) and abs(ph - 58.448348604037) < 1e-9
    assert label(400, ph, 77.7, MB, S, 0.5, 0, 1) == ADM


def test_T2_closure1_verification_point_first_root():
    ph = ejc_phi(5740.3, *V, -1, 20)
    assert isinstance(ph, float) and abs(ph - 53.3178012) < 1e-6
    chk = ejc_check(5740.3, ph, *V, -1, 20)
    assert accept_text(chk, 1e-10, 1e-8) == "REJECTED: EJC"
    assert chk.ejc == pytest.approx(-8.882, abs=1e-3)      # VBA Validation!E5 records -8.882E+00 MPa


def test_T3_near_tensile_limit_admissible_or_controlled_error():
    sn = -98.00411947453206
    ph = chile(sn)
    if not isinstance(ph, ErrorCode):
        assert label(sn, ph, 77.7, MB, S, 0.5, 0, 1) == ADM
        assert ph == pytest.approx(89.6614, abs=1e-4)


def test_T4_valid_hb_point_near_89_95():
    sn = -98.0041250683277
    ph = chile(sn)
    assert isinstance(ph, float) and abs(ph - 89.95) < 0.01
    assert label(sn, ph, 77.7, MB, S, 0.5, 0, 1) == ADM


def test_T5_thesis_closure_chile_has_no_root():
    ph = chile(400, p=-1, q=1 / 0.271)
    assert ph == ErrorCode.NA and ph == "#N/A"


def test_T6_invalid_sigma_ci():
    assert ejc_phi(400, 0, MB, S, 0.5, 0, 1) == ErrorCode.VALUE


def test_T7_isolated_counterexample_has_admissible_root():
    hit = False
    for i in range(1, 5):
        ph = ejc_phi(T7SN, *V, -1, 20, 1e-10, i)
        if isinstance(ph, float) and abs(ph - 39.89156042767222) < 1e-6:
            hit = label(T7SN, ph, *V, -1, 20) == ADM
    assert hit


@pytest.fixture(scope="module")
def t1_state():
    ph = chile(400)
    chk = ejc_check(400, ph, 77.7, MB, S, 0.5, 0, 1)
    c = 1000 * ((chk.s1 - chk.s3) / 2 * math.cos(ph * 3.14159265358979 / 180) - 0.4 * math.tan(ph * 3.14159265358979 / 180))
    return ph, c


def test_T8_force_used_state_of_T1(t1_state):
    ph, c = t1_state
    assert ejc_accept_state(400, c, ph, 77.7, MB, S, 0.5, 0, 1) == ADM


def test_T9_perturbed_force_used_state(t1_state):
    ph, c = t1_state
    assert ejc_accept_state(400, c + 5, ph, 77.7, MB, S, 0.5, 0, 1) != ADM


def test_T10_root_index_zero():
    assert ejc_phi(T7SN, *V, -1, 20, 1e-10, 0) == ErrorCode.VALUE


def test_T11_root_index_negative():
    assert ejc_phi(T7SN, *V, -1, 20, 1e-10, -1) == ErrorCode.VALUE


def test_T12_root2_is_rejected():
    ph = ejc_phi(T7SN, *V, -1, 20, 1e-10, 2)
    assert isinstance(ph, float) and abs(ph - 60.4228863398) < 1e-6
    assert label(T7SN, ph, *V, -1, 20) != ADM


def test_T13_T14_tolerance_zero_or_negative():
    assert ejc_phi(T7SN, *V, -1, 20, 0, 1) == ErrorCode.VALUE
    assert ejc_phi(T7SN, *V, -1, 20, -1e-10, 1) == ErrorCode.VALUE


def test_T15_negative_friction_branch_state():
    lbl = ejc_accept_state(29819.0191175968, 25010.6868294027, -0.012480264998, *V, -1, 20)
    assert lbl in ("REJECTED: DOMAIN", "REJECTED: BRANCH")


def test_T16_candidate_at_minus_45():
    assert ejc_accept(T7SN, -45, *V, -1, 20) == "REJECTED: DOMAIN"


def test_T17_candidate_at_90():
    assert ejc_accept(T7SN, 90, *V, -1, 20) == "REJECTED: DOMAIN"


@pytest.mark.parametrize("ta,te", [(1e-10, 0), (1e-10, -1e-8), (0, 1e-8), (-1e-10, 1e-8)])
def test_T18_to_T21_invalid_acceptance_tolerances(t1_state, ta, te):
    ph, c = t1_state
    assert ejc_accept_state(400, c, ph, 77.7, MB, S, 0.5, 0, 1, ta, te) == "REJECTED: INPUT"


def test_T22_wet_apex_scan_regression():
    sn = -49.768342
    ph = ejc_phi(sn, 77.7, MB, S, 0.509923, 0, 1)
    assert isinstance(ph, float) and abs(ph - 71.27787820469955) < 1e-7
    assert label(sn, ph, 77.7, MB, S, 0.509923, 0, 1) == ADM


def test_non_array_check_label_follows_vba():
    assert accept_text(None, 1e-10, 1e-8) == "REJECTED: INPUT/DOMAIN"   # VBA text; the Python mirror says "DOMAIN"
