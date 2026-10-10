"""The numpy accelerations must agree with the scalar VBA port elementwise."""
import math

import numpy as np
import pytest

from slope3d import Config, strength as st
from slope3d._vector import check_state_vec, fastphi_closure3_vec, gstate_vec, phi_b_vec

C = Config()
MB, S = C.mb, C.s


@pytest.mark.parametrize("a", [0.5, 0.509923])
def test_fastphi_vector_matches_scalar(a):
    rng = np.random.default_rng(1)
    sn = np.concatenate([[-97.9, -92.0, -80.0, -40.0, 0.0, 20.0, 100.0, 400.0, 800.0, 5000.0],
                         rng.uniform(-98.0, 3000.0, 400)])
    guess = np.where(rng.random(sn.size) < 0.5, -1.0, rng.uniform(5, 85, sn.size))
    vec = fastphi_closure3_vec(sn, 77.7, MB, S, a, guess, 1e-12)
    for i, x in enumerate(sn):
        ref = st.fastphi(float(x), 77.7, MB, S, a, 0, 1, float(guess[i]), 1e-12)
        if isinstance(ref, st.ErrorCode):
            assert np.isnan(vec[i]), (x, ref, vec[i])
        else:
            assert vec[i] == pytest.approx(ref, abs=1e-9), (x, guess[i])


def test_gstate_and_phib_vector_match_scalar():
    rng = np.random.default_rng(2)
    sn = rng.uniform(-0.098, 3.0, 200)
    phi = rng.uniform(0.05, 1.55, 200)
    for p, q in ((0, 1), (-1, 1 / 0.271), (0.271, 0.271)):
        g, dg, m, ok = gstate_vec(phi, sn, 77.7, MB, S, 0.5, p, q)
        for i in range(sn.size):
            ref = st.gstate(phi[i], sn[i], 77.7, MB, S, 0.5, p, q)
            assert ok[i] == (ref is not None)
            if ref is not None:
                assert (g[i], dg[i], m[i]) == pytest.approx(ref, rel=1e-12)
    pb, okb = phi_b_vec(sn, 77.7, MB, S, 0.5)
    for i in range(sn.size):
        ref = st.phi_b(sn[i], 77.7, MB, S, 0.5)
        assert okb[i] == (ref is not None)
        if ref is not None:
            assert pb[i] == pytest.approx(ref, abs=1e-13)


@pytest.mark.parametrize("p,q", [(0, 1), (-1, 1 / 0.271), (0.271, 0.271)])
def test_check_state_vector_matches_scalar_labels_and_residuals(p, q):
    rng = np.random.default_rng(3)
    n = 500
    sn = rng.uniform(-100, 2000, n)
    c = rng.uniform(0, 3000, n)
    phi = np.concatenate([rng.uniform(-5, 95, n - 3), [0.0, 89.9999, 90.0]])
    vec = check_state_vec(sn, c, phi, 77.7, MB, S, 0.5, p, q, 1e-10, 1e-8)
    for i in range(n):
        chk = st.ejc_check_state(sn[i], c[i], phi[i], 77.7, MB, S, 0.5, p, q)
        assert vec.label[i] == st.accept_text(chk, 1e-10, 1e-8)
        assert vec.domain[i] == chk.domain and vec.branch[i] == chk.branch and vec.ordering[i] == chk.ordering
        if chk.domain == 1:
            assert vec.angular[i] == pytest.approx(chk.angular, rel=1e-9, abs=1e-12)
            assert vec.ejc[i] == pytest.approx(chk.ejc, rel=1e-9, abs=1e-12)


def test_check_state_vector_on_a_real_admissible_state():
    ph = st.ejc_phi(400, 77.7, MB, S, 0.5, 0, 1)
    chk = st.ejc_check(400, ph, 77.7, MB, S, 0.5, 0, 1)
    c = 1000 * ((chk.s1 - chk.s3) / 2 * math.cos(math.radians(ph)) - 0.4 * math.tan(math.radians(ph)))
    v = check_state_vec(np.array([400.0]), np.array([c]), np.array([ph]), 77.7, MB, S, 0.5, 0, 1, 1e-10, 1e-8)
    assert v.label[0] == "ADMISSIBLE"
