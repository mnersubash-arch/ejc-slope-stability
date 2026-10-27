"""Hungr 3D columns and Bishop 2D slices (brief §1.4-1.6) and the §3 FoS targets."""
import math

import numpy as np
import pytest

import engine as E
from cases import MAIN, WET, r2, r3


def rounds_to(value, target, decimals):
    return value is not None and abs(value - target) <= 0.5 * 10.0 ** (-decimals) * (1 + 1e-12)


# ---- mesh ----------------------------------------------------------------------------------------------------------
def test_active_counts():
    assert E.build_mesh(MAIN, 3, 20).count == 76
    assert E.build_mesh(MAIN, 2, 40).count == 30


def test_column_quantities():
    m = E.build_mesh(MAIN, 3, 20)
    x0 = MAIN.xc - math.sqrt(MAIN.R ** 2 - MAIN.zc ** 2)
    assert m.x0 == pytest.approx(x0) and m.x1 == pytest.approx(MAIN.xc - MAIN.R)
    assert np.allclose(m.x, x0 - (m.i - 0.5) * m.dx) and np.allclose(m.y, (m.j - 0.5) * MAIN.R / 20)
    assert np.all(m.x <= 0) and np.all(m.h > 0)
    assert np.allclose(m.area * m.g, m.ah)
    assert np.allclose(m.w, 27.0 * m.h * m.ah)
    assert np.allclose(m.kap, (MAIN.xc - m.x) / np.sqrt(MAIN.R ** 2 - m.y ** 2))
    assert np.all(m.u == 0)
    w = E.build_mesh(WET, 3, 20)
    zw = np.minimum(360.0, -w.x * math.tan(math.radians(50.0)))
    assert np.allclose(w.u, 9.81 * np.maximum(0.0, zw - w.zb))


def test_2d_slices():
    m = E.build_mesh(MAIN, 2, 40)
    alpha = np.arcsin(m.kap)
    assert np.allclose(m.g, np.cos(alpha)) and np.allclose(m.rho, MAIN.R) and np.allclose(m.ah, m.dx)


# ---- main targets ----------------------------------------------------------------------------------------------
def test_dry_a05():
    a3, a2 = r3(MAIN), r2(MAIN)
    assert a3.accepted and a2.accepted, (a3.reason, a2.reason)
    assert abs(a3.fos - 2.706024342) <= 1e-6
    assert abs(a2.fos - 2.447910895) <= 1e-6


def test_wet_a05():
    a3, a2 = r3(WET), r2(WET)
    assert abs(a3.fos - 1.272331405) <= 1e-6
    assert abs(a2.fos - 1.098069496) <= 1e-6


def test_native_v15_path_iteration_counts():
    """Iteration counts 34 / 63 are checked on the native v15 path only (the default path shows its count)."""
    d, w = r3(MAIN.with_(iteration="native_v15")), r3(WET.with_(iteration="native_v15"))
    assert (d.iterations, w.iterations) == (34, 63)
    assert abs(d.fos - 2.706024342) <= 1e-6 and abs(w.fos - 1.272331405) <= 1e-6
    assert d.projections > 0 and w.projections > 0          # the path passes below the admissible edge


def test_native_v15_path_refuses_c0_t0_and_fm():
    m = E.build_mesh(MAIN, 3, 20)
    for kw in (dict(tension="C0"), dict(tension="T0"), dict(form="FM")):
        with pytest.raises(ValueError):
            E.solve(m, MAIN.strength(), MAIN.xc, MAIN.gamma, iteration="native_v15", **kw)


def test_projection_rule():
    st = MAIN.strength()
    tensile = st.sigma3_tensile * 1000
    edge = E.admissible_sn_min(st)
    assert E.projection_limit(st) == max(tensile, edge) == edge         # Moscarda edge lies above the tensile limit
    assert E.projection_limit(MAIN.with_(model="ghb").strength()) == tensile
    r = r3(MAIN.with_(iteration="native_v15"))
    assert r.projection_value == pytest.approx(edge + max(1e-9, 0.001 * abs(edge)), rel=1e-15)
    assert E.solve_strength([r.projection_value], st).accepted[0]
    names = [c.name for c in r.checks]
    assert "Projected bases at convergence" in names and all(c.passed for c in r.checks)


@pytest.mark.parametrize("wet,a,t3,t2,g3,g2", [(False, 0.5, None, None, 2.554, 2.310),
                                               (True, 0.5, None, None, 1.200, 1.036),
                                               (False, 0.509923, 2.624, 2.376, 2.476, 2.242),
                                               (True, 0.509923, 1.228, 1.061, 1.158, 1.000)])
def test_rounded_targets(wet, a, t3, t2, g3, g2):
    inp = MAIN.with_(wet=wet, a_override=a)
    if t3 is not None:
        assert rounds_to(r3(inp).fos, t3, 3) and rounds_to(r2(inp).fos, t2, 3)
    g = inp.with_(model="ghb")
    assert rounds_to(r3(g).fos, g3, 3) and rounds_to(r2(g).fos, g2, 3)


@pytest.mark.parametrize("wet", [False, True])
@pytest.mark.parametrize("dim", [3, 2])
def test_sigma2_eq_sigma3_reproduces_ghb_fos(wet, dim):
    fn = r3 if dim == 3 else r2
    e = fn(MAIN.with_(wet=wet, model="ejc_sigma2_eq_sigma3")).fos
    g = fn(MAIN.with_(wet=wet, model="ghb")).fos
    assert abs(e - g) / g <= 1e-11


def test_final_checks_all_pass_and_are_independent():
    a = r3(MAIN)
    assert all(c.passed for c in a.checks), [c for c in a.checks if not c.passed]
    assert np.max(a.local_res) <= 1e-9 and a.global_res <= 1e-8
    # vertical force balance of each column at the reported F
    m, F = a.mesh, a.fos
    S = a.tau_used * m.area / F
    assert np.allclose(m.w, a.sn * m.area * m.g + m.u * m.area * m.g + S * m.kap, rtol=1e-9)


def test_2d_matches_bishop_expression():
    a = r2(MAIN)
    m, F = a.mesh, a.fos
    alpha = np.arcsin(m.kap)
    t = np.tan(a.strength.phi)
    bishop = np.sum((a.strength.c * m.dx + (m.w - m.u * m.dx) * t) / (np.cos(alpha) + np.sin(alpha) * t / F)) / \
        np.sum(m.w * np.sin(alpha))
    assert abs(bishop - F) / F <= 1e-8


# ---- §1.4 moment form -------------------------------------------------------------------------------------------
def test_moment_form():
    d, w = r3(MAIN.with_(form="FM")), r3(WET.with_(form="FM"))
    assert rounds_to(d.fos, 2.7022, 4) and rounds_to(w.fos, 1.2695, 4)
    # re-solved with its own iteration, not evaluated at the F_H state
    assert d.fos != pytest.approx(r3(MAIN).fos, rel=1e-4)


# ---- §1.6 tension sensitivity ----------------------------------------------------------------------------------
@pytest.mark.parametrize("wet,mode,t3,t2", [(False, "C0", 2.637, 2.430), (False, "T0", 2.597, 2.429),
                                            (True, "C0", 0.964, 0.978), (True, "T0", None, None)])
def test_tension_sensitivity(wet, mode, t3, t2):
    inp = MAIN.with_(wet=wet, tension=mode)
    a3, a2 = r3(inp), r2(inp)
    if t3 is None:
        assert a3.fos is None and a2.fos is None
    else:
        assert rounds_to(a3.fos, t3, 3) and rounds_to(a2.fos, t2, 3)
        assert "not a calibrated contact law" in a3.reason


# ---- original closure: never a FoS ---------------------------------------------------------------------------------
def test_original_closure_never_gives_fos():
    a = r3(MAIN.with_(model="ejc_original"))
    assert a.status == "NOT APPLICABLE" and a.fos is None


def test_gravity_cap_and_threshold_count():
    md, mw = E.build_mesh(MAIN, 3, 20), E.build_mesh(WET, 3, 20)
    assert rounds_to(E.gravity_cap(md, 27.0).max() / 1000, 1.61, 2)
    assert rounds_to(E.gravity_cap(mw, 27.0).max() / 1000, 1.02, 2)
    mb, s, _ = MAIN.hb
    thr = E.original_closure_minimum(77.7, mb, s, 0.5, 0.05).sn_min * 1000
    assert int(np.sum(E.gravity_cap(md, 27.0) < thr)) == 43


# ---- previous studies ----------------------------------------------------------------------------------------------
CT = MAIN.with_(xc=207.28, zc=586.53, R=622.08)


def test_carranza_torres_default_meshes():
    assert rounds_to(r2(CT.with_(model="ghb")).fos, 2.0039, 4)
    assert rounds_to(r2(CT.with_(model="ghb")).fos, 2.003873, 6)
    assert rounds_to(r3(CT.with_(model="ghb")).fos, 2.142, 3)
    assert rounds_to(r3(CT.with_(model="ghb")).fos, 2.141937, 6)
    assert rounds_to(r2(CT).fos, 2.122598, 6)
    assert rounds_to(r3(CT).fos, 2.269068, 6)


@pytest.mark.slow
def test_carranza_torres_fine_meshes():
    assert rounds_to(r2(CT.with_(n2=5120)).fos, 2.121420, 6)
    assert rounds_to(r3(CT.with_(model="ghb", n3=320)).fos, 2.140776, 6)
    assert rounds_to(r3(CT.with_(n3=320)).fos, 2.267835, 6)


# ---- rejection is explicit ---------------------------------------------------------------------------------------
def test_rejected_run_withholds_fos():
    """Without the projection (diagnostic), a crest-end slice's trial sigma_n falls below the Moscarda admissible range
    on the fine wet 2D mesh: the run is REJECTED and F withheld, never replaced."""
    a = r2(WET.with_(n2=320, safeguard="none"))
    assert a.status == "REJECTED" and a.fos is None and a.f_trial is not None
    assert "no admissible state" in a.reason


def test_projection_only_changes_the_path():
    a = r2(WET.with_(n2=320))
    assert a.accepted and a.projections > 0 and a.history[-1][5] == 0
    b = r3(MAIN.with_(safeguard="none"))
    assert b.fos == pytest.approx(r3(MAIN).fos, rel=1e-12)


# ---- finest mesh (slow) ----------------------------------------------------------------------------------------
@pytest.mark.slow
def test_finest_dry():
    assert rounds_to(r3(MAIN.with_(n3=240)).fos, 2.7112, 4)
    assert rounds_to(r2(MAIN.with_(n2=5120)).fos, 2.4422, 4)


@pytest.mark.slow
def test_finest_wet():
    a3, a2 = r3(WET.with_(n3=240)), r2(WET.with_(n2=5120))
    assert rounds_to(a3.fos, 1.2776, 4) and rounds_to(a2.fos, 1.0916, 4)
    assert a3.projections > 0 and a2.projections > 0         # projected on the path, none at convergence


@pytest.mark.slow
def test_coarse_mesh_underestimates():
    for inp, pct in ((MAIN, 2.6), (WET, 4.4)):
        fine, coarse = r3(inp.with_(n3=240)).fos, r3(inp.with_(n3=10)).fos
        assert rounds_to((fine - coarse) / fine * 100, pct, 1)


def test_mesh_error_formula():
    assert E.mesh_error([2.0, 2.1, 2.05], 2.05) == pytest.approx(0.1 / 2.05 * 100)
