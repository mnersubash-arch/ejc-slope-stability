"""Built-in regression checks behind ``slope3d validate`` (works without pytest, e.g. in the packaged .exe).

The same expectations are asserted in detail by the pytest suites in ``tests/``.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from . import strength as st
from .engine import run_case, solve_model, validate_final
from .mesh_study import compare_with_records, exact_ejc_threshold, run_mesh_study
from .models import Config, Status
from .presets import PRESETS

REL = 1e-9
MB = 10 * math.exp((42 - 100) / 28)
S = math.exp((42 - 100) / 9)
V = (120, 3.5, 0.002, 0.505)
T7SN = 13265.01174607056
ADM = "ADMISSIBLE"


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str = ""
    seconds: float = 0.0


def _lab(sn, ph, *args):
    return st.accept_text(st.ejc_check(sn, ph, *args), 1e-10, 1e-8)


def _approx(x: Optional[float], y: float, rel: float = REL) -> bool:
    return x is not None and abs(x - y) <= rel * abs(y)


# ---- A. strength fixtures T1-T22 (EJC_RunValidation) ----------------------------------------------
def _strength_checks() -> list[tuple[str, Callable[[], str]]]:
    out: list[tuple[str, Callable[[], str]]] = []

    def add(name):
        def deco(fn):
            out.append((name, fn))
            return fn
        return deco

    ph1 = st.ejc_phi(400, 77.7, MB, S, 0.5, 0, 1)
    chk1 = st.ejc_check(400, ph1, 77.7, MB, S, 0.5, 0, 1) if isinstance(ph1, float) else None
    c1 = 0.0
    if chk1 is not None and chk1.s1 is not None:
        r = ph1 * st.PI_ / 180
        c1 = 1000 * ((chk1.s1 - chk1.s3) / 2 * math.cos(r) - 0.4 * math.tan(r))

    @add("A T1 closure 3, Chile, 400 kPa -> 58.448348604037, ADMISSIBLE")
    def _():
        assert isinstance(ph1, float) and abs(ph1 - 58.448348604037) < 1e-9 and _lab(400, ph1, 77.7, MB, S, 0.5, 0, 1) == ADM
        return f"{ph1:.12f}"

    @add("A T2 closure 1 verification point -> 53.3178012, REJECTED: EJC")
    def _():
        ph = st.ejc_phi(5740.3, *V, -1, 20)
        chk = st.ejc_check(5740.3, ph, *V, -1, 20)
        assert abs(ph - 53.3178012) < 1e-6 and st.accept_text(chk, 1e-10, 1e-8) == "REJECTED: EJC"
        return f"{ph:.7f}, EJC residual {chk.ejc:.4f} MPa"

    @add("A T3 near tensile limit -> admissible or controlled error")
    def _():
        sn = -98.00411947453206
        ph = st.ejc_phi(sn, 77.7, MB, S, 0.5, 0, 1)
        assert isinstance(ph, st.ErrorCode) or _lab(sn, ph, 77.7, MB, S, 0.5, 0, 1) == ADM
        return str(ph)

    @add("A T4 valid HB point near 89.95 deg")
    def _():
        sn = -98.0041250683277
        ph = st.ejc_phi(sn, 77.7, MB, S, 0.5, 0, 1)
        assert isinstance(ph, float) and abs(ph - 89.95) < 0.01 and _lab(sn, ph, 77.7, MB, S, 0.5, 0, 1) == ADM
        return f"{ph:.6f}"

    @add("A T5 thesis closure, Chile, nu = 0.271 -> #N/A")
    def _():
        ph = st.ejc_phi(400, 77.7, MB, S, 0.5, -1, 1 / 0.271)
        assert ph == st.ErrorCode.NA
        return str(ph.value)

    @add("A T6 sigma_ci = 0 -> #VALUE!")
    def _():
        assert st.ejc_phi(400, 0, MB, S, 0.5, 0, 1) == st.ErrorCode.VALUE
        return "#VALUE!"

    @add("A T7 isolated counterexample 13.265 MPa -> admissible root at 39.8915604")
    def _():
        hit = False
        for i in range(1, 5):
            ph = st.ejc_phi(T7SN, *V, -1, 20, 1e-10, i)
            if isinstance(ph, float) and abs(ph - 39.89156042767222) < 1e-6:
                hit = _lab(T7SN, ph, *V, -1, 20) == ADM
        assert hit
        return "found"

    @add("A T8 force-used state of T1 -> ADMISSIBLE")
    def _():
        assert st.ejc_accept_state(400, c1, ph1, 77.7, MB, S, 0.5, 0, 1) == ADM
        return f"c = {c1:.6f} kPa"

    @add("A T9 force-used state, c + 5 kPa -> REJECTED")
    def _():
        lab = st.ejc_accept_state(400, c1 + 5, ph1, 77.7, MB, S, 0.5, 0, 1)
        assert lab != ADM
        return lab

    @add("A T10/T11 root index 0 and -1 -> #VALUE!")
    def _():
        assert st.ejc_phi(T7SN, *V, -1, 20, 1e-10, 0) == st.ErrorCode.VALUE
        assert st.ejc_phi(T7SN, *V, -1, 20, 1e-10, -1) == st.ErrorCode.VALUE
        return "#VALUE!"

    @add("A T12 root 2 at T7 inputs -> 60.4228863398, rejected")
    def _():
        ph = st.ejc_phi(T7SN, *V, -1, 20, 1e-10, 2)
        lab = _lab(T7SN, ph, *V, -1, 20)
        assert abs(ph - 60.4228863398) < 1e-6 and lab != ADM
        return f"{ph:.10f} {lab}"

    @add("A T13/T14 tolerance 0 and negative -> #VALUE!")
    def _():
        assert st.ejc_phi(T7SN, *V, -1, 20, 0, 1) == st.ErrorCode.VALUE
        assert st.ejc_phi(T7SN, *V, -1, 20, -1e-10, 1) == st.ErrorCode.VALUE
        return "#VALUE!"

    @add("A T15 negative-friction branch state -> DOMAIN or BRANCH")
    def _():
        lab = st.ejc_accept_state(29819.0191175968, 25010.6868294027, -0.012480264998, *V, -1, 20)
        assert lab in ("REJECTED: DOMAIN", "REJECTED: BRANCH")
        return lab

    @add("A T16/T17 candidate at -45 and 90 deg -> REJECTED: DOMAIN")
    def _():
        assert st.ejc_accept(T7SN, -45, *V, -1, 20) == "REJECTED: DOMAIN"
        assert st.ejc_accept(T7SN, 90, *V, -1, 20) == "REJECTED: DOMAIN"
        return "REJECTED: DOMAIN"

    @add("A T18-T21 invalid acceptance tolerances -> REJECTED: INPUT")
    def _():
        for ta, te in ((1e-10, 0), (1e-10, -1e-8), (0, 1e-8), (-1e-10, 1e-8)):
            assert st.ejc_accept_state(400, c1, ph1, 77.7, MB, S, 0.5, 0, 1, ta, te) == "REJECTED: INPUT"
        return "REJECTED: INPUT"

    @add("A T22 wet a=0.509923, -49.768342 kPa -> 71.27787820469955, ADMISSIBLE")
    def _():
        ph = st.ejc_phi(-49.768342, 77.7, MB, S, 0.509923, 0, 1)
        assert abs(ph - 71.27787820469955) < 1e-7 and _lab(-49.768342, ph, 77.7, MB, S, 0.509923, 0, 1) == ADM
        return f"{ph:.12f}"

    return out


# ---- B. engine regressions ---------------------------------------------------------------------------
def _engine_checks() -> list[tuple[str, Callable[[], str]]]:
    out: list[tuple[str, Callable[[], str]]] = []

    def bench(name):
        def fn():
            p = PRESETS[name]
            r = run_case(p.config, 20, 40)
            for run, e in ((r.run3d, p.expected_3d), (r.run2d, p.expected_2d)):
                assert run.result.status is Status.ACCEPTED, run.result.reason
                assert _approx(run.result.fos, e), (run.result.fos, e)
            return f"3D {r.run3d.result.fos:.12f}  2D {r.run2d.result.fos:.12f}"
        return fn

    for name in ("Dry, a = 0.5", "Wet, a = 0.5", "Dry, a = 0.509923", "Wet, a = 0.509923",
                 "Published circle, dry, a = 0.5", "Mohr-Coulomb c = 100 kPa, phi = 0"):
        out.append((f"B {name}: F 3D 20x20 / 2D 40 slices", bench(name)))

    def thesis():
        r = run_case(Config(), 20, 40)
        for run, n in ((r.run3d, 76), (r.run2d, 30)):
            assert run.result.status is Status.NOT_ACCEPTED and run.result.fos is None
            assert run.result.active == n
        return "NOT ACCEPTED, FoS withheld, 76 columns / 30 slices"
    out.append(("B Thesis closure (1), Chile defaults: both withheld", thesis))

    def moscarda():
        cfg = Config(closure=2, a_override=0.509923, water_angle=50.0, water_height=360.0)
        r = run_case(cfg, 20, 40)
        assert all(x.result.status is Status.NOT_ACCEPTED and x.result.fos is None for x in (r.run3d, r.run2d))
        return "NOT ACCEPTED, FoS withheld"
    out.append(("B Moscarda closure (2), wet, a=0.509923: NOT ACCEPTED", moscarda))

    def checks_off():
        r = run_case(Config(closure=3, checks=0), 20, 40)
        for run, f in ((r.run3d, 2.5536684200376), (r.run2d, 2.31025502757599)):
            assert run.result.status is Status.DIAGNOSTIC_ONLY and run.result.fos is None
            assert _approx(run.result.internal_f, f)
        return "DIAGNOSTIC ONLY, FoS withheld"
    def diagnostic():
        cfg = Config(closure=2, a_override=0.509923, water_angle=50.0, water_height=360.0, checks=0)
        r = run_case(cfg, 20, 40)
        for run, f in ((r.run3d, 1.04327747791523), (r.run2d, 0.90176054371878)):
            assert run.result.status is Status.DIAGNOSTIC_ONLY and run.result.fos is None
            assert _approx(run.result.internal_f, f)
        return "DIAGNOSTIC ONLY, internal F 1.0432774779 / 0.9017605437 (not shown to users)"
    out.append(("B Moscarda wet, checks = 0 (diagnostic mode): internal F matches the worksheet engine", diagnostic))
    out.append(("B checks = 0 on a valid run: DIAGNOSTIC ONLY, FoS withheld", checks_off))

    def perturb():
        cfg = PRESETS["Wet, a = 0.509923"].config
        run = solve_model(cfg, 3, 20)
        res = run.result
        assert res.accepted
        run.elements.final_n[0] += 1.0
        assert not validate_final(cfg, run.elements, res.internal_f, res, rebuild=False)
        return "rejected after +1 kPa on one column"
    out.append(("B +1 kPa on one final normal stress: validation rejects", perturb))

    def refuse(over):
        def fn():
            r = run_case(Config(closure=3, **over), 20, 40)
            for run in (r.run3d, r.run2d):
                x = run.result
                assert x.fos is None and x.active == 0 and x.reason.startswith("UNSUPPORTED GEOMETRY:")
            return r.run3d.result.reason
        return fn
    out.append(("B original MC geometry refused", refuse(dict(xc=-1.45, zc=3.24, radius=4.46, height=2.1,
                slope=1 / 3.5, gamma=19.0, strength=0, c_mc=7.2, phi_mc=30.0))))
    out.append(("B cropped geometry (zc = 300) refused", refuse(dict(xc=207.28, zc=300.0, radius=622.08))))
    return out


# ---- C. slow checks -----------------------------------------------------------------------------------
def _slow_checks(ref_dir: Optional[Path]) -> list[tuple[str, Callable[[], str]]]:
    out: list[tuple[str, Callable[[], str]]] = []

    def thresholds():
        a = exact_ejc_threshold(0.271, 0.5)[0]
        b = exact_ejc_threshold(0.271, 0.509923)[0]
        assert abs(a - 52.966376) < 5e-7 and abs(b - 52.283035) < 5e-7, (a, b)
        return f"a=0.5: {a:.6f} MPa; a=0.509923: {b:.6f} MPa"
    out.append(("C exact-EJC admissibility thresholds at nu = 0.271", thresholds))

    def mesh():
        rec = (ref_dir / "mesh_study_records.csv") if ref_dir else None
        if rec is None or not rec.exists():
            raise FileNotFoundError("reference/mesh_study/mesh_study_records.csv not found (use --reference DIR)")
        problems = compare_with_records(run_mesh_study(), rec)
        assert not problems, "; ".join(problems[:3])
        return "95 records reproduced to 1e-9 (3D to 320x320, 2D to 2560 slices)"
    out.append(("C mesh-study records (95)", mesh))
    return out


def run_all(slow: bool = False, reference_dir: Optional[str | Path] = None,
            echo: Optional[Callable[[str], None]] = None) -> list[CheckResult]:
    ref = Path(reference_dir) if reference_dir else None
    checks = _strength_checks() + _engine_checks() + (_slow_checks(ref) if slow else [])
    results: list[CheckResult] = []
    for name, fn in checks:
        t0 = time.time()
        try:
            detail, ok = fn(), True
        except Exception as exc:  # noqa: BLE001 - a failing check is reported, not raised
            detail, ok = f"{type(exc).__name__}: {exc}", False
        res = CheckResult(name, ok, detail, time.time() - t0)
        results.append(res)
        if echo:
            echo(f"{'PASS' if ok else 'FAIL'}  {name}  [{detail}]")
    return results
