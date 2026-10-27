"""Verification targets of the brief (section 3) and a runner that reports target, computed value, difference and
PASS / FAIL / INFO. No Streamlit imports; used by the tests, the Verification page and ``python verification.py``.

Tolerances (brief section 3): FoS given to 9 decimals must match to 1e-6; values given to 3-4 decimals must match to
rounding (|computed - target| <= half a unit in the last decimal shown); fixtures to the decimals shown.
INFO rows are reference or "about" values that are not acceptance targets.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

import engine as E

MAIN = E.Inputs()                                     # main sphere, Moscarda closure, EJC, dry, a = 0.5, 20x20 / 40
CARRANZA = MAIN.with_(xc=207.28, zc=586.53, R=622.08)


@dataclass
class Row:
    group: str
    item: str
    target: str
    computed: str
    difference: str
    tolerance: str
    status: str                 # PASS / FAIL / INFO
    note: str = ""


class Runner:
    """Runs and caches the analyses needed by the targets under one iteration path / safeguard setting."""

    def __init__(self, iteration: str = "spec", safeguard: str = "none"):
        self.iteration = iteration
        self.safeguard = safeguard
        self._cache: dict = {}

    def inputs(self, base: E.Inputs = MAIN, **kw) -> E.Inputs:
        kw.setdefault("iteration", self.iteration)
        kw.setdefault("safeguard", self.safeguard)
        return base.with_(**kw)

    def r3(self, base: E.Inputs = MAIN, **kw) -> E.LemResult:
        inp = self.inputs(base, **kw)
        key = ("3", inp)
        if key not in self._cache:
            self._cache[key] = E.run_3d(inp)
        return self._cache[key]

    def r2(self, base: E.Inputs = MAIN, **kw) -> E.LemResult:
        inp = self.inputs(base, **kw)
        key = ("2", inp.with_(form="FH", n3=0))
        if key not in self._cache:
            self._cache[key] = E.run_2d(inp)
        return self._cache[key]


def _fmt(v, d=None) -> str:
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "none"
    if isinstance(v, str):
        return v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return f"{v:.{d}f}" if d is not None else f"{v:.10g}"


def by_rounding(group, item, computed, target, decimals, note="") -> Row:
    tol = 0.5 * 10.0 ** (-decimals)
    if computed is None or not math.isfinite(computed):
        return Row(group, item, f"{target:.{decimals}f}", "none", "n/a", f"rounding to {decimals} dp", "FAIL", note)
    diff = computed - target
    ok = abs(diff) <= tol * (1 + 1e-12)
    return Row(group, item, f"{target:.{decimals}f}", _fmt(computed, decimals + 3), f"{diff:+.2e}",
               f"rounding to {decimals} dp", "PASS" if ok else "FAIL", note)


def by_abs(group, item, computed, target, tol, decimals=9, note="") -> Row:
    if computed is None or not math.isfinite(computed):
        return Row(group, item, f"{target:.{decimals}f}", "none", "n/a", f"abs <= {tol:g}", "FAIL", note)
    diff = computed - target
    return Row(group, item, f"{target:.{decimals}f}", _fmt(computed, 12), f"{diff:+.2e}", f"abs <= {tol:g}",
               "PASS" if abs(diff) <= tol else "FAIL", note)


def by_rel(group, item, computed, reference, tol, note="") -> Row:
    if computed is None or reference is None:
        return Row(group, item, "relative difference", "none", "n/a", f"rel <= {tol:g}", "FAIL", note)
    rel = abs(computed - reference) / abs(reference)
    return Row(group, item, _fmt(reference, 12), _fmt(computed, 12), f"{rel:.2e} (rel)", f"rel <= {tol:g}",
               "PASS" if rel <= tol else "FAIL", note)


def by_exact(group, item, computed, target, note="") -> Row:
    ok = computed == target
    diff = "0" if ok else (f"{computed - target:+d}" if isinstance(computed, int) and isinstance(target, int) else "differs")
    return Row(group, item, _fmt(target), _fmt(computed), diff, "exact", "PASS" if ok else "FAIL", note)


def no_fos(group, item, r: E.LemResult, note="") -> Row:
    ok = r.fos is None
    return Row(group, item, "no accepted FoS", r.status if ok else _fmt(r.fos, 6), "n/a", "status",
               "PASS" if ok else "FAIL", note or r.reason[:160])


def info(group, item, computed, target="", note="") -> Row:
    return Row(group, item, target, computed, "", "reference only", "INFO", note)


# ---------------------------------------------------------------------------------------------------------------
def verification_rows(runner: Optional[Runner] = None, slow: bool = False,
                      progress: Optional[Callable[[str], None]] = None) -> list[Row]:
    R = runner or Runner()
    say = progress or (lambda s: None)
    rows: list[Row] = []

    # §1.1 rock mass
    g = "1.1 Rock mass (GSI 42, m_i 10, D 0)"
    mb, s, a_gsi = E.hoek_brown_parameters(42.0, 10.0, 0.0)
    rows += [by_rounding(g, "m_b", mb, 1.260056, 6), by_rounding(g, "s", s, 0.0015893, 7),
             by_rounding(g, "a from GSI", a_gsi, 0.509923, 6)]

    # strength fixtures
    say("strength fixtures")
    g = "3 Strength fixtures (Moscarda, a = 0.5, nu = 0.271, Route 2)"
    st = MAIN.strength()
    fx = E.solve_strength([50.0, 400.0, 1600.0], st)
    for k, (sn, phi, c, tau) in enumerate([(50, 67.45, 362.5, 482.9), (400, 59.90, 487.2, 1177.2),
                                           (1600, 50.21, 935.3, 2856.0)]):
        rows += [by_rounding(g, f"sigma_n = {sn} kPa: phi_i (deg)", float(np.degrees(fx.phi[k])), phi, 2),
                 by_rounding(g, f"sigma_n = {sn} kPa: c_i (kPa)", float(fx.c[k]), c, 1),
                 by_rounding(g, f"sigma_n = {sn} kPa: tau_s (kPa)", float(fx.tau[k]), tau, 1)]
    gh = E.solve_strength([400.0], MAIN.with_(model="ghb").strength())
    rows += [by_rounding(g, "GHB, sigma_n = 400 kPa: phi_i (deg)", float(np.degrees(gh.phi[0])), 58.45, 2),
             by_rounding(g, "GHB, sigma_n = 400 kPa: c_i (kPa)", float(gh.c[0]), 460.0, 0),
             by_rounding(g, "GHB, sigma_n = 400 kPa: tau_s (kPa)", float(gh.tau[0]), 1111.0, 0)]
    xo = E.solve_strength([5740.3], E.Strength(120.0, 3.5, 0.002, 0.505, "ejc_original", 0.05))
    xm = E.solve_strength([5740.3], E.Strength(120.0, 3.5, 0.002, 0.505, "ejc_moscarda", 0.05))
    note = "matched with the ORIGINAL closure sigma2 = sigma3/nu - sigma1"
    rows += [by_rounding(g, "Extra fixture (sigma_ci 120, m_b 3.5, s 0.002, a 0.505, nu 0.05, 5.7403 MPa): phi_i",
                         float(np.degrees(xo.phi[0])), 55.0720, 4, note),
             by_rounding(g, "Extra fixture: c_i (kPa)", float(xo.c[0]), 2079.5, 1, note),
             info(g, "Extra fixture with the Moscarda closure", xm.reason(0), "",
                  "the Moscarda closure has no admissible state here")]

    # EJC at sigma2 = sigma3 reproduces GHB
    say("sigma2 = sigma3 equivalence")
    g = "3 EJC route at sigma2 = sigma3 reproduces GHB"
    for wet in (False, True):
        for dim, fn in (("3D", R.r3), ("2D", R.r2)):
            e = fn(model="ejc_sigma2_eq_sigma3", wet=wet)
            h = fn(model="ghb", wet=wet)
            rows.append(by_rel(g, f"{dim} {'wet' if wet else 'dry'}, a = 0.5: F(EJC, sigma2 = sigma3) vs F(GHB)",
                               e.fos, h.fos, 1e-11))

    # main cases
    say("main cases")
    g = "3 Main sphere, Moscarda closure, EJC (20x20, 40 slices)"
    d3 = R.r3()
    d2 = R.r2()
    rows += [by_exact(g, "Active half-columns (3D 20x20)", d3.mesh.count, 76),
             by_exact(g, "Active slices (2D 40)", d2.mesh.count, 30)]
    rows += [by_abs(g, "Dry, a = 0.5: F 3D", d3.fos, 2.706024342, 1e-6),
             by_exact(g, "Dry, a = 0.5: 3D iterations", d3.iterations, 34),
             by_abs(g, "Dry, a = 0.5: F 2D", d2.fos, 2.447910895, 1e-6)]
    w3 = R.r3(wet=True)
    w2 = R.r2(wet=True)
    rows += [by_abs(g, "Wet, a = 0.5: F 3D", w3.fos, 1.272331405, 1e-6),
             by_exact(g, "Wet, a = 0.5: 3D iterations", w3.iterations, 63),
             by_abs(g, "Wet, a = 0.5: F 2D", w2.fos, 1.098069496, 1e-6)]
    for wet, a, t3, t2, g3, g2 in ((False, 0.5, None, None, 2.554, 2.310), (True, 0.5, None, None, 1.200, 1.036),
                                    (False, 0.509923, 2.624, 2.376, 2.476, 2.242),
                                    (True, 0.509923, 1.228, 1.061, 1.158, 1.000)):
        lab = f"{'Wet' if wet else 'Dry'}, a = {a:g}"
        if t3 is not None:
            rows += [by_rounding(g, f"{lab}: F 3D", R.r3(wet=wet, a_override=a).fos, t3, 3),
                     by_rounding(g, f"{lab}: F 2D", R.r2(wet=wet, a_override=a).fos, t2, 3)]
        rows += [by_rounding(g, f"{lab}: GHB F 3D", R.r3(wet=wet, a_override=a, model="ghb").fos, g3, 3),
                 by_rounding(g, f"{lab}: GHB F 2D", R.r2(wet=wet, a_override=a, model="ghb").fos, g2, 3)]

    # moment form
    say("moment form")
    g = "3 Moment form F_M (20x20, own iteration)"
    fm_d = R.r3(form="FM")
    fm_w = R.r3(form="FM", wet=True)
    rows += [by_rounding(g, "Dry: F_M 3D", fm_d.fos, 2.7022, 4), by_rounding(g, "Wet: F_M 3D", fm_w.fos, 1.2695, 4)]

    # tension sensitivity
    say("tension sensitivity")
    g = "3 Tension sensitivity (20x20, a = 0.5; sensitivity, not a calibrated contact law)"
    rows += [by_rounding(g, "Dry: C0 3D", R.r3(tension="C0").fos, 2.637, 3),
             by_rounding(g, "Dry: T0 3D", R.r3(tension="T0").fos, 2.597, 3),
             by_rounding(g, "Dry: C0 2D", R.r2(tension="C0").fos, 2.430, 3),
             by_rounding(g, "Dry: T0 2D", R.r2(tension="T0").fos, 2.429, 3),
             by_rounding(g, "Wet: C0 3D", R.r3(tension="C0", wet=True).fos, 0.964, 3),
             no_fos(g, "Wet: T0 3D", R.r3(tension="T0", wet=True)),
             by_rounding(g, "Wet: C0 2D", R.r2(tension="C0", wet=True).fos, 0.978, 3),
             no_fos(g, "Wet: T0 2D", R.r2(tension="T0", wet=True))]

    # original closure
    say("original closure")
    g = "1.2(c) / Table 5.2 Original closure (admissibility check only, a = 0.5)"
    mb, s, _ = MAIN.hb
    m0 = E.original_closure_minimum(77.7, mb, s, 0.5, 0.271)
    rows += [by_rounding(g, "nu = 0.271: sigma3,min (MPa)", m0.sigma3_min, 34.376, 3),
             by_rounding(g, "nu = 0.271: sigma_n,min (MPa)", m0.sn_min, 52.97, 2)]
    md = E.build_mesh(MAIN, 3, 20)
    mw = E.build_mesh(MAIN.with_(wet=True), 3, 20)
    cap_d = E.gravity_cap(md, MAIN.gamma)
    cap_w = E.gravity_cap(mw, MAIN.gamma)
    rows += [by_rounding(g, "Gravity cap max(gamma h - u), dry (MPa, 3D 20x20)", float(cap_d.max()) / 1000, 1.61, 2),
             by_rounding(g, "Gravity cap max(gamma h - u), wet (MPa, 3D 20x20)", float(cap_w.max()) / 1000, 1.02, 2)]
    for nu, tgt, d in ((0.05, 0.92, 2), (0.10, 3.43, 2), (0.20, 19.23, 2), (0.271, 52.97, 2), (0.35, 177.8, 1),
                       (0.45, 2188.0, 0)):
        mn = m0 if nu == 0.271 else E.original_closure_minimum(77.7, mb, s, 0.5, nu)
        rows.append(by_rounding(g, f"Table 5.2: min sigma_n at nu = {nu:g} (MPa)", mn.sn_min, tgt, d))
    thr = E.original_closure_minimum(77.7, mb, s, 0.5, 0.05).sn_min * 1000.0
    rows.append(by_exact(g, "nu = 0.05: dry columns with gamma h - u below the threshold (of 76)",
                         int(np.sum(cap_d < thr)), 43))
    orig = R.r3(model="ejc_original")
    rows.append(Row(g, "Original closure never produces a FoS", "NOT APPLICABLE", orig.status, "n/a", "status",
                    "PASS" if orig.fos is None and orig.status == "NOT APPLICABLE" else "FAIL"))

    # previous studies
    say("Carranza-Torres & Hormazabal (2020)")
    g = "3 Carranza-Torres & Hormazabal (2020), published circle, dry, a = 0.5"
    ct2 = R.r2(CARRANZA, model="ghb")
    ct3 = R.r3(CARRANZA, model="ghb")
    ce2 = R.r2(CARRANZA)
    ce3 = R.r3(CARRANZA)
    rows += [by_rounding(g, "2D GHB, 40 slices", ct2.fos, 2.0039, 4),
             info(g, "Published FoS (their method)", "", "2.01", "reference only"),
             info(g, "3D GHB, 20x20", _fmt(ct3.fos, 4), "about 2.14", "mesh to be confirmed"),
             info(g, "2D EJC (Moscarda), 40 slices", _fmt(ce2.fos, 4), "2.1214", "recorded comparison; mesh/inputs to be confirmed"),
             info(g, "3D EJC (Moscarda), 20x20", _fmt(ce3.fos, 4), "2.2678", "recorded comparison; mesh/inputs to be confirmed")]

    # statements on the Limitations page
    g = "6 Limitations statements (about)"
    rows.append(info(g, "F_H vs F_M difference, dry / wet (%)",
                     f"{(d3.fos - fm_d.fos) / d3.fos * 100:.3f} / {(w3.fos - fm_w.fos) / w3.fos * 100:.3f}"
                     if None not in (d3.fos, fm_d.fos, w3.fos, fm_w.fos) else "none", "0.14-0.27"))
    c0w = R.r3(tension="C0", wet=True).fos
    rows.append(info(g, "Wet 3D FoS drop with C0 (%)",
                     f"{(w3.fos - c0w) / w3.fos * 100:.1f}" if None not in (w3.fos, c0w) else "none", "about 25"))

    # finest mesh
    if slow:
        say("finest mesh (slow)")
        g = "3 Finest mesh (240x240, 5120 slices; slow)"
        f3d, f2d = R.r3(n3=240), R.r2(n2=5120)
        f3w, f2w = R.r3(n3=240, wet=True), R.r2(n2=5120, wet=True)
        rows += [by_rounding(g, "Dry: F 3D 240x240", f3d.fos, 2.7112, 4),
                 by_rounding(g, "Dry: F 2D 5120", f2d.fos, 2.4422, 4),
                 by_rounding(g, "Wet: F 3D 240x240", f3w.fos, 1.2776, 4,
                             "" if f3w.fos else f3w.reason[:160]),
                 by_rounding(g, "Wet: F 2D 5120", f2w.fos, 1.0916, 4, "" if f2w.fos else f2w.reason[:160])]
        g6 = "6 Limitations statements (about)"
        c10d, c10w = R.r3(n3=10), R.r3(n3=10, wet=True)
        for lab, c10, fin, tgt in (("dry", c10d, f3d, "about 2.6"), ("wet", c10w, f3w, "about 4.4")):
            val = (f"{(fin.fos - c10.fos) / fin.fos * 100:.2f}" if None not in (fin.fos, c10.fos) else "none")
            rows.append(info(g6, f"10x10 below 240x240, {lab} (%)", val, tgt))
    return rows


def summary(rows: list[Row]) -> dict:
    return {k: sum(r.status == k for r in rows) for k in ("PASS", "FAIL", "INFO")}


def as_markdown(rows: list[Row]) -> str:
    out = ["| Group | Item | Target | Computed | Difference | Tolerance | Status | Note |", "|" + "---|" * 8]
    for r in rows:
        out.append(f"| {r.group} | {r.item} | {r.target} | {r.computed} | {r.difference} | {r.tolerance} | "
                   f"**{r.status}** | {r.note} |")
    return "\n".join(out)


if __name__ == "__main__":
    args = sys.argv[1:]
    it = "workbook" if "--workbook" in args else "spec"
    sg = "project" if "--safeguard" in args else "none"
    rs = verification_rows(Runner(it, sg), slow="--slow" in args, progress=lambda s: print("...", s, file=sys.stderr))
    print(f"iteration = {it}, safeguard = {sg}\n")
    print(as_markdown(rs))
    print("\n", summary(rs))
