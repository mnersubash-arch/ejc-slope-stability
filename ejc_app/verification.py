"""Verification targets of the brief (section 3) and a runner that reports target, computed value, difference and
PASS / FAIL / INFO. No Streamlit imports; used by the tests, the Verification page and ``python verification.py``.

Tolerances (brief section 3): FoS given to 9 decimals must match to 1e-6; values given to 3-4 decimals must match to
rounding (|computed - target| <= half a unit in the last decimal shown); fixtures to the decimals shown.
INFO rows are reference or "about" values that are not acceptance targets. Without ``slow``, the fine-mesh rows are
listed as SKIPPED (their values are in VERIFICATION.md), so every run lists the same rows.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

import engine as E
import presets as PR

MAIN = E.Inputs()           # main sphere, Moscarda closure, EJC, dry, a = 0.5, 20x20 / 40, default path + projection
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
    """Runs and caches the analyses needed by the targets (default path with the Section 4.6 projection unless a row
    asks for the native v15 path)."""

    def __init__(self):
        self._cache: dict = {}

    def inputs(self, base: E.Inputs = MAIN, **kw) -> E.Inputs:
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


SKIP_TEXT = "skipped (slow), see VERIFICATION.md"


def skipped(group, item, target, tolerance) -> Row:
    return Row(group, item, target, SKIP_TEXT, "", tolerance, "SKIPPED")


# slow rows: (item, target, decimals, run). Listed as SKIPPED when the slow targets are not run.
CT_SLOW = [("2D EJC, 5120 slices (slow)", 2.121420, 6, lambda R: R.r2(CARRANZA, n2=5120)),
           ("3D GHB, 320x320 (slow)", 2.140776, 6, lambda R: R.r3(CARRANZA, model="ghb", n3=320)),
           ("3D EJC, 320x320 (slow)", 2.267835, 6, lambda R: R.r3(CARRANZA, n3=320))]
FINEST = [("Dry: F 3D 240x240", 2.7112, 4, lambda R: R.r3(n3=240)),
          ("Dry: F 2D 5120", 2.4422, 4, lambda R: R.r2(n2=5120)),
          ("Wet: F 3D 240x240", 1.2776, 4, lambda R: R.r3(n3=240, wet=True)),
          ("Wet: F 2D 5120", 1.0916, 4, lambda R: R.r2(n2=5120, wet=True))]
COARSE_INFO = [("dry", False, "about 2.6"), ("wet", True, "about 4.4")]


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
    note = "sigma_ci 120 MPa, m_b 3.5, s 0.002, a 0.505, sigma_n 5.7403 MPa"
    rows += [by_rounding(g, "Extra fixture, original closure (nu = 0.05): phi_i (deg)",
                         float(np.degrees(xo.phi[0])), 55.0720, 4, note),
             by_rounding(g, "Extra fixture, original closure (nu = 0.05): c_i (kPa)", float(xo.c[0]), 2079.5, 1, note),
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
    w3 = R.r3(wet=True)
    w2 = R.r2(wet=True)
    rows += [by_abs(g, "Dry, a = 0.5: F 3D", d3.fos, 2.706024342, 1e-6),
             by_abs(g, "Dry, a = 0.5: F 2D", d2.fos, 2.447910895, 1e-6),
             by_abs(g, "Wet, a = 0.5: F 3D", w3.fos, 1.272331405, 1e-6),
             by_abs(g, "Wet, a = 0.5: F 2D", w2.fos, 1.098069496, 1e-6),
             info(g, "Dry / wet, a = 0.5: 3D iterations and projections (default path)",
                  f"{d3.iterations} / {w3.iterations} iterations; {d3.projections} / {w3.projections} projections",
                  "shown, not tested")]
    gn = "3 Native v15 iteration path (start F = 1, Hungr explicit form, projection)"
    n3d, n3w = R.r3(iteration="native_v15"), R.r3(iteration="native_v15", wet=True)
    rows += [by_abs(gn, "Dry, a = 0.5: F 3D", n3d.fos, 2.706024342, 1e-6),
             by_exact(gn, "Dry, a = 0.5: 3D iterations", n3d.iterations, 34),
             by_abs(gn, "Wet, a = 0.5: F 3D", n3w.fos, 1.272331405, 1e-6),
             by_exact(gn, "Wet, a = 0.5: 3D iterations", n3w.iterations, 63),
             info(gn, "Projections dry / wet", f"{n3d.projections} / {n3w.projections}", "")]
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
                       (0.45, 2188.0, 0), (0.49, 59952.0, 0)):
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
    ct2, ct3 = R.r2(CARRANZA, model="ghb"), R.r3(CARRANZA, model="ghb")
    ce2, ce3 = R.r2(CARRANZA), R.r3(CARRANZA)
    rows += [by_rounding(g, "2D GHB, 40 slices (brief)", ct2.fos, 2.0039, 4),
             by_rounding(g, "2D GHB, 40 slices (thesis record)", ct2.fos, 2.003873, 6),
             by_rounding(g, "3D GHB, 20x20 (thesis value)", ct3.fos, 2.142, 3),
             by_rounding(g, "3D GHB, 20x20 (thesis record)", ct3.fos, 2.141937, 6),
             by_rounding(g, "2D EJC, 40 slices", ce2.fos, 2.122598, 6),
             by_rounding(g, "3D EJC, 20x20", ce3.fos, 2.269068, 6),
             info(g, "Published FoS (their method)", "", "2.01", "reference only")]
    if slow:
        say("Carranza-Torres fine meshes (slow)")
    for item, tgt, d, fn in CT_SLOW:
        rows.append(by_rounding(g, item, fn(R).fos, tgt, d) if slow else
                    skipped(g, item, f"{tgt:.{d}f}", f"rounding to {d} dp"))
    g = "3 Moscarda (2025), main sphere, nu = 0.271 (reference values, not targets)"
    for case, wet, a, rep_, sheet in ((1, False, 0.5, 2.24, 2.245), (2, True, 0.5, 1.57, 1.567),
                                      (3, False, 0.509923, 2.04, 2.041), (4, True, 0.509923, 1.13, 1.132)):
        ours = R.r2(wet=wet, a_override=a).fos
        rows.append(info(g, f"Case {case}: {'wet' if wet else 'dry'}, a = {a:g}, 2D FoS",
                         f"this engine (Route 2, 40 slices): {_fmt(ours, 4)}",
                         f"reported {rep_:.2f} (sheet at nu = 0.271: {sheet:.3f})",
                         "his own slice model and the inherited route"))

    # Law (2024): GHB 2D against the thesis reference code (lem_ref.py); reported values for reference only
    say("Law (2024)")
    g = "3 Cross-check: Law (2024) presets, GHB 2D against the reference code lem_ref.py"
    for p in PR.law_presets():
        for key, value in p.data["engine_check"].items():
            n2 = int(key.rsplit("_", 1)[1])
            rows.append(by_rounding(g, f"{p.data['id']}: GHB 2D, {n2} slices", R.r2(p.inputs, n2=n2).fos, value, 4,
                                    "read to 4 decimals (1.953 = 1.9530)" if value == 1.953 else ""))
    g = "3 Law (2024), reported values (reference only, not targets)"
    for p in PR.law_presets():
        rows.append(info(g, f"{p.data['id']}: 2D / 3D FoS",
                         f"this engine (GHB, 40 slices / 20x20): {_fmt(R.r2(p.inputs).fos, 4)} / "
                         f"{_fmt(R.r3(p.inputs).fos, 4)}", PR.values_text(p.data["reported"]), "her own models"))

    # Mohr-Coulomb option against the slope3d package's record (VBA engine, Appendix E sphere, c' = 100 kPa, phi' = 0)
    say("Mohr-Coulomb")
    g = "3 Cross-check: Mohr-Coulomb option (c' = 100 kPa, phi' = 0, main sphere, dry) against the slope3d record"
    mc = dict(model="mc", c_mc=100.0, phi_mc=0.0)
    rows += [by_abs(g, "F 3D, 20x20", R.r3(**mc).fos, 0.255354533635107, 1e-6, 15),
             by_abs(g, "F 2D, 40 slices", R.r2(**mc).fos, 0.190263799071471, 1e-6, 15)]

    # D'Esterre (2024): Mohr-Coulomb, r_u, base circles; reported values for reference only
    say("D'Esterre (2024)")
    g = "3 D'Esterre (2024): this engine against his reported values (reference only, not targets)"
    for p in PR.desterre_presets():
        rep = p.data["reported"]
        rows.append(info(g, p.data["id"],
                         f"2D {_fmt(R.r2(p.inputs).fos, 3)} (40 slices); 3D {_fmt(R.r3(p.inputs).fos, 3)} (20x20)",
                         f"Slide2 2D {rep['Slide2_2D']}; spreadsheet 3D {rep['DEsterre_spreadsheet_3D_Hungr']}; "
                         f"Slide3 3D {rep['Slide3_3D']}",
                         f"needs_floor_columns {str(p.data['needs_floor_columns']).lower()}"))

    # statements on the Limitations page
    g = "6 Limitations statements (about)"
    rows.append(info(g, "F_H vs F_M difference, dry / wet (%)",
                     f"{(d3.fos - fm_d.fos) / d3.fos * 100:.3f} / {(w3.fos - fm_w.fos) / w3.fos * 100:.3f}"
                     if None not in (d3.fos, fm_d.fos, w3.fos, fm_w.fos) else "none", "0.14-0.27"))
    c0w = R.r3(tension="C0", wet=True).fos
    rows.append(info(g, "Wet 3D FoS drop with C0 (%)",
                     f"{(w3.fos - c0w) / w3.fos * 100:.1f}" if None not in (w3.fos, c0w) else "none", "about 25"))

    # finest mesh
    g = "3 Finest mesh (240x240, 5120 slices; slow)"
    g6 = "6 Limitations statements (about)"
    if slow:
        say("finest mesh (slow)")
    for item, tgt, d, fn in FINEST:
        if slow:
            r = fn(R)
            note = f"{r.projections} projections on the path" if r.fos else r.reason[:160]
            rows.append(by_rounding(g, item, r.fos, tgt, d, note))
        else:
            rows.append(skipped(g, item, f"{tgt:.{d}f}", f"rounding to {d} dp"))
    for lab, wet, tgt in COARSE_INFO:
        item = f"10x10 below 240x240, {lab} (%)"
        if slow:
            fin, c10 = R.r3(n3=240, wet=wet), R.r3(n3=10, wet=wet)
            val = (f"{(fin.fos - c10.fos) / fin.fos * 100:.2f}" if None not in (fin.fos, c10.fos) else "none")
            rows.append(info(g6, item, val, tgt))
        else:
            rows.append(skipped(g6, item, tgt, "reference only"))
    return rows


def summary(rows: list[Row]) -> dict:
    return {k: sum(r.status == k for r in rows) for k in ("PASS", "FAIL", "SKIPPED", "INFO")}


def as_markdown(rows: list[Row]) -> str:
    out = ["| Group | Item | Target | Computed | Difference | Tolerance | Status | Note |", "|" + "---|" * 8]
    for r in rows:
        out.append(f"| {r.group} | {r.item} | {r.target} | {r.computed} | {r.difference} | {r.tolerance} | "
                   f"**{r.status}** | {r.note} |")
    return "\n".join(out)


if __name__ == "__main__":
    args = sys.argv[1:]
    rs = verification_rows(Runner(), slow="--slow" in args, progress=lambda s: print("...", s, file=sys.stderr))
    print(as_markdown(rs))
    print("\n", summary(rs))
