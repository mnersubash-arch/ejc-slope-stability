"""Thesis claims that need their own calculation.

1. Table 5.2 at nu = 0.49: the minimum admissible sigma_n of the original closure, at 50 digits, from the joint
   triaxial boundary (sigma2 = sigma3) where the admissible branch starts.
2. "Tangent shear-strength gain 5.9-6.6%" (Section 5.4): tau_EJC / tau_GHB - 1 at equal sigma_n over the stress range.
3. Root order: can both roots of the squared EJC quadratic be admissible? (Appendix D tries the minus root first, the
   engine the plus root.)
4. Route 1 (Table 5.1, inherited angular calculation): an independent rebuild from eqs 4.8 and 4.9, to check the R1
   columns and the unsquared-EJC residuals -0.63 / -1.09 / -1.80 MPa.
"""
import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ejc_app"))
import engine as PY  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"                                    # results (git-ignored)
OUT.mkdir(exist_ok=True)
mp.mp.dps = 50
SCI, MB, S = mp.mpf("77.7"), 10 * mp.e ** (mp.mpf(-58) / 28), mp.e ** (mp.mpf(-58) / 9)
out = {}

# ---------------------------------------------------------------------------------------------------------------
# 1. original closure at nu = 0.49
# ---------------------------------------------------------------------------------------------------------------
def orig_state(s3, nu, a=mp.mpf("0.5")):
    p, q = mp.mpf(-1), 1 / nu
    Q = MB * s3 / SCI + S
    P_, R_, S_ = p * p - p + 1, 2 * p * q - p - q - 1, q * q - q + 1
    K = SCI ** 2 * Q ** (2 * a)
    disc = R_ ** 2 * s3 ** 2 - 4 * P_ * (S_ * s3 ** 2 - K)
    if disc < 0:
        return None
    best = None
    for sign in (1, -1):
        s1 = (-R_ * s3 + sign * mp.sqrt(disc)) / (2 * P_)
        s2 = p * s1 + q * s3
        if not (s1 >= s2 >= s3) or 2 * s1 - s2 - s3 < 0:
            continue
        M = (2 * a * MB * SCI * Q ** (2 * a - 1) - R_ * s1 - 2 * S_ * s3) / (2 * P_ * s1 + R_ * s3)
        if M < 1:
            continue
        sn = s3 + (s1 - s3) / (M + 1)
        best = dict(s1=s1, s2=s2, s3=s3, M=M, sn=sn)
        break
    return best


for nu_s in ("0.271", "0.45", "0.49"):
    nu = mp.mpf(nu_s)
    k = (1 - 2 * nu) / nu
    # joint triaxial boundary: sigma2 = sigma3 and the EJC (= GHB there): k sigma3 = sigma_ci Q^0.5  (a = 0.5)
    s3b = (SCI * MB + mp.sqrt((SCI * MB) ** 2 + 4 * k ** 2 * S * SCI ** 2)) / (2 * k ** 2)
    inside = orig_state(s3b * (1 + mp.mpf(10) ** -30), nu)
    # check that sigma_n increases along the admissible branch from the boundary
    st2 = orig_state(s3b * (1 + mp.mpf("1e-6")), nu)
    eng = PY.original_closure_minimum(77.7, float(MB), float(S), 0.5, float(nu))
    out[f"orig_{nu_s}"] = dict(boundary_s3=float(s3b), boundary_sn=float(inside["sn"]) if inside else None,
                               sn_increasing=bool(st2 and inside and st2["sn"] > inside["sn"]),
                               engine_s3=eng.sigma3_min, engine_sn=eng.sn_min)
    print(f"nu {nu_s}: boundary sigma3 {mp.nstr(s3b, 12)} MPa, sigma_n there {mp.nstr(inside['sn'], 12) if inside else None} MPa, "
          f"sigma_n increases inward: {out[f'orig_{nu_s}']['sn_increasing']}; engine minimum {eng.sn_min:.6f} MPa at sigma3 {eng.sigma3_min:.6f}")

# ---------------------------------------------------------------------------------------------------------------
# 2. tangent shear-strength gain over sigma_n
# ---------------------------------------------------------------------------------------------------------------
ste = PY.Inputs().strength()
sth = PY.Inputs(model="ghb").strength()
sn = np.concatenate([np.linspace(-97.9, 0, 40), np.geomspace(1, 5000, 200)])
te, th = PY.solve_strength(sn, ste), PY.solve_strength(sn, sth)
ok = te.accepted & th.accepted
gain = te.tau[ok] / th.tau[ok] - 1
out["tau_gain"] = dict(sn=sn[ok].tolist(), gain=gain.tolist())
for lo, hi in ((-97.9, 0), (0, 100), (100, 800), (800, 1610), (1610, 5000)):
    sel = (sn[ok] >= lo) & (sn[ok] <= hi)
    if sel.any():
        print(f"tau_EJC/tau_GHB - 1 for sigma_n {lo}..{hi} kPa: {gain[sel].min() * 100:.2f}% to {gain[sel].max() * 100:.2f}%")

# ---------------------------------------------------------------------------------------------------------------
# 3. both roots admissible?
# ---------------------------------------------------------------------------------------------------------------
both = 0
tested = 0
rng = np.random.default_rng(7)
for _ in range(400):
    gsi, mi, D = rng.uniform(10, 100), rng.uniform(1, 40), rng.choice([0, .5, .7, 1])
    mb, s, a = PY.hoek_brown_parameters(gsi, mi, D)
    sci = rng.uniform(1, 300)
    for model in ("ejc_moscarda", "ejc_original"):
        nu = rng.uniform(0.01, 0.49)
        st = PY.Strength(sci, mb, s, a if rng.random() < .5 else .5, model, nu)
        p, q = st.pq
        cP, cR, cS = PY.closure_constants(p, q)
        s3 = st.sigma3_tensile + (10 * sci) * np.geomspace(1e-12, 1, 400)
        Q = mb * s3 / sci + s
        K = sci ** 2 * Q ** (2 * st.a)
        disc = cR * cR * s3 * s3 - 4 * cP * (cS * s3 * s3 - K)
        good = disc >= 0
        sq = np.sqrt(np.where(good, disc, 0))
        adm = []
        for sign in (1, -1):
            s1 = (-cR * s3 + sign * sq) / (2 * cP)
            s2 = p * s1 + q * s3
            f = PY.ejc_residual(s1, s2, s3, Q, sci, st.a)
            adm.append(good & np.isfinite(f) & (np.abs(f) < 1e-8) & (s1 >= s2) & (s2 >= s3))
        tested += good.sum()
        both += int(np.sum(adm[0] & adm[1]))
out["both_roots_admissible"] = dict(points=int(tested), both=both)
print(f"root order: {both} of {tested} scanned states have both roots admissible")

# ---------------------------------------------------------------------------------------------------------------
# 4. Route 1 rebuilt from eqs 4.8 and 4.9 (Moscarda closure, a = 0.5, nu = 0.271)
# ---------------------------------------------------------------------------------------------------------------
sci, mb, s, a, nu = 77.7, float(MB), float(S), 0.5, 0.271
C, B = 1 - nu + nu * nu, 1 + 2 * nu - 2 * nu * nu


def route1(sn_mpa):
    def parts(phi):
        tau = sci * math.cos(phi) / 2 * ((mb * sn_mpa / sci + s) / (1 + math.sin(phi) / a)) ** a      # eq. 4.8
        t = (1 + math.sin(phi)) / math.cos(phi)
        s1, s3 = sn_mpa + tau * t, sn_mpa - tau / t
        s2 = nu * (s1 + s3)
        Q = mb * s3 / sci + s
        if Q <= 0:
            return None
        M = (2 * a * mb * sci * Q ** (2 * a - 1) + B * s1 - 2 * C * s3) / (2 * C * s1 - B * s3)         # eq. 4.5
        return tau, s1, s2, s3, Q, M

    def g(phi):                                                                                         # eq. 4.9
        p = parts(phi)
        return None if p is None else math.sin(phi) - (p[5] - 1) / (p[5] + 1)
    grid = np.linspace(math.radians(1), math.radians(89), 8801)
    roots = []
    prev = None
    for phi in grid:
        v = g(phi)
        if v is not None and prev is not None and prev[1] * v <= 0:
            lo, hi = prev[0], phi
            for _ in range(200):
                m = 0.5 * (lo + hi)
                gm = g(m)
                if (g(lo) > 0) == (gm > 0):
                    lo = m
                else:
                    hi = m
            roots.append(0.5 * (lo + hi))
        prev = (phi, v) if v is not None else None
    res = []
    for phi in roots:
        tau, s1, s2, s3, Q, M = parts(phi)
        f = 2 * s1 - s2 - s3 - math.sqrt(4 * sci ** 2 * Q ** (2 * a) - 3 * (s2 - s3) ** 2)
        res.append(dict(phi=math.degrees(phi), c=(tau - sn_mpa * math.tan(phi)) * 1000, tau=tau * 1000, ejc_residual=f))
    return res


out["route1"] = {}
for sn_k in (50, 400, 1600):
    r = route1(sn_k / 1000)
    out["route1"][sn_k] = r
    print(f"Route 1 at sigma_n {sn_k} kPa: " + "; ".join(f"phi {x['phi']:.2f} deg, c {x['c']:.1f} kPa, tau {x['tau']:.1f} kPa, "
                                                       f"EJC residual {x['ejc_residual']:.2f} MPa" for x in r))
(OUT / "qa_extra.json").write_text(json.dumps(out, indent=1, default=float))
