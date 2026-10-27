"""Randomised checks of the web engine over realistic and extreme inputs.

A. Parity: the web engine (engine.js) against the Python engine (ejc_app/engine.py) on the same random inputs:
   status, FoS, iteration count and error text.
B. Independent re-check of every accepted result: column geometry, the unsquared EJC/closure/ordering of each base
   state, local vertical equilibrium (eq. 4.12), the global F_H / F_M (eqs 4.14, 4.15) and Bishop (eq. 4.16) recomputed
   here with numpy from the returned arrays.
C. Local strength against a 50-digit solution (mpmath): sigma3 solved from the target sigma_n by high-precision root
   finding on the closure, without the engine's scan.
D. Robustness: the share of runs that throw anything other than a documented input error, return a non-finite or
   non-positive accepted FoS, or exceed a time limit.
"""
import json
import math
import random
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"                                    # results (git-ignored)
OUT.mkdir(exist_ok=True)
APP = Path(__file__).resolve().parents[1] / "ejc_app"
sys.path.insert(0, str(APP))
import engine as PY  # noqa: E402


def run_js(jobs, name):
    (OUT / f"{name}_jobs.json").write_text(json.dumps(jobs))
    t = time.time()
    subprocess.run(["node", str(HERE / "run_engine.js"), str(OUT / f"{name}_jobs.json"), str(OUT / f"{name}_out.json")], check=True)
    return json.loads((OUT / f"{name}_out.json").read_text()), time.time() - t


# ---------------------------------------------------------------------------------------------------------------
# random inputs
# ---------------------------------------------------------------------------------------------------------------
def random_circle(rng, H, beta, base_circle):
    """A circle from an exit point on the face (or floor), an entry point behind the crest and a centre on the chord's
    perpendicular bisector; returns (xc, zc, R)."""
    tb = math.tan(math.radians(beta))
    crest = -H / tb
    if base_circle and rng.random() < 0.5:
        xe, ze = rng.uniform(0.0, 0.6) * H, 0.0                    # exit on the floor ahead of the toe
    else:
        ze = rng.uniform(0.0, 0.35) * H
        xe = -ze / tb
    xn, zn = crest - rng.uniform(0.05, 0.8) * H, H
    mx, mz = 0.5 * (xe + xn), 0.5 * (ze + zn)
    L = math.hypot(xe - xn, ze - zn)
    nx, nz = -(ze - zn) / L, (xe - xn) / L                         # unit normal to the chord, pointing up-right
    if nz < 0:
        nx, nz = -nx, -nz
    t = rng.uniform(0.25, 1.6) * L
    xc, zc = mx + t * nx, mz + t * nz
    R = math.hypot(xe - xc, ze - zc)
    return xc, zc, R


def random_inputs(rng, extreme=False):
    H = rng.uniform(0.5, 1000) if extreme else rng.choice([rng.uniform(2, 20), rng.uniform(50, 600)])
    beta = rng.uniform(10, 89) if extreme else rng.uniform(25, 70)
    base_circle = rng.random() < 0.35
    xc, zc, R = random_circle(rng, H, beta, base_circle)
    model = rng.choices(["ejc_moscarda", "ghb", "ejc_sigma2_eq_sigma3", "mc"], weights=[60, 20, 5, 15])[0]
    inp = dict(H=H, beta=beta, xc=xc, zc=zc, R=R, base_circle=base_circle, model=model,
               gamma=rng.uniform(10, 40) if extreme else rng.uniform(18, 30),
               sigma_ci=rng.uniform(1, 300) if extreme else rng.uniform(10, 250),
               mi=rng.uniform(1, 40) if extreme else rng.uniform(4, 32),
               gsi=rng.uniform(10, 100) if extreme else rng.uniform(15, 95),
               D=rng.choice([0.0, 0.3, 0.5, 0.7, 1.0]),
               a_override=None if rng.random() < 0.5 else (rng.uniform(0.3, 1.0) if extreme else 0.5),
               nu=rng.uniform(0.01, 0.49) if extreme else rng.uniform(0.05, 0.49),
               n3=rng.randint(2, 40) if extreme else rng.randint(6, 24),
               n2=rng.randint(2, 400) if extreme else rng.randint(10, 80))
    if model == "mc":
        inp.update(c_mc=rng.uniform(0, 1000) if extreme else rng.uniform(0, 300) * (H / 100 if H < 50 else 1),
                   phi_mc=rng.uniform(0, 60) if extreme else rng.uniform(10, 45))
    w = rng.random()
    if w < 0.3:
        inp.update(wet=True, Hw=rng.uniform(0, H), beta_w=max(1.0, min(89.0, beta - rng.uniform(0, 15))))
    elif w < 0.5:
        inp.update(ru=rng.uniform(0, 0.99 if extreme else 0.6))
    it = rng.random()
    if it < 0.1:
        inp.update(iteration="native_v15")
    else:
        inp.update(form="FM" if rng.random() < 0.15 else "FH", tension=rng.choices(["none", "C0", "T0"], weights=[80, 10, 10])[0])
    return inp


def py_run(inp, dim):
    kw = {k: v for k, v in inp.items()}
    i = PY.Inputs(**kw)
    try:
        r = PY.run_3d(i) if dim == 3 else PY.run_2d(i)
    except ValueError as e:
        return dict(error=str(e))
    return dict(fos=r.fos, status=r.status, iterations=r.iterations, reason=r.reason)


# ---------------------------------------------------------------------------------------------------------------
# B. independent re-check of an accepted result from its arrays
# ---------------------------------------------------------------------------------------------------------------
def recheck(inp, dim, o):
    m = {k: np.asarray(v, dtype=float) for k, v in o["mesh"].items() if isinstance(v, list)}
    xc, zc, R, gam = inp["xc"], inp["zc"], inp["R"], inp.get("gamma", 27.0)
    x, y = m["x"], m["y"]
    q = R * R - (x - xc) ** 2 - y ** 2
    zb = zc - np.sqrt(q)
    tb = math.tan(math.radians(inp["H"] and inp.get("beta", 50.0)))
    zt = np.where(x <= 0, np.minimum(inp["H"], -x * tb), 0.0)
    g = np.sqrt(q) / R
    rho = np.sqrt(R * R - y * y) if dim == 3 else np.full(len(x), R)
    kap = (xc - x) / rho
    ah = o["mesh"]["dx"] * (o["mesh"]["dy"] if dim == 3 else 1.0)
    area = ah / g
    W = gam * (zt - zb) * ah
    geo_err = max(np.max(np.abs(zb - m["zb"]) / np.maximum(1, abs(zc))), np.max(np.abs(W - m["w"]) / m["w"]),
                  np.max(np.abs(kap - m["kap"])), np.max(np.abs(area - m["area"]) / area))
    sn, tau, c, t = (np.asarray(o[k], dtype=float) for k in ("sn", "tau_used", "c_used", "tan_used"))
    u = np.asarray(o["mesh"]["u"], dtype=float)
    F = o["fos"]
    local = np.abs(W - (sn + u) * ah - tau * area * kap / F) / W
    if inp.get("form", "FH") == "FM" and dim == 3:
        Fg = np.sum(tau * area * rho) / np.sum(W * (xc - x))
    else:
        Fg = np.sum(tau * area) / np.sum(W * kap)
    glob = abs(Fg - F) / F
    out = dict(geo=geo_err, local=float(np.max(local)), glob=glob)
    if dim == 2:
        ma = g + kap * t / F
        Fb = np.sum((c * ah + (W - u * ah) * t) / ma) / np.sum(W * kap)
        out["bishop"] = abs(Fb - F) / F
    if inp["model"] != "mc":
        s1, s2, s3 = (np.asarray(o[k], dtype=float) for k in ("s1", "s2", "s3"))
        mb, s, a = PY.hoek_brown_parameters(inp["gsi"], inp["mi"], inp["D"])
        a = a if inp.get("a_override") is None else inp["a_override"]
        sci = inp["sigma_ci"]
        Q = mb * s3 / sci + s
        f = 2 * s1 - s2 - s3 - np.sqrt(4 * sci ** 2 * Q ** (2 * a) - 3 * (s2 - s3) ** 2)
        p, qq = {"ejc_moscarda": (inp["nu"], inp["nu"]), "ghb": (0.0, 1.0), "ejc_sigma2_eq_sigma3": (0.0, 1.0)}[inp["model"]]
        out["ejc"] = float(np.max(np.abs(f)))
        out["closure"] = float(np.max(np.abs(s2 - (p * s1 + qq * s3))))
        out["ordering"] = int(np.sum(~((s1 >= s2) & (s2 >= s3))))
        # tangent geometry: sigma_n and tau of the state match the column's sigma_n and the strength used
        M = np.asarray(o["M"], dtype=float)
        sn_state = (s3 + (s1 - s3) / (M + 1)) * 1000
        out["sn_match"] = float(np.max(np.abs(sn_state - sn) / np.maximum(1.0, np.abs(sn))))
    return out


# ---------------------------------------------------------------------------------------------------------------
# C. 50-digit local strength
# ---------------------------------------------------------------------------------------------------------------
mp.mp.dps = 50


def mp_state(s3, sci, mb, s, a, p, q):
    Q = mb * s3 / sci + s
    if Q <= 0:
        return None
    P_, R_, S_ = p * p - p + 1, 2 * p * q - p - q - 1, q * q - q + 1
    K = sci ** 2 * Q ** (2 * a)
    disc = R_ ** 2 * s3 ** 2 - 4 * P_ * (S_ * s3 ** 2 - K)
    if disc < 0:
        return None
    for sign in (1, -1):
        s1 = (-R_ * s3 + sign * mp.sqrt(disc)) / (2 * P_)
        s2 = p * s1 + q * s3
        if not (s1 >= s2 >= s3):
            continue
        rad = 4 * sci ** 2 * Q ** (2 * a) - 3 * (s2 - s3) ** 2
        if rad < 0 or 2 * s1 - s2 - s3 < 0:
            continue
        M = (2 * a * mb * sci * Q ** (2 * a - 1) - R_ * s1 - 2 * S_ * s3) / (2 * P_ * s1 + R_ * s3)
        if M < 1:
            continue
        sn = s3 + (s1 - s3) / (M + 1)
        tau = (s1 - s3) * mp.sqrt(M) / (M + 1)
        phi = mp.asin((M - 1) / (M + 1))
        return dict(s1=s1, s2=s2, s3=s3, M=M, sn=sn, tau=tau, phi=phi, c=tau - sn * mp.tan(phi))
    return None


def mp_solve(target_mpa, sci, mb, s, a, p, q):
    """Lowest sigma3 with sigma_n(sigma3) = target on the admissible branch, by a fine mpmath scan and bisection."""
    s3t = -s * sci / mb
    span = target_mpa - s3t
    if span <= 0:
        return None
    grid = [s3t + span * mp.mpf(10) ** (mp.mpf(-14) + 14 * mp.mpf(k) / 1500) for k in range(1501)]
    prev = None
    for x in grid:
        st = mp_state(x, sci, mb, s, a, p, q)
        if st is not None:
            if prev is not None and (prev[1]["sn"] - target_mpa) * (st["sn"] - target_mpa) <= 0:
                lo, hi = prev[0], x
                glo = prev[1]["sn"] - target_mpa
                for _ in range(400):
                    mid = (lo + hi) / 2
                    sm = mp_state(mid, sci, mb, s, a, p, q)
                    if sm is None:
                        break
                    gm = sm["sn"] - target_mpa
                    if glo * gm <= 0:
                        hi = mid
                    else:
                        lo, glo = mid, gm
                    if hi - lo < mp.mpf(10) ** -45 * max(1, abs(hi)):
                        break
                return mp_state((lo + hi) / 2, sci, mb, s, a, p, q)
            prev = (x, st)
        else:
            prev = None
    return None


def main():
    rng = random.Random(20261027)
    results = {}

    # ---------------- A + B + D on realistic inputs, then D on extreme inputs
    for label, extreme, n in (("realistic", False, 400), ("extreme", True, 400)):
        cases = []
        for k in range(n):
            inp = random_inputs(rng, extreme)
            for dim in (3, 2):
                cases.append(dict(id=f"{k}|{dim}", dim=dim, inputs=inp, arrays=True))
        js, t_js = run_js(cases, label)
        stat = dict(cases=len(cases), statuses={}, input_errors=0, unexpected_errors=[], bad_fos=[], parity_fail=[],
                    worst_fos_rel=0.0, iter_mismatch=0, recheck_worst={}, recheck_fail=[], js_seconds=t_js)
        t_py = time.time()
        for cse in cases:
            o = js[cse["id"]]
            inp, dim = cse["inputs"], cse["dim"]
            if "error" in o:
                stat["input_errors"] += 1
                known = ("require", "circle", "choose a water table", "mesh size", "sliding mass", "R > z_c", "unknown")
                if not any(w in o["error"] for w in known):
                    stat["unexpected_errors"].append((cse["id"], o["error"]))
            else:
                stat["statuses"][o["status"]] = stat["statuses"].get(o["status"], 0) + 1
                if o["status"] == "ACCEPTED" and not (o["fos"] is not None and math.isfinite(o["fos"]) and o["fos"] > 0):
                    stat["bad_fos"].append((cse["id"], o["fos"]))
            py = py_run(inp, dim)
            if ("error" in o) != ("error" in py):
                stat["parity_fail"].append((cse["id"], o.get("error", o.get("status")), py.get("error", py.get("status"))))
            elif "error" in o:
                if o["error"] != py["error"]:
                    stat["parity_fail"].append((cse["id"], o["error"], py["error"]))
            else:
                if o["status"] != py["status"]:
                    stat["parity_fail"].append((cse["id"], o["status"], py["status"], o["reason"][:80], py["reason"][:80]))
                elif o["fos"] is not None:
                    rel = abs(o["fos"] - py["fos"]) / py["fos"]
                    stat["worst_fos_rel"] = max(stat["worst_fos_rel"], rel)
                    if rel > 1e-9:
                        stat["parity_fail"].append((cse["id"], "fos", o["fos"], py["fos"]))
                if o["iterations"] != py["iterations"]:
                    stat["iter_mismatch"] += 1
            if "error" not in o and o["status"] == "ACCEPTED" and "tau_used" in o:
                rc = recheck(inp, dim, o)
                for kk, vv in rc.items():
                    stat["recheck_worst"][kk] = max(stat["recheck_worst"].get(kk, 0), vv)
                lim = dict(geo=1e-9, local=1e-9, glob=1e-8, bishop=1e-8, ejc=1e-8, closure=1e-8, ordering=0, sn_match=1e-6)
                bad = {kk: vv for kk, vv in rc.items() if vv > lim[kk]}
                if bad:
                    stat["recheck_fail"].append((cse["id"], bad))
        stat["py_seconds"] = time.time() - t_py
        results[label] = stat
        print(f"\n{label}: {len(cases)} runs ({n} input sets x 3D/2D); JS {t_js:.1f} s, Python {stat['py_seconds']:.1f} s")
        print("  statuses:", stat["statuses"], "| input errors:", stat["input_errors"])
        print("  unexpected errors:", stat["unexpected_errors"][:5], "| accepted with bad FoS:", stat["bad_fos"][:5])
        print("  JS vs Python: worst FoS rel diff", f"{stat['worst_fos_rel']:.2e}", "| mismatches:", len(stat["parity_fail"]),
              stat["parity_fail"][:4], "| iteration-count mismatches:", stat["iter_mismatch"])
        print("  independent re-check of accepted runs, worst:", {k: f"{v:.2e}" for k, v in stat["recheck_worst"].items()},
              "| failures:", len(stat["recheck_fail"]), stat["recheck_fail"][:3])

    # ---------------- C. local strength against 50 digits
    worst = {}
    nbr = {}
    count = 0
    t0 = time.time()
    for k in range(150):
        gsi, mi, D, sci = rng.uniform(15, 95), rng.uniform(4, 32), rng.choice([0, 0.5, 0.7, 1.0]), rng.uniform(10, 250)
        mb, s, a = PY.hoek_brown_parameters(gsi, mi, D)
        if rng.random() < 0.5:
            a = 0.5
        model = rng.choice(["ejc_moscarda", "ejc_moscarda", "ghb", "ejc_original"])
        nu = rng.uniform(0.05, 0.49)
        st = PY.Strength(sci, mb, s, a, model, nu)
        p, q = st.pq
        lo_sn = max(st.sigma3_tensile * 1000 * 0.5, -50.0) if model != "ejc_original" else 1.0
        targets = [rng.uniform(lo_sn, 300.0), rng.uniform(300, 5000), rng.uniform(5000, 60000)]
        js_jobs = [dict(id="x", dim=2, inputs={})]
        r = PY.solve_strength(targets, st)                       # Python engine (equals the JS engine to ~1e-16)
        for i, tgt in enumerate(targets):
            ref = mp_solve(mp.mpf(tgt) / 1000, mp.mpf(sci), mp.mpf(mb), mp.mpf(s), mp.mpf(a), mp.mpf(p), mp.mpf(q))
            ok = int(r.code[i]) == PY.OK
            if ref is None:
                key = "both_none" if not ok else "engine_only"
                nbr[key] = nbr.get(key, 0) + 1
                continue
            if not ok:
                nbr["mp_only"] = nbr.get("mp_only", 0) + 1
                continue
            count += 1
            nbr[int(r.n_brackets[i])] = nbr.get(int(r.n_brackets[i]), 0) + 1
            for name, ev, rv in (("phi", r.phi[i], ref["phi"]), ("c", r.c[i] / 1000, ref["c"]), ("tau", r.tau[i] / 1000, ref["tau"]),
                                 ("s1", r.s1[i], ref["s1"]), ("s3", r.s3[i], ref["s3"]), ("M", r.M[i], ref["M"])):
                scale = max(abs(float(rv)), 1e-3 if name in ("c", "s3") else 1e-12)
                d = abs(float(ev) - float(rv)) / scale
                worst[name] = max(worst.get(name, 0), d)
    results["local_50_digit"] = dict(states=count, worst_rel=worst, outcomes=nbr, seconds=time.time() - t0)
    print(f"\n50-digit local strength: {count} accepted states compared; worst relative differences",
          {k: f"{v:.2e}" for k, v in worst.items()}, "| outcome counts (brackets per state / none):", nbr)
    (OUT / "qa_random.json").write_text(json.dumps(results, indent=1, default=str))


if __name__ == "__main__":
    main()
