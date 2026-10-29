"""Baskoyak barite pit: a real slope failure, checked with the apps' engines.

The case (ejc_app/data/presets_baskoyak.json, presets "Baskoyak (real failure)"): a circular failure through one bench
in heavily jointed, dry schist (Sonmez & Ulusay 1999, Case 2). Li, Lo & Cassidy (2016) re-analysed it in SLIDE (Bishop,
the program's critical circle on the real section): F = 1.007 with GSI 16 and D 0.68. The section and the failure
circle are traced from Sonmez & Ulusay (1999) Fig. 7; the apps use a planar face at the overall angle.

1. An independent Bishop / Generalized Hoek-Brown code, written here from the textbook equations (no engine code),
   against the Python and web engines on the app's planar section: the failure circle and three other circles.
2. The three presets in both engines on the failure circle, and the critical-surface search (Python engine, 2D and 3D).
3. Against SLIDE: the critical surface with the SLIDE inputs.
4. The real section (bench, ground rising behind the crest) with the independent code: the failure circle and a critical
   circle search, which size the planar simplification.
5. Sonmez & Ulusay (1999) Table 4: their own parameter sets on the failure circle, by Bishop's method and by the
   ordinary method of slices, against the factors of safety they report from their program HOBRSLP.
6. The EJC (Moscarda closure) for this rock mass: the derivation of nu and the highest admissible sigma_n for each nu.
7. Back-analysis: the disturbance factor D that gives F = 1, GHB and EJC, 2D and 3D, critical surface and failure circle.

Results go to out/qa_baskoyak.json. About 10 minutes on four cores.
"""
import json
import math
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import optimize

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"                                    # results (git-ignored)
OUT.mkdir(exist_ok=True)
APP = HERE.parent / "ejc_app"
sys.path.insert(0, str(APP))
import engine as E  # noqa: E402
import presets as PR  # noqa: E402

REC = PR.BASKOYAK
KEYS = ("baskoyak_slide_ghb", "baskoyak_site_ghb", "baskoyak_site_ejc")
SLIDE_F = REC["reported"]["LLC2016_SLIDE"]["F"]


# ---------------------------------------------------------------------------------------------------------------
# independent Bishop / GHB code (figure coordinates: x into the slope from the toe, z up; the mass moves toward -x)
# ---------------------------------------------------------------------------------------------------------------
def ghb_state(s3, sci, mb, s, a):
    """sigma_n, tau and k = d sigma_1 / d sigma_3 on the GHB envelope at sigma_3 (Balmer's parametric form, kPa)."""
    q = mb * s3 / sci + s
    s1 = s3 + sci * q ** a
    k = 1.0 + a * mb * q ** (a - 1.0)
    return 0.5 * (s1 + s3) - 0.5 * (s1 - s3) * (k - 1.0) / (k + 1.0), (s1 - s3) * np.sqrt(k) / (k + 1.0), k


def ghb_tangent(sn, sci, mb, s, a):
    """Tangent cohesion and tan(phi) of the GHB envelope at sigma_n (kPa), sigma_3 found by bisection in
    log(sigma_3 - sigma_3t). At or below the tensile limit the base has no strength (counted)."""
    sn = np.asarray(sn, dtype=float)
    s3t = -s * sci / mb
    lo = np.full(sn.shape, math.log(1e-14 * sci))
    hi = np.full(sn.shape, math.log(10.0 * sci))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        up = ghb_state(s3t + np.exp(mid), sci, mb, s, a)[0] > sn
        hi, lo = np.where(up, mid, hi), np.where(up, lo, mid)
    v, tau, k = ghb_state(s3t + np.exp(0.5 * (lo + hi)), sci, mb, s, a)
    tanp = (k - 1.0) / (2.0 * np.sqrt(k))
    tens = sn <= s3t
    return np.where(tens, 0.0, tau - v * tanp), np.where(tens, 0.0, tanp), int(tens.sum())


class Section:
    """Ground surface through the points (x, z), flat beyond the end points."""

    def __init__(self, pts):
        p = np.asarray(pts, dtype=float)
        self.x, self.z = p[:, 0], p[:, 1]

    def __call__(self, x):
        return np.interp(x, self.x, self.z)


def mass_interval(sec, xc, zc, R, n=200001):
    xs = np.linspace(xc - R, xc + R, n)[1:-1]
    above = sec(xs) > zc - np.sqrt(np.maximum(R * R - (xs - xc) ** 2, 0.0))
    if not above.any():
        raise ValueError("the circle does not cut the ground")
    i0, i1 = int(np.argmax(above)), len(xs) - 1 - int(np.argmax(above[::-1]))
    if not above[i0:i1 + 1].all():
        raise ValueError("more than one mass")

    def edge(out, inn):
        for _ in range(100):
            m = 0.5 * (out + inn)
            if sec(m) > zc - math.sqrt(max(R * R - (m - xc) ** 2, 0.0)):
                inn = m
            else:
                out = m
        return 0.5 * (out + inn)

    return edge(xs[i0 - 1], xs[i0]), edge(xs[i1 + 1], xs[i1])


def slices(sec, xc, zc, R, gamma, n):
    x0, x1 = mass_interval(sec, xc, zc, R)
    b = (x1 - x0) / n
    x = x0 + (np.arange(n) + 0.5) * b
    root = np.sqrt(R * R - (x - xc) ** 2)
    h = np.maximum(sec(x) - (zc - root), 0.0)
    return b, gamma * h * b, (x - xc) / R, root / R, h


def bishop(sec, circle, gamma, sci, mb, s, a, n=2000):
    """Bishop's simplified method, strength from the GHB tangent at each base's sigma_n, repeated until F and sigma_n
    settle: F = sum[(c b + W tan phi)/m_a] / sum(W sin a), m_a = cos a + sin a tan phi / F,
    sigma_n = (W - c b tan a / F) cos a / (b m_a). Dry."""
    b, W, sa, ca, h = slices(sec, *circle, gamma, n)
    drive = float(np.sum(W * sa))
    sn, F = gamma * h * ca ** 2, 1.0
    for it in range(5000):
        c, tp, nt = ghb_tangent(sn, sci, mb, s, a)
        for _ in range(500):
            Fn = float(np.sum((c * b + W * tp) / (ca + sa * tp / F)) / drive)
            done = abs(Fn - F) < 1e-15 * Fn
            F = Fn
            if done:
                break
        new = (W - c * b * sa / ca / F) * ca / (b * (ca + sa * tp / F))
        move = float(np.max(np.abs(new - sn)))
        sn = new
        if it > 2 and move < 1e-12 * max(1.0, float(np.max(np.abs(sn)))):
            c, tp, nt = ghb_tangent(sn, sci, mb, s, a)
            return float(np.sum((c * b + W * tp) / (ca + sa * tp / F)) / drive), nt
    raise RuntimeError("Bishop iteration did not settle")


def ordinary(sec, circle, gamma, sci, mb, s, a, n=2000):
    """Ordinary method of slices (Fellenius): sigma_n = W cos a / l, F = sum(tau l)/sum(W sin a). Dry."""
    b, W, sa, ca, h = slices(sec, *circle, gamma, n)
    sn = W * ca * ca / b
    c, tp, _ = ghb_tangent(sn, sci, mb, s, a)
    return float(np.sum((c + sn * tp) * b / ca) / np.sum(W * sa))


def circle_from_exits(pf, pb, t):
    """Circle through the front exit pf and the back exit pb, centre on the chord's up-slope normal; t = 1 makes the arc
    vertical at the back exit, smaller t flatter circles."""
    (xf, zf), (xb, zb) = pf, pb
    L = math.hypot(xb - xf, zb - zf)
    nx, nz = -(zb - zf) / L, (xb - xf) / L
    mx, mz = 0.5 * (xf + xb), 0.5 * (zf + zb)
    d = (zb - mz) / nz + (1.0 / t - 1.0) * 0.5 * L
    xc, zc = mx + d * nx, mz + d * nz
    return xc, zc, math.hypot(xf - xc, zf - zc)


def critical_circle(sec, crest_x, gamma, sci, mb, s, a, n=400):
    """Lowest Bishop F over circles through a front exit (floor or face) and a back exit (face or ground behind the
    crest): a grid of starts, then Nelder-Mead from the four best; the minimum re-solved with 2000 slices."""
    best = {"F": math.inf}

    def f(p):
        xf, xb, t = p
        if not (-15.0 <= xf <= 0.6 * crest_x and 0.3 * crest_x <= xb <= crest_x + 40.0 and 0.03 <= t <= 1.0
                and xb > xf + 3.0):
            return 10.0
        try:
            circ = circle_from_exits((xf, float(sec(xf))), (xb, float(sec(xb))), t)
            F, _ = bishop(sec, circ, gamma, sci, mb, s, a, n)
        except (ValueError, RuntimeError):
            return 10.0
        if math.isfinite(F) and 0 < F < best["F"]:
            best.update(F=F, circle=circ, exits=(xf, xb))
        return F if math.isfinite(F) and F > 0 else 10.0

    grid = [(xf, xb, t) for xf in np.linspace(-8, 0.4 * crest_x, 6)
            for xb in np.linspace(0.6 * crest_x, crest_x + 30, 8) for t in (0.15, 0.35, 0.6, 0.9)]
    for _, p in sorted((f(p), p) for p in grid)[:4]:
        optimize.minimize(f, p, method="Nelder-Mead",
                          options={"xatol": 1e-3, "fatol": 1e-7, "maxfev": 600,
                                   "initial_simplex": [p, (p[0] + 2, p[1], p[2]), (p[0], p[1] + 3, p[2]),
                                                       (p[0], p[1], min(1.0, p[2] + 0.1))]})
    best["F2000"] = bishop(sec, best["circle"], gamma, sci, mb, s, a, 2000)[0]
    return best


# ---------------------------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------------------------
def planar_section(inp):
    """The app's section in figure coordinates: pit floor, planar face from the toe to the crest edge, flat crest."""
    crest = inp.H / math.tan(math.radians(inp.beta))
    return Section([(-200.0, 0.0), (0.0, 0.0), (crest, inp.H), (400.0, inp.H)])


def real_section():
    g = REC["geometry"]
    pts = [(-60.0, 0.0), tuple(g["toe"]), *map(tuple, g["bench"]), tuple(g["crest"]), tuple(g["ground_behind_crest"][1])]
    (x1, z1), (x2, z2) = g["ground_behind_crest"]
    pts.append((x2 + 200.0, z2 + 200.0 * (z2 - z1) / (x2 - x1)))
    return Section(pts)


def fig_circle(inp):
    return -inp.xc, inp.zc, inp.R


def run_js(jobs):
    (OUT / "baskoyak_jobs.json").write_text(json.dumps(jobs))
    subprocess.run(["node", str(HERE / "run_engine.js"), str(OUT / "baskoyak_jobs.json"),
                    str(OUT / "baskoyak_jobs_out.json")], check=True)
    return json.loads((OUT / "baskoyak_jobs_out.json").read_text())


def inputs_json(inp):
    return {k: (v if not isinstance(v, float) or math.isfinite(v) else None) for k, v in vars(inp).items()}


def search_job(key):
    inp = PR.BY_KEY[key].inputs
    res = E.critical_search(inp, dims=(2, 3))
    r2, r3 = res["r2"], res["r3"]
    return key, {"fos2": r2["fos"], "fos2_160": r2["check"]["fos"], "fos3": r3["fos"], "fos3_40": r3["check"]["fos"],
                 "circle2_fig": [-r2["xc"], r2["zc"], r2["R"]], "front2": r2["front"], "back2": r2["back"],
                 "withheld2": r2["counts"]["withheld"]}


def back_job(args):
    key, D = args
    inp = PR.BY_KEY[key].inputs.with_(D=D)
    res = E.critical_search(inp, dims=(2, 3))
    return key, D, res["r2"]["fos"], res["r3"]["fos"], E.run_2d(inp.with_(n2=160)).fos, E.run_3d(inp.with_(n3=40)).fos


def d_for_one(Ds, F):
    roots = np.roots(np.polyfit(Ds, np.asarray(F) - 1.0, 3))
    ok = [r.real for r in roots if abs(r.imag) < 1e-9 and min(Ds) - 0.1 <= r.real <= max(Ds) + 0.1]
    return ok[0] if ok else None


def main():
    report = {}
    P = {k: PR.BY_KEY[k] for k in KEYS}
    slide_in, site_in = P["baskoyak_slide_ghb"].inputs, P["baskoyak_site_ghb"].inputs

    # 1. independent Bishop code against the engines on the app's planar section
    print("1. Independent Bishop/GHB code against the engines, app's planar section (GHB)")
    sec = planar_section(slide_in)
    circles = [fig_circle(slide_in), (2.74, 41.47, 41.21), (0.0, 45.0, 45.5), (5.0, 35.0, 36.0)]
    worst, jobs, rows = 0.0, [], []
    for inp in (slide_in, site_in):
        mb, s, a = inp.hb
        for i, circ in enumerate(circles):
            for n in (40, 640, 2000):
                q = inp.with_(xc=-circ[0], zc=circ[1], R=circ[2], n2=n)
                ind, _ = bishop(sec, circ, inp.gamma, inp.sigma_ci * 1000.0, mb, s, a, n)
                py = E.run_2d(q).fos
                worst = max(worst, abs(py / ind - 1))
                rows.append((inp.gsi, circ, n, ind, py))
                jobs.append({"id": f"{len(jobs)}", "dim": 2, "inputs": inputs_json(q)})
    js = run_js(jobs)
    worst_js = max(abs(js[j["id"]]["fos"] / r[4] - 1) for j, r in zip(jobs, rows))
    for gsi, circ, n, ind, py in rows[::3]:
        print(f"   GSI {gsi:4.1f} circle ({circ[0]:5.2f}, {circ[1]:5.2f}) R {circ[2]:5.2f}, {n:4d} slices: "
              f"independent {ind:.10f}, Python engine {py:.10f}")
    print(f"   largest relative difference: independent vs Python {worst:.1e}; web vs Python {worst_js:.1e} "
          f"({len(rows)} analyses)")
    report["independent_vs_engines"] = {"rel_python": worst, "rel_web": worst_js, "analyses": len(rows)}

    # 2. the presets in both engines, failure circle and critical surface
    print("\n2. Presets: failure circle (Python and web) and critical-surface search (Python, 2D and 3D)")
    jobs = []
    for k, p in P.items():
        for dim, n in ((2, 40), (2, 160), (3, 20), (3, 40)):
            q = p.inputs.with_(**({"n2": n} if dim == 2 else {"n3": n}))
            jobs.append({"id": f"{k}|{dim}|{n}", "dim": dim, "inputs": inputs_json(q)})
    js = run_js(jobs)
    with Pool(3) as pool:
        crit = dict(pool.map(search_job, KEYS))
    report["presets"] = {}
    for k, p in P.items():
        f = {}
        for dim, n in ((2, 40), (2, 160), (3, 20), (3, 40)):
            q = p.inputs.with_(**({"n2": n} if dim == 2 else {"n3": n}))
            py = (E.run_2d(q) if dim == 2 else E.run_3d(q)).fos
            web = js[f"{k}|{dim}|{n}"]["fos"]
            f[f"{dim}D_{n}"] = py
            f[f"web_rel_{dim}D_{n}"] = abs(web / py - 1)
        c = crit[k]
        report["presets"][k] = {"failure_circle": f, "critical": c}
        print(f"   {k}: failure circle 2D {f['2D_40']:.4f} (160 slices {f['2D_160']:.4f}), 3D {f['3D_20']:.4f} "
              f"(40x40 {f['3D_40']:.4f}); critical 2D {c['fos2']:.4f} (160 slices {c['fos2_160']:.4f}), 3D "
              f"{c['fos3']:.4f} (40x40 {c['fos3_40']:.4f}); web/Python <= "
              f"{max(v for kk, v in f.items() if kk.startswith('web')):.1e}")

    # 3. SLIDE
    c = report["presets"]["baskoyak_slide_ghb"]["critical"]
    print(f"\n3. SLIDE (Li et al. 2016, Bishop, real section): F = {SLIDE_F}; this engine, planar section, critical "
          f"surface: {c['fos2']:.4f} (40 slices), {c['fos2_160']:.4f} (160 slices): "
          f"{(c['fos2_160'] / SLIDE_F - 1) * 100:+.2f} %")
    report["slide"] = {"published": SLIDE_F, "engine_critical_160": c["fos2_160"]}

    # 4. the real section with the independent code
    print("\n4. Real section (bench, ground rising behind the crest), independent code")
    real, crest_x = real_section(), REC["geometry"]["crest"][0]
    report["real_section"] = {}
    for k in ("baskoyak_slide_ghb", "baskoyak_site_ghb"):
        inp = P[k].inputs
        mb, s, a = inp.hb
        on = bishop(real, fig_circle(inp), inp.gamma, inp.sigma_ci * 1000.0, mb, s, a)[0]
        cr = critical_circle(real, crest_x, inp.gamma, inp.sigma_ci * 1000.0, mb, s, a)
        planar = report["presets"][k]
        report["real_section"][k] = {"failure_circle": on, "critical": cr["F2000"], "critical_circle": cr["circle"]}
        print(f"   {k}: failure circle {on:.4f} (planar {planar['failure_circle']['2D_160']:.4f}); critical "
              f"{cr['F2000']:.4f} (planar {planar['critical']['fos2_160']:.4f}); critical circle centre "
              f"({cr['circle'][0]:.2f}, {cr['circle'][1]:.2f}) R {cr['circle'][2]:.2f}")

    # 5. Sonmez & Ulusay (1999) Table 4
    print("\n5. Sonmez & Ulusay (1999), Case 2 on the failure circle, real section: their sets, Bishop and the ordinary "
          "method of slices (OMS)")
    t4 = REC["reported"]["SU1999_Table4_Case2"]
    g, mi, df = REC["site"]["GSI"]["SU1999"], REC["site"]["m_i"], t4["df"]
    bm, bs = t4["approach_4_limit_equilibrium"]["b_m"], t4["approach_4_limit_equilibrium"]["b_s"]
    sets = {"approach 1": (mi * math.exp((g - 100) / 28), 0.0, 0.65 - g / 200, t4["approach_1"]),
            "approach 2": (mi * math.exp((g * df - 100) / 28), 0.0, 0.65 - g * df / 200, t4["approach_2"]),
            "approach 3": (mi * math.exp((g - 100) / (28 * df)), 0.0, 0.65 - g / 200, t4["approach_3"]),
            "approach 4 (limit eq.)": (mi * math.exp((g - 100) / bm), math.exp((g - 100) / bs), 0.65 - g / 200, 1.0)}
    circ = fig_circle(slide_in)
    report["su1999"] = {}
    for name, (mb, s, a, pub) in sets.items():
        fb = bishop(real, circ, 22.2, 5200.0, mb, s, a)[0]
        fo = ordinary(real, circ, 22.2, 5200.0, mb, s, a)
        report["su1999"][name] = {"mb": mb, "s": s, "a": a, "published": pub, "bishop": fb, "ordinary": fo}
        print(f"   {name:22s} mb {mb:.4f} s {s:.2e} a {a:.3f}: published {pub:.2f}; Bishop {fb:.3f} "
              f"({(fb / pub - 1) * 100:+.1f} %); OMS {fo:.3f} ({(fo / pub - 1) * 100:+.1f} %)")

    # 6. EJC: nu and admissibility
    print("\n6. EJC (Moscarda closure), site rock mass (GSI 16.5, D 0.7)")
    nu_i = 0.25 * (1 + math.exp(-0.2 * site_in.sigma_ci))
    rmr = site_in.gsi + 5
    rmqr = 100 * rmr / (rmr + 1.1 * (100 - rmr))
    factor = 2.5 - 1.5 * rmqr / 100
    print(f"   nu: intact 0.25(1 + e^(-0.2 x {site_in.sigma_ci})) = {nu_i:.3f}; RMR ~ {rmr:.1f}, RMQR = {rmqr:.1f}, factor "
          f"{factor:.3f}; rock mass {nu_i * factor:.3f} (above 0.5); with an intact 0.15-0.20: "
          f"{0.15 * factor:.2f}-{0.20 * factor:.2f}. Preset: {P['baskoyak_site_ejc'].inputs.nu}")
    mb, s, a = site_in.hb
    limits = {}
    for nu in (0.271, 0.30, 0.33, 0.36, 0.40, 0.45, 0.49):
        st = E.Strength(site_in.sigma_ci, mb, s, a, "ejc_moscarda", nu)
        grid = np.geomspace(1.0, 5000.0, 160)
        ok = E.solve_strength(grid, st).accepted
        if ok.all():
            limits[nu] = None
            continue
        i = int(np.argmin(ok))
        lo, hi = (grid[i - 1], grid[i]) if i > 0 else (0.0, grid[0])
        for _ in range(40):
            m = 0.5 * (lo + hi)
            lo, hi = (m, hi) if E.solve_strength([m], st).accepted[0] else (lo, m)
        limits[nu] = lo
    print("   highest admissible sigma_n (kPa): " +
          ", ".join(f"nu {nu}: {'all of 1-5000' if v is None else f'{v:.0f}'}" for nu, v in limits.items()))
    r = E.run_2d(P["baskoyak_site_ghb"].inputs.with_(
        **dict(zip(("xc", "zc", "R"), (-report["presets"]["baskoyak_site_ghb"]["critical"]["circle2_fig"][0],
                                       *report["presets"]["baskoyak_site_ghb"]["critical"]["circle2_fig"][1:]))),
        n2=160))
    print(f"   largest sigma_n on the critical surface (GHB, 160 slices): {float(np.max(r.sn)):.0f} kPa")
    report["ejc"] = {"nu_intact": nu_i, "rmqr": rmqr, "factor": factor, "limits": limits,
                     "critical_sn_max": float(np.max(r.sn))}

    # 7. back-analysis of D
    print("\n7. Back-analysis: D for F = 1 (site inputs)")
    Ds = (0.6, 0.7, 0.8, 0.9)
    with Pool(4) as pool:
        backs = pool.map(back_job, [(k, D) for k in ("baskoyak_site_ghb", "baskoyak_site_ejc") for D in Ds])
    report["back_analysis"] = {}
    for k in ("baskoyak_site_ghb", "baskoyak_site_ejc"):
        rows = sorted(b for b in backs if b[0] == k)
        out = {name: d_for_one(Ds, [row[i] for row in rows])
               for i, name in ((2, "critical_2D"), (3, "critical_3D"), (4, "failure_circle_2D"), (5, "failure_circle_3D"))}
        report["back_analysis"][k] = {"D": out, "F": [list(row[1:]) for row in rows]}
        print(f"   {k}: " + ", ".join(f"{name.replace('_', ' ')} {v:.3f}" for name, v in out.items()))

    (OUT / "qa_baskoyak.json").write_text(json.dumps(report, indent=1, default=float))
    print(f"\nwritten {OUT / 'qa_baskoyak.json'}")


if __name__ == "__main__":
    main()
