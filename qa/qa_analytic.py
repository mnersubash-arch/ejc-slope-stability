"""Independent checks of the web engine that need no external data.

1. Exact phi = 0 solutions. For a purely cohesive material on a circle (2D) Bishop's method is exact moment
   equilibrium, F = c R L / (gamma * int (x_c - x) dA). On a sphere rotating about the y axis the moment form is exact,
   F_M = c * int rho dA_s / (gamma * int (x_c - x) dV), and the simplified Hungr form has the continuum limit
   F_H = c * int dA_s / (gamma * int (x_c - x)/rho dV). The integrals are done here by adaptive quadrature, written from
   the geometry alone (no engine code), and compared with the engine as the mesh is refined.
2. The critical toe circle of a phi = 0 slope (Taylor's problem), found by minimising the exact solution, then solved by
   the engine on the same circle.
3. Metamorphic checks: invariance under geometric scaling with strength scaled alike, independence of pore pressure
   when phi = 0, linearity in c when phi = 0.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"                                    # results (git-ignored)
OUT.mkdir(exist_ok=True)


def run_engine(jobs):
    (OUT / "jobs.json").write_text(json.dumps(jobs))
    subprocess.run(["node", str(HERE / "run_engine.js"), str(OUT / "jobs.json"), str(OUT / "jobs_out.json")], check=True)
    return json.loads((OUT / "jobs_out.json").read_text())


# ---------------------------------------------------------------------------------------------------------------
# geometry written from the thesis text (Section 4.4.1), not from the engine
# ---------------------------------------------------------------------------------------------------------------
class Geo:
    def __init__(self, H, beta, xc, zc, R, base_circle=False):
        self.H, self.tb, self.xc, self.zc, self.R, self.base = H, math.tan(math.radians(beta)), xc, zc, R, base_circle
        self.crest = -H / self.tb

    def zt(self, x):
        return min(self.H, -x * self.tb) if x <= 0 else 0.0

    def zb(self, x, y):
        q = self.R ** 2 - (x - self.xc) ** 2 - y ** 2
        return self.zc - math.sqrt(q) if q > 0 else math.inf

    def h(self, x, y):
        return self.zt(x) - self.zb(x, y)

    def interval(self, y):
        """x-interval of the sliding mass on the line y (None if empty); the mass is one interval for these shapes."""
        xlo, xhi = self.xc - math.sqrt(max(self.R ** 2 - y ** 2, 0.0)), self.xc + math.sqrt(max(self.R ** 2 - y ** 2, 0.0))
        if not self.base:
            xhi = min(xhi, 0.0)
        xs = np.linspace(xlo, xhi, 4001)
        hv = np.array([self.h(x, y) for x in xs])
        pos = np.flatnonzero(hv > 0)
        if len(pos) == 0:
            return None
        i0, i1 = pos[0], pos[-1]

        def edge(a, b):                               # h(a) <= 0 < h(b) or the reverse
            fa = self.h(a, y)
            for _ in range(200):
                m = 0.5 * (a + b)
                if m in (a, b):
                    break
                if (self.h(m, y) > 0) == (fa > 0):
                    a, fa = m, self.h(m, y)
                else:
                    b = m
            return 0.5 * (a + b)
        x_a = xs[0] if i0 == 0 else edge(xs[i0 - 1], xs[i0])
        x_b = xs[-1] if i1 == len(xs) - 1 else edge(xs[i1 + 1], xs[i1])
        return x_a, x_b

    def line_integrals(self, y):
        """(int dA_s/dy, int rho dA_s/dy, int h (x_c-x) dx, int h (x_c-x)/rho dx) on the line y."""
        iv = self.interval(y)
        if iv is None:
            return np.zeros(4)
        a, b = iv
        rho = math.sqrt(self.R ** 2 - y ** 2)
        cuts = sorted({a, b} | ({self.crest} if a < self.crest < b else set()) | ({0.0} if a < 0.0 < b else set()))
        tot = np.zeros(4)
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            def g_inv(x):
                return self.R / math.sqrt(self.R ** 2 - (x - self.xc) ** 2 - y ** 2)
            tot[0] += integrate.quad(g_inv, lo, hi, epsabs=0, epsrel=1e-13, limit=200)[0]
            tot[2] += integrate.quad(lambda x: self.h(x, y) * (self.xc - x), lo, hi, epsabs=0, epsrel=1e-13, limit=200)[0]
        tot[1] = tot[0] * rho
        tot[3] = tot[2] / rho
        return tot

    def exact_3d(self, c, gamma):
        # y extent: largest y with a non-empty line
        lo, hi = 0.0, self.R
        for _ in range(100):
            m = 0.5 * (lo + hi)
            if self.interval(m) is None:
                hi = m
            else:
                lo = m
        ymax = lo
        # the line integrals vanish like (ymax - y)^p near ymax: integrate over t with y = ymax (1 - t^2)
        def comp(k):
            return integrate.quad(lambda t: self.line_integrals(ymax * (1 - t * t))[k] * 2 * t * ymax, 0, 1,
                                  epsabs=0, epsrel=1e-11, limit=200)[0]
        As, Arho, Mw, Mw_rho = comp(0), comp(1), comp(2), comp(3)
        return dict(F_M=c * Arho / (gamma * Mw), F_H=c * As / (gamma * Mw_rho), ymax=ymax)

    def exact_2d(self, c, gamma):
        L, Mw, _, _ = self.line_integrals(0.0)[[0, 2, 1, 3]][[0, 1, 2, 3]] if False else (None, None, None, None)
        t = self.line_integrals(0.0)
        # on y = 0: int dA_s = arc length L, int h (xc - x) dx = area moment; F = c L R / (gamma * moment)
        return c * t[0] * self.R / (gamma * t[2])


def main():
    report = {}
    c, gamma = 100.0, 27.0
    cases = {
        "main sphere (thesis), default grid": dict(H=360, beta=50, xc=355.7423115, zc=644.4192949, R=725.2669975),
        "Carranza-Torres circle (through the toe)": dict(H=360, beta=50, xc=207.28, zc=586.53, R=622.08),
        "D'Esterre DE-3.5-c4 circle, base-circle grid": dict(H=2.1, beta=15.945, xc=-3.46, zc=5.69, R=8.65, base_circle=True),
    }
    jobs = []
    N2, N3 = [40, 160, 640, 2560, 10000], [20, 40, 80, 160, 320]
    for name, g in cases.items():
        cc = c if g["H"] > 10 else 10.0
        inp = dict(g, gamma=gamma, model="mc", c_mc=cc, phi_mc=0.0)
        for n in N2:
            jobs.append(dict(id=f"{name}|2|{n}", dim=2, inputs=dict(inp, n2=n)))
        for n in N3:
            jobs.append(dict(id=f"{name}|3H|{n}", dim=3, inputs=dict(inp, n3=n, form="FH")))
            jobs.append(dict(id=f"{name}|3M|{n}", dim=3, inputs=dict(inp, n3=n, form="FM")))
    res = run_engine(jobs)
    for name, g in cases.items():
        cc = c if g["H"] > 10 else 10.0
        geo = Geo(**g)
        ex2 = geo.exact_2d(cc, gamma)
        ex3 = geo.exact_3d(cc, gamma)
        rows = dict(exact_2d=ex2, exact_3d_FM=ex3["F_M"], exact_3d_FH=ex3["F_H"],
                    engine_2d={n: res[f"{name}|2|{n}"]["fos"] for n in N2},
                    engine_3d_FM={n: res[f"{name}|3M|{n}"]["fos"] for n in N3},
                    engine_3d_FH={n: res[f"{name}|3H|{n}"]["fos"] for n in N3})
        report[name] = rows
        print(f"\n{name}")
        print(f"  2D exact {ex2:.10f}; engine: " + ", ".join(f"N={n}: {v:.10f} ({(v - ex2) / ex2 * 100:+.4f}%)" for n, v in rows['engine_2d'].items()))
        print(f"  3D F_M exact {ex3['F_M']:.10f}; engine: " + ", ".join(f"N={n}: {v:.8f} ({(v - ex3['F_M']) / ex3['F_M'] * 100:+.4f}%)" for n, v in rows['engine_3d_FM'].items()))
        print(f"  3D F_H limit {ex3['F_H']:.10f}; engine: " + ", ".join(f"N={n}: {v:.8f} ({(v - ex3['F_H']) / ex3['F_H'] * 100:+.4f}%)" for n, v in rows['engine_3d_FH'].items()))

    # ---- critical toe circles of phi = 0 slopes (Taylor's problem); stability number N = c / (gamma H F_min)
    print("\nCritical toe circles, phi = 0 (exact minimum vs engine on the same circle)")
    taylor = {}
    jobs = []
    for beta in (90.0, 75.0, 60.0):
        H = 10.0
        def F_of(p, beta=beta, H=H):
            xc, zc = p
            R = math.hypot(xc, zc)
            if R <= zc or zc <= H:
                return 1e9
            geo = Geo(H, beta - 1e-9 if beta == 90 else beta, xc, zc, R)
            try:
                return geo.exact_2d(c, gamma)
            except Exception:
                return 1e9
        best = optimize.minimize(F_of, x0=[0.3 * H, 1.6 * H], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-12, maxiter=4000))
        xc, zc = best.x
        R = math.hypot(xc, zc)
        Fmin = best.fun
        Nst = c / (gamma * H * Fmin)
        taylor[beta] = dict(xc=xc, zc=zc, R=R, F_exact=Fmin, N=Nst)
        for n in (40, 640, 5120):
            jobs.append(dict(id=f"taylor|{beta}|{n}", dim=2, inputs=dict(H=H, beta=min(beta, 89.9999), xc=xc, zc=zc, R=R * (1 + 1e-12), gamma=gamma,
                                                                    model="mc", c_mc=c, phi_mc=0.0, n2=n)))
    res = run_engine(jobs)
    for beta, t in taylor.items():
        eng = {n: res[f"taylor|{beta}|{n}"].get("fos") for n in (40, 640, 5120)}
        t["engine"] = eng
        print(f"  beta {beta:5.1f}: centre ({t['xc']:.4f}, {t['zc']:.4f}) m, R {t['R']:.4f} m; exact F {t['F_exact']:.6f}, "
              f"N = c/(gamma H F) = {t['N']:.4f}; engine " + ", ".join(f"N2={n}: {v:.6f}" for n, v in eng.items() if v))
    report["taylor"] = taylor

    # ---- metamorphic checks
    print("\nMetamorphic checks")
    base = dict(model="ejc_moscarda")
    jobs = []
    for k in (0.1, 1.0, 3.0):                  # lengths x k and sigma_ci x k: the EJC/GHB strength is homogeneous of degree 1
        inp = dict(H=360 * k, xc=355.7423115 * k, zc=644.4192949 * k, R=725.2669975 * k, Hw=360 * k, sigma_ci=77.7 * k)
        for wet in (False, True):
            for model in ("ejc_moscarda", "ghb"):
                jobs.append(dict(id=f"scale|{k}|{wet}|{model}|3", dim=3, inputs=dict(inp, wet=wet, model=model)))
                jobs.append(dict(id=f"scale|{k}|{wet}|{model}|2", dim=2, inputs=dict(inp, wet=wet, model=model)))
    for water in ("dry", "wet", "ru"):          # phi = 0: pore pressure cannot change the FoS
        kw = dict(model="mc", c_mc=150.0, phi_mc=0.0)
        if water == "wet":
            kw["wet"] = True
        if water == "ru":
            kw["ru"] = 0.4
        jobs.append(dict(id=f"u|{water}|3", dim=3, inputs=kw))
        jobs.append(dict(id=f"u|{water}|2", dim=2, inputs=kw))
    for cc in (50.0, 100.0, 400.0):             # phi = 0: F proportional to c
        jobs.append(dict(id=f"c|{cc}|3", dim=3, inputs=dict(model="mc", c_mc=cc, phi_mc=0.0)))
    res = run_engine(jobs)
    meta = {}
    for wet in (False, True):
        for model in ("ejc_moscarda", "ghb"):
            for dim in (3, 2):
                vals = [res[f"scale|{k}|{wet}|{model}|{dim}"]["fos"] for k in (0.1, 1.0, 3.0)]
                spread = (max(vals) - min(vals)) / vals[1]
                meta[f"scale {model} {'wet' if wet else 'dry'} {dim}D"] = dict(values=vals, rel_spread=spread)
                print(f"  scale x0.1/x1/x3, {model}, {'wet' if wet else 'dry'}, {dim}D: {vals} -> relative spread {spread:.2e}")
    for dim in (3, 2):
        vals = [res[f"u|{w}|{dim}"]["fos"] for w in ("dry", "wet", "ru")]
        meta[f"phi0 water {dim}D"] = dict(values=vals, rel_spread=(max(vals) - min(vals)) / vals[0])
        print(f"  phi = 0, dry / water table / r_u 0.4, {dim}D: {vals}")
    vals = [res[f"c|{cc}|3"]["fos"] for cc in (50.0, 100.0, 400.0)]
    meta["phi0 linear in c"] = dict(values=vals, ratios=[vals[1] / vals[0], vals[2] / vals[1]])
    print(f"  phi = 0, c = 50/100/400 kPa, 3D: {vals}; ratios {vals[1] / vals[0]:.15f}, {vals[2] / vals[1]:.15f}")
    report["metamorphic"] = meta
    (OUT / "qa_analytic.json").write_text(json.dumps(report, indent=1, default=float))


if __name__ == "__main__":
    main()
