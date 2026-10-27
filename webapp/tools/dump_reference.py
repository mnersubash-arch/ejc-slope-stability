"""Reference results from the Python engine (ejc_app/engine.py) for checking the JavaScript port.

Usage (from the repository root): python webapp/tools/dump_reference.py [path to ejc_app] [output dir]
Defaults: ejc_app/ next to webapp/, output webapp/tools/ref/. Writes reference.json (analyses, strength fixtures,
original-closure minima, slip extents) for verify_engine.js and verification_rows.json (the verification.py table,
without the slow rows) for verify_checks.js. Takes about 1.5 minutes.
"""
import dataclasses
import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parents[1] / "ejc_app"
OUT = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else HERE / "ref"
sys.path.insert(0, str(APP))

import numpy as np  # noqa: E402

import engine as E  # noqa: E402
import presets as PR  # noqa: E402


def num(v):
    if v is None:
        return None
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    v = float(v)
    if math.isnan(v):
        return "nan"
    if math.isinf(v):
        return "inf" if v > 0 else "-inf"
    return v


def arr(a):
    return None if a is None else [num(x) for x in np.asarray(a).ravel()]


def inputs_dict(inp):
    return {k: num(v) if not isinstance(v, str) else v for k, v in dataclasses.asdict(inp).items()}


def mass(inp, r):
    """Whole-mass forces (E.mass_forces) of a run that reached a final state, else None."""
    if r.sn is None or r.tau_used is None or (r.fos is None and r.f_trial is None):
        return None
    mf = E.mass_forces(inp, r)
    return dict(weight=num(mf["weight"]), vertical_max_rel=num(mf["vertical_max_rel"]),
                totals={k: {c: num(v[c]) for c in ("x", "z")} for k, v in mf["totals"].items()},
                moments={k: num(v) for k, v in mf["moments"].items()},
                method={k: v if isinstance(v, str) else num(v) for k, v in mf["method"].items()})


def lem(r):
    return dict(status=r.status, reason=r.reason, fos=num(r.fos), f_trial=num(r.f_trial), iterations=r.iterations,
                projections=r.projections, count=r.mesh.count, volume=num(r.mesh.volume),
                x0=num(r.mesh.x0), x1=num(r.mesh.x1),
                checks=[dict(name=c.name, value=num(c.value), passed=bool(c.passed)) for c in r.checks],
                history=[[num(v) for v in row] for row in r.history],
                sn=arr(r.sn), tau_used=arr(r.tau_used), projection_limit=num(r.projection_limit),
                global_res=num(r.global_res), bishop_res=num(r.bishop_res))


def run_case(name, inp, dims=(3, 2)):
    out = dict(name=name, inputs=inputs_dict(inp))
    for dim in dims:
        try:
            r = E.run_3d(inp) if dim == 3 else E.run_2d(inp)
            out[f"r{dim}"] = lem(r)
            out[f"m{dim}"] = mass(inp, r)
        except ValueError as e:
            out[f"r{dim}"] = dict(error=str(e))
    return out


MAIN = E.Inputs()
WET = MAIN.with_(wet=True)
CT = MAIN.with_(xc=207.28, zc=586.53, R=622.08)

cases = []
add = lambda name, inp, dims=(3, 2): cases.append(run_case(name, inp, dims))  # noqa: E731
t0 = time.time()
for wet in (False, True):
    w = "wet" if wet else "dry"
    base = MAIN.with_(wet=wet)
    add(f"main {w}", base)
    add(f"main {w} a=0.509923", base.with_(a_override=0.509923))
    add(f"main {w} a from GSI", base.with_(a_override=None))
    add(f"ghb {w}", base.with_(model="ghb"))
    add(f"ghb {w} a=0.509923", base.with_(model="ghb", a_override=0.509923))
    add(f"s2=s3 {w}", base.with_(model="ejc_sigma2_eq_sigma3"))
    add(f"FM {w}", base.with_(form="FM"), (3,))
    add(f"C0 {w}", base.with_(tension="C0"))
    add(f"T0 {w}", base.with_(tension="T0"))
    add(f"native {w}", base.with_(iteration="native_v15"), (3, 2))
    add(f"safeguard none {w}", base.with_(safeguard="none"))
    for n3, n2 in ((10, 80), (30, 160), (40, 20)):
        add(f"mesh {n3}/{n2} {w}", base.with_(n3=n3, n2=n2))
add("original closure", MAIN.with_(model="ejc_original"))
add("ct ghb", CT.with_(model="ghb"))
add("ct ejc", CT)
add("mc c100 phi0", MAIN.with_(model="mc", c_mc=100.0, phi_mc=0.0))
add("mc c50 phi35 wet", WET.with_(model="mc", c_mc=50.0, phi_mc=35.0))
add("ru 0.3 ejc", MAIN.with_(ru=0.3))
add("base circle main", MAIN.with_(base_circle=True))
for nu in (0.05, 0.15, 0.35, 0.49):
    add(f"nu {nu}", MAIN.with_(nu=nu))
add("gsi 60 D 0.5", MAIN.with_(gsi=60.0, D=0.5, a_override=None))
add("sigma_ci 30 mi 20", MAIN.with_(sigma_ci=30.0, mi=20.0))
add("circle behind toe", MAIN.with_(zc=600.0, R=640.0))
add("R <= zc", MAIN.with_(R=600.0))
FACE = MAIN.with_(H=20.0, beta=60.0, xc=2.0, zc=30.0, R=math.hypot(-5.0 / math.tan(math.radians(60.0)) - 2.0, -25.0),
                  gamma=25.0, model="mc", c_mc=30.0, phi_mc=35.0)       # face circle, lowest point above the toe level
add("face circle, default grid", FACE)
add("face circle, base circle", FACE.with_(base_circle=True))
add("ru and wet", WET.with_(ru=0.2))
add("N = 2", MAIN.with_(n3=2, n2=2))
for p in PR.law_presets():
    add(f"preset {p.key}", p.inputs)
    add(f"preset {p.key} 640", p.inputs.with_(n2=640), (2,))
for p in PR.desterre_presets():
    add(f"preset {p.key}", p.inputs)
print(f"{len(cases)} cases in {time.time() - t0:.1f} s", file=sys.stderr)

# strength fixtures and principal states
fixtures = []
targets = [-1.5, -0.5, 0.0, 1e-6, 0.5, 3.0, 50.0, 400.0, 1600.0, 5740.3, 20000.0, 1.0e5]
mb, s, a = MAIN.hb
strengths = [("moscarda", E.Strength(77.7, mb, s, 0.5, "ejc_moscarda", 0.271)),
             ("moscarda a gsi", E.Strength(77.7, mb, s, a, "ejc_moscarda", 0.271)),
             ("original 0.271", E.Strength(77.7, mb, s, 0.5, "ejc_original", 0.271)),
             ("original 0.05", E.Strength(77.7, mb, s, 0.5, "ejc_original", 0.05)),
             ("ghb", E.Strength(77.7, mb, s, 0.5, "ghb")),
             ("s2=s3", E.Strength(77.7, mb, s, 0.5, "ejc_sigma2_eq_sigma3")),
             ("extra original", E.Strength(120.0, 3.5, 0.002, 0.505, "ejc_original", 0.05)),
             ("extra moscarda", E.Strength(120.0, 3.5, 0.002, 0.505, "ejc_moscarda", 0.05)),
             ("mc", E.Strength(0, 0, 0, 0, "mc", c_mc=21.6, phi_mc=20.0))]
for name, st in strengths:
    r = E.solve_strength(targets, st)
    s3s = [st.sigma3_tensile * 0.999 if st.model != "mc" else -0.1, 0.0, 0.01, 1.0, 10.0, 34.376, 100.0]
    ps = E.principal_states(s3s, st) if st.model != "mc" else None
    fixtures.append(dict(
        name=name, strength=dict(sigma_ci=st.sigma_ci, mb=st.mb, s=st.s, a=st.a, model=st.model, nu=st.nu,
                                 c_mc=st.c_mc, phi_mc=st.phi_mc),
        targets=targets, sn=arr(r.sn), phi=arr(r.phi), c=arr(r.c), tau=arr(r.tau), M=arr(r.M), s1=arr(r.s1),
        s3=arr(r.s3), code=arr(r.code), n_brackets=arr(r.n_brackets),
        projection_limit=num(E.projection_limit(st)),
        ps=None if ps is None else dict(s3=s3s, s1=arr(ps.s1), M=arr(ps.M), sn=arr(ps.sn), tau=arr(ps.tau),
                                        code=arr(ps.code))))

ocm = []
for nu in (0.05, 0.10, 0.20, 0.271, 0.35, 0.45, 0.49):
    m0 = E.original_closure_minimum(77.7, mb, s, 0.5, nu)
    ocm.append(dict(nu=nu, sigma3_min=num(m0.sigma3_min), sn_min=num(m0.sn_min), phi_deg=num(m0.phi_deg),
                    found=m0.found))

extents = []
for p in PR.desterre_presets():
    extents.append(dict(key=p.key, inputs=inputs_dict(p.inputs), extent=[num(v) for v in E.slip_extent(p.inputs)]))

ref = dict(hb=[num(v) for v in E.hoek_brown_parameters(42.0, 10.0, 0.0)], cases=cases, fixtures=fixtures,
           original_closure=ocm, extents=extents)
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "reference.json").write_text(json.dumps(ref))

# the verification table (fast rows; slow rows are listed as SKIPPED), for comparison with checks.js
import verification as V  # noqa: E402
rows = V.verification_rows(V.Runner(), slow=False)
(OUT / "verification_rows.json").write_text(json.dumps([dataclasses.asdict(r) for r in rows], ensure_ascii=False))
print(V.summary(rows), f"total {time.time() - t0:.1f} s", file=sys.stderr)
