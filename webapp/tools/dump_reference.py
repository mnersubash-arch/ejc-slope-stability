"""Reference results from the Python engine (ejc_app/engine.py) for checking the JavaScript port.

Usage (from the repository root): python webapp/tools/dump_reference.py [path to ejc_app] [output dir]
Defaults: ejc_app/ next to webapp/, output webapp/tools/ref/. Writes reference.json (analyses, strength fixtures,
original-closure minima, slip extents, critical-surface searches) for verify_engine.js and verification_rows.json (the
verification.py table, without the slow rows) for verify_checks.js. Takes about 5 minutes.
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
    return {k: v if isinstance(v, str) else [[num(a) for a in q] for q in v] if isinstance(v, tuple) else num(v)
            for k, v in dataclasses.asdict(inp).items()}


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
                x0=num(r.mesh.x0), x1=num(r.mesh.x1), n=r.mesh.n, ny=r.mesh.ny, dx=num(r.mesh.dx), dy=num(r.mesh.dy),
                checks=[dict(name=c.name, value=num(c.value), passed=bool(c.passed)) for c in r.checks],
                history=[[num(v) for v in row] for row in r.history],
                sn=arr(r.sn), tau_used=arr(r.tau_used), projection_limit=num(r.projection_limit),
                global_res=num(r.global_res), bishop_res=num(r.bishop_res),
                projection_limits=None if r.projection_limits is None else [num(v) for v in r.projection_limits],
                zone=None if r.mesh.zone is None else int(np.sum(r.mesh.zone)),
                crack=None if r.mesh.crack is None else {k: num(r.mesh.crack[k]) for k in ("x", "level", "P", "MH", "MM")},
                ext=None if r.mesh.ext is None else {k: num(r.mesh.ext[k]) for k in ("Ex", "EH", "EM", "H", "M")})


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
for p in PR.baskoyak_presets():
    add(f"preset {p.key}", p.inputs)
# 3D columns N x N_y (n3y; None is the N x N grid)
add("columns 20 x 20 given", MAIN.with_(n3y=20), (3,))
add("columns 20 x 7", MAIN.with_(n3y=7))
add("columns 20 x 1", MAIN.with_(n3y=1), (3,))
add("columns 20 x 40 wet", WET.with_(n3y=40), (3,))
add("columns 13 x 50 FM", MAIN.with_(n3=13, n3y=50, form="FM"), (3,))
add("columns 30 x 11 ghb", MAIN.with_(model="ghb", n3=30, n3y=11), (3,))
add("columns 20 x 37 native", MAIN.with_(n3y=37, iteration="native_v15"), (3,))
add("columns 20 x 37 C0 wet", WET.with_(n3y=37, tension="C0"), (3,))
add("columns 7 x 90 ru", MAIN.with_(ru=0.3, n3=7, n3y=90), (3,))
add("columns baskoyak ejc 20 x 45", PR.BY_KEY["baskoyak_site_ejc"].inputs.with_(n3y=45), (3,))
add("columns desterre 24 x 9", PR.BY_KEY["desterre_DE-1.5-c3-ru0.5"].inputs.with_(n3=24, n3y=9), (3,))
add("columns n3y = 0", MAIN.with_(n3y=0))
add("columns n3y = 2.5", MAIN.with_(n3y=2.5), (3,))
print(f"{len(cases)} cases in {time.time() - t0:.1f} s", file=sys.stderr)

# geometry, water and loading additions
BAS = PR.BY_KEY["baskoyak_site_ghb"].inputs
REAL = ((0.0, 0.0), (-11.9, 8.5), (-14.1, 8.5), (-30.48, 18.99))          # Baskoyak real section (traced)
CREST = math.degrees(math.atan((20.68 - 18.99) / (49.3 - 30.48)))
PIEZO = ((50.0, 0.0), (0.0, 0.0), (-150.0, 150.0), (-600.0, 300.0))
add("add: profile baskoyak real section", BAS.with_(profile=REAL, crest_slope=CREST))
add("add: profile baskoyak, 2000 slices", BAS.with_(profile=REAL, crest_slope=CREST, n2=2000), (2,))
add("add: profile ejc", PR.BY_KEY["baskoyak_site_ejc"].inputs.with_(profile=REAL, crest_slope=CREST))
add("add: crest slope 10", MAIN.with_(crest_slope=10.0))
add("add: crest slope 25 base circle", MAIN.with_(crest_slope=25.0, base_circle=True))
add("add: concave 800", MAIN.with_(plan_radius=800.0))
add("add: concave 3000 wet", WET.with_(plan_radius=3000.0), (3,))
add("add: convex -1500", MAIN.with_(plan_radius=-1500.0), (3,))
add("add: convex too small", MAIN.with_(plan_radius=-200.0))
add("add: elong 2", MAIN.with_(elong=2.0, n3y=30), (3,))
add("add: elong 0.5", MAIN.with_(elong=0.5), (3,))
add("add: mid 1000", MAIN.with_(mid_length=1000.0, n3y=40), (3,))
add("add: elong 1.5 mid 500 FM", MAIN.with_(elong=1.5, mid_length=500.0, form="FM"), (3,))
add("add: zone 10 ghb D 0.7", MAIN.with_(model="ghb", D=0.7, a_override=None, d_zone=10.0))
add("add: zone 10 ejc D 0.7", MAIN.with_(D=0.7, d_zone=10.0))
add("add: zone 0 ghb", MAIN.with_(model="ghb", D=0.7, a_override=None, d_zone=0.0))
add("add: zone mc ignored", MAIN.with_(model="mc", c_mc=50.0, phi_mc=35.0, d_zone=10.0), (3,))
add("add: piezo", MAIN.with_(wet=True, piezo=PIEZO))
add("add: piezo concave", MAIN.with_(wet=True, piezo=PIEZO, plan_radius=1500.0), (3,))
add("add: crack 20 dry", MAIN.with_(crack_depth=20.0))
add("add: crack 20 water 0.5", MAIN.with_(crack_depth=20.0, crack_water=0.5))
add("add: crack 20 water 1 FM", MAIN.with_(crack_depth=20.0, crack_water=1.0, form="FM"), (3,))
add("add: crack too deep", MAIN.with_(crack_depth=60.0))
add("add: kh 0.1", MAIN.with_(kh=0.1))
add("add: kh 0.1 FM wet", WET.with_(kh=0.1, form="FM"), (3,))
add("add: kh 0.1 native", MAIN.with_(kh=0.1, iteration="native_v15"))
add("add: kh 0.2 ghb C0", MAIN.with_(kh=0.2, model="ghb", tension="C0"))
add("add: all, baskoyak", BAS.with_(profile=REAL, crest_slope=5.0, plan_radius=300.0, elong=1.5, mid_length=20.0,
                                     d_zone=3.0, kh=0.05, crack_depth=1.0, crack_water=0.5, n3y=25))
add("add: all, thesis wet", WET.with_(plan_radius=2000.0, elong=1.2, d_zone=25.0, D=0.5, kh=0.05, crack_depth=10.0,
                                       crack_water=0.5, crest_slope=5.0))
for bad, kw in (("profile not from the toe", {"profile": ((1.0, 0.0), (-10.0, 5.0))}),
                ("profile back and forth", {"profile": ((0.0, 0.0), (-10.0, 5.0), (-5.0, 8.0))}),
                ("profile descends", {"profile": ((0.0, 0.0), (-10.0, 5.0), (-20.0, 4.0))}),
                ("profile one point", {"profile": ((0.0, 0.0),)}),
                ("crest slope 90", {"crest_slope": 90.0}), ("plan radius 0", {"plan_radius": 0.0}),
                ("elong 0", {"elong": 0.0}), ("mid -1", {"mid_length": -1.0}), ("zone -1", {"d_zone": -1.0}),
                ("piezo forward", {"piezo": ((0.0, 0.0), (10.0, 5.0))}), ("crack 0", {"crack_depth": 0.0}),
                ("crack water 2", {"crack_water": 2.0}), ("kh 1", {"kh": 1.0}), ("kh -0.1", {"kh": -0.1})):
    add(f"add error: {bad}", MAIN.with_(**kw))
print(f"{len(cases)} cases with the additions", file=sys.stderr)

# column helpers: columns across y for other N (mesh studies, the search's finer mesh), plan size, counts for a size
columns = []
for inp, n in ((MAIN, 40), (MAIN.with_(n3y=7), 40), (MAIN.with_(n3y=7), 10), (MAIN.with_(n3y=7), 3),
               (MAIN.with_(n3y=37), 320), (MAIN.with_(n3=13, n3y=50), 7), (MAIN.with_(n3y=1), 2),
               (PR.BY_KEY["baskoyak_site_ejc"].inputs.with_(n3y=45), 9)):
    columns.append(dict(inputs=inputs_dict(inp), n=n, across=E.columns_across(inp, n),
                        size=[num(v) for v in E.column_size(inp)]))
for inp, dx, dy in ((MAIN, 20.0, 20.0), (MAIN, 10.0, 10.0), (MAIN, 1000.0, 1000.0), (MAIN, *E.column_size(MAIN)),
                    (MAIN.with_(n3y=37), *E.column_size(MAIN.with_(n3y=37))),
                    (PR.BY_KEY["baskoyak_site_ejc"].inputs, 1.0, 0.5),
                    (PR.BY_KEY["desterre_DE-1.5-c3-ru0.5"].inputs, 0.25, 0.25),
                    (MAIN.with_(elong=1.5, mid_length=300.0), 20.0, 20.0)):
    columns.append(dict(inputs=inputs_dict(inp), dx=num(dx), dy=num(dy), counts=list(E.counts_for_size(inp, dx, dy))))
for inp, n in ((MAIN.with_(elong=1.5, mid_length=300.0, n3y=30), 20),):
    columns.append(dict(inputs=inputs_dict(inp), n=n, across=E.columns_across(inp, n),
                        size=[num(v) for v in E.column_size(inp)]))

# studies: EJC against GHB, sweeps, back-analysis, Monte Carlo, the random stream
studies = {"compare": [], "sweep": [], "back": [], "mc": [], "rng": []}
for inp in (MAIN, WET.with_(kh=0.05), BAS):
    studies["compare"].append(dict(inputs=inputs_dict(inp), result=E.compare_models(inp)))
for inp, name, vals, models, dims in ((MAIN, "gsi", [30.0, 42.0, 60.0], ("ejc_moscarda", "ghb"), (3, 2)),
                                      (MAIN, "kh", [0.0, 0.1, 0.2], ("ghb",), (2,)),
                                      (BAS, "D", [0.5, 0.7, 0.9], ("ghb",), (3, 2))):
    studies["sweep"].append(dict(inputs=inputs_dict(inp), name=name, values=vals, models=list(models), dims=list(dims),
                                 rows=E.parameter_sweep(inp, name, vals, models, dims)))
for inp, name, lo, hi, target, dim in ((BAS, "D", 0.5, 0.9, 1.0, 2), (MAIN.with_(model="ghb"), "kh", 0.0, 0.5, 1.5, 3),
                                       (MAIN, "gsi", 20.0, 30.0, 1.0, 2)):
    studies["back"].append(dict(inputs=inputs_dict(inp), name=name, lo=lo, hi=hi, target=target, dim=dim,
                                result=E.back_analysis(inp, name, lo, hi, target, dim)))
mc_vars = [("gsi", 42.0, 5.0, 10.0, 100.0), ("sigma_ci", 77.7, 15.0, 1.0, 300.0), ("D", 0.1, 0.1, 0.0, 1.0)]
studies["mc"].append(dict(inputs=inputs_dict(MAIN.with_(model="ghb")), variables=[list(v) for v in mc_vars], n=12, seed=7,
                          dims=[2], result=E.monte_carlo(MAIN.with_(model="ghb"), mc_vars, 12, seed=7, dims=(2,))))
for seed in (0, 1, 7, 123456789, 4294967295):
    r = E.mulberry32(seed)
    studies["rng"].append(dict(seed=seed, values=[r() for _ in range(20)]))

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



def search_dict(res):
    """Every field of a critical_search result (the clouds included)."""
    def nums(a):
        return None if a is None else [num(v) for v in a]

    def one(r):
        if r is None:
            return None
        out = {k: num(r.get(k)) for k in ("fos", "xc", "zc", "R", "coverage", "evaluations")}
        out.update(p=nums(r.get("p")), front=nums(r.get("front")), back=nums(r.get("back")), counts=r["counts"],
                   skipped_why=r["skipped_why"], edges=r["edges"], message=r["message"],
                   check=None if r["check"] is None else {k: v if isinstance(v, str) else num(v)
                                                          for k, v in r["check"].items()},
                   runs=[dict(p0=nums(q["p0"]), f0=num(q["f0"]), p=nums(q["p"]), f=num(q["f"])) for q in r["runs"]],
                   cloud={k: nums(r["cloud"][k]) for k in ("xc", "zc", "R", "fos")})
        return out
    return dict(prescribed={k: v if isinstance(v, str) else num(v) for k, v in res["prescribed"].items()},
                settings=res["settings"], domain=res["domain"], r2=one(res["r2"]), r3=one(res["r3"]))


searches = []
t1 = time.time()
TAYLOR75 = E.Inputs(H=10.0, beta=75.0, gamma=27.0, model="mc", c_mc=100.0, phi_mc=0.0, n2=640, xc=0.40, zc=7.20, R=7.21)
for name, inp, dims, opts in (
        ("taylor 75, vertical-back prescribed circle", TAYLOR75, (2,), {}),
        ("taylor 75, back limit 1", TAYLOR75, (2,), {"back": 1.0}),
        ("taylor 60", TAYLOR75.with_(beta=60.0), (2,), {}),
        ("desterre DE-1.5-c1-ru0, 2D and 3D", PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs, (2, 3), {}),
        ("desterre DE-2.5-c4-ru0.5, front limit 1", PR.BY_KEY["desterre_DE-2.5-c4-ru0.5"].inputs, (2,), {"front": 1.0}),
        ("main ghb dry, 2D and 3D", MAIN.with_(model="ghb"), (2, 3), {}),
        ("baskoyak site GHB, base circle, 2D and 3D", PR.BY_KEY["baskoyak_site_ghb"].inputs, (2, 3), {}),
        ("desterre DE-1.5-c1-ru0, 20 x 9 columns, 2D and 3D",
         PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs.with_(n3y=9), (2, 3), {}),
        ("desterre DE-1.5-c1-ru0, 8 x 170 columns (finer mesh at the 320 cap), 2D and 3D",
         PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs.with_(n3=8, n3y=170), (2, 3), {}),
        ("baskoyak site GHB on the real section (profile, crest slope), 2D and 3D",
         BAS.with_(profile=REAL, crest_slope=CREST), (2, 3), {}),
        ("desterre DE-1.5-c1-ru0 with k_h 0.1 and a tension crack 0.5 m half full, 2D",
         PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs.with_(kh=0.1, crack_depth=0.5, crack_water=0.5), (2,), {})):
    searches.append(dict(name=name, inputs=inputs_dict(inp), dims=list(dims), opts=opts,
                         result=search_dict(E.critical_search(inp, dims=dims, **opts))))
print(f"{len(searches)} searches in {time.time() - t1:.1f} s", file=sys.stderr)

faces = []
DE = PR.BY_KEY["desterre_DE-1.5-c1-ru0"].inputs
for inp in (MAIN, MAIN.with_(xc=25.0, zc=250.0, R=251.0), MAIN.with_(xc=500.0, zc=-10.0, R=400.0), TAYLOR75, DE,
            DE.with_(xc=-2.0, zc=-0.5, R=3.0), MAIN.with_(base_circle=True),
            MAIN.with_(xc=25.0, zc=250.0, R=251.0, crack_depth=5.0), MAIN.with_(xc=25.0, zc=250.0, R=251.0, crest_slope=8.0),
            BAS.with_(profile=REAL, crest_slope=CREST)):
    faces.append(dict(inputs=inputs_dict(inp), faces={k: num(v) for k, v in E.vertical_faces(inp).items()}))

ref = dict(hb=[num(v) for v in E.hoek_brown_parameters(42.0, 10.0, 0.0)], cases=cases, fixtures=fixtures,
           original_closure=ocm, extents=extents, searches=searches, faces=faces, columns=columns, studies=studies)
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "reference.json").write_text(json.dumps(ref))

# the verification table (fast rows; slow rows are listed as SKIPPED), for comparison with checks.js
import verification as V  # noqa: E402
rows = V.verification_rows(V.Runner(), slow=False)
(OUT / "verification_rows.json").write_text(json.dumps([dataclasses.asdict(r) for r in rows], ensure_ascii=False))
print(V.summary(rows), f"total {time.time() - t0:.1f} s", file=sys.stderr)
