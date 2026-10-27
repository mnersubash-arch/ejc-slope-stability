"""Replay qa_random.py's random stream and print the one local state the engine accepted but the 50-digit scan did not."""
import random
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import qa_random as QR  # noqa: E402

PY = QR.PY
rng = random.Random(20261027)
for extreme in (False, True):
    for k in range(400):
        QR.random_inputs(rng, extreme)
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
    r = PY.solve_strength(targets, st)
    for i, tgt in enumerate(targets):
        ok = int(r.code[i]) == PY.OK
        if not ok:
            continue
        ref = QR.mp_solve(mp.mpf(tgt) / 1000, mp.mpf(sci), mp.mpf(mb), mp.mpf(s), mp.mpf(a), mp.mpf(p), mp.mpf(q))
        if ref is None:
            print(dict(k=k, gsi=gsi, mi=mi, D=D, sci=sci, mb=mb, s=s, a=a, model=model, nu=nu, p=p, q=q, target_kPa=tgt))
            print("engine state:", dict(s1=float(r.s1[i]), s3=float(r.s3[i]), M=float(r.M[i]), phi=float(r.phi[i]),
                                        c=float(r.c[i]), tau=float(r.tau[i]), brackets=int(r.n_brackets[i])))
            # re-evaluate the engine's state at 50 digits
            stm = QR.mp_state(mp.mpf(float(r.s3[i])), mp.mpf(sci), mp.mpf(mb), mp.mpf(s), mp.mpf(a), mp.mpf(p), mp.mpf(q))
            print("50-digit state at the engine's sigma3:", None if stm is None else {kk: mp.nstr(vv, 15) for kk, vv in stm.items()})
            print("sigma3_tensile", st.sigma3_tensile)
