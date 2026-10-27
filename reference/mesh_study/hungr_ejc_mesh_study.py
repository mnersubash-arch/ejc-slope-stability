"""Independent numerical reimplementation used for the mesh-convergence study of
GENG5512_Thesis_Calculation_Package.xlsx (conventional 2D Hoek-Brown strength, closure 3).
External numerical REIMPLEMENTATION (separate code, same equations); not a result of the thesis EJC closure.
Same equations as the workbook, separate code. Run:  python hungr_ejc_mesh_study.py
Outputs: mesh_study_records.csv (one row per case/dimension/mesh) and reproduction_check.txt."""
import math, csv, sys
import numpy as np
SCI, MI, GSI, D = 77.7, 10.0, 42.0, 0.0
MB = MI*math.exp((GSI-100)/(28-14*D)); S = math.exp((GSI-100)/(9-3*D)); T = S*SCI/MB          # T = HB tensile limit magnitude (MPa)
H, L, GAM, GW = 360.0, 1/math.tan(math.radians(50)), 27.0, 9.81
# acceptance thresholds = the workbook's default Interface tolerances (P16 angle, P17 stress, P18 F); edit both together
TOL_ANG, TOL_MPA, TOL_F, TOL_EQ = 1e-10, 1e-8, 1e-8, 1e-9
SURFACES = {"AppendixE": (355.742311546406, 644.419294864537, 725.2669975187362),   # Moscarda (2025) App. E, mirrored to workbook axes
            "Published": (207.28, 586.53, 622.08)}                                 # Moscarda (2025) Table 4.3 / Carranza-Torres & Hormazabal (2020)
def hb_table(a):
    """exact 2D Hoek-Brown tangent states, parametric in sigma3' (MPa): sorted sigma_n' -> phi (rad), c (MPa)"""
    s3 = -T + np.logspace(-9, math.log10(T + 400), 400000)
    Q = MB*s3/SCI + S; s1 = s3 + SCI*Q**a; M = 1 + a*MB*Q**(a-1); p = (M-1)/(M+1); phi = np.arcsin(p)
    sn = (s1+s3)/2 - (s1-s3)/2*p; tau = (s1-s3)/2*np.cos(phi); o = np.argsort(sn)
    return sn[o], phi[o], (tau - sn*np.tan(phi))[o], s3[o], a
def exact_state(sn_MPa, tab):
    """exact 2D Hoek-Brown tangent state at sigma_n' (MPa): Newton on sigma3' from the table start; returns phi (rad), c (MPa)"""
    snT, phT, cT, s3T, a = tab
    s3 = np.interp(sn_MPa, snT, s3T)
    def snof(x):
        Q = MB*x/SCI + S; s1 = x + SCI*Q**a; M = 1 + a*MB*Q**(a-1); p = (M-1)/(M+1); return (s1+x)/2 - (s1-x)/2*p
    for _ in range(6):
        h = 1e-7*np.maximum(1.0, np.abs(s3)); f = snof(s3) - sn_MPa; df = (snof(s3+h) - snof(s3-h))/(2*h)
        s3 = np.maximum(s3 - f/df, -T*(1-1e-12))
    Q = MB*s3/SCI + S; s1 = s3 + SCI*Q**a; M = 1 + a*MB*Q**(a-1); p = (M-1)/(M+1); phi = np.arcsin(p)
    snx = (s1+s3)/2 - (s1-s3)/2*p; tau = (s1-s3)/2*np.cos(phi); c = tau - snx*np.tan(phi)
    ok = np.isfinite(phi) & np.isfinite(c)
    return np.where(ok, phi, np.interp(sn_MPa, snT, phT)), np.where(ok, c, np.interp(sn_MPa, snT, cT))
def elements(dim, n, surf, wet):
    XC, ZC, R = SURFACES[surf]; x0 = XC - math.sqrt(R*R - ZC*ZC); x1 = XC - R          # toe-side exit -> back of the sphere
    if dim == 3:
        dx, dy = (x0-x1)/n, R/n
        X, Y = np.meshgrid(x0-(np.arange(n)+0.5)*dx, (np.arange(n)+0.5)*dy, indexing="ij"); X, Y = X.ravel(), Y.ravel()
        q = R*R-(X-XC)**2-Y*Y
    else:
        dx = (x0-x1)/n; X = x0-(np.arange(n)+0.5)*dx; Y = np.zeros_like(X); q = R*R-(X-XC)**2
    k = q > 0; X, Y, q = X[k], Y[k], q[k]; zb = ZC-np.sqrt(q)
    top = np.where(X >= 0, 0.0, np.minimum(H, -X/L)); h = top-zb; k = h > 0; X, Y, q, zb, h = X[k], Y[k], q[k], zb[k], h[k]
    if dim == 3:
        g = np.sqrt(q)/R; kap = (XC-X)/np.sqrt(R*R-Y*Y); A = dx*dy/g; W = GAM*h*dx*dy; Ah = A*g             # Hungr column
    else:
        g = np.sqrt(q)/R; kap = (XC-X)/R; A = dx/g; W = GAM*h*dx; Ah = np.full_like(X, dx)                    # Bishop slice (unit width)
    zw = np.where(X >= 0, 0.0, np.minimum(H, -X*math.tan(math.radians(50)))) if wet else np.full_like(X, -1e9)
    u = GW*np.maximum(0.0, zw-zb)
    return dict(W=W, A=A, g=g, kap=kap, u=u, Ah=Ah, vol=float((h*(dx*dy if dim == 3 else dx)).sum()))
def solve(e, tab, tol=1e-12, maxit=50000):
    """fixed point F <-> sigma_n'; trial stresses beyond the tensile limit are moved just inside it (path safeguard only)"""
    snT, phT, cT, s3T, a_ = tab; lim = snT[0]*1000*(1 - 1e-9); W, A, g, kap, u, Ah = (e[k] for k in ("W", "A", "g", "kap", "u", "Ah"))
    sn = (W/Ah - u)*g*g; F = 1.0; nproj = 0
    for it in range(maxit):
        phi, c = exact_state(sn/1000, tab); c = 1000*c; t = np.tan(phi)
        SF = (c*Ah + (W - u*Ah)*t)/(g + kap*t/F); Fn = SF.sum()/(W*kap).sum()
        raw = (SF - c*A)/(A*t); proj = raw <= lim; nproj += int(proj.sum()); new = np.where(proj, -T*1000*(1 - 1e-3), raw)
        if abs(Fn-F) <= tol*abs(Fn) and np.max(np.abs(new-sn)) < 1e-9:
            # ---- final-state acceptance at the reported F with the strength actually used ----
            SFf = (c*Ah + (W - u*Ah)*t)/(g + kap*t/Fn); snf = (SFf - c*A)/(A*t)                  # equilibrium sigma_n' at final F (kPa)
            inside = bool(np.all(snf > snT[0]*1000)); finite = bool(np.all(np.isfinite(snf)) and np.all(np.isfinite(phi)))
            phi_chk, _ = exact_state(snf/1000, tab); cons = float(np.max(np.abs(phi_chk - phi)))                   # strength consistency (rad)
            tau = c + snf*t; p_ = snf + tau*t; r_ = tau/np.cos(phi); s1 = (p_ + r_)/1000; s3 = (p_ - r_)/1000          # actual Mohr state (MPa)
            Q = MB*s3/SCI + S; okQ = bool(np.all(Q > 0)); Qs = np.where(Q > 0, Q, np.nan)
            hbres = float(np.max(np.abs(s1 - s3 - SCI*Qs**a_)))                                                       # HB (sigma2' = sigma3') residual
            Mhb = 1 + a_*MB*Qs**(a_-1); tanres = float(np.max(np.abs(phi - np.arcsin((Mhb-1)/(Mhb+1)))))           # tangency on the Mohr state
            order_fail = int(np.sum(s1 < s3 - 1e-12))
            Sm = (c*A + snf*A*t)/Fn; eqres = float(np.max(np.abs(W - ((snf+u)*A*g + Sm*kap))/W))                    # column vertical equilibrium
            dF = abs(Fn-F)/abs(Fn)
            accepted = (inside and finite and okQ and int(proj.sum()) == 0 and cons <= TOL_ANG and tanres <= TOL_ANG and hbres <= TOL_MPA
                        and order_fail == 0 and eqres <= TOL_EQ and dF <= TOL_F)
            return dict(F=Fn, iters=it, proj_path=nproj, proj_end=int(proj.sum()), min_sn=float(snf.min()), n_tension=int((snf < 0).sum()),
                        inside_limit=inside, finite=finite, Q_positive=okQ, phi_consistency_rad=cons, tangency_residual_rad=tanres, hb_residual_MPa=hbres,
                        ordering_failures=order_fail, col_equilibrium_rel=eqres, dF_rel=dF, accepted=accepted)
        F, sn = Fn, new
    return dict(F=F, iters=-1, proj_path=nproj, proj_end=-1, min_sn=float(sn.min()), n_tension=-1, inside_limit=False, finite=False, Q_positive=False,
                phi_consistency_rad=float("nan"), tangency_residual_rad=float("nan"), hb_residual_MPa=float("nan"), ordering_failures=-1,
                col_equilibrium_rel=float("nan"), dF_rel=float("nan"), accepted=False)
CASES = [("Dry, a = 0.500000", "AppendixE", 0.5, False), ("Wet, a = 0.500000", "AppendixE", 0.5, True),
         ("Dry, a = 0.509923", "AppendixE", 0.509923, False), ("Wet, a = 0.509923", "AppendixE", 0.509923, True),
         ("Dry, a = 0.500000, published circle", "Published", 0.5, False)]
MESH = {3: (10, 15, 20, 30, 40, 60, 80, 120, 160, 240, 320), 2: (20, 40, 80, 160, 320, 640, 1280, 2560)}
if __name__ == "__main__":
    rows = []
    for lab, surf, a, wet in CASES:
        tab = hb_table(a)
        for dim in (3, 2):
            prevF = None
            for n in MESH[dim]:
                e = elements(dim, n, surf, wet); r = solve(e, tab)
                dpct = None if prevF is None else abs(r["F"]-prevF)/r["F"]*100; prevF = r["F"]
                rows.append(dict(case=lab, surface=surf, a=a, wet=wet, dim=f"{dim}D", mesh=(f"{n}x{n}" if dim == 3 else f"{n} slices"), elements=len(e["W"]),
                                 volume_m3=e["vol"]*(2 if dim == 3 else 1), change_from_previous_pct=dpct, **r))
                print(f"{lab:38s} {dim}D {n:5d}: F={r['F']:.6f} accepted={r['accepted']} min sn={r['min_sn']:.2f} proj(path/end)={r['proj_path']}/{r['proj_end']}", flush=True)
    with open("mesh_study_records.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    # reproduction of the workbook at its own meshes (3D 20x20, 2D 40 slices)
    import json, os
    wbvals = json.load(open("workbook_values.json"))     # recorded coarse workbook results (3D 20x20, 2D 40 slices)
    with open("reproduction_check.txt", "w") as f:
        for key, v in wbvals.items():
            lab, dim = key.rsplit("|", 1)
            r = [x for x in rows if x["case"] == lab and x["dim"] == dim and x["mesh"] in ("20x20", "40 slices")][0]
            f.write(f"{lab} {dim} workbook mesh: reimplementation {r['F']:.12f} | workbook {v:.12f} | diff {abs(r['F']-v):.1e}\n")
    print(open("reproduction_check.txt").read())
