"""Reference implementation of EJC_Eq421.bas v6 solver/state logic and 22 shared fixtures.
Python additionally rejects nonfinite acceptance tolerances and Boolean type misuse.
The VBA was statically checked for case-insensitive duplicate declarations, not executed here."""
import math
HI_DEG, NSCAN = 89.9999, 1800
SCAN_OFFSET_DEG = 0.001  # avoid an artificial sign change from roundoff at the tensile apex
def phi_b(sn, sci, mb, s, a):
    """angle above which the eq. 4.21 state is inside the HB tensile limit (Q > 0): closed form (a = 0.5) + Newton"""
    T = s*sci/mb; K = sn + T
    if K <= 0: return None
    Z = 2*K/(sci*(mb*K/sci)**a)
    if Z >= 1: return 0.0
    p = (1+Z*Z) - math.sqrt((1+Z*Z)**2 - (1-Z*Z))
    for _ in range(8):
        f = (1-p)*(1+p/a)**(-a) - Z; fp = -(1+p/a)**(-a)*(1 + (1-p)/(1+p/a)); p = min(max(p - f/fp, 0.0), 1-1e-15)
    return math.asin(p)
def valid_inputs(sci, mb, s, a, tol):
    return all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in (sci,mb,s,a,tol)) and sci > 0 and mb > 0 and s >= 0 and 0 < a <= 1 and 0 < tol <= 1e-6
def gstate(phi, sn, sci, mb, s, a, pC, qC):
    """returns (g, dg) of g = sin(phi) - (M-1)/(M+1) = sin(phi) - 1 + 2D/(N+D), or None outside the domain"""
    P, R, S = pC*pC-pC+1, 2*pC*qC-pC-qC-1, qC*qC-qC+1; Sx, Cx = math.sin(phi), math.cos(phi); B = mb*sn/sci + s
    if B <= 0: return None
    H = (B/(1+Sx/a))**a; dH = -H*Cx/(1+Sx/a)
    s3 = sn - sci/2*H*(1-Sx); ds3 = -sci/2*(dH*(1-Sx) - H*Cx); s1 = sn + sci/2*H*(1+Sx); ds1 = sci/2*(dH*(1+Sx) + H*Cx)
    hbQ = mb*s3/sci + s; dQ = mb*ds3/sci
    if hbQ <= 0: return None
    N = 2*a*mb*sci*hbQ**(2*a-1) - R*s1 - 2*S*s3; dN = 2*a*mb*sci*(2*a-1)*hbQ**(2*a-2)*dQ - R*ds1 - 2*S*ds3
    D = 2*P*s1 + R*s3; dD = 2*P*ds1 + R*ds3
    if N + D == 0 or D == 0: return None
    return Sx - 1 + 2*D/(N+D), Cx - 2*(dN*D - N*dD)/(N+D)**2, N/D
def ang_res(phi, sn, sci, mb, s, a, pC, qC):
    """unclipped eq. 4.19 residual phi - asin((M-1)/(M+1)); None if M+1 = 0 or the argument is outside [-1, 1]"""
    r = gstate(phi, sn, sci, mb, s, a, pC, qC)
    if r is None or r[2] + 1 == 0: return None
    z = (r[2]-1)/(r[2]+1)
    return None if abs(z) > 1 else phi - math.asin(z)
def brackets(sn, sci, mb, s, a, pC, qC):
    """candidate brackets with sign change between defined points, poles removed (|g| grows on refinement)"""
    out = []; prev = None; pb = phi_b(sn, sci, mb, s, a)
    if pb is None: return out
    x0 = max(pb + math.radians(SCAN_OFFSET_DEG), math.radians(0.01)); x1 = math.radians(HI_DEG)
    if x0 >= x1: return out
    for k in range(NSCAN + 1):
        x = x0 + (x1 - x0)*k/NSCAN; r = gstate(x, sn, sci, mb, s, a, pC, qC)
        if r and prev and prev[1]*r[0] <= 0:
            lo, hi, glo, ghi = prev[0], x, prev[1], r[0]; pole = False
            for _ in range(30):                                     # bisect the sign change; a pole makes |g| grow
                m = (lo+hi)/2; rm = gstate(m, sn, sci, mb, s, a, pC, qC)
                if rm is None: pole = True; break
                if glo*rm[0] <= 0: hi, ghi = m, rm[0]
                else: lo, glo = m, rm[0]
            if not pole and min(abs(glo), abs(ghi)) > 1e-3: pole = True
            if not pole: out.append((prev[0], x))
        prev = (x, r[0]) if r else None
    return out
def ejc_phi(sn_kPa, sci, mb, s, a, pC, qC, tol=1e-10, root_index=1):
    if not valid_inputs(sci, mb, s, a, tol) or not all(map(math.isfinite, (sn_kPa, pC, qC))) or type(root_index) is not int or root_index < 1: return "#VALUE!"
    sn = sn_kPa/1000; br = brackets(sn, sci, mb, s, a, pC, qC)
    if len(br) < root_index: return "#N/A"
    lo, hi = br[root_index-1]; glo = gstate(lo, sn, sci, mb, s, a, pC, qC)[0]; x = (lo+hi)/2
    for it in range(200):
        r = gstate(x, sn, sci, mb, s, a, pC, qC)
        if r is None: return "#NUM!"
        g, dg, _ = r
        if glo*g <= 0: hi = x
        else: lo, glo = x, g
        xn = x - g/dg if dg != 0 else (lo+hi)/2
        if xn <= lo or xn >= hi: xn = (lo+hi)/2
        step = abs(xn - x); x = xn; rr = ang_res(x, sn, sci, mb, s, a, pC, qC)
        if rr is not None and abs(rr) <= tol and step <= tol: return math.degrees(x)
    return "#NUM!"                                                   # exhausted: no converged root
DOMAIN_FAIL = lambda pC, qC: dict(angular=None, ejc=None, closure=None, ordering=0, s1=None, s3=None, domain=0, branch=0)
def ejc_check(sn_kPa, phi_deg, sci, mb, s, a, pC, qC, tol=1e-10):
    """checks of the eq. 4.21 (C1) candidate state at (sigma_n', phi'i); controlled domain failure, never an exception"""
    if not valid_inputs(sci, mb, s, a, tol) or not all(map(math.isfinite, (sn_kPa, phi_deg, pC, qC))): return DOMAIN_FAIL(pC, qC)
    if not (0 <= phi_deg < 90): return DOMAIN_FAIL(pC, qC)                       # adopted branch: 0 <= phi < 90 deg
    sn = sn_kPa/1000; phi = math.radians(phi_deg); P, R, S = pC*pC-pC+1, 2*pC*qC-pC-qC-1, qC*qC-qC+1
    base = mb*sn/sci + s; den = 1 + math.sin(phi)/a
    if base <= 0 or den <= 0: return DOMAIN_FAIL(pC, qC)                         # guard before the fractional power
    H = (base/den)**a; s3 = sn - sci/2*H*(1-math.sin(phi)); s1 = sn + sci/2*H*(1+math.sin(phi)); s2 = pC*s1 + qC*s3
    return state_checks(phi, s1, s2, s3, sci, mb, s, a, P, R, S, pC, qC, tol)
def ejc_check_state(sn_kPa, c_kPa, phi_deg, sci, mb, s, a, pC, qC, tol=1e-10):
    """checks of the FORCE-USED state (sigma_n', c'i, phi'i): Mohr circle tau = c + sigma_n tan(phi)"""
    if not valid_inputs(sci, mb, s, a, tol) or not all(map(math.isfinite, (sn_kPa, c_kPa, phi_deg, pC, qC))): return DOMAIN_FAIL(pC, qC)
    if not (0 <= phi_deg < 90): return DOMAIN_FAIL(pC, qC)
    phi = math.radians(phi_deg); tau = c_kPa + sn_kPa*math.tan(phi); p_ = sn_kPa + tau*math.tan(phi); r_ = tau/math.cos(phi)
    s1, s3 = (p_ + r_)/1000, (p_ - r_)/1000; P, R, S = pC*pC-pC+1, 2*pC*qC-pC-qC-1, qC*qC-qC+1
    return state_checks(phi, s1, pC*s1 + qC*s3, s3, sci, mb, s, a, P, R, S, pC, qC, tol)
def state_checks(phi, s1, s2, s3, sci, mb, s, a, P, R, S, pC, qC, tol=1e-10):
    hbQ = mb*s3/sci + s; out = dict(angular=None, ejc=None, closure=s2-(pC*s1+qC*s3), ordering=int(s1 >= s2-1e-9 and s2 >= s3-1e-9), s1=s1, s3=s3, domain=0, branch=0)
    if hbQ <= 0: return out
    den = 2*P*s1 + R*s3
    if den == 0: return out
    M = (2*a*mb*sci*hbQ**(2*a-1) - R*s1 - 2*S*s3)/den
    if M + 1 == 0 or abs((M-1)/(M+1)) > 1: return out
    rad = 4*sci**2*hbQ**(2*a) - 3*(s2-s3)**2
    if rad < 0: return out
    z = (M-1)/(M+1)
    out.update(angular=phi - math.asin(z), ejc=s1 - (s2+s3+math.sqrt(rad))/2, domain=1,
               branch=int(phi >= 0 and phi < math.pi/2 and math.asin(z) >= -tol))      # adopted branch: tangent angle >= 0 (M >= 1) within tol
    return out
def ejc_accept(chk, tol_ang=1e-10, tol_ejc=1e-8):
    if not all(isinstance(t, (int, float)) and not isinstance(t, bool) and math.isfinite(t) and t > 0 for t in (tol_ang, tol_ejc)) or tol_ang > 1e-6 or tol_ejc > 1e-4:
        return "REJECTED: INPUT"
    if not isinstance(chk, dict) or chk.get("domain") != 1: return "REJECTED: DOMAIN"
    if not all(isinstance(chk.get(k), (int, float)) and not isinstance(chk[k], bool) and math.isfinite(chk[k]) for k in ("angular", "ejc", "closure", "s1", "s3")):
        return "REJECTED: DOMAIN"
    if chk.get("branch") != 1: return "REJECTED: BRANCH"
    if abs(chk["angular"]) > tol_ang: return "REJECTED: TANGENCY"
    if abs(chk["ejc"]) > tol_ejc: return "REJECTED: EJC"
    if abs(chk["closure"]) > tol_ejc: return "REJECTED: CLOSURE"
    if chk.get("ordering") != 1: return "REJECTED: ORDERING"
    return "ADMISSIBLE"
if __name__ == "__main__":
    mbC, sC = 10*math.exp((42-100)/28), math.exp((42-100)/9); V = (120, 3.5, 0.002, 0.505)
    tests = []
    ph = ejc_phi(400, 77.7, mbC, sC, 0.5, 0, 1); tests.append(("T1 closure 3, Chile, 400 kPa", ph, isinstance(ph, float) and abs(ph-58.448348604037) < 1e-9 and ejc_accept(ejc_check(400, ph, 77.7, mbC, sC, 0.5, 0, 1)) == "ADMISSIBLE"))
    ph = ejc_phi(5740.3, *V, -1, 20); a_ = ejc_accept(ejc_check(5740.3, ph, *V, -1, 20)); tests.append(("T2 closure 1, verification point (first root)", f"{ph} -> {a_}", abs(ph-53.3178012) < 1e-6 and a_ == "REJECTED: EJC"))
    ph = ejc_phi(-98.00411947453206, 77.7, mbC, sC, 0.5, 0, 1); rr = None if not isinstance(ph, float) else ang_res(math.radians(ph), -0.09800411947453206, 77.7, mbC, sC, 0.5, 0, 1)
    tests.append(("T3 near tensile limit (-98.0041 kPa): converged with |angular residual| <= 1e-10, or a controlled failure", f"{ph} (residual {rr})", (isinstance(ph, str)) or (rr is not None and abs(rr) <= 1e-10)))
    ph = ejc_phi(-98.0041250683277, 77.7, mbC, sC, 0.5, 0, 1); tests.append(("T4 valid HB point with phi near 89.95 deg", ph, isinstance(ph, float) and abs(ph-89.95) < 0.01 and ejc_accept(ejc_check(-98.0041250683277, ph, 77.7, mbC, sC, 0.5, 0, 1)) == "ADMISSIBLE"))
    ph = ejc_phi(400, 77.7, mbC, sC, 0.5, -1, 1/0.271); tests.append(("T5 closure 1, Chile, nu = 0.271: no root", ph, ph == "#N/A"))
    ph = ejc_phi(400, 0, mbC, sC, 0.5, 0, 1); tests.append(("T6 invalid input sigma_ci = 0", ph, ph == "#VALUE!"))
    nb = len(brackets(13.26501174607056, *V, -1, 20)); found = [ejc_phi(13265.01174607056, *V, -1, 20, root_index=i) for i in range(1, nb+1)]
    hit = [f for f in found if isinstance(f, float) and abs(f-39.89156042767222) < 1e-6]
    tests.append(("T7 isolated eq. 4.21 counterexample (13.265 MPa): an ADMISSIBLE root exists", f"brackets {nb}, roots {found}", bool(hit) and ejc_accept(ejc_check(13265.01174607056, hit[0], *V, -1, 20)) == "ADMISSIBLE"))
    ph = ejc_phi(400, 77.7, mbC, sC, 0.5, 0, 1); chk = ejc_check(400, ph, 77.7, mbC, sC, 0.5, 0, 1); H_ = math.radians(ph)
    c_kPa = 1000*((chk["s1"]-chk["s3"])/2*math.cos(H_) - 0.4*math.tan(H_)); tests.append(("T8 force-used state check (sigma_n, c, phi) of T1", ejc_accept(ejc_check_state(400, c_kPa, ph, 77.7, mbC, sC, 0.5, 0, 1)), ejc_accept(ejc_check_state(400, c_kPa, ph, 77.7, mbC, sC, 0.5, 0, 1)) == "ADMISSIBLE"))
    tests.append(("T9 perturbed force-used state (c + 5 kPa) must be rejected", ejc_accept(ejc_check_state(400, c_kPa+5, ph, 77.7, mbC, sC, 0.5, 0, 1)), ejc_accept(ejc_check_state(400, c_kPa+5, ph, 77.7, mbC, sC, 0.5, 0, 1)) != "ADMISSIBLE"))
    T7 = (13265.01174607056, *V, -1, 20)
    tests.append(("T10 root index 0 -> invalid", ejc_phi(*T7, root_index=0), ejc_phi(*T7, root_index=0) == "#VALUE!"))
    tests.append(("T11 root index -1 -> invalid", ejc_phi(*T7, root_index=-1), ejc_phi(*T7, root_index=-1) == "#VALUE!"))
    p2 = ejc_phi(*T7, root_index=2); l2 = ejc_accept(ejc_check(13265.01174607056, p2, *V, -1, 20)) if isinstance(p2, float) else None
    tests.append(("T12 root index 2 at T7 inputs", f"{p2} -> {l2}", isinstance(p2, float) and abs(p2-60.4228863398) < 1e-6 and l2 is not None and l2 != "ADMISSIBLE"))
    tests.append(("T13 tolerance 0 -> invalid", ejc_phi(*T7, tol=0), ejc_phi(*T7, tol=0) == "#VALUE!"))
    tests.append(("T14 tolerance negative -> invalid", ejc_phi(*T7, tol=-1e-10), ejc_phi(*T7, tol=-1e-10) == "#VALUE!"))
    nb = ejc_accept(ejc_check_state(29819.0191175968, 25010.6868294027, -0.012480264998, *V, -1, 20))
    tests.append(("T15 negative-friction branch state (QA example) -> outside the adopted branch", nb, nb in ("REJECTED: DOMAIN", "REJECTED: BRANCH")))
    m45 = ejc_accept(ejc_check(13265.01174607056, -45, *V, -1, 20)); tests.append(("T16 candidate at -45 deg -> controlled invalid domain", m45, m45 == "REJECTED: DOMAIN"))
    p90 = ejc_accept(ejc_check(13265.01174607056, 90, *V, -1, 20)); tests.append(("T17 candidate at 90 deg -> rejected (finite-angle policy)", p90, p90 == "REJECTED: DOMAIN"))
    ck = ejc_check_state(400, c_kPa, ph, 77.7, mbC, sC, 0.5, 0, 1)
    for number, label, kwargs in [
        (18, "acceptance EJC tolerance zero", {"tol_ejc":0}),
        (19, "acceptance EJC tolerance negative", {"tol_ejc":-1e-8}),
        (20, "acceptance angular tolerance zero", {"tol_ang":0}),
        (21, "acceptance angular tolerance negative", {"tol_ang":-1e-10}),
    ]:
        result = ejc_accept(ck, **kwargs)
        tests.append((f"T{number} {label}", result, result == "REJECTED: INPUT"))
    ph = ejc_phi(-49.768342, 77.7, mbC, sC, 0.509923, 0, 1)
    apex_ok = isinstance(ph, float) and abs(ph-71.27787820469955) < 1e-7
    if apex_ok:
        apex_ok = ejc_accept(ejc_check(-49.768342, ph, 77.7, mbC, sC, 0.509923, 0, 1)) == "ADMISSIBLE"
    tests.append(("T22 wet a=.509923 tensile-apex scan regression", ph, apex_ok))
    for name, val, ok in tests: print(("PASS " if ok else "FAIL ") + name + " | " + str(val)[:110])
    print("ALL PASS" if all(t[2] for t in tests) else "SOME FAILED")
    raise SystemExit(0 if all(t[2] for t in tests) else 1)
