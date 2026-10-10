import math
sci, mi, GSI, D = 77.7, 10, 42, 0
mb = mi*math.exp((GSI-100)/(28-14*D)); s = math.exp((GSI-100)/(9-3*D))
def state(s3, nu, a):
    Q = mb*s3/sci + s; k = 1/(2*nu) - 1; rad = sci**2*Q**(2*a) - k*k*s3*s3
    if rad < 0: return None
    s1 = s3/(2*nu) + math.sqrt(rad)/math.sqrt(3); s2 = s3/nu - s1
    if not (s1 >= s2 >= s3): return None
    M = (3*s1/nu - 2*(1-1/nu+1/nu**2)*s3 + 2*a*mb*sci*Q**(2*a-1))/(6*s1 - 3*s3/nu)
    if M < 1: return None
    p = (M-1)/(M+1); return (s1+s3)/2 - (s1-s3)/2*p, math.degrees(math.asin(p))
def min_sn(nu, a, lo=1e-6, hi=1e6, N=200000):
    best = None; L0, L1 = math.log10(lo), math.log10(hi)
    for i in range(N + 1):
        s3 = 10**(L0 + (L1-L0)*i/N); st = state(s3, nu, a)
        if st and (best is None or st[0] < best[0]): best = (st[0], st[1], s3)
    if best is None: return None
    a_, b_ = best[2]/1.01, best[2]*1.01                       # refine (golden section on admissible interval)
    g = (math.sqrt(5)-1)/2
    f = lambda x: (state(x, nu, a) or (1e300,))[0]
    for _ in range(200):
        c1 = b_-g*(b_-a_); c2 = a_+g*(b_-a_)
        if f(c1) < f(c2): b_ = c2
        else: a_ = c1
    x = (a_+b_)/2; st = state(x, nu, a); return (st[0], st[1], x) if st and st[0] <= best[0] else best
print("search: sigma3' from 1e-6 to 1e6 MPa, 200001 log-spaced points + golden-section refinement")
for nu in (0.05, 0.10, 0.15, 0.20, 0.25, 0.271, 0.30, 0.35, 0.40, 0.45):
    r = [min_sn(nu, a) for a in (0.5, 0.509923)]
    f = lambda b: "none in range" if b is None else f"{b[0]:12.6f} MPa (phi {b[1]:6.3f} deg, sigma3' {b[2]:.4g})"
    print(f"nu={nu:5.3f} | a=0.5: {f(r[0])} | a=0.509923: {f(r[1])}")
