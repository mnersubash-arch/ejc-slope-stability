"""Calculation engine: 3D Hungr columns and 2D Bishop slices with the Extended Jiang Criterion (EJC).

No Streamlit imports: the module is tested on its own (``pytest`` in this folder).

Conventions: compression positive, effective stresses. Strength is evaluated in MPa; column forces in kN, stresses
in kPa, areas in m^2 (1 MPa = 1000 kPa). The toe is at the origin and the slope rises towards negative x.
Equation numbers (``eq. 4.12`` and so on) are those of the submitted thesis; ``brief §n`` is the app brief.

Integrity rules: a local state that fails any acceptance test is REJECTED with its reason. GHB is never
substituted for the EJC, and a factor of safety is reported only for a run that passes every check.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, fields, replace
from functools import lru_cache
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------------------------------------------
# constants and tolerances (brief §1.3, §1.4)
# ---------------------------------------------------------------------------------------------------------------
GAMMA_W = 9.81             # kN/m^3
TOL_EJC = 1e-8             # |f| of the unsquared EJC [MPa]
TOL_CLOSURE = 1e-8         # |sigma2 - closure(sigma1, sigma3)| [MPa]
TOL_DF = 1e-12             # |dF|/F
TOL_DSN = 1e-9             # max |d sigma_n| [kPa]
TOL_LOCAL = 1e-9           # relative local (vertical) equilibrium residual
TOL_GLOBAL = 1e-8          # relative global equilibrium residual
MAX_ITER = 2000
N_SCAN = 4000              # sigma3 scan points (geometric in sigma3 - sigma3_t)
SCAN_DECADES = 14          # the scan starts 1e-14 of the span above the tensile limit
N_BISECT = 200             # bisection cap; it stops earlier once the bracket cannot shrink

MODELS = {
    "ejc_moscarda": "EJC, Moscarda closure sigma2 = nu (sigma1 + sigma3) (adopted)",
    "ejc_original": "EJC, original closure sigma2 = sigma3/nu - sigma1 (admissibility check only)",
    "ghb": "Conventional generalised Hoek-Brown (sigma2 = sigma3)",
    "ejc_sigma2_eq_sigma3": "EJC with the closure sigma2 = sigma3 (must reproduce GHB)",
    "mc": "Mohr-Coulomb, constant c' and phi' (D'Esterre presets)",
}

# local-state reason codes (brief §1.3 acceptance); 0 = accepted
OK, DOMAIN, EJC_RES, CLOSURE_RES, ORDERING, GRADIENT, SHEAR, NO_STATE, BISECTION = range(9)
REASONS = {
    OK: "ACCEPTED",
    DOMAIN: "REJECTED: DOMAIN (Q <= 0, complex root or non-finite state)",
    EJC_RES: "REJECTED: EJC residual |f| >= 1e-8 MPa",
    CLOSURE_RES: "REJECTED: closure residual >= 1e-8 MPa",
    ORDERING: "REJECTED: ORDERING (sigma1 >= sigma2 >= sigma3 fails)",
    GRADIENT: "REJECTED: GRADIENT (M < 1 or not finite)",
    SHEAR: "REJECTED: SHEAR (tau_s < 0)",
    NO_STATE: "REJECTED: no admissible state at this sigma_n (outside the admissible range of the closure)",
    BISECTION: "REJECTED: bisection left the admissible interval",
}


# ---------------------------------------------------------------------------------------------------------------
# eq. 2.5 rock-mass parameters (Hoek et al. 2002)
# ---------------------------------------------------------------------------------------------------------------
def hoek_brown_parameters(gsi: float, mi: float, d: float) -> tuple[float, float, float]:
    """Eq. 2.5: return (m_b, s, a) from GSI, m_i and D."""
    mb = mi * math.exp((gsi - 100.0) / (28.0 - 14.0 * d))
    s = math.exp((gsi - 100.0) / (9.0 - 3.0 * d))
    a = 0.5 + (math.exp(-gsi / 15.0) - math.exp(-20.0 / 3.0)) / 6.0
    return mb, s, a


@dataclass(frozen=True)
class Strength:
    """Strength model and constants. ``sigma_ci`` in MPa; Mohr-Coulomb ``c_mc`` in kPa, ``phi_mc`` in degrees."""

    sigma_ci: float
    mb: float
    s: float
    a: float
    model: str = "ejc_moscarda"
    nu: float = 0.271
    c_mc: float = 0.0
    phi_mc: float = 0.0

    def __post_init__(self):
        if self.model not in MODELS:
            raise ValueError(f"unknown strength model {self.model!r}")
        if self.model == "mc":
            if not (self.c_mc >= 0 and 0 <= self.phi_mc < 90):
                raise ValueError("require c' >= 0 and 0 <= phi' < 90 deg")
            return
        if not (self.sigma_ci > 0 and self.mb > 0 and self.s > 0 and 0 < self.a <= 1):
            raise ValueError("require sigma_ci > 0, m_b > 0, s > 0 and 0 < a <= 1")
        if self.model in ("ejc_moscarda", "ejc_original") and not 0 < self.nu < 0.5:
            raise ValueError("require 0 < nu < 0.5")

    @property
    def pq(self) -> tuple[float, float]:
        """(p, q) of the closure sigma2 = p sigma1 + q sigma3."""
        if self.model == "ejc_moscarda":
            return self.nu, self.nu
        if self.model == "ejc_original":
            return -1.0, 1.0 / self.nu
        return 0.0, 1.0                      # ghb and ejc_sigma2_eq_sigma3

    @property
    def sigma3_tensile(self) -> float:
        """sigma3 at Q = 0, the Hoek-Brown tensile limit [MPa]."""
        return -self.s * self.sigma_ci / self.mb


def closure_constants(p: float, q: float) -> tuple[float, float, float]:
    """Squared EJC on sigma2 = p sigma1 + q sigma3:  P sigma1^2 + R sigma1 sigma3 + S sigma3^2 = sigma_ci^2 Q^(2a).

    P = p^2 - p + 1, R = 2pq - p - q - 1, S = q^2 - q + 1. Moscarda closure (eq. 4.3, p = q = nu): P = S = C and
    R = -B (eq. 4.4); original closure (eq. 4.1, p = -1, q = 1/nu): P = 3, R = -3/nu, S = 1 - 1/nu + 1/nu^2 (eq. 4.2).
    """
    return p * p - p + 1.0, 2.0 * p * q - p - q - 1.0, q * q - q + 1.0


# ---------------------------------------------------------------------------------------------------------------
# principal states along the closure (eqs 2.7, 2.8, 4.1-4.5) and tangent parameters (eqs 2.6, 4.6, 4.7)
# ---------------------------------------------------------------------------------------------------------------
@dataclass
class States:
    """Principal states at given sigma3 (MPa), vectorised. ``code`` uses the reason codes above."""

    s1: np.ndarray
    s2: np.ndarray
    s3: np.ndarray
    Q: np.ndarray
    M: np.ndarray
    f: np.ndarray              # unsquared EJC residual [MPa]
    closure: np.ndarray        # sigma2 - (p sigma1 + q sigma3) [MPa]
    sn: np.ndarray             # sigma_n of the tangent point [MPa]
    tau: np.ndarray            # tau_s [MPa]
    phi: np.ndarray            # phi_i [rad]
    c: np.ndarray              # c_i [MPa]
    code: np.ndarray


def ejc_residual(s1, s2, s3, Q, sigma_ci, a):
    """Eq. 2.8 (unsquared EJC): f = 2 sigma1 - sigma2 - sigma3 - sqrt(4 sigma_ci^2 Q^(2a) - 3 (sigma2 - sigma3)^2); NaN if not real."""
    with np.errstate(invalid="ignore"):
        rad = 4.0 * sigma_ci ** 2 * np.power(Q, 2.0 * a) - 3.0 * (s2 - s3) ** 2
        return np.where(rad >= 0, 2.0 * s1 - s2 - s3 - np.sqrt(np.where(rad >= 0, rad, 0.0)), np.nan)


def principal_states(s3, st: Strength) -> States:
    """Principal state, gradient and tangent parameters on the selected model at each sigma3 [MPa].

    EJC models: sigma1 is the root of the squared form (eq. 4.2 or 4.4) that satisfies the unsquared EJC (eq. 2.8),
    Q > 0, a real radicand and sigma1 >= sigma2 >= sigma3 (the + root is tried first); M is the gradient along the
    closure (eq. 4.5). GHB (eq. 2.7): sigma1 = sigma3 + sigma_ci Q^a, M = 1 + a m_b Q^(a-1).
    Tangent parameters (eqs 2.6, 4.6, 4.7): sin(phi_i) = (M-1)/(M+1), sigma_n = sigma3 + (sigma1-sigma3)/(M+1),
    tau_s = (sigma1-sigma3) sqrt(M)/(M+1), c_i = tau_s - sigma_n tan(phi_i).
    """
    s3 = np.asarray(s3, dtype=float)
    sci, mb, s, a = st.sigma_ci, st.mb, st.s, st.a
    p, q = st.pq
    with np.errstate(all="ignore"):
        Q = mb * s3 / sci + s
        qpos = Q > 0
        Qp = np.where(qpos, Q, np.nan)
        if st.model == "ghb":
            s1 = s3 + sci * np.power(Qp, a)
            s2 = s3.copy()
            M = 1.0 + a * mb * np.power(Qp, a - 1.0)
            real = qpos
        else:
            cP, cR, cS = closure_constants(p, q)
            K = sci ** 2 * np.power(Qp, 2.0 * a)
            disc = cR * cR * s3 * s3 - 4.0 * cP * (cS * s3 * s3 - K)
            sq = np.sqrt(np.where(disc >= 0, disc, np.nan))

            def branch(sign):
                b1 = (-cR * s3 + sign * sq) / (2.0 * cP)
                b2 = p * b1 + q * s3
                bf = ejc_residual(b1, b2, s3, Qp, sci, a)
                good = qpos & (disc >= 0) & np.isfinite(bf) & (np.abs(bf) < TOL_EJC) & (b1 >= b2) & (b2 >= s3)
                return b1, good

            s1p, okp = branch(+1.0)
            s1m, okm = branch(-1.0)
            s1 = np.where(okp | ~okm, s1p, s1m)
            s2 = p * s1 + q * s3
            real = qpos & (disc >= 0)
            M = (2.0 * a * mb * sci * np.power(Qp, 2.0 * a - 1.0) - cR * s1 - 2.0 * cS * s3) / (2.0 * cP * s1 + cR * s3)
        f = ejc_residual(s1, s2, s3, Qp, sci, a)
        clo = s2 - (p * s1 + q * s3)
        sinp = (M - 1.0) / (M + 1.0)
        phi = np.arcsin(sinp)
        sn = s3 + (s1 - s3) / (M + 1.0)
        tau = (s1 - s3) * np.sqrt(M) / (M + 1.0)
        c = tau - sn * np.tan(phi)
    code = np.full(s3.shape, OK, dtype=int)
    real = real & np.isfinite(f) & np.isfinite(s1)
    code[~real] = DOMAIN
    code[(code == OK) & ~(np.abs(f) < TOL_EJC)] = EJC_RES
    code[(code == OK) & ~(np.abs(clo) < TOL_CLOSURE)] = CLOSURE_RES
    code[(code == OK) & ~((s1 >= s2) & (s2 >= s3))] = ORDERING
    code[(code == OK) & ~(np.isfinite(M) & (M >= 1.0))] = GRADIENT
    code[(code == OK) & ~(tau >= 0)] = SHEAR
    return States(s1, s2, s3, Q, M, f, clo, sn, tau, phi, c, code)


# ---------------------------------------------------------------------------------------------------------------
# Route 2 (brief §1.3): sigma_n target -> sigma3 by scan + bisection
# ---------------------------------------------------------------------------------------------------------------
@dataclass
class StrengthResult:
    """Local strength at target normal stresses. Stresses in kPa except the principal stresses (MPa)."""

    sn_target: np.ndarray      # kPa
    sn: np.ndarray             # kPa, sigma_n of the accepted state
    phi: np.ndarray            # rad
    c: np.ndarray              # kPa
    tau: np.ndarray            # kPa
    M: np.ndarray
    s1: np.ndarray             # MPa
    s2: np.ndarray
    s3: np.ndarray
    f: np.ndarray              # MPa
    closure: np.ndarray        # MPa
    code: np.ndarray
    n_brackets: np.ndarray     # sign changes seen by the scan (1 = unique on the scan grid)

    @property
    def accepted(self) -> np.ndarray:
        return self.code == OK

    def reason(self, k: int) -> str:
        return REASONS[int(self.code[k])]


def _monotone_pieces(ok: np.ndarray, v: np.ndarray) -> list[tuple[int, int]]:
    """Index ranges [i0, i1] of consecutive admissible scan points on which v is monotone."""
    pieces = []
    n = len(v)
    i = 0
    while i < n - 1:
        if not (ok[i] and ok[i + 1]):
            i += 1
            continue
        i0 = i
        direction = np.sign(v[i + 1] - v[i])
        while i < n - 1 and ok[i + 1] and np.sign(v[i + 1] - v[i]) in (direction, 0):
            i += 1
        pieces.append((i0, i))
    return pieces


def scan_sigma3(st: Strength, top: float, n_scan: int = N_SCAN) -> tuple[np.ndarray, States]:
    """Scan sigma3 over (sigma3_t, top] MPa: points geometric in sigma3 - sigma3_t (1e-14 of the span to the span),
    plus the ends of every admissible interval located by bisection on admissibility (so the first and last
    admissible states are resolved to machine precision, not to the scan spacing)."""
    s3t = st.sigma3_tensile
    grid = s3t + (top - s3t) * np.power(10.0, np.linspace(-SCAN_DECADES, 0.0, n_scan + 1))
    gs = principal_states(grid, st)
    ok = gs.code == OK
    k = np.flatnonzero(ok[:-1] != ok[1:])
    if len(k):
        a, b = grid[k].copy(), grid[k + 1].copy()          # ok[a] != ok[b]
        a_ok = ok[k].copy()
        for _ in range(80):
            mid = 0.5 * (a + b)
            if not np.any((mid > np.minimum(a, b)) & (mid < np.maximum(a, b))):
                break
            m_ok = principal_states(mid, st).code == OK
            same = m_ok == a_ok
            a = np.where(same, mid, a)
            b = np.where(same, b, mid)
        edge = np.where(a_ok, a, b)                         # the admissible end of each transition
        grid = np.unique(np.concatenate([grid, edge]))
        gs = principal_states(grid, st)
    return grid, gs


@lru_cache(maxsize=64)
def _cached_scan(st: Strength, span_exp: int, n_scan: int) -> tuple[np.ndarray, States]:
    """The scan for spans up to 2**span_exp MPa, reused across iterations (results differ only in the bracket)."""
    return scan_sigma3(st, st.sigma3_tensile + 2.0 ** span_exp, n_scan)


@lru_cache(maxsize=64)
def _cached_pieces(st: Strength, span_exp: int, n_scan: int) -> tuple[tuple[int, int], ...]:
    """The monotone pieces of that scan (``_monotone_pieces``), reused with it."""
    _, gs = _cached_scan(st, span_exp, n_scan)
    return tuple(_monotone_pieces(gs.code == OK, gs.sn))


def admissible_sn_min(st: Strength, top: Optional[float] = None) -> float:
    """Lowest sigma_n [kPa] of any admissible state found by the sigma3 scan (NaN if none)."""
    top = 10.0 * st.sigma_ci if top is None else top
    _, gs = scan_sigma3(st, top)
    ok = gs.code == OK
    return float(np.min(gs.sn[ok]) * 1000.0) if ok.any() else math.nan


def _mohr_coulomb(sn_kpa, st: Strength) -> StrengthResult:
    """Mohr-Coulomb with constant c' and phi': tau_s = c' + sigma_n tan(phi'). The principal stresses are those of the
    Mohr circle tangent to that line at sigma_n (centre sigma_n + tau tan(phi'), radius tau/cos(phi'), sigma2 = sigma3)
    and M = (1 + sin phi')/(1 - sin phi'). The EJC residuals do not apply (NaN). A state is accepted when tau_s >= 0."""
    t = np.atleast_1d(np.asarray(sn_kpa, dtype=float))
    n = len(t)
    phi = np.full(n, math.radians(st.phi_mc))
    c = np.full(n, float(st.c_mc))
    tan_p = np.tan(phi)
    tau = c + t * tan_p
    with np.errstate(invalid="ignore"):
        centre = t + tau * tan_p
        radius = tau / np.cos(phi)
    s3 = (centre - radius) / 1000.0
    code = np.where(~np.isfinite(tau), DOMAIN, np.where(tau >= 0, OK, SHEAR))
    return StrengthResult(sn_target=t, sn=t.copy(), phi=phi, c=c, tau=tau, M=(1 + np.sin(phi)) / (1 - np.sin(phi)),
                          s1=(centre + radius) / 1000.0, s2=s3.copy(), s3=s3, f=np.full(n, np.nan),
                          closure=np.full(n, np.nan), code=code, n_brackets=np.ones(n, dtype=int))


def solve_strength(sn_kpa, st: Strength, n_scan: int = N_SCAN) -> StrengthResult:
    """Route 2 (brief §1.3): for each target sigma_n, scan sigma3 for the admissible interval, bracket
    sigma_n(sigma3) - sigma_n,target and bisect. The first bracket (lowest sigma3) is used; ``n_brackets`` counts all.

    The scan runs over sigma3 in (sigma3_t, max target], since sigma3 <= sigma_n for every admissible state; the
    points are geometric in sigma3 - sigma3_t so the tensile end is resolved. One scan serves all targets.
    Mohr-Coulomb (``model="mc"``) needs no scan (see ``_mohr_coulomb``).
    """
    if st.model == "mc":
        return _mohr_coulomb(sn_kpa, st)
    target_kpa = np.atleast_1d(np.asarray(sn_kpa, dtype=float))
    t = target_kpa / 1000.0
    n = len(t)
    s3t = st.sigma3_tensile
    code = np.full(n, NO_STATE, dtype=int)
    nbr = np.zeros(n, dtype=int)
    lo = np.full(n, np.nan)
    hi = np.full(n, np.nan)
    glo = np.full(n, np.nan)
    cand = np.isfinite(t) & (t > s3t)
    if cand.any():
        span_exp = int(math.ceil(math.log2(float(t[cand].max()) - s3t)))
        grid, gs = _cached_scan(st, span_exp, n_scan)
        for i0, i1 in _cached_pieces(st, span_exp, n_scan):
            v = gs.sn[i0:i1 + 1]
            inc = v[-1] >= v[0]
            w = v if inc else v[::-1]
            inside = cand & (t >= w[0]) & (t <= w[-1])
            if not inside.any():
                continue
            k = np.clip(np.searchsorted(w, t[inside], side="left"), 1, len(w) - 1)
            if inc:
                ja, jb = i0 + k - 1, i0 + k
            else:
                ja, jb = i1 - k, i1 - k + 1
            idx = np.flatnonzero(inside)
            nbr[idx] += 1
            first = idx[np.isnan(lo[idx])]
            sel = np.isnan(lo[idx])
            lo[first] = grid[ja[sel]]
            hi[first] = grid[jb[sel]]
            glo[first] = gs.sn[ja[sel]] - t[first]
    have = np.isfinite(lo)
    code[have] = OK
    # bisection on sigma3 (vectorised over all targets that have a bracket)
    for _ in range(N_BISECT):
        act = have & (code == OK)
        if not act.any():
            break
        mid = 0.5 * (lo + hi)
        moving = act & (mid > lo) & (mid < hi)
        if not moving.any():
            break
        ms = principal_states(mid[moving], st)
        idx = np.flatnonzero(moving)
        bad = ms.code != OK
        code[idx[bad]] = BISECTION
        good = idx[~bad]
        gm = ms.sn[~bad] - t[good]
        left = glo[good] * gm <= 0
        hi[good[left]] = mid[good[left]]
        lo[good[~left]] = mid[good[~left]]
        glo[good[~left]] = gm[~left]
    s3 = np.where(have, lo, np.nan)
    use_hi = have & (code == OK)
    if use_hi.any():                                   # take the bracket end closer to the target
        e_lo = principal_states(lo[use_hi], st)
        e_hi = principal_states(hi[use_hi], st)
        pick_hi = np.abs(e_hi.sn - t[use_hi]) < np.abs(e_lo.sn - t[use_hi])
        s3[use_hi] = np.where(pick_hi, hi[use_hi], lo[use_hi])
    fs = principal_states(np.where(have, s3, s3t), st)
    final_code = np.where(code == OK, fs.code, code)
    nanf = lambda arr: np.where(have, arr, np.nan)       # noqa: E731
    return StrengthResult(
        sn_target=target_kpa, sn=nanf(fs.sn * 1000.0), phi=nanf(fs.phi), c=nanf(fs.c * 1000.0),
        tau=nanf(fs.tau * 1000.0), M=nanf(fs.M), s1=nanf(fs.s1), s2=nanf(fs.s2), s3=nanf(fs.s3),
        f=nanf(fs.f), closure=nanf(fs.closure), code=final_code, n_brackets=nbr)


# ---------------------------------------------------------------------------------------------------------------
# original closure (eqs 4.1, 4.2): minimum admissible sigma_n (admissibility check only, never a FoS; Table 5.2)
# ---------------------------------------------------------------------------------------------------------------
@dataclass
class OriginalClosureMinimum:
    nu: float
    a: float
    sigma3_min: float          # MPa, sigma3 at the minimum admissible sigma_n
    sn_min: float              # MPa
    phi_deg: float             # phi_i at that state
    found: bool


def original_closure_minimum(sigma_ci: float, mb: float, s: float, a: float, nu: float,
                             lo: float = 1e-6, hi: float = 1e6, n: int = 200_000) -> OriginalClosureMinimum:
    """Minimum sigma_n over all admissible states of the EJC on sigma2 = sigma3/nu - sigma1.

    Admissible: sigma3/(2 nu) <= sigma1 <= sigma3 (1 - nu)/nu (ordering), sigma1 - sigma3 >= sigma_ci Q^a, |f| < 1e-8 MPa,
    real state, M >= 1. sigma3 is searched on a log grid in [lo, hi] MPa (ordering needs sigma3 >= 0), then refined by a
    golden-section search around the best grid point.
    """
    st = Strength(sigma_ci, mb, s, a, "ejc_original", nu)

    def evaluate(x):
        ps = principal_states(np.atleast_1d(np.asarray(x, dtype=float)), st)
        ok = (ps.code == OK) & (ps.s1 - ps.s3 >= sigma_ci * np.power(np.where(ps.Q > 0, ps.Q, np.nan), a))
        ok &= ps.s1 >= ps.s3 / (2 * nu)
        ok &= ps.s1 <= ps.s3 * (1 - nu) / nu
        return ps, ok

    grid = np.logspace(math.log10(lo), math.log10(hi), n + 1)
    ps, ok = evaluate(grid)
    if not ok.any():
        return OriginalClosureMinimum(nu, a, math.nan, math.nan, math.nan, False)
    vals = np.where(ok, ps.sn, np.inf)
    k = int(np.argmin(vals))
    best = (float(ps.sn[k]), float(grid[k]), float(ps.phi[k]))
    a_, b_ = grid[max(k - 1, 0)], grid[min(k + 1, n)]
    gr = (math.sqrt(5.0) - 1.0) / 2.0

    def fval(x):
        p1, o1 = evaluate(x)
        return float(p1.sn[0]) if o1[0] else math.inf

    for _ in range(200):
        c1 = b_ - gr * (b_ - a_)
        c2 = a_ + gr * (b_ - a_)
        if fval(c1) < fval(c2):
            b_ = c2
        else:
            a_ = c1
    # The minimum can lie on the edge of the admissible set (sigma2 = sigma3), where the bracket closes on the edge and
    # its midpoint can fall just outside it: keep the best admissible point among the midpoint and both bracket ends.
    for x in (0.5 * (a_ + b_), a_, b_):
        p1, o1 = evaluate(x)
        if o1[0] and float(p1.sn[0]) <= best[0]:
            best = (float(p1.sn[0]), float(x), float(p1.phi[0]))
    return OriginalClosureMinimum(nu, a, best[1], best[0], math.degrees(best[2]), True)


# ---------------------------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Inputs:
    """Every input of one analysis. Lengths m, angles deg, gamma kN/m^3, sigma_ci MPa."""

    H: float = 360.0
    beta: float = 50.0
    xc: float = 355.7423115
    zc: float = 644.4192949
    R: float = 725.2669975
    gamma: float = 27.0
    sigma_ci: float = 77.7
    mi: float = 10.0
    gsi: float = 42.0
    D: float = 0.0
    a_override: Optional[float] = 0.5      # None = a from GSI (eq. 2.5)
    nu: float = 0.271
    model: str = "ejc_moscarda"
    wet: bool = False
    Hw: float = 360.0
    beta_w: float = 50.0
    n3: int = 20
    n2: int = 40
    form: str = "FH"                        # "FH" (governing) or "FM" (moment form, 3D only)
    tension: str = "none"                   # "none", "C0" or "T0" (brief §1.6 sensitivities, Table 5.7)
    iteration: str = "default"              # "default" (brief §1.4) or "native_v15" (see solve())
    safeguard: str = "project"              # "project" (thesis Section 4.6) or "none" (diagnostic only)
    # additions for the D'Esterre (2024) presets; all off by default
    c_mc: float = 0.0                       # Mohr-Coulomb c' (kPa), model "mc"
    phi_mc: float = 0.0                     # Mohr-Coulomb phi' (deg), model "mc"
    ru: Optional[float] = None              # pore-pressure ratio: u = ru gamma h (None = dry or water table)
    base_circle: bool = False               # grid from exit to entry, floor columns allowed (see build_mesh)
    # column shape (addition, off by default): 3D columns across y; None = N, the brief's N x N grid
    n3y: Optional[int] = None
    # geometry, water and loading (additions, all off by default; see check_additions and the README)
    profile: Optional[tuple] = None         # face from the toe to the crest edge, ((0, 0), (x, z), ...); None = planar
    crest_slope: float = 0.0                # ground behind the crest edge rises at this angle (deg)
    plan_radius: Optional[float] = None     # toe line in plan: radius (m), > 0 concave wall, < 0 convex; None = straight
    elong: float = 1.0                      # end caps of the slip surface: half-width across y / R (1 = sphere)
    mid_length: float = 0.0                 # cylindrical middle part of the slip surface, length along y (m)
    d_zone: Optional[float] = None          # D applies within this distance (m) of the face; None = everywhere
    piezo: Optional[tuple] = None           # piezometric line ((x, z), ...), x decreasing; replaces H_w, beta_w (wet)
    crack_depth: Optional[float] = None     # tension crack at the back of the mass, depth below the ground (m)
    crack_water: float = 0.0                # share of the crack depth filled with water
    kh: float = 0.0                         # pseudo-static horizontal seismic coefficient (towards the pit)

    @property
    def hb(self) -> tuple[float, float, float]:
        mb, s, a = hoek_brown_parameters(self.gsi, self.mi, self.D)
        return mb, s, (a if self.a_override is None else float(self.a_override))

    def strength(self) -> Strength:
        mb, s, a = self.hb
        return Strength(self.sigma_ci, mb, s, a, self.model, self.nu, self.c_mc, self.phi_mc)

    def with_(self, **kw) -> "Inputs":
        return replace(self, **kw)


# ---------------------------------------------------------------------------------------------------------------
# mesh and column quantities (eqs 4.10, 4.11)
# ---------------------------------------------------------------------------------------------------------------
def ground_z(x, H: float, beta: float):
    """z_t = min(H, -x tan beta) for x <= 0; 0 on the pit floor (x > 0)."""
    x = np.asarray(x, dtype=float)
    return np.where(x <= 0, np.minimum(H, -x * math.tan(math.radians(beta))), 0.0)


def water_z(x, Hw: float, beta_w: float):
    """z_w = min(H_w, -x tan beta_w) for x <= 0; 0 on the pit floor."""
    x = np.asarray(x, dtype=float)
    return np.where(x <= 0, np.minimum(Hw, -x * math.tan(math.radians(beta_w))), 0.0)


def sphere_z(x, y, xc: float, zc: float, R: float):
    """z_b = z_c - sqrt(R^2 - (x - x_c)^2 - y^2); NaN where not real."""
    with np.errstate(invalid="ignore"):
        q = R * R - (np.asarray(x, dtype=float) - xc) ** 2 - np.asarray(y, dtype=float) ** 2
        return np.where(q >= 0, zc - np.sqrt(np.where(q >= 0, q, 0.0)), np.nan)


# ---- additions: slope profile, crest slope, walls curved in plan, piezometric line, tension crack -----------------
# All off by default; without them ground() and water() are ground_z() and water_z(), the thesis surfaces.
def _points(pts, what: str) -> np.ndarray:
    """The polyline ((x, z), ...) as an (n, 2) array: at least two finite points, x strictly decreasing."""
    if not (isinstance(pts, (tuple, list)) and len(pts) >= 2
            and all(isinstance(q, (tuple, list)) and len(q) == 2 for q in pts)):
        raise ValueError(f"the {what} needs at least two points (x, z)")
    a = np.asarray(pts, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"the {what} has a point that is not a number")
    if not np.all(a[1:, 0] < a[:-1, 0]):
        raise ValueError(f"the {what} runs back into the slope: x must decrease from point to point")
    return a


def check_additions(inp: "Inputs") -> None:
    """Validate the additions (each off by default); raise ValueError with the reason. Same as checkAdditions."""
    if inp.profile is not None:
        a = _points(inp.profile, "slope profile")
        if a[0, 0] != 0 or a[0, 1] != 0:
            raise ValueError("the slope profile starts at the toe, (0, 0)")
        if not np.all(a[1:, 1] >= a[:-1, 1]):
            raise ValueError("the slope profile cannot descend: z must not decrease from point to point")
        if not a[-1, 1] > 0:
            raise ValueError("the slope profile must end above the toe (crest edge z > 0)")
    if not (math.isfinite(inp.crest_slope) and 0 <= inp.crest_slope < 90):
        raise ValueError("require 0 <= crest slope < 90 deg")
    if inp.plan_radius is not None and not (math.isfinite(inp.plan_radius) and inp.plan_radius != 0):
        raise ValueError("the plan radius must be a number other than 0 (positive: concave wall; negative: convex)")
    if not (math.isfinite(inp.elong) and inp.elong > 0):
        raise ValueError("the slip-surface elongation must be positive")
    if not (math.isfinite(inp.mid_length) and inp.mid_length >= 0):
        raise ValueError("the length of the cylindrical middle part must be zero or positive")
    if inp.d_zone is not None and not (math.isfinite(inp.d_zone) and inp.d_zone >= 0):
        raise ValueError("the disturbed-zone thickness T must be zero or positive")
    if inp.piezo is not None:
        _points(inp.piezo, "piezometric line")
    if inp.crack_depth is not None and not (math.isfinite(inp.crack_depth) and inp.crack_depth > 0):
        raise ValueError("the tension crack depth must be positive")
    if not (math.isfinite(inp.crack_water) and 0 <= inp.crack_water <= 1):
        raise ValueError("the crack water fill must be from 0 to 1")
    if not (math.isfinite(inp.kh) and 0 <= inp.kh < 1):
        raise ValueError("require 0 <= k_h < 1")


def face_points(inp: "Inputs") -> tuple[np.ndarray, np.ndarray]:
    """The excavated face from the toe (0, 0) to the crest edge, as (x, z) arrays: the profile, or the planar face."""
    if inp.profile is None:
        return np.array([0.0, -inp.H / math.tan(math.radians(inp.beta))]), np.array([0.0, inp.H])
    a = np.asarray(inp.profile, dtype=float)
    return a[:, 0].copy(), a[:, 1].copy()


def face_height(inp: "Inputs") -> float:
    """Height of the crest edge above the toe: H, or the last point of the profile."""
    return float(inp.profile[-1][1]) if inp.profile is not None else inp.H


def _interp_desc(x, px: np.ndarray, pz: np.ndarray):
    """Piecewise-linear z through the points (px, pz), px strictly decreasing; outside [px[-1], px[0]] the end segments
    are extended. Same as interpDesc in webapp/engine.js."""
    x = np.asarray(x, dtype=float)
    k = np.clip(np.searchsorted(-px, -x, side="right") - 1, 0, len(px) - 2)
    return pz[k] + (x - px[k]) * (pz[k + 1] - pz[k]) / (px[k + 1] - px[k])


def plan_x(inp: "Inputs", x, y):
    """Section coordinate of the point (x, y): its signed horizontal distance from the toe line, negative behind it.
    Straight wall: x. Wall curved in plan (``plan_radius`` R_p): the toe line is a circle of radius |R_p| about a
    vertical axis at (R_p, 0), the pit inside it (R_p > 0, concave) or outside it (R_p < 0, convex nose), and the
    coordinate is |R_p| - r (concave) or r - |R_p| (convex), r = sqrt((x - R_p)^2 + y^2). Same as planX."""
    x = np.asarray(x, dtype=float)
    rp = inp.plan_radius
    if rp is None:
        return x
    y = np.asarray(y, dtype=float)
    r = np.sqrt((x - rp) * (x - rp) + y * y)
    return (rp - r) if rp > 0 else (r + rp)


def ground(inp: "Inputs", x, y=None):
    """Ground surface z_t: the pit floor (z = 0) ahead of the toe, the face (planar at beta, or the profile) and the
    crest behind the crest edge (flat, or rising at crest_slope). With y (3D) on a wall curved in plan, the section is
    swept around the plan axis (plan_x). Without additions it is ground_z, the thesis ground. Same as ground in
    webapp/engine.js."""
    xs = np.asarray(x, dtype=float) if y is None else plan_x(inp, x, y)
    if inp.profile is None and inp.crest_slope == 0:
        return ground_z(xs, inp.H, inp.beta)
    px, pz = face_points(inp)
    xe, he = px[-1], pz[-1]
    tc = math.tan(math.radians(inp.crest_slope))
    face = -xs * math.tan(math.radians(inp.beta)) if inp.profile is None else _interp_desc(xs, px, pz)
    return np.where(xs > 0, 0.0, np.where(xs >= xe, face, he + (xe - xs) * tc))


def water(inp: "Inputs", x, y=None):
    """Water table z_w: min(H_w, -x tan beta_w) behind the toe (water_z), or the piezometric line (``piezo``, constant
    beyond its end points). With y on a wall curved in plan, swept as the ground. Same as water in webapp/engine.js."""
    xs = np.asarray(x, dtype=float) if y is None else plan_x(inp, x, y)
    if inp.piezo is None:
        return water_z(xs, inp.Hw, inp.beta_w)
    a = np.asarray(inp.piezo, dtype=float)
    px, pz = a[:, 0], a[:, 1]
    return np.where(xs >= px[0], pz[0], np.where(xs <= px[-1], pz[-1], _interp_desc(xs, px, pz)))


def face_distance(inp: "Inputs", x, z):
    """Shortest distance (m) from the section points (x, z) to the excavated face (toe to crest edge). Same as
    faceDistance in webapp/engine.js."""
    px, pz = face_points(inp)
    x, z = np.asarray(x, dtype=float), np.asarray(z, dtype=float)
    best = np.full(np.broadcast(x, z).shape, np.inf)
    for k in range(len(px) - 1):
        ax, az = px[k], pz[k]
        dx, dz = px[k + 1] - ax, pz[k + 1] - az
        t = np.minimum(np.maximum(((x - ax) * dx + (z - az) * dz) / (dx * dx + dz * dz), 0.0), 1.0)
        ex, ez = x - (ax + t * dx), z - (az + t * dz)
        best = np.minimum(best, np.sqrt(ex * ex + ez * ez))
    return best


@dataclass
class Mesh:
    """Active columns (3D half-model) or slices (2D, unit thickness). Index ``i``, ``j`` are 1-based cell numbers."""

    dim: int
    n: int
    x0: float
    x1: float
    dx: float
    dy: float
    i: np.ndarray
    j: np.ndarray
    x: np.ndarray
    y: np.ndarray
    zb: np.ndarray
    zt: np.ndarray
    h: np.ndarray
    ah: np.ndarray
    g: np.ndarray              # g_z (cos alpha in 2D)
    area: np.ndarray           # A = A_h / g_z
    rho: np.ndarray
    kap: np.ndarray            # kappa (sin alpha in 2D)
    w: np.ndarray              # kN
    u: np.ndarray              # kPa
    zw: np.ndarray
    ny: int = 1                # cells across y (3D: columns_across; 2D: 1)
    n_x: Optional[np.ndarray] = None    # unit base normal into the mass, x and y components (z component: g)
    n_y: Optional[np.ndarray] = None
    xsec: Optional[np.ndarray] = None   # section coordinate (plan_x; x on a straight wall)
    zone: Optional[np.ndarray] = None   # base within d_zone of the face (disturbance zone), else None
    crack: Optional[dict] = None        # tension crack: x and the water force on each row (crack_forces)
    ext: Optional[dict] = None          # driving terms of crack water and seismic load (external_terms), else None

    @property
    def count(self) -> int:
        return len(self.x)

    @property
    def volume(self) -> float:
        """3D: full mass (twice the half-model); 2D: m^3 per m."""
        v = float(np.sum(self.h * self.ah))
        return 2.0 * v if self.dim == 3 else v


def slip_extent(inp: Inputs, n: int = 20001) -> tuple[float, float]:
    """Entry and exit of the slip circle on the central section y = 0: the ends of the x-interval where the ground
    (pit floor z_t = 0 ahead of the toe included) is above the circle. Sampled on n points over [x_c - R, x_c + R],
    each end refined by bisection. Raises ValueError if the circle does not cut the ground."""
    check_additions(inp)
    xc, zc, R = inp.xc, inp.zc, inp.R

    def above(x):
        return bool(ground(inp, x) - sphere_z(x, 0.0, xc, zc, R) > 0)

    xs = np.linspace(xc - R, xc + R, n)
    with np.errstate(invalid="ignore"):
        d = ground(inp, xs) - sphere_z(xs, 0.0, xc, zc, R)
    pos = np.isfinite(d) & (d > 0)
    if not pos.any():
        raise ValueError("the slip circle does not cut the ground: there is no sliding mass")

    def refine(a, b):
        """Boundary between a (outside the mass) and b (inside)."""
        for _ in range(200):
            mid = 0.5 * (a + b)
            if mid in (a, b):
                break
            if above(mid):
                b = mid
            else:
                a = mid
        return 0.5 * (a + b)

    first, last = int(np.argmax(pos)), n - 1 - int(np.argmax(pos[::-1]))
    x_entry = float(xs[0]) if first == 0 else refine(float(xs[first - 1]), float(xs[first]))
    x_exit = float(xs[-1]) if last == n - 1 else refine(float(xs[last + 1]), float(xs[last]))
    return x_entry, x_exit


def crack_x(inp: Inputs, n: int = 20001) -> Optional[float]:
    """Position of the tension crack, the vertical plane x = x_tc (addition, ``crack_depth``): on the central section
    y = 0, the first point from the back of the mass (slip_extent) at which the slip circle lies crack_depth below the
    ground; the back end itself if the mass already ends in a vertical face that high. Sampled on n points, refined by
    bisection. None without a crack. Same as crackX in webapp/engine.js."""
    if inp.crack_depth is None:
        return None
    dc = inp.crack_depth
    x_entry, x_exit = slip_extent(inp)
    if not inp.base_circle:                          # the default grid's mass ends at the toe
        x_exit = min(x_exit, 0.0)
    xc, zc, R = inp.xc, inp.zc, inp.R

    def depth(x):
        return ground(inp, x) - sphere_z(x, 0.0, xc, zc, R)

    xs = np.linspace(x_entry, x_exit, n)
    with np.errstate(invalid="ignore"):
        d = depth(xs)
    deep = np.isfinite(d) & (d >= dc)
    if not deep.any():
        raise ValueError(f"the tension crack ({dc:g} m) is deeper than the sliding mass on the central section "
                         f"(at most {float(np.nanmax(d)):.4g} m deep)")
    k = int(np.argmax(deep))
    if k == 0:
        return float(xs[0])
    a, b = float(xs[k - 1]), float(xs[k])           # shallower than dc at a, at least dc at b
    for _ in range(200):
        mid = 0.5 * (a + b)
        if mid in (a, b):
            break
        if float(depth(mid)) >= dc:
            b = mid
        else:
            a = mid
    return 0.5 * (a + b)


def vertical_faces(inp: Inputs) -> dict:
    """Heights of the vertical faces that bound the sliding mass where the slip circle on y = 0 ends below the ground:
    at its back extreme x_c - R (z = z_c), and with the base-circle grid also at its front extreme x_c + R (the default
    grid stops at the toe). No force acts on such a face: at the back it works as an open tension crack. Heights below
    1e-9 H (rounding) count as zero. Same as verticalFaces in webapp/engine.js."""
    check_additions(inp)
    xc, zc, R, tol = inp.xc, inp.zc, inp.R, 1e-9 * face_height(inp)
    back = float(ground(inp, xc - R)) - zc if (inp.base_circle or xc - R <= 0) and inp.crack_depth is None else 0.0
    front = float(ground(inp, xc + R)) - zc if inp.base_circle else 0.0
    return {"back": back if back > tol else 0.0, "front": front if front > tol else 0.0}


def grid_x_range(inp: Inputs) -> tuple[float, float]:
    """The grid's x range on the central section, (x0, x1) with x0 > x1: brief §1.4, x0 = x_c - sqrt(R^2 - z_c^2) (the
    toe-level crossing) and x1 = x_c - R; with the base-circle option, the circle's exit and entry (slip_extent).
    Raises ValueError for a circle the default grid refuses (see build_mesh)."""
    xc, zc, R = inp.xc, inp.zc, inp.R
    if inp.base_circle:
        x1, x0 = slip_extent(inp)
        return x0, x1
    if not (R > zc):
        raise ValueError("require R > z_c (the sphere must cut the toe level); for a face circle, with its lowest point "
                         "above the toe level, use the base-circle option")
    x0 = xc - math.sqrt(R * R - zc * zc)
    if x0 < -1e-5 * R:
        raise ValueError(f"the circle crosses the toe level {-x0:.4g} m behind the toe, so the default grid would cut "
                         f"off the toe region; use the base-circle option")
    return x0, xc - R


def columns_across(inp: Inputs, n: int) -> int:
    """3D columns across y (y from 0 to R) for a mesh with n columns along x: n for the brief's N x N grid (n3y None);
    otherwise n3y scaled by n/n3 and rounded half up (at least 1), so a mesh study or the search's finer mesh keeps
    the column shape. At n = n3 it is n3y. Same as columnsAcross in webapp/engine.js."""
    if inp.n3y is None:
        return n
    if not (isinstance(inp.n3y, (int, np.integer)) and inp.n3y >= 1 and isinstance(inp.n3, (int, np.integer))
            and inp.n3 >= 1):
        raise ValueError("the number of columns across y must be a positive integer")
    return max(1, (2 * int(inp.n3y) * int(n) + int(inp.n3)) // (2 * int(inp.n3)))


def with_n3(inp: Inputs, n: int) -> Inputs:
    """The inputs with n columns along x and the columns across y scaled to keep their shape (mesh studies)."""
    return inp.with_(n3=n, n3y=None if inp.n3y is None else columns_across(inp, n))


def y_extent(inp: Inputs) -> float:
    """Half-width of the slip surface across y (m), the 3D grid's y range: R for the sphere, lam R + L/2 with the
    ellipsoidal caps (elong lam) and the cylindrical middle part (mid_length L)."""
    half = 0.5 * inp.mid_length
    return inp.elong * inp.R + half if half > 0 else inp.elong * inp.R


def column_size(inp: Inputs) -> tuple[float, float]:
    """Plan size of the 3D columns, (dx, dy) in m: dx = (x0 - x1)/N along x, dy = y_extent/N_y across y (R/N_y for
    the sphere)."""
    x0, x1 = grid_x_range(inp)
    return (x0 - x1) / inp.n3, y_extent(inp) / columns_across(inp, inp.n3)


def counts_for_size(inp: Inputs, dx: float, dy: float) -> tuple[int, int]:
    """Columns along x and across y whose plan size is at most dx by dy (m) on the current geometry:
    N = ceil((x0 - x1)/dx), N_y = ceil(y_extent/dy), each at least 1 (a relative 1e-9 below a whole number counts as
    it)."""
    if not (dx > 0 and dy > 0):
        raise ValueError("the column length and width must be positive")
    x0, x1 = grid_x_range(inp)
    return (max(1, math.ceil((x0 - x1) / dx * (1 - 1e-9))), max(1, math.ceil(y_extent(inp) / dy * (1 - 1e-9))))


def build_mesh(inp: Inputs, dim: int, n: int) -> Mesh:
    """Mesh (brief §1.4): x from x0 = x_c - sqrt(R^2 - z_c^2) to x1 = x_c - R, dx = (x0 - x1)/N, x_i = x0 - (i - 1/2) dx;
    3D: y from 0 to R, dy = R/N, y_j = (j - 1/2) dy. 2D (brief §1.5): y = 0, unit thickness, rho = R.
    Active: real z_b (eq. 4.10; R^2 - (x - x_c)^2 - y^2 > 0, so g_z > 0), x <= 0 and h = z_t - z_b > 0.
    Column quantities (eq. 4.11): A_h, g_z, A = A_h/g_z, W = gamma h A_h, rho = sqrt(R^2 - y^2), kappa = (x_c - x)/rho.

    The default grid needs the circle to cross the toe level at or ahead of the toe (x0 >= -1e-5 R, the tolerance of
    the thesis workbook, which admits the published circle's 2.6 mm rounding); otherwise it would cut off the toe
    region, so it is refused.

    ``base_circle=True`` (off by default): x runs from the circle's exit x0 to its entry x1 on the central section
    (``slip_extent``), with the same cell numbering, and floor columns are allowed: x > 0 with z_t = 0 and h = -z_b > 0.
    Pore pressure: ``ru`` gives u = ru gamma h for every column or slice; otherwise the water table (``wet``) or dry.
    ``n3y`` (off by default): N_y columns across y instead of N (``columns_across``), dy = R/N_y; dx is unchanged.
    """
    if dim not in (2, 3):
        raise ValueError("dim must be 2 or 3")
    if not (isinstance(n, (int, np.integer)) and n >= 1):
        raise ValueError("mesh size must be a positive integer")
    if inp.ru is not None:
        if inp.wet:
            raise ValueError("choose a water table (wet) or r_u, not both")
        if not 0 <= inp.ru < 1:
            raise ValueError("require 0 <= r_u < 1")
    check_additions(inp)
    xc, zc, R = inp.xc, inp.zc, inp.R
    x0, x1 = grid_x_range(inp)
    if dim == 3 and inp.plan_radius is not None and inp.plan_radius < 0 and not x1 > inp.plan_radius:
        raise ValueError(f"a convex wall needs its plan axis behind the sliding mass: the plan radius must exceed "
                         f"{-x1:.4g} m, the distance from the toe to the back of the grid")
    dx = (x0 - x1) / n
    ny = columns_across(inp, n) if dim == 3 else 1
    lam = inp.elong if dim == 3 else 1.0             # 2D: the central section is the circle for any elongation
    half = 0.5 * inp.mid_length if dim == 3 else 0.0  # half-length of the cylindrical middle part
    x_tc = crack_x(inp)
    if dim == 3:
        dy = y_extent(inp) / ny
        ii, jj = np.meshgrid(np.arange(1, n + 1), np.arange(1, ny + 1), indexing="ij")
        ii, jj = ii.ravel(), jj.ravel()
        y = (jj - 0.5) * dy
    else:
        dy = 1.0
        ii = np.arange(1, n + 1)
        jj = np.ones(n, dtype=int)
        y = np.zeros(n)
    x = x0 - (ii - 0.5) * dx
    yc = np.maximum(y - half, 0.0) if half > 0 else y    # distance into the end cap (y beyond the middle part)
    q = R * R - (x - xc) ** 2 - (yc / lam) ** 2
    xs = plan_x(inp, x, y) if dim == 3 else x
    keep = (q > 0) if inp.base_circle else (q > 0) & (xs <= 0)
    if x_tc is not None:
        keep &= x >= x_tc                            # the mass ends at the tension crack
    ii, jj, x, y, yc, q, xs = ii[keep], jj[keep], x[keep], y[keep], yc[keep], q[keep], xs[keep]
    zb = zc - np.sqrt(q)
    zt = ground(inp, x, y) if dim == 3 else ground(inp, x)
    h = zt - zb
    keep = h > 0
    ii, jj, x, y, yc, q, xs, zb, zt, h = (v[keep] for v in (ii, jj, x, y, yc, q, xs, zb, zt, h))
    ah = np.full(len(x), dx * dy)
    if lam == 1.0 and half == 0:
        g = np.sqrt(q) / R
        n_x, n_y = (xc - x) / R, -y / R
        rho = np.sqrt(R * R - y * y) if dim == 3 else np.full(len(x), R)
    else:                     # ellipsoidal caps: normal along (x_c - x, -y_c/lam^2, z_c - z_b), y_c = max(y - L/2, 0)
        yl = yc / (lam * lam)
        nrm = np.sqrt((x - xc) ** 2 + yl * yl + q)
        g = np.sqrt(q) / nrm
        n_x, n_y = (xc - x) / nrm, -yl / nrm
        rho = np.sqrt(R * R - (yc / lam) ** 2)
    kap = (xc - x) / rho
    if inp.ru is not None:
        zw = np.full(len(x), np.nan)
        u = inp.ru * inp.gamma * h
    elif inp.wet:
        zw = water(inp, x, y) if dim == 3 else water(inp, x)
        u = GAMMA_W * np.maximum(0.0, zw - zb)
    else:
        zw = np.full(len(x), np.nan)
        u = np.zeros(len(x))
    zone = face_distance(inp, xs, zb) <= inp.d_zone if inp.d_zone is not None else None
    m = Mesh(dim, n, x0, x1, dx, dy, ii, jj, x, y, zb, zt, h, ah, g, ah / g, rho, kap,
             inp.gamma * h * ah, u, zw, ny, n_x, n_y, xs, zone)
    if x_tc is not None:
        m.crack = crack_forces(inp, m, x_tc)
    if inp.kh > 0 or (m.crack is not None and m.crack["P"] > 0):
        m.ext = external_terms(inp, m)
    return m


def crack_forces(inp: Inputs, m: Mesh, x_tc: float) -> dict:
    """Water in the tension crack (addition): a horizontal force on the crack face x = x_tc, towards the pit. On the
    central section the water stands crack_water times the crack depth above the crack bottom, and that level holds
    across y. On each row (2D: the section, per metre; 3D: each strip dy wide of the half-model) the force is
    P = gamma_w h_w^2 / 2 (times dy) at h_w/3 above the crack bottom, with the lever arm e = z_c - z_P about the
    rotation axis; rho is the row's radius. Same as crackForces in webapp/engine.js."""
    xc, zc, R = inp.xc, inp.zc, inp.R
    lam = inp.elong if m.dim == 3 else 1.0
    half = 0.5 * inp.mid_length if m.dim == 3 else 0.0
    d0 = x_tc - xc
    zb0 = zc - math.sqrt(R * R - d0 * d0)
    level = zb0 + inp.crack_water * (float(ground(inp, x_tc)) - zb0)
    ys = [(j - 0.5) * m.dy for j in range(1, m.ny + 1)] if m.dim == 3 else [0.0]
    width = m.dy if m.dim == 3 else 1.0
    rows = []
    for y in ys:
        yl = (max(y - half, 0.0) if half > 0 else y) / lam
        qq = R * R - d0 * d0 - yl * yl
        if not qq > 0:
            continue
        zbot = zc - math.sqrt(qq)
        ztop = float(ground(inp, x_tc, y)) if m.dim == 3 else float(ground(inp, x_tc))
        if not ztop > zbot:
            continue
        hw = min(max(level - zbot, 0.0), ztop - zbot)
        zp = zbot + hw / 3.0
        rows.append({"y": y, "zbot": zbot, "ztop": ztop, "hw": hw, "P": 0.5 * GAMMA_W * hw * hw * width, "zP": zp,
                     "e": zc - zp, "rho": math.sqrt(R * R - yl * yl) if m.dim == 3 else R})
    tot = mh = mm = 0.0
    for r in rows:
        tot += r["P"]
        mh += r["P"] * r["e"] / r["rho"]
        mm += r["P"] * r["e"]
    return {"x": x_tc, "level": level, "rows": rows, "P": tot, "MH": mh, "MM": mm}


def external_terms(inp: Inputs, m: Mesh) -> dict:
    """Driving terms of the additions in the global equation (eqs 4.14 and 4.15, extended): the seismic force k_h W
    of each column or slice, horizontal towards the pit at its centroid (z_b + z_t)/2, and the crack water force
    (crack_forces). A horizontal force Q at height z_Q has the moment Q (z_c - z_Q) about the rotation axis, in the
    sense of the weight's (W (x_c - x)); F_M adds these moments to sum W (x_c - x), and F_H adds each one divided by the
    lever arm rho of its column or row, as W kappa = W (x_c - x)/rho. Local vertical balance (eq. 4.12) is unchanged,
    since the forces are horizontal. Same as externalTerms in webapp/engine.js."""
    e = inp.zc - 0.5 * (m.zb + m.zt)
    E = inp.kh * m.w
    out = {"kh": inp.kh, "E": E, "e": e, "Ex": float(np.sum(E)), "EH": float(np.sum(E * e / m.rho)),
           "EM": float(np.sum(E * e)), "P": 0.0, "PH": 0.0, "PM": 0.0}
    if m.crack is not None:
        out.update(P=m.crack["P"], PH=m.crack["MH"], PM=m.crack["MM"])
    out["H"] = out["EH"] + out["PH"]
    out["M"] = out["EM"] + out["PM"]
    return out


# ---------------------------------------------------------------------------------------------------------------
# iteration (eqs 4.12-4.16) and final checks
# ---------------------------------------------------------------------------------------------------------------
@dataclass
class Check:
    name: str
    value: float
    limit: str
    passed: bool


@dataclass
class LemResult:
    """Result of one 3D or 2D solution. ``fos`` is set only when every check passes."""

    dim: int
    status: str                       # ACCEPTED / REJECTED / NOT CONVERGED / NOT APPLICABLE
    reason: str
    fos: Optional[float]
    f_trial: Optional[float]          # last trial F (diagnostic, never reported as a FoS)
    iterations: int
    form: str
    tension: str
    iteration: str
    mesh: Mesh
    strength: Optional[StrengthResult] = None
    sn: Optional[np.ndarray] = None   # final sigma_n [kPa]
    c_used: Optional[np.ndarray] = None
    tan_used: Optional[np.ndarray] = None
    tau_used: Optional[np.ndarray] = None
    tensile: Optional[np.ndarray] = None
    m_alpha: Optional[np.ndarray] = None
    local_res: Optional[np.ndarray] = None
    global_res: float = math.nan
    bishop_res: float = math.nan
    checks: list = field(default_factory=list)
    history: list = field(default_factory=list)   # (iteration, F, |dF|/F, max |d sigma_n|, tensile count, projections)
    projections: int = 0                          # trial sigma_n projected over the whole iteration (Section 4.6)
    rejected_index: Optional[np.ndarray] = None
    safeguard: str = "project"
    projection_limit: Optional[float] = None      # kPa
    projection_value: Optional[float] = None      # kPa
    projection_limits: Optional[tuple] = None     # (disturbed zone, beyond) with the disturbance zone, kPa

    @property
    def accepted(self) -> bool:
        return self.status == "ACCEPTED"


def _sensitivity(c, t, tau, sn, tension, tens):
    """Tension sensitivities (brief §1.6, Table 5.7) on the tensile set ``tens``: C0 removes c_i, T0 removes all
    basal shear.

    The set is cumulative: a base joins it the first time its sigma_n is negative during the iteration and stays in it
    (with c_i removed its sigma_n turns positive, so a set rebuilt from the current sigma_n alone would oscillate). The
    thesis reference grew the set across full re-solves until no new base went tensile; this gives the same values.
    """
    c_e, t_e, tau_e = c.copy(), t.copy(), tau.copy()
    if tension == "C0":
        c_e[tens] = 0.0
        tau_e[tens] = sn[tens] * t[tens]
    elif tension == "T0":
        c_e[tens] = 0.0
        t_e[tens] = 0.0
        tau_e[tens] = 0.0
    return c_e, t_e, tau_e


def _drive(m: Mesh, xc: float, form: str):
    """Denominator of the global equation: sum(W kappa) for F_H, sum(W (x_c - x)) for F_M, plus the driving terms of
    the additions (external_terms) when there are any."""
    base = np.sum(m.w * (xc - m.x)) if form == "FM" else np.sum(m.w * m.kap)
    if m.ext is None:
        return base
    return base + (m.ext["M"] if form == "FM" else m.ext["H"])


def _global_f(m: Mesh, tau_e, xc: float, form: str) -> float:
    """F_H = sum(tau_s A)/sum(W kappa) (eq. 4.14, governing); F_M = sum(tau_s A rho)/sum(W (x_c - x)) (eq. 4.15).
    With seismic load or crack water the denominators gain their driving terms (_drive)."""
    if form == "FM":
        return float(np.sum(tau_e * m.area * m.rho) / _drive(m, xc, form))
    return float(np.sum(tau_e * m.area) / _drive(m, xc, form))


def _local_sn(m: Mesh, c_e, t_e, F):
    """Local vertical equilibrium (eq. 4.12) and the sigma_n update (eq. 4.13): m_alpha = g_z + kappa tan(phi)/F,
    sigma_n = (W - u A_h - c A kappa/F)/(A m_alpha)."""
    ma = m.g + m.kap * t_e / F
    return (m.w - m.u * m.ah - c_e * m.area * m.kap / F) / (m.area * ma), ma


def projection_limit(st: Strength) -> Optional[float]:
    """Projection limit of thesis Section 4.6 [kPa]: max(-1000 s sigma_ci/m_b, sigma_n at the admissible edge of the
    closure). For GHB and the sigma2 = sigma3 closure the admissible range starts at the tensile limit itself.
    None for Mohr-Coulomb, whose strength exists at every sigma_n (no projection)."""
    if st.model == "mc":
        return None
    tensile = st.sigma3_tensile * 1000.0
    if st.model in ("ghb", "ejc_sigma2_eq_sigma3"):
        return tensile
    edge = admissible_sn_min(st)
    return max(tensile, edge) if math.isfinite(edge) else tensile


@dataclass(frozen=True, eq=False)
class ZonedStrength:
    """Disturbance zone (addition, ``d_zone``): the strength ``inner`` (D as given) on the bases in ``zone`` and
    ``outer`` (D = 0) on the others. Same as zonedStrength in webapp/engine.js."""

    inner: Strength
    outer: Strength
    zone: np.ndarray

    @property
    def model(self) -> str:
        return self.inner.model


def strength_at(sn_kpa, st) -> StrengthResult:
    """solve_strength for a Strength, or for a ZonedStrength each zone with its own strength."""
    if not isinstance(st, ZonedStrength):
        return solve_strength(sn_kpa, st)
    z = st.zone
    if z.all():
        return solve_strength(sn_kpa, st.inner)
    if not z.any():
        return solve_strength(sn_kpa, st.outer)
    t = np.atleast_1d(np.asarray(sn_kpa, dtype=float))
    a, b = solve_strength(t[z], st.inner), solve_strength(t[~z], st.outer)
    out = {}
    for f in fields(StrengthResult):
        va, vb = getattr(a, f.name), getattr(b, f.name)
        arr = np.empty(len(t), dtype=np.result_type(va, vb))
        arr[z], arr[~z] = va, vb
        out[f.name] = arr
    return StrengthResult(**out)


def strength_for(inp: "Inputs", m: Mesh):
    """The strength of the run: inp.strength(), or with the disturbance zone a ZonedStrength (D within d_zone of the
    face, D = 0 beyond). Mohr-Coulomb has no D, so the zone does not apply to it."""
    st = inp.strength()
    if m.zone is None or st.model == "mc":
        return st
    return ZonedStrength(st, inp.with_(D=0.0).strength(), m.zone)


ITERATIONS = {"default": "Default (brief §1.4): F from F_H/F_M at the current strength, then sigma_n",
              "native_v15": "Native v15 iteration path (start F = 1, Hungr explicit form, projection)"}


def solve(m: Mesh, st: Strength, xc: float, gamma: float, form: str = "FH", tension: str = "none",
          iteration: str = "default", safeguard: str = "project", max_iter: int = MAX_ITER) -> LemResult:
    """Solve F for one mesh.

    ``iteration="default"`` (brief §1.4): sigma_n0 = (gamma h - u) g_z^2; each iteration evaluates the strength at
    sigma_n, then F from the selected global form (F_H eq. 4.14, or F_M eq. 4.15 with its own iteration), then sigma_n
    from local vertical equilibrium at that F (eqs 4.12, 4.13). It stops when |dF|/F < 1e-12 and
    max |d sigma_n| < 1e-9 kPa (2000 iterations at most).

    ``iteration="native_v15"``: the native v15 iteration path of the thesis workbook. F starts at 1, sigma_n is updated
    with the previous F, and the new F is the Hungr explicit form sum[(c A_h + (W - u A_h) tan phi)/m_alpha]/sum(W kappa)
    (equal to F_H at the updated sigma_n). Projection is always on. F_H only, no C0/T0 (those values come from the
    default path, Table 5.7).

    Projection (thesis Section 4.6, ``safeguard="project"``, the default): limit = max(-1000 s sigma_ci/m_b,
    sigma_n at the admissible edge of the closure); a trial sigma_n <= limit is set to limit + max(1e-9, 0.001 |limit|)
    and counted. A solution is accepted only if no base is projected at convergence. ``safeguard="none"`` turns it off
    (diagnostic only: such a base is then REJECTED).

    In 2D (brief §1.5) the same equations hold with rho = R, g_z = cos(alpha), kappa = sin(alpha); the result is also
    checked against Bishop's expression (eq. 4.16).

    Mohr-Coulomb (``model="mc"``): no projection, and tau_s < 0 is allowed on the iteration path (as in the thesis
    workbook); tau_s >= 0 is required of the final state.
    """
    if form not in ("FH", "FM"):
        raise ValueError("form must be FH or FM")
    if tension not in ("none", "C0", "T0"):
        raise ValueError("tension must be none, C0 or T0")
    if iteration not in ITERATIONS:
        raise ValueError("iteration must be default or native_v15")
    if safeguard not in ("none", "project"):
        raise ValueError("safeguard must be none or project")
    if iteration == "native_v15":
        if form != "FH" or tension != "none":
            raise ValueError("the native v15 iteration path supports F_H only, without C0/T0")
        safeguard = "project"
    res = LemResult(m.dim, "REJECTED", "", None, None, 0, form, tension, iteration, m, safeguard=safeguard)
    if st.model == "ejc_original":
        res.status = "NOT APPLICABLE"
        res.reason = ("The original closure is an admissibility check only and never produces a FoS "
                      "(see original_closure_minimum).")
        return res
    if m.count == 0:
        res.reason = "No active columns."
        return res
    drive = _drive(m, xc, form)
    if not drive > 0:
        res.reason = "Non-positive driving sum."
        return res
    sn = (gamma * m.h - m.u) * m.g ** 2
    if safeguard != "project":
        limit = None
    elif isinstance(st, ZonedStrength):              # disturbance zone: each base with its zone's limit
        li, lo = projection_limit(st.inner), projection_limit(st.outer)
        limit = np.where(st.zone, li, lo)
        value = limit + np.maximum(1e-9, 0.001 * np.abs(limit))
        res.projection_limit, res.projection_value = li, li + max(1e-9, 0.001 * abs(li))
        res.projection_limits = (li, lo)
    else:
        limit = projection_limit(st)
    if limit is not None:
        if not isinstance(st, ZonedStrength):
            value = limit + max(1e-9, 0.001 * abs(limit))
            res.projection_limit, res.projection_value = limit, value

        def project(v):
            p = v <= limit
            return np.where(p, value, v), int(p.sum())
    else:
        def project(v):
            return v, 0
    sn, res.projections = project(sn)
    f_old = 1.0 if iteration == "native_v15" else None
    tens = np.zeros(m.count, dtype=bool)
    for it in range(1, max_iter + 1):
        stg = strength_at(sn, st)
        res.iterations = it
        failed = ~stg.accepted
        if st.model == "mc":
            failed &= stg.code != SHEAR             # tau_s < 0 on the path: checked at the final state
        if failed.any():
            bad = np.flatnonzero(failed)
            res.strength, res.sn, res.rejected_index = stg, sn, bad
            res.reason = (f"{len(bad)} of {m.count} bases have no accepted local state at iteration {it}; "
                          f"first: column {int(bad[0])} at sigma_n = {sn[bad[0]]:.6g} kPa, {stg.reason(int(bad[0]))}")
            return res
        t = np.tan(stg.phi)
        if tension != "none":
            tens = tens | (sn < 0)
        c_e, t_e, tau_e = _sensitivity(stg.c, t, stg.tau, sn, tension, tens)
        if iteration == "default":
            F = _global_f(m, tau_e, xc, form)
            if not (math.isfinite(F) and F > 0):
                res.reason = f"Non-positive or non-finite trial F at iteration {it}."
                return res
            new, ma = _local_sn(m, c_e, t_e, F)
            new, n_proj = project(new)
        else:
            new, ma = _local_sn(m, c_e, t_e, f_old)
            new, n_proj = project(new)
            F = _global_f(m, c_e + new * t_e, xc, form)
            if not (math.isfinite(F) and F > 0):
                res.reason = f"Non-positive or non-finite trial F at iteration {it}."
                return res
        res.projections += n_proj
        if (ma <= 0).any():
            res.reason = f"Non-positive m_alpha at iteration {it}."
            return res
        prev = f_old if iteration == "native_v15" else (res.history[-1][1] if res.history else None)
        dF = abs(F - prev) / F if prev is not None else math.inf
        move = float(np.max(np.abs(new - sn)))
        res.history.append((it, F, dF, move, int(tens.sum()), n_proj))
        res.f_trial = F
        sn = new
        f_old = F
        if dF < TOL_DF and move < TOL_DSN:
            if iteration == "native_v15":       # the final state is rebuilt at the converged F
                sn, _ = _local_sn(m, c_e, t_e, F)
            return _final_checks(res, m, st, xc, F, sn, tens, dF, move, limit)
    res.status = "NOT CONVERGED"
    res.reason = f"No convergence in {max_iter} iterations; no accepted FoS."
    res.sn = sn
    return res


def _final_checks(res: LemResult, m: Mesh, st, xc: float, F: float, sn: np.ndarray,
                  tens_used: np.ndarray, dF: float, move: float, limit=None) -> LemResult:
    """Independent check of the converged state: fresh strength at the final sigma_n, then the local (vertical)
    and global residuals at the reported F (brief §1.4 step 4)."""
    stg = strength_at(sn, st)
    t = np.tan(stg.phi)
    tens = tens_used | (sn < 0) if res.tension != "none" else tens_used
    c_e, t_e, tau_e = _sensitivity(stg.c, t, stg.tau, sn, res.tension, tens_used)
    ma = m.g + m.kap * t_e / F
    with np.errstate(invalid="ignore"):
        local = np.abs(m.w - (sn + m.u) * m.ah - tau_e * m.area * m.kap / F) / m.w
    f_glob = _global_f(m, tau_e, xc, res.form)
    glob = abs(f_glob - F) / F
    bishop = float(np.sum((c_e * m.ah + (m.w - m.u * m.ah) * t_e) / ma) / _drive(m, xc, "FH"))
    bishop_res = abs(bishop - F) / F
    acc = stg.accepted
    finite = np.isfinite(stg.f)
    checks = [
        Check("Convergence |dF|/F", dF, "< 1e-12", dF < TOL_DF),
        Check("Convergence max |d sigma_n| (kPa)", move, "< 1e-9", move < TOL_DSN),
        Check("Local states accepted (count rejected)", float((~acc).sum()), "0", bool(acc.all())),
    ]
    if st.model != "mc":                     # EJC residuals (Mohr-Coulomb has none)
        checks += [
            Check("Max |f| unsquared EJC (MPa)", float(np.nanmax(np.abs(stg.f))) if finite.any() else math.nan,
                  "< 1e-8", bool(finite.all() and np.all(np.abs(stg.f) < TOL_EJC))),
            Check("Max closure residual (MPa)", float(np.nanmax(np.abs(stg.closure))), "< 1e-8",
                  bool(np.all(np.abs(stg.closure) < TOL_CLOSURE))),
        ]
    checks += [
        Check("Ordering sigma1 >= sigma2 >= sigma3 (failures)",
              float(np.sum(~((stg.s1 >= stg.s2) & (stg.s2 >= stg.s3)))), "0",
              bool(np.all((stg.s1 >= stg.s2) & (stg.s2 >= stg.s3)))),
        Check("Min M (gradient)", float(np.nanmin(stg.M)), ">= 1, finite",
              bool(np.all(np.isfinite(stg.M) & (stg.M >= 1)))),
        Check("Min tau_s (kPa)", float(np.nanmin(stg.tau)), ">= 0", bool(np.all(stg.tau >= 0))),
        Check("Min m_alpha", float(np.min(ma)), "> 0", bool(np.all(ma > 0))),
        Check("Max relative local residual", float(np.nanmax(local)), "<= 1e-9",
              bool(np.all(local <= TOL_LOCAL))),
        Check(f"Relative global residual ({res.form})", glob, "<= 1e-8", glob <= TOL_GLOBAL),
    ]
    if m.dim == 2:
        checks.append(Check("Relative Bishop residual", bishop_res, "<= 1e-8", bishop_res <= TOL_GLOBAL))
    if res.tension in ("C0", "T0"):
        same = bool(np.array_equal(tens, tens_used))
        checks.append(Check("New tensile bases at the final state", float(np.sum(tens != tens_used)), "0", same))
    if res.projection_limit is not None:
        last_proj = res.history[-1][5] if res.history else 0
        checks.append(Check("Projected bases at convergence", float(last_proj), "0", last_proj == 0))
        if res.projection_limits is None:
            checks.append(Check("Min sigma_n above the projection limit (kPa)",
                                float(np.min(sn) - res.projection_limit), "> 0",
                                bool(np.all(sn > res.projection_limit))))
        else:                                    # disturbance zone: each base against its zone's limit
            checks.append(Check("Min sigma_n above the projection limit (kPa)", float(np.min(sn - limit)), "> 0",
                                bool(np.all(sn > limit))))
    res.strength, res.sn, res.c_used, res.tan_used, res.tau_used = stg, sn, c_e, t_e, tau_e
    res.tensile, res.m_alpha, res.local_res = tens_used, ma, local
    res.global_res, res.bishop_res, res.checks = glob, bishop_res, checks
    res.f_trial = F
    if all(c.passed for c in checks):
        res.status, res.fos = "ACCEPTED", F
        res.reason = "Every local state, local and global equilibrium check passed."
        if res.tension != "none":
            res.reason += f" {res.tension}: sensitivity, not a calibrated contact law."
    else:
        failed = [c.name for c in checks if not c.passed]
        res.status = "REJECTED"
        res.reason = "Final-state check failed: " + "; ".join(failed)
    return res


# ---------------------------------------------------------------------------------------------------------------
# one case = 3D + 2D
# ---------------------------------------------------------------------------------------------------------------
@dataclass
class Analysis:
    inputs: Inputs
    r3: LemResult
    r2: LemResult

    @property
    def ratio(self) -> Optional[float]:
        if self.r3.fos is None or self.r2.fos is None:
            return None
        return self.r3.fos / self.r2.fos


def run_3d(inp: Inputs) -> LemResult:
    m = build_mesh(inp, 3, inp.n3)
    return solve(m, strength_for(inp, m), inp.xc, inp.gamma, inp.form, inp.tension, inp.iteration, inp.safeguard)


def run_2d(inp: Inputs) -> LemResult:
    """Bishop on the central section (brief §1.5, eq. 4.16); F_M = F_H there (rho = R), so F_H is used."""
    m = build_mesh(inp, 2, inp.n2)
    return solve(m, strength_for(inp, m), inp.xc, inp.gamma, "FH", inp.tension, inp.iteration, inp.safeguard)


def run(inp: Inputs) -> Analysis:
    return Analysis(inp, run_3d(inp), run_2d(inp))


def mesh_error(fos_values: list[float], fos_finest: float) -> float:
    """Eq. 4.17: epsilon_mesh = (F_max - F_min)/F_finest x 100 % (0.5 % screening threshold)."""
    return (max(fos_values) - min(fos_values)) / fos_finest * 100.0


def gravity_cap(m: Mesh, gamma: float) -> np.ndarray:
    """Gravity cap sigma_n <= gamma h - u for each column [kPa]."""
    return gamma * m.h - m.u


def columns_table(r: LemResult) -> dict[str, np.ndarray]:
    """Per-column results as named columns (for CSV export)."""
    m = r.mesh
    out = {"i": m.i, "j": m.j, "x_m": m.x, "y_m": m.y, "z_base_m": m.zb, "z_top_m": m.zt, "h_m": m.h,
           "dx_m": np.full(m.count, m.dx), "dy_m": np.full(m.count, m.dy), "A_h_m2": m.ah, "A_m2": m.area,
           "g_z": m.g, "rho_m": m.rho, "kappa": m.kap, "W_kN": m.w, "u_kPa": m.u}
    if r.strength is not None and r.sn is not None:
        s = r.strength
        out.update({"sigma_n_kPa": r.sn, "phi_i_deg": np.degrees(s.phi), "c_i_kPa": s.c, "tau_s_kPa": s.tau,
                    "M": s.M, "sigma1_MPa": s.s1, "sigma2_MPa": s.s2, "sigma3_MPa": s.s3, "f_MPa": s.f,
                    "closure_res_MPa": s.closure, "state": np.array([REASONS[int(k)] for k in s.code]),
                    "n_brackets": s.n_brackets})
        if r.tau_used is not None:
            F = r.f_trial
            out.update({"c_used_kPa": r.c_used, "tan_phi_used": r.tan_used, "tau_used_kPa": r.tau_used,
                        "tensile": r.tensile, "m_alpha": r.m_alpha,
                        "N_prime_kN": r.sn * m.area, "U_kN": m.u * m.area,
                        "S_mobilised_kN": r.tau_used * m.area / F,
                        "resistance_kN": r.tau_used * m.area, "local_residual_rel": r.local_res})
    return out


def mass_forces(inp: Inputs, r: LemResult) -> dict:
    """All forces on the sliding mass at the reported F (free-body diagram of the whole mass).

    For every column or slice: W (vertical), N' = sigma_n A and U = u A along the base normal n, S = tau_s A / F along
    s, as plots.column_forces; their totals for the full mass (both halves in 3D, per metre in 2D) and their moments
    about the rotation axis (parallel to y through (x_c, z_c)), positive in the sense in which the weight drives. N' and
    U act along n, through the sphere centre, so their moments are zero; S has lever arm rho. Vertical balance holds
    element by element (eq. 4.12); horizontal force balance is not part of the method. The method's own global
    equation is sum S = sum W kappa for F_H (eq. 4.14) and sum S rho = sum W (x_c - x) for F_M (eq. 4.15).
    With the additions: the seismic forces E = k_h W and the crack water force P (horizontal, towards the pit) join
    the totals and the moments, and the method's equation gains their driving terms (external_terms).
    Same as massForces in webapp/engine.js.
    """
    m = r.mesh
    F = r.fos if r.fos is not None else r.f_trial
    if r.sn is None or r.tau_used is None or F is None:
        raise ValueError("no converged state for this run")
    mult = 2.0 if m.dim == 3 else 1.0
    nx, ny, nz = m.n_x, m.n_y, m.g.copy()
    sx, sz = -(inp.zc - m.zb) / m.rho, (inp.xc - m.x) / m.rho
    W, N, U = m.w.copy(), r.sn * m.area, m.u * m.area
    S = r.tau_used * m.area / F

    def tot(f, d):
        return float(mult * np.sum(f * d))

    vert = W - (N + U) * nz - S * sz
    totals = {"W": {"x": 0.0, "z": float(-mult * np.sum(W))}, "N": {"x": tot(N, nx), "z": tot(N, nz)},
              "U": {"x": tot(U, nx), "z": tot(U, nz)}, "S": {"x": tot(S, sx), "z": tot(S, sz)},
              "net": {"x": float(mult * np.sum((N + U) * nx + S * sx)), "z": float(-mult * np.sum(vert))}}
    m_w = float(mult * np.sum(W * (inp.xc - m.x)))
    m_s = float(-mult * np.sum(S * m.rho))
    fm = r.form == "FM"
    lhs = -m_s if fm else float(mult * np.sum(S))
    rhs = m_w if fm else float(mult * np.sum(W * m.kap))
    moments = {"W": m_w, "N": 0.0, "U": 0.0, "S": m_s, "net": m_w + m_s}
    el = {"x": m.x, "y": m.y, "zb": m.zb, "zm": (m.zb + m.zt) / 2, "W": W, "N": N, "U": U, "S": S, "nx": nx,
          "ny": ny, "nz": nz, "sx": sx, "sz": sz}
    ext = m.ext
    if ext is not None:                          # seismic load and crack water (additions)
        totals["E"] = {"x": float(mult * ext["Ex"]), "z": 0.0}
        totals["P"] = {"x": float(mult * ext["P"]), "z": 0.0}
        totals["net"]["x"] += totals["E"]["x"] + totals["P"]["x"]
        moments["E"], moments["P"] = float(mult * ext["EM"]), float(mult * ext["PM"])
        moments["net"] += moments["E"] + moments["P"]
        rhs += float(mult * (ext["M"] if fm else ext["H"]))
        el["E"] = ext["E"]
    return {"F": F, "dim": m.dim, "mult": mult, "count": m.count, "form": r.form, "weight": -totals["W"]["z"],
            "totals": totals, "moments": moments,
            "method": {"form": r.form, "lhs": lhs, "rhs": rhs, "residual": (lhs - rhs) / rhs},
            "vertical_max_rel": float(np.max(np.abs(vert) / W)) if m.count else 0.0, "el": el,
            "crack": m.crack}


# ---------------------------------------------------------------------------------------------------------------
# critical-surface search (QA/QC finding R1)
# ---------------------------------------------------------------------------------------------------------------
SEARCH = {"back": 2.0, "front": 2.0, "grid": (4, 9, 4), "grid_base": (6, 9, 4), "starts": 2, "tol": 1e-3,
          "max_eval": 100, "restarts": 2, "coverage": 0.5}
"""Settings of critical_search: limits of the exits (multiples of H behind the crest edge, and ahead of the toe for
base circles), grid of starting surfaces (front exit x back exit x curvature; grid_base with the base-circle option),
Nelder-Mead starts, tolerance on the normalised parameters, evaluations per run, restarts, and the least share of the
default grid's span that the mass must cover. Same as SEARCH in webapp/engine.js."""


def ground_at(s: float, H: float, tan_b: float) -> tuple[float, float]:
    """Point of the ground surface at s (m): s <= 0 on the pit floor ahead of the toe (x = -s), 0 < s <= H up the face
    by height (z = s), s > H on the crest at s - H behind the crest edge."""
    if s <= 0:
        return (-s, 0.0)
    if s <= H:
        return (-s / tan_b, s)
    return (-H / tan_b - (s - H), H)


def ground_point(inp: Inputs, s: float) -> tuple[float, float]:
    """Point of the ground at the search parameter s (m), as ground_at: the pit floor for s <= 0 (x = -s), the face
    for 0 < s <= h_e (h_e the crest edge height) and the crest behind the edge for s > h_e, at s - h_e behind it
    (rising at crest_slope). The planar face is followed by height (ground_at); a profile by its length, scaled so
    that s = h_e at the crest edge (so the search reaches the berms). Same as groundPoint in webapp/engine.js."""
    if inp.profile is None and inp.crest_slope == 0:
        return ground_at(s, inp.H, math.tan(math.radians(inp.beta)))
    he = face_height(inp)
    if s <= 0:
        return (-s, 0.0)
    px, pz = face_points(inp)
    if s <= he:
        if inp.profile is None:
            return (-s / math.tan(math.radians(inp.beta)), s)
        seg = [math.sqrt((px[k + 1] - px[k]) * (px[k + 1] - px[k]) + (pz[k + 1] - pz[k]) * (pz[k + 1] - pz[k]))
               for k in range(len(px) - 1)]
        total = 0.0
        for v in seg:
            total += v
        target, run = s / he * total, 0.0
        for k, v in enumerate(seg):
            if target <= run + v or k == len(seg) - 1:
                t = min(max((target - run) / v, 0.0), 1.0)
                return (float(px[k] + t * (px[k + 1] - px[k])), float(pz[k] + t * (pz[k + 1] - pz[k])))
            run += v
    return (float(px[-1]) - (s - he), he + (s - he) * math.tan(math.radians(inp.crest_slope)))


def circle_from_exits(pf, pb, t: float) -> dict:
    """Circle through the front exit pf and the back exit pb (x_b < x_f, z_b >= z_f) whose lower arc joins them. The
    centre lies on the chord's perpendicular bisector, d = (c/2)(dz/dx + 1/t - 1) above the chord's midpoint, so that
    z_c - z_b = (1/t - 1) dx/2: t = 1 makes the arc vertical at pb, t -> 0 flattens it towards the chord. For t <= 1
    the arc meets the ground at pb before its back extreme x_c - R, so the mass has no vertical back face."""
    dx = pf[0] - pb[0]
    dz = pb[1] - pf[1]
    c = math.sqrt(dx * dx + dz * dz)
    h = c / 2
    d = h * (dz / dx + 1 / t - 1)
    return {"xc": (pf[0] + pb[0]) / 2 + d * dz / c, "zc": (pf[1] + pb[1]) / 2 + d * dx / c, "R": math.sqrt(h * h + d * d)}


def search_surface(inp: Inputs, p) -> dict:
    """The surface of search parameters p = (a, b, t): front exit at s = aH, back exit at s = bH (ground_point; H the
    crest edge height), curvature t (circle_from_exits)."""
    H = face_height(inp)
    front, back = ground_point(inp, p[0] * H), ground_point(inp, p[1] * H)
    return {**circle_from_exits(front, back, p[2]), "front": front, "back": back}


def single_mass(inp: Inputs) -> bool:
    """True if the ground is above the circle on y = 0 over one interval that ends inside (x_c - R, x_c + R), sampled
    as slip_extent samples it: one sliding mass, with no vertical face at the circle's extremes (base-circle grid)."""
    xs = np.linspace(inp.xc - inp.R, inp.xc + inp.R, 20001)
    with np.errstate(invalid="ignore"):
        d = ground(inp, xs) - sphere_z(xs, 0.0, inp.xc, inp.zc, inp.R)
    on = np.isfinite(d) & (d > 0)
    runs = int(on[0]) + int(np.count_nonzero(on[1:] & ~on[:-1]))
    return runs == 1 and not on[0] and not on[-1]


def _nelder_mead(f, x0, step, tol: float, max_eval: int):
    """Nelder-Mead (reflection 1, expansion 2, contractions 1/2, shrink 1/2) from x0 and x0 + step_i e_i, stopped when
    every vertex lies within tol of the best in every parameter or after max_eval calls; ties go to the lower index.
    Same as nelderMead in webapp/engine.js."""
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] = x0[i] + step[i]
        pts.append(p)
    vals = [f(p) for p in pts]
    evals = n + 1
    while evals < max_eval:
        idx = sorted(range(n + 1), key=lambda i: (vals[i], i))
        pts = [pts[i] for i in idx]
        vals = [vals[i] for i in idx]
        spread = 0.0
        for i in range(1, n + 1):
            for j in range(n):
                spread = max(spread, abs(pts[i][j] - pts[0][j]))
        if spread < tol:
            break
        c = []
        for j in range(n):
            s = 0.0
            for i in range(n):
                s += pts[i][j]
            c.append(s / n)
        worst = pts[n]

        def at(t, c=c, worst=worst):
            return [c[j] + t * (worst[j] - c[j]) for j in range(n)]

        xr = at(-1)
        fr = f(xr)
        evals += 1
        if fr < vals[0]:
            xe = at(-2)
            fe = f(xe)
            evals += 1
            if fe < fr:
                pts[n], vals[n] = xe, fe
            else:
                pts[n], vals[n] = xr, fr
        elif fr < vals[n - 1]:
            pts[n], vals[n] = xr, fr
        else:
            outside = fr < vals[n]
            xk = at(-0.5 if outside else 0.5)
            fk = f(xk)
            evals += 1
            if (fk <= fr) if outside else (fk < vals[n]):
                pts[n], vals[n] = xk, fk
            else:
                for i in range(1, n + 1):
                    pts[i] = [pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]) for j in range(n)]
                    vals[i] = f(pts[i])
                    evals += 1
    b = 0
    for i in range(1, n + 1):
        if vals[i] < vals[b]:
            b = i
    return pts[b], vals[b]


def critical_search(inp: Inputs, dims=(2,), progress=None, **opts) -> dict:
    """Critical-surface search: the lowest accepted FoS over circles (2D, Bishop) and spheres centred on y = 0 (3D,
    Hungr) that leave the ground at a front exit (the toe or the face; with the base-circle option also the pit
    floor, up to front H ahead of the toe) and a back exit (face or crest, up to back H behind the crest edge), with
    the lower arc between them (circle_from_exits, so no vertical back face). Each surface is solved with the current
    inputs, mesh and grid option; it is skipped when the grid option cannot analyse it (the default grid refuses it,
    or its mass covers less than ``coverage`` of the default grid's span from x_c - R to the toe-level crossing, so
    that few slices or columns resolve it; with the base-circle option, the ground is not above it over one
    interval), and a result that is not ACCEPTED is withheld. 2D: a grid of starting surfaces (front exit x back exit
    x curvature), then Nelder-Mead from the ``starts`` best grid surfaces that are not neighbours, each restarted from
    its minimum with half the initial step until it stops improving. 3D: Nelder-Mead from the 2D minima. Then each
    critical surface is probed one step 2 tol along every parameter (to name the limits it lies on) and re-solved on a
    finer mesh. A local search: the lowest FoS among the surfaces tried, not a proof of the global minimum.

    ``progress(done, total, label)`` is called as the search advances. Same as criticalSearch in webapp/engine.js;
    the result has the same keys.
    """
    import time

    o = {**SEARCH, **opts}
    dims = list(dims)
    t0 = time.perf_counter()
    base = bool(inp.base_circle)
    na, nb, nt = o["grid_base"] if base else o["grid"]
    lo0 = -o["front"] if base else 0
    hi1 = 1 + o["back"]
    sp = [(1 - lo0) / na, hi1 / nb, 1 / nt]
    step = [v / 2 for v in sp]
    local = o["max_eval"] * (1 + o["restarts"])
    budget = {2: na * nb * nt + o["starts"] * local + 8, 3: o["starts"] * local + 8}
    total = sum(budget[d] for d in dims)
    state = {"done": 0, "stage": ""}

    def say(label):
        if progress:
            progress(min(state["done"], total), total, label)

    faces = vertical_faces(inp)
    out = {"settings": {"back": o["back"], "front": o["front"], "grid": [na, nb, nt], "starts": o["starts"],
                        "tol": o["tol"], "max_eval": o["max_eval"], "restarts": o["restarts"],
                        "coverage": o["coverage"], "dims": list(dims), "base_circle": base},
           "domain": {"front": [lo0, 1], "back": [0, hi1], "curvature": [0, 1]},
           "prescribed": {"xc": inp.xc, "zc": inp.zc, "R": inp.R, "vertical_back": faces["back"],
                          "vertical_front": faces["front"]},
           "r2": None, "r3": None, "ms": 0}

    def solve_at(dim, **kw):
        try:
            r = (run_3d if dim == 3 else run_2d)(inp.with_(**kw))
            return {"fos": r.fos, "status": r.status, "reason": "" if r.fos is not None else r.reason}
        except Exception as e:                    # as the web app's jobs, which report any error as an input error
            return {"fos": None, "status": "INPUT ERROR", "reason": str(e)}

    for dim in dims:
        pr = solve_at(dim)
        out["prescribed"][f"fos{dim}"] = pr["fos"]
        out["prescribed"][f"status{dim}"] = pr["status"]
        out["prescribed"][f"reason{dim}"] = pr["reason"]
    stage_end = {2: budget[2], 3: (budget[2] if 2 in dims else 0) + budget[3]}
    seeds = None
    inf = math.inf
    for dim in dims:
        td = time.perf_counter()
        cache = {}
        counts = {"solved": 0, "withheld": 0, "skipped": 0}
        why = {"coverage": 0, "masses": 0, "grid": 0}
        cloud = {"xc": [], "zc": [], "R": [], "fos": []}

        def evaluate(p, dim=dim, counts=counts, why=why, cloud=cloud):
            a, b, t = p
            if not (lo0 <= a < 1 and max(a, 0) < b <= hi1 and 0 < t <= 1):
                return {"v": inf, "kind": "outside"}
            s = search_surface(inp, p)
            kw = {"xc": s["xc"], "zc": s["zc"], "R": s["R"]}
            skip = None
            if base:
                if not single_mass(inp.with_(**kw)):
                    skip = "masses"
            elif s["R"] > s["zc"]:
                span = (s["xc"] - math.sqrt(s["R"] * s["R"] - s["zc"] * s["zc"])) - (s["xc"] - s["R"])
                if s["front"][0] - s["back"][0] < o["coverage"] * span:
                    skip = "coverage"
            if skip:
                e = {"v": inf, "kind": "skipped", "why": skip}
            else:
                r = solve_at(dim, **kw)
                if r["status"] == "INPUT ERROR":
                    e = {"v": inf, "kind": "skipped", "why": "grid"}
                elif r["fos"] is None:
                    e = {"v": inf, "kind": "withheld"}
                else:
                    e = {"v": r["fos"], "kind": "solved"}
            counts[e["kind"]] += 1
            if e["kind"] == "skipped":
                why[e["why"]] += 1
            for k in ("xc", "zc", "R"):
                cloud[k].append(s[k])
            cloud["fos"].append(e["v"] if e["kind"] == "solved" else None)
            state["done"] += 1
            if progress and state["done"] % 5 == 0:
                say(state["stage"])
            return e

        def look(p, cache=cache, evaluate=evaluate):
            key = tuple(p)
            e = cache.get(key)
            if e is None:
                e = evaluate(p)
                cache[key] = e
            return e

        def f(p, look=look):
            return look(p)["v"]

        runs = []
        res = {"dim": dim, "fos": None, "message": "", "runs": runs, "counts": counts, "skipped_why": why,
               "evaluations": 0, "edges": [], "check": None, "cloud": cloud}
        starts = []
        if dim == 2:
            state["stage"] = "2D: grid of starting surfaces"
            cells = []
            for i in range(na):
                for j in range(1, nb + 1):
                    for k in range(nt):
                        p = [lo0 + sp[0] * i, sp[1] * j, sp[2] * (k + 0.5)]
                        cells.append({"p": p, "ijk": (i, j, k), "v": f(p)})
            for i in sorted(range(len(cells)), key=lambda q: (cells[q]["v"], q)):
                if not cells[i]["v"] < inf or len(starts) >= o["starts"]:
                    break
                ijk = cells[i]["ijk"]
                if any(max(abs(c["ijk"][0] - ijk[0]), abs(c["ijk"][1] - ijk[1]), abs(c["ijk"][2] - ijk[2])) < 2
                       for c in starts):
                    continue
                starts.append({"p": cells[i]["p"], "ijk": ijk, "v": cells[i]["v"]})
            state["done"] = max(state["done"], na * nb * nt)
            if not starts:
                res["message"] = "No surface of the starting grid gave an accepted FoS."
        elif seeds:
            starts = [{"p": p, "v": f(p)} for p in seeds]
        else:
            res["message"] = "The 2D search found no critical circle to start the 3D search from."
        best = None
        for q, c in enumerate(starts):
            state["stage"] = f"{dim}D: Nelder–Mead from start {q + 1} of {len(starts)}"
            say(state["stage"])
            start_done = state["done"]
            x, fx = _nelder_mead(f, c["p"], step, o["tol"], o["max_eval"])
            for _ in range(o["restarts"]):
                x2, f2 = _nelder_mead(f, x, [v / 2 for v in step], o["tol"], o["max_eval"])
                if not f2 < fx:
                    break
                x, fx = x2, f2
            state["done"] = max(state["done"], start_done + local)
            runs.append({"p0": c["p"], "f0": c["v"] if c["v"] < inf else None, "p": x, "f": fx if fx < inf else None})
            if fx < inf and (best is None or fx < best[1]):
                best = (x, fx)
        if dim == 2:
            seeds = [r["p"] for r in runs if r["f"] is not None]
        if best is not None:
            state["stage"] = f"{dim}D: limits and mesh check"
            say(state["stage"])
            names = ("front", "back", "curvature")
            top = best
            for i in range(3):
                for sg in (-1, 1):
                    q = list(best[0])
                    q[i] = best[0][i] + sg * 2 * o["tol"]
                    e = look(q)
                    if e["kind"] == "solved":
                        if e["v"] < top[1]:
                            top = (q, e["v"])
                        continue
                    w = e["why"] if e["kind"] == "skipped" else e["kind"]
                    if e["kind"] == "outside":
                        w = (("front limit" if base else "toe") if q[0] < lo0 else "crest" if q[0] >= 1
                             else "exits meet" if q[1] <= max(q[0], 0) else "back limit" if q[1] > hi1
                             else "vertical back" if q[2] > 1 else "flat")
                    res["edges"].append({"param": names[i], "side": sg, "why": w})
            best = top
            s = search_surface(inp, best[0])
            res.update({"fos": best[1], "p": best[0], "xc": s["xc"], "zc": s["zc"], "R": s["R"], "front": s["front"],
                        "back": s["back"]})
            if not base and s["R"] > s["zc"]:
                res["coverage"] = (s["front"][0] - s["back"][0]) / (
                    (s["xc"] - math.sqrt(s["R"] * s["R"] - s["zc"] * s["zc"])) - (s["xc"] - s["R"]))
            n = min(2 * inp.n3, 320) if dim == 3 else min(4 * inp.n2, 10000)
            kw = {"xc": s["xc"], "zc": s["zc"], "R": s["R"], ("n3" if dim == 3 else "n2"): n}
            ny = None
            if dim == 3 and inp.n3y is not None:     # columns across y in the same ratio, at most 320 as along x
                ny = min(columns_across(inp, n), 320)
                kw["n3y"] = ny
            ck = solve_at(dim, **kw)
            res["check"] = {"n": n, "ny": ny if ny is not None else (n if dim == 3 else 1), "fos": ck["fos"],
                            "status": ck["status"], "reason": ck["reason"],
                            "change": (ck["fos"] - best[1]) / ck["fos"] if ck["fos"] is not None else None}
        res["evaluations"] = counts["solved"] + counts["withheld"] + counts["skipped"]
        res["ms"] = (time.perf_counter() - td) * 1000.0
        out[f"r{dim}"] = res
        state["done"] = stage_end[dim]
        say(f"{dim}D: done")
    out["ms"] = (time.perf_counter() - t0) * 1000.0
    return out


# ---------------------------------------------------------------------------------------------------------------
# studies (additions): EJC against GHB, parameter sweeps, back-analysis, Monte Carlo
# ---------------------------------------------------------------------------------------------------------------
def _attempt(inp: Inputs, dim: int) -> dict:
    """FoS, status and reason of one run; an input error is reported as status INPUT ERROR."""
    try:
        r = run_3d(inp) if dim == 3 else run_2d(inp)
        return {"fos": r.fos, "status": r.status, "reason": "" if r.fos is not None else r.reason}
    except ValueError as e:
        return {"fos": None, "status": "INPUT ERROR", "reason": str(e)}


COMPARE = (("ejc", "ejc_moscarda"), ("ghb", "ghb"))


def compare_models(inp: Inputs) -> dict:
    """EJC (Moscarda closure, at the given nu) and GHB on the same surface, mesh and inputs, 3D and 2D: FoS, status and
    reason of each, the 3D/2D ratio of each model and the EJC/GHB ratio in 3D and 2D (None where a FoS is withheld).
    Same as compareModels in webapp/engine.js."""
    out = {"models": {}, "ratio": {}}
    for key, model in COMPARE:
        b = inp.with_(model=model)
        row = {}
        for dim in (3, 2):
            a = _attempt(b, dim)
            row[f"fos{dim}"], row[f"status{dim}"], row[f"reason{dim}"] = a["fos"], a["status"], a["reason"]
        row["ratio"] = row["fos3"] / row["fos2"] if row["fos3"] is not None and row["fos2"] is not None else None
        out["models"][key] = row
    e, g = out["models"]["ejc"], out["models"]["ghb"]
    for dim in (3, 2):
        out["ratio"][f"ejc_ghb{dim}"] = (e[f"fos{dim}"] / g[f"fos{dim}"]
                                         if e[f"fos{dim}"] is not None and g[f"fos{dim}"] is not None else None)
    return out


def parameter_sweep(inp: Inputs, name: str, values, models=("ejc_moscarda", "ghb"), dims=(3, 2), progress=None) -> list:
    """The FoS for each value of the input ``name`` (any numeric input), for each strength model and dimension:
    rows {"value", "<model>": {"fos3", "status3", "reason3", "fos2", ...}}. Same as parameterSweep."""
    if not hasattr(inp, name):
        raise ValueError(f"unknown input {name!r}")
    rows = []
    for q, v in enumerate(values):
        row = {"value": float(v)}
        for model in models:
            b = inp.with_(**{name: float(v), "model": model})
            cell = {}
            for dim in dims:
                a = _attempt(b, dim)
                cell[f"fos{dim}"], cell[f"status{dim}"], cell[f"reason{dim}"] = a["fos"], a["status"], a["reason"]
            row[model] = cell
        rows.append(row)
        if progress:
            progress(q + 1, len(values), f"{name} = {float(v):g}")
    return rows


def back_analysis(inp: Inputs, name: str, lo: float, hi: float, target: float = 1.0, dim: int = 2,
                  rel_tol: float = 1e-6, max_steps: int = 60, progress=None) -> dict:
    """The value of the input ``name`` in [lo, hi] at which the FoS (dim 3 or 2) equals ``target``, by bisection.
    F - target must change sign between lo and hi, with an accepted FoS at both; each step solves at the midpoint and
    keeps the half with the sign change, until the bracket is narrower than rel_tol (hi - lo). The value reported is
    the linear interpolation of F across the final bracket, solved once more. If the FoS steps (a column entering the
    mass), the result is the step. Same as backAnalysis in webapp/engine.js."""
    if not hasattr(inp, name):
        raise ValueError(f"unknown input {name!r}")
    if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
        raise ValueError("require lo < hi")
    steps = []

    def f(v):
        a = _attempt(inp.with_(**{name: float(v)}), dim)
        steps.append({"value": float(v), "fos": a["fos"], "status": a["status"]})
        if progress:
            progress(len(steps), max_steps + 3, f"{name} = {float(v):.6g}")
        return a

    out = {"name": name, "target": target, "dim": dim, "value": None, "fos": None, "bracket": [lo, hi],
           "steps": steps, "message": ""}
    fa, fb = f(lo), f(hi)
    for v, a in ((lo, fa), (hi, fb)):
        if a["fos"] is None:
            out["message"] = f"no accepted FoS at {name} = {v:g} ({a['status']}): narrow the range"
            return out
    ga, gb = fa["fos"] - target, fb["fos"] - target
    if ga * gb > 0:
        out["message"] = (f"the FoS does not cross {target:g} between {name} = {lo:g} (F = {fa['fos']:.6f}) and "
                          f"{hi:g} (F = {fb['fos']:.6f}): widen the range")
        return out
    a, b, width = lo, hi, rel_tol * (hi - lo)
    for _ in range(max_steps):
        if b - a <= width:
            break
        mid = 0.5 * (a + b)
        fm = f(mid)
        if fm["fos"] is None:
            out["message"] = f"no accepted FoS at {name} = {mid:.6g} ({fm['status']})"
            out["bracket"] = [a, b]
            return out
        gm = fm["fos"] - target
        if ga * gm <= 0:
            b, gb = mid, gm
        else:
            a, ga = mid, gm
    v = a + (0.0 - ga) * (b - a) / (gb - ga) if gb != ga else 0.5 * (a + b)
    fv = f(v)
    out.update(value=v, fos=fv["fos"], bracket=[a, b])
    if fv["fos"] is None:
        out["message"] = f"no accepted FoS at the result ({fv['status']})"
    return out


def mulberry32(seed: int):
    """Seeded uniform random numbers in [0, 1) (Mulberry32, 32-bit integer arithmetic), the same stream as mulberry32
    in webapp/engine.js."""
    a = int(seed) & 0xFFFFFFFF

    def imul(x, y):
        return (x * y) & 0xFFFFFFFF

    def nxt():
        nonlocal a
        a = (a + 0x6D2B79F5) & 0xFFFFFFFF
        t = imul(a ^ (a >> 15), 1 | a)
        t = ((t + imul(t ^ (t >> 7), 61 | t)) & 0xFFFFFFFF) ^ t
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296.0
    return nxt


def monte_carlo(inp: Inputs, variables, n: int, seed: int = 1, models=None, dims=(3, 2), progress=None) -> dict:
    """Monte Carlo analysis (addition): ``variables`` = [(name, mean, sd, lo, hi), ...], each sampled from a normal
    distribution truncated to [lo, hi] (Box-Muller from the seeded Mulberry32 stream, rejection outside the limits; the
    variables in the given order, sample by sample), and the FoS of each sample for each model and dimension. Returns
    the samples and, per model and dimension, the accepted count, mean, standard deviation, reliability index
    (mean - 1)/sd, probability of failure P(F < 1) among accepted results with its standard error, and the 5 %, 50 %
    and 95 % values. Same as monteCarlo in webapp/engine.js."""
    if not (isinstance(n, (int, np.integer)) and n >= 1):
        raise ValueError("the number of samples must be a positive integer")
    models = models or ((inp.model,) if inp.model != "ejc_original" else ("ejc_moscarda",))
    rnd = mulberry32(seed)

    def normal():
        u1 = 1.0 - rnd()
        u2 = rnd()
        return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)

    for name, mean, sd, lo, hi in variables:
        if not hasattr(inp, name):
            raise ValueError(f"unknown input {name!r}")
        if not (sd >= 0 and lo <= mean <= hi):
            raise ValueError(f"{name}: require sd >= 0 and lo <= mean <= hi")
    samples = []
    for k in range(n):
        vals = {}
        for name, mean, sd, lo, hi in variables:
            v = mean
            if sd > 0:
                for _ in range(1000):
                    v = mean + sd * normal()
                    if lo <= v <= hi:
                        break
                else:
                    raise ValueError(f"{name}: no sample inside [{lo:g}, {hi:g}] in 1000 draws")
            vals[name] = v
        b = inp.with_(**vals)
        row = {"values": vals}
        for model in models:
            for dim in dims:
                a = _attempt(b.with_(model=model), dim)
                row[f"{model}{dim}"] = a["fos"]
        samples.append(row)
        if progress:
            progress(k + 1, n, f"sample {k + 1} of {n}")
    stats = {}
    for model in models:
        for dim in dims:
            f = np.array([s[f"{model}{dim}"] for s in samples if s[f"{model}{dim}"] is not None], dtype=float)
            st = {"n": n, "accepted": int(len(f)), "mean": None, "sd": None, "beta": None, "pof": None, "pof_se": None,
                  "p05": None, "p50": None, "p95": None}
            if len(f):
                st["mean"] = float(np.mean(f))
                st["sd"] = float(np.std(f, ddof=1)) if len(f) > 1 else 0.0
                st["beta"] = (st["mean"] - 1.0) / st["sd"] if st["sd"] > 0 else None
                p = float(np.mean(f < 1.0))
                st["pof"], st["pof_se"] = p, math.sqrt(p * (1 - p) / len(f))
                fs = np.sort(f)
                for key, qq in (("p05", 0.05), ("p50", 0.5), ("p95", 0.95)):
                    st[key] = float(fs[min(len(fs) - 1, int(math.floor(qq * len(fs))))])
            stats[f"{model}{dim}"] = st
    return {"samples": samples, "stats": stats, "models": list(models), "dims": list(dims), "seed": seed,
            "variables": [list(v) for v in variables]}
