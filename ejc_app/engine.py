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
from dataclasses import dataclass, field, replace
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
    xc, zc, R = inp.xc, inp.zc, inp.R

    def above(x):
        return bool(ground_z(x, inp.H, inp.beta) - sphere_z(x, 0.0, xc, zc, R) > 0)

    xs = np.linspace(xc - R, xc + R, n)
    with np.errstate(invalid="ignore"):
        d = ground_z(xs, inp.H, inp.beta) - sphere_z(xs, 0.0, xc, zc, R)
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


def vertical_faces(inp: Inputs) -> dict:
    """Heights of the vertical faces that bound the sliding mass where the slip circle on y = 0 ends below the ground:
    at its back extreme x_c - R (z = z_c), and with the base-circle grid also at its front extreme x_c + R (the default
    grid stops at the toe). No force acts on such a face: at the back it works as an open tension crack. Heights below
    1e-9 H (rounding) count as zero. Same as verticalFaces in webapp/engine.js."""
    xc, zc, R, tol = inp.xc, inp.zc, inp.R, 1e-9 * inp.H
    back = float(ground_z(xc - R, inp.H, inp.beta)) - zc if inp.base_circle or xc - R <= 0 else 0.0
    front = float(ground_z(xc + R, inp.H, inp.beta)) - zc if inp.base_circle else 0.0
    return {"back": back if back > tol else 0.0, "front": front if front > tol else 0.0}


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
    xc, zc, R = inp.xc, inp.zc, inp.R
    if inp.base_circle:
        x1, x0 = slip_extent(inp)
    else:
        if not (R > zc):
            raise ValueError("require R > z_c (the sphere must cut the toe level); for a face circle, with its lowest point "
                             "above the toe level, use the base-circle option")
        x0 = xc - math.sqrt(R * R - zc * zc)
        x1 = xc - R
        if x0 < -1e-5 * R:
            raise ValueError(f"the circle crosses the toe level {-x0:.4g} m behind the toe, so the default grid would cut "
                             f"off the toe region; use the base-circle option")
    dx = (x0 - x1) / n
    if dim == 3:
        dy = R / n
        ii, jj = np.meshgrid(np.arange(1, n + 1), np.arange(1, n + 1), indexing="ij")
        ii, jj = ii.ravel(), jj.ravel()
        y = (jj - 0.5) * dy
    else:
        dy = 1.0
        ii = np.arange(1, n + 1)
        jj = np.ones(n, dtype=int)
        y = np.zeros(n)
    x = x0 - (ii - 0.5) * dx
    q = R * R - (x - xc) ** 2 - y ** 2
    keep = (q > 0) if inp.base_circle else (q > 0) & (x <= 0)
    ii, jj, x, y, q = ii[keep], jj[keep], x[keep], y[keep], q[keep]
    zb = zc - np.sqrt(q)
    zt = ground_z(x, inp.H, inp.beta)
    h = zt - zb
    keep = h > 0
    ii, jj, x, y, q, zb, zt, h = (v[keep] for v in (ii, jj, x, y, q, zb, zt, h))
    ah = np.full(len(x), dx * dy)
    g = np.sqrt(q) / R
    rho = np.sqrt(R * R - y * y) if dim == 3 else np.full(len(x), R)
    kap = (xc - x) / rho
    if inp.ru is not None:
        zw = np.full(len(x), np.nan)
        u = inp.ru * inp.gamma * h
    elif inp.wet:
        zw = water_z(x, inp.Hw, inp.beta_w)
        u = GAMMA_W * np.maximum(0.0, zw - zb)
    else:
        zw = np.full(len(x), np.nan)
        u = np.zeros(len(x))
    return Mesh(dim, n, x0, x1, dx, dy, ii, jj, x, y, zb, zt, h, ah, g, ah / g, rho, kap,
                inp.gamma * h * ah, u, zw)


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


def _global_f(m: Mesh, tau_e, xc: float, form: str) -> float:
    """F_H = sum(tau_s A)/sum(W kappa) (eq. 4.14, governing); F_M = sum(tau_s A rho)/sum(W (x_c - x)) (eq. 4.15)."""
    if form == "FM":
        return float(np.sum(tau_e * m.area * m.rho) / np.sum(m.w * (xc - m.x)))
    return float(np.sum(tau_e * m.area) / np.sum(m.w * m.kap))


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
    drive = np.sum(m.w * (xc - m.x)) if form == "FM" else np.sum(m.w * m.kap)
    if not drive > 0:
        res.reason = "Non-positive driving sum."
        return res
    sn = (gamma * m.h - m.u) * m.g ** 2
    limit = projection_limit(st) if safeguard == "project" else None
    if limit is not None:
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
        stg = solve_strength(sn, st)
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
            return _final_checks(res, m, st, xc, F, sn, tens, dF, move)
    res.status = "NOT CONVERGED"
    res.reason = f"No convergence in {max_iter} iterations; no accepted FoS."
    res.sn = sn
    return res


def _final_checks(res: LemResult, m: Mesh, st: Strength, xc: float, F: float, sn: np.ndarray,
                  tens_used: np.ndarray, dF: float, move: float) -> LemResult:
    """Independent check of the converged state: fresh strength at the final sigma_n, then the local (vertical)
    and global residuals at the reported F (brief §1.4 step 4)."""
    stg = solve_strength(sn, st)
    t = np.tan(stg.phi)
    tens = tens_used | (sn < 0) if res.tension != "none" else tens_used
    c_e, t_e, tau_e = _sensitivity(stg.c, t, stg.tau, sn, res.tension, tens_used)
    ma = m.g + m.kap * t_e / F
    with np.errstate(invalid="ignore"):
        local = np.abs(m.w - (sn + m.u) * m.ah - tau_e * m.area * m.kap / F) / m.w
    f_glob = _global_f(m, tau_e, xc, res.form)
    glob = abs(f_glob - F) / F
    bishop = float(np.sum((c_e * m.ah + (m.w - m.u * m.ah) * t_e) / ma) / np.sum(m.w * m.kap))
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
        checks.append(Check("Min sigma_n above the projection limit (kPa)", float(np.min(sn) - res.projection_limit),
                            "> 0", bool(np.all(sn > res.projection_limit))))
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
    return solve(build_mesh(inp, 3, inp.n3), inp.strength(), inp.xc, inp.gamma, inp.form, inp.tension, inp.iteration,
                 inp.safeguard)


def run_2d(inp: Inputs) -> LemResult:
    """Bishop on the central section (brief §1.5, eq. 4.16); F_M = F_H there (rho = R), so F_H is used."""
    return solve(build_mesh(inp, 2, inp.n2), inp.strength(), inp.xc, inp.gamma, "FH", inp.tension, inp.iteration,
                 inp.safeguard)


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
           "A_h_m2": m.ah, "A_m2": m.area, "g_z": m.g, "rho_m": m.rho, "kappa": m.kap, "W_kN": m.w, "u_kPa": m.u}
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
    Same as massForces in webapp/engine.js.
    """
    m = r.mesh
    F = r.fos if r.fos is not None else r.f_trial
    if r.sn is None or r.tau_used is None or F is None:
        raise ValueError("no converged state for this run")
    mult = 2.0 if m.dim == 3 else 1.0
    nx, ny, nz = (inp.xc - m.x) / inp.R, -m.y / inp.R, m.g.copy()
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
    return {"F": F, "dim": m.dim, "mult": mult, "count": m.count, "form": r.form, "weight": -totals["W"]["z"],
            "totals": totals, "moments": {"W": m_w, "N": 0.0, "U": 0.0, "S": m_s, "net": m_w + m_s},
            "method": {"form": r.form, "lhs": lhs, "rhs": rhs, "residual": (lhs - rhs) / rhs},
            "vertical_max_rel": float(np.max(np.abs(vert) / W)) if m.count else 0.0,
            "el": {"x": m.x, "y": m.y, "zb": m.zb, "zm": (m.zb + m.zt) / 2, "W": W, "N": N, "U": U, "S": S, "nx": nx,
                   "ny": ny, "nz": nz, "sx": sx, "sz": sz}}


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
    """The surface of search parameters p = (a, b, t): front exit at s = aH, back exit at s = bH (ground_at),
    curvature t (circle_from_exits)."""
    H, tan_b = inp.H, math.tan(math.radians(inp.beta))
    front, back = ground_at(p[0] * H, H, tan_b), ground_at(p[1] * H, H, tan_b)
    return {**circle_from_exits(front, back, p[2]), "front": front, "back": back}


def single_mass(inp: Inputs) -> bool:
    """True if the ground is above the circle on y = 0 over one interval that ends inside (x_c - R, x_c + R), sampled
    as slip_extent samples it: one sliding mass, with no vertical face at the circle's extremes (base-circle grid)."""
    xs = np.linspace(inp.xc - inp.R, inp.xc + inp.R, 20001)
    with np.errstate(invalid="ignore"):
        d = ground_z(xs, inp.H, inp.beta) - sphere_z(xs, 0.0, inp.xc, inp.zc, inp.R)
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
    H, base = inp.H, bool(inp.base_circle)
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
            ck = solve_at(dim, **kw)
            res["check"] = {"n": n, "fos": ck["fos"], "status": ck["status"], "reason": ck["reason"],
                            "change": (ck["fos"] - best[1]) / ck["fos"] if ck["fos"] is not None else None}
        res["evaluations"] = counts["solved"] + counts["withheld"] + counts["skipped"]
        res["ms"] = (time.perf_counter() - td) * 1000.0
        out[f"r{dim}"] = res
        state["done"] = stage_end[dim]
        say(f"{dim}D: done")
    out["ms"] = (time.perf_counter() - t0) * 1000.0
    return out
