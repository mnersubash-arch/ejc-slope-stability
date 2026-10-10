"""Instantaneous friction angle from eq. 4.21 and the acceptance checks of methodology Section 4.2.

Port of ``reference/EJC_Eq421.bas`` (library revision 6). The VBA is the specification; where the
Python mirror in ``reference/code`` differs, this module follows the VBA (see docs/SPEC.md, section 0).

Units: sigma_n' and c'i in kPa; sigma_ci in MPa; angles in degrees at the public interface, radians
internally. Closure: sigma2' = p*sigma1' + q*sigma3'.

Errors are *values*, not exceptions: ``ejc_phi`` returns an :class:`ErrorCode` (``#N/A``, ``#VALUE!``,
``#NUM!``) exactly where the VBA returns ``CVErr``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Union

import numpy as np

# EJC_Eq421.bas line 21: the strength library uses a 15-digit pi (not math.pi). Kept for fidelity.
PI_ = 3.14159265358979
HI_DEG = 89.9999
NSCAN = 1800


class ErrorCode(str, Enum):
    """Controlled error results (the VBA returns CVErr values; compare with ``==``)."""

    NA = "#N/A"        # no (further) root bracket
    VALUE = "#VALUE!"  # invalid input
    NUM = "#NUM!"      # domain failure / not converged

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"ErrorCode({self.value})"


PhiResult = Union[float, ErrorCode]


@dataclass(frozen=True)
class StateCheck:
    """Result of ``StateChecks`` (VBA 1x8 array). Stresses in MPa, angles in rad.

    ``angular`` and ``ejc`` are ``None`` where the VBA holds ``#NUM!``.
    """

    angular: Optional[float]   # phi - asin((M-1)/(M+1)), unclipped [rad]
    ejc: Optional[float]       # original, unsquared EJC residual [MPa]
    closure: Optional[float]   # sigma2 - (p sigma1 + q sigma3) [MPa]
    ordering: int              # 1 = sigma1 >= sigma2 >= sigma3 (1e-9 slack)
    s1: Optional[float]        # sigma1' [MPa]
    s3: Optional[float]        # sigma3' [MPa]
    domain: int                # 1 = derivative/radicand conditions valid
    branch: int                # 1 = adopted branch (0 <= phi < 90 deg, tangent angle >= -tol)


def closure_constants(p: float, q: float) -> tuple[float, float, float]:
    """P = p^2 - p + 1, R = 2pq - p - q - 1, S = q^2 - q + 1  (VBA GState line 37)."""
    return p * p - p + 1, 2 * p * q - p - q - 1, q * q - q + 1


def valid_inputs(sci: float, mb: float, s: float, a: float, tol_a: float) -> bool:
    """VBA ``ValidInputs`` (line 26), plus a finiteness test the VBA lacks (NaN would compare False anyway)."""
    vals = (sci, mb, s, a, tol_a)
    if any(isinstance(v, bool) or not isinstance(v, (int, float, np.floating, np.integer)) for v in vals):
        return False
    if not all(math.isfinite(v) for v in vals):
        return False
    return sci > 0 and mb > 0 and s >= 0 and 0 < a <= 1 and 0 < tol_a <= 1e-6


# --------------------------------------------------------------------------------------------------
# g(phi) and its derivative
# --------------------------------------------------------------------------------------------------
def gstate(phi: float, sn: float, sci: float, mb: float, s: float, a: float, p: float, q: float
           ) -> Optional[tuple[float, float, float]]:
    """``GState``: (g, dg/dphi, M) with g = sin(phi) - 1 + 2D/(N+D) = sin(phi) - (M-1)/(M+1); ``None`` if undefined.

    ``sn`` is in MPa here (the public functions convert from kPa).
    """
    cP, cR, cS = closure_constants(p, q)
    sin_p, cos_p = math.sin(phi), math.cos(phi)
    b_term = mb * sn / sci + s
    if b_term <= 0:
        return None
    h = (b_term / (1 + sin_p / a)) ** a
    dh = -h * cos_p / (1 + sin_p / a)
    sig3 = sn - sci / 2 * h * (1 - sin_p)
    dsig3 = -sci / 2 * (dh * (1 - sin_p) - h * cos_p)
    sig1 = sn + sci / 2 * h * (1 + sin_p)
    dsig1 = sci / 2 * (dh * (1 + sin_p) + h * cos_p)
    hb_q = mb * sig3 / sci + s
    dhb_q = mb * dsig3 / sci
    if hb_q <= 0:                     # beyond the Hoek-Brown tensile limit
        return None
    n_t = 2 * a * mb * sci * hb_q ** (2 * a - 1) - cR * sig1 - 2 * cS * sig3
    dn_t = (2 * a * mb * sci * (2 * a - 1) * hb_q ** (2 * a - 2) * dhb_q - cR * dsig1 - 2 * cS * dsig3)
    d_t = 2 * cP * sig1 + cR * sig3
    dd_t = 2 * cP * dsig1 + cR * dsig3
    if n_t + d_t == 0 or d_t == 0:
        return None
    g = sin_p - 1 + 2 * d_t / (n_t + d_t)
    dg = cos_p - 2 * (dn_t * d_t - n_t * dd_t) / (n_t + d_t) ** 2
    return g, dg, n_t / d_t


def ang_res(phi: float, sn: float, sci: float, mb: float, s: float, a: float, p: float, q: float) -> Optional[float]:
    """``AngRes``: UNCLIPPED residual phi - asin((M-1)/(M+1)); ``None`` if M+1 = 0 or |argument| > 1."""
    r = gstate(phi, sn, sci, mb, s, a, p, q)
    if r is None:
        return None
    m = r[2]
    if m + 1 == 0:
        return None
    z = (m - 1) / (m + 1)
    if abs(z) > 1:
        return None
    return phi - math.asin(z)


def phi_b(sn: float, sci: float, mb: float, s: float, a: float) -> Optional[float]:
    """``PhiB``: angle above which the eq. 4.21 state is inside the tensile limit (Q > 0). ``None`` = not ok."""
    t_lim = s * sci / mb
    k = sn + t_lim
    if k <= 0:
        return None
    z = 2 * k / (sci * (mb * k / sci) ** a)
    if z >= 1:
        return 0.0
    pv = (1 + z * z) - math.sqrt((1 + z * z) ** 2 - (1 - z * z))
    for _ in range(8):
        f = (1 - pv) * (1 + pv / a) ** (-a) - z
        fp = -(1 + pv / a) ** (-a) * (1 + (1 - pv) / (1 + pv / a))
        if fp == 0:
            break
        pv = pv - f / fp
        if pv < 0:
            pv = 0.0
        if pv > 1 - 1e-15:
            pv = 1 - 1e-15
    return math.asin(pv)


# --------------------------------------------------------------------------------------------------
# bracket scan and root
# --------------------------------------------------------------------------------------------------
def brackets(sn: float, sci: float, mb: float, s: float, a: float, p: float, q: float
             ) -> list[tuple[float, float]]:
    """``FindBracket``: candidate brackets observed by the finite scan, after the preliminary pole filter.

    ``sn`` in MPa. The count is not proof that every root has been found.
    """
    out: list[tuple[float, float]] = []
    pb = phi_b(sn, sci, mb, s, a)
    if pb is None:
        return out
    # A point only 1e-12 rad above the tensile boundary can create a spurious first bracket,
    # so the scan starts 0.001 degree inside the domain (VBA lines 99-101).
    x0 = pb + 0.001 * PI_ / 180
    if x0 < 0.01 * PI_ / 180:
        x0 = 0.01 * PI_ / 180
    x1 = HI_DEG * PI_ / 180
    if x0 >= x1:
        return out
    have_prev = False
    x_prev = g_prev = 0.0
    for k in range(NSCAN + 1):
        x = x0 + (x1 - x0) * k / NSCAN
        st = gstate(x, sn, sci, mb, s, a, p, q)
        if st is not None:
            g_v = st[0]
            if have_prev and g_prev * g_v <= 0:
                lo, hi, glo, ghi = x_prev, x, g_prev, g_v
                is_pole = False
                for _ in range(30):                      # bisect; |g| that stays large means a pole
                    xm = (lo + hi) / 2
                    sm = gstate(xm, sn, sci, mb, s, a, p, q)
                    if sm is None:
                        is_pole = True
                        break
                    if glo * sm[0] <= 0:
                        hi, ghi = xm, sm[0]
                    else:
                        lo, glo = xm, sm[0]
                if not is_pole and abs(glo) > 0.001 and abs(ghi) > 0.001:
                    is_pole = True
                if not is_pole:
                    out.append((x_prev, x))
            have_prev, x_prev, g_prev = True, x, g_v
        else:
            have_prev = False
    return out


def ejc_phi(sn_kpa: float, sci: float, mb: float, s: float, a: float, p: float, q: float,
            tol_a: float = 1e-10, root_index: int = 1) -> PhiResult:
    """``EJC_PHI``: CANDIDATE root only (use :func:`ejc_accept` for acceptance). Degrees, or an :class:`ErrorCode`."""
    if (not valid_inputs(sci, mb, s, a, tol_a) or not math.isfinite(sn_kpa)
            or not (math.isfinite(p) and math.isfinite(q))
            or isinstance(root_index, bool) or not isinstance(root_index, (int, np.integer)) or root_index < 1):
        return ErrorCode.VALUE
    sn = sn_kpa / 1000.0
    br = brackets(sn, sci, mb, s, a, p, q)
    if len(br) < root_index:
        return ErrorCode.NA
    lo, hi = br[root_index - 1]
    st = gstate(lo, sn, sci, mb, s, a, p, q)
    if st is None:
        return ErrorCode.NUM
    glo = st[0]
    x = (lo + hi) / 2
    for _ in range(200):
        r = gstate(x, sn, sci, mb, s, a, p, q)
        if r is None:
            return ErrorCode.NUM
        g_v, dg_v, _m = r
        if glo * g_v <= 0:
            hi = x
        else:
            lo, glo = x, g_v
        xn = x - g_v / dg_v if dg_v != 0 else (lo + hi) / 2
        if xn <= lo or xn >= hi:
            xn = (lo + hi) / 2
        step = abs(xn - x)
        x = xn
        rv = ang_res(x, sn, sci, mb, s, a, p, q)
        if rv is not None and abs(rv) <= tol_a and step <= tol_a:
            return x * 180 / PI_
    return ErrorCode.NUM                                   # iterations exhausted: not converged


def fastphi(sn_kpa: float, sci: float, mb: float, s: float, a: float, p: float, q: float,
            guess_deg: float, tol_a: float = 1e-12) -> PhiResult:
    """``EJC_FASTPHI``: the conventional-HB branch (p = 0, q = 1) has one tangent root, so the tensile-domain bound
    gives a single bracket and the previous phi seeds Newton. Other closures use the scanned FIRST candidate root.
    No root switching, no fallback strength.
    """
    if not valid_inputs(sci, mb, s, a, tol_a):
        return ErrorCode.VALUE
    if p != 0 or q != 1:
        return ejc_phi(sn_kpa, sci, mb, s, a, p, q, tol_a)
    sn = sn_kpa / 1000.0
    pb = phi_b(sn, sci, mb, s, a)
    if pb is None:
        return ErrorCode.NUM
    lo = pb + 1e-9
    hi = HI_DEG * PI_ / 180
    if lo >= hi:
        return ErrorCode.NUM
    sl = gstate(lo, sn, sci, mb, s, a, p, q)
    sh = gstate(hi, sn, sci, mb, s, a, p, q)
    if sl is None or sh is None or sl[0] * sh[0] > 0:
        return ejc_phi(sn_kpa, sci, mb, s, a, p, q, tol_a)         # FullScan
    glo = sl[0]
    x = (lo + hi) / 2
    g_rad = guess_deg * PI_ / 180
    if lo < g_rad < hi:
        x = g_rad
    for _ in range(200):
        st = gstate(x, sn, sci, mb, s, a, p, q)
        if st is None:
            break
        g_v, dg_v, _m = st
        res = ang_res(x, sn, sci, mb, s, a, p, q)
        if res is not None and dg_v != 0:
            if abs(res) <= tol_a and abs(g_v / dg_v) <= tol_a:
                return x * 180 / PI_
        if glo * g_v <= 0:
            hi = x
        else:
            lo, glo = x, g_v
        xn = x - g_v / dg_v if dg_v != 0 else (lo + hi) / 2
        if xn <= lo or xn >= hi:
            xn = (lo + hi) / 2
        step = abs(xn - x)
        x = xn
        res = ang_res(x, sn, sci, mb, s, a, p, q)
        if res is not None and abs(res) <= tol_a and step <= tol_a:
            return x * 180 / PI_
    return ejc_phi(sn_kpa, sci, mb, s, a, p, q, tol_a)


def n_brackets(sn_kpa: float, sci: float, mb: float, s: float, a: float, p: float, q: float) -> Union[int, ErrorCode]:
    """``EJC_NBRACKETS``: number of brackets OBSERVED by the finite scan."""
    if not valid_inputs(sci, mb, s, a, 1e-10):
        return ErrorCode.VALUE
    return len(brackets(sn_kpa / 1000.0, sci, mb, s, a, p, q))


# --------------------------------------------------------------------------------------------------
# state checks and acceptance
# --------------------------------------------------------------------------------------------------
def _domain_fail() -> StateCheck:
    """``DomainFail``: every entry #NUM! except ordering = domain = branch = 0."""
    return StateCheck(None, None, None, 0, None, None, 0, 0)


def state_checks(phi: float, sig1: float, sig2: float, sig3: float, sci: float, mb: float, s: float, a: float,
                 p: float, q: float, tol_a: float) -> StateCheck:
    """``StateChecks`` (VBA line 224): no clipping anywhere. Stresses in MPa, ``phi`` in rad."""
    cP, cR, cS = closure_constants(p, q)
    closure = sig2 - (p * sig1 + q * sig3)
    ordering = 1 if (sig1 >= sig2 - 1e-9 and sig2 >= sig3 - 1e-9) else 0
    angular = ejc = None
    domain = branch = 0
    hb_q = mb * sig3 / sci + s
    if hb_q > 0:
        den_m = 2 * cP * sig1 + cR * sig3
        if den_m != 0:
            m = (2 * a * mb * sci * hb_q ** (2 * a - 1) - cR * sig1 - 2 * cS * sig3) / den_m
            if m + 1 != 0:
                z = (m - 1) / (m + 1)
                rad = 4 * sci ** 2 * hb_q ** (2 * a) - 3 * (sig2 - sig3) ** 2
                if abs(z) <= 1 and rad >= 0:
                    t_ang = math.asin(z)
                    angular = phi - t_ang
                    ejc = sig1 - (sig2 + sig3 + math.sqrt(rad)) / 2
                    domain = 1
                    if phi >= 0 and phi < PI_ / 2 and t_ang >= -tol_a:
                        branch = 1
    return StateCheck(angular, ejc, closure, ordering, sig1, sig3, domain, branch)


def ejc_check(sn_kpa: float, phi_deg: float, sci: float, mb: float, s: float, a: float, p: float, q: float,
              tol_a: float = 1e-10) -> StateCheck:
    """``EJC_CHECK``: checks of the eq. 4.21 CANDIDATE state at (sigma_n', phi'i); controlled domain failure, never raises."""
    if not valid_inputs(sci, mb, s, a, tol_a) or not (math.isfinite(sn_kpa) and math.isfinite(phi_deg)
                                                       and math.isfinite(p) and math.isfinite(q)):
        return _domain_fail()
    if phi_deg < 0 or phi_deg >= 90:                              # adopted branch 0 <= phi < 90 deg
        return _domain_fail()
    sn = sn_kpa / 1000.0
    phi = phi_deg * PI_ / 180
    b_term = mb * sn / sci + s
    d_term = 1 + math.sin(phi) / a
    if b_term <= 0 or d_term <= 0:                                # guard before the fractional power
        return _domain_fail()
    h = (b_term / d_term) ** a
    sig3 = sn - sci / 2 * h * (1 - math.sin(phi))
    sig1 = sn + sci / 2 * h * (1 + math.sin(phi))
    return state_checks(phi, sig1, p * sig1 + q * sig3, sig3, sci, mb, s, a, p, q, tol_a)


def ejc_check_state(sn_kpa: float, c_kpa: float, phi_deg: float, sci: float, mb: float, s: float, a: float,
                    p: float, q: float, tol_a: float = 1e-10) -> StateCheck:
    """``EJC_CHECK_STATE``: checks of the FORCE-USED state (sigma_n', c'i, phi'i).

    tau = c + sigma_n tan(phi); centre = sigma_n + tau tan(phi); radius = tau / cos(phi); sigma1/3 = centre +/- radius.
    """
    if not valid_inputs(sci, mb, s, a, tol_a) or not all(math.isfinite(v) for v in (sn_kpa, c_kpa, phi_deg, p, q)):
        return _domain_fail()
    if phi_deg < 0 or phi_deg >= 90:
        return _domain_fail()
    phi = phi_deg * PI_ / 180
    tau = c_kpa + sn_kpa * math.tan(phi)
    cen = sn_kpa + tau * math.tan(phi)
    radius = tau / math.cos(phi)
    sig1 = (cen + radius) / 1000.0
    sig3 = (cen - radius) / 1000.0
    return state_checks(phi, sig1, p * sig1 + q * sig3, sig3, sci, mb, s, a, p, q, tol_a)


def accept_text(chk: Optional[StateCheck], tol_a: float, tol_e: float) -> str:
    """``AcceptText``: the first failing check, in the order INPUT, DOMAIN, BRANCH, TANGENCY, EJC, CLOSURE, ORDERING."""
    ok_tol = (not isinstance(tol_a, bool) and not isinstance(tol_e, bool)
              and math.isfinite(tol_a) and math.isfinite(tol_e)
              and 0 < tol_a <= 1e-6 and 0 < tol_e <= 1e-4)
    if not ok_tol:
        return "REJECTED: INPUT"
    if not isinstance(chk, StateCheck):
        return "REJECTED: INPUT/DOMAIN"          # VBA: ``Not IsArray(chk)`` (the Python mirror says DOMAIN)
    if chk.domain != 1:
        return "REJECTED: DOMAIN"
    if chk.branch != 1:
        return "REJECTED: BRANCH"
    if abs(chk.angular) > tol_a:                  # type: ignore[arg-type]  (domain == 1 guarantees a number)
        return "REJECTED: TANGENCY"
    if abs(chk.ejc) > tol_e:                      # type: ignore[arg-type]
        return "REJECTED: EJC"
    if abs(chk.closure) > tol_e:                  # type: ignore[arg-type]
        return "REJECTED: CLOSURE"
    if chk.ordering != 1:
        return "REJECTED: ORDERING"
    return "ADMISSIBLE"


def ejc_accept(sn_kpa: float, phi_deg: float, sci: float, mb: float, s: float, a: float, p: float, q: float,
               tol_a: float = 1e-10, tol_e: float = 1e-8) -> str:
    """``EJC_ACCEPT``: acceptance label of the candidate state."""
    return accept_text(ejc_check(sn_kpa, phi_deg, sci, mb, s, a, p, q, tol_a), tol_a, tol_e)


def ejc_accept_state(sn_kpa: float, c_kpa: float, phi_deg: float, sci: float, mb: float, s: float, a: float,
                     p: float, q: float, tol_a: float = 1e-10, tol_e: float = 1e-8) -> str:
    """``EJC_ACCEPT_STATE``: acceptance label of the force-used state."""
    return accept_text(ejc_check_state(sn_kpa, c_kpa, phi_deg, sci, mb, s, a, p, q, tol_a), tol_a, tol_e)
