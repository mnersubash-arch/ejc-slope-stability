"""Vectorised accelerations of ``strength.py`` (numpy), used by the engine for large meshes.

Nothing here changes the science: each function repeats the scalar algorithm of ``strength.py`` elementwise,
with the same branch structure, constants and tolerances. Any element the vector path cannot finish
(undefined g, no bracket, iteration limit) is handed to the scalar routine, which then decides. The tests compare
the two paths on random and engine states.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import strength as st
from .strength import HI_DEG, PI_, ErrorCode, closure_constants, valid_inputs


def gstate_vec(phi, sn, sci, mb, s, a, p, q):
    """Elementwise ``GState``. Returns (g, dg, M, ok); entries with ok == False are undefined (NaN)."""
    cP, cR, cS = closure_constants(p, q)
    with np.errstate(all="ignore"):
        sin_p, cos_p = np.sin(phi), np.cos(phi)
        b_term = mb * sn / sci + s
        ok = b_term > 0
        h = (b_term / (1 + sin_p / a)) ** a
        dh = -h * cos_p / (1 + sin_p / a)
        sig3 = sn - sci / 2 * h * (1 - sin_p)
        dsig3 = -sci / 2 * (dh * (1 - sin_p) - h * cos_p)
        sig1 = sn + sci / 2 * h * (1 + sin_p)
        dsig1 = sci / 2 * (dh * (1 + sin_p) + h * cos_p)
        hb_q = mb * sig3 / sci + s
        dhb_q = mb * dsig3 / sci
        ok &= hb_q > 0
        n_t = 2 * a * mb * sci * hb_q ** (2 * a - 1) - cR * sig1 - 2 * cS * sig3
        dn_t = 2 * a * mb * sci * (2 * a - 1) * hb_q ** (2 * a - 2) * dhb_q - cR * dsig1 - 2 * cS * dsig3
        d_t = 2 * cP * sig1 + cR * sig3
        dd_t = 2 * cP * dsig1 + cR * dsig3
        ok &= (n_t + d_t != 0) & (d_t != 0) & np.isfinite(n_t) & np.isfinite(d_t)
        g = sin_p - 1 + 2 * d_t / (n_t + d_t)
        dg = cos_p - 2 * (dn_t * d_t - n_t * dd_t) / (n_t + d_t) ** 2
        m = n_t / d_t
    ok &= np.isfinite(g) & np.isfinite(dg)
    return g, dg, m, ok


def angres_from_m(phi, m, ok):
    """Elementwise ``AngRes`` given M: (residual, ok2) with the same rejection rules (M+1 = 0, |z| > 1)."""
    with np.errstate(all="ignore"):
        ok2 = ok & (m + 1 != 0)
        z = (m - 1) / (m + 1)
        ok2 &= np.abs(z) <= 1
        res = phi - np.arcsin(z)
    return res, ok2


def phi_b_vec(sn, sci, mb, s, a):
    """Elementwise ``PhiB``. Returns (phi_b, ok); ok == False where the scalar returns ``None``."""
    with np.errstate(all="ignore"):
        k = sn + s * sci / mb
        ok = k > 0
        z = 2 * k / (sci * (mb * k / sci) ** a)
        zero = ok & (z >= 1)
        pv = (1 + z * z) - np.sqrt((1 + z * z) ** 2 - (1 - z * z))
        alive = ok & ~zero
        for _ in range(8):
            f = (1 - pv) * (1 + pv / a) ** (-a) - z
            fp = -(1 + pv / a) ** (-a) * (1 + (1 - pv) / (1 + pv / a))
            alive &= fp != 0
            new = pv - f / fp
            new = np.where(new < 0, 0.0, new)
            new = np.where(new > 1 - 1e-15, 1 - 1e-15, new)
            pv = np.where(alive, new, pv)
        out = np.arcsin(pv)
    out = np.where(zero, 0.0, out)
    return out, ok


def fastphi_closure3_vec(sn_kpa, sci, mb, s, a, guess_deg, tol_a):
    """Vector ``EJC_FASTPHI`` for p = 0, q = 1. Returns phi [deg] with NaN where the scalar returns an error value.

    Elements that leave the warm-start path go to the scalar ``fastphi`` (which falls back to the scanned
    ``ejc_phi``), so the semantics are those of the VBA by construction.
    """
    p, q = 0.0, 1.0
    n = len(sn_kpa)
    out = np.full(n, np.nan)
    if not valid_inputs(sci, mb, s, a, tol_a):
        return out
    sn = sn_kpa / 1000.0
    pb, okb = phi_b_vec(sn, sci, mb, s, a)
    hi_c = HI_DEG * PI_ / 180
    lo = pb + 1e-9
    good = okb & (lo < hi_c)                       # else: #NUM! error value (NaN) as in the scalar
    fallback = np.zeros(n, dtype=bool)
    idx = np.flatnonzero(good)
    if idx.size:
        gl, _, _, okl = gstate_vec(lo[idx], sn[idx], sci, mb, s, a, p, q)
        gh, _, _, okh = gstate_vec(np.full(idx.size, hi_c), sn[idx], sci, mb, s, a, p, q)
        start_ok = okl & okh & ~(gl * gh > 0)
        fallback[idx[~start_ok]] = True
        idx = idx[start_ok]
        lo_a = lo[idx].copy()
        hi_a = np.full(idx.size, hi_c)
        glo_a = gl[start_ok].copy()
        x_a = (lo_a + hi_a) / 2
        g_rad = guess_deg[idx] * PI_ / 180
        use_guess = (g_rad > lo_a) & (g_rad < hi_a)
        x_a = np.where(use_guess, g_rad, x_a)
        sn_a = sn[idx]
        for _ in range(200):
            if idx.size == 0:
                break
            g, dg, m, ok = gstate_vec(x_a, sn_a, sci, mb, s, a, p, q)
            fallback[idx[~ok]] = True                       # GoTo FullScan
            res, okA = angres_from_m(x_a, m, ok)
            with np.errstate(all="ignore"):
                conv1 = ok & okA & (dg != 0) & (np.abs(res) <= tol_a) & (np.abs(g / dg) <= tol_a)
            out[idx[conv1]] = x_a[conv1] * 180 / PI_
            cont = ok & ~conv1
            upd_hi = cont & (glo_a * g <= 0)
            upd_lo = cont & ~upd_hi
            hi_a = np.where(upd_hi, x_a, hi_a)
            lo_a = np.where(upd_lo, x_a, lo_a)
            glo_a = np.where(upd_lo, g, glo_a)
            with np.errstate(all="ignore"):
                xn = np.where(dg != 0, x_a - g / dg, (lo_a + hi_a) / 2)
            xn = np.where((xn <= lo_a) | (xn >= hi_a), (lo_a + hi_a) / 2, xn)
            step = np.abs(xn - x_a)
            g2, dg2, m2, ok2 = gstate_vec(xn, sn_a, sci, mb, s, a, p, q)
            res2, okA2 = angres_from_m(xn, m2, ok2)
            conv2 = cont & okA2 & (np.abs(res2) <= tol_a) & (step <= tol_a)
            out[idx[conv2]] = xn[conv2] * 180 / PI_
            keep = cont & ~conv2
            idx, lo_a, hi_a, glo_a, sn_a = idx[keep], lo_a[keep], hi_a[keep], glo_a[keep], sn_a[keep]
            x_a = xn[keep]
        fallback[idx] = True                                 # iterations exhausted -> FullScan
    for i in np.flatnonzero(fallback):
        r = st.fastphi(float(sn_kpa[i]), sci, mb, s, a, p, q, float(guess_deg[i]), tol_a)
        out[i] = np.nan if isinstance(r, ErrorCode) else r
    return out


# --------------------------------------------------------------------------------------------------
# state checks
# --------------------------------------------------------------------------------------------------
@dataclass
class StateChecksVec:
    angular: np.ndarray
    ejc: np.ndarray
    closure: np.ndarray
    ordering: np.ndarray
    s1: np.ndarray
    s3: np.ndarray
    domain: np.ndarray
    branch: np.ndarray
    label: np.ndarray


def check_state_vec(sn_kpa, c_kpa, phi_deg, sci, mb, s, a, p, q, tol_a, tol_e) -> StateChecksVec:
    """Elementwise ``EJC_CHECK_STATE`` followed by ``AcceptText`` (identical order of checks)."""
    n = len(sn_kpa)
    cP, cR, cS = closure_constants(p, q)
    nan = np.full(n, np.nan)
    inputs_ok = valid_inputs(sci, mb, s, a, tol_a) and math.isfinite(p) and math.isfinite(q)
    in_range = np.isfinite(sn_kpa) & np.isfinite(c_kpa) & np.isfinite(phi_deg) & (phi_deg >= 0) & (phi_deg < 90)
    with np.errstate(all="ignore"):
        phi = phi_deg * PI_ / 180
        tau = c_kpa + sn_kpa * np.tan(phi)
        cen = sn_kpa + tau * np.tan(phi)
        radius = tau / np.cos(phi)
        sig1 = (cen + radius) / 1000.0
        sig3 = (cen - radius) / 1000.0
        sig2 = p * sig1 + q * sig3
        closure = sig2 - (p * sig1 + q * sig3)
        ordering = (sig1 >= sig2 - 1e-9) & (sig2 >= sig3 - 1e-9)
        hb_q = mb * sig3 / sci + s
        dom = hb_q > 0
        den_m = 2 * cP * sig1 + cR * sig3
        dom &= den_m != 0
        m = (2 * a * mb * sci * hb_q ** (2 * a - 1) - cR * sig1 - 2 * cS * sig3) / den_m
        dom &= m + 1 != 0
        z = (m - 1) / (m + 1)
        rad = 4 * sci ** 2 * hb_q ** (2 * a) - 3 * (sig2 - sig3) ** 2
        dom &= (np.abs(z) <= 1) & (rad >= 0)
        t_ang = np.arcsin(z)
        angular = phi - t_ang
        ejc = sig1 - (sig2 + sig3 + np.sqrt(rad)) / 2
        branch = dom & (phi >= 0) & (phi < PI_ / 2) & (t_ang >= -tol_a)
    domain = dom & in_range & inputs_ok
    branch = branch & in_range & inputs_ok
    ordering = ordering & in_range & inputs_ok
    ang_out = np.where(domain, angular, nan)
    ejc_out = np.where(domain, ejc, nan)
    clo_out = np.where(in_range & inputs_ok, closure, nan)
    s1_out = np.where(in_range & inputs_ok, sig1, nan)
    s3_out = np.where(in_range & inputs_ok, sig3, nan)
    tol_ok = (math.isfinite(tol_a) and math.isfinite(tol_e) and 0 < tol_a <= 1e-6 and 0 < tol_e <= 1e-4)
    label = np.full(n, "ADMISSIBLE", dtype=object)
    with np.errstate(invalid="ignore"):
        label[~ordering] = "REJECTED: ORDERING"
        label[domain & ~(np.abs(clo_out) <= tol_e)] = "REJECTED: CLOSURE"
        label[domain & ~(np.abs(ejc_out) <= tol_e)] = "REJECTED: EJC"
        label[domain & ~(np.abs(ang_out) <= tol_a)] = "REJECTED: TANGENCY"
    label[domain & ~branch] = "REJECTED: BRANCH"
    label[~domain] = "REJECTED: DOMAIN"
    if not tol_ok:
        label[:] = "REJECTED: INPUT"
    return StateChecksVec(ang_out, ejc_out, clo_out, ordering.astype(int), s1_out, s3_out,
                          domain.astype(int), branch.astype(int), label)
