"""Hungr (1987) simplified 3D columns and 2D Bishop slices: the iteration, final-state validation and result objects.

Port of ``Thesis_Columns.bas`` (``SolveStrength``, ``ValidateFinal``, ``SolveModel``). The VBA is the specification.

Result rules (non-negotiable): a factor of safety exists only for an ACCEPTED dimension; there is no fallback strength,
no clipped asin and no substitute C1. Where the strength step fails, the run is NOT ACCEPTED and F is withheld.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

from . import strength as st
from ._vector import check_state_vec, fastphi_closure3_vec
from .geometry import PI_C, GeometryError, build_elements, geometry_issue, mesh_issue
from .models import (CaseResult, Config, Elements, IterationRecord, ModelResult, ModelRun, Status)

LABEL_OK = "ADMISSIBLE"


def seqsum(a: np.ndarray) -> float:
    """Sequential (left-to-right) sum, as the VBA accumulates; avoids numpy's pairwise-summation differences."""
    return float(np.cumsum(a)[-1]) if len(a) else 0.0


def solve_strength(cfg: Config, els: Elements, p: float, q: float) -> int:
    """``SolveStrength`` for every element. Updates phi, cohesion, label; returns the number of failures."""
    n = els.n
    if cfg.strength == 0:                                   # Mohr-Coulomb: fixed c and phi, label left as is
        els.phi[:] = cfg.phi_mc
        els.cohesion[:] = cfg.c_mc
        return 0
    sci, mb, s, a = cfg.sci, cfg.mb, cfg.s, cfg.a
    solver_tol = cfg.tol_a
    closure3 = (p == 0 and q == 1)
    if closure3:
        solver_tol = min(1e-12, cfg.tol_a / 10)
        ph = fastphi_closure3_vec(els.sn, sci, mb, s, a, els.phi, solver_tol)
    else:
        ph = np.full(n, np.nan)
        for i in range(n):
            r = st.fastphi(float(els.sn[i]), sci, mb, s, a, p, q, float(els.phi[i]), solver_tol)
            if not isinstance(r, st.ErrorCode):
                ph[i] = r
    conv = ~np.isnan(ph)
    els.label[~conv] = "NO CONVERGED EQ.4.21 ROOT"
    with np.errstate(all="ignore"):
        b = mb * els.sn / (1000 * sci) + s
        outside = conv & (b <= 0)
        els.label[outside] = "OUTSIDE HB DOMAIN"
        good = conv & ~outside
        phi_rad = ph * PI_C / 180
        h_t = (b / (1 + np.sin(phi_rad) / a)) ** a
        coh = 1000 * sci * h_t * np.cos(phi_rad) / 2 - els.sn * np.tan(phi_rad)
    if good.any():
        els.phi[good] = ph[good]
        els.cohesion[good] = coh[good]
        chk = check_state_vec(els.sn[good], coh[good], ph[good], sci, mb, s, a, p, q, cfg.tol_a, cfg.tol_e)
        els.label[good] = chk.label
    if cfg.strength == 1 and cfg.checks == 0:
        # DIAGNOSTIC MODE (Interface!R19, P19 = 0): the first converged eq. 4.21 root is used with NO admissibility
        # test. Only a missing root or a state outside the HB domain still counts as a failure. The acceptance labels
        # stay in the table for information.
        return int(np.count_nonzero((els.label == "NO CONVERGED EQ.4.21 ROOT") | (els.label == "OUTSIDE HB DOMAIN")))
    return int(np.count_nonzero(els.label != LABEL_OK))


def validate_final(cfg: Config, els: Elements, fos: float, result: ModelResult, rebuild: bool = True) -> bool:
    """``ValidateFinal`` on the ACTUAL final state at ``fos``.

    Per element: equilibrium |W - ((sn + u) A_h + tau A kappa / F)| / W <= 1e-9, tau >= 0, EJC acceptance (HB);
    globally |sum(tau A)/sum(W kappa) - F|/F <= tol_f. With ``rebuild=False`` the stored ``final_n`` is used as it is
    (this is how the +1 kPa perturbation test works).
    """
    p, q = cfg.closure_pq
    n = els.n
    result.admissible = 0
    result.evaluated = True
    result.max_angle = result.max_ejc = result.max_eq = 0.0
    result.tension = 0
    result.min_normal = 1e300
    with np.errstate(all="ignore"):
        t = np.tan(els.phi * PI_C / 180)
        ma = els.g + els.kap * t / fos
    bad = (ma <= 0) | (els.g <= 0) | (els.weight <= 0)
    if bad.any():
        els.label[int(np.argmax(bad))] = "INVALID EQUILIBRIUM"
        result.evaluated = False
        return False
    if rebuild:
        els.final_n = (els.weight - els.pressure * els.ah - els.cohesion * els.area * els.kap / fos) / (els.area * ma)
    normal = els.final_n
    els.tau = els.cohesion + normal * t
    els.equilibrium = np.abs(els.weight - ((normal + els.pressure) * els.ah + els.tau * els.area * els.kap / fos)) / els.weight
    result.max_eq = float(els.equilibrium.max())
    result.min_normal = float(normal.min())
    result.tension = int(np.count_nonzero(normal < 0))
    label = np.full(n, LABEL_OK, dtype=object)
    diagnostic = cfg.strength == 1 and cfg.checks == 0
    if diagnostic:
        label[els.equilibrium > 1e-9] = "REJECTED: EQUILIBRIUM/SHEAR"      # local equilibrium only
    else:
        label[(els.tau < 0) | (els.equilibrium > 1e-9)] = "REJECTED: EQUILIBRIUM/SHEAR"
    if cfg.strength == 1:
        chk = check_state_vec(normal, els.cohesion, els.phi, cfg.sci, cfg.mb, cfg.s, cfg.a, p, q, cfg.tol_a, cfg.tol_e)
        if not diagnostic:
            label = np.where(chk.label != LABEL_OK, chk.label, label)
        dom = chk.domain == 1
        els.angular = np.where(dom, chk.angular, els.angular)
        els.ejc = np.where(dom, chk.ejc, els.ejc)
        els.closure_res = np.where(dom, chk.closure, els.closure_res)
        els.s1 = np.where(dom, chk.s1, els.s1)
        els.s3 = np.where(dom, chk.s3, els.s3)
        els.s2 = np.where(dom, p * els.s1 + q * els.s3, els.s2)
        if dom.any():
            result.max_angle = float(np.abs(els.angular[dom]).max())
            result.max_ejc = float(np.abs(els.ejc[dom]).max())
        if not dom.all() and not diagnostic:
            result.evaluated = False
    els.label = label
    result.admissible = int(np.count_nonzero(label == LABEL_OK))
    resist = seqsum(els.tau * els.area)
    drive = seqsum(els.weight * els.kap)
    moment_r = seqsum(els.tau * els.area * els.rho)
    moment_d = seqsum(els.weight * els.kap * els.rho)
    if drive <= 0 or moment_d <= 0:
        return False
    result.global_res = abs(resist / drive - fos) / fos
    result.moment = moment_r / moment_d
    return bool(result.evaluated and result.admissible == n and result.global_res <= cfg.tol_f)


def _refused(result: ModelResult, reason: str) -> ModelRun:
    result.reason = reason
    return ModelRun(result, None)


def solve_model(cfg: Config, dim: int, mesh: int) -> ModelRun:
    """``SolveModel``: iterate F and the element normal stresses, then validate the final state.

    Returns the :class:`ModelResult` (status, reason, F only if ACCEPTED) and the final element states.
    """
    result = ModelResult(dim=dim, mesh=mesh)
    issue = geometry_issue(cfg)
    if issue:
        return _refused(result, issue)
    issue = cfg.input_issue()
    if issue:
        return _refused(result, "INVALID INPUT: " + issue)
    issue = mesh_issue(dim, mesh)
    if issue:
        return _refused(result, "INVALID MESH: " + issue)
    try:
        els, result.volume = build_elements(cfg, dim, mesh)
    except GeometryError as exc:
        return _refused(result, str(exc))
    n = els.n
    result.active = n
    p, q = cfg.closure_pq
    hb = cfg.strength == 1
    limit = -1000 * cfg.s * cfg.sci / cfg.mb               # tensile-limit stress [kPa]
    inside = limit + max(1e-9, 0.001 * abs(limit))
    fos = 1.0
    if hb:
        proj = els.sn <= limit
        els.sn[proj] = inside
        result.projections += int(proj.sum())
    diagnostic = hb and cfg.checks == 0
    tol_conv = min(1e-12, cfg.tol_f / 10)
    if diagnostic:
        # Unchecked roots carry the root-finder noise (tol_a, ~1e-10 rad), so the 1e-12 / 1e-9 kPa criteria are
        # unreachable. Diagnostic mode converges on F alone: |dF|/F <= tol_f/100 (1e-10 by default).
        tol_conv = cfg.tol_f / 100
    for it in range(1, cfg.max_iter + 1):
        failures = solve_strength(cfg, els, p, q)
        drive = seqsum(els.weight * els.kap)
        if failures > 0:
            result.iterations = it
            result.reason = f"{failures} of {n} elements have no admissible eq.4.21 strength; no FoS."
            return ModelRun(result, els)
        if drive <= 0:
            result.reason = "Nonpositive driving sum."
            return ModelRun(result, els)
        with np.errstate(all="ignore"):
            t = np.tan(els.phi * PI_C / 180)
            ma = els.g + els.kap * t / fos
            if (ma <= 0).any():
                result.reason = "Nonpositive m_alpha."
                return ModelRun(result, els)
            sf = (els.cohesion * els.ah + (els.weight - els.pressure * els.ah) * t) / ma
            resist = seqsum(sf)
            raw = (els.weight - els.pressure * els.ah - els.cohesion * els.area * els.kap / fos) / (els.area * ma)
        if hb:
            pr = raw <= limit
            raw = np.where(pr, inside, raw)
            result.projections += int(pr.sum())
        max_move = float(np.abs(raw - els.sn).max())
        els.sn = raw
        fn = resist / drive
        if not fn > 0:
            result.reason = "Nonpositive trial FoS."
            return ModelRun(result, els)
        result.df_rel = abs(fn - fos) / fn
        result.iterations = it
        result.internal_f = fn
        result.history.append(IterationRecord(it, fn, result.df_rel, max_move, result.projections))
        if result.df_rel <= tol_conv and (diagnostic or max_move <= 1e-9):
            if validate_final(cfg, els, fn, result):
                if cfg.checks == 1:
                    result.status = Status.ACCEPTED
                    result.reason = "Every final element and global equilibrium passed."
                else:
                    result.status = Status.DIAGNOSTIC_ONLY
                    result.reason = ("Acceptance switch is off; the diagnostic F is not an accepted result."
                                     if not hb else
                                     "Diagnostic mode (checks = 0): unchecked eq. 4.21 roots; only local equilibrium "
                                     "and convergence of the trial value were tested, so it is never an accepted result.")
                return ModelRun(result, els)
        fos = fn
    result.reason = "Iteration limit reached; no accepted FoS."
    return ModelRun(result, els)


def run_case(cfg: Config, mesh3d: int = 20, mesh2d: int = 40) -> CaseResult:
    """Run the 3D half-model (``mesh3d`` x ``mesh3d`` columns) and the 2D central section (``mesh2d`` slices)."""
    r3 = solve_model(cfg, 3, mesh3d)
    r2 = solve_model(cfg, 2, mesh2d)
    return CaseResult(cfg, mesh3d, mesh2d, r3, r2)
