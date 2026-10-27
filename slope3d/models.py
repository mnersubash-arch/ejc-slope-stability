"""Dataclasses for the configuration, elements and results.

Units: lengths m; unit weight kN/m^3; stresses and pore pressure kPa inside the engine; sigma_ci MPa;
angles degrees at the interface (radians internally); forces kN (per unit width in 2D).
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field, fields
from enum import Enum
from typing import Any, Optional

import numpy as np

CLOSURE_NAMES = {
    1: "Thesis closure (sigma2' = sigma3'/nu - sigma1')",
    2: "Moscarda closure (sigma2' = nu (sigma1' + sigma3'))",
    3: "Conventional HB comparator (3D geometry, 2D strength)",
}
MESH3D_RANGE = (2, 400)
MESH2D_RANGE = (2, 20000)


DIAGNOSTIC_F_LABEL = "Diagnostic F (unchecked roots, not an accepted result)"


class Status(str, Enum):
    """Run status. A factor of safety is reported only for ``ACCEPTED``."""

    ACCEPTED = "ACCEPTED"
    NOT_ACCEPTED = "NOT ACCEPTED"
    DIAGNOSTIC_ONLY = "DIAGNOSTIC ONLY"


@dataclass
class Config:
    """All inputs of one case (the Interface sheet of the workbook). Defaults: Appendix E surface, thesis closure.

    Geometry [m]: ``xc``, ``zc`` slip-sphere centre (y_c = 0 is fixed by the symmetry model), ``radius`` R,
    ``height`` slope height H. ``slope`` = tan(face angle) (-); use :meth:`with_face_angle` to set it in degrees.
    ``gamma`` unit weight [kN/m^3].

    Water: ``water_angle`` [deg] and ``water_height`` [m] describe a planar water table (height 0 = dry);
    ``u_extra`` is an additional uniform pore pressure [kPa].

    Strength: ``strength`` 1 = Hoek-Brown / EJC eq. 4.21, 0 = Mohr-Coulomb (``c_mc`` [kPa], ``phi_mc`` [deg]).
    Rock mass: ``sci`` sigma_ci [MPa], ``mi``, ``gsi``, ``d``; ``a_override`` replaces the GSI-derived ``a`` when set.
    ``closure``: 1 thesis (default), 2 Moscarda, 3 conventional 2D HB comparator; ``nu`` Poisson's ratio (-).

    Tolerances: ``tol_a`` [rad] tangency, ``tol_e`` [MPa] EJC/closure residual, ``tol_f`` relative F;
    ``checks`` 1 = acceptance checks on (0 = result is DIAGNOSTIC ONLY and F is withheld); ``max_iter``.
    """

    # geometry and slope
    xc: float = 355.7423115464064
    zc: float = 644.4192948645371
    radius: float = 725.2669975187362
    height: float = 360.0
    slope: float = math.tan(math.radians(50.0))
    gamma: float = 27.0
    # water
    water_angle: float = 0.0
    water_height: float = 0.0
    u_extra: float = 0.0
    # strength model
    strength: int = 1
    c_mc: float = 7.2
    phi_mc: float = 30.0
    # rock mass
    sci: float = 77.7
    mi: float = 10.0
    gsi: float = 42.0
    d: float = 0.0
    a_override: Optional[float] = 0.5
    # closure
    closure: int = 1
    nu: float = 0.271
    # tolerances and run controls
    tol_a: float = 1e-10
    tol_e: float = 1e-8
    tol_f: float = 1e-8
    checks: int = 1
    max_iter: int = 2000

    # ---- derived rock-mass parameters ---------------------------------------------------------
    @property
    def mb(self) -> float:
        """mb = mi exp((GSI-100)/(28-14D))."""
        return self.mi * math.exp((self.gsi - 100) / (28 - 14 * self.d))

    @property
    def s(self) -> float:
        """s = exp((GSI-100)/(9-3D))."""
        return math.exp((self.gsi - 100) / (9 - 3 * self.d))

    @property
    def a_gsi(self) -> float:
        """a = 1/2 + (1/6)(exp(-GSI/15) - exp(-20/3))."""
        return 0.5 + (math.exp(-self.gsi / 15) - math.exp(-20.0 / 3)) / 6

    @property
    def a(self) -> float:
        """The exponent used: ``a_override`` if given, else the GSI value."""
        return self.a_gsi if self.a_override is None else float(self.a_override)

    @property
    def closure_pq(self) -> tuple[float, float]:
        """(p, q) of sigma2' = p sigma1' + q sigma3'."""
        if self.closure == 1:
            return -1.0, 1.0 / self.nu
        if self.closure == 2:
            return self.nu, self.nu
        if self.closure == 3:
            return 0.0, 1.0
        raise ValueError("closure must be 1, 2 or 3")

    @property
    def face_angle_deg(self) -> float:
        return math.degrees(math.atan(self.slope))

    @property
    def is_thesis_closure(self) -> bool:
        """A thesis-closure result needs HB on, closure 1 and checks on (and an ACCEPTED dimension)."""
        return self.strength == 1 and self.closure == 1 and self.checks == 1

    @property
    def closure_label(self) -> str:
        return "Mohr-Coulomb (verification model)" if self.strength == 0 else CLOSURE_NAMES.get(self.closure, "?")

    def with_face_angle(self, degrees: float) -> "Config":
        """Return a copy with ``slope = tan(degrees)``."""
        return self.replace(slope=math.tan(math.radians(degrees)))

    def replace(self, **kw: Any) -> "Config":
        d = asdict(self)
        unknown = set(kw) - set(d)
        if unknown:
            raise TypeError(f"unknown Config fields: {sorted(unknown)}")
        d.update(kw)
        return Config(**d)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Config":
        data = dict(data)
        face = data.pop("face_angle_deg", None)
        if face is not None and "slope" not in data:
            data["slope"] = math.tan(math.radians(float(face)))
        names = {f.name for f in fields(cls)}
        unknown = set(data) - names
        if unknown:
            raise ValueError(f"unknown input(s): {sorted(unknown)}")
        return cls(**data)

    def input_issue(self) -> str:
        """Non-geometry input checks of ``ReadConfig`` (VBA lines 136-141). Empty string = fine."""
        vals = [self.xc, self.zc, self.radius, self.height, self.slope, self.gamma, self.water_angle,
                self.water_height, self.u_extra, self.c_mc, self.phi_mc, self.sci, self.mi, self.gsi, self.d,
                self.nu, self.tol_a, self.tol_e, self.tol_f]
        if self.a_override is not None:
            vals.append(self.a_override)
        if not all(isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool) and math.isfinite(v)
                   for v in vals):
            return "every input must be a finite number"
        if self.checks not in (0, 1):
            return "acceptance switch must be 0 or 1"
        if self.strength not in (0, 1) or self.closure not in (1, 2, 3):
            return "unsupported model/closure"
        if 28 - 14 * self.d == 0 or 9 - 3 * self.d == 0:
            return "disturbance factor D out of range"
        mb, s, a = self.mb, self.s, self.a
        if self.sci <= 0 or mb <= 0 or s < 0 or a <= 0 or a > 1 or self.nu <= 0 or self.nu >= 0.5:
            return "invalid rock/Poisson input (need sigma_ci > 0, mb > 0, 0 < a <= 1, 0 < nu < 0.5)"
        if (self.c_mc < 0 or self.phi_mc < 0 or self.phi_mc >= 90 or self.u_extra < 0 or self.water_height < 0
                or self.water_angle < 0 or self.water_angle >= 90):
            return "invalid strength/water input"
        if (self.tol_a <= 0 or self.tol_a > 1e-10 or self.tol_e <= 0 or self.tol_e > 1e-8
                or self.tol_f <= 0 or self.tol_f > 1e-8):
            return "invalid acceptance tolerances (angle <= 1e-10 rad, EJC <= 1e-8 MPa, F <= 1e-8)"
        if not isinstance(self.max_iter, (int, np.integer)) or isinstance(self.max_iter, bool) or not 1 <= self.max_iter <= 50000:
            return "iteration limit must be an integer 1-50000"
        return ""


@dataclass
class Elements:
    """Struct-of-arrays element set (3D Hungr columns or 2D Bishop slices). Same fields as the VBA ``ModelElement``.

    x, y, zb, h [m]; ah, area [m^2] (plan area, base area); g, kap, rho (-, -, m); weight [kN];
    pressure u [kPa]; sn trial normal stress [kPa]; phi [deg]; cohesion [kPa]; final_n, tau [kPa];
    angular [rad]; ejc, closure_res [MPa]; equilibrium relative residual |W - ...|/W; s1, s2, s3 [MPa].
    """

    x: np.ndarray
    y: np.ndarray
    zb: np.ndarray
    h: np.ndarray
    ah: np.ndarray
    g: np.ndarray
    area: np.ndarray
    rho: np.ndarray
    kap: np.ndarray
    weight: np.ndarray
    pressure: np.ndarray
    sn: np.ndarray
    ix: np.ndarray          # x-index in the n x n (3D) or n (2D) grid
    iy: np.ndarray          # y-index (0 in 2D)
    phi: np.ndarray = field(default=None)       # type: ignore[assignment]
    cohesion: np.ndarray = field(default=None)  # type: ignore[assignment]
    final_n: np.ndarray = field(default=None)   # type: ignore[assignment]
    tau: np.ndarray = field(default=None)       # type: ignore[assignment]
    angular: np.ndarray = field(default=None)   # type: ignore[assignment]
    ejc: np.ndarray = field(default=None)       # type: ignore[assignment]
    closure_res: np.ndarray = field(default=None)  # type: ignore[assignment]
    equilibrium: np.ndarray = field(default=None)  # type: ignore[assignment]
    s1: np.ndarray = field(default=None)        # type: ignore[assignment]
    s2: np.ndarray = field(default=None)        # type: ignore[assignment]
    s3: np.ndarray = field(default=None)        # type: ignore[assignment]
    label: np.ndarray = field(default=None)     # type: ignore[assignment]

    def __post_init__(self) -> None:
        n = len(self.x)
        z = lambda v: np.full(n, v, dtype=float)   # noqa: E731
        if self.phi is None:
            self.phi = z(-1.0)                      # VBA: .phi = -1 (also the "no previous phi" seed)
        for name in ("cohesion", "final_n", "tau", "angular", "ejc", "closure_res", "equilibrium", "s1", "s2", "s3"):
            if getattr(self, name) is None:
                setattr(self, name, z(0.0))
        if self.label is None:
            self.label = np.full(n, "NOT EVALUATED", dtype=object)

    @property
    def n(self) -> int:
        return len(self.x)


@dataclass
class IterationRecord:
    """One row of the VBA ``VBA Iterations`` log."""

    iteration: int
    f: float
    df_rel: float
    max_move_kpa: float
    projections: int


@dataclass
class ModelResult:
    """Result of one dimension (3D or 2D). ``fos`` is ``None`` unless ``status`` is ACCEPTED.

    ``internal_f`` holds the last trial F *for diagnostics and tests only*; user interfaces must never display it
    for a run that is not ACCEPTED ("FoS withheld"). A DIAGNOSTIC ONLY run additionally exposes ``diagnostic_f``, which is
    shown only under its own label (never as a factor of safety).
    """

    dim: int = 3
    mesh: int = 0
    status: Status = Status.NOT_ACCEPTED
    reason: str = "NOT RUN"
    internal_f: Optional[float] = None
    moment: Optional[float] = None       # spherical moment-form diagnostic sum(tau A rho)/sum(W kappa rho)
    iterations: int = 0
    active: int = 0
    admissible: int = 0
    evaluated: bool = False
    df_rel: float = 0.0
    global_res: float = 0.0
    max_angle: float = 0.0
    max_ejc: float = 0.0
    max_eq: float = 0.0
    min_normal: float = 1e300
    tension: int = 0
    projections: int = 0
    volume: float = 0.0
    history: list[IterationRecord] = field(default_factory=list)

    @property
    def accepted(self) -> bool:
        return self.status == Status.ACCEPTED

    @property
    def fos(self) -> Optional[float]:
        """The factor of safety, or ``None`` (withheld) unless the run is ACCEPTED."""
        return self.internal_f if self.status == Status.ACCEPTED else None

    @property
    def diagnostic_f(self) -> Optional[float]:
        """The converged trial F of a DIAGNOSTIC ONLY run. Never a factor of safety; always show it with
        :data:`DIAGNOSTIC_F_LABEL` and never in the ACCEPTED style. ``None`` for every other status."""
        return self.internal_f if self.status == Status.DIAGNOSTIC_ONLY else None

    @property
    def f_display(self) -> str:
        """The one-line text every interface uses: ``F = x`` (ACCEPTED), the labelled diagnostic F, or ``FoS withheld``."""
        if self.fos is not None:
            return f"F = {self.fos:.6f}"
        if self.diagnostic_f is not None:
            return f"{DIAGNOSTIC_F_LABEL}: {self.diagnostic_f:.6f}"
        return "FoS withheld"

    @property
    def fos_text(self) -> str:
        return f"{self.fos:.6f}" if self.fos is not None else "WITHHELD"

    @property
    def moment_diagnostic(self) -> Optional[float]:
        """Moment-form ratio; shown only for an accepted run and always labelled 'diagnostic only'."""
        return self.moment if self.accepted else None


@dataclass
class ModelRun:
    """A solved dimension: its result and the final element states (``None`` if the run was refused)."""

    result: ModelResult
    elements: Optional[Elements]


@dataclass
class CaseResult:
    """A full case: 3D columns and 2D central-section slices for the same configuration."""

    config: Config
    mesh3d: int
    mesh2d: int
    run3d: ModelRun
    run2d: ModelRun

    @property
    def ratio(self) -> Optional[float]:
        """F(3D)/F(2D), only when both are ACCEPTED."""
        f3, f2 = self.run3d.result.fos, self.run2d.result.fos
        return f3 / f2 if (f3 is not None and f2 is not None and f2 != 0) else None

    @property
    def label(self) -> str:
        c = self.config
        if c.strength == 0:
            return "Mohr-Coulomb verification model"
        if c.is_thesis_closure:
            return "Thesis closure (3D EJC)"
        if c.closure == 3:
            return CLOSURE_NAMES[3]
        return CLOSURE_NAMES[c.closure] + (" - acceptance checks OFF" if c.checks == 0 else "")
