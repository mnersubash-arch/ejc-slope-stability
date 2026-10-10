"""Benchmark presets (the six regression cases of the brief, plus the thesis-closure default)."""
from __future__ import annotations

from dataclasses import dataclass

from .models import Config

PUBLISHED_CIRCLE = dict(xc=207.28, zc=586.53, radius=622.08)
WET = dict(water_angle=50.0, water_height=360.0)


@dataclass(frozen=True)
class Preset:
    name: str
    config: Config
    mesh3d: int = 20
    mesh2d: int = 40
    expected_3d: float | None = None     # recorded F at 3D 20x20 (None = F must be withheld)
    expected_2d: float | None = None     # recorded F at 2D 40 slices
    note: str = ""


def _hb(**kw) -> Config:
    return Config(closure=3, **kw)


PRESETS: dict[str, Preset] = {p.name: p for p in (
    Preset("Dry, a = 0.5", _hb(), expected_3d=2.5536684200376, expected_2d=2.31025502757599,
           note="Conventional HB comparator (3D geometry, 2D strength)."),
    Preset("Wet, a = 0.5", _hb(**WET), expected_3d=1.20012184467095, expected_2d=1.03587083245705,
           note="Conventional HB comparator; water table at the slope surface."),
    Preset("Dry, a = 0.509923", _hb(a_override=0.509923), expected_3d=2.47603808311863, expected_2d=2.24206326532995,
           note="Conventional HB comparator; a from GSI."),
    Preset("Wet, a = 0.509923", _hb(a_override=0.509923, **WET), expected_3d=1.15842745018062,
           expected_2d=1.00046615780691, note="Conventional HB comparator."),
    Preset("Published circle, dry, a = 0.5", _hb(**PUBLISHED_CIRCLE), expected_3d=2.14193739499583,
           expected_2d=2.00387341089056, note="Carranza-Torres & Hormazabal (2020) circle."),
    Preset("Mohr-Coulomb c = 100 kPa, phi = 0", Config(strength=0, closure=3, c_mc=100.0, phi_mc=0.0),
           expected_3d=0.255354533635107, expected_2d=0.190263799071471,
           note="Numerical verification only; not a calibrated material."),
    Preset("Thesis closure (Chile, default)", Config(),
           note="Closure 1, nu = 0.271: no admissible eq. 4.21 strength, so F is withheld (the thesis finding)."),
)}
BENCHMARK_NAMES = [n for n, p in PRESETS.items() if p.expected_3d is not None]
