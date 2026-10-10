"""slope3d - standalone 3D limit-equilibrium rock-slope stability with a 3D Hoek-Brown (EJC eq. 4.21) strength model.

A port of the verified GENG5512 thesis workbook (VBA: EJC_Eq421.bas, Thesis_Columns.bas). The science is unchanged.
"""
from .models import CaseResult, Config, Elements, ModelResult, ModelRun, Status
from .engine import run_case, solve_model, validate_final
from .strength import ErrorCode

__version__ = "0.1.0"
__all__ = ["CaseResult", "Config", "Elements", "ErrorCode", "ModelResult", "ModelRun", "Status",
           "run_case", "solve_model", "validate_final", "__version__"]
