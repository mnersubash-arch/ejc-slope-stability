"""Streamlit Community Cloud entry point: ``streamlit run streamlit_app.py``.

Runs the slope3d GUI (``slope3d/gui/app.py``). The file is executed afresh on every rerun (``runpy``), so the app
behaves exactly as under ``slope3d gui``.
"""
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:          # works even if the package was not pip-installed into the environment
    sys.path.insert(0, str(ROOT))

runpy.run_path(str(ROOT / "slope3d" / "gui" / "app.py"), run_name="__main__")
