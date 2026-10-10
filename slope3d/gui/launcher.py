"""Start the Streamlit app from the command line or from a frozen (PyInstaller) executable."""
from __future__ import annotations

import sys
from pathlib import Path

APP = Path(__file__).with_name("app.py")


def launch(port: int = 8501, headless: bool = False) -> int:
    """Run the GUI in-process (works for ``slope3d gui`` and for the packaged .exe)."""
    try:
        from streamlit.web import bootstrap
    except ImportError:
        print("The GUI needs streamlit and plotly: pip install \"slope3d[gui]\"", file=sys.stderr)
        return 2
    flags = {"server.port": port, "server.headless": headless, "browser.gatherUsageStats": False,
             "global.developmentMode": False}
    bootstrap.load_config_options(flag_options=flags)
    bootstrap.run(str(APP), False, [], flags)
    return 0


if __name__ == "__main__":      # pragma: no cover
    sys.exit(launch())
