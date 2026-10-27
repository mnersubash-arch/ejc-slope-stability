import ast
from pathlib import Path


def test_engine_has_no_streamlit_import():
    for name in ("engine.py", "verification.py"):
        tree = ast.parse((Path(__file__).resolve().parents[1] / name).read_text())
        mods = [a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names]
        mods += [n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
        assert not any(m and m.split(".")[0] == "streamlit" for m in mods)
