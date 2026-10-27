# PyInstaller spec for a single-file slope3d.exe.   Build (from the repo root, Windows):
#     pip install -e ".[exe]"
#     pyinstaller packaging\slope3d.spec --noconfirm
# Result: dist\slope3d.exe  (double-click = GUI; "slope3d.exe run case.yaml" and "slope3d.exe validate" also work).
# NOTE: written on Linux and NOT test-built; if the build misses a Streamlit resource, add it to `datas` below.
from PyInstaller.utils.hooks import collect_all, collect_submodules, copy_metadata

datas, binaries, hiddenimports = [], [], []
for pkg in ("streamlit", "plotly", "matplotlib", "openpyxl", "pandas"):
    d, b, h = collect_all(pkg)
    datas += d
    binaries += b
    hiddenimports += h
for pkg in ("streamlit", "plotly", "numpy", "pandas", "pyarrow", "altair", "pyyaml", "openpyxl", "matplotlib"):
    try:
        datas += copy_metadata(pkg)
    except Exception:
        pass
hiddenimports += collect_submodules("slope3d")
# streamlit executes app.py as a script, so it must exist on disk inside the bundle
datas += [("../slope3d/gui/app.py", "slope3d/gui")]
# the case files and the reference data are optional extras for 'validate --slow'
datas += [("../cases", "cases")]

a = Analysis(["slope3d_app.py"], pathex=[".."], binaries=binaries, datas=datas, hiddenimports=hiddenimports,
             excludes=["tkinter", "pytest"], noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, a.binaries, a.datas, [], name="slope3d", debug=False, strip=False, upx=False,
          console=True, runtime_tmpdir=None)
