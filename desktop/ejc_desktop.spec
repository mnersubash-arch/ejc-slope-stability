# PyInstaller spec for EJC-Slope-Stability(.exe): one file that serves the web app (webapp/) on this computer and
# opens it in the browser. Build from the repository root:
#     python desktop/build_site.py
#     pyinstaller desktop/ejc_desktop.spec --noconfirm --distpath desktop/dist --workpath desktop/build/pyinstaller
# The Windows .exe is built by .github/workflows/desktop.yml on a Windows runner.
from pathlib import Path

HERE = Path(SPECPATH)
a = Analysis([str(HERE / "ejc_desktop.py")], datas=[(str(HERE / "build" / "site"), "site")], excludes=["tkinter"],
             noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, a.binaries, a.datas, [], name="EJC-Slope-Stability", console=True, debug=False,
          strip=False, upx=False, runtime_tmpdir=None)
