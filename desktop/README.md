# desktop: EJC Slope Stability as a desktop program

The web app (`webapp/`) as one program, `EJC-Slope-Stability.exe` on Windows. It serves the app's files on this
computer at `http://127.0.0.1:8611/` (the next free port if that one is busy) and opens the default browser there.
The analysis runs in the browser with the same `engine.js` as the published page, checked against
`ejc_app/engine.py` (`webapp/README.md`). Nothing needs installing, and it works offline: Plotly.js 2.35.2 is bundled.

Research tool for the prescribed-surface comparison in Mahat (2026) thesis; not for design.

## Get the .exe

GitHub builds it on a Windows machine whenever `desktop/`, `webapp/` or the workflow changes on this branch
(`.github/workflows/desktop.yml`):

1. On GitHub, open the repository's **Actions** tab and the workflow **Desktop app for Windows**.
2. Open the latest green run and download **EJC-Slope-Stability-Windows** under **Artifacts** (sign in to GitHub).
3. Unzip it and double-click `EJC-Slope-Stability.exe`; `HOW-TO-RUN.txt` in the zip has the details.

Downloads of a run are kept for 90 days. To build again without a change, use **Run workflow** on that page; it
appears once the workflow is on the repository's default branch.

The program is not signed, so the first start shows Windows SmartScreen ("Windows protected your PC": More info, Run
anyway). Some antivirus programs flag PyInstaller programs they have not seen before; the program only serves files
to `127.0.0.1`.

## What differs from the claude.ai page

- The shared case list needs claude.ai, so it is hidden. Save case, Open case file, case codes and the report work:
  the browser saves the files as downloads.
- Inputs are remembered by the browser for `http://127.0.0.1:8611`. Another port (if 8611 is busy) starts fresh.
- Fonts come from Google Fonts when online; offline the browser's own fonts are used.

## Build it yourself

From the repository root, with Python 3.9 or newer:

```bash
python desktop/build_site.py                 # desktop/build/site: the page, its scripts, Plotly.js (SHA-256 checked)
python desktop/ejc_desktop.py                # run it without building a program (any OS)
pip install pyinstaller==6.11.1
pyinstaller desktop/ejc_desktop.spec --noconfirm --distpath desktop/dist --workpath desktop/build/pyinstaller
desktop/dist/EJC-Slope-Stability --check     # serves every file once and reports
```

PyInstaller builds for the system it runs on: Windows gives the `.exe`, macOS and Linux give a program for
themselves. `build_site.py` downloads Plotly.js from jsDelivr (or the npm registry) unless you pass
`--plotly path/to/plotly.min.js`; a file whose SHA-256 is not the pinned one is refused.

## Files

| File | What it is |
|---|---|
| `build_site.py` | Assembles `build/site`: `webapp/index.html` in the publish skeleton with Plotly loaded from the local file first, the five scripts, and Plotly.js |
| `ejc_desktop.py` | The program: serves `site` on `127.0.0.1` with the right file types, opens the browser, finds a copy already running; `--check` for builds |
| `ejc_desktop.spec` | PyInstaller recipe for the single-file program |
| `HOW-TO-RUN.txt` | The instructions that go with the `.exe` |

Checked on Linux: `--check` passes for the script and for the PyInstaller program, and `webapp/tools/ui_test.js`
passes all 46 checks against the page as this program serves it (with the Web Worker running). The Windows build runs
`--check` on the `.exe` before keeping it.
