"""Assemble the desktop app's files in desktop/build/site from webapp/: the page in the publish skeleton (as
webapp/tools/make_preview.py does), its five scripts, and Plotly.js 2.35.2 as a local file, so the charts work without
an internet connection. The Plotly file is checked against its SHA-256 before it is used.

Usage, from the repository root:
    python desktop/build_site.py                          # downloads Plotly.js 2.35.2 (jsDelivr, else npm)
    python desktop/build_site.py --plotly plotly.min.js   # or uses a copy you already have
"""
from __future__ import annotations

import argparse
import hashlib
import io
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / "webapp"
SITE = Path(__file__).resolve().parent / "build" / "site"
SCRIPTS = ("engine.js", "presets.js", "checks.js", "jobs.js", "worker.js")
PLOTLY_SHA256 = "6d21266ce1bd7d9e5ab4e115989c70c20de0382fd973a8f26ab58619eba4d603"   # plotly.js-dist-min 2.35.2
PLOTLY_FILE = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"
PLOTLY_TGZ = "https://registry.npmjs.org/plotly.js-dist-min/-/plotly.js-dist-min-2.35.2.tgz"
CDN = "https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.35.2/plotly.min.js"   # what the page loads first online
SKELETON = """<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<style>:root{color-scheme:light;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}
body{margin:0;font:14px system-ui,-apple-system,"Segoe UI",sans-serif}img{max-width:100%}[hidden]{display:none!important}</style></head><body>"""


def _get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def plotly_bytes(path: str | None) -> bytes:
    """Plotly.js 2.35.2 from a local file or the network, only if its SHA-256 is the pinned one."""
    tries = []
    if path:
        tries.append(("the file " + path, lambda: Path(path).read_bytes()))
    else:
        tries.append((PLOTLY_FILE, lambda: _get(PLOTLY_FILE)))
        tries.append((PLOTLY_TGZ, lambda: tarfile.open(fileobj=io.BytesIO(_get(PLOTLY_TGZ)), mode="r:gz")
                      .extractfile("package/plotly.min.js").read()))
    for where, fetch in tries:
        try:
            data = fetch()
        except Exception as exc:              # network or file errors: try the next source
            print(f"  could not read Plotly.js from {where}: {exc}", file=sys.stderr)
            continue
        digest = hashlib.sha256(data).hexdigest()
        if digest == PLOTLY_SHA256:
            print(f"  Plotly.js 2.35.2 from {where} ({len(data):,} bytes, SHA-256 checked)")
            return data
        print(f"  Plotly.js from {where} has SHA-256 {digest}, not the pinned {PLOTLY_SHA256}: not used",
              file=sys.stderr)
    sys.exit("No verified copy of Plotly.js 2.35.2: pass --plotly with the file from plotly.js-dist-min 2.35.2.")


def build(plotly: str | None) -> Path:
    if SITE.exists():
        shutil.rmtree(SITE)
    SITE.mkdir(parents=True)
    page = (WEB / "index.html").read_text(encoding="utf-8")
    if page.count(CDN) != 1:
        sys.exit(f"webapp/index.html no longer loads Plotly from {CDN}: update desktop/build_site.py")
    page = page.replace(CDN, "plotly.min.js")                  # the local copy first; jsDelivr stays the fallback
    (SITE / "index.html").write_text(SKELETON + page + "</body></html>", encoding="utf-8")
    for name in SCRIPTS:
        shutil.copyfile(WEB / name, SITE / name)
    (SITE / "plotly.min.js").write_bytes(plotly_bytes(plotly))
    return SITE


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--plotly", help="a local plotly.min.js from plotly.js-dist-min 2.35.2 (else it is downloaded)")
    out = build(ap.parse_args().plotly)
    print(f"Wrote {out.relative_to(ROOT)}: " + ", ".join(sorted(p.name for p in out.iterdir())))
