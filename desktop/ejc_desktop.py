"""EJC Slope Stability on the desktop.

Serves the web app (webapp/, bundled in this program) on this computer at http://127.0.0.1:<port>/ and opens it in
the default browser. The analysis runs in the browser with the same engine as the published page; nothing is sent
anywhere (the page fetches its fonts from Google Fonts when online and works without them). Close this window to quit.
"""
from __future__ import annotations

import argparse
import functools
import http.server
import os
import sys
import threading
import urllib.request
import webbrowser

APP = "EJC Slope Stability"
PORTS = range(8611, 8631)          # a fixed port keeps the browser's remembered inputs from one run to the next
PING = "/__ejc_desktop"
MIME = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8", ".json": "application/json",
        ".css": "text/css; charset=utf-8", ".svg": "image/svg+xml", ".png": "image/png", ".ico": "image/x-icon",
        ".txt": "text/plain; charset=utf-8"}
NOTE = "Research tool for the prescribed-surface comparison in Mahat (2026) thesis; not for design."
LOCAL = urllib.request.build_opener(urllib.request.ProxyHandler({}))   # this computer: never through a proxy


def site_dir() -> str:
    """The app's files: inside this program (PyInstaller), or desktop/build/site when run from the repository."""
    here = os.path.dirname(os.path.abspath(__file__))
    for base in (getattr(sys, "_MEIPASS", None), here):
        for cand in ((os.path.join(base, "site"), os.path.join(base, "build", "site")) if base else ()):
            if os.path.isfile(os.path.join(cand, "index.html")):
                return cand
    sys.exit("The app's files are missing: run  python desktop/build_site.py  first.")


class Handler(http.server.SimpleHTTPRequestHandler):
    # Windows can map .js to text/plain in its registry, and a worker refuses a script served that way
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **MIME}

    def do_GET(self):
        if self.path == PING:                 # lets a second start find this one instead of starting another
            body = APP.encode()
            self.send_response(200)
            self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        super().end_headers()

    def log_message(self, fmt, *args):        # the window shows only what the user needs
        pass


class Server(http.server.ThreadingHTTPServer):
    # on Windows, address reuse would let two programs listen on one port; elsewhere it only skips TIME_WAIT
    allow_reuse_address = sys.platform != "win32"


def is_ours(port: int) -> bool:
    try:
        with LOCAL.open(f"http://127.0.0.1:{port}{PING}", timeout=2) as r:
            return r.read().decode() == APP
    except OSError:
        return False


def open_server(ports) -> tuple:
    """(server, port) on the first free port; (None, port) if this app already runs there."""
    handler = functools.partial(Handler, directory=site_dir())
    for port in ports:
        try:
            httpd = Server(("127.0.0.1", port), handler)
        except OSError:
            if port and is_ours(port):
                return None, port
            continue
        return httpd, httpd.server_address[1]          # the port bound (port 0 asks for any free one)
    httpd = Server(("127.0.0.1", 0), handler)
    return httpd, httpd.server_address[1]


def check() -> int:
    """Start on a free port, fetch every file the page needs, report, and stop (used by the build)."""
    httpd, port = open_server([0])
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    want = {"": ("text/html", "<title>EJC Slope Stability</title>"), "engine.js": ("text/javascript", "EJC"),
            "presets.js": ("text/javascript", "EJC_PRESETS"), "checks.js": ("text/javascript", "EJC_CHECKS"),
            "jobs.js": ("text/javascript", "EJC_JOBS"), "worker.js": ("text/javascript", "importScripts"),
            "plotly.min.js": ("text/javascript", "plotly.js v2.35.2")}
    ok = True
    for path, (ctype, text) in want.items():
        with LOCAL.open(f"http://127.0.0.1:{port}/{path}", timeout=10) as r:
            body, got = r.read().decode("utf-8", "replace"), r.headers.get("Content-Type", "")
        good = r.status == 200 and got.startswith(ctype) and text in body
        ok &= good
        print(f"{'ok  ' if good else 'FAIL'} /{path} ({got}, {len(body):,} characters)")
    with LOCAL.open(f"http://127.0.0.1:{port}/", timeout=10) as r:
        local = 'loadScript("plotly.min.js")' in r.read().decode()
    ok &= local
    print(f"{'ok  ' if local else 'FAIL'} the page loads the bundled Plotly.js first")
    httpd.shutdown()
    print("All files served." if ok else "Some files are missing or wrong.")
    return 0 if ok else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=f"{APP}: runs the app on this computer and opens it in your browser.")
    ap.add_argument("--port", type=int, help="the port to use (default: the first free one from 8611)")
    ap.add_argument("--no-browser", action="store_true", help="do not open the browser")
    ap.add_argument("--check", action="store_true", help="serve, fetch every file, report and exit (for builds)")
    a = ap.parse_args(argv)
    if a.check:
        return check()
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.kernel32.SetConsoleTitleW(APP)
        except Exception:
            pass
    httpd, port = open_server([a.port] if a.port else PORTS)
    url = f"http://127.0.0.1:{port}/"
    if httpd is None:
        print(f"{APP} is already running at {url}; opening it in your browser.")
        webbrowser.open(url)
        return 0
    print(f"{APP}\n\nRunning on this computer at {url}\n"
          "Your browser opens it now; if it does not, open that address.\n"
          "Keep this window open while you use the app. Close it to quit.\n\n" + NOTE)
    if not a.no_browser:
        threading.Timer(0.6, webbrowser.open, args=(url,)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
