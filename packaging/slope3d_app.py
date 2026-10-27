"""Entry point of the packaged executable: double-click = GUI; arguments pass through to the CLI
(``slope3d.exe run case.yaml``, ``slope3d.exe validate`` ...)."""
import multiprocessing
import sys

from slope3d.cli import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main(sys.argv[1:] or ["gui"]))
