"""Command-line interface: ``slope3d run | validate | mesh-study | presets | gui``.

Exit status of ``run``: 0 both dimensions ACCEPTED; 1 at least one NOT ACCEPTED / DIAGNOSTIC ONLY (a valid outcome:
the factor of safety is withheld); 2 usage or input error.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .engine import run_case
from .io import (export_csv, export_xlsx, load_case, result_to_dict)
from .geometry import slip_surface_ends
from .models import DIAGNOSTIC_F_LABEL, CaseResult, Config, Status
from .presets import PRESETS


def format_result(case: CaseResult) -> str:
    """Plain-text result block. A withheld FoS is printed as 'FoS withheld'; a DIAGNOSTIC ONLY F only under its own label."""
    lines = [f"Result identity : {case.label}"]
    for tag, run, mesh, unit in (("3D", case.run3d, case.mesh3d, "columns"), ("2D", case.run2d, case.mesh2d, "slices")):
        r = run.result
        lines.append(f"{tag} ({mesh}{'x' + str(mesh) if tag == '3D' else ''} mesh, {r.active} active {unit})")
        lines.append(f"  status  : {r.status.value}")
        if r.fos is not None:
            lines.append(f"  F       : {r.fos:.12f}")
        elif r.diagnostic_f is not None:
            lines.append(f"  {DIAGNOSTIC_F_LABEL}: {r.diagnostic_f:.12f}")
        else:
            lines.append("  F       : FoS withheld")
        lines.append(f"  reason  : {r.reason}")
        if r.iterations:
            lines.append(f"  iterations: {r.iterations}")
        if r.evaluated:
            lines.append(f"  elements with sigma_n' < 0: {r.tension} (minimum sigma_n' = {r.min_normal:.2f} kPa)")
        if r.accepted and r.moment_diagnostic is not None and tag == "3D":
            lines.append(f"  moment-form ratio (diagnostic only): {r.moment_diagnostic:.6f}")
    ends = slip_surface_ends(case.config)
    if ends.get("supported"):
        lines.append(f"Slip surface    : enters at x = {ends['entry_x']:.2f} m (z = {ends['entry_z']:.2f} m); exits the face at "
                     f"x = {ends['exit_x']:.2f} m, {ends['exit_z']:.2f} m above the toe ({ends['exit_behind_toe']:.2f} m behind the toe)"
                     if ends["exit_on_face"] else
                     f"Slip surface    : enters at x = {ends['entry_x']:.2f} m; reaches the toe plane (x = 0) at z = {ends['exit_z']:.3f} m")
    lines.append(f"3D/2D ratio     : {case.ratio:.6f}" if case.ratio is not None else "3D/2D ratio     : n.a. (needs both ACCEPTED)")
    return "\n".join(lines)


def _cmd_run(args: argparse.Namespace) -> int:
    try:
        if args.preset:
            if args.preset not in PRESETS:
                print(f"unknown preset '{args.preset}'. Available: {', '.join(PRESETS)}", file=sys.stderr)
                return 2
            p = PRESETS[args.preset]
            cfg, m3, m2 = p.config, p.mesh3d, p.mesh2d
        elif args.case:
            cfg, m3, m2 = load_case(args.case)
        else:
            print("give a case file or --preset NAME", file=sys.stderr)
            return 2
    except (OSError, ValueError, TypeError) as exc:
        print(f"cannot read case: {exc}", file=sys.stderr)
        return 2
    if args.mesh3d is not None:
        m3 = args.mesh3d
    if args.mesh2d is not None:
        m2 = args.mesh2d
    t0 = time.time()
    case = run_case(cfg, m3, m2)
    print(format_result(case))
    print(f"(computed in {time.time() - t0:.2f} s)")
    if args.out:
        out = Path(args.out)
        stem = Path(args.case).stem if args.case else "slope3d"
        written = []
        fmts = {f.strip().lower() for f in args.format.split(",")}
        if "csv" in fmts:
            written += export_csv(case, out, stem)
        if "xlsx" in fmts:
            written.append(export_xlsx(case, out / f"{stem}.xlsx"))
        if "json" in fmts:
            import json
            out.mkdir(parents=True, exist_ok=True)
            (out / f"{stem}.json").write_text(json.dumps(result_to_dict(case), indent=2))
            written.append(out / f"{stem}.json")
        if "report" in fmts or "pdf" in fmts or "html" in fmts:
            from .report import write_report   # imported lazily: needs matplotlib
            written += write_report(case, out, stem, pdf=("report" in fmts or "pdf" in fmts),
                                    html=("report" in fmts or "html" in fmts))
        for w in written:
            print(f"wrote {w}")
    ok = case.run3d.result.status is Status.ACCEPTED and case.run2d.result.status is Status.ACCEPTED
    return 0 if ok else 1


def _cmd_validate(args: argparse.Namespace) -> int:
    from .selftest import run_all
    ref = args.reference
    if args.slow and ref is None:
        for cand in (Path("reference/mesh_study"), Path(__file__).resolve().parents[1] / "reference" / "mesh_study"):
            if (cand / "mesh_study_records.csv").exists():
                ref = cand
                break
    results = run_all(slow=args.slow, reference_dir=ref, echo=print)
    n_ok = sum(r.ok for r in results)
    print(f"\n{n_ok} of {len(results)} checks passed" + ("" if args.slow else "  (add --slow for the mesh study and nu thresholds)"))
    print("ALL PASS" if n_ok == len(results) else "FAILURES - do not use these numbers")
    return 0 if n_ok == len(results) else 1


def _cmd_mesh_study(args: argparse.Namespace) -> int:
    from .mesh_study import MESH, compare_with_records, run_mesh_study, write_rows
    meshes3, meshes2 = (MESH[3][:6], MESH[2][:5]) if args.quick else (MESH[3], MESH[2])
    t0 = time.time()
    rows = run_mesh_study(meshes3, meshes2, progress=lambda r: print(
        f"{r['case']:38s} {r['dim']} {r['mesh']:>12s}  F = {r['F'] if r['F'] is None else format(r['F'], '.9f')}  {r['status']}", flush=True))
    out = Path(args.out)
    write_rows(rows, out)
    print(f"wrote {out} ({len(rows)} records, {time.time() - t0:.1f} s)")
    if args.compare:
        problems = compare_with_records(rows, args.compare)
        print("matches the reference records to 1e-9" if not problems else "\n".join(problems))
        return 0 if not problems else 1
    return 0


def _cmd_presets(_: argparse.Namespace) -> int:
    for name, p in PRESETS.items():
        exp = "F withheld" if p.expected_3d is None else f"F3D {p.expected_3d}, F2D {p.expected_2d}"
        print(f"{name:40s} {exp}")
    return 0


def _cmd_gui(args: argparse.Namespace) -> int:
    from .gui.launcher import launch
    return launch(port=args.port, headless=args.headless)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="slope3d", description="3D limit-equilibrium rock-slope stability with a 3D "
                                 "Hoek-Brown (EJC eq. 4.21) strength model.")
    ap.add_argument("--version", action="version", version=f"slope3d {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run one case (3D columns and 2D central section)")
    r.add_argument("case", nargs="?", help="case file (.yaml/.yml/.json)")
    r.add_argument("--preset", help="use a built-in benchmark instead of a case file (see 'slope3d presets')")
    r.add_argument("--mesh3d", type=int, help="3D n x n half-model, 2-400 (default 20 or the case file's value)")
    r.add_argument("--mesh2d", type=int, help="2D slices, 2-20000 (default 40 or the case file's value)")
    r.add_argument("--out", help="directory for exports")
    r.add_argument("--format", default="csv,xlsx,json", help="comma list of csv, xlsx, json, report (pdf+html), pdf, html")
    r.set_defaults(func=_cmd_run)

    v = sub.add_parser("validate", help="run every regression test and print PASS/FAIL")
    v.add_argument("--slow", action="store_true", help="include the 95-record mesh study and the nu thresholds (about 30 s)")
    v.add_argument("--reference", help="directory with mesh_study_records.csv (default: reference/mesh_study)")
    v.set_defaults(func=_cmd_validate)

    m = sub.add_parser("mesh-study", help="mesh-convergence study of the five conventional-HB comparator cases")
    m.add_argument("--out", default="mesh_study_slope3d.csv")
    m.add_argument("--quick", action="store_true", help="coarser meshes only")
    m.add_argument("--compare", help="reference mesh_study_records.csv to compare against")
    m.set_defaults(func=_cmd_mesh_study)

    pr = sub.add_parser("presets", help="list the built-in benchmark cases")
    pr.set_defaults(func=_cmd_presets)

    g = sub.add_parser("gui", help="launch the Streamlit interface")
    g.add_argument("--port", type=int, default=8501)
    g.add_argument("--headless", action="store_true", help="do not open a browser window")
    g.set_defaults(func=_cmd_gui)
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":      # pragma: no cover
    sys.exit(main())
