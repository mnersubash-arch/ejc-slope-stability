#!/usr/bin/env python3
"""Cross-check slope3d against the verified workbook's recorded results (suite D).

Reads the QA sheets of ``GENG5512_Thesis_v14_Final.xlsm`` with openpyxl (CACHED values, no Excel needed) and confirms that
slope3d reproduces every recorded value:

* ``VBA Validation``      the 22 strength fixtures (returned angle, acceptance label, residuals)
* ``Native Workbook QA``  recorded F (expected and the worksheet engine's actual) for each of the 9 cases
* ``Native Engine QA``    recorded F (expected and the VBA engine's actual) for each of the 11 cases
* ``Refined Cases``       finest-mesh F of the 95-record mesh study (3D 320x320, 2D 2,560 slices)
* ``VBA Columns/Slices``  the element tables of the last primary VBA run (Chile, thesis closure, NOT ACCEPTED)

Usage:  python tools/compare_with_workbook.py [WORKBOOK.xlsm] [--tol 1e-9] [--skip-mesh] [--strict]

Exit status 0 only if every comparison passes. Items the workbook records but the VBA engine cannot produce are listed as
OPEN (they never count as a pass); ``--strict`` makes them fail the run.
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from slope3d import Config, Status, solve_model, strength as st, validate_final   # noqa: E402
from slope3d.io import ELEMENT_COLUMNS, elements_table                             # noqa: E402
from slope3d.presets import PRESETS                                                 # noqa: E402

DEFAULT_WB = ROOT / "reference" / "GENG5512_Thesis_v14_Final.xlsm"
MB = 10 * math.exp((42 - 100) / 28)
S = math.exp((42 - 100) / 9)
V = (120, 3.5, 0.002, 0.505)
T7SN = 13265.01174607056


@dataclass
class Report:
    tol: float
    rows: list = field(default_factory=list)
    rels: list = field(default_factory=list)     # (source, what, relative difference) of every numeric comparison

    def add(self, source: str, what: str, recorded, ours, ok: bool, kind: str = "") -> None:
        self.rows.append((source, what, recorded, ours, "PASS" if ok else ("OPEN" if kind == "open" else "FAIL"), kind))

    def num(self, source: str, what: str, recorded, ours, tol=None) -> None:
        tol = self.tol if tol is None else tol
        ok = (isinstance(recorded, (int, float)) and isinstance(ours, (int, float)) and
              abs(recorded - ours) <= tol * max(abs(recorded), 1e-300))
        if isinstance(recorded, (int, float)) and isinstance(ours, (int, float)):
            self.rels.append((source, what, abs(recorded - ours) / max(abs(recorded), 1e-300)))
        self.add(source, what, recorded, ours, ok)

    def text(self, source: str, what: str, recorded, ours) -> None:
        self.add(source, what, recorded, ours, str(recorded) == str(ours))


def fnum(text) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------------------------------------
def check_strength_fixtures(ws, rep: Report) -> None:
    """VBA Validation rows 4-24 and 28 (T1-T22)."""
    def lab(sn, ph, *a):
        return st.accept_text(st.ejc_check(sn, ph, *a), 1e-10, 1e-8)

    ph1 = st.ejc_phi(400, 77.7, MB, S, 0.5, 0, 1)
    chk1 = st.ejc_check(400, ph1, 77.7, MB, S, 0.5, 0, 1)
    c1 = 1000 * ((chk1.s1 - chk1.s3) / 2 * math.cos(ph1 * math.pi / 180) - 0.4 * math.tan(ph1 * math.pi / 180))

    def residuals(chk):
        return (chk.angular, chk.ejc) if chk and chk.domain == 1 else (None, None)

    ph = {}
    cases = {}   # row -> (returned value, label, angular, ejc)
    def put(row, value, chk, label):
        a, e = residuals(chk)
        cases[row] = (value, label, a, e)

    put(4, ph1, chk1, lab(400, ph1, 77.7, MB, S, 0.5, 0, 1))
    ph = st.ejc_phi(5740.3, *V, -1, 20); chk = st.ejc_check(5740.3, ph, *V, -1, 20); put(5, ph, chk, st.accept_text(chk, 1e-10, 1e-8))
    sn = -98.00411947453206; ph = st.ejc_phi(sn, 77.7, MB, S, 0.5, 0, 1); chk = st.ejc_check(sn, ph, 77.7, MB, S, 0.5, 0, 1); put(6, ph, chk, st.accept_text(chk, 1e-10, 1e-8))
    sn = -98.0041250683277; ph = st.ejc_phi(sn, 77.7, MB, S, 0.5, 0, 1); chk = st.ejc_check(sn, ph, 77.7, MB, S, 0.5, 0, 1); put(7, ph, chk, st.accept_text(chk, 1e-10, 1e-8))
    put(8, st.ejc_phi(400, 77.7, MB, S, 0.5, -1, 1 / 0.271), None, "-")
    put(9, st.ejc_phi(400, 0, MB, S, 0.5, 0, 1), None, "-")
    ph = st.ejc_phi(T7SN, *V, -1, 20, 1e-10, 1); chk = st.ejc_check(T7SN, ph, *V, -1, 20); put(10, "found: True", chk, st.accept_text(chk, 1e-10, 1e-8))
    chk = st.ejc_check_state(400, c1, ph1, 77.7, MB, S, 0.5, 0, 1); put(11, st.accept_text(chk, 1e-10, 1e-8), chk, st.accept_text(chk, 1e-10, 1e-8))
    chk = st.ejc_check_state(400, c1 + 5, ph1, 77.7, MB, S, 0.5, 0, 1); put(12, st.accept_text(chk, 1e-10, 1e-8), chk, st.accept_text(chk, 1e-10, 1e-8))
    put(13, st.ejc_phi(T7SN, *V, -1, 20, 1e-10, 0), None, "-")
    put(14, st.ejc_phi(T7SN, *V, -1, 20, 1e-10, -1), None, "-")
    ph = st.ejc_phi(T7SN, *V, -1, 20, 1e-10, 2); chk = st.ejc_check(T7SN, ph, *V, -1, 20); put(15, ph, chk, st.accept_text(chk, 1e-10, 1e-8))
    put(16, st.ejc_phi(T7SN, *V, -1, 20, 0, 1), None, "-")
    put(17, st.ejc_phi(T7SN, *V, -1, 20, -1e-10, 1), None, "-")
    chk = st.ejc_check_state(29819.0191175968, 25010.6868294027, -0.012480264998, *V, -1, 20); lbl = st.accept_text(chk, 1e-10, 1e-8); put(18, lbl, chk, lbl)
    chk = st.ejc_check(T7SN, -45, *V, -1, 20); lbl = st.accept_text(chk, 1e-10, 1e-8); put(19, lbl, chk, lbl)
    chk = st.ejc_check(T7SN, 90, *V, -1, 20); lbl = st.accept_text(chk, 1e-10, 1e-8); put(20, lbl, chk, lbl)
    for row, (ta, te) in zip((21, 22, 23, 24), ((1e-10, 0), (1e-10, -1e-8), (0, 1e-8), (-1e-10, 1e-8))):
        lbl = st.ejc_accept_state(400, c1, ph1, 77.7, MB, S, 0.5, 0, 1, ta, te); put(row, lbl, None, "-")
    ph = st.ejc_phi(-49.768342, 77.7, MB, S, 0.509923, 0, 1); chk = st.ejc_check(-49.768342, ph, 77.7, MB, S, 0.509923, 0, 1); put(28, ph, chk, st.accept_text(chk, 1e-10, 1e-8))

    for row, (value, label, ang, ejc) in cases.items():
        name = str(ws.cell(row, 1).value).split()[0]
        src = f"VBA Validation!C{row}"
        rec = ws.cell(row, 3).value
        if isinstance(value, float):
            r = fnum(rec)
            rep.num(src, f"{name} returned phi [deg]", r, value, tol=1e-9)
        else:
            rep.text(src, f"{name} returned", rec, getattr(value, "value", value))
        rep.text(f"VBA Validation!F{row}", f"{name} acceptance label", ws.cell(row, 6).value, label)
        for col, ours, nm in ((4, ang, "angular residual"), (5, ejc, "EJC residual")):
            rec_t = ws.cell(row, col).value
            r = fnum(rec_t)
            if r is None:
                if ours is not None:
                    rep.add(f"VBA Validation!{'DE'[col - 4]}{row}", f"{name} {nm}", rec_t, ours, False)
                continue
            # residuals are recorded as 4-significant-digit text; tiny values only need to be tiny
            ok = abs(ours - r) <= 1e-9 if abs(r) < 1e-6 else abs(ours - r) <= 6e-4 * abs(r)
            rep.add(f"VBA Validation!{'DE'[col - 4]}{row}", f"{name} {nm}", r, ours, ok)


# ------------------------------------------------------------------------------------------------
def run(cfg: Config, dim: int, n: int):
    return solve_model(cfg, dim, n)


def f_or_withheld(run_) -> object:
    return run_.result.fos if run_.result.fos is not None else "WITHHELD"


def case_table():
    """name used in the sheets -> (Config, extra) for the engine-QA style cases."""
    c = lambda **k: Config(closure=3, **k)   # noqa: E731
    wet = dict(water_angle=50.0, water_height=360.0)
    return {
        "thesis": Config(),
        "dry5": c(), "wet5": c(**wet), "dry509": c(a_override=0.509923), "wet509": c(a_override=0.509923, **wet),
        "pub": c(xc=207.28, zc=586.53, radius=622.08), "mc": c(strength=0, c_mc=100.0, phi_mc=0.0),
        "mosc_dry": Config(closure=2),
        "mosc_wet_off": Config(closure=2, a_override=0.509923, checks=0, **wet),
        "mcgeom": c(xc=-1.45, zc=3.24, radius=4.46, height=2.1, slope=1 / 3.5, gamma=19.0, strength=0, c_mc=7.2, phi_mc=30.0),
        "cropped": c(xc=207.28, zc=300.0, radius=622.08),
    }


def _compare_pair(rep, src_row, label, wsname, recorded_cells, runs, expected_status):
    """recorded_cells: dict col -> value from (B,C,D,E,F,G); runs: (run3d, run2d)."""
    r3, r2 = runs
    for tag, run_, ce, ca, cs in (("3D", r3, "B", "C", "F"), ("2D", r2, "D", "E", "G")):
        res = run_.result
        exp, act, status = recorded_cells[ce], recorded_cells[ca], recorded_cells[cs]
        src = f"{wsname}!{ce}{src_row}"
        if isinstance(exp, (int, float)) and expected_status[tag] == "ACCEPTED":
            rep.num(src, f"{label} {tag} F vs recorded expected", exp, res.fos)
            rep.num(f"{wsname}!{ca}{src_row}", f"{label} {tag} F vs recorded actual", act, res.fos)
            rep.text(f"{wsname}!{cs}{src_row}", f"{label} {tag} status", status, res.status.value)
        else:
            rep.text(src, f"{label} {tag} F withheld", str(exp).upper(), "WITHHELD" if res.fos is None else res.fos)
            rep.text(f"{wsname}!{cs}{src_row}", f"{label} {tag} status", status, res.status.value)


def check_engine_qa(ws, rep: Report) -> None:
    t = case_table()
    plan = [  # row, key, 3D expected status, 2D expected status
        (4, "thesis"), (5, "dry5"), (6, "wet5"), (7, "dry509"), (8, "wet509"), (9, "pub"), (10, "mc"),
        (11, "mosc_dry"), (12, "wet509"), (13, "mcgeom"), (14, "cropped")]
    for row, key in plan:
        cells = {c: ws[f"{c}{row}"].value for c in "ABCDEFGH"}
        cfg = t[key]
        r3, r2 = run(cfg, 3, 20), run(cfg, 2, 40)
        label = str(cells["A"])
        if row == 12:      # perturbation: +1 kPa on the first column's actual final normal stress
            r3.elements.final_n[0] += 1.0
            still = validate_final(cfg, r3.elements, r3.result.internal_f, r3.result, rebuild=False)
            rep.add(f"{ws.title}!F{row}", f"{label} 3D rejected", "NOT ACCEPTED (PERTURBED)", "rejected" if not still else "ACCEPTED", not still)
            rep.num(f"{ws.title}!D{row}", f"{label} 2D F (expected)", cells["D"], r2.result.fos)
            rep.num(f"{ws.title}!E{row}", f"{label} 2D F (actual)", cells["E"], r2.result.fos)
            rep.text(f"{ws.title}!G{row}", f"{label} 2D status", cells["G"], r2.result.status.value)
            continue
        e3 = "ACCEPTED" if isinstance(cells["B"], (int, float)) else "NOT"
        e2 = "ACCEPTED" if isinstance(cells["D"], (int, float)) else "NOT"
        _compare_pair(rep, row, label, ws.title, cells, (r3, r2), {"3D": e3, "2D": e2})
        if row >= 13:
            ok = r3.result.active == 0 and r2.result.active == 0 and r3.result.reason.startswith("UNSUPPORTED GEOMETRY:")
            rep.add(f"{ws.title}!A{row}", f"{label} refused with zero active elements", "0 elements", r3.result.active, ok)


def check_workbook_qa(ws, rep: Report) -> None:
    t = case_table()
    plan = [(4, "thesis"), (5, "dry5"), (6, "wet5"), (7, "dry509"), (8, "wet509"), (9, "pub"), (10, "mc"),
            (11, "mosc_wet_off"), (12, "wet509")]
    for row, key in plan:
        cells = {c: ws[f"{c}{row}"].value for c in "ABCDEFGL"}
        cfg = t[key]
        label = str(cells["A"])
        r3, r2 = run(cfg, 3, 20), run(cfg, 2, 40)
        if row == 11:
            # Diagnostic mode (checks = 0): internal F only, status DIAGNOSTIC ONLY, never a displayed FoS.
            for tag, ce, ca, cs, run_ in (("3D", "B", "C", "F", r3), ("2D", "D", "E", "G", r2)):
                res = run_.result
                rep.num(f"{ws.title}!{ce}{row}", f"{label} {tag} internal F vs expected", cells[ce], res.internal_f)
                rep.num(f"{ws.title}!{ca}{row}", f"{label} {tag} internal F vs worksheet actual", cells[ca], res.internal_f)
                rep.add(f"{ws.title}!{cs}{row}", f"{label} {tag} status", cells[cs],
                        res.status.value, str(cells[cs]).startswith(res.status.value) and res.fos is None)
            continue
        if row == 12:
            r3.elements.final_n[0] += 1.0
            still = validate_final(cfg, r3.elements, r3.result.internal_f, r3.result, rebuild=False)
            rep.add(f"{ws.title}!F{row}", f"{label} 3D rejected", cells["F"], "rejected" if not still else "ACCEPTED", not still)
            rep.num(f"{ws.title}!D{row}", f"{label} 2D F (expected)", cells["D"], r2.result.fos)
            rep.num(f"{ws.title}!E{row}", f"{label} 2D F (worksheet engine actual)", cells["E"], r2.result.fos)
            continue
        e3 = "ACCEPTED" if isinstance(cells["B"], (int, float)) else "NOT"
        e2 = "ACCEPTED" if isinstance(cells["D"], (int, float)) else "NOT"
        if row == 4:
            rep.text(f"{ws.title}!F{row}", f"{label} 3D status", cells["F"], r3.result.status.value)
            rep.text(f"{ws.title}!G{row}", f"{label} 2D status", cells["G"], r2.result.status.value)
            rep.add(f"{ws.title}!B{row}", f"{label} F withheld", "withheld", (r3.result.fos, r2.result.fos),
                    r3.result.fos is None and r2.result.fos is None)
            continue
        _compare_pair(rep, row, label, ws.title, cells, (r3, r2), {"3D": e3, "2D": e2})


def check_last_run_tables(wb, rep: Report) -> None:
    """VBA Columns / VBA Slices: the element states of the last primary run (Chile, closure 1, NOT ACCEPTED)."""
    cfg = Config()
    for sheet, dim, n in (("VBA Columns", 3, 20), ("VBA Slices", 2, 40)):
        ws = wb[sheet]
        run_ = solve_model(cfg, dim, n)
        tab = elements_table(run_, cfg)
        rep.text(f"{sheet}!A2", "status line", ws["A2"].value, f"{run_.result.status.value}: {run_.result.reason}")
        header = [ws.cell(4, c).value for c in range(1, 24)]
        rep.text(f"{sheet}!A4:W4", "column headers", header, list(tab.columns))
        rows = [r for r in range(5, ws.max_row + 1) if ws.cell(r, 1).value is not None]
        rep.add(f"{sheet}!A5:A{rows[-1] if rows else 4}", "number of elements", len(rows), len(tab), len(rows) == len(tab))
        worst = 0.0
        bad = None
        for i, r in enumerate(rows[:len(tab)]):
            for c, col in enumerate(ELEMENT_COLUMNS, start=1):
                rec, ours = ws.cell(r, c).value, tab.iloc[i, c - 1]
                if isinstance(rec, (int, float)):
                    d = abs(rec - ours) / max(abs(rec), 1e-12)
                    worst = max(worst, d)
                    if d > rep.tol and bad is None:
                        bad = (r, col, rec, ours)
                elif rec is not None and str(rec) != str(ours) and not (rec == "NOT EVALUATED" and ours == "NOT EVALUATED"):
                    bad = bad or (r, col, rec, ours)
        rep.rels.append((f"{sheet}!A5:W", "element table cells", worst))
        rep.add(f"{sheet}!A5:W{rows[-1] if rows else 4}", f"all {len(rows)} x 23 cells (worst relative difference {worst:.1e})",
                "recorded", "recomputed", bad is None and worst <= rep.tol)
        if bad:
            rep.add(f"{sheet}!{bad[0]}", f"first mismatch in column '{bad[1]}'", bad[2], bad[3], False)
    run_ws = wb["VBA Run"]
    for cell, expected in (("B12", 76), ("B19", 30)):
        rep.text(f"VBA Run!{cell}", "active elements", run_ws[cell].value, expected)
    rep.text("VBA Run!B13", "3D explanation", run_ws["B13"].value, solve_model(cfg, 3, 20).result.reason)
    rep.text("VBA Run!B20", "2D explanation", run_ws["B20"].value, solve_model(cfg, 2, 40).result.reason)
    # known workbook anomaly, reported not hidden: B18 caches the number 3 although D18 and the status say WITHHELD
    b18, d18 = run_ws["B18"].value, run_ws["D18"].value
    if b18 != d18:
        rep.add("VBA Run!B18", "2D FoS display cell (workbook cache anomaly)", f"{b18!r} (D18 = {d18!r})", "WITHHELD", False, "open")


def check_slip_exit(wb, rep: Report) -> None:
    """Interface!F39:G39 holds the workbook's own chart helper for the slip-surface exit on the face."""
    from slope3d.geometry import slip_surface_ends
    ws = wb["Interface"]
    e = slip_surface_ends(Config())
    rep.num("Interface!F39", "slip-surface exit x [m]", ws["F39"].value, e["exit_x"])
    rep.num("Interface!G39", "slip-surface exit z above toe [m]", ws["G39"].value, e["exit_z"])


def check_refined_cases(ws, rep: Report) -> None:
    from slope3d.mesh_study import CASES, case_config
    names = {c[0]: c for c in CASES}
    for r in range(5, 15):
        cond, dim, mesh, f_rec = ws.cell(r, 1).value, ws.cell(r, 2).value, ws.cell(r, 3).value, ws.cell(r, 4).value
        _, surface, a, wet = names[cond]
        n = int(str(mesh).split("x")[0].split()[0])
        res = solve_model(case_config(surface, a, wet), 3 if dim == "3D" else 2, n).result
        rep.num(f"Refined Cases!D{r}", f"{cond} {dim} {mesh} finest F", f_rec, res.fos)


def print_summary(rep: Report) -> None:
    """Markdown table: per workbook sheet, comparisons / reproduced / failed / open and the worst relative difference."""
    groups: dict = {}
    for source, what, rec, ours, verdict, _ in rep.rows:
        g = source.split("!")[0]
        d = groups.setdefault(g, dict(n=0, ok=0, bad=0, open=0, worst=None, where=""))
        d["n"] += 1
        d["ok"] += verdict == "PASS"
        d["bad"] += verdict == "FAIL"
        d["open"] += verdict == "OPEN"
    for source, what, rel in rep.rels:
        d = groups[source.split("!")[0]]
        if d["worst"] is None or rel > d["worst"]:
            d["worst"], d["where"] = rel, f"{source} {what}"
    print("| Workbook sheet | Comparisons | Reproduced | Failed | Open | Worst relative difference (numeric) | Where |")
    print("|---|---|---|---|---|---|---|")
    for g, d in groups.items():
        w = "n.a." if d["worst"] is None else f"{d['worst']:.1e}"
        print(f"| {g} | {d['n']} | {d['ok']} | {d['bad']} | {d['open']} | {w} | {d['where']} |")
    tot = {k: sum(d[k] for d in groups.values()) for k in ("n", "ok", "bad", "open")}
    worst = max((r[2] for r in rep.rels), default=0.0)
    print(f"| **All** | **{tot['n']}** | **{tot['ok']}** | **{tot['bad']}** | **{tot['open']}** | **{worst:.1e}** | |")


# ------------------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workbook", nargs="?", default=str(DEFAULT_WB))
    ap.add_argument("--tol", type=float, default=1e-9, help="relative tolerance (default 1e-9)")
    ap.add_argument("--skip-mesh", action="store_true", help="skip the 3D 320x320 / 2D 2560-slice comparison (~20 s)")
    ap.add_argument("--strict", action="store_true", help="count OPEN items as failures")
    ap.add_argument("--summary", action="store_true", help="print a Markdown table of counts and worst relative error per sheet")
    ap.add_argument("--quiet", action="store_true", help="print only failures, open items and the summary")
    args = ap.parse_args(argv)

    wb = openpyxl.load_workbook(args.workbook, data_only=True)
    rep = Report(args.tol)
    check_strength_fixtures(wb["VBA Validation"], rep)
    check_workbook_qa(wb["Native Workbook QA"], rep)
    check_engine_qa(wb["Native Engine QA"], rep)
    check_last_run_tables(wb, rep)
    check_slip_exit(wb, rep)
    if not args.skip_mesh:
        check_refined_cases(wb["Refined Cases"], rep)

    if args.summary:
        print_summary(rep)
    n_pass = sum(r[4] == "PASS" for r in rep.rows)
    n_fail = sum(r[4] == "FAIL" for r in rep.rows)
    n_open = sum(r[4] == "OPEN" for r in rep.rows)
    for source, what, rec, ours, verdict, _ in rep.rows:
        if args.quiet and verdict == "PASS":
            continue
        fmt = lambda v: f"{v:.15g}" if isinstance(v, float) else str(v)   # noqa: E731
        print(f"{verdict:4s}  {source:28s} {what:62s} workbook={fmt(rec)}  slope3d={fmt(ours)}")
    print(f"\n{n_pass} of {len(rep.rows)} comparisons reproduced (tolerance {args.tol:g}); {n_fail} failed; {n_open} open")
    if n_open:
        print("OPEN items are workbook records that slope3d does not reproduce; see docs/SPEC.md section 0.")
    bad = n_fail + (n_open if args.strict else 0)
    print("ALL REPRODUCED" if bad == 0 and n_open == 0 else ("NO FAILURES (open items remain)" if bad == 0 else "MISMATCHES FOUND"))
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
