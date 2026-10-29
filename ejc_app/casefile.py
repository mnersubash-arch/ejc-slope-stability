"""Case files: the inputs of one analysis as JSON, shared by this app and the web app (format ``ejc-slope-case/1``).

A case holds the engine inputs (``inputs``, every field of engine.Inputs), the FoS when the results are for those
inputs (``results``), its name, the preset it started from, and optionally the page state of the app that wrote it
(``state`` for the web app, ``streamlit_state`` for this one), which keeps choices such as benches that the inputs
hold only as points. Each app reads its own page state when it is there and the engine inputs otherwise, so a case
saved in one app opens in the other with the same inputs. No Streamlit here, so the module is testable on its own.
"""
from __future__ import annotations

import base64
import json
import math
import re
import zlib
from dataclasses import asdict, fields
from datetime import datetime, timezone
from typing import Optional

import engine as E

CASE_FORMAT = "ejc-slope-case/1"
CHOICE = re.compile(r"^[A-Za-z0-9_.-]{1,40}$")
NULLABLE = {"a_override", "ru", "n3y", "plan_radius", "d_zone", "crack_depth"}
INTS = {"n3", "n2", "n3y"}
POINTS = {"profile", "piezo"}


def inputs_dict(inp: E.Inputs) -> dict:
    """The engine inputs as plain JSON values (point tuples become lists)."""
    d = asdict(inp)
    for k in POINTS:
        if d[k] is not None:
            d[k] = [[float(x), float(z)] for x, z in d[k]]
    return d


def make_case(inp: E.Inputs, *, name: str = "", preset: str = "", ui_state: Optional[dict] = None,
              results: Optional[dict] = None, saved: Optional[str] = None) -> dict:
    """A case for the inputs. ``results`` = {"r3": {"fos", "status"}, "r2": {...}} when the results are for these
    inputs; ``ui_state`` = this app's page state (st.session_state.inp)."""
    c = {"format": CASE_FORMAT, "name": name or "", "saved": saved or datetime.now(timezone.utc).isoformat(),
         "preset": preset or "", "inputs": inputs_dict(inp), "results": results}
    if ui_state is not None:
        c["streamlit_state"] = {k: v for k, v in ui_state.items() if isinstance(v, (bool, int, float, str)) or v is None}
    return json.loads(json.dumps(c, allow_nan=False))


def results_of(A) -> dict:
    """The FoS and status of an analysis, as stored in a case."""
    return {f"r{d}": {"fos": r.fos, "status": r.status} for d, r in ((3, A.r3), (2, A.r2))}


def _number(k: str, v):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ValueError(f"{k} must be a number")
    if k in INTS:
        if float(v) != int(v):
            raise ValueError(f"{k} must be a whole number")
        return int(v)
    return float(v)


def clean_inputs(raw: dict) -> dict:
    """Engine inputs from a case: known fields of the right type only (as cleanInputs in webapp/index.html).
    Raises ValueError naming the first field that is not."""
    if not isinstance(raw, dict):
        raise ValueError("the case holds no inputs")
    defaults = E.Inputs()
    out = {}
    for f in fields(E.Inputs):
        k = f.name
        if k not in raw:
            continue
        v, d = raw[k], getattr(defaults, k)
        if k in POINTS:
            if v is None:
                out[k] = None
            elif (isinstance(v, list) and len(v) <= 500 and all(isinstance(q, list) and len(q) == 2 for q in v)):
                out[k] = tuple((_number(k, q[0]), _number(k, q[1])) for q in v)
            else:
                raise ValueError(f"{k} must be a list of [x, z] points")
        elif isinstance(d, bool):
            if not isinstance(v, bool):
                raise ValueError(f"{k} must be true or false")
            out[k] = v
        elif isinstance(d, str):
            if not (isinstance(v, str) and CHOICE.match(v)):
                raise ValueError(f"{k} is not one of the app's choices")
            out[k] = v
        elif v is None and (d is None or k in NULLABLE):
            out[k] = None
        else:
            out[k] = _number(k, v)
    return out


def read_case(text) -> dict:
    """A case from JSON text (bytes or str). Raises ValueError if it is not one."""
    try:
        c = json.loads(text)
    except (ValueError, TypeError) as exc:
        raise ValueError("it is not valid JSON") from exc
    if not (isinstance(c, dict) and c.get("format") == CASE_FORMAT):
        raise ValueError(f"it is not a case saved by this app (format {CASE_FORMAT})")
    if not isinstance(c.get("inputs"), dict):
        raise ValueError("the case holds no inputs")
    return c


def inputs_from_case(c: dict) -> E.Inputs:
    return E.Inputs(**clean_inputs(c["inputs"]))


def case_json(c: dict) -> str:
    return json.dumps(c, indent=1, allow_nan=False)


def file_name(c: dict, kind: str = "case", ext: str = "json") -> str:
    """ejc-<kind>-<name>-<date>.<ext>, as the web app names its files."""
    slug = re.sub(r"[^a-z0-9]+", "-", str(c.get("name") or c.get("preset") or "case").lower()).strip("-")[:60] or "case"
    return f"ejc-{kind}-{slug}-{str(c.get('saved', ''))[:10]}.{ext}"


# a case in a link: the inputs, name and preset, compressed (the page state stays out to keep the link short)
def encode_link(c: dict) -> str:
    small = {k: c.get(k) for k in ("format", "name", "preset", "inputs")}
    raw = zlib.compress(json.dumps(small, separators=(",", ":"), allow_nan=False).encode(), 9)
    return base64.urlsafe_b64encode(raw).decode().rstrip("=")


def decode_link(code: str) -> dict:
    try:
        z = zlib.decompressobj()
        raw = z.decompress(base64.urlsafe_b64decode(code + "=" * (-len(code) % 4)), 1 << 20)   # at most 1 MiB
    except (ValueError, zlib.error) as exc:
        raise ValueError("the link does not hold a case") from exc
    if z.unconsumed_tail:
        raise ValueError("the link does not hold a case")
    return read_case(raw)
