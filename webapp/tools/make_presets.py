"""Write webapp/presets.js from ejc_app/presets.py (every preset with its inputs, source, targets and reference values).

Usage (from the repository root): python webapp/tools/make_presets.py [path to ejc_app]
Run it after changing ejc_app/presets.py or data/presets_law_desterre.json. The Law and D'Esterre presets carry their
source record (`data`); the file's `_notes` are not exported.
"""
import dataclasses
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parents[1] / "ejc_app"
sys.path.insert(0, str(APP))

import presets as PR  # noqa: E402


def num(v):
    if v is None or isinstance(v, (bool, str)):
        return v
    if isinstance(v, int):
        return v
    v = float(v)
    return v if math.isfinite(v) else None


presets = [dict(key=p.key, label=p.label, group=p.group,
                inputs={k: num(v) for k, v in dataclasses.asdict(p.inputs).items()}, source=p.source,
                targets=p.targets, reference=p.reference, slow=p.slow, action=p.action, info=p.info,
                data=p.data or None)
           for p in PR.PRESETS]
body = json.dumps(presets, ensure_ascii=False, separators=(",", ":"))
js = ("/* Presets of ejc_app/presets.py (this thesis, previous studies, Law (2024), D'Esterre (2024)), exported as data. */\n"
      "(function (root) {\n  \"use strict\";\n  const PRESETS = " + body + ";\n"
      "  if (typeof module !== \"undefined\" && module.exports) module.exports = PRESETS;\n"
      "  root.EJC_PRESETS = PRESETS;\n})(typeof window !== \"undefined\" ? window : globalThis);\n")
(HERE.parent / "presets.js").write_text(js)
print(f"{len(presets)} presets -> {HERE.parent / 'presets.js'}")
