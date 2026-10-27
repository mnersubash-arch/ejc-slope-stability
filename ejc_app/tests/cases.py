"""Shared, cached analyses (the engine is deterministic, so each case is solved once per session)."""
from functools import lru_cache

import engine as E


@lru_cache(maxsize=None)
def r3(inp: E.Inputs) -> E.LemResult:
    return E.run_3d(inp)


@lru_cache(maxsize=None)
def r2(inp: E.Inputs) -> E.LemResult:
    return E.run_2d(inp)


MAIN = E.Inputs()
WET = MAIN.with_(wet=True)
