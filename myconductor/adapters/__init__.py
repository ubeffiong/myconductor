"""Adapters for external AMR engines — the orchestration layer proper.

Myconductor's defensible contribution is not predicting resistance; it is
running validated engines, reconciling what they say, and reporting where they
disagree. These adapters are that boundary.

Validation status
-----------------
Every parser here is written from each tool's **documented output schema**. None
has been validated against real output from an installed tool, because this
repository ships no bioinformatics dependencies. Consequently each adapter
fails loudly on an unrecognised schema (``AdapterSchemaError``) rather than
silently mis-parsing — a wrong resistance call is far more costly than a crash.

Before relying on any adapter, run it against output from your own pinned
version of the tool and add a golden-file test.
"""
from __future__ import annotations

from .base import (
    AdapterSchemaError,
    ConcordanceSummary,
    EngineAdapter,
    EngineReport,
    concordance,
    cross_engine_discordance,
    merge,
)

__all__ = [
    "AdapterSchemaError",
    "ConcordanceSummary",
    "EngineAdapter",
    "EngineReport",
    "concordance",
    "cross_engine_discordance",
    "merge",
]
