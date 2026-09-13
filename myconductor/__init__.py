"""Myconductor -- an orchestration-and-reasoning layer for TB AMR.

Not another variant profiler: a conductor that ingests from the existing tools,
routes every finding to the right specialist lane (catalogue lookup, VUS ML
classifier, or efflux/regulatory model), synthesises a regimen, and -- when the
strain defeats the whole arsenal -- automatically escalates to a drug-discovery
loop. A federated living catalogue lets it learn without moving patient genomes.

Quick start
-----------
    from myconductor import Myconductor
    from myconductor.reporting.render import render_text

    report = Myconductor().analyze("myconductor/data/example_input.vcf")
    print(render_text(report))
"""
from __future__ import annotations

__version__ = "0.1.0"

from .core.pipeline import Myconductor  # noqa: E402

__all__ = ["Myconductor", "__version__"]
