"""Myconductor — an evidence-orchestration layer for AMR genomics.

Not a resistance predictor. Myconductor runs validated AMR engines, reconciles
what they say, keeps uncertainty as a first-class output, and governs how new
evidence changes an interpretation.

The rule everything else follows: absence of evidence is not susceptibility. A
drug is reported susceptible only when its required loci were shown to be
callable by evidence independent of the variant list, and only a
coverage-backed susceptible call counts toward regimen eligibility.

Status
------
Research scaffold. No component has been clinically validated; the bundled
catalogue is a 14-entry illustrative subset, not the WHO catalogue; the engine
adapters are written from documented schemas and have not been run against real
tool output. See README.md for what is and is not implemented.

Quick start
-----------
    from myconductor import Myconductor
    from myconductor.reporting.render import render_text

    report = Myconductor().analyze(
        "myconductor/data/example_input.vcf",
        mask_path="myconductor/data/example_callable.tsv",
    )
    print(render_text(report))
"""
from __future__ import annotations

__version__ = "0.2.0"

from .core.pipeline import Myconductor  # noqa: E402

__all__ = ["Myconductor", "__version__"]
