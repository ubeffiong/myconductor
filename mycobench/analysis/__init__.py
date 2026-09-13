"""Analysis layer: the Class A methods that need only data already collected.

Scope, stated once so it is not mistaken
----------------------------------------
This package implements the parts of the gap-closure plan whose blocker was
analysis rather than laboratory work:

``mic``        MIC on the doubling-dilution scale, with censoring that stays
               censored rather than being silently read as a measurement.
``strata``     Resistance-background stratification — the move that separates
               a causal variant from one merely linked to it.
``effects``    Background-conditional effect sizes with intervals, replacing
               marginal association.
``selective``  Risk-coverage evaluation, so a predictor that declines to answer
               is scored on what it actually claims.
``mechanism``  Mechanism class inferred from MIC distribution shape and
               cross-drug covariance, rather than from unmeasured expression.
``barrier``    Target ranking by resistance barrier rather than essentiality.

What this package deliberately does not contain
-----------------------------------------------
A trained resistance classifier. Building one before the evaluation harness
exists is the failure this whole codebase was rebuilt to remove: prediction
machinery outrunning any means of checking it. ``selective`` and ``effects``
are the means of checking; a model is welcome afterwards and is not welcome
before.

Nothing here needs a laboratory, and nothing here can substitute for one. The
gaps that need expression measurement, dilution series or paired African DST
are untouched by this package and remain open.
"""
from __future__ import annotations

__all__ = ["mic", "strata", "effects", "selective", "mechanism", "barrier",
           "stats"]
