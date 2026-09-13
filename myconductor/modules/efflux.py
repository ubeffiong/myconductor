"""Efflux / regulatory lane.

Two blind spots a catalogue lookup shares: (1) efflux-pump overexpression, where
the mere *presence* of a gene like Rv0678 tells you nothing about whether it is
actively overexpressed and driving low-to-moderate multi-drug resistance; and
(2) promoter / intergenic variants that modulate expression of a drug target.

This lane flags both. In production the "is it overexpressed" question is
answered with an expression-inference model (ideally RNA evidence, or a
regulatory-variant model). Here we apply a transparent rule over the regulatory
context so the lane is exercised end to end.
"""
from __future__ import annotations

from typing import Optional

from ..core.models import Call, DrugEvidence, Region, Route, Variant
from .base import VariantModule

# Loss-of-function in these regulators de-represses the MmpS5-MmpL5 efflux pump,
# raising bedaquiline and clofazimine MICs together.
_EFFLUX_GENE_DRUGS = {
    "Rv0678": ["bedaquiline", "clofazimine"],
    "mmpR5": ["bedaquiline", "clofazimine"],
    "pepQ": ["bedaquiline", "clofazimine"],
}

_PROMOTER_DRUGS = {
    "inhA": "isoniazid",
    "eis": "amikacin",
    "ahpC": "isoniazid",
    "fabG1": "isoniazid",
}


class EffluxRegulatoryModule(VariantModule):
    name = "efflux_regulatory"

    def evaluate(self, variant: Variant) -> Optional[DrugEvidence]:
        # Efflux regulator: any non-silent change is a de-repression candidate.
        if variant.gene in _EFFLUX_GENE_DRUGS and not variant.silent:
            drugs = _EFFLUX_GENE_DRUGS[variant.gene]
            return DrugEvidence(
                drug=drugs[0],
                call=Call.PREDICTED_RESISTANT,
                confidence=0.55,
                route=Route.EFFLUX_REGULATORY,
                variant_key=variant.key(),
                rationale=(
                    f"Change in efflux regulator {variant.gene}: likely MmpS5-MmpL5 "
                    f"de-repression -> low/moderate resistance across {', '.join(drugs)}."
                ),
            )

        # Promoter / intergenic variant that up-regulates a target.
        if variant.region in (Region.PROMOTER, Region.INTERGENIC):
            drug = _PROMOTER_DRUGS.get(variant.gene)
            if drug:
                return DrugEvidence(
                    drug=drug,
                    call=Call.PREDICTED_RESISTANT,
                    confidence=0.6,
                    route=Route.EFFLUX_REGULATORY,
                    variant_key=variant.key(),
                    rationale=(
                        f"Regulatory variant in {variant.gene} promoter: modulates "
                        f"expression affecting {drug}."
                    ),
                )
        return None
