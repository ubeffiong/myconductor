"""Efflux / regulatory lane.

Two blind spots a catalogue lookup shares: efflux-pump de-repression, where the
mere presence of a variant in a regulator says nothing about whether the pump is
actually overexpressed; and promoter variants that modulate expression of a drug
target.

This lane reports **INDETERMINATE**, never RESISTANT
----------------------------------------------------
The original implementation returned ``PREDICTED_RESISTANT`` for any non-silent
change in ``Rv0678``. That overclaims: a hand-written rule over gene membership
is not graded evidence, and without expression data (RNA, or a validated
regulatory model) de-repression is a hypothesis rather than a finding.

``INDETERMINATE`` is both honest and clinically protective. It withholds
susceptibility — so the drug cannot count toward a regimen — without asserting
a mechanism nobody measured. ``Tier.INFERRED`` is structurally barred from
emitting RESISTANT by ``DrugEvidence.__post_init__``.

And it emits one piece of evidence **per affected drug**: ``Rv0678`` raises
bedaquiline and clofazimine MICs together, and reporting only the first drug in
the list silently loses the second.
"""
from __future__ import annotations

from typing import Optional

from ..catalogue.profile import OrganismProfile, load_profile
from ..core.models import (
    Call,
    Consequence,
    DrugEvidence,
    EngineRef,
    Lane,
    Region,
    Tier,
    Variant,
)
from .base import VariantModule

_ENGINE = EngineRef(
    name="myconductor-efflux-rules",
    version="0.2.0",
    database="drug_loci.json (efflux_regulators, promoter_loci)",
)

_NO_EXPRESSION_DATA = (
    "no expression evidence: de-repression is inferred from gene membership "
    "alone, not measured"
)
_RULE_BASED = (
    "rule-based inference, not graded catalogue evidence and not a validated "
    "predictive model"
)
_LOF_NOT_GRADED = (
    "loss of function is the strongest inference this lane can make, but it is "
    "still not a graded catalogue entry; an ingested WHO catalogue grades many "
    "regulator variants directly and should be preferred where it does"
)
_MISSENSE_WEAKER = (
    "a missense change in a regulator carries materially weaker evidence than "
    "loss of function; prioritise it for validation rather than acting on it"
)


class EffluxRegulatoryModule(VariantModule):
    name = "efflux_regulatory"
    lane = Lane.EFFLUX_REGULATORY

    def __init__(self, profile: Optional[OrganismProfile] = None):
        self.profile = profile or load_profile()

    def applies_to(self, variant: Variant) -> bool:
        if variant.silent:
            return False
        if variant.gene in self.profile.efflux_genes:
            return True
        return (variant.region in (Region.PROMOTER, Region.INTERGENIC)
                and bool(self.profile.drugs_for_promoter(variant.gene)))

    def evaluate(self, variant: Variant) -> list[DrugEvidence]:
        if variant.silent:
            return []

        evidence: list[DrugEvidence] = []

        # Efflux regulator. Loss of function and missense are separated
        # because the published evidence for them is not comparable: repeated,
        # independent truncation of mmpR5/Rv0678 accounts for the large
        # majority of described bedaquiline resistance, whereas a missense
        # change there is a variant of unknown significance like any other.
        # Treating both identically — as this lane previously did — discards
        # the strongest signal available for the drug that matters most.
        drugs = self.profile.drugs_for_efflux_regulator(variant.gene)
        if drugs:
            truncating = variant.consequence.is_truncating
            for drug in drugs:
                if truncating:
                    limitations = (_NO_EXPRESSION_DATA, _RULE_BASED,
                                   _LOF_NOT_GRADED)
                    rationale = (
                        f"Loss of function ({variant.consequence.value}) in "
                        f"efflux regulator {variant.gene}: MmpS5-MmpL5 "
                        f"de-repression would raise {drug} MIC, and repeated "
                        f"independent truncation of this regulator is the "
                        f"commonest described route to bedaquiline resistance. "
                        f"Susceptibility withheld; resistance NOT asserted, "
                        f"because this lane holds no graded evidence."
                    )
                else:
                    limitations = (_NO_EXPRESSION_DATA, _RULE_BASED,
                                   _MISSENSE_WEAKER)
                    rationale = (
                        f"Missense change in efflux regulator {variant.gene}: "
                        f"de-repression is plausible but unquantified, and "
                        f"carries materially weaker evidence than loss of "
                        f"function. Susceptibility withheld pending "
                        f"phenotypic testing."
                    )
                evidence.append(DrugEvidence(
                    drug=drug,
                    call=Call.INDETERMINATE,
                    tier=Tier.INFERRED,
                    lane=Lane.EFFLUX_REGULATORY,
                    confidence=None,
                    variant=variant.identity,
                    engine=_ENGINE,
                    limitations=limitations,
                    rationale=rationale,
                ))

        # Promoter / intergenic variant modulating target expression.
        if variant.region in (Region.PROMOTER, Region.INTERGENIC):
            for drug in self.profile.drugs_for_promoter(variant.gene):
                evidence.append(DrugEvidence(
                    drug=drug,
                    call=Call.INDETERMINATE,
                    tier=Tier.INFERRED,
                    lane=Lane.EFFLUX_REGULATORY,
                    confidence=None,
                    variant=variant.identity,
                    engine=_ENGINE,
                    limitations=(_NO_EXPRESSION_DATA, _RULE_BASED),
                    rationale=(
                        f"Regulatory variant upstream of {variant.gene} "
                        f"({variant.change}): may modulate expression affecting "
                        f"{drug}. Effect size unmeasured — susceptibility "
                        f"withheld."
                    ),
                ))

        return _dedupe(evidence)


def _dedupe(evidence: list[DrugEvidence]) -> list[DrugEvidence]:
    """A gene can be both a regulator and a promoter locus (e.g. pepQ).

    Keep one piece of evidence per drug so the report does not double-count the
    same inference.
    """
    seen: set[str] = set()
    out: list[DrugEvidence] = []
    for ev in evidence:
        if ev.drug in seen:
            continue
        seen.add(ev.drug)
        out.append(ev)
    return out
