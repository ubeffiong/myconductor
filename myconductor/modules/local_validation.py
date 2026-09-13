"""Local validation lane — site-confirmed VUS results as evidence.

Reads from a ``federated.vus_feedback.LocalValidationStore``. A variant this
site has already validated is no longer just a laboratory ranking; it is a
real, local, phenotype-derived answer, so it is evidence at ``Tier.PHENOTYPIC``
— the tier a laboratory DST result is entitled to. Two things keep this from
overreaching:

* a **conflicting** local verdict (some isolates resistant, some susceptible)
  is reported as ``INDETERMINATE``, never resolved by picking a side — the
  same "disagreement is a finding" rule every other lane follows;
* this is deliberately **not** the global catalogue. It carries a limitation
  saying so on every piece of evidence it emits, and nothing here writes to
  ``federated.catalogue_update``'s review queue — that promotion path, with
  its confounding checks and expert curation, is unchanged and unbypassed.
"""
from __future__ import annotations

from ..core.models import Call, DrugEvidence, Lane, Tier, Variant
from ..federated.vus_feedback import LocalValidationStore
from .base import VariantModule


class LocalValidationModule(VariantModule):
    name = "local_validation"
    lane = Lane.PHENOTYPE

    def __init__(self, store: LocalValidationStore):
        self.store = store

    def applies_to(self, variant: Variant) -> bool:
        return bool(self.store.drugs_for(variant.key()))

    def evaluate(self, variant: Variant) -> list[DrugEvidence]:
        out: list[DrugEvidence] = []
        for drug in self.store.drugs_for(variant.key()):
            verdict = self.store.verdict(variant.key(), drug)
            if verdict is None:
                continue
            sites = sorted({r.site_id for r in verdict.records})
            base_limitations = (
                "site-local validation, not a global catalogue entry: not "
                "yet promoted through federated.catalogue_update's review "
                "queue",
                f"{verdict.describe()}",
            )
            if verdict.consistent:
                out.append(DrugEvidence(
                    drug=drug,
                    call=verdict.call,
                    tier=Tier.PHENOTYPIC,
                    lane=Lane.PHENOTYPE,
                    confidence=None,
                    variant=variant.identity,
                    rule_id=f"local-validation:{'+'.join(sites)}:"
                            f"{variant.key()}:{drug}",
                    limitations=base_limitations,
                    rationale=(
                        f"Locally validated {verdict.call.value} for {drug}: "
                        f"{verdict.describe()}."
                    ),
                ))
            else:
                out.append(DrugEvidence(
                    drug=drug,
                    call=Call.INDETERMINATE,
                    tier=Tier.PHENOTYPIC,
                    lane=Lane.PHENOTYPE,
                    confidence=None,
                    variant=variant.identity,
                    rule_id=f"local-validation:{'+'.join(sites)}:"
                            f"{variant.key()}:{drug}",
                    limitations=base_limitations + (
                        "conflicting local results are not resolved by "
                        "picking a side",
                    ),
                    rationale=(
                        f"Conflicting local validation results for {drug}: "
                        f"{verdict.describe()}. Susceptibility withheld; "
                        f"resistance not asserted."
                    ),
                ))
        return out
