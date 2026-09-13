"""Historical variant associations are research evidence, never another isolate's phenotype."""
from ..core.models import Call, DrugEvidence, Lane, Tier
from .base import VariantModule


class LocalValidationModule(VariantModule):
    name = "local_validation"
    lane = Lane.VUS

    def __init__(self, store):
        self.store = store

    def applies_to(self, variant):
        return bool(self.store.drugs_for(variant.key()))

    def evaluate(self, variant):
        out = []
        for drug in self.store.drugs_for(variant.key()):
            verdict = self.store.verdict(variant.key(), drug)
            if verdict is None:
                continue
            out.append(DrugEvidence(
                drug=drug, call=Call.INDETERMINATE, tier=Tier.INFERRED,
                lane=Lane.VUS, variant=variant.identity, scope="historical",
                rule_id=f"historical-association:{variant.key()}:{drug}",
                rationale=f"Historical observations: {verdict.describe()}. "
                          "An isolate phenotype does not establish variant causality.",
                limitations=("site-local validation is research context only; "
                             "not a phenotype measurement of the current isolate",),
                metadata={"records": [r.to_dict() for r in verdict.records]},
            ))
        return out
