"""The lane router.

What changed
------------
The original router assigned each variant to exactly one lane, in priority
order, and only that lane was allowed to speak. That discards evidence by
construction: a variant can simultaneously be catalogued, sit in an efflux
regulator, be structurally interesting and be present as a minority allele.
Those are not competing classifications — they are facts that co-exist, and a
reconciliation layer needs all of them.

So the router now asks every lane whether it *applies*, runs all that do, and
collects everything they return. It decides which analyses are applicable, not
which single analysis is permitted to speak.

The router still owns no biology. It dispatches and counts.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

from ..modules.base import VariantModule
from ..modules.catalogue import CatalogueModule
from ..modules.efflux import EffluxRegulatoryModule
from .models import DrugEvidence, Lane, Variant


@dataclass
class RoutingDecision:
    """Which lanes examined one variant, and what they produced."""

    variant: Variant
    lanes: tuple[Lane, ...]
    evidence: list[DrugEvidence] = field(default_factory=list)

    @property
    def examined(self) -> bool:
        return bool(self.lanes)

    @property
    def spoke(self) -> tuple[Lane, ...]:
        """Lanes that actually returned evidence, as opposed to abstaining."""
        return tuple(dict.fromkeys(ev.lane for ev in self.evidence))


@dataclass
class RoutingOutcome:
    decisions: list[RoutingDecision] = field(default_factory=list)

    @property
    def evidence(self) -> list[DrugEvidence]:
        return [ev for d in self.decisions for ev in d.evidence]

    @property
    def unexamined(self) -> list[Variant]:
        """Variants no lane could speak to — reported, never silently dropped."""
        return [d.variant for d in self.decisions if not d.examined]

    def lane_counts(self) -> dict[str, int]:
        """How many variants each lane *examined*.

        Sums to more than the variant count when lanes overlap, which is the
        intended behaviour and the point of multi-lane routing.
        """
        counts = {lane.value: 0 for lane in Lane if lane is not Lane.NONE}
        counts["unexamined"] = 0
        for d in self.decisions:
            if not d.examined:
                counts["unexamined"] += 1
            for lane in d.lanes:
                counts[lane.value] = counts.get(lane.value, 0) + 1
        return counts


class TriageRouter:
    def __init__(
        self,
        catalogue: Optional[CatalogueModule] = None,
        efflux: Optional[EffluxRegulatoryModule] = None,
        extra_lanes: Sequence[VariantModule] = (),
    ):
        self.catalogue = catalogue or CatalogueModule()
        self.efflux = efflux or EffluxRegulatoryModule()
        #: Ordered for reporting stability only; all applicable lanes run.
        self.lanes: list[VariantModule] = [self.catalogue, self.efflux,
                                           *extra_lanes]

    def applicable_lanes(self, variant: Variant) -> list[VariantModule]:
        """Every lane that has something to say about this variant.

        Silent variants are excluded wholesale: a synonymous change is the one
        case where no lane applies for a structural reason rather than a lack
        of data.
        """
        if variant.silent:
            return []
        return [lane for lane in self.lanes if lane.applies_to(variant)]

    def route(self, variants: list[Variant]) -> RoutingOutcome:
        outcome = RoutingOutcome()
        for variant in variants:
            lanes = self.applicable_lanes(variant)
            evidence: list[DrugEvidence] = []
            for lane in lanes:
                evidence.extend(lane.evaluate(variant))
            outcome.decisions.append(RoutingDecision(
                variant=variant,
                lanes=tuple(lane.lane for lane in lanes),
                evidence=evidence,
            ))
        return outcome
