"""The triage router -- the brain of Myconductor.

This is the component that does not exist in today's tools. For every variant it
decides *which lane* should judge it, then collects that lane's evidence. The
routing rule, in order:

    1. Silent / off-target          -> IGNORED
    2. Known to the catalogue       -> CATALOGUE lane
    3. Efflux gene or promoter/non-coding -> EFFLUX_REGULATORY lane
    4. Coding, in a modelled gene   -> VUS_ML lane
    5. Otherwise                    -> IGNORED (nothing can speak to it yet)

The router owns no biology of its own; it only dispatches. That keeps the
decision logic auditable and lets each lane evolve independently.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from ..modules.catalogue import CatalogueModule
from ..modules.efflux import EffluxRegulatoryModule
from ..modules.vus_classifier import VUSClassifierModule
from .models import DrugEvidence, Region, Route, Variant


@dataclass
class RoutingDecision:
    variant: Variant
    route: Route
    evidence: Optional[DrugEvidence]


@dataclass
class RoutingOutcome:
    decisions: list[RoutingDecision] = field(default_factory=list)

    @property
    def evidence(self) -> list[DrugEvidence]:
        return [d.evidence for d in self.decisions if d.evidence is not None]

    def counts(self) -> dict[str, int]:
        c: dict[str, int] = {r.value: 0 for r in Route}
        for d in self.decisions:
            c[d.route.value] += 1
        return c


class TriageRouter:
    def __init__(
        self,
        catalogue: Optional[CatalogueModule] = None,
        vus: Optional[VUSClassifierModule] = None,
        efflux: Optional[EffluxRegulatoryModule] = None,
    ):
        self.catalogue = catalogue or CatalogueModule()
        self.vus = vus or VUSClassifierModule()
        self.efflux = efflux or EffluxRegulatoryModule()

    def _route_of(self, variant: Variant) -> Route:
        if variant.silent:
            return Route.IGNORED
        if self.catalogue.knows(variant):
            return Route.CATALOGUE
        if variant.gene in self.catalogue.efflux_genes:
            return Route.EFFLUX_REGULATORY
        if variant.region in (Region.PROMOTER, Region.INTERGENIC):
            return Route.EFFLUX_REGULATORY
        return Route.VUS_ML

    def _dispatch(self, variant: Variant, route: Route) -> Optional[DrugEvidence]:
        if route == Route.CATALOGUE:
            return self.catalogue.evaluate(variant)
        if route == Route.EFFLUX_REGULATORY:
            return self.efflux.evaluate(variant)
        if route == Route.VUS_ML:
            return self.vus.evaluate(variant)
        return None

    def route(self, variants: list[Variant]) -> RoutingOutcome:
        outcome = RoutingOutcome()
        for v in variants:
            route = self._route_of(v)
            evidence = self._dispatch(v, route)
            # If a lane abstains, the variant is effectively unrouteable.
            actual = route if (evidence is not None or route == Route.IGNORED) else Route.IGNORED
            outcome.decisions.append(RoutingDecision(v, actual, evidence))
        return outcome
