"""Lane interfaces.

Two changes from the original contract, both forced by biology:

``evaluate`` returns a **list**
    A single variant can affect several drugs. Loss of function in ``Rv0678``
    de-represses the MmpS5-MmpL5 pump and raises bedaquiline *and* clofazimine
    MICs together; an interface capped at one ``DrugEvidence`` silently drops
    the second drug. Returning a list makes that impossible.

``applies_to`` is separate from ``evaluate``
    The router asks every lane whether it *applies* to a variant, and a variant
    may be examined by several lanes at once — catalogued, sitting in an efflux
    regulator, and structurally interesting are not mutually exclusive facts.
    Exclusive routing throws away evidence by construction.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

from ..core.models import DrugEvidence, Lane, Variant


class VariantModule(ABC):
    """A lane that turns one variant into zero or more pieces of evidence."""

    name: str = "unnamed"
    lane: Lane = Lane.NONE

    def applies_to(self, variant: Variant) -> bool:
        """Should this lane look at this variant at all?

        Default: yes. Lanes that are cheap to run and abstain on their own can
        leave this alone; lanes gated on gene membership should override it so
        the router's lane counts reflect what was actually considered.
        """
        return True

    @abstractmethod
    def evaluate(self, variant: Variant) -> list[DrugEvidence]:
        """Return every piece of evidence this lane has, or an empty list.

        An empty list means abstention. It must never mean susceptibility —
        only the coverage-backed assessment in ``modules.synthesis`` may
        conclude that.
        """
        raise NotImplementedError
