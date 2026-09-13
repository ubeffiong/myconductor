"""Module interfaces.

Every specialist plugs into the router through one of these. Swapping the
demo VUS scorer for a trained XGBoost model, or the mini catalogue for the
full WHO catalogue, means implementing the same tiny interface -- nothing
else in the system changes.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional

from ..core.models import DrugEvidence, Variant


class VariantModule(ABC):
    """A module that turns a single variant into at most one DrugEvidence."""

    name: str = "unnamed"

    @abstractmethod
    def evaluate(self, variant: Variant) -> Optional[DrugEvidence]:
        """Return evidence for this variant, or None if this module abstains."""
        raise NotImplementedError
