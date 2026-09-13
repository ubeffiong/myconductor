"""Known-variant module: the wrapped, catalogue-based lookup.

This is the lane that stands in for TB-Profiler / Mykrobe / MAST -- a static,
graded catalogue of known resistance mutations. In a real deployment you would
back this with the full WHO catalogue (or delegate to one of those tools and
parse its output); here we load a small illustrative JSON.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ..core.models import Call, DrugEvidence, Route, Variant
from .base import VariantModule

_CATALOGUE_PATH = Path(__file__).resolve().parent.parent / "catalogue" / "mtb_amr_catalogue.json"


class CatalogueModule(VariantModule):
    name = "catalogue"

    def __init__(self, path: Path = _CATALOGUE_PATH):
        data = json.loads(Path(path).read_text())
        self.version: str = data["catalogue_version"]
        self.efflux_genes: set[str] = set(data.get("efflux_genes", []))
        self.regulatory_regions: set[str] = set(data.get("regulatory_regions", []))
        self._by_key: dict[str, dict] = {
            f"{v['gene']}_{v['change']}": v for v in data["variants"]
        }

    def knows(self, variant: Variant) -> bool:
        return variant.key() in self._by_key

    def evaluate(self, variant: Variant) -> Optional[DrugEvidence]:
        entry = self._by_key.get(variant.key())
        if entry is None:
            return None
        return DrugEvidence(
            drug=entry["drug"],
            call=Call(entry["call"]),
            confidence=float(entry["confidence"]),
            route=Route.CATALOGUE,
            variant_key=variant.key(),
            who_grade=entry.get("who_grade"),
            rationale=f"Catalogued mutation ({entry.get('who_grade', 'graded')}).",
        )
