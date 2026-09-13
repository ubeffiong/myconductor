"""Catalogue lane — graded, known resistance mutations.

This is the only lane (alongside laboratory phenotype) permitted to establish
``RESISTANT``, because it is the only one backed by graded evidence.

The bundled JSON is a small illustrative subset. It is not the WHO catalogue,
and the module says so in its provenance rather than leaving the reader to
assume otherwise: ``EngineRef.database_version`` carries the demo tag through
to the report and the FHIR bundle. ``adapters/who_catalogue.py`` is the
ingester for the real thing.

Lookup is by ``VariantIdentity``: coordinate key first, falling back to the
gene/label pair for catalogue entries and inputs that carry no coordinates.
Matching on labels alone is recorded as a limitation on the evidence, because
annotator spellings differ and a label match is weaker than a coordinate match.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ..core.models import (
    Call,
    DrugEvidence,
    EngineRef,
    Lane,
    Tier,
    Variant,
)
from .base import VariantModule

_CATALOGUE_PATH = (Path(__file__).resolve().parent.parent
                   / "catalogue" / "mtb_amr_catalogue.json")

#: Bundled illustrative catalogues, by organism name. See
#: ``catalogue/profile.py::BUNDLED_ORGANISMS`` for the matching drug/loci
#: profile of each; the two registries are kept in step so
#: ``organism="mabscessus"`` picks a consistent pair everywhere.
BUNDLED_CATALOGUES: dict[str, Path] = {
    "mtbc": _CATALOGUE_PATH,
    "mabscessus": (Path(__file__).resolve().parent.parent / "catalogue"
                  / "organisms" / "mabscessus" / "catalogue.json"),
}

_LABEL_MATCH_LIMITATION = (
    "matched on gene/label, not coordinates; annotator spellings differ, so "
    "this match is weaker than a coordinate match"
)


class CatalogueModule(VariantModule):
    name = "catalogue"
    lane = Lane.CATALOGUE

    def __init__(self, path: Optional[Path] = None,
                organism: Optional[str] = None):
        if organism is not None and path is None:
            if organism not in BUNDLED_CATALOGUES:
                raise ValueError(
                    f"unknown organism {organism!r}; bundled catalogues are "
                    f"{sorted(BUNDLED_CATALOGUES)}. Ship your own by "
                    f"passing path= directly."
                )
            path = BUNDLED_CATALOGUES[organism]
        path = path or _CATALOGUE_PATH
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        self.version: str = data["catalogue_version"]
        self.is_illustrative: bool = bool(data.get("illustrative", True))
        self.efflux_genes: set[str] = set(data.get("efflux_genes", []))

        self.engine = EngineRef(
            name="myconductor-bundled-catalogue",
            version="0.2.0",
            database="mtb_amr_catalogue.json",
            database_version=self.version,
        )

        self._by_coord: dict[str, dict] = {}
        self._by_label: dict[str, dict] = {}
        for entry in data["variants"]:
            label = f"{entry['gene']}_{entry['change']}"
            self._by_label[label] = entry
            coord = entry.get("coordinate_key")
            if coord:
                self._by_coord[coord] = entry

    # -- lookup -----------------------------------------------------------
    def _lookup(self, variant: Variant) -> tuple[Optional[dict], bool]:
        """Return ``(entry, matched_on_coordinates)``."""
        key = variant.identity.key()
        if key in self._by_coord:
            return self._by_coord[key], True
        return self._by_label.get(variant.label()), False

    @property
    def labels(self) -> set[str]:
        """Every gene/label the catalogue knows.

        Used to stop the VUS workbench re-examining variants the catalogue has
        already graded.
        """
        return set(self._by_label)

    def knows(self, variant: Variant) -> bool:
        entry, _ = self._lookup(variant)
        return entry is not None

    def applies_to(self, variant: Variant) -> bool:
        return self.knows(variant)

    # -- evaluation -------------------------------------------------------
    def evaluate(self, variant: Variant) -> list[DrugEvidence]:
        entry, by_coord = self._lookup(variant)
        if entry is None:
            return []

        call = _CALL_MAP.get(entry["call"])
        if call is None:
            raise ValueError(
                f"catalogue entry {variant.label()!r} has unrecognised call "
                f"{entry['call']!r}"
            )

        limitations: list[str] = []
        if not by_coord:
            limitations.append(_LABEL_MATCH_LIMITATION)
        if self.is_illustrative:
            limitations.append(
                "from the bundled illustrative catalogue, not the WHO catalogue"
            )

        # Catalogue entries may name several drugs.
        drugs = entry.get("drugs") or [entry["drug"]]
        grade = entry.get("who_grade")
        return [
            DrugEvidence(
                drug=drug,
                call=call,
                tier=Tier.CATALOGUED,
                lane=Lane.CATALOGUE,
                confidence=float(entry["confidence"]),
                variant=variant.identity,
                who_grade=grade,
                engine=self.engine,
                limitations=tuple(limitations),
                rationale=(
                    f"Catalogued {call.value} for {drug} "
                    f"(grade: {grade or 'ungraded'})."
                ),
                # Traceability for multi-site governance (see
                # federated/catalogue_governance.py): exactly which catalogue
                # version and which entry produced this call, so two sites
                # disagreeing because they run different catalogue versions
                # can be told apart from two sites disagreeing despite
                # running the same one.
                catalogue_version=self.version,
                rule_id=f"{self.version}:{variant.label()}:{drug}:{call.value}",
            )
            for drug in drugs
        ]


_CALL_MAP = {
    "resistant": Call.RESISTANT,
    "susceptible": Call.SUSCEPTIBLE,
    "indeterminate": Call.INDETERMINATE,
    # "Not assoc w R" entries are evidence that a variant does not confer
    # resistance. They do not by themselves make the drug usable -- coverage
    # still has to be shown -- but they are a susceptible-leaning catalogue
    # statement about this variant.
    "not_associated": Call.SUSCEPTIBLE,
}
