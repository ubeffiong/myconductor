"""Organism profile loading.

A profile is the versioned bundle that makes the engine organism-agnostic in
practice rather than in intention: reference assembly, the loci that must be
callable per drug, efflux regulators, promoter loci, the drug list and the
regimen definitions. Swapping organisms means shipping another profile, not
editing the engine.

What a profile is *not* is a validation claim. ``ships_validated`` is False for
every profile in this repository, and the renderer says so on every report.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

_HERE = Path(__file__).resolve().parent
DRUG_LOCI_PATH = _HERE / "drug_loci.json"
DRUGS_PATH = _HERE / "drugs.json"


@dataclass
class RegimenDefinition:
    name: str
    drugs: list[str]
    min_required: int
    note: str = ""


@dataclass
class OrganismProfile:
    """Everything organism-specific, loaded from versioned JSON."""

    name: str
    version: str
    reference_assembly: str
    drugs: list[str]
    drug_loci: dict[str, dict[str, list[str]]]
    loci: dict[str, dict]
    efflux_regulators: dict[str, list[str]]
    promoter_loci: dict[str, list[str]]
    regimens: list[RegimenDefinition] = field(default_factory=list)
    ships_validated: bool = False

    # -- lookups ----------------------------------------------------------
    def required_loci(self, drug: str, include_tier2: bool = False) -> list[str]:
        """Loci that must be callable before ``drug`` may be called susceptible.

        Tier 2 loci are excluded by default: requiring them would make almost
        every drug NOT_ASSESSED on targeted assays, which is unhelpfully strict.
        Deployments wanting the stricter gate can opt in.
        """
        entry = self.drug_loci.get(drug)
        if not entry:
            return []
        loci = list(entry.get("tier1", []))
        if include_tier2:
            loci += list(entry.get("tier2", []))
        return loci

    def supports_drug(self, drug: str) -> bool:
        return drug in self.drugs

    def drugs_for_efflux_regulator(self, gene: str) -> list[str]:
        return list(self.efflux_regulators.get(gene, []))

    def drugs_for_promoter(self, gene: str) -> list[str]:
        return list(self.promoter_loci.get(gene, []))

    @property
    def efflux_genes(self) -> set[str]:
        return set(self.efflux_regulators)

    def locus_lengths(self) -> dict[str, int]:
        """Locus spans, for callers that compute callable fractions from BED.

        Empty for every bundled profile — the repository ships no reference
        annotation, and will not invent one.
        """
        return {
            locus: meta["length_bp"]
            for locus, meta in self.loci.items()
            if meta.get("length_bp")
        }

    @property
    def loci_with_unknown_length(self) -> list[str]:
        return sorted(l for l, m in self.loci.items() if not m.get("length_bp"))


def load_profile(drug_loci_path: Optional[Path] = None,
                 drugs_path: Optional[Path] = None) -> OrganismProfile:
    loci_data = json.loads(Path(drug_loci_path or DRUG_LOCI_PATH).read_text(encoding="utf-8"))
    drug_data = json.loads(Path(drugs_path or DRUGS_PATH).read_text(encoding="utf-8"))

    regimens = [
        RegimenDefinition(
            name=r["name"],
            drugs=list(r["drugs"]),
            min_required=int(r["min_effective"]),
            note=r.get("note", ""),
        )
        for r in drug_data.get("regimens", [])
    ]

    return OrganismProfile(
        name=loci_data.get("profile", "unknown"),
        version=loci_data.get("profile_version", "unknown"),
        reference_assembly=loci_data.get("reference_assembly", "unknown"),
        drugs=list(drug_data.get("drugs", [])),
        drug_loci=loci_data.get("drugs", {}),
        loci=loci_data.get("loci", {}),
        efflux_regulators=loci_data.get("efflux_regulators", {}),
        promoter_loci=loci_data.get("promoter_loci", {}),
        regimens=regimens,
        ships_validated=False,
    )
