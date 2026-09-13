"""Versioned regulatory-region imports. Region membership is not a resistance grade."""
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional
from ..core.models import Region


@dataclass(frozen=True)
class RegulatoryRegion:
    id: str
    name: str
    organism: str
    reference_assembly: str
    type: str
    target_genes: list[str]
    drug_associations: list[str]
    evidence_tier: str
    source: str
    coordinates: Optional[dict] = None
    variant_labels: tuple[str, ...] = ()
    status: str = "active"

    def __post_init__(self):
        for name in ("id", "name", "organism", "reference_assembly", "source"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"regulatory region requires {name}")
        if self.type not in {"PROMOTER", "INTERGENIC", "OPERATOR", "TERMINATOR"}:
            raise ValueError("unknown regulatory region type")
        if self.evidence_tier not in {"CATALOGUED", "CANDIDATE", "UNDER_INVESTIGATION"}:
            raise ValueError("unknown regulatory evidence tier")
        if self.status not in {"active", "retracted"}:
            raise ValueError("unknown regulatory entry status")
        for name in ("target_genes", "drug_associations"):
            if not isinstance(getattr(self, name), (list, tuple)) or not getattr(self, name) or not all(
                    isinstance(x, str) and x for x in getattr(self, name)):
                raise ValueError(f"regulatory region requires {name}")
        if self.coordinates is not None:
            c = self.coordinates
            if not c.get("chrom") or c.get("strand") not in {"+", "-"} or not all(
                    isinstance(c.get(k), int) and not isinstance(c[k], bool) for k in ("start", "end")):
                raise ValueError("regulatory coordinates require chrom, integer start/end, and strand")
            if not 1 <= c["start"] <= c["end"]:
                raise ValueError("regulatory coordinates use 1-based closed intervals")


class RegulatoryCatalogue:
    def __init__(self, entries=(), version="empty-1", illustrative=True):
        self.entries = tuple(entries)
        self.version, self.illustrative = version, illustrative
        if not isinstance(version, str) or not version or not isinstance(illustrative, bool):
            raise ValueError("regulatory catalogue requires version and boolean illustrative")
        if len({e.id for e in self.entries}) != len(self.entries):
            raise ValueError("duplicate regulatory region ID")

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        if d.get("schema") != "myconductor.regulatory.v1":
            raise ValueError("unknown regulatory schema")
        return cls([RegulatoryRegion(**e) for e in d["regions"]], d["version"], d["illustrative"])

    def to_dict(self):
        return {"schema": "myconductor.regulatory.v1", "version": self.version,
                "illustrative": self.illustrative, "regions": [asdict(e) for e in self.entries]}

    def matches(self, variant, organism):
        matches = []
        for entry in self.entries:
            if entry.status != "active" or entry.organism != organism or entry.reference_assembly != variant.identity.assembly:
                continue
            if entry.variant_labels:
                matched = variant.label() in entry.variant_labels
            elif entry.coordinates:
                c, identity = entry.coordinates, variant.identity
                matched = identity.chrom == c["chrom"] and identity.pos is not None and c["start"] <= identity.pos <= c["end"]
            else:
                matched = variant.gene in entry.target_genes and variant.region in {Region.PROMOTER, Region.INTERGENIC}
            if matched:
                matches.append(entry)
        return matches

    def evidence_tier_for(self, gene, change):
        # Return every attributed match; never pick a stronger tier by vote.
        return [(e.id, e.evidence_tier) for e in self.entries if e.status == "active"
                and gene in e.target_genes and f"{gene}_{change}" in e.variant_labels]
