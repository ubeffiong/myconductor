"""Attributed co-observation notes. No interaction prediction or confidence boost."""
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ..core.models import Call, Tier


@dataclass(frozen=True)
class EpistasisRule:
    rule_id: str
    drug: str
    primary_variant_key_or_gene: str
    partner_variant_key_or_gene: str
    interaction: str
    note: str
    source: str
    organism: str
    reference_assembly: str
    reviewed_by: str
    reviewed_at: str
    confidence_adjustment: float = 0
    status: str = "active"

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name != "confidence_adjustment" and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"interaction rule requires {name}")
        if self.interaction not in {"compensatory", "antagonistic"}:
            raise ValueError("unsupported interaction annotation")
        if self.status not in {"active", "retracted"}:
            raise ValueError("unknown interaction rule status")
        if self.confidence_adjustment != 0:
            raise ValueError("uncalibrated interaction notes cannot adjust resistance confidence")


class EpistasisTable:
    def __init__(self, rules=(), version="empty-1", illustrative=False):
        self.rules = tuple(rules)
        self.version = version
        self.illustrative = illustrative
        if not isinstance(version, str) or not version.strip():
            raise ValueError("interaction table requires version")
        if not isinstance(illustrative, bool):
            raise ValueError("illustrative must be boolean")
        if len({r.rule_id for r in self.rules}) != len(self.rules):
            raise ValueError("duplicate interaction rule ID")

    @classmethod
    def load(cls, path):
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data.get("schema") != "myconductor.epistasis.v1":
            raise ValueError("unknown interaction-table schema")
        return cls([EpistasisRule(**r) for r in data["rules"]], data["version"], data["illustrative"])

    def to_dict(self):
        return {"schema": "myconductor.epistasis.v1", "version": self.version,
                "illustrative": self.illustrative, "rules": [asdict(r) for r in self.rules]}

    def annotate(self, results, variants, organism, reference_assembly):
        notes = []
        for rule in self.rules:
            if rule.status != "active" or rule.organism != organism or rule.reference_assembly != reference_assembly:
                continue
            result = next((r for r in results if r.drug == rule.drug), None)
            if result is None or result.call is not Call.RESISTANT or result.tier is not Tier.CATALOGUED:
                continue
            def matches(selector):
                return [v for v in variants if v.passed_filters and v.identity.assembly == reference_assembly
                        and selector in {v.key(), v.gene}]
            primaries, partners = matches(rule.primary_variant_key_or_gene), matches(rule.partner_variant_key_or_gene)
            pairs = sorted({(a.key(), b.key()) for a in primaries for b in partners if a.key() != b.key()})
            if pairs:
                notes.append(dict(asdict(rule), matched_pairs=[list(p) for p in pairs],
                                  table_version=self.version, illustrative=self.illustrative,
                                  effect="annotation_only", tier=Tier.INFERRED.value,
                                  interpretation="Co-observed variants match a supplied rule. Interaction, cellular linkage and fitness effects are not established here."))
        return notes
