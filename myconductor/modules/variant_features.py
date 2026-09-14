"""Imported feature vectors: what a model was shown, recorded so it can be audited.

Myconductor does not compute features. Conservation scores, structural distances
and stability changes come from tools with their own licences, databases and
failure modes, and recomputing them here would add a second unvalidated
component to explain. What this module does is record the feature vector a model
was given, alongside the prediction, so a reviewer can ask the question that
actually matters about a surprising prediction: *what did it see?*

A missing feature is missing, never zero
---------------------------------------
This is the whole point. A model handed ``delta_delta_g = 0.0`` for a variant
with no structure cannot distinguish "this substitution is energetically
neutral" from "nobody folded this protein". Imputing a default at the feature
level reintroduces, one layer down, exactly the confusion between absence of
evidence and evidence of absence that the call states exist to prevent.

So availability is not a parallel flag that can drift out of step with the
values — that failure mode is real and this codebase has already been bitten by
it once. A feature is in ``values`` or it is in ``unavailable`` with a stated
reason, never both and never neither. The two dictionaries are validated against
each other on construction.

Nothing here ranks or calls
---------------------------
Features attach to a VUS as a ``Dimension`` after ranking, exactly as imported
structural annotations do. They carry no tier, establish nothing, and are
excluded from scores. A feature vector is context for a human reading a
prediction, not an input to any decision this system makes.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

from ..core.models import Dimension

#: Features this schema knows how to describe. An unknown name is refused
#: rather than silently carried, so a typo cannot become a phantom feature that
#: a reviewer reads as meaningful.
KNOWN_FEATURES = {
    # Conservation
    "sift_score", "polyphen_score", "phylop_score", "gerp_score",
    # Structural — require a structure to have been solved or predicted
    "residue_plddt", "binding_pocket_distance", "delta_delta_g", "sasa_change",
    # Population
    "global_prevalence", "lineage_prevalence",
    # Regulatory
    "distance_to_tss",
}

#: Reasons a feature may legitimately be unavailable. Free text is refused: a
#: reason a reader cannot act on is not a reason.
UNAVAILABLE_REASONS = {
    "no_structure",            # no solved or predicted structure for this protein
    "low_confidence_region",   # structure exists but the residue is unreliable
    "non_coding",              # feature is meaningless for a non-coding variant
    "not_in_reference",        # coordinate absent from the reference used
    "tool_failed",             # the computing tool errored on this variant
    "not_computed",            # simply never run
    "licence_restricted",      # the computing tool is not available to this site
}


@dataclass(frozen=True)
class VariantFeatures:
    """One variant's feature vector, with absences stated rather than filled."""

    variant_key: str
    organism: str
    reference_assembly: str
    source: str
    source_version: str
    timestamp: str
    #: Features that were computed. Names must be in ``KNOWN_FEATURES``.
    values: dict[str, float] = field(default_factory=dict)
    #: Features that were not, mapped to a reason from ``UNAVAILABLE_REASONS``.
    unavailable: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ("variant_key", "organism", "reference_assembly",
                     "source", "source_version", "timestamp"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"variant features require {name}")

        for name, value in self.values.items():
            if name not in KNOWN_FEATURES:
                raise ValueError(f"unknown feature {name!r}")
            if not isinstance(value, (int, float)) or isinstance(value, bool) \
                    or not math.isfinite(value):
                raise ValueError(f"feature {name} must be a finite number")

        for name, reason in self.unavailable.items():
            if name not in KNOWN_FEATURES:
                raise ValueError(f"unknown feature {name!r}")
            if reason not in UNAVAILABLE_REASONS:
                raise ValueError(
                    f"feature {name} needs a reason from "
                    f"{sorted(UNAVAILABLE_REASONS)}, got {reason!r}")

        # The invariant that keeps availability from drifting: a feature is
        # present or absent, never both. Two structures that each claim to know
        # a feature's status will eventually disagree, and the disagreement is
        # silent — so it is refused at construction instead.
        both = set(self.values) & set(self.unavailable)
        if both:
            raise ValueError(
                f"feature(s) {sorted(both)} are both valued and unavailable")

    def available(self, name: str) -> bool:
        """Derived, never stored: there is exactly one source of truth."""
        if name not in KNOWN_FEATURES:
            raise ValueError(f"unknown feature {name!r}")
        return name in self.values

    def get(self, name: str) -> Optional[float]:
        """The value, or ``None``. ``None`` means unknown, never zero."""
        if name not in KNOWN_FEATURES:
            raise ValueError(f"unknown feature {name!r}")
        return self.values.get(name)

    def reason(self, name: str) -> Optional[str]:
        return self.unavailable.get(name)

    @property
    def described(self) -> set[str]:
        return set(self.values) | set(self.unavailable)

    @property
    def undeclared(self) -> set[str]:
        """Features neither valued nor explained. Reported, not assumed."""
        return KNOWN_FEATURES - self.described

    def describe(self) -> str:
        parts = [f"{len(self.values)} feature(s) from {self.source} "
                 f"({self.source_version})"]
        if self.unavailable:
            grouped: dict[str, list[str]] = {}
            for name, reason in sorted(self.unavailable.items()):
                grouped.setdefault(reason, []).append(name)
            parts.append("; ".join(f"{reason}: {', '.join(names)}"
                                   for reason, names in sorted(grouped.items())))
        if self.undeclared:
            parts.append(f"{len(self.undeclared)} feature(s) not declared "
                         f"either way")
        return " | ".join(parts)


class VariantFeatureTable:
    """Feature vectors keyed by variant, organism and source version."""

    def __init__(self, records=()):
        self.records = tuple(records)
        keys = [(r.variant_key, r.organism, r.reference_assembly,
                 r.source, r.source_version) for r in self.records]
        if len(keys) != len(set(keys)):
            raise ValueError(
                "duplicate feature vector; use separate source versions")

    @classmethod
    def load(cls, path) -> "VariantFeatureTable":
        rows = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(rows, list):
            raise ValueError("variant features must be a JSON list")
        return cls(VariantFeatures(**row) for row in rows)

    def to_dict(self) -> list[dict]:
        return [asdict(r) for r in self.records]

    def for_variant(self, variant_key: str, organism: str,
                    reference_assembly: str) -> list[VariantFeatures]:
        return [r for r in self.records
                if r.variant_key == variant_key and r.organism == organism
                and r.reference_assembly == reference_assembly]

    def attach(self, priorities, organism: str,
               reference_assembly: str) -> list[dict]:
        """Attach feature vectors as dimensions. Ranking is untouched.

        Distinct sources are preserved rather than merged: two tools that
        disagree about a variant's conservation are reporting a disagreement,
        and averaging it away destroys the only signal that matters.
        """
        attached: list[dict] = []
        for priority in priorities:
            records = self.for_variant(priority.variant_key, organism,
                                       reference_assembly)
            if not records:
                continue
            payload = [{
                "source": r.source, "source_version": r.source_version,
                "values": dict(r.values), "unavailable": dict(r.unavailable),
                "undeclared": sorted(r.undeclared),
            } for r in records]
            dimension = Dimension(
                "variant_features", True, payload,
                source="; ".join(r.source for r in records),
                note=("Imported feature vectors as supplied to a model; "
                      "missing features are reported unavailable, never "
                      "imputed, and none of this enters the ranking."))
            priority.dimensions = [d for d in priority.dimensions
                                   if d.name != "variant_features"]
            priority.dimensions.append(dimension)
            attached.extend(payload)
        return attached
