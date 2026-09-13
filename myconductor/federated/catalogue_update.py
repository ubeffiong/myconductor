"""Federated living catalogue -- the feedback loop.

The design goal: local sites keep patient genomes on-premise (privacy laws) and
share only *aggregated evidence* -- never raw sequence. New confirmed
VUS->phenotype pairs accumulate and, once enough sites agree, promote a variant
into the catalogue. A drift monitor watches predictive performance so a decaying
model is caught rather than trusted.

This module is a local-node scaffold: it shows the shape of the share/aggregate/
promote contract without any network transport.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from ..core.models import AnalysisReport


@dataclass
class VariantObservation:
    """The only thing that leaves a site: a variant key + phenotype tally.

    Deliberately carries NO patient identifiers and NO raw sequence.
    """

    variant_key: str
    drug: str
    resistant_count: int = 0
    susceptible_count: int = 0

    def merge(self, other: "VariantObservation") -> None:
        self.resistant_count += other.resistant_count
        self.susceptible_count += other.susceptible_count

    @property
    def total(self) -> int:
        return self.resistant_count + self.susceptible_count

    @property
    def resistant_fraction(self) -> float:
        return self.resistant_count / self.total if self.total else 0.0


def extract_shareable(report: AnalysisReport,
                      confirmed_phenotypes: dict[str, bool] | None = None
                      ) -> list[VariantObservation]:
    """Turn a local report into privacy-preserving observations to share.

    ``confirmed_phenotypes`` maps variant_key -> was_resistant (e.g. from later
    phenotypic DST). Only confirmed pairs are shared, keeping the signal clean.
    """
    confirmed_phenotypes = confirmed_phenotypes or {}
    obs: list[VariantObservation] = []
    for r in report.drug_results:
        for ev in r.evidence:
            if ev.variant_key not in confirmed_phenotypes:
                continue
            resistant = confirmed_phenotypes[ev.variant_key]
            obs.append(
                VariantObservation(
                    variant_key=ev.variant_key,
                    drug=ev.drug,
                    resistant_count=int(resistant),
                    susceptible_count=int(not resistant),
                )
            )
    return obs


def aggregate(site_batches: list[list[VariantObservation]]) -> dict[str, VariantObservation]:
    """Federated aggregation across sites (last step a coordinator runs)."""
    merged: dict[str, VariantObservation] = {}
    for batch in site_batches:
        for o in batch:
            key = f"{o.variant_key}|{o.drug}"
            if key in merged:
                merged[key].merge(o)
            else:
                merged[key] = VariantObservation(
                    o.variant_key, o.drug, o.resistant_count, o.susceptible_count
                )
    return merged


def promote_candidates(aggregated: dict[str, VariantObservation],
                       min_total: int = 20,
                       min_fraction: float = 0.75) -> list[dict]:
    """Variants that now clear the evidence bar to enter the catalogue."""
    promoted = []
    for o in aggregated.values():
        if o.total >= min_total and o.resistant_fraction >= min_fraction:
            promoted.append(
                {
                    "gene_change": o.variant_key,
                    "drug": o.drug,
                    "call": "resistant",
                    "evidence_n": o.total,
                    "resistant_fraction": round(o.resistant_fraction, 3),
                    "who_grade": "Assoc w R - provisional (federated)",
                }
            )
    return promoted


def drift_alarm(recent_auroc: list[float], baseline: float = 0.9,
                tolerance: float = 0.1) -> bool:
    """Flag model drift if recent AUROC falls a tolerance below baseline."""
    if not recent_auroc:
        return False
    return (baseline - (sum(recent_auroc) / len(recent_auroc))) >= tolerance
