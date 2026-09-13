"""Continuous, lineage-stratified accuracy monitoring across cohort runs.

Why this is separate from a one-time validation
-------------------------------------------------
The README's phase-3 gate asks for "per-drug error rates against
pre-registered thresholds on independent isolates — with African lineages
represented, where published tools are weakest." As written, that is a single
measurement from a single run. But lineage generalisability is not something
one cohort can establish: a panel dominated by lineage 4 says nothing about
lineage 1 or 3, and the next cohort ingested might shift that balance
entirely. So this module turns the one-time measurement into a cumulative
one: every ``mycobench run`` that supplies paired phenotypes adds its
isolates to a persistent, per-``(drug, lineage)`` confusion-matrix state
(:class:`MonitoringState`), and :func:`lineage_summary` reports where the
evidence is strong and where it is thin — visible in the HTML report every
run produces (``mycobench/report.py``), which is the closest thing this
offline, no-server tool has to a shared dashboard: every participating site
can run the same command over its own results and get the same rendering.

Idempotence
-----------
Re-running the same cohort must not inflate its own contribution.
:func:`ingest_cohort` refuses (with a note, not silently) to add a
``cohort_id`` it has already recorded — the fix for accidentally doubling a
cohort's isolates just by validating it twice.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from .metrics import AccuracyResult
from .thresholds import MIN_LINEAGES_FOR_GENERALIZABILITY, target_for

_CALLED = ("resistant", "susceptible")
_COUNT_KEYS = ("tp", "fp", "tn", "fn", "abstained")


@dataclass
class MonitoringState:
    """Cumulative per-``(drug, lineage)`` confusion counts across cohorts.

    A plain accumulator, not a hash-chained ledger like
    ``federated.catalogue_update.EvidenceLedger`` — there is nothing here to
    govern, curate or roll back, only counts to add up and a set of cohort
    ids already counted.
    """

    counts: dict[tuple[str, str], dict[str, int]] = field(default_factory=dict)
    ingested_cohort_ids: set[str] = field(default_factory=set)

    def to_json(self) -> dict:
        return {
            "counts": {f"{drug}␟{lineage}": dict(c)
                      for (drug, lineage), c in sorted(self.counts.items())},
            "ingested_cohort_ids": sorted(self.ingested_cohort_ids),
        }

    @classmethod
    def from_json(cls, data: dict) -> "MonitoringState":
        counts = {}
        for key, c in data.get("counts", {}).items():
            drug, _, lineage = key.partition("␟")
            counts[(drug, lineage)] = {k: int(c.get(k, 0)) for k in _COUNT_KEYS}
        return cls(counts=counts,
                  ingested_cohort_ids=set(data.get("ingested_cohort_ids", [])))

    def save(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_json(), indent=2) + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "MonitoringState":
        p = Path(path)
        if not p.is_file():
            return cls()
        return cls.from_json(json.loads(p.read_text(encoding="utf-8")))


def ingest_cohort(
    state: MonitoringState,
    cohort_id: str,
    rows: Iterable[tuple[str, str, str, str]],
) -> tuple[MonitoringState, list[str]]:
    """Fold one cohort's ``(drug, lineage, predicted, phenotype)`` rows in.

    Returns ``(state, notes)``. If ``cohort_id`` was already ingested, this is
    a no-op and the note says so explicitly — silently re-adding the same
    isolates would inflate the apparent evidence base without a single new
    observation.
    """
    if cohort_id in state.ingested_cohort_ids:
        return state, [f"cohort {cohort_id!r} was already ingested; skipped "
                       f"to avoid double-counting its isolates"]

    for drug, lineage, predicted, phenotype in rows:
        phenotype_u = (phenotype or "").strip().upper()
        if phenotype_u not in ("R", "S"):
            continue
        key = (drug, (lineage or "unknown").strip() or "unknown")
        bucket = state.counts.setdefault(
            key, {k: 0 for k in _COUNT_KEYS})
        predicted_l = (predicted or "").strip().lower()
        if predicted_l not in _CALLED:
            bucket["abstained"] += 1
            continue
        if predicted_l == "resistant":
            bucket["tp" if phenotype_u == "R" else "fp"] += 1
        else:
            bucket["fn" if phenotype_u == "R" else "tn"] += 1

    state.ingested_cohort_ids.add(cohort_id)
    return state, []


@dataclass
class LineageAccuracy:
    drug: str
    lineage: str
    accuracy: AccuracyResult


@dataclass
class DrugGeneralizability:
    drug: str
    n_lineages: int
    lineages: list[str]
    generalizable: bool
    detail: str


def lineage_summary(state: MonitoringState) -> dict[str, list[LineageAccuracy]]:
    """Per-drug, per-lineage :class:`~mycobench.metrics.AccuracyResult`."""
    by_drug: dict[str, list[LineageAccuracy]] = {}
    for (drug, lineage), c in state.counts.items():
        n_called = c["tp"] + c["fp"] + c["tn"] + c["fn"]
        result = AccuracyResult(
            drug=drug, n_evaluable=n_called + c["abstained"],
            n_called=n_called, n_abstained=c["abstained"],
            true_positive=c["tp"], false_positive=c["fp"],
            true_negative=c["tn"], false_negative=c["fn"],
            target=target_for(drug))
        by_drug.setdefault(drug, []).append(
            LineageAccuracy(drug=drug, lineage=lineage, accuracy=result))
    for entries in by_drug.values():
        entries.sort(key=lambda la: la.lineage)
    return by_drug


def generalizability(state: MonitoringState) -> dict[str, DrugGeneralizability]:
    """Per-drug: how many independent lineages has this ever been measured in?

    Distinct from ``AccuracyResult.powered`` (which asks "enough isolates
    overall") — a drug can clear the isolate-count bar entirely within one
    lineage and still say nothing about another. Mirrors
    ``federated.catalogue_update``'s own ``min_lineages=2`` bar for the same
    underlying reason: an effect seen in exactly one lineage cannot be told
    apart from a lineage-specific effect.
    """
    lineages_by_drug: dict[str, set[str]] = {}
    for drug, lineage in state.counts:
        if lineage == "unknown":
            continue
        lineages_by_drug.setdefault(drug, set()).add(lineage)

    out: dict[str, DrugGeneralizability] = {}
    for drug, lineages in lineages_by_drug.items():
        n = len(lineages)
        ok = n >= MIN_LINEAGES_FOR_GENERALIZABILITY
        out[drug] = DrugGeneralizability(
            drug=drug, n_lineages=n, lineages=sorted(lineages),
            generalizable=ok,
            detail=(f"measured across {n} lineage(s) "
                    f"(need {MIN_LINEAGES_FOR_GENERALIZABILITY}); "
                    + ("no lineage-specific-effect concern raised by this "
                       "alone" if ok else "not yet established as "
                       "lineage-generalisable")))
    return out
