"""Governed multi-site catalogue learning.

The statistical defect this rewrite fixes
-----------------------------------------
The previous module counted **variants**. It promoted a variant into the
catalogue once 20 observations reached 75% resistance. That denominator is
wrong, and the consequence is not subtle.

Resistance is a property of an *isolate–drug observation*, not of a variant. A
resistant isolate typically carries several variants: the causal one, lineage
markers it inherited, and hitchhikers in linkage with it. Counting per variant
credits every one of them with the phenotype. Promote on that basis and you
enshrine lineage markers as resistance determinants, then distribute them to
other sites in a signed catalogue — where they will cause susceptible isolates
to be reported resistant, and patients to lose usable drugs.

So this module works from ``IsolateObservation`` records, counts each isolate
once per drug, and refuses promotion unless the association survives three
confounding checks:

* **lineage stratification** — the association must hold within at least two
  lineages, so a marker of a resistant-enriched lineage cannot pass;
* **co-occurrence** — a variant that almost always appears alongside a known
  resistance determinant cannot be credited with the phenotype;
* **site diversity** — evidence from one laboratory is one laboratory's
  systematic error.

Promotion is still never automatic. A candidate that clears all three enters a
review queue for expert curation, and every state change is written to an
append-only hash-chained ledger so a bad release can be traced and reversed.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Iterable, Optional

# -- the unit of evidence -------------------------------------------------
@dataclass(frozen=True)
class IsolateObservation:
    """One isolate, one drug, one phenotype — the correct denominator.

    Carries no patient identifier and no sequence. It does carry a site and a
    lineage, which is why ``transport.k_anonymity_violations`` exists: a rare
    variant plus a site plus a phenotype can re-identify a patient even with no
    name attached.
    """

    isolate_id: str
    site_id: str
    drug: str
    phenotype: Optional[bool]          # True = resistant, False = susceptible
    variant_keys: tuple[str, ...] = ()
    lineage: Optional[str] = None
    dst_method: Optional[str] = None
    critical_concentration: Optional[float] = None
    mic: Optional[float] = None
    mic_unit: str = "mg/L"
    lab_quality: str = "unknown"       # accredited | participating | unknown

    @property
    def usable(self) -> bool:
        """Only a completed DST with a known phenotype contributes."""
        return self.phenotype is not None and bool(self.variant_keys)

    @property
    def weight(self) -> float:
        """Laboratory-quality weighting, applied to reporting only.

        Deliberately not applied to the promotion counts: down-weighting a
        laboratory's results is a governance decision, not a silent
        statistical one.
        """
        return {"accredited": 1.0, "participating": 0.75}.get(self.lab_quality, 0.5)


@dataclass
class VariantDrugAssociation:
    """Isolate-level evidence for one variant–drug pair."""

    variant_key: str
    drug: str
    resistant_isolates: set[str] = field(default_factory=set)
    susceptible_isolates: set[str] = field(default_factory=set)
    sites: set[str] = field(default_factory=set)
    by_lineage: dict[str, list[int]] = field(default_factory=dict)  # [R, S]
    cooccurring: dict[str, int] = field(default_factory=dict)
    mics: list[float] = field(default_factory=list)
    dst_methods: set[str] = field(default_factory=set)
    lab_qualities: set[str] = field(default_factory=set)
    #: The drug's whole tested cohort, isolate-level -- including isolates that
    #: do NOT carry this variant. Without it there is no comparison group, and
    #: a variant carried by every tested isolate looks perfectly associated
    #: with resistance while saying nothing at all.
    cohort_resistant: int = 0
    cohort_susceptible: int = 0

    @property
    def n_resistant(self) -> int:
        return len(self.resistant_isolates)

    @property
    def n_susceptible(self) -> int:
        return len(self.susceptible_isolates)

    @property
    def n_isolates(self) -> int:
        return self.n_resistant + self.n_susceptible

    @property
    def n_sites(self) -> int:
        return len(self.sites)

    @property
    def resistant_fraction(self) -> float:
        return self.n_resistant / self.n_isolates if self.n_isolates else 0.0

    def supporting_lineages(self, min_isolates: int,
                            min_fraction: float) -> list[str]:
        """Lineages in which the association independently holds."""
        out = []
        for lineage, (r, s) in self.by_lineage.items():
            total = r + s
            if total >= min_isolates and total and (r / total) >= min_fraction:
                out.append(lineage)
        return sorted(out)

    def top_cooccurring(self, n: int = 5) -> list[tuple[str, int]]:
        return sorted(self.cooccurring.items(), key=lambda kv: -kv[1])[:n]

    # -- the 2x2 the association is actually judged on --------------------
    @property
    def cohort_size(self) -> int:
        return self.cohort_resistant + self.cohort_susceptible

    @property
    def absent_resistant(self) -> int:
        """Resistant isolates tested for this drug that lack the variant."""
        return max(0, self.cohort_resistant - self.n_resistant)

    @property
    def absent_susceptible(self) -> int:
        """Susceptible isolates tested for this drug that lack the variant."""
        return max(0, self.cohort_susceptible - self.n_susceptible)

    @property
    def has_variant_absent_controls(self) -> bool:
        """Is there any comparison group without this variant at all?"""
        return (self.absent_resistant + self.absent_susceptible) > 0


def aggregate(batches: Iterable[Iterable[IsolateObservation]]
              ) -> dict[tuple[str, str], VariantDrugAssociation]:
    """Aggregate across sites with isolate-level denominators.

    Two passes, because a variant's association can only be judged against the
    isolates that *lack* it. The first builds each drug's tested cohort; the
    second builds the per-variant table within it.

    Requiring susceptible isolates that *carry* the variant would be the wrong
    control and would penalise precisely the strongest determinants, which are
    rarely seen in susceptible isolates. The right control is the
    variant-absent group.

    An isolate contributes at most once to any variant–drug cell, however many
    times it appears in the input.
    """
    batches = [list(batch) for batch in batches]

    # Pass 1: each drug's tested cohort, isolate-level.
    cohort: dict[str, dict[str, bool]] = {}
    for batch in batches:
        for obs in batch:
            if not obs.usable:
                continue
            cohort.setdefault(obs.drug, {}).setdefault(
                obs.isolate_id, bool(obs.phenotype))

    # Pass 2: the per variant–drug tables.
    merged: dict[tuple[str, str], VariantDrugAssociation] = {}
    for batch in batches:
        for obs in batch:
            if not obs.usable:
                continue
            for variant_key in set(obs.variant_keys):
                key = (variant_key, obs.drug)
                assoc = merged.get(key)
                if assoc is None:
                    assoc = VariantDrugAssociation(variant_key, obs.drug)
                    merged[key] = assoc

                # Isolate counted once, in exactly one bucket.
                if obs.isolate_id in assoc.resistant_isolates or \
                        obs.isolate_id in assoc.susceptible_isolates:
                    continue

                if obs.phenotype:
                    assoc.resistant_isolates.add(obs.isolate_id)
                else:
                    assoc.susceptible_isolates.add(obs.isolate_id)

                assoc.sites.add(obs.site_id)
                lineage = obs.lineage or "unknown"
                bucket = assoc.by_lineage.setdefault(lineage, [0, 0])
                bucket[0 if obs.phenotype else 1] += 1

                for other in set(obs.variant_keys):
                    if other != variant_key:
                        assoc.cooccurring[other] = \
                            assoc.cooccurring.get(other, 0) + 1
                if obs.mic is not None:
                    assoc.mics.append(obs.mic)
                if obs.dst_method:
                    assoc.dst_methods.add(obs.dst_method)
                assoc.lab_qualities.add(obs.lab_quality)

    # Attach the cohort each association is judged against.
    for (_variant_key, drug), assoc in merged.items():
        phenotypes = cohort.get(drug, {})
        assoc.cohort_resistant = sum(1 for r in phenotypes.values() if r)
        assoc.cohort_susceptible = sum(1 for r in phenotypes.values() if not r)
    return merged


# -- promotion ------------------------------------------------------------
@dataclass
class PromotionThresholds:
    min_isolates: int = 50
    min_resistant_fraction: float = 0.80
    min_sites: int = 3
    #: Susceptible isolates tested for the drug anywhere in the cohort -- NOT
    #: susceptible isolates carrying the variant. See ``aggregate``.
    min_cohort_susceptible: int = 10
    #: There must be some isolate tested for this drug that lacks the variant.
    require_variant_absent_controls: bool = True
    min_lineages: int = 2
    min_isolates_per_lineage: int = 10
    max_cooccurrence_fraction: float = 0.80
    require_mic: bool = False


@dataclass
class PromotionCandidate:
    variant_key: str
    drug: str
    n_isolates: int
    n_resistant: int
    n_susceptible: int
    n_sites: int
    resistant_fraction: float
    supporting_lineages: list[str]
    blocked_reasons: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    mic_range: Optional[tuple[float, float]] = None
    cohort_susceptible: int = 0
    variant_absent_isolates: int = 0

    @property
    def eligible(self) -> bool:
        return not self.blocked_reasons

    def to_catalogue_entry(self) -> dict:
        if not self.eligible:
            raise ValueError(
                f"{self.variant_key}/{self.drug} is not eligible for promotion: "
                + "; ".join(self.blocked_reasons)
            )
        gene, _, change = self.variant_key.partition("_")
        return {
            "gene": gene or self.variant_key,
            "change": change or self.variant_key,
            "drugs": [self.drug],
            "call": "resistant",
            "who_grade": "provisional (federated, pending expert curation)",
            "confidence": round(self.resistant_fraction, 3),
            "evidence_n_isolates": self.n_isolates,
            "evidence_n_sites": self.n_sites,
            "supporting_lineages": self.supporting_lineages,
        }


def evaluate_candidate(
    assoc: VariantDrugAssociation,
    known_resistance_variants: Optional[set[str]] = None,
    thresholds: Optional[PromotionThresholds] = None,
) -> PromotionCandidate:
    """Assess one association, recording every reason it falls short."""
    t = thresholds or PromotionThresholds()
    known = known_resistance_variants or set()

    lineages = assoc.supporting_lineages(t.min_isolates_per_lineage,
                                         t.min_resistant_fraction)
    candidate = PromotionCandidate(
        variant_key=assoc.variant_key,
        drug=assoc.drug,
        n_isolates=assoc.n_isolates,
        n_resistant=assoc.n_resistant,
        n_susceptible=assoc.n_susceptible,
        n_sites=assoc.n_sites,
        resistant_fraction=round(assoc.resistant_fraction, 3),
        supporting_lineages=lineages,
        mic_range=((min(assoc.mics), max(assoc.mics)) if assoc.mics else None),
        cohort_susceptible=assoc.cohort_susceptible,
        variant_absent_isolates=assoc.absent_resistant + assoc.absent_susceptible,
    )

    if assoc.n_isolates < t.min_isolates:
        candidate.blocked_reasons.append(
            f"{assoc.n_isolates} isolates (need {t.min_isolates})")
    if assoc.resistant_fraction < t.min_resistant_fraction:
        candidate.blocked_reasons.append(
            f"resistant fraction {assoc.resistant_fraction:.2f} "
            f"(need {t.min_resistant_fraction:.2f})")
    if assoc.n_sites < t.min_sites:
        candidate.blocked_reasons.append(
            f"{assoc.n_sites} site(s) (need {t.min_sites}); single-site "
            f"evidence cannot be separated from that laboratory's systematic "
            f"error")
    if assoc.cohort_susceptible < t.min_cohort_susceptible:
        candidate.blocked_reasons.append(
            f"only {assoc.cohort_susceptible} susceptible control(s) tested "
            f"for {assoc.drug} in the whole cohort (need "
            f"{t.min_cohort_susceptible}); without susceptible controls the "
            f"specificity of the association is unknown")
    if t.require_variant_absent_controls and not assoc.has_variant_absent_controls:
        candidate.blocked_reasons.append(
            f"every isolate tested for {assoc.drug} carries this variant, so "
            f"there is no variant-absent comparison group; the association "
            f"cannot be separated from the cohort's baseline resistance")

    # Lineage confounding.
    if len(lineages) < t.min_lineages:
        observed = len([l for l, (r, s) in assoc.by_lineage.items() if r + s])
        candidate.blocked_reasons.append(
            f"association holds in {len(lineages)} lineage(s) of {observed} "
            f"observed (need {t.min_lineages}); a variant that tracks one "
            f"lineage may be a lineage marker rather than a determinant")

    # Co-occurrence confounding.
    for other, count in assoc.top_cooccurring():
        if other not in known:
            continue
        fraction = count / assoc.n_resistant if assoc.n_resistant else 0.0
        if fraction >= t.max_cooccurrence_fraction:
            candidate.blocked_reasons.append(
                f"co-occurs with the known determinant {other} in "
                f"{fraction:.0%} of resistant isolates; causality cannot be "
                f"attributed to this variant")

    if t.require_mic and not assoc.mics:
        candidate.blocked_reasons.append(
            "no MIC measurements; binary DST alone was required to carry MICs")

    if len(assoc.dst_methods) > 1:
        candidate.notes.append(
            f"mixed DST methods ({', '.join(sorted(assoc.dst_methods))}); "
            f"critical concentrations are not interchangeable")
    if assoc.lab_qualities and assoc.lab_qualities <= {"unknown"}:
        candidate.notes.append(
            "no laboratory accreditation status reported for any contributing "
            "site")
    if candidate.eligible:
        candidate.notes.append(
            "clears the statistical bar; requires expert curation before "
            "release — promotion is never automatic")
    return candidate


def evaluate_all(
    aggregated: dict[tuple[str, str], VariantDrugAssociation],
    known_resistance_variants: Optional[set[str]] = None,
    thresholds: Optional[PromotionThresholds] = None,
) -> list[PromotionCandidate]:
    return [
        evaluate_candidate(assoc, known_resistance_variants, thresholds)
        for _, assoc in sorted(aggregated.items())
    ]


# -- governance -----------------------------------------------------------
class ReviewState(str, Enum):
    PENDING = "pending"
    APPROVED = "approved"
    REJECTED = "rejected"
    ROLLED_BACK = "rolled_back"


@dataclass
class LedgerEntry:
    seq: int
    timestamp_utc: str
    action: str
    payload: dict
    prev_hash: str
    entry_hash: str


class EvidenceLedger:
    """Append-only, hash-chained record of every catalogue decision.

    Not a blockchain and not a substitute for access control — a tamper-evident
    audit trail, so a bad release can be found and reversed rather than
    silently overwritten.
    """

    GENESIS = "0" * 64

    def __init__(self) -> None:
        self.entries: list[LedgerEntry] = []

    def append(self, action: str, payload: dict) -> LedgerEntry:
        prev = self.entries[-1].entry_hash if self.entries else self.GENESIS
        seq = len(self.entries)
        timestamp = datetime.now(timezone.utc).isoformat()
        body = json.dumps(
            {"seq": seq, "timestamp_utc": timestamp, "action": action,
             "payload": payload, "prev_hash": prev},
            sort_keys=True, separators=(",", ":"),
        )
        entry = LedgerEntry(
            seq=seq, timestamp_utc=timestamp, action=action, payload=payload,
            prev_hash=prev,
            entry_hash=hashlib.sha256(body.encode()).hexdigest(),
        )
        self.entries.append(entry)
        return entry

    def verify(self) -> tuple[bool, Optional[str]]:
        prev = self.GENESIS
        for entry in self.entries:
            if entry.prev_hash != prev:
                return False, f"chain broken at seq {entry.seq}"
            body = json.dumps(
                {"seq": entry.seq, "timestamp_utc": entry.timestamp_utc,
                 "action": entry.action, "payload": entry.payload,
                 "prev_hash": entry.prev_hash},
                sort_keys=True, separators=(",", ":"),
            )
            if hashlib.sha256(body.encode()).hexdigest() != entry.entry_hash:
                return False, f"entry {entry.seq} has been altered"
            prev = entry.entry_hash
        return True, None


@dataclass
class ReviewItem:
    candidate: PromotionCandidate
    state: ReviewState = ReviewState.PENDING
    reviewer: Optional[str] = None
    rationale: Optional[str] = None


class ReviewQueue:
    """Expert curation, with every transition written to the ledger."""

    def __init__(self, ledger: Optional[EvidenceLedger] = None):
        self.ledger = ledger or EvidenceLedger()
        self.items: dict[tuple[str, str], ReviewItem] = {}

    def submit(self, candidate: PromotionCandidate) -> ReviewItem:
        if not candidate.eligible:
            raise ValueError(
                f"refusing to queue an ineligible candidate "
                f"({'; '.join(candidate.blocked_reasons)})"
            )
        key = (candidate.variant_key, candidate.drug)
        item = ReviewItem(candidate=candidate)
        self.items[key] = item
        self.ledger.append("submitted", asdict(candidate))
        return item

    def _transition(self, variant_key: str, drug: str, state: ReviewState,
                    reviewer: str, rationale: str) -> ReviewItem:
        key = (variant_key, drug)
        item = self.items.get(key)
        if item is None:
            raise KeyError(f"{variant_key}/{drug} is not in the review queue")
        item.state, item.reviewer, item.rationale = state, reviewer, rationale
        self.ledger.append(state.value, {
            "variant_key": variant_key, "drug": drug,
            "reviewer": reviewer, "rationale": rationale,
        })
        return item

    def approve(self, variant_key: str, drug: str, reviewer: str,
                rationale: str) -> ReviewItem:
        return self._transition(variant_key, drug, ReviewState.APPROVED,
                                reviewer, rationale)

    def reject(self, variant_key: str, drug: str, reviewer: str,
               rationale: str) -> ReviewItem:
        return self._transition(variant_key, drug, ReviewState.REJECTED,
                                reviewer, rationale)

    def roll_back(self, variant_key: str, drug: str, reviewer: str,
                  rationale: str) -> ReviewItem:
        """Reverse an approval after release."""
        return self._transition(variant_key, drug, ReviewState.ROLLED_BACK,
                                reviewer, rationale)

    def approved_entries(self) -> list[dict]:
        return [item.candidate.to_catalogue_entry()
                for item in self.items.values()
                if item.state is ReviewState.APPROVED]


# -- drift ----------------------------------------------------------------
@dataclass
class DriftReport:
    assessed: bool
    alarm: bool
    detail: str
    mean: Optional[float] = None


def assess_drift(recent_auroc: list[float], baseline: float = 0.90,
                 tolerance: float = 0.05, min_points: int = 5) -> DriftReport:
    """Flag predictive decay, and refuse to judge on too few points.

    The previous version returned False for an empty list, which reads as "no
    drift" when it means "not measured".
    """
    if len(recent_auroc) < min_points:
        return DriftReport(
            assessed=False, alarm=False,
            detail=(f"{len(recent_auroc)} measurement(s); need {min_points} "
                    f"before drift can be assessed. Not assessed is not the "
                    f"same as no drift."),
        )
    mean = sum(recent_auroc) / len(recent_auroc)
    drop = baseline - mean
    if drop >= tolerance:
        return DriftReport(
            assessed=True, alarm=True, mean=round(mean, 4),
            detail=(f"mean AUROC {mean:.3f} is {drop:.3f} below the "
                    f"{baseline:.2f} baseline (tolerance {tolerance:.2f})"),
        )
    return DriftReport(
        assessed=True, alarm=False, mean=round(mean, 4),
        detail=f"mean AUROC {mean:.3f} within tolerance of {baseline:.2f}",
    )
