"""Metrics for the two tracks, with abstention treated as a first-class outcome.

The design decision that shapes this module
-------------------------------------------
Myconductor can decline to call a drug. A naive confusion matrix has nowhere to
put that, and the usual workaround — folding "not assessed" into "susceptible"
— is the exact defect the interpretation engine was rebuilt to remove. Putting
it back in the scoring would hide the failure it was designed to surface.

So every accuracy figure here is **conditional on a call being made**, and is
reported beside the **call rate**. Neither number means anything alone: a tool
that abstains on every hard isolate posts perfect accuracy on the easy ones,
and a tool that calls everything confidently posts a high call rate with bad
errors. The pre-registered targets require both.

Error names follow clinical method-comparison usage rather than machine
learning usage, because that is what the numbers will be read against:

``VME`` very major error — predicted susceptible, phenotypically resistant.
``ME``  major error — predicted resistant, phenotypically susceptible.

Intervals
---------
Wilson score intervals, which behave sensibly at the proportions that actually
occur here (0 of 12, 19 of 19). A normal-approximation interval gives
impossible bounds at those values. Where a denominator is below the
pre-registered minimum, no interval is produced at all — an interval from four
isolates implies a precision the data do not have.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Optional

from .thresholds import DrugTarget, target_for

#: Myconductor verdicts that constitute a call. Everything else is abstention.
CALLED = ("resistant", "susceptible")


def wilson_interval(successes: int, total: int,
                    z: float = 1.959963985) -> Optional[tuple[float, float]]:
    """95% Wilson score interval, or None when there is nothing to bound."""
    if total <= 0:
        return None
    p = successes / total
    denominator = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    margin = (z * math.sqrt(p * (1 - p) / total
                            + z * z / (4 * total * total)) / denominator)
    return max(0.0, centre - margin), min(1.0, centre + margin)


@dataclass
class AccuracyResult:
    """Accuracy for one drug against a laboratory phenotype."""

    drug: str
    n_evaluable: int             # isolates with an accepted-quality phenotype
    n_called: int
    n_abstained: int
    true_positive: int = 0       # predicted R, phenotype R
    false_positive: int = 0      # predicted R, phenotype S  -> ME
    true_negative: int = 0       # predicted S, phenotype S
    false_negative: int = 0      # predicted S, phenotype R  -> VME
    abstained_resistant: int = 0
    abstained_susceptible: int = 0
    target: Optional[DrugTarget] = None
    notes: list[str] = field(default_factory=list)

    # -- rates, all conditional on a call being made ---------------------
    @property
    def call_rate(self) -> Optional[float]:
        return self.n_called / self.n_evaluable if self.n_evaluable else None

    @property
    def n_phenotype_resistant(self) -> int:
        return self.true_positive + self.false_negative

    @property
    def n_phenotype_susceptible(self) -> int:
        return self.true_negative + self.false_positive

    @property
    def sensitivity(self) -> Optional[float]:
        d = self.n_phenotype_resistant
        return self.true_positive / d if d else None

    @property
    def specificity(self) -> Optional[float]:
        d = self.n_phenotype_susceptible
        return self.true_negative / d if d else None

    @property
    def ppv(self) -> Optional[float]:
        d = self.true_positive + self.false_positive
        return self.true_positive / d if d else None

    @property
    def npv(self) -> Optional[float]:
        d = self.true_negative + self.false_negative
        return self.true_negative / d if d else None

    @property
    def vme_rate(self) -> Optional[float]:
        """Predicted susceptible among phenotypically resistant isolates."""
        d = self.n_phenotype_resistant
        return self.false_negative / d if d else None

    @property
    def me_rate(self) -> Optional[float]:
        """Predicted resistant among phenotypically susceptible isolates."""
        d = self.n_phenotype_susceptible
        return self.false_positive / d if d else None

    def interval(self, name: str) -> Optional[tuple[float, float]]:
        pairs = {
            "sensitivity": (self.true_positive, self.n_phenotype_resistant),
            "specificity": (self.true_negative, self.n_phenotype_susceptible),
            "vme_rate": (self.false_negative, self.n_phenotype_resistant),
            "me_rate": (self.false_positive, self.n_phenotype_susceptible),
            "call_rate": (self.n_called, self.n_evaluable),
        }
        if name not in pairs:
            raise KeyError(f"no interval defined for {name!r}")
        successes, total = pairs[name]
        if not self.powered:
            return None
        return wilson_interval(successes, total)

    # -- verdict against the pre-registration ----------------------------
    @property
    def powered(self) -> bool:
        """Are there enough evaluable isolates to make any claim at all?"""
        minimum = self.target.min_evaluable if self.target else 20
        return self.n_evaluable >= minimum

    def verdict(self) -> tuple[str, list[str]]:
        """``(outcome, reasons)`` where outcome is pass / fail / underpowered.

        Underpowered is deliberately distinct from fail: not having measured
        something is not the same as having measured it and found it wanting.
        """
        if self.target is None:
            return "not-targeted", [
                f"{self.drug} has no pre-registered target; reported for "
                f"information only"]
        if not self.powered:
            return "underpowered", [
                f"{self.n_evaluable} evaluable isolate(s), below the "
                f"pre-registered minimum of {self.target.min_evaluable}; no "
                f"pass or fail is claimed"]

        # VME is checked, not sensitivity, because they are the same quantity
        # (VME = 1 - sensitivity) and checking both would report one breach
        # twice. Likewise ME rather than specificity. The derived ceilings come
        # from the declared sensitivity/specificity targets.
        reasons = []
        checks = (
            ("VME", self.vme_rate, self.target.max_vme, "above"),
            ("ME", self.me_rate, self.target.max_me, "above"),
            ("call rate", self.call_rate, self.target.min_call_rate, "below"),
        )
        for name, observed, bound, direction in checks:
            if observed is None:
                reasons.append(f"{name} not estimable (no isolates in its "
                               f"denominator)")
                continue
            breached = (observed > bound if direction == "above"
                        else observed < bound)
            if breached:
                reasons.append(
                    f"{name} {observed:.3f} is {direction} the pre-registered "
                    f"bound {bound:.3f}")
        if any("is above" in r or "is below" in r for r in reasons):
            return "fail", reasons
        return "pass", reasons or ["every pre-registered bound met"]


def score_accuracy(observations: Iterable[tuple[str, str]], drug: str,
                   n_evaluable: Optional[int] = None) -> AccuracyResult:
    """Score one drug from ``(predicted_call, phenotype)`` pairs.

    ``predicted_call`` is a Myconductor ``Call`` value; ``phenotype`` is ``R``
    or ``S``. Pairs whose phenotype is neither are not evaluable and must be
    filtered out before calling this.
    """
    observations = list(observations)
    result = AccuracyResult(
        drug=drug,
        n_evaluable=n_evaluable if n_evaluable is not None else len(observations),
        n_called=0, n_abstained=0, target=target_for(drug))

    for predicted, phenotype in observations:
        phenotype = (phenotype or "").strip().upper()
        predicted = (predicted or "").strip().lower()
        if phenotype not in ("R", "S"):
            continue
        if predicted not in CALLED:
            result.n_abstained += 1
            if phenotype == "R":
                result.abstained_resistant += 1
            else:
                result.abstained_susceptible += 1
            continue
        result.n_called += 1
        if predicted == "resistant":
            if phenotype == "R":
                result.true_positive += 1
            else:
                result.false_positive += 1
        else:
            if phenotype == "R":
                result.false_negative += 1
            else:
                result.true_negative += 1

    if result.abstained_resistant and result.n_phenotype_resistant:
        result.notes.append(
            f"{result.abstained_resistant} phenotypically resistant isolate(s) "
            f"were not called and are excluded from sensitivity; folding them "
            f"into 'susceptible' would understate VME")
    if result.n_abstained and not result.n_called:
        result.notes.append(
            "every evaluable isolate was abstained on; no accuracy figure "
            "exists for this drug")
    return result


@dataclass
class ConcordanceResult:
    """Agreement between two engines on one drug. Not an accuracy measure."""

    drug: str
    n_compared: int
    agree: int = 0
    disagree: int = 0
    only_a_called: int = 0
    only_b_called: int = 0
    neither_called: int = 0
    disagreements: list[tuple[str, str, str]] = field(default_factory=list)

    @property
    def agreement_rate(self) -> Optional[float]:
        both = self.agree + self.disagree
        return self.agree / both if both else None

    def describe(self) -> str:
        if self.agree + self.disagree == 0:
            return (f"{self.drug}: no isolate had a call from both engines; "
                    f"nothing to reconcile")
        return (f"{self.drug}: {self.agree}/{self.agree + self.disagree} "
                f"concordant where both called "
                f"({self.only_a_called + self.only_b_called} single-engine, "
                f"{self.neither_called} neither)")


def score_concordance(observations: Iterable[tuple[str, str, str]],
                      drug: str) -> ConcordanceResult:
    """Score one drug from ``(sample_id, call_a, call_b)`` triples.

    Agreement is computed only over isolates both engines called. An engine
    abstaining is recorded separately, never counted as agreement — two tools
    both declining to answer have not agreed about anything.
    """
    observations = list(observations)
    result = ConcordanceResult(drug=drug, n_compared=len(observations))
    for sample_id, call_a, call_b in observations:
        a = (call_a or "").strip().lower()
        b = (call_b or "").strip().lower()
        a_called, b_called = a in CALLED, b in CALLED
        if a_called and b_called:
            if a == b:
                result.agree += 1
            else:
                result.disagree += 1
                result.disagreements.append((sample_id, a, b))
        elif a_called:
            result.only_a_called += 1
        elif b_called:
            result.only_b_called += 1
        else:
            result.neither_called += 1
    return result


@dataclass
class SpeciesControlResult:
    """Did the pipeline correctly refuse a non-MTBC genome?"""

    sample_id: str
    organism: str
    refused: bool
    offending_calls: list[tuple[str, str]] = field(default_factory=list)
    detail: str = ""

    @property
    def passed(self) -> bool:
        return self.refused and not self.offending_calls


def score_species_control(sample_id: str, organism: str,
                          drug_calls: Iterable[tuple[str, str]]
                          ) -> SpeciesControlResult:
    """A drug call on a non-MTBC isolate is a failure, not a result."""
    offending = [(drug, call) for drug, call in drug_calls
                 if (call or "").strip().lower() in CALLED]
    if offending:
        return SpeciesControlResult(
            sample_id, organism, refused=False, offending_calls=offending,
            detail=(f"{len(offending)} drug call(s) produced for a "
                    f"{organism} genome against the MTBC profile: "
                    + ", ".join(f"{d}={c}" for d, c in offending[:4])))
    return SpeciesControlResult(
        sample_id, organism, refused=True,
        detail=f"no drug call produced for this {organism} genome, as required")
