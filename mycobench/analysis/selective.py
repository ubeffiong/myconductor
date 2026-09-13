"""Risk-coverage evaluation: scoring a predictor that is allowed to decline.

The mismatch this fixes
-----------------------
Every published TB-AMR predictor answers every query and is reported as a
single pooled AUC or accuracy. The system these predictions feed into declines
to answer when its evidence is insufficient. Training and deployment therefore
disagree about what the task is, and the reported metric cannot see the
disagreement.

A predictor allowed to abstain has two numbers, not one:

**coverage** — the fraction of queries it answered.
**selective risk** — the error rate *among the queries it answered*.

Neither is interpretable alone. Abstaining on everything difficult drives
selective risk to zero at trivial coverage; answering everything confidently
drives coverage to one at whatever risk the data imposes. The honest summary is
the curve traced as the abstention threshold moves, and the honest comparison
is **at matched coverage**.

The falsification test this module implements
---------------------------------------------
``beats_baseline_at_matched_coverage`` is the check that decides whether a
model earns its place: at the coverage the existing catalogue achieves, does
the model make fewer errors than the catalogue? If not, it has added
complexity and no information, however good its AUC looked.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from . import stats

#: Below this many answered queries a selective risk is not reported: an error
#: rate over three predictions is noise with a decimal point.
MIN_COVERED_FOR_RISK = 10


@dataclass(frozen=True)
class Prediction:
    """One prediction, its truth, and how confident the model was.

    ``confidence`` orders abstention: the threshold sweep answers the most
    confident queries first. ``predicted`` may be ``None`` for a model that
    abstains unconditionally on a query, which is then never covered at any
    threshold.
    """

    key: str
    predicted: Optional[str]
    truth: str
    confidence: float

    @property
    def answerable(self) -> bool:
        return self.predicted is not None

    @property
    def correct(self) -> bool:
        return self.answerable and self.predicted == self.truth


@dataclass(frozen=True)
class RiskCoveragePoint:
    threshold: float
    n_total: int
    n_covered: int
    n_errors: int

    @property
    def coverage(self) -> float:
        return self.n_covered / self.n_total if self.n_total else 0.0

    @property
    def selective_risk(self) -> Optional[float]:
        if self.n_covered < MIN_COVERED_FOR_RISK:
            return None
        return self.n_errors / self.n_covered

    def describe(self) -> str:
        risk = ("not reported (too few answered)"
                if self.selective_risk is None
                else f"{self.selective_risk:.3f}")
        return (f"threshold {self.threshold:.3f}: coverage "
                f"{self.coverage:.3f} ({self.n_covered}/{self.n_total}), "
                f"selective risk {risk}")


@dataclass
class RiskCoverageCurve:
    points: list[RiskCoveragePoint] = field(default_factory=list)
    n_total: int = 0
    n_unanswerable: int = 0

    @property
    def full_coverage_risk(self) -> Optional[float]:
        """Error rate if the model is forced to answer everything it can."""
        if not self.points:
            return None
        return max(self.points, key=lambda p: p.n_covered).selective_risk

    def risk_at_coverage(self, target: float) -> Optional[RiskCoveragePoint]:
        """The least-abstaining point that still reaches ``target`` coverage.

        Returns ``None`` when the model cannot reach that coverage at all —
        which is itself the answer, not a missing value.
        """
        eligible = [p for p in self.points
                    if p.coverage >= target and p.selective_risk is not None]
        return min(eligible, key=lambda p: p.coverage) if eligible else None

    def area_under_curve(self) -> Optional[float]:
        """Mean selective risk over coverage — lower is better.

        A single number for convenience only; it hides where on the curve the
        model is good, which is usually the thing that matters clinically.
        """
        usable = [p for p in self.points if p.selective_risk is not None]
        if len(usable) < 2:
            return None
        usable.sort(key=lambda p: p.coverage)
        area = 0.0
        span = 0.0
        for left, right in zip(usable, usable[1:]):
            width = right.coverage - left.coverage
            if width <= 0:
                continue
            area += width * (left.selective_risk + right.selective_risk) / 2
            span += width
        return area / span if span > 0 else None

    def describe(self) -> str:
        lines = [f"{self.n_total} query/queries, "
                 f"{self.n_unanswerable} never answerable"]
        area = self.area_under_curve()
        lines.append("area under risk-coverage: "
                     + (f"{area:.4f}" if area is not None else "not computable"))
        return "\n".join(lines)


def risk_coverage_curve(predictions: Iterable[Prediction],
                        steps: int = 50) -> RiskCoverageCurve:
    """Sweep the abstention threshold and record coverage against risk."""
    predictions = list(predictions)
    answerable = [p for p in predictions if p.answerable]
    curve = RiskCoverageCurve(
        n_total=len(predictions),
        n_unanswerable=len(predictions) - len(answerable))
    if not predictions:
        return curve

    confidences = sorted({p.confidence for p in answerable})
    if not confidences:
        return curve

    # Sweep from the most permissive threshold to the strictest, so coverage
    # falls monotonically and the curve reads left to right.
    thresholds = [confidences[0]]
    if len(confidences) > 1:
        stride = max(1, len(confidences) // max(1, steps - 1))
        thresholds.extend(confidences[i] for i in range(0, len(confidences), stride))
        thresholds.append(confidences[-1])

    seen: set[float] = set()
    for threshold in thresholds:
        if threshold in seen:
            continue
        seen.add(threshold)
        covered = [p for p in answerable if p.confidence >= threshold]
        errors = sum(1 for p in covered if not p.correct)
        curve.points.append(RiskCoveragePoint(
            threshold=threshold, n_total=len(predictions),
            n_covered=len(covered), n_errors=errors))
    curve.points.sort(key=lambda p: -p.coverage)
    return curve


@dataclass
class BaselineComparison:
    """Does the model beat a fixed-behaviour baseline at matched coverage?"""

    baseline_coverage: float
    baseline_risk: float
    matched: Optional[RiskCoveragePoint]
    verdict: str
    reasons: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.verdict == "beats-baseline"


def beats_baseline_at_matched_coverage(
        curve: RiskCoverageCurve,
        baseline_coverage: float,
        baseline_risk: float,
        margin: float = 0.0) -> BaselineComparison:
    """The falsification test for a VUS model.

    ``baseline_coverage`` and ``baseline_risk`` describe the incumbent — for
    TB-AMR that is the catalogue, which answers where it has a graded entry and
    abstains elsewhere. A model is only worth deploying if, at the coverage the
    catalogue already achieves, it makes fewer errors.
    """
    matched = curve.risk_at_coverage(baseline_coverage)
    if matched is None:
        return BaselineComparison(
            baseline_coverage, baseline_risk, None, "cannot-match-coverage",
            [f"the model never reaches {baseline_coverage:.0%} coverage with "
             f"enough answered queries to estimate a risk, so it cannot be "
             f"compared against the baseline on the baseline's own terms"])

    risk = matched.selective_risk
    if risk is None:
        return BaselineComparison(
            baseline_coverage, baseline_risk, matched, "underpowered",
            [f"only {matched.n_covered} answered query/queries at matched "
             f"coverage; below the {MIN_COVERED_FOR_RISK} needed to report a "
             f"risk"])

    if risk + margin < baseline_risk:
        return BaselineComparison(
            baseline_coverage, baseline_risk, matched, "beats-baseline",
            [f"selective risk {risk:.3f} at coverage {matched.coverage:.3f} "
             f"versus baseline {baseline_risk:.3f} at "
             f"{baseline_coverage:.3f}"])

    return BaselineComparison(
        baseline_coverage, baseline_risk, matched, "no-better-than-baseline",
        [f"selective risk {risk:.3f} at matched coverage is not below the "
         f"baseline's {baseline_risk:.3f}; the model adds complexity without "
         f"adding information"])


def catalogue_baseline(predictions: Iterable[Prediction]) -> tuple[float, float]:
    """Coverage and risk of a catalogue-style predictor over the same queries.

    A catalogue answers where it has an entry and abstains otherwise, with no
    threshold to tune — so it is a single point rather than a curve. Passing
    ``Prediction`` objects whose ``predicted`` is ``None`` for uncatalogued
    variants yields exactly that point.
    """
    predictions = list(predictions)
    if not predictions:
        return 0.0, 0.0
    answered = [p for p in predictions if p.answerable]
    coverage = len(answered) / len(predictions)
    if len(answered) < MIN_COVERED_FOR_RISK:
        return coverage, 0.0
    risk = sum(1 for p in answered if not p.correct) / len(answered)
    return coverage, risk


def bootstrap_risk_interval(predictions: Sequence[Prediction],
                            threshold: float) -> Optional[tuple[float, float]]:
    """Interval on selective risk at one threshold, over answered queries."""
    covered = [p for p in predictions
               if p.answerable and p.confidence >= threshold]
    if len(covered) < MIN_COVERED_FOR_RISK:
        return None
    errors = [0.0 if p.correct else 1.0 for p in covered]
    return stats.bootstrap_ci(errors, stats.mean,
                              min_n=MIN_COVERED_FOR_RISK)
