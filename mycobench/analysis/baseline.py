"""What the catalogue itself achieves — the number every later claim needs.

Nothing in this project could previously say how good the incumbent is. That
gap sits underneath several others: ``ModelRegistry`` refuses to approve a model
that cannot beat a named baseline, ``selective.beats_baseline_at_matched_coverage``
compares a risk-coverage curve against a baseline point, and neither had a
measured one to use. This module measures it.

The incumbent, stated precisely
-------------------------------
The WHO catalogue is not a classifier with a threshold to tune. It answers where
it holds a graded entry and abstains elsewhere, so it is a single
(coverage, error-rate) point rather than a curve. Concretely, for one drug and
one isolate:

* carries a variant the catalogue grades **resistant** for this drug → ``R``;
* carries none, **and every such variant's coordinate was examined** → ``S``;
* otherwise → **abstain**.

The middle rule is the one that matters and it is deliberately strict. Reporting
susceptibility because no resistant variant was *found* would repeat the defect
this codebase was rebuilt to remove: not looking is not the same as looking and
seeing nothing. ``assessed_variants`` carries which coordinates were actually
examined — evidence drawn from the reference calls in the VCF, not from the
variant list — so "we checked every place resistance could hide" is a claim with
evidence behind it. ``ASSESSED_FRACTION_FOR_SUSCEPTIBLE`` relaxes "every" only
because a single uncallable base in one rarely-graded locus should not force
abstention on an otherwise complete examination; it is reported so a reader can
see how much was relaxed.

Abstention is counted, not hidden
---------------------------------
Every isolate the catalogue declines is recorded with the reason it declined.
A baseline that reports only its accuracy on the questions it chose to answer
flatters itself; coverage and error rate are meaningless apart.

What this is not
----------------
This measures the catalogue against CRyPTIC binary phenotypes on the isolates
supplied. CRyPTIC and the WHO catalogue share underlying isolates, so this is
**not** an independent evaluation — it is the incumbent's home-ground
performance, which is the right bar for a challenger to clear and the wrong
number to quote as accuracy. Every result says so.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from .selective import MIN_COVERED_FOR_RISK, Prediction, catalogue_baseline
from .strata import DeterminantIndex, Isolate

#: Share of a drug's graded resistant coordinates that must have been examined
#: before "no resistant variant found" is allowed to mean susceptible. Below
#: this the catalogue abstains rather than guess.
ASSESSED_FRACTION_FOR_SUSCEPTIBLE = 0.95

#: Truth labels accepted from the reuse table. Anything else is not a phenotype.
TRUTH_LABELS = {"R": "R", "S": "S"}

#: Reasons the catalogue declined to answer, counted per drug.
ABSTAIN_INCOMPLETE = "loci not fully examined"
ABSTAIN_NO_DETERMINANTS = "no graded determinant for this drug"


@dataclass
class DrugBaseline:
    """The incumbent's single operating point for one drug."""

    drug: str
    n_isolates: int = 0
    n_answered: int = 0
    n_errors: int = 0
    n_resistant_truth: int = 0
    n_susceptible_truth: int = 0
    n_determinants: int = 0
    coverage: float = 0.0
    error_rate: float = 0.0
    abstained: dict[str, int] = field(default_factory=dict)
    #: Errors split by direction. ``called_S_truth_R`` is the one that matters:
    #: a drug wrongly cleared for use, which is the failure mode this whole
    #: codebase is built to avoid.
    errors_by_kind: dict[str, int] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    @property
    def estimable(self) -> bool:
        """Enough answered queries for the error rate to mean anything."""
        return self.n_answered >= MIN_COVERED_FOR_RISK

    @property
    def n_false_susceptible(self) -> int:
        return self.errors_by_kind.get("called_S_truth_R", 0)

    def as_registry_baseline(self, source: str) -> Optional[dict]:
        """The block ``ModelRegistry`` requires on a performance row.

        ``None`` when the point is not estimable, so an unmeasurable baseline
        cannot be quoted as though it were a measured one.
        """
        if not self.estimable:
            return None
        return {"source": source,
                "coverage": round(self.coverage, 4),
                "error_rate": round(self.error_rate, 4)}

    def describe(self) -> str:
        if not self.estimable:
            return (f"{self.drug}: answered {self.n_answered} of "
                    f"{self.n_isolates}; fewer than {MIN_COVERED_FOR_RISK} "
                    f"answers, so no error rate is reportable")
        false_s = self.n_false_susceptible
        return (f"{self.drug}: answers {self.coverage:.0%} of "
                f"{self.n_isolates} isolate(s), errs on {self.error_rate:.1%} "
                f"of those ({self.n_errors}/{self.n_answered}; "
                f"{false_s} called susceptible but resistant)")


def _truth_for(row: dict, code: str) -> Optional[str]:
    return TRUTH_LABELS.get((row.get(f"{code}_BINARY_PHENOTYPE") or "").strip().upper())


def catalogue_predictions(
        isolates: Sequence[Isolate],
        truths: dict[str, str],
        determinants: Iterable[str],
        assessed_fraction: float = ASSESSED_FRACTION_FOR_SUSCEPTIBLE,
) -> tuple[list[Prediction], dict[str, int]]:
    """Catalogue answers for one drug, plus why it abstained where it did.

    ``truths`` maps isolate id to ``R``/``S``. ``determinants`` are the variant
    labels the catalogue grades resistant for this drug.
    """
    determinants = set(determinants)
    predictions: list[Prediction] = []
    abstained: dict[str, int] = {}

    for isolate in isolates:
        truth = truths.get(isolate.isolate_id)
        if truth is None:
            continue  # no phenotype: not a question the catalogue was asked
        if not determinants:
            abstained[ABSTAIN_NO_DETERMINANTS] = \
                abstained.get(ABSTAIN_NO_DETERMINANTS, 0) + 1
            predictions.append(Prediction(isolate.isolate_id, None, truth, 1.0))
            continue

        carried = determinants & set(isolate.genotype)
        if carried:
            # A graded resistant variant is present. The catalogue answers.
            predictions.append(Prediction(isolate.isolate_id, "R", truth, 1.0))
            continue

        # No resistant variant found. That is only susceptibility if we looked
        # everywhere one could have been.
        assessed = isolate.assessed_variants
        if assessed is None:
            examined = 0.0
        else:
            examined = len(determinants & set(assessed)) / len(determinants)
        if examined >= assessed_fraction:
            predictions.append(Prediction(isolate.isolate_id, "S", truth, 1.0))
        else:
            abstained[ABSTAIN_INCOMPLETE] = \
                abstained.get(ABSTAIN_INCOMPLETE, 0) + 1
            predictions.append(Prediction(isolate.isolate_id, None, truth, 1.0))

    return predictions, abstained


def measure_drug(drug: str, isolates: Sequence[Isolate],
                 truths: dict[str, str], determinants: Iterable[str],
                 assessed_fraction: float = ASSESSED_FRACTION_FOR_SUSCEPTIBLE
                 ) -> DrugBaseline:
    determinants = set(determinants)
    predictions, abstained = catalogue_predictions(
        isolates, truths, determinants, assessed_fraction)
    coverage, error_rate = catalogue_baseline(predictions)

    answered = [p for p in predictions if p.answerable]
    kinds: dict[str, int] = {}
    for p in answered:
        if not p.correct:
            kinds[f"called_{p.predicted}_truth_{p.truth}"] = \
                kinds.get(f"called_{p.predicted}_truth_{p.truth}", 0) + 1

    baseline = DrugBaseline(
        drug=drug,
        n_isolates=len(predictions),
        n_answered=len(answered),
        n_errors=sum(1 for p in answered if not p.correct),
        n_resistant_truth=sum(1 for p in predictions if p.truth == "R"),
        n_susceptible_truth=sum(1 for p in predictions if p.truth == "S"),
        n_determinants=len(determinants),
        coverage=coverage,
        error_rate=error_rate if len(answered) >= MIN_COVERED_FOR_RISK else 0.0,
        abstained=abstained,
        errors_by_kind=kinds,
    )
    if not baseline.estimable:
        baseline.notes.append(
            f"only {len(answered)} answered quer(ies); an error rate over "
            f"fewer than {MIN_COVERED_FOR_RISK} describes the sample, not the "
            f"catalogue")
    if baseline.n_resistant_truth == 0:
        baseline.notes.append(
            "no resistant isolate in this sample, so the error rate says "
            "nothing about sensitivity")
    return baseline


def measure(isolates: Sequence[Isolate], rows: Iterable[dict],
            index: DeterminantIndex, drugs: Optional[Iterable[str]] = None,
            assessed_fraction: float = ASSESSED_FRACTION_FOR_SUSCEPTIBLE
            ) -> dict[str, DrugBaseline]:
    """Measure the catalogue across every drug with a binary phenotype."""
    from ..phenotypes import DRUG_CODES

    wanted = set(drugs) if drugs else set(DRUG_CODES.values())
    by_run = {f"cr_{(r.get('ENA_RUN') or '').strip()}": r for r in rows}
    results: dict[str, DrugBaseline] = {}

    for code, drug in sorted(DRUG_CODES.items(), key=lambda kv: kv[1]):
        if drug not in wanted:
            continue
        truths = {}
        for isolate in isolates:
            row = by_run.get(isolate.isolate_id)
            if row is None:
                continue
            truth = _truth_for(row, code)
            if truth is not None:
                truths[isolate.isolate_id] = truth
        determinants = (index.by_drug_label.get(drug, set())
                        | index.by_drug_coordinate.get(drug, set()))
        results[drug] = measure_drug(drug, isolates, truths, determinants,
                                     assessed_fraction)
    return results


BASELINE_COLUMNS = (
    "drug", "n_isolates", "n_answered", "coverage", "n_errors", "error_rate",
    "n_false_susceptible", "n_resistant_truth", "n_susceptible_truth",
    "n_determinants", "estimable", "abstained", "notes",
)


def baseline_rows(results: dict[str, DrugBaseline]) -> list[dict]:
    rows = []
    for drug, b in sorted(results.items()):
        rows.append({
            "drug": drug, "n_isolates": b.n_isolates,
            "n_answered": b.n_answered,
            "coverage": f"{b.coverage:.4f}",
            "n_errors": b.n_errors,
            "error_rate": (f"{b.error_rate:.4f}" if b.estimable else ""),
            "n_false_susceptible": b.n_false_susceptible,
            "n_resistant_truth": b.n_resistant_truth,
            "n_susceptible_truth": b.n_susceptible_truth,
            "n_determinants": b.n_determinants,
            "estimable": "yes" if b.estimable else "no",
            "abstained": "; ".join(f"{k}: {v}" for k, v in sorted(b.abstained.items())),
            "notes": " | ".join(b.notes),
        })
    return rows


CAVEAT = (
    "CRyPTIC and the WHO catalogue share underlying isolates, so this is the "
    "incumbent's home-ground performance: the bar a challenger must clear, not "
    "an independent accuracy estimate for the catalogue."
)
