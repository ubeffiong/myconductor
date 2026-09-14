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

from ..metrics import AccuracyResult, score_accuracy, score_accuracy_stratified
from .selective import MIN_COVERED_FOR_RISK, Prediction, catalogue_baseline
from .strata import DeterminantIndex, Isolate

#: Myconductor call values the catalogue's R/S answers correspond to, so the
#: baseline can be scored by the same code as any engine. Reusing
#: ``metrics.score_accuracy`` rather than recomputing here is deliberate: it
#: already excludes abstained-resistant isolates from sensitivity and says so,
#: which is the difference between an honest sensitivity and one inflated by
#: declining the hard cases.
_CALL_FOR = {"R": "resistant", "S": "susceptible"}
_ABSTAINED = "not_assessed"

#: Share of a drug's graded resistant coordinates that must have been examined
#: before "no resistant variant found" is allowed to mean susceptible. Below
#: this the catalogue abstains rather than guess.
ASSESSED_FRACTION_FOR_SUSCEPTIBLE = 0.95

#: Truth labels accepted from the reuse table. Anything else is not a phenotype.
TRUTH_LABELS = {"R": "R", "S": "S"}

#: Below this many isolates on one side of the truth, the point estimate that
#: side supports is not readable on its own. Sensitivity rests only on the
#: resistant isolates and specificity only on the susceptible ones, so a large
#: cohort says nothing about whether either is well determined.
MIN_FOR_DIRECTIONAL_CLAIM = 10

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
    #: Sensitivity, specificity, PPV, NPV, VME and ME with Wilson intervals,
    #: scored by the same code that scores any engine. A single composite error
    #: rate cannot distinguish a tool that misses resistance from one that
    #: over-calls it, and those have opposite clinical consequences.
    accuracy: Optional[AccuracyResult] = None
    #: The same, split by lineage. Empty when no isolate carries one — the
    #: CRyPTIC reuse table has no lineage column, so this stays empty unless
    #: lineages are supplied from outside.
    by_lineage: dict[str, AccuracyResult] = field(default_factory=dict)
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


def locus_of(label: str) -> str:
    """The gene a determinant label names, e.g. ``rpoB_p.Ser450Leu`` -> ``rpoB``.

    Labels are ``gene_change``; the change may itself contain underscores, so
    only the first segment is the locus.
    """
    return label.split("_", 1)[0] if label else ""


def _locus_fraction(determinants: set, assessed: Optional[frozenset]) -> float:
    """Share of a drug's catalogued **loci** that were examined.

    The denominator is the set of genes the drug's determinants sit in, not the
    number of determinants, so a drug with 136 graded variants in one gene is
    not held to a stricter standard than one with eleven across two.
    """
    loci = {locus_of(label) for label in determinants if locus_of(label)}
    if not loci:
        return 0.0
    if not assessed:
        return 0.0
    seen = {locus_of(label) for label in assessed if locus_of(label)}
    return len(loci & seen) / len(loci)


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
        # everywhere one could have been — and "everywhere" is a question about
        # loci, not about enumerating catalogued variants.
        #
        # Counting variant coordinates instead made the rule silently stricter
        # for better-studied drugs: rifampicin has 136 graded determinants and
        # isoniazid 143, against eleven for amikacin, so a 95% variant-level
        # threshold was essentially unreachable for the first two. The
        # catalogue then never answered "susceptible" for them at all, every
        # answered isolate was an R call, and specificity came out at exactly
        # zero — an artifact of the denominator, not a property of the
        # catalogue. Locus coverage is what the callable mask actually
        # establishes, and it is what ``profile.required_loci`` means.
        assessed = isolate.assessed_variants
        examined = _locus_fraction(determinants, assessed)
        if examined >= assessed_fraction:
            predictions.append(Prediction(isolate.isolate_id, "S", truth, 1.0))
        else:
            abstained[ABSTAIN_INCOMPLETE] = \
                abstained.get(ABSTAIN_INCOMPLETE, 0) + 1
            predictions.append(Prediction(isolate.isolate_id, None, truth, 1.0))

    return predictions, abstained


def _scored(drug: str, predictions: Sequence[Prediction],
            lineages: Optional[dict[str, str]] = None
            ) -> tuple[AccuracyResult, dict[str, AccuracyResult]]:
    """Score the catalogue's answers as if it were any other predictor."""
    pairs = [(_CALL_FOR.get(p.predicted, _ABSTAINED), p.truth)
             for p in predictions]
    pooled = score_accuracy(pairs, drug)
    if not lineages:
        return pooled, {}
    triples = [(_CALL_FOR.get(p.predicted, _ABSTAINED), p.truth,
                lineages.get(p.key, "unknown")) for p in predictions]
    stratified = score_accuracy_stratified(triples, drug)
    # A single "unknown" bucket is not stratification, it is the absence of it.
    if set(stratified) <= {"unknown"}:
        return pooled, {}
    return pooled, stratified


def measure_drug(drug: str, isolates: Sequence[Isolate],
                 truths: dict[str, str], determinants: Iterable[str],
                 assessed_fraction: float = ASSESSED_FRACTION_FOR_SUSCEPTIBLE,
                 lineages: Optional[dict[str, str]] = None
                 ) -> DrugBaseline:
    determinants = set(determinants)
    predictions, abstained = catalogue_predictions(
        isolates, truths, determinants, assessed_fraction)
    coverage, error_rate = catalogue_baseline(predictions)
    accuracy, by_lineage = _scored(drug, predictions, lineages)

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
        accuracy=accuracy,
        by_lineage=by_lineage,
    )
    baseline.notes.extend(accuracy.notes)

    # ``powered`` in metrics.py asks whether the *cohort* is large enough. It
    # does not ask whether the resistant isolates within it are, and sensitivity
    # rests only on those. Bedaquiline here has four: a sensitivity of 0% over
    # four isolates carries a 95% interval reaching 49%, which is a different
    # statement from the same 0% over twenty-one. Both would otherwise read
    # "powered: yes" beside a headline figure.
    if 0 < accuracy.n_phenotype_resistant < MIN_FOR_DIRECTIONAL_CLAIM:
        baseline.notes.append(
            f"sensitivity and VME rest on {accuracy.n_phenotype_resistant} "
            f"resistant isolate(s); read the interval, not the point estimate")
    if 0 < accuracy.n_phenotype_susceptible < MIN_FOR_DIRECTIONAL_CLAIM:
        baseline.notes.append(
            f"specificity and ME rest on "
            f"{accuracy.n_phenotype_susceptible} susceptible isolate(s); read "
            f"the interval, not the point estimate")
    if not baseline.estimable:
        baseline.notes.append(
            f"only {len(answered)} answered quer(ies); an error rate over "
            f"fewer than {MIN_COVERED_FOR_RISK} describes the sample, not the "
            f"catalogue")
    if baseline.n_resistant_truth == 0:
        baseline.notes.append(
            "no resistant isolate in this sample, so the error rate says "
            "nothing about sensitivity")

    # A namespace mismatch abstains on everything and looks exactly like a
    # cohort nobody sequenced deeply. It is not: it means the determinant set
    # and the isolate data are keyed differently — labels against coordinate
    # keys, say — so nothing could ever match. That failed silently once and
    # produced a whole baseline table of plausible-looking wrong numbers.
    if determinants and isolates and baseline.n_answered == 0:
        seen = set()
        for isolate in isolates:
            seen |= set(isolate.genotype)
            if isolate.assessed_variants:
                seen |= set(isolate.assessed_variants)
        if seen and not (seen & determinants):
            baseline.notes.append(
                f"no isolate carried or assessed ANY of the {len(determinants)} "
                f"determinant(s) for this drug; the determinant set and the "
                f"isolate data appear to be keyed differently (e.g. gene "
                f"labels against coordinate keys) rather than the cohort "
                f"genuinely lacking coverage")
    return baseline


def measure(isolates: Sequence[Isolate], rows: Iterable[dict],
            index: DeterminantIndex, drugs: Optional[Iterable[str]] = None,
            assessed_fraction: float = ASSESSED_FRACTION_FOR_SUSCEPTIBLE,
            lineages: Optional[dict[str, str]] = None
            ) -> dict[str, DrugBaseline]:
    """Measure the catalogue across every drug with a binary phenotype."""
    from ..phenotypes import DRUG_CODES

    wanted = set(drugs) if drugs else set(DRUG_CODES.values())
    if lineages is None:
        lineages = {i.isolate_id: i.lineage for i in isolates if i.lineage}
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
        # Labels only, and deliberately not unioned with the coordinate keys.
        # ``DeterminantIndex`` carries both spellings of every determinant, but
        # ``load_genotypes`` produces label-keyed genotypes and assessed sets,
        # so mixing the namespaces inflates the denominator without adding a
        # single matchable entry: isoniazid holds 143 labels against 3,517
        # coordinate keys, which would cap the achievable assessed fraction at
        # about 4% and abstain on every isolate. The determinant set must be in
        # the same namespace as the isolate data being measured.
        determinants = index.by_drug_label.get(drug, set())
        results[drug] = measure_drug(drug, isolates, truths, determinants,
                                     assessed_fraction, lineages)
    return results


BASELINE_COLUMNS = (
    "drug", "n_isolates", "n_answered", "coverage", "n_errors", "error_rate",
    "sensitivity", "sensitivity_ci", "specificity", "specificity_ci",
    "ppv", "npv", "vme_rate", "me_rate",
    "n_false_susceptible", "n_false_resistant",
    "n_resistant_truth", "n_susceptible_truth", "resistance_prevalence",
    "abstained_resistant", "n_determinants", "estimable", "cohort_powered",
    # Sensitivity rests only on resistant isolates and specificity only on
    # susceptible ones, so cohort size says nothing about whether either is
    # well determined. Both counts are shown beside the rates they support.
    "sensitivity_basis_n", "specificity_basis_n",
    "abstained", "notes",
)


def _rate(value: Optional[float]) -> str:
    return "" if value is None else f"{value:.4f}"


def _interval(result: Optional[AccuracyResult], name: str) -> str:
    if result is None:
        return ""
    span = result.interval(name)
    return "" if span is None else f"{span[0]:.4f}-{span[1]:.4f}"


def baseline_rows(results: dict[str, DrugBaseline]) -> list[dict]:
    rows = []
    for drug, b in sorted(results.items()):
        accuracy = b.accuracy
        prevalence = (b.n_resistant_truth / b.n_isolates
                      if b.n_isolates else None)
        rows.append({
            "drug": drug, "n_isolates": b.n_isolates,
            "n_answered": b.n_answered,
            "coverage": f"{b.coverage:.4f}",
            "n_errors": b.n_errors,
            "error_rate": (f"{b.error_rate:.4f}" if b.estimable else ""),
            "sensitivity": _rate(accuracy.sensitivity if accuracy else None),
            "sensitivity_ci": _interval(accuracy, "sensitivity"),
            "specificity": _rate(accuracy.specificity if accuracy else None),
            "specificity_ci": _interval(accuracy, "specificity"),
            # Prevalence-dependent, and this cohort is not a clinical
            # population; the manifest carries the caveat alongside.
            "ppv": _rate(accuracy.ppv if accuracy else None),
            "npv": _rate(accuracy.npv if accuracy else None),
            "vme_rate": _rate(accuracy.vme_rate if accuracy else None),
            "me_rate": _rate(accuracy.me_rate if accuracy else None),
            "n_false_susceptible": b.n_false_susceptible,
            "n_false_resistant": (accuracy.false_positive if accuracy else 0),
            "n_resistant_truth": b.n_resistant_truth,
            "n_susceptible_truth": b.n_susceptible_truth,
            "resistance_prevalence": _rate(prevalence),
            # Resistant isolates the catalogue declined. Excluded from
            # sensitivity by construction, so shown beside it: a predictor can
            # always look sensitive by abstaining on the hard ones.
            "abstained_resistant": (accuracy.abstained_resistant
                                    if accuracy else 0),
            "n_determinants": b.n_determinants,
            "estimable": "yes" if b.estimable else "no",
            # Renamed from "powered": it asks about the cohort, not about the
            # side of the truth each rate actually rests on.
            "cohort_powered": ("yes" if accuracy and accuracy.powered else "no"),
            "sensitivity_basis_n": (accuracy.n_phenotype_resistant
                                    if accuracy else 0),
            "specificity_basis_n": (accuracy.n_phenotype_susceptible
                                    if accuracy else 0),
            "abstained": "; ".join(f"{k}: {v}" for k, v in sorted(b.abstained.items())),
            "notes": " | ".join(b.notes),
        })
    return rows


LINEAGE_COLUMNS = (
    "drug", "lineage", "n_evaluable", "n_called", "sensitivity",
    "specificity", "vme_rate", "me_rate", "n_phenotype_resistant",
    "n_phenotype_susceptible", "powered",
)


def lineage_rows(results: dict[str, DrugBaseline]) -> list[dict]:
    """Per-lineage performance, where lineages were supplied.

    A pooled figure can be carried entirely by one lineage, and a tool that
    works on Lineage 2 and fails on Lineage 4 is a different product in each
    setting. Empty when no lineage is known, which is itself the finding: the
    CRyPTIC reuse table carries no lineage column, so this stays empty unless
    lineages are supplied from a metadata file.
    """
    rows = []
    for drug, baseline in sorted(results.items()):
        for lineage, score in sorted(baseline.by_lineage.items()):
            rows.append({
                "drug": drug, "lineage": lineage,
                "n_evaluable": score.n_evaluable, "n_called": score.n_called,
                "sensitivity": _rate(score.sensitivity),
                "specificity": _rate(score.specificity),
                "vme_rate": _rate(score.vme_rate),
                "me_rate": _rate(score.me_rate),
                "n_phenotype_resistant": score.n_phenotype_resistant,
                "n_phenotype_susceptible": score.n_phenotype_susceptible,
                "powered": "yes" if score.powered else "no",
            })
    return rows


PPV_NPV_CAVEAT = (
    "PPV and NPV depend on the prevalence of resistance in the population "
    "tested, and this cohort is not a clinical population: CRyPTIC was "
    "assembled to contain resistance, so its prevalence is far above what a "
    "routine diagnostic service sees. Sensitivity and specificity transfer "
    "between populations; PPV and NPV do not. The prevalence each pair was "
    "computed at is in resistance_prevalence, and they must be recomputed at "
    "a setting's own prevalence before they mean anything there."
)


#: Drugs worth reporting on even when this pairing cannot measure them. A
#: benchmark that silently lists only what it happens to cover reads as a
#: complete panel; naming the absences, and which side each is missing from,
#: is the difference between a gap and an oversight. BPaL/M components are
#: here because they are the backbone of modern drug-resistant TB treatment
#: and are exactly where genomic prediction is weakest.
DRUGS_OF_INTEREST = (
    "bedaquiline", "pretomanid", "linezolid", "moxifloxacin",   # BPaL/M
    "clofazimine", "delamanid",
    "rifampicin", "isoniazid", "ethambutol", "pyrazinamide",
    "levofloxacin", "amikacin", "kanamycin", "streptomycin",
    "capreomycin", "ethionamide", "cycloserine", "rifabutin",
)

MEASURABILITY_COLUMNS = ("drug", "in_catalogue", "in_phenotype_source",
                         "measurable", "missing_side", "consequence")


def measurability(catalogue_drugs: Iterable[str],
                  phenotype_drugs: Iterable[str],
                  drugs: Iterable[str] = DRUGS_OF_INTEREST) -> list[dict]:
    """Which drugs this pairing can evaluate at all, and which side is missing.

    Measuring the catalogue needs both a **genotypic** side (graded
    determinants to make a call from) and a **phenotypic** side (a laboratory
    result to score it against). A drug missing either cannot be evaluated,
    and the two absences have different remedies: no catalogue entry means
    nothing can be predicted; no phenotype means a prediction cannot be
    checked. Reporting them apart tells a reader which problem to go and solve.

    Pretomanid is the case that motivates this. It is absent from the WHO
    catalogue's fifteen drugs *and* from CRyPTIC's thirteen MIC columns, so no
    pairing of these two sources can say anything about it — yet it is a
    quarter of the BPaL regimen. Omitting it silently reads as an oversight;
    naming it reads as the infrastructure gap it is.
    """
    catalogue_drugs = {d.lower() for d in catalogue_drugs}
    phenotype_drugs = {d.lower() for d in phenotype_drugs}
    rows = []
    for drug in sorted(set(drugs)):
        genotypic = drug.lower() in catalogue_drugs
        phenotypic = drug.lower() in phenotype_drugs
        if genotypic and phenotypic:
            missing, consequence = "", "measurable from these two sources"
        elif genotypic:
            missing, consequence = (
                "phenotype",
                "a call can be made but not scored; supply a phenotype source")
        elif phenotypic:
            missing, consequence = (
                "catalogue",
                "a phenotype exists but nothing predicts it; supply graded "
                "determinants")
        else:
            missing, consequence = (
                "both",
                "no pairing of these sources can evaluate this drug at all")
        rows.append({
            "drug": drug,
            "in_catalogue": "yes" if genotypic else "no",
            "in_phenotype_source": "yes" if phenotypic else "no",
            "measurable": "yes" if (genotypic and phenotypic) else "no",
            "missing_side": missing,
            "consequence": consequence,
        })
    return rows


def isolates_needed(baseline: "DrugBaseline",
                    min_resistant: int = MIN_COVERED_FOR_RISK) -> Optional[int]:
    """Roughly how many more isolates would make this drug estimable.

    Extrapolated from the resistance prevalence observed in this sample, so it
    is an order-of-magnitude planning figure, not a power calculation. Returned
    rather than a bare "not estimable" because "collect about 1,200 more" is
    something a site can act on.
    """
    if baseline.estimable or not baseline.n_isolates:
        return None
    prevalence = baseline.n_resistant_truth / baseline.n_isolates
    if prevalence <= 0:
        return None
    answered_rate = baseline.n_answered / baseline.n_isolates
    if answered_rate <= 0:
        return None
    needed = (min_resistant / prevalence) / answered_rate
    return max(0, int(needed) - baseline.n_isolates)


CAVEAT = (
    "CRyPTIC and the WHO catalogue share underlying isolates, so this is the "
    "incumbent's home-ground performance: the bar a challenger must clear, not "
    "an independent accuracy estimate for the catalogue."
)
