"""Background-conditional effect sizes, replacing marginal association.

What this computes, and why it is not a correlation
---------------------------------------------------
For a candidate variant and a drug, compare the MICs of carriers against
non-carriers **within each resistance background**, then pool those
within-background comparisons. Holding the background fixed is what removes the
hitchhiker: a variant that only ever travels with ``rpoB S450L`` has no
variation left to explain once that determinant is held constant, and its
pooled effect collapses toward zero.

Five outcome states, deliberately distinct
------------------------------------------
``evidence-of-effect``   a pooled shift whose interval excludes zero
``no-evidence``          comparable, and the shift is not distinguishable from zero
``confounded``           every carrier sits in a stratum with no non-carrier, so
                         the variant never varies at fixed background — the
                         hitchhiker signature, and a finding in its own right
``underpowered``         too few carriers, or too few informative strata
``not-estimable``        censoring leaves too few decidable pairs

Collapsing ``confounded`` and ``underpowered`` into ``no-evidence`` is the
mistake that lets a lineage marker look merely unproven rather than
structurally unevaluable.

Scans are multiple-testing corrected. A sweep over the 20,843 uncertain entries
in the real WHO catalogue will manufacture findings otherwise: Holm for a
confirmatory claim about one variant, Benjamini-Hochberg for a laboratory
follow-up queue, and both are reported because they answer different questions.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from . import mic, stats
from .strata import Isolate, Stratification, Stratum

#: A stratum needs this many on each side before it contributes.
MIN_PER_SIDE = 5
#: Fewer carriers than this across all informative strata and no claim is made.
MIN_CARRIERS = 5
#: Co-occurrence with an established determinant above this fraction is
#: reported as a confounding warning even when strata are informative.
COOCCURRENCE_WARN = 0.80


@dataclass
class StratumEffect:
    """The comparison inside one resistance background."""

    background: frozenset[str]
    label: str
    shift: mic.StochasticShift
    n_carriers: int
    n_non_carriers: int

    @property
    def usable(self) -> bool:
        return self.shift.available

    @property
    def weight(self) -> int:
        """Decidable pairs — the information the stratum actually contributes."""
        return self.shift.n_decisive if self.shift.available else 0


@dataclass
class VariantEffect:
    """Pooled, background-conditional effect of one variant on one drug."""

    variant: str
    drug: str
    verdict: str
    delta: Optional[float] = None
    interval: Optional[tuple[float, float]] = None
    pvalue: Optional[float] = None
    pvalue_holm: Optional[float] = None
    survives_fdr: Optional[bool] = None
    n_carriers: int = 0
    n_informative_strata: int = 0
    #: None when no carrier has a recorded lineage: unknown, not 1.0.
    lineage_diversity: Optional[float] = None
    strata: list[StratumEffect] = field(default_factory=list)
    cooccurrence: dict[str, float] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)

    @property
    def significant(self) -> bool:
        return self.verdict == "evidence-of-effect"

    @property
    def top_cooccurrence(self) -> Optional[tuple[str, float]]:
        if not self.cooccurrence:
            return None
        determinant, fraction = max(self.cooccurrence.items(),
                                    key=lambda kv: kv[1])
        return determinant, fraction

    def describe(self) -> str:
        head = f"{self.variant} / {self.drug}: {self.verdict}"
        if self.delta is not None:
            head += f" delta={self.delta:+.3f}"
            if self.interval:
                head += f" (CI {self.interval[0]:+.3f} to {self.interval[1]:+.3f})"
        if self.reasons:
            head += " — " + "; ".join(self.reasons)
        return head


def _pooled_delta(effects: Sequence[StratumEffect]) -> Optional[float]:
    """Weight each stratum by its decidable pairs.

    Weighting by information rather than by isolate count matters under heavy
    censoring: a large stratum in which almost nothing is decidable should not
    dominate a smaller one where the comparison is clean.
    """
    usable = [e for e in effects if e.usable and e.weight > 0]
    if not usable:
        return None
    total = sum(e.weight for e in usable)
    return sum(e.shift.delta * e.weight for e in usable) / total


#: Above this many carrier x non-carrier pairs in one stratum, the groups are
#: subsampled for the interval. The point estimate still uses every isolate;
#: only the resampling is bounded, because a bootstrap that never finishes
#: yields no interval at all — which is strictly worse than a slightly wider one.
MAX_BOOTSTRAP_PAIRS = 40_000
#: Total pairwise comparisons budgeted across the whole resampling.
BOOTSTRAP_PAIR_BUDGET = 8_000_000
MIN_BOOTSTRAP_ITERATIONS = 200
MAX_BOOTSTRAP_ITERATIONS = 1000


def _comparison_matrix(group_a: Sequence[mic.Observation],
                       group_b: Sequence[mic.Observation]
                       ) -> list[tuple[int, ...]]:
    """Pairwise order, computed once. ``None`` becomes 0 and is not counted.

    Building this up front is what makes the bootstrap affordable: the
    censoring rules are evaluated ``n_a x n_b`` times in total rather than once
    per resample, and each iteration then does integer lookups.
    """
    matrix = []
    for a in group_a:
        row = []
        for b in group_b:
            order = mic.compare(a, b)
            # 2 marks "undecidable" so it can be excluded from the denominator
            # without confusing it with a genuine tie at 0.
            row.append(2 if order is None else order)
        matrix.append(tuple(row))
    return matrix


def _delta_from_matrix(matrix: Sequence[Sequence[int]],
                       rows: Sequence[int],
                       cols: Sequence[int]) -> Optional[float]:
    greater = lesser = decisive = 0
    for i in rows:
        row = matrix[i]
        for j in cols:
            order = row[j]
            if order == 2:
                continue
            decisive += 1
            if order > 0:
                greater += 1
            elif order < 0:
                lesser += 1
    if decisive == 0:
        return None
    return (greater - lesser) / decisive


def _stratified_bootstrap(stratification: Stratification,
                          panel: mic.DrugPanel,
                          seed: int = stats.BOOTSTRAP_SEED
                          ) -> Optional[tuple[float, float]]:
    """Resample within strata, preserving the stratification.

    Resampling across strata would break exactly the conditioning the estimate
    depends on, so each stratum is resampled independently and the per-stratum
    deltas are recombined by their decidable-pair weight.

    The iteration count adapts to the size of the comparison, because the
    statistic is quadratic in group size: a fixed 1000 iterations over a
    200-by-200 stratum is 40 million censoring decisions per variant, which on
    a real cohort does not terminate in useful time.
    """
    informative = [s for s in stratification.informative_strata
                   if len(s.carriers) >= MIN_PER_SIDE
                   and len(s.non_carriers) >= MIN_PER_SIDE]
    if not informative:
        return None

    rng = random.Random(seed)
    prepared = []
    total_pairs = 0
    for stratum in informative:
        carriers = [i.isolate_id for i in stratum.carriers]
        non_carriers = [i.isolate_id for i in stratum.non_carriers]
        a_obs = panel.subset(carriers)
        b_obs = panel.subset(non_carriers)
        if len(a_obs) < MIN_PER_SIDE or len(b_obs) < MIN_PER_SIDE:
            continue
        # Bound the resampling cost without touching the point estimate.
        if len(a_obs) * len(b_obs) > MAX_BOOTSTRAP_PAIRS:
            cap = max(MIN_PER_SIDE, int(MAX_BOOTSTRAP_PAIRS ** 0.5))
            a_obs = rng.sample(a_obs, min(len(a_obs), cap))
            b_obs = rng.sample(b_obs, min(len(b_obs), cap))
        prepared.append((_comparison_matrix(a_obs, b_obs),
                         len(a_obs), len(b_obs)))
        total_pairs += len(a_obs) * len(b_obs)

    if not prepared or total_pairs == 0:
        return None

    iterations = max(MIN_BOOTSTRAP_ITERATIONS,
                     min(MAX_BOOTSTRAP_ITERATIONS,
                         BOOTSTRAP_PAIR_BUDGET // total_pairs))

    draws: list[float] = []
    for _ in range(iterations):
        weighted_sum = 0.0
        weight_total = 0
        for matrix, n_a, n_b in prepared:
            rows = [rng.randrange(n_a) for _ in range(n_a)]
            cols = [rng.randrange(n_b) for _ in range(n_b)]
            delta = _delta_from_matrix(matrix, rows, cols)
            if delta is None:
                continue
            weight = sum(matrix[i][j] != 2 for i in rows for j in cols)
            weighted_sum += delta * weight
            weight_total += weight
        if weight_total > 0:
            draws.append(weighted_sum / weight_total)

    if len(draws) < iterations // 4:
        return None
    draws.sort()
    return (draws[int(0.025 * len(draws))],
            draws[min(len(draws) - 1, int(0.975 * len(draws)))])


def _stratified_permutation(stratification, panel, observed, min_per_side=MIN_PER_SIDE):
    """Label-permutation p-value on the pooled effect size itself.

    The statistic permuted here is the same censoring-aware delta reported as
    the effect, not a difference in medians. That is not a stylistic choice.
    MIC values lie on a discrete doubling-dilution series, so a median takes
    only a handful of distinct values and a median-difference permutation test
    loses almost all of its resolution: on a perfectly separated 8-versus-8
    comparison it returns p near 0.6, because any split placing most of the
    high values on one side reproduces the observed difference exactly. Delta
    varies continuously with the split, so it discriminates.

    Labels are permuted only within a stratum, which is what keeps the null
    conditional on background, site and lineage rather than destroying the
    stratification the point estimate depends on.
    """
    prepared = []
    pairs = 0
    for stratum in stratification.informative_strata:
        a = panel.subset(i.isolate_id for i in stratum.carriers)
        b = panel.subset(i.isolate_id for i in stratum.non_carriers)
        if min(len(a), len(b)) < min_per_side:
            continue
        pooled = a + b
        pairs += len(a) * len(b)
        if pairs > 100000 or len(pooled) ** 2 > 1000000:
            return None
        # Matrix preserves censoring when the labels change. The statistic
        # cannot be computed by treating censored bounds as exact MICs.
        prepared.append((_comparison_matrix(pooled, pooled), len(a), len(b)))
    if not prepared or observed is None:
        return None
    # Refuse excessive work rather than silently subsampling the hypothesis.
    if pairs > 100000:
        return None
    iterations = max(199, min(1999, 2000000 // max(1, pairs)))
    rng = random.Random(stats.PERMUTATION_SEED)
    extreme = valid = 0
    for _ in range(iterations):
        numerator = denominator = 0
        for matrix, na, nb in prepared:
            indices = list(range(na + nb))
            rng.shuffle(indices)
            for i in indices[:na]:
                for j in indices[na:]:
                    order = matrix[i][j]
                    if order != 2:
                        numerator += order
                        denominator += 1
        if denominator:
            valid += 1
            extreme += abs(numerator / denominator) >= abs(observed) - 1e-12
    return (extreme + 1) / (valid + 1) if valid >= iterations // 2 else None


def variant_effect(stratification: Stratification, panel: mic.DrugPanel,
                   cooccurrence: Optional[dict[str, float]] = None,
                   min_carriers: int = MIN_CARRIERS,
                   min_per_side: int = MIN_PER_SIDE) -> VariantEffect:
    """Estimate one variant's background-conditional effect on one drug."""
    effect = VariantEffect(variant=stratification.variant,
                           drug=stratification.drug,
                           verdict="underpowered",
                           n_carriers=stratification.n_carriers,
                           cooccurrence=dict(cooccurrence or {}))

    if stratification.fully_confounded:
        effect.verdict = "confounded"
        top = effect.top_cooccurrence
        detail = (f"; always co-occurs with {top[0]} in {top[1]:.0%} of "
                  f"carriers" if top else "")
        effect.reasons.append(
            f"all {stratification.n_carriers} carrier(s) sit in strata with no "
            f"non-carrier, so the variant never varies at fixed background"
            + detail)
        return effect

    if stratification.n_carriers < min_carriers:
        effect.reasons.append(
            f"{stratification.n_carriers} carrier(s), below the "
            f"{min_carriers} minimum")
        return effect

    for stratum in stratification.informative_strata:
        if (len(stratum.carriers) < min_per_side
                or len(stratum.non_carriers) < min_per_side):
            continue
        shift = mic.stochastic_shift(
            panel.subset(i.isolate_id for i in stratum.carriers),
            panel.subset(i.isolate_id for i in stratum.non_carriers),
            min_per_group=min_per_side)
        effect.strata.append(StratumEffect(
            background=stratum.background, label=stratum.label(), shift=shift,
            n_carriers=len(stratum.carriers),
            n_non_carriers=len(stratum.non_carriers)))

    usable = [s for s in effect.strata if s.usable]
    effect.n_informative_strata = len(usable)
    effect.lineage_diversity = stratification.lineage_diversity()

    if not usable:
        effect.verdict = "not-estimable"
        if effect.strata:
            effect.reasons.append(
                "no stratum yielded a decidable comparison: "
                + (effect.strata[0].shift.refused_because or "censoring"))
        else:
            effect.reasons.append(
                f"no stratum had at least {min_per_side} isolates on both sides")
        return effect

    effect.delta = _pooled_delta(usable)
    effect.interval = _stratified_bootstrap(stratification, panel)

    # Test the same censor-aware pooled ordering statistic as the effect,
    # permuting labels only within the same background/site/lineage strata.
    # A median test on one stratum answers a different question and can have
    # very low power even when every carrier exceeds every non-carrier.
    effect.pvalue = _stratified_permutation(stratification, panel, effect.delta, min_per_side)

    top = effect.top_cooccurrence
    if top and top[1] >= COOCCURRENCE_WARN:
        effect.warnings.append(
            f"co-occurs with the established determinant {top[0]} in "
            f"{top[1]:.0%} of carriers; the conditional estimate rests on the "
            f"minority of carriers that lack it")
    if effect.lineage_diversity is None:
        # Not a diversity of 1.0. Saying "these carriers span one lineage"
        # about isolates that were never typed would report an absence of
        # evidence as evidence of restriction — the exact substitution this
        # codebase exists to prevent.
        effect.warnings.append(
            "lineage is not recorded for any carrier, so lineage confounding "
            "was NOT assessed; this is unknown, not absent")
    elif effect.lineage_diversity < 1.5 and effect.n_carriers >= min_carriers:
        effect.warnings.append(
            f"carriers span an effective {effect.lineage_diversity:.1f} "
            f"lineage(s); the effect may be lineage-specific rather than causal")

    if effect.pvalue is None:
        effect.verdict = "underpowered"
        effect.reasons.append("no valid pooled permutation test; effect is descriptive only")
    elif effect.interval is not None and (effect.interval[0] > 0
                                        or effect.interval[1] < 0):
        effect.verdict = "evidence-of-effect"
        effect.reasons.append(
            f"pooled delta {effect.delta:+.3f} with an interval excluding zero "
            f"across {effect.n_informative_strata} stratum/strata")
    elif effect.interval is None:
        effect.verdict = "underpowered"
        effect.reasons.append(
            "pooled estimate available but no stable interval; too few "
            "informative strata to resample")
    else:
        effect.verdict = "no-evidence"
        effect.reasons.append(
            f"pooled delta {effect.delta:+.3f} with an interval spanning zero")
    return effect


@dataclass
class EffectScan:
    """A multiple-testing-corrected sweep over many variant-drug pairs."""

    effects: list[VariantEffect] = field(default_factory=list)
    fdr_alpha: float = 0.05

    def finalise(self) -> "EffectScan":
        """Apply Holm and Benjamini-Hochberg across the whole scan."""
        pvalues = [e.pvalue for e in self.effects]
        for effect, adjusted in zip(self.effects, stats.holm_adjust(pvalues)):
            effect.pvalue_holm = adjusted
        for effect, survives in zip(self.effects,
                                    stats.benjamini_hochberg(pvalues,
                                                             self.fdr_alpha)):
            effect.survives_fdr = survives
        # A nominally significant pooled interval that does not survive
        # family-wise correction is demoted, so a scan cannot launder a single
        # uncorrected finding into a claim.
        for effect in self.effects:
            if (effect.verdict == "evidence-of-effect"
                    and effect.pvalue_holm is not None
                    and effect.pvalue_holm > self.fdr_alpha):
                effect.verdict = "no-evidence"
                effect.reasons.append(
                    f"interval excluded zero but Holm-adjusted p is "
                    f"{effect.pvalue_holm:.3f} across {len(self.effects)} "
                    f"tested pair(s); demoted")
        return self

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for effect in self.effects:
            out[effect.verdict] = out.get(effect.verdict, 0) + 1
        return out

    def validation_queue(self, limit: int = 50) -> list[VariantEffect]:
        """Variants worth laboratory follow-up, most informative first.

        Ranked on FDR survival then effect magnitude — a discovery queue wants
        false-discovery control, not the family-wise control a single
        confirmatory claim needs.
        """
        candidates = [e for e in self.effects
                      if e.verdict in ("evidence-of-effect", "no-evidence")
                      and e.delta is not None]
        candidates.sort(key=lambda e: (not bool(e.survives_fdr),
                                       -abs(e.delta)))
        return candidates[:limit]

    def confounded(self) -> list[VariantEffect]:
        """Variants structurally unevaluable — the hitchhiker list."""
        return [e for e in self.effects if e.verdict == "confounded"]

    def describe(self) -> str:
        counts = self.counts()
        lines = [f"{len(self.effects)} variant-drug pair(s) tested"]
        for verdict, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {n:>6}  {verdict}")
        return "\n".join(lines)


def scan(stratifications: Iterable[tuple[Stratification, mic.DrugPanel,
                                         dict[str, float]]],
         fdr_alpha: float = 0.05) -> EffectScan:
    """Estimate effects for many pairs and correct across the family."""
    result = EffectScan(fdr_alpha=fdr_alpha)
    for stratification, panel, cooccur in stratifications:
        result.effects.append(
            variant_effect(stratification, panel, cooccur))
    return result.finalise()
