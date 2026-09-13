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
    lineage_diversity: float = 0.0
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


def _stratified_bootstrap(stratification: Stratification,
                          panel: mic.DrugPanel,
                          iterations: int = 1000,
                          seed: int = stats.BOOTSTRAP_SEED
                          ) -> Optional[tuple[float, float]]:
    """Resample within strata, preserving the stratification.

    Resampling across strata would break exactly the conditioning the estimate
    depends on.
    """
    informative = [s for s in stratification.informative_strata
                   if len(s.carriers) >= MIN_PER_SIDE
                   and len(s.non_carriers) >= MIN_PER_SIDE]
    if not informative:
        return None

    rng = random.Random(seed)
    draws: list[float] = []
    for _ in range(iterations):
        weighted_sum = 0.0
        weight_total = 0
        for stratum in informative:
            carriers = [i.isolate_id for i in
                        rng.choices(stratum.carriers, k=len(stratum.carriers))]
            non_carriers = [i.isolate_id for i in
                            rng.choices(stratum.non_carriers,
                                        k=len(stratum.non_carriers))]
            shift = mic.stochastic_shift(panel.subset(carriers),
                                         panel.subset(non_carriers),
                                         min_per_group=MIN_PER_SIDE)
            if shift.available and shift.n_decisive > 0:
                weighted_sum += shift.delta * shift.n_decisive
                weight_total += shift.n_decisive
        if weight_total > 0:
            draws.append(weighted_sum / weight_total)

    if len(draws) < iterations // 4:
        return None
    draws.sort()
    return (draws[int(0.025 * len(draws))],
            draws[min(len(draws) - 1, int(0.975 * len(draws)))])


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

    # A permutation test on the largest informative stratum. Pooling p-values
    # across strata would need a combination procedure whose assumptions are
    # harder to defend than simply reporting the best-powered stratum.
    largest = max(usable, key=lambda s: s.weight)
    stratum = next((s for s in stratification.informative_strata
                    if s.background == largest.background), None)
    if stratum is not None:
        carrier_values = [o.log2_bound for o in
                          panel.subset(i.isolate_id for i in stratum.carriers)]
        other_values = [o.log2_bound for o in
                        panel.subset(i.isolate_id for i in stratum.non_carriers)]
        effect.pvalue = stats.permutation_pvalue(carrier_values, other_values,
                                                 min_n=min_per_side)

    top = effect.top_cooccurrence
    if top and top[1] >= COOCCURRENCE_WARN:
        effect.warnings.append(
            f"co-occurs with the established determinant {top[0]} in "
            f"{top[1]:.0%} of carriers; the conditional estimate rests on the "
            f"minority of carriers that lack it")
    if effect.lineage_diversity < 1.5 and effect.n_carriers >= min_carriers:
        effect.warnings.append(
            f"carriers span an effective {effect.lineage_diversity:.1f} "
            f"lineage(s); the effect may be lineage-specific rather than causal")

    if effect.interval is not None and (effect.interval[0] > 0
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
