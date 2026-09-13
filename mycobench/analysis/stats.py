"""Statistical primitives. Standard library only, deterministic by construction.

Two rules run through this module.

**An interval from too few observations is a lie about precision.** Every
interval function refuses below a minimum n and returns ``None`` rather than a
number that would be plotted as though it meant something.

**Resampling is seeded.** A confidence interval that moves between two runs of
the same analysis is not reportable. The seeds are fixed constants, not clock
reads, so a figure in a report can be reproduced exactly.
"""
from __future__ import annotations

import math
import random
import statistics
from typing import Callable, Optional, Sequence

#: Fixed so an interval is reproducible across runs and machines.
BOOTSTRAP_SEED = 20260913
PERMUTATION_SEED = 20260914

#: Below this, a resampled interval describes the resampling, not the data.
MIN_N_FOR_INTERVAL = 8
#: Below this, a two-group comparison is not attempted at all.
MIN_N_PER_GROUP = 5


def mean(values: Sequence[float]) -> Optional[float]:
    return statistics.fmean(values) if values else None


def median(values: Sequence[float]) -> Optional[float]:
    return statistics.median(values) if values else None


def quantile(values: Sequence[float], q: float) -> Optional[float]:
    """Type-7 quantile, defined for a single observation."""
    if not values:
        return None
    if not 0.0 <= q <= 1.0:
        raise ValueError(f"q must be in [0, 1], got {q}")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = q * (len(ordered) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def bootstrap_ci(values: Sequence[float],
                 statistic: Callable[[Sequence[float]], Optional[float]] = median,
                 iterations: int = 2000,
                 level: float = 0.95,
                 min_n: int = MIN_N_FOR_INTERVAL,
                 seed: int = BOOTSTRAP_SEED
                 ) -> Optional[tuple[float, float]]:
    """Percentile bootstrap interval, or ``None`` when n is too small.

    Descriptive, not a hypothesis test: it says how much the statistic moves
    under resampling, and nothing about a null.
    """
    values = list(values)
    if len(values) < min_n:
        return None
    rng = random.Random(seed)
    draws = []
    for _ in range(iterations):
        sample = rng.choices(values, k=len(values))
        value = statistic(sample)
        if value is not None:
            draws.append(value)
    if len(draws) < iterations // 2:
        return None
    draws.sort()
    tail = (1.0 - level) / 2.0
    low = draws[int(tail * len(draws))]
    high = draws[min(len(draws) - 1, int((1.0 - tail) * len(draws)))]
    return low, high


def difference_ci(group_a: Sequence[float], group_b: Sequence[float],
                  statistic: Callable[[Sequence[float]], Optional[float]] = median,
                  iterations: int = 2000, level: float = 0.95,
                  min_n: int = MIN_N_PER_GROUP,
                  seed: int = BOOTSTRAP_SEED
                  ) -> Optional[tuple[float, float]]:
    """Bootstrap interval for ``statistic(a) - statistic(b)``.

    Both groups are resampled independently, which is the right model when the
    two are different isolates rather than paired measurements.
    """
    group_a, group_b = list(group_a), list(group_b)
    if len(group_a) < min_n or len(group_b) < min_n:
        return None
    rng = random.Random(seed)
    draws = []
    for _ in range(iterations):
        a = statistic(rng.choices(group_a, k=len(group_a)))
        b = statistic(rng.choices(group_b, k=len(group_b)))
        if a is not None and b is not None:
            draws.append(a - b)
    if len(draws) < iterations // 2:
        return None
    draws.sort()
    tail = (1.0 - level) / 2.0
    return (draws[int(tail * len(draws))],
            draws[min(len(draws) - 1, int((1.0 - tail) * len(draws)))])


def permutation_pvalue(group_a: Sequence[float], group_b: Sequence[float],
                       statistic: Callable[[Sequence[float]], Optional[float]] = median,
                       iterations: int = 10000,
                       min_n: int = MIN_N_PER_GROUP,
                       seed: int = PERMUTATION_SEED) -> Optional[float]:
    """Two-sided label-permutation p-value, or ``None`` when underpowered.

    Deliberately a permutation test rather than a parametric one: log-MIC
    distributions are censored and multimodal, so a t-test's assumptions do not
    hold here.

    **Do not use the default median statistic on MIC data.** Values lie on a
    discrete doubling-dilution series, so a median takes only a few distinct
    values and the test loses nearly all resolution: on a perfectly separated
    eight-versus-eight comparison it returns p near 0.6, because any split
    placing most high values on one side reproduces the observed difference
    exactly. ``analysis.effects._permutation_delta_pvalue`` permutes the
    censoring-aware effect size instead, which varies continuously with the
    split. This function remains for continuous statistics.
    """
    group_a, group_b = list(group_a), list(group_b)
    if len(group_a) < min_n or len(group_b) < min_n:
        return None
    observed_a, observed_b = statistic(group_a), statistic(group_b)
    if observed_a is None or observed_b is None:
        return None
    observed = abs(observed_a - observed_b)

    pooled = group_a + group_b
    cut = len(group_a)
    rng = random.Random(seed)
    extreme = 0
    for _ in range(iterations):
        rng.shuffle(pooled)
        a, b = statistic(pooled[:cut]), statistic(pooled[cut:])
        if a is None or b is None:
            continue
        extreme += abs(a - b) >= observed
    # Add-one correction: a p-value of exactly zero is not a finding, it is a
    # statement that the resampling did not go far enough.
    return (extreme + 1) / (iterations + 1)


def holm_adjust(pvalues: Sequence[Optional[float]]) -> list[Optional[float]]:
    """Holm-Bonferroni adjustment, preserving input order and ``None`` gaps.

    Needed because these analyses test many variant-drug pairs at once, and an
    unadjusted scan over thousands of pairs manufactures findings.
    """
    indexed = sorted(((p, i) for i, p in enumerate(pvalues) if p is not None),
                     key=lambda item: item[0])
    adjusted: list[Optional[float]] = [None] * len(pvalues)
    running = 0.0
    total = len(indexed)
    for rank, (pvalue, index) in enumerate(indexed):
        value = min(1.0, max(running, pvalue * (total - rank)))
        adjusted[index] = value
        running = value
    return adjusted


def benjamini_hochberg(pvalues: Sequence[Optional[float]],
                       alpha: float = 0.05) -> list[bool]:
    """Which hypotheses survive at a given false-discovery rate.

    Offered alongside Holm because these scans have two different purposes: a
    confirmatory claim about one variant wants Holm's family-wise control, and
    a discovery queue for laboratory follow-up wants FDR.
    """
    indexed = sorted(((p, i) for i, p in enumerate(pvalues) if p is not None),
                     key=lambda item: item[0])
    total = len(indexed)
    survives = [False] * len(pvalues)
    threshold_rank = -1
    for rank, (pvalue, _) in enumerate(indexed, start=1):
        if pvalue <= alpha * rank / total:
            threshold_rank = rank
    for rank, (_, index) in enumerate(indexed, start=1):
        if rank <= threshold_rank:
            survives[index] = True
    return survives


def cliffs_delta(group_a: Sequence[float],
                 group_b: Sequence[float]) -> Optional[float]:
    """Non-parametric effect size in [-1, 1]; robust to censoring.

    Reported beside a difference in medians because a median shift of one
    doubling dilution means something very different when the distributions
    barely overlap than when they almost coincide.
    """
    group_a, group_b = list(group_a), list(group_b)
    if not group_a or not group_b:
        return None
    greater = sum(1 for a in group_a for b in group_b if a > b)
    lesser = sum(1 for a in group_a for b in group_b if a < b)
    return (greater - lesser) / (len(group_a) * len(group_b))


def inverse_simpson(counts: Sequence[int]) -> float:
    """Effective number of groups. Used for study and lineage clustering."""
    total = sum(counts)
    if total <= 0:
        return 0.0
    return 1.0 / sum((count / total) ** 2 for count in counts if count > 0)
