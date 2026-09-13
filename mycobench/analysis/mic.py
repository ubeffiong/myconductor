"""MIC on the doubling-dilution scale, with censoring that stays censored.

Why this module is not three lines of ``log2()``
-----------------------------------------------
A microtitre plate tests a finite dilution series. A result of ``<=0.25`` does
not mean 0.25; it means the true value lies somewhere below the lowest tested
well. Reading it as 0.25 converts a bound into a measurement, and the bias is
not small — measured on the real CRyPTIC release:

    amikacin      69% left-censored
    delamanid     77% left-censored
    clofazimine   56% left-censored
    rifabutin     65% left-censored + 25% right-censored
    rifampicin    30% left-censored + 37% right-censored
    ethambutol     1% left-censored +  7% right-censored
    linezolid      1% left-censored +  1% right-censored

For the first four, **more than half the observations sit outside the tested
range, so the median is not identifiable.** A "median log2 MIC" for amikacin
computed by treating ``<=0.25`` as 0.25 is not a median of anything; it is the
lowest well of the plate, reported with false precision.

So this module does two things differently. It **refuses** to return a location
statistic when censoring makes it unidentifiable, naming the reason rather than
returning a number. And it provides a comparison that remains valid under heavy
censoring, by only counting a pairwise comparison when the censoring bounds
make its direction unambiguous. The effect is conditional on decidable pairs;
it is not a lower bound on the uncensored population effect.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from ..phenotypes import MIC, parse_mic
from . import stats

#: Censoring fraction above which a one-sided order statistic is unidentifiable.
#: At exactly 0.5 the median sits on the censoring bound, which is already
#: uninformative, so the check is inclusive.
CENSORING_LIMIT = 0.5

#: Minimum fraction of pairwise comparisons that must be decidable before a
#: stochastic-ordering effect size is reported. Below this the two groups are
#: mostly mutually censored and the estimate describes the plate, not the
#: biology.
MIN_DECISIVE_FRACTION = 0.30

NONE, LEFT, RIGHT = "none", "left", "right"


@dataclass(frozen=True)
class Observation:
    """One MIC on the log2 dilution scale, with its censoring direction kept.

    ``log2_bound`` is the log2 of the reported number. For an exact value that
    is the measurement; for a censored value it is the *bound*, and the true
    value lies strictly below (left) or above (right) it.
    """

    log2_bound: float
    censoring: str = NONE
    raw: str = ""

    @property
    def exact(self) -> bool:
        return self.censoring == NONE

    def describe(self) -> str:
        if self.censoring == LEFT:
            return f"<={self.log2_bound:+.1f} log2"
        if self.censoring == RIGHT:
            return f">{self.log2_bound:+.1f} log2"
        return f"{self.log2_bound:+.1f} log2"


def to_observation(value: MIC | str | None) -> Optional[Observation]:
    """Convert a parsed or raw MIC onto the log2 scale.

    Returns ``None`` for anything unusable — an absent value, a non-numeric
    one, or a non-positive one (log2 of zero is not a dilution).
    """
    if value is None:
        return None
    mic = parse_mic(value) if isinstance(value, str) else value
    if mic is None or mic.value is None or not math.isfinite(mic.value) or mic.value <= 0:
        return None
    return Observation(log2_bound=math.log2(mic.value),
                       censoring=mic.censoring, raw=mic.raw)


def observations(values: Iterable[MIC | str | None]) -> list[Observation]:
    return [obs for obs in (to_observation(v) for v in values) if obs is not None]


# -- censoring diagnostics ------------------------------------------------
@dataclass
class CensoringProfile:
    n: int
    n_exact: int
    n_left: int
    n_right: int

    @property
    def left_fraction(self) -> float:
        return self.n_left / self.n if self.n else 0.0

    @property
    def right_fraction(self) -> float:
        return self.n_right / self.n if self.n else 0.0

    @property
    def median_identifiable(self) -> bool:
        """Is a median defined, given how much of the sample is out of range?"""
        return (self.n > 0
                and self.left_fraction < CENSORING_LIMIT
                and self.right_fraction < CENSORING_LIMIT)

    def why_not_identifiable(self) -> Optional[str]:
        if self.n == 0:
            return "no usable MIC observations"
        if self.left_fraction >= CENSORING_LIMIT:
            return (f"{self.left_fraction:.0%} of observations are below the "
                    f"lowest tested dilution, so the median is the censoring "
                    f"bound rather than a measurement")
        if self.right_fraction >= CENSORING_LIMIT:
            return (f"{self.right_fraction:.0%} of observations are above the "
                    f"highest tested dilution, so the median is the censoring "
                    f"bound rather than a measurement")
        return None

    def describe(self) -> str:
        return (f"n={self.n} exact={self.n_exact} "
                f"left={self.left_fraction:.0%} right={self.right_fraction:.0%}")


def censoring_profile(obs: Sequence[Observation]) -> CensoringProfile:
    return CensoringProfile(
        n=len(obs),
        n_exact=sum(1 for o in obs if o.censoring == NONE),
        n_left=sum(1 for o in obs if o.censoring == LEFT),
        n_right=sum(1 for o in obs if o.censoring == RIGHT))


@dataclass
class Location:
    """A location statistic, or an explicit refusal to provide one."""

    value: Optional[float]
    interval: Optional[tuple[float, float]] = None
    profile: Optional[CensoringProfile] = None
    refused_because: Optional[str] = None

    @property
    def available(self) -> bool:
        return self.value is not None

    def describe(self) -> str:
        if not self.available:
            return f"not identifiable: {self.refused_because}"
        text = f"{self.value:+.2f} log2"
        if self.interval:
            text += f" (95% CI {self.interval[0]:+.2f} to {self.interval[1]:+.2f})"
        return text


def median_log2(obs: Sequence[Observation],
                with_interval: bool = True) -> Location:
    """Median log2 MIC, refused when censoring makes it undefined.

    When identifiable, censored observations are used at their bounds: with
    fewer than half the sample censored on either side the median falls inside
    the observed range, so a bound cannot move it past a neighbouring exact
    value. That is why the identifiability check comes first.
    """
    profile = censoring_profile(obs)
    reason = profile.why_not_identifiable()
    if reason:
        return Location(None, profile=profile, refused_because=reason)
    def bounds(sample):
        lower = stats.median([float("-inf") if o.censoring == LEFT else o.log2_bound for o in sample])
        upper = stats.median([float("inf") if o.censoring == RIGHT else o.log2_bound for o in sample])
        return lower, upper
    low, high = bounds(obs)
    if low != high or not math.isfinite(low):
        return Location(None, profile=profile,
                        refused_because="censoring intervals do not identify a unique median")
    interval = None
    if with_interval:
        import random
        rng = random.Random(stats.BOOTSTRAP_SEED)
        draws = [bounds([obs[rng.randrange(len(obs))] for _ in obs]) for _ in range(1000)]
        lower = sorted(x[0] for x in draws)[25]
        upper = sorted(x[1] for x in draws)[975]
        if math.isfinite(lower) and math.isfinite(upper):
            interval = (lower, upper)
    return Location(low, interval, profile)


# -- comparison that survives censoring ----------------------------------
def compare(a: Observation, b: Observation) -> Optional[int]:
    """Order two observations, or ``None`` when censoring makes it undecidable.

    The rules follow from what a bound actually asserts:

    * two exact values compare directly;
    * ``<=x`` versus an exact ``y > x`` is decidable (the censored one is
      smaller), because the true value is below ``x`` which is below ``y``;
    * ``<=x`` versus an exact ``y <= x`` is **not** decidable — the true value
      could be either side of ``y``;
    * two left-censored values are never decidable against each other, however
      different their bounds, because both true values lie in the same open
      region below the plate;
    * left versus right censored is decidable only when the bounds do not overlap.

    Only decidable pairs enter the denominator. The effect is conditional on
    those pairs; undecidable pairs do not dilute it toward zero.
    """
    if a.exact and b.exact:
        if a.log2_bound > b.log2_bound:
            return 1
        if a.log2_bound < b.log2_bound:
            return -1
        return 0

    if a.censoring == LEFT and b.censoring == LEFT:
        return None
    if a.censoring == RIGHT and b.censoring == RIGHT:
        return None
    if a.censoring == LEFT and b.censoring == RIGHT:
        return -1 if a.log2_bound <= b.log2_bound else None
    if a.censoring == RIGHT and b.censoring == LEFT:
        return 1 if b.log2_bound <= a.log2_bound else None

    # One censored, one exact.
    if a.censoring == LEFT:
        return -1 if b.log2_bound > a.log2_bound else None
    if a.censoring == RIGHT:
        return 1 if b.log2_bound <= a.log2_bound else None
    if b.censoring == LEFT:
        return 1 if a.log2_bound > b.log2_bound else None
    # b right-censored, a exact
    return -1 if a.log2_bound <= b.log2_bound else None


@dataclass
class StochasticShift:
    """A censoring-robust comparison of two MIC groups.

    ``delta`` is Cliff's delta over decidable pairs only: the probability a
    carrier exceeds a non-carrier minus the reverse. It is bounded in
    [-1, 1] and is conditional on the decidable subset. Excluding undecidable
    pairs does not make it a bound on the full population effect.
    """

    delta: Optional[float]
    n_a: int
    n_b: int
    n_pairs: int
    n_decisive: int
    median_shift: Optional[float] = None
    median_shift_interval: Optional[tuple[float, float]] = None
    refused_because: Optional[str] = None
    profile_a: Optional[CensoringProfile] = None
    profile_b: Optional[CensoringProfile] = None

    @property
    def available(self) -> bool:
        return self.delta is not None

    @property
    def decisive_fraction(self) -> float:
        return self.n_decisive / self.n_pairs if self.n_pairs else 0.0

    def describe(self) -> str:
        if not self.available:
            return f"not estimable: {self.refused_because}"
        text = (f"delta={self.delta:+.3f} "
                f"({self.n_decisive}/{self.n_pairs} pairs decidable)")
        if self.median_shift is not None:
            text += f", median shift {self.median_shift:+.2f} log2"
            if self.median_shift_interval:
                low, high = self.median_shift_interval
                text += f" (CI {low:+.2f} to {high:+.2f})"
        else:
            text += ", median shift not identifiable"
        return text


def stochastic_shift(group_a: Sequence[Observation],
                     group_b: Sequence[Observation],
                     min_per_group: int = stats.MIN_N_PER_GROUP,
                     min_decisive: float = MIN_DECISIVE_FRACTION
                     ) -> StochasticShift:
    """Compare two MIC groups without assuming the censoring away.

    ``group_a`` is conventionally the carriers of the variant under test and
    ``group_b`` the non-carriers, so a positive delta means carriers have
    higher MICs.
    """
    profile_a = censoring_profile(group_a)
    profile_b = censoring_profile(group_b)
    n_pairs = len(group_a) * len(group_b)
    base = StochasticShift(None, len(group_a), len(group_b), n_pairs, 0,
                           profile_a=profile_a, profile_b=profile_b)

    if len(group_a) < min_per_group or len(group_b) < min_per_group:
        base.refused_because = (
            f"{len(group_a)} vs {len(group_b)} observations; both groups need "
            f"at least {min_per_group}")
        return base

    greater = lesser = decisive = 0
    for a in group_a:
        for b in group_b:
            order = compare(a, b)
            if order is None:
                continue
            decisive += 1
            if order > 0:
                greater += 1
            elif order < 0:
                lesser += 1
    base.n_decisive = decisive

    if decisive == 0 or base.decisive_fraction < min_decisive:
        base.refused_because = (
            f"only {base.decisive_fraction:.0%} of pairs are decidable "
            f"({decisive}/{n_pairs}); the two groups are largely mutually "
            f"censored, so any effect size would describe the dilution range "
            f"rather than the isolates")
        return base

    base.delta = (greater - lesser) / decisive

    # The median shift is a secondary, more interpretable statistic — reported
    # only when both groups' medians are themselves identifiable.
    location_a = median_log2(group_a, with_interval=False)
    location_b = median_log2(group_b, with_interval=False)
    if location_a.available and location_b.available:
        base.median_shift = location_a.value - location_b.value
        base.median_shift_interval = stats.difference_ci(
            [o.log2_bound for o in group_a], [o.log2_bound for o in group_b],
            stats.median, min_n=min_per_group)
    return base


# -- per-drug panel ------------------------------------------------------
@dataclass
class DrugPanel:
    """Every usable MIC for one drug, with its censoring diagnostics."""

    drug: str
    by_isolate: dict[str, Observation] = field(default_factory=dict)

    def add(self, isolate: str, value: MIC | str | None) -> bool:
        obs = to_observation(value)
        if obs is None:
            return False
        self.by_isolate[isolate] = obs
        return True

    def subset(self, isolates: Iterable[str]) -> list[Observation]:
        return [self.by_isolate[i] for i in isolates if i in self.by_isolate]

    @property
    def profile(self) -> CensoringProfile:
        return censoring_profile(list(self.by_isolate.values()))

    def usable_for_location(self) -> bool:
        return self.profile.median_identifiable

    def describe(self) -> str:
        profile = self.profile
        note = "" if profile.median_identifiable else \
            f" — location statistics refused: {profile.why_not_identifiable()}"
        return f"{self.drug}: {profile.describe()}{note}"
