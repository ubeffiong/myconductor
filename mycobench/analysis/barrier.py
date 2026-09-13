"""Target ranking by resistance barrier, rather than by essentiality.

Why the usual pipeline returns the same answer every time
---------------------------------------------------------
Subtractive genomics — drop human homologs, keep essential genes, keep
druggable ones — is run repeatedly and converges on roughly the same few
hundred proteins, because those three filters are properties of the reference
genome and the reference genome does not change. None of the three says
anything about how quickly resistance will arrive, which is the property that
decided the fate of every TB drug so far.

The alternative criterion
-------------------------
Rank candidate targets by how *hard* it is for the organism to escape them.
A target where escape requires a rare, fitness-costly mutation is worth more
than one that is merely essential. That is estimable from surveillance data
already being collected:

* **escape route count** — how many distinct variants in that gene have been
  observed conferring resistance. Many independent routes means a low barrier.
* **lineage breadth** — whether escape variants arise across the phylogeny
  (tolerated, low barrier) or only in one lineage (possibly compensated or
  lineage-specific).
* **truncation tolerance** — whether loss of function is a viable escape.
  A gene the organism can simply switch off is a poor target; one where escape
  needs a specific substitution that preserves function is a better one.
* **observed frequency** — how common escape is among treated isolates.

The falsification test is built in
----------------------------------
``validate_against_known`` requires the criterion to reproduce a contrast the
bedaquiline literature already settles: ``mmpR5``/``Rv0678`` must rank as a
**low** barrier, because truncating it is common and tolerated and accounts for
most bedaquiline resistance; ``atpE`` must rank as a **high** barrier, because
escape there is rare. A criterion that cannot recover that ordering is not
measuring resistance barrier and should not be used to choose targets.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from . import stats

#: Consequence classes that constitute loss of function.
TRUNCATING = ("frameshift", "nonsense", "deletion", "insertion",
              "feature_ablation", "stop_gained", "frameshift_variant")

#: Fewer observed escape events than this and the barrier is not estimated:
#: absence of observed escape in a small sample is not evidence of a high
#: barrier, it is absence of data.
MIN_ESCAPE_EVENTS = 5
#: Isolates screened, below which nothing is claimed either way.
MIN_SCREENED = 200


@dataclass
class EscapeVariant:
    """One observed resistance variant in a candidate target gene."""

    variant: str
    consequence: str = ""
    n_isolates: int = 0
    lineages: tuple[str, ...] = ()

    @property
    def truncating(self) -> bool:
        value = (self.consequence or "").strip().lower()
        return any(token in value for token in TRUNCATING)


@dataclass
class TargetBarrier:
    """A candidate target, scored by how easily the organism escapes it."""

    gene: str
    drug: str
    escape_variants: list[EscapeVariant] = field(default_factory=list)
    n_screened: int = 0
    essential: Optional[bool] = None
    human_homolog: Optional[bool] = None
    score: Optional[float] = None
    band: str = "not-estimated"
    reasons: list[str] = field(default_factory=list)

    # -- observable quantities ------------------------------------------
    @property
    def n_routes(self) -> int:
        """Distinct variants observed conferring resistance."""
        return len(self.escape_variants)

    @property
    def n_escape_isolates(self) -> int:
        return sum(v.n_isolates for v in self.escape_variants)

    @property
    def escape_frequency(self) -> Optional[float]:
        if not self.n_screened:
            return None
        return self.n_escape_isolates / self.n_screened

    @property
    def truncating_fraction(self) -> Optional[float]:
        if not self.escape_variants:
            return None
        truncating = sum(v.n_isolates for v in self.escape_variants
                         if v.truncating)
        return (truncating / self.n_escape_isolates
                if self.n_escape_isolates else None)

    @property
    def lineage_breadth(self) -> float:
        counts: dict[str, int] = {}
        for variant in self.escape_variants:
            for lineage in variant.lineages:
                counts[lineage] = counts.get(lineage, 0) + 1
        return stats.inverse_simpson(list(counts.values()))

    def describe(self) -> str:
        if self.score is None:
            return f"{self.gene}/{self.drug}: {self.band} — " + \
                   "; ".join(self.reasons)
        return (f"{self.gene}/{self.drug}: {self.band} barrier "
                f"(score {self.score:.3f}) — {self.n_routes} escape route(s), "
                f"frequency {self.escape_frequency:.3%}, "
                f"truncating {self.truncating_fraction:.0%}, "
                f"lineage breadth {self.lineage_breadth:.1f}")


def score_barrier(target: TargetBarrier,
                  min_escape: int = MIN_ESCAPE_EVENTS,
                  min_screened: int = MIN_SCREENED) -> TargetBarrier:
    """Score one target. Higher score means a higher barrier.

    The score is a deliberately simple, unweighted combination of four
    observable quantities. A weighted combination would imply a calibration
    against drug-development outcomes that nobody has, and this is a ranking
    criterion rather than a probability.
    """
    if target.n_screened < min_screened:
        target.band = "not-estimated"
        target.reasons.append(
            f"{target.n_screened} isolate(s) screened, below the "
            f"{min_screened} minimum; absence of observed escape in a small "
            f"sample is not evidence of a high barrier")
        return target

    if target.n_escape_isolates < min_escape:
        target.band = "no-observed-escape"
        target.reasons.append(
            f"{target.n_escape_isolates} escape event(s) across "
            f"{target.n_screened} isolate(s). Consistent with a high barrier, "
            f"but equally with a drug nobody has used long enough — this is "
            f"reported as a distinct state rather than scored as high")
        return target

    # Each component is mapped so that 1.0 means "hard to escape".
    frequency = target.escape_frequency or 0.0
    route_penalty = min(1.0, target.n_routes / 20.0)
    breadth_penalty = min(1.0, max(0.0, (target.lineage_breadth - 1.0) / 3.0))
    truncation_penalty = target.truncating_fraction or 0.0
    frequency_penalty = min(1.0, frequency / 0.10)

    components = {
        "few escape routes": 1.0 - route_penalty,
        "lineage-restricted escape": 1.0 - breadth_penalty,
        "loss of function not viable": 1.0 - truncation_penalty,
        "escape is rare": 1.0 - frequency_penalty,
    }
    target.score = sum(components.values()) / len(components)
    target.band = ("high" if target.score >= 0.66
                   else "moderate" if target.score >= 0.40 else "low")

    weakest = min(components.items(), key=lambda kv: kv[1])
    target.reasons.append(
        f"lowest component: {weakest[0]} at {weakest[1]:.2f}")
    if target.truncating_fraction and target.truncating_fraction >= 0.5:
        target.reasons.append(
            f"{target.truncating_fraction:.0%} of escape is loss of function; "
            f"a gene the organism can switch off is a poor target regardless "
            f"of how essential it appears")
    return target


@dataclass
class BarrierRanking:
    targets: list[TargetBarrier] = field(default_factory=list)

    def ranked(self) -> list[TargetBarrier]:
        """Highest barrier first; unscored states sort to the end."""
        scored = [t for t in self.targets if t.score is not None]
        unscored = [t for t in self.targets if t.score is None]
        scored.sort(key=lambda t: -t.score)
        return scored + unscored

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for target in self.targets:
            out[target.band] = out.get(target.band, 0) + 1
        return out

    def describe(self) -> str:
        lines = [f"{len(self.targets)} candidate target(s)"]
        for band, n in sorted(self.counts().items(), key=lambda kv: -kv[1]):
            lines.append(f"  {n:>4}  {band}")
        return "\n".join(lines)


def rank(targets: Iterable[TargetBarrier]) -> BarrierRanking:
    return BarrierRanking([score_barrier(t) for t in targets])


# -- the experiment that would kill this ---------------------------------
@dataclass
class BarrierValidation:
    """Does the criterion reproduce a contrast the literature already settles?"""

    low_barrier_gene: str
    high_barrier_gene: str
    low_band: Optional[str] = None
    high_band: Optional[str] = None
    low_score: Optional[float] = None
    high_score: Optional[float] = None
    verdict: str = "not-run"
    reasons: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.verdict == "reproduces-known-contrast"


def validate_against_known(ranking: BarrierRanking,
                           low_barrier_gene: str = "Rv0678",
                           high_barrier_gene: str = "atpE",
                           aliases: Optional[dict[str, str]] = None
                           ) -> BarrierValidation:
    """Require the ranking to place ``mmpR5`` below ``atpE`` for bedaquiline.

    This contrast is not a guess. Truncation of ``mmpR5``/``Rv0678`` accounts
    for most bedaquiline resistance and is clearly tolerated, while escape via
    ``atpE`` is rare — so any criterion claiming to measure resistance barrier
    must rank the first low and the second high. Failing this is decisive.
    """
    aliases = aliases or {"mmpR5": "Rv0678", "Rv0678": "Rv0678"}
    result = BarrierValidation(low_barrier_gene=low_barrier_gene,
                               high_barrier_gene=high_barrier_gene)

    def find(gene: str) -> Optional[TargetBarrier]:
        canonical = aliases.get(gene, gene)
        for target in ranking.targets:
            if aliases.get(target.gene, target.gene) == canonical:
                return target
        return None

    low = find(low_barrier_gene)
    high = find(high_barrier_gene)
    if low is None or high is None:
        result.reasons.append(
            f"both {low_barrier_gene} and {high_barrier_gene} must be in the "
            f"ranking to run this check; found "
            f"{'both' if low and high else low_barrier_gene if low else high_barrier_gene if high else 'neither'}")
        return result

    result.low_band, result.high_band = low.band, high.band
    result.low_score, result.high_score = low.score, high.score

    if low.score is None or high.score is None:
        result.verdict = "underpowered"
        result.reasons.append(
            f"{low_barrier_gene} is '{low.band}' and {high_barrier_gene} is "
            f"'{high.band}'; at least one was not scored, so the contrast "
            f"cannot be evaluated")
        return result

    if low.score < high.score:
        result.verdict = "reproduces-known-contrast"
        result.reasons.append(
            f"{low_barrier_gene} scored {low.score:.3f} ({low.band}) below "
            f"{high_barrier_gene} at {high.score:.3f} ({high.band}), matching "
            f"the published contrast")
    else:
        result.verdict = "fails-known-contrast"
        result.reasons.append(
            f"{low_barrier_gene} scored {low.score:.3f} at or above "
            f"{high_barrier_gene} at {high.score:.3f}. The criterion does not "
            f"reproduce a contrast the literature settles, so it is not "
            f"measuring resistance barrier and must not be used to select "
            f"targets")
    return result
