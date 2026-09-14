"""Calibrating model predictions — to raise a barrier, never to lower one.

A model's raw probability is not a probability of anything until it has been
checked against outcomes. Calibration does that check, and this module performs
it against site-local phenotypic results. What it is allowed to *do* with the
answer is deliberately one-directional.

Calibration may raise priority. It may not grant a call.
--------------------------------------------------------
A calibrated high probability of resistance is a good reason to move a variant
up the laboratory validation queue: the lab confirms it, the confirmation
enters as PHENOTYPIC evidence, and that can establish resistance. Calibration
earns its keep there.

What it must never do is license a SUSCEPTIBLE call. The tempting argument —
"a low calibrated probability asserts the *absence* of resistance, which is the
safe direction" — has the asymmetry backwards in this system. ``SUSCEPTIBLE``
is the only call that admits a drug to a regimen (see ``modules/synthesis``),
so it is the *actionable* one. A wrong RESISTANT costs a usable drug; a wrong
SUSCEPTIBLE puts a patient on a failing regimen, which is how amplified
resistance and onward transmission of a resistant strain are manufactured.
``INDETERMINATE`` already says "unknown" and withholds the drug safely.

There is a structural reason too. ``SUSCEPTIBLE`` requires positive evidence
that the drug's loci were *callable*, drawn from outside the variant list. A
model probability is not coverage evidence: if a locus was never callable, the
model has no information about it either — it is scoring features of a variant
nobody could observe. Promotion on a probability would route around the
callable mask, which is the gate the whole design rests on.

Why a count of positives, not a count of observations
-----------------------------------------------------
The usual minimum — "calibrate on at least 30 paired observations" — fails
silently on exactly the drugs that matter most here. Resistance to bedaquiline
appears in roughly 0.8% of the CRyPTIC compendium, so thirty pairs contain a
quarter of one resistant isolate on average. No curve can be fitted to that,
and one fitted anyway would be describing noise. The binding constraint is the
number of *resistant* observations, not the total, so that is what is required
— and a site is told how many more it needs rather than handed a number that
looks calibrated.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

#: Paired observations required before calibration is attempted at all.
MIN_PAIRED = 30
#: Resistant observations required. This, not MIN_PAIRED, is what usually
#: binds: at 0.8% prevalence, 30 pairs hold a quarter of one resistant isolate.
MIN_POSITIVES = 10
#: Susceptible observations required, for the same reason in the other
#: direction: a set that is all-resistant cannot show a false-positive rate.
MIN_NEGATIVES = 10
#: Bins used for the reliability estimate. Few enough to be populated at the
#: sample sizes a single site can reach.
CALIBRATION_BINS = 5
#: Above this expected calibration error the model's probabilities are not
#: usable for ordering a queue, never mind anything else.
MAX_USABLE_ECE = 0.15
#: A calibrated probability at or above this is worth laboratory time.
PRIORITY_PROBABILITY = 0.50


@dataclass(frozen=True)
class PairedObservation:
    """One model probability set against the phenotype that followed."""

    variant_key: str
    drug: str
    probability: float
    resistant: bool
    site_id: str
    lineage: Optional[str] = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.probability) or not 0 <= self.probability <= 1:
            raise ValueError("paired probability must be in [0,1]")
        for name in ("variant_key", "drug", "site_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"paired observation requires {name}")


@dataclass
class CalibrationResult:
    """What the pairing supports, and what it explicitly does not."""

    drug: str
    site_id: str
    n_paired: int = 0
    n_resistant: int = 0
    n_susceptible: int = 0
    ece: Optional[float] = None
    usable_for_priority: bool = False
    refusals: list[str] = field(default_factory=list)
    lineages: dict[str, int] = field(default_factory=dict)

    #: Fixed, and not configurable through this result. Stated as data so a
    #: report can show the reader that no path to SUSCEPTIBLE exists here.
    may_establish_resistance: bool = False
    may_grant_susceptibility: bool = False

    @property
    def calibrated(self) -> bool:
        return self.ece is not None and not self.refusals

    def describe(self) -> str:
        if self.refusals:
            return (f"{self.drug} at {self.site_id}: not calibrated — "
                    + "; ".join(self.refusals))
        return (f"{self.drug} at {self.site_id}: calibrated on "
                f"{self.n_paired} pair(s) ({self.n_resistant} resistant), "
                f"expected calibration error {self.ece:.3f}"
                + ("" if self.usable_for_priority else
                   f"; above {MAX_USABLE_ECE} so probabilities are not used "
                   f"for ordering"))


def expected_calibration_error(observations: Sequence[PairedObservation],
                               bins: int = CALIBRATION_BINS) -> Optional[float]:
    """Mean gap between predicted probability and observed frequency.

    Bin-count weighted, over populated bins only. An empty bin contributes
    nothing rather than a zero, which would flatter the estimate.
    """
    observations = list(observations)
    if not observations:
        return None
    buckets: dict[int, list[PairedObservation]] = {}
    for observation in observations:
        index = min(bins - 1, int(observation.probability * bins))
        buckets.setdefault(index, []).append(observation)

    total = len(observations)
    error = 0.0
    for group in buckets.values():
        predicted = sum(o.probability for o in group) / len(group)
        observed = sum(1 for o in group if o.resistant) / len(group)
        error += (len(group) / total) * abs(predicted - observed)
    return error


def calibrate(observations: Sequence[PairedObservation], drug: str,
              site_id: str, min_paired: int = MIN_PAIRED,
              min_positives: int = MIN_POSITIVES,
              min_negatives: int = MIN_NEGATIVES) -> CalibrationResult:
    """Assess whether this site's pairing supports using the probabilities.

    Refusals are returned rather than raised, and they say what is missing and
    how much more of it is needed, because "collect 7 more resistant isolates"
    is actionable where "insufficient data" is not.
    """
    relevant = [o for o in observations
                if o.drug == drug and o.site_id == site_id]
    resistant = [o for o in relevant if o.resistant]
    susceptible = [o for o in relevant if not o.resistant]

    result = CalibrationResult(
        drug=drug, site_id=site_id, n_paired=len(relevant),
        n_resistant=len(resistant), n_susceptible=len(susceptible))
    for observation in relevant:
        key = observation.lineage or "untyped"
        result.lineages[key] = result.lineages.get(key, 0) + 1

    if len(relevant) < min_paired:
        result.refusals.append(
            f"{len(relevant)} paired observation(s); {min_paired} needed "
            f"({min_paired - len(relevant)} more)")
    if len(resistant) < min_positives:
        result.refusals.append(
            f"{len(resistant)} resistant observation(s); {min_positives} "
            f"needed ({min_positives - len(resistant)} more). This is usually "
            f"the binding constraint: at low resistance prevalence a large "
            f"paired set can still contain almost no resistant isolates")
    if len(susceptible) < min_negatives:
        result.refusals.append(
            f"{len(susceptible)} susceptible observation(s); {min_negatives} "
            f"needed ({min_negatives - len(susceptible)} more)")
    if result.refusals:
        return result

    result.ece = expected_calibration_error(relevant)
    result.usable_for_priority = result.ece <= MAX_USABLE_ECE
    return result


def priority_signal(probability: float,
                    calibration: CalibrationResult) -> Optional[str]:
    """A reason to move a variant *up* the laboratory queue, or ``None``.

    Returns a sentence for the workbench, never a call and never a score
    adjustment that could be mistaken for evidence. A low probability returns
    ``None`` — it does not return a reason to deprioritise, because a model
    being unexcited about a variant is not evidence the variant is harmless,
    and letting it push things *down* the queue would be the model quietly
    deciding what never gets tested.
    """
    if not calibration.calibrated or not calibration.usable_for_priority:
        return None
    if probability < PRIORITY_PROBABILITY:
        return None
    return (f"calibrated probability of resistance {probability:.0%} at "
            f"{calibration.site_id} (expected calibration error "
            f"{calibration.ece:.3f} over {calibration.n_paired} pair(s), "
            f"{calibration.n_resistant} resistant); prioritise phenotypic "
            f"validation. This does not establish resistance.")


class PromotionRefused(RuntimeError):
    """Raised when something tries to turn a prediction into a call."""


def refuse_promotion(target_call: str,
                     calibration: Optional[CalibrationResult] = None) -> None:
    """The explicit guard. There is no calibration that unlocks a call.

    This exists so the refusal is a named, tested behaviour rather than the
    absence of a feature — a future contributor looking for "how do I promote a
    calibrated prediction" finds this and its reasoning instead of an
    unguarded gap.
    """
    raise PromotionRefused(
        f"a calibrated prediction cannot be promoted to {target_call}. "
        f"Calibration may raise laboratory priority; it cannot establish "
        f"resistance, and it cannot grant susceptibility because SUSCEPTIBLE "
        f"is the call that admits a drug to a regimen and requires positive "
        f"callable-locus evidence, which a probability is not."
        + ("" if calibration is None else f" ({calibration.describe()})"))
