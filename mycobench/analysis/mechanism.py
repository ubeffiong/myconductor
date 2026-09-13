"""Mechanism class from MIC evidence, rather than from unmeasured expression.

The problem
-----------
Efflux-mediated resistance cannot be read off DNA: the presence of a pump
regulator says nothing about whether the pump is overexpressed. The established
way to settle it is a wet-lab efflux-inhibitor assay — a verapamil MIC shift.

Two routes that need no new data
--------------------------------
**Distribution shape.** An expression-mediated effect should move MIC modestly
and continuously; a target-site knockout should move it far and discretely. So
the *size* of the shift carries mechanism information. Implemented as a
two-component fit on log2 MIC, and **refused outright** when censoring makes the
distribution's shape unrecoverable — which on the real CRyPTIC release rules
this out for amikacin, delamanid, clofazimine and rifabutin.

**Cross-drug covariance.** Pump de-repression raises a *correlated set* of MICs:
loss of ``mmpR5`` raises bedaquiline and clofazimine together, because both are
substrates of MmpS5-MmpL5. A target-site mutation in ``atpE`` raises bedaquiline
alone. CRyPTIC measures all 13 drugs on one plate, so the joint structure is
already collected.

This second route is the one part of this plan I could not find in the
literature, so it ships with the experiment that would kill it:
``validate_discriminator`` runs the known contrast — ``Rv0678`` carriers versus
``atpE`` carriers on the bedaquiline/clofazimine pair — and reports whether the
discriminator separates them. If it does not, the idea is wrong and this module
should be deleted rather than tuned.
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

from . import mic, stats

#: Drug sets known to share an efflux route, so a correlated shift across the
#: set is the efflux signature. Sourced from the pump's substrate range, not
#: inferred from data, so the test below is not circular.
SHARED_EFFLUX_SETS: dict[str, tuple[str, ...]] = {
    "mmpS5_mmpL5": ("bedaquiline", "clofazimine"),
}

#: Genes whose variants act at the drug target rather than through a pump.
#: Used only as the control arm of the validation experiment.
TARGET_SITE_GENES: dict[str, tuple[str, ...]] = {
    "bedaquiline": ("atpE",),
    "rifampicin": ("rpoB",),
    "isoniazid": ("katG",),
}

#: Efflux regulators, for the same purpose.
EFFLUX_REGULATOR_GENES = ("Rv0678", "mmpR5", "pepQ", "Rv1979c")

#: A two-component fit needs enough exact observations to have a shape at all.
MIN_EXACT_FOR_SHAPE = 40
#: And censoring below this on both sides, or the modes are artefacts of the
#: plate's boundaries rather than of biology.
MAX_CENSORING_FOR_SHAPE = 0.25


# -- distribution shape ---------------------------------------------------
@dataclass
class TwoComponentFit:
    """A two-component Gaussian mixture on log2 MIC, or a refusal."""

    means: Optional[tuple[float, float]] = None
    sds: Optional[tuple[float, float]] = None
    weights: Optional[tuple[float, float]] = None
    n_exact: int = 0
    iterations: int = 0
    refused_because: Optional[str] = None

    @property
    def available(self) -> bool:
        return self.means is not None

    @property
    def separation(self) -> Optional[float]:
        """Distance between component means, in doubling dilutions."""
        if not self.available:
            return None
        return abs(self.means[1] - self.means[0])

    def describe(self) -> str:
        if not self.available:
            return f"not fitted: {self.refused_because}"
        return (f"components at {self.means[0]:+.2f} and {self.means[1]:+.2f} "
                f"log2 (separation {self.separation:.2f} dilutions), weights "
                f"{self.weights[0]:.2f}/{self.weights[1]:.2f}, "
                f"n_exact={self.n_exact}")


def fit_two_component(observations: Sequence[mic.Observation],
                      iterations: int = 200,
                      tolerance: float = 1e-6) -> TwoComponentFit:
    """Fit two Gaussians to the exact log2 MICs by expectation-maximisation.

    Only exact observations are used. A censored value carries no information
    about *where* in the tail it sits, so including it at its bound would pile
    mass onto the plate boundary and invent a component there.
    """
    exact = [o.log2_bound for o in observations if o.exact]
    profile = mic.censoring_profile(observations)

    if len(exact) < MIN_EXACT_FOR_SHAPE:
        return TwoComponentFit(
            n_exact=len(exact),
            refused_because=(f"{len(exact)} exact observation(s), below the "
                             f"{MIN_EXACT_FOR_SHAPE} needed to describe a "
                             f"distribution's shape"))
    if (profile.left_fraction > MAX_CENSORING_FOR_SHAPE
            or profile.right_fraction > MAX_CENSORING_FOR_SHAPE):
        return TwoComponentFit(
            n_exact=len(exact),
            refused_because=(f"censoring is {profile.left_fraction:.0%} left / "
                             f"{profile.right_fraction:.0%} right; above "
                             f"{MAX_CENSORING_FOR_SHAPE:.0%} the apparent modes "
                             f"are the plate's boundaries, not the isolates'"))

    low, high = min(exact), max(exact)
    if math.isclose(low, high):
        return TwoComponentFit(
            n_exact=len(exact),
            refused_because="every exact observation is the same value")

    # Initialise the components at the quartiles: deterministic, so the fit is
    # reproducible, and far enough apart to avoid collapsing immediately.
    means = [stats.quantile(exact, 0.25), stats.quantile(exact, 0.75)]
    spread = max(statistics.pstdev(exact), 0.25)
    sds = [spread, spread]
    weights = [0.5, 0.5]

    def density(x: float, mu: float, sd: float) -> float:
        sd = max(sd, 1e-3)
        return (1.0 / (sd * math.sqrt(2 * math.pi))
                * math.exp(-0.5 * ((x - mu) / sd) ** 2))

    previous = None
    performed = 0
    for performed in range(1, iterations + 1):
        responsibilities = []
        log_likelihood = 0.0
        for x in exact:
            a = weights[0] * density(x, means[0], sds[0])
            b = weights[1] * density(x, means[1], sds[1])
            total = a + b
            if total <= 0:
                responsibilities.append(0.5)
                continue
            responsibilities.append(a / total)
            log_likelihood += math.log(total)

        mass_a = sum(responsibilities)
        mass_b = len(exact) - mass_a
        if mass_a < 1.0 or mass_b < 1.0:
            return TwoComponentFit(
                n_exact=len(exact), iterations=performed,
                refused_because="one component collapsed; the data support a "
                                "single mode rather than two")
        means = [
            sum(r * x for r, x in zip(responsibilities, exact)) / mass_a,
            sum((1 - r) * x for r, x in zip(responsibilities, exact)) / mass_b,
        ]
        sds = [
            math.sqrt(max(1e-6, sum(r * (x - means[0]) ** 2
                                    for r, x in zip(responsibilities, exact))
                          / mass_a)),
            math.sqrt(max(1e-6, sum((1 - r) * (x - means[1]) ** 2
                                    for r, x in zip(responsibilities, exact))
                          / mass_b)),
        ]
        weights = [mass_a / len(exact), mass_b / len(exact)]
        if previous is not None and abs(log_likelihood - previous) < tolerance:
            break
        previous = log_likelihood

    order = sorted(range(2), key=lambda i: means[i])
    return TwoComponentFit(
        means=(means[order[0]], means[order[1]]),
        sds=(sds[order[0]], sds[order[1]]),
        weights=(weights[order[0]], weights[order[1]]),
        n_exact=len(exact), iterations=performed)


# -- cross-drug covariance ------------------------------------------------
@dataclass
class JointShift:
    """How a variant moves MICs across a set of drugs sharing an efflux route."""

    variant: str
    drugs: tuple[str, ...]
    shifts: dict[str, mic.StochasticShift] = field(default_factory=dict)

    @property
    def estimable(self) -> list[str]:
        return [d for d, s in self.shifts.items() if s.available]

    @property
    def deltas(self) -> dict[str, float]:
        return {d: s.delta for d, s in self.shifts.items() if s.available}

    @property
    def concordant(self) -> Optional[bool]:
        """Did every estimable drug move in the same direction?"""
        deltas = list(self.deltas.values())
        if len(deltas) < 2:
            return None
        return all(d > 0 for d in deltas) or all(d < 0 for d in deltas)

    @property
    def min_magnitude(self) -> Optional[float]:
        deltas = [abs(d) for d in self.deltas.values()]
        return min(deltas) if deltas else None

    def describe(self) -> str:
        if len(self.estimable) < 2:
            return (f"{self.variant}: fewer than two drugs estimable "
                    f"({', '.join(self.estimable) or 'none'})")
        parts = ", ".join(f"{d}={v:+.3f}" for d, v in sorted(self.deltas.items()))
        return (f"{self.variant}: {parts}; "
                f"{'concordant' if self.concordant else 'discordant'}")


def joint_shift(variant: str, drugs: Sequence[str],
                panels: dict[str, mic.DrugPanel],
                carriers: Iterable[str],
                non_carriers: Iterable[str]) -> JointShift:
    """Shift for one variant across several drugs, on the same isolates."""
    carriers = list(carriers)
    non_carriers = list(non_carriers)
    result = JointShift(variant=variant, drugs=tuple(drugs))
    for drug in drugs:
        panel = panels.get(drug)
        if panel is None:
            continue
        result.shifts[drug] = mic.stochastic_shift(
            panel.subset(carriers), panel.subset(non_carriers))
    return result


@dataclass
class MechanismCall:
    """An exploratory mechanism class. Never a resistance call."""

    variant: str
    mechanism: str          # efflux-like | target-like | indeterminate
    confidence: str         # exploratory, always
    joint: Optional[JointShift] = None
    reasons: list[str] = field(default_factory=list)

    def describe(self) -> str:
        return (f"{self.variant}: {self.mechanism} ({self.confidence})"
                + ("; " + "; ".join(self.reasons) if self.reasons else ""))


def classify(variant: str, panels: dict[str, mic.DrugPanel],
             carriers: Iterable[str], non_carriers: Iterable[str],
             efflux_set: str = "mmpS5_mmpL5",
             min_delta: float = 0.15) -> MechanismCall:
    """Classify a variant's mechanism from its cross-drug MIC signature.

    Exploratory by construction: the output is a mechanism hypothesis for
    laboratory triage, and it never contributes to a drug's resistance verdict.
    """
    drugs = SHARED_EFFLUX_SETS.get(efflux_set)
    if not drugs:
        raise KeyError(f"unknown efflux set {efflux_set!r}")

    carriers, non_carriers = list(carriers), list(non_carriers)
    shift = joint_shift(variant, drugs, panels, carriers, non_carriers)
    call = MechanismCall(variant=variant, mechanism="indeterminate",
                         confidence="exploratory", joint=shift)

    if len(shift.estimable) < 2:
        call.reasons.append(
            f"only {len(shift.estimable)} of {len(drugs)} drug(s) estimable; "
            f"a shared-route signature needs both")
        return call

    magnitude = shift.min_magnitude or 0.0
    if shift.concordant and magnitude >= min_delta:
        call.mechanism = "efflux-like"
        call.reasons.append(
            f"both {' and '.join(drugs)} shift in the same direction with "
            f"magnitude at least {magnitude:.3f}, consistent with a shared "
            f"efflux route rather than two independent target changes")
        return call

    moved = [d for d, v in shift.deltas.items() if abs(v) >= min_delta]
    if len(moved) == 1:
        call.mechanism = "target-like"
        call.reasons.append(
            f"only {moved[0]} shifted materially; a shared pump would move "
            f"both substrates")
        return call

    call.reasons.append(
        f"no drug shifted by at least {min_delta:.2f}, or the shifts "
        f"disagreed in direction")
    return call


# -- the experiment that would kill this ---------------------------------
@dataclass
class DiscriminatorResult:
    """Does the cross-drug signature separate a known pump from a known target?"""

    efflux_variants: list[str] = field(default_factory=list)
    target_variants: list[str] = field(default_factory=list)
    efflux_calls: dict[str, str] = field(default_factory=dict)
    target_calls: dict[str, str] = field(default_factory=dict)
    verdict: str = "not-run"
    reasons: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.verdict == "discriminates"


#: A reference group smaller than this cannot support a comparison.
MIN_REFERENCE = 5


def validate_discriminator(panels: dict[str, mic.DrugPanel],
                           carriers_by_variant: dict[str, list[str]],
                           all_isolates: Sequence[str],
                           efflux_genes: Sequence[str] = EFFLUX_REGULATOR_GENES,
                           target_gene: str = "atpE",
                           min_agreement: float = 0.60,
                           min_reference: int = MIN_REFERENCE
                           ) -> DiscriminatorResult:
    """Run the known contrast, and say plainly whether the method works.

    ``mmpR5``/``Rv0678`` variants are the positive control for the efflux class
    and ``atpE`` variants for the target class, because the bedaquiline
    literature settles what each one is. A discriminator that cannot reproduce
    that separation is not measuring mechanism, and the honest response is to
    abandon it rather than adjust the threshold until it agrees.

    The reference group is every isolate carrying **none** of the variants
    under test. Using "everyone who lacks this variant" instead is wrong and
    quietly inverts the result: the atpE carriers would then be compared
    against Rv0678 carriers, whose clofazimine MICs are raised, so atpE would
    appear to *lower* clofazimine and be classified indeterminate. A reference
    group that is itself resistant is not a reference group.
    """
    result = DiscriminatorResult()

    tested: set[str] = set()
    for carriers in carriers_by_variant.values():
        tested |= set(carriers)
    reference = [i for i in all_isolates if i not in tested]

    if len(reference) < min_reference:
        result.verdict = "not-run"
        result.reasons.append(
            f"{len(reference)} isolate(s) carry none of the variants under "
            f"test; at least {min_reference} are needed as a reference group "
            f"that is not itself resistant")
        return result

    for variant, carriers in carriers_by_variant.items():
        gene = variant.split("_", 1)[0]
        if gene in efflux_genes:
            result.efflux_variants.append(variant)
            result.efflux_calls[variant] = classify(
                variant, panels, carriers, reference).mechanism
        elif gene == target_gene:
            result.target_variants.append(variant)
            result.target_calls[variant] = classify(
                variant, panels, carriers, reference).mechanism

    if not result.efflux_calls or not result.target_calls:
        result.verdict = "not-run"
        result.reasons.append(
            f"need carriers of both an efflux regulator and {target_gene}; "
            f"found {len(result.efflux_calls)} and "
            f"{len(result.target_calls)}")
        return result

    result.reasons.append(
        f"compared against {len(reference)} isolate(s) carrying none of the "
        f"variants under test")

    efflux_right = sum(1 for m in result.efflux_calls.values()
                       if m == "efflux-like") / len(result.efflux_calls)
    target_right = sum(1 for m in result.target_calls.values()
                       if m == "target-like") / len(result.target_calls)

    result.reasons.append(
        f"efflux regulators classified efflux-like: {efflux_right:.0%} "
        f"({len(result.efflux_calls)} variant(s)); {target_gene} classified "
        f"target-like: {target_right:.0%} ({len(result.target_calls)})")

    if efflux_right >= min_agreement and target_right >= min_agreement:
        result.verdict = "discriminates"
    else:
        result.verdict = "does-not-discriminate"
        result.reasons.append(
            "the cross-drug signature does not reproduce a contrast the "
            "literature already settles, so it carries no mechanism "
            "information and should be dropped rather than retuned")
    return result
