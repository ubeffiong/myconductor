"""Evidence reconciliation and guideline eligibility.

This module contains the fix for the defect that invalidated every regimen the
tool previously produced. The old rule was:

    if r is None:
        return True          # no evidence of resistance -> presumed usable

That treats "we never looked at this drug's loci" as "this drug will work".
A drug can have no resistance result because its locus was not covered, the
caller failed, the gene was absent from the input, the input carried variants
but no reference-call information, the profile does not support the drug, or the
mechanism is simply unknown. None of those are susceptibility.

The rule is now inverted and coverage-gated. A drug reaches ``SUSCEPTIBLE``
only when *positive* evidence exists that its required loci were callable, and
that evidence has to come from outside the variant list (see
``io.callable_mask``). Everything else is ``NOT_ASSESSED``, ``INDETERMINATE``,
``NO_CALL`` or ``UNSUPPORTED`` — and only ``SUSCEPTIBLE`` counts toward a
regimen.

The second change is what this module outputs. It no longer proposes a regimen.
Regimen construction depends on treatment history, disease site, age,
pregnancy, comorbidity, drug interactions, toxicity, drug availability,
baseline ECG and laboratory findings, and national policy — none of which this
tool sees. What it can honestly produce is a per-regimen **eligibility
assessment** against the genotypic evidence, flagged for clinical review.
"""
from __future__ import annotations

from typing import Iterable, Optional

from ..catalogue.profile import OrganismProfile, load_profile
from ..core.models import (
    Call,
    Discordance,
    DrugEvidence,
    DrugResult,
    EligibilityReport,
    LocusCoverage,
    RegimenAssessment,
    Tier,
)
from ..io.callable_mask import CallableMask


def _best_tier(evs: Iterable[DrugEvidence]) -> Tier:
    tiers = [ev.tier for ev in evs]
    return max(tiers, key=lambda t: t.rank) if tiers else Tier.NONE


def _confidence(evs: Iterable[DrugEvidence]) -> Optional[float]:
    vals = [ev.confidence for ev in evs if ev.confidence is not None]
    return max(vals) if vals else None


def _discordance(drug: str, evs: list[DrugEvidence]) -> Optional[Discordance]:
    """Flag materially conflicting verdicts rather than collapsing them.

    Two sources disagreeing about a drug is among the most clinically useful
    things this platform can report, so it is retained as a finding instead of
    being resolved away by picking a winner.
    """
    resistant = [ev for ev in evs if ev.call.is_resistant]
    susceptible = [ev for ev in evs if ev.call.is_susceptible]
    if not (resistant and susceptible):
        return None
    return Discordance(
        drug=drug,
        calls=tuple(dict.fromkeys(ev.call.value for ev in evs)),
        sources=tuple(dict.fromkeys(ev.source_name for ev in evs)),
        note=(
            f"{len(resistant)} source(s) call {drug} resistant and "
            f"{len(susceptible)} call it susceptible. Resistance is carried "
            f"forward as the safer position; the conflict is unresolved and "
            f"needs phenotypic testing."
        ),
    )


class EvidenceReconciler:
    """Turns lane evidence plus coverage evidence into per-drug results."""

    def __init__(
        self,
        profile: Optional[OrganismProfile] = None,
        depth_floor: int = 10,
        callable_fraction_floor: float = 0.95,
        include_tier2_loci: bool = False,
    ):
        self.profile = profile or load_profile()
        self.depth_floor = depth_floor
        self.callable_fraction_floor = callable_fraction_floor
        self.include_tier2_loci = include_tier2_loci

    def reconcile(
        self,
        evidence: list[DrugEvidence],
        mask: Optional[CallableMask] = None,
        drugs: Optional[Iterable[str]] = None,
    ) -> list[DrugResult]:
        """Produce one result per drug in the profile.

        Every profile drug gets a result, including drugs nothing was found
        for: a drug missing from the output is indistinguishable from a drug
        that was assessed and cleared, which is the ambiguity this whole module
        exists to remove.
        """
        mask = mask or CallableMask.absent()
        by_drug: dict[str, list[DrugEvidence]] = {}
        for ev in evidence:
            by_drug.setdefault(ev.drug, []).append(ev)

        universe = list(drugs) if drugs is not None else list(self.profile.drugs)
        # Include any drug evidence arrived for, even if outside the profile,
        # so nothing is silently discarded.
        for drug in by_drug:
            if drug not in universe:
                universe.append(drug)

        return [self._resolve(drug, by_drug.get(drug, []), mask)
                for drug in sorted(universe)]

    def _resolve(self, drug: str, evs: list[DrugEvidence],
                 mask: CallableMask) -> DrugResult:
        if not self.profile.supports_drug(drug):
            return DrugResult(
                drug=drug, call=Call.UNSUPPORTED, tier=Tier.NONE,
                evidence=evs,
                reason=(f"{drug} is not in the {self.profile.name} profile "
                        f"v{self.profile.version}"),
            )

        discordance = _discordance(drug, evs)

        # 1. Established resistance wins, and only graded or phenotypic
        #    evidence can establish it (enforced in DrugEvidence).
        resistant = [ev for ev in evs if ev.call.is_resistant]
        if resistant:
            return DrugResult(
                drug=drug, call=Call.RESISTANT, tier=_best_tier(resistant),
                confidence=_confidence(resistant), evidence=evs,
                discordance=discordance,
                reason=resistant[0].rationale,
            )

        # 2. A lane withheld susceptibility without asserting resistance.
        indeterminate = [ev for ev in evs if ev.call is Call.INDETERMINATE]
        if indeterminate:
            return DrugResult(
                drug=drug, call=Call.INDETERMINATE,
                tier=_best_tier(indeterminate), evidence=evs,
                discordance=discordance,
                reason=indeterminate[0].rationale,
            )

        # 3. An unreliable genotype at the locus.
        no_call = [ev for ev in evs if ev.call is Call.NO_CALL]
        if no_call:
            return DrugResult(
                drug=drug, call=Call.NO_CALL, tier=_best_tier(no_call),
                evidence=evs, discordance=discordance,
                reason=no_call[0].rationale,
            )

        # 4. No resistance evidence. This is the branch that used to return
        #    "usable". Susceptibility now has to be earned with coverage.
        loci = self.profile.required_loci(drug, self.include_tier2_loci)
        callable_ok, coverages, reasons = mask.assess(
            loci,
            depth_floor=self.depth_floor,
            fraction_floor=self.callable_fraction_floor,
        )

        # 4a. An external engine may have done its own coverage assessment.
        #     Honour it, and attribute it, rather than discarding a validated
        #     tool's susceptible call for want of our own mask.
        asserted = [ev for ev in evs
                    if ev.asserts_coverage and ev.call.is_susceptible]
        if not callable_ok and asserted:
            sources = ", ".join(dict.fromkeys(ev.source_name for ev in asserted))
            return DrugResult(
                drug=drug, call=Call.SUSCEPTIBLE, tier=_best_tier(asserted),
                confidence=_confidence(asserted), evidence=evs,
                coverage=coverages, discordance=discordance,
                reason=(f"susceptibility asserted by {sources}, which performed "
                        f"its own callable-locus assessment; Myconductor's own "
                        f"mask did not cover it ({'; '.join(reasons)})"),
            )

        if not callable_ok:
            return DrugResult(
                drug=drug, call=Call.NOT_ASSESSED, tier=Tier.NONE,
                evidence=evs, coverage=coverages, discordance=discordance,
                reason=("cannot establish susceptibility: " + "; ".join(reasons)),
            )

        return DrugResult(
            drug=drug,
            call=Call.SUSCEPTIBLE,
            tier=_best_tier(evs) if evs else Tier.CATALOGUED,
            confidence=_confidence(evs),
            evidence=evs,
            coverage=coverages,
            discordance=discordance,
            reason=(
                f"no resistance-associated variant found across "
                f"{len(loci)} callable locus/loci "
                f"({', '.join(sorted(loci))}); "
                f"coverage source: {mask.source}"
            ),
        )


class EligibilityAssessor:
    """Assesses guideline eligibility. Does not prescribe.

    The output names which regimens the genotypic evidence is *consistent
    with*, and which drugs fall short and why. It deliberately does not emit a
    drug list to administer.
    """

    def __init__(self, profile: Optional[OrganismProfile] = None):
        self.profile = profile or load_profile()

    def assess(self, results: list[DrugResult]) -> EligibilityReport:
        by_drug = {r.drug: r for r in results}
        assessments: list[RegimenAssessment] = []

        for regimen in self.profile.regimens:
            usable, resistant, unestablished = [], [], {}
            for drug in regimen.drugs:
                r = by_drug.get(drug)
                if r is None:
                    unestablished[drug] = "no result produced for this drug"
                elif r.permits_use:
                    usable.append(drug)
                elif r.call.is_resistant:
                    resistant.append(drug)
                else:
                    unestablished[drug] = (
                        r.reason or r.call.reason_unestablished or "unknown"
                    )

            shortfall = max(0, regimen.min_required - len(usable))
            eligible = shortfall == 0
            if eligible:
                note = (f"{len(usable)}/{len(regimen.drugs)} companion drugs have "
                        f"established, coverage-backed susceptibility "
                        f"(>= {regimen.min_required} required).")
            else:
                parts = [f"{len(usable)}/{regimen.min_required} required drugs have "
                         f"established susceptibility"]
                if resistant:
                    parts.append(f"resistant: {', '.join(sorted(resistant))}")
                if unestablished:
                    parts.append(
                        "not established: " + ", ".join(sorted(unestablished))
                    )
                note = "; ".join(parts) + "."

            assessments.append(RegimenAssessment(
                name=regimen.name, drugs=list(regimen.drugs),
                min_required=regimen.min_required, usable=sorted(usable),
                resistant=sorted(resistant), unestablished=unestablished,
                eligible=eligible, shortfall=shortfall, note=note,
            ))

        eligible_names = [a.name for a in assessments if a.eligible]
        resistant_drugs = sorted(r.drug for r in results if r.call.is_resistant)
        unestablished_drugs = {
            r.drug: (r.reason or r.call.reason_unestablished or "unknown")
            for r in results if not r.call.is_established
        }

        if eligible_names:
            summary = (
                f"Genotypic evidence is consistent with: {', '.join(eligible_names)}. "
                f"Regimen selection requires clinical assessment this tool does "
                f"not perform (treatment history, site of disease, comorbidity, "
                f"interactions, toxicity, availability, national policy)."
            )
        else:
            summary = (
                "No regimen reaches its required count of drugs with "
                "established, coverage-backed susceptibility. This is a "
                "statement about the available evidence, not a conclusion that "
                "no regimen exists"
                + (f" — {len(unestablished_drugs)} drug(s) could not be "
                   f"assessed at all." if unestablished_drugs else ".")
            )

        return EligibilityReport(
            assessments=assessments,
            any_eligible=bool(eligible_names),
            resistant_drugs=resistant_drugs,
            unestablished_drugs=unestablished_drugs,
            summary=summary,
            requires_clinical_review=True,
        )


def collect_discordances(results: list[DrugResult]) -> list[Discordance]:
    return [r.discordance for r in results if r.discordance is not None]
