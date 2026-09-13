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
from ..core.context import InterpretationPolicy, SampleContext


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
    susceptible = [ev for ev in evs if ev.call.is_susceptible and ev.scope == "isolate"]
    if not (resistant and susceptible):
        return None
    return Discordance(
        drug=drug,
        calls=tuple(dict.fromkeys(ev.call.value for ev in evs)),
        sources=tuple(dict.fromkeys(ev.source_name for ev in evs)),
        note=(
            f"{len(resistant)} source(s) call {drug} resistant and "
            f"{len(susceptible)} call it susceptible. The conclusion is "
            f"indeterminate; the conflict is unresolved and "
            f"needs phenotypic testing."
        ),
    )


class EvidenceReconciler:
    """Separate analytical adequacy, genomic inference and measured phenotype.

    An empty variant list plus coverage is only a negative genomic screen.
    Susceptibility additionally requires an externally reviewed validation scope.
    Direct, quality-reviewed phenotype applies only to the tested isolate.
    """

    def __init__(self, profile=None, depth_floor=10, callable_fraction_floor=0.95,
                 include_tier2_loci=False, policy=None, context=None,
                 catalogue_engine=None, catalogue_sha256=None, illustrative=True):
        if depth_floor <= 0 or not 0 < callable_fraction_floor <= 1:
            raise ValueError("invalid coverage thresholds")
        self.profile = profile or load_profile()
        self.depth_floor = depth_floor
        self.callable_fraction_floor = callable_fraction_floor
        self.include_tier2_loci = include_tier2_loci
        self.policy = policy or InterpretationPolicy()
        self.context = context
        self.catalogue_engine = catalogue_engine
        self.catalogue_sha256 = catalogue_sha256
        self.illustrative = illustrative

    def reconcile(self, evidence, mask=None, drugs=None):
        mask = mask or CallableMask.absent()
        by_drug = {}
        for ev in evidence:
            if ev.tier is Tier.PHENOTYPIC:
                if (self.context is None or ev.sample_id != self.context.sample_id
                        or not ev.observation_id):
                    raise ValueError("phenotypic evidence requires the current sample and observation ID")
            by_drug.setdefault(ev.drug, []).append(ev)
        universe = set(drugs if drugs is not None else self.profile.drugs) | set(by_drug)
        return [self._resolve(drug, by_drug.get(drug, []), mask) for drug in sorted(universe)]

    def _resolve(self, drug, evs, mask):
        genomic = [e for e in evs if e.tier is not Tier.PHENOTYPIC]
        phenotypes = [e for e in evs if e.tier is Tier.PHENOTYPIC]
        result = self._genomic(drug, genomic, mask)
        result.genomic_call = result.call
        result.evidence = evs
        result.discordance = _discordance(drug, evs)
        if not self.profile.supports_drug(drug):
            return result
        if phenotypes:
            calls = {e.call for e in phenotypes}
            measured = next(iter(calls)) if len(calls) == 1 else Call.INDETERMINATE
            result.phenotypic_call = measured
            if measured.is_established:
                result.call = measured
                result.tier = Tier.PHENOTYPIC
                result.reason = "Measured phenotype for this isolate: " + phenotypes[0].rationale
            else:
                result.call = Call.INDETERMINATE
                result.reason = "Phenotype measurements conflict or lack adequate quality/context."
        if result.discordance:
            result.call = Call.INDETERMINATE
            result.confidence = None
            result.reason = result.discordance.note
        return result

    def _genomic(self, drug, evs, mask):
        if not self.profile.supports_drug(drug):
            return DrugResult(drug, Call.UNSUPPORTED, Tier.NONE, evidence=evs,
                              reason=f"{drug} is outside the {self.profile.name} profile")
        loci = self.profile.required_loci(drug, self.include_tier2_loci)
        ok, coverage, reasons = mask.assess(loci, depth_floor=self.depth_floor,
                                          fraction_floor=self.callable_fraction_floor)
        status = "adequate" if ok else "insufficient"
        def result(call, tier, reason):
            return DrugResult(drug, call, tier, evidence=evs, coverage=coverage,
                              reason=reason, assay_status=status)
        discordance = _discordance(drug, evs)
        if discordance:
            return result(Call.INDETERMINATE, _best_tier(evs), discordance.note)
        for call in (Call.RESISTANT, Call.INDETERMINATE, Call.NO_CALL):
            matched = [e for e in evs if e.call is call]
            if matched:
                return result(call, _best_tier(matched), matched[0].rationale)
        # Every NTM macrolide negative call needs the relevant reference state,
        # including functional wild-type alleles invisible in a variant-only VCF.
        if self.profile.name == "mabscessus" and drug == "clarithromycin":
            c = self.context
            if c is None or not c.subspecies or c.erm41_status != "nonfunctional" or c.rrl_status != "wild_type":
                return result(Call.INDETERMINATE, Tier.NONE,
                              "M. abscessus macrolide interpretation requires subspecies, "
                              "nonfunctional erm(41) evidence and assessed rrl status; "
                              "inducible resistance cannot be excluded from a variant list.")
        asserted = [e for e in evs if e.call is Call.SUSCEPTIBLE
                    and e.scope == "isolate" and e.asserts_coverage]
        for ev in asserted:
            if self.policy.permits(drug, self.context, ev.engine):
                r = result(Call.SUSCEPTIBLE, ev.tier,
                           f"Susceptibility asserted by {ev.source_name} within the "
                           f"externally reviewed scope of policy {self.policy.version}.")
                r.assay_status = "engine_assessed"
                return r
        if not ok and not asserted:
            return result(Call.NOT_ASSESSED, Tier.NONE,
                          "cannot establish susceptibility: " + "; ".join(reasons))
        if (ok and not self.illustrative and self.policy.permits(
                drug, self.context, self.catalogue_engine, self.catalogue_sha256)):
            return result(Call.SUSCEPTIBLE, Tier.CATALOGUED,
                          "No resistance-associated variant detected in callable loci; "
                          f"interpretation covered by policy {self.policy.version}.")
        return result(Call.INDETERMINATE, Tier.NONE,
                      "No resistance-associated variant detected, but susceptibility "
                      "is not established: no matching validated organism/drug/assay/"
                      "engine/database scope with required QC. Coverage alone is insufficient.")


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
