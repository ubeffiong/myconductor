import tempfile
import unittest
from pathlib import Path

from myconductor.catalogue.profile import load_profile
from myconductor.core.context import SampleContext, InterpretationPolicy, ValidationScope
from myconductor.core.models import QCFinding
from myconductor.core.models import (
    Call,
    DrugEvidence,
    EngineRef,
    Lane,
    Tier,
    VariantIdentity,
)
from myconductor.io.callable_mask import CallableMask
from myconductor.modules.synthesis import (
    EligibilityAssessor,
    EvidenceReconciler,
    collect_discordances,
)


def mask_for(*loci, fraction=0.99, depth=80):
    rows = ["locus\tmean_depth\tcallable_fraction"]
    rows += [f"{l}\t{depth}\t{fraction}" for l in loci]
    path = Path(tempfile.mkdtemp()) / "m.tsv"
    path.write_text("\n".join(rows) + "\n")
    return CallableMask.from_tsv(path)


def catalogued(drug, call=Call.RESISTANT, gene="rpoB", change="S450L"):
    return DrugEvidence(
        drug=drug, call=call, tier=Tier.CATALOGUED, lane=Lane.CATALOGUE,
        confidence=0.98, variant=VariantIdentity(gene=gene, hgvs_p=change),
        engine=EngineRef(name="test-catalogue", version="1"),
        rationale=f"catalogued {call.value}",
    )


def validated_reconciler(*drugs):
    context = SampleContext("s", "i", "site", "mtbc", assay="WGS")
    context.qc = [QCFinding(k, "pass", "reviewed fixture") for k in InterpretationPolicy().required_qc]
    engine = EngineRef("test", "1", "db", "v1")
    policy = InterpretationPolicy("fixture", [
        ValidationScope("mtbc", "WGS", d, "test", "1", "v1", "fixture validation",
                        "reviewer", "2026-01-01", "digest") for d in drugs])
    return EvidenceReconciler(context=context, policy=policy, catalogue_engine=engine,
                              catalogue_sha256="digest", illustrative=False)


class CoverageGateTests(unittest.TestCase):
    def setUp(self):
        self.profile = load_profile()
        self.reconciler = EvidenceReconciler(profile=self.profile)

    def test_no_evidence_and_no_mask_is_not_assessed(self):
        results = self.reconciler.reconcile([], CallableMask.absent())
        for r in results:
            self.assertIs(r.call, Call.NOT_ASSESSED, r.drug)
            self.assertFalse(r.permits_use)

    def test_coverage_without_validation_is_indeterminate(self):
        mask = mask_for("rplC", "rrl")
        results = self.reconciler.reconcile([], mask)
        linezolid = next(r for r in results if r.drug == "linezolid")
        self.assertIs(linezolid.call, Call.INDETERMINATE)
        self.assertFalse(linezolid.permits_use)

    def test_matching_validation_and_coverage_permit_susceptibility(self):
        results = validated_reconciler("linezolid").reconcile([], mask_for("rplC", "rrl"))
        self.assertIs(next(r for r in results if r.drug == "linezolid").call, Call.SUSCEPTIBLE)

    def test_partial_locus_coverage_is_not_enough(self):
        # linezolid needs rplC AND rrl.
        results = self.reconciler.reconcile([], mask_for("rplC"))
        linezolid = next(r for r in results if r.drug == "linezolid")
        self.assertIs(linezolid.call, Call.NOT_ASSESSED)
        self.assertIn("rrl", linezolid.reason)

    def test_every_profile_drug_gets_a_result(self):
        results = self.reconciler.reconcile([], CallableMask.absent())
        self.assertEqual({r.drug for r in results}, set(self.profile.drugs))

    def test_drug_outside_the_profile_is_unsupported(self):
        evidence = [catalogued("ceftazidime")]
        results = self.reconciler.reconcile(evidence, CallableMask.absent())
        drug = next(r for r in results if r.drug == "ceftazidime")
        self.assertIs(drug.call, Call.UNSUPPORTED)
        self.assertFalse(drug.permits_use)

    def test_resistance_wins_over_coverage(self):
        results = self.reconciler.reconcile(
            [catalogued("rifampicin")], mask_for("rpoB"))
        rif = next(r for r in results if r.drug == "rifampicin")
        self.assertIs(rif.call, Call.RESISTANT)

    def test_indeterminate_blocks_susceptibility_despite_coverage(self):
        evidence = [DrugEvidence(
            drug="bedaquiline", call=Call.INDETERMINATE, tier=Tier.INFERRED,
            lane=Lane.EFFLUX_REGULATORY,
            variant=VariantIdentity(gene="Rv0678", hgvs_p="L117R"),
            rationale="efflux inference",
        )]
        results = self.reconciler.reconcile(evidence, mask_for("atpE", "Rv0678"))
        bdq = next(r for r in results if r.drug == "bedaquiline")
        self.assertIs(bdq.call, Call.INDETERMINATE)
        self.assertFalse(bdq.permits_use)


class EngineCoverageAssertionTests(unittest.TestCase):
    def test_engine_assertion_without_validation_is_indeterminate(self):
        evidence = [DrugEvidence(
            drug="linezolid", call=Call.SUSCEPTIBLE, tier=Tier.CATALOGUED,
            lane=Lane.ENGINE,
            engine=EngineRef(name="mykrobe", version="0.13"),
            asserts_coverage=True, rationale="Mykrobe predicts S",
        )]
        results = EvidenceReconciler().reconcile(evidence, CallableMask.absent())
        lzd = next(r for r in results if r.drug == "linezolid")
        self.assertIs(lzd.call, Call.INDETERMINATE)
        self.assertIn("validated", lzd.reason)

    def test_susceptible_without_the_assertion_does_not_license_it(self):
        evidence = [DrugEvidence(
            drug="linezolid", call=Call.SUSCEPTIBLE, tier=Tier.CATALOGUED,
            lane=Lane.ENGINE,
            engine=EngineRef(name="some-tool", version="1"),
            asserts_coverage=False, rationale="reports susceptible",
        )]
        results = EvidenceReconciler().reconcile(evidence, CallableMask.absent())
        lzd = next(r for r in results if r.drug == "linezolid")
        self.assertIs(lzd.call, Call.NOT_ASSESSED)


class DiscordanceTests(unittest.TestCase):
    def test_conflicting_calls_are_retained_as_a_finding(self):
        evidence = [
            catalogued("rifampicin", Call.RESISTANT),
            DrugEvidence(drug="rifampicin", call=Call.SUSCEPTIBLE,
                         tier=Tier.CATALOGUED, lane=Lane.ENGINE,
                         engine=EngineRef(name="other-tool", version="2"),
                         rationale="reports susceptible"),
        ]
        results = EvidenceReconciler().reconcile(evidence, mask_for("rpoB"))
        rif = next(r for r in results if r.drug == "rifampicin")
        self.assertIsNotNone(rif.discordance)
        self.assertIs(rif.call, Call.INDETERMINATE, "unresolved conflicts require adjudication")
        self.assertEqual(len(rif.evidence), 2, "no evidence is discarded")
        self.assertEqual(len(collect_discordances(results)), 1)

    def test_agreement_produces_no_discordance(self):
        evidence = [catalogued("rifampicin"), catalogued("rifampicin")]
        results = EvidenceReconciler().reconcile(evidence, mask_for("rpoB"))
        rif = next(r for r in results if r.drug == "rifampicin")
        self.assertIsNone(rif.discordance)


class EligibilityTests(unittest.TestCase):
    def setUp(self):
        self.profile = load_profile()
        self.assessor = EligibilityAssessor(self.profile)

    def test_nothing_eligible_without_coverage(self):
        results = EvidenceReconciler().reconcile([], CallableMask.absent())
        report = self.assessor.assess(results)
        self.assertFalse(report.any_eligible)
        self.assertTrue(report.requires_clinical_review)

    def test_eligibility_requires_established_susceptibility(self):
        mask = mask_for("atpE", "Rv0678", "ddn", "fbiA", "fbiB", "fbiC",
                        "fgd1", "rplC", "rrl", "gyrA", "gyrB")
        results = validated_reconciler("bedaquiline", "pretomanid", "linezolid", "moxifloxacin").reconcile([], mask)
        report = self.assessor.assess(results)
        bpalm = next(a for a in report.assessments if a.name == "BPaLM")
        self.assertTrue(bpalm.eligible)
        self.assertGreaterEqual(len(bpalm.usable), bpalm.min_required)

    def test_unestablished_drugs_carry_a_reason(self):
        results = EvidenceReconciler().reconcile([], CallableMask.absent())
        report = self.assessor.assess(results)
        for assessment in report.assessments:
            for drug, reason in assessment.unestablished.items():
                self.assertTrue(reason, f"{drug} has no stated reason")

    def test_summary_does_not_claim_no_regimen_exists(self):
        results = EvidenceReconciler().reconcile([], CallableMask.absent())
        summary = self.assessor.assess(results).summary
        self.assertIn("not a conclusion that", summary)


if __name__ == "__main__":
    unittest.main()
