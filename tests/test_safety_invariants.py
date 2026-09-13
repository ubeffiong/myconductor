"""The invariants that make this codebase safe to run.

Each test here corresponds to a defect in the previous implementation. They are
kept in one file, separate from the unit tests, because a failure in any of them
is a clinical-safety regression rather than a functional bug.
"""
import unittest
from pathlib import Path

from myconductor import Myconductor
from myconductor.core.models import (
    Call,
    DrugEvidence,
    Lane,
    Tier,
    Variant,
)
from myconductor.modules.features import DemoAnnotator

DATA = Path(__file__).resolve().parent.parent / "myconductor" / "data"
VCF = DATA / "example_input.vcf"
MASK = DATA / "example_callable.tsv"


class MissingEvidenceIsNotSusceptibility(unittest.TestCase):
    """The defect: `if r is None: return True  # presumed usable`."""

    def test_no_mask_means_nothing_is_susceptible(self):
        report = Myconductor(platform="illumina").analyze(VCF)
        self.assertEqual(
            [r.drug for r in report.drug_results if r.permits_use], [],
            "no drug may be susceptible without callable-locus evidence",
        )

    def test_no_mask_means_nothing_is_eligible(self):
        report = Myconductor(platform="illumina").analyze(VCF)
        self.assertFalse(report.eligibility.any_eligible)
        for assessment in report.eligibility.assessments:
            self.assertEqual(assessment.usable, [])

    def test_drug_with_no_evidence_is_not_assessed(self):
        report = Myconductor(platform="illumina").analyze(VCF)
        linezolid = report.result_for("linezolid")
        self.assertIsNotNone(linezolid)
        self.assertIs(linezolid.call, Call.NOT_ASSESSED)
        self.assertIn("callable", (linezolid.reason or "").lower())

    def test_coverage_evidence_is_what_licenses_susceptibility(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        linezolid = report.result_for("linezolid")
        self.assertIs(linezolid.call, Call.SUSCEPTIBLE)
        self.assertTrue(linezolid.permits_use)

    def test_locus_absent_from_mask_stays_unassessed(self):
        # rrs is deliberately absent from the demo mask.
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        amikacin = report.result_for("amikacin")
        self.assertFalse(amikacin.permits_use)


class OnlyGradedEvidenceEstablishesResistance(unittest.TestCase):
    """Rule-based and model lanes must not assert resistance."""

    def test_inferred_tier_cannot_emit_resistant(self):
        with self.assertRaises(ValueError) as ctx:
            DrugEvidence(drug="bedaquiline", call=Call.RESISTANT,
                         tier=Tier.INFERRED, lane=Lane.EFFLUX_REGULATORY)
        self.assertIn("only catalogued or phenotypic", str(ctx.exception))

    def test_predicted_tier_cannot_emit_resistant(self):
        with self.assertRaises(ValueError):
            DrugEvidence(drug="isoniazid", call=Call.RESISTANT,
                         tier=Tier.PREDICTED, lane=Lane.VUS)

    def test_catalogued_tier_may_emit_resistant(self):
        ev = DrugEvidence(drug="isoniazid", call=Call.RESISTANT,
                          tier=Tier.CATALOGUED, lane=Lane.CATALOGUE)
        self.assertIs(ev.call, Call.RESISTANT)

    def test_efflux_withholds_rather_than_asserts(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        bedaquiline = report.result_for("bedaquiline")
        # Rv0678 L117R is present; the old code called this predicted-resistant.
        self.assertIs(bedaquiline.call, Call.INDETERMINATE)
        self.assertFalse(bedaquiline.permits_use)


class MultiDrugEffectsSurvive(unittest.TestCase):
    """The defect: efflux returned `drugs[0]`, dropping clofazimine."""

    def test_rv0678_reaches_both_bedaquiline_and_clofazimine(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        for drug in ("bedaquiline", "clofazimine"):
            result = report.result_for(drug)
            self.assertIsNotNone(result, f"{drug} missing from results")
            lanes = {ev.lane for ev in result.evidence}
            self.assertIn(Lane.EFFLUX_REGULATORY, lanes,
                          f"{drug} lost its efflux evidence")


class SyntheticValuesCannotReachAReport(unittest.TestCase):
    """The defect: SHA-256 digests presented as conservation scores."""

    def test_demo_annotator_requires_acknowledgement(self):
        with self.assertRaises(ValueError):
            DemoAnnotator()

    def test_pipeline_refuses_synthetic_annotator_without_demo_mode(self):
        with self.assertRaises(ValueError) as ctx:
            Myconductor(annotator=DemoAnnotator(acknowledged=True))
        self.assertIn("synthetic", str(ctx.exception).lower())

    def test_demo_mode_stamps_the_report(self):
        conductor = Myconductor(annotator=DemoAnnotator(acknowledged=True),
                                demo_mode=True, platform="illumina")
        report = conductor.analyze(VCF, mask_path=MASK)
        self.assertTrue(report.demo_mode)
        from myconductor.reporting.render import render_text
        self.assertIn("SYNTHETIC", render_text(report))

    def test_default_pipeline_is_not_synthetic(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        self.assertFalse(report.demo_mode)


class UninterpretedVariantsBlockSusceptibility(unittest.TestCase):
    def test_vus_in_a_required_locus_prevents_a_clean_susceptible(self):
        # pncA D12A is not in the catalogue, and pncA is pyrazinamide's tier-1
        # locus. Coverage alone must not make pyrazinamide susceptible.
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        pyrazinamide = report.result_for("pyrazinamide")
        self.assertIs(pyrazinamide.call, Call.INDETERMINATE)
        self.assertFalse(pyrazinamide.permits_use)

    def test_vus_lane_never_emits_a_resistance_call(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        for result in report.drug_results:
            for ev in result.evidence:
                if ev.lane is Lane.VUS:
                    self.assertIsNot(ev.call, Call.RESISTANT)
                    self.assertIsNot(ev.call, Call.SUSCEPTIBLE)


class EligibilityIsNotAPrescription(unittest.TestCase):
    def test_report_carries_no_proposed_drug_list(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        self.assertFalse(hasattr(report.eligibility, "proposed"))
        for assessment in report.eligibility.assessments:
            self.assertFalse(hasattr(assessment, "proposed"))

    def test_clinical_review_is_always_required(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        self.assertTrue(report.eligibility.requires_clinical_review)


class DiscoveryIsDecoupled(unittest.TestCase):
    def test_pipeline_module_does_not_import_discovery(self):
        # Checked against the import graph rather than the raw text: the
        # module's docstring legitimately explains why discovery is absent.
        import ast

        source = (Path(__file__).resolve().parent.parent
                  / "myconductor" / "core" / "pipeline.py").read_text()
        imported: list[str] = []
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.Import):
                imported += [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                imported.append(module)
                imported += [f"{module}.{alias.name}" for alias in node.names]
        offenders = [name for name in imported if "discovery" in name]
        self.assertEqual(
            offenders, [],
            "core/pipeline.py must not import the discovery workflow; "
            f"found {offenders}",
        )

    def test_pipeline_module_exposes_no_discovery_attribute(self):
        from myconductor.core import pipeline as pipeline_module
        self.assertFalse(hasattr(pipeline_module, "discovery"))

    def test_report_has_no_discovery_field(self):
        report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
        self.assertFalse(hasattr(report, "discovery"))

    def test_discovery_refuses_a_single_isolate(self):
        from myconductor.modules.discovery import (
            CohortSignal,
            DiscoveryInputError,
            run,
        )

        class NoProteome:
            def candidates(self, genes):
                return []

        signal = CohortSignal(
            cohort_id="C1", n_isolates=1, n_sites=1,
            unexplained_resistant_isolates=1, drug="bedaquiline",
            convergent_genes=["atpE"],
        )
        with self.assertRaises(DiscoveryInputError):
            run(signal, NoProteome())


class NoLaneMayInventCoverage(unittest.TestCase):
    def test_only_adapters_can_assert_coverage(self):
        v = Variant.of("katG", "S315T", depth=80)
        conductor = Myconductor(platform="illumina")
        outcome = conductor.router.route([v])
        for ev in outcome.evidence:
            self.assertFalse(
                ev.asserts_coverage,
                f"{ev.lane.value} asserted coverage; only engine adapters may",
            )


class ConflictingLocalValidationNeverPicksASide(unittest.TestCase):
    """The defect this guards against: a lane resolving a conflicting local
    laboratory result by choosing the more recent or more frequent answer,
    instead of reporting the conflict itself."""

    def test_conflicting_verdict_is_indeterminate_not_a_choice(self):
        from myconductor.core.models import Call as _Call
        from myconductor.federated.vus_feedback import (
            LocalValidationStore,
            ValidationMethod,
            VUSValidationRecord,
        )
        from myconductor.modules.local_validation import LocalValidationModule

        v = Variant.of("Rv0678", "L114R")
        store = LocalValidationStore()
        for i, result in enumerate([_Call.RESISTANT, _Call.SUSCEPTIBLE,
                                    _Call.RESISTANT, _Call.RESISTANT]):
            store.ingest(VUSValidationRecord(
                variant_key=v.key(), variant_label=v.label(), gene=v.gene,
                drug="bedaquiline", isolate_id=f"iso{i}", site_id="site-A",
                method=ValidationMethod.MIC, result=result))
        evidence = LocalValidationModule(store).evaluate(v)
        self.assertEqual(len(evidence), 1)
        # 3 of 4 records say RESISTANT; a majority-vote lane would report
        # that. This lane must not, because disagreement is a finding.
        self.assertEqual(evidence[0].call, _Call.INDETERMINATE)


class CalibratedPosteriorCannotBecomeAResistanceCall(unittest.TestCase):
    """The defect this guards against: a statistical posterior, however
    close to 1.0, being treated as strong enough evidence to assert
    RESISTANT. Only Tier.CATALOGUED/PHENOTYPIC may do that
    (DrugEvidence.__post_init__), and this lane is structurally Tier.PREDICTED."""

    def test_near_certain_posterior_still_yields_indeterminate_predicted(self):
        from myconductor.modules.calibration import (
            CalibrationTable,
            DetectionCalibration,
            FederatedPriorSource,
        )
        from myconductor.modules.heteroresistance import (
            HeteroresistanceAssessor,
            heteroresistance_evidence,
        )

        v = Variant.of("gyrA", "D94G", vaf=0.5, depth=200, alt_depth=100,
                       platform="illumina")
        table = CalibrationTable()
        table.add(DetectionCalibration(
            platform="illumina", drug="moxifloxacin", lineage="*",
            lod50=0.02, width=0.01, background_fp_rate=0.001,
            n_replicates=200, source="test"))

        class _AlmostCertain:
            n_isolates = 500
            resistant_fraction = 0.999

        prior_source = FederatedPriorSource(
            associations={(v.key(), "moxifloxacin"): _AlmostCertain()})
        assessor = HeteroresistanceAssessor(
            depth_floor=10, calibration=table, prior_source=prior_source)
        ev = DrugEvidence(drug="moxifloxacin", call=Call.RESISTANT,
                          tier=Tier.CATALOGUED, lane=Lane.CATALOGUE,
                          variant=v.identity, rationale="catalogued")
        findings = assessor.assess([v], [ev])
        self.assertGreater(findings[0].posterior_resistance_probability, 0.9)

        evidence = heteroresistance_evidence(findings, [v])
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].tier, Tier.PREDICTED)
        self.assertEqual(evidence[0].call, Call.INDETERMINATE)


if __name__ == "__main__":
    unittest.main()
