import unittest

from myconductor.core.models import Call, DrugEvidence, Lane, Tier, Variant
from myconductor.modules.calibration import (
    CalibrationTable,
    DetectionCalibration,
    FederatedPriorSource,
    NullPriorSource,
    posterior_resistance_probability,
)
from myconductor.modules.heteroresistance import (
    HeteroresistanceAssessor,
    heteroresistance_evidence,
)


def variant(**kwargs):
    kwargs.setdefault("gene", "gyrA")
    kwargs.setdefault("change", "D94G")
    return Variant.of(kwargs.pop("gene"), kwargs.pop("change"), **kwargs)


def evidence_for(v, drug="moxifloxacin"):
    return DrugEvidence(drug=drug, call=Call.RESISTANT, tier=Tier.CATALOGUED,
                        lane=Lane.CATALOGUE, variant=v.identity,
                        rationale="catalogued resistant")


def calibrated_table():
    table = CalibrationTable()
    table.add(DetectionCalibration(
        platform="illumina", drug="moxifloxacin", lineage="*",
        lod50=0.05, width=0.04, background_fp_rate=0.02, n_replicates=30,
        source="in-house dilution series"))
    return table


class AssociationStub:
    def __init__(self, n_isolates, resistant_fraction):
        self.n_isolates = n_isolates
        self.resistant_fraction = resistant_fraction


class PosteriorMathTests(unittest.TestCase):
    def test_missing_calibration_yields_no_posterior(self):
        result = posterior_resistance_probability(
            variant(vaf=0.15, platform="illumina"), "moxifloxacin",
            CalibrationTable(), NullPriorSource())
        self.assertIsNone(result)

    def test_missing_prior_yields_no_posterior(self):
        result = posterior_resistance_probability(
            variant(vaf=0.15, platform="illumina"), "moxifloxacin",
            calibrated_table(), NullPriorSource())
        self.assertIsNone(result)

    def test_configured_pair_produces_a_bounded_probability(self):
        v = variant(vaf=0.15, platform="illumina")
        prior_source = FederatedPriorSource(
            associations={(v.key(), "moxifloxacin"): AssociationStub(40, 0.9)})
        result = posterior_resistance_probability(
            v, "moxifloxacin", calibrated_table(), prior_source)
        self.assertIsNotNone(result)
        self.assertGreaterEqual(result.probability, 0.0)
        self.assertLessEqual(result.probability, 1.0)
        lo, hi = result.credible_interval
        self.assertLessEqual(lo, result.probability)
        self.assertGreaterEqual(hi, result.probability)

    def test_thin_prior_evidence_is_not_trusted(self):
        v = variant(vaf=0.15, platform="illumina")
        prior_source = FederatedPriorSource(
            associations={(v.key(), "moxifloxacin"): AssociationStub(3, 0.9)},
            min_isolates=10)
        result = posterior_resistance_probability(
            v, "moxifloxacin", calibrated_table(), prior_source)
        self.assertIsNone(result)

    def test_lineage_specific_calibration_is_preferred_over_wildcard(self):
        table = CalibrationTable()
        table.add(DetectionCalibration(
            platform="illumina", drug="moxifloxacin", lineage="*",
            lod50=0.10, width=0.04, background_fp_rate=0.05, n_replicates=10,
            source="generic"))
        table.add(DetectionCalibration(
            platform="illumina", drug="moxifloxacin", lineage="lineage4",
            lod50=0.03, width=0.02, background_fp_rate=0.01, n_replicates=50,
            source="lineage4-specific"))
        hit = table.get("illumina", "moxifloxacin", "lineage4")
        self.assertEqual(hit.source, "lineage4-specific")
        fallback = table.get("illumina", "moxifloxacin", "lineage2")
        self.assertEqual(fallback.source, "generic")


class HeteroresistanceIntegrationTests(unittest.TestCase):
    """The posterior must never let this lane cross into RESISTANT."""

    def setUp(self):
        probe = variant(vaf=0.15, depth=64, alt_depth=10, platform="illumina")
        self.prior_source = FederatedPriorSource(
            associations={(probe.key(), "moxifloxacin"): AssociationStub(50, 0.98)})
        self.assessor = HeteroresistanceAssessor(
            depth_floor=10, calibration=calibrated_table(),
            prior_source=self.prior_source)

    def test_default_assessor_computes_no_posterior(self):
        default = HeteroresistanceAssessor(depth_floor=10)
        v = variant(vaf=0.15, depth=64, alt_depth=10, platform="illumina")
        finding = default.assess([v], [evidence_for(v)])[0]
        self.assertIsNone(finding.posterior_resistance_probability)
        self.assertEqual(finding.calibration_source, "uncalibrated default")

    def test_configured_assessor_computes_a_high_posterior(self):
        v = variant(vaf=0.15, depth=64, alt_depth=10, platform="illumina")
        finding = self.assessor.assess([v], [evidence_for(v)])[0]
        self.assertIsNotNone(finding.posterior_resistance_probability)
        self.assertGreater(finding.posterior_resistance_probability, 0.5)

    def test_evidence_from_a_high_posterior_is_still_indeterminate_predicted(self):
        v = variant(vaf=0.15, depth=64, alt_depth=10, platform="illumina")
        findings = self.assessor.assess([v], [evidence_for(v)])
        evidence = heteroresistance_evidence(findings, [v])
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].call, Call.INDETERMINATE)
        self.assertEqual(evidence[0].tier, Tier.PREDICTED)
        self.assertGreater(evidence[0].confidence, 0.5)

    def test_no_evidence_emitted_when_posterior_is_absent(self):
        v = variant(vaf=0.15, depth=64, alt_depth=10, platform="illumina")
        default = HeteroresistanceAssessor(depth_floor=10)
        findings = default.assess([v], [evidence_for(v)])
        self.assertEqual(heteroresistance_evidence(findings, [v]), [])


if __name__ == "__main__":
    unittest.main()
