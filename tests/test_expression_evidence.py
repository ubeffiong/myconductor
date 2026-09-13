import unittest
from myconductor.core.models import Call, Tier, Lane, DrugEvidence, DrugResult, Variant, ExpressionEvidence
from myconductor.core.context import SampleContext
from myconductor.modules.expression_evidence import reconcile_expression


class ExpressionTests(unittest.TestCase):
    def obs(self, **kwargs):
        data = dict(gene="testGene", measurement_type="rna_seq", expression_level=10, reference_level=1,
                    fold_change=10, provenance="lab/1", sample_id="s", observation_id="o", isolate_id="i",
                    site_id="site", organism="test-organism", unit="normalized units", timestamp="2026-09-01")
        data.update(kwargs)
        return ExpressionEvidence(**data)

    def test_large_fold_change_does_not_change_call_or_confidence(self):
        variant = Variant.of("testGene", "X1Y")
        result = DrugResult("test-drug", Call.INDETERMINATE, Tier.INFERRED, confidence=None,
                            evidence=[DrugEvidence("test-drug", Call.INDETERMINATE, Tier.INFERRED,
                                      Lane.EFFLUX_REGULATORY, variant=variant.identity)])
        findings = reconcile_expression([result], [self.obs(fold_change=1e10)],
                                       SampleContext("s", "i", "site", "test-organism"))
        self.assertEqual(result.call, Call.INDETERMINATE)
        self.assertIsNone(result.confidence)
        self.assertEqual(result.evidence[-1].tier, Tier.PREDICTED)
        self.assertEqual(findings[0]["reported_conclusion"], "uncertain")

    def test_wrong_sample_and_unsupported_link_refused(self):
        with self.assertRaises(ValueError):
            reconcile_expression([], [self.obs()], SampleContext("other", "i", "site", "test-organism"))
        with self.assertRaises(ValueError):
            self.obs(linked_variant_keys=("v",))

    def test_linear_abundance_and_raw_metric_are_distinguished(self):
        with self.assertRaises(ValueError):
            self.obs(fold_change=2, measurement_scale="linear_abundance")
        self.obs(fold_change=2, measurement_scale="reported_metric")

    def test_no_change_does_not_establish_susceptibility(self):
        findings = reconcile_expression([], [self.obs(fold_change=1)], SampleContext("s", "i", "site", "test-organism"))
        self.assertEqual(findings[0]["matching_status"], "unlinked")
