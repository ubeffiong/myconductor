import unittest
from dataclasses import replace
from myconductor.core.context import SampleContext
from myconductor.core.models import Call, Tier, DrugResult
from myconductor.modules.mic_evidence import mic_call, MICPrediction, reconcile_predictions
from myconductor.io.phenotypes import PhenotypeObservation


def prediction(**kwargs):
    data = dict(prediction_id="pred1", sample_id="s", isolate_id="i", site_id="site",
                organism="test-organism", drug="test-drug", value=1, unit="mg/L",
                method="external model", model_version="1", source="model-output/1",
                timestamp="2026-09-01", critical_concentration=2, breakpoint_reference="test-SOP-v1")
    data.update(kwargs)
    return MICPrediction(**data)


class MICTests(unittest.TestCase):
    def test_censoring_and_equality_preserve_ambiguity(self):
        self.assertEqual(mic_call(1, "mg/L", 2), Call.SUSCEPTIBLE)
        self.assertEqual(mic_call(4, "mg/L", 2), Call.RESISTANT)
        self.assertIsNone(mic_call(2, "mg/L", 2))
        self.assertEqual(mic_call(2, "mg/L", 2, susceptible_inclusive=True), Call.SUSCEPTIBLE)
        self.assertEqual(mic_call(2, "mg/L", 2, "right"), Call.RESISTANT)
        self.assertIsNone(mic_call(4, "mg/L", 2, "left"))
        self.assertIsNone(mic_call(1, "mg/L", 2, "right"))
        self.assertIsNone(mic_call(2, "mg/L", 2, interval=(1, 4)))

    def test_invalid_units_bounds_and_missing_interval_rejected(self):
        for args in ((1, "mmol/L", 2), (float("nan"), "mg/L", 2),
                     (1, "mg/L", -2), (1, "mg/L", 2, "interval")):
            with self.assertRaises(ValueError):
                mic_call(*args)
        with self.assertRaises(ValueError):
            prediction(interval=(2, 1))

    def test_prediction_conflict_abstains_but_keeps_genomic_call(self):
        result = DrugResult("test-drug", Call.RESISTANT, Tier.CATALOGUED, genomic_call=Call.RESISTANT)
        findings = reconcile_predictions([result], [prediction()], SampleContext("s", "i", "site", "test-organism"))
        self.assertEqual(result.call, Call.INDETERMINATE)
        self.assertEqual(result.genomic_call, Call.RESISTANT)
        self.assertEqual(result.evidence[-1].tier, Tier.PREDICTED)
        self.assertEqual(result.evidence[-1].call, Call.INDETERMINATE)
        self.assertTrue(findings[0]["conflict"])

    def test_predictions_never_establish_resistance(self):
        for base in (Call.SUSCEPTIBLE, Call.NOT_ASSESSED, Call.INDETERMINATE):
            result = DrugResult("test-drug", base, Tier.NONE, genomic_call=base)
            reconcile_predictions([result], [prediction(value=4)], SampleContext("s", "i", "site", "test-organism"))
            self.assertNotEqual(result.call, Call.RESISTANT)

    def test_ambiguous_prediction_is_not_definite_disagreement(self):
        result = DrugResult("test-drug", Call.RESISTANT, Tier.CATALOGUED, genomic_call=Call.RESISTANT)
        findings = reconcile_predictions([result], [prediction(interval=(0.5, 4))],
                                        SampleContext("s", "i", "site", "test-organism"))
        self.assertFalse(findings[0]["conflict"])
        self.assertEqual(result.call, Call.RESISTANT)

    def test_matched_phenotype_not_overridden_by_prediction(self):
        result = DrugResult("test-drug", Call.RESISTANT, Tier.PHENOTYPIC,
                            genomic_call=Call.RESISTANT, phenotypic_call=Call.RESISTANT)
        reconcile_predictions([result], [prediction()], SampleContext("s", "i", "site", "test-organism"))
        self.assertEqual(result.call, Call.RESISTANT)
        self.assertEqual(result.tier, Tier.PHENOTYPIC)

    def test_wrong_isolate_prediction_rejected(self):
        with self.assertRaises(ValueError):
            reconcile_predictions([], [prediction()], SampleContext("s", "different", "site", "test-organism"))

    def test_measured_mic_mismatch_is_a_note_not_silent_override(self):
        obs = PhenotypeObservation("dst", "s", "i", "site", "test-organism", "test-drug", "resistant",
              "lab DST", "2026-09-01", "lab", "pass", "SOP", "1", "lab/1",
              mic=1, mic_unit="mg/L", critical_concentration=2)
        evidence = obs.evidence(SampleContext("s", "i", "site", "test-organism"))
        self.assertEqual(evidence.call, Call.RESISTANT)
        self.assertEqual(evidence.metadata["mic_comparison"], "susceptible")
        self.assertTrue(any("mismatch" in x for x in evidence.limitations))
