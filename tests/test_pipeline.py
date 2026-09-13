import unittest
from pathlib import Path

from myconductor import Myconductor
from myconductor.reporting.fhir import to_fhir_bundle
from myconductor.reporting.render import render_text

DATA = Path(__file__).resolve().parent.parent / "myconductor" / "data"


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.report = Myconductor().analyze(DATA / "example_input.vcf")

    def test_sample_id_propagates(self):
        self.assertEqual(self.report.sample_id, "MTB-DEMO-001")

    def test_resistance_detected_for_first_line(self):
        calls = {r.drug: r.call for r in self.report.drug_results}
        self.assertTrue(calls["isoniazid"].is_resistant)
        self.assertTrue(calls["rifampicin"].is_resistant)

    def test_heteroresistance_flagged_for_minority_variant(self):
        drugs = {h.drug for h in self.report.heteroresistance}
        # gyrA D94G is present at 15% VAF in the example.
        self.assertIn("moxifloxacin", drugs)

    def test_low_depth_locus_dropped(self):
        self.assertTrue(any("rrs" in w for w in self.report.qc_warnings))

    def test_discovery_triggers_when_regimen_inadequate(self):
        self.assertFalse(self.report.regimen.adequate)
        self.assertTrue(self.report.discovery.triggered)
        self.assertTrue(self.report.discovery.targets)

    def test_render_and_fhir_produce_output(self):
        self.assertIn("MYCONDUCTOR REPORT", render_text(self.report))
        bundle = to_fhir_bundle(self.report)
        self.assertEqual(bundle["resourceType"], "Bundle")

    def test_predicted_calls_are_marked_preliminary_in_fhir(self):
        bundle = to_fhir_bundle(self.report)
        obs = [e["resource"] for e in bundle["entry"]
               if e["resource"]["resourceType"] == "Observation"]
        predicted = [o for o in obs
                     if any(x["url"] == "predicted" and x["valueBoolean"]
                            for x in o["extension"])]
        for o in predicted:
            self.assertEqual(o["status"], "preliminary")


if __name__ == "__main__":
    unittest.main()
