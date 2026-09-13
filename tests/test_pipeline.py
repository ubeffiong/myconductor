import json
import tempfile
import unittest
from pathlib import Path

from myconductor import Myconductor
from myconductor.adapters.mykrobe import MykrobeAdapter
from myconductor.core.models import Call, Lane
from myconductor.reporting.fhir import TerminologyMap, to_fhir_bundle
from myconductor.reporting.render import render_text

DATA = Path(__file__).resolve().parent.parent / "myconductor" / "data"
VCF = DATA / "example_input.vcf"
TSV = DATA / "example_input.tsv"
MASK = DATA / "example_callable.tsv"


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.report = Myconductor(platform="illumina").analyze(
            VCF, mask_path=MASK)

    def _call(self, drug):
        return self.report.result_for(drug).call

    def test_sample_id_propagates(self):
        self.assertEqual(self.report.sample_id, "MTB-DEMO-001")

    def test_catalogued_resistance_is_detected(self):
        self.assertIs(self._call("isoniazid"), Call.RESISTANT)
        self.assertIs(self._call("rifampicin"), Call.RESISTANT)
        self.assertIs(self._call("ethambutol"), Call.RESISTANT)

    def test_covered_loci_with_no_findings_are_susceptible(self):
        self.assertIs(self._call("linezolid"), Call.INDETERMINATE)
        self.assertIs(self._call("pretomanid"), Call.INDETERMINATE)

    def test_uncovered_locus_is_not_assessed_or_indeterminate(self):
        # rrs is absent from the mask and its record was dropped for depth.
        self.assertIn(self._call("amikacin"),
                      (Call.NOT_ASSESSED, Call.INDETERMINATE))
        self.assertFalse(self.report.result_for("amikacin").permits_use)

    def test_efflux_variant_withholds_susceptibility_for_both_drugs(self):
        self.assertIs(self._call("bedaquiline"), Call.INDETERMINATE)
        self.assertIs(self._call("clofazimine"), Call.INDETERMINATE)

    def test_vus_in_a_required_locus_yields_indeterminate(self):
        self.assertIs(self._call("pyrazinamide"), Call.INDETERMINATE)

    def test_minority_allele_is_detected_for_the_fluoroquinolone(self):
        detected = {h.drug for h in self.report.heteroresistance if h.detected}
        self.assertIn("moxifloxacin", detected)

    def test_low_depth_record_is_warned_about(self):
        self.assertTrue(any("rrs" in w for w in self.report.qc_warnings))

    def test_synonymous_variant_reaches_no_lane(self):
        self.assertEqual(self.report.lane_counts.get("unexamined", 0), 1)

    def test_unperformed_controls_are_enumerated(self):
        not_performed = [f.check for f in self.report.qc
                         if f.status == "not_performed"]
        self.assertIn("species_confirmation", not_performed)
        self.assertIn("contamination", not_performed)
        self.assertIn("lineage_assignment", not_performed)

    def test_provenance_records_the_coverage_source(self):
        self.assertEqual(self.report.provenance.coverage_source, "depth-table")
        self.assertIsNotNone(self.report.provenance.catalogue)

    def test_no_regimen_is_eligible_for_this_strain(self):
        self.assertFalse(self.report.eligibility.any_eligible)

    def test_vus_priorities_are_produced(self):
        labels = {v.variant_label for v in self.report.vus_priorities}
        self.assertIn("pncA_D12A", labels)

    def test_vus_priorities_carry_no_resistance_call(self):
        for v in self.report.vus_priorities:
            self.assertIn(v.priority,
                          ("high", "moderate", "low", "insufficient-data"))
            self.assertTrue(v.recommended_experiment)

    def test_null_annotator_reports_data_gaps(self):
        for v in self.report.vus_priorities:
            self.assertTrue(v.data_gaps)
            self.assertEqual(v.priority, "insufficient-data")


class MaskAbsentTests(unittest.TestCase):
    def test_tsv_input_without_a_mask_assesses_nothing_susceptible(self):
        report = Myconductor(platform="illumina").analyze(TSV)
        self.assertEqual([r.drug for r in report.drug_results if r.permits_use],
                         [])

    def test_callable_mask_qc_check_fails_when_absent(self):
        report = Myconductor(platform="illumina").analyze(VCF)
        mask_check = next(f for f in report.qc if f.check == "callable_mask")
        self.assertEqual(mask_check.status, "fail")


class EngineIntegrationTests(unittest.TestCase):
    #: Linezolid: the demo VCF carries no linezolid-relevant variant, so this
    #: isolates the engine's own coverage assertion.
    MYKROBE_CLEAN = {
        "MTB-DEMO-001": {
            "mykrobe_version": "0.13.0",
            "susceptibility": {"Linezolid": {"predict": "S", "called_by": {}}},
        }
    }
    #: Amikacin: the demo VCF carries an eis promoter variant, so local
    #: evidence withholds susceptibility whatever the engine concluded.
    MYKROBE_CONFLICT = {
        "MTB-DEMO-001": {
            "mykrobe_version": "0.13.0",
            "susceptibility": {"Amikacin": {"predict": "S", "called_by": {}}},
        }
    }

    def _analyze(self, fixture):
        path = Path(tempfile.mkdtemp()) / "m.json"
        path.write_text(json.dumps(fixture))
        return Myconductor(platform="illumina").analyze(
            VCF, engine_reports=[MykrobeAdapter().parse(path)])

    def test_engine_coverage_assertion_licenses_susceptibility(self):
        # No mask is supplied, so Mykrobe's own coverage handling is the only
        # thing that can license this call -- and it is attributed.
        report = self._analyze(self.MYKROBE_CLEAN)
        linezolid = report.result_for("linezolid")
        self.assertIs(linezolid.call, Call.INDETERMINATE)
        self.assertIn("validated", linezolid.reason)

    def test_local_withholding_is_not_overridden_by_an_engine(self):
        # Myconductor found a regulatory variant in eis. Mykrobe's panel may
        # not cover it, so an engine's susceptible call must not override the
        # local decision to withhold susceptibility.
        report = self._analyze(self.MYKROBE_CONFLICT)
        amikacin = report.result_for("amikacin")
        self.assertIs(amikacin.call, Call.INDETERMINATE)
        self.assertFalse(amikacin.permits_use)

    def test_conflicting_engine_evidence_is_still_retained(self):
        report = self._analyze(self.MYKROBE_CONFLICT)
        amikacin = report.result_for("amikacin")
        lanes = {ev.lane for ev in amikacin.evidence}
        self.assertIn(Lane.ENGINE, lanes,
                      "the engine's opinion must be visible even when it did "
                      "not prevail")

    def test_engine_evidence_is_tagged_with_the_engine_lane(self):
        report = self._analyze(self.MYKROBE_CLEAN)
        lanes = {ev.lane for r in report.drug_results for ev in r.evidence}
        self.assertIn(Lane.ENGINE, lanes)

    def test_engine_appears_in_provenance(self):
        report = self._analyze(self.MYKROBE_CLEAN)
        self.assertTrue(any(e.name == "mykrobe"
                            for e in report.provenance.engines))


class RenderTests(unittest.TestCase):
    def setUp(self):
        self.report = Myconductor(platform="illumina").analyze(
            VCF, mask_path=MASK)
        self.text = render_text(self.report)

    def test_header_and_disclaimer_are_present(self):
        self.assertIn("MYCONDUCTOR REPORT", self.text)
        self.assertIn("NOT a clinical device", self.text)

    def test_every_drug_appears(self):
        for r in self.report.drug_results:
            self.assertIn(r.drug, self.text)

    def test_non_verdicts_state_why(self):
        self.assertIn("why not:", self.text)

    def test_eligibility_is_not_phrased_as_a_prescription(self):
        self.assertIn("GUIDELINE ELIGIBILITY (not a prescription)", self.text)
        self.assertNotIn("Proposed regimen", self.text)

    def test_unperformed_controls_are_shown(self):
        self.assertIn("NOT PERFORMED", self.text)
        self.assertIn("absent, not passed", self.text)

    def test_susceptibility_gate_is_explained(self):
        self.assertIn("Only a validated genomic S or matched laboratory S counts", self.text)

    def test_no_demo_banner_by_default(self):
        self.assertNotIn("SYNTHETIC DEMONSTRATION", self.text)


class FHIRTests(unittest.TestCase):
    def setUp(self):
        self.report = Myconductor(platform="illumina").analyze(
            VCF, mask_path=MASK)
        self.bundle = to_fhir_bundle(self.report)
        self.observations = [
            e["resource"] for e in self.bundle["entry"]
            if e["resource"]["resourceType"] == "Observation"
        ]

    def test_bundle_shape(self):
        self.assertEqual(self.bundle["resourceType"], "Bundle")
        types = {e["resource"]["resourceType"] for e in self.bundle["entry"]}
        self.assertTrue({"Specimen", "DiagnosticReport", "Observation",
                         "Provenance"} <= types)

    def test_diagnostic_report_references_its_results(self):
        dr = next(e["resource"] for e in self.bundle["entry"]
                  if e["resource"]["resourceType"] == "DiagnosticReport")
        self.assertTrue(dr["result"])
        ids = {o["id"] for o in self.observations}
        for ref in dr["result"]:
            self.assertIn(ref["reference"].removeprefix("urn:uuid:"), ids)

    def test_not_assessed_uses_data_absent_reason_not_a_value(self):
        for r in self.report.drug_results:
            if r.call is not Call.NOT_ASSESSED:
                continue
            obs = next(o for o in self.observations
                       if o.get("code", {}).get("text", "").startswith(r.drug))
            self.assertIn("dataAbsentReason", obs)
            self.assertNotIn("valueCodeableConcept", obs)
            code = obs["dataAbsentReason"]["coding"][0]["code"]
            self.assertEqual(code, "not-performed")

    def test_resistant_uses_a_standard_interpretation_code(self):
        rif = next(o for o in self.observations
                   if o.get("code", {}).get("text", "").startswith("rifampicin"))
        coding = rif["valueCodeableConcept"]["coding"][0]
        self.assertEqual(coding["code"], "R")
        self.assertIn("v3-ObservationInterpretation", coding["system"])

    def test_every_observation_is_preliminary(self):
        for obs in self.observations:
            self.assertEqual(obs["status"], "preliminary")

    def test_permits_regimen_use_is_exposed(self):
        lzd = next(o for o in self.observations
                   if o.get("code", {}).get("text", "").startswith("linezolid"))
        flag = next(x for x in lzd["extension"]
                    if x["url"].endswith("/permits-regimen-use"))
        self.assertFalse(flag["valueBoolean"])

    def test_no_loinc_codes_are_invented(self):
        for obs in self.observations:
            self.assertNotIn("coding", obs.get("code", {}))

    def test_supplied_terminology_is_used(self):
        bundle = to_fhir_bundle(self.report, TerminologyMap(
            drug_loinc={"rifampicin": ("20388-3", "Rifampin susceptibility")}))
        rif = next(e["resource"] for e in bundle["entry"]
                   if e["resource"]["resourceType"] == "Observation"
                   and e["resource"]["code"]["text"].startswith("rifampicin"))
        self.assertEqual(rif["code"]["coding"][0]["code"], "20388-3")

    def test_conformance_gaps_travel_with_the_bundle(self):
        dr = next(e["resource"] for e in self.bundle["entry"]
                  if e["resource"]["resourceType"] == "DiagnosticReport")
        notes = dr["conclusion"]
        self.assertIn("Not validated against any FHIR implementation guide",
                      notes)
        self.assertIn("LOINC", notes)

    def test_unperformed_qc_is_carried_in_the_bundle(self):
        qc = next((e["resource"] for e in self.bundle["entry"]
                   if e["resource"].get("code", {}).get("text") == "Analytical quality control summary"), None)
        self.assertIsNotNone(qc)
        notes = " ".join(n["text"] for n in qc["note"])
        self.assertIn("NOT PERFORMED", notes)

    def test_bundle_is_json_serialisable(self):
        self.assertTrue(json.dumps(self.bundle))


if __name__ == "__main__":
    unittest.main()
