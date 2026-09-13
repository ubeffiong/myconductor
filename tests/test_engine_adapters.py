import json
import tempfile
import unittest
from pathlib import Path

from myconductor.adapters import concordance
from myconductor.adapters.amrfinderplus import AMRFinderPlusAdapter
from myconductor.adapters.base import AdapterSchemaError, EngineReport
from myconductor.adapters.mykrobe import MykrobeAdapter
from myconductor.adapters.tbprofiler import TBProfilerAdapter
from myconductor.adapters.who_catalogue import ingest_csv
from myconductor.core.models import Call, DrugEvidence, EngineRef, Lane, Tier


def write(name, text):
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return path


class TBProfilerTests(unittest.TestCase):
    RESULTS = {
        "id": "SAMPLE-1",
        "tbprofiler_version": "6.2.0",
        "main_lin": "lineage4",
        "sublin": "lineage4.3.4.2",
        "dr_variants": [
            {"gene": "rpoB", "change": "p.Ser450Leu", "freq": 0.98,
             "genome_pos": 761155, "ref": "C", "alt": "T",
             "type": "missense_variant",
             "drugs": [{"drug": "rifampicin", "confidence": "Assoc w R"}]},
            {"gene": "pncA", "change": "p.Asp12Ala", "freq": 1.0,
             "type": "missense_variant",
             "drugs": [{"drug": "pyrazinamide",
                        "confidence": "Uncertain significance"}]},
        ],
        "other_variants": [
            {"gene": "rpoC", "change": "p.Ala542Ala",
             "type": "synonymous_variant"},
        ],
        "qc": {
            "median_coverage": 65,
            "percent_reads_mapped": 98.2,
            "gene_coverage": [
                {"gene": "rpoB", "fraction": 99.5, "median_depth": 70},
                {"gene": "pncA", "fraction": 97.0, "median_depth": 62},
            ],
        },
    }

    def test_parses_resistance_and_lineage(self):
        path = write("r.json", json.dumps(self.RESULTS))
        report = TBProfilerAdapter().parse(path)
        self.assertEqual(report.sample_id, "SAMPLE-1")
        self.assertEqual(report.lineage, "lineage4.3.4.2")
        self.assertEqual(report.engine.version, "6.2.0")

    def test_graded_resistance_becomes_a_catalogued_resistant_call(self):
        report = TBProfilerAdapter().parse(write("r.json", json.dumps(self.RESULTS)))
        rif = [ev for ev in report.evidence if ev.drug == "rifampicin"]
        self.assertEqual(len(rif), 1)
        self.assertIs(rif[0].call, Call.RESISTANT)
        self.assertIs(rif[0].tier, Tier.CATALOGUED)

    def test_uncertain_grading_becomes_indeterminate(self):
        report = TBProfilerAdapter().parse(write("r.json", json.dumps(self.RESULTS)))
        pza = [ev for ev in report.evidence if ev.drug == "pyrazinamide"]
        self.assertIs(pza[0].call, Call.INDETERMINATE)

    def test_gene_coverage_becomes_a_callable_mask(self):
        report = TBProfilerAdapter().parse(write("r.json", json.dumps(self.RESULTS)))
        self.assertIsNotNone(report.mask)
        ok, _, reasons = report.mask.assess(
            ["rpoB"], depth_floor=10, fraction_floor=0.95)
        self.assertTrue(ok, reasons)

    def test_percentages_are_converted_to_fractions(self):
        report = TBProfilerAdapter().parse(write("r.json", json.dumps(self.RESULTS)))
        cov = report.mask.coverage_for("rpoB")
        self.assertAlmostEqual(cov.callable_fraction, 0.995)

    def test_minority_allele_is_recorded_as_a_limitation(self):
        data = json.loads(json.dumps(self.RESULTS))
        data["dr_variants"][0]["freq"] = 0.2
        report = TBProfilerAdapter().parse(write("r.json", json.dumps(data)))
        rif = [ev for ev in report.evidence if ev.drug == "rifampicin"][0]
        self.assertTrue(any("allele fraction" in l for l in rif.limitations))

    def test_unknown_schema_raises_rather_than_mis_parsing(self):
        with self.assertRaises(AdapterSchemaError) as ctx:
            TBProfilerAdapter().parse(write("r.json", json.dumps({"foo": []})))
        self.assertIn("dr_variants", str(ctx.exception))

    def test_invalid_json_raises(self):
        with self.assertRaises(AdapterSchemaError):
            TBProfilerAdapter().parse(write("r.json", "not json"))


class MykrobeTests(unittest.TestCase):
    RESULTS = {
        "SAMPLE-1": {
            "mykrobe_version": "0.13.0",
            "phylogenetics": {
                "species": {"Mycobacterium_tuberculosis": {}},
                "lineage": {"lineage": ["lineage4", "lineage4.3"]},
            },
            "susceptibility": {
                "Isoniazid": {"predict": "R",
                              "called_by": {"katG_S315T-GC2155168GA": {}}},
                "Linezolid": {"predict": "S", "called_by": {}},
                "Amikacin": {"predict": "N", "called_by": {}},
                "Bedaquiline": {"predict": "U", "called_by": {}},
                "Moxifloxacin": {"predict": "r",
                                 "called_by": {"gyrA_D94G-AG7570GA": {}}},
            },
        }
    }

    def setUp(self):
        self.report = MykrobeAdapter().parse(
            write("m.json", json.dumps(self.RESULTS)))

    def _call(self, drug):
        return [ev for ev in self.report.evidence if ev.drug == drug][0]

    def test_species_and_lineage_are_captured(self):
        self.assertEqual(self.report.species, "Mycobacterium_tuberculosis")
        self.assertEqual(self.report.lineage, "lineage4.3")

    def test_resistant_call_is_imported(self):
        self.assertIs(self._call("isoniazid").call, Call.RESISTANT)

    def test_susceptible_call_asserts_coverage(self):
        lzd = self._call("linezolid")
        self.assertIs(lzd.call, Call.SUSCEPTIBLE)
        self.assertTrue(lzd.asserts_coverage)

    def test_no_coverage_becomes_not_assessed_and_asserts_nothing(self):
        amk = self._call("amikacin")
        self.assertIs(amk.call, Call.NOT_ASSESSED)
        self.assertFalse(amk.asserts_coverage)

    def test_unknown_becomes_indeterminate(self):
        self.assertIs(self._call("bedaquiline").call, Call.INDETERMINATE)

    def test_minority_call_is_resistant_with_a_limitation(self):
        mfx = self._call("moxifloxacin")
        self.assertIs(mfx.call, Call.RESISTANT)
        self.assertTrue(any("minority" in l for l in mfx.limitations))

    def test_called_by_key_is_parsed_into_a_variant(self):
        inh = self._call("isoniazid")
        self.assertEqual(inh.variant_label, "katG_S315T")

    def test_unrecognised_predict_code_is_skipped_with_a_warning(self):
        data = {"S": {"susceptibility": {"Isoniazid": {"predict": "Z"}}}}
        report = MykrobeAdapter().parse(write("m.json", json.dumps(data)))
        self.assertEqual(report.evidence, [])
        self.assertTrue(any("unrecognised" in w for w in report.warnings))

    def test_missing_susceptibility_raises(self):
        with self.assertRaises(AdapterSchemaError):
            MykrobeAdapter().parse(
                write("m.json", json.dumps({"S": {"phylogenetics": {}}})))


class AMRFinderPlusTests(unittest.TestCase):
    TSV = (
        "Gene symbol\tSequence name\tScope\tElement type\tElement subtype\t"
        "Class\tSubclass\tMethod\t% Coverage of reference sequence\t"
        "% Identity to reference sequence\tContig id\tStart\n"
        "blaKPC-2\tKPC family carbapenemase\tcore\tAMR\tAMR\tBETA-LACTAM\t"
        "CARBAPENEM\tALLELEX\t100.00\t100.00\tcontig1\t1000\n"
        "aac(6')-Ib\taminoglycoside acetyltransferase\tcore\tAMR\tAMR\t"
        "AMINOGLYCOSIDE\tAMIKACIN\tBLASTX\t65.00\t99.00\tcontig2\t2000\n"
        "fimH\ttype 1 fimbrin\tplus\tVIRULENCE\tVIRULENCE\t\t\tBLASTX\t"
        "100.00\t99.00\tcontig3\t3000\n"
    )

    def setUp(self):
        self.report = AMRFinderPlusAdapter().parse(write("a.tsv", self.TSV))

    def test_only_amr_elements_become_evidence(self):
        self.assertEqual(len(self.report.evidence), 2)

    def test_drug_is_the_class_not_an_agent(self):
        drugs = {ev.drug for ev in self.report.evidence}
        self.assertEqual(drugs, {"carbapenem", "amikacin"})
        for ev in self.report.evidence:
            self.assertTrue(any("drug class" in l for l in ev.limitations))

    def test_full_length_high_identity_hit_is_resistant(self):
        kpc = [ev for ev in self.report.evidence if ev.drug == "carbapenem"][0]
        self.assertIs(kpc.call, Call.RESISTANT)

    def test_partial_coverage_downgrades_to_indeterminate(self):
        partial = [ev for ev in self.report.evidence if ev.drug == "amikacin"][0]
        self.assertIs(partial.call, Call.INDETERMINATE)
        self.assertTrue(any("coverage" in l for l in partial.limitations))

    def test_missing_required_column_raises(self):
        with self.assertRaises(AdapterSchemaError):
            AMRFinderPlusAdapter().parse(write("a.tsv", "Foo\tBar\n1\t2\n"))


class WHOCatalogueIngestTests(unittest.TestCase):
    CSV = (
        "variant,drug,FINAL CONFIDENCE GRADING,Genome position\n"
        "rpoB_S450L,Rifampicin,1) Assoc w R,761155\n"
        "Rv0678_R94Q,Bedaquiline,2) Assoc w R - interim,4247429\n"
        "Rv0678_R94Q,Clofazimine,2) Assoc w R - interim,4247429\n"
        "katG_R463L,Isoniazid,5) Not assoc w R,2154724\n"
        "pncA_D12A,Pyrazinamide,3) Uncertain significance,2288725\n"
        "bad_row,,,\n"
    )

    def setUp(self):
        self.result = ingest_csv(write("who.csv", self.CSV),
                                 catalogue_version="test-v2")

    def test_entries_are_merged_per_variant(self):
        rv0678 = [e for e in self.result.entries if e["gene"] == "Rv0678"]
        self.assertEqual(len(rv0678), 1)
        self.assertEqual(sorted(rv0678[0]["drugs"]),
                         ["bedaquiline", "clofazimine"])

    def test_grading_maps_to_calls(self):
        by_gene = {e["gene"]: e for e in self.result.entries}
        self.assertEqual(by_gene["rpoB"]["call"], "resistant")
        self.assertEqual(by_gene["katG"]["call"], "not_associated")
        self.assertEqual(by_gene["pncA"]["call"], "indeterminate")

    def test_confidence_is_labelled_as_ordinal_not_probability(self):
        entry = self.result.entries[0]
        self.assertIn("not a probability", entry["confidence_note"])

    def test_unusable_rows_are_skipped_and_reported(self):
        self.assertTrue(self.result.skipped)

    def test_bundled_json_is_not_marked_illustrative(self):
        bundle = self.result.to_bundled_json()
        self.assertFalse(bundle["illustrative"])
        self.assertEqual(bundle["catalogue_version"], "test-v2")

    def test_missing_required_column_raises_with_guidance(self):
        with self.assertRaises(AdapterSchemaError) as ctx:
            ingest_csv(write("who.csv", "foo,bar\n1,2\n"))
        self.assertIn("column", str(ctx.exception))


class ConcordanceTests(unittest.TestCase):
    def _report(self, name, drug, call):
        engine = EngineRef(name=name, version="1")
        return EngineReport(engine=engine, evidence=[DrugEvidence(
            drug=drug, call=call,
            tier=Tier.CATALOGUED if call is Call.RESISTANT else Tier.NONE,
            lane=Lane.ENGINE, engine=engine, rationale="test",
        )])

    def test_agreement_is_reported(self):
        summary = concordance([
            self._report("a", "rifampicin", Call.RESISTANT),
            self._report("b", "rifampicin", Call.RESISTANT),
        ])
        self.assertEqual(summary.agreed, {"rifampicin": "resistant"})
        self.assertEqual(summary.agreement_rate, 1.0)

    def test_disagreement_is_reported_not_resolved(self):
        summary = concordance([
            self._report("a", "rifampicin", Call.RESISTANT),
            self._report("b", "rifampicin", Call.INDETERMINATE),
        ])
        self.assertEqual(len(summary.disagreed), 1)
        self.assertIn("Not resolved by vote", summary.disagreed[0].note)

    def test_single_source_is_flagged_as_such(self):
        summary = concordance([self._report("a", "rifampicin", Call.RESISTANT)])
        self.assertIn("rifampicin", summary.single_source)
        self.assertIsNone(summary.agreement_rate)

    def test_nothing_to_compare_is_stated_plainly(self):
        summary = concordance([])
        self.assertIn("nothing to reconcile", summary.describe())


if __name__ == "__main__":
    unittest.main()
