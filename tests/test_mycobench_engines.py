"""TB-Profiler parsers.

``COLLATE_FIXTURE`` is a verbatim excerpt of ``tests/example_collate.txt`` from
the TB-Profiler repository — a real artifact, not one composed for this test.
That is what makes ``parse_collate`` a validated parser rather than a
schema-derived guess, and these assertions pin its exact value format.
"""
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.engines import (
    COLLATE_DRUGS,
    EngineError,
    parse_collate,
    parse_results_json,
    to_variant_tsv,
)

#: Verbatim from jodyphelan/TBProfiler tests/example_collate.txt.
COLLATE_FIXTURE = (
    "sample\tmain_lineage\tsub_lineage\tspoligotype\tdrtype\t"
    "target_median_depth\tpct_reads_mapped\tnum_reads_mapped\t"
    "num_dr_variants\tnum_other_variants\trifampicin\trifapentine\t"
    "isoniazid\tethambutol\tpyrazinamide\tmoxifloxacin\tlevofloxacin\t"
    "bedaquiline\tdelamanid\tpretomanid\tlinezolid\tstreptomycin\tamikacin\t"
    "kanamycin\tcapreomycin\tclofazimine\tethionamide\tprothionamide\n"
    "por5A_bcftools\tlineage4\tlineage4.3.4.2\t-\tMDR-TB\t59.0\t99.9\t54238\t"
    "4\t24\trpoB p.Ser450Leu (1.00)\trpoB p.Ser450Leu (1.00)\t"
    "inhA c.-777C>T (1.00)\tembB p.Met306Val (1.00)\t"
    "pncA p.Val125Gly (1.00)\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t"
    "inhA c.-777C>T (1.00)\tinhA c.-777C>T (1.00)\n"
)


def write(name, text):
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text, encoding="utf-8")
    return path


class CollateTests(unittest.TestCase):
    def setUp(self):
        self.samples = parse_collate(write("collate.txt", COLLATE_FIXTURE))
        self.sample = self.samples[0]

    def test_parses_the_real_fixture(self):
        self.assertEqual(len(self.samples), 1)
        self.assertEqual(self.sample.sample, "por5A_bcftools")
        self.assertEqual(self.sample.main_lineage, "lineage4")
        self.assertEqual(self.sample.sub_lineage, "lineage4.3.4.2")
        self.assertEqual(self.sample.drtype, "MDR-TB")
        self.assertAlmostEqual(self.sample.median_depth, 59.0)
        self.assertAlmostEqual(self.sample.pct_reads_mapped, 99.9)

    def test_parses_the_variant_value_format(self):
        rifampicin = self.sample.findings["rifampicin"]
        self.assertEqual(len(rifampicin.variants), 1)
        variant = rifampicin.variants[0]
        self.assertEqual(variant.gene, "rpoB")
        self.assertEqual(variant.change, "p.Ser450Leu")
        self.assertAlmostEqual(variant.frequency, 1.0)
        self.assertEqual(variant.label(), "rpoB_p.Ser450Leu")

    def test_dash_means_no_variant(self):
        self.assertEqual(self.sample.findings["moxifloxacin"].variants, [])

    def test_parser_is_marked_validated(self):
        self.assertEqual(self.sample.parser, "collate")
        self.assertTrue(self.sample.validated_parser)

    def test_uncovered_drugs_are_recorded_not_dropped(self):
        # ethionamide has a variant but is not in the mycobacterial profile.
        self.assertIn("ethionamide", self.sample.uncompared_drugs)
        self.assertNotIn("ethionamide", self.sample.findings)

    def test_unmapped_drug_columns_map_to_none(self):
        self.assertIsNone(COLLATE_DRUGS["levofloxacin"])
        self.assertEqual(COLLATE_DRUGS["rifampicin"], "rifampicin")

    def test_missing_column_raises_rather_than_guessing(self):
        with self.assertRaises(EngineError):
            parse_collate(write("bad.txt", "foo\tbar\n1\t2\n"))

    def test_empty_file_raises(self):
        with self.assertRaises(EngineError):
            parse_collate(write("empty.txt",
                                "sample\tmain_lineage\tdrtype\n"))


class NoSusceptibleCallTests(unittest.TestCase):
    """TB-Profiler's '-' is not a coverage-backed susceptible call."""

    def setUp(self):
        self.sample = parse_collate(write("c.txt", COLLATE_FIXTURE))[0]

    def test_absence_of_a_variant_is_not_assessed_not_susceptible(self):
        calls = self.sample.resistance_calls()
        self.assertEqual(calls["moxifloxacin"], "not_assessed")
        self.assertNotIn("susceptible", set(calls.values()))

    def test_detected_resistance_is_reported_as_resistant(self):
        self.assertEqual(self.sample.resistance_calls()["rifampicin"],
                         "resistant")

    def test_coverage_rows_leave_the_fraction_empty(self):
        # A genome-wide median depth says nothing about any locus's breadth,
        # so the mask must not receive a fabricated callable fraction.
        rows = self.sample.coverage_rows(["rpoB", "katG"])
        self.assertEqual(len(rows), 2)
        for row in rows:
            self.assertEqual(row["callable_fraction"], "")
            self.assertEqual(row["mean_depth"], 59)


class VariantExportTests(unittest.TestCase):
    def test_variants_export_in_the_adapter_tsv_shape(self):
        sample = parse_collate(write("c.txt", COLLATE_FIXTURE))[0]
        rows = to_variant_tsv(sample)
        labels = {f"{r['gene']}_{r['change']}" for r in rows}
        self.assertIn("rpoB_p.Ser450Leu", labels)
        self.assertIn("inhA_c.-777C>T", labels)
        for row in rows:
            self.assertEqual(row["platform"], "illumina")

    def test_variants_are_deduplicated_across_drugs(self):
        # rpoB p.Ser450Leu appears under both rifampicin and rifapentine.
        sample = parse_collate(write("c.txt", COLLATE_FIXTURE))[0]
        labels = [v.label() for v in sample.all_variants()]
        self.assertEqual(len(labels), len(set(labels)))


class ResultsJsonTests(unittest.TestCase):
    RESULTS = {
        "id": "SAMPLE-1",
        "main_lin": "lineage4",
        "sublin": "lineage4.3.4.2",
        "dr_variants": [
            {"gene": "rpoB", "change": "p.Ser450Leu", "freq": 0.98,
             "drugs": [{"drug": "rifampicin", "confidence": "Assoc w R"}]},
        ],
        "qc": {"median_coverage": 65, "percent_reads_mapped": 98.2},
    }

    def test_parses_and_marks_itself_unvalidated(self):
        path = write("r.json", json.dumps(self.RESULTS))
        sample = parse_results_json(path)
        self.assertEqual(sample.sample, "SAMPLE-1")
        self.assertEqual(sample.parser, "results-json")
        self.assertFalse(sample.validated_parser)
        self.assertTrue(any("not been validated" in w for w in sample.warnings))

    def test_drug_with_no_variant_is_still_present_as_not_assessed(self):
        sample = parse_results_json(write("r.json", json.dumps(self.RESULTS)))
        calls = sample.resistance_calls()
        self.assertEqual(calls["rifampicin"], "resistant")
        self.assertEqual(calls["isoniazid"], "not_assessed")

    def test_unknown_schema_raises(self):
        with self.assertRaises(EngineError):
            parse_results_json(write("r.json", json.dumps({"foo": 1})))

    def test_invalid_json_raises(self):
        with self.assertRaises(EngineError):
            parse_results_json(write("r.json", "not json"))


if __name__ == "__main__":
    unittest.main()
