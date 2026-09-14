"""The HTML report: self-contained, and never showing a number it does not have.

Two properties carry this file.

**Self-contained.** The output must open on a machine that has never had a
network — no CDN, no webfont, no ``<script src>``. This tool is aimed at
settings where a report that silently fails to fetch a charting library is a
report that cannot be read.

**Nothing generated.** The design this was ported from hard-coded seventeen
datasets and manufactured four more with a seeded PRNG: call matrices, coverage
grids, lineage strata, sequence bases. Rendered against a one-sample run it
showed "326 / 400 isolates" and "17 very major errors" directly beneath a
header that correctly said one isolate — two numbers on one screen, disagreeing,
both looking authoritative. Every such figure is now derived, and a panel with
no data renders an empty state. These tests hold that line, because it is the
kind of thing that creeps back one convenient default at a time.
"""
from __future__ import annotations

import json
import re
import unittest

from myconductor.core.pipeline import Myconductor
from myconductor.reporting import html_data as shape
from myconductor.reporting.html_report import render_html

INPUT = "myconductor/data/example_input.vcf"
MASK = "myconductor/data/example_callable.tsv"


def demo_report():
    return Myconductor(platform="illumina", demo_mode=True).analyze(
        INPUT, mask_path=MASK)


def payload_of(html: str) -> dict:
    match = re.search(r"window\.__MYCONDUCTOR__ = (.*?);</script>", html, re.S)
    assert match, "no injected payload"
    return json.loads(match.group(1).replace("\\u003c", "<"))


class SelfContainedTests(unittest.TestCase):
    def setUp(self):
        self.html = render_html(demo_report())

    def test_no_external_resource_is_referenced(self):
        for pattern in (r'src="http', r'href="http', r'src="//', r'href="//',
                        r"@import", r"cdnjs", r"jsdelivr", r"unpkg",
                        r"fonts\.googleapis", r"fonts\.gstatic"):
            self.assertNotRegex(self.html, pattern, pattern)

    def test_every_script_is_inline(self):
        for tag in re.findall(r"<script[^>]*>", self.html):
            self.assertNotIn("src=", tag)

    def test_every_stylesheet_is_inline(self):
        self.assertNotIn("<link", self.html)
        self.assertIn("<style>", self.html)

    def test_it_is_one_file(self):
        self.assertTrue(self.html.startswith("<!doctype html>"))
        self.assertIn("</html>", self.html)

    def test_no_template_placeholder_survives(self):
        leftovers = {t for t in re.findall(r"__[A-Z_]+__", self.html)
                     if t != "__MYCONDUCTOR__"}
        self.assertEqual(leftovers, set())


class NothingFabricatedTests(unittest.TestCase):
    """The mock's numbers must not survive anywhere."""

    def setUp(self):
        self.html = render_html(demo_report())

    def test_none_of_the_mock_figures_appear(self):
        for ghost in ("326 / 400", "81.5%", "4,800", "n = 128 lab",
                      "42 variants", "k-anonymity threshold = 5",
                      "mbx-2026-09-14-0042", "1,847"):
            self.assertNotIn(ghost, self.html, ghost)

    def test_the_runtime_carries_no_random_generator(self):
        from pathlib import Path
        import myconductor.reporting.html_report as module
        runtime = (Path(module.__file__).parent / "assets" / "report.js"
                   ).read_text(encoding="utf-8")
        self.assertNotIn("Math.random", runtime)
        # The seeded PRNG that manufactured every heatmap and sequence.
        self.assertNotRegex(runtime, r"\brng\(")

    def test_absent_measurements_render_as_a_dash_not_a_number(self):
        """A null coverage once interpolated as the literal text "null%".

        Asserted against the markup with the inline stylesheet and runtime
        stripped out: the runtime's own source legitimately discusses the
        defect it guards against, and matching that comment would make this
        test pass or fail on prose rather than on output.
        """
        markup = re.sub(r"<(script|style)\b.*?</\1>", "", self.html, flags=re.S)
        for bad in ("null%", ">null<", "NaN", "undefined"):
            self.assertNotIn(bad, markup, bad)

    def test_a_single_sample_run_claims_no_accuracy(self):
        html = self.html
        self.assertIn("No phenotypic reference was supplied", html)
        payload = payload_of(html)
        for drug in payload["drugs"]:
            self.assertFalse(drug["estimable"], drug["id"])
            self.assertIsNone(drug["error"])
            self.assertIsNone(drug["sensitivity"])

    def test_sections_with_no_data_are_empty_not_invented(self):
        payload = payload_of(self.html)
        for key in ("audit", "federated_sites", "validation_outcomes",
                    "validation_timeline", "error_trend", "lineages"):
            self.assertEqual(payload[key], [], key)


class RealDataTests(unittest.TestCase):
    def setUp(self):
        self.html = render_html(demo_report())
        self.payload = payload_of(self.html)

    def test_the_call_distribution_matches_the_run(self):
        total = sum(entry["v"] for entry in self.payload["call_distribution"])
        self.assertEqual(total, len(self.payload["drugs"]))

    def test_the_alignment_shows_only_resolved_coordinates(self):
        for locus in self.payload["alignment_loci"].values():
            self.assertTrue(locus["positions"])
            for record in locus["positions"]:
                self.assertIsInstance(record["pos"], int)
                self.assertTrue(record["ref"])
                self.assertTrue(record["alt"])

    def test_vus_carry_the_workbench_band_not_an_invented_score(self):
        for item in self.payload["vus"]:
            self.assertIn(item["priority"],
                          ("high", "moderate", "low", "insufficient-data"))

    def test_coverage_loci_come_from_the_mask(self):
        self.assertTrue(self.payload["coverage_loci"])
        for locus in self.payload["coverage_loci"]:
            self.assertTrue(locus["id"])


class PreservedContentTests(unittest.TestCase):
    """The previous report's sections must all still be here."""

    def setUp(self):
        self.html = render_html(demo_report())

    def test_the_per_drug_review_table_survives(self):
        self.assertIn('id="review"', self.html)
        for header in ("Conclusion", "Genomic evidence", "Phenotype",
                       "Assay adequacy", "Reason"):
            self.assertIn(header, self.html, header)

    def test_evidence_chains_are_still_expandable(self):
        self.assertIn("<details>", self.html)
        self.assertIn("Evidence chain", self.html)

    def test_the_investigation_queue_survives(self):
        self.assertIn('id="followup"', self.html)

    def test_quality_control_survives(self):
        self.assertIn('id="qc"', self.html)

    def test_provenance_survives(self):
        self.assertIn("Analysis fingerprint", self.html)
        self.assertIn("Tool version", self.html)

    def test_the_research_scaffold_warning_survives(self):
        self.assertIn("not for clinical use", self.html.lower())

    def test_demo_mode_is_stamped(self):
        self.assertIn("DEMO MODE", self.html.upper())

    def test_a_non_demo_run_is_not_stamped_as_demo(self):
        html = render_html(Myconductor(platform="illumina").analyze(
            INPUT, mask_path=MASK))
        self.assertNotIn("DEMO MODE", html.upper())


class CohortInputTests(unittest.TestCase):
    """Optional cohort context flows through; omitting it never fabricates."""

    BASELINE = [{
        "drug": "rifampicin", "coverage": "1.0000", "error_rate": "0.0960",
        "estimable": "yes", "n_isolates": "391", "n_answered": "391",
        "n_false_susceptible": "9", "n_false_resistant": "24",
        "sensitivity": "0.9360", "specificity": "0.9040",
        "sensitivity_ci": "0.8831-0.9661", "ppv": "0.8460", "npv": "0.9620",
        "resistance_prevalence": "0.3610", "sensitivity_basis_n": "141",
        "specificity_basis_n": "250", "notes": "",
    }]
    LINEAGE = [
        {"drug": "rifampicin", "lineage": "L2", "n_called": "120",
         "sensitivity": "0.95", "specificity": "0.91", "vme_rate": "0.05",
         "me_rate": "0.09", "powered": "yes"},
        {"drug": "rifampicin", "lineage": "L4", "n_called": "80",
         "sensitivity": "0.88", "specificity": "0.93", "vme_rate": "0.12",
         "me_rate": "0.07", "powered": "yes"},
    ]

    def test_a_benchmark_makes_the_drug_estimable(self):
        html = render_html(demo_report(), baseline_rows=self.BASELINE)
        drugs = {d["id"]: d for d in payload_of(html)["drugs"]}
        self.assertTrue(drugs["RIF"]["estimable"])
        self.assertEqual(drugs["RIF"]["error"], 9.6)
        self.assertEqual(drugs["RIF"]["vme"], 9)

    def test_lineage_rows_populate_the_stratified_view(self):
        html = render_html(demo_report(), baseline_rows=self.BASELINE,
                           lineage_rows=self.LINEAGE)
        payload = payload_of(html)
        self.assertEqual(payload["lineages"], ["L2", "L4"])
        self.assertIn("RIF", payload["lineage"])

    def test_an_unknown_only_lineage_is_not_stratification(self):
        rows = [dict(self.LINEAGE[0], lineage="unknown")]
        html = render_html(demo_report(), lineage_rows=rows)
        payload = payload_of(html)
        self.assertEqual(payload["lineages"], [])
        self.assertIn("Lineage was not assessed", html)

    def test_measurability_reaches_the_interpretation(self):
        html = render_html(demo_report(), measurability=[
            {"drug": "pretomanid", "measurable": "no", "missing_side": "both",
             "consequence": "no pairing of these sources can evaluate it"}])
        self.assertIn("pretomanid", html)
        self.assertIn("absent from both", html)

    def test_audit_entries_render_when_supplied(self):
        html = render_html(demo_report(), audit_entries=[
            {"timestamp": "2026-09-14T10:22:14Z", "action": "call_promoted",
             "actor": "myconductor", "target": "x", "hash": "a3f9c2e1dead"}])
        payload = payload_of(html)
        self.assertEqual(payload["audit"][0]["event"], "CALL_PROMOTED")
        self.assertEqual(payload["audit"][0]["hash"], "a3f9c2e1")


if __name__ == "__main__":
    unittest.main()
