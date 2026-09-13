"""The HTML report, especially what it says when a measurement is absent."""
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.cohort import write_rows
from mycobench.report import GLOSSARY, build_report

COHORT = [
    {"sample_id": "ng_A", "sra_run": "SRR1", "biosample": "SAMN1",
     "bioproject": "PRJNA111", "organism": "Mycobacterium tuberculosis",
     "geo_loc_name": "Nigeria", "platform": "ILLUMINA", "layout": "PAIRED",
     "bases": "400000000", "cohort_tier": "B", "phenotype_source": "",
     "source_study": "pending"},
    {"sample_id": "ng_B", "sra_run": "SRR2", "biosample": "SAMN2",
     "bioproject": "PRJNA111", "organism": "Mycobacterium tuberculosis",
     "geo_loc_name": "Nigeria:BENUE", "platform": "ILLUMINA",
     "layout": "PAIRED", "bases": "300000000", "cohort_tier": "B",
     "phenotype_source": "", "source_study": "pending"},
    {"sample_id": "ng_C", "sra_run": "SRR3", "biosample": "SAMN3",
     "bioproject": "PRJNA222", "organism": "Mycobacterium tuberculosis",
     "geo_loc_name": "Nigeria", "platform": "ILLUMINA", "layout": "PAIRED",
     "bases": "300000000", "cohort_tier": "B", "phenotype_source": "",
     "source_study": "pending"},
]
COHORT_COLUMNS = list(COHORT[0])


def build(results: dict, cohort=COHORT) -> tuple[str, Path]:
    root = Path(tempfile.mkdtemp())
    results_dir = root / "results"
    results_dir.mkdir(parents=True)
    sheet = root / "panel.tsv"
    write_rows(sheet, cohort, COHORT_COLUMNS)
    for name, (rows, columns) in results.items():
        if name.endswith(".json"):
            (results_dir / name).write_text(json.dumps(rows), encoding="utf-8")
        else:
            write_rows(results_dir / name, rows, columns)
    out = build_report(results_dir, sheet, results_dir / "report.html")
    return out.read_text(encoding="utf-8"), out


MANIFEST = {
    "mycobench_version": "0.1.0", "myconductor_version": "0.2.0",
    "generated_utc": "2026-09-13T00:00:00Z",
    "cohort": {"path": "panel.tsv"},
    "catalogue": {"source": "bundled illustrative subset",
                  "warning": "not the WHO catalogue"},
    "phenotypes": {"source": None, "note": "no paired phenotypes"},
    "pre_registration": {"description": "pre-registration v1.1.0",
                         "sha256": "abc123"},
}

CONCORDANCE = ([
    {"drug": "rifampicin", "compared": "3", "agree": "2", "disagree": "1",
     "agreement_rate": "0.6667", "only_myconductor": "0", "only_engine": "1",
     "neither": "0"},
], ("drug", "compared", "agree", "disagree", "agreement_rate",
    "only_myconductor", "only_engine", "neither"))


class NoAccuracyTests(unittest.TestCase):
    """A concordance-only run must not read as an accuracy measurement."""

    def setUp(self):
        self.html, self.path = build({
            "run_manifest.json": (MANIFEST, None),
            "concordance.tsv": CONCORDANCE,
        })

    def test_headline_states_accuracy_was_not_computed(self):
        self.assertIn("No accuracy figure was computed", self.html)

    def test_headline_warns_agreement_is_not_correctness(self):
        self.assertIn("Agreement is not correctness", self.html)

    def test_accuracy_section_explains_how_to_get_one(self):
        self.assertIn("Not measured in this run", self.html)
        self.assertIn("fetch-phenotypes", self.html)

    def test_no_zero_is_presented_as_an_accuracy_value(self):
        self.assertNotIn("<td class='n'>0.000</td>", self.html)

    def test_concordance_is_rendered(self):
        self.assertIn("rifampicin", self.html)
        self.assertIn("66.7%", self.html)

    def test_catalogue_warning_surfaces(self):
        self.assertIn("not the WHO catalogue", self.html)

    def test_report_is_self_contained(self):
        for marker in ("http://", "https://", "<link"):
            self.assertNotIn(marker, self.html,
                             "the report must render offline")

    def test_boundaries_section_is_present(self):
        self.assertIn("What this report does not establish", self.html)
        self.assertIn("has not been clinically evaluated", self.html)


class ClusteringTests(unittest.TestCase):
    def test_dominant_bioproject_is_called_out_with_effective_n(self):
        cohort = [dict(COHORT[0], sample_id=f"ng_{i}", sra_run=f"SRR{i}",
                       biosample=f"SAMN{i}", bioproject="PRJNA111")
                  for i in range(9)]
        cohort.append(dict(COHORT[2], sample_id="ng_x", sra_run="SRR99",
                           biosample="SAMN99", bioproject="PRJNA222"))
        html, _ = build({"run_manifest.json": (MANIFEST, None)}, cohort)
        self.assertIn("effective n", html)
        self.assertIn("90%", html)
        self.assertIn("Treat intervals", html)


class AccuracyRenderingTests(unittest.TestCase):
    ACCURACY = ([
        {"drug": "rifampicin", "evaluable": "60", "called": "58",
         "call_rate": "0.9667", "abstained": "2", "tp": "28", "fp": "1",
         "tn": "28", "fn": "1", "sensitivity": "0.9655",
         "specificity": "0.9655", "ppv": "0.9655", "npv": "0.9655",
         "vme_rate": "0.0345", "me_rate": "0.0345", "verdict": "pass",
         "reasons": "every pre-registered bound met"},
        {"drug": "bedaquiline", "evaluable": "4", "called": "4",
         "call_rate": "1.0", "abstained": "0", "tp": "0", "fp": "0",
         "tn": "4", "fn": "0", "sensitivity": "", "specificity": "1.0",
         "ppv": "", "npv": "1.0", "vme_rate": "", "me_rate": "0.0",
         "verdict": "underpowered",
         "reasons": "4 evaluable isolate(s), below the pre-registered minimum"},
    ], ("drug", "evaluable", "called", "call_rate", "abstained", "tp", "fp",
        "tn", "fn", "sensitivity", "specificity", "ppv", "npv", "vme_rate",
        "me_rate", "verdict", "reasons"))

    def setUp(self):
        self.html, _ = build({
            "run_manifest.json": (MANIFEST, None),
            "concordance.tsv": CONCORDANCE,
            "accuracy.tsv": self.ACCURACY,
        })

    def test_headline_reports_the_verdict_mix(self):
        self.assertIn("Accuracy was measured", self.html)
        self.assertIn("1 met", self.html)

    def test_underpowered_is_shown_as_its_own_state(self):
        self.assertIn("underpowered", self.html)
        self.assertIn("v-under", self.html)

    def test_unmeasurable_sensitivity_is_not_rendered_as_zero(self):
        self.assertIn("not measured", self.html)

    def test_call_rate_is_shown_beside_the_error_rates(self):
        self.assertIn("call rate", self.html)
        self.assertIn("96.7%", self.html)

    def test_pre_registration_hash_is_printed(self):
        self.assertIn("pre-registration v1.1.0", self.html)


class SpeciesControlRenderingTests(unittest.TestCase):
    def test_a_failed_control_is_reported_prominently(self):
        controls = ([
            {"sample_id": "ntm_1", "organism": "Mycobacterium avium",
             "refused": "no", "passed": "no",
             "detail": "1 drug call(s) produced for a Mycobacterium avium genome"},
        ], ("sample_id", "organism", "refused", "passed", "detail"))
        html, _ = build({"run_manifest.json": (MANIFEST, None),
                         "species_control.tsv": controls})
        self.assertIn("species control(s) FAILED", html)
        self.assertIn("v-fail", html)

    def test_all_refused_reads_as_a_pass(self):
        controls = ([
            {"sample_id": "ntm_1", "organism": "Mycobacterium avium",
             "refused": "yes", "passed": "yes", "detail": "no drug call"},
        ], ("sample_id", "organism", "refused", "passed", "detail"))
        html, _ = build({"run_manifest.json": (MANIFEST, None),
                         "species_control.tsv": controls})
        self.assertIn("correctly refused", html)


class StageStatusTests(unittest.TestCase):
    def test_failures_remain_visible_with_their_reasons(self):
        status = ([
            {"sample_id": "ng_A", "status": "ok", "reason": ""},
            {"sample_id": "ng_B", "status": "failed",
             "reason": "prefetch exited 3"},
            {"sample_id": "ng_C", "status": "skipped", "reason": "reads absent"},
        ], ("sample_id", "status", "reason"))
        html, _ = build({"run_manifest.json": (MANIFEST, None),
                         "download_status.tsv": status})
        self.assertIn("1 failed", html)
        self.assertIn("prefetch exited 3", html)
        self.assertIn("never silently dropped", html)


class GlossaryTests(unittest.TestCase):
    def test_every_term_has_a_what_and_a_why(self):
        for term, body in GLOSSARY.items():
            with self.subTest(term=term):
                self.assertTrue(body["what"])
                self.assertTrue(body["why"])

    def test_glossary_is_embedded_for_tooltips(self):
        html, _ = build({"run_manifest.json": (MANIFEST, None)})
        self.assertIn('id="gl-data"', html)
        self.assertIn("Concordance", html)


if __name__ == "__main__":
    unittest.main()
