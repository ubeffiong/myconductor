"""Cohort schema, lock, review and estimate — plus the shipped manifests."""
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.cohort import (
    CohortError,
    RowEvidence,
    derive_tier,
    is_mtbc,
    read_rows,
    schema_errors,
    summarise,
    verify_lock,
    write_lock,
    write_rows,
)
from mycobench.estimate import estimate, human
from mycobench.review import review, study_dependence

REPO = Path(__file__).resolve().parent.parent
COHORTS = REPO / "cohorts"

GOOD_ROW = {
    "sample_id": "ng_SRR1", "sra_run": "SRR1234567",
    "biosample": "SAMN12345678", "bioproject": "PRJNA123456",
    "organism": "Mycobacterium tuberculosis", "platform": "ILLUMINA",
    "layout": "PAIRED", "cohort_tier": "B", "bases": "400000000",
}
FIELDS = list(GOOD_ROW) + ["geo_loc_name", "phenotype_source",
                           "source_study", "expected_outcome"]


def sheet(rows, fields=None):
    path = Path(tempfile.mkdtemp()) / "cohort.tsv"
    write_rows(path, rows, fields or FIELDS)
    return path


class SchemaTests(unittest.TestCase):
    def test_a_good_row_passes(self):
        self.assertEqual(schema_errors([GOOD_ROW], FIELDS), [])

    def test_missing_required_column_is_reported_once(self):
        errors = schema_errors([GOOD_ROW], ["sample_id"])
        self.assertEqual(len(errors), 1)
        self.assertIn("must include", errors[0])

    def test_bad_accessions_are_caught(self):
        row = {**GOOD_ROW, "sra_run": "NOTARUN", "biosample": "X",
               "bioproject": "Y"}
        errors = schema_errors([row], FIELDS)
        self.assertTrue(any("SRA run" in e for e in errors))
        self.assertTrue(any("BioSample" in e for e in errors))
        self.assertTrue(any("BioProject" in e for e in errors))

    def test_duplicate_run_is_caught(self):
        rows = [GOOD_ROW, {**GOOD_ROW, "sample_id": "other"}]
        self.assertTrue(any("duplicate sra_run" in e
                            for e in schema_errors(rows, FIELDS)))

    def test_single_end_is_rejected(self):
        errors = schema_errors([{**GOOD_ROW, "layout": "SINGLE"}], FIELDS)
        self.assertTrue(any("PAIRED" in e for e in errors))

    def test_non_mtbc_row_needs_a_control_declaration(self):
        row = {**GOOD_ROW, "organism": "Mycobacterium avium"}
        errors = schema_errors([row], FIELDS)
        self.assertTrue(any("expected_outcome" in e for e in errors), errors)

    def test_non_mtbc_row_with_the_declaration_passes(self):
        row = {**GOOD_ROW, "organism": "Mycobacterium avium",
               "expected_outcome": "species_mismatch_refused"}
        self.assertEqual(schema_errors([row], FIELDS), [])

    def test_mtbc_row_cannot_declare_itself_a_species_control(self):
        row = {**GOOD_ROW, "expected_outcome": "species_mismatch_refused"}
        errors = schema_errors([row], FIELDS)
        self.assertTrue(any("is MTBC but declares" in e for e in errors))

    def test_is_mtbc_recognises_the_complex(self):
        self.assertTrue(is_mtbc("Mycobacterium tuberculosis"))
        self.assertTrue(is_mtbc("Mycobacterium tuberculosis variant bovis"))
        self.assertFalse(is_mtbc("Mycobacteroides abscessus"))


class TierTests(unittest.TestCase):
    def test_reviewed_study_derives_tier_a(self):
        self.assertEqual(derive_tier([], "Smith_2021_cohort"), "A")

    def test_absent_or_placeholder_study_derives_tier_b(self):
        self.assertEqual(derive_tier([], ""), "B")
        self.assertEqual(
            derive_tier([], "NCBI_discovery_pending_publication_review"), "B")

    def test_a_review_word_inside_a_real_citation_still_grades_a(self):
        self.assertEqual(derive_tier([], "Smith_2021_systematic_review"), "A")

    def test_evidence_errors_derive_no_tier(self):
        self.assertIsNone(derive_tier(["mismatch"], "Smith_2021"))


class LockTests(unittest.TestCase):
    def _lock(self, path, evidence=None):
        evidence = evidence or [RowEvidence(
            sample_id="ng_SRR1", run={"run": "SRR1234567", "bases": "400000000"},
            biosample={"accession": "SAMN12345678"}, derived_tier="B")]
        lock_path = path.with_suffix(".lock.json")
        return write_lock(lock_path, path, evidence, "Nigeria")

    def test_lock_verifies_the_exact_sheet(self):
        path = sheet([GOOD_ROW])
        lock = self._lock(path)
        self.assertEqual(verify_lock(lock, path), 1)

    def test_edited_sheet_fails_its_lock(self):
        path = sheet([GOOD_ROW])
        lock = self._lock(path)
        path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        with self.assertRaises(CohortError) as ctx:
            verify_lock(lock, path)
        self.assertIn("checksum", str(ctx.exception))

    def test_lock_with_recorded_errors_fails(self):
        path = sheet([GOOD_ROW])
        lock = self._lock(path, [RowEvidence(
            sample_id="ng_SRR1", errors=["SRA BioSample does not match"])])
        with self.assertRaises(CohortError) as ctx:
            verify_lock(lock, path)
        self.assertIn("evidence errors", str(ctx.exception))

    def test_stale_schema_version_is_rejected(self):
        path = sheet([GOOD_ROW])
        lock = self._lock(path)
        payload = json.loads(lock.read_text(encoding="utf-8"))
        payload["schema_version"] = "0.1"
        lock.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaises(CohortError) as ctx:
            verify_lock(lock, path)
        self.assertIn("schema", str(ctx.exception))

    def test_missing_lock_is_an_error_not_a_pass(self):
        path = sheet([GOOD_ROW])
        with self.assertRaises(CohortError):
            verify_lock(path.with_suffix(".absent.json"), path)


class ReviewTests(unittest.TestCase):
    def _rows(self, counts):
        rows = []
        for project, n in counts.items():
            for i in range(n):
                rows.append({**GOOD_ROW, "sample_id": f"{project}_{i}",
                             "sra_run": f"SRR{project[-3:]}{i:04d}",
                             "bioproject": project,
                             "bases": str(1_000_000 * (n - i))})
        return rows

    def test_effective_n_reflects_clustering_not_row_count(self):
        dependence = study_dependence(self._rows({"PRJNA111": 90,
                                                  "PRJNA222": 10}))
        self.assertEqual(dependence.n_rows, 100)
        self.assertLess(dependence.effective_n, 2.0,
                        "one dominant study must not look like 100 independent "
                        "observations")

    def test_even_spread_has_effective_n_near_group_count(self):
        dependence = study_dependence(self._rows({f"PRJNA{i:03d}": 10
                                                  for i in range(1, 6)}))
        self.assertAlmostEqual(dependence.effective_n, 5.0, places=6)

    def test_cap_holds_rows_back_rather_than_deleting_them(self):
        rows = self._rows({"PRJNA111": 10, "PRJNA222": 3})
        result = review(rows, max_per_bioproject=4)
        self.assertEqual(len(result.shortlist), 7)
        self.assertEqual(len(result.held_back), 6)
        self.assertEqual(len(result.shortlist) + len(result.held_back),
                         len(rows))

    def test_cap_improves_effective_n(self):
        rows = self._rows({"PRJNA111": 90, "PRJNA222": 10})
        result = review(rows, max_per_bioproject=5)
        before = result.dependence["bioproject"].effective_n
        after = study_dependence(result.shortlist).effective_n
        self.assertGreater(after, before)

    def test_cap_is_deterministic_and_prefers_depth(self):
        rows = self._rows({"PRJNA111": 5})
        first = [r["sample_id"] for r in review(rows, max_per_bioproject=2).shortlist]
        second = [r["sample_id"] for r in review(rows, max_per_bioproject=2).shortlist]
        self.assertEqual(first, second)
        self.assertEqual(first, ["PRJNA111_0", "PRJNA111_1"])

    def test_every_held_back_row_has_a_reason(self):
        result = review(self._rows({"PRJNA111": 6}), max_per_bioproject=2)
        self.assertEqual(len(result.reasons), len(result.held_back))


class EstimateTests(unittest.TestCase):
    def test_human_readable_sizes(self):
        self.assertEqual(human(0), "0 B")
        self.assertIn("GB", human(5 * 1024 ** 3))

    def test_estimate_uses_the_bases_column(self):
        path = sheet([GOOD_ROW])
        result = estimate(path, Path(tempfile.mkdtemp()))
        self.assertEqual(result.pending, 1)
        self.assertEqual(result.known_bases, 400_000_000)
        self.assertGreater(result.mid, 0)
        self.assertLess(result.low, result.mid)
        self.assertGreater(result.high, result.mid)

    def test_missing_bases_are_named_not_invented(self):
        path = sheet([{**GOOD_ROW, "bases": ""}])
        result = estimate(path, Path(tempfile.mkdtemp()))
        self.assertEqual(result.known_bases, 0)
        self.assertEqual(result.unknown, ["ng_SRR1"])
        text = result.render(path)
        self.assertIn("size unknown", text)
        self.assertIn("NOT included in the total", text)

    def test_present_reads_are_skipped(self):
        data_dir = Path(tempfile.mkdtemp())
        directory = data_dir / "ng_SRR1"
        directory.mkdir(parents=True)
        (directory / "SRR1234567_1.fastq.gz").write_bytes(b"x")
        (directory / "SRR1234567_2.fastq.gz").write_bytes(b"x")
        result = estimate(sheet([GOOD_ROW]), data_dir)
        self.assertEqual(result.pending, 0)
        self.assertEqual(result.already_present, 1)
        self.assertTrue(result.nothing_to_fetch)


class ShippedManifestTests(unittest.TestCase):
    """The committed panels must satisfy the schema they are validated by."""

    def _rows(self, name):
        return read_rows(COHORTS / name)

    def test_every_shipped_panel_passes_schema_validation(self):
        for name in ("nigeria-v1.tsv", "nigeria-full.tsv", "ntm-controls.tsv"):
            with self.subTest(panel=name):
                rows, fields = self._rows(name)
                self.assertTrue(rows, f"{name} has no rows")
                self.assertEqual(schema_errors(rows, fields), [], name)

    def test_nigeria_v1_is_capped_at_four_per_bioproject(self):
        rows, _ = self._rows("nigeria-v1.tsv")
        counts = study_dependence(rows).counts
        self.assertTrue(all(n <= 4 for n in counts.values()), counts)

    def test_every_nigerian_row_confirms_nigeria(self):
        for name in ("nigeria-v1.tsv", "nigeria-full.tsv"):
            rows, _ = self._rows(name)
            for row in rows:
                self.assertTrue(
                    row["geo_loc_name"].lower().startswith("nigeria"),
                    f"{name}: {row['sample_id']} has geo_loc_name "
                    f"{row['geo_loc_name']!r}")

    def test_no_nigerian_row_claims_a_phenotype(self):
        rows, _ = self._rows("nigeria-full.tsv")
        for row in rows:
            self.assertEqual(row.get("phenotype_source", ""), "",
                             "the Nigerian panels carry no phenotypes")

    def test_ntm_panel_is_entirely_non_mtbc_controls(self):
        rows, _ = self._rows("ntm-controls.tsv")
        for row in rows:
            self.assertFalse(is_mtbc(row["organism"]), row["sample_id"])
            self.assertEqual(row["expected_outcome"],
                             "species_mismatch_refused")

    def test_full_panel_clustering_is_declared_by_the_summary(self):
        rows, _ = self._rows("nigeria-full.tsv")
        stats = summarise(rows)
        self.assertGreater(stats["largest_bioproject_share"], 0.30,
                           "this panel IS clustered; the summary must show it")

    def test_panels_do_not_share_runs(self):
        seen = {}
        for name in ("nigeria-v1.tsv", "ntm-controls.tsv"):
            rows, _ = self._rows(name)
            for row in rows:
                self.assertNotIn(row["sra_run"], seen)
                seen[row["sra_run"]] = name

    def test_v1_is_a_subset_of_full(self):
        v1, _ = self._rows("nigeria-v1.tsv")
        full, _ = self._rows("nigeria-full.tsv")
        full_runs = {r["sra_run"] for r in full}
        for row in v1:
            self.assertIn(row["sra_run"], full_runs)


if __name__ == "__main__":
    unittest.main()
