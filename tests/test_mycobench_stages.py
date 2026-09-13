"""Stage orchestration, exercised without sra-tools or TB-Profiler installed.

``RecordingRunner`` is why this is testable: the external commands are injected
rather than hard-coded, so the contracts that matter — a failing isolate is
recorded and skipped, a stage aborts only when nothing succeeded — can be
checked on a machine with no bioinformatics toolchain at all.
"""
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.cohort import write_rows
from mycobench.stages import (
    RecordingRunner,
    RunContext,
    StageError,
    StageOutcome,
    reads_for,
    stage_download,
    write_manifest,
)

ROWS = [
    {"sample_id": "ng_A", "sra_run": "SRR1", "biosample": "SAMN1",
     "bioproject": "PRJNA111", "organism": "Mycobacterium tuberculosis",
     "platform": "ILLUMINA", "layout": "PAIRED", "cohort_tier": "B",
     "bases": "400000000", "expected_outcome": ""},
    {"sample_id": "ng_B", "sra_run": "SRR2", "biosample": "SAMN2",
     "bioproject": "PRJNA222", "organism": "Mycobacterium tuberculosis",
     "platform": "ILLUMINA", "layout": "PAIRED", "cohort_tier": "B",
     "bases": "300000000", "expected_outcome": ""},
    {"sample_id": "ntm_C", "sra_run": "SRR3", "biosample": "SAMN3",
     "bioproject": "PRJNA333", "organism": "Mycobacterium avium",
     "platform": "ILLUMINA", "layout": "PAIRED", "cohort_tier": "B",
     "bases": "500000000", "expected_outcome": "species_mismatch_refused"},
]
COLUMNS = list(ROWS[0])


class FakeRunner(RecordingRunner):
    """Records commands and can fail for chosen accessions."""

    def __init__(self, fail_for=(), missing_tools=()):
        super().__init__()
        self.fail_for = set(fail_for)
        self.present = tuple(t for t in ("prefetch", "fasterq-dump",
                                          "tb-profiler", "gzip")
                             if t not in missing_tools)
        self.touch = True

    def run(self, command, log_path, timeout=None):
        self.commands.append(list(command))
        if any(part in self.fail_for for part in command):
            return 3
        # Simulate the tool's side effect so the stage's own file checks run.
        if self.touch and command[0] == "fasterq-dump":
            out = Path(command[command.index("-O") + 1])
            run = command[-1]
            out.mkdir(parents=True, exist_ok=True)
            for suffix in ("_1", "_2"):
                (out / f"{run}{suffix}.fastq.gz").write_bytes(b"x")
        return 0


def context(rows=ROWS, runner=None) -> RunContext:
    root = Path(tempfile.mkdtemp())
    sheet = root / "panel.tsv"
    write_rows(sheet, rows, COLUMNS)
    return RunContext.create(sheet, root / "data", root / "results",
                             root / "logs", threads=2,
                             runner=runner or FakeRunner())


class ContextTests(unittest.TestCase):
    def test_controls_are_separated_from_benchmark_rows(self):
        ctx = context()
        self.assertEqual([r["sample_id"] for r in ctx.benchmark_rows],
                         ["ng_A", "ng_B"])
        self.assertEqual([r["sample_id"] for r in ctx.control_rows], ["ntm_C"])

    def test_directories_are_created(self):
        ctx = context()
        for directory in (ctx.data_dir, ctx.results_dir, ctx.log_dir):
            self.assertTrue(directory.is_dir())


class DownloadStageTests(unittest.TestCase):
    def test_each_isolate_gets_prefetch_then_fasterq_dump(self):
        runner = FakeRunner()
        ctx = context(runner=runner)
        outcome = stage_download(ctx)
        self.assertEqual(outcome.n_ok, 3)
        tools = [c[0] for c in runner.commands]
        self.assertEqual(tools.count("prefetch"), 3)
        self.assertEqual(tools.count("fasterq-dump"), 3)

    def test_threads_are_passed_through(self):
        runner = FakeRunner()
        stage_download(context(runner=runner))
        fasterq = next(c for c in runner.commands if c[0] == "fasterq-dump")
        self.assertIn("2", fasterq)

    def test_a_failing_isolate_is_recorded_and_the_others_continue(self):
        runner = FakeRunner(fail_for={"SRR2"})
        ctx = context(runner=runner)
        outcome = stage_download(ctx)
        self.assertEqual(outcome.n_ok, 2)
        self.assertEqual([s for s, _ in outcome.failed], ["ng_B"])
        self.assertIn("prefetch exited 3", outcome.failed[0][1])

    def test_status_table_is_written_with_reasons(self):
        ctx = context(runner=FakeRunner(fail_for={"SRR2"}))
        stage_download(ctx)
        status = ctx.results_dir / "download_status.tsv"
        self.assertTrue(status.is_file())
        text = status.read_text(encoding="utf-8")
        self.assertIn("ng_B", text)
        self.assertIn("failed", text)

    def test_present_reads_are_skipped_not_refetched(self):
        ctx = context()
        row = ctx.rows[0]
        r1, r2 = reads_for(ctx, row)
        r1.parent.mkdir(parents=True, exist_ok=True)
        r1.write_bytes(b"x")
        r2.write_bytes(b"x")
        outcome = stage_download(ctx)
        self.assertIn(("ng_A", "reads already present"), outcome.skipped)
        self.assertEqual(outcome.n_ok, 2)

    def test_total_failure_aborts_because_nothing_can_be_benchmarked(self):
        runner = FakeRunner(fail_for={"SRR1", "SRR2", "SRR3"})
        with self.assertRaises(StageError) as ctx_mgr:
            stage_download(context(runner=runner))
        self.assertIn("nothing to benchmark", str(ctx_mgr.exception))

    def test_missing_tool_is_a_clear_error_not_a_crash(self):
        runner = FakeRunner(missing_tools={"prefetch"})
        with self.assertRaises(StageError) as ctx_mgr:
            stage_download(context(runner=runner))
        message = str(ctx_mgr.exception)
        self.assertIn("prefetch", message)
        self.assertIn("container", message)

    def test_unpaired_output_is_a_recorded_failure(self):
        runner = FakeRunner()
        runner.touch = False          # tool "succeeds" but writes nothing
        outcome_ctx = context(runner=runner)
        with self.assertRaises(StageError):
            stage_download(outcome_ctx)
        text = (outcome_ctx.results_dir / "download_status.tsv").read_text(
            encoding="utf-8")
        self.assertIn("paired FASTQ not produced", text)


class OutcomeTests(unittest.TestCase):
    def test_status_rows_cover_every_state(self):
        outcome = StageOutcome("download", ok=["a"], failed=[("b", "boom")],
                               skipped=[("c", "present")])
        states = {r["status"] for r in outcome.status_rows()}
        self.assertEqual(states, {"ok", "failed", "skipped"})

    def test_summary_reports_all_three_counts(self):
        outcome = StageOutcome("profile", ok=["a", "b"],
                               failed=[("c", "x")], skipped=[])
        self.assertIn("2 ok", outcome.summary())
        self.assertIn("1 failed", outcome.summary())


class ManifestTests(unittest.TestCase):
    def test_manifest_records_versions_and_the_pre_registration_hash(self):
        ctx = context()
        outcome = StageOutcome("download", ok=["ng_A"])
        path = write_manifest(ctx, [outcome])
        payload = json.loads(path.read_text(encoding="utf-8"))
        self.assertIn("mycobench_version", payload)
        self.assertIn("myconductor_version", payload)
        self.assertEqual(len(payload["pre_registration"]["sha256"]), 64)
        self.assertEqual(payload["cohort"]["n_benchmark"], 2)
        self.assertEqual(payload["cohort"]["n_controls"], 1)

    def test_manifest_warns_when_the_bundled_catalogue_was_used(self):
        ctx = context()
        payload = json.loads(
            write_manifest(ctx, []).read_text(encoding="utf-8"))
        self.assertIn("not the WHO catalogue",
                      payload["catalogue"]["warning"])

    def test_manifest_states_accuracy_was_not_computed_without_phenotypes(self):
        ctx = context()
        payload = json.loads(
            write_manifest(ctx, []).read_text(encoding="utf-8"))
        self.assertIn("accuracy not computed", payload["phenotypes"]["note"])

    def test_stage_counts_are_recorded(self):
        ctx = context()
        outcomes = [StageOutcome("download", ok=["a", "b"],
                                 failed=[("c", "x")])]
        payload = json.loads(
            write_manifest(ctx, outcomes).read_text(encoding="utf-8"))
        self.assertEqual(payload["stages"][0],
                         {"stage": "download", "ok": 2, "failed": 1,
                          "skipped": 0})


if __name__ == "__main__":
    unittest.main()
