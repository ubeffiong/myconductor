import json
import tempfile
import unittest
from pathlib import Path

from mycobench.metrics import score_accuracy, score_accuracy_stratified
from mycobench.monitoring import MonitoringState, generalizability, ingest_cohort, lineage_summary
from mycobench.thresholds import MIN_LINEAGES_FOR_GENERALIZABILITY


def rows_lineage2(n_resistant=8, n_susceptible=8):
    rows = [("resistant", "R", "lineage2")] * n_resistant
    rows += [("susceptible", "S", "lineage2")] * n_susceptible
    return [("moxifloxacin", lineage, predicted, phenotype)
            for predicted, phenotype, lineage in rows]


def rows_lineage4(n_resistant=6, n_susceptible=6):
    rows = [("resistant", "R", "lineage4")] * n_resistant
    rows += [("susceptible", "S", "lineage4")] * n_susceptible
    return [("moxifloxacin", lineage, predicted, phenotype)
            for predicted, phenotype, lineage in rows]


class StratifiedScoringTests(unittest.TestCase):
    def test_stratification_matches_pooled_counts(self):
        pairs = [("resistant", "R"), ("susceptible", "S")] * 5
        triples = [(p, ph, "lineage2") for p, ph in pairs]
        pooled = score_accuracy(pairs, "moxifloxacin")
        by_lineage = score_accuracy_stratified(triples, "moxifloxacin")
        self.assertEqual(set(by_lineage), {"lineage2"})
        self.assertEqual(by_lineage["lineage2"].n_called, pooled.n_called)

    def test_unknown_lineage_is_a_visible_bucket_not_a_drop(self):
        triples = [("resistant", "R", ""), ("susceptible", "S", None)]
        by_lineage = score_accuracy_stratified(triples, "moxifloxacin")
        self.assertIn("unknown", by_lineage)
        self.assertEqual(by_lineage["unknown"].n_called, 2)


class MonitoringStateTests(unittest.TestCase):
    def test_ingesting_a_cohort_updates_counts(self):
        state = MonitoringState()
        state, notes = ingest_cohort(state, "cohort-1", rows_lineage2())
        self.assertEqual(notes, [])
        summary = lineage_summary(state)
        self.assertIn("moxifloxacin", summary)
        entry = summary["moxifloxacin"][0]
        self.assertEqual(entry.lineage, "lineage2")
        self.assertEqual(entry.accuracy.n_called, 16)

    def test_reingesting_the_same_cohort_is_a_noop(self):
        state = MonitoringState()
        state, _ = ingest_cohort(state, "cohort-1", rows_lineage2())
        before = json.dumps(state.to_json(), sort_keys=True)
        state, notes = ingest_cohort(state, "cohort-1", rows_lineage2())
        after = json.dumps(state.to_json(), sort_keys=True)
        self.assertEqual(before, after)
        self.assertTrue(notes)
        self.assertIn("already ingested", notes[0])

    def test_a_new_cohort_accumulates_on_top(self):
        state = MonitoringState()
        state, _ = ingest_cohort(state, "cohort-1", rows_lineage2(8, 8))
        state, _ = ingest_cohort(state, "cohort-2", rows_lineage2(2, 2))
        summary = lineage_summary(state)
        entry = summary["moxifloxacin"][0]
        self.assertEqual(entry.accuracy.n_called, 20)

    def test_persistence_round_trips(self):
        state = MonitoringState()
        state, _ = ingest_cohort(state, "cohort-1", rows_lineage2())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "monitoring_state.json"
            state.save(path)
            restored = MonitoringState.load(path)
        self.assertEqual(restored.ingested_cohort_ids, {"cohort-1"})
        self.assertEqual(restored.counts, state.counts)

    def test_missing_state_file_loads_empty(self):
        state = MonitoringState.load("/nonexistent/path/monitoring_state.json")
        self.assertEqual(state.counts, {})


class GeneralizabilityTests(unittest.TestCase):
    def test_one_lineage_is_not_generalizable(self):
        state = MonitoringState()
        state, _ = ingest_cohort(state, "cohort-1", rows_lineage2())
        gen = generalizability(state)["moxifloxacin"]
        self.assertEqual(gen.n_lineages, 1)
        self.assertFalse(gen.generalizable)

    def test_two_lineages_clears_the_bar(self):
        state = MonitoringState()
        state, _ = ingest_cohort(state, "cohort-1", rows_lineage2())
        state, _ = ingest_cohort(state, "cohort-2", rows_lineage4())
        gen = generalizability(state)["moxifloxacin"]
        self.assertEqual(gen.n_lineages, MIN_LINEAGES_FOR_GENERALIZABILITY)
        self.assertTrue(gen.generalizable)

    def test_unknown_lineage_does_not_count_toward_generalizability(self):
        state = MonitoringState()
        rows = [("moxifloxacin", "unknown", "resistant", "R")] * 30
        state, _ = ingest_cohort(state, "cohort-1", rows)
        gen = generalizability(state)
        self.assertNotIn("moxifloxacin", gen)


if __name__ == "__main__":
    unittest.main()
