"""The bring-your-own-model benchmarking bridge.

Covers the join contract (a row with no matching accepted phenotype is
dropped with a reason, never silently scored), the matched-coverage
comparison against a real catalogue baseline, the fix for the drug that
kills the whole run when its baseline is missing (skipped-and-reported
instead of fatal), and that the resulting registry payload is actually
accepted by ``ModelRegistry``/``RegisteredModel`` unmodified.
"""
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.analysis import external_model as em
from myconductor.federated.model_registry import ModelRegistry, RegisteredModel


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


class LoadPredictionsTests(unittest.TestCase):
    def test_parses_a_well_formed_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "preds.tsv",
                        "isolate_id\tdrug\tpredicted\tconfidence\n"
                        "iso1\trifampicin\tresistant\t0.9\n"
                        "iso2\trifampicin\tsusceptible\t0.6\n")
            rows = em.load_predictions(path)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0].predicted, "R")
        self.assertEqual(rows[1].confidence, 0.6)

    def test_abstain_values_become_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "preds.tsv",
                        "isolate_id\tdrug\tpredicted\n"
                        "iso1\trifampicin\tabstain\n")
            rows = em.load_predictions(path)
        self.assertIsNone(rows[0].predicted)

    def test_unknown_predicted_value_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "preds.tsv",
                        "isolate_id\tdrug\tpredicted\n"
                        "iso1\trifampicin\tmaybe\n")
            with self.assertRaises(ValueError):
                em.load_predictions(path)

    def test_missing_column_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "preds.tsv", "isolate_id\tdrug\n iso1\trifampicin\n")
            with self.assertRaises(ValueError):
                em.load_predictions(path)

    def test_duplicate_isolate_drug_row_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "preds.tsv",
                        "isolate_id\tdrug\tpredicted\n"
                        "iso1\trifampicin\tresistant\n"
                        "iso1\trifampicin\tsusceptible\n")
            with self.assertRaises(ValueError):
                em.load_predictions(path)

    def test_out_of_range_confidence_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "preds.tsv",
                        "isolate_id\tdrug\tpredicted\tconfidence\n"
                        "iso1\trifampicin\tresistant\t1.5\n")
            with self.assertRaises(ValueError):
                em.load_predictions(path)


class LoadTruthsTests(unittest.TestCase):
    def test_generic_sample_drug_phenotype_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "truths.tsv",
                        "sample_id\tdrug\tphenotype\tquality\n"
                        "iso1\trifampicin\tR\tHIGH\n"
                        "iso2\trifampicin\tS\tLOW\n")
            truths = em.load_truths(path)
        self.assertEqual(truths[("iso1", "rifampicin")], "R")
        # LOW quality is not in ACCEPTED_QUALITY -> dropped, not defaulted.
        self.assertNotIn(("iso2", "rifampicin"), truths)

    def test_cryptic_style_reuse_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "reuse.csv",
                        "ENA_RUN,RIF_BINARY_PHENOTYPE,RIF_PHENOTYPE_QUALITY\n"
                        "ERR1,R,HIGH\n")
            truths = em.load_truths(path)
        self.assertEqual(truths[("cr_ERR1", "rifampicin")], "R")

    def test_unrecognised_schema_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "truths.tsv", "foo\tbar\n1\t2\n")
            with self.assertRaises(ValueError):
                em.load_truths(path)


class LoadCatalogueBaselinesTests(unittest.TestCase):
    def test_json_baseline_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "baseline.json"
            path.write_text(json.dumps({
                "provenance": {"catalogue": "who-v2"},
                "baselines": {"rifampicin": {"source": "who-v2", "coverage": 0.9,
                                             "error_rate": 0.05}},
            }), encoding="utf-8")
            baselines, source = em.load_catalogue_baselines(path)
        self.assertEqual(baselines["rifampicin"]["coverage"], 0.9)
        self.assertEqual(source, "who-v2")

    def test_csv_baseline_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write(Path(tmp) / "baseline.csv",
                        "drug,coverage,error_rate\nrifampicin,0.9,0.05\n")
            baselines, _source = em.load_catalogue_baselines(path)
        self.assertEqual(baselines["rifampicin"]["error_rate"], 0.05)


class MeasureExternalModelTests(unittest.TestCase):
    def _predictions(self, calls):
        return [em.ExternalPrediction(isolate, "rifampicin", predicted, 0.9)
                for isolate, predicted in calls]

    def test_join_drops_rows_with_no_accepted_truth_rather_than_scoring_them(self):
        predictions = self._predictions([("iso1", "R"), ("ghost", "S")])
        truths = {("iso1", "rifampicin"): "R"}
        grouped, dropped = em.to_selective_predictions(predictions, truths)
        self.assertEqual(len(grouped["rifampicin"]), 1)
        self.assertEqual(len(dropped), 1)
        self.assertEqual(dropped[0].isolate_id, "ghost")

    def test_a_drug_with_no_baseline_is_skipped_not_fatal(self):
        # Regression test: measure_external_model used to raise ValueError
        # for the whole run when ANY drug lacked a baseline. One
        # unmeasurable drug (e.g. pretomanid) must not discard every other
        # drug's result.
        predictions = (self._predictions([("iso1", "R"), ("iso2", "S")])
                      + [em.ExternalPrediction("iso3", "pretomanid", "R", 0.9)])
        truths = {("iso1", "rifampicin"): "R", ("iso2", "rifampicin"): "S",
                  ("iso3", "pretomanid"): "R"}
        baselines = {"rifampicin": {"source": "who-v2", "coverage": 0.9,
                                    "error_rate": 0.05}}
        results, skipped = em.measure_external_model(predictions, truths, baselines)
        self.assertIn("rifampicin", results)
        self.assertNotIn("pretomanid", results)
        self.assertTrue(any("pretomanid" in note for note in skipped))

    def test_registry_payload_is_accepted_by_registeredmodel_unmodified(self):
        # 12 isolates: 6 resistant correctly called, 6 susceptible correctly
        # called -> perfect accuracy, comfortably beating a lenient baseline.
        calls = [(f"r{i}", "R") for i in range(6)] + [(f"s{i}", "S") for i in range(6)]
        predictions = self._predictions(calls)
        truths = {(isolate, "rifampicin"): call for isolate, call in calls}
        baselines = {"rifampicin": {"source": "who-v2", "coverage": 0.5,
                                    "error_rate": 0.5}}
        results, skipped = em.measure_external_model(
            predictions, truths, baselines,
            model_id="test-model", model_version="1.0")
        self.assertEqual(skipped, [])
        self.assertIn("rifampicin", results)
        block = em.registry_payload(
            results, model_id="test-model", model_version="1.0",
            training_data_provenance="unit-test-fixture")
        # Must round-trip into RegisteredModel with no field massaging.
        model = RegisteredModel(**block)
        self.assertEqual(model.model_id, "test-model")
        self.assertTrue(model.performance)

    def test_a_model_no_better_than_baseline_is_measured_but_refused_on_approval(self):
        # registry_payload() reports what was measured, whether or not it
        # beats the baseline -- ModelRegistry's "approved" transition is
        # where beating the baseline is actually enforced (baseline_verdict).
        calls = [(f"r{i}", "S") for i in range(6)] + [(f"s{i}", "R") for i in range(6)]
        predictions = self._predictions(calls)  # every call flipped -> all wrong
        truths = {(isolate, "rifampicin"): ("R" if call == "S" else "S")
                 for isolate, call in calls}
        baselines = {"rifampicin": {"source": "who-v2", "coverage": 0.5,
                                    "error_rate": 0.1}}
        results, _skipped = em.measure_external_model(predictions, truths, baselines)
        self.assertFalse(results["rifampicin"].comparison.passed)
        block = em.registry_payload(
            results, model_id="bad-model", model_version="1.0",
            training_data_provenance="unit-test-fixture")
        self.assertEqual(len(block["performance"]), 1)

        model = RegisteredModel(**block)
        registry = ModelRegistry()
        registry.submit(model, reviewer="dr-a", rationale="submission")
        registry.transition(model.model_id, model.version, "reviewed",
                           reviewer="dr-b", rationale="reviewed")
        with self.assertRaises(ValueError):
            registry.transition(model.model_id, model.version, "approved",
                               reviewer="dr-b", rationale="attempted approval")

    def test_end_to_end_through_model_registry_approval(self):
        calls = [(f"r{i}", "R") for i in range(10)] + [(f"s{i}", "S") for i in range(10)]
        predictions = self._predictions(calls)
        truths = {(isolate, "rifampicin"): call for isolate, call in calls}
        baselines = {"rifampicin": {"source": "who-v2", "coverage": 0.5,
                                    "error_rate": 0.5}}
        results, _skipped = em.measure_external_model(
            predictions, truths, baselines,
            model_id="e2e-model", model_version="1.0")
        block = em.registry_payload(
            results, model_id="e2e-model", model_version="1.0",
            training_data_provenance="unit-test-fixture")
        model = RegisteredModel(**block)

        registry = ModelRegistry()
        registry.submit(model, reviewer="dr-a", rationale="initial submission")
        registry.transition(model.model_id, model.version, "reviewed",
                           reviewer="dr-b", rationale="looks complete")
        registry.transition(model.model_id, model.version, "approved",
                           reviewer="dr-b", rationale="beats catalogue baseline")
        state = registry.state()[(model.model_id, model.version)]
        self.assertEqual(state["status"], "approved")


if __name__ == "__main__":
    unittest.main()
