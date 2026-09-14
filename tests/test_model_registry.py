import tempfile
import unittest
from pathlib import Path
from myconductor.core.models import InSilicoPrediction
from myconductor.federated.model_registry import (
    ModelRegistry, RegisteredModel, baseline_verdict,
)


def performance(call_rate=0.95, error_rate=0.04,
                base_coverage=0.90, base_error=0.08, **kwargs):
    """One evaluated scope. Defaults describe a model that beats its baseline."""
    row = dict(drug="test-drug", lineage="L1", cohort="external-cohort",
               source="evaluation/1", independent=True,
               sensitivity=0.8, specificity=0.8,
               call_rate=call_rate, error_rate=error_rate,
               baseline=dict(source="WHO catalogue", coverage=base_coverage,
                             error_rate=base_error))
    row.update(kwargs)
    return row


class RegistryTests(unittest.TestCase):
    def model(self, rows=None):
        return RegisteredModel("test-model", "1", "training/1", "test-organism",
                               ["external-cohort"], rows or [performance()])

    def prediction(self):
        return InSilicoPrediction("v", "test-drug", "test-model", "1", "training/1",
                "uncertain", "predictions/1", "s", "i", "site", "test-organism", "2026-09-01")

    def test_approval_and_rollback_are_replayed(self):
        registry = ModelRegistry()
        registry.submit(self.model(), "submitter", "new version")
        with self.assertRaises(ValueError):
            registry.require_approved(self.prediction(), "L1")
        registry.transition("test-model", "1", "reviewed", "reviewer", "reviewed evidence")
        registry.transition("test-model", "1", "approved", "reviewer", "external evaluation accepted")
        registry.require_approved(self.prediction(), "L1")
        with self.assertRaises(ValueError):
            registry.require_approved(self.prediction(), "L2")
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "registry.json"
            registry.save(path)
            loaded = ModelRegistry.load(path)
            self.assertEqual(registry.state(), loaded.state())
        registry.transition("test-model", "1", "rolled_back", "reviewer", "validation withdrawn")
        with self.assertRaises(ValueError):
            registry.require_approved(self.prediction(), "L1")

    def test_approval_cannot_skip_review_or_evaluation(self):
        registry = ModelRegistry()
        registry.submit(self.model(), "submitter", "new")
        with self.assertRaises(ValueError):
            registry.transition("test-model", "1", "approved", "reviewer", "skip")
        self.assertEqual(len(registry.ledger.entries), 1)


class BaselineFalsificationTests(unittest.TestCase):
    """A reported number is not evidence of improvement.

    Without an incumbent to compare against, sensitivity 0.40 and 0.95 are both
    "independently evaluated performance" and both pass a paperwork check. For
    TB-AMR the incumbent is never nothing: it is the catalogue, which answers
    where it holds a graded entry and abstains elsewhere. These tests hold the
    registry to that comparison rather than to the presence of a number.
    """

    def _model(self, rows):
        return RegisteredModel("m", "1", "training/1", "test-organism",
                               ["external-cohort"], rows)

    def _approve(self, rows):
        registry = ModelRegistry()
        registry.submit(self._model(rows), "submitter", "new")
        registry.transition("m", "1", "reviewed", "reviewer", "reviewed")
        registry.transition("m", "1", "approved", "reviewer", "accepted")
        return registry

    def test_a_model_that_beats_the_catalogue_is_approved(self):
        registry = self._approve([performance()])
        self.assertEqual(registry.state()[("m", "1")]["status"], "approved")

    def test_a_worse_error_rate_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self._approve([performance(error_rate=0.12, base_error=0.08)])
        self.assertIn("no improvement to deploy", str(caught.exception))

    def test_abstaining_more_than_the_incumbent_is_refused(self):
        """The abstention loophole: a model can always look accurate by
        answering less. That shifts work back to the clinician."""
        with self.assertRaises(ValueError) as caught:
            self._approve([performance(call_rate=0.55, error_rate=0.001,
                                       base_coverage=0.90)])
        self.assertIn("abstains more than the incumbent", str(caught.exception))

    def test_one_failing_scope_refuses_the_whole_approval(self):
        good = performance(drug="drug-a")
        bad = performance(drug="drug-b", error_rate=0.30)
        with self.assertRaises(ValueError) as caught:
            self._approve([good, bad])
        self.assertIn("drug-b", str(caught.exception))

    def test_performance_without_a_baseline_cannot_be_registered(self):
        row = performance()
        del row["baseline"]
        with self.assertRaises(ValueError) as caught:
            self._model([row])
        self.assertIn("requires a baseline", str(caught.exception))

    def test_baseline_needs_a_named_source(self):
        row = performance()
        row["baseline"] = dict(coverage=0.9, error_rate=0.08)
        with self.assertRaises(ValueError):
            self._model([row])

    def test_verdict_explains_itself_in_both_directions(self):
        beat, why = baseline_verdict(performance())
        self.assertTrue(beat)
        self.assertIn("WHO catalogue", why)
        beat, why = baseline_verdict(performance(error_rate=0.5))
        self.assertFalse(beat)
        self.assertIn("WHO catalogue", why)

    def test_evidence_for_states_the_basis_not_just_the_verdict(self):
        registry = self._approve([performance()])
        prediction = InSilicoPrediction(
            "v", "test-drug", "m", "1", "training/1", "uncertain",
            "predictions/1", "s", "i", "site", "test-organism", "2026-09-01")
        evidence = registry.evidence_for(prediction, "L1")
        self.assertEqual(evidence["baseline"], "WHO catalogue")
        self.assertEqual(evidence["cohort"], "external-cohort")
        self.assertIn("95%", evidence["basis"])
        # Beating a baseline is still not calibration or clinical validation.
        self.assertIn("cannot establish resistance", evidence["caveat"])
