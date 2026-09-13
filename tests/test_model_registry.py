import tempfile
import unittest
from pathlib import Path
from myconductor.core.models import InSilicoPrediction
from myconductor.federated.model_registry import ModelRegistry, RegisteredModel


class RegistryTests(unittest.TestCase):
    def model(self):
        return RegisteredModel("test-model", "1", "training/1", "test-organism", ["external-cohort"],
            [dict(drug="test-drug", lineage="L1", cohort="external-cohort", source="evaluation/1",
                  independent=True, sensitivity=0.8, specificity=0.8, call_rate=0.9)])

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
