import unittest

from myconductor.core.models import Call, Tier, Variant
from myconductor.federated.vus_feedback import (
    LocalValidationStore,
    VUSFeedbackError,
    ValidationMethod,
    VUSValidationRecord,
)
from myconductor.modules.local_validation import LocalValidationModule


def variant(**kwargs):
    kwargs.setdefault("gene", "Rv0678")
    kwargs.setdefault("change", "L114R")
    return Variant.of(kwargs.pop("gene"), kwargs.pop("change"), **kwargs)


def record(v, result, isolate_id="iso1", site_id="site-A", **kwargs):
    return VUSValidationRecord(
        variant_key=v.key(), variant_label=v.label(), gene=v.gene,
        drug="bedaquiline", isolate_id=isolate_id, site_id=site_id,
        method=ValidationMethod.MIC, result=result,
        submitted_by="lab-tech", rationale="MIC shift assay", **kwargs)


class RecordValidationTests(unittest.TestCase):
    def test_an_inconclusive_result_is_refused(self):
        v = variant()
        with self.assertRaises(VUSFeedbackError):
            VUSValidationRecord(
                variant_key=v.key(), variant_label=v.label(), gene=v.gene,
                drug="bedaquiline", isolate_id="iso1", site_id="site-A",
                method=ValidationMethod.MIC, result=Call.INDETERMINATE)


class StoreVerdictTests(unittest.TestCase):
    def setUp(self):
        self.store = LocalValidationStore()
        self.v = variant()

    def test_no_records_yields_no_verdict(self):
        self.assertIsNone(self.store.verdict(self.v.key(), "bedaquiline"))

    def test_consistent_records_yield_an_established_call(self):
        self.store.ingest(record(self.v, Call.RESISTANT, isolate_id="iso1"))
        self.store.ingest(record(self.v, Call.RESISTANT, isolate_id="iso2"))
        verdict = self.store.verdict(self.v.key(), "bedaquiline")
        self.assertEqual(verdict.call, Call.RESISTANT)
        self.assertTrue(verdict.consistent)
        self.assertEqual(verdict.n_isolates, 2)

    def test_conflicting_records_yield_indeterminate(self):
        self.store.ingest(record(self.v, Call.RESISTANT, isolate_id="iso1"))
        self.store.ingest(record(self.v, Call.SUSCEPTIBLE, isolate_id="iso2"))
        verdict = self.store.verdict(self.v.key(), "bedaquiline")
        self.assertEqual(verdict.call, Call.INDETERMINATE)
        self.assertFalse(verdict.consistent)

    def test_retraction_is_reversible_and_ledger_audited(self):
        self.store.ingest(record(self.v, Call.SUSCEPTIBLE, isolate_id="iso1"))
        self.store.ingest(record(self.v, Call.RESISTANT, isolate_id="iso2"))
        verdict = self.store.verdict(self.v.key(), "bedaquiline")
        self.assertEqual(verdict.call, Call.INDETERMINATE)

        self.store.retract(self.v.key(), "bedaquiline", 0, "reviewer-1",
                           "iso1's culture was later found mixed")
        verdict = self.store.verdict(self.v.key(), "bedaquiline")
        self.assertEqual(verdict.call, Call.RESISTANT)
        actions = [e.action for e in self.store.ledger.entries]
        self.assertEqual(actions, ["validated", "validated", "retracted"])
        ok, _ = self.store.ledger_verified()
        self.assertTrue(ok)

    def test_double_retraction_is_refused(self):
        self.store.ingest(record(self.v, Call.RESISTANT))
        self.store.retract(self.v.key(), "bedaquiline", 0, "r", "why")
        with self.assertRaises(VUSFeedbackError):
            self.store.retract(self.v.key(), "bedaquiline", 0, "r", "again")

    def test_persistence_round_trips(self):
        self.store.ingest(record(self.v, Call.RESISTANT))
        data = self.store.to_json()
        restored = LocalValidationStore.from_json(data)
        verdict = restored.verdict(self.v.key(), "bedaquiline")
        self.assertEqual(verdict.call, Call.RESISTANT)
        ok, _ = restored.ledger_verified()
        self.assertTrue(ok)


class LocalValidationModuleTests(unittest.TestCase):
    def setUp(self):
        self.store = LocalValidationStore()
        self.v = variant()
        self.lane = LocalValidationModule(self.store)

    def test_variant_with_no_record_does_not_apply(self):
        self.assertFalse(self.lane.applies_to(self.v))
        self.assertEqual(self.lane.evaluate(self.v), [])

    def test_historical_association_does_not_become_current_phenotype(self):
        self.store.ingest(record(self.v, Call.RESISTANT))
        self.assertTrue(self.lane.applies_to(self.v))
        evidence = self.lane.evaluate(self.v)
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].call, Call.INDETERMINATE)
        self.assertEqual(evidence[0].tier, Tier.INFERRED)
        self.assertEqual(evidence[0].scope, "historical")

    def test_conflicting_verdict_never_asserts_resistance(self):
        self.store.ingest(record(self.v, Call.RESISTANT, isolate_id="iso1"))
        self.store.ingest(record(self.v, Call.SUSCEPTIBLE, isolate_id="iso2"))
        evidence = self.lane.evaluate(self.v)
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].call, Call.INDETERMINATE)

    def test_validated_labels_excludes_the_vus_workbench(self):
        self.store.ingest(record(self.v, Call.SUSCEPTIBLE))
        self.assertIn(self.v.label(), self.store.validated_labels())

    def test_retracted_variant_no_longer_applies(self):
        self.store.ingest(record(self.v, Call.RESISTANT))
        self.store.retract(self.v.key(), "bedaquiline", 0, "r", "mislabelled")
        self.assertFalse(self.lane.applies_to(self.v))


if __name__ == "__main__":
    unittest.main()
