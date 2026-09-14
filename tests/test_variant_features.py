"""Imported feature vectors: a missing feature must stay missing.

A model handed ``delta_delta_g = 0.0`` for a variant nobody folded cannot tell
"energetically neutral" apart from "no structure exists". Imputing at the
feature level reintroduces the absence-of-evidence confusion one layer below the
call states, where it is harder to see.

The second concern is drift. An ``available`` flag stored beside the values is a
second source of truth about the same fact, and the two eventually disagree
silently — this codebase has already been bitten by exactly that with a position
index maintained beside the data it indexed. So availability here is derived,
and a feature valued *and* declared unavailable is refused at construction.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from myconductor.core.models import VUSPriority
from myconductor.modules.variant_features import (
    KNOWN_FEATURES, UNAVAILABLE_REASONS, VariantFeatureTable, VariantFeatures,
)

ORGANISM = "Mycobacterium tuberculosis"
ASSEMBLY = "NC_000962.3"
KEY = "NC_000962.3:761155:C>T"


def features(**kwargs) -> VariantFeatures:
    base = dict(variant_key=KEY, organism=ORGANISM,
                reference_assembly=ASSEMBLY, source="feature-pipeline",
                source_version="1.4.0", timestamp="2026-09-14T00:00:00Z",
                values={"sift_score": 0.01, "residue_plddt": 92.4},
                unavailable={"delta_delta_g": "no_structure"})
    base.update(kwargs)
    return VariantFeatures(**base)


def priority(**kwargs) -> VUSPriority:
    base = dict(variant_label="rpoB_p.Ser450Leu", variant_key=KEY,
                gene="rpoB", drug="rifampicin", priority="high")
    base.update(kwargs)
    return VUSPriority(**base)


class MissingIsNotZeroTests(unittest.TestCase):
    def test_an_unavailable_feature_reads_as_none_not_zero(self):
        record = features()
        self.assertIsNone(record.get("delta_delta_g"))
        self.assertFalse(record.available("delta_delta_g"))
        self.assertEqual(record.reason("delta_delta_g"), "no_structure")

    def test_a_present_feature_reads_as_its_value(self):
        record = features()
        self.assertEqual(record.get("sift_score"), 0.01)
        self.assertTrue(record.available("sift_score"))
        self.assertIsNone(record.reason("sift_score"))

    def test_a_genuine_zero_is_distinguishable_from_absence(self):
        """The case the whole design exists for."""
        neutral = features(values={"delta_delta_g": 0.0}, unavailable={})
        absent = features(values={}, unavailable={"delta_delta_g": "no_structure"})
        self.assertEqual(neutral.get("delta_delta_g"), 0.0)
        self.assertTrue(neutral.available("delta_delta_g"))
        self.assertIsNone(absent.get("delta_delta_g"))
        self.assertFalse(absent.available("delta_delta_g"))


class NoDriftTests(unittest.TestCase):
    """Availability is derived, so it cannot disagree with the values."""

    def test_valued_and_unavailable_at_once_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            features(values={"sift_score": 0.01},
                     unavailable={"sift_score": "tool_failed"})
        self.assertIn("both valued and unavailable", str(caught.exception))

    def test_features_declared_neither_way_are_reported_not_assumed(self):
        record = features(values={"sift_score": 0.01}, unavailable={})
        self.assertIn("gerp_score", record.undeclared)
        self.assertIn("not declared either way", record.describe())

    def test_availability_tracks_the_values_with_no_second_flag(self):
        record = features(values={"gerp_score": 4.2}, unavailable={})
        self.assertTrue(record.available("gerp_score"))
        self.assertFalse(record.available("sift_score"))


class SchemaTests(unittest.TestCase):
    def test_an_unknown_feature_name_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            features(values={"vibes": 1.0})
        self.assertIn("unknown feature", str(caught.exception))

    def test_an_unknown_feature_name_is_refused_when_unavailable_too(self):
        with self.assertRaises(ValueError):
            features(unavailable={"vibes": "no_structure"})

    def test_a_free_text_reason_is_refused(self):
        """A reason a reader cannot act on is not a reason."""
        with self.assertRaises(ValueError) as caught:
            features(unavailable={"delta_delta_g": "dunno"})
        self.assertIn("needs a reason from", str(caught.exception))

    def test_non_finite_values_are_refused(self):
        for bad in (float("nan"), float("inf"), float("-inf")):
            with self.assertRaises(ValueError):
                features(values={"sift_score": bad})

    def test_booleans_are_not_numbers(self):
        with self.assertRaises(ValueError):
            features(values={"sift_score": True})

    def test_provenance_is_mandatory(self):
        for field_name in ("source", "source_version", "timestamp",
                           "organism", "reference_assembly", "variant_key"):
            with self.assertRaises(ValueError):
                features(**{field_name: "  "})

    def test_every_known_feature_has_a_usable_name(self):
        self.assertIn("sift_score", KNOWN_FEATURES)
        self.assertIn("binding_pocket_distance", KNOWN_FEATURES)
        self.assertIn("licence_restricted", UNAVAILABLE_REASONS)


class TableTests(unittest.TestCase):
    def test_duplicate_vectors_from_one_source_version_are_refused(self):
        with self.assertRaises(ValueError):
            VariantFeatureTable([features(), features()])

    def test_two_source_versions_coexist(self):
        table = VariantFeatureTable(
            [features(source_version="1.4.0"), features(source_version="2.0.0")])
        self.assertEqual(len(table.for_variant(KEY, ORGANISM, ASSEMBLY)), 2)

    def test_round_trips_through_json(self):
        table = VariantFeatureTable([features()])
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "features.json"
            path.write_text(json.dumps(table.to_dict()), encoding="utf-8")
            loaded = VariantFeatureTable.load(path)
        self.assertEqual(loaded.to_dict(), table.to_dict())

    def test_a_vector_for_another_assembly_does_not_match(self):
        table = VariantFeatureTable([features()])
        self.assertEqual(table.for_variant(KEY, ORGANISM, "GCF_other"), [])


class AttachmentTests(unittest.TestCase):
    def test_features_attach_without_touching_the_score(self):
        item = priority()
        item.score = 0.75
        VariantFeatureTable([features()]).attach([item], ORGANISM, ASSEMBLY)
        self.assertEqual(item.score, 0.75)
        dimension = item.dimensions[-1]
        self.assertEqual(dimension.name, "variant_features")
        self.assertIn("never imputed", dimension.note)
        self.assertIn("none of this enters the ranking", dimension.note)

    def test_the_attached_payload_keeps_absences_visible(self):
        item = priority()
        VariantFeatureTable([features()]).attach([item], ORGANISM, ASSEMBLY)
        payload = item.dimensions[-1].value[0]
        self.assertEqual(payload["unavailable"], {"delta_delta_g": "no_structure"})
        self.assertNotIn("delta_delta_g", payload["values"])
        self.assertIn("gerp_score", payload["undeclared"])

    def test_disagreeing_sources_are_preserved_not_averaged(self):
        """Two tools disagreeing is the signal, not noise to smooth away."""
        item = priority()
        table = VariantFeatureTable([
            features(source="tool-a", values={"sift_score": 0.01}, unavailable={}),
            features(source="tool-b", values={"sift_score": 0.98}, unavailable={}),
        ])
        table.attach([item], ORGANISM, ASSEMBLY)
        payload = item.dimensions[-1].value
        self.assertEqual(len(payload), 2)
        self.assertEqual({p["values"]["sift_score"] for p in payload},
                         {0.01, 0.98})

    def test_attaching_twice_replaces_rather_than_duplicates(self):
        item = priority()
        table = VariantFeatureTable([features()])
        table.attach([item], ORGANISM, ASSEMBLY)
        table.attach([item], ORGANISM, ASSEMBLY)
        names = [d.name for d in item.dimensions]
        self.assertEqual(names.count("variant_features"), 1)

    def test_a_variant_with_no_vector_is_left_alone(self):
        item = priority(variant_key="NC_000962.3:999999:A>G")
        VariantFeatureTable([features()]).attach([item], ORGANISM, ASSEMBLY)
        self.assertEqual(item.dimensions, [])


if __name__ == "__main__":
    unittest.main()
