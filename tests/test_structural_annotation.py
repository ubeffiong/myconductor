import unittest
from dataclasses import replace
from myconductor.core.models import Variant
from myconductor.modules.structural_annotation import StructuralAnnotation, StructuralAnnotationTable
from myconductor.modules.vus_workbench import VUSWorkbench


class StructuralTests(unittest.TestCase):
    def annotation(self, variant, **kwargs):
        data = dict(variant_key=variant.key(), location="unknown", predicted_effect="External uncertain annotation",
                    source="test archive", organism="mtbc", reference_assembly="NC_000962.3",
                    source_version="1", timestamp="2026-09-01")
        data.update(kwargs)
        return StructuralAnnotation(**data)

    def test_annotation_does_not_change_rank_or_invent_ligand_distance(self):
        variant = Variant.of("pncA", "X99Y")
        priorities = VUSWorkbench().priorities([variant])
        before = [(p.score, p.priority) for p in priorities]
        records, hypotheses = StructuralAnnotationTable([self.annotation(variant)]).attach(priorities, "mtbc", "NC_000962.3")
        self.assertEqual(before, [(p.score, p.priority) for p in priorities])
        self.assertTrue(records)
        self.assertTrue(hypotheses)
        self.assertFalse(next(d for d in priorities[0].dimensions if d.name == "ligand_distance").available)

    def test_cross_organism_and_wrong_identity_do_not_attach(self):
        variant = Variant.of("pncA", "X99Y")
        table = StructuralAnnotationTable([self.annotation(variant, organism="mabscessus")])
        self.assertEqual(table.attach(VUSWorkbench().priorities([variant]), "mtbc", "NC_000962.3"), ([], []))

    def test_distance_requires_units_and_source(self):
        variant = Variant.of("testGene", "X99Y")
        for kwargs in ({"ligand_distance": 3}, {"confidence": float("nan")}, {"source": ""}):
            with self.assertRaises(ValueError):
                self.annotation(variant, **kwargs)

    def test_conflicting_sources_remain_separate(self):
        variant = Variant.of("pncA", "X99Y")
        annotations = [self.annotation(variant), self.annotation(variant, source="other", location="distal")]
        records, _ = StructuralAnnotationTable(annotations).attach(VUSWorkbench().priorities([variant]), "mtbc", "NC_000962.3")
        self.assertEqual(len(records), 2)
