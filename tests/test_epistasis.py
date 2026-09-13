import unittest
from dataclasses import replace
from myconductor.core.models import Variant, DrugResult, Call, Tier
from myconductor.modules.epistasis import EpistasisRule, EpistasisTable


class InteractionTests(unittest.TestCase):
    def rule(self, **kwargs):
        args = dict(rule_id="rule1", drug="test-drug", primary_variant_key_or_gene="geneA",
                    partner_variant_key_or_gene="geneB", interaction="compensatory", note="external annotation",
                    source="test reference", organism="test-organism", reference_assembly="NC_000962.3",
                    reviewed_by="curator", reviewed_at="2026-09-01")
        args.update(kwargs)
        return EpistasisRule(**args)

    def test_annotations_never_change_calls_tiers_or_confidence(self):
        variants = [Variant.of("geneA", "X1Y"), Variant.of("geneB", "X2Y")]
        for interaction in ("compensatory", "antagonistic"):
            result = DrugResult("test-drug", Call.RESISTANT, Tier.CATALOGUED, confidence=0.7)
            notes = EpistasisTable([self.rule(interaction=interaction)]).annotate(
                [result], variants, "test-organism", "NC_000962.3")
            self.assertEqual(len(notes), 1)
            self.assertEqual((result.call, result.tier, result.confidence), (Call.RESISTANT, Tier.CATALOGUED, 0.7))

    def test_wrong_organism_retracted_and_unresolved_have_no_annotation(self):
        variants = [Variant.of("geneA", "X1Y"), Variant.of("geneB", "X2Y")]
        result = DrugResult("test-drug", Call.INDETERMINATE, Tier.CATALOGUED)
        self.assertEqual(EpistasisTable([self.rule()]).annotate([result], variants, "test-organism", "NC_000962.3"), [])
        result.call = Call.RESISTANT
        for table, organism in ((EpistasisTable([self.rule(status="retracted")]), "test-organism"),
                                 (EpistasisTable([self.rule()]), "other")):
            self.assertEqual(table.annotate([result], variants, organism, "NC_000962.3"), [])

    def test_uncalibrated_boost_and_missing_source_rejected(self):
        for kwargs in ({"confidence_adjustment": 0.1}, {"source": ""}):
            with self.assertRaises(ValueError):
                self.rule(**kwargs)

    def test_same_variant_cannot_be_its_own_partner(self):
        table = EpistasisTable([self.rule(partner_variant_key_or_gene="geneA")])
        self.assertEqual(table.annotate([DrugResult("test-drug", Call.RESISTANT, Tier.CATALOGUED)],
                         [Variant.of("geneA", "X1Y")], "test-organism", "NC_000962.3"), [])
