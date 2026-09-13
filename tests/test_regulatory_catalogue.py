import unittest
from myconductor.catalogue.regulatory import RegulatoryCatalogue, RegulatoryRegion
from myconductor.core.models import Variant, Region


class RegulatoryTests(unittest.TestCase):
    def entry(self, **kwargs):
        data = dict(id="r1", name="test region", organism="test-organism", reference_assembly="NC_000962.3",
                    type="PROMOTER", target_genes=["testGene"], drug_associations=["test-drug"],
                    evidence_tier="CANDIDATE", source="test source")
        data.update(kwargs)
        return RegulatoryRegion(**data)

    def test_unknown_coordinates_are_not_invented(self):
        entry = self.entry()
        self.assertIsNone(entry.coordinates)
        variant = Variant.of("testGene", "c.-1A>T", region=Region.PROMOTER)
        self.assertEqual(RegulatoryCatalogue([entry]).matches(variant, "test-organism"), [entry])
        self.assertEqual(RegulatoryCatalogue([entry]).matches(variant, "other"), [])

    def test_coordinate_boundaries_and_assembly(self):
        entry = self.entry(coordinates={"chrom":"c", "start":10, "end":20, "strand":"+"})
        table = RegulatoryCatalogue([entry])
        for pos, match in ((9, False), (10, True), (20, True), (21, False)):
            variant = Variant.of("testGene", "X", chrom="c", pos=pos, ref="A", alt="G")
            self.assertEqual(bool(table.matches(variant, "test-organism")), match)

    def test_invalid_coordinates_and_tiers_rejected(self):
        for args in ({"coordinates":{"chrom":"c","start":0,"end":2,"strand":"+"}}, {"evidence_tier":"RESISTANT"}):
            with self.assertRaises(ValueError):
                self.entry(**args)
