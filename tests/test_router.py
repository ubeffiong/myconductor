import unittest

from myconductor.core.models import Consequence, Lane, Region, Variant
from myconductor.core.router import TriageRouter
from myconductor.modules.catalogue import CatalogueModule
from myconductor.modules.vus_workbench import VUSWorkbench


class RouterTests(unittest.TestCase):
    def setUp(self):
        catalogue = CatalogueModule()
        self.router = TriageRouter(
            catalogue=catalogue,
            extra_lanes=[VUSWorkbench(known=catalogue.labels)],
        )

    def _lanes(self, variant):
        return set(self.router.route([variant]).decisions[0].lanes)

    def test_catalogued_variant_reaches_the_catalogue_lane(self):
        v = Variant.of("katG", "S315T", depth=80)
        self.assertIn(Lane.CATALOGUE, self._lanes(v))

    def test_efflux_gene_variant_reaches_the_efflux_lane(self):
        v = Variant.of("Rv0678", "L117R", depth=50,
                       consequence=Consequence.MISSENSE)
        self.assertIn(Lane.EFFLUX_REGULATORY, self._lanes(v))

    def test_promoter_variant_reaches_the_efflux_lane(self):
        v = Variant.of("eis", "c.-14C>T", region=Region.PROMOTER, depth=60,
                       consequence=Consequence.UPSTREAM)
        self.assertIn(Lane.EFFLUX_REGULATORY, self._lanes(v))

    def test_uncatalogued_coding_variant_reaches_the_vus_lane(self):
        v = Variant.of("pncA", "D12A", depth=70,
                       consequence=Consequence.MISSENSE)
        self.assertIn(Lane.VUS, self._lanes(v))

    def test_catalogued_variant_is_not_sent_to_the_vus_lane(self):
        v = Variant.of("katG", "S315T", depth=80)
        self.assertNotIn(Lane.VUS, self._lanes(v))

    def test_silent_variant_reaches_no_lane(self):
        v = Variant.of("rpoC", "A542A", depth=70,
                       consequence=Consequence.SYNONYMOUS)
        self.assertEqual(self._lanes(v), set())

    def test_unmodelled_gene_reaches_no_lane_but_is_reported(self):
        v = Variant.of("Rv9999", "X1Y", depth=70,
                       consequence=Consequence.MISSENSE)
        outcome = self.router.route([v])
        self.assertEqual(set(outcome.decisions[0].lanes), set())
        self.assertEqual(len(outcome.unexamined), 1,
                         "an uninterpretable variant must be reported, "
                         "not silently dropped")


class MultiLaneTests(unittest.TestCase):
    """The defect: exclusive routing let one lane silence the others."""

    def setUp(self):
        catalogue = CatalogueModule()
        self.router = TriageRouter(
            catalogue=catalogue,
            extra_lanes=[VUSWorkbench(known=catalogue.labels)],
        )

    def test_a_catalogued_efflux_variant_reaches_both_lanes(self):
        # Rv0678 R94Q is both catalogued and in an efflux regulator. The old
        # router picked the catalogue and discarded the efflux inference.
        v = Variant.of("Rv0678", "R94Q", depth=60,
                       consequence=Consequence.MISSENSE)
        lanes = set(self.router.route([v]).decisions[0].lanes)
        self.assertIn(Lane.CATALOGUE, lanes)
        self.assertIn(Lane.EFFLUX_REGULATORY, lanes)

    def test_evidence_from_several_lanes_is_all_retained(self):
        v = Variant.of("Rv0678", "R94Q", depth=60,
                       consequence=Consequence.MISSENSE)
        evidence = self.router.route([v]).evidence
        lanes = {ev.lane for ev in evidence}
        self.assertGreaterEqual(len(lanes), 2)

    def test_catalogue_evidence_covers_every_named_drug(self):
        v = Variant.of("Rv0678", "R94Q", depth=60)
        catalogue_evidence = [
            ev for ev in self.router.route([v]).evidence
            if ev.lane is Lane.CATALOGUE
        ]
        drugs = {ev.drug for ev in catalogue_evidence}
        self.assertEqual(drugs, {"bedaquiline", "clofazimine"})

    def test_lane_counts_may_exceed_the_variant_count(self):
        v = Variant.of("Rv0678", "R94Q", depth=60,
                       consequence=Consequence.MISSENSE)
        counts = self.router.route([v]).lane_counts()
        self.assertGreaterEqual(
            counts[Lane.CATALOGUE.value] + counts[Lane.EFFLUX_REGULATORY.value],
            2,
        )


if __name__ == "__main__":
    unittest.main()
