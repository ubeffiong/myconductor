import unittest

from myconductor.core.models import Region, Route, Variant
from myconductor.core.router import TriageRouter


class RouterTests(unittest.TestCase):
    def setUp(self):
        self.router = TriageRouter()

    def _route(self, variant):
        return self.router.route([variant]).decisions[0].route

    def test_catalogued_variant_goes_to_catalogue(self):
        v = Variant(gene="katG", change="S315T", depth=80)
        self.assertEqual(self._route(v), Route.CATALOGUE)

    def test_efflux_gene_variant_goes_to_efflux(self):
        v = Variant(gene="Rv0678", change="L117R", depth=50)
        self.assertEqual(self._route(v), Route.EFFLUX_REGULATORY)

    def test_promoter_variant_goes_to_efflux(self):
        v = Variant(gene="eis", change="c.-14C>T", region=Region.PROMOTER, depth=60)
        self.assertEqual(self._route(v), Route.EFFLUX_REGULATORY)

    def test_unknown_coding_variant_goes_to_vus(self):
        v = Variant(gene="pncA", change="D12A", depth=70)
        self.assertEqual(self._route(v), Route.VUS_ML)

    def test_silent_variant_is_ignored(self):
        v = Variant(gene="rpoC", change="A542A", silent=True, depth=70)
        self.assertEqual(self._route(v), Route.IGNORED)

    def test_unmodelled_gene_is_ignored(self):
        v = Variant(gene="Rv9999", change="X1Y", depth=70)
        self.assertEqual(self._route(v), Route.IGNORED)

    def test_evidence_is_produced_for_catalogue_hit(self):
        v = Variant(gene="rpoB", change="S450L", depth=80)
        ev = self.router.route([v]).evidence
        self.assertEqual(len(ev), 1)
        self.assertEqual(ev[0].drug, "rifampicin")
        self.assertTrue(ev[0].call.is_resistant)


if __name__ == "__main__":
    unittest.main()
