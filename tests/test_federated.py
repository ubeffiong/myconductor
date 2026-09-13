import unittest

from myconductor.federated.catalogue_update import (
    VariantObservation,
    aggregate,
    drift_alarm,
    promote_candidates,
)


class FederatedTests(unittest.TestCase):
    def test_aggregation_merges_across_sites(self):
        site_a = [VariantObservation("pncA_D12A", "pyrazinamide", 8, 1)]
        site_b = [VariantObservation("pncA_D12A", "pyrazinamide", 9, 2)]
        merged = aggregate([site_a, site_b])
        obs = merged["pncA_D12A|pyrazinamide"]
        self.assertEqual(obs.resistant_count, 17)
        self.assertEqual(obs.total, 20)

    def test_promotion_requires_evidence_bar(self):
        weak = aggregate([[VariantObservation("x_A1B", "rifampicin", 3, 0)]])
        self.assertEqual(promote_candidates(weak), [])
        strong = aggregate([[VariantObservation("x_A1B", "rifampicin", 19, 2)]])
        promoted = promote_candidates(strong)
        self.assertEqual(len(promoted), 1)
        self.assertEqual(promoted[0]["drug"], "rifampicin")

    def test_drift_alarm(self):
        self.assertFalse(drift_alarm([0.9, 0.88, 0.91]))
        self.assertTrue(drift_alarm([0.75, 0.72, 0.70]))


if __name__ == "__main__":
    unittest.main()
