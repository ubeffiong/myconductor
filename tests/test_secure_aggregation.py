import unittest

from myconductor.federated.secure_aggregation import (
    PairwiseMaskingAggregator,
    SecureAggregationError,
    derive_test_secrets,
    secrets_for_site,
    unmask_sum,
)

SEED = b"test-coordinator-seed-not-a-real-secret"


def build_aggregators(site_ids, vectors_by_site, seed=SEED):
    pairwise = derive_test_secrets(site_ids, seed)
    aggregators = {}
    for site_id in site_ids:
        others = [s for s in site_ids if s != site_id]
        aggregators[site_id] = PairwiseMaskingAggregator(
            site_id=site_id, other_site_ids=others,
            pairwise_secrets=secrets_for_site(site_id, pairwise),
            acknowledged=True,
        )
    return aggregators


class RoundTripTests(unittest.TestCase):
    def test_masked_sum_recovers_the_true_aggregate(self):
        site_ids = ["site-A", "site-B", "site-C"]
        vectors = {
            "site-A": [5, 0, 12],
            "site-B": [1, 3, 0],
            "site-C": [2, 2, 4],
        }
        aggregators = build_aggregators(site_ids, vectors)
        masked = [aggregators[s].mask(vectors[s]) for s in site_ids]
        result = unmask_sum(masked)
        self.assertEqual(result, [8, 5, 16])

    def test_masked_vectors_do_not_equal_the_raw_values(self):
        site_ids = ["a", "b"]
        vectors = {"a": [10, 20], "b": [1, 2]}
        aggregators = build_aggregators(site_ids, vectors)
        masked_a = aggregators["a"].mask(vectors["a"])
        self.assertNotEqual(masked_a, vectors["a"])

    def test_two_site_round_trips(self):
        site_ids = ["north", "south"]
        vectors = {"north": [100], "south": [7]}
        aggregators = build_aggregators(site_ids, vectors)
        masked = [aggregators[s].mask(vectors[s]) for s in site_ids]
        self.assertEqual(unmask_sum(masked), [107])

    def test_a_dropped_site_produces_a_wrong_silent_result(self):
        """Documents the known limitation (module docstring point 3): this
        is not a feature to celebrate, it is the reason a real deployment
        needs the full protocol's dropout recovery path."""
        site_ids = ["a", "b", "c"]
        vectors = {"a": [10], "b": [10], "c": [10]}
        aggregators = build_aggregators(site_ids, vectors)
        masked = [aggregators["a"].mask(vectors["a"]),
                 aggregators["b"].mask(vectors["b"])]  # "c" drops out
        result = unmask_sum(masked)
        self.assertNotEqual(result, [20])  # true sum of a+b alone


class GuardrailTests(unittest.TestCase):
    def test_construction_without_acknowledgement_is_refused(self):
        with self.assertRaises(SecureAggregationError):
            PairwiseMaskingAggregator(
                site_id="a", other_site_ids=["b"],
                pairwise_secrets={"b": b"x" * 32})

    def test_missing_pairwise_secret_is_refused(self):
        with self.assertRaises(SecureAggregationError):
            PairwiseMaskingAggregator(
                site_id="a", other_site_ids=["b", "c"],
                pairwise_secrets={"b": b"x" * 32}, acknowledged=True)

    def test_mismatched_vector_lengths_are_refused(self):
        with self.assertRaises(SecureAggregationError):
            unmask_sum([[1, 2], [1, 2, 3]])

    def test_empty_round_is_refused(self):
        with self.assertRaises(SecureAggregationError):
            unmask_sum([])


class SecretDerivationTests(unittest.TestCase):
    def test_derived_secrets_are_symmetric_between_the_two_sites(self):
        pairwise = derive_test_secrets(["a", "b"], SEED)
        self.assertEqual(len(pairwise), 1)
        a_secrets = secrets_for_site("a", pairwise)
        b_secrets = secrets_for_site("b", pairwise)
        self.assertEqual(a_secrets["b"], b_secrets["a"])

    def test_a_site_never_receives_a_secret_it_is_not_party_to(self):
        pairwise = derive_test_secrets(["a", "b", "c"], SEED)
        a_secrets = secrets_for_site("a", pairwise)
        self.assertEqual(set(a_secrets), {"b", "c"})


if __name__ == "__main__":
    unittest.main()
