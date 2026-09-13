import unittest

from myconductor.federated.catalogue_governance import (
    DiscordanceCategory,
    SiteCallRecord,
    reconcile_sites,
)


def rec(site_id, call, version, variant_key="rpoB_S450L", drug="rifampicin",
       rule_id=None):
    return SiteCallRecord(
        site_id=site_id, isolate_id=f"{site_id}-iso1", variant_key=variant_key,
        gene="rpoB", drug=drug, catalogue_version=version,
        rule_id=rule_id or f"{version}:{variant_key}:{drug}:{call}",
        call=call, tier="catalogued", timestamp_utc="2026-01-01T00:00:00+00:00",
    )


class ReconciliationTests(unittest.TestCase):
    def test_agreement_produces_no_discordance(self):
        records = [rec("site-A", "resistant", "v2"),
                  rec("site-B", "resistant", "v2")]
        self.assertEqual(reconcile_sites(records), [])

    def test_different_catalogue_versions_are_classified_as_such(self):
        records = [rec("site-A", "indeterminate", "v2"),
                  rec("site-B", "resistant", "v3")]
        discordances = reconcile_sites(records)
        self.assertEqual(len(discordances), 1)
        self.assertEqual(discordances[0].category,
                         DiscordanceCategory.CATALOGUE_VERSION)
        self.assertIn("catalogue version", discordances[0].note)

    def test_same_catalogue_version_disagreement_is_a_rule_problem(self):
        records = [rec("site-A", "resistant", "v3"),
                  rec("site-B", "susceptible", "v3")]
        discordances = reconcile_sites(records)
        self.assertEqual(len(discordances), 1)
        self.assertEqual(discordances[0].category, DiscordanceCategory.RULE)
        self.assertIn("pipeline/rule", discordances[0].note)

    def test_missing_catalogue_version_does_not_establish_a_cause(self):
        records = [rec("site-A", "resistant", None),
                  rec("site-B", "susceptible", None)]
        discordances = reconcile_sites(records)
        self.assertEqual(discordances[0].category, DiscordanceCategory.UNKNOWN)

    def test_never_resolves_by_vote(self):
        records = [rec("site-A", "resistant", "v3"),
                  rec("site-B", "resistant", "v3"),
                  rec("site-C", "susceptible", "v3")]
        discordances = reconcile_sites(records)
        self.assertEqual(len(discordances), 1)
        self.assertEqual(set(discordances[0].calls_by_site.values()),
                         {"resistant", "susceptible"})

    def test_independent_variant_drug_pairs_are_reported_separately(self):
        records = [rec("site-A", "resistant", "v2", variant_key="rpoB_S450L"),
                  rec("site-B", "indeterminate", "v3", variant_key="rpoB_S450L"),
                  rec("site-A", "resistant", "v3", drug="isoniazid",
                      variant_key="katG_S315T"),
                  rec("site-B", "resistant", "v3", drug="isoniazid",
                      variant_key="katG_S315T")]
        discordances = reconcile_sites(records)
        self.assertEqual(len(discordances), 1)
        self.assertEqual(discordances[0].drug, "rifampicin")


if __name__ == "__main__":
    unittest.main()
