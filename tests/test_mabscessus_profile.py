"""A second bundled organism profile, exercised end to end.

Not a claim that M. abscessus resistance calling has been validated —
``OrganismProfile.ships_validated`` is False for this profile exactly as it
is for the bundled MTBC one. These tests only confirm the profile loads,
wires correctly, and that swapping organisms really does route through
profile/catalogue selection rather than editing engine code.
"""
import unittest

from myconductor import Myconductor
from myconductor.catalogue.profile import load_profile
from myconductor.core.models import Call, Tier, Variant
from myconductor.modules.catalogue import CatalogueModule


class ProfileLoadingTests(unittest.TestCase):
    def test_mabscessus_profile_loads_with_its_own_drugs(self):
        profile = load_profile(organism="mabscessus")
        self.assertEqual(profile.name, "mabscessus")
        self.assertIn("clarithromycin", profile.drugs)
        self.assertIn("amikacin", profile.drugs)
        self.assertNotIn("rifampicin", profile.drugs)
        self.assertFalse(profile.ships_validated)

    def test_mtbc_profile_is_unaffected_and_still_the_default(self):
        default_profile = load_profile()
        named_profile = load_profile(organism="mtbc")
        self.assertEqual(default_profile.drugs, named_profile.drugs)
        self.assertIn("rifampicin", default_profile.drugs)

    def test_unknown_organism_raises_a_clear_error(self):
        with self.assertRaises(ValueError):
            load_profile(organism="not-a-real-organism")

    def test_mabscessus_has_no_efflux_regulators_or_mtbc_genes(self):
        profile = load_profile(organism="mabscessus")
        self.assertEqual(profile.efflux_regulators, {})
        self.assertNotIn("Rv0678", profile.loci)


class CatalogueSelectionTests(unittest.TestCase):
    def test_mabscessus_catalogue_grades_its_own_variants(self):
        catalogue = CatalogueModule(organism="mabscessus")
        self.assertTrue(catalogue.is_illustrative)
        self.assertIn("erm41_T28", catalogue.labels)
        self.assertIn("rrs_1408A>G", catalogue.labels)
        self.assertNotIn("rpoB_S450L", catalogue.labels)

    def test_unknown_organism_catalogue_raises(self):
        with self.assertRaises(ValueError):
            CatalogueModule(organism="not-a-real-organism")


class EndToEndAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.conductor = Myconductor(organism="mabscessus", platform="illumina")

    def test_erm41_functional_sequevar_calls_clarithromycin_resistant(self):
        v = Variant.of("erm41", "T28")
        outcome = self.conductor.router.route([v])
        evidence = [e for e in outcome.evidence if e.drug == "clarithromycin"]
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].call, Call.RESISTANT)
        self.assertEqual(evidence[0].tier, Tier.CATALOGUED)

    def test_rrs_a1408g_calls_amikacin_resistant(self):
        v = Variant.of("rrs", "1408A>G", region="rrna")
        outcome = self.conductor.router.route([v])
        evidence = [e for e in outcome.evidence if e.drug == "amikacin"]
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].call, Call.RESISTANT)

    def test_erm41_nonfunctional_sequevar_is_susceptible_leaning(self):
        v = Variant.of("erm41", "C28")
        outcome = self.conductor.router.route([v])
        evidence = [e for e in outcome.evidence if e.drug == "clarithromycin"]
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0].call, Call.SUSCEPTIBLE)

    def test_mtbc_gene_does_not_apply_under_the_mabscessus_profile(self):
        v = Variant.of("rpoB", "S450L")
        outcome = self.conductor.router.route([v])
        self.assertEqual(outcome.evidence, [])

    def test_organism_does_not_change_the_default_mtbc_pipeline(self):
        default_conductor = Myconductor(platform="illumina")
        self.assertIn("rifampicin", default_conductor.profile.drugs)
        self.assertNotIn("clarithromycin", default_conductor.profile.drugs)


if __name__ == "__main__":
    unittest.main()
