import unittest
from types import SimpleNamespace
from myconductor.core.models import Subpopulation, Variant
from myconductor.modules.population_structure import cluster_subpopulations, population_structure


class PopulationTests(unittest.TestCase):
    def variants(self, values):
        return [Variant.of("testGene", f"X{i}Y", vaf=v) for i, v in enumerate(values)]

    def findings(self, variants):
        return [SimpleNamespace(variant_key=v.key(), assessable=True) for v in variants]

    def test_proximity_is_not_a_clone_or_microevolution(self):
        variants = self.variants([0.2, 0.21, 0.6, 0.61, None])
        result = population_structure("s", variants, self.findings(variants))
        self.assertEqual(len(result.subpopulations), 2)
        self.assertEqual(result.classification, "indeterminate")
        self.assertEqual(result.unclustered_variant_keys, [variants[-1].key()])
        self.assertTrue(result.data_gaps)

    def test_total_span_prevents_chain_merging(self):
        variants = self.variants([0.1, 0.14, 0.18])
        groups = cluster_subpopulations(variants, self.findings(variants))
        self.assertEqual(groups[0].variant_keys, {variants[0].key(), variants[1].key()})

    def test_linked_different_lineages_only_support_possible_mixture(self):
        variants = self.variants([0.2, 0.21, 0.6, 0.61])
        groups = cluster_subpopulations(variants, self.findings(variants))
        hints = [dict(variant_keys=sorted(g.variant_keys), lineage=f"L{i}",
                      source="external report", linkage_source="external phased report") for i,g in enumerate(groups)]
        result = population_structure("s", variants, self.findings(variants), lineage_hints=hints)
        self.assertEqual(result.classification, "possible_coinfection")
        hints[1]["lineage"] = "L0"
        self.assertEqual(population_structure("s", variants, self.findings(variants), lineage_hints=hints).classification, "indeterminate")

    def test_invalid_fraction_and_unmatched_lineage_refused(self):
        with self.assertRaises(ValueError):
            Subpopulation(frozenset({"x"}), float("nan"))
        variants = self.variants([0.2, 0.21])
        with self.assertRaises(ValueError):
            population_structure("s", variants, self.findings(variants), lineage_hints=[dict(
                variant_keys=["missing"], lineage="L", source="x", linkage_source="x")])

    def test_single_finding_does_not_create_population(self):
        variants = self.variants([0.2])
        self.assertIsNone(population_structure("s", variants, self.findings(variants)))
