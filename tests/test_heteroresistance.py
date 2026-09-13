import unittest

from myconductor.core.models import (
    Call,
    DrugEvidence,
    Lane,
    Tier,
    Variant,
    VariantIdentity,
)
from myconductor.modules.heteroresistance import (
    PLATFORM_LOD,
    HeteroresistanceAssessor,
)


def evidence_for(variant, drug="moxifloxacin"):
    return DrugEvidence(
        drug=drug, call=Call.RESISTANT, tier=Tier.CATALOGUED,
        lane=Lane.CATALOGUE, variant=variant.identity,
        rationale="catalogued resistant",
    )


def variant(**kwargs):
    kwargs.setdefault("gene", "gyrA")
    kwargs.setdefault("change", "D94G")
    return Variant.of(kwargs.pop("gene"), kwargs.pop("change"), **kwargs)


class AssessabilityTests(unittest.TestCase):
    """"Could not assess" must never render as "nothing found"."""

    def setUp(self):
        self.assessor = HeteroresistanceAssessor(depth_floor=10)

    def _one(self, v):
        return self.assessor.assess([v], [evidence_for(v)])[0]

    def test_unknown_platform_is_not_assessable(self):
        f = self._one(variant(vaf=0.15, depth=64, alt_depth=10))
        self.assertFalse(f.assessable)
        self.assertFalse(f.detected)
        self.assertIn("platform unknown", f.note)

    def test_missing_vaf_is_not_assessable(self):
        f = self._one(variant(depth=64, platform="illumina"))
        self.assertFalse(f.assessable)
        self.assertIn("no allele fraction", f.note)

    def test_low_depth_is_not_assessable(self):
        f = self._one(variant(vaf=0.15, depth=5, alt_depth=1,
                              platform="illumina"))
        self.assertFalse(f.assessable)
        self.assertIn("below floor", f.note)

    def test_too_few_alternate_reads_is_not_assessable(self):
        f = self._one(variant(vaf=0.05, depth=40, alt_depth=2,
                              platform="illumina"))
        self.assertFalse(f.assessable)
        self.assertIn("alternate read", f.note)

    def test_unconfigured_platform_is_not_assessable(self):
        f = self._one(variant(vaf=0.15, depth=64, alt_depth=10,
                              platform="some-new-sequencer"))
        self.assertFalse(f.assessable)
        self.assertIn("no limit of detection", f.note)


class DetectionTests(unittest.TestCase):
    def setUp(self):
        self.assessor = HeteroresistanceAssessor(depth_floor=10)

    def _one(self, v):
        return self.assessor.assess([v], [evidence_for(v)])[0]

    def test_minority_allele_above_lod_is_detected(self):
        f = self._one(variant(vaf=0.15, depth=64, alt_depth=10,
                              platform="illumina"))
        self.assertTrue(f.assessable)
        self.assertTrue(f.detected)
        self.assertIn("Minority allele at 15%", f.note)
        self.assertIn("Read-level confirmation", f.note)

    def test_fraction_below_platform_lod_is_not_detected(self):
        f = self._one(variant(vaf=0.06, depth=400, alt_depth=24,
                              platform="nanopore"))
        self.assertTrue(f.assessable)
        self.assertFalse(f.detected)
        self.assertIn("below the nanopore limit of detection", f.note)

    def test_same_fraction_differs_by_platform(self):
        illumina = self._one(variant(vaf=0.08, depth=400, alt_depth=32,
                                     platform="illumina"))
        nanopore = self._one(variant(vaf=0.08, depth=400, alt_depth=32,
                                     platform="nanopore"))
        self.assertTrue(illumina.detected,
                        "8% is above the Illumina limit of detection")
        self.assertFalse(nanopore.detected,
                         "8% is below the Nanopore limit of detection")

    def test_majority_allele_is_not_a_minority_population(self):
        f = self._one(variant(vaf=0.98, depth=80, alt_depth=78,
                              platform="illumina"))
        self.assertTrue(f.assessable)
        self.assertFalse(f.detected)
        self.assertIn("majority population", f.note)

    def test_lod_is_recorded_on_the_finding(self):
        f = self._one(variant(vaf=0.15, depth=64, alt_depth=10,
                              platform="illumina"))
        self.assertAlmostEqual(f.limit_of_detection, PLATFORM_LOD["illumina"])


class ScopeTests(unittest.TestCase):
    def test_indeterminate_evidence_is_also_assessed(self):
        v = variant(gene="Rv0678", change="L117R", vaf=0.2, depth=60,
                    alt_depth=12, platform="illumina")
        ev = DrugEvidence(
            drug="bedaquiline", call=Call.INDETERMINATE, tier=Tier.INFERRED,
            lane=Lane.EFFLUX_REGULATORY, variant=v.identity,
            rationale="efflux inference",
        )
        findings = HeteroresistanceAssessor().assess([v], [ev])
        self.assertEqual(len(findings), 1)
        self.assertTrue(findings[0].detected)

    def test_susceptible_evidence_is_not_assessed(self):
        v = variant(vaf=0.5, depth=60, alt_depth=30, platform="illumina")
        ev = DrugEvidence(
            drug="isoniazid", call=Call.SUSCEPTIBLE, tier=Tier.CATALOGUED,
            lane=Lane.CATALOGUE, variant=v.identity, rationale="not associated",
        )
        self.assertEqual(HeteroresistanceAssessor().assess([v], [ev]), [])

    def test_one_finding_per_variant_drug_pair(self):
        v = variant(vaf=0.2, depth=60, alt_depth=12, platform="illumina")
        findings = HeteroresistanceAssessor().assess(
            [v], [evidence_for(v), evidence_for(v)])
        self.assertEqual(len(findings), 1)

    def test_default_platform_is_applied(self):
        v = variant(vaf=0.15, depth=64, alt_depth=10)
        assessor = HeteroresistanceAssessor(default_platform="illumina")
        f = assessor.assess([v], [evidence_for(v)])[0]
        self.assertTrue(f.assessable)
        self.assertEqual(f.platform, "illumina")


if __name__ == "__main__":
    unittest.main()
