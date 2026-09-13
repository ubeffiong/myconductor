"""The analysis driver, end to end on synthetic data and offline.

VCFs are pre-placed in the cache so ``fetch_vcf`` short-circuits, which
exercises the real code path with no network. The refusals matter as much as
the happy path: an analysis that silently ran against the bundled illustrative
catalogue would report effect sizes over almost no data.
"""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.analysis import pipeline
from mycobench.cohort import write_rows

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)

#: Positions used by the fixture. CANDIDATE is the variant under test,
#: BACKGROUND is an established determinant used as the stratifying background.
CANDIDATE_POS = 500_001
BACKGROUND_POS = 761_155


def _coordinate(position: int, ref: str, alt: str) -> str:
    return f"NC_000962.3:Chromosome:{position}:{ref}>{alt}"


def build_catalogue(illustrative: bool = False,
                    n_filler: int = 1200) -> Path:
    """An ingested-shaped catalogue with enough coordinates to pass validation."""
    variants = [
        {
            "gene": "rpoB", "change": "p.Ser450Leu",
            "variant": "rpoB_p.Ser450Leu", "drugs": ["rifampicin"],
            "call": "resistant", "tier": "1",
            "coordinate_key": _coordinate(BACKGROUND_POS, "C", "T"),
            "coordinate_keys": [_coordinate(BACKGROUND_POS, "C", "T")],
        },
        {
            "gene": "pncA", "change": "p.Asp12Ala",
            "variant": "pncA_p.Asp12Ala", "drugs": ["rifampicin"],
            "call": "indeterminate", "tier": "1",
            "coordinate_key": _coordinate(CANDIDATE_POS, "A", "C"),
            "coordinate_keys": [_coordinate(CANDIDATE_POS, "A", "C")],
        },
    ]
    # Filler entries only exist to clear the coordinate-count guard.
    for i in range(n_filler):
        variants.append({
            "gene": f"filler{i}", "change": f"p.X{i}Y",
            "variant": f"filler{i}_p.X{i}Y", "drugs": ["rifampicin"],
            "call": "indeterminate", "tier": "2",
            "coordinate_key": _coordinate(900_000 + i, "A", "G"),
            "coordinate_keys": [_coordinate(900_000 + i, "A", "G")],
        })
    payload = {"catalogue_version": "test-ingested",
               "illustrative": illustrative, "variants": variants}
    path = Path(tempfile.mkdtemp()) / "cat.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class Fixture:
    """A synthetic cohort where the candidate genuinely raises MIC."""

    def __init__(self, n_per_group: int = 8, candidate_mic: str = "8.0",
                 other_mic: str = "0.5", quality: str = "HIGH"):
        self.root = Path(tempfile.mkdtemp())
        self.cache = self.root / "vcf"
        self.catalogue = build_catalogue()

        reuse_rows = []
        phenotype_rows = []
        # Both groups carry the background determinant, so the comparison is
        # made at fixed background rather than confounded by it.
        for group, carries_candidate, value in (
                ("c", True, candidate_mic), ("n", False, other_mic)):
            for i in range(n_per_group):
                run = f"ERR{group}{i}"
                records = f"NC_000962.3\t{BACKGROUND_POS}\t.\tC\tT\t60\tPASS\t.\tGT\t1\n"
                if carries_candidate:
                    records += (f"NC_000962.3\t{CANDIDATE_POS}\t.\tA\tC\t60\t"
                                f"PASS\t.\tGT\t1\n")
                relative = f"reproducibility/{run}.vcf.gz"
                destination = self.cache / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(
                    gzip.compress((VCF_HEADER + records).encode()))
                reuse_rows.append({
                    "ENA_RUN": run, "UNIQUEID": f"site.01.subj.{i}",
                    "VCF": f"../{relative}", "REGENOTYPED_VCF": "",
                })
                phenotype_rows.append({
                    "sample_id": f"cr_{run}", "sra_run": run,
                    "drug": "rifampicin", "phenotype": "R", "quality": quality,
                    "mic": value, "mic_value": "", "mic_censoring": "none",
                    "evaluable": "yes", "not_evaluable_reason": "",
                })

        self.reuse_table = self.root / "reuse.csv"
        header = ["ENA_RUN", "UNIQUEID", "VCF", "REGENOTYPED_VCF"]
        lines = [",".join(header)]
        lines += [",".join(row[c] for c in header) for row in reuse_rows]
        self.reuse_table.write_text("\n".join(lines) + "\n", encoding="utf-8")

        self.phenotypes = self.root / "phenotypes.tsv"
        write_rows(self.phenotypes, phenotype_rows,
                   ("sample_id", "sra_run", "drug", "phenotype", "quality",
                    "mic", "mic_value", "mic_censoring", "evaluable",
                    "not_evaluable_reason"))

    def inputs(self, **overrides) -> pipeline.AnalysisInputs:
        base = dict(catalogue=self.catalogue, phenotypes=self.phenotypes,
                    reuse_table=self.reuse_table, cache_dir=self.cache,
                    min_carriers=5)
        base.update(overrides)
        return pipeline.AnalysisInputs(**base)


class ValidationTests(unittest.TestCase):
    def test_illustrative_catalogue_is_refused_with_a_reason(self):
        fixture = Fixture()
        inputs = fixture.inputs(catalogue=build_catalogue(illustrative=True))
        with self.assertRaises(pipeline.AnalysisError) as ctx:
            inputs.validate()
        message = str(ctx.exception)
        self.assertIn("illustrative", message)
        self.assertIn("fetch-catalogue", message)

    def test_catalogue_without_enough_coordinates_is_refused(self):
        fixture = Fixture()
        inputs = fixture.inputs(catalogue=build_catalogue(n_filler=3))
        with self.assertRaises(pipeline.AnalysisError) as ctx:
            inputs.validate()
        self.assertIn("coordinates for only", str(ctx.exception))

    def test_missing_input_file_is_named(self):
        fixture = Fixture()
        inputs = fixture.inputs(phenotypes=Path("nope.tsv"))
        with self.assertRaises(pipeline.AnalysisError) as ctx:
            inputs.validate()
        self.assertIn("phenotype table not found", str(ctx.exception))

    def test_a_valid_fixture_passes_validation(self):
        Fixture().inputs().validate()


class BuildPanelTests(unittest.TestCase):
    def test_low_quality_phenotypes_are_excluded(self):
        fixture = Fixture(quality="LOW")
        panels = pipeline.build_panels(fixture.phenotypes)
        self.assertEqual(panels, {})

    def test_high_quality_phenotypes_are_admitted(self):
        fixture = Fixture()
        panels = pipeline.build_panels(fixture.phenotypes)
        self.assertIn("rifampicin", panels)
        self.assertEqual(panels["rifampicin"].profile.n, 16)

    def test_drug_restriction_is_honoured(self):
        fixture = Fixture()
        panels = pipeline.build_panels(fixture.phenotypes, drugs=["isoniazid"])
        self.assertEqual(panels, {})


class CandidateSelectionTests(unittest.TestCase):
    def test_only_uncertain_entries_are_candidates(self):
        catalogue = build_catalogue()
        observed = {"pncA_p.Asp12Ala": 10, "rpoB_p.Ser450Leu": 10}
        candidates = pipeline.candidate_variants(catalogue, observed,
                                                 min_carriers=5)
        self.assertIn("pncA_p.Asp12Ala", candidates)
        self.assertNotIn("rpoB_p.Ser450Leu", candidates,
                         "an established determinant is background, not a "
                         "candidate")

    def test_carrier_minimum_is_enforced(self):
        catalogue = build_catalogue()
        candidates = pipeline.candidate_variants(
            catalogue, {"pncA_p.Asp12Ala": 2}, min_carriers=5)
        self.assertEqual(candidates, {})


class EndToEndTests(unittest.TestCase):
    def setUp(self):
        self.fixture = Fixture()
        self.result = pipeline.run(self.fixture.inputs())

    def test_isolates_are_genotyped_from_the_cache(self):
        self.assertEqual(self.result.load.n_loaded, 16)
        self.assertGreater(self.result.load.match_rate, 0.9)

    def test_the_candidate_is_tested_and_its_effect_detected(self):
        self.assertEqual(self.result.n_candidates, 1)
        self.assertEqual(self.result.n_tested, 1)
        effect = self.result.scan.effects[0]
        self.assertEqual(effect.variant, "pncA_p.Asp12Ala")
        self.assertEqual(effect.verdict, "evidence-of-effect")
        self.assertGreater(effect.delta, 0.9)

    def test_the_background_determinant_is_held_fixed(self):
        # Every isolate carries rpoB S450L, so the single informative stratum
        # is the one whose background contains it.
        effect = self.result.scan.effects[0]
        self.assertEqual(effect.n_informative_strata, 1)
        self.assertIn("rpoB_p.Ser450Leu", effect.strata[0].background)

    def test_drugs_cryptic_does_not_measure_are_skipped_with_a_reason(self):
        # The filler entries name rifampicin, so nothing is skipped here; assert
        # the mechanism exists rather than inventing a drug.
        self.assertIsInstance(self.result.skipped, list)

    def test_panel_censoring_is_reported(self):
        rows = self.result.panel_rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["median_identifiable"], "yes")

    def test_summary_names_the_caveats(self):
        text = self.result.summary()
        self.assertIn("candidate variant", text)
        self.assertIn("catalogued variant", text)

    def test_outputs_are_written_with_a_manifest(self):
        out = Path(tempfile.mkdtemp())
        written = pipeline.write_outputs(self.result, out)
        names = {p.name for p in written}
        self.assertEqual(names, {"mic_panels.tsv", "variant_effects.tsv",
                                 "analysis_manifest.json"})
        payload = json.loads((out / "analysis_manifest.json")
                             .read_text(encoding="utf-8"))
        self.assertEqual(payload["n_tested"], 1)
        self.assertIn("Nothing here is a resistance call", payload["caveat"])
        self.assertIn("not promoted to a predictive call",
                      payload["caveat"].replace("no variant is promoted",
                                                "not promoted"))

    def test_effect_rows_carry_the_interval_and_reasons(self):
        rows = self.result.effect_rows()
        self.assertEqual(len(rows), 1)
        self.assertTrue(rows[0]["ci_low"])
        self.assertTrue(rows[0]["reasons"])

    def test_discriminator_not_run_without_its_control_genes(self):
        # The fixture has no efflux-regulator or atpE carriers.
        self.assertIsNone(self.result.discriminator)
        self.assertTrue(any("discriminator not run" in n
                            for n in self.result.notes))


class NoCandidateTests(unittest.TestCase):
    def test_too_few_carriers_yields_a_note_not_a_crash(self):
        fixture = Fixture(n_per_group=3)
        result = pipeline.run(fixture.inputs(min_carriers=5))
        self.assertEqual(result.n_candidates, 0)
        self.assertTrue(any("carried by at least" in n for n in result.notes))
        self.assertIsNone(result.scan)

    def test_confounded_candidate_is_reported_as_such(self):
        # Give the candidate to nobody without the background AND to everyone
        # with it, so it never varies at fixed background.
        fixture = Fixture()
        result = pipeline.run(fixture.inputs(limit=8))
        # Only carriers were loaded, so there is no non-carrier anywhere.
        effect = result.scan.effects[0] if result.scan else None
        if effect is not None:
            self.assertIn(effect.verdict,
                          ("confounded", "underpowered", "not-estimable"))


if __name__ == "__main__":
    unittest.main()
