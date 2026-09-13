"""WHO catalogue ingest.

The fixtures use the **real column names** from the published files — the
master file's 114-column header and the coordinates file's five — so the
name-based lookup in ``catalogue.py`` is exercised against the spelling it will
actually meet. A subset of columns is enough precisely because the ingester
looks columns up by name rather than by position.
"""
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.catalogue import (
    CatalogueError,
    PINNED_COMMIT,
    REFERENCE_ASSEMBLY,
    coordinate_keys,
    ingest,
    read_coordinates,
    write_outputs,
)

POINTER = '(see "Genomic_coordinates" sheet)'

MASTER = "\t".join([
    "drug", "gene", "mutation", "variant", "tier", "effect",
    "genomic position", "FINAL CONFIDENCE GRADING", "Silent mutation",
    "Comment"]) + "\n" + "\n".join([
    "\t".join(["Rifampicin", "rpoB", "p.Ser450Leu", "rpoB_p.Ser450Leu", "1",
               "missense_variant", POINTER, "1) Assoc w R", "", "borderline"]),
    "\t".join(["Bedaquiline", "Rv0678", "p.Arg94Gln", "Rv0678_p.Arg94Gln", "1",
               "missense_variant", POINTER, "2) Assoc w R - interim", "", ""]),
    "\t".join(["Clofazimine", "Rv0678", "p.Arg94Gln", "Rv0678_p.Arg94Gln", "1",
               "missense_variant", POINTER, "2) Assoc w R - interim", "", ""]),
    "\t".join(["Isoniazid", "katG", "p.Arg463Leu", "katG_p.Arg463Leu", "1",
               "missense_variant", POINTER, "5) Not assoc w R", "", ""]),
    "\t".join(["Amikacin", "bacA", "c.102G>A", "bacA_c.102G>A", "2",
               "synonymous_variant", POINTER, "3) Uncertain significance",
               "True", ""]),
    "\t".join(["Rifampicin", "rpoB", "p.Zzz9Zzz", "rpoB_p.Zzz9Zzz", "1",
               "missense_variant", POINTER, "not a WHO grading", "", ""]),
]) + "\n"

COORDINATES = "\t".join([
    "variant", "chromosome", "position", "reference_nucleotide",
    "alternative_nucleotide"]) + "\n" + "\n".join([
    "rpoB_p.Ser450Leu\tNC_000962.3\t761155\tC\tT",
    # The same change also has an MNV spelling; both must be retained.
    "rpoB_p.Ser450Leu\tNC_000962.3\t761155\tCT\tTG",
    "Rv0678_p.Arg94Gln\tNC_000962.3\t779010\tG\tA",
    "rpoB_p.Bad\tNC_000962.3\tnot-a-number\tC\tT",
]) + "\n"


def write(name, text):
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text, encoding="utf-8")
    return path


class CoordinateTests(unittest.TestCase):
    def setUp(self):
        self.coords = read_coordinates(write("coords.txt", COORDINATES))

    def test_all_spellings_of_one_variant_are_retained(self):
        forms = self.coords["rpoB_p.Ser450Leu"]
        self.assertEqual(len(forms), 2,
                         "keeping only one spelling would make matching depend "
                         "on which form an annotator emitted")

    def test_non_numeric_position_is_dropped(self):
        self.assertNotIn("rpoB_p.Bad", self.coords)

    def test_coordinate_keys_are_myconductor_identities(self):
        keys = coordinate_keys(self.coords["rpoB_p.Ser450Leu"])
        self.assertIn("NC_000962.3:Chromosome:761155:C>T", keys)
        self.assertEqual(len(keys), 2)

    def test_missing_column_raises_with_guidance(self):
        with self.assertRaises(CatalogueError) as ctx:
            read_coordinates(write("bad.txt", "foo\tbar\n1\t2\n"))
        self.assertIn("missing column", str(ctx.exception))


class IngestTests(unittest.TestCase):
    def setUp(self):
        self.result = ingest(write("master.txt", MASTER),
                             write("coords.txt", COORDINATES),
                             provenance={"commit": PINNED_COMMIT})
        self.by_gene = {e["gene"]: e for e in self.result.entries}

    def test_one_entry_per_variant_with_drugs_merged(self):
        # Rv0678 p.Arg94Gln is graded for bedaquiline AND clofazimine.
        entry = self.by_gene["Rv0678"]
        self.assertEqual(sorted(entry["drugs"]),
                         ["bedaquiline", "clofazimine"])

    def test_grading_groups_map_to_calls(self):
        self.assertEqual(self.by_gene["rpoB"]["call"], "resistant")
        self.assertEqual(self.by_gene["Rv0678"]["call"], "resistant")
        self.assertEqual(self.by_gene["katG"]["call"], "not_associated")
        self.assertEqual(self.by_gene["bacA"]["call"], "indeterminate")

    def test_unrecognised_grading_is_skipped_and_reported(self):
        self.assertNotIn("rpoB_p.Zzz9Zzz",
                         {e["variant"] for e in self.result.entries})
        self.assertTrue(any("not one of WHO's five groups" in s
                            for s in self.result.skipped), self.result.skipped)

    def test_coordinates_are_attached_and_counted(self):
        entry = self.by_gene["rpoB"]
        self.assertEqual(entry["coordinate_key"],
                         "NC_000962.3:Chromosome:761155:C>T")
        self.assertEqual(len(entry["coordinate_keys"]), 2)
        # rpoB and Rv0678 have coordinates; katG and bacA do not.
        self.assertEqual(self.result.n_with_coordinates, 2)

    def test_entry_without_coordinates_has_no_key_rather_than_a_fake_one(self):
        self.assertNotIn("coordinate_key", self.by_gene["katG"])

    def test_confidence_is_labelled_as_ordinal(self):
        self.assertIn("not a probability",
                      self.by_gene["rpoB"]["confidence_note"])

    def test_silent_flag_is_carried(self):
        self.assertTrue(self.by_gene["bacA"]["silent"])
        self.assertFalse(self.by_gene["rpoB"]["silent"])

    def test_comment_is_preserved_when_present(self):
        self.assertEqual(self.by_gene["rpoB"]["comment"], "borderline")

    def test_drug_tiers_come_from_whos_own_tiering(self):
        self.assertIn("rpoB", self.result.drug_tiers["rifampicin"]["1"])
        self.assertIn("bacA", self.result.drug_tiers["amikacin"]["2"])

    def test_missing_master_column_raises_by_name(self):
        broken = MASTER.replace("FINAL CONFIDENCE GRADING", "SOMETHING ELSE")
        with self.assertRaises(CatalogueError) as ctx:
            ingest(write("m.txt", broken), write("c.txt", COORDINATES))
        self.assertIn("FINAL CONFIDENCE GRADING", str(ctx.exception))

    def test_all_rows_unusable_raises(self):
        header, _, _ = MASTER.partition("\n")
        only_bad = header + "\n" + "\t".join([
            "Rifampicin", "rpoB", "p.Q", "rpoB_p.Q", "1", "missense",
            POINTER, "nonsense", "", ""]) + "\n"
        with self.assertRaises(CatalogueError):
            ingest(write("m.txt", only_bad), write("c.txt", COORDINATES))


class OutputTests(unittest.TestCase):
    def setUp(self):
        self.result = ingest(write("master.txt", MASTER),
                             write("coords.txt", COORDINATES))

    def test_catalogue_json_is_not_marked_illustrative(self):
        payload = self.result.to_catalogue_json()
        self.assertFalse(payload["illustrative"])
        self.assertEqual(payload["reference_assembly"], REFERENCE_ASSEMBLY)
        self.assertIn("carry genomic coordinates", payload["note"])

    def test_drug_loci_json_keeps_lengths_null(self):
        # The catalogue supplies variant coordinates, not gene spans; inventing
        # a locus length would fabricate the evidence the mask demands.
        payload = self.result.to_drug_loci_json()
        self.assertTrue(payload["loci"])
        for locus, meta in payload["loci"].items():
            self.assertIsNone(meta["length_bp"], locus)

    def test_written_catalogue_loads_in_myconductor(self):
        from myconductor.modules.catalogue import CatalogueModule

        out = Path(tempfile.mkdtemp()) / "cat.json"
        write_outputs(self.result, out)
        module = CatalogueModule(out)
        self.assertFalse(module.is_illustrative)
        self.assertIn("rpoB_p.Ser450Leu", module.labels)

    def test_written_catalogue_resolves_a_coordinate_match(self):
        from myconductor.core.models import Consequence, Variant
        from myconductor.modules.catalogue import CatalogueModule

        out = Path(tempfile.mkdtemp()) / "cat.json"
        write_outputs(self.result, out)
        module = CatalogueModule(out)

        variant = Variant.of(
            "rpoB", "p.Ser450Leu", chrom="Chromosome", pos=761155,
            ref="C", alt="T", consequence=Consequence.MISSENSE)
        evidence = module.evaluate(variant)
        self.assertTrue(evidence)
        self.assertEqual(evidence[0].drug, "rifampicin")
        # A coordinate match must NOT carry the weaker label-match limitation.
        self.assertFalse(
            any("matched on gene/label" in limitation
                for limitation in evidence[0].limitations),
            evidence[0].limitations)

    def test_label_only_match_records_the_weaker_evidence(self):
        from myconductor.core.models import Variant
        from myconductor.modules.catalogue import CatalogueModule

        out = Path(tempfile.mkdtemp()) / "cat.json"
        write_outputs(self.result, out)
        module = CatalogueModule(out)

        evidence = module.evaluate(Variant.of("rpoB", "p.Ser450Leu"))
        self.assertTrue(evidence)
        self.assertTrue(
            any("matched on gene/label" in limitation
                for limitation in evidence[0].limitations))


if __name__ == "__main__":
    unittest.main()
