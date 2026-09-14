"""Importing precomputed AlphaFold models, without touching the network.

Two refusals carry this module, and both are grounded in measurements against
the live database rather than in caution.

**Numbering.** The WHO catalogue's rpoB Ser450 is residue *456* of UniProt
P9WGY9; position 450 there is a threonine, and the offset is a consistent +6
across Leu430, Ser431, Asp435, His445 and Ser450. Nothing about that mismatch
announces itself — every number in range resolves to some residue with some
pLDDT — so an unchecked mapping produces a confident, well-formatted annotation
of the wrong amino acid.

**Ligands.** AlphaFold DB entries are apo. There is no rifampicin in the rpoB
model, so a "distance to the rifampicin pocket" computed from these coordinates
would be a distance to something absent from the file.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from mycobench.alphafold import (
    NO_LIGAND_REASON, PLDDT_UNRELIABLE, AlphaFoldError, annotation_rows,
    fetch_structure, ligand_distance, model_url, parse_pdb, write_annotations,
)

ORGANISM = "Mycobacterium tuberculosis"
ASSEMBLY = "NC_000962.3"
WHEN = "2026-09-14T00:00:00Z"


def atom(number: int, name: str, x: float, plddt: float) -> str:
    """One alpha-carbon record in PDB column format."""
    return (f"ATOM  {number * 4:5d}  CA  {name:>3s} A{number:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{plddt:6.2f}           C")


def model(records=None) -> str:
    records = records or [
        atom(1, "MET", 0.0, 95.0),
        atom(2, "THR", 3.8, 97.3),
        atom(3, "SER", 7.6, 96.9),
        atom(4, "GLY", 11.4, 40.5),   # unreliable region
    ]
    return "\n".join(["HEADER    PREDICTED MODEL", *records, "END"])


class ParsingTests(unittest.TestCase):
    def test_alpha_carbons_and_plddt_are_read(self):
        structure = parse_pdb(model(), "P9WGY9", "v6")
        self.assertEqual(len(structure.residues), 4)
        self.assertEqual(structure.residue(2).name, "THR")
        self.assertAlmostEqual(structure.residue(2).plddt, 97.3)

    def test_confidence_bands_are_applied(self):
        structure = parse_pdb(model(), "P9WGY9", "v6")
        self.assertTrue(structure.residue(2).confident)
        self.assertFalse(structure.residue(4).confident)
        self.assertTrue(structure.residue(4).plddt < PLDDT_UNRELIABLE)

    def test_distance_is_euclidean_over_alpha_carbons(self):
        structure = parse_pdb(model(), "P9WGY9", "v6")
        self.assertAlmostEqual(
            structure.residue(1).distance_to(structure.residue(2)), 3.8, places=3)

    def test_mean_plddt_summarises_the_model(self):
        structure = parse_pdb(model(), "P9WGY9", "v6")
        self.assertAlmostEqual(structure.mean_plddt, (95.0 + 97.3 + 96.9 + 40.5) / 4)
        self.assertIn("mean pLDDT", structure.describe())

    def test_a_file_with_no_alpha_carbons_is_refused(self):
        with self.assertRaises(AlphaFoldError) as caught:
            parse_pdb("HEADER\nEND", "P9WGY9", "v6")
        self.assertIn("no alpha-carbon records", str(caught.exception))

    def test_a_malformed_record_is_refused_not_skipped(self):
        broken = "ATOM      1  CA  MET A   1     x.xxx   0.000   0.000  1.00 95.00"
        with self.assertRaises(AlphaFoldError):
            parse_pdb(broken, "P9WGY9", "v6")

    def test_empty_mean_for_an_empty_structure(self):
        from mycobench.alphafold import Structure
        self.assertIsNone(Structure("P0", "v6").mean_plddt)


class NumberingGuardTests(unittest.TestCase):
    """The +6 rpoB offset, reproduced in miniature."""

    def setUp(self):
        self.structure = parse_pdb(model(), "P9WGY9", "v6")

    def _rows(self, **variant):
        base = dict(variant_key="NC_000962.3:761155:C>T")
        base.update(variant)
        return annotation_rows(self.structure, [base], ORGANISM, ASSEMBLY, WHEN)

    def test_a_residue_that_is_not_what_the_caller_expected_is_refused(self):
        rows = self._rows(residue=2, expected_residue="Ser")
        self.assertIn("NUMBERING MISMATCH", rows[0]["predicted_effect"])
        self.assertIn("is THR, not the expected SER", rows[0]["predicted_effect"])
        self.assertIsNone(rows[0]["confidence"])

    def test_the_correct_residue_is_annotated(self):
        rows = self._rows(residue=3, expected_residue="Ser")
        self.assertNotIn("MISMATCH", rows[0]["predicted_effect"])
        self.assertAlmostEqual(rows[0]["confidence"], 0.969, places=3)

    def test_one_letter_and_three_letter_codes_both_work(self):
        for spelling in ("S", "SER", "Ser", "ser"):
            rows = self._rows(residue=3, expected_residue=spelling)
            self.assertNotIn("MISMATCH", rows[0]["predicted_effect"])

    def test_omitting_the_expectation_is_refused_outright(self):
        with self.assertRaises(AlphaFoldError) as caught:
            self._rows(residue=3)
        self.assertIn("expected_residue is required", str(caught.exception))
        self.assertIn("+6", str(caught.exception))

    def test_a_residue_beyond_the_model_is_reported_absent(self):
        rows = self._rows(residue=9999, expected_residue="Ser")
        self.assertIn("absent from the predicted model",
                      rows[0]["predicted_effect"])

    def test_a_variant_without_a_key_is_refused(self):
        with self.assertRaises(AlphaFoldError):
            annotation_rows(self.structure, [{"residue": 3}], ORGANISM,
                            ASSEMBLY, WHEN)


class LowConfidenceTests(unittest.TestCase):
    def test_an_unreliable_region_yields_no_interpretation(self):
        structure = parse_pdb(model(), "P9WGY9", "v6")
        rows = annotation_rows(
            structure, [{"variant_key": "k", "residue": 4,
                         "expected_residue": "Gly"}],
            ORGANISM, ASSEMBLY, WHEN)
        self.assertIn("does not support a structural interpretation",
                      rows[0]["predicted_effect"])
        # Reported as low confidence, not as a low score to be used anyway.
        self.assertAlmostEqual(rows[0]["confidence"], 0.405, places=3)


class LigandRefusalTests(unittest.TestCase):
    def test_ligand_distance_always_refuses_with_its_reason(self):
        with self.assertRaises(AlphaFoldError) as caught:
            ligand_distance()
        self.assertIn("predicted apo", str(caught.exception))
        self.assertEqual(str(caught.exception), NO_LIGAND_REASON)

    def test_emitted_rows_never_claim_a_ligand_distance(self):
        structure = parse_pdb(model(), "P9WGY9", "v6")
        rows = annotation_rows(
            structure, [{"variant_key": "k", "residue": 3,
                         "expected_residue": "Ser"}],
            ORGANISM, ASSEMBLY, WHEN)
        self.assertIsNone(rows[0]["ligand_distance"])
        self.assertIsNone(rows[0]["ligand_reference"])
        # And never claim a pocket, which an apo model cannot place.
        self.assertEqual(rows[0]["location"], "unknown")

    def test_caller_cited_reference_residues_are_attributed_to_the_caller(self):
        structure = parse_pdb(model(), "P9WGY9", "v6")
        rows = annotation_rows(
            structure, [{"variant_key": "k", "residue": 3,
                         "expected_residue": "Ser"}],
            ORGANISM, ASSEMBLY, WHEN,
            reference_residues={1: "Met1 (caller-cited)"})
        self.assertIn("as cited by the caller", rows[0]["predicted_effect"])
        self.assertIn("7.6 A", rows[0]["predicted_effect"])


class CacheTests(unittest.TestCase):
    def test_a_cached_model_is_read_without_the_network(self):
        with tempfile.TemporaryDirectory() as root:
            cache = Path(root)
            (cache / "AF-P9WGY9-F1-model_v6.pdb").write_text(model(),
                                                             encoding="utf-8")
            structure = fetch_structure("P9WGY9", cache, version="v6",
                                        cached_only=True)
        self.assertEqual(len(structure.residues), 4)

    def test_cached_only_without_a_version_is_refused(self):
        """Resolving the latest version needs the network, so say so."""
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaises(AlphaFoldError) as caught:
                fetch_structure("P9WGY9", root, cached_only=True)
        self.assertIn("requires the network", str(caught.exception))

    def test_a_missing_cached_model_is_refused(self):
        with tempfile.TemporaryDirectory() as root:
            with self.assertRaises(AlphaFoldError):
                fetch_structure("P00000", root, version="v6", cached_only=True)

    def test_the_url_carries_accession_and_version(self):
        url = model_url("p9wgy9", "v6")
        self.assertIn("AF-P9WGY9-F1-model_v6.pdb", url)

    def test_an_empty_accession_is_refused(self):
        with self.assertRaises(AlphaFoldError):
            model_url("  ", "v6")

    def test_annotations_write_as_loadable_json(self):
        from myconductor.modules.structural_annotation import (
            StructuralAnnotationTable,
        )
        structure = parse_pdb(model(), "P9WGY9", "v6")
        rows = annotation_rows(
            structure, [{"variant_key": "k", "residue": 3,
                         "expected_residue": "Ser"}],
            ORGANISM, ASSEMBLY, WHEN)
        with tempfile.TemporaryDirectory() as root:
            path = write_annotations(rows, Path(root) / "annotations.json")
            # The core loads it like any other attributed import.
            table = StructuralAnnotationTable.load(path)
        self.assertEqual(len(table.annotations), 1)
        self.assertEqual(table.annotations[0].source, "alphafold-db:P9WGY9")


if __name__ == "__main__":
    unittest.main()
