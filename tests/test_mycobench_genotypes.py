"""Genotyping from CRyPTIC VCFs, tested without touching the network.

``fetch_vcf`` returns early when the cache already holds the file, so
pre-populating the cache exercises the real parsing path with no download.
"""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from mycobench.analysis import genotypes

CATALOGUE = {
    "catalogue_version": "test-2023.7",
    "illustrative": False,
    "variants": [
        {
            "gene": "rpoB", "change": "p.Ser450Leu",
            "variant": "rpoB_p.Ser450Leu", "drugs": ["rifampicin"],
            "call": "resistant",
            # WHO publishes both an MNV and an SNV spelling of this change.
            "coordinate_key": "NC_000962.3:Chromosome:761155:C>T",
            "coordinate_keys": ["NC_000962.3:Chromosome:761155:C>T",
                                "NC_000962.3:Chromosome:761155:CT>TG"],
        },
        {
            "gene": "Rv0678", "change": "p.Arg94Gln",
            "variant": "Rv0678_p.Arg94Gln",
            "drugs": ["bedaquiline", "clofazimine"], "call": "resistant",
            "coordinate_key": "NC_000962.3:Chromosome:779010:G>A",
            "coordinate_keys": ["NC_000962.3:Chromosome:779010:G>A"],
        },
        {
            # No coordinates: must be skipped rather than indexed on a guess.
            "gene": "katG", "change": "p.Arg463Leu",
            "variant": "katG_p.Arg463Leu", "drugs": ["isoniazid"],
            "call": "not_associated",
        },
    ],
}


def catalogue_file() -> Path:
    path = Path(tempfile.mkdtemp()) / "cat.json"
    path.write_text(json.dumps(CATALOGUE), encoding="utf-8")
    return path


VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def write_vcf(records: str, gz: bool = True, name: str = "s.vcf.gz") -> Path:
    path = Path(tempfile.mkdtemp()) / name
    payload = (VCF_HEADER + records).encode()
    if gz:
        path.write_bytes(gzip.compress(payload))
    else:
        path.write_bytes(payload)
    return path


class TrimTests(unittest.TestCase):
    def test_substitution_unchanged(self):
        self.assertEqual(genotypes._trim(100, "G", "C"), (100, "G", "C"))

    def test_equivalent_spellings_converge(self):
        self.assertEqual(genotypes._trim(100, "CTT", "C"),
                         genotypes._trim(100, "CTTT", "CT"))

    def test_shared_prefix_advances_position(self):
        self.assertEqual(genotypes._trim(100, "AG", "AT"), (101, "G", "T"))


class CoordinateIndexTests(unittest.TestCase):
    def setUp(self):
        self.index = genotypes.CoordinateIndex.from_catalogue(catalogue_file())

    def test_only_variants_with_coordinates_are_indexed(self):
        self.assertEqual(self.index.n_variants, 2)
        self.assertIn("test-2023.7", self.index.describe())

    def test_both_published_spellings_resolve_to_one_label(self):
        snv = self.index.lookup(761155, "C", "T")
        mnv = self.index.lookup(761155, "CT", "TG")
        self.assertEqual(snv, {"rpoB_p.Ser450Leu"})
        self.assertEqual(mnv, {"rpoB_p.Ser450Leu"})

    def test_lookup_tries_the_trimmed_spelling(self):
        # An isolate reporting the MNV form must still match, and the trimmed
        # form of CT>TG is C>T at the same position.
        self.assertIn("rpoB_p.Ser450Leu", self.index.lookup(761155, "CT", "TG"))

    def test_unknown_allele_returns_empty(self):
        self.assertEqual(self.index.lookup(999999, "A", "G"), set())

    def test_lowercase_alleles_match(self):
        self.assertEqual(self.index.lookup(761155, "c", "t"),
                         {"rpoB_p.Ser450Leu"})

    def test_malformed_key_is_skipped_not_fatal(self):
        payload = {"catalogue_version": "x", "variants": [
            {"gene": "g", "change": "c", "drugs": ["d"], "call": "resistant",
             "coordinate_keys": ["not-a-key"]}]}
        path = Path(tempfile.mkdtemp()) / "bad.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        index = genotypes.CoordinateIndex.from_catalogue(path)
        self.assertEqual(index.n_keys, 0)


class URLTests(unittest.TestCase):
    def test_relative_prefix_is_stripped(self):
        url = genotypes.vcf_url("../reproducibility/00/01/x.vcf.gz")
        self.assertTrue(url.startswith("https://"))
        self.assertIn("reproducibility/00/01/x.vcf.gz", url)
        self.assertNotIn("../", url)

    def test_multiple_relative_prefixes_are_stripped(self):
        self.assertNotIn("..", genotypes.vcf_url("../../a/b.vcf.gz"))


class ParseVCFTests(unittest.TestCase):
    def setUp(self):
        self.index = genotypes.CoordinateIndex.from_catalogue(catalogue_file())

    def test_matches_a_catalogued_variant(self):
        path = write_vcf(
            "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        result = genotypes.parse_vcf(path, self.index)
        self.assertEqual(result.variants, frozenset({"rpoB_p.Ser450Leu"}))
        self.assertEqual(result.n_matched, 1)

    def test_contig_naming_does_not_prevent_a_match(self):
        # The catalogue writes NC_000962.3; a VCF may write anything.
        path = write_vcf("Chromosome\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        result = genotypes.parse_vcf(path, self.index)
        self.assertIn("rpoB_p.Ser450Leu", result.variants)
        self.assertEqual(result.contigs, frozenset({"Chromosome"}))

    def test_uncatalogued_variant_is_counted_but_not_reported(self):
        path = write_vcf("NC_000962.3\t5\t.\tA\tG\t60\tPASS\t.\tGT\t1\n")
        result = genotypes.parse_vcf(path, self.index)
        self.assertEqual(result.variants, frozenset())
        self.assertEqual(result.n_records, 1)
        self.assertEqual(result.n_matched, 0)
        self.assertEqual(result.match_rate, 0.0)

    def test_filtered_records_are_excluded(self):
        path = write_vcf(
            "NC_000962.3\t761155\t.\tC\tT\t60\tMASKED\t.\tGT\t1\n")
        result = genotypes.parse_vcf(path, self.index)
        self.assertEqual(result.variants, frozenset())
        self.assertEqual(result.n_filtered_out, 1)

    def test_dot_filter_is_accepted(self):
        path = write_vcf("NC_000962.3\t761155\t.\tC\tT\t60\t.\t.\tGT\t1\n")
        self.assertIn("rpoB_p.Ser450Leu",
                      genotypes.parse_vcf(path, self.index).variants)

    def test_symbolic_and_absent_alts_are_skipped(self):
        path = write_vcf(
            "NC_000962.3\t761155\t.\tC\t<NON_REF>\t60\tPASS\t.\tGT\t1\n"
            "NC_000962.3\t761156\t.\tC\t.\t60\tPASS\t.\tGT\t1\n"
            "NC_000962.3\t761157\t.\tC\t*\t60\tPASS\t.\tGT\t1\n")
        result = genotypes.parse_vcf(path, self.index)
        self.assertEqual(result.n_records, 0)

    def test_only_the_genotyped_multiallelic_alt_is_present(self):
        path = write_vcf(
            "NC_000962.3\t761155\t.\tC\tA,T\t60\tPASS\t.\tGT\t1\n")
        result = genotypes.parse_vcf(path, self.index)
        self.assertEqual(result.n_records, 1)
        self.assertNotIn("rpoB_p.Ser450Leu", result.variants)
        self.assertIn("rpoB_p.Ser450Leu", result.assessed_variants)

    def test_multi_drug_variant_resolves_once(self):
        path = write_vcf("NC_000962.3\t779010\t.\tG\tA\t60\tPASS\t.\tGT\t1\n")
        result = genotypes.parse_vcf(path, self.index)
        self.assertEqual(result.variants, frozenset({"Rv0678_p.Arg94Gln"}))

    def test_plain_uncompressed_vcf_is_read(self):
        path = write_vcf("NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n",
                         gz=False, name="s.vcf")
        self.assertIn("rpoB_p.Ser450Leu",
                      genotypes.parse_vcf(path, self.index).variants)

    def test_non_numeric_position_is_skipped(self):
        path = write_vcf("NC_000962.3\tXYZ\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        self.assertEqual(genotypes.parse_vcf(path, self.index).n_records, 0)

    def test_corrupt_gzip_raises_a_clear_error(self):
        path = Path(tempfile.mkdtemp()) / "bad.vcf.gz"
        path.write_bytes(b"not gzip at all")
        with self.assertRaises(genotypes.GenotypeError):
            genotypes.parse_vcf(path, self.index)


class LoadGenotypesTests(unittest.TestCase):
    """Exercised through the cache, so no network is touched."""

    def setUp(self):
        self.index = genotypes.CoordinateIndex.from_catalogue(catalogue_file())
        self.cache = Path(tempfile.mkdtemp())

    def _cache(self, relative: str, records: str) -> None:
        destination = self.cache / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(gzip.compress((VCF_HEADER + records).encode()))

    def test_cached_file_is_used_without_a_download(self):
        self._cache("reproducibility/a.vcf.gz",
                    "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        path = genotypes.fetch_vcf("../reproducibility/a.vcf.gz", self.cache)
        self.assertTrue(path.is_file())

    def test_loads_isolates_and_keys_them_for_the_phenotype_table(self):
        self._cache("reproducibility/a.vcf.gz",
                    "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        rows = [{"ENA_RUN": "ERR1", "VCF": "../reproducibility/a.vcf.gz",
                 "REGENOTYPED_VCF": "", "UNIQUEID": "site.02.subj.0001"}]
        load = genotypes.load_genotypes(rows, self.index, self.cache,
                                        progress_every=0)
        self.assertEqual(load.n_loaded, 1)
        # Must match phenotypes.phenotype_table's sample_id convention.
        self.assertEqual(load.isolates["ERR1"].isolate_id, "cr_ERR1")
        self.assertEqual(load.isolates["ERR1"].genotype,
                         frozenset({"rpoB_p.Ser450Leu"}))
        self.assertEqual(load.isolates["ERR1"].site, "02")

    def test_regenotyped_vcf_is_preferred(self):
        self._cache("reproducibility/plain.vcf.gz",
                    "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        self._cache("reproducibility/regeno.vcf.gz",
                    "NC_000962.3\t779010\t.\tG\tA\t60\tPASS\t.\tGT\t1\n")
        rows = [{"ENA_RUN": "ERR1", "VCF": "../reproducibility/plain.vcf.gz",
                 "REGENOTYPED_VCF": "../reproducibility/regeno.vcf.gz"}]
        load = genotypes.load_genotypes(rows, self.index, self.cache,
                                        progress_every=0)
        self.assertEqual(load.isolates["ERR1"].genotype,
                         frozenset({"Rv0678_p.Arg94Gln"}))

    def test_a_missing_vcf_is_recorded_and_skipped_not_fatal(self):
        self._cache("reproducibility/a.vcf.gz",
                    "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        rows = [
            {"ENA_RUN": "ERR1", "VCF": "../reproducibility/a.vcf.gz"},
            {"ENA_RUN": "ERR2", "VCF": ""},
        ]
        load = genotypes.load_genotypes(rows, self.index, self.cache,
                                        progress_every=0)
        self.assertEqual(load.n_loaded, 1)
        self.assertEqual(len(load.failures), 1)
        self.assertIn("no VCF path", load.failures[0])

    def test_rows_without_a_run_accession_are_ignored(self):
        rows = [{"ENA_RUN": "None", "VCF": "../x.vcf.gz"},
                {"ENA_RUN": "", "VCF": "../y.vcf.gz"}]
        load = genotypes.load_genotypes(rows, self.index, self.cache,
                                        progress_every=0)
        self.assertEqual(load.n_requested, 0)

    def test_limit_caps_the_load(self):
        for name in ("a", "b", "c"):
            self._cache(f"reproducibility/{name}.vcf.gz",
                        "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        rows = [{"ENA_RUN": f"ERR{i}",
                 "VCF": f"../reproducibility/{n}.vcf.gz"}
                for i, n in enumerate("abc")]
        load = genotypes.load_genotypes(rows, self.index, self.cache, limit=2,
                                        progress_every=0)
        self.assertEqual(load.n_loaded, 2)

    def test_carriers_per_variant_is_counted(self):
        for name in ("a", "b"):
            self._cache(f"reproducibility/{name}.vcf.gz",
                        "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\t.\tGT\t1\n")
        rows = [{"ENA_RUN": f"ERR{i}",
                 "VCF": f"../reproducibility/{n}.vcf.gz"}
                for i, n in enumerate("ab")]
        load = genotypes.load_genotypes(rows, self.index, self.cache,
                                        progress_every=0)
        self.assertEqual(load.carriers_per_variant["rpoB_p.Ser450Leu"], 2)
        self.assertIn("distinct catalogued variant", load.describe())


if __name__ == "__main__":
    unittest.main()
