import tempfile
import unittest
from pathlib import Path

from myconductor.core.models import Consequence, Region
from myconductor.io.adapter import AdapterError, load, normalise_alleles


class NormalisationTests(unittest.TestCase):
    def test_substitution_is_unchanged(self):
        self.assertEqual(normalise_alleles(100, "G", "C"), (100, "G", "C"))

    def test_shared_suffix_is_trimmed(self):
        self.assertEqual(normalise_alleles(100, "CTT", "CT"), (100, "CT", "C"))

    def test_shared_prefix_advances_position(self):
        self.assertEqual(normalise_alleles(100, "AG", "AT"), (101, "G", "T"))

    def test_equivalent_deletions_converge(self):
        a = normalise_alleles(100, "CTT", "C")
        b = normalise_alleles(100, "CTTT", "CT")
        self.assertEqual(a, b)

    def test_at_least_one_base_remains(self):
        pos, ref, alt = normalise_alleles(100, "AAA", "AAA")
        self.assertTrue(ref and alt)


class VCFTests(unittest.TestCase):
    def _write(self, name, text):
        path = Path(tempfile.mkdtemp()) / name
        path.write_text(text)
        return path

    HEADER = (
        "##fileformat=VCFv4.2\n"
        "##sampleID=S1\n"
        "##reference=NC_000962.3\n"
        "##platform=illumina\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )

    def test_parses_info_format_and_platform(self):
        vcf = self.HEADER + (
            "NC_000962.3\t100\t.\tG\tC\t60\tPASS\t"
            "GENE=katG;AACHANGE=S315T;REGION=coding;EFFECT=missense\t"
            "AF:DP:AD\t0.9:80:8,72\n"
        )
        res = load(self._write("x.vcf", vcf))
        self.assertEqual(res.sample_id, "S1")
        self.assertEqual(res.platform, "illumina")
        self.assertEqual(len(res.variants), 1)
        v = res.variants[0]
        self.assertEqual(v.gene, "katG")
        self.assertAlmostEqual(v.vaf, 0.9)
        self.assertEqual(v.depth, 80)
        self.assertEqual(v.alt_depth, 72)
        self.assertIs(v.consequence, Consequence.MISSENSE)
        self.assertTrue(v.identity.is_coordinate_resolved)

    def test_multiallelic_record_becomes_two_variants(self):
        vcf = self.HEADER + (
            "NC_000962.3\t100\t.\tG\tC,T\t60\tPASS\tGENE=rpoB\t"
            "AF:DP\t0.6,0.3:90\n"
        )
        res = load(self._write("x.vcf", vcf))
        self.assertEqual(len(res.variants), 2)
        self.assertAlmostEqual(res.variants[0].vaf, 0.6)
        self.assertAlmostEqual(res.variants[1].vaf, 0.3)
        self.assertNotEqual(res.variants[0].key(), res.variants[1].key())

    def test_multi_sample_vcf_refuses_to_guess(self):
        vcf = (
            "##sampleID=S1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chr\t100\t.\tG\tC\t60\tPASS\tGENE=katG\tAF:DP\t1.0:80\t1.0:80\n"
        )
        with self.assertRaises(AdapterError) as ctx:
            load(self._write("x.vcf", vcf))
        self.assertIn("multi-sample", str(ctx.exception))

    def test_named_sample_is_selected(self):
        vcf = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chr\t100\t.\tG\tC\t60\tPASS\tGENE=katG\tAF:DP\t0.5:80\t0.9:90\n"
        )
        res = load(self._write("x.vcf", vcf), sample="S2")
        self.assertEqual(res.sample_id, "S2")
        self.assertAlmostEqual(res.variants[0].vaf, 0.9)

    def test_symbolic_alt_is_rejected_not_coerced(self):
        vcf = self.HEADER + (
            "NC_000962.3\t100\t.\tG\t<DEL>\t60\tPASS\tGENE=katG\tAF:DP\t1.0:80\n"
        )
        res = load(self._write("x.vcf", vcf))
        self.assertEqual(len(res.variants), 0)
        self.assertEqual(len(res.rejected), 1)
        self.assertIn("symbolic", res.rejected[0].reason)

    def test_filter_failing_record_is_excluded(self):
        vcf = self.HEADER + (
            "NC_000962.3\t100\t.\tG\tC\t60\tLowQual\tGENE=katG\tAF:DP\t1.0:80\n"
        )
        res = load(self._write("x.vcf", vcf))
        self.assertEqual(len(res.variants), 0)
        self.assertTrue(any("FILTER" in w for w in res.warnings))

    def test_low_depth_warning_says_it_is_not_susceptibility(self):
        vcf = self.HEADER + (
            "NC_000962.3\t1\t.\tA\tG\t60\tPASS\tGENE=rrs;REGION=rrna\tAF:DP\t1.0:4\n"
        )
        res = load(self._write("x.vcf", vcf), depth_floor=10)
        self.assertEqual(len(res.variants), 0)
        self.assertTrue(any("NOT_ASSESSED" in w for w in res.warnings))

    def test_synonymy_is_derived_not_trusted(self):
        vcf = self.HEADER + (
            "NC_000962.3\t100\t.\tC\tT\t60\tPASS\t"
            "GENE=rpoC;AACHANGE=A542A;REGION=coding\tAF:DP\t1.0:70\n"
        )
        res = load(self._write("x.vcf", vcf))
        self.assertTrue(res.variants[0].silent)

    def test_promoter_region_is_derived_from_hgvs_c(self):
        vcf = self.HEADER + (
            "NC_000962.3\t100\t.\tC\tT\t60\tPASS\t"
            "GENE=eis;AACHANGE=c.-14C>T\tAF:DP\t1.0:66\n"
        )
        res = load(self._write("x.vcf", vcf))
        self.assertIs(res.variants[0].region, Region.PROMOTER)

    def test_missing_chrom_header_is_an_error(self):
        with self.assertRaises(AdapterError):
            load(self._write("x.vcf", "##fileformat=VCFv4.2\n"))

    def test_platform_override_applies_to_variants(self):
        vcf = self.HEADER + (
            "NC_000962.3\t100\t.\tG\tC\t60\tPASS\tGENE=katG\tAF:DP\t1.0:80\n"
        )
        res = load(self._write("x.vcf", vcf), platform="nanopore")
        self.assertEqual(res.variants[0].platform, "nanopore")


class TSVJSONTests(unittest.TestCase):
    def _write(self, name, text):
        path = Path(tempfile.mkdtemp()) / name
        path.write_text(text)
        return path

    def test_tsv_parsing(self):
        tsv = ("gene\tchange\tvaf\tdepth\tplatform\n"
               "katG\tS315T\t1.0\t50\tillumina\n")
        res = load(self._write("x.tsv", tsv))
        self.assertEqual(res.variants[0].gene, "katG")
        self.assertEqual(res.variants[0].platform, "illumina")

    def test_tsv_missing_required_column(self):
        with self.assertRaises(AdapterError):
            load(self._write("x.tsv", "gene\tvaf\nkatG\t1.0\n"))

    def test_json_parsing_keeps_coordinates(self):
        payload = (
            '{"sample_id": "S9", "platform": "illumina", "variants": ['
            '{"gene": "rpoB", "change": "S450L", "chrom": "NC_000962.3",'
            ' "pos": 761155, "ref": "C", "alt": "T", "depth": 70}]}'
        )
        res = load(self._write("x.json", payload))
        self.assertEqual(res.sample_id, "S9")
        self.assertTrue(res.variants[0].identity.is_coordinate_resolved)

    def test_unsupported_extension(self):
        with self.assertRaises(AdapterError):
            load(self._write("x.bam", "binary"))


if __name__ == "__main__":
    unittest.main()
