import tempfile
import unittest
from pathlib import Path

from myconductor.io.adapter import load


class AdapterTests(unittest.TestCase):
    def _write(self, name, text):
        d = Path(tempfile.mkdtemp())
        p = d / name
        p.write_text(text)
        return p

    def test_vcf_parsing_reads_info_and_format(self):
        vcf = (
            "##fileformat=VCFv4.2\n"
            "##sampleID=S1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr\t100\t.\tG\tC\t60\tPASS\tGENE=katG;AACHANGE=S315T;REGION=coding\tAF:DP\t0.9:80\n"
        )
        res = load(self._write("x.vcf", vcf))
        self.assertEqual(res.sample_id, "S1")
        self.assertEqual(len(res.variants), 1)
        v = res.variants[0]
        self.assertEqual(v.gene, "katG")
        self.assertEqual(v.change, "S315T")
        self.assertAlmostEqual(v.vaf, 0.9)
        self.assertEqual(v.depth, 80)

    def test_depth_floor_drops_low_coverage(self):
        vcf = (
            "##sampleID=S1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr\t1\t.\tA\tG\t60\tPASS\tGENE=rrs;AACHANGE=a1401g;REGION=rrna\tAF:DP\t1.0:4\n"
        )
        res = load(self._write("x.vcf", vcf), depth_floor=10)
        self.assertEqual(len(res.variants), 0)
        self.assertEqual(len(res.warnings), 1)

    def test_tsv_parsing(self):
        tsv = "gene\tchange\tvaf\tdepth\nkatG\tS315T\t1.0\t50\n"
        res = load(self._write("x.tsv", tsv))
        self.assertEqual(res.variants[0].gene, "katG")


if __name__ == "__main__":
    unittest.main()
