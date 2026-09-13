import tempfile
import unittest
from pathlib import Path

from myconductor.io.callable_mask import CallableMask, CallableMaskError


def write(name, text):
    path = Path(tempfile.mkdtemp()) / name
    path.write_text(text)
    return path


class AbsentMaskTests(unittest.TestCase):
    def test_absent_mask_is_not_present(self):
        mask = CallableMask.absent()
        self.assertFalse(mask.is_present)

    def test_absent_mask_never_passes_a_locus(self):
        ok, _, reasons = CallableMask.absent().assess(
            ["katG"], depth_floor=10, fraction_floor=0.95)
        self.assertFalse(ok)
        self.assertIn("no callable-locus evidence", reasons[0])

    def test_empty_locus_list_is_not_a_pass(self):
        mask = CallableMask.from_tsv(write(
            "m.tsv", "locus\tmean_depth\tcallable_fraction\nkatG\t80\t1.0\n"))
        ok, _, reasons = mask.assess([], depth_floor=10, fraction_floor=0.95)
        self.assertFalse(ok, "a drug with no declared loci must not pass")
        self.assertIn("no loci declared", reasons[0])


class TSVTests(unittest.TestCase):
    def test_explicit_fraction(self):
        mask = CallableMask.from_tsv(write(
            "m.tsv", "locus\tmean_depth\tcallable_fraction\nkatG\t82\t0.99\n"))
        ok, coverages, reasons = mask.assess(
            ["katG"], depth_floor=10, fraction_floor=0.95)
        self.assertTrue(ok, reasons)
        self.assertAlmostEqual(coverages[0].callable_fraction, 0.99)

    def test_fraction_computed_from_bp(self):
        mask = CallableMask.from_tsv(write(
            "m.tsv", "locus\tcovered_bp\ttotal_bp\nkatG\t2200\t2223\n"))
        cov = mask.coverage_for("katG")
        self.assertAlmostEqual(cov.callable_fraction, 2200 / 2223, places=4)

    def test_unknown_fraction_is_not_callable(self):
        mask = CallableMask.from_tsv(write(
            "m.tsv", "locus\tmean_depth\nkatG\t82\n"))
        ok, _, reasons = mask.assess(
            ["katG"], depth_floor=10, fraction_floor=0.95)
        self.assertFalse(ok, "unknown must never be treated as callable")
        self.assertIn("unknown", reasons[0])

    def test_below_fraction_floor_fails_with_the_number(self):
        mask = CallableMask.from_tsv(write(
            "m.tsv", "locus\tmean_depth\tcallable_fraction\nrrs\t40\t0.60\n"))
        ok, _, reasons = mask.assess(
            ["rrs"], depth_floor=10, fraction_floor=0.95)
        self.assertFalse(ok)
        self.assertIn("60%", reasons[0])

    def test_below_depth_floor_fails(self):
        mask = CallableMask.from_tsv(write(
            "m.tsv", "locus\tmean_depth\tcallable_fraction\nrrs\t4\t1.0\n"))
        ok, _, reasons = mask.assess(
            ["rrs"], depth_floor=10, fraction_floor=0.95)
        self.assertFalse(ok)
        self.assertIn("depth", reasons[0])

    def test_locus_absent_from_mask_is_reported(self):
        mask = CallableMask.from_tsv(write(
            "m.tsv", "locus\tmean_depth\tcallable_fraction\nkatG\t82\t1.0\n"))
        ok, _, reasons = mask.assess(
            ["rrs"], depth_floor=10, fraction_floor=0.95)
        self.assertFalse(ok)
        self.assertIn("absent", reasons[0])

    def test_missing_locus_column_is_an_error(self):
        with self.assertRaises(CallableMaskError):
            CallableMask.from_tsv(write("m.tsv", "gene\tdepth\nkatG\t80\n"))

    def test_out_of_range_fraction_is_an_error(self):
        with self.assertRaises(CallableMaskError):
            CallableMask.from_tsv(write(
                "m.tsv", "locus\tcallable_fraction\nkatG\t1.5\n"))


class BEDTests(unittest.TestCase):
    def test_bed_without_lengths_leaves_fraction_unknown(self):
        mask = CallableMask.from_bed(write(
            "m.bed", "NC_000962.3\t100\t200\tkatG\n"))
        cov = mask.coverage_for("katG")
        self.assertIsNone(cov.callable_fraction)
        ok, _, _ = mask.assess(["katG"], depth_floor=10, fraction_floor=0.95)
        self.assertFalse(ok)

    def test_bed_with_lengths_computes_fraction(self):
        mask = CallableMask.from_bed(
            write("m.bed", "NC_000962.3\t0\t99\tkatG\n"),
            locus_lengths={"katG": 100})
        self.assertAlmostEqual(mask.coverage_for("katG").callable_fraction, 0.99)

    def test_intervals_accumulate(self):
        mask = CallableMask.from_bed(
            write("m.bed", "c\t0\t50\tkatG\nc\t50\t100\tkatG\n"),
            locus_lengths={"katG": 100})
        self.assertAlmostEqual(mask.coverage_for("katG").callable_fraction, 1.0)

    def test_too_few_columns_is_an_error(self):
        with self.assertRaises(CallableMaskError):
            CallableMask.from_bed(write("m.bed", "c\t0\t100\n"))

    def test_inverted_interval_is_an_error(self):
        with self.assertRaises(CallableMaskError):
            CallableMask.from_bed(write("m.bed", "c\t200\t100\tkatG\n"))


class GVCFTests(unittest.TestCase):
    GVCF = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "c\t1\t.\tA\t<NON_REF>\t.\t.\tEND=100\tMIN_DP\t50\n"
    )

    def test_requires_locus_spans(self):
        with self.assertRaises(CallableMaskError) as ctx:
            CallableMask.from_gvcf(write("m.gvcf", self.GVCF), locus_spans={})
        self.assertIn("locus_spans", str(ctx.exception))

    def test_computes_fraction_over_a_span(self):
        mask = CallableMask.from_gvcf(
            write("m.gvcf", self.GVCF),
            locus_spans={"katG": (1, 100)}, depth_floor=10)
        self.assertAlmostEqual(mask.coverage_for("katG").callable_fraction, 1.0)

    def test_blocks_below_depth_do_not_count(self):
        gvcf = self.GVCF.replace("MIN_DP\t50", "MIN_DP\t5")
        mask = CallableMask.from_gvcf(
            write("m.gvcf", gvcf),
            locus_spans={"katG": (1, 100)}, depth_floor=10)
        self.assertAlmostEqual(mask.coverage_for("katG").callable_fraction, 0.0)


class MergeTests(unittest.TestCase):
    def test_known_fraction_beats_unknown(self):
        unknown = CallableMask.from_tsv(write(
            "a.tsv", "locus\tmean_depth\nkatG\t80\n"))
        known = CallableMask.from_tsv(write(
            "b.tsv", "locus\tcallable_fraction\nkatG\t0.99\n"))
        merged = CallableMask.merge(unknown, known)
        self.assertAlmostEqual(merged.coverage_for("katG").callable_fraction, 0.99)

    def test_higher_fraction_wins(self):
        low = CallableMask.from_tsv(write(
            "a.tsv", "locus\tcallable_fraction\nkatG\t0.70\n"))
        high = CallableMask.from_tsv(write(
            "b.tsv", "locus\tcallable_fraction\nkatG\t0.98\n"))
        merged = CallableMask.merge(low, high)
        self.assertAlmostEqual(merged.coverage_for("katG").callable_fraction, 0.98)

    def test_merging_absent_masks_yields_absent(self):
        merged = CallableMask.merge(CallableMask.absent(), CallableMask.absent())
        self.assertFalse(merged.is_present)

    def test_loci_from_both_masks_are_retained(self):
        a = CallableMask.from_tsv(write(
            "a.tsv", "locus\tcallable_fraction\nkatG\t0.99\n"))
        b = CallableMask.from_tsv(write(
            "b.tsv", "locus\tcallable_fraction\nrpoB\t0.98\n"))
        merged = CallableMask.merge(a, b)
        ok, _, reasons = merged.assess(
            ["katG", "rpoB"], depth_floor=10, fraction_floor=0.95)
        self.assertTrue(ok, reasons)


if __name__ == "__main__":
    unittest.main()
