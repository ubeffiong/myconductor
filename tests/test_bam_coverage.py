import gzip
import tempfile
import unittest
from pathlib import Path

from myconductor.io.bam_coverage import mask_from_bam
from myconductor.io.callable_mask import CallableMaskError


class FakeMosdepthRunner:
    """Simulates mosdepth by writing canned output files, per its documented
    ``--by``/``--thresholds`` schema. Lets the parser be tested without a real
    mosdepth binary."""

    def __init__(self, available=("mosdepth",), exit_code=0, depth_floor=10):
        self._available = set(available)
        self.exit_code = exit_code
        self.depth_floor = depth_floor

    def available(self, executable):
        return executable in self._available

    def run(self, command, timeout=None):
        if self.exit_code != 0:
            return self.exit_code
        prefix = Path(command[-2])
        with gzip.open(f"{prefix}.regions.bed.gz", "wt") as f:
            f.write("Chromosome\t0\t100\trpoB\t42.5\n")
            f.write("Chromosome\t200\t260\tkatG\t5.0\n")
        with gzip.open(f"{prefix}.thresholds.bed.gz", "wt") as f:
            f.write(f"chrom\tstart\tend\tregion\t{self.depth_floor}X\n")
            f.write("Chromosome\t0\t100\trpoB\t95\n")
            f.write("Chromosome\t200\t260\tkatG\t10\n")
        return 0

    def run_capturing_stdout(self, command, out_path, timeout=None):
        raise NotImplementedError


class FakeSamtoolsRunner:
    def __init__(self, available=("samtools",), exit_code=0, lines=()):
        self._available = set(available)
        self.exit_code = exit_code
        self.lines = lines

    def available(self, executable):
        return executable in self._available

    def run(self, command, timeout=None):
        raise NotImplementedError

    def run_capturing_stdout(self, command, out_path, timeout=None):
        if self.exit_code != 0:
            return self.exit_code
        Path(out_path).write_text("\n".join(self.lines) + ("\n" if self.lines else ""),
                                  encoding="utf-8")
        return 0


class BamCoverageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.bed = Path(self.tmp.name) / "targets.bed"
        self.bed.write_text(
            "Chromosome\t0\t100\trpoB\nChromosome\t200\t260\tkatG\n",
            encoding="utf-8")
        self.bam = Path(self.tmp.name) / "sample.bam"
        self.bam.write_text("not a real bam", encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def test_mosdepth_path_produces_a_usable_mask(self):
        mask = mask_from_bam(self.bam, self.bed, depth_floor=10,
                             runner=FakeMosdepthRunner(depth_floor=10))
        self.assertTrue(mask.is_present)
        callable_ok, coverages, reasons = mask.assess(
            ["rpoB"], depth_floor=10, fraction_floor=0.90)
        self.assertTrue(callable_ok, reasons)
        cov = mask.coverage_for("rpoB")
        self.assertEqual(cov.mean_depth, 42)
        self.assertAlmostEqual(cov.callable_fraction, 0.95)

    def test_mosdepth_low_threshold_locus_is_not_callable(self):
        mask = mask_from_bam(self.bam, self.bed, depth_floor=10,
                             runner=FakeMosdepthRunner(depth_floor=10))
        callable_ok, _coverages, reasons = mask.assess(
            ["katG"], depth_floor=10, fraction_floor=0.90)
        self.assertFalse(callable_ok)
        self.assertTrue(reasons)

    def test_samtools_fallback_is_used_when_mosdepth_absent(self):
        lines = [f"Chromosome\t{pos}\t20" for pos in range(1, 101)]
        mask = mask_from_bam(
            self.bam, self.bed, depth_floor=10,
            runner=FakeSamtoolsRunner(available=("samtools",), lines=lines))
        cov = mask.coverage_for("rpoB")
        self.assertEqual(cov.source, "bam:samtools")
        self.assertEqual(cov.callable_fraction, 1.0)

    def test_neither_tool_available_raises_actionably(self):
        with self.assertRaises(CallableMaskError) as ctx:
            mask_from_bam(self.bam, self.bed,
                         runner=FakeMosdepthRunner(available=()))
        self.assertIn("mosdepth", str(ctx.exception))

    def test_tool_failure_raises_rather_than_returning_an_absent_mask(self):
        with self.assertRaises(CallableMaskError):
            mask_from_bam(self.bam, self.bed,
                         runner=FakeMosdepthRunner(exit_code=1))

    def test_explicit_tool_choice_is_honoured(self):
        with self.assertRaises(CallableMaskError):
            mask_from_bam(self.bam, self.bed, tool="samtools",
                         runner=FakeMosdepthRunner(available=("mosdepth",)))


if __name__ == "__main__":
    unittest.main()
