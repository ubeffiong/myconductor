import tempfile
import unittest
from pathlib import Path

from myconductor.io.read_calling import (
    ReadCallingError,
    align_reads,
    call_variants,
    reads_to_vcf,
)


class FakeReadRunner:
    """Simulates every external tool this module shells out to, by writing a
    plausible output file wherever the real tool would, so the orchestration
    logic (argument construction, failure propagation, output verification)
    can be tested without any of these tools installed."""

    def __init__(self, available=(), exit_code=0):
        self._available = set(available)
        self.exit_code = exit_code
        self.commands: list[list[str]] = []

    def available(self, executable):
        return executable in self._available

    def _out_after(self, command, flag):
        return Path(command[command.index(flag) + 1])

    def run(self, command, timeout=None):
        self.commands.append(list(command))
        if self.exit_code != 0:
            return self.exit_code
        prog, sub = command[0], command[1] if len(command) > 1 else ""
        if prog == "bwa-mem2" and sub == "index":
            Path(f"{command[2]}.bwt.2bit.64").write_bytes(b"fake-index")
        elif prog == "samtools" and sub == "faidx":
            Path(f"{command[2]}.fai").write_text("chr\t100\n", encoding="utf-8")
        elif prog == "samtools" and sub == "sort":
            self._out_after(command, "-o").write_bytes(b"fake-bam")
        elif prog == "samtools" and sub == "index":
            Path(f"{command[2]}.bai").write_bytes(b"fake-bai")
        elif prog == "bcftools" and sub == "mpileup":
            self._out_after(command, "-o").write_bytes(b"fake-bcf")
        elif prog == "bcftools" and sub == "call":
            self._out_after(command, "-o").write_bytes(b"fake-vcf")
        elif prog == "gatk" and sub == "HaplotypeCaller":
            self._out_after(command, "-O").write_bytes(b"fake-vcf")
        return 0

    def run_capturing_stdout(self, command, out_path, timeout=None):
        self.commands.append(list(command))
        if self.exit_code != 0:
            return self.exit_code
        Path(out_path).write_text("@HD\tVN:1.6\nFAKE\n", encoding="utf-8")
        return 0


class AlignReadsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.r1, self.r2 = d / "r1.fastq", d / "r2.fastq"
        self.r1.write_text("@x\nACGT\n+\nIIII\n", encoding="utf-8")
        self.r2.write_text("@x\nACGT\n+\nIIII\n", encoding="utf-8")
        self.ref = d / "ref.fa"
        self.ref.write_text(">chr\nACGTACGT\n", encoding="utf-8")
        self.out_bam = d / "out.bam"

    def tearDown(self):
        self.tmp.cleanup()

    def test_minimap2_is_preferred_when_available(self):
        runner = FakeReadRunner(available=("minimap2", "samtools"))
        result = align_reads(self.r1, self.r2, self.ref, self.out_bam,
                             runner=runner)
        self.assertEqual(result, self.out_bam)
        self.assertTrue(self.out_bam.is_file())
        self.assertTrue(Path(f"{self.out_bam}.bai").is_file())
        self.assertEqual(runner.commands[0][0], "minimap2")

    def test_bwa_mem2_is_used_when_minimap2_absent(self):
        runner = FakeReadRunner(available=("bwa-mem2", "samtools"))
        align_reads(self.r1, self.r2, self.ref, self.out_bam, runner=runner)
        self.assertTrue(Path(f"{self.ref}.bwt.2bit.64").is_file())
        progs = [c[0] for c in runner.commands]
        self.assertIn("bwa-mem2", progs)

    def test_bwa_index_is_not_rebuilt_when_already_present(self):
        Path(f"{self.ref}.bwt.2bit.64").write_bytes(b"already-indexed")
        runner = FakeReadRunner(available=("bwa-mem2", "samtools"))
        align_reads(self.r1, self.r2, self.ref, self.out_bam, runner=runner)
        index_commands = [c for c in runner.commands
                          if c[0] == "bwa-mem2" and c[1] == "index"]
        self.assertEqual(index_commands, [])

    def test_no_aligner_available_raises(self):
        runner = FakeReadRunner(available=("samtools",))
        with self.assertRaises(ReadCallingError):
            align_reads(self.r1, self.r2, self.ref, self.out_bam, runner=runner)

    def test_aligner_failure_raises_rather_than_a_partial_bam(self):
        runner = FakeReadRunner(available=("minimap2", "samtools"), exit_code=1)
        with self.assertRaises(ReadCallingError):
            align_reads(self.r1, self.r2, self.ref, self.out_bam, runner=runner)
        self.assertFalse(self.out_bam.is_file())


class CallVariantsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        self.ref = d / "ref.fa"
        self.ref.write_text(">chr\nACGTACGT\n", encoding="utf-8")
        self.bam = d / "sample.bam"
        self.bam.write_bytes(b"fake-bam")
        self.out_vcf = d / "out.vcf.gz"

    def tearDown(self):
        self.tmp.cleanup()

    def test_bcftools_is_preferred_when_available(self):
        runner = FakeReadRunner(available=("bcftools", "samtools"))
        result = call_variants(self.bam, self.ref, self.out_vcf, runner=runner)
        self.assertEqual(result, self.out_vcf)
        self.assertTrue(self.out_vcf.is_file())
        self.assertTrue(Path(f"{self.ref}.fai").is_file())

    def test_gatk_fallback_when_bcftools_absent(self):
        runner = FakeReadRunner(available=("gatk", "samtools"))
        call_variants(self.bam, self.ref, self.out_vcf, runner=runner)
        progs = [c[0] for c in runner.commands]
        self.assertIn("gatk", progs)

    def test_no_caller_available_raises(self):
        runner = FakeReadRunner(available=("samtools",))
        with self.assertRaises(ReadCallingError):
            call_variants(self.bam, self.ref, self.out_vcf, runner=runner)

    def test_caller_failure_raises(self):
        runner = FakeReadRunner(available=("bcftools", "samtools"), exit_code=1)
        with self.assertRaises(ReadCallingError):
            call_variants(self.bam, self.ref, self.out_vcf, runner=runner)


class ReadsToVcfTests(unittest.TestCase):
    def test_end_to_end_convenience_wrapper(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            r1, r2 = d / "r1.fastq", d / "r2.fastq"
            r1.write_text("@x\nACGT\n+\nIIII\n", encoding="utf-8")
            r2.write_text("@x\nACGT\n+\nIIII\n", encoding="utf-8")
            ref = d / "ref.fa"
            ref.write_text(">chr\nACGTACGT\n", encoding="utf-8")
            out_vcf = d / "out.vcf.gz"
            runner = FakeReadRunner(available=("minimap2", "samtools", "bcftools"))
            result = reads_to_vcf(r1, r2, ref, out_vcf, runner=runner)
            self.assertEqual(result, out_vcf)
            self.assertTrue(out_vcf.is_file())


if __name__ == "__main__":
    unittest.main()
