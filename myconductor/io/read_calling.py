"""FASTQ -> BAM -> VCF, via bwa-mem2/minimap2 + bcftools/GATK.

Scope and honesty
------------------
This is the remaining half of "FASTQ/BAM ingestion" — the README previously
named it explicitly as needing "a validated aligner+caller choice" that this
project would not pick on a user's behalf. It still doesn't: this module is a
thin, swappable orchestration layer over established external tools, not a
new aligner or caller, and not a validated pipeline. Every default below
(minimap2's short-read preset, bcftools over GATK) is a common, defensible
choice for bacterial resequencing — not this project's own measured
recommendation. The output is a standard VCF that flows into
``io/adapter.py`` exactly like a VCF from anywhere else; Myconductor does not
treat it as more or less reliable than a VCF from a different pipeline
choice, and nothing downstream is aware this module produced it.

Written from each tool's documented CLI contract, in the same spirit as
``adapters/tbprofiler.py`` and ``io/bam_coverage.py`` — **not run against real
sequencing reads in this environment.** Before trusting its output for
anything beyond exercising the interface, validate it against a truth set
(paired FASTQ + an independently-called VCF for the same isolate, or a
spiked-in synthetic mutation) for your own aligner/caller/version
combination.

What this deliberately does not do
-----------------------------------
No base-quality recalibration, no duplicate marking, no joint genotyping
across samples, no filtering beyond each caller's own default. Those are
real pipeline decisions with their own accuracy trade-offs; adding silent
defaults for them here would be exactly the kind of unvalidated,
unattributable choice this project refuses to make on a deployment's behalf.
Pass ``extra_align_args``/``extra_call_args`` to add them explicitly, and
record what you passed in your own provenance — this module does not.
"""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional, Sequence

from .bam_coverage import SubprocessToolRunner, ToolRunner
from .callable_mask import CallableMaskError


class ReadCallingError(CallableMaskError):
    """Alignment or variant calling could not be completed."""


def _ensure_faidx(reference_fasta: Path, runner: ToolRunner) -> None:
    if Path(f"{reference_fasta}.fai").is_file():
        return
    if not runner.available("samtools"):
        raise ReadCallingError(
            "samtools is not on PATH and is needed to index the reference "
            f"({reference_fasta}.fai is missing)")
    code = runner.run(["samtools", "faidx", str(reference_fasta)])
    if code != 0:
        raise ReadCallingError(
            f"samtools faidx exited {code} on {reference_fasta}")


def _ensure_bwa_index(reference_fasta: Path, runner: ToolRunner) -> None:
    marker = Path(f"{reference_fasta}.bwt.2bit.64")
    if marker.is_file():
        return
    code = runner.run(["bwa-mem2", "index", str(reference_fasta)])
    if code != 0:
        raise ReadCallingError(
            f"bwa-mem2 index exited {code} on {reference_fasta}")


def align_reads(
    fastq1: str | Path,
    fastq2: str | Path,
    reference_fasta: str | Path,
    out_bam: str | Path,
    aligner: str = "auto",
    runner: Optional[ToolRunner] = None,
    threads: int = 4,
    read_group: Optional[str] = None,
    extra_align_args: Sequence[str] = (),
) -> Path:
    """Align paired-end reads and produce a coordinate-sorted, indexed BAM.

    ``aligner``: ``"auto"`` prefers minimap2 (its ``-ax sr`` preset needs no
    separate index step, so there is less to get wrong without being able to
    test against real reads) and falls back to bwa-mem2; pass either name to
    require it specifically. Raises :class:`ReadCallingError` — never
    silently produces an empty or partial BAM — on any tool failure.
    """
    fastq1, fastq2 = Path(fastq1), Path(fastq2)
    reference_fasta, out_bam = Path(reference_fasta), Path(out_bam)
    runner = runner or SubprocessToolRunner()
    read_group = read_group or f"@RG\\tID:sample\\tSM:sample\\tPL:ILLUMINA"

    chosen = aligner
    if aligner == "auto":
        if runner.available("minimap2"):
            chosen = "minimap2"
        elif runner.available("bwa-mem2"):
            chosen = "bwa-mem2"
        else:
            raise ReadCallingError(
                "neither minimap2 nor bwa-mem2 is on PATH. Install one of "
                "them, or align by hand and pass the resulting BAM to "
                "io.bam_coverage.mask_from_bam / your variant caller of "
                "choice directly.")
    if chosen not in ("minimap2", "bwa-mem2"):
        raise ReadCallingError(f"unknown aligner {chosen!r}; choose from "
                               f"auto, minimap2, bwa-mem2")
    if not runner.available(chosen):
        raise ReadCallingError(f"{chosen} is not on PATH")

    with tempfile.TemporaryDirectory(prefix="myconductor-align-") as tmp:
        sam_path = Path(tmp) / "aligned.sam"
        if chosen == "minimap2":
            command = ["minimap2", "-ax", "sr", "-t", str(threads),
                      "-R", read_group, *extra_align_args,
                      str(reference_fasta), str(fastq1), str(fastq2)]
        else:
            _ensure_bwa_index(reference_fasta, runner)
            command = ["bwa-mem2", "mem", "-t", str(threads),
                      "-R", read_group, *extra_align_args,
                      str(reference_fasta), str(fastq1), str(fastq2)]
        code = runner.run_capturing_stdout(command, sam_path)
        if code != 0 or not sam_path.is_file() or not sam_path.stat().st_size:
            raise ReadCallingError(
                f"{chosen} exited {code} aligning {fastq1}/{fastq2} against "
                f"{reference_fasta}; run it by hand to see stderr")

        if not runner.available("samtools"):
            raise ReadCallingError(
                "samtools is not on PATH (needed to sort and index the "
                "alignment)")
        out_bam.parent.mkdir(parents=True, exist_ok=True)
        code = runner.run(["samtools", "sort", "-@", str(threads),
                          "-o", str(out_bam), str(sam_path)])
        if code != 0:
            raise ReadCallingError(f"samtools sort exited {code}")
        code = runner.run(["samtools", "index", str(out_bam)])
        if code != 0:
            raise ReadCallingError(f"samtools index exited {code}")

    if not out_bam.is_file() or not out_bam.stat().st_size:
        raise ReadCallingError(f"no BAM produced at {out_bam}")
    return out_bam


def call_variants(
    bam_path: str | Path,
    reference_fasta: str | Path,
    out_vcf: str | Path,
    caller: str = "auto",
    runner: Optional[ToolRunner] = None,
    region_bed: Optional[str | Path] = None,
    extra_call_args: Sequence[str] = (),
) -> Path:
    """Call variants from an aligned, indexed BAM.

    ``caller``: ``"auto"`` prefers ``bcftools`` (a single, dependency-light
    binary; no Java runtime, no reference dictionary/index beyond ``.fai``)
    and falls back to GATK's ``HaplotypeCaller``. ``region_bed`` restricts
    calling to the given regions (the same BED a coverage mask would use);
    without it, the whole reference is called, which is slower and not
    necessary when only specific resistance loci matter.
    """
    bam_path = Path(bam_path)
    reference_fasta, out_vcf = Path(reference_fasta), Path(out_vcf)
    runner = runner or SubprocessToolRunner()
    _ensure_faidx(reference_fasta, runner)

    chosen = caller
    if caller == "auto":
        if runner.available("bcftools"):
            chosen = "bcftools"
        elif runner.available("gatk"):
            chosen = "gatk"
        else:
            raise ReadCallingError(
                "neither bcftools nor gatk is on PATH. Install one of "
                "them, or call variants by hand and feed the resulting VCF "
                "to myconductor analyze directly.")
    if chosen not in ("bcftools", "gatk"):
        raise ReadCallingError(f"unknown caller {chosen!r}; choose from "
                               f"auto, bcftools, gatk")
    if not runner.available(chosen):
        raise ReadCallingError(f"{chosen} is not on PATH")

    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="myconductor-call-") as tmp:
        if chosen == "bcftools":
            pileup_bcf = Path(tmp) / "pileup.bcf"
            mpileup_cmd = ["bcftools", "mpileup", "-f", str(reference_fasta)]
            if region_bed:
                mpileup_cmd += ["-R", str(region_bed)]
            mpileup_cmd += ["-Ob", "-o", str(pileup_bcf), str(bam_path)]
            code = runner.run(mpileup_cmd)
            if code != 0:
                raise ReadCallingError(f"bcftools mpileup exited {code}")

            call_cmd = ["bcftools", "call", "-mv", "-Oz",
                       "-o", str(out_vcf), *extra_call_args, str(pileup_bcf)]
            code = runner.run(call_cmd)
            if code != 0:
                raise ReadCallingError(f"bcftools call exited {code}")
        else:
            command = ["gatk", "HaplotypeCaller", "-R", str(reference_fasta),
                      "-I", str(bam_path), "-O", str(out_vcf)]
            if region_bed:
                command += ["-L", str(region_bed)]
            command += list(extra_call_args)
            code = runner.run(command)
            if code != 0:
                raise ReadCallingError(f"gatk HaplotypeCaller exited {code}")

    if not out_vcf.is_file() or not out_vcf.stat().st_size:
        raise ReadCallingError(f"no VCF produced at {out_vcf}")
    return out_vcf


def reads_to_vcf(
    fastq1: str | Path,
    fastq2: str | Path,
    reference_fasta: str | Path,
    out_vcf: str | Path,
    aligner: str = "auto",
    caller: str = "auto",
    runner: Optional[ToolRunner] = None,
    threads: int = 4,
    region_bed: Optional[str | Path] = None,
    read_group: Optional[str] = None,
) -> Path:
    """Convenience: FASTQ pair straight to a VCF.

    Equivalent to :func:`align_reads` then :func:`call_variants`, keeping the
    intermediate BAM in a temporary directory. Call the two functions
    separately (as ``mycobench``'s stages do) when the BAM itself is also
    wanted — e.g. to derive a coverage mask with
    ``io.bam_coverage.mask_from_bam`` over the same alignment.
    """
    runner = runner or SubprocessToolRunner()
    with tempfile.TemporaryDirectory(prefix="myconductor-reads2vcf-") as tmp:
        bam_path = Path(tmp) / "aligned.bam"
        align_reads(fastq1, fastq2, reference_fasta, bam_path,
                   aligner=aligner, runner=runner, threads=threads,
                   read_group=read_group)
        return call_variants(bam_path, reference_fasta, out_vcf,
                            caller=caller, runner=runner,
                            region_bed=region_bed)
