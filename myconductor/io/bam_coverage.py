"""Derive a callable-locus mask directly from a BAM, via mosdepth/samtools.

Scope, precisely
-----------------
This closes the *coverage* half of "FASTQ/BAM ingestion" — nothing else. The
README already documents the exact manual command
(``mosdepth --by targets.bed --thresholds 10``) a user runs today to produce
the depth table ``io/callable_mask.py::CallableMask.from_tsv`` reads. This
module only automates running that command (or a ``samtools depth``
fallback) and parsing its output, so a BAM plus a BED of target loci can go
straight to a :class:`~myconductor.io.callable_mask.CallableMask` without a
manual step in between.

It does **not** perform alignment or variant calling. FASTQ/BAM to genotype
remains an external, per-deployment choice of aligner (bwa-mem2, minimap2)
and variant caller (GATK, bcftools, DeepVariant) that needs its own
validation; that caller's VCF output is the input ``io/adapter.py`` already
knows how to normalise. Conflating the two would mean Myconductor silently
picking (and thereby implicitly endorsing) an alignment/calling pipeline it
has not evaluated — exactly the kind of unvalidated default this codebase
refuses elsewhere.

The BED of target loci is supplied by the caller, not invented here — the
bundled organism profile deliberately ships no genomic coordinates (see
``catalogue/profile.py``), and a locus BED is exactly what mosdepth's
``--by`` flag already needs.

Unvalidated against real tool output
-------------------------------------
The mosdepth parser below is written from mosdepth's documented
``--thresholds``/``--by`` output schema, in the same spirit as
``adapters/tbprofiler.py`` and ``adapters/mykrobe.py`` — it has not been run
against a real mosdepth invocation in this environment. Add a golden-file
test against your own pinned mosdepth version before relying on it in
production; the ``samtools depth`` fallback is simpler and closer to
``samtools``'s stable plain-text contract.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Protocol

from ..core.models import LocusCoverage
from .callable_mask import CallableMask, CallableMaskError


class ToolRunner(Protocol):
    """A seam so this module is unit-testable without mosdepth/samtools
    installed. Deliberately not shared with ``mycobench.stages.ToolRunner``:
    myconductor has zero runtime dependencies and must not depend on
    mycobench.
    """

    def run(self, command: list[str], timeout: Optional[int] = None) -> int:
        """Run a command that writes its own output files (mosdepth)."""
        ...

    def run_capturing_stdout(self, command: list[str], out_path: Path,
                             timeout: Optional[int] = None) -> int:
        """Run a command, writing its stdout to ``out_path`` (samtools)."""
        ...

    def available(self, executable: str) -> bool:
        ...


class SubprocessToolRunner:
    def run(self, command: list[str], timeout: Optional[int] = None) -> int:
        try:
            completed = subprocess.run(
                command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                timeout=timeout, check=False)
            return completed.returncode
        except FileNotFoundError:
            return 127
        except subprocess.TimeoutExpired:
            return 124

    def run_capturing_stdout(self, command: list[str], out_path: Path,
                             timeout: Optional[int] = None) -> int:
        try:
            with open(out_path, "wb") as handle:
                completed = subprocess.run(
                    command, stdout=handle, stderr=subprocess.PIPE,
                    timeout=timeout, check=False)
                return completed.returncode
        except FileNotFoundError:
            return 127
        except subprocess.TimeoutExpired:
            return 124

    def available(self, executable: str) -> bool:
        return shutil.which(executable) is not None


def _read_bed_loci(bed_path: str | Path) -> dict[str, tuple[str, int, int]]:
    """chrom/start/end per locus name, from column 4 of a BED."""
    loci: dict[str, tuple[str, int, int]] = {}
    for raw in Path(bed_path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 4:
            raise CallableMaskError(
                f"{bed_path}: BED needs at least 4 columns "
                f"(chrom, start, end, locus); got {len(f)}")
        loci[f[3].strip()] = (f[0], int(f[1]), int(f[2]))
    if not loci:
        raise CallableMaskError(f"{bed_path}: no locus regions found")
    return loci


def _open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _mosdepth_mask(bam_path: Path, bed_path: Path, depth_floor: int,
                   runner: ToolRunner) -> CallableMask:
    """Parse mosdepth's ``--by``/``--thresholds`` output.

    Expected schema (mosdepth >= 0.3): ``<prefix>.regions.bed.gz`` with columns
    ``chrom start end name mean_depth``, and ``<prefix>.thresholds.bed.gz``
    with a header ``chrom start end region <T>X ...`` and one count column per
    requested threshold (bases at/above that depth).
    """
    with tempfile.TemporaryDirectory(prefix="myconductor-mosdepth-") as tmp:
        prefix = Path(tmp) / "sample"
        code = runner.run([
            "mosdepth", "--by", str(bed_path), "--thresholds", str(depth_floor),
            "--no-per-base", str(prefix), str(bam_path),
        ])
        if code != 0:
            raise CallableMaskError(
                f"mosdepth exited {code} on {bam_path}; see its stderr for "
                f"detail, or run it manually: mosdepth --by {bed_path} "
                f"--thresholds {depth_floor} <prefix> {bam_path}")

        thresholds_path = Path(f"{prefix}.thresholds.bed.gz")
        regions_path = Path(f"{prefix}.regions.bed.gz")
        if not thresholds_path.is_file() or not regions_path.is_file():
            raise CallableMaskError(
                f"mosdepth did not produce the expected output files under "
                f"{prefix}; got exit code 0 but no thresholds/regions BED")

        mean_depth: dict[str, float] = {}
        with _open_maybe_gzip(regions_path) as handle:
            for line in handle:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5:
                    mean_depth[f[3]] = float(f[4])

        loci: dict[str, LocusCoverage] = {}
        with _open_maybe_gzip(thresholds_path) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            wanted_col = None
            for i, name in enumerate(header):
                if name.upper() == f"{depth_floor}X":
                    wanted_col = i
                    break
            if wanted_col is None:
                raise CallableMaskError(
                    f"mosdepth's thresholds output has no {depth_floor}X "
                    f"column (header: {header}); it may not have honoured "
                    f"--thresholds {depth_floor}")
            for line in handle:
                f = line.rstrip("\n").split("\t")
                if len(f) <= wanted_col:
                    continue
                locus = f[3]
                start, end = int(f[1]), int(f[2])
                total = max(1, end - start)
                bases_at_floor = int(f[wanted_col])
                loci[locus] = LocusCoverage(
                    locus=locus,
                    mean_depth=(round(mean_depth[locus])
                               if locus in mean_depth else None),
                    callable_fraction=min(1.0, bases_at_floor / total),
                    source="bam:mosdepth",
                )
        if not loci:
            raise CallableMaskError(
                f"mosdepth produced no per-locus thresholds for {bam_path}")
        return CallableMask(loci, source="bam:mosdepth")


def _samtools_depth_mask(bam_path: Path, bed_path: Path, depth_floor: int,
                         runner: ToolRunner) -> CallableMask:
    regions = _read_bed_loci(bed_path)
    by_chrom: dict[str, list[tuple[str, int, int]]] = {}
    for locus, (chrom, start, end) in regions.items():
        by_chrom.setdefault(chrom, []).append((locus, start, end))

    with tempfile.TemporaryDirectory(prefix="myconductor-samtools-") as tmp:
        out_path = Path(tmp) / "depth.tsv"
        code = runner.run_capturing_stdout(
            ["samtools", "depth", "-a", "-b", str(bed_path), str(bam_path)],
            out_path)
        if code != 0:
            raise CallableMaskError(
                f"samtools depth exited {code} on {bam_path}; run it "
                f"manually to see why: samtools depth -a -b {bed_path} "
                f"{bam_path}")

        callable_bp = {locus: 0 for locus in regions}
        depth_sum = {locus: 0 for locus in regions}
        depth_n = {locus: 0 for locus in regions}

        if out_path.is_file() and out_path.stat().st_size:
            for raw in out_path.read_text(encoding="utf-8").splitlines():
                f = raw.split("\t")
                if len(f) < 3:
                    continue
                chrom, pos, depth = f[0], int(f[1]) - 1, int(f[2])
                for locus, start, end in by_chrom.get(chrom, []):
                    if start <= pos < end:
                        depth_sum[locus] += depth
                        depth_n[locus] += 1
                        if depth >= depth_floor:
                            callable_bp[locus] += 1

        loci: dict[str, LocusCoverage] = {}
        for locus, (chrom, start, end) in regions.items():
            total = max(1, end - start)
            loci[locus] = LocusCoverage(
                locus=locus,
                mean_depth=(round(depth_sum[locus] / depth_n[locus])
                           if depth_n[locus] else None),
                callable_fraction=min(1.0, callable_bp[locus] / total),
                source="bam:samtools",
            )
        return CallableMask(loci, source="bam:samtools")


def mask_from_bam(
    bam_path: str | Path,
    bed_path: str | Path,
    depth_floor: int = 10,
    tool: str = "auto",
    runner: Optional[ToolRunner] = None,
) -> CallableMask:
    """Build a :class:`CallableMask` directly from a BAM and a locus BED.

    ``tool``: ``"auto"`` prefers ``mosdepth`` (matches the README's documented
    manual workflow) and falls back to ``samtools depth``; pass ``"mosdepth"``
    or ``"samtools"`` to require one specifically. Raises
    :class:`CallableMaskError` — never silently returns an absent mask — when
    neither tool is available or the chosen tool fails, so a missing
    dependency is never mistaken for "no coverage evidence was supplied."
    """
    bam_path, bed_path = Path(bam_path), Path(bed_path)
    runner = runner or SubprocessToolRunner()

    chosen = tool
    if tool == "auto":
        if runner.available("mosdepth"):
            chosen = "mosdepth"
        elif runner.available("samtools"):
            chosen = "samtools"
        else:
            raise CallableMaskError(
                "neither mosdepth nor samtools is on PATH. Install one of "
                "them, or run the equivalent by hand and pass the result via "
                "--mask: mosdepth --by targets.bed --thresholds "
                f"{depth_floor} <prefix> {bam_path}")

    if chosen == "mosdepth":
        if not runner.available("mosdepth"):
            raise CallableMaskError("mosdepth is not on PATH")
        return _mosdepth_mask(bam_path, bed_path, depth_floor, runner)
    if chosen == "samtools":
        if not runner.available("samtools"):
            raise CallableMaskError("samtools is not on PATH")
        return _samtools_depth_mask(bam_path, bed_path, depth_floor, runner)
    raise CallableMaskError(f"unknown coverage tool {chosen!r}; choose from "
                            f"auto, mosdepth, samtools")
