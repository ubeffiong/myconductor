"""Say what a run is about to download, before it downloads it.

A cohort run fetches tens of gigabytes over several hours. A user is entitled to
see the number first rather than discover it from a disk that filled up.

The estimate is labelled as one. Compressed FASTQ size is derived from the base
count NCBI reports, using a bytes-per-base factor; where no base count is
available the script says so rather than inventing a total, and the samples with
unknown size are named.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .cohort import lock_bases, read_rows

#: Compressed-FASTQ bytes per sequenced base. Illumina short-read FASTQ
#: compresses to roughly 0.5-0.7 bytes per base depending on read length and
#: quality-score entropy; the midpoint is used for the headline figure and the
#: bounds are reported alongside it.
BYTES_PER_BASE = 0.615
BYTES_PER_BASE_LOW = 0.50
BYTES_PER_BASE_HIGH = 0.68

#: TB-Profiler writes a BAM, a VCF and a results bundle per isolate. Measured
#: against nothing here, so it is declared as a rule-of-thumb multiplier on the
#: download rather than a per-sample constant.
WORKING_SPACE_MULTIPLIER = 3


def human(num_bytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(num_bytes) < 1024 or unit == "TB":
            return (f"{num_bytes:.0f} B" if unit == "B"
                    else f"{num_bytes:.1f} {unit}")
        num_bytes /= 1024
    return f"{num_bytes:.1f} TB"


@dataclass
class Estimate:
    n_rows: int = 0
    pending: int = 0
    already_present: int = 0
    known_bases: int = 0
    unknown: list[str] = field(default_factory=list)

    @property
    def mid(self) -> float:
        return self.known_bases * BYTES_PER_BASE

    @property
    def low(self) -> float:
        return self.known_bases * BYTES_PER_BASE_LOW

    @property
    def high(self) -> float:
        return self.known_bases * BYTES_PER_BASE_HIGH

    @property
    def nothing_to_fetch(self) -> bool:
        return self.pending == 0

    def render(self, sheet: Path) -> str:
        if self.n_rows == 0:
            return "No samples in the sheet; nothing to download."
        if self.nothing_to_fetch:
            return (f"All {self.n_rows} sample(s) are already downloaded; "
                    f"nothing to fetch.")

        lines = ["About to download from NCBI:"]
        if self.known_bases:
            lines.append(
                f"  read sets        : {self.pending:>4}  ~{human(self.mid)}"
                f"  ({human(self.low)} - {human(self.high)})")
            lines.append(
                f"  ESTIMATED TOTAL  :       ~{human(self.mid)}")
            lines.append(
                f"  working space    :       ~"
                f"{human(self.mid * WORKING_SPACE_MULTIPLIER)} including "
                f"alignments and tool output")
        else:
            lines.append(f"  read sets        : {self.pending:>4}  size unknown")

        if self.unknown:
            shown = ", ".join(self.unknown[:4])
            lines.append(
                f"  NOTE: no base count for {len(self.unknown)} sample(s) "
                f"({shown}{' ...' if len(self.unknown) > 4 else ''});")
            lines.append(
                f"        their reads are NOT included in the total above.")
            lines.append(
                f"        Run 'mycobench validate-cohort --samples {sheet} "
                f"--online --write-lock {sheet.with_suffix('.lock.json')}' "
                f"for an exact figure.")
        if self.already_present:
            lines.append(f"  ({self.already_present} sample(s) already "
                         f"downloaded and will be skipped)")
        lines.append("")
        lines.append("  Sizes are estimated from NCBI base counts; the real")
        lines.append("  total will differ somewhat.")
        return "\n".join(lines)


def reads_present(data_dir: Path, sample_id: str, run: str) -> bool:
    directory = data_dir / sample_id
    if (directory / f"{run}_1.fastq.gz").is_file() and \
            (directory / f"{run}_2.fastq.gz").is_file():
        return True
    # A single-file dump is accepted as present too, so a completed download in
    # either layout is not re-fetched.
    return (directory / f"{run}.fastq.gz").is_file()


def estimate(samples: str | Path, data_dir: str | Path) -> Estimate:
    sheet = Path(samples)
    data_dir = Path(data_dir)
    rows, _ = read_rows(sheet)
    bases = lock_bases(sheet)

    result = Estimate(n_rows=len(rows))
    for row in rows:
        sample_id = (row.get("sample_id") or "").strip()
        run = (row.get("sra_run") or "").strip()
        if not sample_id or not run:
            continue
        if reads_present(data_dir, sample_id, run):
            result.already_present += 1
            continue
        result.pending += 1
        # Prefer the lock's figure: it is what NCBI reported at verification
        # time, whereas the sheet column can have been edited by hand.
        count = bases.get(sample_id)
        if not count:
            try:
                count = int(row.get("bases") or 0)
            except ValueError:
                count = 0
        if count:
            result.known_bases += count
        else:
            result.unknown.append(sample_id)
    return result


def confirm(estimate_text: str, assume_yes: bool = False,
            interactive: bool = True) -> bool:
    """Show the estimate and get consent. Non-interactive callers proceed."""
    print(estimate_text)
    if "nothing to fetch" in estimate_text or "nothing to download" in estimate_text:
        return True
    if assume_yes:
        return True
    if not interactive:
        print("Not an interactive terminal; proceeding with the download.")
        return True
    try:
        reply = input("Download these files now? [y/N] ").strip().lower()
    except EOFError:
        return True
    return reply in ("y", "yes")
