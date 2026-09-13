"""Callable-locus evidence — the module that breaks a circularity.

A variant-only VCF proves that a variant was *seen*. It proves nothing about
the positions where no variant was called: those may be wild type, or they may
never have been sequenced at all. Inferring coverage from the variant list is
circular, and it is the mechanism by which a pipeline reports a drug
susceptible whose target locus was never read.

So susceptibility requires evidence from *outside* the variant list. Three
sources are accepted:

    depth table   TSV of per-locus mean depth and callable fraction — the
                  easiest to produce (``mosdepth --by``, ``samtools depth``
                  thresholded, GATK DepthOfCoverage)
    BED mask      callable regions with the locus named in column 4, e.g.
                  GATK CallableLoci or ``mosdepth --quantize`` output
    gVCF          non-variant reference blocks carrying END and MIN_DP

Unknown is not callable
-----------------------
Where the callable fraction cannot be established, ``LocusCoverage`` reports
``callable_fraction=None`` and ``is_callable`` returns False. That is the whole
point: the failure mode this module exists to prevent is treating an unknown as
a pass. A missing mask therefore yields NOT_ASSESSED for every drug, which is
the correct answer to "is this drug usable?" when nobody checked.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

from ..core.models import LocusCoverage


class CallableMaskError(ValueError):
    """The mask input could not be parsed or is internally inconsistent."""


class CallableMask:
    """Per-locus callable evidence, independent of the variant list."""

    def __init__(self, loci: Optional[dict[str, LocusCoverage]] = None,
                 source: str = "absent"):
        self._loci: dict[str, LocusCoverage] = dict(loci or {})
        self.source = source

    # -- construction -----------------------------------------------------
    @classmethod
    def absent(cls) -> "CallableMask":
        """The honest default when no coverage evidence was supplied."""
        return cls({}, source="absent")

    @classmethod
    def merge(cls, *masks: "CallableMask") -> "CallableMask":
        """Combine masks, preferring the better-evidenced coverage per locus.

        Used when an engine adapter supplies coverage (e.g. TB-Profiler's
        per-gene QC) alongside a locally generated mask. A known fraction always
        beats an unknown one; between two known fractions the higher wins, and
        the source is recorded so a report can say who supplied it.
        """
        combined: dict[str, LocusCoverage] = {}
        sources: list[str] = []
        for mask in masks:
            if not mask.is_present:
                continue
            sources.append(mask.source)
            for locus, cov in mask._loci.items():
                existing = combined.get(locus)
                if existing is None:
                    combined[locus] = cov
                    continue
                if existing.callable_fraction is None:
                    combined[locus] = cov
                elif (cov.callable_fraction is not None
                      and cov.callable_fraction > existing.callable_fraction):
                    combined[locus] = cov
        if not combined:
            return cls.absent()
        return cls(combined, source=" + ".join(dict.fromkeys(sources)))

    @classmethod
    def from_tsv(cls, path: str | Path) -> "CallableMask":
        """Load a per-locus depth table.

        Required column: ``locus``. Then either ``callable_fraction``, or both
        ``covered_bp`` and ``total_bp`` so the fraction can be computed.
        Optional: ``mean_depth``.
        """
        text = Path(path).read_text(encoding="utf-8")
        rows = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
        if not rows:
            raise CallableMaskError(f"{path}: no data rows")

        header = [h.strip().lower() for h in rows[0].split("\t")]
        if "locus" not in header:
            raise CallableMaskError(f"{path}: missing required 'locus' column")

        loci: dict[str, LocusCoverage] = {}
        for ln in rows[1:]:
            cols = dict(zip(header, [c.strip() for c in ln.split("\t")]))
            locus = cols.get("locus", "")
            if not locus:
                continue
            fraction = _fraction_from(cols, path)
            loci[locus] = LocusCoverage(
                locus=locus,
                mean_depth=_opt_int(cols.get("mean_depth")),
                callable_fraction=fraction,
                source="depth-table",
            )
        if not loci:
            raise CallableMaskError(f"{path}: parsed no loci")
        return cls(loci, source="depth-table")

    @classmethod
    def from_bed(cls, path: str | Path,
                 locus_lengths: Optional[dict[str, int]] = None) -> "CallableMask":
        """Load a callable-regions BED whose column 4 names the locus.

        ``locus_lengths`` supplies each locus's total span so a fraction can be
        computed. Without it, covered base pairs are recorded but the fraction
        stays unknown — and unknown is not callable. Lengths come from the
        reference annotation, which Myconductor does not ship; the organism
        profile carries them once a catalogue has been ingested.
        """
        locus_lengths = locus_lengths or {}
        covered: dict[str, int] = {}
        depths: dict[str, list[int]] = {}

        for raw in Path(path).read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 4:
                raise CallableMaskError(
                    f"{path}: BED needs at least 4 columns "
                    f"(chrom, start, end, locus); got {len(f)}"
                )
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError as exc:
                raise CallableMaskError(f"{path}: bad BED coordinates in {line!r}") from exc
            if end < start:
                raise CallableMaskError(f"{path}: BED end < start in {line!r}")
            locus = f[3].strip()
            covered[locus] = covered.get(locus, 0) + (end - start)
            if len(f) >= 5:
                d = _opt_int(f[4])
                if d is not None:
                    depths.setdefault(locus, []).append(d)

        if not covered:
            raise CallableMaskError(f"{path}: parsed no callable intervals")

        loci: dict[str, LocusCoverage] = {}
        for locus, bp in covered.items():
            total = locus_lengths.get(locus)
            fraction = (min(1.0, bp / total) if total else None)
            mean_depth = (round(sum(depths[locus]) / len(depths[locus]))
                          if depths.get(locus) else None)
            loci[locus] = LocusCoverage(
                locus=locus,
                mean_depth=mean_depth,
                callable_fraction=fraction,
                source="bed" if total else "bed (fraction unknown: no locus length)",
            )
        return cls(loci, source="bed")

    @classmethod
    def from_gvcf(cls, path: str | Path,
                  locus_spans: dict[str, tuple[int, int]],
                  depth_floor: int = 10) -> "CallableMask":
        """Derive callable fractions from gVCF non-variant reference blocks.

        ``locus_spans`` maps locus -> (start, end) in assembly coordinates and
        is required: a gVCF says which positions are covered, but only the
        reference annotation says which positions belong to which gene.
        """
        if not locus_spans:
            raise CallableMaskError(
                "from_gvcf requires locus_spans; a gVCF alone cannot say which "
                "coordinates belong to which locus. Supply spans from the "
                "reference annotation, or use from_tsv."
            )

        blocks: list[tuple[int, int, Optional[int]]] = []
        for raw in Path(path).read_text(encoding="utf-8").splitlines():
            if raw.startswith("#") or not raw.strip():
                continue
            f = raw.split("\t")
            if len(f) < 8:
                continue
            try:
                pos = int(f[1])
            except ValueError:
                continue
            info = _info_dict(f[7])
            end = _opt_int(info.get("END")) or pos
            min_dp = _opt_int(info.get("MIN_DP"))
            if min_dp is None and len(f) >= 10:
                fmt = dict(zip(f[8].split(":"), f[9].split(":")))
                min_dp = _opt_int(fmt.get("MIN_DP") or fmt.get("DP"))
            blocks.append((pos, end, min_dp))

        if not blocks:
            raise CallableMaskError(f"{path}: no reference blocks found")

        loci: dict[str, LocusCoverage] = {}
        for locus, (lo, hi) in locus_spans.items():
            total = max(0, hi - lo + 1)
            if total == 0:
                continue
            callable_bp, depth_sum, depth_bp = 0, 0, 0
            for start, end, min_dp in blocks:
                overlap = min(end, hi) - max(start, lo) + 1
                if overlap <= 0:
                    continue
                if min_dp is not None and min_dp >= depth_floor:
                    callable_bp += overlap
                if min_dp is not None:
                    depth_sum += min_dp * overlap
                    depth_bp += overlap
            loci[locus] = LocusCoverage(
                locus=locus,
                mean_depth=(round(depth_sum / depth_bp) if depth_bp else None),
                callable_fraction=callable_bp / total,
                source="gvcf",
            )
        return cls(loci, source="gvcf")

    # -- use --------------------------------------------------------------
    @property
    def is_present(self) -> bool:
        return bool(self._loci)

    def coverage_for(self, locus: str) -> LocusCoverage:
        return self._loci.get(locus, LocusCoverage(locus=locus, source="absent"))

    def assess(self, loci: Iterable[str], *, depth_floor: int,
               fraction_floor: float) -> tuple[bool, list[LocusCoverage], list[str]]:
        """Are all of these loci callable?

        Returns ``(all_callable, per_locus_coverage, reasons_not_callable)``.
        An empty locus list is *not* a pass: a drug with no declared loci
        cannot be assessed.
        """
        loci = list(loci)
        if not loci:
            return False, [], ["no loci declared for this drug in the profile"]
        if not self.is_present:
            return False, [], ["no callable-locus evidence supplied "
                               "(pass a depth table, BED mask or gVCF)"]

        coverages, reasons = [], []
        for locus in loci:
            cov = self.coverage_for(locus)
            coverages.append(cov)
            if cov.source == "absent":
                reasons.append(f"{locus}: absent from the callable mask")
            elif cov.callable_fraction is None:
                reasons.append(f"{locus}: callable fraction unknown")
            elif cov.callable_fraction < fraction_floor:
                reasons.append(
                    f"{locus}: only {cov.callable_fraction:.0%} callable "
                    f"(floor {fraction_floor:.0%})"
                )
            elif cov.mean_depth is not None and cov.mean_depth < depth_floor:
                reasons.append(
                    f"{locus}: mean depth {cov.mean_depth}x below floor {depth_floor}x"
                )
        return (not reasons), coverages, reasons

    def __repr__(self) -> str:  # pragma: no cover - diagnostics only
        return f"CallableMask(source={self.source!r}, loci={len(self._loci)})"


# -- helpers --------------------------------------------------------------
def _opt_int(value: Optional[str]) -> Optional[int]:
    if value is None or value == "" or value == ".":
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _opt_float(value: Optional[str]) -> Optional[float]:
    if value is None or value == "" or value == ".":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _fraction_from(cols: dict[str, str], path: str | Path) -> Optional[float]:
    explicit = _opt_float(cols.get("callable_fraction"))
    if explicit is not None:
        if not 0.0 <= explicit <= 1.0:
            raise CallableMaskError(
                f"{path}: callable_fraction {explicit} outside 0..1 for "
                f"locus {cols.get('locus')!r}"
            )
        return explicit
    covered, total = _opt_int(cols.get("covered_bp")), _opt_int(cols.get("total_bp"))
    if covered is not None and total:
        return min(1.0, covered / total)
    return None


def _info_dict(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            out[k.strip().upper()] = v.strip()
        elif token.strip():
            out[token.strip().upper()] = "true"
    return out
