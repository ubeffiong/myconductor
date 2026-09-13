"""Minority-allele (heteroresistance) assessment.

What this module can and cannot do
----------------------------------
The previous implementation compared a VAF *supplied in the input* against a
fixed 3% floor and called the result heteroresistance detection. Two problems.

First, a fixed 3% floor is not universally valid. The limit below which an
alternate-allele fraction is indistinguishable from sequencing error depends on
the platform: Illumina, Nanopore and targeted amplicon sequencing have very
different error profiles, and a floor valid for one produces false positives on
another.

Second, and more fundamentally, real minority-variant detection needs read-level
evidence this module never sees — base and mapping qualities, strand bias, local
realignment around indels, platform error models, and contamination/mixed-
infection separation. Thresholding a number someone else computed is not
detection.

So this module reports an **assessment**, with ``assessable`` carried separately
from ``detected``. When the platform is unknown, alternate-read support is
missing, or depth is too low for the claimed fraction to be meaningful, the
finding says so rather than reporting "nothing found". A negative result and an
unperformed test are different things.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..core.models import (
    Call,
    DrugEvidence,
    HeteroresistanceFinding,
    Variant,
)

#: Approximate lower limits of reliable alternate-allele detection by platform.
#: These are conservative defaults for triage, NOT validated limits of
#: detection: a deployment must establish its own per-platform, per-locus LOD
#: against a dilution series and override these.
PLATFORM_LOD: dict[str, float] = {
    "illumina": 0.05,
    "targeted-amplicon": 0.03,
    "pacbio": 0.10,
    "iontorrent": 0.10,
    "nanopore": 0.20,
}

#: Below this many alternate reads, a fraction is not interpretable regardless
#: of what the VAF says.
MIN_ALT_READS = 5

#: Above this fraction the allele is the majority population, so "minority
#: subpopulation" does not apply.
MAJORITY_FLOOR = 0.90


@dataclass
class HeteroresistanceAssessor:
    depth_floor: int = 10
    min_alt_reads: int = MIN_ALT_READS
    majority_floor: float = MAJORITY_FLOOR
    default_platform: Optional[str] = None
    lod_table: Optional[dict[str, float]] = None

    def limit_of_detection(self, platform: Optional[str]) -> Optional[float]:
        if not platform:
            return None
        table = self.lod_table or PLATFORM_LOD
        return table.get(platform.strip().lower())

    def assess(
        self,
        variants: list[Variant],
        evidence: list[DrugEvidence],
    ) -> list[HeteroresistanceFinding]:
        """Assess every variant that carries resistance-relevant evidence.

        Runs across all lanes: a minority catalogued ``katG`` variant and a
        minority efflux-regulator change are both worth surfacing.
        """
        by_key = {v.key(): v for v in variants}
        findings: list[HeteroresistanceFinding] = []
        seen: set[tuple[str, str]] = set()

        for ev in evidence:
            # Anything that establishes resistance or withholds susceptibility
            # is worth checking for subclonality.
            if ev.call not in (Call.RESISTANT, Call.INDETERMINATE):
                continue
            if ev.variant_key is None:
                continue
            variant = by_key.get(ev.variant_key)
            if variant is None:
                continue
            dedupe_key = (ev.variant_key, ev.drug)
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            findings.append(self._assess_one(variant, ev))

        return findings

    def _assess_one(self, variant: Variant,
                    ev: DrugEvidence) -> HeteroresistanceFinding:
        platform = variant.platform or self.default_platform
        lod = self.limit_of_detection(platform)

        def finding(assessable: bool, detected: bool,
                    note: str) -> HeteroresistanceFinding:
            return HeteroresistanceFinding(
                variant_label=variant.label(),
                variant_key=variant.key(),
                drug=ev.drug,
                assessable=assessable,
                detected=detected,
                note=note,
                vaf=variant.vaf,
                depth=variant.depth,
                alt_depth=variant.alt_depth,
                limit_of_detection=lod,
                platform=platform,
            )

        # -- the conditions under which we decline to assess ---------------
        if variant.vaf is None:
            return finding(False, False,
                           "no allele fraction in the input; minority-allele "
                           "status cannot be assessed")
        if platform is None:
            return finding(False, False,
                           "sequencing platform unknown, so no limit of "
                           "detection applies; a fixed fraction floor is not "
                           "valid across platforms")
        if lod is None:
            return finding(False, False,
                           f"no limit of detection configured for platform "
                           f"{platform!r}")
        if variant.depth is None or variant.depth < self.depth_floor:
            return finding(False, False,
                           f"depth {variant.depth or 'unknown'} below floor "
                           f"{self.depth_floor}x; fraction not interpretable")

        alt_reads = variant.alt_depth
        if alt_reads is None and variant.depth is not None:
            alt_reads = int(round(variant.vaf * variant.depth))
        if alt_reads is not None and alt_reads < self.min_alt_reads:
            return finding(False, False,
                           f"only {alt_reads} alternate read(s), below the "
                           f"{self.min_alt_reads}-read minimum; fraction not "
                           f"interpretable")

        # -- assessable ----------------------------------------------------
        if variant.vaf >= self.majority_floor:
            return finding(True, False,
                           f"allele is the majority population at "
                           f"{variant.vaf:.0%}; not a minority subpopulation")
        if variant.vaf < lod:
            return finding(True, False,
                           f"allele fraction {variant.vaf:.1%} is below the "
                           f"{platform} limit of detection ({lod:.0%}); "
                           f"indistinguishable from sequencing error")

        return finding(
            True, True,
            f"Minority allele at {variant.vaf:.0%} "
            f"({alt_reads}/{variant.depth} reads), above the {platform} limit "
            f"of detection ({lod:.0%}). A consensus-only caller may not report "
            f"this. Read-level confirmation (strand bias, mapping quality, "
            f"contamination check) is required before acting on it."
        )


def assess(variants: list[Variant], evidence: list[DrugEvidence],
           depth_floor: int = 10,
           default_platform: Optional[str] = None) -> list[HeteroresistanceFinding]:
    """Convenience wrapper around :class:`HeteroresistanceAssessor`."""
    return HeteroresistanceAssessor(
        depth_floor=depth_floor, default_platform=default_platform
    ).assess(variants, evidence)
