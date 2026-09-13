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
    EngineRef,
    HeteroresistanceFinding,
    Lane,
    Tier,
    Variant,
)
from .calibration import CalibrationTable, PriorSource, posterior_resistance_probability

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
    #: Calibrated per-platform/drug/lineage detection curves and a real prior
    #: source. Both default to ``None``: without them, behaviour is identical
    #: to a deployment with no calibration configured at all, and no posterior
    #: is ever computed. See ``modules.calibration``.
    calibration: Optional[CalibrationTable] = None
    prior_source: Optional[PriorSource] = None

    def limit_of_detection(self, platform: Optional[str]) -> Optional[float]:
        if not platform:
            return None
        table = self.lod_table or PLATFORM_LOD
        return table.get(platform.strip().lower())

    def assess(
        self,
        variants: list[Variant],
        evidence: list[DrugEvidence],
        lineage: Optional[str] = None,
    ) -> list[HeteroresistanceFinding]:
        """Assess every variant that carries resistance-relevant evidence.

        Runs across all lanes: a minority catalogued ``katG`` variant and a
        minority efflux-regulator change are both worth surfacing. ``lineage``
        is the sample's lineage, when known, and is used only to look up a
        lineage-specific calibration entry; it never changes the
        assessable/detected logic below.
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
            findings.append(self._assess_one(variant, ev, lineage))

        return findings

    def _assess_one(self, variant: Variant, ev: DrugEvidence,
                    lineage: Optional[str] = None) -> HeteroresistanceFinding:
        platform = variant.platform or self.default_platform
        lod = self.limit_of_detection(platform)

        def finding(assessable: bool, detected: bool, note: str,
                    posterior=None, interval=None,
                    calibration_source: str = "uncalibrated default"
                    ) -> HeteroresistanceFinding:
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
                posterior_resistance_probability=posterior,
                posterior_interval=interval,
                calibration_source=calibration_source,
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

        posterior_estimate = None
        if self.calibration is not None and self.prior_source is not None:
            posterior_estimate = posterior_resistance_probability(
                variant, ev.drug, self.calibration, self.prior_source, lineage)

        note = (
            f"Minority allele at {variant.vaf:.0%} "
            f"({alt_reads}/{variant.depth} reads), above the {platform} limit "
            f"of detection ({lod:.0%}). A consensus-only caller may not report "
            f"this. Read-level confirmation (strand bias, mapping quality, "
            f"contamination check) is required before acting on it."
        )
        if posterior_estimate is not None:
            note += (
                f" Calibrated posterior probability of a true "
                f"resistance-conferring subpopulation: "
                f"{posterior_estimate.probability:.0%} "
                f"(range {posterior_estimate.credible_interval[0]:.0%}-"
                f"{posterior_estimate.credible_interval[1]:.0%})."
            )
        return finding(
            True, True, note,
            posterior=(posterior_estimate.probability if posterior_estimate else None),
            interval=(posterior_estimate.credible_interval if posterior_estimate else None),
            calibration_source=(posterior_estimate.basis if posterior_estimate
                                else "uncalibrated default"),
        )


def assess(variants: list[Variant], evidence: list[DrugEvidence],
           depth_floor: int = 10,
           default_platform: Optional[str] = None) -> list[HeteroresistanceFinding]:
    """Convenience wrapper around :class:`HeteroresistanceAssessor`."""
    return HeteroresistanceAssessor(
        depth_floor=depth_floor, default_platform=default_platform
    ).assess(variants, evidence)


_CALIBRATION_ENGINE = EngineRef(
    name="myconductor-heteroresistance-calibration", version="0.1.0")


def heteroresistance_evidence(
    findings: list[HeteroresistanceFinding],
    variants: list[Variant],
) -> list[DrugEvidence]:
    """Turn a calibrated posterior into evidence — never a resistance call.

    Only emitted for findings that actually carry a computed posterior (i.e.
    both a calibration curve and a prior were configured and resolved); with
    neither configured this returns an empty list and pipeline behaviour is
    unchanged. ``Tier.PREDICTED`` structurally forbids ``Call.RESISTANT``
    (``DrugEvidence.__post_init__``), so however high the posterior, the most
    this can do is make ``INDETERMINATE`` — withholding susceptibility —
    quantitatively better justified.
    """
    by_key = {v.key(): v for v in variants}
    out: list[DrugEvidence] = []
    for f in findings:
        if f.posterior_resistance_probability is None:
            continue
        variant = by_key.get(f.variant_key)
        interval = f.posterior_interval or (None, None)
        out.append(DrugEvidence(
            drug=f.drug,
            call=Call.INDETERMINATE,
            tier=Tier.PREDICTED,
            lane=Lane.VUS,
            confidence=f.posterior_resistance_probability,
            variant=variant.identity if variant else None,
            engine=_CALIBRATION_ENGINE,
            limitations=(
                "posterior probability from an explicit measurement model "
                "(calibrated detection curve + isolate-level prior), not a "
                "graded catalogue entry or a laboratory phenotype",
                "cannot, and structurally does not, establish resistance on "
                "its own",
                f"basis: {f.calibration_source}",
            ),
            rationale=(
                f"Minority allele {f.variant_label} at {f.vaf:.0%} VAF "
                f"({f.platform}): calibrated posterior probability of a true "
                f"resistance-conferring subpopulation is "
                f"{f.posterior_resistance_probability:.0%}"
                + (f" (range {interval[0]:.0%}-{interval[1]:.0%})"
                   if interval[0] is not None else "")
                + f". Susceptibility withheld for {f.drug}; resistance not "
                  f"asserted."
            ),
        ))
    return out
