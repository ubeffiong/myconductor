"""Upstream quality assessment — including what was *not* assessed.

Myconductor's analysis begins after variant calling, but its reports look like
complete genomic analyses. That gap is dangerous: a reader cannot tell from a
drug-susceptibility table whether the species was confirmed, whether the sample
was contaminated, or whether two patients' samples were swapped.

This module closes the gap in the only honest way available to a tool that
receives a variant file: by enumerating every control that matters and
reporting ``not_performed`` for each one it cannot run. A report that says
"species confirmation: not performed — requires reads or an assembly" is far
safer than a report that silently omits the question.

Checks marked ``not_performed`` here become ``pass``/``fail`` once Phase 3
lands FASTQ/BAM ingestion or integration with an upstream pipeline.
"""
from __future__ import annotations

from typing import Optional

from ..catalogue.profile import OrganismProfile
from ..core.models import QCFinding
from .adapter import AdapterResult
from .callable_mask import CallableMask

#: Controls that require read-level or assembly-level data. Each entry is
#: (check name, what it needs, why it matters).
_READ_LEVEL_CONTROLS = [
    ("species_confirmation", "reads or an assembly",
     "an NTM misidentified as MTBC would be interpreted against the wrong "
     "catalogue entirely"),
    ("mtbc_vs_ntm", "reads or an assembly",
     "NTM species have different intrinsic resistance and different targets"),
    ("read_quality", "FASTQ",
     "low-quality bases inflate false-positive minority alleles"),
    ("contamination", "reads",
     "cross-contamination produces spurious minority resistance alleles"),
    ("coverage_breadth", "BAM/CRAM or a callable mask",
     "breadth is what licenses a susceptible call; see callable_mask"),
    ("mapping_quality", "BAM/CRAM",
     "repetitive and paralogous regions produce mismapped false variants"),
    ("mixed_infection", "reads",
     "two strains in one sample invalidate a single consensus interpretation"),
    ("lineage_assignment", "reads or a genome-wide SNP set",
     "lineage confounds variant-phenotype association and catalogue transfer"),
    ("sample_swap_detection", "genotype fingerprint across samples",
     "a swapped sample attaches one patient's resistance to another's record"),
    ("reference_bias", "BAM/CRAM and the reference used",
     "reference-biased alignment loses divergent alleles"),
]


def assess(
    adapted: AdapterResult,
    mask: CallableMask,
    profile: OrganismProfile,
    depth_floor: int = 10,
) -> list[QCFinding]:
    """Build the QC block for a report, including unperformed controls."""
    findings: list[QCFinding] = list(adapted.qc)

    # -- what we can actually check from a variant file -------------------
    if mask.is_present:
        findings.append(QCFinding(
            "callable_mask", "pass",
            f"callable-locus evidence supplied from {mask.source}"))
    else:
        findings.append(QCFinding(
            "callable_mask", "fail",
            "no callable-locus evidence supplied, so no drug can be reported "
            "susceptible. Supply a depth table, BED mask or gVCF."))

    declared = adapted.assembly
    if profile.reference_assembly and declared != profile.reference_assembly:
        findings.append(QCFinding(
            "reference_assembly", "fail",
            f"input declares assembly {declared!r} but the "
            f"{profile.name} profile is written against "
            f"{profile.reference_assembly!r}. Coordinates are not comparable "
            f"and liftover is not implemented."))
    else:
        findings.append(QCFinding(
            "reference_assembly", "pass",
            f"input matches the profile assembly ({declared})"))

    if profile.loci_with_unknown_length:
        findings.append(QCFinding(
            "profile_annotation", "warn",
            f"{len(profile.loci_with_unknown_length)} locus/loci in the "
            f"{profile.name} profile have no declared length, so callable "
            f"fractions must be supplied explicitly by the input. Reference "
            f"annotation arrives with an ingested catalogue."))

    if not profile.ships_validated:
        findings.append(QCFinding(
            "profile_validation", "warn",
            f"the {profile.name} profile v{profile.version} ships with no "
            f"validation data. No sensitivity, specificity or error rate has "
            f"been established for it."))

    # -- what we cannot check, named explicitly ---------------------------
    for check, needs, why in _READ_LEVEL_CONTROLS:
        findings.append(QCFinding(
            check, "not_performed",
            f"requires {needs}; not available from a variant file. Matters "
            f"because {why}."))

    return findings


def blocking(findings: list[QCFinding]) -> list[QCFinding]:
    return [f for f in findings if f.blocking]


def summarise(findings: list[QCFinding]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for f in findings:
        counts[f.status] = counts.get(f.status, 0) + 1
    return counts
