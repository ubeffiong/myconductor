"""Heteroresistance detection.

A single patient can carry a mix of susceptible and resistant subpopulations.
Standard consensus calling misses low-frequency resistant minorities, producing
false-negative resistance reports. This module scans variant allele frequencies
and flags resistance-associated variants present below the majority threshold
but above the sequencing noise floor.

It runs across *all* resistance evidence, independent of which lane produced it,
so a minority katG S315T (catalogued) and a minority Rv0678 change (efflux lane)
are both surfaced.
"""
from __future__ import annotations

from ..core.models import DrugEvidence, HeteroresistanceFlag, Variant


def detect(
    variants: list[Variant],
    resistance_evidence: list[DrugEvidence],
    majority_floor: float = 0.9,
    noise_floor: float = 0.03,
) -> list[HeteroresistanceFlag]:
    by_key = {v.key(): v for v in variants}
    flags: list[HeteroresistanceFlag] = []
    for ev in resistance_evidence:
        if not ev.call.is_resistant:
            continue
        var = by_key.get(ev.variant_key)
        if var is None:
            continue
        if noise_floor <= var.vaf < majority_floor:
            flags.append(
                HeteroresistanceFlag(
                    variant_key=var.key(),
                    drug=ev.drug,
                    vaf=var.vaf,
                    depth=var.depth,
                    note=(
                        f"Minority resistant subpopulation at {var.vaf:.0%} "
                        f"(depth {var.depth}). A consensus caller would likely "
                        f"miss this and report {ev.drug} susceptible."
                    ),
                )
            )
    return flags
