"""Descriptive allele-frequency groups; never infer clones from VAF proximity."""
import math
from dataclasses import replace
from statistics import median

from ..core.models import PopulationStructure, Subpopulation


def cluster_subpopulations(variants, heteroresistance_findings, tolerance=0.05):
    if not math.isfinite(tolerance) or not 0 < tolerance <= 0.25:
        raise ValueError("VAF grouping tolerance must be in (0,0.25]")
    assessed = {h.variant_key for h in heteroresistance_findings if h.assessable}
    unique = {}
    for v in variants:
        if v.key() in unique and unique[v.key()].vaf != v.vaf:
            raise ValueError("conflicting VAF observations for one variant")
        unique[v.key()] = v
    usable = sorted((v for v in unique.values() if v.key() in assessed
                     and v.vaf is not None and math.isfinite(v.vaf)
                     and 0 <= v.vaf <= 1 and v.passed_filters),
                    key=lambda v: (v.vaf, v.key()))
    groups = []
    for v in usable:
        # Bound total span, preventing transitive chaining of distant VAFs.
        if not groups or v.vaf - groups[-1][0].vaf > tolerance:
            groups.append([])
        groups[-1].append(v)
    return [Subpopulation(frozenset(v.key() for v in group), median(v.vaf for v in group),
                          note="VAF proximity only; shared cellular origin is unproven.")
            for group in groups if len(group) >= 2]


def classify_structure(subpopulations, lineage_hints=()):
    """Different externally linked lineages support a possible mixture only."""
    hints = {frozenset(h["variant_keys"]): h for h in lineage_hints}
    supported = [hints.get(p.variant_keys) for p in subpopulations]
    supported = [h for h in supported if h and h.get("lineage") and h.get("source") and h.get("linkage_source")]
    if len({h["lineage"] for h in supported}) >= 2:
        return "possible_coinfection", "Distinct lineage assignments with supplied linkage provenance; verify contamination and specimen identity."
    return "indeterminate", "VAF groups and shared lineage cannot distinguish clones, same-lineage infection, or within-host evolution."


def population_structure(sample_id, variants, findings, tolerance=0.05, lineage_hints=()):
    if len({h.variant_key for h in findings}) < 2:
        return None
    groups = cluster_subpopulations(variants, findings, tolerance)
    seen = set()
    for hint in lineage_hints:
        for name in ("variant_keys", "lineage", "source", "linkage_source"):
            if not hint.get(name):
                raise ValueError(f"population lineage hint requires {name}")
        key = frozenset(hint["variant_keys"])
        if key in seen or key not in {g.variant_keys for g in groups}:
            raise ValueError("lineage hints must match one unique frequency group exactly")
        seen.add(key)
        groups = [replace(g, lineage=hint["lineage"], lineage_source=hint["source"],
                          linkage_source=hint["linkage_source"]) if g.variant_keys == key else g for g in groups]
    classification, basis = classify_structure(groups, lineage_hints)
    grouped = set().union(*(g.variant_keys for g in groups)) if groups else set()
    gaps = ["Frequency proximity does not establish co-occurrence or cellular linkage.",
            "Fractions are median allele frequencies, not estimated clone proportions."]
    if not lineage_hints:
        gaps.append("No per-subpopulation lineage/linkage evidence supplied; classification requires independent evidence.")
    return PopulationStructure(sample_id, groups, sorted({v.key() for v in variants} - grouped),
                               classification, basis, gaps, tolerance=tolerance)
