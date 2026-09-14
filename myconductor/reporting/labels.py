"""Turning internal identifiers into something a clinician can read.

Field names, enum values and check identifiers are written for code:
``not_assessed``, ``mtbc_vs_ntm``, ``_ranking_basis``, ``NOT_ASSESSED``. They
leak into a report because the renderer has them to hand and they are not
obviously wrong — they are merely unreadable, and a reader who cannot parse a
label stops reading the value beside it.

Sentence case is the default, but a blind ``token.replace("_", " ").capitalize()``
mangles the terms that matter most: ``mic_association`` becomes "Mic
association", ``mtbc_vs_ntm`` becomes "Mtbc vs ntm", ``ppv`` becomes "Ppv".
So the acronyms and domain terms are named explicitly and everything else
falls through to the mechanical rule.

The internal token is never discarded — ``humanise`` is for display, and the
exports and the JSON payload keep the original spelling so anything reading
this report programmatically still joins on a stable key.
"""
from __future__ import annotations

import re

#: Terms the mechanical rule gets wrong. Keys are the internal token.
EXPLICIT = {
    # call states
    "resistant": "Resistant",
    "susceptible": "Susceptible",
    "indeterminate": "Indeterminate",
    "not_assessed": "Not assessed",
    "no_call": "No call",
    "unsupported": "Unsupported",
    # tiers
    "catalogued": "Catalogued",
    "phenotypic": "Phenotypic",
    "inferred": "Inferred",
    "predicted": "Predicted",
    "none": "None",
    # statuses
    "not_performed": "Not performed",
    "pass": "Pass",
    "warn": "Warning",
    "fail": "Fail",
    "adequate": "Adequate",
    "insufficient": "Insufficient",
    "insufficient-data": "Insufficient data",
    "high": "High",
    "moderate": "Moderate",
    "low": "Low",
    # QC checks
    "mtbc_vs_ntm": "MTBC vs NTM",
    "callable_mask": "Callable mask",
    "coverage_breadth": "Coverage breadth",
    "input_parsed": "Input parsed",
    "lineage_assignment": "Lineage assignment",
    "mapping_quality": "Mapping quality",
    "mixed_infection": "Mixed infection",
    "profile_annotation": "Profile annotation",
    "profile_validation": "Profile validation",
    "read_quality": "Read quality",
    "record_depth": "Record depth",
    "reference_assembly": "Reference assembly",
    "reference_bias": "Reference bias",
    "sample_swap_detection": "Sample-swap detection",
    "sequencing_platform": "Sequencing platform",
    "species_confirmation": "Species confirmation",
    "contamination": "Contamination",
    "availability_unknown": "Availability unknown",
    # investigation categories
    "coverage_gap": "Coverage gap",
    "uncertain_interpretation": "Uncertain interpretation",
    "validation_gap": "Validation gap",
    # VUS dimensions
    "_ranking_basis": "Ranking basis",
    "in_resistance_region": "In a resistance region",
    "consequence": "Consequence",
    "conservation": "Conservation",
    "cooccurring_known_variants": "Co-occurring known variants",
    "homoplasy": "Homoplasy",
    "ligand_distance": "Ligand distance",
    "lineage_distribution": "Lineage distribution",
    "literature_support": "Literature support",
    "mic_association": "MIC association",
    "population_prevalence": "Population prevalence",
    "structural_impact": "Structural impact",
    "substitution_severity": "Substitution severity",
    # metrics
    "vme": "Very major errors",
    "me": "Major errors",
    "vme_rate": "Very major error rate",
    "me_rate": "Major error rate",
    "ppv": "PPV",
    "npv": "NPV",
    "sensitivity": "Sensitivity",
    "specificity": "Specificity",
    "error_rate": "Error rate",
    "call_rate": "Call rate",
    # evidence lanes
    "efflux_regulatory": "Efflux / regulatory",
    "vus": "VUS workbench",
    "engine": "Engine",
    "catalogue": "Catalogue",
    "phenotype": "Phenotype",
    "coverage": "Coverage",
    # assay kinds
    "targeted-amplicon": "Targeted amplicon",
    "targeted-capture": "Targeted capture",
    "wgs": "Whole genome",
}

#: Acronyms to upper-case wherever the mechanical rule produces them.
ACRONYMS = {"mic", "vme", "me", "ppv", "npv", "qc", "dst", "vcf", "bed",
            "who", "ntm", "mtbc", "tss", "hgvs", "snv", "mnv", "id", "url",
            "json", "tsv", "csv", "ci", "ml", "rna", "dna", "atp", "hth"}


def humanise(token) -> str:
    """``some_internal_name`` -> ``Some internal name``.

    Falls through to a mechanical rule only for tokens not named in
    ``EXPLICIT``, so a new identifier renders readably by default and can be
    corrected by adding one line rather than by changing the renderer.
    """
    if token is None:
        return ""
    raw = str(token).strip()
    if not raw:
        return ""
    if raw in EXPLICIT:
        return EXPLICIT[raw]
    lowered = raw.lower()
    if lowered in EXPLICIT:
        return EXPLICIT[lowered]

    # camelCase -> spaced, then split on the usual separators.
    spaced = re.sub(r"(?<=[a-z0-9])(?=[A-Z])", " ", raw)
    words = [w for w in re.split(r"[\s_\-./]+", spaced) if w]
    if not words:
        return raw

    out = []
    for index, word in enumerate(words):
        if word.lower() in ACRONYMS:
            out.append(word.upper())
        elif word.isupper() and len(word) <= 5:
            out.append(word)                  # already an acronym or a code
        elif index == 0:
            out.append(word[0].upper() + word[1:].lower()
                       if word.isupper() else word[0].upper() + word[1:])
        else:
            out.append(word if any(c.isupper() for c in word[1:])
                       else word.lower())
    return " ".join(out)


def humanise_all(tokens) -> str:
    """A comma-joined, readable list. Empty in, empty out."""
    items = [humanise(t) for t in (tokens or [])]
    return ", ".join(i for i in items if i)


def variant_name(label) -> str:
    """``Rv0678_L117R`` -> ``Rv0678 L117R``.

    Gene symbols are case-sensitive identifiers and are left exactly as
    written; only the separator between gene and change is opened up.
    """
    if not label:
        return ""
    text = str(label)
    return text.replace("_", " ", 1) if "_" in text else text
