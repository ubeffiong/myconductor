"""Feature engineering for the VUS classifier.

The whole point of the classifier lane is that you do NOT feed raw genomic
letters to the model -- you convert each variant into biological features. This
module defines that feature contract. The demo computes features
deterministically from the variant so the pipeline runs offline; each function
is a labelled hook where a real annotator (SIFT/PolyPhen, an AlphaFold structure
service, a Grantham matrix) would be wired in.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass

from ..core.models import Region, Variant

# Genes with well-characterised resistance-determining regions. A substitution
# landing in one of these is prior-weighted as more likely functional.
_RESISTANCE_REGION_GENES = {
    "katG", "rpoB", "gyrA", "gyrB", "pncA", "embB", "rrs", "atpE",
    "Rv0678", "ddn", "rplC", "inhA",
}


@dataclass
class FeatureVector:
    conservation: float           # 0..1, higher = more conserved
    structural_impact: float      # 0..1, higher = more destabilising
    substitution_severity: float  # 0..1, higher = more chemically radical
    in_resistance_region: float   # 0 or 1


def _stable_unit(seed: str) -> float:
    """Deterministic pseudo-value in [0,1] from a string.

    Stands in for a real score so the demo is reproducible without network or
    model files. Replace the body with a real annotator call.
    """
    digest = hashlib.sha256(seed.encode()).hexdigest()
    return int(digest[:8], 16) / 0xFFFFFFFF


def conservation_score(variant: Variant) -> float:
    # Hook for SIFT / PolyPhen-2 / phyloP over the codon.
    return _stable_unit("cons:" + variant.key())


def structural_impact_score(variant: Variant) -> float:
    # Hook for AlphaFold-derived binding-pocket / stability (ddG) delta.
    return _stable_unit("struct:" + variant.key())


def substitution_severity(variant: Variant) -> float:
    # Hook for a Grantham/BLOSUM-style physicochemical distance.
    # Nonsense / frameshift changes are maximally severe.
    change = variant.change.lower()
    if "stop" in change or "*" in change or "fs" in change or "del" in change:
        return 1.0
    return _stable_unit("subst:" + variant.key())


def extract_features(variant: Variant) -> FeatureVector:
    in_region = 1.0 if variant.gene in _RESISTANCE_REGION_GENES else 0.0
    # Non-coding coding-region features are damped: a promoter SNP is handled by
    # the regulatory lane, so if one reaches here we lower its structural weight.
    struct = structural_impact_score(variant)
    if variant.region in (Region.PROMOTER, Region.INTERGENIC):
        struct *= 0.3
    return FeatureVector(
        conservation=conservation_score(variant),
        structural_impact=struct,
        substitution_severity=substitution_severity(variant),
        in_resistance_region=in_region,
    )
