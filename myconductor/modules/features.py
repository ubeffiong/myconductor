"""Feature annotation for the VUS workbench.

What changed and why
--------------------
The previous implementation computed every "feature" as a SHA-256 digest of the
variant's own name:

    conservation      = sha256("cons:katG_S315T")  / 0xFFFFFFFF
    structural_impact = sha256("struct:katG_S315T") / 0xFFFFFFFF

Those numbers are reproducible, which made them look like scores. They carry no
biological information whatsoever, and feeding them to a scorer produced a
probability with no meaning — then explaining that probability with SHAP-style
contributions manufactured confidence in noise. That is worse than reporting
nothing.

So the default annotator is now ``NullAnnotator``, which reports every dimension
as **unavailable**. The workbench renders "no data source configured" instead of
a number. Unavailable is a legitimate, useful output; an invented score is not.

``DemoAnnotator`` retains the hash behaviour for interface demonstrations, but
it must be constructed explicitly, it stamps ``synthetic=True`` on everything it
returns, and the pipeline refuses to use it unless ``demo_mode`` is set. It can
never reach a report that is not labelled synthetic.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Optional, Protocol

from ..core.models import Consequence, Dimension, Variant

#: Genes with well-characterised resistance-determining regions. Membership is
#: a real, citable fact about the gene, so it is the one dimension the null
#: annotator can populate without a data source.
RESISTANCE_REGION_GENES = frozenset({
    "katG", "rpoB", "gyrA", "gyrB", "pncA", "embB", "rrs", "rrl", "atpE",
    "Rv0678", "ddn", "rplC", "inhA", "fabG1", "eis", "fbiC", "fgd1",
})


@dataclass
class FeatureSet:
    """Annotation dimensions for one variant, each with availability attached."""

    dimensions: list[Dimension] = field(default_factory=list)
    synthetic: bool = False

    def get(self, name: str) -> Optional[Dimension]:
        for d in self.dimensions:
            if d.name == name:
                return d
        return None

    @property
    def available_names(self) -> list[str]:
        return [d.name for d in self.dimensions if d.available]

    @property
    def gaps(self) -> list[str]:
        return [d.name for d in self.dimensions if not d.available]


class AnnotatorProtocol(Protocol):
    """Anything that can describe a variant's biology.

    A real implementation wires in SIFT/PolyPhen or phyloP for conservation, an
    AlphaFold or PDB structure service for stability and ligand distance, a
    Grantham/BLOSUM matrix for substitution severity, and a population database
    for prevalence, lineage distribution and homoplasy.
    """

    synthetic: bool

    def annotate(self, variant: Variant) -> FeatureSet:
        ...


#: The dimensions the workbench asks for. Named here so a partial annotator
#: still produces a complete, gap-labelled report rather than a short one.
DIMENSION_NAMES = (
    "in_resistance_region",
    "consequence",
    "conservation",
    "structural_impact",
    "ligand_distance",
    "substitution_severity",
    "population_prevalence",
    "lineage_distribution",
    "homoplasy",
    "mic_association",
    "cooccurring_known_variants",
    "literature_support",
)


class NullAnnotator:
    """The honest default: reports what is knowable without a data source.

    Two dimensions are derivable from the variant itself — whether its gene has
    a characterised resistance-determining region, and its consequence class.
    Everything else requires an external database, and is reported unavailable.
    """

    synthetic = False
    name = "null"

    def annotate(self, variant: Variant) -> FeatureSet:
        dims = [
            Dimension(
                name="in_resistance_region",
                available=True,
                value=variant.gene in RESISTANCE_REGION_GENES,
                source="drug_loci.json gene list",
            ),
            Dimension(
                name="consequence",
                available=variant.consequence is not Consequence.UNKNOWN,
                value=variant.consequence.value,
                source="input annotation",
                note=("input carried no consequence annotation"
                      if variant.consequence is Consequence.UNKNOWN else ""),
            ),
        ]
        for name in DIMENSION_NAMES:
            if name in ("in_resistance_region", "consequence"):
                continue
            dims.append(Dimension(
                name=name,
                available=False,
                source="none",
                note=_REQUIRES[name],
            ))
        return FeatureSet(dimensions=dims, synthetic=False)


_REQUIRES = {
    "conservation": "requires a conservation service (SIFT, PolyPhen-2, phyloP)",
    "structural_impact": "requires a structure source (PDB or AlphaFold) and a ddG model",
    "ligand_distance": "requires a structure with the ligand or active site annotated",
    "substitution_severity": "requires a Grantham or BLOSUM matrix and a resolved substitution",
    "population_prevalence": "requires an isolate population database",
    "lineage_distribution": "requires lineage-typed isolate counts",
    "homoplasy": "requires a phylogeny and ancestral-state reconstruction",
    "mic_association": "requires paired MIC measurements",
    "cooccurring_known_variants": "requires isolate-level co-occurrence counts",
    "literature_support": "requires a curated literature index",
}


class DemoAnnotator:
    """Hash-derived values, for exercising the interface only.

    Every dimension it returns is marked synthetic, and ``FeatureSet.synthetic``
    is True, which the pipeline checks before allowing this into a report.
    These numbers are not scores. They are deterministic noise.
    """

    synthetic = True
    name = "demo-hash"

    def __init__(self, acknowledged: bool = False):
        if not acknowledged:
            raise ValueError(
                "DemoAnnotator produces hash-derived values with no biological "
                "meaning. Construct it as DemoAnnotator(acknowledged=True) to "
                "confirm the output is for interface demonstration only."
            )

    def annotate(self, variant: Variant) -> FeatureSet:
        base = NullAnnotator().annotate(variant)
        dims = []
        for d in base.dimensions:
            if d.available:
                dims.append(d)
                continue
            dims.append(Dimension(
                name=d.name,
                available=True,
                value=round(_stable_unit(f"{d.name}:{variant.key()}"), 3),
                source="demo-hash (SYNTHETIC — not a measurement)",
                note="deterministic hash of the variant name; no biological meaning",
            ))
        return FeatureSet(dimensions=dims, synthetic=True)


def _stable_unit(seed: str) -> float:
    """Deterministic pseudo-value in [0, 1). Not a score. See DemoAnnotator."""
    digest = hashlib.sha256(seed.encode()).hexdigest()
    return int(digest[:8], 16) / 0x100000000
