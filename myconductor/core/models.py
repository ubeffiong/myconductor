"""Core domain model for Myconductor.

Everything downstream (router, modules, reporting) speaks in these types, so a
new module only has to accept a ``Variant`` and return a ``DrugEvidence``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class Call(str, Enum):
    """Per-drug resistance verdict.

    The distinction between a *catalogued* call and a *predicted* call is
    deliberate and load-bearing: a clinician must always be able to tell an
    evidence-graded WHO catalogue hit from a model prediction.
    """

    RESISTANT = "resistant"
    SUSCEPTIBLE = "susceptible"
    PREDICTED_RESISTANT = "predicted_resistant"
    PREDICTED_SUSCEPTIBLE = "predicted_susceptible"
    INDETERMINATE = "indeterminate"

    @property
    def is_resistant(self) -> bool:
        return self in (Call.RESISTANT, Call.PREDICTED_RESISTANT)

    @property
    def is_predicted(self) -> bool:
        return self in (Call.PREDICTED_RESISTANT, Call.PREDICTED_SUSCEPTIBLE)


class Route(str, Enum):
    """Which lane of the router handled a variant."""

    CATALOGUE = "catalogue"          # known WHO-graded mutation
    VUS_ML = "vus_ml"                # coding variant of unknown significance
    EFFLUX_REGULATORY = "efflux"     # efflux locus or non-coding/promoter
    IGNORED = "ignored"              # silent / synonymous / off-target


class Region(str, Enum):
    CODING = "coding"
    PROMOTER = "promoter"
    INTERGENIC = "intergenic"
    RRNA = "rrna"


@dataclass
class Variant:
    """A canonical variant, normalised from any input format."""

    gene: str
    change: str                       # e.g. "S315T", "c.-15C>T", "a1401g"
    genomic_pos: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    vaf: float = 1.0                  # variant allele frequency, 0..1
    depth: int = 0                    # read depth at the locus
    region: Region = Region.CODING
    silent: bool = False

    def key(self) -> str:
        return f"{self.gene}_{self.change}"


@dataclass
class DrugEvidence:
    """One module's opinion about one drug, backed by one variant."""

    drug: str
    call: Call
    confidence: float                 # 0..1
    route: Route
    variant_key: str
    rationale: str
    who_grade: Optional[str] = None   # e.g. "Assoc w R", "Not assoc w R"


@dataclass
class DrugResult:
    """The reconciled verdict for a drug after all evidence is weighed."""

    drug: str
    call: Call
    confidence: float
    evidence: list[DrugEvidence] = field(default_factory=list)

    @property
    def predicted(self) -> bool:
        return self.call.is_predicted


@dataclass
class HeteroresistanceFlag:
    variant_key: str
    drug: str
    vaf: float
    depth: int
    note: str


@dataclass
class Regimen:
    proposed: list[str]
    effective_drugs: list[str]
    ineffective_drugs: list[str]
    adequate: bool
    rationale: str


@dataclass
class DiscoveryTarget:
    protein: str
    gene: str
    essential: bool
    human_homolog: bool
    druggable: bool
    best_docking_kcal_mol: Optional[float] = None
    top_compound: Optional[str] = None


@dataclass
class DiscoveryResult:
    triggered: bool
    reason: str
    targets: list[DiscoveryTarget] = field(default_factory=list)


@dataclass
class Provenance:
    tool: str
    version: str
    catalogue_version: str
    vus_model: str
    coverage_min: int
    het_vaf_floor: float


@dataclass
class AnalysisReport:
    sample_id: str
    drug_results: list[DrugResult]
    heteroresistance: list[HeteroresistanceFlag]
    regimen: Regimen
    discovery: DiscoveryResult
    provenance: Provenance
    qc_warnings: list[str] = field(default_factory=list)
    routed_counts: dict[str, int] = field(default_factory=dict)
