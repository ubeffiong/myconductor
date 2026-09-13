"""Target-discovery research workflow — cohort-level, never per-patient.

Why this is no longer wired into the pipeline
---------------------------------------------
The previous design fired a drug-discovery funnel whenever one patient had no
adequate regimen, and printed docking scores in that patient's report. Three
things were wrong with it.

It was scientifically weak: subtractive genomics and docking operate at the
level of a pathogen, a lineage or a resistance mechanism, not a single isolate.
One patient's unexplained resistance is a sample size of one.

It conflated evidence with speculation on the same page. A clinician reading a
susceptibility table should not encounter candidate compounds beside it.

And the docking scores were fabricated — ``-11.0 + 4.0 * hash(gene)`` — while
being formatted to two decimal places in kcal/mol, which is the most
persuasive possible presentation of a number that means nothing. Docking scores
alone establish neither inhibition, permeability, whole-cell activity,
selectivity, toxicity, nor resistance barrier, even when they are real.

So discovery is now a **separate workflow** that consumes aggregated cohort
evidence and emits a laboratory prioritisation brief. It is not imported by
``core.pipeline``, and it cannot appear in a patient report. Docking requires an
explicitly injected real backend; there is no built-in scorer.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Protocol

from ..core.models import EngineRef


class DiscoveryInputError(ValueError):
    """The cohort evidence supplied is too thin to justify a discovery run."""


@dataclass
class CohortSignal:
    """The aggregated evidence that legitimately motivates target discovery."""

    cohort_id: str
    n_isolates: int
    n_sites: int
    unexplained_resistant_isolates: int
    drug: str
    convergent_genes: list[str] = field(default_factory=list)
    lineages: list[str] = field(default_factory=list)
    note: str = ""

    @property
    def unexplained_fraction(self) -> float:
        return (self.unexplained_resistant_isolates / self.n_isolates
                if self.n_isolates else 0.0)


@dataclass
class DockingResult:
    best_kcal_mol: float
    top_compound: str
    engine: EngineRef
    caveats: tuple[str, ...] = ()


@dataclass
class DiscoveryTarget:
    protein: str
    gene: str
    essential: Optional[bool] = None
    human_homolog: Optional[bool] = None
    druggable: Optional[bool] = None
    evidence_sources: list[str] = field(default_factory=list)
    docking: Optional[DockingResult] = None

    @property
    def passes_subtractive_filter(self) -> Optional[bool]:
        """None where any input annotation is unknown — not False.

        An unknown essentiality is not evidence of dispensability, and
        collapsing the two would quietly drop real candidates or admit bad ones.
        """
        if None in (self.essential, self.human_homolog, self.druggable):
            return None
        return bool(self.essential and not self.human_homolog and self.druggable)


@dataclass
class DiscoveryBrief:
    """A laboratory prioritisation brief. Not a clinical document."""

    cohort_id: str
    signal: CohortSignal
    targets: list[DiscoveryTarget] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    audience: str = "laboratory prioritisation — not for clinical use"

    @property
    def ranked_targets(self) -> list[DiscoveryTarget]:
        return [t for t in self.targets if t.passes_subtractive_filter]

    @property
    def indeterminate_targets(self) -> list[DiscoveryTarget]:
        return [t for t in self.targets if t.passes_subtractive_filter is None]


class DockingBackend(Protocol):
    """A real docking engine. There is no built-in implementation."""

    engine: EngineRef

    def dock(self, gene: str) -> DockingResult:
        ...


class ProteomeSource(Protocol):
    """Supplies candidate proteins with essentiality and homology annotation.

    A real implementation cross-references the Database of Essential Genes,
    BLASTP against the human proteome, and pocket detection over PDB or
    AlphaFold structures.
    """

    def candidates(self, genes: list[str]) -> list[DiscoveryTarget]:
        ...


#: Minimum cohort evidence before a discovery run is scientifically meaningful.
MIN_ISOLATES = 20
MIN_SITES = 2
MIN_UNEXPLAINED = 5


def run(
    signal: CohortSignal,
    proteome: ProteomeSource,
    docking: Optional[DockingBackend] = None,
    min_isolates: int = MIN_ISOLATES,
    min_sites: int = MIN_SITES,
    min_unexplained: int = MIN_UNEXPLAINED,
) -> DiscoveryBrief:
    """Build a discovery brief from cohort-level evidence.

    Raises rather than producing a thin brief: a discovery run motivated by one
    isolate is exactly the failure this module was rewritten to prevent.
    """
    problems = []
    if signal.n_isolates < min_isolates:
        problems.append(f"{signal.n_isolates} isolates (need {min_isolates})")
    if signal.n_sites < min_sites:
        problems.append(f"{signal.n_sites} site(s) (need {min_sites})")
    if signal.unexplained_resistant_isolates < min_unexplained:
        problems.append(
            f"{signal.unexplained_resistant_isolates} unexplained resistant "
            f"isolate(s) (need {min_unexplained})")
    if problems:
        raise DiscoveryInputError(
            "cohort evidence is insufficient to motivate target discovery: "
            + "; ".join(problems)
            + ". Target discovery is a pathogen- or mechanism-level activity, "
              "not a per-isolate one."
        )

    if not signal.convergent_genes:
        raise DiscoveryInputError(
            "no convergent genes supplied; without recurrent, independently "
            "arising variants there is no signal to prioritise against."
        )

    targets = proteome.candidates(signal.convergent_genes)

    caveats = [
        "Candidates are hypothesis-generating only and have no therapeutic "
        "standing.",
        "Subtractive-genomics filters reflect the supplied annotation; "
        "unknown essentiality or homology yields an indeterminate verdict, "
        "not a pass.",
    ]

    if docking is not None:
        for t in targets:
            if t.passes_subtractive_filter:
                t.docking = docking.dock(t.gene)
        caveats.append(
            f"Docking performed by {docking.engine}. A binding score "
            f"establishes none of: inhibition, permeability, whole-cell "
            f"activity, selectivity, toxicity, or resistance barrier."
        )
    else:
        caveats.append(
            "No docking backend was supplied, so no affinity is reported. "
            "There is deliberately no built-in scorer: the previous "
            "hash-derived values were fabricated."
        )

    return DiscoveryBrief(
        cohort_id=signal.cohort_id,
        signal=signal,
        targets=targets,
        caveats=caveats,
    )
