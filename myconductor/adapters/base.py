"""The common evidence schema, and cross-engine reconciliation.

Every adapter normalises its tool's output into ``list[DrugEvidence]`` carrying
assertion, drug, determinant, source engine and version, database and version,
evidence grade, confidence, and the limitations that apply to that particular
assertion. That shared shape is what makes reconciliation possible at all.

Reporting disagreement is the point
-----------------------------------
Two engines differing on rifampicin — because they read different catalogue
versions, or handle a minority allele differently — is among the most useful
things this platform can surface. ``concordance`` reports agreement and
disagreement per drug; it does not resolve conflicts by voting, because the
majority of three tools sharing one catalogue is not three pieces of evidence.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from ..core.models import (
    Call,
    Discordance,
    DrugEvidence,
    EngineRef,
    QCFinding,
    Variant,
)
from ..io.callable_mask import CallableMask


class AdapterSchemaError(ValueError):
    """The tool's output did not match the schema this adapter understands.

    Raised rather than guessing. Schema drift between tool versions is the
    normal case, and a mis-parse produces a wrong clinical call silently.
    """


@dataclass
class EngineReport:
    """One external engine's contribution, in Myconductor's own vocabulary."""

    engine: EngineRef
    sample_id: Optional[str] = None
    evidence: list[DrugEvidence] = field(default_factory=list)
    variants: list[Variant] = field(default_factory=list)
    lineage: Optional[str] = None
    species: Optional[str] = None
    qc: list[QCFinding] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    #: Coverage evidence the engine itself computed, where it exposes any.
    mask: Optional[CallableMask] = None

    @property
    def drugs(self) -> set[str]:
        return {ev.drug for ev in self.evidence}

    def call_for(self, drug: str) -> Optional[Call]:
        """This engine's strongest position on one drug."""
        calls = [ev.call for ev in self.evidence if ev.drug == drug]
        if not calls:
            return None
        for preferred in (Call.RESISTANT, Call.INDETERMINATE, Call.NO_CALL,
                          Call.SUSCEPTIBLE, Call.NOT_ASSESSED):
            if preferred in calls:
                return preferred
        return calls[0]


class EngineAdapter(ABC):
    """Parses one tool's output into an ``EngineReport``."""

    engine: EngineRef

    @abstractmethod
    def parse(self, path: str | Path) -> EngineReport:
        raise NotImplementedError

    def __str__(self) -> str:  # pragma: no cover - diagnostics only
        return str(self.engine)


# -- reconciliation -------------------------------------------------------
@dataclass
class ConcordanceSummary:
    """Per-drug agreement across engines."""

    agreed: dict[str, str] = field(default_factory=dict)
    disagreed: list[Discordance] = field(default_factory=list)
    single_source: dict[str, str] = field(default_factory=dict)
    unassessed: list[str] = field(default_factory=list)

    @property
    def n_compared(self) -> int:
        return len(self.agreed) + len(self.disagreed)

    @property
    def agreement_rate(self) -> Optional[float]:
        """Only defined where at least two engines assessed the same drug."""
        if self.n_compared == 0:
            return None
        return round(len(self.agreed) / self.n_compared, 3)

    def describe(self) -> str:
        if self.n_compared == 0:
            return ("no drug was assessed by more than one engine; nothing to "
                    "reconcile")
        return (f"{len(self.agreed)}/{self.n_compared} drug(s) concordant "
                f"across engines; {len(self.disagreed)} discordant")


def concordance(reports: Iterable[EngineReport],
                drugs: Optional[Iterable[str]] = None) -> ConcordanceSummary:
    """Compare engines drug by drug, without resolving conflicts by vote."""
    reports = list(reports)
    universe = set(drugs) if drugs else {d for r in reports for d in r.drugs}
    summary = ConcordanceSummary()

    for drug in sorted(universe):
        positions = {
            r.engine.name: r.call_for(drug)
            for r in reports
            if r.call_for(drug) is not None
        }
        if not positions:
            summary.unassessed.append(drug)
        elif len(positions) == 1:
            name, call = next(iter(positions.items()))
            summary.single_source[drug] = f"{call.value} ({name} only)"
        elif len(set(positions.values())) == 1:
            summary.agreed[drug] = next(iter(positions.values())).value
        else:
            summary.disagreed.append(Discordance(
                drug=drug,
                calls=tuple(dict.fromkeys(c.value for c in positions.values())),
                sources=tuple(positions),
                note="; ".join(f"{n} says {c.value}"
                               for n, c in sorted(positions.items()))
                + ". Not resolved by vote: engines sharing a catalogue are not "
                  "independent evidence.",
            ))
    return summary


def cross_engine_discordance(reports: Iterable[EngineReport]) -> list[Discordance]:
    return concordance(reports).disagreed


def merge(reports: Iterable[EngineReport]) -> list[DrugEvidence]:
    """Flatten every engine's evidence into one list for reconciliation."""
    return [ev for r in reports for ev in r.evidence]
