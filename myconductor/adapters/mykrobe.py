"""Mykrobe adapter.

Written from Mykrobe's documented JSON predict schema. **Not validated against
real output.**

Mykrobe is the one engine whose susceptible call this adapter imports directly,
because Mykrobe distinguishes ``S`` (assessed, susceptible) from ``N`` (not
enough coverage to call). That distinction is exactly the one Myconductor
enforces everywhere else, so an ``S`` from Mykrobe is a coverage-backed claim
and is imported with ``asserts_coverage=True`` — attributed to Mykrobe in the
report, so a reader can see whose coverage logic licensed it.

``r`` (lower case) is Mykrobe's minority-allele resistance call. It is imported
as resistance with the minority basis recorded as a limitation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ..core.models import (
    Call,
    DrugEvidence,
    EngineRef,
    Lane,
    QCFinding,
    Tier,
    Variant,
    VariantIdentity,
)
from .base import AdapterSchemaError, EngineAdapter, EngineReport

#: Mykrobe's prediction codes.
_PREDICT = {
    "R": (Call.RESISTANT, "resistant"),
    "r": (Call.RESISTANT, "resistant from a minority allele"),
    "S": (Call.SUSCEPTIBLE, "susceptible with adequate coverage"),
    "N": (Call.NOT_ASSESSED, "insufficient coverage to call"),
    "U": (Call.INDETERMINATE, "unknown / uncalled"),
}


class MykrobeAdapter(EngineAdapter):
    def __init__(self, version: str = "unknown", panel: str = "unknown"):
        self.engine = EngineRef(
            name="mykrobe", version=version,
            database="mykrobe-panel", database_version=panel,
        )

    def parse(self, path: str | Path) -> EngineReport:
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise AdapterSchemaError(f"{path}: not valid JSON") from exc

        if not isinstance(data, dict) or not data:
            raise AdapterSchemaError(
                f"{path}: expected a JSON object keyed by sample name")

        # Mykrobe nests everything under the sample name.
        sample_id, payload = next(iter(data.items()))
        if not isinstance(payload, dict):
            raise AdapterSchemaError(
                f"{path}: top-level value for {sample_id!r} is not an object")

        susceptibility = payload.get("susceptibility")
        if not isinstance(susceptibility, dict):
            raise AdapterSchemaError(
                f"{path}: no 'susceptibility' object under {sample_id!r}. "
                f"Keys present: {sorted(payload)[:12]}. This adapter "
                f"understands `mykrobe predict --format json` output."
            )

        engine = EngineRef(
            name=self.engine.name,
            version=str(payload.get("mykrobe_version", self.engine.version)),
            database=self.engine.database,
            database_version=str(payload.get("probe_sets",
                                             self.engine.database_version)),
        )

        report = EngineReport(
            engine=engine,
            sample_id=str(sample_id),
            lineage=_lineage_of(payload),
            species=_species_of(payload),
        )

        for drug_name, entry in susceptibility.items():
            if not isinstance(entry, dict):
                report.warnings.append(f"{drug_name}: unexpected entry shape")
                continue
            code = str(entry.get("predict", "")).strip()
            mapping = _PREDICT.get(code)
            if mapping is None:
                report.warnings.append(
                    f"{drug_name}: unrecognised predict code {code!r}; skipped "
                    f"rather than guessed")
                continue
            call, description = mapping
            drug = drug_name.strip().lower()

            called_by = entry.get("called_by") or {}
            variants = [_variant_of(name) for name in called_by]
            variants = [v for v in variants if v is not None]
            report.variants.extend(variants)

            limitations = [f"imported from {engine}"]
            if code == "r":
                limitations.append(
                    "minority-allele call; read-level confirmation required")
            if code == "S":
                limitations.append(
                    "susceptibility asserted by Mykrobe's own coverage "
                    "handling (it reports 'N' where coverage is insufficient), "
                    "not by a Myconductor callable mask")

            report.evidence.append(DrugEvidence(
                drug=drug,
                call=call,
                tier=Tier.CATALOGUED if code in ("R", "r", "S") else Tier.NONE,
                lane=Lane.ENGINE,
                confidence=None,
                variant=variants[0].identity if variants else None,
                engine=engine,
                limitations=tuple(limitations),
                asserts_coverage=(code == "S"),
                rationale=(
                    f"Mykrobe predicts {code} ({description}) for {drug}"
                    + (f", called by {', '.join(sorted(called_by))}"
                       if called_by else "")
                    + "."
                ),
            ))

        report.qc.extend(_qc_of(payload, report))
        return report


def _variant_of(called_by: str) -> Optional[Variant]:
    """Parse a Mykrobe called_by key such as ``katG_S315T-GC2155168GA``."""
    head = called_by.split("-", 1)[0]
    if "_" not in head:
        return None
    gene, change = head.split("_", 1)
    if not gene or not change:
        return None
    return Variant(identity=VariantIdentity(gene=gene, hgvs_p=change))


def _lineage_of(payload: dict) -> Optional[str]:
    phylo = payload.get("phylogenetics")
    if not isinstance(phylo, dict):
        return None
    lineage = phylo.get("lineage")
    if isinstance(lineage, dict):
        names = lineage.get("lineage")
        if isinstance(names, list) and names:
            return str(names[-1])
        return None
    if isinstance(lineage, str):
        return lineage
    return None


def _species_of(payload: dict) -> Optional[str]:
    phylo = payload.get("phylogenetics")
    if not isinstance(phylo, dict):
        return None
    species = phylo.get("species")
    if isinstance(species, dict) and species:
        return str(next(iter(species)))
    return None


def _qc_of(payload: dict, report: EngineReport) -> list[QCFinding]:
    findings = []
    if report.species:
        findings.append(QCFinding(
            "engine_species", "pass",
            f"Mykrobe identified species: {report.species}"))
    if report.lineage:
        findings.append(QCFinding(
            "engine_lineage", "pass",
            f"Mykrobe lineage: {report.lineage}"))
    asserted = [ev.drug for ev in report.evidence if ev.asserts_coverage]
    if asserted:
        findings.append(QCFinding(
            "engine_coverage_assertion", "pass",
            f"Mykrobe asserts adequate coverage for: {', '.join(sorted(asserted))}"))
    return findings
