"""NCBI AMRFinderPlus adapter — acquired determinants, for non-TB profiles.

Written from AMRFinderPlus's documented TSV output. **Not validated against
real output.**

The honesty constraint that shapes this adapter
-----------------------------------------------
AMRFinderPlus reports a drug **class** or subclass (``BETA-LACTAM``,
``AMINOGLYCOSIDE``, ``QUINOLONE``), not a specific agent. Turning "aminoglycoside
resistance gene detected" into "amikacin resistant" is an interpretation step
with real clinical consequences, and it is not one this adapter performs. Every
piece of evidence it emits therefore names the class as its drug and carries a
limitation saying so; mapping class to agent belongs to an organism profile with
its own validation.

It also matters that acquired determinants are a different kind of object from
chromosomal point mutations: gene presence, allele identity, coverage of the
reference sequence and mobile context all bear on interpretation. Those fields
are preserved on the evidence rather than flattened away.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Optional

from ..core.models import (
    Call,
    Consequence,
    DrugEvidence,
    EngineRef,
    Lane,
    QCFinding,
    Tier,
    Variant,
    VariantIdentity,
)
from .base import AdapterSchemaError, EngineAdapter, EngineReport

_REQUIRED_COLUMNS = ("gene symbol", "element type")

_CLASS_LIMITATION = (
    "AMRFinderPlus reports a drug class, not a specific agent; mapping class to "
    "agent requires a validated organism profile and is not done here"
)


class AMRFinderPlusAdapter(EngineAdapter):
    def __init__(self, version: str = "unknown", database_version: str = "unknown",
                 min_identity: float = 90.0, min_coverage: float = 90.0):
        self.engine = EngineRef(
            name="amrfinderplus", version=version,
            database="NCBI Reference Gene Catalog",
            database_version=database_version,
        )
        self.min_identity = min_identity
        self.min_coverage = min_coverage

    def parse(self, path: str | Path) -> EngineReport:
        path = Path(path)
        text = path.read_text(encoding="utf-8-sig")
        reader = csv.DictReader(text.splitlines(), delimiter="\t")
        if not reader.fieldnames:
            raise AdapterSchemaError(f"{path}: no header row")

        lowered = {h.strip().lower(): h for h in reader.fieldnames}
        missing = [c for c in _REQUIRED_COLUMNS if c not in lowered]
        if missing:
            raise AdapterSchemaError(
                f"{path}: missing required column(s) {missing}. "
                f"Columns found: {sorted(lowered)[:20]}. This adapter "
                f"understands AMRFinderPlus TSV output."
            )

        report = EngineReport(engine=self.engine, source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        n_rows = 0
        n_virulence = 0

        for row in reader:
            n_rows += 1
            element_type = _get(row, lowered, "element type", "").upper()
            if element_type != "AMR":
                n_virulence += 1
                continue

            gene = _get(row, lowered, "gene symbol", "").strip()
            if not gene:
                report.warnings.append("row with empty Gene symbol skipped")
                continue

            drug_class = (_get(row, lowered, "subclass", "")
                          or _get(row, lowered, "class", "")).strip().lower()
            if not drug_class:
                report.warnings.append(
                    f"{gene}: no Class/Subclass; cannot attribute to a drug class")
                continue

            identity_pct = _as_float(_get(row, lowered,
                                          "% identity to reference sequence"))
            coverage_pct = _as_float(_get(row, lowered,
                                          "% coverage of reference sequence"))
            method = _get(row, lowered, "method", "").strip()
            subtype = _get(row, lowered, "element subtype", "").strip().upper()

            limitations = [f"imported from {self.engine}", _CLASS_LIMITATION]
            call = Call.INDETERMINATE
            tier = Tier.CATALOGUED

            if identity_pct is not None and identity_pct < self.min_identity:
                call, tier = Call.INDETERMINATE, Tier.NONE
                limitations.append(
                    f"{identity_pct:.1f}% identity to the reference allele, "
                    f"below the {self.min_identity:.0f}% threshold")
            if coverage_pct is not None and coverage_pct < self.min_coverage:
                call, tier = Call.INDETERMINATE, Tier.NONE
                limitations.append(
                    f"{coverage_pct:.1f}% coverage of the reference sequence, "
                    f"below the {self.min_coverage:.0f}% threshold; a partial "
                    f"gene may not be functional")
            if subtype == "POINT":
                limitations.append("point mutation, not an acquired gene")
            else:
                limitations.append(
                    "acquired gene: presence alone does not establish "
                    "expression, and mobile context is not assessed")

            variant = Variant(identity=VariantIdentity(
                gene=gene,
                chrom=_get(row, lowered, "contig id") or None,
                pos=_as_int(_get(row, lowered, "start")),
                hgvs_p=_get(row, lowered, "sequence name") or gene,
                consequence=(Consequence.MISSENSE if subtype == "POINT"
                             else Consequence.UNKNOWN),
            ))
            report.variants.append(variant)

            report.evidence.append(DrugEvidence(
                drug=drug_class,
                call=call,
                tier=tier,
                lane=Lane.ENGINE, scope="determinant",
                confidence=None,
                variant=variant.identity,
                engine=self.engine,
                limitations=tuple(limitations),
                rationale=(
                    f"AMRFinderPlus detected {gene} "
                    f"({subtype or element_type}) associated with {drug_class}; "
                    f"phenotypic resistance is not predicted"
                    + (f"; {identity_pct:.1f}% identity" if identity_pct else "")
                    + (f", {coverage_pct:.1f}% coverage" if coverage_pct else "")
                    + f", method {method or 'unreported'}."
                ),
            ))

        report.qc.append(QCFinding(
            "engine_rows", "pass",
            f"AMRFinderPlus: {n_rows} row(s), {len(report.evidence)} AMR "
            f"determinant(s), {n_virulence} non-AMR element(s) ignored"))
        report.qc.append(QCFinding(
            "drug_class_granularity", "warn",
            "evidence is at drug-class granularity; agent-level calls require "
            "a validated organism profile mapping"))
        return report


def _get(row: dict, lowered: dict[str, str], key: str,
         default: str = "") -> str:
    column = lowered.get(key)
    if column is None:
        return default
    value = row.get(column)
    return default if value is None else str(value)


def _as_float(value: str) -> Optional[float]:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _as_int(value: str) -> Optional[int]:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None
