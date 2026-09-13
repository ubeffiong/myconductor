"""Explicit specimen identity, assay context and deployment validation scope.

No validation approvals are bundled. A deployment supplies its own evidence
references for a specific organism, assay, drug, engine and database version.
This records an external validation decision; it does not perform or certify it.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Optional

from .models import EngineRef, QCFinding


@dataclass
class SampleContext:
    sample_id: str
    isolate_id: str
    site_id: str
    organism: str
    assay: str = "unknown"
    specimen_id: Optional[str] = None
    collected_at: Optional[str] = None
    lineage: Optional[str] = None
    subspecies: Optional[str] = None
    # Exact aliases supplied by the laboratory, never inferred from filenames.
    engine_sample_ids: dict[str, str] = field(default_factory=dict)
    qc: list[QCFinding] = field(default_factory=list)
    erm41_status: str = "unknown"  # functional | nonfunctional | unknown
    erm41_source: Optional[str] = None
    rrl_status: str = "unknown"    # wild_type | resistance_variant | unknown
    rrl_source: Optional[str] = None
    rrs_status: str = "unknown"
    rrs_source: Optional[str] = None

    def __post_init__(self) -> None:
        for key in ("sample_id", "isolate_id", "site_id", "organism", "assay"):
            if not isinstance(getattr(self, key), str) or not getattr(self, key).strip():
                raise ValueError(f"context {key} must be a nonempty string")
        for name, allowed in (
            ("erm41", {"functional", "nonfunctional", "unknown"}),
            ("rrl", {"wild_type", "resistance_variant", "unknown"}),
            ("rrs", {"wild_type", "resistance_variant", "unknown"}),
        ):
            status = getattr(self, name + "_status")
            if status not in allowed:
                raise ValueError(f"invalid {name} status: {status}")
            if status != "unknown" and not getattr(self, name + "_source"):
                raise ValueError(f"{name} status requires an evidence source")
        names = [q.check for q in self.qc]
        if len(names) != len(set(names)):
            raise ValueError("context contains duplicate QC checks")
        if any(q.status not in {"pass", "fail", "warn", "not_performed"} for q in self.qc):
            raise ValueError("invalid QC status")

    @classmethod
    def from_dict(cls, data: dict) -> "SampleContext":
        data = dict(data)
        data["qc"] = [QCFinding(**q) for q in data.get("qc", [])]
        return cls(**data)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class ValidationScope:
    organism: str
    assay: str
    drug: str
    engine: str
    engine_version: str
    database_version: str
    evidence_reference: str
    reviewed_by: str
    reviewed_at: str
    # Local catalogue predictions must additionally match the catalogue digest.
    catalogue_sha256: Optional[str] = None

    def __post_init__(self) -> None:
        for key in ("organism", "assay", "drug", "engine", "engine_version",
                    "database_version", "evidence_reference", "reviewed_by", "reviewed_at"):
            value = getattr(self, key)
            if not isinstance(value, str) or value.strip().lower() in ("", "unknown", "*"):
                raise ValueError(f"validation scope requires a concrete {key}")


@dataclass
class InterpretationPolicy:
    version: str = "research-default"
    scopes: list[ValidationScope] = field(default_factory=list)
    required_qc: tuple[str, ...] = (
        "species_confirmation", "contamination", "mapping_quality", "mixed_infection",
    )

    @classmethod
    def from_dict(cls, data: dict) -> "InterpretationPolicy":
        if set(data) - {"version", "scopes", "required_qc"}:
            raise ValueError("unknown interpretation policy fields")
        if not data.get("version"):
            raise ValueError("policy requires a version")
        if not isinstance(data.get("scopes", []), list):
            raise ValueError("policy scopes must be a list")
        required = tuple(data.get("required_qc", cls().required_qc))
        if not {"species_confirmation", "contamination"} <= set(required):
            raise ValueError("policy must require species_confirmation and contamination")
        return cls(data["version"], [ValidationScope(**s) for s in data.get("scopes", [])],
                   required)

    def permits(self, drug: str, context: Optional[SampleContext],
                engine: Optional[EngineRef], catalogue_sha256: Optional[str] = None) -> bool:
        if context is None or engine is None:
            return False
        statuses = {q.check: q.status for q in context.qc}
        if any(statuses.get(q) != "pass" for q in self.required_qc):
            return False
        return any(
            s.organism == context.organism and s.assay == context.assay
            and s.drug == drug and s.engine == engine.name
            and s.engine_version == engine.version
            and s.database_version == engine.database_version
            and (catalogue_sha256 is None or s.catalogue_sha256 == catalogue_sha256)
            for s in self.scopes
        )
