"""Stable, versioned research report serialization and atomic local writes."""
from __future__ import annotations

import dataclasses
import enum
import json
import os
import tempfile
from pathlib import Path


def primitive(value):
    if isinstance(value, enum.Enum):
        return value.value
    if dataclasses.is_dataclass(value):
        return {f.name: primitive(getattr(value, f.name)) for f in dataclasses.fields(value)}
    if isinstance(value, dict):
        return {str(k): primitive(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [primitive(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return [primitive(v) for v in sorted(value)]
    if hasattr(value, "_loci") and hasattr(value, "source"):
        return {"source": value.source, "loci": primitive(value._loci)}
    return value


def report_dict(report):
    import hashlib
    from ..modules.investigation import evidence_id
    data = primitive(report)
    for raw, result in zip(data["drug_results"], report.drug_results):
        for record, evidence in zip(raw["evidence"], result.evidence):
            record["evidence_id"] = evidence_id(evidence)
    data["report_sha256"] = hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return data


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def load_report(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_report(data)
    return data


def validate_report(data):
    import hashlib
    content = {k: v for k, v in data.items() if k != "report_sha256"}
    digest = hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if data.get("report_sha256") != digest:
        raise ValueError("report content checksum mismatch")
    if data.get("report_schema") != "myconductor.report.v1":
        raise ValueError("unsupported report schema")
    if not data.get("sample_id") or not isinstance(data.get("drug_results"), list):
        raise ValueError("invalid report")
    manifest = data.get("analysis_manifest")
    if not isinstance(manifest, dict) or not manifest:
        raise ValueError("report requires its frozen analysis manifest")
    expected = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    if expected != data.get("provenance", {}).get("analysis_fingerprint"):
        raise ValueError("report analysis fingerprint does not match its manifest")
    if manifest.get("context") != data.get("context"):
        raise ValueError("report context differs from frozen manifest")
    if manifest.get("sources") != data["provenance"].get("source_hashes"):
        raise ValueError("report sources differ from frozen manifest")
    seen = set()
    for result in data["drug_results"]:
        if result["drug"] in seen:
            raise ValueError("report has duplicate drug results")
        seen.add(result["drug"])
        for evidence in result["evidence"]:
            raw = {k: v for k, v in evidence.items() if k != "evidence_id"}
            digest = hashlib.sha256(json.dumps(raw, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
            if digest != evidence.get("evidence_id"):
                raise ValueError("report evidence fingerprint mismatch")
