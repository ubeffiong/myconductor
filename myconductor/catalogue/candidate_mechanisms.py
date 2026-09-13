"""Research-note archive, deliberately independent of the reconciliation engine."""
import json
from pathlib import Path


class CandidateMechanismCatalogue:
    def __init__(self, entries=(), version="empty-1"):
        self.entries = list(entries)
        self.version = version
        if not version:
            raise ValueError("candidate mechanism catalogue requires a version")
        if len({e["id"] for e in self.entries}) != len(self.entries):
            raise ValueError("duplicate candidate mechanism ID")
        for entry in self.entries:
            if entry.get("tier") not in {"HYPOTHESIS", "EMERGING", "ESTABLISHED"}:
                raise ValueError("unknown research mechanism tier")
            if not all(entry.get(k) for k in ("id", "source", "note", "reviewed_by", "reviewed_at")):
                raise ValueError("candidate mechanism record requires provenance and review")

    @classmethod
    def load(cls, path):
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data.get("schema") != "myconductor.candidate-mechanisms.v1":
            raise ValueError("unknown candidate-mechanism schema")
        return cls(data["entries"], data["version"])
