"""Shared local persistence for append-only governance records."""
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from .catalogue_update import EvidenceLedger, LedgerEntry
from ..reporting.json_report import atomic_json


class LedgerStore:
    schema = "myconductor.ledger.v1"

    def __init__(self):
        self.ledger = EvidenceLedger()
        self.file_hash = None

    def validate(self):
        ok, reason = self.ledger.verify()
        if not ok:
            raise ValueError(f"ledger integrity failed: {reason}")
        self.state()

    def state(self):
        raise NotImplementedError

    def to_dict(self):
        return {"schema": self.schema, "ledger": [asdict(e) for e in self.ledger.entries]}

    @classmethod
    def load(cls, path):
        obj = cls()
        path = Path(path)
        if path.exists():
            raw = path.read_bytes()
            data = json.loads(raw)
            if data.get("schema") != cls.schema:
                raise ValueError("unknown ledger schema")
            obj.ledger.entries = [LedgerEntry(**e) for e in data["ledger"]]
            obj.validate()
            obj.file_hash = hashlib.sha256(raw).hexdigest()
        return obj

    def append(self, action, payload):
        import copy
        self.ledger.append(action, copy.deepcopy(payload))
        try:
            self.validate()
        except Exception:
            self.ledger.entries.pop()
            raise

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        lock = path.with_name(path.name + ".lock")
        try:
            stream = lock.open("x")
        except FileExistsError as exc:
            raise ValueError("another writer holds the ledger lock") from exc
        try:
            current = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
            if current != self.file_hash:
                raise ValueError("ledger changed; reload before saving")
            self.validate()
            atomic_json(path, self.to_dict())
            self.file_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        finally:
            stream.close()
            lock.unlink()
