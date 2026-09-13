"""Local case review ledger. Review outcomes never silently amend reports.

Hash chaining detects modification relative to a retained ledger head; it is
not authentication or protection against wholesale replacement by an attacker.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from .catalogue_update import EvidenceLedger, LedgerEntry
from ..reporting.json_report import atomic_json, validate_report


def _evidence_ids(report):
    return sorted({e["evidence_id"] for r in report["drug_results"] for e in r["evidence"]}
                  | {"report:" + report["provenance"]["analysis_fingerprint"]})


class CaseReviewStore:
    def __init__(self):
        self.ledger = EvidenceLedger()
        self.path = None
        self.file_hash = None

    def cases(self):
        cases = {}
        for e in self.ledger.entries:
            p = e.payload
            if e.action == "opened":
                cases[p["case_id"]] = dict(p, state="open")
            elif e.action == "attached":
                case = cases[p["case_id"]]
                case["evidence_ids"] = sorted(set(case["evidence_ids"]) | set(_evidence_ids(p["report"])))
                case["latest_report"] = p["report"]
            else:
                cases[p["case_id"]].update(p, state=e.action)
        return cases

    def open(self, report, reviewer):
        validate_report(report)
        if not reviewer.strip():
            raise ValueError("reviewer required")
        fingerprint = report["provenance"].get("analysis_fingerprint")
        if not fingerprint:
            raise ValueError("case requires a fingerprinted report")
        case_id = hashlib.sha256((report["sample_id"] + "|" + fingerprint).encode()).hexdigest()[:24]
        if case_id not in self.cases():
            self.ledger.append("opened", {
                "case_id": case_id, "sample_id": report["sample_id"],
                "analysis_fingerprint": fingerprint, "reviewer": reviewer,
                "evidence_ids": _evidence_ids(report),
                "report": report,
            })
        return case_id

    def attach(self, case_id, report, reviewer, rationale):
        validate_report(report)
        case = self.cases().get(case_id)
        if not case or case["state"] == "resolved":
            raise ValueError("attach requires an open case; reopen resolved cases first")
        if not reviewer.strip() or not rationale.strip():
            raise ValueError("reviewer and rationale required")
        previous = case["report"]
        for key in ("sample_id", "isolate_id", "site_id", "organism"):
            if previous["context"].get(key) != report["context"].get(key):
                raise ValueError(f"attached report has different {key}")
        self.ledger.append("attached", {"case_id": case_id, "report": report,
                                       "reviewer": reviewer, "rationale": rationale})

    def transition(self, case_id, state, reviewer, rationale, evidence_ids=(),
                   resolution=None, reviewer_minutes=None):
        transitions = {"open": {"in_review"}, "in_review": {"resolved"},
                       "resolved": {"reopened"}, "reopened": {"in_review"}}
        case = self.cases().get(case_id)
        if not case or state not in transitions[case["state"]]:
            raise ValueError("invalid case transition")
        if not reviewer.strip() or not rationale.strip():
            raise ValueError("reviewer and rationale required")
        evidence_ids = sorted(set(evidence_ids))
        if set(evidence_ids) - set(case["evidence_ids"]):
            raise ValueError("review cites evidence not present in the frozen report")
        if state == "resolved" and (not evidence_ids or not resolution):
            raise ValueError("resolution requires cited evidence and a written disposition")
        if reviewer_minutes is not None:
            import math
            if not math.isfinite(reviewer_minutes) or reviewer_minutes < 0:
                raise ValueError("reviewer minutes must be finite and nonnegative")
        self.ledger.append(state, {
            "case_id": case_id, "reviewer": reviewer, "rationale": rationale,
            "cited_evidence_ids": evidence_ids, "resolution": resolution,
            "reviewer_minutes": reviewer_minutes,
        })

    @classmethod
    def load(cls, path):
        store = cls()
        store.path = Path(path)
        if not store.path.exists():
            return store
        raw = store.path.read_bytes()
        data = json.loads(raw)
        if data.get("schema") != "myconductor.case-review.v1":
            raise ValueError("unknown case-review schema")
        store.ledger.entries = [LedgerEntry(**e) for e in data["ledger"]]
        ok, reason = store.ledger.verify()
        if not ok:
            raise ValueError(f"case ledger integrity failed: {reason}")
        # Replay transitions, not a mutable cached case table.
        replay = cls()
        for i, e in enumerate(store.ledger.entries):
            if e.seq != i:
                raise ValueError("case ledger sequence is not contiguous")
            p = e.payload
            if e.action == "opened":
                if p["case_id"] in replay.cases():
                    raise ValueError("duplicate case opening")
                expected = replay.open(p["report"], p["reviewer"])
                if expected != p["case_id"]:
                    raise ValueError("case ID does not match report fingerprint")
                if replay.ledger.entries[-1].payload != p:
                    raise ValueError("case evidence index differs from frozen report")
            elif e.action == "attached":
                replay.attach(p["case_id"], p["report"], p["reviewer"], p["rationale"])
            else:
                replay.transition(p["case_id"], e.action, p["reviewer"], p["rationale"],
                                  p["cited_evidence_ids"], p["resolution"], p["reviewer_minutes"])
        store.file_hash = hashlib.sha256(raw).hexdigest()
        return store

    def save(self, path):
        path = Path(path)
        lock = path.with_name(path.name + ".lock")
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            stream = lock.open("x")
        except FileExistsError as exc:
            raise ValueError("another writer holds the case-store lock") from exc
        try:
            current = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
            if current != self.file_hash:
                raise ValueError("case store changed; reload before saving")
            atomic_json(path, {"schema": "myconductor.case-review.v1",
                               "ledger": [asdict(e) for e in self.ledger.entries]})
            self.file_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        finally:
            stream.close()
            lock.unlink()
