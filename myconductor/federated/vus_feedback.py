"""Closed-loop VUS validation — the laboratory feedback path.

The problem
-----------
The WHO catalogue grades many variants group 3 (uncertain significance) or
group 4 (not associated, interim). Published surveillance of bedaquiline,
clofazimine, delamanid and pretomanid genes across Central and West Africa has
found that most group-3/4 variants there show no phenotypic resistance —
meaning the uncertainty is not just a knowledge gap, it is a live source of
overtreatment risk while it persists. ``modules.vus_workbench`` already ranks
these variants for validation instead of guessing at them. What was missing
is the other half of the loop: when a laboratory actually validates one, that
result has nowhere to go. It is not global catalogue evidence — one isolate's
DST does not get to reclassify a WHO grade — but it is real, local knowledge
that should stop the same variant, at the same site, being re-flagged as
unresolved the next time it is seen.

What this module is, and is not
--------------------------------
This is a **local evidence tier**, not a catalogue promotion path.
``federated/catalogue_update.py`` already has the correct, much stricter
machinery for that (isolate-level denominators, lineage/co-occurrence/
site-diversity confounding checks, expert review, hash-chained ledger) and
this module reuses its ``EvidenceLedger`` directly rather than inventing a
second audit trail. A validation recorded here never edits the global
catalogue and never bypasses that review queue; it only informs Myconductor's
own site of what its own laboratory has already confirmed, with a retraction
path if a result turns out to be wrong.

Reversibility
-------------
Every ingestion and retraction is a ledger entry. ``retract`` does not delete
a record — a validated result standing today may be found wrong tomorrow (an
isolate mislabelled, a mixed culture, a second DST that disagrees), and the
audit trail must show that it changed rather than silently vanish.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Optional

from ..core.models import Call
from .catalogue_update import EvidenceLedger


class ValidationMethod(str, Enum):
    PHENOTYPIC_DST = "phenotypic_dst"
    MIC = "mic"
    EFFLUX_ASSAY = "efflux_assay"
    ALLELIC_EXCHANGE = "allelic_exchange"
    RNA_EXPRESSION = "rna_expression"
    OTHER = "other"


class VUSFeedbackError(ValueError):
    """A validation record or a store operation is invalid."""


@dataclass(frozen=True)
class VUSValidationRecord:
    """One laboratory's validation of one variant, for one drug, on one isolate.

    ``result`` must resolve to an established call (``RESISTANT`` or
    ``SUSCEPTIBLE``) — a validation experiment that itself came back
    inconclusive is not evidence to feed back in, and constructing this with
    anything else raises rather than silently accepting an unresolved result.
    """

    variant_key: str
    #: Display label (``gene_change``), e.g. ``"Rv0678_L114R"``. Needed
    #: separately from ``variant_key`` (which is coordinate-based where
    #: possible) because ``modules.vus_workbench`` excludes VUS candidates by
    #: label, not by coordinate key — see ``VUSWorkbench.known``.
    variant_label: str
    gene: str
    drug: str
    isolate_id: str
    site_id: str
    method: ValidationMethod
    result: Call
    lineage: Optional[str] = None
    mic: Optional[float] = None
    submitted_by: str = "unknown"
    rationale: str = ""
    timestamp_utc: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def __post_init__(self) -> None:
        if not self.result.is_established:
            raise VUSFeedbackError(
                f"a validation record must resolve to RESISTANT or "
                f"SUSCEPTIBLE; got {self.result.value!r} for "
                f"{self.variant_key}/{self.drug}. An inconclusive experiment "
                f"is not feedback to ingest — repeat it, or record it "
                f"elsewhere as a lab note.")

    def to_dict(self) -> dict:
        return {
            "variant_key": self.variant_key, "variant_label": self.variant_label,
            "gene": self.gene,
            "drug": self.drug, "isolate_id": self.isolate_id,
            "site_id": self.site_id, "method": self.method.value,
            "result": self.result.value, "lineage": self.lineage,
            "mic": self.mic, "submitted_by": self.submitted_by,
            "rationale": self.rationale, "timestamp_utc": self.timestamp_utc,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "VUSValidationRecord":
        return cls(
            variant_key=data["variant_key"],
            variant_label=data.get("variant_label", data["variant_key"]),
            gene=data["gene"],
            drug=data["drug"], isolate_id=data["isolate_id"],
            site_id=data["site_id"], method=ValidationMethod(data["method"]),
            result=Call(data["result"]), lineage=data.get("lineage"),
            mic=data.get("mic"), submitted_by=data.get("submitted_by", "unknown"),
            rationale=data.get("rationale", ""),
            timestamp_utc=data.get("timestamp_utc", ""),
        )


@dataclass
class LocalValidationVerdict:
    """The site's current position on one variant/drug, from its own records."""

    variant_key: str
    drug: str
    call: Call                      # RESISTANT/SUSCEPTIBLE if consistent, else INDETERMINATE
    records: list[VUSValidationRecord]

    @property
    def n_isolates(self) -> int:
        return len({r.isolate_id for r in self.records})

    @property
    def consistent(self) -> bool:
        return self.call is not Call.INDETERMINATE

    @property
    def lineages(self) -> list[str]:
        return sorted({r.lineage for r in self.records if r.lineage})

    def describe(self) -> str:
        if self.consistent:
            return (f"{self.n_isolates} isolate(s) at "
                    f"{len({r.site_id for r in self.records})} site(s) "
                    f"consistently validated {self.call.value}")
        calls = sorted({r.result.value for r in self.records})
        return (f"{self.n_isolates} isolate(s) give conflicting results "
                f"({'/'.join(calls)}); held as INDETERMINATE pending "
                f"resolution")


class LocalValidationStore:
    """Site-local ledger of validated VUS results.

    Not a catalogue and not a substitute for ``federated.catalogue_update``'s
    review queue — a smaller, purely local, immediately-reversible record of
    what this site's own laboratory has confirmed about variants it has seen.
    """

    def __init__(self, ledger: Optional[EvidenceLedger] = None):
        self.ledger = ledger or EvidenceLedger()
        self._records: dict[tuple[str, str], list[VUSValidationRecord]] = {}
        self._retracted: dict[tuple[str, str], list[int]] = {}

    # -- ingestion ----------------------------------------------------------
    def ingest(self, record: VUSValidationRecord) -> dict:
        key = (record.variant_key, record.drug)
        self._records.setdefault(key, []).append(record)
        entry = self.ledger.append("validated", record.to_dict())
        return {"entry_hash": entry.entry_hash, "seq": entry.seq}

    def retract(self, variant_key: str, drug: str, record_index: int,
                reviewer: str, rationale: str) -> dict:
        """Reverse one contributing record without deleting its history.

        ``record_index`` is the position of the record within
        ``records_for(variant_key, drug)`` (the order they were ingested).
        A retracted record no longer contributes to :meth:`verdict`, but stays
        in ``self._records`` and in the ledger, so a reader can see that a
        result was withdrawn and why, not just that it is now absent.
        """
        key = (variant_key, drug)
        records = self._records.get(key, [])
        if not records or record_index < 0 or record_index >= len(records):
            raise VUSFeedbackError(
                f"no record at index {record_index} for {variant_key}/{drug}")
        if record_index in self._retracted.get(key, []):
            raise VUSFeedbackError(
                f"record {record_index} for {variant_key}/{drug} is already "
                f"retracted")
        self._retracted.setdefault(key, []).append(record_index)
        entry = self.ledger.append("retracted", {
            "variant_key": variant_key, "drug": drug,
            "record_index": record_index, "reviewer": reviewer,
            "rationale": rationale,
        })
        return {"entry_hash": entry.entry_hash, "seq": entry.seq}

    # -- reading --------------------------------------------------------------
    def records_for(self, variant_key: str, drug: str) -> list[VUSValidationRecord]:
        key = (variant_key, drug)
        retracted = set(self._retracted.get(key, []))
        return [r for i, r in enumerate(self._records.get(key, []))
                if i not in retracted]

    def verdict(self, variant_key: str, drug: str) -> Optional[LocalValidationVerdict]:
        records = self.records_for(variant_key, drug)
        if not records:
            return None
        calls = {r.result for r in records}
        call = calls.pop() if len(calls) == 1 else Call.INDETERMINATE
        return LocalValidationVerdict(variant_key=variant_key, drug=drug,
                                      call=call, records=records)

    def drugs_for(self, variant_key: str) -> list[str]:
        """Every drug with at least one non-retracted record for this variant."""
        return sorted({drug for (vk, drug) in self._records
                      if vk == variant_key and self.records_for(vk, drug)})

    def validated_labels(self) -> set[str]:
        """Every variant display label with at least one non-retracted record.

        Wired into ``VUSWorkbench.known`` (which excludes by label, not
        coordinate key) so a locally-resolved VUS drops out of the "needs
        validation" ranking — it already has an answer here.
        """
        out: set[str] = set()
        for (variant_key, drug), records in self._records.items():
            if self.records_for(variant_key, drug):
                out.update(r.variant_label for r in records)
        return out

    def ledger_verified(self) -> tuple[bool, Optional[str]]:
        return self.ledger.verify()

    # -- persistence ----------------------------------------------------------
    def to_json(self) -> dict:
        return {
            "records": {
                f"{vk}␟{drug}": [r.to_dict() for r in records]
                for (vk, drug), records in self._records.items()
            },
            "retracted": {
                f"{vk}␟{drug}": indices
                for (vk, drug), indices in self._retracted.items()
            },
            "ledger": [
                {"seq": e.seq, "timestamp_utc": e.timestamp_utc,
                 "action": e.action, "payload": e.payload,
                 "prev_hash": e.prev_hash, "entry_hash": e.entry_hash}
                for e in self.ledger.entries
            ],
        }

    def save(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_json(), indent=2) + "\n", encoding="utf-8")

    @classmethod
    def from_json(cls, data: dict) -> "LocalValidationStore":
        store = cls()
        for compound_key, rows in data.get("records", {}).items():
            vk, _, drug = compound_key.partition("␟")
            store._records[(vk, drug)] = [VUSValidationRecord.from_dict(r)
                                          for r in rows]
        for compound_key, indices in data.get("retracted", {}).items():
            vk, _, drug = compound_key.partition("␟")
            store._retracted[(vk, drug)] = list(indices)
        from .catalogue_update import LedgerEntry
        store.ledger.entries = [
            LedgerEntry(seq=e["seq"], timestamp_utc=e["timestamp_utc"],
                       action=e["action"], payload=e["payload"],
                       prev_hash=e["prev_hash"], entry_hash=e["entry_hash"])
            for e in data.get("ledger", [])
        ]
        return store

    @classmethod
    def load(cls, path: str | Path) -> "LocalValidationStore":
        return cls.from_json(json.loads(Path(path).read_text(encoding="utf-8")))
