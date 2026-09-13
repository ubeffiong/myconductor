"""Local discordance tickets and disclosure-gated count summaries."""
import hashlib
from dataclasses import asdict, dataclass, field
from .ledger_store import LedgerStore
from .transport import DisclosureThresholds


@dataclass(frozen=True)
class DiscordanceTicket:
    ticket_id: str
    sample_id: str
    isolate_id: str
    site_id: str
    organism: str
    drug: str
    genomic_call: str
    phenotypic_call: str
    source: str
    timestamp: str
    candidate_mechanisms: list[str] = field(default_factory=list)
    evidence_gaps: list[str] = field(default_factory=list)
    status: str = "open"
    resolution: str = ""

    def __post_init__(self):
        for key in ("ticket_id", "sample_id", "isolate_id", "site_id", "organism", "drug", "source", "timestamp"):
            if not isinstance(getattr(self, key), str) or not getattr(self, key).strip():
                raise ValueError(f"discordance ticket requires {key}")
        if self.genomic_call not in {"resistant", "susceptible"} or self.phenotypic_call not in {"resistant", "susceptible"} or self.genomic_call == self.phenotypic_call:
            raise ValueError("ticket requires opposing categorical genomic/phenotypic calls")
        if self.status != "open" or self.resolution:
            raise ValueError("new ticket must be open; use ledger transitions")


class WatchlistStore(LedgerStore):
    schema = "myconductor.watchlist.v1"

    def state(self):
        state = {}
        for entry in self.ledger.entries:
            p = entry.payload
            if entry.action == "opened":
                ticket = DiscordanceTicket(**p)
                if ticket.ticket_id in state:
                    raise ValueError("duplicate watchlist ticket")
                state[ticket.ticket_id] = dict(p)
            else:
                key = p["ticket_id"]
                allowed = {"open": {"under_investigation"}, "under_investigation": {"resolved"},
                           "resolved": {"under_investigation"}}
                if key not in state or entry.action not in allowed[state[key]["status"]]:
                    raise ValueError("invalid watchlist transition")
                if not p.get("reviewer") or not p.get("rationale"):
                    raise ValueError("ticket transition requires reviewer and rationale")
                if entry.action == "resolved" and not p.get("resolution"):
                    raise ValueError("ticket resolution must be explicit")
                state[key].update(p, status=entry.action, timestamp=entry.timestamp_utc)
        return state

    def capture(self, report):
        tickets = []
        for result in report.drug_results:
            a, b = result.genomic_call, result.phenotypic_call
            if a is None or b is None or a.value not in {"resistant", "susceptible"} or b.value not in {"resistant", "susceptible"} or a == b:
                continue
            context = report.context
            source = report.provenance.analysis_fingerprint
            key = "|".join((context["site_id"], context["isolate_id"], result.drug, source))
            ticket_id = hashlib.sha256(key.encode()).hexdigest()
            ticket = DiscordanceTicket(ticket_id, report.sample_id, context["isolate_id"], context["site_id"],
                context["organism"], result.drug, a.value, b.value, source, report.provenance.generated_utc,
                evidence_gaps=["Review assay identity, QC, method and categorical disagreement; no mechanism is inferred."])
            if ticket_id not in self.state():
                self.append("opened", asdict(ticket))
            tickets.append(self.state()[ticket_id])
        return tickets

    def transition(self, ticket_id, action, reviewer, rationale, resolution=""):
        self.append(action, dict(ticket_id=ticket_id, reviewer=reviewer, rationale=rationale, resolution=resolution))


def aggregate(tickets_by_site, thresholds=None):
    """Release count cells only, never an isolate-level list under a k label."""
    t = thresholds or DisclosureThresholds()
    if t.min_isolates_per_cell < 2 or t.min_sites_per_cell < 2:
        raise ValueError("watchlist disclosure requires at least two isolates and two sites")
    cells = {}
    for site, tickets in tickets_by_site.items():
        for ticket in tickets:
            if ticket["site_id"] != site:
                raise ValueError("ticket site differs from aggregation envelope")
            if ticket["status"] == "resolved":
                continue
            cells.setdefault(ticket["drug"], set()).add((site, ticket["isolate_id"]))
    released, suppressed = [], False
    for drug, identities in sorted(cells.items()):
        sites = {site for site, _ in identities}
        if len(identities) < t.min_isolates_per_cell or len(sites) < t.min_sites_per_cell:
            suppressed = True
            continue
        released.append({"drug": drug, "unresolved_isolates": len(identities), "contributing_sites": len(sites)})
    return {"schema": "myconductor.watchlist-aggregate.v1", "released": released,
            "has_suppressed_cells": suppressed,
            "limitations": "Count thresholds are a disclosure floor, not anonymity. No variant, patient, site identifier or ticket text is exported."}
