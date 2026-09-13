"""Explain unresolved results and prioritize laboratory review.

Priorities are explicit workflow rules, not calibrated probabilities or
claims that a particular test will resolve a case. Costs/availability come
only from the laboratory's own menu. No test is ordered by this module.
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict

from ..core.models import Call, Tier


def evidence_id(evidence) -> str:
    body = json.dumps(asdict(evidence), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def investigate(report) -> list[dict]:
    findings = []

    def add(drug, category, certainty, explanation, actions, evidence=(), priority=2):
        key = f"{report.sample_id}|{drug}|{category}"
        findings.append({
            "id": hashlib.sha256(key.encode()).hexdigest()[:24],
            "drug": drug, "category": category, "certainty": certainty,
            "explanation": explanation, "candidate_actions": list(actions),
            "evidence_ids": [evidence_id(e) for e in evidence],
            "priority": priority, "rule_version": "investigation-1",
        })

    failures = [q.check for q in report.qc if q.status == "fail" and q.check != "callable_mask"]
    if failures:
        add("*", "analytical_failure", "observed", "; ".join(failures),
            ["review_sample_and_qc"], priority=0)
    for result in report.drug_results:
        evs = result.evidence
        if result.discordance:
            phenotypes = [e for e in evs if e.tier is Tier.PHENOTYPIC]
            if phenotypes:
                add(result.drug, "phenotype_discordance", "observed",
                    "Conflicting reported conclusions; neither mechanism nor measurement error "
                    "is established as the cause. Review assay methods, timing and QC.",
                    ["review_phenotype", "confirm_phenotype"], evs, 0)
            else:
                versions = {e.catalogue_version or (e.engine.database_version if e.engine else None)
                            for e in evs}
                versions.discard(None)
                versions.discard("unknown")
                add(result.drug, "engine_discordance", "observed", result.discordance.note,
                    ["review_engine_evidence"], evs, 0)
                if len(versions) > 1:
                    add(result.drug, "catalogue_drift", "possible",
                        "Different catalogue versions coexist. Replay identical input against "
                        "both catalogues before attributing the disagreement to an update.",
                        ["replay_catalogues"], evs, 1)
        if result.assay_status == "insufficient" and result.phenotypic_call is None:
            add(result.drug, "coverage_gap", "observed",
                "Required genomic loci are not independently shown callable.",
                ["review_coverage", "repeat_sequencing"], evs, 1)
        if any(e.scope == "historical" for e in evs):
            add(result.drug, "historical_association", "observed",
                "Prior isolate measurements are association evidence, not this isolate's phenotype.",
                ["review_association", "confirm_phenotype"], evs, 2)
        uncertain = [e for e in evs if e.call is Call.INDETERMINATE
                     and e.scope != "historical" and e.tier is not Tier.PHENOTYPIC]
        if uncertain:
            add(result.drug, "uncertain_interpretation", "observed",
                "Uninterpreted genomic evidence requires expert assessment; resistance is not inferred.",
                ["review_variant", "confirm_phenotype"], uncertain, 2)
        if result.genomic_call is Call.INDETERMINATE and not uncertain and not result.discordance:
            add(result.drug, "validation_gap", "observed",
                result.reason or "No matching externally validated interpretation scope.",
                ["review_validation_scope"], evs, 2)
        if (report.provenance.organism_profile == "mabscessus"
                and result.drug == "clarithromycin"
                and not result.call.is_established):
            add(result.drug, "ntm_induction_context", "observed",
                "Review subspecies, functional erm(41), rrl evidence, and the laboratory's "
                "incubation/inducible-resistance assessment. An early susceptible result "
                "does not exclude inducible resistance.",
                ["review_ntm_context", "assess_inducible_resistance"], evs, 1)
    return sorted(findings, key=lambda f: (f["priority"], f["drug"], f["category"]))


def plan_follow_up(findings, menu=(), budget=None):
    """Deterministic priority-first allocation; not an optimal test strategy.

    At most one available action per finding is selected. A shared action on
    the same drug is charged once. Unknown costs never become zero.
    """
    if budget is not None and (not math.isfinite(budget) or budget < 0):
        raise ValueError("follow-up budget must be finite and nonnegative")
    options = []
    currencies = set()
    ids = set()
    allowed = {a for f in findings for a in f["candidate_actions"]}
    for option in menu:
        required = {"id", "action", "available"}
        if not required <= set(option) or not isinstance(option["available"], bool):
            raise ValueError("test menu requires id, action and boolean available")
        if option["id"] in ids:
            raise ValueError("duplicate test-menu id")
        ids.add(option["id"])
        for field in ("cost", "turnaround_hours"):
            value = option.get(field)
            if value is not None and (not math.isfinite(value) or value < 0):
                raise ValueError(f"menu {field} must be finite and nonnegative")
        if option.get("cost") is not None:
            if not option.get("currency"):
                raise ValueError("a cost requires a currency")
            currencies.add(option["currency"])
        if option["action"] in allowed:
            options.append(dict(option))
    if budget is not None and len(currencies) > 1:
        raise ValueError("budget allocation requires one currency; no currency conversion is inferred")
    remaining = budget
    selected = set()
    out = []
    for finding in sorted(findings, key=lambda f: (f["priority"], f["drug"], f["id"])):
        matches = [o for o in options if o["action"] in finding["candidate_actions"]
                   and (not o.get("drugs") or finding["drug"] in o["drugs"])]
        available = [o for o in matches if o["available"]]
        available.sort(key=lambda o: (
            (finding["drug"], o["id"]) not in selected,
            o.get("cost") is None, o.get("cost", 0) or 0,
            o.get("turnaround_hours") is None, o.get("turnaround_hours", 0) or 0, o["id"],
        ))
        choice = None
        for option in available:
            reused = (finding["drug"], option["id"]) in selected
            if budget is None or reused or (option.get("cost") is not None
                                           and option["cost"] <= remaining):
                choice = option
                break
        status = "recommended"
        charged = None
        if choice:
            key = (finding["drug"], choice["id"])
            charged = 0 if key in selected else choice.get("cost")
            selected.add(key)
            if remaining is not None:
                remaining -= charged
        else:
            status = ("deferred_budget" if available and any(o.get("cost") is not None for o in available)
                      else "cost_unknown" if available
                      else "unavailable" if matches else "availability_unknown")
        out.append({
            "investigation_id": finding["id"], "drug": finding["drug"],
            "priority": finding["priority"], "status": status, "selected": choice,
            "alternatives": matches, "charged_cost": charged,
            "remaining_budget": remaining,
            "basis": "priority-first heuristic; no measured probability of resolution",
            "requires_laboratory_review": True,
        })
    return out
