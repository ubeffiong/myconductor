"""JSON-LD export: the evidence graph, with its semantics attached.

A report that leaves this system is read by something that did not read the
domain model. The risk is not serialisation — it is that every downstream
consumer wants two states, and this system has six. A consumer that maps
``INDETERMINATE`` to "not resistant", or drops ``NOT_ASSESSED`` rows as empty,
has reconstructed the presumption this codebase exists to remove, outside the
reach of any guard in it.

So the export does three things beyond writing fields:

**It names the states in the payload, not only in a schema.** Every call
carries ``establishesUse``, ``isEstablished`` and a human-readable
``reasonUnestablished``. A consumer filtering on "may this drug be used?" never
has to interpret the enum, and one that ignores the flags has to ignore
something explicit rather than merely fail to notice something implicit.

**It refuses to emit a binary field.** There is no ``resistant: true/false``
anywhere, at any nesting level, and a test asserts it. Offering one would be
the whole problem: the field would be consumed and the other four states
discarded.

**It carries the tier beside the call.** ``RESISTANT`` from a graded catalogue
entry and ``INDETERMINATE`` from a model are different claims about the world,
and a flat export that keeps only the verdict throws away the axis that says
how much to trust it.

Why JSON-LD and not RDF
-----------------------
A ``@context`` gives the fields stable, dereferenceable meaning and costs one
dictionary. A triple store, a query language and an ontology would be real
infrastructure, and the questions this system currently answers are per-sample
and answered by the report object directly. Build the graph when a question
arrives that cannot be answered without one; until then this is the
interoperability that was actually missing.

The vocabulary is local and self-describing. It deliberately does **not** claim
a mapping to SNOMED, LOINC or NCIT: those mappings are a curation task with
clinical consequences, and inventing codes here would be the same failure as
inventing terminology codes in the FHIR bundle, which this project already
refuses.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from ..core.models import AnalysisReport, Call, DrugResult, Tier

#: Local vocabulary namespace. Not a resolvable URL today; it is a stable
#: identifier so two exports can be compared and merged without guessing
#: whether two fields named "call" mean the same thing.
VOCAB = "https://github.com/ubeffiong/myconductor/ns/v1#"

SCHEMA = "myconductor.semantic.v1"

#: Terms whose presence in an export would invite exactly the collapse this
#: module exists to prevent. Asserted against every emitted payload.
FORBIDDEN_TERMS = ("resistant", "susceptible", "isResistant", "isSusceptible",
                   "binaryCall", "sir")


def context() -> dict[str, Any]:
    """The ``@context``: what each term means, and what it does not."""
    return {
        "@vocab": VOCAB,
        "myco": VOCAB,
        "id": "@id",
        "type": "@type",
        "drug": {"@id": f"{VOCAB}drug"},
        "call": {"@id": f"{VOCAB}call"},
        "tier": {"@id": f"{VOCAB}tier"},
        "evidence": {"@id": f"{VOCAB}evidence", "@container": "@set"},
        "coverage": {"@id": f"{VOCAB}coverage", "@container": "@set"},
        "discordance": {"@id": f"{VOCAB}discordance"},
        # Stated in the context so a consumer reading only this knows the
        # shape of the vocabulary before parsing a single record.
        "myco:note": (
            "Six call states, not two. A drug may be used only where "
            "establishesUse is true. INDETERMINATE, NOT_ASSESSED, NO_CALL and "
            "UNSUPPORTED are four distinct reasons no verdict was reached and "
            "must not be merged with each other or with SUSCEPTIBLE."),
    }


def call_term(call: Call) -> dict[str, Any]:
    """One call, with its semantics inline rather than implied."""
    return {
        "type": "Call",
        "id": f"{VOCAB}call/{call.value}",
        "value": call.value,
        # The only question most consumers actually want answered, answered
        # explicitly so they never derive it from the enum themselves.
        "establishesUse": call is Call.SUSCEPTIBLE,
        "isEstablished": call.is_established,
        "reasonUnestablished": call.reason_unestablished,
    }


def tier_term(tier: Tier) -> dict[str, Any]:
    return {
        "type": "Tier",
        "id": f"{VOCAB}tier/{tier.value}",
        "value": tier.value,
        "rank": tier.rank,
        "mayEstablishResistance": tier.may_establish_resistance,
    }


def _coverage(result: DrugResult) -> list[dict[str, Any]]:
    return [{
        "type": "LocusCoverage",
        "locus": cov.locus,
        "meanDepth": cov.mean_depth,
        "callableFraction": cov.callable_fraction,
        "source": cov.source,
        # "unknown" and "not callable" are different, and the second is what
        # a None fraction means here.
        "callableEvidencePresent": cov.callable_fraction is not None,
    } for cov in (result.coverage or [])]


def _evidence(result: DrugResult) -> list[dict[str, Any]]:
    out = []
    for item in (result.evidence or []):
        out.append({
            "type": "Evidence",
            "call": call_term(item.call),
            "tier": tier_term(item.tier),
            "lane": getattr(item.lane, "value", str(item.lane)),
            "source": item.source_name,
            "variantKey": item.variant_key,
            "variantLabel": item.variant_label,
            "scope": getattr(item, "scope", None),
            "whoGrade": item.who_grade,
            "rationale": item.rationale,
            "confidence": item.confidence,
            "limitations": list(item.limitations or []),
        })
    return out


def drug_node(result: DrugResult) -> dict[str, Any]:
    node = {
        "type": "DrugResult",
        "drug": result.drug,
        "call": call_term(result.call),
        "tier": tier_term(result.tier),
        "permitsUse": result.permits_use,
        "reason": result.reason,
        "assayStatus": result.assay_status,
        "confidence": result.confidence,
        "evidence": _evidence(result),
        "coverage": _coverage(result),
    }
    if result.discordance is not None:
        node["discordance"] = {
            "type": "Discordance",
            "calls": list(result.discordance.calls),
            "sources": list(result.discordance.sources),
            "note": result.discordance.note,
            "context": list(result.discordance.context),
            "resolved": False,
        }
    return node


def to_jsonld(report: AnalysisReport) -> dict[str, Any]:
    """The report as a linked-data document."""
    document: dict[str, Any] = {
        "@context": context(),
        "type": "AnalysisReport",
        "id": f"urn:myconductor:report:{report.sample_id}",
        "schema": SCHEMA,
        "sampleId": report.sample_id,
        "demoMode": report.demo_mode,
        "drugResults": [drug_node(r) for r in report.drug_results],
        "callStates": [call_term(state) for state in Call],
        "tiers": [tier_term(tier) for tier in Tier],
        "interpretation": {
            "type": "InterpretationGuidance",
            "note": (
                "A drug may be used only where establishesUse is true. The "
                "four unestablished states are distinct and must not be "
                "merged: INDETERMINATE means evidence conflicted or a lane "
                "abstained, NOT_ASSESSED means the required loci were not "
                "shown to be callable, NO_CALL means the genotype is "
                "unreliable, UNSUPPORTED means the drug is outside this "
                "organism profile. None of them means susceptible."),
            "noBinaryField": (
                "This document deliberately carries no resistant/susceptible "
                "boolean. Deriving one discards four states and reintroduces "
                "the presumption that absent evidence implies susceptibility."),
            "terminology": (
                "The vocabulary is local. No SNOMED CT, LOINC or NCIT codes "
                "are asserted, because mapping them is a curation task with "
                "clinical consequences and inventing them here would be the "
                "same failure as inventing terminology codes in a FHIR "
                "bundle."),
        },
    }
    if report.panel:
        document["assayPanel"] = {
            "type": "AssayPanel",
            "name": report.panel.get("name"),
            "version": report.panel.get("version"),
            "assay": report.panel.get("assay"),
            "locusListVerified": report.panel.get("verified", False),
            "unsupportedDrugs": report.panel.get("unsupported_drugs", {}),
            "note": report.panel.get("note", ""),
        }
    if report.provenance is not None:
        provenance = report.provenance
        document["provenance"] = {
            "type": "Provenance",
            "tool": provenance.tool,
            "version": provenance.version,
            "organismProfile": provenance.organism_profile,
            "profileVersion": provenance.profile_version,
            "referenceAssembly": provenance.reference_assembly,
            "depthFloor": provenance.depth_floor,
            "callableFractionFloor": provenance.callable_fraction_floor,
            "coverageSource": provenance.coverage_source,
            "demoMode": provenance.demo_mode,
            "generatedUtc": provenance.generated_utc,
        }
    return document


def audit(document: dict[str, Any]) -> list[str]:
    """Problems that would let a consumer collapse the six states.

    Run over every export. Returns the offending JSON paths, so a regression
    names where it is rather than only that it happened.
    """
    problems: list[str] = []

    def walk(node: Any, path: str) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                # A boolean under a forbidden name is the collapse hazard; the
                # same word as a *string value* (a call's own name) is fine.
                if isinstance(value, bool) and key in FORBIDDEN_TERMS:
                    problems.append(f"{path}.{key}: binary call field")
                walk(value, f"{path}.{key}")
        elif isinstance(node, list):
            for index, item in enumerate(node):
                walk(item, f"{path}[{index}]")

    walk(document, "$")
    if "@context" not in document:
        problems.append("$: no @context, so terms have no stable meaning")
    for drug in document.get("drugResults", []):
        call = drug.get("call") or {}
        if "establishesUse" not in call:
            problems.append(
                f"$.drugResults[{drug.get('drug')}].call: no establishesUse, "
                f"so a consumer must interpret the enum itself")
        if "tier" not in drug:
            problems.append(
                f"$.drugResults[{drug.get('drug')}]: no tier, so the basis "
                f"for the call is lost")
    return problems


def write_jsonld(report: AnalysisReport, path: str | Path,
                 strict: bool = True) -> Path:
    """Write the document, refusing by default to emit a collapsible one."""
    document = to_jsonld(report)
    problems = audit(document)
    if problems and strict:
        raise ValueError(
            "refusing to write a semantic export that a consumer could "
            "collapse to two states: " + "; ".join(problems))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return path


def summarise(document: dict[str, Any]) -> str:
    counts: dict[str, int] = {}
    for drug in document.get("drugResults", []):
        value = (drug.get("call") or {}).get("value", "?")
        counts[value] = counts.get(value, 0) + 1
    usable = sum(1 for drug in document.get("drugResults", [])
                 if (drug.get("call") or {}).get("establishesUse"))
    parts = ", ".join(f"{count} {state}" for state, count in sorted(counts.items()))
    return (f"{document.get('sampleId')}: {parts}; {usable} drug(s) with "
            f"establishesUse true")
