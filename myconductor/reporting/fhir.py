"""FHIR R4 output.

What is and is not conformant
-----------------------------
The structure here is real FHIR R4: a ``Bundle`` of ``Specimen``,
``DiagnosticReport`` with populated ``result`` references, one ``Observation``
per drug, and ``Provenance``. Three codings are used that are genuinely
standard and that this module is confident about:

* ``http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation`` —
  ``R``, ``S``, ``IND``;
* ``http://terminology.hl7.org/CodeSystem/data-absent-reason`` —
  ``not-performed``, ``not-applicable``, ``error``, ``unknown``;
* ``http://terminology.hl7.org/CodeSystem/observation-category`` —
  ``laboratory``.

What is deliberately **not** here: invented LOINC and SNOMED CT codes. A drug
susceptibility observation should carry a LOINC code for that specific
drug/method pair and a SNOMED CT code for the organism, and this module does
not know them. Fabricating plausible-looking codes would be worse than omitting
them, because a receiving system would key on them. So each Observation
carries ``code.text`` plus, optionally, real codes injected through
``TerminologyMap`` by a deployment whose terminology service supplies them.

``conformance_notes()`` lists every gap, and the bundle itself carries that
list — a receiving system is told what is missing rather than left to discover
it. Validate with the official HL7 validator against a declared implementation
guide before any integration.

The modelling choice that matters most
--------------------------------------
``NOT_ASSESSED``, ``UNSUPPORTED`` and ``NO_CALL`` become an Observation with
``dataAbsentReason`` and **no** value — not a value of "susceptible", and not
an omitted Observation. That is the correct FHIR representation of "we did not
establish this", and it is the interoperability form of this whole project's
central rule.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Optional

from ..core.models import AnalysisReport, Call, DrugResult

_EXT = "https://myconductor.org/fhir/StructureDefinition"
_INTERP_SYSTEM = "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation"
_ABSENT_SYSTEM = "http://terminology.hl7.org/CodeSystem/data-absent-reason"
_CATEGORY_SYSTEM = "http://terminology.hl7.org/CodeSystem/observation-category"

#: Call -> (interpretation code, display) for the calls that have a value.
_INTERPRETATION = {
    Call.RESISTANT: ("R", "Resistant"),
    Call.SUSCEPTIBLE: ("S", "Susceptible"),
    Call.INDETERMINATE: ("IND", "Indeterminate"),
}

#: Call -> data-absent-reason code for the calls that have no value.
_ABSENT_REASON = {
    Call.NOT_ASSESSED: ("not-performed", "Not Performed"),
    Call.UNSUPPORTED: ("not-applicable", "Not Applicable"),
    Call.NO_CALL: ("error", "Error"),
}


@dataclass
class TerminologyMap:
    """Injection point for real codes from a deployment's terminology service.

    ``drug_loinc`` maps a drug name to a ``(code, display)`` LOINC pair;
    ``organism_snomed`` is a ``(code, display)`` SNOMED CT pair for the
    organism. Anything not supplied is simply omitted — never guessed.
    """

    drug_loinc: dict[str, tuple[str, str]] = field(default_factory=dict)
    organism_snomed: Optional[tuple[str, str]] = None
    method_loinc: Optional[tuple[str, str]] = None

    @property
    def is_empty(self) -> bool:
        return not (self.drug_loinc or self.organism_snomed or self.method_loinc)


def conformance_notes(report: AnalysisReport,
                      terminology: Optional[TerminologyMap] = None) -> list[str]:
    """Every reason this bundle is not yet ready for clinical exchange."""
    terminology = terminology or TerminologyMap()
    notes = [
        "Not validated against any FHIR implementation guide. Run the HL7 "
        "validator against a declared IG before integration.",
        "No Patient resource: this bundle carries a Specimen and results only. "
        "A deployment must link its own Patient reference.",
    ]
    missing_loinc = sorted({r.drug for r in report.drug_results
                            if r.drug not in terminology.drug_loinc})
    if missing_loinc:
        notes.append(
            f"No LOINC code for {len(missing_loinc)} drug observation(s) "
            f"({', '.join(missing_loinc)}); code.text only. Supply codes via "
            f"TerminologyMap — none are invented here.")
    if terminology.organism_snomed is None:
        notes.append("No SNOMED CT organism code; organism given as text only.")
    if terminology.method_loinc is None:
        notes.append("No method coding: the susceptibility method (genotypic "
                     "prediction vs. phenotypic DST) is described in text only.")
    if report.demo_mode:
        notes.append("SYNTHETIC DEMONSTRATION OUTPUT: values carry no "
                     "biological meaning.")
    if not report.provenance.catalogue:
        notes.append("No catalogue provenance recorded.")
    return notes


def _observation(index: int, r: DrugResult, report: AnalysisReport,
                 terminology: TerminologyMap) -> dict:
    obs: dict = {
        "resourceType": "Observation",
        "id": f"obs-{index}",
        # Genotypic prediction is never 'final'. A catalogued call is as
        # settled as this tool gets, and that is still preliminary.
        "status": "preliminary",
        "category": [{"coding": [{
            "system": _CATEGORY_SYSTEM,
            "code": "laboratory",
            "display": "Laboratory",
        }]}],
        "code": _drug_code(r.drug, terminology),
        "specimen": {"reference": f"Specimen/specimen-{report.sample_id}"},
        "extension": [
            {"url": f"{_EXT}/call", "valueCode": r.call.value},
            {"url": f"{_EXT}/evidence-tier", "valueCode": r.tier.value},
            {"url": f"{_EXT}/permits-regimen-use",
             "valueBoolean": r.permits_use},
        ],
    }

    if r.confidence is not None:
        obs["extension"].append(
            {"url": f"{_EXT}/confidence", "valueDecimal": r.confidence})

    interpretation = _INTERPRETATION.get(r.call)
    if interpretation:
        code, display = interpretation
        obs["valueCodeableConcept"] = {"coding": [{
            "system": _INTERP_SYSTEM, "code": code, "display": display,
        }], "text": r.call.value}
        obs["interpretation"] = [{"coding": [{
            "system": _INTERP_SYSTEM, "code": code, "display": display,
        }]}]
    else:
        code, display = _ABSENT_REASON.get(r.call, ("unknown", "Unknown"))
        obs["dataAbsentReason"] = {"coding": [{
            "system": _ABSENT_SYSTEM, "code": code, "display": display,
        }], "text": r.reason or r.call.reason_unestablished or r.call.value}

    notes = []
    if r.reason:
        notes.append({"text": r.reason})
    for ev in r.evidence:
        detail = f"[{ev.lane.value}/{ev.tier.value}] {ev.rationale}"
        if ev.limitations:
            detail += " Limitations: " + "; ".join(ev.limitations)
        notes.append({"text": detail})
        if ev.variant_label:
            obs["extension"].append({
                "url": f"{_EXT}/determinant",
                "extension": [
                    {"url": "label", "valueString": ev.variant_label},
                    {"url": "identity", "valueString": ev.variant_key or ""},
                    {"url": "source", "valueString": str(ev.engine or ev.lane.value)},
                ] + ([{"url": "whoGrade", "valueString": ev.who_grade}]
                     if ev.who_grade else []),
            })
    for cov in r.coverage:
        notes.append({"text": f"coverage — {cov.describe()}"})
    if r.discordance:
        notes.append({"text": f"DISCORDANT: {r.discordance.note}"})
    if notes:
        obs["note"] = notes

    return obs


def _drug_code(drug: str, terminology: TerminologyMap) -> dict:
    code = {"text": f"{drug} susceptibility"}
    loinc = terminology.drug_loinc.get(drug)
    if loinc:
        code["coding"] = [{
            "system": "http://loinc.org", "code": loinc[0], "display": loinc[1],
        }]
    return code


def to_fhir_bundle(report: AnalysisReport,
                   terminology: Optional[TerminologyMap] = None) -> dict:
    terminology = terminology or TerminologyMap()
    specimen_id = f"specimen-{report.sample_id}"

    observations = [
        _observation(i, r, report, terminology)
        for i, r in enumerate(report.drug_results)
    ]

    organism = {"text": "Mycobacterium tuberculosis complex"}
    if terminology.organism_snomed:
        organism["coding"] = [{
            "system": "http://snomed.info/sct",
            "code": terminology.organism_snomed[0],
            "display": terminology.organism_snomed[1],
        }]

    p = report.provenance
    diagnostic_report = {
        "resourceType": "DiagnosticReport",
        "id": "report-1",
        "status": "preliminary",
        "category": [{"coding": [{
            "system": _CATEGORY_SYSTEM, "code": "laboratory",
            "display": "Laboratory",
        }]}],
        "code": {"text": f"{organism['text']} drug-resistance profile "
                         f"(genotypic prediction)"},
        "specimen": [{"reference": f"Specimen/{specimen_id}"}],
        "result": [{"reference": f"Observation/{o['id']}"} for o in observations],
        "conclusion": report.eligibility.summary,
        "extension": [
            {"url": f"{_EXT}/organism", "valueCodeableConcept": organism},
            {"url": f"{_EXT}/requires-clinical-review", "valueBoolean": True},
            {"url": f"{_EXT}/synthetic-demonstration",
             "valueBoolean": report.demo_mode},
        ],
        "note": [{"text": n} for n in conformance_notes(report, terminology)],
    }

    entries = [
        {"fullUrl": f"urn:uuid:{specimen_id}", "resource": {
            "resourceType": "Specimen",
            "id": specimen_id,
            "identifier": [{"value": report.sample_id}],
            "note": [{"text": f"reference assembly {p.reference_assembly}"}],
        }},
        {"fullUrl": "urn:uuid:report-1", "resource": diagnostic_report},
        *[{"fullUrl": f"urn:uuid:{o['id']}", "resource": o} for o in observations],
        {"fullUrl": "urn:uuid:provenance-1", "resource": {
            "resourceType": "Provenance",
            "id": "provenance-1",
            "target": [{"reference": "DiagnosticReport/report-1"}],
            "recorded": p.generated_utc,
            "agent": [{
                "type": {"text": "assembler"},
                "who": {"display": f"{p.tool} v{p.version}"},
            }],
            "extension": [
                {"url": f"{_EXT}/organismProfile",
                 "valueString": f"{p.organism_profile} v{p.profile_version}"},
                {"url": f"{_EXT}/catalogue", "valueString": str(p.catalogue or "none")},
                {"url": f"{_EXT}/engines",
                 "valueString": ", ".join(str(e) for e in p.engines) or "none"},
                {"url": f"{_EXT}/coverageSource", "valueString": p.coverage_source},
                {"url": f"{_EXT}/depthFloor", "valueInteger": p.depth_floor},
                {"url": f"{_EXT}/callableFractionFloor",
                 "valueDecimal": p.callable_fraction_floor},
            ],
        }},
    ]

    # Unperformed QC controls travel with the bundle: a receiving system should
    # be able to see which checks were never run.
    not_performed = [f for f in report.qc if f.status == "not_performed"]
    failed = [f for f in report.qc if f.status == "fail"]
    if not_performed or failed:
        entries.append({"fullUrl": "urn:uuid:qc-1", "resource": {
            "resourceType": "Observation",
            "id": "qc-1",
            "status": "preliminary",
            "category": [{"coding": [{
                "system": _CATEGORY_SYSTEM, "code": "laboratory",
                "display": "Laboratory",
            }]}],
            "code": {"text": "Analytical quality control summary"},
            "specimen": {"reference": f"Specimen/{specimen_id}"},
            "note": (
                [{"text": f"FAILED — {f.check}: {f.detail}"} for f in failed]
                + [{"text": f"NOT PERFORMED — {f.check}: {f.detail}"}
                   for f in not_performed]
            ),
        }})
        diagnostic_report["result"].append({"reference": "Observation/qc-1"})

    return {
        "resourceType": "Bundle",
        "type": "collection",
        "timestamp": p.generated_utc,
        "entry": entries,
    }


def to_fhir_json(report: AnalysisReport, indent: int = 2,
                 terminology: Optional[TerminologyMap] = None) -> str:
    return json.dumps(to_fhir_bundle(report, terminology), indent=indent)
