"""Emit a FHIR-style DiagnosticReport bundle (JSON).

Simplified but shaped like an HL7 FHIR R4 Bundle so downstream LIMS/EHR
integration is a real path, not an afterthought. Each drug becomes an
Observation; provenance is preserved for auditability.
"""
from __future__ import annotations

import json

from ..core.models import AnalysisReport, Call

_INTERP = {
    Call.RESISTANT: "R",
    Call.PREDICTED_RESISTANT: "R",
    Call.SUSCEPTIBLE: "S",
    Call.PREDICTED_SUSCEPTIBLE: "S",
    Call.INDETERMINATE: "IND",
}


def to_fhir_bundle(report: AnalysisReport) -> dict:
    observations = []
    for i, r in enumerate(report.drug_results):
        observations.append(
            {
                "resource": {
                    "resourceType": "Observation",
                    "id": f"obs-{i}",
                    "status": "preliminary" if r.predicted else "final",
                    "category": [{"text": "microbiology-susceptibility"}],
                    "code": {"text": f"{r.drug} susceptibility"},
                    "interpretation": [{"text": _INTERP[r.call]}],
                    "note": [{"text": e.rationale} for e in r.evidence],
                    "extension": [
                        {"url": "confidence", "valueDecimal": r.confidence},
                        {"url": "predicted", "valueBoolean": r.predicted},
                    ],
                }
            }
        )

    p = report.provenance
    bundle = {
        "resourceType": "Bundle",
        "type": "collection",
        "entry": [
            {
                "resource": {
                    "resourceType": "DiagnosticReport",
                    "status": "preliminary",
                    "code": {"text": "M. tuberculosis drug-resistance profile"},
                    "subject": {"reference": f"Sample/{report.sample_id}"},
                    "conclusion": report.regimen.rationale,
                }
            },
            *observations,
            {
                "resource": {
                    "resourceType": "Provenance",
                    "agent": [{"who": {"display": f"{p.tool} v{p.version}"}}],
                    "extension": [
                        {"url": "catalogueVersion", "valueString": p.catalogue_version},
                        {"url": "vusModel", "valueString": p.vus_model},
                        {"url": "coverageMin", "valueInteger": p.coverage_min},
                    ],
                }
            },
        ],
    }
    return bundle


def to_fhir_json(report: AnalysisReport, indent: int = 2) -> str:
    return json.dumps(to_fhir_bundle(report), indent=indent)
