"""Render an AnalysisReport as human-readable text.

Predicted calls are visually marked (*) and separated from catalogued calls, so
a reader never mistakes a model prediction for a graded catalogue result.
"""
from __future__ import annotations

from ..core.models import AnalysisReport, Call

_MARK = {
    Call.RESISTANT: "R ",
    Call.SUSCEPTIBLE: "S ",
    Call.PREDICTED_RESISTANT: "R*",
    Call.PREDICTED_SUSCEPTIBLE: "S*",
    Call.INDETERMINATE: "? ",
}


def render_text(report: AnalysisReport) -> str:
    lines: list[str] = []
    p = report.provenance
    lines.append("=" * 64)
    lines.append(f" MYCONDUCTOR REPORT  ·  sample: {report.sample_id}")
    lines.append("=" * 64)
    lines.append(f" tool {p.tool} v{p.version} | catalogue {p.catalogue_version}")
    lines.append(f" VUS model: {p.vus_model} | depth floor: {p.coverage_min}x")
    lines.append("")

    lines.append("-- Drug susceptibility --------------------------------")
    lines.append("   (R/S = catalogued, R*/S* = predicted, ? = indeterminate)")
    for r in sorted(report.drug_results, key=lambda x: x.drug):
        route = r.evidence[0].route.value if r.evidence else "-"
        lines.append(
            f"   {_MARK[r.call]} {r.drug:<14} conf={r.confidence:<5.2f} via {route}"
        )
        for ev in r.evidence:
            lines.append(f"        - {ev.variant_key}: {ev.rationale}")

    if report.heteroresistance:
        lines.append("")
        lines.append("-- Heteroresistance warnings --------------------------")
        for h in report.heteroresistance:
            lines.append(f"   ! {h.drug} ({h.variant_key}): {h.note}")

    lines.append("")
    lines.append("-- Proposed regimen -----------------------------------")
    reg = report.regimen
    status = "ADEQUATE" if reg.adequate else "INADEQUATE"
    lines.append(f"   status: {status}")
    lines.append(f"   {reg.rationale}")
    if reg.proposed:
        lines.append(f"   drugs: {', '.join(reg.proposed)}")
    if reg.ineffective_drugs:
        lines.append(f"   dropped: {', '.join(reg.ineffective_drugs)}")

    if report.discovery.triggered:
        lines.append("")
        lines.append("-- Discovery loop (novel-target proposals) ------------")
        lines.append("   NOTE: hypothesis-generating only; NOT validated therapy.")
        for t in report.discovery.targets:
            if t.best_docking_kcal_mol is None:
                continue
            lines.append(
                f"   > {t.gene:<8} {t.protein[:34]:<34} "
                f"dock {t.best_docking_kcal_mol:>6} kcal/mol | {t.top_compound}"
            )

    if report.qc_warnings:
        lines.append("")
        lines.append("-- QC warnings ----------------------------------------")
        for w in report.qc_warnings:
            lines.append(f"   ~ {w}")

    lines.append("")
    lines.append("-- Routing ---------------------------------------------")
    lines.append(f"   {report.routed_counts}")
    lines.append("=" * 64)
    lines.append(" Research/decision-support scaffold. Not for clinical use.")
    lines.append("=" * 64)
    return "\n".join(lines)
