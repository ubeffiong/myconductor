"""Render an AnalysisReport as human-readable text.

The rendering rules exist to stop a reader drawing a conclusion the evidence
does not support:

* every drug appears, including drugs nothing was found for — a drug missing
  from a table is indistinguishable from a drug that was assessed and cleared;
* a non-verdict states *why* it is a non-verdict, on the same line;
* the marker legend separates catalogued, inferred and predicted bases;
* the minority-allele section distinguishes "assessed, nothing found" from
  "could not assess";
* the QC block lists controls that were **not performed**, so an absent check
  is visible rather than merely absent;
* eligibility is phrased as eligibility, never as a recommended regimen.
"""
from __future__ import annotations

from ..core.models import AnalysisReport, Call

_WIDTH = 74
_RULE = "=" * _WIDTH
_THIN = "-" * _WIDTH

_CALL_ORDER = {
    Call.RESISTANT: 0,
    Call.INDETERMINATE: 1,
    Call.NO_CALL: 2,
    Call.NOT_ASSESSED: 3,
    Call.SUSCEPTIBLE: 4,
    Call.UNSUPPORTED: 5,
}


def render_text(report: AnalysisReport) -> str:
    lines: list[str] = []
    p = report.provenance

    lines.append(_RULE)
    lines.append(f" MYCONDUCTOR REPORT  ·  sample: {report.sample_id}")
    lines.append(_RULE)

    if report.demo_mode:
        lines += [
            "",
            " !! SYNTHETIC DEMONSTRATION OUTPUT !!",
            " Values in this report were produced by a synthetic annotator and",
            " carry no biological meaning. Not usable for any purpose beyond",
            " exercising the software interfaces.",
            "",
        ]

    lines.append(f" tool      : {p.tool} v{p.version}")
    lines.append(f" profile   : {p.organism_profile} v{p.profile_version} "
                 f"(assembly {p.reference_assembly})")
    lines.append(f" catalogue : {p.catalogue or 'none'}")
    if p.engines:
        lines.append(f" engines   : {', '.join(str(e) for e in p.engines)}")
    lines.append(f" coverage  : {p.coverage_source} "
                 f"(depth floor {p.depth_floor}x, callable floor "
                 f"{p.callable_fraction_floor:.0%})")
    if p.generated_utc:
        lines.append(f" generated : {p.generated_utc}")

    # -- drug table --------------------------------------------------------
    lines += ["", " DRUG SUSCEPTIBILITY", _THIN]
    lines.append("   legend: R/S catalogued or phenotypic · R*/S* predicted · "
                 "?† inferred")
    lines.append("           ?  indeterminate · -  not assessed · "
                 "x  no call · n/a unsupported")
    lines.append("")
    lines.append("   Only a coverage-backed S counts toward regimen eligibility.")
    lines.append("")

    for r in sorted(report.drug_results,
                    key=lambda x: (_CALL_ORDER.get(x.call, 9), x.drug)):
        conf = f" conf={r.confidence:.2f}" if r.confidence is not None else ""
        lines.append(f"   {r.marker:<4} {r.drug:<14}{conf}")
        if r.call.is_established:
            if r.reason:
                lines.append(f"          {_wrap(r.reason, 10)}")
        else:
            lines.append(f"          why not: "
                         f"{_wrap(r.reason or r.call.reason_unestablished or '', 19)}")
        for ev in r.evidence:
            label = ev.variant_label or "—"
            lines.append(f"          · [{ev.lane.value}/{ev.tier.value}] "
                         f"{label}: {_wrap(ev.rationale, 12)}")
            for limitation in ev.limitations:
                lines.append(f"              ! {_wrap(limitation, 16)}")
        if r.coverage:
            for cov in r.coverage:
                lines.append(f"          ~ {cov.describe()}")

    # -- discordance -------------------------------------------------------
    if report.discordances:
        lines += ["", " DISCORDANT EVIDENCE", _THIN]
        lines.append("   Conflicts are reported, not resolved by vote.")
        lines.append("")
        for d in report.discordances:
            lines.append(f"   ! {d.drug}: {'/'.join(d.calls)} "
                         f"from {', '.join(d.sources)}")
            lines.append(f"       {_wrap(d.note, 7)}")

    # -- minority alleles --------------------------------------------------
    if report.heteroresistance:
        detected = [h for h in report.heteroresistance if h.detected]
        assessed = [h for h in report.heteroresistance
                    if h.assessable and not h.detected]
        blocked = [h for h in report.heteroresistance if not h.assessable]

        lines += ["", " MINORITY-ALLELE ASSESSMENT", _THIN]
        if detected:
            lines.append("   Detected:")
            for h in detected:
                lines.append(f"     ! {h.drug} ({h.variant_label})")
                lines.append(f"         {_wrap(h.note, 9)}")
                if h.posterior_resistance_probability is not None:
                    lo, hi = h.posterior_interval or (None, None)
                    rng = f" [{lo:.0%}-{hi:.0%}]" if lo is not None else ""
                    lines.append(
                        f"         posterior P(true resistance-conferring "
                        f"subpopulation) = {h.posterior_resistance_probability:.0%}"
                        f"{rng}  ({h.calibration_source})")
        if assessed:
            lines.append("   Assessed, no minority population found:")
            for h in assessed:
                lines.append(f"     · {h.drug} ({h.variant_label}): "
                             f"{_wrap(h.note, 7)}")
        if blocked:
            lines.append("   COULD NOT ASSESS (this is not a negative result):")
            for h in blocked:
                lines.append(f"     ? {h.drug} ({h.variant_label}): "
                             f"{_wrap(h.note, 7)}")

    # -- eligibility -------------------------------------------------------
    e = report.eligibility
    lines += ["", " GUIDELINE ELIGIBILITY (not a prescription)", _THIN]
    lines.append(f"   {_wrap(e.summary, 3)}")
    lines.append("")
    for a in e.assessments:
        status = "ELIGIBLE" if a.eligible else "NOT ELIGIBLE"
        lines.append(f"   [{status}] {a.name} "
                     f"({len(a.usable)}/{a.min_required} required)")
        lines.append(f"       {_wrap(a.note, 7)}")
        if a.usable:
            lines.append(f"       established susceptible: {', '.join(a.usable)}")
        if a.resistant:
            lines.append(f"       resistant: {', '.join(a.resistant)}")
        for drug, reason in sorted(a.unestablished.items()):
            lines.append(f"       not established — {drug}: {_wrap(reason, 9)}")
    if e.requires_clinical_review:
        lines += ["",
                  "   Regimen construction requires treatment history, site of",
                  "   disease, age, pregnancy, comorbidity, interactions,",
                  "   toxicity, availability and national policy. None of those",
                  "   are visible to this tool."]

    # -- VUS ---------------------------------------------------------------
    if report.vus_priorities:
        lines += ["", " VARIANTS OF UNKNOWN SIGNIFICANCE", _THIN]
        lines.append("   Ranked for laboratory validation. No resistance call is")
        lines.append("   made for any variant below.")
        lines.append("")
        for v in report.vus_priorities:
            score = f" score={v.score:.2f}" if v.score is not None else ""
            lines.append(f"   > {v.variant_label:<22} [{v.priority}]{score}"
                         f"  context: {v.drug or 'unknown'}")
            lines.append(f"       experiment: {_wrap(v.recommended_experiment, 19)}")
            available = [d for d in v.dimensions if d.available]
            for d in available:
                lines.append(f"       {d.name}: {d.value}  ({d.source})")
            if v.data_gaps:
                lines.append(f"       no data for {len(v.data_gaps)} dimension(s): "
                             f"{_wrap(', '.join(v.data_gaps), 19)}")

    # -- mechanism research queue --------------------------------------------
    if report.mechanism_queue:
        lines += ["", " MECHANISM RESEARCH QUEUE", _THIN]
        lines.append("   INDETERMINATE calls with a named hypothesis and the")
        lines.append("   specific evidence gaps that would resolve them. Not a")
        lines.append("   resistance call — a laboratory work list.")
        lines.append("")
        for m in report.mechanism_queue:
            lines.append(f"   > {m.variant_label:<22} context: {m.drug}")
            lines.append(f"       hypothesis: {_wrap(m.hypothesis, 19)}")
            if m.evidence_gaps:
                lines.append(f"       evidence gaps: "
                             f"{_wrap('; '.join(m.evidence_gaps), 23)}")
            if m.resolving_experiments:
                lines.append(f"       would resolve it: "
                             f"{_wrap('; '.join(m.resolving_experiments), 26)}")

    # -- QC ----------------------------------------------------------------
    if report.qc:
        by_status: dict[str, list] = {}
        for f in report.qc:
            by_status.setdefault(f.status, []).append(f)
        lines += ["", " QUALITY CONTROL", _THIN]
        for status in ("fail", "warn", "pass"):
            for f in by_status.get(status, []):
                lines.append(f"   [{status.upper():<4}] {f.check}: "
                             f"{_wrap(f.detail, 12)}")
        not_performed = by_status.get("not_performed", [])
        if not_performed:
            lines.append("")
            lines.append(f"   NOT PERFORMED ({len(not_performed)} control(s)) — "
                         f"absent, not passed:")
            for f in not_performed:
                lines.append(f"     · {f.check}: {_wrap(f.detail, 7)}")

    if report.qc_warnings:
        lines += ["", " INPUT WARNINGS", _THIN]
        for w in report.qc_warnings:
            lines.append(f"   ~ {_wrap(w, 5)}")

    # -- routing -----------------------------------------------------------
    lines += ["", " LANE ROUTING", _THIN]
    lines.append("   Counts exceed the variant total where lanes overlap; a")
    lines.append("   variant may be examined by several lanes.")
    for lane, n in sorted(report.lane_counts.items()):
        if n:
            lines.append(f"     {lane:<20} {n}")

    lines += ["", _RULE]
    lines.append(" Research and decision-support scaffold. NOT a clinical device.")
    lines.append(" No component has been clinically validated. Every call above")
    lines.append(" requires review by a qualified clinician against phenotypic")
    lines.append(" testing and current treatment guidelines.")
    lines.append(_RULE)
    return "\n".join(lines)


def _wrap(text: str, indent: int, width: int = _WIDTH) -> str:
    """Soft-wrap a paragraph to the report width, continuation lines indented."""
    if not text:
        return ""
    limit = max(24, width - indent)
    words, out, line = text.split(), [], ""
    for word in words:
        if line and len(line) + 1 + len(word) > limit:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        out.append(line)
    pad = " " * indent
    return ("\n" + pad).join(out)
