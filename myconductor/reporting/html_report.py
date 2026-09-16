"""Offline review report: escaped content, no external assets or scripts.

One self-contained file
-----------------------
The stylesheet and the runtime are package data, inlined at render time, so the
output opens from a USB stick on a machine that has never had a network. That
is not a stylistic preference: this tool is aimed at settings where a report
that silently fails to load a font or a charting library is a report that
cannot be read. There is no CDN, no webfont, no ``<script src>``, and a test
asserts it stays that way.

Nothing is generated that was not measured
------------------------------------------
The design this is ported from carried a seeded PRNG that produced call
matrices, coverage grids, lineage strata and sequence bases, plus seventeen
hard-coded datasets. All of it is gone. Every panel reads what the run
produced, and a panel with nothing to show renders an empty state saying so.
A report that fabricates plausible numbers is worse than one with visible
gaps: it is indistinguishable from a real result to everyone downstream.

The previous report's content is preserved whole
------------------------------------------------
The per-drug review table with its evidence chains, the investigation and
follow-up queue, the quality-control list and the provenance block are all
still here, as sections R1-R3 and in the methods grid. The new sections are
additions, not replacements.
"""
from __future__ import annotations

import json
import re
from html import escape
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

from . import html_data as shape
from .labels import EXPLICIT, humanise, variant_name
from .json_report import report_dict

ASSETS = Path(__file__).resolve().parent / "assets"

#: Status classes for the QC list.
_QC_CLASS = {"pass": "qc-pass", "warn": "qc-warn", "fail": "qc-fail"}


def _asset(name: str) -> str:
    return (ASSETS / name).read_text(encoding="utf-8")


def _esc(value: Any) -> str:
    return escape(str(value if value is not None else ""))


# -- preserved sections ---------------------------------------------------
def _review_rows(data: dict) -> str:
    rows = []
    for result in data.get("drug_results", []):
        cells = [
            result.get("drug", "").capitalize(),
            humanise(result.get("call", "")),
            humanise(result.get("genomic_call")) or "Unavailable",
            humanise(result.get("phenotypic_call")) or "Not measured",
            humanise(result.get("assay_status", "")),
            _prose(result.get("reason") or ""),
        ]
        evidence = "".join(
            "<li>" + _esc(_prose(item.get("rationale", ""))) + " <small>"
            + _esc((item.get("evidence_id") or "")[:16]) + "</small></li>"
            for item in result.get("evidence", []))
        call_class = shape.CALL_CLASS.get(result.get("call"), "na")
        body = "".join(f"<td>{_esc(c)}</td>" for c in cells[2:])
        rows.append(
            f'<tr><td style="font-family:var(--mono)">{_esc(cells[0])}</td>'
            f'<td><span class="pill call-{call_class}">{_esc(cells[1])}</span></td>'
            + body + "</tr>"
            + "<tr><td colspan='6'><details><summary>Evidence chain</summary><ul>"
            + (evidence or "<li>No evidence recorded for this drug.</li>")
            + "</ul></details></td></tr>")
    return "".join(rows) or (
        "<tr><td colspan='6' style='padding:22px;color:var(--text-2)'>"
        "No drug results in this run.</td></tr>")


def _followup(data: dict) -> str:
    follow = {f["investigation_id"]: f for f in data.get("follow_up", [])}
    blocks = []
    for finding in data.get("investigations", []):
        action = follow.get(finding.get("id"), {})
        chosen = action.get("selected") or {}
        status = chosen.get("action", action.get("status", "unknown"))
        blocks.append(
            '<div class="interp-block">'
            f'<h3>{_esc(finding.get("drug", "").capitalize())}: '
            f'{_esc(humanise(finding.get("category", "")))}</h3>'
            f'<p>{_esc(humanise(finding.get("certainty", "")))} — '
            f'{_esc(_prose(finding.get("explanation", "")))}</p>'
            f'<p><strong>Follow-up:</strong> {_esc(humanise(status))}</p>'
            f'<p>Priority {_esc(finding.get("priority", ""))}; estimated cost '
            f'{_esc(chosen.get("cost", "unknown"))} '
            f'{_esc(chosen.get("currency", ""))}</p>'
            "</div>")
    return "".join(blocks) or (
        '<div class="empty-state"><div class="empty-title">No unresolved findings'
        '</div><div class="empty-detail">Nothing in this run was queued for '
        'investigation.</div></div>')


def _qc(data: dict) -> str:
    items = []
    for check in data.get("qc", []):
        status = (check.get("status") or "").lower()
        cls = _QC_CLASS.get(status, "qc-none")
        items.append(
            f'<li><span class="qc-status {cls}">'
            f'{_esc(humanise(status or "not_performed"))}</span>'
            f'<span><strong>{_esc(humanise(check.get("check", "")))}</strong> — '
            f'{_esc(_prose(check.get("detail", "")))}</span></li>')
    return "".join(items) or (
        '<li><span class="qc-status qc-none">none</span>'
        '<span>No quality-control findings were recorded. That is not the same '
        'as every control passing.</span></li>')


def _prose(text) -> str:
    """Replace known internal identifiers inside upstream prose.

    Only tokens present in the label map are touched, so this cannot reword a
    clinical explanation — it just stops ``callable_mask`` and
    ``availability_unknown`` appearing mid-sentence as though they were words.
    """
    if not text:
        return ""
    out = str(text)
    for token in sorted(EXPLICIT, key=len, reverse=True):
        if "_" in token and token in out:
            out = out.replace(token, EXPLICIT[token].lower())
    out = re.sub(r"\b[A-Za-z][A-Za-z0-9]*_[A-Za-z0-9cC.>\-]+",
                 lambda match: variant_name(match.group(0)), out)
    return out

# -- generated commentary -------------------------------------------------
def _interpretation(data: dict, drugs: Sequence[dict],
                    lineages: Sequence[str],
                    unmeasurable: Sequence[dict]) -> str:
    """Commentary derived from this run, never a fixed narrative.

    The mock's interpretation section was written prose about one illustrative
    cohort. Shipping that text beside different numbers would make the report
    argue for findings it does not contain, so each block here is emitted only
    when the run supports it, and quotes the figure it is about.
    """
    blocks: list[str] = []

    thin = [d for d in drugs if d["estimable"] and d["coverage"] is not None
            and d["coverage"] < 60]
    if thin:
        named = ", ".join(f"{d['name'].lower()} at {d['coverage']:.0f}%"
                          for d in thin[:3])
        blocks.append(_interp(
            "I.1", "Coverage is the first-order finding",
            f"The most important numbers here are not error rates — they are "
            f"{_esc(named)}. On this cohort those drugs cannot be scored for "
            f"most isolates, and every error rate reported for them is "
            f"computed on the selected subset that happened to have the "
            f"relevant loci callable.",
            "bad",
            "Consequence: those error rates are conditional on the locus being "
            "sequenced. Reporting them without the coverage figure beside them "
            "would overstate operational readiness."))

    vme = [(d["name"], d["vme"]) for d in drugs
           if d["estimable"] and (d["vme"] or 0) > 0]
    if vme:
        vme.sort(key=lambda item: -item[1])
        named = ", ".join(f"{n.lower()} ({v})" for n, v in vme[:4])
        blocks.append(_interp(
            "I.2", "The false-susceptible asymmetry",
            f"Very major errors — a drug wrongly cleared for use — are not "
            f"evenly distributed. They concentrate in {_esc(named)}. A false "
            f"susceptible call puts a patient on a failing regimen, which is "
            f"how amplified resistance and onward transmission are made; a "
            f"false resistant call costs a usable drug. The two are not "
            f"interchangeable and this report never averages them into one "
            f"figure without also showing them apart.",
            "",
            "Hypothesis to test: whether these cluster in one lineage or in "
            "isolates carrying non-canonical variants. The lineage section "
            "addresses the first; the alignment view allows direct inspection."))

    if not lineages:
        blocks.append(_interp(
            "I.3", "Lineage was not assessed",
            "No isolate in this run carried a recorded lineage, so the "
            "lineage-stratified view is empty. That is reported as unknown "
            "rather than as a single lineage: folding untyped isolates into "
            "one bucket produces a diversity of exactly one, which is "
            "indistinguishable from a genuinely lineage-restricted cohort.",
            "",
            "An aggregate error rate can be carried entirely by one lineage. "
            "Until lineages are supplied, this report cannot say whether it is."))

    not_estimable = [d["name"] for d in drugs if not d["estimable"]]
    if not_estimable:
        blocks.append(_interp(
            "I.4", "Not-estimable is a finding, not a gap",
            f"{len(not_estimable)} drug(s) — {_esc(', '.join(not_estimable))} — "
            f"cannot be scored on this cohort. That is a measurement of the "
            f"field rather than a defect of the benchmark: there are too few "
            f"phenotypically characterised resistant isolates for these drugs "
            f"to support an error-rate estimate. Any tool claiming accuracy "
            f"for them from a cohort like this is overstating its evidence.",
            "",
            "For these drugs the report shows INDETERMINATE with a coverage "
            "note, never a susceptibility call. That is enforced structurally, "
            "not by convention."))

    if unmeasurable:
        both = [r["drug"] for r in unmeasurable if r.get("missing_side") == "both"]
        detail = "; ".join(
            f"{_esc(r['drug'])} (missing {_esc(r.get('missing_side', ''))})"
            for r in unmeasurable[:6])
        extra = ""
        if both:
            extra = (f" {_esc(', '.join(both))} "
                     f"{'is' if len(both) == 1 else 'are'} absent from both the "
                     f"catalogue and the phenotype source, so no pairing of "
                     f"these two can evaluate "
                     f"{'it' if len(both) == 1 else 'them'} at all.")
        blocks.append(_interp(
            "I.5", "Drugs this pairing cannot evaluate",
            f"Evaluating a call needs a genotypic side to make it from and a "
            f"phenotypic side to score it against. These lack one or both: "
            f"{detail}.{extra} Listing only the drugs that happen to be "
            f"covered would read as a complete panel.",
            "",
            "The two absences have different remedies: no catalogue entry means "
            "nothing can be predicted; no phenotype means a prediction cannot "
            "be checked."))

    blocks.append(_interp(
        "I.6", "What this report does not claim",
        "No component here has been clinically evaluated, and no sensitivity, "
        "specificity or error rate in this document should be read as a "
        "validated performance characteristic of a diagnostic. Where a "
        "benchmark figure is shown it is measured against the reference "
        "phenotype named in the methods grid, on the cohort described there, "
        "and nowhere else.",
        "",
        "Every panel with no data renders an empty state rather than a "
        "generated one. A blank section in this report means the run produced "
        "nothing for it."))
    return "\n".join(blocks)


def _interp(number: str, title: str, body: str, tone: str, callout: str) -> str:
    tone_class = f" {tone}" if tone else ""
    return (
        f'<div class="interp-block"><h3><span class="num">{number}</span>'
        f'{title}</h3><p>{body}</p>'
        f'<div class="callout{tone_class}">{callout}</div></div>')


def _methods(data: dict, drugs: Sequence[dict], n_samples: int,
             reference: Optional[str]) -> str:
    provenance = data.get("provenance", {})
    items = [
        ("Cohort", f"{n_samples} sample(s) analysed in this run. "
                   f"{len(drugs)} drug(s) reported."),
        ("Reference phenotype",
         _esc(reference) if reference else
         "No phenotypic reference was supplied, so no accuracy figure is "
         "computed. The benchmark section is empty rather than estimated."),
        ("Call ontology",
         "Six states: resistant, susceptible, indeterminate, not assessed, "
         "no call, unsupported. Only catalogued and phenotypic evidence may "
         "establish RESISTANT; coverage is a precondition for SUSCEPTIBLE."),
        ("Coverage",
         f"Callable-locus evidence from outside the variant list. Depth floor "
         f"{_esc(provenance.get('depth_floor'))}x, callable-fraction floor "
         f"{_esc(provenance.get('callable_fraction_floor'))}. Source: "
         f"{_esc(provenance.get('coverage_source', 'absent'))}."),
        ("Organism profile",
         f"{_esc(provenance.get('organism_profile'))} "
         f"{_esc(provenance.get('profile_version'))}, against "
         f"{_esc(provenance.get('reference_assembly'))}."),
        ("Tool version",
         f"{_esc(provenance.get('tool'))} {_esc(provenance.get('version'))}. "
         f"Generated {_esc(provenance.get('generated_utc') or 'unknown')}."),
        ("Analysis fingerprint",
         "<code>" + _esc(provenance.get("analysis_fingerprint") or "unavailable")
         + "</code>. Policy: " + _esc(provenance.get("policy_version", "unknown"))),
        ("What this is not",
         "A research and decision-support scaffold, not a clinical device. "
         "Every call requires review by a qualified clinician against "
         "phenotypic testing and current treatment guidelines."),
    ]
    if data.get("panel"):
        panel = data["panel"]
        verified = "verified" if panel.get("verified") else "NOT verified"
        items.insert(1, (
            "Assay panel",
            f"{_esc(panel.get('name'))} {_esc(panel.get('version'))} "
            f"({_esc(panel.get('assay'))}), {panel.get('n_loci', 0)} locus/loci; "
            f"locus list {verified} against the manufacturer's design. A locus "
            f"the panel does not target is never callable, whatever a coverage "
            f"file claims."))
    return "".join(
        f'<div class="method-item"><h4>{title}</h4><p>{body}</p></div>'
        for title, body in items)



#: The export bar. Restoring this mattered more than it looks: the runtime
#: always carried ``downloadReport``, but the bar was lost when the design's
#: header was replaced, so nothing could reach it and a reader had no way to
#: get the data back out of the page.
DOWNLOADS = (
    ("json", "Full report (JSON)", True),
    ("calls", "Call matrix (CSV)", False),
    ("bench", "Benchmark (CSV)", False),
    ("disc", "Discordance (TSV)", False),
    ("vcf", "Variants (VCF)", False),
    ("coverage", "Coverage (TSV)", False),
    ("audit", "Audit ledger (JSONL)", False),
    ("html", "This report (HTML)", False),
)

_DL_ICON = (
    '<svg class="ico" viewBox="0 0 16 16" aria-hidden="true">'
    '<path d="M8 1v9m0 0L4.5 6.5M8 10l3.5-3.5M2 12v1.5A1.5 1.5 0 003.5 15h9a1.5'
    ' 1.5 0 001.5-1.5V12" fill="none" stroke="currentColor" stroke-width="1.5"'
    ' stroke-linecap="round" stroke-linejoin="round"/></svg>')


def _download_bar() -> str:
    buttons = []
    for kind, label, primary in DOWNLOADS:
        cls = "dl-btn primary" if primary else "dl-btn"
        onclick = "downloadReport('" + kind + "')"
        buttons.append('<button class="' + cls + '" onclick="' + onclick
                       + '">' + _DL_ICON + _esc(label) + "</button>")
    return ('<div class="download-bar" id="downloadBar">'
            + "".join(buttons) + "</div>")

# -- summary figures ------------------------------------------------------
def _kpi(label: str, value: Any, sub: str, tone: str = "") -> str:
    tone_class = f" {tone}" if tone else ""
    value_class = f" {tone}" if tone else " accent"
    return (f'<div class="kpi{tone_class}"><div class="label">{_esc(label)}</div>'
            f'<div class="value{value_class}">{_esc(value)}</div>'
            f'<div class="sub">{_esc(sub)}</div></div>')


def _summary(payload: dict, reports: Sequence[dict],
             drugs: Sequence[dict]) -> dict:
    """The executive strip, computed. Never a figure the run did not produce.

    The design this is ported from hard-coded every one of these — 81.5%
    callable, 17 very major errors, 24 discordances — and a single-sample run
    rendered them unchanged beside a header that correctly said one isolate.
    Two numbers on one screen, disagreeing, both looking authoritative. Every
    tile here is derived, and a tile with nothing behind it says so instead of
    showing a zero, because zero is a measurement.
    """
    n_samples = len(reports)
    estimable = [d for d in drugs if d["estimable"]]
    coverages = [d["coverage"] for d in estimable if d["coverage"] is not None]
    vme = sum(d["vme"] or 0 for d in estimable)
    discordances = len(payload["discordance"])
    vus = len(payload["vus"])

    calls = [r for report in reports for r in report.get("drug_results", [])]
    usable = sum(1 for c in calls if c.get("call") == "susceptible")
    unassessed = sum(1 for c in calls if c.get("call") == "not_assessed")

    tiles = [
        _kpi("Drugs reported", len(drugs),
             f"{len(estimable)} with a measured benchmark"),
        _kpi("Coverage-backed susceptible", usable,
             f"of {len(calls)} call(s) — the only state that admits a drug"),
        _kpi("Not assessed", unassessed,
             "loci not shown callable; not susceptibility",
             "warn" if unassessed else ""),
    ]
    if coverages:
        coverages = sorted(coverages)
        median = coverages[len(coverages) // 2]
        tiles.append(_kpi(
            f"Median coverage ({len(coverages)} drug(s))", f"{median:.0f}%",
            f"range {min(coverages):.0f}% – {max(coverages):.0f}%",
            "warn" if median < 85 else ""))
    if estimable:
        tiles.append(_kpi(
            "Very major errors", vme,
            f"false-susceptible across {len(estimable)} drug(s)",
            "bad" if vme else ""))
    not_estimable = [d["id"] for d in drugs if not d["estimable"]]
    if not_estimable:
        tiles.append(_kpi(
            "Not estimable", len(not_estimable),
            " · ".join(not_estimable[:6]), "warn"))
    tiles.append(_kpi("Discordances", discordances,
                      "retained, never resolved by vote"))
    tiles.append(_kpi("VUS queued", vus, "ranked for laboratory validation"))

    if estimable:
        desc = (f"Measured on {n_samples} isolate(s) against the phenotypic "
                f"reference named in the methods grid, with coverage-gated "
                f"susceptibility enforced throughout. "
                f"{len(estimable)} of {len(drugs)} drug(s) carry a benchmark; "
                f"the rest are reported as not estimable rather than scored.")
    else:
        desc = (f"This run analysed {n_samples} sample(s). No phenotypic "
                f"reference was supplied, so no accuracy figure is computed "
                f"anywhere in this report — the benchmark panels render empty "
                f"rather than estimated. What follows is the evidence the run "
                f"produced and the reasoning applied to it.")

    n_calls = len(calls)
    title = ("Per-drug performance against the reference" if estimable
             else "What this run established, and what it did not")
    return {
        "__SUMMARY_TITLE__": title,
        "__KPIS__": "".join(tiles),
        "__SUMMARY_DESC__": desc,
        "__CALLS_NOTE__": (f"n = {n_calls} call(s) ({n_samples} sample(s) "
                           f"× {len(drugs)} drug(s))"),
        "__TREND_NOTE__": ("ordered error series from the cohort dataset"
                           if payload["error_trend"] else "no series in this run"),
        "__TREND_PROSE__": (
            "An error trend needs an ordered cohort series from the "
            "report-context dataset. This run did not supply one, so the "
            "panel is empty rather than smoothed."
            if not payload["error_trend"] else
            "Each point is one ordered cohort-window value supplied to the "
            "renderer. Interpret movement using the source dataset's window "
            "size, ordering, and denominator; none is inferred from the values."),
        "__SCATTER_PROSE__": (
            "Each point pairs a drug's coverage with its error rate. A drug in "
            "the low-coverage corner is not necessarily a tool failure: its "
            "error rate is conditional on the loci having been sequenced at "
            "all, and the coverage figure is what says how often that held."
            if estimable else
            "No drug in this run carries both a coverage and an error figure, "
            "so there is nothing to plot."),
        "__HEADLINE_NOTE__": (
            f"n = {n_samples} sample(s) · phenotypic reference"
            if estimable else "no phenotypic reference supplied"),
        "__MATRIX_NOTE__": (
            f"{len(payload['isolates'])} sample(s) × {len(drugs)} drug(s)"
            if payload["isolates"] else "no call matrix in this run"),
        "__VUS_NOTE__": (f"{vus} variant(s) ranked for validation" if vus
                         else "nothing queued"),
        "__VALIDATION_NOTE__": (
            f"n = {sum(o.get('v', 0) for o in payload['validation_outcomes'])} "
            f"laboratory validation(s)"
            if payload["validation_outcomes"] else "no validations supplied"),
        "__TIMELINE_NOTE__": ("site-local validation activity"
                              if payload["validation_timeline"]
                              else "no activity recorded"),
        "__FEDERATED_NOTE__": ("aggregated statistics only; isolate-level data "
                               "never leaves a site"),
        "__FEDERATED_KPIS__": _federated_kpis(payload["federated_sites"]),
    }


def _federated_kpis(sites: Sequence[dict]) -> str:
    if not sites:
        return ('<div class="empty-state" style="grid-column:1/-1">'
                '<div class="empty-title">No federated contributions</div>'
                '<div class="empty-detail">No site submitted aggregated '
                'statistics for this run.</div></div>')
    total = sum(_int(s.get("n")) for s in sites)
    ks = [_int(s.get("k")) for s in sites if s.get("k") is not None]
    return (_kpi("Participating sites", len(sites), "aggregated submissions")
            + _kpi("Total submissions", total, "signed and replay-protected")
            + _kpi("Minimum k", min(ks) if ks else "—",
                   "k-anonymity floor across sites"))


def _int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


# -- header ---------------------------------------------------------------
def _meta(data: dict, drugs: Sequence[dict], n_samples: int,
          lineages: Sequence[str]) -> str:
    provenance = data.get("provenance", {})
    estimable = sum(1 for d in drugs if d["estimable"])
    catalogue = provenance.get("catalogue") or {}
    primary = data.get("sample_id", "—")
    pairs = []
    if n_samples > 1:
        pairs.extend([
            ("Report scope", f"Cohort report · {n_samples} isolates"),
            ("Primary isolate", f"{primary} (first record used for interoperability exports)"),
        ])
    else:
        pairs.append(("Report scope", f"Single-isolate report · {primary}"))
    pairs.extend([
        ("Drugs reported",
         f"{estimable} / {len(drugs)} estimable" if drugs else "—"),
        ("Catalogue", catalogue.get("version") or catalogue.get("name") or "bundled"),
        ("Engines", ", ".join(e.get("name", "") for e in provenance.get("engines", []))
         or "none supplied"),
        ("Lineages", " · ".join(lineages) if lineages else "not recorded"),
    ])
    return "".join(
        f'<div class="meta-item"><div class="k">{_esc(k)}</div>'
        f'<div class="v">{_esc(v)}</div></div>' for k, v in pairs)


# -- entry point ----------------------------------------------------------
def render_html(report,
                baseline_rows: Optional[Iterable[dict]] = None,
                lineage_rows: Optional[Iterable[dict]] = None,
                measurability: Optional[Iterable[dict]] = None,
                audit_entries: Optional[Iterable[dict]] = None,
                validation: Optional[dict] = None,
                federated_sites: Optional[Iterable[dict]] = None,
                error_trend: Optional[Sequence[float]] = None,
                prevalence_rows: Optional[Iterable[dict]] = None,
                target_rows: Optional[Iterable[dict]] = None,
                watchlist_rows: Optional[Iterable[dict]] = None,
                external_benchmark_rows: Optional[Iterable[dict]] = None,
                extra_reports: Sequence[Any] = (),
                reference_method: Optional[str] = None) -> str:
    """Render one self-contained HTML report.

    Every argument after ``report`` is optional cohort context. Omitting one
    renders its section's empty state; none of them is ever substituted with a
    generated stand-in.
    """
    data = report_dict(report)
    others = [report_dict(r) for r in extra_reports]
    reports = [data] + others

    drugs = shape.drugs_from_baseline(baseline_rows or [])
    if not drugs:
        # No benchmark: build cards from the run's own drugs so the call
        # matrix and evidence views still have an axis, with every accuracy
        # field left None so nothing reads as measured.
        seen: list[str] = []
        for item in reports:
            for result in item.get("drug_results", []):
                name = result.get("drug", "")
                if name and name not in seen:
                    seen.append(name)
        drugs = [{
            "id": shape.code_for(name), "name": name.capitalize(),
            "cls": shape.DRUG_CLASS.get(name.lower(), "Other"),
            "coverage": None, "error": None, "vme": None, "me": None,
            "n": len(reports), "evaluable": 0, "estimable": False,
            "sensitivity": None, "specificity": None, "sensitivity_ci": "",
            "specificity_ci": "", "ppv": None, "npv": None, "prevalence": None,
            "sens_basis": 0, "spec_basis": 0,
            "notes": ["no phenotypic reference supplied in this run"],
        } for name in sorted(seen)]

    lineages, strata = shape.lineage_strata(lineage_rows or [])
    provenance = data.get("provenance", {})
    alignment = shape.alignment_loci(reports)
    vus = shape.vus_items(reports)
    workflows = shape.implemented_workflows(reports)
    payload = {
        # The run's real identity, so an exported file names the analysis it
        # came from. The mock shipped a fixed string here.
        "run": {
            "sample_id": data.get("sample_id"),
            "scope": "cohort" if len(reports) > 1 else "single_isolate",
            "label": (f"cohort-{len(reports)}-isolates" if len(reports) > 1
                      else data.get("sample_id")),
            "primary_sample_id": data.get("sample_id"),
            "fingerprint": provenance.get("analysis_fingerprint"),
            "generated_utc": provenance.get("generated_utc"),
            "tool": provenance.get("tool"),
            "version": provenance.get("version"),
            "organism_profile": provenance.get("organism_profile"),
            "demo_mode": bool(data.get("demo_mode")),
        },
        "drugs": drugs,
        "lineages": lineages,
        "lineage": strata,
        "coverage_loci": shape.coverage_loci(reports),
        "coverage_matrix": shape.coverage_matrix(reports),
        "call_distribution": shape.call_distribution(reports),
        "isolates": shape.isolate_rows(reports, drugs),
        "tiers": shape.tier_composition(reports, drugs),
        "discordance": shape.discordance_rows(reports),
        "discordance_tickets": shape.discordance_tickets(reports),
        "sample_details": shape.sample_details(reports),
        "vus": vus,
        "mechanisms": shape.mechanism_cards(reports),
        "alignment_loci": alignment,
        "mutation_index": shape.mutation_index(alignment, vus),
        "prevalence": shape.prevalence_series(prevalence_rows or []),
        "targets": shape.target_evidence(target_rows or []),
        "epistasis": shape.epistasis_notes(reports),
        "watchlist": shape.watchlist_rows(watchlist_rows or []),
        "external_benchmarks": shape.external_model_benchmarks(external_benchmark_rows or []),
        "workflows": workflows,
        "audit": shape.audit_events(audit_entries or []),
        "validation_outcomes": (validation or {}).get("outcomes", []),
        "validation_timeline": [{
            "date": row.get("date") or row.get("timestamp") or row.get("period") or "Date not supplied",
            "text": row.get("text") or " — ".join(filter(None, (row.get("title"), row.get("detail")))) or "Validation event",
        } for row in (validation or {}).get("timeline", [])],
        "federated_sites": [{
            "id": row.get("id") or row.get("site") or row.get("site_id") or "Unnamed site",
            "n": row.get("n") if row.get("n") is not None else row.get("submissions", 0),
            "k": row.get("k") if row.get("k") is not None else row.get("anonymity_threshold", 0),
            "status": humanise(row.get("status") or "Not supplied"),
        } for row in (federated_sites or [])],
        "error_trend": list(error_trend or []),
        "measurability": list(measurability or []),
        # One label source for both sides, so the runtime and the
        # renderer cannot drift into two spellings of one term.
        "labels": dict(EXPLICIT),
    }

    body = _asset("body.html")
    for token, value in _summary(payload, reports, drugs).items():
        body = body.replace(token, value)
    body = body.replace("__REVIEW_ROWS__", _review_rows(data))
    body = body.replace("__FOLLOWUP__", _followup(data))
    body = body.replace("__QC__", _qc(data))
    body = body.replace("__INTERPRETATION__", _interpretation(
        data, drugs, lineages, payload["measurability"]))
    body = body.replace("__METHODS__", _methods(
        data, drugs, len(reports), reference_method))

    # The demo warning is additive. Demo mode does not stop this being a
    # research scaffold, and replacing the clinical-use warning with it
    # removed the one line that must appear on every report.
    banner = "Research scaffold · not for clinical use"
    if data.get("demo_mode"):
        banner = ("Synthetic annotator — demo mode. Values are illustrative "
                  "and carry no biological meaning. " + banner)

    # json.dumps with no raw "<" so the payload cannot close the script early.
    encoded = json.dumps(payload, default=str).replace("<", "\\u003c")

    report_label = (f"cohort report, {len(reports)} isolates" if len(reports) > 1
                    else f"isolate {data.get('sample_id')}")

    return (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\n"
        f"<title>Myconductor evidence review — {_esc(report_label)}</title>\n"
        "<style>\n" + _asset("report.css") + "\n</style>\n</head>\n<body>\n"
        '<header class="report-header">\n'
        f'  <div class="scaffold-badge">{_esc(banner)}</div>\n'
        '  <div class="brand-row"><span class="brand-dot"></span>'
        '<span class="brand">Myconductor</span></div>\n'
        '  <h1 class="report-title">Evidence &amp; Benchmark Report</h1>\n'
        '  <p class="report-sub">Six-state call ontology · coverage-gated '
        'susceptibility · governed evidence. Nothing in this document is a '
        'clinical decision.</p>\n'
        f'  <div class="meta-grid">{_meta(data, drugs, len(reports), lineages)}</div>\n'
        + '  ' + _download_bar() + '\n'
        '</header>\n'
        + body +
        '\n<div class="tooltip" id="tooltip"></div>\n'
        '<aside class="drill-panel" id="drillPanel">\n'
        '  <div class="drill-head"><h3 id="drillTitle">—</h3>'
        '<button class="drill-close" id="drillClose">✕</button></div>\n'
        '  <div class="drill-body" id="drillBody"></div>\n'
        '</aside>\n'
        "<script>window.__MYCONDUCTOR__ = " + encoded + ";</script>\n"
        "<script>\n" + _asset("report.js") + "\n</script>\n"
        "</body>\n</html>\n"
    )
