"""Shaping a real analysis into the report's view model.

The design mock hard-coded seventeen datasets and manufactured four more with a
seeded PRNG — call matrices, coverage grids, lineage strata, sequence bases.
This module replaces all of that with what a run actually produced, and where a
run produced nothing it emits nothing so the page renders an empty state.

That rule is the whole point. A report full of plausible generated numbers is
worse than one with visible gaps: it looks identical to a real result, and
nothing downstream — a reader, a reviewer, an importer — can tell them apart.
It is the same failure as the hash-derived conservation scores this project
removed, wearing a nicer typeface.

Cohort inputs are optional. A single-sample analysis has no benchmark, no
lineage stratification and no rolling error curve, and says so rather than
borrowing a shape from somewhere else.
"""
from __future__ import annotations

from typing import Any, Iterable, Optional, Sequence

from .labels import humanise, humanise_all, variant_name

#: Call value -> the CSS class the mock's palette uses.
CALL_CLASS = {
    "resistant": "res", "susceptible": "sus", "indeterminate": "ind",
    "not_assessed": "na", "no_call": "nc", "unsupported": "un",
}

#: The six states in the order the legend lists them, with their palette token.
#: The six states in legend order: (display label, enum value, palette token).
#: The label is what a reader sees; the enum value is what the payload and
#: every export carry, so a consumer still joins on a stable key.
CALL_ORDER = (
    ("Resistant", "resistant", "var(--res)"),
    ("Susceptible", "susceptible", "var(--sus)"),
    ("Indeterminate", "indeterminate", "var(--ind)"),
    ("Not assessed", "not_assessed", "var(--na)"),
    ("No call", "no_call", "var(--nc)"),
    ("Unsupported", "unsupported", "var(--un)"),
)

TIER_COLOUR = {
    "catalogued": "var(--res)", "phenotypic": "#c17a86",
    "inferred": "var(--ind)", "predicted": "#8d9b8d", "none": "var(--na)",
}

#: Short codes for the drugs the bundled profile covers. A drug without one
#: falls back to its first three letters rather than being dropped.
DRUG_CODES = {
    "amikacin": "AMK", "bedaquiline": "BDQ", "capreomycin": "CAP",
    "clofazimine": "CFZ", "cycloserine": "CYC", "delamanid": "DLM",
    "ethambutol": "EMB", "ethionamide": "ETH", "isoniazid": "INH",
    "kanamycin": "KAN", "levofloxacin": "LFX", "linezolid": "LZD",
    "moxifloxacin": "MFX", "pretomanid": "PMD", "pyrazinamide": "PZA",
    "rifabutin": "RFB", "rifampicin": "RIF", "streptomycin": "STR",
}

DRUG_CLASS = {
    "isoniazid": "First-line", "rifampicin": "First-line",
    "ethambutol": "First-line", "pyrazinamide": "First-line",
    "rifabutin": "First-line",
    "levofloxacin": "Fluoroquinolone", "moxifloxacin": "Fluoroquinolone",
    "amikacin": "Injectable", "kanamycin": "Injectable",
    "capreomycin": "Injectable", "streptomycin": "Injectable",
    "bedaquiline": "BPaL/M", "pretomanid": "BPaL/M", "linezolid": "BPaL/M",
    "clofazimine": "BPaL/M", "delamanid": "BPaL/M",
    "ethionamide": "Second-line", "cycloserine": "Second-line",
}


def code_for(drug: str) -> str:
    return DRUG_CODES.get(drug.lower(), drug[:3].upper())


def _float(value: Any) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


# -- benchmark ------------------------------------------------------------
def drugs_from_baseline(rows: Iterable[dict]) -> list[dict]:
    """Per-drug benchmark cards from ``mycobench baseline`` rows.

    ``estimable`` is carried through rather than recomputed: the baseline
    already decided, on its own stated rule, whether a drug's point can be
    reported, and a second opinion here could disagree with the table it sits
    beside.
    """
    out = []
    for row in rows or []:
        drug = row.get("drug", "")
        if not drug:
            continue
        estimable = str(row.get("estimable", "")).lower() == "yes"
        error = _float(row.get("error_rate"))
        coverage = _float(row.get("coverage"))
        out.append({
            "id": code_for(drug),
            "name": drug.capitalize(),
            "cls": DRUG_CLASS.get(drug.lower(), "Other"),
            "coverage": round(coverage * 100, 1) if coverage is not None else None,
            "error": round(error * 100, 1) if (error is not None and estimable) else None,
            "vme": _int(row.get("n_false_susceptible")) if estimable else None,
            "me": _int(row.get("n_false_resistant")) if estimable else None,
            "n": _int(row.get("n_isolates")),
            "evaluable": _int(row.get("n_answered")) if estimable else 0,
            "estimable": estimable,
            # Carried for the drill panel; the mock had nowhere to put these.
            "sensitivity": _float(row.get("sensitivity")),
            "specificity": _float(row.get("specificity")),
            "sensitivity_ci": row.get("sensitivity_ci") or "",
            "specificity_ci": row.get("specificity_ci") or "",
            "ppv": _float(row.get("ppv")),
            "npv": _float(row.get("npv")),
            "prevalence": _float(row.get("resistance_prevalence")),
            "sens_basis": _int(row.get("sensitivity_basis_n")),
            "spec_basis": _int(row.get("specificity_basis_n")),
            "notes": [n for n in (row.get("notes") or "").split(" | ") if n],
        })
    order = {"estimable": 0}
    out.sort(key=lambda d: (not d["estimable"],
                            d["error"] if d["error"] is not None else 999))
    return out


def lineage_strata(rows: Iterable[dict]) -> tuple[list[str], dict]:
    """``(lineages, {drug_code: {lineage: {error, n}}})`` from lineage rows.

    Returns an empty lineage list when nothing was stratified, which the page
    renders as "lineage was not recorded" rather than as a single bucket. A
    lone ``unknown`` group is indistinguishable from a genuinely
    lineage-restricted cohort, and reporting it as the latter would be a claim
    nobody measured.
    """
    strata: dict[str, dict[str, dict]] = {}
    lineages: list[str] = []
    for row in rows or []:
        lineage = (row.get("lineage") or "").strip()
        drug = row.get("drug", "")
        if not lineage or lineage == "unknown" or not drug:
            continue
        if lineage not in lineages:
            lineages.append(lineage)
        sensitivity = _float(row.get("sensitivity"))
        specificity = _float(row.get("specificity"))
        vme = _float(row.get("vme_rate"))
        me = _float(row.get("me_rate"))
        # One number for the heatmap cell: the error rate this lineage saw.
        error = None
        if vme is not None and me is not None:
            error = round((vme + me) / 2 * 100, 1)
        elif vme is not None:
            error = round(vme * 100, 1)
        strata.setdefault(code_for(drug), {})[lineage] = {
            "error": error if error is not None else 0.0,
            "n": _int(row.get("n_called")),
            "sensitivity": sensitivity, "specificity": specificity,
            "powered": str(row.get("powered", "")).lower() == "yes",
        }
    return sorted(lineages), strata


def measurability_rows(rows: Iterable[dict]) -> list[dict]:
    return [r for r in (rows or []) if str(r.get("measurable", "")).lower() == "no"]


# -- single-sample --------------------------------------------------------
def call_distribution(reports: Sequence[dict]) -> list[dict]:
    """Counts per call state across every drug result supplied."""
    counts = {value: 0 for _label, value, _c in CALL_ORDER}
    for report in reports:
        for result in report.get("drug_results", []):
            call = result.get("call")
            if call in counts:
                counts[call] += 1
    return [{"k": label, "v": counts[value], "c": colour, "enum": value}
            for label, value, colour in CALL_ORDER if counts[value]]


def isolate_rows(reports: Sequence[dict], drugs: Sequence[dict]) -> list[dict]:
    """One call-matrix row per analysed sample."""
    codes = [d["id"] for d in drugs]
    rows = []
    for report in reports:
        calls = {}
        for result in report.get("drug_results", []):
            code = code_for(result.get("drug", ""))
            if code in codes:
                calls[code] = CALL_CLASS.get(result.get("call"), "na")
        if not calls:
            continue
        rows.append({
            "id": report.get("sample_id", "sample"),
            "lin": (report.get("context") or {}).get("lineage") or "untyped",
            "calls": calls,
        })
    return rows


def coverage_loci(reports: Sequence[dict]) -> list[dict]:
    """Per-locus callability, aggregated over the supplied reports."""
    seen: dict[str, dict] = {}
    for report in reports:
        for result in report.get("drug_results", []):
            drug = code_for(result.get("drug", ""))
            for coverage in result.get("coverage", []) or []:
                locus = coverage.get("locus")
                if not locus:
                    continue
                entry = seen.setdefault(
                    locus, {"id": locus, "drugs": set(), "callable": [],
                            "unknown": 0})
                entry["drugs"].add(drug)
                fraction = _float(coverage.get("callable_fraction"))
                if fraction is None:
                    entry["unknown"] += 1
                else:
                    entry["callable"].append(fraction)
    out = []
    for locus, entry in sorted(seen.items()):
        values = entry["callable"]
        # Share of observations at or above the 95% floor. With no known
        # fraction this stays None and the bar renders as unknown, not zero.
        callable_pct = (round(100 * sum(1 for v in values if v >= 0.95) / len(values))
                        if values else None)
        out.append({"id": locus, "drug": "/".join(sorted(entry["drugs"])),
                    "callable": callable_pct})
    return out


def coverage_matrix(reports: Sequence[dict]) -> list[dict]:
    """Per-sample, per-locus callable fraction for the heatmap."""
    rows = []
    for report in reports:
        loci: dict[str, Optional[float]] = {}
        for result in report.get("drug_results", []):
            for coverage in result.get("coverage", []) or []:
                locus = coverage.get("locus")
                if locus:
                    loci[locus] = _float(coverage.get("callable_fraction"))
        if loci:
            rows.append({"id": report.get("sample_id", "sample"), "loci": loci})
    return rows


def tier_composition(reports: Sequence[dict], drugs: Sequence[dict]) -> dict:
    """``{drug_code: {tier: share}}`` over the evidence actually attached."""
    codes = {d["id"] for d in drugs}
    counts: dict[str, dict[str, int]] = {}
    for report in reports:
        for result in report.get("drug_results", []):
            code = code_for(result.get("drug", ""))
            if code not in codes:
                continue
            bucket = counts.setdefault(code, {})
            evidence = result.get("evidence") or []
            if not evidence:
                bucket["none"] = bucket.get("none", 0) + 1
                continue
            for item in evidence:
                tier = item.get("tier", "none")
                bucket[tier] = bucket.get(tier, 0) + 1
    out = {}
    for code, bucket in counts.items():
        total = sum(bucket.values()) or 1
        out[code] = {tier: round(100 * n / total, 1)
                     for tier, n in bucket.items()}
    return out


def discordance_rows(reports: Sequence[dict]) -> list[dict]:
    """Conflicts, with the reconciled call and the context beside them."""
    rows = []
    for report in reports:
        for item in report.get("discordances", []) or []:
            sources = list(item.get("sources") or [])
            calls = list(item.get("calls") or [])
            rows.append({
                "var": item.get("drug", ""),
                "drug": code_for(item.get("drug", "")),
                "prof": humanise(calls[0]) if calls else "—",
                "myk": humanise(calls[1]) if len(calls) > 1 else "—",
                "rec": humanise("indeterminate"),
                "reason": item.get("note", ""),
                "sources": ", ".join(sources),
                # Populated by modules/discordance_context: model and
                # structural signals shown beside the conflict, never
                # resolving it.
                "context": list(item.get("context") or []),
            })
    return rows


def vus_items(reports: Sequence[dict]) -> list[dict]:
    out = []
    for report in reports:
        for rank, item in enumerate(report.get("vus_priorities", []) or [], 1):
            dimensions = {d.get("name"): d for d in (item.get("dimensions") or [])}
            available = [n for n, d in dimensions.items() if d.get("available")]
            out.append({
                "rank": rank,
                "variant": variant_name(item.get("variant_label", "")),
                "gene": item.get("gene", ""),
                "drug": code_for(item.get("drug") or ""),
                # The workbench's own word, not a number invented here. A
                # score of None means it declined to rank, and the card says so.
                "priority": humanise(item.get("priority", "insufficient-data")),
                "score": item.get("score"),
                "features": {
                    "dimensions": humanise_all(available) or "None available",
                    # The available dimensions were humanised in the first
                    # pass but the gaps were not, so the drill panel still
                    # listed eight raw dimension names. The page-level scan
                    # could not see it: the panel only exists once opened.
                    "gaps": ("; ".join(humanise(g) for g in
                                       (item.get("data_gaps") or []))
                             or "None recorded"),
                    "experiment": item.get("recommended_experiment") or "not specified",
                },
            })
    return out


def mechanism_cards(reports: Sequence[dict]) -> list[dict]:
    cards = []
    for report in reports:
        for item in report.get("mechanism_queue", []) or []:
            tags = [{"t": "ind", "l": humanise("indeterminate")}]
            if item.get("gene"):
                tags.append({"t": "accent", "l": item["gene"]})
            if item.get("evidence_gaps"):
                tags.append({"t": "na", "l": "evidence gap"})
            body = item.get("hypothesis", "")
            if item.get("evidence_gaps"):
                body += " Missing: " + "; ".join(item["evidence_gaps"])
            if item.get("resolving_experiments"):
                body += " Resolved by: " + "; ".join(item["resolving_experiments"])
            cards.append({
                "title": (f"{variant_name(item.get('variant_label', ''))} — "
                          f"{humanise(item.get('drug', ''))}"),
                "body": body, "tags": tags,
            })
    return cards


def alignment_loci(reports: Sequence[dict]) -> dict:
    """Real called variants, grouped by gene, for the alignment view.

    Only coordinate-resolved records appear. A variant whose position could not
    be resolved has no place on a positional axis, and guessing one would be
    the invention this module exists to prevent.
    """
    loci: dict[str, dict] = {}
    for report in reports:
        sample = report.get("sample_id", "sample")
        lineage = (report.get("context") or {}).get("lineage") or "untyped"
        carried: dict[str, list[str]] = {}
        for result in report.get("drug_results", []):
            drug = result.get("drug", "")
            for item in result.get("evidence") or []:
                variant = item.get("variant") or {}
                gene = variant.get("gene")
                pos = variant.get("pos")
                ref, alt = variant.get("ref"), variant.get("alt")
                if not gene or pos is None or not ref or not alt:
                    continue
                entry = loci.setdefault(gene, {
                    "name": gene, "assembly": variant.get("assembly", "reference"),
                    "positions": {}, "isolates": {},
                })
                record = entry["positions"].setdefault(int(pos), {
                    "pos": int(pos), "ref": ref, "alt": alt,
                    "label": item.get("variant_label") or "",
                    "catalogued": item.get("tier") == "catalogued",
                    "drugs": set(),
                })
                record["drugs"].add(drug)
                carried.setdefault(gene, []).append(str(int(pos)))
        for gene, positions in carried.items():
            loci[gene]["isolates"][sample] = {
                "id": sample, "lineage": lineage, "carried": sorted(set(positions))}

    out = {}
    for gene, entry in sorted(loci.items()):
        positions = sorted(entry["positions"].values(), key=lambda p: p["pos"])
        for record in positions:
            record["drugs"] = sorted(record["drugs"])
        out[gene] = {
            "name": entry["name"], "assembly": entry["assembly"],
            "positions": positions,
            "isolates": list(entry["isolates"].values()),
        }
    return out


def audit_events(entries: Iterable[dict]) -> list[dict]:
    out = []
    for entry in entries or []:
        out.append({
            "ts": (entry.get("timestamp_utc") or entry.get("timestamp")
                   or entry.get("ts") or ""),
            "event": humanise(entry.get("action") or entry.get("event") or ""),
            "actor": entry.get("actor") or entry.get("reviewer") or "unknown",
            "target": entry.get("target") or entry.get("summary") or "",
            "hash": (entry.get("hash") or entry.get("entry_hash") or "")[:8],
        })
    return out
