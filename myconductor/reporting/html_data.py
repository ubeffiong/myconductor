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

import re
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


def display_variant(value: Any) -> str:
    """Human label for either a gene/change token or a coordinate identity."""
    raw = str(value or "")
    match = re.search(r":(\d+):([^:>]+)>([^:>]+)$", raw)
    if match:
        pos, ref, alt = match.groups()
        return f"{ref} to {alt} at genomic position {int(pos):,}"
    return variant_name(raw)


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


def _optional_int(value: Any) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


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
                "drug": humanise(item.get("drug") or ""),
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
    samples = [{
        "id": report.get("sample_id", "sample"),
        "lineage": (report.get("context") or {}).get("lineage") or "untyped",
        "site": (report.get("context") or {}).get("site_id") or "site not recorded",
    } for report in reports]

    for report in reports:
        sample = report.get("sample_id", "sample")
        carried: dict[str, dict[str, dict]] = {}
        for result in report.get("drug_results", []):
            drug = result.get("drug", "")
            for item in result.get("evidence") or []:
                variant = item.get("variant") or {}
                gene = variant.get("gene")
                pos = variant.get("pos")
                ref, alt = variant.get("ref"), variant.get("alt")
                if not gene or pos is None or not ref or not alt:
                    continue
                position = int(pos)
                change = (item.get("variant_label") or variant.get("hgvs_p")
                          or variant.get("hgvs_c") or f"{ref}{position}{alt}")
                entry = loci.setdefault(gene, {
                    "name": gene, "assembly": variant.get("assembly", "reference"),
                    "positions": {}, "samples": {},
                })
                record = entry["positions"].setdefault(position, {
                    "pos": position, "ref": ref, "alt": alt,
                    "label": change, "change": change,
                    "display": variant_name(f"{gene}_{change}"),
                    "consequence": humanise(variant.get("consequence") or "not supplied"),
                    "region": humanise(variant.get("region") or "not supplied"),
                    "catalogued": False, "drugs": set(), "tiers": set(),
                    "lanes": set(), "calls": set(), "rationales": set(),
                    "limitations": set(), "evidence_ids": set(), "carriers": set(),
                    "depths": [], "vafs": [],
                })
                record["catalogued"] = record["catalogued"] or item.get("tier") == "catalogued"
                record["drugs"].add(humanise(drug))
                record["tiers"].add(humanise(item.get("tier") or "unknown"))
                record["lanes"].add(humanise(item.get("lane") or "unknown"))
                record["calls"].add(humanise(item.get("call") or "unknown"))
                if item.get("rationale"):
                    record["rationales"].add(str(item["rationale"]))
                for limitation in item.get("limitations") or []:
                    record["limitations"].add(humanise(limitation))
                if item.get("evidence_id"):
                    record["evidence_ids"].add(str(item["evidence_id"])[:12])
                if variant.get("depth") is not None:
                    record["depths"].append(variant.get("depth"))
                if variant.get("vaf") is not None:
                    record["vafs"].append(variant.get("vaf"))
                record["carriers"].add(sample)
                carried.setdefault(gene, {})[str(position)] = {
                    "pos": str(position), "alt": alt, "depth": variant.get("depth"),
                    "vaf": variant.get("vaf"), "tier": humanise(item.get("tier") or "unknown"),
                    "lane": humanise(item.get("lane") or "unknown"),
                    "call": humanise(item.get("call") or "unknown"),
                }
        for gene, positions in carried.items():
            entry = loci.setdefault(gene, {"name": gene, "assembly": "reference", "positions": {}, "samples": {}})
            sample_entry = entry["samples"].setdefault(sample, {"carried": {}})
            sample_entry["carried"].update(positions)

    out = {}
    n_samples = max(1, len(samples))
    for gene, entry in sorted(loci.items()):
        positions = sorted(entry["positions"].values(), key=lambda p: p["pos"])
        for record in positions:
            depths = [d for d in record.pop("depths") if d is not None]
            vafs = [v for v in record.pop("vafs") if v is not None]
            carriers = sorted(record.pop("carriers"))
            record["drugs"] = sorted(record["drugs"])
            record["tiers"] = sorted(record["tiers"])
            record["lanes"] = sorted(record["lanes"])
            record["calls"] = sorted(record["calls"])
            record["rationales"] = sorted(record["rationales"])[:3]
            record["limitations"] = sorted(record["limitations"])[:4]
            record["evidence_ids"] = sorted(record["evidence_ids"])[:4]
            record["carrier_count"] = len(carriers)
            record["carrier_fraction"] = len(carriers) / n_samples
            record["median_depth"] = sorted(depths)[len(depths)//2] if depths else None
            record["mean_vaf"] = round(sum(vafs) / len(vafs), 3) if vafs else None
        isolates = []
        for sample in samples:
            sample_entry = entry.get("samples", {}).get(sample["id"], {})
            carried_map = sample_entry.get("carried", {})
            isolates.append({
                "id": sample["id"], "lineage": sample["lineage"], "site": sample["site"],
                "carried": sorted(carried_map), "variants": carried_map,
            })
        out[gene] = {
            "name": entry["name"], "assembly": entry["assembly"],
            "positions": positions, "isolates": isolates,
        }
    return out



def external_model_benchmarks(rows: Iterable[dict]) -> list[dict]:
    """Normalise external-model benchmark rows for the HTML report."""
    out = []
    for row in rows or []:
        verdict = str(row.get("verdict") or row.get("baseline_verdict") or "not supplied")
        out.append({
            "model": str(row.get("model") or row.get("model_id") or "External model"),
            "version": str(row.get("version") or row.get("model_version") or "not supplied"),
            "drug": humanise(row.get("drug") or "not supplied"),
            "lineage": humanise(row.get("lineage") or "All lineages"),
            "cohort": str(row.get("cohort") or "not supplied"),
            "source": str(row.get("source") or "not supplied"),
            "call_rate": _percent_value(row.get("call_rate")),
            "error_rate": _percent_value(row.get("error_rate")),
            "baseline_coverage": _percent_value(row.get("baseline_coverage") or row.get("baseline_call_rate")),
            "baseline_error_rate": _percent_value(row.get("baseline_error_rate") or row.get("baseline_risk")),
            "matched_coverage": _percent_value(row.get("matched_coverage")),
            "matched_error_rate": _percent_value(row.get("matched_error_rate")),
            "sensitivity": _percent_value(row.get("sensitivity")),
            "specificity": _percent_value(row.get("specificity")),
            "n_evaluable": _int_value(row.get("n_evaluable") or row.get("n_predictions")),
            "n_called": _int_value(row.get("n_called")),
            "n_dropped": _int_value(row.get("n_dropped")),
            "verdict": humanise(verdict),
            "verdict_key": verdict.lower().replace("_", "-").replace(" ", "-"),
            "registry_ready": str(row.get("registry_ready") or "").lower() in {"yes", "true", "1"},
            "notes": str(row.get("notes") or row.get("reason") or "No benchmark note supplied"),
        })
    return out


def _percent_value(value) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return round(numeric * 100, 2) if numeric <= 1 else round(numeric, 2)


def _int_value(value) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0

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


def mutation_index(alignment: dict, vus: Sequence[dict]) -> list[dict]:
    """Searchable, deduplicated mutation records from measured report data."""
    records: dict[tuple, dict] = {}
    for gene, locus in alignment.items():
        for item in locus.get("positions", []):
            key = (gene, item.get("pos"), item.get("ref"), item.get("alt"))
            records[key] = {
                "gene": gene,
                "variant": (variant_name(item.get("label")) if item.get("label")
                            else (f"{item.get('ref', '')} to {item.get('alt', '')} "
                                  f"at genomic position {item.get('pos', '')}")),
                "position": item.get("pos"),
                "reference": item.get("ref"),
                "alternate": item.get("alt"),
                "drugs": [humanise(d) for d in item.get("drugs", [])],
                "classification": ("Catalogued resistance-associated" if
                                   item.get("catalogued") else
                                   "Not established by the catalogue"),
                "source": "Called variant evidence",
            }
    for item in vus:
        key = (item.get("gene"), None, item.get("variant"), None)
        if key not in records:
            records[key] = {
                "gene": item.get("gene") or "Unresolved locus",
                "variant": item.get("variant") or "Unlabelled variant",
                "position": None, "reference": None, "alternate": None,
                "drugs": [item.get("drug")] if item.get("drug") else [],
                "classification": "Variant of uncertain significance",
                "source": "Laboratory validation queue",
            }
    return sorted(records.values(), key=lambda r: (
        str(r["gene"]).lower(), r["position"] is None,
        r["position"] or 0, str(r["variant"])))


def prevalence_series(rows: Iterable[dict]) -> list[dict]:
    """Normalise supplied cohort prevalence; never infer a denominator."""
    out = []
    for row in rows or []:
        raw_total = row.get("total") if row.get("total") is not None else row.get("n_tested")
        raw_resistant = (row.get("resistant") if row.get("resistant") is not None
                         else row.get("n_resistant"))
        total = _optional_int(raw_total)
        resistant = _optional_int(raw_resistant)
        rate = _float(row.get("prevalence") or row.get("rate"))
        if rate is not None and rate <= 1:
            rate *= 100
        if rate is None and total and resistant is not None:
            rate = 100 * resistant / total
        out.append({
            "period": str(row.get("period") or row.get("date") or "Unspecified"),
            "drug": humanise(row.get("drug") or "All reported drugs"),
            "lineage": humanise(row.get("lineage") or "All lineages"),
            "geography": str(row.get("geography") or row.get("site") or
                             "All participating sites"),
            "resistant": resistant, "total": total, "rate": rate,
            "source": str(row.get("source") or "Report context dataset"),
        })
    return out


def target_evidence(rows: Iterable[dict]) -> list[dict]:
    """Normalise externally supplied target-liability evidence for review."""
    out = []
    for row in rows or []:
        name = row.get("target") or row.get("gene")
        if not name:
            continue
        out.append({
            "target": str(name),
            "essentiality": humanise(row.get("essentiality") or "Not supplied"),
            "druggability": humanise(row.get("druggability") or "Not supplied"),
            "human_homology": humanise(row.get("human_homology") or "Not supplied"),
            "resistance_liability": humanise(
                row.get("resistance_liability") or "Not supplied"),
            "evidence": str(row.get("evidence") or row.get("rationale") or
                            "No supporting evidence supplied"),
            "source": str(row.get("source") or "Report context dataset"),
        })
    return out


def epistasis_notes(reports: Sequence[dict]) -> list[dict]:
    """Render supplied co-observation rules without changing call confidence."""
    out = []
    for report in reports:
        for item in report.get("epistasis_notes", []) or []:
            pairs = item.get("matched_pairs") or []
            variants = []
            for pair in pairs:
                variants.append(" + ".join(display_variant(v) for v in pair))
            out.append({
                "drug": humanise(item.get("drug") or "Drug not recorded"),
                "interaction": humanise(item.get("interaction") or "annotation only"),
                "primary": variant_name(item.get("primary_variant_key_or_gene") or ""),
                "partner": variant_name(item.get("partner_variant_key_or_gene") or ""),
                "matched_pairs": variants,
                "note": str(item.get("note") or item.get("interpretation") or
                            "Co-observed variants matched a supplied rule."),
                "source": str(item.get("source") or "Supplied epistasis table"),
                "rule_id": str(item.get("rule_id") or ""),
                "table_version": str(item.get("table_version") or ""),
                "effect": "Annotation only",
            })
    return out


def watchlist_rows(rows: Iterable[dict]) -> list[dict]:
    """Normalise privacy-gated collaborative watch-list aggregates."""
    out = []
    for row in rows or []:
        drug = row.get("drug")
        if not drug:
            continue
        out.append({
            "drug": humanise(drug),
            "variant": variant_name(row.get("variant") or
                                    row.get("variant_label") or "Aggregate drug-level signal"),
            "unresolved_isolates": _optional_int(row.get("unresolved_isolates") or
                                                 row.get("isolate_count") or
                                                 row.get("n_isolates")),
            "contributing_sites": _optional_int(row.get("contributing_sites") or
                                                row.get("site_count") or
                                                row.get("n_sites")),
            "status": humanise(row.get("status") or "Under investigation"),
            "evidence_gaps": humanise_all(row.get("evidence_gaps") or []) or
                             str(row.get("evidence_gap") or "Not supplied"),
            "source": str(row.get("source") or "Privacy-gated report context"),
            "limitations": str(row.get("limitations") or
                               "Aggregate counts only; no isolate, patient or site identifier is exported."),
        })
    return out


def discordance_tickets(reports: Sequence[dict]) -> list[dict]:
    """Per-isolate genomic/phenotypic disagreement tickets from this run.

    Distinct from ``discordance_rows`` (cross-tool disagreement inside one
    report) and from ``watchlist_rows`` (externally supplied, privacy-gated
    cross-site aggregates): this reads ``AnalysisReport.discordance_tickets``
    directly — the tickets *this run itself* opened when its own genomic and
    phenotypic calls for one drug disagreed (see
    ``federated/watchlist.py::WatchlistStore.capture``).
    """
    out = []
    for report in reports:
        for item in report.get("discordance_tickets", []) or []:
            status = str(item.get("status") or "open")
            out.append({
                "ticket_id": str(item.get("ticket_id") or "")[:12],
                "sample": item.get("sample_id") or report.get("sample_id") or "Unknown sample",
                "drug": humanise(item.get("drug") or ""),
                "genomic_call": humanise(item.get("genomic_call") or ""),
                "genomic_call_key": CALL_CLASS.get(item.get("genomic_call"), "na"),
                "phenotypic_call": humanise(item.get("phenotypic_call") or ""),
                "phenotypic_call_key": CALL_CLASS.get(item.get("phenotypic_call"), "na"),
                "status": humanise(status),
                "status_key": status.lower().replace(" ", "_"),
                "candidate_mechanisms": humanise_all(
                    item.get("candidate_mechanisms") or []) or "None recorded",
                "evidence_gaps": humanise_all(
                    item.get("evidence_gaps") or []) or "Not supplied",
                "resolution": str(item.get("resolution") or ""),
                "source": str(item.get("source") or "Not supplied"),
            })
    return out


def sample_details(reports: Sequence[dict]) -> dict[str, dict]:
    """Per-sample drilldown: site, lineage, drug calls, called variants and
    locus coverage — read from the same per-report payload every cohort-level
    view already reads. Selecting one sample never re-derives or re-fetches
    anything; it only narrows the same data to one ``sample_id``.
    """
    out: dict[str, dict] = {}
    for report in reports:
        sample_id = report.get("sample_id") or "sample"
        context = report.get("context") or {}
        drug_calls = []
        coverage_by_locus: dict[str, dict] = {}
        variants_seen: dict[str, dict] = {}
        for result in report.get("drug_results", []):
            drug_calls.append({
                "drug": humanise(result.get("drug") or ""),
                "call": humanise(result.get("call") or ""),
                "call_key": CALL_CLASS.get(result.get("call"), "na"),
                "tier": humanise(result.get("tier") or ""),
                "reason": str(result.get("reason") or ""),
            })
            for cov in result.get("coverage", []) or []:
                locus = cov.get("locus")
                if not locus:
                    continue
                coverage_by_locus[locus] = {
                    "locus": locus,
                    "callable_fraction": _float(cov.get("callable_fraction")),
                    "mean_depth": cov.get("mean_depth"),
                }
            for item in result.get("evidence", []) or []:
                variant = item.get("variant") or {}
                key = item.get("variant_label") or variant.get("gene")
                if not key or key in variants_seen:
                    continue
                variants_seen[key] = {
                    "variant": display_variant(item.get("variant_label") or key),
                    "gene": variant.get("gene") or "Not supplied",
                    "consequence": humanise(variant.get("consequence") or "not supplied"),
                    "drug": humanise(result.get("drug") or ""),
                    "call": humanise(item.get("call") or ""),
                    "tier": humanise(item.get("tier") or ""),
                }
        out[sample_id] = {
            "sample_id": sample_id,
            "site": context.get("site_id") or "Not recorded",
            "lineage": context.get("lineage") or "Untyped",
            "organism": context.get("organism") or "Not recorded",
            "assay": humanise(context.get("assay") or "unknown"),
            "drug_calls": sorted(drug_calls, key=lambda d: d["drug"]),
            "variants": sorted(variants_seen.values(), key=lambda v: v["gene"]),
            "coverage": sorted(coverage_by_locus.values(), key=lambda c: c["locus"]),
        }
    return out


def implemented_workflows(reports: Sequence[dict]) -> dict[str, list[dict]]:
    """Shape advanced evidence produced by the analysis pipeline.

    These records used to survive only in the JSON report.  Keeping the
    transformation here gives the HTML a small, human-readable view model and
    prevents internal field names from leaking into visible tables.
    """
    shaped = {key: [] for key in (
        "population", "mic", "structural", "regulatory", "expression",
        "models", "panels")}
    gene_list = lambda values: ", ".join(str(item) for item in (values or []))
    for report in reports:
        population = report.get("population_structure")
        if population:
            shaped["population"].append({
                "sample": report.get("sample_id") or "Unknown sample",
                "classification": humanise(population.get("classification") or "indeterminate"),
                "basis": str(population.get("classification_basis") or "No classification basis supplied"),
                "method": str(population.get("method") or "Not supplied"),
                "groups": [{
                    "fraction": _float(group.get("estimated_fraction")),
                    "lineage": humanise(group.get("lineage") or "Not assigned"),
                    "variants": [display_variant(v) for v in group.get("variant_keys", [])],
                    "linkage": str(group.get("linkage_source") or "No molecular linkage supplied"),
                    "note": str(group.get("note") or "Frequency group only; not a reconstructed clone."),
                } for group in population.get("subpopulations", [])],
                "unclustered": [display_variant(v) for v in population.get("unclustered_variant_keys", [])],
                "gaps": humanise_all(population.get("data_gaps", [])),
            })
        for row in report.get("quantitative_findings", []) or []:
            interval = row.get("interval")
            shaped["mic"].append({
                "drug": humanise(row.get("drug")), "value": row.get("value"),
                "unit": row.get("unit"), "interval": interval,
                "comparison": humanise(row.get("comparison")),
                "critical_concentration": row.get("critical_concentration"),
                "method": str(row.get("method") or "Not supplied"),
                "source": str(row.get("source") or "Not supplied"),
                "conflict": bool(row.get("conflict")),
                "interpretation": str(row.get("interpretation") or ""),
            })
        for row in report.get("structural_annotations", []) or []:
            distance = None
            if row.get("ligand_distance") is not None:
                distance = f"{row['ligand_distance']} {row.get('distance_unit') or ''}".strip()
            shaped["structural"].append({
                "variant": display_variant(row.get("variant_key") or ""),
                "gene": str(row.get("gene") or "Not supplied"),
                "location": humanise(row.get("location") or "unknown"),
                "effect": str(row.get("predicted_effect") or "Not supplied"),
                "distance": distance or "Not measured",
                "status": humanise(row.get("validation_status") or row.get("status") or "unknown"),
                "source": str(row.get("source") or "Not supplied"),
                "ranking_effect": humanise(row.get("ranking_effect") or "none"),
            })
        for row in report.get("regulatory_findings", []) or []:
            shaped["regulatory"].append({
                "name": str(row.get("name") or "Unnamed region"),
                "variant": display_variant(row.get("variant_key") or ""),
                "type": humanise(row.get("type") or "unknown"),
                "targets": humanise_all(row.get("target_genes", [])),
                "drugs": humanise_all(row.get("drug_associations", [])),
                "tier": humanise(row.get("evidence_tier") or "unknown"),
                "source": str(row.get("source") or "Not supplied"),
                "interpretation": "Region membership is contextual evidence and does not establish expression or resistance.",
            })
        for row in report.get("expression_findings", []) or []:
            shaped["expression"].append({
                "gene": str(row.get("gene") or "Not supplied"),
                "measurement": humanise(row.get("measurement_type") or "unknown"),
                "fold_change": row.get("fold_change"), "unit": row.get("unit"),
                "conclusion": humanise(row.get("reported_conclusion") or "uncertain"),
                "matching": humanise(row.get("matching_status") or "unlinked"),
                "source": str(row.get("provenance") or "Not supplied"),
                "interpretation": str(row.get("interpretation") or "Expression evidence is contextual and does not independently establish a drug call."),
            })
        for row in report.get("in_silico_findings", []) or []:
            basis = row.get("baseline_basis") or {}
            shaped["models"].append({
                "variant": display_variant(row.get("variant_key") or ""),
                "drug": humanise(row.get("drug")),
                "prediction": humanise(row.get("prediction") or "uncertain"),
                "confidence": row.get("confidence"),
                "model": f"{row.get('model_id', 'Unknown model')} {row.get('model_version', '')}".strip(),
                "cohort": str(basis.get("cohort") or "Not supplied"),
                "baseline": str(basis.get("baseline") or "Not supplied"),
                "basis": str(basis.get("basis") or "No governed comparison supplied"),
                "effect": "No call or ranking effect",
                "interpretation": str(row.get("interpretation") or ""),
            })
        panel = report.get("panel") or {}
        if panel:
            shaped["panels"].append({
                "name": str(panel.get("name") or "Unnamed panel"),
                "version": str(panel.get("version") or "Not supplied"),
                "assay": humanise(panel.get("assay") or "unknown"),
                "verified": bool(panel.get("verified")),
                "loci": panel.get("n_loci"),
                "kept": gene_list(panel.get("kept", [])),
                "discarded": gene_list(panel.get("discarded", [])),
                "unsupported": [{"drug": humanise(drug), "missing": gene_list(missing)}
                                for drug, missing in sorted((panel.get("unsupported_drugs") or {}).items())],
                "off_panel": ", ".join(display_variant(v) for v in panel.get("off_panel_variants", [])),
                "note": str(panel.get("note") or ""),
                "source": str(panel.get("source") or "Not supplied"),
            })
    return shaped
