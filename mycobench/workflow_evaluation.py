"""Paired evaluation of baseline versus MyConductor-assisted review.

The unit is site/isolate/drug. Unresolved results remain in the denominator.
Operational measurements are compared only when both arms have observations.
No before/after causal or clinical benefit is asserted by this evaluator.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from html import escape
from pathlib import Path

from myconductor.core.models import Call
from myconductor.reporting.json_report import atomic_json

CALLS = {c.value for c in Call}
MEASURES = ("reviewer_minutes", "turnaround_hours", "additional_tests", "cost")


def _number(row, field):
    value = row.get(field)
    if value in ("", None):
        return None
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{field} must be finite and nonnegative")
    return value


def validate(rows):
    seen = set()
    groups = {}
    for row in rows:
        for name in ("site_id", "isolate_id", "drug", "truth", "truth_quality",
                     "truth_method", "truth_source", "training_overlap",
                     "independent_adjudication", "baseline_call", "intervention_call"):
            if not row.get(name):
                raise ValueError(f"evaluation row requires {name}")
        key = (row["site_id"], row["isolate_id"], row["drug"])
        if key in seen:
            raise ValueError(f"duplicate isolate-drug observation: {key}")
        seen.add(key)
        if row.get("partition") not in (None, "", "development", "evaluation"):
            raise ValueError("partition must be development or evaluation")
        if row["truth"] not in ("resistant", "susceptible"):
            raise ValueError("truth must be a reviewed binary phenotype")
        if row["truth_quality"] not in ("pass", "fail", "unknown"):
            raise ValueError("unknown truth quality")
        if row["training_overlap"] not in ("yes", "no", "unknown"):
            raise ValueError("training overlap must be yes/no/unknown")
        if row["independent_adjudication"] not in ("yes", "no", "unknown"):
            raise ValueError("independent adjudication must be yes/no/unknown")
        for arm in ("baseline", "intervention"):
            if row[arm + "_call"] not in CALLS:
                raise ValueError(f"unknown {arm} call")
            for measure in MEASURES:
                _number(row, arm + "_" + measure)
        if any(_number(row, arm + "_cost") is not None for arm in ("baseline", "intervention")) and not row.get("currency"):
            raise ValueError("cost observations require currency")
        # A declared genetic cluster is global; patient IDs are namespaced to
        # their site. No relatedness is guessed when these fields are missing.
        for field in ("patient_id", "cluster_id"):
            if row.get(field) and row.get("partition"):
                group = (field, row["site_id"] if field == "patient_id" else "", row[field])
                if group in groups and groups[group] != row["partition"]:
                    raise ValueError(f"{field} crosses development/evaluation partitions")
                groups[group] = row["partition"]


def _rate(numerator, denominator):
    return numerator / denominator if denominator else None


def accuracy(rows, arm):
    counts = dict(tp=0, tn=0, fp=0, fn=0, abstained_resistant=0, abstained_susceptible=0)
    for row in rows:
        predicted = row[arm + "_call"]
        resistant = row["truth"] == "resistant"
        if predicted not in ("resistant", "susceptible"):
            counts["abstained_resistant" if resistant else "abstained_susceptible"] += 1
        else:
            counts[("tp" if resistant else "fp") if predicted == "resistant"
                   else ("fn" if resistant else "tn")] += 1
    tp, tn, fp, fn = (counts[k] for k in ("tp", "tn", "fp", "fn"))
    nr = tp + fn + counts["abstained_resistant"]
    ns = tn + fp + counts["abstained_susceptible"]
    called = tp + tn + fp + fn
    return dict(counts, n=len(rows), n_called=called,
                call_rate=_rate(called, len(rows)),
                conditional_accuracy=_rate(tp + tn, called),
                conditional_sensitivity=_rate(tp, tp + fn),
                conditional_specificity=_rate(tn, tn + fp),
                resistant_detection_yield=_rate(tp, nr),
                susceptible_detection_yield=_rate(tn, ns),
                false_susceptible_per_resistant=_rate(fn, nr),
                false_resistant_per_susceptible=_rate(fp, ns),
                unresolved_rate=_rate(len(rows) - called, len(rows)))


def paired_operational(rows):
    out = {}
    for measure in MEASURES:
        groups = defaultdict(list)
        for row in rows:
            before = _number(row, "baseline_" + measure)
            after = _number(row, "intervention_" + measure)
            if before is not None and after is not None:
                groups[row.get("currency", "") if measure == "cost" else ""].append(after - before)
        out[measure] = {
            unit or "native_unit": {"n_pairs": len(values), "mean_change": statistics.fmean(values),
                                   "median_change": statistics.median(values)}
            for unit, values in groups.items()
        }
        # Empty maps mean unavailable; never report an unmeasured saving of 0.
    return out


def evaluate(rows):
    rows = list(rows)
    validate(rows)
    quality = [r for r in rows if r["truth_quality"] == "pass"]
    independent = [r for r in quality if r["training_overlap"] == "no"
                   and r["independent_adjudication"] == "yes"
                   and r.get("partition") == "evaluation"]
    strata = {}
    for field in ("drug", "lineage", "site_id", "country"):
        groups = defaultdict(list)
        for row in independent:
            groups[row.get(field) or "unknown"].append(row)
        strata[field] = {name: {a: accuracy(group, a) for a in ("baseline", "intervention")}
                         for name, group in sorted(groups.items())}
    discordant_pairs = {"improved": 0, "worsened": 0, "unchanged": 0}
    for row in independent:
        before = row["baseline_call"] == row["truth"]
        after = row["intervention_call"] == row["truth"]
        discordant_pairs["improved" if after and not before else
                         "worsened" if before and not after else "unchanged"] += 1
    return {
        "schema": "mycobench.workflow-evaluation.v1",
        "n_submitted": len(rows), "n_quality_reviewed": len(quality),
        "n_independent_evaluation": len(independent),
        "excluded_quality": len(rows) - len(quality),
        "excluded_overlap_adjudication_or_partition": len(quality) - len(independent),
        "primary": {a: accuracy(independent, a) for a in ("baseline", "intervention")},
        "descriptive_all_quality_pass": {a: accuracy(quality, a) for a in ("baseline", "intervention")},
        "strata": strata, "paired_correctness": discordant_pairs,
        "operational": paired_operational(independent),
        "limitations": [
            "Paired descriptive comparison, not evidence of a causal intervention benefit.",
            "No deployment prevalence or clinical outcome is inferred from selected cohorts.",
            "Training/catalogue overlap and independent adjudication are externally supplied declarations.",
            "Relatedness-aware confidence intervals require a prespecified independent cluster design; none are fabricated.",
            f"{sum(not r.get('cluster_id') for r in independent)} primary rows lack genetic-cluster metadata.",
            "Abstentions are retained; conditional accuracy must be read with call rate and detection yield.",
            "Operational observations are per isolate-drug review; do not duplicate whole-case costs across drugs.",
        ],
    }


def run(args):
    path = Path(args.input)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="," if path.suffix == ".csv" else "\t"))
    result = evaluate(rows)
    result["input_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    atomic_json(args.out, result)
    if args.html:
        body = "<h1>Paired laboratory workflow evaluation</h1>"
        body += "<p>Independent evaluable isolate–drug pairs: " + str(result["n_independent_evaluation"]) + "</p>"
        body += "<table><tr><th>Metric</th><th>Baseline</th><th>MyConductor-assisted</th></tr>"
        for name in result["primary"]["baseline"]:
            body += "<tr><td>" + escape(name) + "</td>" + "".join(
                "<td>" + escape(str(result["primary"][a][name])) + "</td>"
                for a in ("baseline", "intervention")) + "</tr>"
        body += "</table><ul>" + "".join("<li>" + escape(x) + "</li>" for x in result["limitations"]) + "</ul>"
        Path(args.html).parent.mkdir(parents=True, exist_ok=True)
        Path(args.html).write_text("<!doctype html><meta charset='utf-8'><title>Workflow evaluation</title>"
                                  "<style>body{font:16px system-ui;max-width:1000px;margin:32px auto}"
                                  "td,th{padding:10px;text-align:left;border-bottom:1px solid #ddd}</style>"
                                  + body, encoding="utf-8")
    print(f"Evaluated {len(rows)} rows; {result['n_independent_evaluation']} independent pairs. Wrote {args.out}")
    return 0
