"""Validated optional datasets for the self-contained HTML report.

The analysis report is the source of isolate-level calls.  This companion
document supplies cohort and programme context that cannot honestly be
inferred from one isolate.  Keeping it separate prevents a missing cohort
measurement from turning into a plausible-looking generated value.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Union


COLLECTIONS = {
    "baseline_rows", "lineage_rows", "measurability", "audit_entries",
    "federated_sites", "error_trend", "prevalence_rows", "target_rows",
    "watchlist_rows", "external_benchmark_rows",
}
ALLOWED = COLLECTIONS | {"validation", "reference_method"}


def load_report_context(path: Union[str, Path]) -> dict[str, Any]:
    """Load and validate a ``myconductor.report-context.v1`` document."""
    source = Path(path)
    data = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("report context must be a JSON object")
    schema = data.pop("schema", None)
    if schema != "myconductor.report-context.v1":
        raise ValueError("report context schema must be myconductor.report-context.v1")
    unknown = sorted(set(data) - ALLOWED)
    if unknown:
        raise ValueError("unknown report context field(s): " + ", ".join(unknown))
    for key in COLLECTIONS:
        if key in data and not isinstance(data[key], list):
            raise ValueError(f"report context field {key} must be an array")
        for index, row in enumerate(data.get(key, [])):
            if key == "error_trend":
                if isinstance(row, bool) or not isinstance(row, (int, float)):
                    raise ValueError(f"error_trend item {index} must be numeric")
                if row < 0 or row > 100:
                    raise ValueError(f"error_trend item {index} must be between 0 and 100")
            elif not isinstance(row, dict):
                raise ValueError(f"{key} item {index} must be an object")
    for index, row in enumerate(data.get("prevalence_rows", [])):
        total = row.get("total", row.get("n_tested"))
        resistant = row.get("resistant", row.get("n_resistant"))
        if total is not None and not _valid_non_negative_int(total):
            raise ValueError(f"prevalence_rows item {index} has invalid tested count")
        if resistant is not None and not _valid_non_negative_int(resistant):
            raise ValueError(f"prevalence_rows item {index} has invalid resistant count")
        if total is not None and resistant is not None and int(resistant) > int(total):
            raise ValueError(f"prevalence_rows item {index} has resistant count above tested count")
    for index, row in enumerate(data.get("target_rows", [])):
        if not (row.get("target") or row.get("gene")):
            raise ValueError(f"target_rows item {index} requires target or gene")
    if "validation" in data:
        validation = data["validation"]
        if not isinstance(validation, dict):
            raise ValueError("report context field validation must be an object")
        unknown_validation = sorted(set(validation) - {"outcomes", "timeline"})
        if unknown_validation:
            raise ValueError("unknown validation field(s): " +
                             ", ".join(unknown_validation))
        for key in ("outcomes", "timeline"):
            if key in validation and not isinstance(validation[key], list):
                raise ValueError(f"validation field {key} must be an array")
    if "reference_method" in data and not isinstance(data["reference_method"], str):
        raise ValueError("report context field reference_method must be text")
    return data


def _valid_non_negative_int(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    try:
        return int(value) >= 0
    except (TypeError, ValueError):
        return False
