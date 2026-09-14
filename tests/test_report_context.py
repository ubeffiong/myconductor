import json

import pytest

from myconductor.reporting.report_context import load_report_context


def write(tmp_path, value):
    path = tmp_path / "context.json"
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_valid_context_loads_all_discovery_collections(tmp_path):
    value = {
        "schema": "myconductor.report-context.v1",
        "prevalence_rows": [{"period": "2026", "total": 10, "resistant": 2}],
        "target_rows": [{"target": "DprE1"}],
        "watchlist_rows": [{"drug": "rifampicin"}],
        "validation": {"outcomes": [], "timeline": []},
        "reference_method": "MGIT 960",
    }
    loaded = load_report_context(write(tmp_path, value))
    assert loaded["prevalence_rows"][0]["total"] == 10
    assert loaded["target_rows"][0]["target"] == "DprE1"
    assert loaded["watchlist_rows"][0]["drug"] == "rifampicin"


@pytest.mark.parametrize("value,message", [
    ({}, "schema"),
    ({"schema": "wrong"}, "schema"),
    ({"schema": "myconductor.report-context.v1", "made_up": []}, "unknown"),
    ({"schema": "myconductor.report-context.v1", "target_rows": {}}, "array"),
    ({"schema": "myconductor.report-context.v1", "validation": []}, "object"),
])
def test_invalid_context_is_rejected(tmp_path, value, message):
    with pytest.raises(ValueError, match=message):
        load_report_context(write(tmp_path, value))
