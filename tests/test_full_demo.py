import json
import re

from myconductor.full_demo import (
    DEMO_COHORT_SIZE,
    build_full_demo,
    build_full_demo_cohort,
    synthetic_report_context,
    write_full_demo,
)
from myconductor.reporting.html_report import render_html


def payload_of(html):
    match = re.search(r"window\.__MYCONDUCTOR__ = (.*?);</script>", html, re.S)
    assert match
    return json.loads(match.group(1).replace("\\u003c", "<"))


def test_full_demo_exercises_advanced_pipeline_workflows():
    report = build_full_demo()
    assert report.demo_mode
    assert report.population_structure is not None
    assert report.epistasis_notes
    assert report.quantitative_findings
    assert report.structural_annotations
    assert report.regulatory_findings
    assert report.expression_findings
    assert report.in_silico_findings
    assert report.panel
    assert report.discordances
    assert report.investigations
    assert report.follow_up


def test_full_demo_html_populates_every_advanced_view():
    reports = build_full_demo_cohort()
    html = render_html(reports[0], extra_reports=reports[1:], **{
        k: v for k, v in synthetic_report_context().items() if k != "schema"
    })
    payload = payload_of(html)
    assert len(payload["isolates"]) == DEMO_COHORT_SIZE
    assert len({row["id"] for row in payload["isolates"]}) == DEMO_COHORT_SIZE
    assert len(payload["drugs"]) == 10
    assert payload["validation_timeline"][0]["text"] != "Validation event"
    assert "undefined" not in json.dumps(payload)
    assert sum(site["n"] for site in payload["federated_sites"]) == 389
    for key in ("population", "mic", "structural", "expression", "models", "panels"):
        assert len(payload["workflows"][key]) == DEMO_COHORT_SIZE, key
    assert len(payload["workflows"]["regulatory"]) >= DEMO_COHORT_SIZE
    assert "NC_000962" not in json.dumps(payload["workflows"])
    assert "atpE" in json.dumps(payload["workflows"]["panels"])
    assert payload["run"]["scope"] == "cohort"
    assert payload["run"]["label"] == f"cohort-{DEMO_COHORT_SIZE}-isolates"
    assert "Report scope" in html
    assert "Cohort report" in html
    assert '<div class="k">Sample</div>' not in html
    rpo_b = payload["alignment_loci"]["rpoB"]
    assert len(rpo_b["positions"]) >= 3
    assert len(rpo_b["isolates"]) == DEMO_COHORT_SIZE
    assert all("carrier_fraction" in row for row in rpo_b["positions"])
    assert "Called variants in genomic context" in html
    assert "av-position-card" in html
    assert "prev-${field}" in html
    assert "aria-label" in html
    for key in ("prevalence", "targets", "epistasis", "watchlist", "audit",
                "validation_outcomes", "validation_timeline", "federated_sites",
                "error_trend", "lineages"):
        assert payload[key], key


def test_full_demo_writer_produces_all_review_artifacts(tmp_path):
    paths = write_full_demo(tmp_path)
    assert set(paths) == {"html", "json", "cohort_json", "fhir", "jsonld", "context"}
    assert all(path.is_file() and path.stat().st_size > 0 for path in paths.values())
    assert "SYNTHETIC" in paths["html"].read_text(encoding="utf-8").upper()
    assert len(json.loads(paths["cohort_json"].read_text(encoding="utf-8"))) == DEMO_COHORT_SIZE
    assert json.loads(paths["context"].read_text(encoding="utf-8"))["schema"] == "myconductor.report-context.v1"

