"""Offline review report: escaped content, no external assets or scripts."""
from html import escape
from .json_report import report_dict


def render_html(report):
    d = report_dict(report)
    rows = []
    for r in d["drug_results"]:
        cells = [r["drug"], r["call"], r.get("genomic_call") or "unavailable",
                 r.get("phenotypic_call") or "not measured", r["assay_status"], r.get("reason") or ""]
        evidence = "".join("<li>" + escape(e["rationale"]) + " <small>" +
                           escape(e["evidence_id"][:16]) + "</small></li>" for e in r["evidence"])
        rows.append("<tr>" + "".join("<td>" + escape(str(c)) + "</td>" for c in cells)
                    + "</tr><tr><td colspan='6'><details><summary>Evidence</summary><ul>"
                    + evidence + "</ul></details></td></tr>")
    follow = {f["investigation_id"]: f for f in d["follow_up"]}
    findings = []
    for f in d["investigations"]:
        action = follow.get(f["id"], {})
        chosen = action.get("selected") or {}
        findings.append("<article><h3>" + escape(f["drug"] + ": " + f["category"]) +
                        "</h3><p>" + escape(f["certainty"] + " — " + f["explanation"]) +
                        "</p><p>Follow-up: " + escape(chosen.get("action", action.get("status", "unknown"))) +
                        "</p><p>Priority " + str(f["priority"]) +
                        "; estimated cost: " + escape(str(chosen.get("cost", "unknown"))) +
                        " " + escape(chosen.get("currency", "")) + "</p></article>")
    qc = "".join("<li>" + escape(q["status"] + " — " + q["check"] + ": " + q["detail"]) +
                 "</li>" for q in d["qc"])
    return """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>MyConductor evidence review</title>
<style>body{font:16px/1.55 system-ui,sans-serif;background:#f4f6f8;color:#172d3c;margin:0}
main{max-width:1250px;margin:auto;padding:32px}h1{font-size:32px}table{border-collapse:collapse;
width:100%;background:white}td,th{text-align:left;padding:12px;border-bottom:1px solid #ccd5dd}
th{background:#173b50;color:white}article{background:white;border-left:4px solid #25817c;
padding:12px 20px;margin:16px 0}.table{overflow:auto}small{color:#52616d}
.notice{padding:16px;background:#fff0d3}code{overflow-wrap:anywhere}</style><main>""" + (
        "<h1>MyConductor evidence review</h1><p>Sample: <strong>" + escape(d["sample_id"]) +
        "</strong> · Organism: " + escape(d["provenance"]["organism_profile"]) +
        "</p><p class='notice'>Research output. Priorities are review rules, not measured "
        "probabilities of resolution. No tests are ordered and no treatment is selected.</p>"
        "<h2>Interpretation</h2><div class='table'><table><thead><tr>" +
        "".join("<th>" + c + "</th>" for c in ("Drug", "Conclusion", "Genomic evidence",
                                               "Phenotype", "Assay adequacy", "Reason")) +
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"
        "<h2>Investigation and follow-up queue</h2>" +
        ("".join(findings) or "<p>No unresolved findings.</p>") +
        "<h2>Quality control</h2><ul>" + qc + "</ul><h2>Provenance</h2><p>Analysis: <code>" +
        escape(d["provenance"].get("analysis_fingerprint") or "unavailable") +
        "</code></p><p>Policy: " + escape(d["provenance"]["policy_version"]) +
        "</p></main></html>"
    )
