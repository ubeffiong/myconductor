"""Compare frozen reports without attributing correlation to catalogue causation."""


def compare_reports(before, after):
    if before["sample_id"] != after["sample_id"]:
        raise ValueError("cannot compare different samples")
    a, b = before["provenance"], after["provenance"]
    if not a.get("input_sha256") or a["input_sha256"] != b.get("input_sha256"):
        raise ValueError("change impact requires identical variant input")
    changed_sources = sorted(k for k in set(a["source_hashes"]) | set(b["source_hashes"])
                             if a["source_hashes"].get(k) != b["source_hashes"].get(k))
    before_by_drug = {r["drug"]: r for r in before["drug_results"]}
    after_by_drug = {r["drug"]: r for r in after["drug_results"]}
    changes = []
    for drug in sorted(set(before_by_drug) | set(after_by_drug)):
        left, right = before_by_drug.get(drug), after_by_drug.get(drug)
        if left != right:
            changes.append({
                "drug": drug,
                "before": left["call"] if left else "absent",
                "after": right["call"] if right else "absent",
                "call_changed": left is None or right is None or left["call"] != right["call"],
                "before_evidence": left["evidence"] if left else [],
                "after_evidence": right["evidence"] if right else [],
            })
    return {
        "schema": "myconductor.change-impact.v1", "sample_id": before["sample_id"],
        "before_fingerprint": a["analysis_fingerprint"], "after_fingerprint": b["analysis_fingerprint"],
        "changed_sources": changed_sources,
        "context_changed": before.get("context") != after.get("context"),
        "policy_changed": before.get("analysis_manifest", {}).get("policy") != after.get("analysis_manifest", {}).get("policy"),
        "manifest_changes": sorted(k for k in set(before.get("analysis_manifest", {})) | set(after.get("analysis_manifest", {}))
                                   if before.get("analysis_manifest", {}).get(k) != after.get("analysis_manifest", {}).get(k)),
        "changes": changes,
        "interpretation": "Observed report differences. A catalogue version difference alone "
                          "does not prove the cause; inspect rule/evidence changes and replay "
                          "with all other inputs held fixed.",
        "requires_review": True,
    }
