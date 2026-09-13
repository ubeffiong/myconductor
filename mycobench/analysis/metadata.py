"""Explicit cohort metadata and conservative removal of related replicates."""
import csv
from dataclasses import replace
from pathlib import Path


def apply_metadata(isolates, path, partition=None, independent_clusters=False):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    index, groups = {}, {}
    for row in rows:
        key = row.get("isolate_id")
        if not key or key in index:
            raise ValueError("metadata requires unique isolate_id")
        if row.get("partition") not in ("development", "evaluation", "", None):
            raise ValueError("unknown metadata partition")
        index[key] = row
        for field in ("cluster_id", "patient_id"):
            if row.get(field) and row.get("partition"):
                group = (field, row.get("site_id", "") if field == "patient_id" else "", row[field])
                if group in groups and groups[group] != row["partition"]:
                    raise ValueError(f"{field} crosses development/evaluation partitions")
                groups[group] = row["partition"]
    selected, seen_clusters, seen_patients = [], set(), set()
    for isolate in sorted(isolates, key=lambda i: i.isolate_id):
        row = index.get(isolate.isolate_id)
        if row is None:
            if partition or independent_clusters:
                continue
            selected.append(isolate)
            continue
        if partition and row.get("partition") != partition:
            continue
        if independent_clusters:
            if not all(row.get(k) for k in ("cluster_id", "patient_id", "site_id")):
                raise ValueError("independent-clusters requires cluster_id, patient_id and site_id")
            patient = (row["site_id"], row["patient_id"])
            if row["cluster_id"] in seen_clusters or patient in seen_patients:
                continue
            seen_clusters.add(row["cluster_id"])
            seen_patients.add(patient)
        selected.append(replace(isolate, lineage=row.get("lineage") or isolate.lineage,
                                site=row.get("site_id") or isolate.site))
    if not selected:
        raise ValueError("no genotyped isolates remain after metadata selection")
    return selected, [
        f"Metadata selection retained {len(selected)}/{len(isolates)} genotyped isolates.",
        "Cluster IDs are externally reviewed genetic relatedness groups, not inferred from lineage.",
        "Independent-clusters selects the first isolate ID per patient/cluster without looking at MIC or genotype."
        if independent_clusters else "No removal of related isolates requested.",
    ]
