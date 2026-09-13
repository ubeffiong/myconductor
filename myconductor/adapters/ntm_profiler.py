"""NTM-Profiler adapter verified against the pinned upstream ERR459870 example.

Drug-class annotations remain uncertain at agent level. Assembly-based QC
depth values are not read-depth evidence. Reference states are preserved.
"""
from dataclasses import replace
import hashlib
import json
from pathlib import Path

from .base import EngineAdapter, EngineReport, AdapterSchemaError
from .tbprofiler import _variant_of
from ..core.models import Call, DrugEvidence, EngineRef, Lane, QCFinding, Tier


class NTMProfilerAdapter(EngineAdapter):
    def __init__(self, version="unknown"):
        self.engine = EngineRef("ntm-profiler", version, "ntmdb", "unknown")

    def parse(self, path):
        raw = Path(path).read_bytes()
        data = json.loads(raw)
        if not isinstance(data, dict) or not data.get("id"):
            raise AdapterSchemaError("NTM-Profiler output requires id")
        for key in ("dr_variants", "dr_genes", "other_variants"):
            if not isinstance(data.get(key), list):
                raise AdapterSchemaError(f"NTM-Profiler output requires list {key}")
        species = data.get("species", {}).get("species", [])
        if len(species) != 1 or not species[0].get("species"):
            raise AdapterSchemaError("NTM species assignment is absent or mixed")
        name = species[0]["species"]
        if "abscessus" not in name.lower():
            raise AdapterSchemaError("this adapter's organism-specific interpretation supports M. abscessus only")
        db = data.get("resistance_db") or {}
        engine = EngineRef("ntm-profiler", str(data.get("version", self.engine.version)),
                           "ntmdb", str(db.get("version", db.get("Date", "unknown"))))
        report = EngineReport(engine, sample_id=data["id"], species=name,
                              source_sha256=hashlib.sha256(raw).hexdigest())
        barcodes = data.get("barcode", [])
        subspecies = [b["id"].removeprefix("subsp. ") for b in barcodes
                      if str(b.get("id", "")).startswith("subsp. ") and b.get("frequency") == 1]
        if len(subspecies) == 1:
            report.context_findings["subspecies"] = subspecies[0]
        for item in data["dr_variants"] + data["other_variants"]:
            v = _variant_of(item)
            if v is None:
                raise AdapterSchemaError("unrecognised NTM variant record")
            v.identity = replace(v.identity, assembly="CU458896")
            if str(item.get("filter", "pass")).lower() != "pass":
                v.filters = (str(item["filter"]),)
            report.variants.append(v)
            for ann in item.get("drugs", []):
                drug = {"macrolides": "clarithromycin", "amikacin": "amikacin"}.get(ann.get("drug"))
                if not drug:
                    report.warnings.append(f"unmapped NTM drug annotation: {ann.get('drug')}")
                    continue
                report.evidence.append(DrugEvidence(
                    drug, Call.INDETERMINATE, Tier.INFERRED, Lane.ENGINE,
                    variant=v.identity, engine=engine, sample_id=data["id"],
                    rationale=f"NTM-Profiler annotates {v.label()} for {ann['drug']}; agent-level review required.",
                    scope="variant", metadata={"upstream_annotation": ann},
                    catalogue_version=engine.database_version,
                ))
                if drug == "clarithromycin" and v.gene == "rrl" and v.passed_filters and ann.get("type") == "drug_resistance":
                    report.context_findings.update(rrl_status="resistance_variant",
                                                  rrl_source=f"ntm-profiler:{report.source_sha256}")
        for gene in data["dr_genes"] + data.get("other_genes", []):
            if gene.get("gene_name") != "erm(41)" or str(gene.get("filter", "")).lower() != "pass":
                continue
            # A generic abnormal annotation does not establish loss of induction.
            status = {"functionally_normal": "functional"}.get(gene.get("type"))
            if status:
                report.context_findings.update(erm41_status=status,
                                              erm41_source=f"ntm-profiler:{report.source_sha256}")
        report.qc.append(QCFinding("species_confirmation", "pass", f"NTM-Profiler identifies {name}"))
        report.qc.append(QCFinding("ntm_reference_states", "warn",
                                  "Gene states imported; negative rrl/rrs calls and read callability are not inferred."))
        return report
