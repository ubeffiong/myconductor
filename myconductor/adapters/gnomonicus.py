"""Discrete gnomonicus JSON adapter, pinned to the upstream saveJSON schema.

Not a clinical validation. A sample binding is required because output schema
versions need not contain a specimen ID. Continuous MIC predictions are refused.
"""
import hashlib
import json
from pathlib import Path

from .base import EngineAdapter, EngineReport, AdapterSchemaError
from ..core.models import Call, DrugEvidence, EngineRef, Lane, Tier, Variant

DRUGS = {"RIF": "rifampicin", "INH": "isoniazid", "PZA": "pyrazinamide",
         "EMB": "ethambutol", "MXF": "moxifloxacin", "LEV": "levofloxacin",
         "BDQ": "bedaquiline", "CFZ": "clofazimine", "DLM": "delamanid",
         "LZD": "linezolid", "AMI": "amikacin", "PMD": "pretomanid"}
CALLS = {"R": Call.RESISTANT, "S": Call.SUSCEPTIBLE, "U": Call.INDETERMINATE,
         "F": Call.NO_CALL, "N": Call.NO_CALL}


class GnomonicusAdapter(EngineAdapter):
    def __init__(self, sample_id=None):
        self.sample_id = sample_id
        self.engine = EngineRef("gnomonicus")

    def parse(self, path):
        raw = Path(path).read_bytes()
        payload = json.loads(raw)
        meta, data = payload.get("meta"), payload.get("data")
        if not isinstance(meta, dict) or not isinstance(data, dict):
            raise AdapterSchemaError("gnomonicus requires meta and data objects")
        catalogue_types = meta.get("catalogue_type")
        catalogue_types = catalogue_types if isinstance(catalogue_types, list) else [catalogue_types]
        if (meta.get("workflow_name") != "gnomonicus"
                or meta.get("status") != "success"
                or not catalogue_types
                or any(not isinstance(t, str) or not t or not set(t) <= set(CALLS) for t in catalogue_types)):
            raise AdapterSchemaError("only gnomonicus discrete resistance predictions are supported")
        if not self.sample_id:
            raise AdapterSchemaError("gnomonicus requires an explicit sample binding")
        if (meta.get("sample_id") and meta["sample_id"] != self.sample_id
                or meta.get("guid") and meta["guid"] != self.sample_id):
            raise AdapterSchemaError("gnomonicus sample binding conflicts with metadata")
        if not isinstance(data.get("antibiogram"), dict):
            raise AdapterSchemaError("gnomonicus antibiogram must be an object")
        engine = EngineRef("gnomonicus", str(meta.get("workflow_version", "unknown")),
                           str(meta.get("catalogue_name", "unknown")),
                           str(meta.get("catalogue_version", "unknown")))
        report = EngineReport(engine, sample_id=self.sample_id,
                              source_sha256=hashlib.sha256(raw).hexdigest())
        if payload.get("errors"):
            from ..core.models import QCFinding
            report.qc.append(QCFinding("engine_mutation_errors", "fail",
                                       "gnomonicus reports mutation-processing errors"))
        assembly = meta.get("reference")
        if assembly not in ("NC_000962.3", "AL123456.3"):
            raise AdapterSchemaError("gnomonicus MTBC reference must be explicitly NC_000962.3 or AL123456.3")
        for name, code in data["antibiogram"].items():
            if code not in CALLS:
                raise AdapterSchemaError(f"unknown gnomonicus discrete call: {code!r}")
            drug = DRUGS.get(name, name.lower())
            report.evidence.append(DrugEvidence(
                drug, CALLS[code], Tier.CATALOGUED, Lane.ENGINE, engine=engine,
                sample_id=self.sample_id, catalogue_version=engine.database_version,
                rationale=f"gnomonicus reports {code} for {name}.",
                # The JSON antibiogram alone does not certify callable loci.
                asserts_coverage=False,
                metadata={"effects": data.get("effects", {}).get(name, [])},
            ))
        for mutation in data.get("mutations", []):
            if mutation.get("gene") and mutation.get("mutation"):
                report.variants.append(Variant.of(mutation["gene"], mutation["mutation"]))
        return report
