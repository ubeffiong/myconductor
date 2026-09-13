"""Governance for externally validated model versions; approval is not calibration."""
import copy
import math
from dataclasses import asdict, dataclass, field
from .ledger_store import LedgerStore


@dataclass(frozen=True)
class RegisteredModel:
    model_id: str
    version: str
    training_data_provenance: str
    organism: str
    validated_cohorts: list[str] = field(default_factory=list)
    performance: list[dict] = field(default_factory=list)

    def __post_init__(self):
        for name in ("model_id", "version", "training_data_provenance", "organism"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"model requires {name}")
        for row in self.performance:
            for key in ("drug", "lineage", "cohort", "source"):
                if not row.get(key):
                    raise ValueError(f"model performance requires {key}")
            for key in ("sensitivity", "specificity", "call_rate"):
                value = row.get(key)
                if value is None or not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f"model {key} must be in [0,1]")
            if row["cohort"] not in self.validated_cohorts or row.get("independent") is not True:
                raise ValueError("performance must name an independently evaluated cohort")


class ModelRegistry(LedgerStore):
    schema = "myconductor.model-registry.v1"

    def state(self):
        state = {}
        transitions = {"submitted": {"reviewed", "rejected"}, "reviewed": {"approved", "rejected"},
                       "approved": {"rolled_back"}, "rejected": set(), "rolled_back": set()}
        for entry in self.ledger.entries:
            p = entry.payload
            if not p.get("reviewer") or not p.get("rationale"):
                raise ValueError("model governance requires reviewer and rationale")
            key = (p["model_id"], p["version"])
            if entry.action == "submitted":
                if key in state:
                    raise ValueError("model version is immutable; submit a new version")
                model = RegisteredModel(**p["model"])
                if key != (model.model_id, model.version):
                    raise ValueError("registry model identity mismatch")
                state[key] = dict(model=asdict(model), status="submitted")
            else:
                if key not in state or entry.action not in transitions[state[key]["status"]]:
                    raise ValueError("invalid model governance transition")
                if entry.action == "approved" and not state[key]["model"]["performance"]:
                    raise ValueError("approval requires independently reported performance")
                state[key]["status"] = entry.action
            state[key]["reviewer"] = p["reviewer"]
        return copy.deepcopy(state)

    def submit(self, model, reviewer, rationale):
        self.append("submitted", dict(model_id=model.model_id, version=model.version,
                    model=asdict(model), reviewer=reviewer, rationale=rationale))

    def transition(self, model_id, version, action, reviewer, rationale):
        self.append(action, dict(model_id=model_id, version=version, reviewer=reviewer, rationale=rationale))

    def require_approved(self, prediction, lineage=None):
        record = self.state().get((prediction.model_id, prediction.model_version))
        if record is None or record["status"] != "approved":
            raise ValueError("in-silico prediction model/version is not approved")
        model = record["model"]
        if model["organism"] != prediction.organism or model["training_data_provenance"] != prediction.training_data_reference:
            raise ValueError("prediction differs from registered model provenance")
        scopes = [s for s in model["performance"] if s["drug"] == prediction.drug and s["lineage"] == (lineage or "unknown")]
        if not scopes:
            raise ValueError("model has no approved performance record for this drug and lineage")
        return record
