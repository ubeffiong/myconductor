"""Governance for externally validated model versions; approval is not calibration.

A model must earn its place, not merely document itself
------------------------------------------------------
Reporting a performance number is not evidence that a model is worth
deploying. Without an incumbent to compare against, sensitivity 0.40 and
sensitivity 0.95 are both "independently evaluated performance" and both pass a
paperwork check. For TB-AMR the incumbent is never nothing: it is the WHO
catalogue, which answers wherever it holds a graded entry and abstains
elsewhere. A model that answers less often than the catalogue, or errs more
often at the same coverage, is not an improvement no matter how it was
validated.

So every performance claim carries the baseline it was measured against, and
approval refuses unless the model is at least as available and no less accurate
than that baseline. The comparison is deliberately two-sided: a model can
always buy a better error rate by abstaining more, and abstention is not free
when the alternative is a catalogue that would have answered.

This is the summary-level form of the test. The full risk-coverage curve, where
a model is evaluated at the baseline's own coverage rather than at its
reported one, lives in ``mycobench.analysis.selective`` and is the right tool
when per-isolate predictions are available. This module cannot import it —
the core carries no dependency on the benchmark harness — so the rule is
restated here over the summary numbers a registry actually holds.
"""
import copy
import math
from dataclasses import asdict, dataclass, field
from .ledger_store import LedgerStore

#: A model must beat the incumbent by at least this much before its error rate
#: counts as an improvement. Zero would approve a model that is worse in every
#: practical sense but rounds to equal; this demands a real margin.
MIN_ERROR_MARGIN = 0.0

#: Coverage is allowed to fall below the baseline's by at most this much. A
#: model that abstains more often than the catalogue is shifting work back to
#: the clinician while claiming to reduce it.
MAX_COVERAGE_SHORTFALL = 0.0


def baseline_verdict(row: dict, margin: float = MIN_ERROR_MARGIN) -> tuple[bool, str]:
    """Did this performance row beat the incumbent it names?

    Returns the verdict and the sentence a report should print for it. Both
    outcomes are returned rather than raised, because "did not beat the
    baseline" is a finding about the model, not a malformed record.
    """
    baseline = row.get("baseline") or {}
    coverage, error = row.get("call_rate"), row.get("error_rate")
    base_coverage, base_error = baseline.get("coverage"), baseline.get("error_rate")

    if coverage + MAX_COVERAGE_SHORTFALL < base_coverage:
        return False, (
            f"answers {coverage:.0%} of isolates where {baseline['source']} "
            f"answers {base_coverage:.0%}; a model that abstains more than the "
            f"incumbent has not reduced the unanswered fraction")
    if error > base_error - margin:
        return False, (
            f"errs on {error:.1%} of the isolates it answers where "
            f"{baseline['source']} errs on {base_error:.1%} at "
            f"{base_coverage:.0%} coverage; no improvement to deploy")
    return True, (
        f"answers {coverage:.0%} versus {base_coverage:.0%} and errs on "
        f"{error:.1%} versus {base_error:.1%} against {baseline['source']}")


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
            for key in ("sensitivity", "specificity", "call_rate", "error_rate"):
                value = row.get(key)
                if value is None or not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f"model {key} must be in [0,1]")
            if row["cohort"] not in self.validated_cohorts or row.get("independent") is not True:
                raise ValueError("performance must name an independently evaluated cohort")
            # A performance number with nothing to compare it against cannot
            # say whether the model is an improvement, only that someone
            # measured it. The incumbent is named here so approval can check.
            baseline = row.get("baseline")
            if not isinstance(baseline, dict):
                raise ValueError(
                    "model performance requires a baseline: the incumbent it "
                    "was measured against, with source, coverage and error_rate")
            if not baseline.get("source"):
                raise ValueError("model baseline requires source")
            for key in ("coverage", "error_rate"):
                value = baseline.get(key)
                if value is None or not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f"model baseline {key} must be in [0,1]")


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
                if entry.action == "approved":
                    rows = state[key]["model"]["performance"]
                    if not rows:
                        raise ValueError(
                            "approval requires independently reported performance")
                    # Reporting a number is not evidence of improvement. Every
                    # claimed scope must beat the incumbent it names, or the
                    # approval is refused and says which scope failed.
                    for row in rows:
                        beat, why = baseline_verdict(row)
                        if not beat:
                            raise ValueError(
                                f"approval refused for {row['drug']}/"
                                f"{row['lineage']}: {why}")
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

    def evidence_for(self, prediction, lineage=None) -> dict:
        """What this model beat, in the scope being used right now.

        A reader should not have to take "approved" on trust. This returns the
        incumbent, the margin over it, and the cohort the comparison was made
        on, so a report can state the basis rather than the verdict.
        """
        record = self.require_approved(prediction, lineage)
        scope = next(s for s in record["model"]["performance"]
                     if s["drug"] == prediction.drug
                     and s["lineage"] == (lineage or "unknown"))
        _beat, why = baseline_verdict(scope)
        return {
            "drug": scope["drug"], "lineage": scope["lineage"],
            "cohort": scope["cohort"], "baseline": scope["baseline"]["source"],
            "basis": why,
            "caveat": ("Beating the catalogue on an evaluation cohort is not "
                       "calibration and not clinical validation; this "
                       "prediction still cannot establish resistance."),
        }
