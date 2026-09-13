"""Quantitative evidence imports; no MIC prediction model or breakpoint table."""
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Protocol

from ..core.models import Call, DrugEvidence, Lane, Tier

UNITS = {"mg/L", "ug/mL", "µg/mL", "μg/mL"}


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def mic_call(mic, mic_unit, critical_concentration, censoring="none",
             interval=None, susceptible_inclusive=None):
    """Compare in one explicit concentration unit. Equality defaults to unknown.

    left means <= bound, right means > bound. An interval is a closed interval.
    A prediction interval is propagated, never reduced to its central estimate.
    """
    _positive(mic, "MIC")
    _positive(critical_concentration, "critical concentration")
    if mic_unit not in UNITS:
        raise ValueError("MIC unit must be mg/L or ug/mL (including micro-symbol spellings)")
    if susceptible_inclusive not in (None, True, False):
        raise ValueError("susceptible_inclusive must be boolean or null")
    if censoring not in {"none", "left", "right", "interval"}:
        raise ValueError("unsupported MIC censoring")
    lower, upper = mic, mic
    if interval is not None:
        if len(interval) != 2:
            raise ValueError("MIC interval must have two bounds")
        lower, upper = interval
        _positive(lower, "MIC lower bound")
        _positive(upper, "MIC upper bound")
        if not lower <= mic <= upper:
            raise ValueError("MIC interval must contain the supplied value")
        if censoring in {"left", "right"}:
            raise ValueError("do not combine a one-sided censoring bound with an interval")
    elif censoring == "interval":
        raise ValueError("interval censoring requires both bounds")
    if censoring == "left":
        lower = 0
    if censoring == "right":
        upper = float("inf")
    cc = critical_concentration
    if lower > cc or (lower == cc and (censoring == "right" or susceptible_inclusive is False)):
        return Call.RESISTANT
    if upper < cc or (upper == cc and susceptible_inclusive is True):
        return Call.SUSCEPTIBLE
    return None


@dataclass(frozen=True)
class MICPrediction:
    prediction_id: str
    sample_id: str
    isolate_id: str
    site_id: str
    organism: str
    drug: str
    value: float
    unit: str
    method: str
    model_version: str
    source: str
    timestamp: str
    critical_concentration: float
    breakpoint_reference: str
    interval: Optional[tuple[float, float]] = None
    censoring: str = "none"
    susceptible_inclusive: Optional[bool] = None
    interval_kind: str = "unspecified"
    interval_level: Optional[float] = None

    def __post_init__(self):
        for name in ("prediction_id", "sample_id", "isolate_id", "site_id", "organism", "drug",
                     "method", "model_version", "source", "timestamp", "breakpoint_reference"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"MIC prediction requires {name}")
        if self.interval is not None:
            object.__setattr__(self, "interval", tuple(self.interval))
        if self.interval_kind not in {"unspecified", "prediction", "confidence", "measurement"}:
            raise ValueError("unknown MIC interval kind")
        if self.interval_level is not None and (self.interval is None or not math.isfinite(self.interval_level) or not 0 < self.interval_level < 1):
            raise ValueError("interval level requires bounds and a level in (0,1)")
        self.comparison()

    def comparison(self):
        return mic_call(self.value, self.unit, self.critical_concentration,
                        self.censoring, self.interval, self.susceptible_inclusive)


class MICPredictorProtocol(Protocol):
    def predict(self, context, evidence) -> list[MICPrediction]:
        ...


def load_predictions(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("MIC predictions must be a list")
    predictions = [MICPrediction(**row) for row in data]
    if len({p.prediction_id for p in predictions}) != len(predictions):
        raise ValueError("duplicate MIC prediction ID")
    return predictions


def reconcile_predictions(results, predictions, context):
    findings = []
    if len({p.prediction_id for p in predictions}) != len(predictions):
        raise ValueError("duplicate MIC prediction ID")
    for prediction in predictions:
        for key in ("sample_id", "isolate_id", "site_id", "organism"):
            if getattr(prediction, key) != getattr(context, key):
                raise ValueError(f"MIC prediction {key} differs from current context")
        result = next((r for r in results if r.drug == prediction.drug), None)
        if result is None or result.call is Call.UNSUPPORTED:
            raise ValueError("MIC prediction drug is outside the current result profile")
        call = prediction.comparison()
        conflict = call is not None and result.genomic_call in {Call.RESISTANT, Call.SUSCEPTIBLE} and call is not result.genomic_call
        finding = dict(asdict(prediction), comparison=call.value if call else "ambiguous",
                       genomic_call=result.genomic_call.value if result.genomic_call else None,
                       conflict=conflict, tier=Tier.PREDICTED.value,
                       interpretation="Supplied quantitative prediction; not a measured phenotype or validated categorical call.")
        findings.append(finding)
        if conflict:
            ev = DrugEvidence(result.drug, Call.INDETERMINATE, Tier.PREDICTED, Lane.ENGINE,
                              sample_id=context.sample_id, observation_id=prediction.prediction_id,
                              rationale="Supplied MIC prediction disagrees with genomic interpretation; laboratory review required.",
                              metadata={"mic_prediction": asdict(prediction)},
                              limitations=("No MIC model is implemented or validated by Myconductor.",))
            result.evidence.append(ev)
            # A supplied model does not overwrite an independently measured phenotype.
            if result.phenotypic_call is None or not result.phenotypic_call.is_established:
                result.call, result.tier = Call.INDETERMINATE, Tier.PREDICTED
                result.reason = ev.rationale
                result.confidence = None
    return findings
