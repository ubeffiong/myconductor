"""Governed model-output display. No model execution or resistance-target ranking."""
import json
from dataclasses import asdict
from pathlib import Path
from ..core.models import Dimension, InSilicoPrediction, Tier


def load_in_silico(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("in-silico predictions must be a list")
    return [InSilicoPrediction(**r) for r in data]


def reconcile_in_silico_predictions(priorities, predictions, context, registry):
    if predictions and registry is None:
        raise ValueError("in-silico predictions require a model registry")
    findings, seen = [], set()
    for prediction in predictions:
        key = (prediction.variant_key, prediction.drug, prediction.model_id, prediction.model_version)
        if key in seen:
            raise ValueError("duplicate model prediction")
        seen.add(key)
        for name in ("sample_id", "isolate_id", "site_id", "organism"):
            if getattr(prediction, name) != getattr(context, name):
                raise ValueError(f"in-silico {name} differs from current context")
        # Not just "approved": the incumbent it beat, on which cohort, by how
        # much. A reader should be able to judge the basis, not the verdict.
        basis = registry.evidence_for(prediction, context.lineage)
        matching = [p for p in priorities if p.variant_key == prediction.variant_key and p.drug == prediction.drug]
        if not matching:
            raise ValueError("prediction does not match a current VUS/drug")
        finding = dict(asdict(prediction), tier=Tier.PREDICTED.value, call_effect="none", ranking_effect="none",
                       baseline_basis=basis,
                       interpretation="Externally supplied prediction and attributions; approval does not calibrate confidence or establish drug response.")
        findings.append(finding)
        for priority in matching:
            priority.dimensions.append(Dimension("in_silico_prediction", True, finding,
                source=prediction.source, note="Imported after ranking; confidence is not pooled with unrelated feature scales."))
    return findings
