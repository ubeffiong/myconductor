"""VUS classifier: predicts resistance for coding variants of unknown significance.

This is where a *catalogue lookup would give up* ("unknown") and where the real
value is. The production model is a supervised classifier (XGBoost or a
hierarchical attentive network) trained on CRyPTIC-style paired genome/DST data,
consuming engineered features:

    * evolutionary conservation   (SIFT / PolyPhen-2 style score)
    * structural impact           (AlphaFold binding-pocket / stability delta)
    * biochemical substitution    (e.g. Grantham distance)

To keep this scaffold dependency-free and runnable offline, ``predict`` here is
a transparent, deterministic scorer over the *same feature vector* the real
model would consume, and it emits a SHAP-style per-feature contribution so the
explainability contract is exercised end to end. Drop a trained estimator in via
``PredictorProtocol`` and nothing else changes.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Protocol

from ..core.models import Call, DrugEvidence, Route, Variant
from .base import VariantModule
from .features import FeatureVector, extract_features

# Map the gene a VUS sits in to the drug whose resistance it would most plausibly
# affect. A production system derives this from gene->drug knowledge graphs.
_GENE_DRUG = {
    "katG": "isoniazid", "inhA": "isoniazid", "fabG1": "isoniazid",
    "rpoB": "rifampicin", "rpoC": "rifampicin",
    "pncA": "pyrazinamide",
    "embB": "ethambutol", "embA": "ethambutol",
    "gyrA": "moxifloxacin", "gyrB": "moxifloxacin",
    "rrs": "amikacin", "eis": "amikacin",
    "atpE": "bedaquiline", "Rv0678": "bedaquiline", "pepQ": "bedaquiline",
    "ddn": "pretomanid", "fgd1": "pretomanid", "fbiA": "pretomanid",
    "rplC": "linezolid", "rrl": "linezolid",
}


class PredictorProtocol(Protocol):
    """Anything with ``score(features) -> (probability, contributions)`` fits."""

    def score(self, features: FeatureVector) -> tuple[float, dict[str, float]]:
        ...


@dataclass
class TransparentScorer:
    """A stand-in for the trained estimator.

    Linear-additive over normalised features so the SHAP-style contributions are
    exact rather than approximate -- fine for a scaffold, and it makes the
    explanation contract testable. Weights are hand-set, NOT learned.
    """

    weights = {
        "conservation": 0.42,
        "structural_impact": 0.34,
        "substitution_severity": 0.16,
        "in_resistance_region": 0.08,
    }
    bias = -0.45

    def score(self, features: FeatureVector) -> tuple[float, dict[str, float]]:
        contributions: dict[str, float] = {}
        z = self.bias
        for name, w in self.weights.items():
            x = getattr(features, name)
            contributions[name] = w * x
            z += w * x
        prob = 1.0 / (1.0 + pow(2.718281828, -6.0 * z))  # sigmoid, steepened
        return prob, contributions


class VUSClassifierModule(VariantModule):
    name = "vus_classifier"

    def __init__(self, predictor: Optional[PredictorProtocol] = None,
                 decision_threshold: float = 0.6,
                 abstain_band: float = 0.1):
        self.predictor = predictor or TransparentScorer()
        self.threshold = decision_threshold
        # If the probability lands within +/- abstain_band of the threshold the
        # model declines to call -- honest uncertainty beats a coin-flip.
        self.abstain_band = abstain_band
        self.model_name = type(self.predictor).__name__

    def evaluate(self, variant: Variant) -> Optional[DrugEvidence]:
        drug = _GENE_DRUG.get(variant.gene)
        if drug is None:
            return None  # not a resistance-associated gene we model

        features = extract_features(variant)
        prob, contributions = self.predictor.score(features)

        if abs(prob - self.threshold) <= self.abstain_band:
            call = Call.INDETERMINATE
            conf = 1.0 - abs(prob - self.threshold) / max(self.abstain_band, 1e-9)
        elif prob >= self.threshold:
            call, conf = Call.PREDICTED_RESISTANT, prob
        else:
            call, conf = Call.PREDICTED_SUSCEPTIBLE, 1.0 - prob

        top = sorted(contributions.items(), key=lambda kv: abs(kv[1]), reverse=True)
        explain = ", ".join(f"{k}={v:+.2f}" for k, v in top[:3])
        return DrugEvidence(
            drug=drug,
            call=call,
            confidence=round(conf, 3),
            route=Route.VUS_ML,
            variant_key=variant.key(),
            rationale=f"VUS predicted by {self.model_name} (p={prob:.2f}); drivers: {explain}.",
        )
