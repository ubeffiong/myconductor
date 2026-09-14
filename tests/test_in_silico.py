"""The machine-learning path, end to end — and where it is required to stop.

Myconductor does not train models, extract features or fold proteins. It
*governs* model output supplied from outside. That is a deliberate boundary,
and these tests hold it, because the value of the boundary is entirely in it
not moving:

* a model reaches the report only through an approved registry entry;
* approval requires beating the incumbent catalogue, not merely reporting a
  number (see ``test_model_registry``);
* and even then the prediction is ``Tier.PREDICTED``, which is structurally
  incapable of establishing resistance.

That last point is what makes the whole arrangement safe rather than merely
documented. The most capable model imaginable, validated on a million isolates
and approved by every reviewer, still cannot flip a drug to RESISTANT through
this path. It can only withhold susceptibility and raise a hypothesis for the
laboratory — which is the clinically protective direction.
"""
from __future__ import annotations

import unittest

from myconductor.core.context import SampleContext
from myconductor.core.models import (
    Call, DrugEvidence, InSilicoPrediction, Lane, Tier, VUSPriority,
)
from myconductor.federated.model_registry import ModelRegistry, RegisteredModel
from myconductor.modules.in_silico import reconcile_in_silico_predictions

ORGANISM = "Mycobacterium tuberculosis"


def context(**kwargs) -> SampleContext:
    base = dict(sample_id="S1", isolate_id="I1", site_id="site-01",
                organism=ORGANISM, lineage="L4")
    base.update(kwargs)
    return SampleContext(**base)


def performance(**kwargs) -> dict:
    row = dict(drug="rifampicin", lineage="L4", cohort="external-cohort",
               source="evaluation/1", independent=True,
               sensitivity=0.93, specificity=0.97,
               call_rate=0.96, error_rate=0.04,
               baseline=dict(source="WHO-UCN-TB-2023.7 catalogue",
                             coverage=0.90, error_rate=0.08))
    row.update(kwargs)
    return row


def approved_registry(rows=None) -> ModelRegistry:
    model = RegisteredModel("vus-classifier", "2.1", "training/cryptic-2022",
                            ORGANISM, ["external-cohort"],
                            rows or [performance()])
    registry = ModelRegistry()
    registry.submit(model, "submitter", "new version")
    registry.transition("vus-classifier", "2.1", "reviewed", "reviewer", "read")
    registry.transition("vus-classifier", "2.1", "approved", "reviewer", "beats catalogue")
    return registry


def prediction(**kwargs) -> InSilicoPrediction:
    base = dict(
        variant_key="NC_000962.3:761155:C>T", drug="rifampicin",
        model_id="vus-classifier", model_version="2.1",
        training_data_reference="training/cryptic-2022",
        prediction="resistant", source="predictions/run-7",
        sample_id="S1", isolate_id="I1", site_id="site-01",
        organism=ORGANISM, timestamp="2026-09-14T00:00:00Z",
        confidence=0.98,
        feature_attributions={"conservation": 0.41, "ligand_distance": -0.22},
    )
    base.update(kwargs)
    return InSilicoPrediction(**base)


def priority(**kwargs) -> VUSPriority:
    base = dict(variant_label="rpoB_p.Ser450Leu",
                variant_key="NC_000962.3:761155:C>T", gene="rpoB",
                drug="rifampicin", priority="high")
    base.update(kwargs)
    return VUSPriority(**base)


class BoundaryTests(unittest.TestCase):
    """The line myconductor does not cross."""

    def test_a_confident_prediction_cannot_establish_resistance(self):
        """The single most important assertion about this feature.

        The prediction below says "resistant" with confidence 0.98 from an
        approved model that beat the catalogue. It still resolves to
        Tier.PREDICTED, and that tier cannot carry a RESISTANT call at all.
        """
        self.assertFalse(Tier.PREDICTED.may_establish_resistance)
        with self.assertRaises(ValueError) as caught:
            DrugEvidence("rifampicin", Call.RESISTANT, Tier.PREDICTED, Lane.ENGINE)
        self.assertIn("only catalogued or phenotypic evidence",
                      str(caught.exception))

    def test_a_model_may_still_withhold_susceptibility(self):
        """Refusing to cross the line is not refusing to be useful.

        INDETERMINATE is permitted, and it is the clinically protective
        action: the drug does not reach a regimen on unestablished evidence.
        """
        evidence = DrugEvidence("rifampicin", Call.INDETERMINATE,
                                Tier.PREDICTED, Lane.ENGINE)
        self.assertIs(evidence.call, Call.INDETERMINATE)

    def test_predictions_neither_call_nor_rank(self):
        findings = reconcile_in_silico_predictions(
            [priority()], [prediction()], context(), approved_registry())
        self.assertEqual(findings[0]["call_effect"], "none")
        self.assertEqual(findings[0]["ranking_effect"], "none")
        self.assertEqual(findings[0]["tier"], Tier.PREDICTED.value)


class GovernanceTests(unittest.TestCase):
    def test_a_prediction_without_a_registry_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            reconcile_in_silico_predictions([priority()], [prediction()],
                                            context(), None)
        self.assertIn("require a model registry", str(caught.exception))

    def test_an_unapproved_model_cannot_reach_the_report(self):
        registry = ModelRegistry()
        registry.submit(
            RegisteredModel("vus-classifier", "2.1", "training/cryptic-2022",
                            ORGANISM, ["external-cohort"], [performance()]),
            "submitter", "new")
        with self.assertRaises(ValueError):
            reconcile_in_silico_predictions([priority()], [prediction()],
                                            context(), registry)

    def test_a_rolled_back_model_stops_being_usable(self):
        registry = approved_registry()
        registry.transition("vus-classifier", "2.1", "rolled_back",
                            "reviewer", "validation withdrawn")
        with self.assertRaises(ValueError):
            reconcile_in_silico_predictions([priority()], [prediction()],
                                            context(), registry)

    def test_a_lineage_outside_the_approved_scope_is_refused(self):
        """Approval is scoped. A model validated on L4 says nothing about L2."""
        with self.assertRaises(ValueError):
            reconcile_in_silico_predictions([priority()], [prediction()],
                                            context(lineage="L2"),
                                            approved_registry())

    def test_a_prediction_for_another_specimen_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            reconcile_in_silico_predictions(
                [priority()], [prediction(isolate_id="SOMEONE-ELSE")],
                context(), approved_registry())
        self.assertIn("differs from current context", str(caught.exception))

    def test_provenance_drift_from_the_registered_model_is_refused(self):
        with self.assertRaises(ValueError):
            reconcile_in_silico_predictions(
                [priority()], [prediction(training_data_reference="other")],
                context(), approved_registry())

    def test_duplicate_predictions_are_refused(self):
        with self.assertRaises(ValueError) as caught:
            reconcile_in_silico_predictions(
                [priority()], [prediction(), prediction()],
                context(), approved_registry())
        self.assertIn("duplicate", str(caught.exception))

    def test_a_prediction_for_a_variant_not_under_review_is_refused(self):
        with self.assertRaises(ValueError):
            reconcile_in_silico_predictions(
                [priority(variant_key="NC_000962.3:999999:A>G")],
                [prediction()], context(), approved_registry())


class ExplainabilityTests(unittest.TestCase):
    """Feature attributions — the SHAP layer — are carried, not computed."""

    def test_attributions_reach_the_finding_intact(self):
        findings = reconcile_in_silico_predictions(
            [priority()], [prediction()], context(), approved_registry())
        self.assertEqual(findings[0]["feature_attributions"],
                         {"conservation": 0.41, "ligand_distance": -0.22})

    def test_attributions_must_be_named_finite_values(self):
        for bad in ({"": 0.5}, {"conservation": float("nan")},
                    {"conservation": float("inf")}):
            with self.assertRaises(ValueError):
                prediction(feature_attributions=bad)

    def test_confidence_outside_zero_to_one_is_refused(self):
        for bad in (-0.1, 1.2, float("nan")):
            with self.assertRaises(ValueError):
                prediction(confidence=bad)

    def test_an_unknown_prediction_label_is_refused(self):
        with self.assertRaises(ValueError):
            prediction(prediction="probably-resistant")

    def test_the_finding_states_what_the_model_beat(self):
        """A reader should be able to judge the basis, not just the verdict."""
        findings = reconcile_in_silico_predictions(
            [priority()], [prediction()], context(), approved_registry())
        basis = findings[0]["baseline_basis"]
        self.assertEqual(basis["baseline"], "WHO-UCN-TB-2023.7 catalogue")
        self.assertEqual(basis["cohort"], "external-cohort")
        self.assertIn("cannot establish resistance", basis["caveat"])

    def test_the_dimension_attached_to_the_vus_does_not_pool_confidence(self):
        item = priority()
        reconcile_in_silico_predictions([item], [prediction()], context(),
                                        approved_registry())
        dimension = item.dimensions[-1]
        self.assertEqual(dimension.name, "in_silico_prediction")
        self.assertIn("not pooled", dimension.note)


if __name__ == "__main__":
    unittest.main()
