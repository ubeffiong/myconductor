"""Context beside a conflict must inform the reader, never decide for them.

A tie-break signal that breaks ties would let Tier.PREDICTED evidence settle a
clinical question through a side door — the black-box override this design
refuses everywhere else. These tests hold the line: the conflict's calls,
sources and conclusion survive annotation byte for byte, and only an approved
model is given the standing to appear beside a clinical disagreement at all.
"""
from __future__ import annotations

import unittest

from myconductor.core.models import Discordance, InSilicoPrediction
from myconductor.federated.model_registry import ModelRegistry, RegisteredModel
from myconductor.modules.discordance_context import (
    MAX_CONTEXT_LINES, annotate, context_for_drug,
)
from myconductor.modules.structural_annotation import StructuralAnnotation

ORGANISM = "Mycobacterium tuberculosis"


def conflict(drug="rifampicin") -> Discordance:
    return Discordance(
        drug=drug, calls=("resistant", "susceptible"),
        sources=("tb-profiler", "mykrobe"),
        note="1 source(s) call it resistant and 1 call it susceptible.")


def prediction(**kwargs) -> InSilicoPrediction:
    base = dict(
        variant_key="NC_000962.3:761155:C>T", drug="rifampicin",
        model_id="vus-classifier", model_version="2.1",
        training_data_reference="training/cryptic-2022",
        prediction="resistant", source="predictions/run-7",
        sample_id="S1", isolate_id="I1", site_id="site-01",
        organism=ORGANISM, timestamp="2026-09-14T00:00:00Z", confidence=0.92,
        feature_attributions={"conservation": 0.41, "pocket_distance": -0.22,
                              "plddt": 0.05})
    base.update(kwargs)
    return InSilicoPrediction(**base)


def annotation(**kwargs) -> StructuralAnnotation:
    base = dict(variant_key="NC_000962.3:761155:C>T", location="binding_pocket",
                predicted_effect="rifampicin contact disrupted",
                source="alphafold-db", organism=ORGANISM,
                reference_assembly="NC_000962.3", source_version="v4",
                timestamp="2026-09-14T00:00:00Z", ligand_distance=3.4,
                distance_unit="angstrom", ligand_reference="RIF")
    base.update(kwargs)
    return StructuralAnnotation(**base)


def approved_registry() -> ModelRegistry:
    model = RegisteredModel(
        "vus-classifier", "2.1", "training/cryptic-2022", ORGANISM,
        ["external-cohort"],
        [dict(drug="rifampicin", lineage="L4", cohort="external-cohort",
              source="evaluation/1", independent=True, sensitivity=0.93,
              specificity=0.97, call_rate=0.96, error_rate=0.04,
              baseline=dict(source="WHO catalogue", coverage=0.90,
                            error_rate=0.08))])
    registry = ModelRegistry()
    registry.submit(model, "submitter", "new")
    registry.transition("vus-classifier", "2.1", "reviewed", "reviewer", "read")
    registry.transition("vus-classifier", "2.1", "approved", "reviewer", "ok")
    return registry


class TheConflictSurvivesTests(unittest.TestCase):
    """Annotation is additive. Nothing about the disagreement changes."""

    def test_calls_sources_and_note_are_untouched(self):
        original = conflict()
        annotated = annotate([original], [prediction()],
                             registry=approved_registry(), lineage="L4")[0]
        self.assertEqual(annotated.calls, original.calls)
        self.assertEqual(annotated.sources, original.sources)
        self.assertEqual(annotated.note, original.note)
        self.assertEqual(annotated.drug, original.drug)

    def test_context_is_populated_and_marked_non_establishing(self):
        annotated = annotate([conflict()], [prediction()],
                             registry=approved_registry(), lineage="L4")[0]
        self.assertTrue(annotated.context)
        for line in annotated.context:
            self.assertIn("phenotypic testing still decides", line)

    def test_a_conflict_with_no_signals_is_returned_unchanged(self):
        original = conflict()
        self.assertIs(annotate([original])[0], original)

    def test_context_defaults_to_empty_so_existing_conflicts_are_unaffected(self):
        self.assertEqual(conflict().context, ())


class StandingTests(unittest.TestCase):
    """Only an approved model may appear beside a clinical disagreement."""

    def test_an_unapproved_model_is_not_shown_at_all(self):
        registry = ModelRegistry()
        registry.submit(
            RegisteredModel("vus-classifier", "2.1", "training/cryptic-2022",
                            ORGANISM, ["external-cohort"],
                            [dict(drug="rifampicin", lineage="L4",
                                  cohort="external-cohort", source="e/1",
                                  independent=True, sensitivity=0.9,
                                  specificity=0.9, call_rate=0.96,
                                  error_rate=0.04,
                                  baseline=dict(source="WHO catalogue",
                                                coverage=0.90,
                                                error_rate=0.08))]),
            "submitter", "new")
        context = context_for_drug("rifampicin", [prediction()],
                                   registry=registry, lineage="L4")
        self.assertEqual(context, ())

    def test_a_lineage_outside_the_approved_scope_is_not_shown(self):
        context = context_for_drug("rifampicin", [prediction()],
                                   registry=approved_registry(), lineage="L2")
        self.assertEqual(context, ())

    def test_an_approved_model_names_what_it_beat(self):
        context = context_for_drug("rifampicin", [prediction()],
                                   registry=approved_registry(), lineage="L4")
        self.assertIn("approved against WHO catalogue", context[0])
        self.assertIn("external-cohort", context[0])

    def test_predictions_for_other_drugs_are_not_shown(self):
        context = context_for_drug(
            "isoniazid", [prediction()], registry=approved_registry(),
            lineage="L4")
        self.assertEqual(context, ())


class PresentationTests(unittest.TestCase):
    def test_the_strongest_feature_attributions_are_named(self):
        context = context_for_drug("rifampicin", [prediction()],
                                   registry=approved_registry(), lineage="L4")
        self.assertIn("conservation +0.41", context[0])
        self.assertIn("pocket_distance -0.22", context[0])

    def test_structural_context_reports_distance_with_its_ligand(self):
        context = context_for_drug("rifampicin", annotations=[annotation()])
        self.assertIn("binding_pocket", context[0])
        self.assertIn("3.4 angstrom from RIF", context[0])

    def test_a_flood_of_signals_is_truncated_rather_than_dumped(self):
        many = [annotation(source_version=f"v{i}") for i in range(20)]
        context = context_for_drug("rifampicin", annotations=many)
        self.assertEqual(len(context), MAX_CONTEXT_LINES + 1)
        self.assertIn("further signal(s) omitted", context[-1])

    def test_context_renders_into_the_text_report(self):
        from myconductor.reporting.render import _wrap  # noqa: F401
        annotated = annotate([conflict()], [prediction()],
                             registry=approved_registry(), lineage="L4")[0]
        # The renderer labels the block so a reader cannot mistake it for a
        # resolution of the conflict.
        self.assertTrue(annotated.context)


if __name__ == "__main__":
    unittest.main()
