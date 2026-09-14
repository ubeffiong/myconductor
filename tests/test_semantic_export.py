"""The semantic export, and the one thing it must never let a consumer do.

Serialisation is the easy half. The risk is that every downstream consumer
wants two states and this system has six, so an export that offers any route to
a boolean will have that route taken — and a consumer mapping INDETERMINATE to
"not resistant", or dropping NOT_ASSESSED rows as empty, has rebuilt the
presumption this codebase removed, somewhere no guard in it can reach.

These tests hold three properties: no binary field exists at any depth, the
usability question is answered explicitly so nobody has to interpret the enum,
and the tier travels with the call so the basis is not lost.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from myconductor.core.models import (
    AnalysisReport, Call, Discordance, DrugEvidence, DrugResult, Lane,
    LocusCoverage, Tier,
)
from myconductor.reporting.semantic import (
    FORBIDDEN_TERMS, SCHEMA, VOCAB, audit, call_term, summarise, tier_term,
    to_jsonld, write_jsonld,
)


def result(drug="rifampicin", call=Call.RESISTANT, tier=Tier.CATALOGUED,
           **kwargs) -> DrugResult:
    base = dict(
        drug=drug, call=call, tier=tier,
        evidence=[DrugEvidence(drug, call, tier, Lane.CATALOGUE,
                               rationale="graded entry",
                               limitations=("illustrative catalogue",))],
        coverage=[LocusCoverage(locus="rpoB", mean_depth=120,
                                callable_fraction=0.99, source="depth-table")],
        reason="a graded catalogue entry fired")
    base.update(kwargs)
    return DrugResult(**base)


def report(results=None, **kwargs) -> AnalysisReport:
    base = dict(sample_id="S1",
                drug_results=results or [
                    result(),
                    result("isoniazid", Call.SUSCEPTIBLE, Tier.CATALOGUED),
                    result("bedaquiline", Call.NOT_ASSESSED, Tier.NONE),
                    result("linezolid", Call.INDETERMINATE, Tier.INFERRED),
                    result("pretomanid", Call.UNSUPPORTED, Tier.NONE),
                    result("amikacin", Call.NO_CALL, Tier.NONE),
                ],
                provenance=None)
    base.update(kwargs)
    return AnalysisReport(**base)


class NoBinaryFieldTests(unittest.TestCase):
    """The collapse hazard, guarded at every depth."""

    def test_no_boolean_is_named_for_a_call_state(self):
        document = to_jsonld(report())
        found: list[str] = []

        def walk(node, path):
            if isinstance(node, dict):
                for key, value in node.items():
                    if isinstance(value, bool) and key in FORBIDDEN_TERMS:
                        found.append(f"{path}.{key}")
                    walk(value, f"{path}.{key}")
            elif isinstance(node, list):
                for index, item in enumerate(node):
                    walk(item, f"{path}[{index}]")

        walk(document, "$")
        self.assertEqual(found, [])

    def test_the_audit_catches_an_injected_binary_field(self):
        document = to_jsonld(report())
        document["drugResults"][0]["resistant"] = True
        problems = audit(document)
        self.assertTrue(any("binary call field" in p for p in problems))

    def test_writing_a_collapsible_document_is_refused(self):
        class Sneaky(dict):
            pass

        document = to_jsonld(report())
        document["susceptible"] = False
        problems = audit(document)
        self.assertTrue(problems)

    def test_a_clean_document_audits_clean(self):
        self.assertEqual(audit(to_jsonld(report())), [])


class SemanticsTravelTests(unittest.TestCase):
    """A consumer must never have to interpret the enum itself."""

    def test_every_call_answers_the_usability_question_explicitly(self):
        for state in Call:
            term = call_term(state)
            self.assertIn("establishesUse", term)
            self.assertEqual(term["establishesUse"], state is Call.SUSCEPTIBLE)

    def test_only_susceptible_establishes_use(self):
        document = to_jsonld(report())
        usable = [d["drug"] for d in document["drugResults"]
                  if d["call"]["establishesUse"]]
        self.assertEqual(usable, ["isoniazid"])

    def test_the_four_unestablished_states_carry_distinct_reasons(self):
        reasons = {
            state.value: call_term(state)["reasonUnestablished"]
            for state in (Call.INDETERMINATE, Call.NOT_ASSESSED,
                          Call.NO_CALL, Call.UNSUPPORTED)
        }
        self.assertEqual(len(set(reasons.values())), 4)
        for value in reasons.values():
            self.assertTrue(value)

    def test_established_states_carry_no_reason(self):
        for state in (Call.RESISTANT, Call.SUSCEPTIBLE):
            self.assertIsNone(call_term(state)["reasonUnestablished"])

    def test_the_tier_travels_with_the_call(self):
        document = to_jsonld(report())
        for drug in document["drugResults"]:
            self.assertIn("tier", drug)
            self.assertIn("mayEstablishResistance", drug["tier"])

    def test_only_catalogued_and_phenotypic_may_establish_resistance(self):
        for tier in Tier:
            term = tier_term(tier)
            self.assertEqual(term["mayEstablishResistance"],
                             tier in (Tier.CATALOGUED, Tier.PHENOTYPIC))

    def test_the_document_enumerates_every_state_it_could_emit(self):
        """A consumer can build a complete mapping without reading our source."""
        document = to_jsonld(report())
        self.assertEqual(len(document["callStates"]), len(list(Call)))
        self.assertEqual(len(document["tiers"]), len(list(Tier)))


class DocumentShapeTests(unittest.TestCase):
    def test_a_context_gives_the_terms_stable_meaning(self):
        document = to_jsonld(report())
        self.assertIn("@context", document)
        self.assertEqual(document["@context"]["@vocab"], VOCAB)
        self.assertIn("must not be merged", document["@context"]["myco:note"])

    def test_the_schema_and_identity_are_stated(self):
        document = to_jsonld(report())
        self.assertEqual(document["schema"], SCHEMA)
        self.assertEqual(document["id"], "urn:myconductor:report:S1")

    def test_guidance_warns_against_deriving_a_boolean(self):
        guidance = to_jsonld(report())["interpretation"]
        self.assertIn("no resistant/susceptible boolean",
                      guidance["noBinaryField"])
        self.assertIn("None of them means susceptible", guidance["note"])

    def test_no_terminology_codes_are_invented(self):
        document = to_jsonld(report())
        self.assertIn("No SNOMED CT, LOINC or NCIT codes",
                      document["interpretation"]["terminology"])
        serialised = json.dumps(document)
        for system in ("snomed", "loinc", "ncit"):
            # Mentioned only in the disclaimer, never as an asserted code.
            self.assertNotIn(f'"{system}', serialised.lower())

    def test_discordance_is_exported_unresolved(self):
        conflict = Discordance(drug="rifampicin",
                               calls=("resistant", "susceptible"),
                               sources=("a", "b"), note="conflict",
                               context=("a model agreed with a",))
        document = to_jsonld(report([result(discordance=conflict)]))
        exported = document["drugResults"][0]["discordance"]
        self.assertFalse(exported["resolved"])
        self.assertEqual(exported["context"], ["a model agreed with a"])

    def test_evidence_limitations_survive_the_export(self):
        document = to_jsonld(report())
        evidence = document["drugResults"][0]["evidence"][0]
        self.assertIn("illustrative catalogue", evidence["limitations"])

    def test_unknown_coverage_is_distinguishable_from_zero(self):
        unknown = result(coverage=[LocusCoverage(locus="rpoB", source="absent")])
        document = to_jsonld(report([unknown]))
        exported = document["drugResults"][0]["coverage"][0]
        self.assertIsNone(exported["callableFraction"])
        self.assertFalse(exported["callableEvidencePresent"])

    def test_a_panel_is_exported_with_its_verification_status(self):
        document = to_jsonld(report(panel={
            "name": "First-line core", "version": "example-0",
            "assay": "targeted-amplicon", "verified": False,
            "unsupported_drugs": {"pyrazinamide": ["pncA"]}, "note": "n"}))
        panel = document["assayPanel"]
        self.assertFalse(panel["locusListVerified"])
        self.assertIn("pyrazinamide", panel["unsupportedDrugs"])

    def test_no_panel_means_no_panel_node(self):
        self.assertNotIn("assayPanel", to_jsonld(report()))


class WriteTests(unittest.TestCase):
    def test_the_document_round_trips_as_json(self):
        with tempfile.TemporaryDirectory() as root:
            path = write_jsonld(report(), Path(root) / "report.jsonld")
            loaded = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(loaded["sampleId"], "S1")
        self.assertIn("@context", loaded)

    def test_writing_refuses_a_document_that_fails_the_audit(self):
        import myconductor.reporting.semantic as semantic

        original = semantic.audit
        semantic.audit = lambda _doc: ["$.injected: binary call field"]
        try:
            with tempfile.TemporaryDirectory() as root:
                with self.assertRaises(ValueError) as caught:
                    write_jsonld(report(), Path(root) / "r.jsonld")
            self.assertIn("could collapse to two states", str(caught.exception))
        finally:
            semantic.audit = original

    def test_the_summary_counts_states_without_collapsing_them(self):
        text = summarise(to_jsonld(report()))
        self.assertIn("1 drug(s) with establishesUse true", text)
        for state in ("resistant", "susceptible", "not_assessed",
                      "indeterminate", "unsupported", "no_call"):
            self.assertIn(state, text)


if __name__ == "__main__":
    unittest.main()
