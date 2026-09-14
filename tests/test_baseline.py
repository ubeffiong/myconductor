"""Measuring the incumbent catalogue.

The load-bearing rule is the one for *susceptibility*. A catalogue that reports
"susceptible" because it found no resistant variant is reporting that it did not
look, which is the defect this codebase was rebuilt to remove. So the catalogue
is only allowed to answer ``S`` when the coordinates where resistance could have
hidden were actually examined — and ``assessed_variants`` carries that, drawn
from the VCF's reference calls rather than from the variant list.

Getting this wrong in the permissive direction would inflate the baseline's
coverage, which in turn lowers the bar a model must clear in
``ModelRegistry``. A wrong baseline is worse than no baseline: it launders a
weak model through the approval gate.
"""
from __future__ import annotations

import unittest

from mycobench.analysis import baseline
from mycobench.analysis.strata import DeterminantIndex, Isolate

DRUG = "rifampicin"
DETERMINANTS = {"rpoB_S450L", "rpoB_H445Y", "rpoB_D435V", "rpoB_S431T"}


def isolate(name, carries=(), assessed=None):
    return Isolate(isolate_id=name, genotype=frozenset(carries),
                   assessed_variants=(None if assessed is None
                                      else frozenset(assessed)))


class SusceptibilityRuleTests(unittest.TestCase):
    """When may the catalogue say "susceptible"?"""

    def test_a_carried_determinant_is_answered_resistant(self):
        predictions, _ = baseline.catalogue_predictions(
            [isolate("a", carries={"rpoB_S450L"})], {"a": "R"}, DETERMINANTS)
        self.assertEqual(predictions[0].predicted, "R")
        self.assertTrue(predictions[0].correct)

    def test_no_determinant_with_every_locus_examined_is_susceptible(self):
        predictions, abstained = baseline.catalogue_predictions(
            [isolate("a", assessed=DETERMINANTS)], {"a": "S"}, DETERMINANTS)
        self.assertEqual(predictions[0].predicted, "S")
        self.assertEqual(abstained, {})

    def test_no_determinant_and_nothing_examined_abstains(self):
        """The rule that matters. Not finding is not the same as not looking."""
        predictions, abstained = baseline.catalogue_predictions(
            [isolate("a", assessed=set())], {"a": "S"}, DETERMINANTS)
        self.assertIsNone(predictions[0].predicted)
        self.assertFalse(predictions[0].answerable)
        self.assertEqual(abstained[baseline.ABSTAIN_INCOMPLETE], 1)

    def test_partial_examination_below_the_threshold_abstains(self):
        # Two of four coordinates examined is 50%, under the 95% requirement.
        predictions, abstained = baseline.catalogue_predictions(
            [isolate("a", assessed={"rpoB_S450L", "rpoB_H445Y"})],
            {"a": "S"}, DETERMINANTS)
        self.assertIsNone(predictions[0].predicted)
        self.assertEqual(abstained[baseline.ABSTAIN_INCOMPLETE], 1)

    def test_unknown_assessed_set_is_treated_as_nothing_examined(self):
        """``None`` means the VCF carried no coverage evidence at all."""
        predictions, _ = baseline.catalogue_predictions(
            [isolate("a", assessed=None)], {"a": "S"}, DETERMINANTS)
        self.assertIsNone(predictions[0].predicted)

    def test_a_drug_with_no_graded_determinant_abstains_everywhere(self):
        predictions, abstained = baseline.catalogue_predictions(
            [isolate("a", assessed=DETERMINANTS)], {"a": "S"}, set())
        self.assertIsNone(predictions[0].predicted)
        self.assertEqual(abstained[baseline.ABSTAIN_NO_DETERMINANTS], 1)

    def test_an_isolate_without_a_phenotype_is_not_a_question(self):
        predictions, _ = baseline.catalogue_predictions(
            [isolate("a", assessed=DETERMINANTS)], {}, DETERMINANTS)
        self.assertEqual(predictions, [])

    def test_a_relaxed_threshold_can_be_requested_and_is_honoured(self):
        predictions, _ = baseline.catalogue_predictions(
            [isolate("a", assessed={"rpoB_S450L", "rpoB_H445Y", "rpoB_D435V"})],
            {"a": "S"}, DETERMINANTS, assessed_fraction=0.5)
        self.assertEqual(predictions[0].predicted, "S")


class MeasurementTests(unittest.TestCase):
    def _isolates(self, n_resistant=8, n_susceptible=8):
        isolates, truths = [], {}
        for i in range(n_resistant):
            isolates.append(isolate(f"r{i}", carries={"rpoB_S450L"},
                                    assessed=DETERMINANTS))
            truths[f"r{i}"] = "R"
        for i in range(n_susceptible):
            isolates.append(isolate(f"s{i}", assessed=DETERMINANTS))
            truths[f"s{i}"] = "S"
        return isolates, truths

    def test_a_perfect_catalogue_covers_everything_and_errs_on_nothing(self):
        isolates, truths = self._isolates()
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.coverage, 1.0)
        self.assertEqual(result.error_rate, 0.0)
        self.assertTrue(result.estimable)
        self.assertEqual(result.n_false_susceptible, 0)

    def test_false_susceptible_is_counted_separately(self):
        """The error that puts a patient on a failing drug is called out."""
        isolates, truths = self._isolates()
        # A resistant isolate whose mechanism the catalogue does not grade.
        isolates.append(isolate("missed", assessed=DETERMINANTS))
        truths["missed"] = "R"
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.n_errors, 1)
        self.assertEqual(result.n_false_susceptible, 1)
        self.assertIn("called susceptible but resistant", result.describe())

    def test_abstention_lowers_coverage_without_flattering_the_error_rate(self):
        isolates, truths = self._isolates()
        for i in range(10):
            isolates.append(isolate(f"u{i}", assessed=set()))
            truths[f"u{i}"] = "R"
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.n_answered, 16)
        self.assertEqual(result.n_isolates, 26)
        self.assertAlmostEqual(result.coverage, 16 / 26)
        # Declining the hard ones must not make the catalogue look accurate.
        self.assertEqual(result.error_rate, 0.0)
        self.assertEqual(result.abstained[baseline.ABSTAIN_INCOMPLETE], 10)

    def test_too_few_answers_is_not_estimable(self):
        isolates, truths = self._isolates(n_resistant=2, n_susceptible=2)
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertFalse(result.estimable)
        self.assertIn("no error rate is reportable", result.describe())
        self.assertTrue(any("describes the sample" in n for n in result.notes))

    def test_a_sample_without_resistance_says_nothing_about_sensitivity(self):
        isolates, truths = self._isolates(n_resistant=0, n_susceptible=20)
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertTrue(any("nothing about sensitivity" in n
                            for n in result.notes))


class RegistryHandoffTests(unittest.TestCase):
    """The measured point must be usable by the approval gate, or withheld."""

    def test_an_estimable_baseline_produces_a_registry_block(self):
        isolates, truths = [], {}
        for i in range(12):
            isolates.append(isolate(f"s{i}", assessed=DETERMINANTS))
            truths[f"s{i}"] = "S"
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        block = result.as_registry_baseline("WHO catalogue")
        self.assertEqual(block["source"], "WHO catalogue")
        self.assertEqual(block["coverage"], 1.0)
        self.assertEqual(block["error_rate"], 0.0)

    def test_an_unmeasurable_baseline_is_withheld_not_guessed(self):
        """A baseline nobody could measure must not be quoted as measured."""
        isolates = [isolate("s0", assessed=DETERMINANTS)]
        result = baseline.measure_drug(DRUG, isolates, {"s0": "S"}, DETERMINANTS)
        self.assertIsNone(result.as_registry_baseline("WHO catalogue"))

    def test_the_block_is_accepted_by_the_model_registry(self):
        """End to end: a measured baseline satisfies the approval gate."""
        from myconductor.federated.model_registry import baseline_verdict
        isolates, truths = [], {}
        for i in range(10):
            isolates.append(isolate(f"r{i}", carries={"rpoB_S450L"},
                                    assessed=DETERMINANTS))
            truths[f"r{i}"] = "R"
        for i in range(10):
            isolates.append(isolate(f"s{i}", assessed=DETERMINANTS))
            truths[f"s{i}"] = "S"
        measured = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        block = measured.as_registry_baseline("WHO catalogue")

        # A model matching the catalogue exactly does not beat it.
        row = dict(call_rate=block["coverage"], error_rate=block["error_rate"],
                   baseline=block)
        beat, _why = baseline_verdict(row)
        self.assertFalse(beat)

    def test_rows_render_for_a_report(self):
        isolates, truths = [], {}
        for i in range(12):
            isolates.append(isolate(f"s{i}", assessed=DETERMINANTS))
            truths[f"s{i}"] = "S"
        rows = baseline.baseline_rows(
            {DRUG: baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)})
        self.assertEqual(rows[0]["drug"], DRUG)
        self.assertEqual(rows[0]["estimable"], "yes")
        for column in baseline.BASELINE_COLUMNS:
            self.assertIn(column, rows[0])


class MeasureAcrossDrugsTests(unittest.TestCase):
    def test_binary_phenotypes_are_read_per_drug_code(self):
        index = DeterminantIndex(by_drug_label={"rifampicin": {"rpoB_S450L"}})
        isolates = [isolate("cr_ERR1", carries={"rpoB_S450L"},
                            assessed={"rpoB_S450L"})]
        rows = [{"ENA_RUN": "ERR1", "RIF_BINARY_PHENOTYPE": "R",
                 "INH_BINARY_PHENOTYPE": "NA"}]
        results = baseline.measure(isolates, rows, index,
                                   drugs=["rifampicin", "isoniazid"])
        self.assertEqual(results["rifampicin"].n_isolates, 1)
        # NA is not a phenotype, so isoniazid was never asked.
        self.assertEqual(results["isoniazid"].n_isolates, 0)

    def test_the_caveat_names_the_shared_isolate_problem(self):
        self.assertIn("share underlying isolates", baseline.CAVEAT)
        self.assertIn("not an independent accuracy estimate", baseline.CAVEAT)


if __name__ == "__main__":
    unittest.main()
