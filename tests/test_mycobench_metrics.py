"""Metrics, with the abstention behaviour that makes them honest."""
import unittest

from mycobench.metrics import (
    score_accuracy,
    score_concordance,
    score_species_control,
    wilson_interval,
)
from mycobench.thresholds import (
    TARGETS,
    registration_hash,
    target_for,
)


class WilsonTests(unittest.TestCase):
    def test_zero_total_has_no_interval(self):
        self.assertIsNone(wilson_interval(0, 0))

    def test_bounds_stay_inside_zero_one_at_extremes(self):
        low, high = wilson_interval(0, 12)
        self.assertGreaterEqual(low, 0.0)
        self.assertLess(high, 1.0)
        low, high = wilson_interval(19, 19)
        self.assertGreater(low, 0.0)
        self.assertLessEqual(high, 1.0)

    def test_interval_narrows_as_n_grows(self):
        narrow = wilson_interval(500, 1000)
        wide = wilson_interval(5, 10)
        self.assertLess(narrow[1] - narrow[0], wide[1] - wide[0])


class AbstentionTests(unittest.TestCase):
    """Abstention must never be scored as susceptible."""

    def test_abstained_resistant_is_not_a_false_negative(self):
        result = score_accuracy(
            [("not_assessed", "R"), ("not_assessed", "R")], "rifampicin")
        self.assertEqual(result.false_negative, 0,
                         "an abstention is not a susceptible call")
        self.assertEqual(result.abstained_resistant, 2)
        self.assertEqual(result.n_called, 0)

    def test_vme_excludes_abstentions_from_its_denominator(self):
        # 1 real VME out of 2 called resistant-phenotype isolates; two more
        # were abstained on and must not dilute the rate.
        result = score_accuracy(
            [("susceptible", "R"), ("resistant", "R"),
             ("not_assessed", "R"), ("indeterminate", "R")], "rifampicin")
        self.assertEqual(result.n_phenotype_resistant, 2)
        self.assertAlmostEqual(result.vme_rate, 0.5)

    def test_call_rate_exposes_the_abstention(self):
        result = score_accuracy(
            [("resistant", "R")] + [("not_assessed", "S")] * 3, "rifampicin")
        self.assertAlmostEqual(result.call_rate, 0.25)

    def test_all_abstained_yields_no_accuracy_and_says_so(self):
        result = score_accuracy([("not_assessed", "R")] * 5, "rifampicin")
        self.assertEqual(result.n_called, 0)
        self.assertIsNone(result.sensitivity)
        self.assertTrue(any("no accuracy figure" in n for n in result.notes))

    def test_indeterminate_counts_as_abstention(self):
        result = score_accuracy([("indeterminate", "S")], "rifampicin")
        self.assertEqual(result.n_abstained, 1)
        self.assertEqual(result.true_negative, 0)

    def test_non_binary_phenotype_is_not_evaluable(self):
        result = score_accuracy(
            [("resistant", "U"), ("resistant", ""), ("resistant", "R")],
            "rifampicin")
        self.assertEqual(result.true_positive, 1)
        self.assertEqual(result.n_called, 1)


class ErrorClassTests(unittest.TestCase):
    def test_vme_is_predicted_susceptible_phenotype_resistant(self):
        result = score_accuracy([("susceptible", "R")], "rifampicin")
        self.assertEqual(result.false_negative, 1)
        self.assertAlmostEqual(result.vme_rate, 1.0)
        self.assertIsNone(result.me_rate, "no susceptible isolates to rate")

    def test_me_is_predicted_resistant_phenotype_susceptible(self):
        result = score_accuracy([("resistant", "S")], "rifampicin")
        self.assertEqual(result.false_positive, 1)
        self.assertAlmostEqual(result.me_rate, 1.0)

    def test_ppv_and_npv(self):
        result = score_accuracy(
            [("resistant", "R"), ("resistant", "S"),
             ("susceptible", "S"), ("susceptible", "R")], "rifampicin")
        self.assertAlmostEqual(result.ppv, 0.5)
        self.assertAlmostEqual(result.npv, 0.5)


class VerdictTests(unittest.TestCase):
    def test_underpowered_is_distinct_from_failure(self):
        result = score_accuracy([("resistant", "R")], "rifampicin")
        verdict, reasons = result.verdict()
        self.assertEqual(verdict, "underpowered")
        self.assertTrue(any("no pass or fail is claimed" in r for r in reasons))

    def test_clean_run_passes_its_bounds(self):
        target = target_for("rifampicin")
        pairs = ([("resistant", "R")] * 40) + ([("susceptible", "S")] * 40)
        result = score_accuracy(pairs, "rifampicin")
        self.assertGreaterEqual(result.n_evaluable, target.min_evaluable)
        verdict, _ = result.verdict()
        self.assertEqual(verdict, "pass")

    def test_excess_vme_fails_and_names_the_bound(self):
        pairs = ([("susceptible", "R")] * 20) + ([("resistant", "R")] * 20) \
            + ([("susceptible", "S")] * 40)
        result = score_accuracy(pairs, "rifampicin")
        verdict, reasons = result.verdict()
        self.assertEqual(verdict, "fail")
        self.assertTrue(any("VME" in r for r in reasons), reasons)

    def test_low_call_rate_fails_even_when_errors_are_zero(self):
        # Perfect on what it called, but it only called a third of them.
        pairs = ([("resistant", "R")] * 10) + ([("susceptible", "S")] * 10) \
            + ([("not_assessed", "S")] * 40)
        result = score_accuracy(pairs, "rifampicin")
        self.assertEqual(result.false_negative, 0)
        self.assertEqual(result.false_positive, 0)
        verdict, reasons = result.verdict()
        self.assertEqual(verdict, "fail")
        self.assertTrue(any("call rate" in r for r in reasons), reasons)

    def test_untargeted_drug_is_reported_not_judged(self):
        result = score_accuracy([("resistant", "R")] * 40, "streptomycin")
        verdict, _ = result.verdict()
        self.assertEqual(verdict, "not-targeted")

    def test_no_interval_when_underpowered(self):
        result = score_accuracy([("resistant", "R")] * 3, "rifampicin")
        self.assertIsNone(result.interval("sensitivity"))


class ThresholdTests(unittest.TestCase):
    def test_registration_hash_is_stable(self):
        first = registration_hash()
        self.assertEqual(first, registration_hash())
        self.assertEqual(len(first), 64)

    def test_error_ceilings_are_derived_not_declared_independently(self):
        # VME rate IS 1 - sensitivity and ME rate IS 1 - specificity. Declaring
        # them separately let the two disagree; deriving them cannot.
        for name, target in TARGETS.items():
            with self.subTest(drug=name):
                self.assertAlmostEqual(target.max_vme,
                                       1.0 - target.min_sensitivity)
                self.assertAlmostEqual(target.max_me,
                                       1.0 - target.min_specificity)

    def test_targets_declare_no_independent_error_fields(self):
        target = TARGETS["rifampicin"]
        self.assertNotIn("max_vme", target.__dataclass_fields__)
        self.assertNotIn("max_me", target.__dataclass_fields__)

    def test_rifampicin_and_isoniazid_carry_the_strictest_bounds(self):
        strictest = min(t.max_vme for t in TARGETS.values())
        self.assertAlmostEqual(TARGETS["rifampicin"].max_vme, strictest)
        self.assertAlmostEqual(TARGETS["isoniazid"].max_vme, strictest)

    def test_every_target_is_internally_sane(self):
        for name, target in TARGETS.items():
            with self.subTest(drug=name):
                self.assertTrue(target.rationale, f"{name} has no rationale")
                self.assertGreater(target.min_evaluable, 0)
                for value in (target.min_sensitivity, target.min_specificity,
                              target.min_call_rate):
                    self.assertGreater(value, 0.0)
                    self.assertLessEqual(value, 1.0)

    def test_verdict_does_not_double_report_one_breach(self):
        # A sensitivity breach and a VME breach are the same breach.
        pairs = ([("susceptible", "R")] * 20) + ([("resistant", "R")] * 20) \
            + ([("susceptible", "S")] * 40)
        _, reasons = score_accuracy(pairs, "rifampicin").verdict()
        breaches = [r for r in reasons if "is above" in r or "is below" in r]
        self.assertEqual(len(breaches), 1, reasons)


class ConcordanceTests(unittest.TestCase):
    def test_agreement_only_counts_isolates_both_engines_called(self):
        result = score_concordance([
            ("s1", "resistant", "resistant"),
            ("s2", "not_assessed", "resistant"),
            ("s3", "not_assessed", "not_assessed"),
        ], "rifampicin")
        self.assertEqual(result.agree, 1)
        self.assertEqual(result.only_b_called, 1)
        self.assertEqual(result.neither_called, 1)
        self.assertAlmostEqual(result.agreement_rate, 1.0)

    def test_mutual_abstention_is_not_agreement(self):
        result = score_concordance(
            [("s1", "not_assessed", "not_assessed")], "rifampicin")
        self.assertEqual(result.agree, 0)
        self.assertIsNone(result.agreement_rate)
        self.assertIn("nothing to reconcile", result.describe())

    def test_disagreement_records_the_pair(self):
        result = score_concordance(
            [("s1", "resistant", "susceptible")], "rifampicin")
        self.assertEqual(result.disagree, 1)
        self.assertEqual(result.disagreements,
                         [("s1", "resistant", "susceptible")])


class SpeciesControlTests(unittest.TestCase):
    def test_no_call_on_an_ntm_genome_passes(self):
        result = score_species_control(
            "ntm_1", "Mycobacterium avium",
            [("rifampicin", "unsupported"), ("isoniazid", "not_assessed")])
        self.assertTrue(result.passed)
        self.assertIn("as required", result.detail)

    def test_a_drug_call_on_an_ntm_genome_fails(self):
        result = score_species_control(
            "ntm_1", "Mycobacterium avium",
            [("rifampicin", "resistant"), ("isoniazid", "not_assessed")])
        self.assertFalse(result.passed)
        self.assertFalse(result.refused)
        self.assertEqual(result.offending_calls, [("rifampicin", "resistant")])

    def test_a_susceptible_call_also_fails(self):
        result = score_species_control(
            "ntm_1", "Mycobacterium abscessus", [("linezolid", "susceptible")])
        self.assertFalse(result.passed)


if __name__ == "__main__":
    unittest.main()
