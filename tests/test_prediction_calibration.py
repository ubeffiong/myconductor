"""Calibration raises a barrier. It never lowers one.

Two properties are load-bearing here.

**Direction.** There is no calibration good enough to turn a prediction into a
call. The tempting argument is that a low calibrated probability asserts the
absence of resistance and is therefore the safe direction — but SUSCEPTIBLE is
the only call that admits a drug to a regimen, so it is the actionable one. A
wrong RESISTANT costs a usable drug; a wrong SUSCEPTIBLE puts a patient on a
failing regimen. INDETERMINATE already withholds safely.

**Arithmetic.** "At least 30 paired observations" fails silently on the drugs
that matter most. Bedaquiline resistance runs near 0.8% in CRyPTIC, so thirty
pairs hold a quarter of one resistant isolate. The binding constraint is the
count of resistant observations, and these tests hold the code to that.
"""
from __future__ import annotations

import unittest

from myconductor.modules.prediction_calibration import (
    MAX_USABLE_ECE, MIN_POSITIVES, CalibrationResult, PairedObservation,
    PromotionRefused, calibrate, expected_calibration_error, priority_signal,
    refuse_promotion,
)

DRUG = "bedaquiline"
SITE = "site-01"


def pair(probability, resistant, i=0, drug=DRUG, site=SITE, lineage="L4"):
    return PairedObservation(f"v{i}", drug, probability, resistant, site,
                             lineage)


def cohort(n_resistant=12, n_susceptible=40, well_calibrated=True):
    """Paired observations whose probabilities do or do not track outcomes."""
    out = []
    for i in range(n_resistant):
        out.append(pair(0.9 if well_calibrated else 0.1, True, i))
    for i in range(n_susceptible):
        out.append(pair(0.05 if well_calibrated else 0.95, False,
                        1000 + i))
    return out


class DirectionTests(unittest.TestCase):
    """No calibration unlocks a call, in either direction."""

    def test_promotion_to_susceptible_is_refused_with_its_reasoning(self):
        result = calibrate(cohort(), DRUG, SITE)
        self.assertTrue(result.calibrated)
        with self.assertRaises(PromotionRefused) as caught:
            refuse_promotion("SUSCEPTIBLE", result)
        message = str(caught.exception)
        self.assertIn("admits a drug to a regimen", message)
        self.assertIn("callable-locus evidence", message)

    def test_promotion_to_resistant_is_refused(self):
        with self.assertRaises(PromotionRefused):
            refuse_promotion("RESISTANT", calibrate(cohort(), DRUG, SITE))

    def test_a_calibrated_result_advertises_that_it_grants_nothing(self):
        result = calibrate(cohort(), DRUG, SITE)
        self.assertFalse(result.may_establish_resistance)
        self.assertFalse(result.may_grant_susceptibility)

    def test_a_low_probability_produces_no_signal_at_all(self):
        """Not a reason to deprioritise.

        A model being unexcited about a variant is not evidence the variant is
        harmless. Letting it push work down the queue would be the model
        quietly deciding what never gets tested.
        """
        result = calibrate(cohort(), DRUG, SITE)
        self.assertIsNone(priority_signal(0.02, result))

    def test_a_high_probability_raises_priority_and_says_what_it_is_not(self):
        result = calibrate(cohort(), DRUG, SITE)
        signal = priority_signal(0.91, result)
        self.assertIn("prioritise phenotypic validation", signal)
        self.assertIn("does not establish resistance", signal)


class ArithmeticTests(unittest.TestCase):
    """Positives bind, not totals."""

    def test_many_pairs_with_too_few_resistant_is_refused(self):
        """The realistic bedaquiline case: plenty of data, no positives."""
        observations = cohort(n_resistant=2, n_susceptible=200)
        result = calibrate(observations, DRUG, SITE)
        self.assertFalse(result.calibrated)
        self.assertEqual(result.n_paired, 202)
        self.assertTrue(any("resistant observation" in r
                            for r in result.refusals))
        self.assertIsNone(result.ece)

    def test_the_refusal_says_how_many_more_are_needed(self):
        result = calibrate(cohort(n_resistant=3, n_susceptible=200), DRUG, SITE)
        self.assertTrue(any(f"{MIN_POSITIVES - 3} more" in r
                            for r in result.refusals))

    def test_an_all_resistant_set_is_refused_too(self):
        result = calibrate(cohort(n_resistant=40, n_susceptible=1), DRUG, SITE)
        self.assertTrue(any("susceptible observation" in r
                            for r in result.refusals))

    def test_too_few_pairs_overall_is_refused(self):
        result = calibrate(cohort(n_resistant=2, n_susceptible=3), DRUG, SITE)
        self.assertTrue(any("paired observation" in r for r in result.refusals))

    def test_a_sufficient_balanced_cohort_calibrates(self):
        result = calibrate(cohort(), DRUG, SITE)
        self.assertTrue(result.calibrated)
        self.assertIsNotNone(result.ece)
        self.assertTrue(result.usable_for_priority)


class ScopeTests(unittest.TestCase):
    def test_observations_from_another_site_are_not_borrowed(self):
        mixed = cohort() + [pair(0.9, True, 9000, site="site-02")]
        result = calibrate(mixed, DRUG, "site-02")
        self.assertFalse(result.calibrated)
        self.assertEqual(result.n_paired, 1)

    def test_observations_for_another_drug_are_not_borrowed(self):
        result = calibrate(cohort(), "linezolid", SITE)
        self.assertEqual(result.n_paired, 0)
        self.assertFalse(result.calibrated)

    def test_lineage_composition_is_reported(self):
        observations = (cohort()
                        + [pair(0.9, True, 5000, lineage="L2") for _ in range(1)])
        result = calibrate(observations, DRUG, SITE)
        self.assertIn("L4", result.lineages)
        self.assertIn("L2", result.lineages)

    def test_untyped_lineages_are_labelled_not_dropped(self):
        observations = [PairedObservation(f"v{i}", DRUG, 0.9, True, SITE, None)
                        for i in range(12)]
        observations += [PairedObservation(f"s{i}", DRUG, 0.05, False, SITE, None)
                         for i in range(40)]
        result = calibrate(observations, DRUG, SITE)
        self.assertEqual(result.lineages, {"untyped": 52})


class ReliabilityTests(unittest.TestCase):
    def test_a_well_calibrated_model_has_low_error(self):
        self.assertLess(expected_calibration_error(cohort()), 0.15)

    def test_a_badly_calibrated_model_has_high_error(self):
        error = expected_calibration_error(cohort(well_calibrated=False))
        self.assertGreater(error, MAX_USABLE_ECE)

    def test_a_badly_calibrated_model_cannot_order_the_queue(self):
        result = calibrate(cohort(well_calibrated=False), DRUG, SITE)
        self.assertTrue(result.calibrated)      # enough data to judge
        self.assertFalse(result.usable_for_priority)  # and the judgement is no
        self.assertIsNone(priority_signal(0.99, result))
        self.assertIn("not used for ordering", result.describe())

    def test_empty_input_returns_no_error_rather_than_zero(self):
        self.assertIsNone(expected_calibration_error([]))

    def test_probabilities_outside_zero_to_one_are_refused(self):
        for bad in (-0.01, 1.01, float("nan")):
            with self.assertRaises(ValueError):
                PairedObservation("v", DRUG, bad, True, SITE)

    def test_an_uncalibrated_result_never_signals(self):
        self.assertIsNone(priority_signal(0.99, CalibrationResult(DRUG, SITE)))


if __name__ == "__main__":
    unittest.main()
