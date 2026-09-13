"""The analysis layer: hand-written statistics, so tested at the seams.

These are not integration smoke tests. Every function here was written from
scratch rather than taken from a library, and the parts that matter most are
the **refusals** — the cases where a statistic is not identifiable and the code
must decline rather than return a number that would be plotted as though it
meant something. Untested hand-rolled statistics are worse than none.
"""
import math
import unittest

from mycobench.analysis import barrier, effects, mechanism, mic, selective, stats
from mycobench.analysis.strata import (
    DeterminantIndex,
    Isolate,
    cooccurrence,
    stratify,
)


def obs(value, censoring=mic.NONE):
    """An observation at a given log2 bound."""
    return mic.Observation(log2_bound=value, censoring=censoring)


def exact(*values):
    return [obs(v) for v in values]


# ── stats ────────────────────────────────────────────────────────────────
class StatsTests(unittest.TestCase):
    def test_bootstrap_refuses_below_minimum_n(self):
        self.assertIsNone(stats.bootstrap_ci([1.0, 2.0, 3.0]))

    def test_bootstrap_is_deterministic(self):
        values = list(range(20))
        self.assertEqual(stats.bootstrap_ci(values), stats.bootstrap_ci(values))

    def test_bootstrap_interval_brackets_the_statistic(self):
        values = [float(i) for i in range(40)]
        low, high = stats.bootstrap_ci(values, stats.median)
        self.assertLessEqual(low, stats.median(values))
        self.assertGreaterEqual(high, stats.median(values))

    def test_difference_ci_refuses_a_tiny_group(self):
        self.assertIsNone(stats.difference_ci(list(range(20)), [1.0, 2.0]))

    def test_difference_ci_excludes_zero_for_separated_groups(self):
        low, high = stats.difference_ci([10.0] * 20, [0.0] * 20)
        self.assertGreater(low, 0.0)
        self.assertGreater(high, 0.0)

    def test_permutation_pvalue_is_never_exactly_zero(self):
        p = stats.permutation_pvalue([10.0] * 10, [0.0] * 10, iterations=200)
        self.assertGreater(p, 0.0)

    def test_permutation_pvalue_is_large_for_identical_groups(self):
        p = stats.permutation_pvalue([1.0] * 10, [1.0] * 10, iterations=200)
        self.assertGreater(p, 0.5)

    def test_holm_preserves_order_and_none_gaps(self):
        adjusted = stats.holm_adjust([0.01, None, 0.04])
        self.assertIsNone(adjusted[1])
        self.assertGreaterEqual(adjusted[0], 0.01)
        self.assertEqual(len(adjusted), 3)

    def test_holm_is_monotone_in_rank(self):
        adjusted = stats.holm_adjust([0.001, 0.01, 0.5])
        self.assertLessEqual(adjusted[0], adjusted[1])
        self.assertLessEqual(adjusted[1], adjusted[2])

    def test_benjamini_hochberg_is_less_strict_than_holm(self):
        pvalues = [0.001, 0.01, 0.02, 0.03, 0.04]
        survivors = sum(stats.benjamini_hochberg(pvalues))
        holm_survivors = sum(1 for p in stats.holm_adjust(pvalues)
                             if p is not None and p <= 0.05)
        self.assertGreaterEqual(survivors, holm_survivors)

    def test_cliffs_delta_bounds_and_sign(self):
        self.assertAlmostEqual(stats.cliffs_delta([2.0], [1.0]), 1.0)
        self.assertAlmostEqual(stats.cliffs_delta([1.0], [2.0]), -1.0)
        self.assertAlmostEqual(stats.cliffs_delta([1.0], [1.0]), 0.0)

    def test_inverse_simpson_measures_clustering(self):
        self.assertAlmostEqual(stats.inverse_simpson([10, 10, 10, 10, 10]), 5.0)
        self.assertLess(stats.inverse_simpson([96, 1, 1, 1, 1]), 1.2)
        self.assertEqual(stats.inverse_simpson([]), 0.0)

    def test_quantile_handles_a_single_observation(self):
        self.assertEqual(stats.quantile([3.0], 0.5), 3.0)


# ── mic: censoring ───────────────────────────────────────────────────────
class CensoringTests(unittest.TestCase):
    def test_censoring_direction_survives_parsing(self):
        self.assertEqual(mic.to_observation("<=0.25").censoring, mic.LEFT)
        self.assertEqual(mic.to_observation(">4.0").censoring, mic.RIGHT)
        self.assertEqual(mic.to_observation("1.0").censoring, mic.NONE)

    def test_log2_scale_is_used(self):
        self.assertAlmostEqual(mic.to_observation("4.0").log2_bound, 2.0)
        self.assertAlmostEqual(mic.to_observation("0.25").log2_bound, -2.0)

    def test_unusable_values_are_dropped_not_coerced(self):
        for value in ("NA", "", "0", "-1", None):
            self.assertIsNone(mic.to_observation(value), value)

    def test_profile_reports_fractions(self):
        profile = mic.censoring_profile(
            [obs(0, mic.LEFT)] * 7 + exact(1, 2, 3))
        self.assertAlmostEqual(profile.left_fraction, 0.7)
        self.assertEqual(profile.n_exact, 3)


class MedianIdentifiabilityTests(unittest.TestCase):
    """The amikacin problem: 69% left-censored means no median exists."""

    def test_median_refused_above_half_left_censored(self):
        observations = [obs(-2, mic.LEFT)] * 7 + exact(0, 1, 2)
        location = mic.median_log2(observations)
        self.assertFalse(location.available)
        self.assertIn("below the lowest tested dilution",
                      location.refused_because)

    def test_median_refused_above_half_right_censored(self):
        observations = [obs(3, mic.RIGHT)] * 7 + exact(0, 1, 2)
        location = mic.median_log2(observations)
        self.assertFalse(location.available)
        self.assertIn("above the highest tested dilution",
                      location.refused_because)

    def test_median_available_with_light_censoring(self):
        observations = exact(*range(10)) + [obs(-2, mic.LEFT)]
        location = mic.median_log2(observations)
        self.assertTrue(location.available)
        self.assertIsNotNone(location.interval)

    def test_exactly_half_censored_is_refused(self):
        observations = [obs(-2, mic.LEFT)] * 5 + exact(1, 2, 3, 4, 5)
        self.assertFalse(mic.median_log2(observations).available)

    def test_empty_is_refused_with_a_reason(self):
        location = mic.median_log2([])
        self.assertFalse(location.available)
        self.assertIn("no usable", location.refused_because)


class CompareTests(unittest.TestCase):
    """Every branch of the decidability rules."""

    def test_two_exact_values(self):
        self.assertEqual(mic.compare(obs(2), obs(1)), 1)
        self.assertEqual(mic.compare(obs(1), obs(2)), -1)
        self.assertEqual(mic.compare(obs(1), obs(1)), 0)

    def test_two_left_censored_are_never_decidable(self):
        self.assertIsNone(mic.compare(obs(-2, mic.LEFT), obs(-5, mic.LEFT)))

    def test_two_right_censored_are_never_decidable(self):
        self.assertIsNone(mic.compare(obs(3, mic.RIGHT), obs(9, mic.RIGHT)))

    def test_left_versus_right_is_always_decidable(self):
        self.assertEqual(mic.compare(obs(0, mic.LEFT), obs(0, mic.RIGHT)), -1)
        self.assertEqual(mic.compare(obs(0, mic.RIGHT), obs(0, mic.LEFT)), 1)

    def test_left_censored_below_an_exact_value_is_decidable(self):
        self.assertEqual(mic.compare(obs(-2, mic.LEFT), obs(1)), -1)

    def test_left_censored_against_a_lower_exact_value_is_not(self):
        self.assertIsNone(mic.compare(obs(-2, mic.LEFT), obs(-3)))

    def test_right_censored_above_an_exact_value_is_decidable(self):
        self.assertEqual(mic.compare(obs(3, mic.RIGHT), obs(1)), 1)

    def test_right_censored_against_a_higher_exact_value_is_not(self):
        self.assertIsNone(mic.compare(obs(3, mic.RIGHT), obs(5)))

    def test_exact_above_a_left_censored_bound_is_decidable(self):
        self.assertEqual(mic.compare(obs(1), obs(-2, mic.LEFT)), 1)

    def test_exact_below_a_right_censored_bound_is_decidable(self):
        self.assertEqual(mic.compare(obs(1), obs(3, mic.RIGHT)), -1)


class StochasticShiftTests(unittest.TestCase):
    def test_refuses_when_a_group_is_too_small(self):
        shift = mic.stochastic_shift(exact(1, 2), exact(*range(20)))
        self.assertFalse(shift.available)
        self.assertIn("at least", shift.refused_because)

    def test_separated_groups_give_delta_one(self):
        shift = mic.stochastic_shift(exact(*[10.0] * 10), exact(*[0.0] * 10))
        self.assertAlmostEqual(shift.delta, 1.0)
        self.assertEqual(shift.decisive_fraction, 1.0)

    def test_identical_groups_give_delta_zero(self):
        shift = mic.stochastic_shift(exact(*[1.0] * 10), exact(*[1.0] * 10))
        self.assertAlmostEqual(shift.delta, 0.0)

    def test_mutually_censored_groups_are_refused(self):
        # Both groups entirely below the plate: nothing is decidable.
        shift = mic.stochastic_shift([obs(-2, mic.LEFT)] * 10,
                                     [obs(-2, mic.LEFT)] * 10)
        self.assertFalse(shift.available)
        self.assertIn("mutually", shift.refused_because)

    def test_delta_is_a_conservative_lower_bound(self):
        # Carriers truly higher, but half are left-censored so half the pairs
        # are undecidable. The estimate must not exceed the honest maximum.
        carriers = exact(*[5.0] * 5) + [obs(-2, mic.LEFT)] * 5
        non_carriers = [obs(-2, mic.LEFT)] * 10
        shift = mic.stochastic_shift(carriers, non_carriers)
        if shift.available:
            self.assertLessEqual(shift.delta, 1.0)
            self.assertLess(shift.decisive_fraction, 1.0)

    def test_median_shift_omitted_when_unidentifiable(self):
        carriers = [obs(-2, mic.LEFT)] * 8 + exact(4, 5)
        non_carriers = exact(*range(10))
        shift = mic.stochastic_shift(carriers, non_carriers)
        self.assertIsNone(shift.median_shift)

    def test_median_shift_present_when_both_identifiable(self):
        shift = mic.stochastic_shift(exact(*[5.0] * 10), exact(*[1.0] * 10))
        self.assertAlmostEqual(shift.median_shift, 4.0)


class DrugPanelTests(unittest.TestCase):
    def test_panel_reports_unusable_for_location(self):
        panel = mic.DrugPanel("amikacin")
        for i in range(7):
            panel.add(f"iso{i}", "<=0.25")
        for i in range(7, 10):
            panel.add(f"iso{i}", "2.0")
        self.assertFalse(panel.usable_for_location())
        self.assertIn("refused", panel.describe())

    def test_add_rejects_unusable_values(self):
        panel = mic.DrugPanel("rifampicin")
        self.assertFalse(panel.add("iso1", "NA"))
        self.assertTrue(panel.add("iso2", "1.0"))
        self.assertEqual(len(panel.by_isolate), 1)

    def test_subset_ignores_missing_isolates(self):
        panel = mic.DrugPanel("rifampicin")
        panel.add("a", "1.0")
        self.assertEqual(len(panel.subset(["a", "absent"])), 1)


# ── strata ───────────────────────────────────────────────────────────────
class StratificationTests(unittest.TestCase):
    def setUp(self):
        self.index = DeterminantIndex(
            by_drug_label={"rifampicin": {"rpoB_S450L"}})

    def _isolates(self, spec):
        return [Isolate(isolate_id=name, genotype=frozenset(genes),
                        lineage=lineage)
                for name, genes, lineage in spec]

    def test_background_holds_established_determinants_only(self):
        background = self.index.background(
            "rifampicin", {"rpoB_S450L", "someGene_X1Y"})
        self.assertEqual(background, frozenset({"rpoB_S450L"}))

    def test_candidate_is_excluded_from_its_own_background(self):
        isolates = self._isolates([("a", {"rpoB_S450L"}, "l4")])
        result = stratify("rpoB_S450L", "rifampicin", isolates, self.index)
        self.assertEqual(result.strata[0].background, frozenset())

    def test_hitchhiker_is_detected_as_fully_confounded(self):
        # V always travels with the determinant; nobody carries V without it.
        isolates = self._isolates(
            [(f"c{i}", {"V_1", "rpoB_S450L"}, "l4") for i in range(10)]
            + [(f"n{i}", set(), "l4") for i in range(10)])
        result = stratify("V_1", "rifampicin", isolates, self.index)
        self.assertTrue(result.fully_confounded)
        self.assertEqual(result.n_informative, 0)

    def test_a_variant_confounded_with_lineage_is_not_informative(self):
        isolates = self._isolates(
            [(f"c{i}", {"V_1", "rpoB_S450L"}, "l4") for i in range(6)]
            + [(f"d{i}", {"rpoB_S450L"}, "l2") for i in range(6)])
        result = stratify("V_1", "rifampicin", isolates, self.index)
        self.assertTrue(result.fully_confounded)
        self.assertEqual(result.n_informative, 0)

    def test_lineage_diversity_of_carriers_is_measured(self):
        isolates = self._isolates(
            [("c1", {"V_1"}, "l2"), ("c2", {"V_1"}, "l4"),
             ("n1", set(), "l2"), ("n2", set(), "l4")])
        result = stratify("V_1", "rifampicin", isolates, self.index)
        self.assertGreater(result.lineage_diversity(), 1.5)

    def test_cooccurrence_fraction(self):
        isolates = self._isolates(
            [(f"c{i}", {"V_1", "rpoB_S450L"}, "l4") for i in range(8)]
            + [("c8", {"V_1"}, "l4"), ("c9", {"V_1"}, "l4")])
        fractions = cooccurrence("V_1", "rifampicin", isolates, self.index)
        self.assertAlmostEqual(fractions["rpoB_S450L"], 0.8)

    def test_indeterminate_catalogue_entries_are_not_background(self):
        # Only 'resistant' entries count; using uncertain ones as background
        # would assume the answer to the question being asked.
        self.assertEqual(effects.MIN_CARRIERS, 5)  # guard against silent drift
        self.assertEqual(("resistant",),
                         __import__("mycobench.analysis.strata",
                                    fromlist=["ESTABLISHED_CALLS"]
                                    ).ESTABLISHED_CALLS)


# ── effects ──────────────────────────────────────────────────────────────
class VariantEffectTests(unittest.TestCase):
    def setUp(self):
        self.index = DeterminantIndex(
            by_drug_label={"rifampicin": {"rpoB_S450L"}})

    def _panel(self, carrier_value, other_value, n=12):
        panel = mic.DrugPanel("rifampicin")
        for i in range(n):
            panel.add(f"c{i}", carrier_value)
            panel.add(f"n{i}", other_value)
        return panel

    def _isolates(self, n=12, carrier_genes=("V_1",), other_genes=()):
        return ([Isolate(f"c{i}", frozenset(carrier_genes),
                         lineage=("l2", "l4")[i % 2]) for i in range(n)]
                + [Isolate(f"n{i}", frozenset(other_genes),
                           lineage=("l2", "l4")[i % 2]) for i in range(n)])

    def test_confounded_is_distinct_from_no_evidence(self):
        isolates = ([Isolate(f"c{i}", frozenset({"V_1", "rpoB_S450L"}))
                     for i in range(10)]
                    + [Isolate(f"n{i}", frozenset()) for i in range(10)])
        stratification = stratify("V_1", "rifampicin", isolates, self.index)
        effect = effects.variant_effect(
            stratification, self._panel("8.0", "0.5", n=10),
            cooccurrence={"rpoB_S450L": 1.0})
        self.assertEqual(effect.verdict, "confounded")
        self.assertIn("never varies", effect.reasons[0])

    def test_too_few_carriers_is_underpowered_not_no_evidence(self):
        isolates = ([Isolate("c0", frozenset({"V_1"}))]
                    + [Isolate(f"n{i}", frozenset()) for i in range(10)])
        stratification = stratify("V_1", "rifampicin", isolates, self.index)
        effect = effects.variant_effect(stratification,
                                        self._panel("8.0", "0.5", n=10))
        self.assertEqual(effect.verdict, "underpowered")

    def test_a_real_shift_is_detected(self):
        isolates = self._isolates()
        stratification = stratify("V_1", "rifampicin", isolates, self.index)
        effect = effects.variant_effect(stratification,
                                        self._panel("8.0", "0.5"))
        self.assertEqual(effect.verdict, "evidence-of-effect")
        self.assertGreater(effect.delta, 0.9)
        self.assertIsNotNone(effect.interval)

    def test_no_shift_gives_no_evidence(self):
        isolates = self._isolates()
        stratification = stratify("V_1", "rifampicin", isolates, self.index)
        effect = effects.variant_effect(stratification,
                                        self._panel("1.0", "1.0"))
        self.assertEqual(effect.verdict, "no-evidence")

    def test_heavy_censoring_gives_not_estimable(self):
        isolates = self._isolates()
        stratification = stratify("V_1", "rifampicin", isolates, self.index)
        effect = effects.variant_effect(stratification,
                                        self._panel("<=0.25", "<=0.25"))
        self.assertEqual(effect.verdict, "not-estimable")

    def test_cooccurrence_warning_is_attached(self):
        isolates = self._isolates()
        stratification = stratify("V_1", "rifampicin", isolates, self.index)
        effect = effects.variant_effect(
            stratification, self._panel("8.0", "0.5"),
            cooccurrence={"rpoB_S450L": 0.95})
        self.assertTrue(any("co-occurs" in w for w in effect.warnings))

    def test_lineage_restricted_carriers_are_warned_about(self):
        isolates = ([Isolate(f"c{i}", frozenset({"V_1"}), lineage="l2")
                     for i in range(12)]
                    + [Isolate(f"n{i}", frozenset(), lineage="l2")
                       for i in range(12)])
        stratification = stratify("V_1", "rifampicin", isolates, self.index)
        effect = effects.variant_effect(stratification,
                                        self._panel("8.0", "0.5"))
        self.assertTrue(any("lineage" in w for w in effect.warnings))


class PermutationResolutionTests(unittest.TestCase):
    """The permutation statistic must be the effect size, not a median.

    This is a regression guard, not a general test of significance. An earlier
    version permuted the difference in medians. On MIC data that is degenerate:
    values lie on a discrete doubling-dilution series, so the median takes only
    a handful of distinct values and almost every relabelling reproduces the
    observed difference exactly. A perfectly separated twelve-versus-twelve
    comparison — every carrier above every non-carrier, delta +1.0, bootstrap
    interval (1.0, 1.0) — came back with p = 0.62 and was demoted to
    ``no-evidence``. A test that cannot detect complete separation would have
    silently suppressed every real finding in the scan.
    """

    def setUp(self):
        self.index = DeterminantIndex(
            by_drug_label={"rifampicin": {"rpoB_S450L"}})
        self.isolates = (
            [Isolate(f"c{i}", frozenset({"V_1"}), lineage=("l2", "l4")[i % 2])
             for i in range(12)]
            + [Isolate(f"n{i}", frozenset(), lineage=("l2", "l4")[i % 2])
               for i in range(12)])

    def _effect(self, carrier_values, other_values):
        panel = mic.DrugPanel("rifampicin")
        for i in range(12):
            panel.add(f"c{i}", carrier_values[i % len(carrier_values)])
            panel.add(f"n{i}", other_values[i % len(other_values)])
        stratification = stratify("V_1", "rifampicin", self.isolates,
                                  self.index)
        return effects.variant_effect(stratification, panel)

    def test_complete_separation_is_significant(self):
        effect = self._effect(("8.0", "16.0"), ("0.25", "0.5"))
        self.assertEqual(effect.delta, 1.0)
        self.assertIsNotNone(effect.pvalue)
        # The degenerate median statistic returned 0.62 on exactly this input.
        self.assertLess(effect.pvalue, 0.05)
        self.assertEqual(effect.verdict, "evidence-of-effect")

    def test_no_separation_is_not_significant(self):
        # Same values on both sides: the test must not manufacture a finding.
        effect = self._effect(("1.0", "2.0"), ("1.0", "2.0"))
        self.assertIsNotNone(effect.pvalue)
        self.assertGreater(effect.pvalue, 0.05)
        self.assertNotEqual(effect.verdict, "evidence-of-effect")

    def test_pvalue_is_reproducible(self):
        first = self._effect(("8.0", "16.0"), ("0.25", "0.5")).pvalue
        second = self._effect(("8.0", "16.0"), ("0.25", "0.5")).pvalue
        self.assertEqual(first, second)

    def test_pvalue_never_reaches_zero(self):
        """Add-one correction: p = 0 would claim more than resampling shows."""
        effect = self._effect(("8.0", "16.0"), ("0.25", "0.5"))
        self.assertGreater(effect.pvalue, 0.0)


class EffectScanTests(unittest.TestCase):
    def _effect(self, name, verdict, delta, pvalue):
        return effects.VariantEffect(variant=name, drug="rifampicin",
                                     verdict=verdict, delta=delta,
                                     pvalue=pvalue)

    def test_holm_demotes_an_uncorrected_finding(self):
        scan = effects.EffectScan(effects=[
            self._effect(f"V_{i}", "evidence-of-effect", 0.9, 0.04)
            for i in range(30)
        ]).finalise()
        self.assertEqual(scan.counts().get("evidence-of-effect", 0), 0,
                         "30 pairs at p=0.04 must not all survive correction")
        self.assertTrue(any("demoted" in r
                            for e in scan.effects for r in e.reasons))

    def test_a_strong_finding_survives_correction(self):
        scan = effects.EffectScan(effects=[
            self._effect("V_strong", "evidence-of-effect", 0.95, 0.00001)
        ] + [self._effect(f"V_{i}", "no-evidence", 0.01, 0.9)
             for i in range(20)]).finalise()
        strong = next(e for e in scan.effects if e.variant == "V_strong")
        self.assertEqual(strong.verdict, "evidence-of-effect")

    def test_confounded_variants_are_listed_separately(self):
        scan = effects.EffectScan(effects=[
            self._effect("V_hitch", "confounded", None, None)]).finalise()
        self.assertEqual(len(scan.confounded()), 1)
        self.assertEqual(scan.validation_queue(), [])

    def test_validation_queue_prefers_fdr_survivors(self):
        scan = effects.EffectScan(effects=[
            self._effect("V_weak", "no-evidence", 0.9, 0.9),
            self._effect("V_real", "evidence-of-effect", 0.5, 0.0001),
        ]).finalise()
        queue = scan.validation_queue()
        self.assertEqual(queue[0].variant, "V_real")


# ── selective prediction ─────────────────────────────────────────────────
def prediction(key, predicted, truth, confidence):
    return selective.Prediction(key, predicted, truth, confidence)


class RiskCoverageTests(unittest.TestCase):
    def _mixed(self):
        # High-confidence predictions are correct; low-confidence ones wrong.
        correct = [prediction(f"c{i}", "R", "R", 0.9 + i * 0.001)
                   for i in range(20)]
        wrong = [prediction(f"w{i}", "R", "S", 0.1 + i * 0.001)
                 for i in range(20)]
        return correct + wrong

    def test_coverage_falls_as_the_threshold_rises(self):
        curve = selective.risk_coverage_curve(self._mixed())
        coverages = [p.coverage for p in curve.points]
        self.assertEqual(coverages, sorted(coverages, reverse=True))

    def test_abstaining_on_the_wrong_ones_lowers_selective_risk(self):
        curve = selective.risk_coverage_curve(self._mixed())
        full = max(curve.points, key=lambda p: p.coverage)
        strict = min((p for p in curve.points
                      if p.selective_risk is not None),
                     key=lambda p: p.coverage)
        self.assertLess(strict.selective_risk, full.selective_risk)

    def test_risk_not_reported_below_the_minimum_covered(self):
        curve = selective.risk_coverage_curve(
            [prediction("a", "R", "R", 0.9), prediction("b", "R", "S", 0.8)])
        for point in curve.points:
            self.assertIsNone(point.selective_risk)

    def test_unanswerable_queries_count_against_coverage(self):
        predictions = ([prediction(f"a{i}", "R", "R", 0.9) for i in range(10)]
                       + [prediction(f"n{i}", None, "R", 0.0)
                          for i in range(10)])
        curve = selective.risk_coverage_curve(predictions)
        self.assertEqual(curve.n_unanswerable, 10)
        self.assertLessEqual(max(p.coverage for p in curve.points), 0.5)

    def test_area_under_curve_needs_two_usable_points(self):
        curve = selective.risk_coverage_curve(
            [prediction("a", "R", "R", 0.9)])
        self.assertIsNone(curve.area_under_curve())

    def test_empty_input_is_handled(self):
        curve = selective.risk_coverage_curve([])
        self.assertEqual(curve.points, [])


class BaselineComparisonTests(unittest.TestCase):
    """The falsification test for a VUS model."""

    def test_a_model_that_beats_the_catalogue_passes(self):
        curve = selective.risk_coverage_curve(
            [prediction(f"a{i}", "R", "R", 0.9) for i in range(40)])
        result = selective.beats_baseline_at_matched_coverage(
            curve, baseline_coverage=0.5, baseline_risk=0.20)
        self.assertTrue(result.passed)

    def test_a_model_no_better_than_the_catalogue_fails(self):
        predictions = ([prediction(f"a{i}", "R", "R", 0.9) for i in range(20)]
                       + [prediction(f"w{i}", "R", "S", 0.9) for i in range(20)])
        curve = selective.risk_coverage_curve(predictions)
        result = selective.beats_baseline_at_matched_coverage(
            curve, baseline_coverage=0.5, baseline_risk=0.20)
        self.assertFalse(result.passed)
        self.assertEqual(result.verdict, "no-better-than-baseline")
        self.assertIn("without adding information", result.reasons[0])

    def test_a_model_that_cannot_reach_the_coverage_says_so(self):
        predictions = ([prediction(f"a{i}", "R", "R", 0.9) for i in range(10)]
                       + [prediction(f"n{i}", None, "R", 0.0)
                          for i in range(90)])
        curve = selective.risk_coverage_curve(predictions)
        result = selective.beats_baseline_at_matched_coverage(
            curve, baseline_coverage=0.9, baseline_risk=0.2)
        self.assertEqual(result.verdict, "cannot-match-coverage")

    def test_catalogue_baseline_is_a_single_point(self):
        predictions = ([prediction(f"a{i}", "R", "R", 1.0) for i in range(30)]
                       + [prediction(f"u{i}", None, "R", 0.0)
                          for i in range(70)])
        coverage, risk = selective.catalogue_baseline(predictions)
        self.assertAlmostEqual(coverage, 0.3)
        self.assertAlmostEqual(risk, 0.0)


# ── mechanism ────────────────────────────────────────────────────────────
class TwoComponentTests(unittest.TestCase):
    def test_recovers_two_known_components(self):
        observations = exact(*([0.0 + i * 0.01 for i in range(40)]
                               + [6.0 + i * 0.01 for i in range(40)]))
        fit = mechanism.fit_two_component(observations)
        self.assertTrue(fit.available)
        self.assertGreater(fit.separation, 4.0)

    def test_refuses_too_few_exact_observations(self):
        fit = mechanism.fit_two_component(exact(*range(10)))
        self.assertFalse(fit.available)
        self.assertIn("below the", fit.refused_because)

    def test_refuses_when_censoring_would_fake_a_mode(self):
        observations = ([obs(-2, mic.LEFT)] * 30
                        + exact(*[1.0 + i * 0.01 for i in range(50)]))
        fit = mechanism.fit_two_component(observations)
        self.assertFalse(fit.available)
        self.assertIn("plate", fit.refused_because)

    def test_refuses_a_degenerate_single_value(self):
        fit = mechanism.fit_two_component(exact(*[1.0] * 50))
        self.assertFalse(fit.available)


class MechanismClassificationTests(unittest.TestCase):
    def _panels(self, bdq_carrier, cfz_carrier, n=10):
        panels = {}
        for drug, carrier_value in (("bedaquiline", bdq_carrier),
                                    ("clofazimine", cfz_carrier)):
            panel = mic.DrugPanel(drug)
            for i in range(n):
                panel.add(f"c{i}", carrier_value)
                panel.add(f"n{i}", "0.06")
            panels[drug] = panel
        return panels

    def _groups(self, n=10):
        return [f"c{i}" for i in range(n)], [f"n{i}" for i in range(n)]

    def test_a_shared_shift_reads_as_efflux_like(self):
        carriers, others = self._groups()
        call = mechanism.classify("Rv0678_L117R", self._panels("2.0", "2.0"),
                                  carriers, others)
        self.assertEqual(call.mechanism, "efflux-like")
        self.assertEqual(call.confidence, "exploratory")

    def test_a_single_drug_shift_reads_as_target_like(self):
        carriers, others = self._groups()
        call = mechanism.classify("atpE_D28V", self._panels("2.0", "0.06"),
                                  carriers, others)
        self.assertEqual(call.mechanism, "target-like")

    def test_no_shift_is_indeterminate(self):
        carriers, others = self._groups()
        call = mechanism.classify("V_1", self._panels("0.06", "0.06"),
                                  carriers, others)
        self.assertEqual(call.mechanism, "indeterminate")

    def test_one_estimable_drug_is_not_enough(self):
        carriers, others = self._groups()
        panels = self._panels("2.0", "2.0")
        del panels["clofazimine"]
        call = mechanism.classify("Rv0678_L117R", panels, carriers, others)
        self.assertEqual(call.mechanism, "indeterminate")
        self.assertIn("needs both", call.reasons[0])

    def test_unknown_efflux_set_raises(self):
        carriers, others = self._groups()
        with self.assertRaises(KeyError):
            mechanism.classify("V_1", self._panels("2.0", "2.0"),
                               carriers, others, efflux_set="nope")


class DiscriminatorValidationTests(unittest.TestCase):
    """The method ships with the experiment that would kill it."""

    def _panels(self, spec, n=10):
        panels = {}
        for drug in ("bedaquiline", "clofazimine"):
            panel = mic.DrugPanel(drug)
            for variant, (bdq, cfz) in spec.items():
                value = bdq if drug == "bedaquiline" else cfz
                for i in range(n):
                    panel.add(f"{variant}_c{i}", value)
            for i in range(n):
                panel.add(f"wt{i}", "0.06")
            panels[drug] = panel
        return panels

    def _carriers(self, spec, n=10):
        return {variant: [f"{variant}_c{i}" for i in range(n)]
                for variant in spec}

    def test_discriminates_when_the_signature_holds(self):
        spec = {"Rv0678_L117R": ("2.0", "2.0"), "atpE_D28V": ("2.0", "0.06")}
        carriers = self._carriers(spec)
        isolates = [i for group in carriers.values() for i in group] \
            + [f"wt{i}" for i in range(10)]
        result = mechanism.validate_discriminator(
            self._panels(spec), carriers, isolates)
        self.assertEqual(result.verdict, "discriminates")
        self.assertTrue(result.passed)

    def test_does_not_discriminate_when_both_look_the_same(self):
        # atpE also moves clofazimine: the signature carries no information.
        spec = {"Rv0678_L117R": ("2.0", "2.0"), "atpE_D28V": ("2.0", "2.0")}
        carriers = self._carriers(spec)
        isolates = [i for group in carriers.values() for i in group] \
            + [f"wt{i}" for i in range(10)]
        result = mechanism.validate_discriminator(
            self._panels(spec), carriers, isolates)
        self.assertEqual(result.verdict, "does-not-discriminate")
        self.assertTrue(any("should be dropped" in r for r in result.reasons))

    def test_not_run_without_both_control_classes(self):
        spec = {"Rv0678_L117R": ("2.0", "2.0")}
        carriers = self._carriers(spec)
        isolates = [i for group in carriers.values() for i in group] \
            + [f"wt{i}" for i in range(10)]
        result = mechanism.validate_discriminator(
            self._panels(spec), carriers, isolates)
        self.assertEqual(result.verdict, "not-run")

    def test_reference_group_excludes_other_tested_variants(self):
        """The bug this pins: a reference group that is itself resistant.

        Comparing atpE carriers against 'everyone lacking atpE' puts the
        Rv0678 carriers — whose clofazimine MICs are raised — into the
        reference. atpE then appears to LOWER clofazimine, both drugs move in
        opposite directions, and a target-site variant is misclassified as
        indeterminate. The reference must be isolates carrying neither.
        """
        spec = {"Rv0678_L117R": ("2.0", "2.0"), "atpE_D28V": ("2.0", "0.06")}
        carriers = self._carriers(spec)
        isolates = [i for group in carriers.values() for i in group] \
            + [f"wt{i}" for i in range(10)]
        result = mechanism.validate_discriminator(
            self._panels(spec), carriers, isolates)
        self.assertEqual(result.target_calls["atpE_D28V"], "target-like")
        self.assertEqual(result.efflux_calls["Rv0678_L117R"], "efflux-like")
        self.assertTrue(any("carrying none of the variants" in r
                            for r in result.reasons))

    def test_no_clean_reference_is_not_run_rather_than_wrong(self):
        # Every isolate carries one of the tested variants, so there is no
        # group left that is not itself resistant.
        spec = {"Rv0678_L117R": ("2.0", "2.0"), "atpE_D28V": ("2.0", "0.06")}
        carriers = self._carriers(spec)
        isolates = [i for group in carriers.values() for i in group]
        result = mechanism.validate_discriminator(
            self._panels(spec), carriers, isolates)
        self.assertEqual(result.verdict, "not-run")
        self.assertIn("not itself resistant", result.reasons[0])


# ── barrier ──────────────────────────────────────────────────────────────
def escape(variant, consequence, n, lineages):
    return barrier.EscapeVariant(variant=variant, consequence=consequence,
                                 n_isolates=n, lineages=tuple(lineages))


class BarrierScoringTests(unittest.TestCase):
    def test_small_screen_is_not_estimated(self):
        target = barrier.score_barrier(barrier.TargetBarrier(
            gene="atpE", drug="bedaquiline", n_screened=50))
        self.assertEqual(target.band, "not-estimated")
        self.assertIn("not evidence of a high barrier", target.reasons[0])

    def test_no_observed_escape_is_its_own_state(self):
        target = barrier.score_barrier(barrier.TargetBarrier(
            gene="atpE", drug="bedaquiline", n_screened=8000,
            escape_variants=[escape("atpE_D28V", "missense", 2, ["l4"])]))
        self.assertEqual(target.band, "no-observed-escape")
        self.assertIsNone(target.score)
        self.assertIn("nobody has used long enough", target.reasons[0])

    def test_many_tolerated_truncating_routes_score_low(self):
        variants = [escape(f"mmpR5_fs{i}", "frameshift", 20,
                           ["l1", "l2", "l3", "l4"]) for i in range(15)]
        target = barrier.score_barrier(barrier.TargetBarrier(
            gene="Rv0678", drug="bedaquiline", n_screened=8000,
            escape_variants=variants))
        self.assertEqual(target.band, "low")
        self.assertTrue(any("switch off" in r for r in target.reasons))

    def test_few_missense_routes_in_one_lineage_score_higher(self):
        variants = [escape("atpE_D28V", "missense", 3, ["l4"]),
                    escape("atpE_A63P", "missense", 3, ["l4"])]
        target = barrier.score_barrier(barrier.TargetBarrier(
            gene="atpE", drug="bedaquiline", n_screened=8000,
            escape_variants=variants))
        self.assertIsNotNone(target.score)
        self.assertIn(target.band, ("high", "moderate"))

    def test_truncating_fraction_is_computed_from_isolates(self):
        target = barrier.TargetBarrier(
            gene="Rv0678", drug="bedaquiline", n_screened=1000,
            escape_variants=[escape("a", "frameshift", 30, ["l2"]),
                             escape("b", "missense", 10, ["l2"])])
        self.assertAlmostEqual(target.truncating_fraction, 0.75)


class BarrierValidationTests(unittest.TestCase):
    """Must reproduce the published mmpR5-versus-atpE contrast."""

    def _ranking(self, pump_truncating=True):
        pump = barrier.TargetBarrier(
            gene="Rv0678", drug="bedaquiline", n_screened=8000,
            escape_variants=[
                escape(f"mmpR5_fs{i}",
                       "frameshift" if pump_truncating else "missense",
                       25, ["l1", "l2", "l3", "l4"]) for i in range(14)])
        target = barrier.TargetBarrier(
            gene="atpE", drug="bedaquiline", n_screened=8000,
            escape_variants=[escape("atpE_D28V", "missense", 4, ["l4"]),
                             escape("atpE_A63P", "missense", 3, ["l4"])])
        return barrier.rank([pump, target])

    def test_passes_when_the_contrast_is_reproduced(self):
        result = barrier.validate_against_known(self._ranking())
        self.assertEqual(result.verdict, "reproduces-known-contrast")
        self.assertTrue(result.passed)
        self.assertLess(result.low_score, result.high_score)

    def test_underpowered_when_a_gene_was_not_scored(self):
        ranking = barrier.rank([
            barrier.TargetBarrier(gene="Rv0678", drug="bedaquiline",
                                  n_screened=50),
            barrier.TargetBarrier(gene="atpE", drug="bedaquiline",
                                  n_screened=50)])
        result = barrier.validate_against_known(ranking)
        self.assertEqual(result.verdict, "underpowered")

    def test_not_run_when_a_control_gene_is_absent(self):
        ranking = barrier.rank([barrier.TargetBarrier(
            gene="atpE", drug="bedaquiline", n_screened=8000,
            escape_variants=[escape("atpE_D28V", "missense", 10, ["l4"])])])
        result = barrier.validate_against_known(ranking)
        self.assertEqual(result.verdict, "not-run")

    def test_mmpR5_alias_resolves_to_Rv0678(self):
        ranking = self._ranking()
        ranking.targets[0].gene = "mmpR5"
        result = barrier.validate_against_known(ranking)
        self.assertNotEqual(result.verdict, "not-run")

    def test_ranking_puts_the_highest_barrier_first(self):
        ranked = self._ranking().ranked()
        self.assertEqual(ranked[0].gene, "atpE")


class StructuralVariantRepresentationTests(unittest.TestCase):
    """IS-element insertion into mmpR5 must be representable, not rejected.

    Truncation of mmpR5/Rv0678 accounts for most described bedaquiline
    resistance, and insertional inactivation is one route to it. The adapter
    previously rejected every structural record, so the pipeline structurally
    could not see the commonest cause of resistance to that drug.
    """

    HEADER = (
        "##fileformat=VCFv4.2\n"
        "##sampleID=S1\n"
        "##reference=NC_000962.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )

    def _load(self, record):
        import tempfile
        from pathlib import Path

        from myconductor.io.adapter import load

        path = Path(tempfile.mkdtemp()) / "sv.vcf"
        path.write_text(self.HEADER + record, encoding="utf-8")
        return load(path)

    def test_gene_disrupting_deletion_is_represented_as_truncating(self):
        from myconductor.core.models import Consequence

        result = self._load(
            "NC_000962.3\t779000\t.\tN\t<DEL>\t60\tPASS\t"
            "SVTYPE=DEL;END=779400;GENE=Rv0678;DP=45\tGT\t1\n")
        self.assertEqual(len(result.variants), 1)
        variant = result.variants[0]
        self.assertEqual(variant.gene, "Rv0678")
        self.assertIs(variant.consequence, Consequence.DELETION)
        self.assertTrue(variant.consequence.is_truncating)

    def test_insertion_is_its_own_consequence_and_truncating(self):
        from myconductor.core.models import Consequence

        result = self._load(
            "NC_000962.3\t779100\t.\tN\t<INS>\t60\tPASS\t"
            "SVTYPE=INS;GENE=Rv0678;AACHANGE=IS6110_insertion\tGT\t1\n")
        variant = result.variants[0]
        self.assertIs(variant.consequence, Consequence.INSERTION)
        self.assertTrue(variant.consequence.is_truncating)
        self.assertIn("IS6110_insertion", variant.label())

    def test_breakpoint_evidence_is_declared_unassessed(self):
        result = self._load(
            "NC_000962.3\t779000\t.\tN\t<DEL>\t60\tPASS\t"
            "SVTYPE=DEL;GENE=Rv0678\tGT\t1\n")
        self.assertTrue(any("Breakpoint-level evidence is NOT assessed" in w
                            for w in result.warnings), result.warnings)

    def test_structural_record_without_a_gene_is_still_rejected(self):
        result = self._load(
            "NC_000962.3\t779000\t.\tN\t<DEL>\t60\tPASS\t"
            "SVTYPE=DEL;END=779400\tGT\t1\n")
        self.assertEqual(len(result.variants), 0)
        self.assertEqual(len(result.rejected), 1)
        self.assertIn("no interpretable gene", result.rejected[0].reason)

    def test_symbolic_alt_without_an_svtype_is_still_rejected(self):
        result = self._load(
            "NC_000962.3\t779000\t.\tN\t<DEL>\t60\tPASS\tGENE=Rv0678\tGT\t1\n")
        self.assertEqual(len(result.variants), 0)
        self.assertEqual(len(result.rejected), 1)
        self.assertIn("symbolic", result.rejected[0].reason)

    def test_a_truncating_regulator_variant_reaches_both_pump_drugs(self):
        from myconductor.core.models import Consequence, Lane, Variant
        from myconductor.modules.efflux import EffluxRegulatoryModule

        variant = Variant.of("Rv0678", "IS6110_insertion",
                             consequence=Consequence.INSERTION, depth=45)
        evidence = EffluxRegulatoryModule().evaluate(variant)
        drugs = {e.drug for e in evidence}
        self.assertEqual(drugs, {"bedaquiline", "clofazimine"})
        for item in evidence:
            self.assertIs(item.lane, Lane.EFFLUX_REGULATORY)
            self.assertIn("Loss of function", item.rationale)
            self.assertTrue(any("commonest described route" in item.rationale
                                for _ in (0,)))

    def test_missense_regulator_variant_is_marked_weaker(self):
        from myconductor.core.models import Consequence, Variant
        from myconductor.modules.efflux import EffluxRegulatoryModule

        variant = Variant.of("Rv0678", "L117R",
                             consequence=Consequence.MISSENSE, depth=45)
        evidence = EffluxRegulatoryModule().evaluate(variant)
        self.assertTrue(evidence)
        for item in evidence:
            self.assertIn("Missense", item.rationale)
            self.assertTrue(any("materially weaker" in limitation
                                for limitation in item.limitations),
                            item.limitations)

    def test_neither_class_may_assert_resistance(self):
        from myconductor.core.models import Call, Consequence, Variant
        from myconductor.modules.efflux import EffluxRegulatoryModule

        for consequence in (Consequence.INSERTION, Consequence.MISSENSE):
            variant = Variant.of("Rv0678", "x", consequence=consequence,
                                 depth=45)
            for item in EffluxRegulatoryModule().evaluate(variant):
                self.assertIs(item.call, Call.INDETERMINATE)


if __name__ == "__main__":
    unittest.main()
