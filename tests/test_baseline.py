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

    def test_partial_locus_examination_abstains(self):
        """Examined katG but not inhA: resistance could be at the unexamined one."""
        multi_gene = {"katG_S315T", "inhA_I21T", "fabG1_c.-15C>T"}
        predictions, abstained = baseline.catalogue_predictions(
            [isolate("a", assessed={"katG_S315T"})], {"a": "S"}, multi_gene)
        self.assertIsNone(predictions[0].predicted)
        self.assertEqual(abstained[baseline.ABSTAIN_INCOMPLETE], 1)

    def test_examining_every_locus_answers_susceptible(self):
        multi_gene = {"katG_S315T", "inhA_I21T", "fabG1_c.-15C>T"}
        predictions, _ = baseline.catalogue_predictions(
            [isolate("a", assessed={"katG_S315T", "inhA_I21T",
                                    "fabG1_c.-15C>T"})],
            {"a": "S"}, multi_gene)
        self.assertEqual(predictions[0].predicted, "S")

    def test_one_variant_per_locus_is_enough_to_examine_that_locus(self):
        """The artifact this rule was changed to remove.

        Counting graded *variants* rather than loci made the threshold silently
        stricter for better-studied drugs. Rifampicin has 136 catalogued
        determinants, almost all in rpoB; requiring 95% of them to be
        individually examined meant the catalogue never answered susceptible
        for it at all, every answered isolate was an R call, and specificity
        came out at exactly zero. That was the denominator, not the catalogue.
        """
        many_in_one_gene = {f"rpoB_v{i}" for i in range(136)}
        predictions, _ = baseline.catalogue_predictions(
            [isolate("a", assessed={"rpoB_v0"})], {"a": "S"}, many_in_one_gene)
        self.assertEqual(predictions[0].predicted, "S")

    def test_locus_is_the_first_label_segment(self):
        self.assertEqual(baseline.locus_of("rpoB_p.Ser450Leu"), "rpoB")
        self.assertEqual(baseline.locus_of("fabG1_c.-15C>T"), "fabG1")
        self.assertEqual(baseline.locus_of(""), "")

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


class DiagnosticMetricTests(unittest.TestCase):
    """A composite error rate cannot say which way a tool is wrong.

    Missing resistance and over-calling it have opposite clinical
    consequences — one puts a patient on a failing drug, the other withholds a
    working one — and a single error rate averages them together. These come
    from ``metrics.score_accuracy``, the same code that scores any engine,
    rather than a second implementation that could drift from it.
    """

    def _cohort(self, tp=8, fn=2, tn=8, fp=2):
        """Carriers are called R, non-carriers with full coverage are called S."""
        isolates, truths = [], {}
        for i in range(tp):
            isolates.append(isolate(f"tp{i}", carries={"rpoB_S450L"},
                                    assessed=DETERMINANTS))
            truths[f"tp{i}"] = "R"
        for i in range(fp):
            isolates.append(isolate(f"fp{i}", carries={"rpoB_S450L"},
                                    assessed=DETERMINANTS))
            truths[f"fp{i}"] = "S"
        for i in range(tn):
            isolates.append(isolate(f"tn{i}", assessed=DETERMINANTS))
            truths[f"tn{i}"] = "S"
        for i in range(fn):
            isolates.append(isolate(f"fn{i}", assessed=DETERMINANTS))
            truths[f"fn{i}"] = "R"
        return isolates, truths

    def test_sensitivity_and_specificity_are_separated(self):
        isolates, truths = self._cohort(tp=8, fn=2, tn=8, fp=2)
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertAlmostEqual(result.accuracy.sensitivity, 0.8)
        self.assertAlmostEqual(result.accuracy.specificity, 0.8)

    def test_vme_and_me_are_derived_not_declared(self):
        """VME is 1 - sensitivity and ME is 1 - specificity, by construction."""
        isolates, truths = self._cohort(tp=9, fn=1, tn=7, fp=3)
        accuracy = baseline.measure_drug(DRUG, isolates, truths,
                                         DETERMINANTS).accuracy
        self.assertAlmostEqual(accuracy.vme_rate, 1 - accuracy.sensitivity)
        self.assertAlmostEqual(accuracy.me_rate, 1 - accuracy.specificity)

    def test_a_false_susceptible_is_a_very_major_error(self):
        isolates, truths = self._cohort(tp=9, fn=1, tn=10, fp=0)
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.n_false_susceptible, 1)
        self.assertEqual(result.accuracy.false_negative, 1)
        self.assertGreater(result.accuracy.vme_rate, 0)

    def test_a_false_resistant_is_a_major_error_and_counted_apart(self):
        isolates, truths = self._cohort(tp=10, fn=0, tn=7, fp=3)
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.accuracy.false_positive, 3)
        self.assertEqual(result.n_false_susceptible, 0)
        self.assertEqual(result.accuracy.vme_rate, 0.0)

    def test_intervals_accompany_the_point_estimates_when_powered(self):
        # rifampicin's pre-registered target requires 30 evaluable isolates.
        isolates, truths = self._cohort(tp=16, fn=4, tn=16, fp=4)
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertTrue(result.accuracy.powered)
        rows = baseline.baseline_rows({DRUG: result})
        self.assertTrue(rows[0]["sensitivity_ci"])
        self.assertIn("-", rows[0]["sensitivity_ci"])
        self.assertEqual(rows[0]["powered"], "yes")

    def test_no_interval_is_offered_on_an_underpowered_sample(self):
        """An interval from too few observations is a lie about precision."""
        isolates, truths = self._cohort(tp=4, fn=1, tn=4, fp=1)
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertFalse(result.accuracy.powered)
        rows = baseline.baseline_rows({DRUG: result})
        self.assertEqual(rows[0]["sensitivity_ci"], "")
        self.assertEqual(rows[0]["powered"], "no")
        # The point estimate is still shown; only the precision claim is not.
        self.assertTrue(rows[0]["sensitivity"])

    def test_abstained_resistant_isolates_are_shown_beside_sensitivity(self):
        """A predictor can always look sensitive by declining the hard ones."""
        isolates, truths = self._cohort()
        for i in range(5):
            isolates.append(isolate(f"skip{i}", assessed=set()))
            truths[f"skip{i}"] = "R"
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.accuracy.abstained_resistant, 5)
        # Excluded from sensitivity, and the exclusion is stated.
        self.assertEqual(result.accuracy.n_phenotype_resistant, 10)
        self.assertTrue(any("excluded from sensitivity" in n
                            for n in result.notes))

    def test_ppv_and_npv_carry_the_prevalence_they_assume(self):
        isolates, truths = self._cohort()
        rows = baseline.baseline_rows(
            {DRUG: baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)})
        self.assertTrue(rows[0]["ppv"])
        self.assertTrue(rows[0]["resistance_prevalence"])
        self.assertIn("do not", baseline.PPV_NPV_CAVEAT)
        self.assertIn("recomputed at", baseline.PPV_NPV_CAVEAT)


class LineageStratificationTests(unittest.TestCase):
    """A pooled figure can be carried entirely by one lineage."""

    def _cohort(self):
        isolates, truths, lineages = [], {}, {}
        for i in range(10):
            name = f"l2r{i}"
            isolates.append(isolate(name, carries={"rpoB_S450L"},
                                    assessed=DETERMINANTS))
            truths[name] = "R"
            lineages[name] = "L2"
        for i in range(10):
            name = f"l4s{i}"
            isolates.append(isolate(name, assessed=DETERMINANTS))
            truths[name] = "S"
            lineages[name] = "L4"
        return isolates, truths, lineages

    def test_lineages_are_reported_separately_when_supplied(self):
        isolates, truths, lineages = self._cohort()
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS,
                                       lineages=lineages)
        self.assertEqual(sorted(result.by_lineage), ["L2", "L4"])
        self.assertEqual(result.by_lineage["L2"].n_phenotype_resistant, 10)
        self.assertEqual(result.by_lineage["L4"].n_phenotype_susceptible, 10)

    def test_a_single_unknown_bucket_is_not_stratification(self):
        """CRyPTIC has no lineage column, so this is the common case."""
        isolates, truths, _ = self._cohort()
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.by_lineage, {})

    def test_lineage_rows_render_only_where_lineages_exist(self):
        isolates, truths, lineages = self._cohort()
        with_lineage = baseline.measure_drug(DRUG, isolates, truths,
                                             DETERMINANTS, lineages=lineages)
        without = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(len(baseline.lineage_rows({DRUG: with_lineage})), 2)
        self.assertEqual(baseline.lineage_rows({DRUG: without}), [])
        for column in baseline.LINEAGE_COLUMNS:
            self.assertIn(column, baseline.lineage_rows({DRUG: with_lineage})[0])


class MeasurabilityTests(unittest.TestCase):
    """Which drugs this pairing cannot evaluate, and which side is missing."""

    CATALOGUE = ("rifampicin", "isoniazid", "capreomycin", "bedaquiline")
    PHENOTYPE = ("rifampicin", "isoniazid", "bedaquiline", "rifabutin")

    def _row(self, drug):
        rows = baseline.measurability(self.CATALOGUE, self.PHENOTYPE)
        return {r["drug"]: r for r in rows}[drug]

    def test_a_drug_in_both_sources_is_measurable(self):
        row = self._row("rifampicin")
        self.assertEqual(row["measurable"], "yes")
        self.assertEqual(row["missing_side"], "")

    def test_a_catalogued_drug_with_no_phenotype_names_that_side(self):
        row = self._row("capreomycin")
        self.assertEqual(row["measurable"], "no")
        self.assertEqual(row["missing_side"], "phenotype")
        self.assertIn("not scored", row["consequence"])

    def test_a_phenotyped_drug_with_no_catalogue_entry_names_that_side(self):
        row = self._row("rifabutin")
        self.assertEqual(row["missing_side"], "catalogue")
        self.assertIn("nothing predicts it", row["consequence"])

    def test_pretomanid_is_absent_from_both_and_says_so(self):
        """A quarter of the BPaL regimen, unmeasurable from these sources."""
        row = self._row("pretomanid")
        self.assertEqual(row["missing_side"], "both")
        self.assertIn("no pairing of these sources", row["consequence"])

    def test_every_bpalm_component_is_reported_on(self):
        drugs = {r["drug"] for r in baseline.measurability(self.CATALOGUE,
                                                           self.PHENOTYPE)}
        for component in ("bedaquiline", "pretomanid", "linezolid",
                          "moxifloxacin", "clofazimine"):
            self.assertIn(component, drugs)


class IsolatesNeededTests(unittest.TestCase):
    """"Collect about N more" beats "not estimable"."""

    def test_an_estimable_drug_needs_none(self):
        isolates = [isolate(f"s{i}", assessed=DETERMINANTS) for i in range(12)]
        truths = {f"s{i}": "S" for i in range(12)}
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertIsNone(baseline.isolates_needed(result))

    def test_a_rare_resistance_drug_reports_a_planning_figure(self):
        isolates, truths = [], {}
        for i in range(8):
            isolates.append(isolate(f"s{i}", assessed=DETERMINANTS))
            truths[f"s{i}"] = "S"
        isolates.append(isolate("r0", carries={"rpoB_S450L"},
                                assessed=DETERMINANTS))
        truths["r0"] = "R"
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertFalse(result.estimable)
        needed = baseline.isolates_needed(result)
        self.assertIsNotNone(needed)
        self.assertGreater(needed, 0)

    def test_a_cohort_with_no_resistance_cannot_be_extrapolated_from(self):
        isolates = [isolate(f"s{i}", assessed=DETERMINANTS) for i in range(5)]
        truths = {f"s{i}": "S" for i in range(5)}
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertIsNone(baseline.isolates_needed(result))


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

    def test_a_namespace_mismatch_is_named_not_reported_as_no_coverage(self):
        """Labels against coordinate keys abstains on everything, silently.

        That looks identical to a cohort nobody sequenced deeply, and it
        produced an entire table of plausible wrong numbers once. The
        determinant set here is coordinate-keyed while the isolates are
        label-keyed, so nothing can ever match.
        """
        coordinate_keys = {"NC_000962.3:Chromosome:761155:C>T"}
        isolates = [isolate(f"s{i}", assessed=DETERMINANTS) for i in range(12)]
        truths = {f"s{i}": "S" for i in range(12)}
        result = baseline.measure_drug(DRUG, isolates, truths, coordinate_keys)
        self.assertEqual(result.n_answered, 0)
        self.assertTrue(any("keyed differently" in n for n in result.notes))

    def test_genuine_absence_of_coverage_is_not_blamed_on_namespaces(self):
        # Determinants and isolates share a namespace; the cohort simply was
        # not examined. No mismatch note should appear.
        isolates = [isolate(f"s{i}", assessed=set()) for i in range(12)]
        truths = {f"s{i}": "S" for i in range(12)}
        result = baseline.measure_drug(DRUG, isolates, truths, DETERMINANTS)
        self.assertEqual(result.n_answered, 0)
        self.assertFalse(any("keyed differently" in n for n in result.notes))

    def test_measure_uses_labels_not_the_coordinate_union(self):
        """The bug this guards: unioning namespaces caps the achievable
        assessed fraction and abstains on everything."""
        index = DeterminantIndex(
            by_drug_label={"rifampicin": {"rpoB_S450L"}},
            by_drug_coordinate={"rifampicin": {
                "NC_000962.3:Chromosome:761155:C>T",
                "NC_000962.3:Chromosome:761140:A>G"}})
        isolates = [isolate(f"cr_ERR{i}", assessed={"rpoB_S450L"})
                    for i in range(12)]
        rows = [{"ENA_RUN": f"ERR{i}", "RIF_BINARY_PHENOTYPE": "S"}
                for i in range(12)]
        results = baseline.measure(isolates, rows, index, drugs=["rifampicin"])
        # With the union the denominator would be 3 and nothing would answer.
        self.assertEqual(results["rifampicin"].n_answered, 12)
        self.assertEqual(results["rifampicin"].coverage, 1.0)

    def test_the_caveat_names_the_shared_isolate_problem(self):
        self.assertIn("share underlying isolates", baseline.CAVEAT)
        self.assertIn("not an independent accuracy estimate", baseline.CAVEAT)


if __name__ == "__main__":
    unittest.main()
