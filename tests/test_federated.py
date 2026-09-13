import unittest
from datetime import timedelta

from myconductor.federated.catalogue_update import (
    EvidenceLedger,
    IsolateObservation,
    PromotionThresholds,
    ReviewQueue,
    ReviewState,
    aggregate,
    assess_drift,
    evaluate_candidate,
)
from myconductor.federated.transport import (
    DisclosureThresholds,
    Site,
    SiteRegistry,
    Submission,
    SubmissionVerifier,
    TransportError,
    build_submission,
    k_anonymity_violations,
    suppress_small_cells,
)


def observations(n, *, variant, resistant, lineage_cycle=("l2", "l4"),
                 sites=3, drug="rifampicin", prefix="ISO", mic=None):
    out = []
    for i in range(n):
        out.append(IsolateObservation(
            isolate_id=f"{prefix}-{i}", site_id=f"SITE-{i % sites}",
            drug=drug, phenotype=resistant, variant_keys=(variant,),
            lineage=lineage_cycle[i % len(lineage_cycle)],
            dst_method="MGIT", lab_quality="accredited", mic=mic,
        ))
    return out


class IsolateLevelDenominatorTests(unittest.TestCase):
    """The defect: counting variants instead of isolates."""

    def test_an_isolate_is_counted_once_per_drug(self):
        duplicate = [
            IsolateObservation("ISO-1", "SITE-A", "rifampicin", True,
                               ("rpoB_S450L",), lineage="l2"),
            IsolateObservation("ISO-1", "SITE-A", "rifampicin", True,
                               ("rpoB_S450L",), lineage="l2"),
        ]
        assoc = aggregate([duplicate])[("rpoB_S450L", "rifampicin")]
        self.assertEqual(assoc.n_isolates, 1)

    def test_co_occurring_variants_are_tracked_not_credited_equally(self):
        obs = [IsolateObservation(
            "ISO-1", "SITE-A", "rifampicin", True,
            ("rpoB_S450L", "Rv1234_A1B"), lineage="l2")]
        aggregated = aggregate([obs])
        causal = aggregated[("rpoB_S450L", "rifampicin")]
        self.assertEqual(causal.cooccurring, {"Rv1234_A1B": 1})

    def test_observations_without_a_phenotype_are_ignored(self):
        obs = [IsolateObservation("ISO-1", "SITE-A", "rifampicin", None,
                                  ("rpoB_S450L",))]
        self.assertEqual(aggregate([obs]), {})

    def test_aggregation_merges_across_sites(self):
        a = observations(10, variant="rpoB_S450L", resistant=True, prefix="A")
        b = observations(10, variant="rpoB_S450L", resistant=True, prefix="B")
        assoc = aggregate([a, b])[("rpoB_S450L", "rifampicin")]
        self.assertEqual(assoc.n_isolates, 20)
        self.assertEqual(assoc.n_sites, 3)


class ConfoundingTests(unittest.TestCase):
    def _assoc(self, obs):
        return aggregate([obs])[("V_1", "rifampicin")]

    def test_lineage_marker_is_refused_despite_perfect_association(self):
        # 60 resistant isolates, 100% association, three sites -- but only
        # ever seen in one lineage.
        obs = observations(60, variant="V_1", resistant=True,
                           lineage_cycle=("l2",))
        obs += observations(20, variant="V_other", resistant=False,
                            lineage_cycle=("l2",), prefix="S")
        candidate = evaluate_candidate(self._assoc(obs))
        self.assertFalse(candidate.eligible)
        self.assertTrue(any("lineage" in r for r in candidate.blocked_reasons),
                        candidate.blocked_reasons)

    def test_association_across_two_lineages_passes_the_lineage_check(self):
        obs = observations(80, variant="V_1", resistant=True,
                           lineage_cycle=("l2", "l4"))
        obs += observations(20, variant="V_wt", resistant=False,
                            lineage_cycle=("l2", "l4"), prefix="S")
        candidate = evaluate_candidate(self._assoc(obs))
        self.assertFalse(any("lineage" in r
                             for r in candidate.blocked_reasons),
                         candidate.blocked_reasons)

    def test_hitchhiker_with_a_known_determinant_is_refused(self):
        obs = []
        for i in range(80):
            obs.append(IsolateObservation(
                isolate_id=f"ISO-{i}", site_id=f"SITE-{i % 3}",
                drug="rifampicin", phenotype=True,
                variant_keys=("V_1", "rpoB_S450L"),
                lineage=("l2", "l4")[i % 2], lab_quality="accredited",
            ))
        for i in range(20):
            obs.append(IsolateObservation(
                isolate_id=f"S-{i}", site_id=f"SITE-{i % 3}",
                drug="rifampicin", phenotype=False, variant_keys=("V_wt",),
                lineage=("l2", "l4")[i % 2], lab_quality="accredited",
            ))
        candidate = evaluate_candidate(
            self._assoc(obs), known_resistance_variants={"rpoB_S450L"})
        self.assertFalse(candidate.eligible)
        self.assertTrue(any("co-occurs" in r for r in candidate.blocked_reasons),
                        candidate.blocked_reasons)

    def test_single_site_evidence_is_refused(self):
        obs = observations(80, variant="V_1", resistant=True, sites=1)
        obs += observations(20, variant="V_wt", resistant=False, sites=1,
                            prefix="S")
        candidate = evaluate_candidate(self._assoc(obs))
        self.assertTrue(any("site" in r for r in candidate.blocked_reasons))

    def test_missing_susceptible_controls_are_refused(self):
        obs = observations(60, variant="V_1", resistant=True)
        candidate = evaluate_candidate(self._assoc(obs))
        self.assertTrue(any("control" in r for r in candidate.blocked_reasons))

    def test_susceptible_controls_need_not_carry_the_variant(self):
        # The control group is the variant-ABSENT isolates. Requiring
        # susceptible isolates that carry the variant would be the wrong
        # control, and would penalise precisely the strongest determinants.
        obs = observations(80, variant="V_1", resistant=True)
        obs += observations(20, variant="V_wt", resistant=False, prefix="S")
        assoc = self._assoc(obs)
        self.assertEqual(assoc.n_susceptible, 0, "no carrier tested susceptible")
        self.assertEqual(assoc.cohort_susceptible, 20)
        self.assertTrue(assoc.has_variant_absent_controls)
        candidate = evaluate_candidate(assoc)
        self.assertTrue(candidate.eligible, candidate.blocked_reasons)

    def test_variant_in_every_tested_isolate_is_refused(self):
        # A 100% association with nothing to compare it against says nothing.
        obs = observations(80, variant="V_1", resistant=True)
        obs += observations(20, variant="V_1", resistant=False, prefix="S")
        candidate = evaluate_candidate(self._assoc(obs))
        self.assertFalse(candidate.eligible)
        self.assertTrue(
            any("variant-absent comparison group" in r
                for r in candidate.blocked_reasons),
            candidate.blocked_reasons,
        )

    def test_cohort_totals_count_isolates_without_the_variant(self):
        obs = observations(80, variant="V_1", resistant=True)
        obs += observations(20, variant="V_wt", resistant=False, prefix="S")
        assoc = self._assoc(obs)
        self.assertEqual(assoc.cohort_resistant, 80)
        self.assertEqual(assoc.cohort_susceptible, 20)
        self.assertEqual(assoc.absent_susceptible, 20)
        self.assertEqual(assoc.cohort_size, 100)

    def test_a_clean_association_is_eligible_but_needs_curation(self):
        obs = observations(80, variant="V_1", resistant=True)
        obs += observations(20, variant="V_wt", resistant=False, prefix="S")
        candidate = evaluate_candidate(self._assoc(obs))
        self.assertTrue(candidate.eligible, candidate.blocked_reasons)
        self.assertTrue(any("curation" in n for n in candidate.notes))

    def test_mic_requirement_can_be_enforced(self):
        obs = observations(80, variant="V_1", resistant=True)
        obs += observations(20, variant="V_wt", resistant=False, prefix="S")
        candidate = evaluate_candidate(
            self._assoc(obs), thresholds=PromotionThresholds(require_mic=True))
        self.assertTrue(any("MIC" in r for r in candidate.blocked_reasons))


class GovernanceTests(unittest.TestCase):
    def _eligible_candidate(self):
        obs = observations(80, variant="V_1", resistant=True)
        obs += observations(20, variant="V_wt", resistant=False, prefix="S")
        return evaluate_candidate(aggregate([obs])[("V_1", "rifampicin")])

    def test_ledger_chain_verifies(self):
        ledger = EvidenceLedger()
        ledger.append("submitted", {"x": 1})
        ledger.append("approved", {"x": 1})
        ok, problem = ledger.verify()
        self.assertTrue(ok, problem)

    def test_ledger_detects_tampering(self):
        ledger = EvidenceLedger()
        ledger.append("submitted", {"x": 1})
        ledger.entries[0].payload["x"] = 2
        ok, problem = ledger.verify()
        self.assertFalse(ok)
        self.assertIn("altered", problem)

    def test_ineligible_candidates_cannot_be_queued(self):
        obs = observations(5, variant="V_1", resistant=True)
        candidate = evaluate_candidate(aggregate([obs])[("V_1", "rifampicin")])
        with self.assertRaises(ValueError):
            ReviewQueue().submit(candidate)

    def test_promotion_is_never_automatic(self):
        queue = ReviewQueue()
        item = queue.submit(self._eligible_candidate())
        self.assertIs(item.state, ReviewState.PENDING)
        self.assertEqual(queue.approved_entries(), [])

    def test_approval_records_the_reviewer_and_releases_the_entry(self):
        queue = ReviewQueue()
        queue.submit(self._eligible_candidate())
        queue.approve("V_1", "rifampicin", reviewer="curator-1",
                      rationale="consistent with published evidence")
        entries = queue.approved_entries()
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]["call"], "resistant")
        self.assertIn("provisional", entries[0]["who_grade"])

    def test_rollback_is_possible_and_audited(self):
        queue = ReviewQueue()
        queue.submit(self._eligible_candidate())
        queue.approve("V_1", "rifampicin", "curator-1", "looks right")
        queue.roll_back("V_1", "rifampicin", "curator-2", "refuted by new DST")
        self.assertEqual(queue.approved_entries(), [])
        actions = [e.action for e in queue.ledger.entries]
        self.assertEqual(actions, ["submitted", "approved", "rolled_back"])
        self.assertTrue(queue.ledger.verify()[0])

    def test_entry_cannot_be_built_from_an_ineligible_candidate(self):
        obs = observations(5, variant="V_1", resistant=True)
        candidate = evaluate_candidate(aggregate([obs])[("V_1", "rifampicin")])
        with self.assertRaises(ValueError):
            candidate.to_catalogue_entry()


class DriftTests(unittest.TestCase):
    def test_too_few_points_is_not_assessed_rather_than_no_drift(self):
        report = assess_drift([0.9, 0.88])
        self.assertFalse(report.assessed)
        self.assertFalse(report.alarm)
        self.assertIn("not the same as no drift", report.detail)

    def test_stable_performance_raises_no_alarm(self):
        report = assess_drift([0.90, 0.91, 0.89, 0.90, 0.92])
        self.assertTrue(report.assessed)
        self.assertFalse(report.alarm)

    def test_decay_raises_an_alarm(self):
        report = assess_drift([0.75, 0.72, 0.70, 0.71, 0.69])
        self.assertTrue(report.assessed)
        self.assertTrue(report.alarm)


class TransportTests(unittest.TestCase):
    def setUp(self):
        self.registry = SiteRegistry()
        self.site = self.registry.register(Site(
            site_id="SITE-0", name="Lab 0", secret=SiteRegistry.new_secret()))
        self.obs = [IsolateObservation(
            "ISO-1", "SITE-0", "rifampicin", True, ("rpoB_S450L",),
            lineage="l2")]

    def test_round_trip_verifies(self):
        submission = build_submission(self.site, self.obs)
        accepted = SubmissionVerifier(self.registry).verify(submission)
        self.assertEqual(len(accepted), 1)
        self.assertEqual(accepted[0].isolate_id, "ISO-1")

    def test_tampered_payload_is_rejected(self):
        submission = build_submission(self.site, self.obs)
        submission.payload[0]["phenotype"] = False
        with self.assertRaises(TransportError) as ctx:
            SubmissionVerifier(self.registry).verify(submission)
        self.assertIn("signature mismatch", str(ctx.exception))

    def test_replay_is_rejected(self):
        submission = build_submission(self.site, self.obs)
        verifier = SubmissionVerifier(self.registry)
        verifier.verify(submission)
        with self.assertRaises(TransportError) as ctx:
            verifier.verify(submission)
        self.assertIn("replay", str(ctx.exception))

    def test_stale_submission_is_rejected(self):
        submission = build_submission(self.site, self.obs)
        verifier = SubmissionVerifier(self.registry,
                                      max_age=timedelta(seconds=0))
        with self.assertRaises(TransportError):
            verifier.verify(submission)

    def test_unknown_site_is_rejected(self):
        submission = Submission("SITE-9", "n", "2024-01-01T00:00:00+00:00")
        with self.assertRaises(TransportError):
            SubmissionVerifier(self.registry).verify(submission)

    def test_site_without_a_secret_cannot_register(self):
        with self.assertRaises(TransportError):
            SiteRegistry().register(Site(site_id="X", name="X"))

    def test_observation_from_another_site_cannot_be_submitted(self):
        foreign = [IsolateObservation("ISO-2", "SITE-9", "rifampicin", True,
                                      ("rpoB_S450L",))]
        with self.assertRaises(TransportError):
            build_submission(self.site, foreign)


class DisclosureTests(unittest.TestCase):
    def test_rare_cell_violates_the_disclosure_floor(self):
        obs = [IsolateObservation("ISO-1", "SITE-A", "rifampicin", True,
                                  ("rare_V1",), lineage="l2")]
        violations = k_anonymity_violations(obs)
        self.assertTrue(any("isolate" in v for v in violations))
        self.assertTrue(any("site" in v for v in violations))

    def test_large_cell_passes(self):
        obs = observations(20, variant="common_V", resistant=True)
        self.assertEqual(k_anonymity_violations(obs), [])

    def test_small_cells_are_suppressed(self):
        obs = observations(20, variant="common_V", resistant=True)
        obs.append(IsolateObservation(
            "ISO-rare", "SITE-0", "rifampicin", True, ("rare_V",),
            lineage="l2"))
        kept, suppressed = suppress_small_cells(obs)
        remaining = {v for o in kept for v in o.variant_keys}
        self.assertNotIn("rare_V", remaining)
        self.assertTrue(suppressed)

    def test_thresholds_are_configurable(self):
        obs = observations(6, variant="V", resistant=True)
        strict = DisclosureThresholds(min_isolates_per_cell=10)
        self.assertTrue(k_anonymity_violations(obs, strict))


if __name__ == "__main__":
    unittest.main()
