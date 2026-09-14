import copy
import hashlib
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from myconductor.core.context import SampleContext, InterpretationPolicy, ValidationScope
from myconductor.core.models import Call, DrugEvidence, EngineRef, Lane, Tier, QCFinding
from myconductor.core.pipeline import Myconductor
from myconductor.io.phenotypes import PhenotypeObservation
from myconductor.modules.synthesis import EvidenceReconciler
from myconductor.modules.investigation import plan_follow_up
from myconductor.federated.case_review import CaseReviewStore
from myconductor.federated.change_impact import compare_reports
from myconductor.reporting.json_report import report_dict, validate_report
from myconductor.reporting.html_report import render_html
from myconductor.reporting.fhir import to_fhir_bundle
from myconductor.adapters.base import EngineReport, AdapterSchemaError
from myconductor.adapters.ntm_profiler import NTMProfilerAdapter
from myconductor.adapters.gnomonicus import GnomonicusAdapter
from mycobench.workflow_evaluation import evaluate
from mycobench.analysis.genotypes import CoordinateIndex, parse_vcf, fetch_vcf, GenotypeError
from mycobench.analysis.strata import Isolate, DeterminantIndex, stratify
from mycobench.analysis.metadata import apply_metadata


def phenotype(**kwargs):
    data = dict(observation_id="dst1", sample_id="s", isolate_id="i", site_id="site",
                organism="mtbc", drug="linezolid", result="susceptible", method="laboratory DST",
                measured_at="2026-09-01", laboratory="reference-lab", quality="pass",
                interpretation_standard="lab-SOP", standard_version="1", source="record/1")
    data.update(kwargs)
    return PhenotypeObservation(**data)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.input = self.root / "variants.json"
        self.input.write_text(json.dumps({"sample_id": "s", "variants": []}))
        self.context = SampleContext("s", "i", "site", "mtbc", assay="WGS")

    def report(self, **kwargs):
        return Myconductor().analyze(self.input, context=kwargs.pop("context", self.context), **kwargs)

    def test_matched_phenotype_can_establish_s_without_genomic_coverage(self):
        result = next(r for r in self.report(phenotypes=[phenotype()]).drug_results if r.drug == "linezolid")
        self.assertEqual(result.call, Call.SUSCEPTIBLE)
        self.assertEqual(result.genomic_call, Call.NOT_ASSESSED)
        self.assertEqual(result.phenotypic_call, Call.SUSCEPTIBLE)
        self.assertEqual(result.assay_status, "insufficient")

    def test_identity_dimensions_are_checked(self):
        for field in ("sample_id", "isolate_id", "site_id", "organism"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.report(phenotypes=[phenotype(**{field: "different"})])

    def test_failed_or_unknown_phenotype_quality_abstains(self):
        for quality in ("fail", "unknown"):
            result = next(r for r in self.report(phenotypes=[phenotype(quality=quality)]).drug_results if r.drug == "linezolid")
            self.assertEqual(result.call, Call.INDETERMINATE)

    def test_repeat_phenotypes_conflict_without_majority_vote(self):
        obs = [phenotype(), phenotype(observation_id="dst2", result="resistant"),
               phenotype(observation_id="dst3")]
        result = next(r for r in self.report(phenotypes=obs).drug_results if r.drug == "linezolid")
        self.assertEqual(result.call, Call.INDETERMINATE)
        self.assertIsNotNone(result.discordance)

    def test_genomic_qc_failure_does_not_erase_matched_phenotype(self):
        context = replace(self.context, qc=[QCFinding("contamination", "fail", "reviewed")])
        report = self.report(context=context, phenotypes=[phenotype()])
        result = next(r for r in report.drug_results if r.drug == "linezolid")
        self.assertEqual(result.call, Call.SUSCEPTIBLE)
        self.assertFalse(result.genomic_call.is_established)
        self.assertTrue(any(f["category"] == "analytical_failure" for f in report.investigations))

    def test_variant_not_associated_is_not_an_isolate_disagreement(self):
        evidence = [DrugEvidence("linezolid", Call.RESISTANT, Tier.CATALOGUED, Lane.CATALOGUE, scope="variant"),
                    DrugEvidence("linezolid", Call.SUSCEPTIBLE, Tier.CATALOGUED, Lane.CATALOGUE, scope="variant")]
        result = next(r for r in EvidenceReconciler().reconcile(evidence) if r.drug == "linezolid")
        self.assertEqual(result.call, Call.RESISTANT)
        self.assertIsNone(result.discordance)

    def test_engine_identity_and_duplicate_report_rejected(self):
        for reports in ([EngineReport(EngineRef("x"), sample_id="other")],
                        [EngineReport(EngineRef("x"), sample_id="s")] * 2):
            with self.assertRaises(ValueError):
                self.report(engine_reports=reports)

    def test_engine_s_requires_exact_scope_and_qc(self):
        engine = EngineRef("mykrobe", "1", "db", "v1")
        ev = DrugEvidence("linezolid", Call.SUSCEPTIBLE, Tier.CATALOGUED, Lane.ENGINE,
                          engine=engine, asserts_coverage=True)
        policy = InterpretationPolicy("v1", [ValidationScope("mtbc", "WGS", "linezolid",
             "mykrobe", "1", "v1", "external-validation", "reviewer", "2026-09-01")])
        context = replace(self.context, qc=[QCFinding(k, "pass", "reviewed") for k in policy.required_qc])
        def call(c, e=ev):
            return next(r.call for r in EvidenceReconciler(policy=policy, context=c).reconcile([e]) if r.drug == "linezolid")
        self.assertEqual(call(context), Call.SUSCEPTIBLE)
        self.assertEqual(call(self.context), Call.INDETERMINATE)
        self.assertEqual(call(context, replace(ev, engine=replace(engine, database_version="v2"))), Call.INDETERMINATE)

    def test_ntm_early_s_requires_induction_context(self):
        context = SampleContext("s", "i", "site", "mabscessus")
        obs = phenotype(organism="mabscessus", drug="clarithromycin", incubation_days=3)
        self.assertEqual(obs.evidence(context).call, Call.INDETERMINATE)
        self.assertEqual(replace(obs, incubation_days=14).evidence(context).call, Call.SUSCEPTIBLE)
        functional = replace(context, erm41_status="functional", erm41_source="reviewed sequence")
        self.assertEqual(replace(obs, incubation_days=14).evidence(functional).call, Call.INDETERMINATE)
        nonfunctional = replace(context, erm41_status="nonfunctional", erm41_source="reviewed sequence")
        self.assertEqual(obs.evidence(nonfunctional).call, Call.SUSCEPTIBLE)

    def test_report_roundtrip_and_tamper_detection(self):
        data = report_dict(self.report(phenotypes=[phenotype()]))
        validate_report(json.loads(json.dumps(data)))
        for mutate in (lambda d: d["context"].update(site_id="other"),
                       lambda d: d["analysis_manifest"]["policy"].update(version="changed"),
                       lambda d: d["drug_results"][-1]["evidence"].append({"evidence_id": "fake"})):
            changed = copy.deepcopy(data)
            mutate(changed)
            with self.assertRaises(ValueError):
                validate_report(changed)

    def test_case_open_attach_resolve_reopen_and_reload(self):
        store_path = self.root / "cases.json"
        store = CaseReviewStore.load(store_path)
        before = report_dict(self.report())
        case = store.open(before, "reviewer")
        store.transition(case, "in_review", "reviewer", "review started")
        after = report_dict(self.report(phenotypes=[phenotype()]))
        store.attach(case, after, "reviewer", "follow-up DST received")
        ev = next(e["evidence_id"] for r in after["drug_results"] for e in r["evidence"])
        store.transition(case, "resolved", "reviewer", "assessed DST", [ev], "laboratory review completed", 8)
        with self.assertRaises(ValueError):
            store.attach(case, after, "reviewer", "must reopen")
        store.transition(case, "reopened", "reviewer", "new review requested")
        store.save(store_path)
        loaded = CaseReviewStore.load(store_path)
        self.assertEqual(loaded.cases()[case]["state"], "reopened")
        self.assertEqual(loaded.cases()[case]["report"], before)
        self.assertEqual(loaded.cases()[case]["latest_report"], after)

    def test_case_rejects_unknown_evidence_and_stale_writer(self):
        path = self.root / "cases.json"
        a, b = CaseReviewStore.load(path), CaseReviewStore.load(path)
        case = a.open(report_dict(self.report()), "reviewer")
        a.transition(case, "in_review", "reviewer", "started")
        with self.assertRaises(ValueError):
            a.transition(case, "resolved", "reviewer", "done", ["invented"], "complete")
        a.save(path)
        with self.assertRaises(ValueError):
            b.save(path)

    def test_policy_content_change_detected_even_with_same_version(self):
        before = report_dict(self.report())
        after = report_dict(self.report(phenotypes=[phenotype()]))
        diff = compare_reports(before, after)
        self.assertIn("phenotypes", diff["manifest_changes"])
        self.assertTrue(any(c["call_changed"] for c in diff["changes"]))

    def test_html_escapes_laboratory_text(self):
        """Laboratory text must never reach the page as markup.

        The report now carries its own inline runtime, so "no <script> tag
        anywhere" is not the property to assert — it would be unsatisfiable
        and would say nothing about injection. What matters is that the
        *supplied* payload appears only in escaped form, on both routes into
        the page: the HTML body, and the JSON blob the runtime reads.
        """
        payload = "<script>alert(1)</script>"
        text = render_html(self.report(phenotypes=[phenotype(method=payload)]))

        # The injected markup never appears as markup.
        self.assertNotIn(payload, text)
        self.assertNotIn("alert(1)</script>", text)

        # The page's own runtime is still there and is the only script source.
        self.assertIn("<script>", text)

        # A "<" inside the injected JSON is unicode-escaped, so the payload
        # cannot close the script element early and become executable.
        start = text.index("window.__MYCONDUCTOR__")
        blob = text[start:text.index("</script>", start)]
        self.assertNotIn("<", blob)

    def test_fhir_references_resolve_to_unique_resources(self):
        bundle = to_fhir_bundle(self.report())
        urls = {e["fullUrl"] for e in bundle["entry"]}
        self.assertEqual(len(urls), len(bundle["entry"]))
        def walk(value):
            if isinstance(value, dict):
                if "reference" in value:
                    self.assertIn(value["reference"], urls)
                for v in value.values():
                    walk(v)
            elif isinstance(value, list):
                for v in value:
                    walk(v)
        walk(bundle)


class FollowUpTests(unittest.TestCase):
    def test_shared_action_charged_once_and_budget_respected(self):
        findings = [dict(id=str(i), drug="linezolid", priority=i, candidate_actions=["confirm"]) for i in range(2)]
        menu = [dict(id="dst", action="confirm", available=True, cost=5, currency="NGN")]
        result = plan_follow_up(findings, menu, 5)
        self.assertEqual([r["charged_cost"] for r in result], [5, 0])
        self.assertEqual(result[-1]["remaining_budget"], 0)
        self.assertEqual(plan_follow_up(findings, menu, 4)[0]["status"], "deferred_budget")

    def test_unknown_cost_is_not_zero(self):
        findings = [dict(id="1", drug="linezolid", priority=0, candidate_actions=["confirm"])]
        menu = [dict(id="dst", action="confirm", available=True)]
        self.assertEqual(plan_follow_up(findings, menu, 5)[0]["status"], "cost_unknown")
        self.assertIsNone(plan_follow_up(findings, menu)[0]["charged_cost"])

    def test_invalid_cost_and_mixed_currency_rejected(self):
        for menu in ([dict(id="1", action="a", available=True, cost=-1, currency="NGN")],
                     [dict(id=str(i), action="a", available=True, cost=1, currency=c) for i,c in enumerate(("NGN","USD"))]):
            with self.assertRaises(ValueError):
                plan_follow_up([], menu, 5)


class EvaluationTests(unittest.TestCase):
    def row(self, **kwargs):
        row = dict(site_id="s", isolate_id="i", drug="linezolid", truth="resistant",
                   truth_quality="pass", truth_method="DST", truth_source="lab-record",
                   training_overlap="no", independent_adjudication="yes", partition="evaluation",
                   baseline_call="susceptible", intervention_call="indeterminate")
        row.update(kwargs)
        return row

    def test_abstention_stays_in_denominator(self):
        result = evaluate([self.row()])
        metrics = result["primary"]["intervention"]
        self.assertEqual(metrics["n"], 1)
        self.assertEqual(metrics["unresolved_rate"], 1)
        self.assertEqual(metrics["resistant_detection_yield"], 0)
        self.assertIsNone(metrics["conditional_accuracy"])
        self.assertEqual(result["operational"]["cost"], {})

    def test_development_unknown_overlap_and_unblinded_truth_excluded(self):
        for key, value in (("partition", "development"), ("partition", ""),
                           ("training_overlap", "unknown"), ("independent_adjudication", "no")):
            self.assertEqual(evaluate([self.row(**{key: value})])["n_independent_evaluation"], 0)

    def test_cluster_leakage_and_duplicate_pairs_rejected(self):
        for rows in ([self.row(), self.row()],
                     [self.row(cluster_id="c"), self.row(isolate_id="j", cluster_id="c", partition="development")]):
            with self.assertRaises(ValueError):
                evaluate(rows)


class UpstreamTests(unittest.TestCase):
    def test_fixtures_match_pinned_hashes(self):
        root = Path(__file__).parent / "fixtures/upstream"
        for fixture in json.loads((root / "manifest.json").read_text()):
            self.assertEqual(hashlib.sha256((root / fixture["file"]).read_bytes()).hexdigest(), fixture["sha256"])

    def test_real_ntm_output_preserves_reference_gene_state(self):
        report = NTMProfilerAdapter().parse(Path(__file__).parent / "fixtures/upstream/ntm-results.json")
        self.assertEqual(report.sample_id, "ERR459870_fasta")
        self.assertEqual(report.context_findings["erm41_status"], "functional")
        self.assertEqual(report.context_findings["subspecies"], "abscessus")
        self.assertIsNone(report.mask)
        self.assertTrue(report.evidence)
        self.assertTrue(all(e.call == Call.INDETERMINATE for e in report.evidence))


class CohortGenotypeTests(unittest.TestCase):
    def test_reference_multiallelic_and_missing_genotypes(self):
        index = CoordinateIndex(by_allele={(10, "A", "C"): {"v1"}, (10, "A", "G"): {"v2"}})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "s.vcf"
            for gt, expected, assessed in (("0/0", set(), {"v1","v2"}), ("2", {"v2"}, {"v1","v2"}), (".", set(), set())):
                path.write_text(f"NC_000962.3\t10\t.\tA\tC,G\t60\tPASS\t.\tGT\t{gt}\n")
                result = parse_vcf(path, index)
                self.assertEqual(result.variants, expected)
                self.assertEqual(result.assessed_variants, assessed)

    def test_unassessed_is_not_a_control(self):
        isolates = [Isolate("carrier", frozenset({"v"}), assessed_variants=frozenset({"v"})),
                    Isolate("missing", assessed_variants=frozenset())]
        result = stratify("v", "linezolid", isolates, DeterminantIndex())
        self.assertEqual(result.n_non_carriers, 0)

    def test_cache_cannot_escape_directory_or_download_in_cached_mode(self):
        with tempfile.TemporaryDirectory() as root:
            for path in ("a/../../outside.vcf", "missing.vcf"):
                with self.assertRaises(GenotypeError):
                    fetch_vcf(path, root, cached_only=True)

    def test_metadata_cluster_selection_and_leakage(self):
        isolates = [Isolate("a"), Isolate("b")]
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "metadata.tsv"
            path.write_text("isolate_id\tsite_id\tpatient_id\tcluster_id\tpartition\n"
                            "a\ts\tp1\tc\tevaluation\nb\ts\tp2\tc\tevaluation\n")
            selected, _ = apply_metadata(isolates, path, "evaluation", True)
            self.assertEqual([i.isolate_id for i in selected], ["a"])
            path.write_text(path.read_text().replace("p2\tc\tevaluation", "p2\tc\tdevelopment"))
            with self.assertRaises(ValueError):
                apply_metadata(isolates, path)
