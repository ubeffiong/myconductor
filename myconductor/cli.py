"""Command-line interface for Myconductor.

    myconductor analyze <input> [--mask PATH | --bam PATH --bed PATH]
                        [--platform NAME] [--sample NAME] [--fhir OUT]
                        [--depth-floor N] [--callable-floor F]
                        [--tbprofiler JSON] [--mykrobe JSON]
                        [--amrfinder TSV] [--local-validation-store PATH]
    myconductor demo
    myconductor ingest-catalogue <csv|tsv|xlsx> --out catalogue.json
    myconductor validate-vus --store PATH --variant-key K --gene G --drug D
                        --isolate-id ID --site-id SITE --method M --result R
    myconductor governance-report --site-calls SITE=path.json [SITE=path.json ...]
    myconductor call-variants --fastq1 R1 --fastq2 R2 --reference ref.fa --out out.vcf.gz
    myconductor version
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .catalogue.profile import BUNDLED_ORGANISMS
from .core.pipeline import Myconductor
from .reporting.fhir import to_fhir_json
from .reporting.render import render_text

_DATA = Path(__file__).resolve().parent / "data"
_DEMO_INPUT = _DATA / "example_input.vcf"
_DEMO_MASK = _DATA / "example_callable.tsv"


def _engine_reports(args: argparse.Namespace) -> list:
    reports = []
    if getattr(args, "tbprofiler", None):
        from .adapters.tbprofiler import TBProfilerAdapter
        reports.append(TBProfilerAdapter().parse(args.tbprofiler))
    if getattr(args, "mykrobe", None):
        from .adapters.mykrobe import MykrobeAdapter
        reports.append(MykrobeAdapter().parse(args.mykrobe))
    if getattr(args, "amrfinder", None):
        from .adapters.amrfinderplus import AMRFinderPlusAdapter
        report = AMRFinderPlusAdapter().parse(args.amrfinder)
        if not args.amrfinder_sample:
            raise ValueError("AMRFinder TSV requires --amrfinder-sample")
        report.sample_id = args.amrfinder_sample
        reports.append(report)
    if getattr(args, "ntm_profiler", None):
        from .adapters.ntm_profiler import NTMProfilerAdapter
        reports.append(NTMProfilerAdapter().parse(args.ntm_profiler))
    if getattr(args, "gnomonicus", None):
        from .adapters.gnomonicus import GnomonicusAdapter
        reports.append(GnomonicusAdapter(sample_id=args.gnomonicus_sample).parse(args.gnomonicus))
    return reports


def _run_analyze(args: argparse.Namespace) -> int:
    if bool(args.compare_catalogue) != bool(args.change_impact):
        raise ValueError("catalogue replay requires both --compare-catalogue and --change-impact")
    if args.case_store and not args.reviewer:
        raise ValueError("--case-store requires --reviewer")
    if args.bam and not args.bed:
        raise ValueError("--bam requires --bed (the locus regions to assess)")
    if args.bed and not args.bam:
        raise ValueError("--bed requires --bam")
    if args.bam and args.mask:
        raise ValueError("pass either --mask or --bam/--bed, not both")

    mask = None
    if args.bam:
        from .io.bam_coverage import mask_from_bam
        mask = mask_from_bam(args.bam, args.bed, depth_floor=args.depth_floor,
                             tool=args.coverage_tool)

    local_validation_store = None
    if args.local_validation_store:
        from .federated.vus_feedback import LocalValidationStore
        local_validation_store = LocalValidationStore.load(
            args.local_validation_store)

    from .core.context import SampleContext, InterpretationPolicy
    from .catalogue.profile import load_profile
    from .io.phenotypes import load_phenotypes
    from .reporting.json_report import report_dict, atomic_json
    context = SampleContext.from_dict(json.loads(Path(args.context).read_text(encoding="utf-8"))) if args.context else None
    policy = InterpretationPolicy.from_dict(json.loads(Path(args.policy).read_text(encoding="utf-8"))) if args.policy else None
    if bool(args.drug_loci) != bool(args.drugs):
        raise ValueError("custom profiles require both --drug-loci and --drugs")
    profile = load_profile(Path(args.drug_loci), Path(args.drugs)) if args.drug_loci else None
    conductor = Myconductor(
        profile=profile, catalogue_path=args.catalogue, policy=policy,
        depth_floor=args.depth_floor,
        callable_fraction_floor=args.callable_floor,
        platform=args.platform,
        local_validation_store=local_validation_store,
        organism=args.organism,
    )
    analysis_args = dict(
        mask=mask,
        mask_path=args.mask,
        panel_path=args.panel,
        sample=args.sample,
        engine_reports=_engine_reports(args), context=context,
        phenotypes=load_phenotypes(args.phenotypes) if args.phenotypes else (),
        test_menu=json.loads(Path(args.test_menu).read_text(encoding="utf-8")) if args.test_menu else (),
        follow_up_budget=args.follow_up_budget,
    )
    report = conductor.analyze(args.input, **analysis_args)
    if args.compare_catalogue:
        from .federated.change_impact import compare_reports
        replay = Myconductor(
            profile=conductor.profile, catalogue_path=args.compare_catalogue, policy=policy,
            depth_floor=args.depth_floor, callable_fraction_floor=args.callable_floor,
            platform=args.platform, local_validation_store=local_validation_store,
        ).analyze(args.input, **analysis_args)
        impact_path = Path(args.change_impact)
        before_path = impact_path.with_name(impact_path.stem + ".before.json")
        after_path = impact_path.with_name(impact_path.stem + ".after.json")
        before_data, after_data = report_dict(report), report_dict(replay)
        atomic_json(before_path, before_data)
        atomic_json(after_path, after_data)
        impact = compare_reports(before_data, after_data)
        impact["replay_reports"] = {"before": str(before_path), "after": str(after_path)}
        atomic_json(impact_path, impact)
    print(render_text(report))

    if args.json:
        atomic_json(args.json, report_dict(report))
    if args.html:
        from .reporting.html_report import render_html
        Path(args.html).parent.mkdir(parents=True, exist_ok=True)
        Path(args.html).write_text(render_html(report), encoding="utf-8")
    if args.case_store:
        if not args.reviewer:
            raise ValueError("--case-store requires --reviewer")
        from .federated.case_review import CaseReviewStore
        store = CaseReviewStore.load(args.case_store)
        case_id = store.open(report_dict(report), args.reviewer)
        store.save(args.case_store)
        print(f"Case opened: {case_id}")

    if args.fhir:
        Path(args.fhir).parent.mkdir(parents=True, exist_ok=True)
        Path(args.fhir).write_text(to_fhir_json(report), encoding="utf-8")
        print(f"\n[FHIR bundle written to {args.fhir}]")

    if args.jsonld:
        from .reporting.semantic import summarise, to_jsonld, write_jsonld

        written = write_jsonld(report, args.jsonld)
        print(f"\n[JSON-LD written to {written}]")
        print(f"  {summarise(to_jsonld(report))}")
        print("  Six call states, not two: a drug may be used only where "
              "establishesUse is true.")
    return 0


def _run_demo(args: argparse.Namespace) -> int:
    print(f"Myconductor demo — {_DEMO_INPUT.name} with {_DEMO_MASK.name}\n")
    print("The bundled input and catalogue are illustrative. This demonstrates "
          "the evidence flow, not a validated analysis.\n")
    report = Myconductor(platform="illumina", demo_mode=True).analyze(
        _DEMO_INPUT, mask_path=_DEMO_MASK)
    print(render_text(report))
    return 0


def _run_ingest(args: argparse.Namespace) -> int:
    path = Path(args.catalogue)
    if path.suffix.lower() == ".xlsx":
        from .adapters.who_catalogue import ingest_xlsx
        result = ingest_xlsx(path, sheet=args.sheet, catalogue_version=args.version)
    else:
        from .adapters.who_catalogue import ingest_csv
        result = ingest_csv(path, catalogue_version=args.version)

    out = result.write(args.out)
    print(f"Ingested {result.n_entries} entries across "
          f"{len(result.drugs)} drug(s) -> {out}")
    if result.skipped:
        print(f"Skipped {len(result.skipped)} row(s). First few:")
        for line in result.skipped[:5]:
            print(f"  - {line}")
    return 0


def _run_validate_vus(args: argparse.Namespace) -> int:
    from .core.models import Call
    from .federated.vus_feedback import (
        LocalValidationStore,
        ValidationMethod,
        VUSValidationRecord,
    )

    store_path = Path(args.store)
    store = LocalValidationStore.load(store_path) if store_path.is_file() \
        else LocalValidationStore()

    record = VUSValidationRecord(
        variant_key=args.variant_key,
        variant_label=args.variant_label or args.variant_key,
        gene=args.gene, drug=args.drug, isolate_id=args.isolate_id,
        site_id=args.site_id, lineage=args.lineage,
        method=ValidationMethod(args.method), result=Call(args.result),
        mic=args.mic, submitted_by=args.submitted_by or "unknown",
        rationale=args.rationale or "",
    )
    ledger_entry = store.ingest(record)
    store.save(store_path)

    verdict = store.verdict(record.variant_key, record.drug)
    print(f"Ingested: {record.variant_key} / {record.drug} = "
         f"{record.result.value} (isolate {record.isolate_id}, site "
         f"{record.site_id})")
    print(f"Ledger entry #{ledger_entry['seq']} "
         f"({ledger_entry['entry_hash'][:16]}...)")
    print(f"Current site verdict: {verdict.describe()}")
    ok, detail = store.ledger_verified()
    print(f"Ledger integrity: {'OK' if ok else f'FAILED — {detail}'}")
    print(f"\nThis is a site-local record, not a global catalogue entry. "
         f"Wrote {store_path}")
    return 0


def _run_call_variants(args: argparse.Namespace) -> int:
    from .io.read_calling import align_reads, call_variants, reads_to_vcf

    if args.bam_out:
        bam = align_reads(args.fastq1, args.fastq2, args.reference,
                          args.bam_out, aligner=args.aligner,
                          threads=args.threads)
        print(f"Wrote alignment: {bam}")
        vcf = call_variants(bam, args.reference, args.out, caller=args.caller,
                            region_bed=args.region_bed)
    else:
        vcf = reads_to_vcf(args.fastq1, args.fastq2, args.reference, args.out,
                           aligner=args.aligner, caller=args.caller,
                           threads=args.threads, region_bed=args.region_bed)
    print(f"Wrote variants: {vcf}")
    print("\nThis VCF was produced by an unvalidated orchestration of "
         "external tools (see myconductor/io/read_calling.py) — treat it "
         "the same as output from any variant caller you have not yet "
         "benchmarked. Feed it to 'myconductor analyze' like any other VCF.")
    return 0


def _run_governance_report(args: argparse.Namespace) -> int:
    from .federated.catalogue_governance import SiteCallRecord, reconcile_sites

    records = []
    for entry in args.site_calls:
        site_id, _, path = entry.partition("=")
        if not path:
            raise ValueError(f"--site-calls expects SITE_ID=path.json; got {entry!r}")
        rows = json.loads(Path(path).read_text(encoding="utf-8"))
        for row in rows:
            row.setdefault("site_id", site_id)
            records.append(SiteCallRecord.from_dict(row))

    discordances = reconcile_sites(records)
    print(f"{len(records)} record(s) from {len(args.site_calls)} site(s); "
         f"{len(discordances)} discordance(s) found.\n")
    for d in discordances:
        print(f"[{d.category.value}] {d.variant_key} / {d.drug}")
        for site, call in sorted(d.calls_by_site.items()):
            version = d.catalogue_versions_by_site.get(site) or "unknown"
            print(f"    {site}: {call} (catalogue {version})")
        print(f"    -> {d.note}\n")
    if not discordances:
        print("Every site agrees on every shared variant/drug pair.")
    return 0


def _run_case_review(args):
    from .federated.case_review import CaseReviewStore
    store = CaseReviewStore.load(args.store)
    if args.action == "list":
        print(json.dumps(store.cases(), indent=2))
        return 0
    if not args.case_id or not args.reviewer or not args.rationale:
        raise ValueError("case transitions require --case-id, --reviewer and --rationale")
    if args.action == "attach":
        from .reporting.json_report import load_report
        if not args.report:
            raise ValueError("attach requires --report")
        store.attach(args.case_id, load_report(args.report), args.reviewer, args.rationale)
    else:
        store.transition(args.case_id, args.action, args.reviewer, args.rationale,
                         args.evidence_id or (), args.resolution, args.reviewer_minutes)
    store.save(args.store)
    print(f"Case {args.case_id}: {args.action}")
    return 0


def _run_change_impact(args):
    from .federated.change_impact import compare_reports
    from .reporting.json_report import load_report, atomic_json
    result = compare_reports(load_report(args.before), load_report(args.after))
    atomic_json(args.out, result)
    print(f"{len(result['changes'])} drug evidence change(s); wrote {args.out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="myconductor",
        description="Evidence-orchestration layer for AMR genomics. "
                    "Research scaffold; not a clinical device.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    a = sub.add_parser("analyze", help="Analyze variants (.vcf/.tsv/.json).")
    a.add_argument("input", help="Path to input variants.")
    a.add_argument("--compare-catalogue", help="Replay identical inputs against this second catalogue.")
    a.add_argument("--change-impact", help="Replay comparison JSON; also freezes before/after reports.")
    a.add_argument("--mask", metavar="PATH",
                   help="Callable-locus evidence (.tsv depth table or .bed "
                        "mask). Without it, no drug can be reported "
                        "susceptible.")
    a.add_argument("--panel", metavar="PATH",
                   help="Targeted-panel declaration (JSON). Coverage claimed "
                        "for a locus the panel does not target is discarded: "
                        "a targeted assay cannot call a locus it never "
                        "amplified, whatever a coverage file says.")
    a.add_argument("--bam", metavar="PATH",
                   help="Derive the callable-locus mask directly from this "
                        "BAM via mosdepth/samtools (requires --bed). "
                        "Coverage only — does not call variants; mutually "
                        "exclusive with --mask.")
    a.add_argument("--bed", metavar="PATH",
                   help="Locus regions (chrom/start/end/name) for --bam.")
    a.add_argument("--coverage-tool", choices=("auto", "mosdepth", "samtools"),
                   default="auto", help="Tool to derive coverage from --bam "
                                       "(default: auto-detect).")
    a.add_argument("--local-validation-store", metavar="PATH",
                   help="A site-local VUS validation store (see "
                        "'validate-vus') to consult during this analysis.")
    a.add_argument("--organism", choices=sorted(BUNDLED_ORGANISMS),
                   help="Bundled organism profile + catalogue to use "
                        "(default: mtbc). Every non-mtbc profile is "
                        "illustrative in the same sense the bundled MTBC "
                        "catalogue is -- see docs/ARCHITECTURE.md.")
    a.add_argument("--platform", help="Sequencing platform (illumina, "
                                      "nanopore, pacbio, targeted-amplicon); "
                                      "required for minority-allele limits of "
                                      "detection.")
    a.add_argument("--sample", help="Sample name, for multi-sample VCFs.")
    a.add_argument("--fhir", metavar="PATH", help="Also write a FHIR bundle.")
    a.add_argument("--jsonld", metavar="PATH",
                   help="Also write a JSON-LD export carrying the six call "
                        "states and their semantics, for consumers outside "
                        "this system. Refuses to emit any binary "
                        "resistant/susceptible field.")
    a.add_argument("--depth-floor", type=int, default=10,
                   help="Minimum read depth to accept a locus (default 10).")
    a.add_argument("--callable-floor", type=float, default=0.95,
                   help="Minimum callable fraction of a locus (default 0.95).")
    a.add_argument("--tbprofiler", metavar="JSON",
                   help="TB-Profiler results JSON to reconcile.")
    a.add_argument("--mykrobe", metavar="JSON",
                   help="Mykrobe predict JSON to reconcile.")
    a.add_argument("--amrfinder", metavar="TSV",
                   help="AMRFinderPlus TSV to reconcile.")
    a.add_argument("--context", help="Specimen/assay context JSON; IDs must match input.")
    a.add_argument("--policy", help="Versioned, externally reviewed validation scopes JSON.")
    a.add_argument("--catalogue", help="Explicit ingested catalogue JSON (default: illustrative).")
    a.add_argument("--drug-loci", help="Custom organism drug-locus profile JSON.")
    a.add_argument("--drugs", help="Custom organism drug definitions JSON.")
    a.add_argument("--phenotypes", help="Matched laboratory phenotype observations JSON.")
    a.add_argument("--test-menu", help="Local follow-up options, availability and costs JSON.")
    a.add_argument("--follow-up-budget", type=float, help="Budget in the menu's single currency.")
    a.add_argument("--json", help="Write versioned evidence/review JSON.")
    a.add_argument("--html", help="Write an offline HTML review report.")
    a.add_argument("--case-store", help="Open an auditable local review case.")
    a.add_argument("--reviewer", help="Reviewer identity for opening a case.")
    a.add_argument("--ntm-profiler", help="NTM-Profiler results JSON.")
    a.add_argument("--gnomonicus", help="gnomonicus JSON output.")
    a.add_argument("--gnomonicus-sample", help="Explicit sample binding for gnomonicus output.")
    a.add_argument("--amrfinder-sample", help="Explicit sample binding for AMRFinder TSV.")
    a.set_defaults(func=_run_analyze)

    d = sub.add_parser("demo", help="Run the bundled illustrative demo.")
    d.set_defaults(func=_run_demo)

    i = sub.add_parser("ingest-catalogue",
                       help="Ingest a WHO catalogue export into JSON.")
    i.add_argument("catalogue", help="CSV, TSV or XLSX export.")
    i.add_argument("--out", required=True, help="Output catalogue JSON path.")
    i.add_argument("--sheet", help="Worksheet name, for .xlsx.")
    i.add_argument("--version", help="Catalogue version label to record.")
    i.set_defaults(func=_run_ingest)

    vv = sub.add_parser(
        "validate-vus",
        help="Ingest a laboratory VUS validation result into a site-local "
             "store (the closed loop back from the VUS workbench).")
    vv.add_argument("--store", required=True, metavar="PATH",
                    help="JSON file for the local validation store; created "
                         "if absent, updated in place otherwise.")
    vv.add_argument("--variant-key", required=True,
                    help="Variant identity, e.g. from Variant.identity.key().")
    vv.add_argument("--variant-label",
                    help="Display label (gene_change); defaults to "
                         "--variant-key if omitted.")
    vv.add_argument("--gene", required=True)
    vv.add_argument("--drug", required=True)
    vv.add_argument("--isolate-id", required=True)
    vv.add_argument("--site-id", required=True)
    vv.add_argument("--lineage")
    vv.add_argument("--method", required=True,
                    choices=("phenotypic_dst", "mic", "efflux_assay",
                            "allelic_exchange", "rna_expression", "other"))
    vv.add_argument("--result", required=True,
                    choices=("resistant", "susceptible"),
                    help="The validation's established outcome. An "
                         "inconclusive result is not accepted here.")
    vv.add_argument("--mic", type=float)
    vv.add_argument("--submitted-by")
    vv.add_argument("--rationale")
    vv.set_defaults(func=_run_validate_vus)

    cv = sub.add_parser(
        "call-variants",
        help="Align FASTQ reads and call variants (unvalidated orchestration "
             "of external tools; see io/read_calling.py). Feed the resulting "
             "VCF to 'analyze' like any other VCF.")
    cv.add_argument("--fastq1", required=True, metavar="R1")
    cv.add_argument("--fastq2", required=True, metavar="R2")
    cv.add_argument("--reference", required=True, metavar="PATH",
                    help="Reference FASTA.")
    cv.add_argument("--out", required=True, metavar="VCF",
                    help="Output VCF path.")
    cv.add_argument("--bam-out", metavar="PATH",
                    help="Also keep the intermediate BAM here (e.g. to "
                         "derive a coverage mask with 'analyze --bam').")
    cv.add_argument("--aligner", choices=("auto", "minimap2", "bwa-mem2"),
                    default="auto")
    cv.add_argument("--caller", choices=("auto", "bcftools", "gatk"),
                    default="auto")
    cv.add_argument("--region-bed", metavar="PATH",
                    help="Restrict calling to these regions.")
    cv.add_argument("--threads", type=int, default=4)
    cv.set_defaults(func=_run_call_variants)

    gr = sub.add_parser(
        "governance-report",
        help="Reconcile per-site calls, distinguishing catalogue-version "
             "discordance from pipeline/rule discordance.")
    gr.add_argument("--site-calls", action="append", required=True,
                    metavar="SITE_ID=path.json",
                    help="A site's exported SiteCallRecord list (see "
                         "federated.catalogue_governance.site_call_records). "
                         "Repeat for each site.")
    gr.set_defaults(func=_run_governance_report)

    cr = sub.add_parser("review-case", help="List or transition a local review case.")
    cr.add_argument("--store", required=True)
    cr.add_argument("--action", choices=("list", "attach", "in_review", "resolved", "reopened"), required=True)
    cr.add_argument("--report", help="Frozen analysis report to attach.")
    cr.add_argument("--case-id")
    cr.add_argument("--reviewer")
    cr.add_argument("--rationale")
    cr.add_argument("--evidence-id", action="append")
    cr.add_argument("--resolution")
    cr.add_argument("--reviewer-minutes", type=float)
    cr.set_defaults(func=_run_case_review)
    ci = sub.add_parser("change-impact", help="Compare frozen reports for one input sample.")
    ci.add_argument("--before", required=True)
    ci.add_argument("--after", required=True)
    ci.add_argument("--out", required=True)
    ci.set_defaults(func=_run_change_impact)

    v = sub.add_parser("version", help="Print version.")
    v.set_defaults(func=lambda _a: (print(f"Myconductor {__version__}"), 0)[1])

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, OSError, TypeError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
