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
        reports.append(AMRFinderPlusAdapter().parse(args.amrfinder))
    return reports


def _run_analyze(args: argparse.Namespace) -> int:
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

    conductor = Myconductor(
        depth_floor=args.depth_floor,
        callable_fraction_floor=args.callable_floor,
        platform=args.platform,
        local_validation_store=local_validation_store,
        organism=args.organism,
    )
    report = conductor.analyze(
        args.input,
        mask=mask,
        mask_path=args.mask,
        sample=args.sample,
        engine_reports=_engine_reports(args),
    )
    print(render_text(report))

    if not args.mask and not args.bam:
        print("\n[note] no --mask or --bam supplied, so no drug can be "
              "reported susceptible. This is intended: susceptibility "
              "requires evidence that the loci were sequenced.",
              file=sys.stderr)

    if args.fhir:
        Path(args.fhir).write_text(to_fhir_json(report))
        print(f"\n[FHIR bundle written to {args.fhir}]")
    return 0


def _run_demo(args: argparse.Namespace) -> int:
    print(f"Myconductor demo — {_DEMO_INPUT.name} with {_DEMO_MASK.name}\n")
    print("The bundled input and catalogue are illustrative. This demonstrates "
          "the evidence flow, not a validated analysis.\n")
    report = Myconductor(platform="illumina").analyze(
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="myconductor",
        description="Evidence-orchestration layer for AMR genomics. "
                    "Research scaffold; not a clinical device.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    a = sub.add_parser("analyze", help="Analyze variants (.vcf/.tsv/.json).")
    a.add_argument("input", help="Path to input variants.")
    a.add_argument("--mask", metavar="PATH",
                   help="Callable-locus evidence (.tsv depth table or .bed "
                        "mask). Without it, no drug can be reported "
                        "susceptible.")
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

    v = sub.add_parser("version", help="Print version.")
    v.set_defaults(func=lambda _a: (print(f"Myconductor {__version__}"), 0)[1])

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
