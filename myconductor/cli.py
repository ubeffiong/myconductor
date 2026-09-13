"""Command-line interface for Myconductor.

    myconductor analyze <input> [--mask PATH] [--platform NAME]
                        [--sample NAME] [--fhir OUT] [--depth-floor N]
                        [--callable-floor F] [--tbprofiler JSON]
                        [--mykrobe JSON] [--amrfinder TSV]
    myconductor demo
    myconductor ingest-catalogue <csv|tsv|xlsx> --out catalogue.json
    myconductor version
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
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
    conductor = Myconductor(
        depth_floor=args.depth_floor,
        callable_fraction_floor=args.callable_floor,
        platform=args.platform,
    )
    report = conductor.analyze(
        args.input,
        mask_path=args.mask,
        sample=args.sample,
        engine_reports=_engine_reports(args),
    )
    print(render_text(report))

    if not args.mask:
        print("\n[note] no --mask supplied, so no drug can be reported "
              "susceptible. This is intended: susceptibility requires "
              "evidence that the loci were sequenced.", file=sys.stderr)

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
