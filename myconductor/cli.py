"""Command-line interface for Myconductor.

    python -m myconductor analyze <input> [--fhir out.json] [--depth-floor 10]
    python -m myconductor demo
    python -m myconductor version
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .core.pipeline import Myconductor
from .reporting.fhir import to_fhir_json
from .reporting.render import render_text

_DEMO_INPUT = Path(__file__).resolve().parent / "data" / "example_input.vcf"


def _run_analyze(args: argparse.Namespace) -> int:
    conductor = Myconductor(depth_floor=args.depth_floor)
    report = conductor.analyze(args.input)
    print(render_text(report))
    if args.fhir:
        Path(args.fhir).write_text(to_fhir_json(report))
        print(f"\n[FHIR bundle written to {args.fhir}]")
    return 0


def _run_demo(_args: argparse.Namespace) -> int:
    print(f"Running Myconductor demo on {_DEMO_INPUT.name}\n")
    conductor = Myconductor()
    report = conductor.analyze(_DEMO_INPUT)
    print(render_text(report))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="myconductor",
        description="Orchestration-and-reasoning layer for TB antimicrobial resistance.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    a = sub.add_parser("analyze", help="Analyze a variant file (.vcf/.tsv/.json).")
    a.add_argument("input", help="Path to input variants.")
    a.add_argument("--fhir", metavar="PATH", help="Also write a FHIR bundle here.")
    a.add_argument("--depth-floor", type=int, default=10,
                   help="Minimum read depth to accept a locus (default 10).")
    a.set_defaults(func=_run_analyze)

    d = sub.add_parser("demo", help="Run the bundled end-to-end demo.")
    d.set_defaults(func=_run_demo)

    v = sub.add_parser("version", help="Print version.")
    v.set_defaults(func=lambda _a: (print(f"Myconductor {__version__}"), 0)[1])

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
