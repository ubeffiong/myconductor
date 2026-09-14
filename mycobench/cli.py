"""Command-line interface for mycobench.

    mycobench discover-cohort   --organism "..." [--country Nigeria] --out-dir DIR
    mycobench validate-cohort   --samples S [--online] [--write-lock L] [--verify-lock L]
    mycobench review-candidates --candidates C --out-dir DIR [--max-per-bioproject N]
    mycobench estimate-download --samples S [--data-dir D]
    mycobench fetch-catalogue   --out-dir DIR [--write-profile]
    mycobench fetch-phenotypes  --out-dir DIR [--cohort-out C]
    mycobench run               --cohort NAME [--stages ...] [--write-script F]
    mycobench report            --results-dir R [--out report.html]
    mycobench version
"""
from __future__ import annotations

import argparse
import os
import shlex
import sys
from pathlib import Path

from . import __version__

COHORTS_DIR = Path("cohorts")


def _resolve_samples(args) -> Path:
    if getattr(args, "cohort", None):
        path = COHORTS_DIR / f"{args.cohort}.tsv"
        if not path.is_file():
            available = sorted(p.stem for p in COHORTS_DIR.glob("*.tsv"))
            raise SystemExit(
                f"error: no cohort named {args.cohort!r} "
                f"({path} not found). Available: {', '.join(available) or 'none'}")
        return path
    if getattr(args, "samples", None):
        return Path(args.samples)
    raise SystemExit("error: pass --cohort NAME or --samples PATH")


# -- discover -------------------------------------------------------------
def _run_discover(args) -> int:
    from .cohort import write_rows
    from .discover import ACCEPTED_COLUMNS, REJECTED_COLUMNS, discover
    from .ncbi import Credentials

    creds = Credentials.from_env(args.email, args.api_key)
    print(f"[discover] {creds.describe()}")
    if args.country:
        print(f"[discover] origin will be CONFIRMED against BioSample "
              f"geo_loc_name, not free text")

    result = discover(
        args.organism, creds, country=args.country, max_runs=args.max_runs,
        prefix=args.prefix, require_mtbc=not args.non_mtbc,
        expected_outcome=("species_mismatch_refused" if args.non_mtbc else ""),
        min_bases=args.min_bases)

    out = Path(args.out_dir)
    write_rows(out / "accepted.tsv", result.accepted, ACCEPTED_COLUMNS)
    write_rows(out / "rejected.tsv", result.rejected, REJECTED_COLUMNS)
    print(result.summary())
    print(f"\nWrote {out / 'accepted.tsv'} and {out / 'rejected.tsv'}")
    print("Acceptance proves deposited metadata linkage only. Review each "
          "isolate's publication before any release claim.")
    return 0


# -- validate -------------------------------------------------------------
def _run_validate(args) -> int:
    from .cohort import (
        CohortError, read_rows, schema_errors, summarise, verify_lock,
        verify_rows, write_lock,
    )
    from .ncbi import Credentials

    samples = _resolve_samples(args)
    rows, fields = read_rows(samples)
    errors = schema_errors(rows, fields)

    if args.verify_lock:
        try:
            n = verify_lock(args.verify_lock, samples)
            print(f"COHORT LOCK VERIFIED: {n} evidence record(s)")
        except CohortError as exc:
            errors.append(str(exc))

    evidence = []
    if args.online and not errors:
        creds = Credentials.from_env(args.email, args.api_key)
        print(f"[validate] {creds.describe()}")
        print(f"[validate] verifying {len(rows)} row(s) against NCBI ...")
        evidence = verify_rows(rows, creds, country=args.country)
        for record in evidence:
            errors.extend(f"{record.sample_id}: {message}"
                          for message in record.errors)

    if errors:
        raise SystemExit("COHORT VALIDATION FAILED\n" + "\n".join(errors))

    if args.write_lock:
        if not args.online:
            raise SystemExit("error: --write-lock requires --online")
        path = write_lock(args.write_lock, samples, evidence, args.country)
        print(f"Wrote cohort verification lock: {path}")

    stats = summarise(rows)
    scope = "NCBI-verified" if args.online else "schema only"
    print(f"COHORT VALIDATION PASSED: {stats['n_rows']} row(s) ({scope})")
    print(f"  BioProjects: {stats['n_bioprojects']}, largest share "
          f"{stats['largest_bioproject_share']:.0%}")
    print(f"  controls: {stats['n_controls']}, with phenotypes: "
          f"{stats['n_phenotyped']}")
    if stats["largest_bioproject_share"] > 0.30:
        print("  NOTE: one BioProject supplies more than 30% of this panel. "
              "Run review-candidates before any interval or significance "
              "claim.")
    if not stats["n_phenotyped"]:
        print("  NOTE: no row carries a phenotype, so this panel supports the "
              "concordance track only; accuracy cannot be measured from it.")
    return 0


# -- review ---------------------------------------------------------------
def _run_review(args) -> int:
    from .cohort import ALL_COLUMNS, write_rows
    from .review import DEPENDENCE_COLUMNS, dependence_rows, review_file

    result = review_file(args.candidates, args.max_per_bioproject,
                         args.max_per_organism, args.max_total)
    out = Path(args.out_dir)
    write_rows(out / "balanced_shortlist.pending_review.tsv",
               result.shortlist, ALL_COLUMNS)
    write_rows(out / "held_back.tsv", result.held_back, ALL_COLUMNS)
    rows = []
    for dependence in result.dependence.values():
        rows.extend(dependence_rows(dependence))
    write_rows(out / "study_dependence.tsv", rows, DEPENDENCE_COLUMNS)
    print(result.summary())
    print(f"\nWrote {out}/ (shortlist, held_back, study_dependence)")
    print("Nothing was deleted: a cap is a statement about independence, not "
          "about data quality.")
    return 0


# -- estimate -------------------------------------------------------------
def _run_estimate(args) -> int:
    from .estimate import estimate

    samples = _resolve_samples(args)
    print(estimate(samples, args.data_dir).render(samples))
    return 0


# -- catalogue ------------------------------------------------------------
def _run_fetch_catalogue(args) -> int:
    from .catalogue import (
        PINNED_COMMIT, REPOSITORY, fetch_all, ingest_from_dir, write_outputs,
    )

    out = Path(args.out_dir)
    print(f"[catalogue] fetching from {REPOSITORY}@{PINNED_COMMIT[:8]} ...")
    fetched = fetch_all(out, include_workbook=args.workbook)
    for f in fetched.values():
        print(f"  {f.describe()}")

    print("[catalogue] ingesting ...")
    result = ingest_from_dir(out, fetched)
    catalogue_out = Path(args.catalogue_out)
    targets = write_outputs(
        result, catalogue_out,
        drug_loci_out=args.profile_out if args.write_profile else None,
        drug_loci_template=(Path("myconductor/catalogue/drug_loci.json")
                            if args.write_profile else None))
    print(f"  {result.n_entries} variant entries across "
          f"{len(result.drugs)} drug(s)")
    print(f"  {result.n_with_coordinates} entries carry genomic coordinates")
    if result.skipped:
        print(f"  skipped {len(result.skipped)} row(s); first: "
              f"{result.skipped[0]}")
    for path in targets:
        print(f"Wrote {path}")
    print("\nPoint Myconductor at it:")
    print(f"  myconductor analyze sample.vcf --mask mask.tsv  "
          f"# after replacing myconductor/catalogue/mtb_amr_catalogue.json")
    return 0


# -- phenotypes -----------------------------------------------------------
def _run_fetch_phenotypes(args) -> int:
    from .cohort import write_rows
    from .ncbi import Credentials, runinfo_for_accessions
    from .phenotypes import (
        COHORT_COLUMNS, EXCLUDED_TABLE, PHENOTYPE_COLUMNS, REUSE_TABLE,
        build_cohort_rows, fetch, phenotype_table, read_excluded, read_table,
        summarise,
    )

    out = Path(args.out_dir)
    print(f"[phenotypes] fetching {REUSE_TABLE} ...")
    provenance = fetch(out)
    for name, meta in provenance["files"].items():
        print(f"  {name} ({meta['size']:,} bytes, "
              f"sha256:{meta['sha256'][:12]})")
    print(f"  cite: {provenance['citation']}")

    excluded = read_excluded(out / EXCLUDED_TABLE)
    print(f"  {len(excluded)} withdrawn sample(s) will be dropped")
    isolates = read_table(out / REUSE_TABLE, excluded)
    stats = summarise(isolates)
    print(f"  {stats['n_isolates']} isolate(s) parsed")
    print("  evaluable phenotypes at accepted quality, per drug:")
    for drug, counts in sorted(stats["per_drug"].items()):
        if counts["evaluable"]:
            print(f"    {drug:<14} {counts['evaluable']:>6} "
                  f"({counts['resistant']} R / {counts['susceptible']} S)")

    write_rows(out / "cryptic_phenotypes.tsv", phenotype_table(isolates),
               PHENOTYPE_COLUMNS)
    print(f"Wrote {out / 'cryptic_phenotypes.tsv'}")

    if args.cohort_out:
        if args.limit:
            isolates = isolates[:args.limit]
        creds = Credentials.from_env(args.email, args.api_key)
        print(f"[phenotypes] resolving {len(isolates)} ENA run(s) at NCBI "
              f"({creds.describe()}) ...")
        info = runinfo_for_accessions([i.ena_run for i in isolates], creds)
        accepted, rejected = build_cohort_rows(isolates, info)
        write_rows(args.cohort_out, accepted, COHORT_COLUMNS)
        print(f"Wrote {args.cohort_out} ({len(accepted)} row(s), "
              f"{len(rejected)} rejected)")
        for reason in rejected[:5]:
            print(f"  - {reason}")
    return 0


# -- run ------------------------------------------------------------------
ALL_STAGES = ("download", "profile", "interpret", "validate", "report")


def _run_pipeline(args) -> int:
    from .stages import (
        RunContext, StageError, stage_download, stage_interpret, stage_profile,
        stage_validate, write_manifest,
    )

    samples = _resolve_samples(args)
    stages = args.stages or list(ALL_STAGES)
    unknown = [s for s in stages if s not in ALL_STAGES]
    if unknown:
        raise SystemExit(f"error: unknown stage(s) {unknown}; choose from "
                         f"{', '.join(ALL_STAGES)}")

    if args.write_script:
        return _write_script(args, samples, stages)

    context = RunContext.create(samples, args.data_dir, args.results_dir,
                               args.log_dir, args.threads)
    context.image_digests = {
        name: os.environ.get("CONTAINER_IMAGE_DIGEST", "")
        for name in ("mycobench",) if os.environ.get("CONTAINER_IMAGE_DIGEST")
    }
    print(f"[run] cohort {samples} ({len(context.rows)} row(s): "
          f"{len(context.benchmark_rows)} benchmark, "
          f"{len(context.control_rows)} control)")

    outcomes = []
    validation = None
    try:
        if "download" in stages:
            from .estimate import confirm, estimate
            if not confirm(estimate(samples, args.data_dir).render(samples),
                           assume_yes=args.yes, interactive=sys.stdin.isatty()):
                print("Download declined. Nothing was fetched.")
                return 1
            outcomes.append(stage_download(context))
        if "profile" in stages:
            outcomes.append(stage_profile(context))
        if "interpret" in stages:
            outcomes.append(stage_interpret(
                context, Path(args.catalogue) if args.catalogue else None))
        if "validate" in stages:
            validation = stage_validate(
                context, Path(args.phenotypes) if args.phenotypes else None)
    except StageError as exc:
        raise SystemExit(f"error: {exc}")

    for outcome in outcomes:
        print(outcome.summary())
    if validation:
        print(validation.summary())
        for note in validation.notes:
            print(f"NOTE: {note}")

    manifest = write_manifest(context, outcomes, validation)
    print(f"Wrote {manifest}")

    if "report" in stages:
        from .report import build_report
        path = build_report(context.results_dir, samples,
                            Path(args.results_dir) / "mycobench_report.html")
        print(f"Wrote {path}")
    return 0


def _write_script(args, samples: Path, stages: list[str]) -> int:
    """Emit the exact invocation as an editable script instead of running it."""
    path = Path(args.write_script)
    command = [sys.executable, "-m", "mycobench", "run",
               "--samples", str(samples), "--threads", str(args.threads),
               "--data-dir", str(args.data_dir),
               "--results-dir", str(args.results_dir),
               "--log-dir", str(args.log_dir), "--stages", *stages]
    if args.catalogue:
        command += ["--catalogue", str(args.catalogue)]
    if args.phenotypes:
        command += ["--phenotypes", str(args.phenotypes)]
    lines = [
        "#!/usr/bin/env bash",
        f"# Generated by 'mycobench run --write-script {path.name}'.",
        "# A plain, editable script, not a hidden CLI internal: change any",
        "# value, drop a stage, or add your own commands, then run it.",
        "#",
        "# The download and profile stages need sra-tools and tb-profiler.",
        "# Run those inside their containers:",
        "#   benchmark/docker_run.sh sra        ...",
        "#   benchmark/docker_run.sh tbprofiler ...",
        "set -euo pipefail",
        "",
        " ".join(shlex.quote(part) for part in command),
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    try:
        path.chmod(path.stat().st_mode | 0o111)
    except OSError:
        pass
    print(f"Wrote {path}\n")
    print("\n".join(lines))
    print(f"\nReview or edit it, then run: bash {path}")
    return 0


# -- report ---------------------------------------------------------------
def _run_report(args) -> int:
    from .report import build_report

    samples = Path(args.samples) if args.samples else None
    path = build_report(Path(args.results_dir), samples, Path(args.out))
    print(f"Wrote {path}")
    return 0


# -- analyse --------------------------------------------------------------

def _run_baseline(args) -> int:
    """Measure what the incumbent catalogue achieves on this cohort."""
    import json
    from datetime import datetime, timezone

    from .analysis import baseline as baseline_module
    from .analysis.genotypes import CoordinateIndex, load_genotypes
    from .analysis.pipeline import _read_reuse_rows
    from .analysis.strata import DeterminantIndex
    from .cohort import write_rows

    catalogue = Path(args.catalogue)
    print(f"[baseline] catalogue  {catalogue}")
    print(f"[baseline] VCF cache  {args.cache_dir}")

    rows = _read_reuse_rows(Path(args.reuse_table))
    coordinates = CoordinateIndex.from_catalogue(catalogue)
    determinants = DeterminantIndex.from_catalogue(catalogue)
    load = load_genotypes(
        rows, coordinates, Path(args.cache_dir), limit=args.limit,
        sample_seed=args.sample_seed, jobs=args.jobs,
        **({"cached_only": True} if args.cached_only else {}))
    if not load.isolates:
        raise SystemExit(
            "error: no isolate was genotyped; nothing can be measured")

    isolates = list(load.isolates.values())
    results = baseline_module.measure(
        isolates, rows, determinants,
        drugs=tuple(args.drug) if args.drug else None,
        assessed_fraction=args.assessed_fraction)

    print()
    print(load.describe())
    print()
    for drug in sorted(results):
        measured = results[drug]
        print(f"  {measured.describe()}")
        for note in measured.notes:
            print(f"      note: {note}")

    out_dir = Path(args.out_dir)
    written = [write_rows(out_dir / "catalogue_baseline.tsv",
                          baseline_module.baseline_rows(results),
                          baseline_module.BASELINE_COLUMNS)]
    payload = {
        "provenance": {
            "mycobench_version": __version__,
            "generated_at": datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"),
            "catalogue": str(catalogue),
            "reuse_table": str(args.reuse_table),
            "limit": args.limit,
            "sample_seed": args.sample_seed,
            "assessed_fraction": args.assessed_fraction,
            "n_isolates_genotyped": load.n_loaded,
        },
        "baselines": {
            drug: measured.as_registry_baseline(args.source)
            for drug, measured in sorted(results.items())
        },
        "caveat": baseline_module.CAVEAT,
    }
    manifest = out_dir / "catalogue_baseline.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    written.append(manifest)
    for path in written:
        print(f"Wrote {path}")

    print()
    print("These blocks are what a model must beat to be approved; paste one")
    print("into a registered model's performance row as its 'baseline'.")
    print(baseline_module.CAVEAT)
    return 0


def _run_analyse(args) -> int:
    from .analysis.pipeline import AnalysisError, AnalysisInputs, run, write_outputs

    inputs = AnalysisInputs(
        catalogue=Path(args.catalogue),
        phenotypes=Path(args.phenotypes),
        reuse_table=Path(args.reuse_table),
        cache_dir=Path(args.cache_dir),
        limit=args.limit,
        min_carriers=args.min_carriers,
        drugs=tuple(args.drug) if args.drug else None,
        cached_only=args.cached_only, metadata=Path(args.metadata) if args.metadata else None,
        partition=args.partition, independent_clusters=args.independent_clusters,
        jobs=args.jobs, sample_seed=args.sample_seed,
    )
    print(f"[analyse] catalogue  {inputs.catalogue}")
    print(f"[analyse] phenotypes {inputs.phenotypes}")
    print(f"[analyse] VCF cache  {inputs.cache_dir}")
    if args.limit:
        print(f"[analyse] limited to the first {args.limit} isolate(s)")
    else:
        print("[analyse] no --limit: every isolate in the reuse table will be "
              "genotyped. Measured, that is ~20 MB per isolate and ~248 GB "
              "for the full compendium, so use --limit unless the storage is "
              "provisioned")

    try:
        result = run(inputs)
    except AnalysisError as exc:
        raise SystemExit(f"error: {exc}")

    print()
    print(result.summary())
    if result.skipped:
        print(f"\nskipped {len(result.skipped)} variant-drug pair(s); first:")
        for line in result.skipped[:5]:
            print(f"  - {line}")

    for path in write_outputs(result, args.out_dir):
        print(f"Wrote {path}")

    print("\nNothing above is a resistance call. Effect sizes are conditional")
    print("on resistance background and restricted to catalogued coordinates;")
    print("no variant is promoted to a predictive call by this analysis.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mycobench",
        description="Real-data validation harness for Myconductor. "
                    "Concordance always; accuracy only with paired phenotypes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  mycobench validate-cohort --cohort nigeria-v1 --online "
               "--country Nigeria --write-lock cohorts/nigeria-v1.lock.json\n"
               "  mycobench estimate-download --cohort nigeria-v1\n"
               "  mycobench fetch-catalogue --out-dir data/who\n"
               "  mycobench run --cohort nigeria-v1 --write-script run_v1.sh\n"
               "  mycobench report --results-dir results\n")
    parser.add_argument("--version", action="version",
                        version=f"mycobench {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    from .workflow_evaluation import run as evaluate_workflow
    evaluation = sub.add_parser("evaluate-workflow", help="Evaluate paired laboratory review against independent phenotypes.")
    evaluation.add_argument("--input", required=True, help="Paired isolate-drug TSV or CSV.")
    evaluation.add_argument("--out", required=True)
    evaluation.add_argument("--html")
    evaluation.set_defaults(func=evaluate_workflow)

    d = sub.add_parser("discover-cohort",
                       help="Find candidate runs at NCBI, confirming origin.")
    d.add_argument("--organism", action="append", required=True)
    d.add_argument("--out-dir", required=True)
    d.add_argument("--country", help="Require BioSample geo_loc_name to "
                                     "confirm this country.")
    d.add_argument("--max-runs", type=int, default=600)
    d.add_argument("--prefix", default="", help="sample_id prefix, e.g. ng")
    d.add_argument("--non-mtbc", action="store_true",
                   help="Build a non-MTBC species-control panel instead.")
    d.add_argument("--min-bases", type=int, default=50_000_000)
    d.add_argument("--email")
    d.add_argument("--api-key")
    d.set_defaults(func=_run_discover)

    v = sub.add_parser("validate-cohort", help="Validate a cohort, or verify its lock.")
    group = v.add_mutually_exclusive_group()
    group.add_argument("--samples")
    group.add_argument("--cohort", help="Shorthand for cohorts/NAME.tsv")
    v.add_argument("--online", action="store_true")
    v.add_argument("--country")
    v.add_argument("--write-lock")
    v.add_argument("--verify-lock")
    v.add_argument("--email")
    v.add_argument("--api-key")
    v.set_defaults(func=_run_validate)

    r = sub.add_parser("review-candidates",
                       help="Measure study clustering and cap a candidate table.")
    r.add_argument("--candidates", required=True)
    r.add_argument("--out-dir", required=True)
    r.add_argument("--max-per-bioproject", type=int)
    r.add_argument("--max-per-organism", type=int)
    r.add_argument("--max-total", type=int)
    r.set_defaults(func=_run_review)

    e = sub.add_parser("estimate-download",
                       help="Say what a run will download, before it does.")
    eg = e.add_mutually_exclusive_group()
    eg.add_argument("--samples")
    eg.add_argument("--cohort")
    e.add_argument("--data-dir", default="data")
    e.set_defaults(func=_run_estimate)

    c = sub.add_parser("fetch-catalogue",
                       help="Fetch and ingest the published WHO catalogue.")
    c.add_argument("--out-dir", default="data/who")
    c.add_argument("--catalogue-out",
                   default="data/who/mtb_amr_catalogue.ingested.json")
    c.add_argument("--write-profile", action="store_true",
                   help="Also emit drug->loci tiers from WHO's own tiering.")
    c.add_argument("--profile-out", default="data/who/drug_loci.ingested.json")
    c.add_argument("--workbook", action="store_true",
                   help="Also download the 31 MB .xlsx (not needed to ingest).")
    c.set_defaults(func=_run_fetch_catalogue)

    p = sub.add_parser("fetch-phenotypes",
                       help="Fetch CRyPTIC paired MIC/DST phenotypes.")
    p.add_argument("--out-dir", default="data/cryptic")
    p.add_argument("--cohort-out",
                   help="Also build a phenotyped cohort TSV (resolves runs at NCBI).")
    p.add_argument("--limit", type=int,
                   help="Cap isolates resolved at NCBI; useful for a first run.")
    p.add_argument("--email")
    p.add_argument("--api-key")
    p.set_defaults(func=_run_fetch_phenotypes)

    run = sub.add_parser("run", help="Run the pipeline stages.")
    rg = run.add_mutually_exclusive_group()
    rg.add_argument("--samples")
    rg.add_argument("--cohort")
    run.add_argument("--stages", nargs="+", metavar="STAGE",
                     help=f"Subset of: {', '.join(ALL_STAGES)}")
    run.add_argument("--threads", type=int, default=4)
    run.add_argument("--data-dir", default=os.environ.get("DATA_DIR", "data"))
    run.add_argument("--results-dir",
                     default=os.environ.get("RESULTS_DIR", "results"))
    run.add_argument("--log-dir", default=os.environ.get("LOG_DIR", "logs"))
    run.add_argument("--catalogue", help="Ingested catalogue JSON to use.")
    run.add_argument("--phenotypes",
                     help="Phenotype TSV; without it NO accuracy is computed.")
    run.add_argument("--write-script", help="Emit the invocation instead of running.")
    run.add_argument("--yes", action="store_true",
                     help="Skip the download confirmation prompt.")
    run.set_defaults(func=_run_pipeline)

    rep = sub.add_parser("report", help="Build the HTML report from results.")
    rep.add_argument("--results-dir", default="results")
    rep.add_argument("--samples", help="Cohort TSV, for composition detail.")
    rep.add_argument("--out", default="results/mycobench_report.html")
    rep.set_defaults(func=_run_report)

    an = sub.add_parser(
        "analyse",
        help="Background-conditional effect sizes for uncertain catalogue "
             "variants, from CRyPTIC genotypes and MICs.")
    an.add_argument("--catalogue",
                    default="data/who/mtb_amr_catalogue.ingested.json",
                    help="INGESTED catalogue JSON; the bundled illustrative "
                         "file is refused.")
    an.add_argument("--phenotypes",
                    default="data/cryptic/cryptic_phenotypes.tsv")
    an.add_argument("--reuse-table",
                    default="data/cryptic/CRyPTIC_reuse_table_20240917.csv")
    an.add_argument("--cache-dir", default="data/cryptic/vcf",
                    help="Where per-isolate VCFs are cached (~25 KB each).")
    an.add_argument("--out-dir", default="results/analysis")
    an.add_argument("--limit", type=int,
                    help="Genotype only the first N isolates; useful for a "
                         "first run.")
    an.add_argument("--min-carriers", type=int, default=5,
                    help="Skip a variant carried by fewer isolates (default 5).")
    an.add_argument("--drug", action="append",
                    help="Restrict to these drugs; repeat as needed.")
    an.add_argument("--cached-only", action="store_true", help="Never download missing VCFs.")
    an.add_argument("--sample-seed", type=int, default=None, metavar="S",
                    help="Draw --limit isolates at random with this seed "
                         "instead of taking the first N, which come from one "
                         "or two sites and leave the confounding checks inert.")
    an.add_argument("--jobs", type=int, default=1, metavar="N",
                    help="Download this many VCFs concurrently before parsing "
                         "(parsing stays serial). Each is ~20 MB.")
    an.add_argument("--metadata", help="TSV with isolate_id, lineage, site_id, patient_id, cluster_id, partition.")
    an.add_argument("--partition", choices=("development", "evaluation"))
    an.add_argument("--independent-clusters", action="store_true",
                    help="Keep one deterministic representative per patient and reviewed genetic cluster.")
    an.set_defaults(func=_run_analyse)

    from .analysis.baseline import ASSESSED_FRACTION_FOR_SUSCEPTIBLE

    bl = sub.add_parser(
        "baseline",
        help="Measure what the WHO catalogue itself achieves on a cohort.")
    bl.add_argument("--catalogue", required=True,
                    help="INGESTED catalogue JSON.")
    bl.add_argument("--reuse-table", required=True)
    bl.add_argument("--cache-dir", required=True,
                    help="Where per-isolate VCFs are cached.")
    bl.add_argument("--out-dir", default="results/baseline")
    bl.add_argument("--limit", type=int, default=None)
    bl.add_argument("--sample-seed", type=int, default=None, metavar="S",
                    help="Draw --limit isolates at random with this seed.")
    bl.add_argument("--jobs", type=int, default=1, metavar="N")
    bl.add_argument("--cached-only", action="store_true")
    bl.add_argument("--drug", action="append")
    bl.add_argument("--assessed-fraction", type=float,
                    default=ASSESSED_FRACTION_FOR_SUSCEPTIBLE,
                    metavar="F",
                    help="Share of a drug's graded coordinates that must have "
                         "been examined before 'no resistant variant found' "
                         "may mean susceptible.")
    bl.add_argument("--source", default="WHO catalogue",
                    help="Name recorded as the baseline's source.")
    bl.set_defaults(func=_run_baseline)

    sub.add_parser("version", help="Print version.").set_defaults(
        func=lambda _a: (print(f"mycobench {__version__}"), 0)[1])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
