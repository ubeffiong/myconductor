"""Pipeline stages: download, profile, interpret, validate.

    download  SRA reads per isolate                      (sra-tools)
    profile   TB-Profiler over each isolate              (tb-profiler)
    interpret Myconductor over TB-Profiler's output      (in-process)
    validate  concordance / accuracy / species controls  (in-process)

Two contracts carried over from the stage design that works
-----------------------------------------------------------
**A failing isolate is recorded and skipped, never fatal.** One transient SRA
error must not discard hours of successful downloads in a 213-isolate cohort.
Every stage writes a per-isolate status table and aborts only when *nothing*
succeeded, because then there is nothing to benchmark.

**External commands are injected.** ``ToolRunner`` is a seam, so the
orchestration is unit-testable without sra-tools or TB-Profiler installed. The
default runner shells out; the recording runner used in tests does not.

Where coverage comes from, and why it matters here
--------------------------------------------------
Myconductor will not report a drug susceptible without independent evidence its
loci were callable. TB-Profiler's *collate* output carries only a genome-wide
median depth, which says nothing about the breadth of any single locus — so
this stage takes the mask from TB-Profiler's per-sample JSON, whose QC block
reports per-gene coverage. Without that block no drug becomes susceptible, and
the run reports that as the reason rather than quietly producing empty results.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Optional, Protocol

from . import __version__
from .cohort import read_rows, write_rows
from .metrics import (
    AccuracyResult,
    ConcordanceResult,
    SpeciesControlResult,
    score_accuracy,
    score_concordance,
    score_species_control,
)
from .monitoring import MonitoringState, generalizability, ingest_cohort, lineage_summary
from .thresholds import (
    MIN_LINEAGES_FOR_GENERALIZABILITY,
    TARGETS,
    describe as describe_thresholds,
    registration_hash,
)

STATUS_COLUMNS = ("sample_id", "status", "reason")
ENGINE_CALL_COLUMNS = ("sample_id", "drug", "engine_call", "variants",
                       "lineage", "drtype", "median_depth", "parser")
INTERPRET_CALL_COLUMNS = ("sample_id", "drug", "call", "tier", "permits_use",
                          "reason", "coverage_source")
CONCORDANCE_COLUMNS = ("drug", "compared", "agree", "disagree",
                       "agreement_rate", "only_myconductor", "only_engine",
                       "neither")
ACCURACY_COLUMNS = ("drug", "evaluable", "called", "call_rate", "abstained",
                    "tp", "fp", "tn", "fn", "sensitivity", "specificity",
                    "ppv", "npv", "vme_rate", "me_rate", "verdict", "reasons")
CONTROL_COLUMNS = ("sample_id", "organism", "refused", "passed", "detail")
LINEAGE_ACCURACY_COLUMNS = ("drug", "lineage", "evaluable", "called",
                           "call_rate", "abstained", "tp", "fp", "tn", "fn",
                           "sensitivity", "specificity", "vme_rate", "me_rate",
                           "verdict", "n_lineages_for_drug", "generalizable")


class StageError(RuntimeError):
    """A stage could not run at all."""


class ToolRunner(Protocol):
    def run(self, command: list[str], log_path: Path,
            timeout: Optional[int] = None) -> int:
        ...

    def available(self, executable: str) -> bool:
        ...


@dataclass
class SubprocessRunner:
    """Shells out, tees output to a per-sample log."""

    def run(self, command: list[str], log_path: Path,
            timeout: Optional[int] = None) -> int:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "ab") as log:
            log.write(f"\n$ {' '.join(command)}\n".encode())
            log.flush()
            try:
                completed = subprocess.run(command, stdout=log, stderr=log,
                                           timeout=timeout, check=False)
                return completed.returncode
            except subprocess.TimeoutExpired:
                log.write(b"\n[timed out]\n")
                return 124
            except FileNotFoundError as exc:
                log.write(f"\n[not found: {exc}]\n".encode())
                return 127

    def available(self, executable: str) -> bool:
        return shutil.which(executable) is not None


@dataclass
class RecordingRunner:
    """Records commands without executing them. For tests and --write-script."""

    commands: list[list[str]] = field(default_factory=list)
    exit_code: int = 0
    present: tuple[str, ...] = ("prefetch", "fasterq-dump", "tb-profiler")

    def run(self, command: list[str], log_path: Path,
            timeout: Optional[int] = None) -> int:
        self.commands.append(list(command))
        return self.exit_code

    def available(self, executable: str) -> bool:
        return executable in self.present


@dataclass
class RunContext:
    cohort_path: Path
    data_dir: Path
    results_dir: Path
    log_dir: Path
    threads: int = 4
    runner: ToolRunner = field(default_factory=SubprocessRunner)
    rows: list[dict] = field(default_factory=list)
    image_digests: dict[str, str] = field(default_factory=dict)

    @classmethod
    def create(cls, cohort_path: str | Path, data_dir: str | Path,
               results_dir: str | Path, log_dir: str | Path,
               threads: int = 4,
               runner: Optional[ToolRunner] = None) -> "RunContext":
        rows, _ = read_rows(cohort_path)
        context = cls(Path(cohort_path), Path(data_dir), Path(results_dir),
                      Path(log_dir), threads,
                      runner or SubprocessRunner(), rows)
        for directory in (context.data_dir, context.results_dir, context.log_dir):
            directory.mkdir(parents=True, exist_ok=True)
        return context

    def sample_dir(self, sample_id: str) -> Path:
        return self.data_dir / sample_id

    def result_dir(self, sample_id: str) -> Path:
        return self.results_dir / sample_id

    def log(self, sample_id: str, suffix: str) -> Path:
        return self.log_dir / f"{sample_id}.{suffix}.log"

    @property
    def benchmark_rows(self) -> list[dict]:
        """Rows scored as benchmark subjects, excluding declared controls."""
        return [r for r in self.rows
                if not (r.get("expected_outcome") or "").strip()]

    @property
    def control_rows(self) -> list[dict]:
        return [r for r in self.rows
                if (r.get("expected_outcome") or "").strip()]


@dataclass
class StageOutcome:
    stage: str
    ok: list[str] = field(default_factory=list)
    failed: list[tuple[str, str]] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)

    @property
    def n_ok(self) -> int:
        return len(self.ok)

    def status_rows(self) -> list[dict]:
        rows = [{"sample_id": s, "status": "ok", "reason": ""} for s in self.ok]
        rows += [{"sample_id": s, "status": "failed", "reason": r}
                 for s, r in self.failed]
        rows += [{"sample_id": s, "status": "skipped", "reason": r}
                 for s, r in self.skipped]
        return rows

    def write(self, results_dir: Path) -> Path:
        return write_rows(results_dir / f"{self.stage}_status.tsv",
                          self.status_rows(), STATUS_COLUMNS)

    def summary(self) -> str:
        return (f"stage {self.stage}: {self.n_ok} ok, "
                f"{len(self.failed)} failed, {len(self.skipped)} skipped")


# -- stage 1: download ----------------------------------------------------
def reads_for(context: RunContext, row: dict) -> tuple[Path, Path]:
    directory = context.sample_dir(row["sample_id"])
    run = row["sra_run"]
    return directory / f"{run}_1.fastq.gz", directory / f"{run}_2.fastq.gz"


def stage_download(context: RunContext) -> StageOutcome:
    outcome = StageOutcome("download")
    for tool in ("prefetch", "fasterq-dump"):
        if not context.runner.available(tool):
            raise StageError(
                f"{tool} is not on PATH. Run this stage inside the sra-tools "
                f"container (benchmark/docker_run.sh), or install sra-tools.")

    for row in context.rows:
        sample_id, run = row["sample_id"], row["sra_run"]
        directory = context.sample_dir(sample_id)
        directory.mkdir(parents=True, exist_ok=True)
        r1, r2 = reads_for(context, row)
        if r1.is_file() and r2.is_file():
            outcome.skipped.append((sample_id, "reads already present"))
            continue

        log = context.log(sample_id, "download")
        code = context.runner.run(
            ["prefetch", "--max-size", "100G", "-O", str(directory), run], log)
        if code != 0:
            outcome.failed.append((sample_id, f"prefetch exited {code}; see {log}"))
            continue
        code = context.runner.run(
            ["fasterq-dump", "--split-files", "--threads", str(context.threads),
             "-O", str(directory), run], log)
        if code != 0:
            outcome.failed.append(
                (sample_id, f"fasterq-dump exited {code}; see {log}"))
            continue

        # fasterq-dump writes uncompressed FASTQ; compress so a 213-isolate
        # cohort does not need several hundred gigabytes of plain text.
        for suffix in ("_1", "_2"):
            plain = directory / f"{run}{suffix}.fastq"
            if plain.is_file():
                context.runner.run(["gzip", "-f", str(plain)], log)
        if not (r1.is_file() and r2.is_file()):
            outcome.failed.append(
                (sample_id, "paired FASTQ not produced (single-end run?)"))
            continue
        outcome.ok.append(sample_id)

    outcome.write(context.results_dir)
    if not outcome.ok and not outcome.skipped:
        raise StageError(
            "no isolate downloaded successfully; nothing to benchmark. See "
            f"{context.results_dir / 'download_status.tsv'}")
    return outcome


# -- stage 2: profile -----------------------------------------------------
def tbprofiler_json(context: RunContext, sample_id: str) -> Path:
    return (context.result_dir(sample_id) / "tbprofiler" / "results"
            / f"{sample_id}.results.json")


def stage_profile(context: RunContext) -> StageOutcome:
    outcome = StageOutcome("profile")
    if not context.runner.available("tb-profiler"):
        raise StageError(
            "tb-profiler is not on PATH. Run this stage inside the "
            "TB-Profiler container (benchmark/docker_run.sh tbprofiler).")

    for row in context.rows:
        sample_id = row["sample_id"]
        r1, r2 = reads_for(context, row)
        if not (r1.is_file() and r2.is_file()):
            outcome.skipped.append((sample_id, "reads absent; download failed"))
            continue
        destination = context.result_dir(sample_id) / "tbprofiler"
        if tbprofiler_json(context, sample_id).is_file():
            outcome.skipped.append((sample_id, "already profiled"))
            continue
        destination.mkdir(parents=True, exist_ok=True)
        log = context.log(sample_id, "tbprofiler")
        code = context.runner.run(
            ["tb-profiler", "profile", "-1", str(r1), "-2", str(r2),
             "-p", sample_id, "-d", str(destination),
             "-t", str(context.threads), "--txt", "--csv"], log)
        if code != 0:
            outcome.failed.append(
                (sample_id, f"tb-profiler exited {code}; see {log}"))
            continue
        outcome.ok.append(sample_id)

    if outcome.ok:
        # One collate table across the cohort: the validated parser's input.
        context.runner.run(
            ["tb-profiler", "collate", "-d",
             str(context.results_dir), "-p",
             str(context.results_dir / "tbprofiler_collate")],
            context.log_dir / "collate.log")

    outcome.write(context.results_dir)
    if not outcome.ok and not outcome.skipped:
        raise StageError("no isolate was profiled; nothing to interpret")
    return outcome


# -- stage 3: interpret ---------------------------------------------------
def stage_interpret(context: RunContext,
                    catalogue_path: Optional[Path] = None) -> StageOutcome:
    """Run Myconductor over each isolate's TB-Profiler output, in process."""
    from myconductor import Myconductor
    from myconductor.adapters.base import AdapterSchemaError
    from myconductor.adapters.tbprofiler import TBProfilerAdapter
    from myconductor.modules.catalogue import CatalogueModule
    from myconductor.reporting.render import render_text

    from .engines import VARIANT_TSV_COLUMNS, parse_results_json, to_variant_tsv

    outcome = StageOutcome("interpret")
    engine_rows: list[dict] = []
    interpret_rows: list[dict] = []

    catalogue = (CatalogueModule(catalogue_path) if catalogue_path
                 else CatalogueModule())
    conductor = Myconductor(platform="illumina")
    conductor.router.catalogue = catalogue

    for row in context.rows:
        sample_id = row["sample_id"]
        results_json = tbprofiler_json(context, sample_id)
        if not results_json.is_file():
            outcome.skipped.append((sample_id, "no TB-Profiler result"))
            continue
        try:
            engine_sample = parse_results_json(results_json)
            engine_report = TBProfilerAdapter().parse(results_json)
        except (AdapterSchemaError, ValueError) as exc:
            outcome.failed.append((sample_id, f"could not parse output: {exc}"))
            continue

        for drug, call in sorted(engine_sample.resistance_calls().items()):
            finding = engine_sample.findings.get(drug)
            engine_rows.append({
                "sample_id": sample_id, "drug": drug, "engine_call": call,
                "variants": ";".join(v.label() for v in finding.variants)
                            if finding else "",
                "lineage": engine_sample.sub_lineage or engine_sample.main_lineage,
                "drtype": engine_sample.drtype,
                "median_depth": engine_sample.median_depth or "",
                "parser": engine_sample.parser,
            })

        # Variants go through Myconductor's own input adapter, so they are
        # normalised the same way any other input is.
        variant_tsv = context.result_dir(sample_id) / "engine_variants.tsv"
        write_rows(variant_tsv, to_variant_tsv(engine_sample),
                   VARIANT_TSV_COLUMNS)
        try:
            report = conductor.analyze(
                variant_tsv, mask=engine_report.mask,
                engine_reports=[engine_report])
        except (ValueError, OSError) as exc:
            outcome.failed.append((sample_id, f"interpretation failed: {exc}"))
            continue

        for result in report.drug_results:
            interpret_rows.append({
                "sample_id": sample_id, "drug": result.drug,
                "call": result.call.value, "tier": result.tier.value,
                "permits_use": "yes" if result.permits_use else "no",
                "reason": result.reason or "",
                "coverage_source": report.provenance.coverage_source,
            })
        (context.result_dir(sample_id) / "myconductor_report.txt").write_text(
            render_text(report), encoding="utf-8")
        outcome.ok.append(sample_id)

    write_rows(context.results_dir / "engine_calls.tsv", engine_rows,
               ENGINE_CALL_COLUMNS)
    write_rows(context.results_dir / "myconductor_calls.tsv", interpret_rows,
               INTERPRET_CALL_COLUMNS)
    outcome.write(context.results_dir)
    if not outcome.ok:
        raise StageError("no isolate was interpreted; nothing to validate")
    return outcome


# -- stage 4: validate ----------------------------------------------------
def _load_calls(path: Path, call_field: str) -> dict[tuple[str, str], str]:
    if not path.is_file():
        return {}
    rows, _ = read_rows(path)
    return {(r["sample_id"], r["drug"]): r.get(call_field, "") for r in rows}


def _load_lineages(path: Path) -> dict[str, str]:
    """Per-sample lineage, from ``engine_calls.tsv``'s ``lineage`` column
    (already resolved in ``stage_interpret`` as
    ``sub_lineage or main_lineage``). First non-empty value wins."""
    if not path.is_file():
        return {}
    rows, _ = read_rows(path)
    out: dict[str, str] = {}
    for r in rows:
        lineage = (r.get("lineage") or "").strip()
        if lineage and r["sample_id"] not in out:
            out[r["sample_id"]] = lineage
    return out


@dataclass
class ValidationOutcome:
    concordance: list[ConcordanceResult] = field(default_factory=list)
    accuracy: list[AccuracyResult] = field(default_factory=list)
    controls: list[SpeciesControlResult] = field(default_factory=list)
    lineage_accuracy: dict = field(default_factory=dict)   # drug -> [LineageAccuracy]
    generalizability: dict = field(default_factory=dict)   # drug -> DrugGeneralizability
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = []
        if self.concordance:
            lines.append(f"concordance computed for {len(self.concordance)} drug(s)")
        if self.accuracy:
            verdicts: dict[str, int] = {}
            for result in self.accuracy:
                verdict, _ = result.verdict()
                verdicts[verdict] = verdicts.get(verdict, 0) + 1

            lines.append("accuracy verdicts: " + ", ".join(
                f"{v}={n}" for v, n in sorted(verdicts.items())))
        else:
            lines.append("accuracy: NOT COMPUTED (no paired phenotypes supplied)")
        if self.controls:
            passed = sum(1 for c in self.controls if c.passed)
            lines.append(f"species controls: {passed}/{len(self.controls)} refused "
                         f"as required")
        if self.generalizability:
            not_yet = sorted(d for d, g in self.generalizability.items()
                             if not g.generalizable)
            lines.append(
                f"lineage generalizability: {len(self.generalizability) - len(not_yet)}/"
                f"{len(self.generalizability)} drug(s) measured across "
                f"≥{MIN_LINEAGES_FOR_GENERALIZABILITY} lineages"
                + (f" (not yet: {', '.join(not_yet)})" if not_yet else ""))
        return "\n".join(lines)


def stage_validate(context: RunContext,
                   phenotypes_path: Optional[Path] = None,
                   monitoring_state_path: Optional[Path] = None) -> ValidationOutcome:
    outcome = ValidationOutcome()
    mine = _load_calls(context.results_dir / "myconductor_calls.tsv", "call")
    theirs = _load_calls(context.results_dir / "engine_calls.tsv", "engine_call")
    lineages = _load_lineages(context.results_dir / "engine_calls.tsv")

    benchmark_ids = {r["sample_id"] for r in context.benchmark_rows}
    drugs = sorted({drug for (sample, drug) in mine if sample in benchmark_ids})

    # -- concordance (always available) --------------------------------
    for drug in drugs:
        triples = [(sample, mine.get((sample, drug), ""),
                    theirs.get((sample, drug), ""))
                   for sample in sorted(benchmark_ids)
                   if (sample, drug) in mine or (sample, drug) in theirs]
        if triples:
            outcome.concordance.append(score_concordance(triples, drug))
    write_rows(context.results_dir / "concordance.tsv",
               [{"drug": c.drug, "compared": c.n_compared, "agree": c.agree,
                 "disagree": c.disagree,
                 "agreement_rate": (f"{c.agreement_rate:.4f}"
                                    if c.agreement_rate is not None else ""),
                 "only_myconductor": c.only_a_called,
                 "only_engine": c.only_b_called,
                 "neither": c.neither_called}
                for c in outcome.concordance], CONCORDANCE_COLUMNS)

    # -- accuracy (only with paired phenotypes) ------------------------
    if phenotypes_path and Path(phenotypes_path).is_file():
        phenotype_rows, _ = read_rows(phenotypes_path)
        by_pair: dict[tuple[str, str], str] = {}
        evaluable: dict[str, int] = {}
        for row in phenotype_rows:
            if (row.get("evaluable") or "").lower() != "yes":
                continue
            key = (row["sample_id"], row["drug"])
            by_pair[key] = (row.get("phenotype") or "").upper()
            evaluable[row["drug"]] = evaluable.get(row["drug"], 0) + 1

        for drug in sorted(set(list(evaluable) + drugs)):
            pairs = [
                (mine.get((sample_id, phenotype_drug), ""), phenotype)
                for (sample_id, phenotype_drug), phenotype in by_pair.items()
                if phenotype_drug == drug and sample_id in benchmark_ids
            ]
            if not pairs:
                continue
            outcome.accuracy.append(
                score_accuracy(pairs, drug, n_evaluable=len(pairs)))
        write_rows(context.results_dir / "accuracy.tsv",
                   [_accuracy_row(a) for a in outcome.accuracy],
                   ACCURACY_COLUMNS)

        # -- lineage-stratified, continuously-monitored accuracy ----------
        if any(lineages.values()):
            state_path = monitoring_state_path or (
                context.results_dir / "monitoring_state.json")
            state = MonitoringState.load(state_path)
            rows = [
                (phenotype_drug, lineages.get(sample_id, "unknown"),
                 mine.get((sample_id, phenotype_drug), ""), phenotype)
                for (sample_id, phenotype_drug), phenotype in by_pair.items()
                if sample_id in benchmark_ids
            ]
            state, notes = ingest_cohort(state, str(context.cohort_path), rows)
            outcome.notes.extend(notes)
            state.save(state_path)

            outcome.lineage_accuracy = lineage_summary(state)
            outcome.generalizability = generalizability(state)
            lineage_rows = []
            for drug, entries in sorted(outcome.lineage_accuracy.items()):
                gen = outcome.generalizability.get(drug)
                for la in entries:
                    row = _accuracy_row(la.accuracy)
                    row["lineage"] = la.lineage
                    row["n_lineages_for_drug"] = gen.n_lineages if gen else 0
                    row["generalizable"] = ("yes" if gen and gen.generalizable
                                            else "no")
                    lineage_rows.append(row)
            write_rows(context.results_dir / "lineage_accuracy.tsv",
                       lineage_rows, LINEAGE_ACCURACY_COLUMNS)
        else:
            outcome.notes.append(
                "No isolate carried a resolved lineage, so lineage-stratified "
                "accuracy was not computed. TB-Profiler's per-sample output "
                "supplies main_lineage/sub_lineage; without it, generalizability "
                "across lineages cannot be assessed.")
    else:
        outcome.notes.append(
            "No paired phenotypes supplied, so NO accuracy figure was "
            "computed. Concordance with TB-Profiler measures agreement, not "
            "correctness. Build a phenotyped cohort with "
            "'mycobench fetch-phenotypes'.")

    # -- species controls ----------------------------------------------
    for row in context.control_rows:
        sample_id = row["sample_id"]
        calls = [(drug, call) for (sample, drug), call in mine.items()
                 if sample == sample_id]
        if not calls:
            outcome.notes.append(
                f"{sample_id}: no interpretation produced, so the species "
                f"control could not be evaluated")
            continue
        outcome.controls.append(
            score_species_control(sample_id, row.get("organism", ""), calls))
    if outcome.controls:
        write_rows(context.results_dir / "species_control.tsv",
                   [{"sample_id": c.sample_id, "organism": c.organism,
                     "refused": "yes" if c.refused else "no",
                     "passed": "yes" if c.passed else "no",
                     "detail": c.detail} for c in outcome.controls],
                   CONTROL_COLUMNS)
    return outcome


def _accuracy_row(a: AccuracyResult) -> dict:
    verdict, reasons = a.verdict()

    def fmt(value: Optional[float]) -> str:
        return f"{value:.4f}" if value is not None else ""

    return {"drug": a.drug, "evaluable": a.n_evaluable, "called": a.n_called,
            "call_rate": fmt(a.call_rate), "abstained": a.n_abstained,
            "tp": a.true_positive, "fp": a.false_positive,
            "tn": a.true_negative, "fn": a.false_negative,
            "sensitivity": fmt(a.sensitivity),
            "specificity": fmt(a.specificity), "ppv": fmt(a.ppv),
            "npv": fmt(a.npv), "vme_rate": fmt(a.vme_rate),
            "me_rate": fmt(a.me_rate), "verdict": verdict,
            "reasons": "; ".join(reasons)}


# -- run manifest ---------------------------------------------------------
def write_manifest(context: RunContext, outcomes: Iterable[StageOutcome],
                   validation: Optional[ValidationOutcome] = None,
                   catalogue_provenance: Optional[dict] = None,
                   phenotype_provenance: Optional[dict] = None) -> Path:
    """Everything needed to say what produced these numbers."""
    import myconductor

    payload = {
        "mycobench_version": __version__,
        "myconductor_version": myconductor.__version__,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "cohort": {
            "path": str(context.cohort_path),
            "n_rows": len(context.rows),
            "n_benchmark": len(context.benchmark_rows),
            "n_controls": len(context.control_rows),
        },
        "container_images": dict(context.image_digests),
        "catalogue": catalogue_provenance or {
            "source": "bundled illustrative subset",
            "warning": "not the WHO catalogue; run 'mycobench fetch-catalogue'",
        },
        "phenotypes": phenotype_provenance or {
            "source": None,
            "note": "no paired phenotypes; accuracy not computed",
        },
        "pre_registration": {
            "description": describe_thresholds(),
            "sha256": registration_hash(),
            "drugs": sorted(TARGETS),
        },
        "stages": [{"stage": o.stage, "ok": o.n_ok,
                    "failed": len(o.failed), "skipped": len(o.skipped)}
                   for o in outcomes],
        "notes": list(validation.notes) if validation else [],
    }
    path = context.results_dir / "run_manifest.json"
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path
