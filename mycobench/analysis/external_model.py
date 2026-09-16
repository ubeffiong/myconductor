"""External-model benchmarking bridge for bring-your-own TB-AMR predictions.

The benchmark harness already knows how to evaluate a predictor that may
abstain: ``selective.risk_coverage_curve`` and
``selective.beats_baseline_at_matched_coverage``. This module supplies the
missing bridge from a plain prediction file into that machinery and then back
out into the performance-row shape expected by Myconductor's model registry.

It intentionally does not import or wrap named research tools. A deployment can
run DeepAMR, xAI-MTBDR, HANN-TT, Treesist-TB, TB-DROP, or any future model by
exporting the same four columns. If a tool is not publicly runnable, no adapter
is invented for it here.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Iterable, Optional, Sequence

from ..metrics import score_accuracy
from .selective import (
    BaselineComparison,
    Prediction,
    RiskCoverageCurve,
    beats_baseline_at_matched_coverage,
    bootstrap_risk_interval,
    risk_coverage_curve,
)

PREDICTION_COLUMNS = ("isolate_id", "drug", "predicted", "confidence")
RESULT_COLUMNS = (
    "model_id", "model_version", "drug", "lineage", "cohort", "source",
    "n_predictions", "n_evaluable", "n_called", "n_dropped", "call_rate",
    "error_rate", "baseline_source", "baseline_coverage", "baseline_error_rate",
    "matched_coverage", "matched_threshold", "matched_error_rate",
    "risk_ci", "verdict", "sensitivity", "specificity", "vme_rate",
    "me_rate", "registry_ready", "notes",
)
TRUTH_VALUES = {"R": "R", "RESISTANT": "R", "S": "S", "SUSCEPTIBLE": "S"}
PREDICTED_VALUES = {
    "R": "R", "RESISTANT": "R", "S": "S", "SUSCEPTIBLE": "S",
    "": None, "ABSTAIN": None, "NA": None, "N/A": None, "NONE": None,
    "NO_CALL": None, "NO CALL": None, "NOT_ASSESSED": None,
    "NOT ASSESSED": None, "INDETERMINATE": None, "UNKNOWN": None,
}
ACCEPTED_QUALITY = {"", "HIGH", "MEDIUM", "PASS", "ACCEPTED"}


@dataclass(frozen=True)
class ExternalPrediction:
    isolate_id: str
    drug: str
    predicted: Optional[str]
    confidence: Optional[float] = None


@dataclass(frozen=True)
class JoinNote:
    isolate_id: str
    drug: str
    reason: str


@dataclass
class ExternalModelResult:
    model_id: str
    model_version: str
    drug: str
    lineage: str
    cohort: str
    source: str
    predictions: list[Prediction]
    curve: RiskCoverageCurve
    comparison: BaselineComparison
    accuracy: object
    dropped: list[JoinNote] = field(default_factory=list)
    risk_interval: Optional[tuple[float, float]] = None

    @property
    def n_predictions(self) -> int:
        return len(self.predictions) + len(self.dropped)

    @property
    def n_evaluable(self) -> int:
        return len(self.predictions)

    @property
    def n_called(self) -> int:
        return sum(1 for p in self.predictions if p.answerable)

    @property
    def call_rate(self) -> float:
        return self.n_called / self.n_evaluable if self.n_evaluable else 0.0

    @property
    def error_rate(self) -> float:
        if not self.n_called:
            return 0.0
        return sum(1 for p in self.predictions if p.answerable and not p.correct) / self.n_called

    @property
    def baseline_source(self) -> str:
        return self.comparison.baseline_source or "catalogue baseline"

    def registry_exclusion(self) -> Optional[str]:
        # The registry compares summary numbers and cannot see how many
        # predictions stand behind them: five correct answers would pass it.
        # Only a scope with a measurable matched-coverage risk and estimable
        # sensitivity and specificity is exported; the rest keep their reason.
        matched = self.comparison.matched
        if matched is None or matched.selective_risk is None:
            return ("no measurable risk at the baseline's coverage ("
                    + self.comparison.verdict + ")")
        missing = [name for name in ("sensitivity", "specificity")
                   if getattr(self.accuracy, name) is None]
        if missing:
            return (" and ".join(missing) + " not estimable: the evaluable set "
                    "lacks phenotypically resistant or susceptible isolates")
        return None

    def registry_performance(self, independent: bool = True) -> dict:
        # Rates are taken at the matched-coverage point the verdict was decided
        # on, so the registry's summary check and this benchmark cannot
        # disagree about the same model.
        matched = self.comparison.matched
        if matched is not None and matched.selective_risk is not None:
            coverage, error = matched.coverage, matched.selective_risk
        else:
            coverage, error = self.call_rate, self.error_rate
        return {
            "drug": self.drug,
            "lineage": self.lineage,
            "cohort": self.cohort,
            "source": self.source,
            "independent": independent,
            "sensitivity": _round(self.accuracy.sensitivity),
            "specificity": _round(self.accuracy.specificity),
            "call_rate": round(coverage, 4),
            "error_rate": round(error, 4),
            "baseline": {
                "source": self.baseline_source,
                "coverage": round(self.comparison.baseline_coverage, 4),
                "error_rate": round(self.comparison.baseline_risk, 4),
            },
        }

    def row(self) -> dict:
        matched = self.comparison.matched
        risk = matched.selective_risk if matched else None
        return {
            "model_id": self.model_id,
            "model_version": self.model_version,
            "drug": self.drug,
            "lineage": self.lineage,
            "cohort": self.cohort,
            "source": self.source,
            "n_predictions": self.n_predictions,
            "n_evaluable": self.n_evaluable,
            "n_called": self.n_called,
            "n_dropped": len(self.dropped),
            "call_rate": f"{self.call_rate:.4f}",
            "error_rate": f"{self.error_rate:.4f}" if self.n_called else "",
            "baseline_source": self.baseline_source,
            "baseline_coverage": f"{self.comparison.baseline_coverage:.4f}",
            "baseline_error_rate": f"{self.comparison.baseline_risk:.4f}",
            "matched_coverage": f"{matched.coverage:.4f}" if matched else "",
            "matched_threshold": f"{matched.threshold:.4f}" if matched else "",
            "matched_error_rate": f"{risk:.4f}" if risk is not None else "",
            "risk_ci": (_interval_text(self.risk_interval) if self.risk_interval else ""),
            "verdict": self.comparison.verdict,
            "sensitivity": _rate(self.accuracy.sensitivity),
            "specificity": _rate(self.accuracy.specificity),
            "vme_rate": _rate(self.accuracy.vme_rate),
            "me_rate": _rate(self.accuracy.me_rate),
            "registry_ready": ("yes" if self.comparison.passed
                               and self.registry_exclusion() is None else "no"),
            "notes": " | ".join(
                self.comparison.reasons
                + ([f"not exported to the registry: {self.registry_exclusion()}"]
                   if self.registry_exclusion() else [])
                + _drop_summary(self.dropped)),
        }


def _round(value: Optional[float]) -> Optional[float]:
    # None stays None: a missing sensitivity rounded to 0.0 reads as measured.
    return round(float(value), 4) if value is not None else None


def _rate(value: Optional[float]) -> str:
    return "" if value is None else f"{value:.4f}"


def _interval_text(value: tuple[float, float]) -> str:
    return f"{value[0]:.4f}-{value[1]:.4f}"


def _drop_summary(dropped: Sequence[JoinNote]) -> list[str]:
    if not dropped:
        return []
    counts: dict[str, int] = {}
    for note in dropped:
        counts[note.reason] = counts.get(note.reason, 0) + 1
    return [f"dropped {n} row(s): {reason}" for reason, n in sorted(counts.items())]


def _dialect(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def _prediction_value(raw: str) -> Optional[str]:
    key = (raw or "").strip().upper().replace("-", "_")
    if key not in PREDICTED_VALUES:
        raise ValueError(f"unknown predicted value {raw!r}; use resistant, susceptible, or abstain/empty")
    return PREDICTED_VALUES[key]


def _truth_value(raw: str) -> Optional[str]:
    return TRUTH_VALUES.get((raw or "").strip().upper())


def load_predictions(path: str | Path) -> list[ExternalPrediction]:
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter=_dialect(path),
        )
        fields = set(reader.fieldnames or [])
        missing = {"isolate_id", "drug", "predicted"} - fields
        if missing:
            raise ValueError(f"{path}: prediction file missing column(s): {', '.join(sorted(missing))}")
        out: list[ExternalPrediction] = []
        seen: set[tuple[str, str]] = set()
        for number, row in enumerate(reader, 2):
            isolate = (row.get("isolate_id") or "").strip()
            drug = (row.get("drug") or "").strip().lower()
            if not isolate or not drug:
                raise ValueError(f"{path}:{number}: isolate_id and drug are required")
            key = (isolate, drug)
            if key in seen:
                raise ValueError(f"{path}:{number}: duplicate prediction for {isolate}/{drug}")
            seen.add(key)
            confidence = None
            raw_conf = (row.get("confidence") or "").strip()
            if raw_conf:
                confidence = float(raw_conf)
                if not 0 <= confidence <= 1:
                    raise ValueError(f"{path}:{number}: confidence must be in [0,1]")
            out.append(ExternalPrediction(isolate, drug, _prediction_value(row.get("predicted", "")), confidence))
    return out


def load_truths(path: str | Path) -> dict[tuple[str, str], str]:
    """Load accepted binary phenotype truths from generic or CRyPTIC-style rows."""
    from ..phenotypes import DRUG_CODES

    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter=_dialect(path),
        )
        fields = set(reader.fieldnames or [])
        truths: dict[tuple[str, str], str] = {}
        if {"sample_id", "drug", "phenotype"} <= fields:
            for row in reader:
                quality = (row.get("quality") or "").strip().upper()
                if quality not in ACCEPTED_QUALITY:
                    continue
                truth = _truth_value(row.get("phenotype", ""))
                if truth:
                    truths[((row.get("sample_id") or "").strip(), (row.get("drug") or "").strip().lower())] = truth
            return truths
        if "ENA_RUN" in fields:
            for row in reader:
                isolate = f"cr_{(row.get('ENA_RUN') or '').strip()}"
                for code, drug in DRUG_CODES.items():
                    truth = _truth_value(row.get(f"{code}_BINARY_PHENOTYPE", ""))
                    quality = (row.get(f"{code}_PHENOTYPE_QUALITY") or "").strip().upper()
                    if truth and quality in ACCEPTED_QUALITY:
                        truths[(isolate, drug)] = truth
            return truths
    raise ValueError(f"{path}: expected either sample_id/drug/phenotype rows or a CRyPTIC reuse table")


def load_catalogue_baselines(path: str | Path) -> tuple[dict[str, dict], str]:
    path = Path(path)
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        baselines = payload.get("baselines") or payload
        out = {str(drug).lower(): block for drug, block in baselines.items() if isinstance(block, dict)}
        source = payload.get("provenance", {}).get("catalogue") or payload.get("source") or str(path)
        return out, str(source)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=_dialect(path))
        out = {}
        for row in reader:
            drug = (row.get("drug") or "").strip().lower()
            if not drug:
                continue
            coverage = _float(row.get("coverage"))
            error = _float(row.get("error_rate"))
            if coverage is not None and error is not None:
                out[drug] = {"source": str(path), "coverage": coverage, "error_rate": error}
        return out, str(path)


def _float(raw) -> Optional[float]:
    if raw is None or raw == "":
        return None
    return float(raw)


def to_selective_predictions(
    external: Iterable[ExternalPrediction],
    truths: dict[tuple[str, str], str],
    weight: Optional[float] = None,
) -> tuple[dict[str, list[Prediction]], list[JoinNote]]:
    grouped: dict[str, list[Prediction]] = {}
    dropped: list[JoinNote] = []
    for row in external:
        truth = truths.get((row.isolate_id, row.drug.lower()))
        if truth is None:
            dropped.append(JoinNote(row.isolate_id, row.drug, "no accepted binary phenotype"))
            continue
        confidence = row.confidence if row.confidence is not None else (weight if weight is not None else 1.0)
        grouped.setdefault(row.drug.lower(), []).append(
            Prediction(row.isolate_id, row.predicted, truth, confidence))
    return grouped, dropped


def measure_external_model(
    predictions: Iterable[ExternalPrediction],
    truths: dict[tuple[str, str], str],
    catalogue_baseline: dict[str, dict],
    drugs: Optional[Iterable[str]] = None,
    model_id: str = "external-model",
    model_version: str = "unspecified",
    cohort: str = "external-evaluation",
    source: str = "external prediction file",
    lineage: str = "unknown",
    baseline_source: Optional[str] = None,
) -> tuple[dict[str, ExternalModelResult], list[str]]:
    """Score ``predictions`` per drug against ``catalogue_baseline``.

    A drug with predictions but no supplied baseline is **skipped, not
    fatal** — one drug the catalogue baseline doesn't cover (e.g.
    pretomanid, unmeasurable from either source; see
    ``baseline.py::measurability``) must not discard every other drug's
    result. The skip is reported back as a note, mirroring this codebase's
    "record and continue" posture everywhere else (stage failures,
    unassessed drugs, dropped predictions).
    """
    wanted = {d.lower() for d in drugs} if drugs else None
    external = [p for p in predictions if wanted is None or p.drug.lower() in wanted]
    grouped, dropped = to_selective_predictions(external, truths)
    dropped_by_drug: dict[str, list[JoinNote]] = {}
    for note in dropped:
        dropped_by_drug.setdefault(note.drug.lower(), []).append(note)
    results: dict[str, ExternalModelResult] = {}
    skipped: list[str] = []
    for drug in sorted(set(grouped) | set(dropped_by_drug)):
        baseline = catalogue_baseline.get(drug)
        if not baseline:
            skipped.append(
                f"{drug}: no catalogue baseline supplied; cannot compare at "
                f"matched coverage")
            continue
        rows = grouped.get(drug, [])
        curve = risk_coverage_curve(rows)
        comparison = beats_baseline_at_matched_coverage(
            curve, float(baseline["coverage"]), float(baseline["error_rate"]))
        comparison = replace(
            comparison,
            baseline_source=baseline.get("source") or baseline_source
            or "catalogue baseline")
        accuracy = score_accuracy(
            [("resistant" if p.predicted == "R" else "susceptible" if p.predicted == "S" else "not_assessed", p.truth) for p in rows], drug)
        interval = None
        if comparison.matched is not None:
            interval = bootstrap_risk_interval(rows, comparison.matched.threshold)
        results[drug] = ExternalModelResult(
            model_id, model_version, drug, lineage, cohort, source, rows,
            curve, comparison, accuracy, dropped_by_drug.get(drug, []), interval)
    return results, skipped


def result_rows(results: dict[str, ExternalModelResult]) -> list[dict]:
    return [result.row() for _, result in sorted(results.items())]


def registry_payload(results: dict[str, ExternalModelResult], *, model_id: str,
                     model_version: str, training_data_provenance: str,
                     organism: str = "mtbc", independent: bool = True) -> dict:
    # Exactly RegisteredModel's fields, so RegisteredModel(**block) accepts it
    # unmodified. Unexportable scopes are reported by excluded_scopes().
    performance = [r.registry_performance(independent=independent)
                   for _, r in sorted(results.items())
                   if r.registry_exclusion() is None]
    return {
        "model_id": model_id,
        "version": model_version,
        "training_data_provenance": training_data_provenance,
        "organism": organism,
        "validated_cohorts": sorted({r["cohort"] for r in performance}),
        "performance": performance,
    }


def excluded_scopes(results: dict[str, ExternalModelResult]) -> list[dict]:
    return [{"drug": r.drug, "lineage": r.lineage, "verdict": r.comparison.verdict,
             "reason": r.registry_exclusion()}
            for _, r in sorted(results.items()) if r.registry_exclusion()]


def payload(results: dict[str, ExternalModelResult], *, model_id: str,
            model_version: str, training_data_provenance: str,
            organism: str = "mtbc") -> dict:
    return {
        "schema": "mycobench.external-model-benchmark.v1",
        "model": registry_payload(
            results, model_id=model_id, model_version=model_version,
            training_data_provenance=training_data_provenance,
            organism=organism),
        "excluded_scopes": excluded_scopes(results),
        "notes": [
            "Generated by mycobench benchmark-model from externally supplied predictions.",
            "Approval still requires ModelRegistry governance review; this payload only supplies measured performance rows.",
        ],
        "results": [asdict_row(r) for _, r in sorted(results.items())],
    }


def asdict_row(result: ExternalModelResult) -> dict:
    row = result.row()
    row["registry_performance"] = result.registry_performance()
    row["comparison_reasons"] = list(result.comparison.reasons)
    row["dropped"] = [asdict(n) for n in result.dropped]
    row["curve"] = [asdict(p) | {"coverage": p.coverage, "selective_risk": p.selective_risk}
                     for p in result.curve.points]
    return row
