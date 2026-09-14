"""The driver: catalogue + genotypes + phenotypes to effect sizes.

    ingested catalogue ──┬─► coordinate index ──► genotypes (from VCFs)
                         └─► determinant index ──┐
    CRyPTIC phenotypes ───► MIC panels ──────────┴─► stratify ─► effect scan
                                                              ─► mechanism check

What it refuses to do
---------------------
Run against the bundled illustrative catalogue. That file holds 14 hand-picked
entries with no coordinates, so it can neither genotype an isolate nor define a
resistance background — and an analysis that silently used it would produce
effect sizes over almost no data and report them as though they meant
something. ``AnalysisError`` says so instead.

Candidate selection
-------------------
The variants worth testing are the ones the catalogue itself grades as
**uncertain**: on the real release that is 20,843 of 30,699 entries. Those are
the VUS gap, they have coordinates, and their whole point is that nobody knows
what they do. Variants already graded resistant are carried only as the
background the candidates are conditioned on.
"""
from __future__ import annotations

import csv
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from ..cohort import write_rows
from ..phenotypes import DRUG_CODES, REUSE_TABLE
from ..thresholds import ACCEPTED_PHENOTYPE_QUALITY
from . import effects, mechanism, mic
from .genotypes import CoordinateIndex, GenotypeLoad, load_genotypes
from .strata import DeterminantIndex, Isolate, cooccurrence, stratify

#: Catalogue calls whose entries are worth testing. These are the VUS.
CANDIDATE_CALLS = ("indeterminate",)

#: A variant observed in fewer isolates than this is not tested at all.
MIN_CARRIERS = 5

EFFECT_COLUMNS = (
    "variant", "drug", "verdict", "delta", "ci_low", "ci_high",
    "pvalue", "pvalue_holm", "survives_fdr", "n_carriers",
    "n_informative_strata", "n_decidable_pairs", "n_possible_pairs",
    "lineage_diversity", "top_cooccurrence",
    "top_cooccurrence_fraction", "warnings", "reasons",
)
PANEL_COLUMNS = ("drug", "n_isolates", "n_exact", "left_censored",
                 "right_censored", "median_identifiable", "note")


class AnalysisError(RuntimeError):
    """The inputs cannot support the analysis being asked for."""


@dataclass
class AnalysisInputs:
    catalogue: Path
    phenotypes: Path
    reuse_table: Path
    cache_dir: Path
    limit: Optional[int] = None
    min_carriers: int = MIN_CARRIERS
    drugs: Optional[tuple[str, ...]] = None
    cached_only: bool = False
    metadata: Optional[Path] = None
    partition: Optional[str] = None
    independent_clusters: bool = False
    jobs: int = 1
    sample_seed: Optional[int] = None

    def validate(self) -> None:
        if self.limit is not None and self.limit <= 0 or self.min_carriers <= 0:
            raise AnalysisError("limit and min_carriers must be positive")
        if self.jobs < 1:
            raise AnalysisError("jobs must be at least 1")
        if (self.partition or self.independent_clusters) and not self.metadata:
            raise AnalysisError("partition/independent-clusters requires metadata")
        for label, path in (("catalogue", self.catalogue),
                            ("phenotype table", self.phenotypes),
                            ("reuse table", self.reuse_table)):
            if not Path(path).is_file():
                raise AnalysisError(f"{label} not found: {path}")

        data = json.loads(Path(self.catalogue).read_text(encoding="utf-8"))
        if data.get("illustrative", True):
            raise AnalysisError(
                f"{self.catalogue} is the bundled illustrative catalogue: 14 "
                f"hand-picked entries with no coordinates. It cannot genotype "
                f"an isolate or define a resistance background. Run "
                f"'mycobench fetch-catalogue' and point --catalogue at the "
                f"ingested file.")
        with_coordinates = sum(1 for e in data.get("variants", [])
                               if e.get("coordinate_key"))
        if with_coordinates < 1000:
            raise AnalysisError(
                f"{self.catalogue} carries coordinates for only "
                f"{with_coordinates} variant(s); genotyping needs the full "
                f"ingested catalogue")


@dataclass
class AnalysisResult:
    load: Optional[GenotypeLoad] = None
    panels: dict[str, mic.DrugPanel] = field(default_factory=dict)
    scan: Optional[effects.EffectScan] = None
    discriminator: Optional[mechanism.DiscriminatorResult] = None
    n_candidates: int = 0
    n_tested: int = 0
    #: The threshold actually applied, not the module default. A summary that
    #: prints the default while a different one was used misdescribes the run.
    min_carriers: int = MIN_CARRIERS
    skipped: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0

    def panel_rows(self) -> list[dict]:
        rows = []
        for drug, panel in sorted(self.panels.items()):
            profile = panel.profile
            rows.append({
                "drug": drug, "n_isolates": profile.n,
                "n_exact": profile.n_exact,
                "left_censored": f"{profile.left_fraction:.4f}",
                "right_censored": f"{profile.right_fraction:.4f}",
                "median_identifiable": ("yes" if profile.median_identifiable
                                        else "no"),
                "note": profile.why_not_identifiable() or "",
            })
        return rows

    def effect_rows(self) -> list[dict]:
        rows = []
        for effect in (self.scan.effects if self.scan else []):
            top = effect.top_cooccurrence
            rows.append({
                "variant": effect.variant, "drug": effect.drug,
                "verdict": effect.verdict,
                "delta": ("" if effect.delta is None
                          else f"{effect.delta:.4f}"),
                "ci_low": ("" if not effect.interval
                           else f"{effect.interval[0]:.4f}"),
                "ci_high": ("" if not effect.interval
                            else f"{effect.interval[1]:.4f}"),
                "pvalue": ("" if effect.pvalue is None
                           else f"{effect.pvalue:.6f}"),
                "pvalue_holm": ("" if effect.pvalue_holm is None
                                else f"{effect.pvalue_holm:.6f}"),
                "survives_fdr": ("" if effect.survives_fdr is None
                                 else "yes" if effect.survives_fdr else "no"),
                "n_carriers": effect.n_carriers,
                "n_informative_strata": effect.n_informative_strata,
                # The information behind the estimate: a delta of -1.000 over
                # 40 decidable pairs is not the same claim as over 4,000.
                "n_decidable_pairs": effect.n_decidable_pairs,
                "n_possible_pairs": effect.n_possible_pairs,
                # Empty, never "1.00", when lineage was never typed.
                "lineage_diversity": ("" if effect.lineage_diversity is None
                                      else f"{effect.lineage_diversity:.2f}"),
                "top_cooccurrence": top[0] if top else "",
                "top_cooccurrence_fraction": (f"{top[1]:.4f}" if top else ""),
                "warnings": " | ".join(effect.warnings),
                "reasons": " | ".join(effect.reasons),
            })
        return rows

    def summary(self) -> str:
        lines = []
        if self.load:
            lines.append(self.load.describe())
        usable = [d for d, p in self.panels.items()
                  if p.profile.median_identifiable]
        lines.append(
            f"{len(self.panels)} drug panel(s); median identifiable for "
            f"{len(usable)} ({', '.join(sorted(usable)) or 'none'})")
        lines.append(f"{self.n_candidates} candidate variant(s) with at least "
                     f"{self.min_carriers} carrier(s); {self.n_tested} "
                     f"variant-drug pair(s) tested")
        if self.scan:
            lines.append(self.scan.describe())
        if self.discriminator:
            lines.append(f"mechanism discriminator: "
                         f"{self.discriminator.verdict}")
            lines.extend(f"  {r}" for r in self.discriminator.reasons)
        lines.extend(f"NOTE: {n}" for n in self.notes)
        lines.append(f"elapsed {self.elapsed_seconds:.0f}s")
        return "\n".join(lines)


def _read_reuse_rows(path: Path) -> list[dict]:
    text = Path(path).read_text(encoding="utf-8-sig")
    return list(csv.DictReader(text.splitlines()))


def build_panels(phenotype_path: Path,
                 drugs: Optional[Iterable[str]] = None
                 ) -> dict[str, mic.DrugPanel]:
    """MIC panels from the long-form phenotype table, evaluable rows only.

    Only rows at accepted phenotype quality are admitted: measuring a variant
    against a low-quality phenotype measures the phenotype.
    """
    wanted = set(drugs) if drugs else set(DRUG_CODES.values())
    panels: dict[str, mic.DrugPanel] = {}
    with open(phenotype_path, newline="", encoding="utf-8") as handle:
        stream = (line for line in handle
                  if line.strip() and not line.lstrip().startswith("#"))
        for row in csv.DictReader(stream, delimiter="\t"):
            drug = (row.get("drug") or "").strip()
            if drug not in wanted:
                continue
            if (row.get("evaluable") or "").lower() != "yes":
                continue
            if (row.get("quality") or "").upper() not in ACCEPTED_PHENOTYPE_QUALITY:
                continue
            panel = panels.setdefault(drug, mic.DrugPanel(drug))
            panel.add((row.get("sample_id") or "").strip(),
                      (row.get("mic") or "").strip())
    return panels


def candidate_variants(catalogue: Path,
                       observed: dict[str, int],
                       min_carriers: int) -> dict[str, list[str]]:
    """``variant -> drugs`` for uncertain entries with enough carriers."""
    data = json.loads(Path(catalogue).read_text(encoding="utf-8"))
    candidates: dict[str, list[str]] = {}
    for entry in data.get("variants", []):
        if entry.get("call") not in CANDIDATE_CALLS:
            continue
        label = f"{entry['gene']}_{entry['change']}"
        if observed.get(label, 0) < min_carriers:
            continue
        drugs = [d for d in entry.get("drugs", [])]
        if drugs:
            candidates[label] = sorted(set(candidates.get(label, [])) | set(drugs))
    return candidates


def run(inputs: AnalysisInputs) -> AnalysisResult:
    """Run the whole analysis. Every refusal is recorded, not silent."""
    started = time.monotonic()
    inputs.validate()
    result = AnalysisResult(min_carriers=inputs.min_carriers)

    coordinate_index = CoordinateIndex.from_catalogue(inputs.catalogue)
    determinant_index = DeterminantIndex.from_catalogue(inputs.catalogue)
    result.notes.append(coordinate_index.describe())
    result.notes.append(
        f"{determinant_index.n_entries} established determinant(s) define the "
        f"resistance backgrounds")

    rows = _read_reuse_rows(inputs.reuse_table)
    result.load = load_genotypes(rows, coordinate_index, inputs.cache_dir,
                                 limit=inputs.limit, jobs=inputs.jobs,
                                 sample_seed=inputs.sample_seed,
                                 **({"cached_only": True} if inputs.cached_only else {}))
    if not result.load.isolates:
        raise AnalysisError(
            "no isolate was genotyped; nothing can be analysed. Check the "
            f"VCF cache at {inputs.cache_dir} and the reuse table.")

    result.panels = build_panels(inputs.phenotypes, inputs.drugs)
    if not result.panels:
        raise AnalysisError(
            f"no MIC panel could be built from {inputs.phenotypes}; the "
            f"accuracy track needs evaluable phenotypes at quality "
            f"{'/'.join(ACCEPTED_PHENOTYPE_QUALITY)}")

    isolates: list[Isolate] = list(result.load.isolates.values())
    if inputs.metadata:
        from .metadata import apply_metadata
        isolates, notes = apply_metadata(isolates, inputs.metadata, inputs.partition, inputs.independent_clusters)
        result.notes.extend(notes)
        result.load.isolates = {i.isolate_id: i for i in isolates}
    result.notes.extend([
        f"{sum(not i.lineage for i in isolates)} isolates lack lineage; {sum(not i.site for i in isolates)} lack site.",
        "Conditional association is a discovery result, not a causal determinant or susceptibility prediction.",
        "CRyPTIC and WHO catalogue samples may overlap; this scan is not independent predictive validation.",
        "Relatedness is not controlled unless independent-clusters is explicitly enabled with reviewed cluster metadata.",
    ])
    observed = result.load.carriers_per_variant
    candidates = candidate_variants(inputs.catalogue, observed,
                                    inputs.min_carriers)
    result.n_candidates = len(candidates)
    if not candidates:
        result.notes.append(
            f"no uncertain catalogue entry was carried by at least "
            f"{inputs.min_carriers} of the {len(isolates)} genotyped "
            f"isolate(s); widen --limit or lower --min-carriers")
        result.elapsed_seconds = time.monotonic() - started
        return result

    jobs = []
    for variant, drugs in sorted(candidates.items()):
        for drug in drugs:
            panel = result.panels.get(drug)
            if panel is None:
                result.skipped.append(
                    f"{variant}/{drug}: no MIC panel (CRyPTIC does not "
                    f"measure this drug)")
                continue
            stratification = stratify(variant, drug, isolates,
                                      determinant_index)
            cooccur = cooccurrence(variant, drug, isolates, determinant_index)
            jobs.append((stratification, panel, cooccur))

    result.n_tested = len(jobs)
    result.scan = effects.scan(jobs)

    # The exploratory mechanism discriminator, with its own control built in.
    carriers_by_variant = {
        variant: [i.isolate_id for i in isolates if i.carries(variant)]
        for variant in observed
        if observed[variant] >= inputs.min_carriers
        and variant.split("_", 1)[0] in
        (mechanism.EFFLUX_REGULATOR_GENES + ("atpE",))
    }
    if carriers_by_variant:
        result.discriminator = mechanism.validate_discriminator(
            result.panels, carriers_by_variant,
            [i.isolate_id for i in isolates])
    else:
        result.notes.append(
            "mechanism discriminator not run: no efflux-regulator or atpE "
            "variant reached the carrier minimum")

    result.elapsed_seconds = time.monotonic() - started
    return result


def write_outputs(result: AnalysisResult, out_dir: str | Path) -> list[Path]:
    out_dir = Path(out_dir)
    written = [
        write_rows(out_dir / "mic_panels.tsv", result.panel_rows(),
                   PANEL_COLUMNS),
        write_rows(out_dir / "variant_effects.tsv", result.effect_rows(),
                   EFFECT_COLUMNS),
    ]
    payload = {
        "reuse_release": REUSE_TABLE,
        "n_isolates_genotyped": result.load.n_loaded if result.load else 0,
        "allele_match_rate": (result.load.match_rate if result.load else 0.0),
        "n_candidates": result.n_candidates,
        "n_tested": result.n_tested,
        "verdicts": result.scan.counts() if result.scan else {},
        "discriminator": (result.discriminator.verdict
                          if result.discriminator else "not-run"),
        "discriminator_reasons": (result.discriminator.reasons
                                  if result.discriminator else []),
        "skipped": result.skipped[:50],
        "notes": result.notes,
        "elapsed_seconds": round(result.elapsed_seconds, 1),
        "caveat": ("Effect sizes are conditional on resistance background and "
                   "restricted to catalogued coordinates. Nothing here is a "
                   "resistance call, and no variant is promoted to a "
                   "predictive call by this analysis."),
    }
    manifest = out_dir / "analysis_manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    written.append(manifest)
    return written
