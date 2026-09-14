"""Deterministic synthetic demonstration of the complete reporting workflow.

The fixture is deliberately rich and conspicuously synthetic.  It exercises
the production analysis and evidence-governance APIs, then adds cohort context
through the same validated report-context contract used by real deployments.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from .adapters.base import EngineReport
from .catalogue.regulatory import RegulatoryCatalogue, RegulatoryRegion
from .core.context import SampleContext
from .core.models import (
    Call, DrugEvidence, EngineRef, ExpressionEvidence, InSilicoPrediction,
    Lane, QCFinding, Tier,
)
from .core.pipeline import Myconductor
from .federated.model_registry import ModelRegistry, RegisteredModel
from .io.panel import Panel
from .modules.epistasis import EpistasisRule, EpistasisTable
from .modules.mic_evidence import MICPrediction
from .modules.structural_annotation import (
    StructuralAnnotation, StructuralAnnotationTable,
)
from .reporting.fhir import to_fhir_json
from .reporting.html_report import render_html
from .reporting.json_report import atomic_json, report_dict
from .reporting.report_context import load_report_context
from .reporting.semantic import write_jsonld

ROOT = Path(__file__).resolve().parent.parent
VCF = ROOT / "myconductor" / "data" / "example_input.vcf"
MASK = ROOT / "myconductor" / "data" / "example_callable.tsv"

PNCA_KEY = "NC_000962.3:NC_000962.3:2288725:A>C"
RV0678_KEY = "NC_000962.3:NC_000962.3:4247429:T>G"
DEMO_COHORT_SIZE = 12
DEMO_LINEAGES = ("Lineage 1", "Lineage 2", "Lineage 3", "Lineage 4")
DEMO_SITES = (
    "SYNTHETIC-SITE-01",
    "SYNTHETIC-SITE-02",
    "SYNTHETIC-SITE-03",
    "SYNTHETIC-SITE-04",
)
DEMO_EXTRA_VARIANTS = (
    ("NC_000962.3", 761110, "A", "G", "rpoB", "L430P", "coding", "missense", 0.84, 62, "10,52"),
    ("NC_000962.3", 761115, "G", "T", "rpoB", "D435V", "coding", "missense", 0.92, 76, "6,70"),
    ("NC_000962.3", 761172, "T", "C", "rpoB", "I491T", "coding", "missense", 0.36, 88, "56,32"),
    ("NC_000962.3", 2155104, "C", "T", "katG", "R463L", "coding", "missense", 0.95, 69, "3,66"),
    ("NC_000962.3", 2155204, "T", "A", "katG", "W328R", "coding", "missense", 0.28, 72, "52,20"),
    ("NC_000962.3", 2288759, "C", "T", "pncA", "H57D", "coding", "missense", 0.81, 54, "10,44"),
    ("NC_000962.3", 4247520, "C", "A", "embB", "Q497R", "coding", "missense", 0.42, 63, "37,26"),
    ("NC_000962.3", 7582, "C", "T", "gyrA", "A90V", "coding", "missense", 0.77, 73, "17,56"),
    ("NC_000962.3", 4247385, "G", "A", "Rv0678", "R94Q", "coding", "missense", 0.69, 58, "18,40"),
    ("NC_000962.3", 4247467, "A", "T", "Rv0678", "Q133*", "coding", "stop_gained", 0.23, 61, "47,14"),
    ("NC_000962.3", 2715328, "G", "A", "eis", "c.-37G>A", "promoter", "upstream", 0.88, 52, "6,46"),
    ("NC_000962.3", 1473246, "A", "G", "rrs", "a1401g", "rrna", "rrna", 1.0, 44, "0,44"),
)


def synthetic_report_context() -> dict:
    """A full, deterministic programme dataset for visual inspection only."""
    drugs = [
        "rifampicin", "isoniazid", "moxifloxacin", "bedaquiline",
        "linezolid", "pyrazinamide", "ethambutol", "amikacin",
        "clofazimine", "pretomanid",
    ]
    baseline = []
    for index, drug in enumerate(drugs):
        n = 360 - index * 17
        resistant = 104 - index * 7
        susceptible = n - resistant
        false_susceptible = index % 4
        false_resistant = (index + 2) % 5
        baseline.append({
            "drug": drug, "estimable": "yes", "n_isolates": str(n),
            "n_answered": str(n - index * 2),
            "coverage": f"{(n-index*2)/n:.4f}",
            "error_rate": f"{(false_susceptible+false_resistant)/(n-index*2):.4f}",
            "n_false_susceptible": str(false_susceptible),
            "n_false_resistant": str(false_resistant),
            "sensitivity": f"{(resistant-false_susceptible)/resistant:.4f}",
            "specificity": f"{(susceptible-false_resistant)/susceptible:.4f}",
            "sensitivity_ci": "synthetic demonstration interval",
            "specificity_ci": "synthetic demonstration interval",
            "ppv": f"{0.91-index*0.012:.4f}", "npv": f"{0.98-index*0.006:.4f}",
            "resistance_prevalence": f"{resistant/n:.4f}",
            "sensitivity_basis_n": str(resistant),
            "specificity_basis_n": str(susceptible),
            "notes": "Synthetic workflow fixture; no clinical performance claim.",
        })
    lineage_rows = []
    for d_index, drug in enumerate(drugs):
        for l_index, lineage in enumerate(DEMO_LINEAGES):
            called = 42 + d_index * 4 + l_index * 7
            lineage_rows.append({
                "drug": drug, "lineage": lineage, "n_called": str(called),
                "sensitivity": f"{0.94-d_index*.012-l_index*.008:.3f}",
                "specificity": f"{0.97-d_index*.008+l_index*.003:.3f}",
                "vme_rate": f"{0.02+d_index*.004+l_index*.003:.3f}",
                "me_rate": f"{0.03+d_index*.003:.3f}", "powered": "yes",
            })
    prevalence = []
    for p_index, period in enumerate(("2025 Q1", "2025 Q2", "2025 Q3", "2025 Q4", "2026 Q1")):
        for drug, base in (("rifampicin", 31), ("isoniazid", 39), ("moxifloxacin", 14)):
            total = 120 + p_index * 9
            prevalence.append({
                "period": period, "drug": drug,
                "lineage": "Lineage 2" if p_index % 2 else "Lineage 4",
                "geography": "North Central synthetic network" if p_index < 3 else "South West synthetic network",
                "resistant": base + p_index * 2, "total": total,
                "source": "Synthetic programme register",
            })
    return {
        "schema": "myconductor.report-context.v1",
        "baseline_rows": baseline, "lineage_rows": lineage_rows,
        "measurability": [
            {"drug": "delamanid", "measurable": "no", "missing_side": "phenotype",
             "consequence": "No paired synthetic phenotype was supplied."},
            {"drug": "cycloserine", "measurable": "no", "missing_side": "both",
             "consequence": "Neither side is represented in this fixture."},
        ],
        "error_trend": [8.4, 7.8, 7.1, 6.7, 6.2, 5.9, 5.5, 5.1],
        "prevalence_rows": prevalence,
        "target_rows": [
            {"target": "DprE1", "essentiality": "Supported in external studies", "druggability": "Biochemical and inhibitor evidence", "human_homology": "Low reported homology", "resistance_liability": "On-target resistance observed", "evidence": "Synthetic summary representing a curated CRISPR-interference and biochemical evidence bundle.", "source": "Synthetic target-evidence registry"},
            {"target": "MmpL3", "essentiality": "Supported", "druggability": "Multiple chemical series", "human_homology": "No close human orthologue reported", "resistance_liability": "Resistance mutations reported", "evidence": "Synthetic target card for workflow testing.", "source": "Synthetic target-evidence registry"},
            {"target": "PknB", "essentiality": "Context dependent evidence", "druggability": "Early-stage", "human_homology": "Kinase selectivity requires review", "resistance_liability": "Insufficient evidence", "evidence": "Included to test incomplete target evidence without converting gaps to scores.", "source": "Synthetic target-evidence registry"},
        ],
        "watchlist_rows": [
            {"drug": "bedaquiline", "variant": "Rv0678_L117R", "unresolved_isolates": 18, "contributing_sites": 4, "status": "Phenotype review pending", "evidence_gaps": ["paired MIC", "lineage-balanced replication"], "limitations": "Synthetic privacy-gated aggregate; no isolate records.", "source": "Synthetic multi-site watch list"},
            {"drug": "pyrazinamide", "variant": "pncA_D12A", "unresolved_isolates": 11, "contributing_sites": 3, "status": "Functional assay proposed", "evidence_gaps": ["enzyme activity", "replicate phenotype"], "source": "Synthetic multi-site watch list"},
            {"drug": "amikacin", "variant": "eis_c.-14C>T", "unresolved_isolates": 7, "contributing_sites": 2, "status": "Catalogue reconciliation", "evidence_gaps": ["catalogue-version alignment"], "source": "Synthetic multi-site watch list"},
        ],
        "external_benchmark_rows": [
            {"model_id": "synthetic-vus-transformer", "model_version": "0.4.0", "drug": "rifampicin", "lineage": "Lineage 2", "cohort": "synthetic-independent-evaluation", "source": "Synthetic external prediction file", "n_predictions": 84, "n_evaluable": 80, "n_called": 76, "n_dropped": 4, "call_rate": 0.95, "error_rate": 0.035, "baseline_coverage": 0.91, "baseline_error_rate": 0.082, "matched_coverage": 0.925, "matched_error_rate": 0.041, "sensitivity": 0.96, "specificity": 0.94, "verdict": "beats-baseline", "registry_ready": "yes", "notes": "Synthetic model beats the catalogue baseline at matched coverage; governance review still required."},
            {"model_id": "synthetic-vus-transformer", "model_version": "0.4.0", "drug": "pyrazinamide", "lineage": "Lineage 4", "cohort": "synthetic-independent-evaluation", "source": "Synthetic external prediction file", "n_predictions": 79, "n_evaluable": 74, "n_called": 65, "n_dropped": 5, "call_rate": 0.878, "error_rate": 0.092, "baseline_coverage": 0.86, "baseline_error_rate": 0.088, "matched_coverage": 0.878, "matched_error_rate": 0.092, "sensitivity": 0.89, "specificity": 0.91, "verdict": "no-better-than-baseline", "registry_ready": "no", "notes": "Synthetic comparison does not beat the catalogue at matched coverage; do not register for approval."},
            {"model_id": "candidate-paper-output", "model_version": "paper-table", "drug": "bedaquiline", "lineage": "All lineages", "cohort": "synthetic-independent-evaluation", "source": "Synthetic paper-derived table", "n_predictions": 9, "n_evaluable": 9, "n_called": 9, "n_dropped": 0, "call_rate": 1.0, "error_rate": "", "baseline_coverage": 0.84, "baseline_error_rate": 0.11, "matched_coverage": 1.0, "matched_error_rate": "", "sensitivity": "", "specificity": "", "verdict": "underpowered", "registry_ready": "no", "notes": "Fewer than the minimum answered queries; shown as evidence gap, not as a failed model."}
        ],
        "validation": {
            "outcomes": [{"k": "Confirmed resistant", "v": 14}, {"k": "Confirmed susceptible", "v": 9}, {"k": "Inconclusive", "v": 4}],
            "timeline": [
                {"date": "2025-10", "title": "Bedaquiline VUS batch reviewed", "detail": "Six synthetic isolates received paired MIC testing."},
                {"date": "2025-12", "title": "pncA functional assay completed", "detail": "Synthetic activity measurements were appended to the local ledger."},
                {"date": "2026-02", "title": "Discordance adjudication", "detail": "Two synthetic engine disagreements were reviewed without majority voting."},
                {"date": "2026-04", "title": "Retraction audit exercised", "detail": "One synthetic validation was retracted and retained in history."},
            ],
        },
        "federated_sites": [
            {"site": "Synthetic reference laboratory A", "submissions": 148, "k": 10, "status": "Accepted"},
            {"site": "Synthetic reference laboratory B", "submissions": 121, "k": 10, "status": "Accepted"},
            {"site": "Synthetic sequencing centre C", "submissions": 86, "k": 8, "status": "Accepted with threshold note"},
            {"site": "Synthetic sentinel site D", "submissions": 34, "k": 5, "status": "Accepted"},
        ],
        "audit_entries": [
            {"timestamp": "2026-04-01T09:10:00Z", "action": "catalogue_imported", "actor": "Synthetic curator", "target": "Catalogue demonstration 2026.1", "hash": "91d1e6e4a2634b88"},
            {"timestamp": "2026-04-02T12:00:00Z", "action": "model_approved", "actor": "Synthetic model review board", "target": "VUS classifier 2.1", "hash": "38b7c4532258c2dd"},
            {"timestamp": "2026-04-05T15:20:00Z", "action": "validation_recorded", "actor": "Synthetic laboratory A", "target": "pncA D12A", "hash": "ac919b0982ce1dd1"},
            {"timestamp": "2026-04-08T08:45:00Z", "action": "validation_retracted", "actor": "Synthetic reviewer", "target": "pncA D12A record 3", "hash": "e46307b8ad8e6f07"},
        ],
        "reference_method": "Synthetic paired MGIT and microdilution demonstration dataset",
    }


def _model_registry() -> ModelRegistry:
    performance = []
    for index, lineage in enumerate(DEMO_LINEAGES):
        performance.append({
            "drug": "pyrazinamide", "lineage": lineage,
            "cohort": f"synthetic-independent-evaluation-{index + 1}",
            "source": "synthetic-evaluation/1",
            "independent": True, "sensitivity": .93 - index * .01,
            "specificity": .96 - index * .005,
            "call_rate": .95 - index * .01, "error_rate": .04 + index * .006,
            "baseline": {
                "source": "Synthetic catalogue baseline",
                "coverage": .91 - index * .015,
                "error_rate": .08 + index * .004,
            },
        })
    registry = ModelRegistry()
    registry.submit(RegisteredModel("vus-classifier", "2.1", "synthetic-training/1", "mtbc",
                                    [row["cohort"] for row in performance], performance),
                    "Synthetic submitter", "Exercise governed model intake")
    registry.transition("vus-classifier", "2.1", "reviewed", "Synthetic reviewer", "Synthetic documentation reviewed")
    registry.transition("vus-classifier", "2.1", "approved", "Synthetic reviewer", "Synthetic fixture beats its declared baseline")
    return registry


def _sample_vcf(sample_id: str, sample_index: int = 0) -> Path:
    """Create a valid per-sample synthetic VCF with cohort-level variety."""
    target_dir = Path(tempfile.gettempdir()) / "myconductor-full-demo-inputs"
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / f"{sample_id}.vcf"
    lines = VCF.read_text(encoding="utf-8").splitlines()
    out: list[str] = []
    for line in lines:
        if line.startswith("##sampleID="):
            out.append(f"##sampleID={sample_id}")
        elif line.startswith("#CHROM"):
            parts = line.split("\t")
            parts[-1] = sample_id
            out.append("\t".join(parts))
        elif line and not line.startswith("#"):
            out.append(line)
        else:
            out.append(line)

    for variant_index, row in enumerate(DEMO_EXTRA_VARIANTS):
        if (sample_index + variant_index) % 4 == 1:
            continue
        chrom, pos, ref, alt, gene, change, region, effect, af, depth, ad = row
        shifted_af = max(0.05, min(1.0, af - ((sample_index + variant_index) % 3) * 0.04))
        info = f"GENE={gene};AACHANGE={change};REGION={region};EFFECT={effect}"
        out.append(
            f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t57\tPASS\t{info}\tAF:DP:AD\t{shifted_af:.2f}:{depth}:{ad}"
        )
    target.write_text("\n".join(out) + "\n", encoding="utf-8")
    return target


def build_full_demo(sample_index: int = 0):
    """Run one synthetic sample through the advanced evidence paths."""
    sample_id = f"MTB-DEMO-{sample_index + 1:03d}"
    isolate_id = f"ISO-DEMO-{sample_index + 1:03d}"
    site_id = DEMO_SITES[sample_index % len(DEMO_SITES)]
    lineage = DEMO_LINEAGES[sample_index % len(DEMO_LINEAGES)]
    collected_day = 1 + sample_index
    context = SampleContext(
        sample_id, isolate_id, site_id, "mtbc",
        assay="targeted-amplicon", lineage=lineage,
        specimen_id=f"SYNTHETIC-SPECIMEN-{sample_index + 1:03d}",
        collected_at=f"2026-03-{collected_day:02d}",
        qc=[QCFinding(name, "pass", "Synthetic pass used to exercise the workflow")
            for name in ("species_confirmation", "contamination", "mapping_quality", "mixed_infection")],
    )
    panel = Panel(
        "Synthetic comprehensive TB panel", "2026.1", "targeted-amplicon",
        "Synthetic manufacturer declaration", frozenset({
            "katG", "fabG1", "inhA", "rpoB", "pncA", "embB", "gyrA", "gyrB",
            "rplC", "rrl", "ddn", "fbiA", "fbiB", "fbiC", "fgd1", "atpE",
            "Rv0678", "rrs", "eis", "rpoC"}), organism="mtbc",
        reference_assembly="NC_000962.3", verified=False,
        note="Synthetic declaration created only to exercise panel restriction and reporting.")
    tb_profiler_call = Call.RESISTANT if sample_index % 3 != 1 else Call.SUSCEPTIBLE
    mykrobe_call = Call.SUSCEPTIBLE if sample_index % 2 == 0 else Call.RESISTANT
    mykrobe_tier = Tier.CATALOGUED if mykrobe_call is Call.RESISTANT else Tier.NONE
    engines = [
        EngineReport(EngineRef("TB-Profiler", "4.4", "Synthetic catalogue", "2026.1"),
                     sample_id=context.sample_id, evidence=[DrugEvidence(
                         "isoniazid", tb_profiler_call, Tier.CATALOGUED, Lane.ENGINE,
                         rationale="Synthetic TB-Profiler resistance assertion.")]),
        EngineReport(EngineRef("Mykrobe", "0.13", "Synthetic panel", "2026.1"),
                     sample_id=context.sample_id, evidence=[DrugEvidence(
                         "isoniazid", mykrobe_call, mykrobe_tier, Lane.ENGINE,
                         rationale="Synthetic Mykrobe assertion used to exercise engine discordance.")]),
    ]
    epistasis = EpistasisTable([EpistasisRule(
        "synthetic-rpob-rpoc-1", "rifampicin", "rpoB", "rpoC", "compensatory",
        "Synthetic co-observation used to exercise annotation-only epistasis reporting.",
        "Synthetic reviewed interaction table", "mtbc", "NC_000962.3",
        "Synthetic curator", "2026-04-01")], "synthetic-2026.1", True)
    structural = StructuralAnnotationTable([StructuralAnnotation(
        PNCA_KEY, "binding_pocket", "Synthetic report of a possible local packing change",
        "Synthetic structural archive", "mtbc", "NC_000962.3", "2026.1",
        "2026-04-01", confidence=.72, ligand_distance=4.8, distance_unit="angstrom",
        ligand_reference="Synthetic pyrazinamide analogue", validation_status="pending")])
    regulatory = RegulatoryCatalogue([RegulatoryRegion(
        "synthetic-eis-promoter-1", "eis upstream regulatory region", "mtbc",
        "NC_000962.3", "PROMOTER", ["eis"], ["amikacin"], "CANDIDATE",
        "Synthetic regulatory catalogue")], "synthetic-2026.1", True)
    expression = [ExpressionEvidence(
        "Rv0678", "rna_seq", 7.2 + sample_index * 0.25, 1.0, 7.2 + sample_index * 0.25,
        f"Synthetic RNA-seq run {11 + sample_index}",
        context.sample_id, f"SYN-EXP-{sample_index + 1:03d}", context.isolate_id, context.site_id,
        context.organism, "normalised counts", "2026-04-01", confidence=.88,
        reported_conclusion="elevated", interpretation_source="Synthetic analysis plan",
        linked_variant_keys=(RV0678_KEY,), linkage_source="Same synthetic isolate")]
    mic_value = 2.0 + sample_index * 0.1
    mic = [MICPrediction(
        f"SYN-MIC-PRED-{sample_index + 1:03d}", context.sample_id, context.isolate_id, context.site_id,
        context.organism, "pyrazinamide", mic_value, "mg/L", "Synthetic external MIC model",
        "1.3", "Synthetic model output", "2026-04-01", 1.0,
        "Synthetic breakpoint table", interval=(round(mic_value - 0.6, 1), round(mic_value + 0.8, 1)), interval_kind="prediction",
        interval_level=.95)]
    predictions = [InSilicoPrediction(
        PNCA_KEY, "pyrazinamide", "vus-classifier", "2.1", "synthetic-training/1",
        "resistant", "Synthetic predictions/run-7", context.sample_id,
        context.isolate_id, context.site_id, context.organism, "2026-04-01",
        confidence=.91, feature_attributions={"conservation": .36, "structure": .24})]
    menu = [
        {"id": "review-v1", "action": "review_variant", "available": True, "cost": 20.0, "currency": "USD", "turnaround_hours": 24},
        {"id": "phenotype-v1", "action": "confirm_phenotype", "available": True, "cost": 45.0, "currency": "USD", "turnaround_hours": 72},
        {"id": "coverage-v1", "action": "review_coverage", "available": True, "cost": 5.0, "currency": "USD", "turnaround_hours": 4},
        {"id": "engine-v1", "action": "review_engine_evidence", "available": True, "cost": 10.0, "currency": "USD", "turnaround_hours": 8},
    ]
    return Myconductor(platform="illumina", demo_mode=True).analyze(
        _sample_vcf(sample_id, sample_index), mask_path=MASK, panel=panel,
        context=context, engine_reports=engines,
        include_population_structure=True, epistasis_table=epistasis,
        mic_predictions=mic, structural_annotations=structural,
        regulatory_catalogue=regulatory, expression_evidence=expression,
        in_silico_predictions=predictions, model_registry=_model_registry(),
        test_menu=menu, follow_up_budget=120.0)


def build_full_demo_cohort(count: int = DEMO_COHORT_SIZE):
    """Run a deterministic multi-sample cohort through the same workflow."""
    if count < 2:
        raise ValueError("full synthetic cohort demo requires at least two samples")
    return [build_full_demo(index) for index in range(count)]


def write_full_demo(output_dir: str | Path) -> dict[str, Path]:
    """Generate the complete offline demo and its interoperability artifacts."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    context_path = output / "synthetic-report-context.json"
    atomic_json(context_path, synthetic_report_context())
    context = load_report_context(context_path)
    reports = build_full_demo_cohort()
    report = reports[0]
    paths = {
        "html": output / "myconductor-full-synthetic-demo.html",
        "json": output / "myconductor-full-synthetic-demo.json",
        "cohort_json": output / "myconductor-full-synthetic-demo.cohort.json",
        "fhir": output / "myconductor-full-synthetic-demo.fhir.json",
        "jsonld": output / "myconductor-full-synthetic-demo.jsonld",
        "context": context_path,
    }
    paths["html"].write_text(
        render_html(report, extra_reports=reports[1:], **context),
        encoding="utf-8",
    )
    atomic_json(paths["json"], report_dict(report))
    atomic_json(paths["cohort_json"], [report_dict(item) for item in reports])
    paths["fhir"].write_text(to_fhir_json(report) + "\n", encoding="utf-8")
    write_jsonld(report, paths["jsonld"])
    manifest = {name: str(path) for name, path in paths.items()}
    (output / "README.txt").write_text(
        "MYCONDUCTOR FULL SYNTHETIC DEMO\n\n"
        "All observations, metrics, identities, sites, dates, and evidence sources in this directory are synthetic.\n"
        "Open myconductor-full-synthetic-demo.html or serve this directory over localhost.\n\n"
        + json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return paths
