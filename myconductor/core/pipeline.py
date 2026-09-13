"""The conductor: wires every stage into one end-to-end analysis.

    input + callable mask + engine reports
        -> QC
        -> multi-lane routing
        -> minority-allele assessment
        -> coverage-gated reconciliation
        -> guideline eligibility
        -> VUS prioritisation
        -> report

Two structural changes from the original pipeline
-------------------------------------------------
**The discovery loop is gone from this file.** It used to fire whenever one
patient had no adequate regimen, putting candidate compounds in a clinical
report. Target discovery is a cohort-level research activity and now lives in
``modules.discovery``, which this module does not import. There is no code path
from a patient's report to a drug-discovery run.

**Synthetic annotators cannot reach a report.** If the VUS annotator reports
itself synthetic, the constructor refuses unless ``demo_mode=True`` — and
``demo_mode`` stamps every rendered report and FHIR bundle.
"""
from __future__ import annotations

from datetime import datetime, timezone
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
from typing import Iterable, Optional, Sequence

from .. import __version__
from ..adapters.base import EngineReport, concordance, merge as merge_evidence
from ..catalogue.profile import OrganismProfile, load_profile
from ..core.models import AnalysisReport, EngineRef, Provenance, Call, DrugEvidence, Lane, Tier
from ..io import qc as qc_module
from ..io.adapter import load
from ..io.callable_mask import CallableMask
from ..modules.catalogue import CatalogueModule
from ..modules.efflux import EffluxRegulatoryModule
from ..modules.features import AnnotatorProtocol
from ..federated.vus_feedback import LocalValidationStore
from ..modules.heteroresistance import HeteroresistanceAssessor, heteroresistance_evidence
from ..modules.local_validation import LocalValidationModule
from ..modules.synthesis import (
    EligibilityAssessor,
    EvidenceReconciler,
    collect_discordances,
)
from ..modules.vus_workbench import VUSWorkbench
from .router import TriageRouter
from .context import SampleContext, InterpretationPolicy
from ..io.phenotypes import PhenotypeObservation


class Myconductor:
    def __init__(
        self,
        profile: Optional[OrganismProfile] = None,
        depth_floor: int = 10,
        callable_fraction_floor: float = 0.95,
        platform: Optional[str] = None,
        demo_mode: bool = False,
        annotator: Optional[AnnotatorProtocol] = None,
        router: Optional[TriageRouter] = None,
        reconciler: Optional[EvidenceReconciler] = None,
        assessor: Optional[EligibilityAssessor] = None,
        het_assessor: Optional[HeteroresistanceAssessor] = None,
        include_tier2_loci: bool = False,
        local_validation_store: Optional[LocalValidationStore] = None,
        organism: Optional[str] = None,
        catalogue_path: Optional[str | Path] = None,
        policy: Optional[InterpretationPolicy] = None,
    ):
        """``organism`` selects a bundled organism profile + catalogue pair
        (see ``catalogue/profile.py::BUNDLED_ORGANISMS`` and
        ``modules/catalogue.py::BUNDLED_CATALOGUES``) — default ``None``
        keeps the original MTBC profile. It is ignored wherever ``profile``
        is passed explicitly, since an explicit profile already says which
        organism this is. Only ``mtbc`` ships with any real-world track
        record here; every other bundled organism (e.g. ``mabscessus``) is
        illustrative in the same sense the bundled MTBC catalogue is —
        ``OrganismProfile.ships_validated`` is ``False`` for all of them.
        """
        self.profile = profile or load_profile(organism=organism)
        self.depth_floor = depth_floor
        self.callable_fraction_floor = callable_fraction_floor
        self.platform = platform
        self.demo_mode = demo_mode
        self.local_validation_store = local_validation_store

        self.policy = policy or InterpretationPolicy()
        catalogue = CatalogueModule(path=catalogue_path, organism=organism or self.profile.name)
        efflux = EffluxRegulatoryModule(self.profile)
        known = set(catalogue.resolved_labels)
        self.workbench = VUSWorkbench(annotator=annotator, known=known,
                                      profile=self.profile)
        if self.workbench.synthetic and not demo_mode:
            raise ValueError(
                f"{self.workbench.model_name} produces synthetic, "
                f"hash-derived values with no biological meaning. Pass "
                f"demo_mode=True to permit it; every report will then be "
                f"stamped as synthetic."
            )

        local_validation = (LocalValidationModule(local_validation_store)
                           if local_validation_store is not None else None)
        self.router = router or TriageRouter(
            catalogue=catalogue, efflux=efflux, extra_lanes=[self.workbench],
            local_validation=local_validation)
        self.catalogue = self.router.catalogue
        self._custom_reconciler = reconciler
        self.reconciler = reconciler or EvidenceReconciler(
            profile=self.profile, depth_floor=depth_floor,
            callable_fraction_floor=callable_fraction_floor,
            include_tier2_loci=include_tier2_loci)
        self.assessor = assessor or EligibilityAssessor(self.profile)
        self.het = het_assessor or HeteroresistanceAssessor(
            depth_floor=depth_floor, default_platform=platform)

    # -- mask construction ------------------------------------------------
    def _mask_from_path(self, mask_path: str | Path) -> CallableMask:
        path = Path(mask_path)
        suffix = path.suffix.lower()
        if suffix == ".bed":
            return CallableMask.from_bed(
                path, locus_lengths=self.profile.locus_lengths())
        if suffix in (".tsv", ".txt", ".csv"):
            return CallableMask.from_tsv(path)
        if suffix in (".gvcf", ".g.vcf", ".vcf"):
            raise ValueError(
                "gVCF masks need locus spans from the reference annotation; "
                "build the mask with CallableMask.from_gvcf(path, locus_spans) "
                "and pass it via mask=, or supply a depth table instead."
            )
        raise ValueError(f"unsupported mask format {suffix!r}")

    # -- the analysis -----------------------------------------------------
    def analyze(
        self,
        input_path: str | Path,
        mask: Optional[CallableMask] = None,
        mask_path: Optional[str | Path] = None,
        sample: Optional[str] = None,
        platform: Optional[str] = None,
        engine_reports: Sequence[EngineReport] = (),
        context: Optional[SampleContext] = None,
        phenotypes: Sequence[PhenotypeObservation] = (),
        test_menu: Sequence[dict] = (),
        follow_up_budget: Optional[float] = None,
        include_population_structure: bool = False,
        population_tolerance: float = 0.05,
        population_lineage_hints: Sequence[dict] = (),
        epistasis_table=None,
        mic_predictions=(),
        structural_annotations=None,
        regulatory_catalogue=None,
        expression_evidence=(),
        in_silico_predictions=(),
        model_registry=None,
        watchlist=None,
    ) -> AnalysisReport:
        if population_lineage_hints and not include_population_structure:
            raise ValueError("population lineage hints require population-structure analysis")
        adapted = load(
            input_path, depth_floor=self.depth_floor, sample=sample,
            platform=platform or self.platform,
        )

        if mask is None and mask_path is not None:
            mask = self._mask_from_path(mask_path)
        mask = mask or CallableMask.absent()

        if context is not None and context.sample_id != adapted.sample_id:
            raise ValueError("context sample_id does not match variant input")
        context = context or SampleContext(adapted.sample_id, adapted.sample_id,
                                          "unspecified-site", self.profile.name)
        if context.organism != self.profile.name:
            raise ValueError("context organism does not match profile")
        if adapted.assembly != self.profile.reference_assembly:
            raise ValueError(f"reference assembly mismatch: {adapted.assembly} vs "
                             f"{self.profile.reference_assembly}; liftover is not implemented")
        engine_reports = list(engine_reports)
        seen_engines = set()
        for report in engine_reports:
            name = report.engine.name
            if name in seen_engines:
                raise ValueError(f"duplicate engine report: {name}")
            seen_engines.add(name)
            expected = context.engine_sample_ids.get(name, context.sample_id)
            if report.sample_id != expected:
                raise ValueError(f"{name}: sample identity mismatch or missing "
                                 f"({report.sample_id!r} vs {expected!r})")
            if report.species and not _species_matches(report.species, self.profile.name):
                raise ValueError(f"{name}: species {report.species!r} conflicts with profile")
            if any(v.identity.assembly != self.profile.reference_assembly
                   for v in report.variants):
                raise ValueError(f"{name}: variant assembly mismatch")
        # Imported organism states are attributed and cannot overwrite a
        # contradictory observation supplied by the laboratory.
        states = {}
        for er in engine_reports:
            for key, value in er.context_findings.items():
                previous = states.get(key, getattr(context, key, None))
                if previous not in (None, "unknown", value):
                    raise ValueError(f"conflicting organism context: {key}")
                states[key] = value
        if states:
            context = replace(context, **states)
        # Preserve both pass and fail findings. A supplied pass cannot erase a
        # failure reported by the engine.
        imported_qc = [q for er in engine_reports for q in er.qc]
        context_qc = {q.check: q for q in context.qc}
        for q in imported_qc:
            if q.status == "fail" or q.check not in context_qc:
                context_qc[q.check] = q
        context = replace(context, qc=list(context_qc.values()))
        if len({p.observation_id for p in phenotypes}) != len(phenotypes):
            raise ValueError("duplicate phenotype observation IDs")
        engine_masks = [r.mask for r in engine_reports if r.mask is not None]
        if engine_masks:
            mask = CallableMask.merge(mask, *engine_masks)

        # -- evidence ------------------------------------------------------
        # External variants must also reach the uncertainty workbench. Deduplicate
        # by identity without inventing equivalence between HGVS dialects.
        variants_by_key = {v.key(): v for v in adapted.variants}
        for report in engine_reports:
            for v in report.variants:
                variants_by_key.setdefault(v.key(), v)
        outcome = self.router.route(list(variants_by_key.values()))
        evidence = list(outcome.evidence) + merge_evidence(engine_reports)
        regulatory_findings = []
        if regulatory_catalogue is not None:
            for variant in variants_by_key.values():
                for entry in regulatory_catalogue.matches(variant, self.profile.name):
                    regulatory_findings.append(dict(asdict(entry), variant_key=variant.key(),
                        catalogue_version=regulatory_catalogue.version, illustrative=regulatory_catalogue.illustrative))
                    if entry.evidence_tier == "CATALOGUED":
                        for drug in entry.drug_associations:
                            if self.profile.supports_drug(drug):
                                evidence.append(DrugEvidence(drug, Call.INDETERMINATE, Tier.INFERRED,
                                    Lane.EFFLUX_REGULATORY, variant=variant.identity,
                                    rationale=f"Regulatory-region record {entry.id} ({entry.evidence_tier}); region membership does not establish expression or resistance.",
                                    scope="variant", catalogue_version=regulatory_catalogue.version,
                                    rule_id=entry.id, metadata={"regulatory_region": asdict(entry)}))

        all_variants = list(adapted.variants)
        for report in engine_reports:
            all_variants.extend(report.variants)

        lineage = next((r.lineage for r in engine_reports if r.lineage), None)
        het = self.het.assess(all_variants, evidence, lineage=lineage)
        if self.het.calibration is not None and self.het.prior_source is not None:
            evidence = evidence + heteroresistance_evidence(het, all_variants)
        evidence += [p.evidence(context) for p in phenotypes]
        qc = qc_module.assess(adapted, mask, self.profile, self.depth_floor)
        overrides = {q.check: q for q in context.qc}
        qc = [overrides.pop(q.check, q) for q in qc] + list(overrides.values())
        for report in engine_reports:
            qc.extend(report.qc)
        # Caller-level rejected records cannot disappear into a negative screen.
        # Failures in identity/contamination/etc invalidate genomic conclusions,
        # but an independent, matched laboratory phenotype remains visible.
        blockers = [q for q in qc if q.status == "fail" and q.check != "callable_mask"]
        if adapted.rejected:
            blockers.append(qc_module.QCFinding("rejected_records", "fail",
                                              "input records rejected; incomplete genomic evidence"))
            qc.append(blockers[-1])
        if blockers:
            evidence = [replace(e, call=Call.INDETERMINATE,
                                limitations=e.limitations + ("genomic QC failed",))
                        if e.tier.value != "phenotypic" else e for e in evidence]
            evidence += [DrugEvidence(drug=d, call=Call.NO_CALL, tier=Tier.NONE,
                                      lane=Lane.NONE, rationale="Genomic QC failed: " +
                                      "; ".join(q.check for q in blockers))
                         for d in self.profile.drugs]
        reconciler = self._custom_reconciler or EvidenceReconciler(
            profile=self.profile, depth_floor=self.depth_floor,
            callable_fraction_floor=self.callable_fraction_floor,
            include_tier2_loci=self.reconciler.include_tier2_loci,
            policy=self.policy, context=context, catalogue_engine=self.catalogue.engine,
            catalogue_sha256=self.catalogue.sha256,
            illustrative=self.catalogue.is_illustrative)
        results = reconciler.reconcile(evidence, mask)
        from ..modules.mic_evidence import reconcile_predictions
        quantitative_findings = reconcile_predictions(results, mic_predictions, context)
        from ..modules.expression_evidence import reconcile_expression
        expression_findings = reconcile_expression(results, expression_evidence, context)
        eligibility = self.assessor.assess(results)
        vus = self.workbench.priorities(list(variants_by_key.values()))
        mechanism_queue = self.router.efflux.mechanism_queue(adapted.variants)
        structural_records = []
        if structural_annotations is not None:
            structural_records, structural_hypotheses = structural_annotations.attach(
                vus, self.profile.name, self.profile.reference_assembly)
            mechanism_queue.extend(structural_hypotheses)
        from ..modules.in_silico import reconcile_in_silico_predictions
        in_silico_findings = reconcile_in_silico_predictions(vus, in_silico_predictions, context, model_registry)

        discordances = list(collect_discordances(results))
        if len(engine_reports) > 1:
            discordances.extend(concordance(engine_reports).disagreed)

        # -- QC ------------------------------------------------------------
        warnings = list(adapted.warnings)
        for report in engine_reports:
            warnings.extend(report.warnings)
        for v in outcome.unexamined:
            warnings.append(
                f"{v.label()}: no lane could interpret this variant "
                f"({v.consequence.value} in {v.gene}); reported, not dropped."
            )

        from ..reporting.json_report import primitive
        source_hashes = {"variants": hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
                         "catalogue": self.catalogue.sha256}
        implementation = hashlib.sha256()
        package_root = Path(__file__).resolve().parents[1]
        for source in sorted(package_root.rglob("*.py")):
            implementation.update(source.relative_to(package_root).as_posix().encode())
            implementation.update(source.read_bytes().replace(b"\r\n", b"\n"))
        source_hashes["implementation"] = implementation.hexdigest()
        if self.local_validation_store is not None:
            source_hashes["historical_associations"] = hashlib.sha256(json.dumps(
                self.local_validation_store.to_json(), sort_keys=True).encode()).hexdigest()
        if mask_path:
            source_hashes["mask"] = hashlib.sha256(Path(mask_path).read_bytes()).hexdigest()
        for er in engine_reports:
            source_hashes[er.engine.name] = er.source_sha256 or hashlib.sha256(
                json.dumps(primitive(er), sort_keys=True).encode()).hexdigest()
        manifest = {"sources": source_hashes, "context": context.to_dict(),
                    "policy": asdict(self.policy), "profile": asdict(self.profile),
                    "phenotypes": [asdict(p) for p in phenotypes],
                    "coverage": {k: asdict(v) for k, v in mask._loci.items()},
                    "depth_floor": self.depth_floor, "callable_floor": self.callable_fraction_floor,
                    "tool_version": __version__, "demo_mode": self.demo_mode,
                    "include_tier2_loci": self.reconciler.include_tier2_loci,
                    "test_menu": test_menu, "follow_up_budget": follow_up_budget}
        manifest["population_options"] = {"enabled": include_population_structure,
                                          "tolerance": population_tolerance,
                                          "lineage_hints": list(population_lineage_hints)}
        manifest["epistasis_table"] = epistasis_table.to_dict() if epistasis_table else None
        manifest["mic_predictions"] = [asdict(p) for p in mic_predictions]
        manifest["structural_annotations"] = structural_annotations.to_dict() if structural_annotations is not None else []
        manifest["regulatory_catalogue"] = regulatory_catalogue.to_dict() if regulatory_catalogue is not None else None
        manifest["expression_evidence"] = [asdict(e) for e in expression_evidence]
        manifest["in_silico_predictions"] = [asdict(p) for p in in_silico_predictions]
        manifest["model_registry"] = model_registry.to_dict() if model_registry is not None else None
        manifest["watchlist_enabled"] = watchlist is not None
        fingerprint = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
        provenance = Provenance(
            input_sha256=source_hashes["variants"], analysis_fingerprint=fingerprint,
            policy_version=self.policy.version, source_hashes=source_hashes,
            tool="Myconductor",
            version=__version__,
            organism_profile=self.profile.name,
            profile_version=self.profile.version,
            reference_assembly=adapted.assembly,
            depth_floor=self.depth_floor,
            callable_fraction_floor=self.callable_fraction_floor,
            catalogue=self.catalogue.engine,
            engines=[r.engine for r in engine_reports],
            coverage_source=mask.source,
            demo_mode=self.demo_mode,
            generated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        )

        report = AnalysisReport(
            expression_findings=expression_findings,
            in_silico_findings=in_silico_findings,
            regulatory_findings=regulatory_findings,
            structural_annotations=structural_records,
            mic_predictions=[asdict(p) for p in mic_predictions],
            quantitative_findings=quantitative_findings,
            analysis_manifest=manifest,
            sample_id=adapted.sample_id,
            context=context.to_dict(),
            drug_results=results,
            provenance=provenance,
            eligibility=eligibility,
            heteroresistance=het,
            vus_priorities=vus,
            mechanism_queue=mechanism_queue,
            discordances=discordances,
            qc=qc,
            qc_warnings=warnings,
            lane_counts=outcome.lane_counts(),
            demo_mode=self.demo_mode,
        )

        if include_population_structure:
            from ..modules.population_structure import population_structure
            report.population_structure = population_structure(
                report.sample_id, list(variants_by_key.values()), het,
                population_tolerance, population_lineage_hints)
        if epistasis_table:
            report.epistasis_notes = epistasis_table.annotate(
                report.drug_results, list(variants_by_key.values()),
                self.profile.name, self.profile.reference_assembly)

        from ..modules.investigation import investigate, plan_follow_up
        report.investigations = investigate(report)
        report.follow_up = plan_follow_up(report.investigations, test_menu, follow_up_budget)
        if watchlist is not None:
            report.discordance_tickets = watchlist.capture(report)
        return report


def _species_matches(species: str, organism: str) -> bool:
    name = species.lower().replace("_", " ").strip()
    if organism == "mtbc":
        return name in {"mtbc", "mycobacterium tuberculosis", "mycobacterium tuberculosis complex",
                        "mycobacterium africanum", "mycobacterium bovis"}
    if organism == "mabscessus":
        return name.startswith(("mycobacterium abscessus", "mycobacteroides abscessus"))
    return name == organism.lower()
