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
from pathlib import Path
from typing import Iterable, Optional, Sequence

from .. import __version__
from ..adapters.base import EngineReport, concordance, merge as merge_evidence
from ..catalogue.profile import OrganismProfile, load_profile
from ..core.models import AnalysisReport, EngineRef, Provenance
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
    ):
        self.profile = profile or load_profile()
        self.depth_floor = depth_floor
        self.callable_fraction_floor = callable_fraction_floor
        self.platform = platform
        self.demo_mode = demo_mode
        self.local_validation_store = local_validation_store

        catalogue = CatalogueModule()
        efflux = EffluxRegulatoryModule(self.profile)
        known = set(catalogue.labels)
        if local_validation_store is not None:
            # A variant this site has already validated has an answer; it
            # should stop being re-ranked as "needs laboratory validation".
            known |= local_validation_store.validated_labels()
        self.workbench = VUSWorkbench(annotator=annotator, known=known)
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
    ) -> AnalysisReport:
        adapted = load(
            input_path, depth_floor=self.depth_floor, sample=sample,
            platform=platform or self.platform,
        )

        if mask is None and mask_path is not None:
            mask = self._mask_from_path(mask_path)
        mask = mask or CallableMask.absent()

        engine_reports = list(engine_reports)
        engine_masks = [r.mask for r in engine_reports if r.mask is not None]
        if engine_masks:
            mask = CallableMask.merge(mask, *engine_masks)

        # -- evidence ------------------------------------------------------
        outcome = self.router.route(adapted.variants)
        evidence = list(outcome.evidence) + merge_evidence(engine_reports)

        all_variants = list(adapted.variants)
        for report in engine_reports:
            all_variants.extend(report.variants)

        lineage = next((r.lineage for r in engine_reports if r.lineage), None)
        het = self.het.assess(all_variants, evidence, lineage=lineage)
        if self.het.calibration is not None and self.het.prior_source is not None:
            evidence = evidence + heteroresistance_evidence(het, all_variants)
        results = self.reconciler.reconcile(evidence, mask)
        eligibility = self.assessor.assess(results)
        vus = self.workbench.priorities(adapted.variants)
        mechanism_queue = self.router.efflux.mechanism_queue(adapted.variants)

        discordances = list(collect_discordances(results))
        if len(engine_reports) > 1:
            discordances.extend(concordance(engine_reports).disagreed)

        # -- QC ------------------------------------------------------------
        qc = qc_module.assess(adapted, mask, self.profile, self.depth_floor)
        for report in engine_reports:
            qc.extend(report.qc)
        warnings = list(adapted.warnings)
        for report in engine_reports:
            warnings.extend(report.warnings)
        for v in outcome.unexamined:
            warnings.append(
                f"{v.label()}: no lane could interpret this variant "
                f"({v.consequence.value} in {v.gene}); reported, not dropped."
            )

        provenance = Provenance(
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

        return AnalysisReport(
            sample_id=adapted.sample_id,
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
