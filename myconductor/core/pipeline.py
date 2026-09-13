"""The conductor: wires every stage into one end-to-end analysis.

    adapter -> router -> heteroresistance -> synthesis -> (discovery) -> report

This is the object a caller instantiates. Each collaborator is injectable, so a
deployment can swap in a real catalogue, a trained VUS model, or a live docking
backend without touching this orchestration code.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from .. import __version__
from ..io.adapter import load
from ..modules import discovery, heteroresistance
from ..modules.synthesis import RegimenSynthesiser, reconcile
from .models import AnalysisReport, Provenance
from .router import TriageRouter


class Myconductor:
    def __init__(
        self,
        router: Optional[TriageRouter] = None,
        synthesiser: Optional[RegimenSynthesiser] = None,
        depth_floor: int = 10,
        het_vaf_floor: float = 0.03,
    ):
        self.router = router or TriageRouter()
        self.synthesiser = synthesiser or RegimenSynthesiser()
        self.depth_floor = depth_floor
        self.het_vaf_floor = het_vaf_floor

    def analyze(self, input_path: str | Path) -> AnalysisReport:
        adapted = load(input_path, depth_floor=self.depth_floor)

        outcome = self.router.route(adapted.variants)
        evidence = outcome.evidence

        het = heteroresistance.detect(
            adapted.variants, evidence, noise_floor=self.het_vaf_floor
        )

        drug_results = reconcile(evidence)
        regimen = self.synthesiser.propose(drug_results)

        disc = discovery.run(
            reason=regimen.rationale, triggered=not regimen.adequate
        )

        provenance = Provenance(
            tool="Myconductor",
            version=__version__,
            catalogue_version=self.router.catalogue.version,
            vus_model=self.router.vus.model_name,
            coverage_min=self.depth_floor,
            het_vaf_floor=self.het_vaf_floor,
        )

        return AnalysisReport(
            sample_id=adapted.sample_id,
            drug_results=drug_results,
            heteroresistance=het,
            regimen=regimen,
            discovery=disc,
            provenance=provenance,
            qc_warnings=adapted.warnings,
            routed_counts=outcome.counts(),
        )
