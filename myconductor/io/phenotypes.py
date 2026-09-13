"""Measured isolate phenotypes, independent of variant-association evidence."""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from ..core.context import SampleContext
from ..core.models import Call, DrugEvidence, Lane, Tier


@dataclass(frozen=True)
class PhenotypeObservation:
    observation_id: str
    sample_id: str
    isolate_id: str
    site_id: str
    organism: str
    drug: str
    result: str
    method: str
    measured_at: str
    laboratory: str
    quality: str
    interpretation_standard: str
    standard_version: str
    source: str
    mic: Optional[float] = None
    mic_unit: Optional[str] = None
    censoring: str = "none"
    critical_concentration: Optional[float] = None
    incubation_days: Optional[float] = None
    replicate_id: Optional[str] = None
    mic_interval: Optional[tuple[float, float]] = None
    susceptible_inclusive: Optional[bool] = None

    def __post_init__(self) -> None:
        for key in ("observation_id", "sample_id", "isolate_id", "site_id", "organism",
                    "drug", "method", "measured_at", "laboratory", "quality",
                    "interpretation_standard", "standard_version", "source"):
            if not isinstance(getattr(self, key), str) or not getattr(self, key).strip():
                raise ValueError(f"phenotype requires {key}")
        if self.result not in ("resistant", "susceptible", "indeterminate"):
            raise ValueError("phenotype result must be resistant, susceptible or indeterminate")
        if self.quality not in ("pass", "fail", "unknown"):
            raise ValueError("phenotype quality must be pass, fail or unknown")
        if self.censoring not in ("none", "left", "right", "interval"):
            raise ValueError("unknown MIC censoring")
        for key in ("mic", "critical_concentration", "incubation_days"):
            value = getattr(self, key)
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError(f"{key} must be finite and positive")
        if self.mic is not None and not self.mic_unit:
            raise ValueError("MIC requires an explicit unit")
        if self.critical_concentration is not None and not self.mic_unit:
            raise ValueError("critical concentration requires an explicit unit")
        if self.censoring == "interval" and self.mic_interval is None:
            raise ValueError("interval-censored MIC requires mic_interval")

    def evidence(self, context: SampleContext) -> DrugEvidence:
        for key in ("sample_id", "isolate_id", "site_id", "organism"):
            if getattr(self, key) != getattr(context, key):
                raise ValueError(f"phenotype {self.observation_id}: {key} does not match context")
        call = Call(self.result)
        limitations = []
        mic_comparison = None
        if self.mic is not None and self.critical_concentration is not None:
            from ..modules.mic_evidence import mic_call
            mic_comparison = mic_call(self.mic, self.mic_unit, self.critical_concentration,
                                      self.censoring, self.mic_interval, self.susceptible_inclusive)
            if mic_comparison is not None and mic_comparison != call:
                limitations.append("MIC versus reported category mismatch; review method and breakpoint convention. Laboratory category retained.")
            elif mic_comparison is None:
                limitations.append("MIC bounds or boundary convention do not establish a categorical comparison.")
        if self.quality != "pass":
            call = Call.INDETERMINATE
            limitations.append("phenotype quality has not passed laboratory review")
        if (context.organism == "mabscessus" and self.drug == "clarithromycin"
                and call is Call.SUSCEPTIBLE):
            # Early susceptibility cannot exclude induction unless erm(41)
            # has been independently shown nonfunctional.
            if context.erm41_status == "functional":
                call = Call.INDETERMINATE
                limitations.append("susceptible phenotype conflicts with functional erm(41)")
            elif (context.erm41_status != "nonfunctional"
                  and (self.incubation_days is None or self.incubation_days < 14)):
                call = Call.INDETERMINATE
                limitations.append("inducible macrolide resistance has not been assessed")
        return DrugEvidence(
            drug=self.drug, call=call, tier=Tier.PHENOTYPIC, lane=Lane.PHENOTYPE,
            sample_id=self.sample_id, observation_id=self.observation_id,
            rationale=f"Laboratory reports {self.result}: {self.method}, {self.measured_at}.",
            limitations=tuple(limitations), metadata=dict(asdict(self),
                mic_comparison=mic_comparison.value if mic_comparison else None),
        )


def load_phenotypes(path: str | Path) -> list[PhenotypeObservation]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("phenotypes must be a JSON list of observations")
    observations = [PhenotypeObservation(**r) for r in data]
    ids = [r.observation_id for r in observations]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate phenotype observation_id")
    return observations
