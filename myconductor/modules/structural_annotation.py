"""Import attributed structural annotations without predicting function or ranking resistance."""
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from ..core.models import Dimension, MechanismHypothesis


@dataclass(frozen=True)
class StructuralAnnotation:
    variant_key: str
    location: str
    predicted_effect: str
    source: str
    organism: str
    reference_assembly: str
    source_version: str
    timestamp: str
    confidence: Optional[float] = None
    ligand_distance: Optional[float] = None
    distance_unit: Optional[str] = None
    ligand_reference: Optional[str] = None
    validation_status: str = "unknown"

    def __post_init__(self):
        for name in ("variant_key", "predicted_effect", "source", "organism", "reference_assembly", "source_version", "timestamp"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"structural annotation requires {name}")
        if self.location not in {"active_site", "binding_pocket", "distal", "unknown"}:
            raise ValueError("unknown structural location")
        if self.confidence is not None and (not math.isfinite(self.confidence) or not 0 <= self.confidence <= 1):
            raise ValueError("annotation confidence must be in [0,1]")
        if self.ligand_distance is not None:
            if not math.isfinite(self.ligand_distance) or self.ligand_distance < 0:
                raise ValueError("ligand distance must be finite and nonnegative")
            if self.distance_unit not in {"angstrom", "nm"} or not self.ligand_reference:
                raise ValueError("ligand distance requires an explicit unit and ligand reference")
        if self.validation_status not in {"unknown", "pending", "reviewed"}:
            raise ValueError("unknown structural validation status")

    @property
    def annotation_id(self):
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()


class StructuralAnnotationTable:
    def __init__(self, annotations=()):
        self.annotations = tuple(annotations)
        keys = [(a.variant_key, a.organism, a.reference_assembly, a.source, a.source_version) for a in self.annotations]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate structural observation; use separate source versions")

    @classmethod
    def load(cls, path):
        path = Path(path)
        if path.suffix.lower() == ".tsv":
            with path.open(encoding="utf-8-sig", newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            for row in rows:
                for name in ("confidence", "ligand_distance"):
                    if name in row:
                        row[name] = float(row[name]) if row[name] else None
        else:
            rows = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(rows, list):
            raise ValueError("structural annotations must be a JSON list or TSV")
        return cls(StructuralAnnotation(**r) for r in rows)

    def to_dict(self):
        return [asdict(a) for a in self.annotations]

    def attach(self, priorities, organism, reference_assembly):
        """Attach dimensions after ranking; imported claims cannot change scores."""
        used, hypotheses = {}, []
        for priority in priorities:
            annotations = [a for a in self.annotations if a.variant_key == priority.variant_key
                           and a.organism == organism and a.reference_assembly == reference_assembly]
            if not annotations:
                continue
            for a in annotations:
                used[a.annotation_id] = dict(asdict(a), annotation_id=a.annotation_id, gene=priority.gene,
                                            status="external_claim", ranking_effect="none")
                hypotheses.append(MechanismHypothesis(
                    priority.variant_label, priority.variant_key, priority.gene, priority.drug,
                    hypothesis="External structural annotation: " + a.predicted_effect,
                    evidence_gaps=[f"Source: {a.source} ({a.source_version}); structural annotation does not establish drug response."],
                    resolving_experiments=["Review the source method and independently measured phenotype."]))
            # Preserve distinct/contradictory source assertions rather than averaging confidence.
            replacements = {
                "structural_impact": Dimension("structural_impact", True,
                    [{"location": a.location, "reported_effect": a.predicted_effect, "source_confidence": a.confidence}
                     for a in annotations], source="; ".join(a.source for a in annotations),
                    note="Imported claims, not a functional score; excluded from ranking."),
            }
            distances = [a for a in annotations if a.ligand_distance is not None]
            if distances:
                replacements["ligand_distance"] = Dimension("ligand_distance", True,
                    [{"value": a.ligand_distance, "unit": a.distance_unit, "ligand": a.ligand_reference} for a in distances],
                    source="; ".join(a.source for a in distances), note="Reported distance only; no inferred binding effect.")
            priority.dimensions = [replacements.get(d.name, d) for d in priority.dimensions]
            priority.data_gaps = [g for g in priority.data_gaps if g not in replacements]
        return [used[k] for k in sorted(used)], hypotheses
