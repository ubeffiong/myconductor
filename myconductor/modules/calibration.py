"""Calibrated, per-drug/per-lineage uncertainty for minority-allele calls.

Why this exists
----------------
``heteroresistance.py`` compares a measured allele fraction against a fixed
per-platform floor (``PLATFORM_LOD``) and reports a binary
assessable/detected. That is deliberately conservative, and the module says so:
those floors are "conservative defaults for triage, NOT validated limits of
detection." But a floor throws away the one thing a clinician actually wants
at a borderline fraction: how likely is this to be real, and actionable,
*here* — for this drug, this platform, this lineage?

What this module adds
----------------------
Two real, measurable inputs, both required, neither fabricated:

1. A **detection-probability curve** fit to a dilution series
   (``DetectionCalibration``): sensitivity at a given VAF, and a background
   false-positive rate from wild-type-only negative controls at the same VAF
   range. This is standard analytical validation, not a model invented here.
2. A **prior** — the real-world prevalence of a true resistance-conferring
   minority subpopulation at this exact variant/drug (``PriorSource``). The
   only non-fabricated source of this in the codebase today is isolate-level
   federated evidence (``FederatedPriorSource``, built from
   ``federated.catalogue_update.VariantDrugAssociation``); the honest default
   (``NullPriorSource``) reports no prior at all.

Both are required to compute anything. Missing either yields ``None`` — never
a guessed number — from ``posterior_resistance_probability``. This ships with
**no bundled calibration data**, for the same reason the WHO catalogue is not
redistributed here: publishing a plausible-looking table with no dilution
series behind it would be worse than shipping nothing.

What this does NOT change
--------------------------
The posterior this module computes can only ever surface as
``Tier.PREDICTED`` / ``Call.INDETERMINATE`` evidence (see
``heteroresistance.heteroresistance_evidence``). A statistical posterior is
not graded catalogue evidence or a laboratory phenotype, so
``DrugEvidence.__post_init__`` structurally forbids it from ever asserting
``RESISTANT`` — however high the number.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Protocol

from ..core.models import Variant

_WILDCARD_LINEAGE = "*"


@dataclass(frozen=True)
class DetectionCalibration:
    """A measured detection curve for one platform/drug/lineage combination.

    Every field must come from a dilution-series study with wild-type-only
    negative controls — never a guess:

    ``lod50``              alternate-allele fraction detected 50% of the time.
    ``width``               logistic scale over which detection rises from
                            unreliable to reliable (roughly the VAF span
                            covering 12%-88% detection probability).
    ``background_fp_rate``  fraction of wild-type-only replicates that show an
                            apparent minority signal in this VAF range —
                            the assay's limit of blank, not zero by assumption.
    ``n_replicates``        how many dilution-series replicates this rests on;
                            used only to size the credible interval, never to
                            silently upgrade a thin study to a confident one.
    """

    platform: str
    drug: str
    lineage: str
    lod50: float
    width: float
    background_fp_rate: float
    n_replicates: int
    source: str

    def p_detect(self, vaf: float) -> float:
        """Probability a truly-present variant at this VAF is detected."""
        scale = max(self.width, 1e-6) / 4.394449  # logistic 12%..88% span == width
        z = (vaf - self.lod50) / scale
        return 1.0 / (1.0 + math.exp(-z))


@dataclass
class CalibrationTable:
    """Loaded, deployment-specific detection calibrations.

    Ships with no entries. A deployment adds its own dilution-series results
    via :meth:`from_json` or :meth:`add`. Lookup falls back from a specific
    lineage to the wildcard lineage ``"*"`` (a platform/drug calibration that
    does not yet vary by lineage), and to ``None`` — never to a fabricated
    default — when nothing is configured for a platform/drug at all. (The flat
    ``PLATFORM_LOD`` floors in ``heteroresistance.py`` remain the fallback
    triage floor when no calibration is configured; this table is additive,
    not a replacement.)
    """

    _entries: dict[tuple[str, str, str], DetectionCalibration] = field(
        default_factory=dict)

    @classmethod
    def from_json(cls, path: str | Path) -> "CalibrationTable":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        table = cls()
        for row in data.get("calibrations", []):
            table.add(DetectionCalibration(
                platform=row["platform"].strip().lower(),
                drug=row["drug"].strip().lower(),
                lineage=(row.get("lineage") or _WILDCARD_LINEAGE).strip(),
                lod50=float(row["lod50"]),
                width=float(row["width"]),
                background_fp_rate=float(row.get("background_fp_rate", 0.0)),
                n_replicates=int(row.get("n_replicates", 0)),
                source=row.get("source", "unspecified"),
            ))
        return table

    def add(self, calibration: DetectionCalibration) -> None:
        key = (calibration.platform, calibration.drug, calibration.lineage)
        self._entries[key] = calibration

    def get(self, platform: Optional[str], drug: Optional[str],
            lineage: Optional[str]) -> Optional[DetectionCalibration]:
        if not platform or not drug:
            return None
        platform, drug = platform.strip().lower(), drug.strip().lower()
        if lineage:
            hit = self._entries.get((platform, drug, lineage.strip()))
            if hit is not None:
                return hit
        return self._entries.get((platform, drug, _WILDCARD_LINEAGE))

    def __len__(self) -> int:
        return len(self._entries)


class PriorSource(Protocol):
    def prior_for(self, variant_key: str, drug: str,
                  lineage: Optional[str]) -> Optional[float]:
        ...


class NullPriorSource:
    """The honest default: no prevalence data, so no prior.

    Matches ``modules.features.NullAnnotator``'s posture — reporting
    "unavailable" is a legitimate output, and is preferable to inventing a
    number.
    """

    def prior_for(self, variant_key: str, drug: str,
                  lineage: Optional[str]) -> Optional[float]:
        return None


@dataclass
class FederatedPriorSource:
    """Derives a prior from real, isolate-level federated evidence.

    Built from ``federated.catalogue_update.VariantDrugAssociation`` objects a
    deployment already has from :func:`federated.catalogue_update.aggregate`
    — the one non-fabricated source of prevalence data this codebase has for
    "how often does this exact variant, at this drug, turn out to be
    phenotypically resistant." Below ``min_isolates`` the association is not
    trusted as a prior and ``None`` is returned rather than a number computed
    from too few observations.
    """

    associations: dict  # (variant_key, drug) -> VariantDrugAssociation
    min_isolates: int = 10

    def prior_for(self, variant_key: str, drug: str,
                  lineage: Optional[str]) -> Optional[float]:
        assoc = self.associations.get((variant_key, drug))
        if assoc is None or assoc.n_isolates < self.min_isolates:
            return None
        return assoc.resistant_fraction


@dataclass(frozen=True)
class PosteriorEstimate:
    probability: float
    credible_interval: tuple[float, float]
    basis: str


def posterior_resistance_probability(
    variant: Variant,
    drug: str,
    calibration: CalibrationTable,
    prior_source: PriorSource,
    lineage: Optional[str] = None,
) -> Optional[PosteriorEstimate]:
    """Bayesian posterior that a detected minority allele is real resistance.

    An explicit, auditable measurement model, not a fixed rule:

        P(true | detected) = sens * prior / (sens * prior + fp_rate * (1 - prior))

    where ``sens`` is the calibration curve's detection probability at this
    VAF and ``fp_rate`` is its measured background false-positive rate. Both
    inputs — the calibration and the prior — must be configured for this
    platform/drug/lineage; this returns ``None``, never a guessed number,
    whenever either is missing. Callers must not treat ``None`` as "low
    probability" — it means "not measured here."

    This never becomes a resistance call by itself. See
    ``heteroresistance.heteroresistance_evidence``, which can only emit
    ``Tier.PREDICTED`` / ``Call.INDETERMINATE`` evidence from this number —
    structurally barred from ``RESISTANT`` regardless of how high the
    posterior is.
    """
    if variant.vaf is None or variant.platform is None:
        return None
    cal = calibration.get(variant.platform, drug, lineage)
    if cal is None:
        return None
    prior = prior_source.prior_for(variant.key(), drug, lineage)
    if prior is None:
        return None
    if not 0.0 <= prior <= 1.0:
        return None

    sensitivity = cal.p_detect(variant.vaf)
    numerator = sensitivity * prior
    denominator = numerator + cal.background_fp_rate * (1.0 - prior)
    if denominator <= 0.0:
        return None
    posterior = min(1.0, max(0.0, numerator / denominator))

    # An honestly-bounded interval, not a substitute for a hierarchical model:
    # it widens as the dilution series thins, and stays wide by default when
    # replicate count is small, rather than reporting a false point estimate.
    spread = 0.5 if cal.n_replicates < 5 else max(0.03, 1.0 / math.sqrt(cal.n_replicates))
    lower, upper = max(0.0, posterior - spread), min(1.0, posterior + spread)

    return PosteriorEstimate(
        probability=round(posterior, 4),
        credible_interval=(round(lower, 4), round(upper, 4)),
        basis=(f"calibration: {cal.source} (n={cal.n_replicates} dilution "
               f"replicate(s), background FP rate {cal.background_fp_rate:.1%}); "
               f"prior: isolate-level federated evidence"),
    )
