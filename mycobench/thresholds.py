"""Pre-registered performance targets, fixed before any data is examined.

Why this file exists separately, and is hashed
----------------------------------------------
A threshold chosen after seeing the results is not a threshold. These targets
are declared here, in version control, with a content hash that the report
prints — so a reader can tell whether the bar was set before or after the
measurement, and a silent edit is visible.

What these numbers are, and are not
-----------------------------------
They are **this project's declared targets**, informed by the expectations set
out in WHO target product profiles for drug-susceptibility tests and by the
error classes used in clinical microbiology method comparison. They are **not
WHO-endorsed thresholds**, and clearing one is not regulatory acceptance.

One bound, two names
--------------------
The very major error rate is *identically* ``1 - sensitivity``: both are
``FN / (TP + FN)``. The major error rate is identically ``1 - specificity``.
An earlier version of this file declared ``max_vme`` and ``min_sensitivity``
as independent fields, and six of eleven drugs ended up with contradictory
values — so one bound silently overrode the other depending on which check ran.

They are now **derived**, not declared. Sensitivity and specificity are the
inputs; VME and ME are the same bounds under the names a clinician reads them
by:

``VME`` very major error — predicted susceptible, phenotypically resistant.
    A patient is given a drug that will not work. The most consequential error
    here, which is why sensitivity targets are set first and the VME ceiling
    follows from them.

``ME``  major error — predicted resistant, phenotypically susceptible.
    A usable drug is withheld. Harmful, but recoverable.

Abstention is not an error
--------------------------
Myconductor can decline to call a drug (``NOT_ASSESSED``, ``INDETERMINATE``).
Those are excluded from the error rates and reported separately as a call rate,
because a tool that abstains on every hard isolate would otherwise post perfect
accuracy. ``min_call_rate`` is what stops that: a drug must be *called* often
enough for its error rates to mean anything.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from typing import Optional

#: Bump when any target below changes, and say why in the changelog.
REGISTRATION_VERSION = "1.1.0"
REGISTRATION_DATE = "2026-09-13"


@dataclass(frozen=True)
class DrugTarget:
    """Pre-registered targets for one drug.

    Only sensitivity, specificity, call rate and the minimum evaluable count
    are declared. The VME and ME ceilings are derived, because they are the
    same two quantities.
    """

    drug: str
    min_sensitivity: float
    min_specificity: float
    min_call_rate: float
    min_evaluable: int
    rationale: str

    @property
    def max_vme(self) -> float:
        """Ceiling on the very major error rate, i.e. ``1 - sensitivity``."""
        return round(1.0 - self.min_sensitivity, 6)

    @property
    def max_me(self) -> float:
        """Ceiling on the major error rate, i.e. ``1 - specificity``."""
        return round(1.0 - self.min_specificity, 6)


#: Drugs are grouped by how much a decision depends on them and how good the
#: phenotypic reference standard is. Rifampicin and isoniazid carry the
#: strictest targets: they drive regimen choice and their DST is the most
#: reliable reference available. Bedaquiline, pretomanid, linezolid and
#: delamanid carry the loosest — a wide target there reflects the reference
#: standard's own uncertainty, not tolerance for error.
TARGETS: dict[str, DrugTarget] = {
    "rifampicin": DrugTarget(
        "rifampicin", min_sensitivity=0.95, min_specificity=0.98,
        min_call_rate=0.95, min_evaluable=30,
        rationale="drives MDR classification; rpoB DST is the most reliable "
                  "phenotypic reference available"),
    "isoniazid": DrugTarget(
        "isoniazid", min_sensitivity=0.95, min_specificity=0.98,
        min_call_rate=0.95, min_evaluable=30,
        rationale="companion to rifampicin for the MDR definition; katG and "
                  "inhA mechanisms are well characterised"),
    "moxifloxacin": DrugTarget(
        "moxifloxacin", min_sensitivity=0.90, min_specificity=0.95,
        min_call_rate=0.90, min_evaluable=30,
        rationale="fluoroquinolone eligibility gates BPaLM; gyrA and gyrB "
                  "mechanisms are known but MIC distributions overlap the "
                  "epidemiological cutoff"),
    "levofloxacin": DrugTarget(
        "levofloxacin", min_sensitivity=0.90, min_specificity=0.95,
        min_call_rate=0.90, min_evaluable=30,
        rationale="as moxifloxacin; reported separately because CRyPTIC "
                  "measures both and cross-resistance is not total"),
    "ethambutol": DrugTarget(
        "ethambutol", min_sensitivity=0.85, min_specificity=0.90,
        min_call_rate=0.90, min_evaluable=30,
        rationale="embB mechanisms are heterogeneous and its phenotypic DST is "
                  "poorly reproducible near the critical concentration"),
    "amikacin": DrugTarget(
        "amikacin", min_sensitivity=0.90, min_specificity=0.95,
        min_call_rate=0.85, min_evaluable=20,
        rationale="rrs a1401g dominates and is well characterised; the call "
                  "rate target is lower because rrs coverage often fails QC"),
    "ethionamide": DrugTarget(
        "ethionamide", min_sensitivity=0.80, min_specificity=0.88,
        min_call_rate=0.85, min_evaluable=20,
        rationale="shares inhA mechanisms with isoniazid; its phenotypic "
                  "reference is comparatively unreliable"),
    "clofazimine": DrugTarget(
        "clofazimine", min_sensitivity=0.80, min_specificity=0.90,
        min_call_rate=0.80, min_evaluable=20,
        rationale="Rv0678 efflux de-repression raises MICs without a clean "
                  "genotype-phenotype boundary; shares mechanisms with "
                  "bedaquiline"),
    "bedaquiline": DrugTarget(
        "bedaquiline", min_sensitivity=0.75, min_specificity=0.90,
        min_call_rate=0.75, min_evaluable=20,
        rationale="the weakest target here, deliberately: known mechanisms "
                  "explain only part of observed phenotypic resistance, and "
                  "the cutoff sits close to the wild-type MIC distribution"),
    "linezolid": DrugTarget(
        "linezolid", min_sensitivity=0.75, min_specificity=0.90,
        min_call_rate=0.75, min_evaluable=20,
        rationale="rplC and rrl mechanisms are rare in public data, so few "
                  "resistant isolates exist to measure sensitivity against"),
    "delamanid": DrugTarget(
        "delamanid", min_sensitivity=0.75, min_specificity=0.90,
        min_call_rate=0.75, min_evaluable=20,
        rationale="as linezolid; nitroimidazole activation mutations are "
                  "diverse and individually rare"),
    "pretomanid": DrugTarget(
        "pretomanid", min_sensitivity=0.75, min_specificity=0.90,
        min_call_rate=0.75, min_evaluable=20,
        rationale="shares the nitroimidazole activation pathway with "
                  "delamanid; CRyPTIC does not measure it, so this target "
                  "applies only where another phenotype source supplies it"),
}

#: Phenotype qualities accepted into the accuracy track. CRyPTIC grades each
#: phenotype HIGH / MEDIUM / LOW; measuring a predictor against a LOW-quality
#: phenotype measures the phenotype.
ACCEPTED_PHENOTYPE_QUALITY = ("HIGH",)

#: Drugs CRyPTIC measures that the mycobacterial profile does not cover.
#: Listed so a report can say "not compared" rather than omitting them.
UNCOMPARED_CRYPTIC_DRUGS = ("rifabutin", "kanamycin")


def registration_hash() -> str:
    """Content hash of every target, so a silent edit is visible in a report."""
    payload = json.dumps(
        {"version": REGISTRATION_VERSION, "date": REGISTRATION_DATE,
         "accepted_phenotype_quality": list(ACCEPTED_PHENOTYPE_QUALITY),
         "targets": {name: asdict(target)
                     for name, target in sorted(TARGETS.items())}},
        sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def target_for(drug: str) -> Optional[DrugTarget]:
    return TARGETS.get(drug.strip().lower())


def describe() -> str:
    return (f"pre-registration v{REGISTRATION_VERSION} ({REGISTRATION_DATE}), "
            f"sha256:{registration_hash()[:12]}, "
            f"{len(TARGETS)} drug target(s), "
            f"phenotype quality accepted: "
            f"{'/'.join(ACCEPTED_PHENOTYPE_QUALITY)}")
