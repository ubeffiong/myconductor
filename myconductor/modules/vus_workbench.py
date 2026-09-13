"""VUS workbench — ranks variants for validation, never calls resistance.

Why this replaced the classifier
--------------------------------
The previous ``vus_classifier`` emitted ``PREDICTED_RESISTANT`` /
``PREDICTED_SUSCEPTIBLE`` with a probability and a SHAP-style explanation. Its
features were SHA-256 digests of the variant name, so the probability carried no
information and the explanation described the structure of a hash. Explaining an
unvalidated model does not make it reliable; it makes it persuasive, which is
worse.

The defensible product at this stage is **prioritisation**. For each variant of
unknown significance, report what is known across each evidence dimension, say
plainly which dimensions have no data source, and name the experiment that would
resolve it. A laboratory can act on that. Nobody can act on a hash.

What this lane may and may not conclude
--------------------------------------
It can never make a drug ``RESISTANT``; ``Tier`` forbids it structurally. Nor
may it make a drug ``SUSCEPTIBLE``. What it does emit is ``INDETERMINATE`` for
the affected drug, which says the honest thing: an uninterpreted variant was
found in a locus that matters for this drug, so susceptibility cannot be
concluded even if the locus was well covered.

That middle position is the whole point. Letting a VUS pass silently would
allow a drug whose target carries an uncharacterised change to be reported
susceptible on coverage alone. Calling it resistant would invent evidence.
Promotion to a genuine predictive call happens only after independent
validation and per-drug calibration — Phase 4's exit gate, not a code path.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from ..core.models import (
    Call,
    Consequence,
    DrugEvidence,
    EngineRef,
    Lane,
    Tier,
    Variant,
    VUSPriority,
)
from .base import VariantModule
from .features import AnnotatorProtocol, FeatureSet, NullAnnotator

if TYPE_CHECKING:
    from ..catalogue.profile import OrganismProfile

#: Gene -> the drug whose resistance a variant there would most plausibly
#: affect. Context for the reader and for experiment design; NOT a prediction.
GENE_DRUG_CONTEXT = {
    "katG": "isoniazid", "inhA": "isoniazid", "fabG1": "isoniazid",
    "ahpC": "isoniazid", "ndh": "isoniazid",
    "rpoB": "rifampicin", "rpoA": "rifampicin", "rpoC": "rifampicin",
    "pncA": "pyrazinamide", "panD": "pyrazinamide", "rpsA": "pyrazinamide",
    "embB": "ethambutol", "embA": "ethambutol", "embC": "ethambutol",
    "gyrA": "moxifloxacin", "gyrB": "moxifloxacin",
    "rrs": "amikacin", "eis": "amikacin", "tlyA": "amikacin",
    "atpE": "bedaquiline", "Rv0678": "bedaquiline", "pepQ": "bedaquiline",
    "ddn": "pretomanid", "fgd1": "pretomanid", "fbiA": "pretomanid",
    "fbiB": "pretomanid", "fbiC": "pretomanid",
    "rplC": "linezolid", "rrl": "linezolid",
}

#: How many dimensions must carry real data before a numeric rank is meaningful.
MIN_DIMENSIONS_FOR_SCORE = 5

_EXPERIMENTS = {
    Consequence.NONSENSE: (
        "Allelic exchange or CRISPRi knockdown of {gene}, then MIC against "
        "{drug} — a truncating change is the most tractable functional test."
    ),
    Consequence.FRAMESHIFT: (
        "Allelic exchange or CRISPRi knockdown of {gene}, then MIC against "
        "{drug} — a truncating change is the most tractable functional test."
    ),
    Consequence.DELETION: (
        "Complementation of {gene} in a deletion background, then MIC against "
        "{drug}."
    ),
    Consequence.MISSENSE: (
        "Site-directed mutagenesis of {gene} {change} in a susceptible "
        "background, then MIC against {drug}; pair with a wild-type control."
    ),
    Consequence.UPSTREAM: (
        "Reporter fusion or RT-qPCR to measure {gene} expression change, then "
        "MIC against {drug} — a regulatory variant needs an expression "
        "readout, not a structural one."
    ),
}
_DEFAULT_EXPERIMENT = (
    "Resolve the consequence first (annotate against the reference), then "
    "site-directed mutagenesis of {gene} {change} and MIC against {drug}."
)


class VUSWorkbench(VariantModule):
    """Ranks unknown variants for functional validation.

    Implements ``VariantModule`` so the router can run it as a lane, but
    ``evaluate`` always returns an empty evidence list by design. The rankings
    are collected separately via ``priorities``.
    """

    name = "vus_workbench"
    lane = Lane.VUS

    def __init__(self, annotator: Optional[AnnotatorProtocol] = None,
                 known: Optional[set[str]] = None,
                 min_dimensions: int = MIN_DIMENSIONS_FOR_SCORE,
                 profile: Optional["OrganismProfile"] = None):
        self.annotator: AnnotatorProtocol = annotator or NullAnnotator()
        self.known = known or set()
        self.min_dimensions = min_dimensions
        self.model_name = type(self.annotator).__name__
        # GENE_DRUG_CONTEXT is written against the bundled MTBC gene set. A
        # profile scopes it to genes that organism's own drug_loci actually
        # declares, so a coincidentally-named gene from another organism's
        # profile (e.g. M. abscessus's own rpoB) is never silently
        # attributed to an MTBC drug it has no validated relationship to.
        # Without a profile (e.g. constructing this lane standalone), every
        # hardcoded gene applies, exactly as before this parameter existed —
        # every gene GENE_DRUG_CONTEXT names is already declared in the
        # bundled MTBC profile's own loci, so this is a no-op for the
        # default pipeline.
        self.gene_drugs = {}
        if profile is not None:
            for drug, tiers in profile.drug_loci.items():
                for genes in tiers.values():
                    for gene in genes:
                        self.gene_drugs.setdefault(gene, set()).add(drug)
        else:
            self.gene_drugs = {g: {d} for g, d in GENE_DRUG_CONTEXT.items()}
        self.applicable_genes = set(self.gene_drugs)

    @property
    def synthetic(self) -> bool:
        return bool(getattr(self.annotator, "synthetic", False))

    def applies_to(self, variant: Variant) -> bool:
        """Uncatalogued variants in a resistance-associated gene."""
        if variant.silent:
            return False
        if variant.label() in self.known:
            return False
        return variant.gene in self.applicable_genes

    def evaluate(self, variant: Variant) -> list[DrugEvidence]:
        """Withhold susceptibility for the affected drug. Never call resistance.

        An uninterpreted variant in a resistance-associated gene means the drug
        cannot be concluded susceptible, however well covered the locus was.
        That is ``INDETERMINATE`` — not a prediction, and not silence.
        """
        if not self.applies_to(variant):
            return []
        drugs = sorted(self.gene_drugs.get(variant.gene, ()))
        if not drugs:
            return []
        return [DrugEvidence(
            drug=drug,
            call=Call.INDETERMINATE,
            tier=Tier.NONE,
            lane=Lane.VUS,
            confidence=None,
            variant=variant.identity,
            engine=EngineRef(name="myconductor-vus-workbench", version="0.2.0",
                             database=self.model_name),
            limitations=(
                "variant of unknown significance: no graded catalogue entry",
                "no predictive call is made; see the VUS priority list for the "
                "evidence dimensions and the recommended validation experiment",
            ),
            rationale=(
                f"Uninterpreted {variant.consequence.value} variant "
                f"{variant.label()} in {variant.gene}, a locus relevant to "
                f"{drug}. Susceptibility cannot be concluded; resistance is not "
                f"asserted."
            ),
        ) for drug in drugs]

    # -- the actual output ------------------------------------------------
    def priorities(self, variants: list[Variant]) -> list[VUSPriority]:
        out = [self._rank(v, drug) for v in variants
               if self.applies_to(v) and v.label() not in self.known
               for drug in sorted(self.gene_drugs.get(v.gene, ()))]
        # Present the most tractable and most likely-functional first. With the
        # null annotator this is a triage ordering over consequence class, not
        # an evidence-weighted ranking -- `priority` says which it is.
        out.sort(key=lambda p: (
            p.score is None,
            -(p.score or 0.0),
            _consequence_order(p),
        ))
        return out

    def _rank(self, variant: Variant, drug: str) -> VUSPriority:
        features: FeatureSet = self.annotator.annotate(variant)
        available = features.available_names
        gaps = features.gaps

        score: Optional[float] = None
        priority = "insufficient-data"
        if len(available) >= self.min_dimensions and not features.synthetic:
            score = _score(features)
            priority = ("high" if score >= 0.66
                        else "moderate" if score >= 0.33 else "low")

        template = _EXPERIMENTS.get(variant.consequence, _DEFAULT_EXPERIMENT)
        experiment = template.format(
            gene=variant.gene, change=variant.change, drug=drug or "the relevant drug"
        )

        note = ""
        if features.synthetic:
            note = ("annotator is SYNTHETIC: values are hash-derived and carry "
                    "no biological meaning")
        elif score is None:
            note = (f"only {len(available)} of {len(features.dimensions)} "
                    f"dimensions have a data source ("
                    f"{self.min_dimensions} needed to rank); ordering below is "
                    f"a consequence-class triage, not an evidence ranking")

        dims = list(features.dimensions)
        if note:
            dims.append(
                type(dims[0])(name="_ranking_basis", available=True,
                              value=note, source=self.model_name)
                if dims else None
            )
        dims = [d for d in dims if d is not None]

        return VUSPriority(
            variant_label=variant.label(),
            variant_key=variant.key(),
            gene=variant.gene,
            drug=drug,
            priority=priority,
            dimensions=dims,
            recommended_experiment=experiment,
            data_gaps=gaps,
            score=score,
        )


def _score(features: FeatureSet) -> float:
    """Mean of the available numeric dimensions, on 0..1.

    Deliberately simple and unweighted: a weighted combination implies a
    calibration that does not exist yet. This is a rank, not a probability, and
    it is only computed when enough real dimensions are present.
    """
    values = []
    for d in features.dimensions:
        if not d.available:
            continue
        if isinstance(d.value, bool):
            values.append(1.0 if d.value else 0.0)
        elif isinstance(d.value, (int, float)):
            values.append(max(0.0, min(1.0, float(d.value))))
    return round(sum(values) / len(values), 3) if values else 0.0


_CONSEQUENCE_TRIAGE_ORDER = {
    Consequence.NONSENSE: 0,
    Consequence.FRAMESHIFT: 1,
    Consequence.DELETION: 2,
    Consequence.MISSENSE: 3,
    Consequence.INFRAME_INDEL: 4,
    Consequence.UPSTREAM: 5,
    Consequence.RRNA: 6,
    Consequence.UNKNOWN: 7,
    Consequence.SYNONYMOUS: 8,
}


def _consequence_order(priority: VUSPriority) -> int:
    for d in priority.dimensions:
        if d.name == "consequence" and d.available:
            try:
                return _CONSEQUENCE_TRIAGE_ORDER[Consequence(d.value)]
            except (ValueError, KeyError):
                return 7
    return 7
