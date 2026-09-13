"""Core domain model for Myconductor.

Four commitments are encoded in these types, because encoding them here is what
keeps every module downstream honest:

1. **Absence of evidence is not susceptibility.** ``Call`` carries four
   distinct "we did not establish susceptibility" states. Only a
   coverage-backed ``SUSCEPTIBLE`` ever counts toward regimen eligibility
   (``DrugResult.permits_use``).

2. **The verdict and its basis are separate axes.** ``Call`` says *what*;
   ``Tier`` says *how we know*. A graded catalogue hit and a model prediction
   can reach the same verdict on very different footings, and a reader must
   always be able to see which. Only ``CATALOGUED`` and ``PHENOTYPIC`` tiers
   are permitted to establish ``RESISTANT``; rule-based and model-based lanes
   report ``INDETERMINATE`` instead of asserting resistance.

3. **A variant's identity is its coordinates, not its label.** Two engines
   cannot be reconciled on ``"katG_S315T"`` — that string is a display label
   whose spelling varies by annotator. ``VariantIdentity`` carries assembly,
   coordinate and alleles so reconciliation is well defined.

4. **Disagreement is a finding, not noise.** ``DrugResult`` retains every
   piece of evidence it received and reports ``Discordance`` rather than
   collapsing to a single winner.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional
import math

#: H37Rv, the MTBC reference assembly the bundled profile is written against.
MTB_ASSEMBLY = "NC_000962.3"


class Call(str, Enum):
    """Per-drug verdict.

    The three-state (R / S / indeterminate) model that most tools expose cannot
    distinguish "we sequenced these loci and found no resistance mechanism"
    from "we never looked". Conflating those two is how a report ends up
    recommending a drug whose target was never covered.
    """

    RESISTANT = "resistant"
    SUSCEPTIBLE = "susceptible"
    INDETERMINATE = "indeterminate"
    NOT_ASSESSED = "not_assessed"
    NO_CALL = "no_call"
    UNSUPPORTED = "unsupported"

    @property
    def is_resistant(self) -> bool:
        return self is Call.RESISTANT

    @property
    def is_susceptible(self) -> bool:
        return self is Call.SUSCEPTIBLE

    @property
    def is_established(self) -> bool:
        """True only where the analysis actually reached a verdict."""
        return self in (Call.RESISTANT, Call.SUSCEPTIBLE)

    @property
    def reason_unestablished(self) -> Optional[str]:
        """Why no verdict was reached, for the states where that applies."""
        return {
            Call.INDETERMINATE: "conflicting evidence, or a lane that abstained",
            Call.NOT_ASSESSED: "required loci were not shown to be callable",
            Call.NO_CALL: "genotype at the locus is not reliable",
            Call.UNSUPPORTED: "drug is not covered by this organism profile",
        }.get(self)


class Tier(str, Enum):
    """The basis for a call — how we know what we claim to know."""

    PHENOTYPIC = "phenotypic"   # laboratory DST / MIC
    CATALOGUED = "catalogued"   # graded catalogue entry
    INFERRED = "inferred"       # rule-based mechanism, not graded
    PREDICTED = "predicted"     # statistical or ML model
    NONE = "none"               # no basis at all

    @property
    def rank(self) -> int:
        return {
            Tier.PHENOTYPIC: 4,
            Tier.CATALOGUED: 3,
            Tier.INFERRED: 2,
            Tier.PREDICTED: 1,
            Tier.NONE: 0,
        }[self]

    @property
    def may_establish_resistance(self) -> bool:
        """Rule-based and model-based lanes must not assert resistance.

        They can withhold susceptibility — which is the clinically protective
        action — without claiming an unvalidated mechanism is established.
        """
        return self in (Tier.PHENOTYPIC, Tier.CATALOGUED)

    @property
    def marker(self) -> str:
        """Single-character suffix used in rendered reports."""
        return {
            Tier.PHENOTYPIC: "",
            Tier.CATALOGUED: "",
            Tier.INFERRED: "†",   # dagger
            Tier.PREDICTED: "*",
            Tier.NONE: "",
        }[self]


class Lane(str, Enum):
    """Which analysis produced a piece of evidence.

    A variant may be examined by several lanes; see ``core.router``. This
    records provenance of evidence, not an exclusive assignment.
    """

    CATALOGUE = "catalogue"
    EFFLUX_REGULATORY = "efflux_regulatory"
    VUS = "vus"
    ENGINE = "engine"          # an external tool via adapters/
    PHENOTYPE = "phenotype"    # laboratory DST
    COVERAGE = "coverage"      # the callable-locus assessment itself
    NONE = "none"


class Region(str, Enum):
    CODING = "coding"
    PROMOTER = "promoter"
    INTERGENIC = "intergenic"
    RRNA = "rrna"


class Consequence(str, Enum):
    """Normalised variant consequence.

    Replaces the old free-floating ``silent`` boolean: synonymy is a property
    of the consequence, not an independent flag an input can contradict.
    """

    MISSENSE = "missense"
    NONSENSE = "nonsense"
    FRAMESHIFT = "frameshift"
    INFRAME_INDEL = "inframe_indel"
    SYNONYMOUS = "synonymous"
    UPSTREAM = "upstream"
    RRNA = "rrna"
    DELETION = "deletion"
    #: Insertional inactivation, including mobile-element insertion. Kept
    #: distinct from DELETION because the detection evidence differs, and
    #: added because IS-element insertion into mmpR5 is a leading route to
    #: bedaquiline resistance — a mechanism the model previously could not
    #: represent at all.
    INSERTION = "insertion"
    UNKNOWN = "unknown"

    @property
    def is_silent(self) -> bool:
        return self is Consequence.SYNONYMOUS

    @property
    def is_truncating(self) -> bool:
        """Does this consequence abolish the product's function?

        Load-bearing for efflux regulators: loss of function in mmpR5 is a
        materially stronger inference than a missense change of unknown
        effect, and is reported as such.
        """
        return self in (Consequence.NONSENSE, Consequence.FRAMESHIFT,
                        Consequence.DELETION, Consequence.INSERTION)


@dataclass(frozen=True)
class VariantIdentity:
    """A reconcilable variant identity.

    ``key()`` is what two engines can be joined on. ``label()`` is what a human
    reads. They are deliberately different methods so that using a display
    label as a join key is a visible choice rather than an accident.
    """

    gene: str
    assembly: str = MTB_ASSEMBLY
    chrom: Optional[str] = None
    pos: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    consequence: Consequence = Consequence.UNKNOWN

    @property
    def is_coordinate_resolved(self) -> bool:
        return self.pos is not None and bool(self.ref) and bool(self.alt)

    def key(self) -> str:
        """Canonical identity. Coordinates where we have them."""
        if self.is_coordinate_resolved:
            return f"{self.assembly}:{self.chrom or '.'}:{self.pos}:{self.ref}>{self.alt}"
        # Unresolved: fall back to a namespaced label, and say so in the key so
        # a downstream consumer can tell these apart from real coordinates.
        return f"unresolved:{self.assembly}:{self.gene}:{self._change()}"

    def label(self) -> str:
        """Display label. Not an identity — never join on this across engines."""
        return f"{self.gene}_{self._change()}"

    def _change(self) -> str:
        if self.hgvs_p:
            return self.hgvs_p
        if self.hgvs_c:
            return self.hgvs_c
        if self.is_coordinate_resolved:
            return f"{self.ref}{self.pos}{self.alt}"
        return "unknown"

    def __str__(self) -> str:  # pragma: no cover - convenience only
        return self.label()


@dataclass
class Variant:
    """An observed variant, normalised from any input format."""

    identity: VariantIdentity
    vaf: Optional[float] = None
    depth: Optional[int] = None
    alt_depth: Optional[int] = None
    region: Region = Region.CODING
    filters: tuple[str, ...] = ()
    qual: Optional[float] = None
    genotype_quality: Optional[int] = None
    platform: Optional[str] = None

    @classmethod
    def of(cls, gene: str, change: str, **kwargs: Any) -> "Variant":
        """Ergonomic constructor for tests, examples and simple adapters.

        ``change`` is sorted into ``hgvs_p`` or ``hgvs_c`` by shape. This is a
        convenience for inputs that carry only a label; it does NOT invent
        coordinates, so the resulting identity reports itself as unresolved.
        """
        identity_fields = {
            "assembly", "chrom", "pos", "ref", "alt", "hgvs_c", "hgvs_p",
            "consequence",
        }
        ident_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in identity_fields}
        if "hgvs_p" not in ident_kwargs and "hgvs_c" not in ident_kwargs:
            if _looks_like_protein_change(change):
                ident_kwargs["hgvs_p"] = change
            else:
                ident_kwargs["hgvs_c"] = change
        identity = VariantIdentity(gene=gene, **ident_kwargs)
        return cls(identity=identity, **kwargs)

    # -- passthroughs, so call sites read naturally -----------------------
    @property
    def gene(self) -> str:
        return self.identity.gene

    @property
    def change(self) -> str:
        return self.identity._change()

    @property
    def consequence(self) -> Consequence:
        return self.identity.consequence

    @property
    def silent(self) -> bool:
        return self.identity.consequence.is_silent

    @property
    def passed_filters(self) -> bool:
        return not self.filters or set(self.filters) <= {"PASS", "."}

    def key(self) -> str:
        return self.identity.key()

    def label(self) -> str:
        return self.identity.label()


def _looks_like_protein_change(change: str) -> bool:
    """Heuristic split between protein and nucleotide notation."""
    c = change.strip()
    if c.startswith(("c.", "n.", "g.", "-")):
        return False
    if c.startswith("p."):
        return True
    lowered = c.lower()
    if any(tok in lowered for tok in ("stop", "fs", "del", "ins", "dup", "*")):
        return True
    # e.g. "S315T", "D94G": upper, digits, upper
    return bool(c) and c[0].isalpha() and c[0].isupper() and any(ch.isdigit() for ch in c)


@dataclass(frozen=True)
class LocusCoverage:
    """Independent evidence that a locus was actually sequenced.

    ``callable_fraction`` is the fraction of the locus at or above the depth
    floor. ``None`` means unknown — which is treated as *not* callable, never
    as callable.
    """

    locus: str
    mean_depth: Optional[int] = None
    callable_fraction: Optional[float] = None
    source: str = "absent"

    def __post_init__(self) -> None:
        if self.callable_fraction is not None and (
            not math.isfinite(self.callable_fraction)
            or not 0 <= self.callable_fraction <= 1
        ):
            raise ValueError("callable_fraction must be finite and in [0, 1]")
        if self.mean_depth is not None and (
            not math.isfinite(self.mean_depth) or self.mean_depth < 0
        ):
            raise ValueError("mean_depth must be finite and nonnegative")

    def is_callable(self, depth_floor: int, fraction_floor: float) -> bool:
        if self.callable_fraction is None:
            return False
        if self.callable_fraction < fraction_floor:
            return False
        if self.mean_depth is not None and self.mean_depth < depth_floor:
            return False
        return True

    def describe(self) -> str:
        frac = ("unknown" if self.callable_fraction is None
                else f"{self.callable_fraction:.0%}")
        depth = "unknown" if self.mean_depth is None else f"{self.mean_depth}x"
        return f"{self.locus}: {frac} callable at {depth} (source: {self.source})"


@dataclass(frozen=True)
class EngineRef:
    """Which tool, at which version, over which database, said this."""

    name: str
    version: str = "unknown"
    database: Optional[str] = None
    database_version: Optional[str] = None

    def __str__(self) -> str:
        s = f"{self.name} {self.version}"
        if self.database:
            s += f" / {self.database} {self.database_version or '?'}"
        return s


@dataclass
class DrugEvidence:
    """One lane's opinion about one drug, with its basis attached."""

    drug: str
    call: Call
    tier: Tier
    lane: Lane
    confidence: Optional[float] = None
    variant: Optional[VariantIdentity] = None
    rationale: str = ""
    who_grade: Optional[str] = None
    engine: Optional[EngineRef] = None
    limitations: tuple[str, ...] = ()
    #: The exact catalogue version that produced this call, and the exact rule
    #: (which entry, which grade) that fired. Needed so multi-site aggregation
    #: can tell "different catalogue versions disagree" apart from "the same
    #: catalogue version was applied differently" — two different governance
    #: problems with two different fixes. ``None`` where a lane is not
    #: catalogue-driven (efflux, VUS, minority-allele).
    catalogue_version: Optional[str] = None
    rule_id: Optional[str] = None
    #: Set only by an adapter whose engine performed its own callable-locus
    #: assessment and asserts the drug's loci were adequately covered (e.g.
    #: Mykrobe distinguishing "S" from "N"). It lets orchestration honour a
    #: validated engine's coverage logic instead of discarding its susceptible
    #: call — and the assertion is attributed to that engine in the report.
    #: Never set by a lane that only inspects a variant list.
    asserts_coverage: bool = False
    sample_id: Optional[str] = None
    observation_id: Optional[str] = None
    # Variant-level non-association is not a susceptible isolate verdict.
    scope: str = "isolate"
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.scope not in ("isolate", "variant", "historical", "determinant"):
            raise ValueError(f"unknown evidence scope: {self.scope}")
        if self.confidence is not None and (
            not math.isfinite(self.confidence) or not 0 <= self.confidence <= 1
        ):
            raise ValueError("confidence must be finite and in [0, 1]")
        # The invariant that makes the rest of the system safe.
        if self.call is Call.RESISTANT and not self.tier.may_establish_resistance:
            raise ValueError(
                f"{self.lane.value} lane emitted RESISTANT on tier "
                f"{self.tier.value}; only catalogued or phenotypic evidence may "
                f"establish resistance. Use INDETERMINATE to withhold "
                f"susceptibility without asserting a mechanism."
            )

    @property
    def variant_key(self) -> Optional[str]:
        return self.variant.key() if self.variant else None

    @property
    def variant_label(self) -> Optional[str]:
        return self.variant.label() if self.variant else None

    @property
    def source_name(self) -> str:
        return self.engine.name if self.engine else self.lane.value


@dataclass
class Discordance:
    """Two or more sources reached materially different verdicts."""

    drug: str
    calls: tuple[str, ...]
    sources: tuple[str, ...]
    note: str


@dataclass
class DrugResult:
    """The reconciled position on one drug, with all evidence retained."""

    drug: str
    call: Call
    tier: Tier
    confidence: Optional[float] = None
    evidence: list[DrugEvidence] = field(default_factory=list)
    discordance: Optional[Discordance] = None
    coverage: list[LocusCoverage] = field(default_factory=list)
    reason: Optional[str] = None
    genomic_call: Optional[Call] = None
    phenotypic_call: Optional[Call] = None
    assay_status: str = "not_assessed"

    @property
    def permits_use(self) -> bool:
        """The single gate for regimen eligibility.

        Deliberately narrow: only an established, coverage-backed susceptible
        call lets a drug count toward a regimen.
        """
        return self.call is Call.SUSCEPTIBLE

    @property
    def predicted(self) -> bool:
        return self.tier is Tier.PREDICTED

    @property
    def marker(self) -> str:
        return {
            Call.RESISTANT: "R",
            Call.SUSCEPTIBLE: "S",
            Call.INDETERMINATE: "?",
            Call.NOT_ASSESSED: "-",
            Call.NO_CALL: "x",
            Call.UNSUPPORTED: "n/a",
        }[self.call] + self.tier.marker


@dataclass
class HeteroresistanceFinding:
    """Result of the minority-allele assessment for one variant.

    Carries ``assessable`` separately from ``detected`` so "we could not look"
    is never rendered as "nothing found".
    """

    variant_label: str
    variant_key: str
    drug: str
    assessable: bool
    detected: bool
    note: str
    vaf: Optional[float] = None
    depth: Optional[int] = None
    alt_depth: Optional[int] = None
    limit_of_detection: Optional[float] = None
    platform: Optional[str] = None
    #: Posterior probability that this minority allele reflects a true
    #: resistance-conferring subpopulation, under an explicit measurement
    #: model (see ``modules.calibration``). ``None`` — not zero, not omitted —
    #: whenever a calibrated detection curve or a prior is not configured for
    #: this platform/drug/lineage: an absent calibration is not evidence of
    #: absence of resistance.
    posterior_resistance_probability: Optional[float] = None
    posterior_interval: Optional[tuple[float, float]] = None
    calibration_source: str = "uncalibrated default"


@dataclass
class MechanismHypothesis:
    """A named, testable hypothesis for why a call was withheld — not a call.

    Generalises the VUS workbench's "research queue, not a dead end" idea to
    every lane that reports ``INDETERMINATE`` without a graded mechanism (today:
    the efflux/regulatory lane). ``evidence_gaps`` names what is missing;
    ``resolving_experiments`` names what would resolve it. Neither field
    upgrades the underlying call — only catalogued or phenotypic evidence may
    do that (see ``Tier.may_establish_resistance``).
    """

    variant_label: str
    variant_key: str
    gene: str
    drug: str
    hypothesis: str
    evidence_gaps: list[str] = field(default_factory=list)
    resolving_experiments: list[str] = field(default_factory=list)


@dataclass
class RegimenAssessment:
    """Guideline *eligibility* for one named regimen.

    Note what this does not contain: a list of drugs to give. Regimen
    construction depends on treatment history, disease site, pregnancy,
    comorbidity, interactions, toxicity, availability, baseline ECG and
    national policy — none of which this tool sees.
    """

    name: str
    drugs: list[str]
    min_required: int
    usable: list[str]
    resistant: list[str]
    unestablished: dict[str, str]
    eligible: bool
    shortfall: int
    note: str


@dataclass
class EligibilityReport:
    assessments: list[RegimenAssessment] = field(default_factory=list)
    any_eligible: bool = False
    resistant_drugs: list[str] = field(default_factory=list)
    unestablished_drugs: dict[str, str] = field(default_factory=dict)
    summary: str = ""
    requires_clinical_review: bool = True


@dataclass
class Dimension:
    """One axis of evidence about a variant of unknown significance.

    ``available`` is mandatory and load-bearing: a dimension with no data
    source reports itself unavailable rather than defaulting to a number.
    """

    name: str
    available: bool
    value: Any = None
    source: str = "none"
    note: str = ""


@dataclass
class VUSPriority:
    """A ranked candidate for functional validation — not a resistance call."""

    variant_label: str
    variant_key: str
    gene: str
    drug: Optional[str]
    priority: str                       # high | moderate | low | insufficient-data
    dimensions: list[Dimension] = field(default_factory=list)
    recommended_experiment: str = ""
    data_gaps: list[str] = field(default_factory=list)
    score: Optional[float] = None

    @property
    def rankable(self) -> bool:
        return self.score is not None


@dataclass
class QCFinding:
    check: str
    status: str        # pass | warn | fail | not_performed
    detail: str

    @property
    def blocking(self) -> bool:
        return self.status == "fail"


@dataclass
class Provenance:
    tool: str
    version: str
    organism_profile: str
    profile_version: str
    reference_assembly: str
    depth_floor: int
    callable_fraction_floor: float
    catalogue: Optional[EngineRef] = None
    engines: list[EngineRef] = field(default_factory=list)
    coverage_source: str = "absent"
    demo_mode: bool = False
    generated_utc: Optional[str] = None
    input_sha256: Optional[str] = None
    analysis_fingerprint: Optional[str] = None
    policy_version: str = "research-default"
    source_hashes: dict[str, str] = field(default_factory=dict)


@dataclass
class AnalysisReport:
    sample_id: str
    drug_results: list[DrugResult]
    provenance: Provenance
    eligibility: EligibilityReport = field(default_factory=EligibilityReport)
    heteroresistance: list[HeteroresistanceFinding] = field(default_factory=list)
    vus_priorities: list[VUSPriority] = field(default_factory=list)
    mechanism_queue: list[MechanismHypothesis] = field(default_factory=list)
    discordances: list[Discordance] = field(default_factory=list)
    qc: list[QCFinding] = field(default_factory=list)
    qc_warnings: list[str] = field(default_factory=list)
    lane_counts: dict[str, int] = field(default_factory=dict)
    demo_mode: bool = False
    context: dict[str, Any] = field(default_factory=dict)
    investigations: list[dict[str, Any]] = field(default_factory=list)
    follow_up: list[dict[str, Any]] = field(default_factory=list)
    report_schema: str = "myconductor.report.v1"
    analysis_manifest: dict = field(default_factory=dict)

    @property
    def resistant_drugs(self) -> list[str]:
        return sorted(r.drug for r in self.drug_results if r.call.is_resistant)

    @property
    def usable_drugs(self) -> list[str]:
        return sorted(r.drug for r in self.drug_results if r.permits_use)

    @property
    def unestablished_drugs(self) -> dict[str, str]:
        return {
            r.drug: (r.reason or r.call.reason_unestablished or "unknown")
            for r in self.drug_results
            if not r.call.is_established
        }

    def result_for(self, drug: str) -> Optional[DrugResult]:
        for r in self.drug_results:
            if r.drug == drug:
                return r
        return None
