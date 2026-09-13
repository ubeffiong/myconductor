# Myconductor architecture

Myconductor is an **evidence-orchestration layer**: it runs validated AMR
engines, normalises what they produce into one schema, reconciles their
disagreements, keeps uncertainty explicit, and governs how new evidence changes
an interpretation.

It is not another profiler, and the case for it does not rest on existing tools
being weak. TB-Profiler, Mykrobe, GenTB, Deeplex Myc-TB, SAM-TB, MTBseq,
Pathogenwatch, AMRFinderPlus and CARD/RGI collectively cover read processing,
graph genotyping, machine learning, minority-variant detection, mixture and
lineage analysis, transmission analysis, surveillance and curated determinant
detection. Their capabilities are strong but uneven, and few expose a unified,
auditable framework for reconciling different evidence types, incorporating
phenotypic feedback, and governing catalogue updates across sovereign sites.
Myconductor sits in that gap.

## Three layers

```
┌─ Layer 1: Core ───────────┐
│ adapters, normalisation,  │
│ QC, coverage evidence     │
└────────────┬──────────────┘
             ▼
┌─ Shared spine ────────────┐
│ evidence graph:           │
│ multi-lane evidence with  │
│ conflicts retained        │
└──────┬─────────────┬──────┘
       ▼             ▼
┌─ Layer 2 ────┐ ┌─ Layer 3 ──────────┐
│ Interpret    │ │ Discover           │
│ coverage-    │ │ cohort-level only  │
│ aware calls  │ │                    │
└──────┬───────┘ └────────┬───────────┘
       ▼                  ▼
 clinician report   laboratory queue
                          │
            governed catalogue feedback
            (MIC / DST, expert curation)
                          │
                          └──► back to the evidence graph
```

The feedback path is the only route by which new evidence changes an
interpretation, and it passes through expert review. Nothing promotes itself.

Layers 2 and 3 are separated because they have different audiences and different
evidence bars. A clinician reading a susceptibility table should not encounter
candidate compounds beside it, and a discovery run motivated by one patient is a
sample size of one. `core/pipeline.py` does not import `modules/discovery.py`,
and CI enforces that.

## Data flow

```
input (.vcf/.tsv/.json)        callable mask            engine reports
  multiallelic split           (TSV / BED / gVCF)       (TB-Profiler,
  indel normalisation                 │                  Mykrobe,
  FILTER + depth gating               │                  AMRFinderPlus)
  consequence derivation              │                        │
         └──────────────┬─────────────┴────────────────────────┘
                        ▼
                  QC (including
                  not-performed
                  controls)
                        ▼
              multi-lane router  ──► catalogue lane
                        │            efflux/regulatory lane
                        │            VUS workbench lane
                        ▼
              minority-allele assessment
              (per-platform limits of detection)
                        ▼
              coverage-gated reconciliation
              (discordance retained, not voted away)
                        ▼
              guideline eligibility  +  VUS priorities
                        ▼
              text report  /  FHIR R4 bundle
```

## The domain vocabulary

Everything hinges on `core/models.py`, which encodes four commitments.

**1. Absence of evidence is not susceptibility.** `Call` has six states, four of
which mean "no verdict was reached", each carrying its own reason.
`DrugResult.permits_use` is True for exactly one state.

**2. The verdict and its basis are separate axes.** `Call` says *what*; `Tier`
says *how we know* — phenotypic, catalogued, inferred, predicted, or none. Only
`PHENOTYPIC` and `CATALOGUED` may establish resistance, and
`DrugEvidence.__post_init__` raises otherwise. A rule-based lane can withhold
susceptibility without claiming an unvalidated mechanism.

**3. A variant's identity is its coordinates, not its label.** Two engines
cannot be joined on `"katG_S315T"` — that spelling varies by annotator.
`VariantIdentity.key()` uses assembly, coordinate and alleles; `label()` is for
display. They are separate methods so that joining on a label is a visible
choice rather than an accident.

**4. Disagreement is a finding.** `DrugResult` retains every piece of evidence
it received and carries a `Discordance` rather than collapsing to a winner.
Conflicts are never resolved by vote, because engines sharing a catalogue are
not independent evidence.

## Modules

| Stage | File | What it does | Real backend to plug in |
|---|---|---|---|
| Adapter | `io/adapter.py` | Normalises input; splits multiallelic records, trims indels, honours FILTER, derives consequence; rejects records the model cannot express | any variant caller |
| Coverage | `io/callable_mask.py` | Independent evidence that loci were sequenced — the module that breaks the circularity | mosdepth, GATK CallableLoci, gVCF |
| QC | `io/qc.py` | What was checked, and every control that was **not** | upstream pipeline integration |
| Profile | `catalogue/profile.py` | Versioned organism bundle: assembly, drug→loci, regulators, regimens | one profile per organism |
| Catalogue lane | `modules/catalogue.py` | Graded lookup; may establish resistance; stamps `catalogue_version`/`rule_id` on every piece of evidence | full WHO catalogue via the ingester |
| Efflux lane | `modules/efflux.py` | Flags de-repression and promoter effects as `INDETERMINATE`, per affected drug; `mechanism_queue()` turns each into a named hypothesis with evidence gaps and resolving experiments | expression-inference model, RNA evidence |
| VUS workbench | `modules/vus_workbench.py` | Ranks for validation; withholds susceptibility; never predicts resistance | real annotation sources, then a calibrated model |
| Local validation lane | `modules/local_validation.py` + `federated/vus_feedback.py` | Reads a site's own laboratory-validated VUS results back as site-local `Tier.PHENOTYPIC` evidence, reversibly (hash-chained ledger); deliberately not the global catalogue | a deployment's own DST/MIC/efflux-assay results |
| Annotation | `modules/features.py` | Dimension contract with availability attached; default annotator reports unavailable | SIFT/PolyPhen, AlphaFold, population DB |
| Minority alleles | `modules/heteroresistance.py` + `modules/calibration.py` | Per-platform assessment; separates "not assessable" from "nothing found"; when a calibration curve and a prior are configured, attaches a Bayesian posterior (still capped at `Tier.PREDICTED`/`INDETERMINATE`) | read-level caller with error modelling; a deployment's own dilution-series calibration |
| Reconciliation | `modules/synthesis.py` | Coverage-gated calls, discordance, guideline eligibility | external, versioned clinical decision rules |
| Router | `core/router.py` | Runs **all** applicable lanes | — |
| Engine adapters | `adapters/` | Normalise external tools into one evidence schema | the tools themselves |
| BAM coverage | `io/bam_coverage.py` | Automates deriving a `CallableMask` from a BAM via mosdepth/samtools; coverage only, no alignment or variant calling | mosdepth or samtools on PATH |
| Read calling | `io/read_calling.py` | Orchestrates FASTQ -> BAM -> VCF via minimap2/bwa-mem2 + bcftools/GATK; the resulting VCF is treated exactly like one from any other source once it reaches `io/adapter.py` | minimap2/bwa-mem2 + bcftools/GATK on PATH |
| Reporting | `reporting/` | Text + FHIR R4 | LIMS / EHR |
| Federated | `federated/` | Isolate-level learning, governance, signed transport, closed-loop VUS validation, catalogue-version reconciliation across sites, an unreviewed secure-aggregation reference implementation (not wired into the transport) | a cryptographer's review before the reference implementation could be trusted, then wiring it in |
| Discovery | `modules/discovery.py` | Cohort-level prioritisation brief | DEG, BLASTP, a real docking backend |

## Extension points

A lane implements `modules/base.py::VariantModule`: `applies_to(variant) -> bool`
and `evaluate(variant) -> list[DrugEvidence]`. Returning a list is required
because one variant can affect several drugs. An empty list means abstention —
it must never mean susceptibility.

An engine adapter implements `adapters/base.py::EngineAdapter`:
`parse(path) -> EngineReport`. An adapter may set `asserts_coverage=True` on
susceptible evidence when its engine performed its own callable-locus assessment
(Mykrobe's `S` versus `N`); that assertion is attributed to the engine in the
report. No lane that merely inspects a variant list may set it, and a test
enforces this.

An organism profile is JSON: `catalogue/drug_loci.json` plus
`catalogue/drugs.json`. Adding an organism means shipping a profile, not editing
the engine — but a profile is only as good as its own evaluation data, and every
bundled profile reports `ships_validated = False`. A second bundled profile
(`catalogue/organisms/mabscessus/`, selected via `Myconductor(organism=...)` or
`myconductor analyze --organism mabscessus`) exercises this claim rather than
just asserting it, and exposed one real gap while doing so:
`modules/vus_workbench.py`'s gene->drug context was a flat, MTBC-only
dictionary with no organism awareness, so a gene from a foreign profile that
happened to share a name with an MTBC gene (M. abscessus's own `rpoB`, say)
would have been silently attributed to an MTBC drug it has no declared
relationship to. `VUSWorkbench` now intersects that dictionary with the active
profile's own `loci` before using it — a pure narrowing, so the bundled MTBC
pipeline's behaviour is unchanged (every gene the dictionary names is already
in the MTBC profile's own loci), while a profile that doesn't declare a gene
can no longer have it silently borrowed from another organism's knowledge.
See `tests/test_safety_invariants.py::OrganismProfilesDoNotLeakIntoEachOther`.

Everything in the pipeline is constructor-injected.

## What is deliberately absent

Each of these was removed or withheld because a plausible-looking placeholder is
more dangerous than a gap:

- **No presumed-usable fallback.** Missing evidence yields `NOT_ASSESSED`.
- **No hash-derived feature scores.** The previous conservation and structural
  "scores" were SHA-256 digests of the variant name. A synthetic annotator
  remains for interface demonstration, but the pipeline refuses it unless
  `demo_mode=True`, which stamps every report.
- **No docking scores.** There is no built-in scorer; a real backend must be
  injected. A binding affinity establishes none of inhibition, permeability,
  whole-cell activity, selectivity, toxicity, or resistance barrier.
- **No invented terminology codes.** The FHIR bundle carries real HL7
  interpretation and data-absent-reason codes and omits LOINC and SNOMED
  entirely, with the gap listed in the bundle itself.
- **No fabricated reference annotation.** Locus lengths are `null` throughout the
  bundled profile; where a length is unknown, the input must supply an explicit
  callable fraction.
- **No claim of secure aggregation or differential privacy in the shipped
  transport.** `federated/transport.py` authenticates sites and gates
  disclosure. It does not hide a site's contribution from the coordinator,
  and it says so. `federated/secure_aggregation.py` is a separate, explicitly
  unreviewed reference implementation of the pairwise-masking arithmetic
  behind Bonawitz et al. 2017 — it exists to make the idea testable, is
  gated behind `acknowledged=True`, does not solve pairwise key agreement or
  dropout tolerance, and is not wired into `transport.py`. See its own module
  docstring before reading its presence here as more than that.
- **No bundled minority-allele calibration or resistance prior.**
  `modules/calibration.py::CalibrationTable` ships empty and
  `NullPriorSource` is the default; a deployment supplies its own
  dilution-series curve and, optionally, a prior built from its own federated
  evidence. Without both, no posterior is computed — the flat, conservative
  `PLATFORM_LOD` floors remain the fallback, exactly as before.

## Federated learning: the statistical core

The previous implementation counted variants and promoted one at 20 observations
with 75% resistance. Resistance belongs to an **isolate–drug observation**, not
to a variant: a resistant isolate carries the causal variant, lineage markers it
inherited, and hitchhikers in linkage. Counting per variant credits all of them,
so lineage markers become "resistance determinants" and are then distributed to
other sites in a signed catalogue — costing patients usable drugs.

`federated/catalogue_update.py` therefore counts each isolate once per drug and
blocks promotion unless the association survives:

- **lineage stratification** — it must hold within at least two lineages;
- **co-occurrence** — a variant that almost always appears with a known
  determinant cannot be credited with the phenotype;
- **site diversity** — single-site evidence is one laboratory's systematic error;
- **susceptible controls** — without susceptible isolates in the drug's tested
  cohort, specificity is unknown;
- **a variant-absent comparison group** — if every isolate tested for the drug
  carries the variant, a 100% association says nothing, because it cannot be
  separated from the cohort's baseline resistance.

The last two are why `aggregate` runs two passes and tracks each drug's whole
tested cohort, not only the isolates carrying the variant. Requiring susceptible
isolates that *carry* the variant would be the wrong control, and would penalise
precisely the strongest determinants — which are rarely seen in susceptible
isolates.

Clearing that bar produces a *candidate*, not a catalogue entry. Candidates enter
a review queue for expert curation; every transition is written to an
append-only hash-chained ledger, and approvals can be rolled back. There is no
automatic promotion path.

## Two more federated concerns: local validation and catalogue-version drift

`federated/vus_feedback.py` is a *smaller, local* counterpart to the module
above — not a replacement for it. A laboratory's validation of one VUS
(a DST, an MIC shift, an efflux-inhibitor assay) is ingested as a
`VUSValidationRecord`, linked to its isolate, lineage and site, and read back
through `modules/local_validation.py` as `Tier.PHENOTYPIC` evidence *at that
site only*. It reuses `EvidenceLedger` for its own audit trail and never
writes to, or bypasses, the review queue above; promoting a variant into the
shared catalogue still requires the full confounding-check and
expert-curation path. A conflicting set of local results is reported
`INDETERMINATE`, never resolved by majority.

`federated/catalogue_governance.py` addresses a different failure mode:
catalogue heterogeneity, not confounding. Every catalogue-lane and
local-validation piece of `DrugEvidence` now carries `catalogue_version` and
`rule_id`. When sites disagree on the same variant/drug, `reconcile_sites`
classifies the disagreement as catalogue-version drift (different sites
running different catalogue versions — the fix is propagation) or a
same-version rule discordance (the fix is in the pipeline). Neither category
is resolved by outvoting the minority site.

## Reproducibility

CI runs the suite on Python 3.9–3.13, executes the demo and the example
end to end, and greps for the specific defects this codebase was rebuilt to
remove. Still outstanding: containerised environments, a workflow language
(Nextflow/WDL/CWL), version-pinned external databases, model cards, golden
benchmark datasets with expected outputs, SBOM and dependency scanning, signed
releases, and an archival DOI.

## Scope

A research and decision-support scaffold, not a clinical device. Nothing here
has been clinically evaluated. The hard part of turning it into a product is not
the code — it is data breadth, external evaluation against phenotypic results
(especially for bedaquiline and pretomanid, and especially on under-represented
lineages), and regulatory approval.
