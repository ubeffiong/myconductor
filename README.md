# Myconductor

**An evidence-orchestration layer for antimicrobial-resistance genomics.**

*myco*bacterium + *conductor* — the layer that conducts the existing tools.

> **Status: research scaffold.** Not a clinical device. No component has been
> clinically evaluated, and no sensitivity, specificity or error rate has been
> established for anything in this repository. The bundled catalogue is a
> 14-entry illustrative subset, not the WHO catalogue. The engine adapters are
> written from published output schemas and have never been run against real
> tool output. Read [What is actually implemented](#what-is-actually-implemented)
> before drawing any conclusion from its output.

## What this is, and what it is not

Myconductor is **not another resistance predictor**, and it does not claim that
existing tools are inadequate at prediction. TB-Profiler, Mykrobe, GenTB,
Deeplex Myc-TB, SAM-TB, MTBseq, Pathogenwatch, AMRFinderPlus and CARD/RGI each
already do at least one thing well: read processing, graph genotyping, machine
learning, minority-variant detection, mixture and lineage analysis, transmission
analysis, surveillance, or curated determinant detection.

What they provide is strong but **uneven**, and few of them expose a unified,
auditable framework for reconciling these different kinds of evidence,
incorporating phenotypic feedback, and governing how catalogue updates propagate
across sovereign sites. That gap is what Myconductor addresses:

- **engine adapters** that normalise several tools into one evidence schema, and
  report where those tools *disagree* rather than resolving conflicts by vote;
- **coverage-gated calls** — a drug is susceptible only where evidence
  independent of the variant list shows its loci were sequenced;
- **uncertainty as a first-class output** — six call states, so "we did not
  look" is never rendered as "this drug will work";
- **a VUS workbench** that ranks variants for laboratory validation instead of
  predicting resistance for them;
- **governed multi-site learning** at isolate level, with lineage and
  co-occurrence confounding checks, expert curation and a reversible audit trail.

### The rule everything follows

> **Absence of evidence is not susceptibility.**

A drug may have no resistance result because its locus was not covered, the
caller failed, the gene was absent from the input, the input carried variants
but no reference-call information, the profile does not support the drug, or the
mechanism is unknown. None of those mean the drug will work. Only a
coverage-backed `SUSCEPTIBLE` counts toward regimen eligibility.

| State | Meaning | Counts toward a regimen? |
|---|---|---|
| `RESISTANT` | Graded catalogue hit or laboratory phenotype | No |
| `SUSCEPTIBLE` | Wild type **and** required loci shown callable | **Yes — only this** |
| `INDETERMINATE` | Conflicting evidence, an uninterpreted variant, or a lane that abstained | No |
| `NOT_ASSESSED` | Loci not shown callable; caller failed; gene absent | No |
| `NO_CALL` | Genotype at the locus is not reliable | No |
| `UNSUPPORTED` | Drug not covered by this organism profile | No |

Only **catalogued** or **phenotypic** evidence may establish `RESISTANT`. The
efflux rule and the VUS workbench report `INDETERMINATE` instead — withholding
susceptibility, which is the clinically protective action, without asserting a
mechanism nobody measured. This is enforced structurally: constructing a
`RESISTANT` evidence object on an inferred or predicted tier raises.

## Install

```bash
pip install -e .          # or run from the repo root with no install
```

The core has **zero runtime dependencies** (Python 3.9+ standard library only),
so it runs offline. `openpyxl` is optional, and only to read the WHO
catalogue's published `.xlsx`.

## Quickstart

```bash
python -m myconductor demo
```

Analyse your own data. Note that `--mask` is what permits susceptible calls:

```bash
python -m myconductor analyze sample.vcf \
    --mask sample.callable.tsv \
    --platform illumina \
    --fhir sample.fhir.json
```

Reconcile external engines:

```bash
python -m myconductor analyze sample.vcf --mask sample.callable.tsv \
    --tbprofiler sample.results.json \
    --mykrobe sample.mykrobe.json
```

Ingest a real catalogue in place of the bundled subset:

```bash
python -m myconductor ingest-catalogue who_catalogue.csv --out catalogue.json
```

Programmatic:

```python
from myconductor import Myconductor
from myconductor.reporting.render import render_text

report = Myconductor(platform="illumina").analyze(
    "myconductor/data/example_input.vcf",
    mask_path="myconductor/data/example_callable.tsv",
)
print(render_text(report))
print(report.usable_drugs)          # only coverage-backed susceptibles
print(report.unestablished_drugs)   # drug -> why no verdict was reached
```

### Providing coverage evidence

Without a mask, **every drug reports `NOT_ASSESSED`**. That is intended, not a
bug. Three accepted sources:

| Source | How to produce it |
|---|---|
| depth table (TSV) | `mosdepth --by targets.bed --thresholds 10`, then one row per locus: `locus`, `mean_depth`, `callable_fraction` |
| BED mask | GATK `CallableLoci` or `mosdepth --quantize`, with the locus in column 4 |
| gVCF | non-variant reference blocks with `END` and `MIN_DP`; needs locus spans from the reference annotation |

## What the demo shows

The bundled synthetic sample is a pre-XDR-style strain. With coverage evidence
supplied it produces:

- **resistant** — isoniazid, rifampicin, ethambutol, moxifloxacin (catalogued);
- **indeterminate** — bedaquiline and clofazimine (an `Rv0678` change: efflux
  de-repression is plausible but unmeasured, and it correctly reaches *both*
  drugs); pyrazinamide (an uninterpreted `pncA` variant in a required locus);
- **susceptible** — linezolid and pretomanid, because their loci were shown
  callable and nothing was found;
- **not assessed** — amikacin, whose `rrs` record failed the depth floor and
  which is absent from the mask;
- a **minority allele** at 15% for the fluoroquinolone, assessed against the
  Illumina limit of detection rather than a fixed threshold;
- **no eligible regimen**, reported as a statement about the available evidence
  rather than a conclusion that no regimen exists.

Run the same input *without* `--mask` and nothing is susceptible. Comparing the
two is the clearest demonstration of what this rebuild changed:

```bash
python examples/run_demo.py
```

## What is actually implemented

Being precise about this is the point of the project.

| Capability | Status |
|---|---|
| Coverage-gated susceptibility | **Implemented and tested** |
| Six-state calls with reasons | **Implemented and tested** |
| Multi-lane evidence, multi-drug effects | **Implemented and tested** |
| Discordance reporting | **Implemented and tested** |
| Variant normalisation (multiallelic, indel trimming, FILTER, consequence) | **Implemented and tested** |
| Isolate-level federated statistics with confounding checks | **Implemented and tested** |
| Signed submissions, replay protection, k-anonymity gate | **Implemented and tested** |
| Append-only audit ledger, review queue, rollback | **Implemented and tested** |
| Platform-specific minority-allele assessment | **Implemented** — the limits of detection are conservative defaults for triage, not measured limits |
| Calibrated heteroresistance posteriors | **Implemented** (`modules/calibration.py`) as an explicit Bayesian measurement model (a dilution-series detection curve + an isolate-level prior), returning a posterior probability instead of a fixed rule. Ships with **no bundled calibration data or prior** — a deployment supplies its own dilution series and, optionally, a `FederatedPriorSource` built from its own federated evidence. Structurally capped at `Tier.PREDICTED`/`INDETERMINATE`; it can never assert `RESISTANT` |
| Engine adapters (TB-Profiler, Mykrobe, AMRFinderPlus) | **Written from published schemas, never run against real output.** They fail loudly on unrecognised schemas rather than mis-parsing. Add a golden-file test against your own pinned version before relying on one |
| WHO catalogue ingester | **Implemented**; the catalogue data is not redistributed here and must be obtained from WHO |
| VUS prioritisation | **Implemented** as ranking with explicit data gaps. With no annotator configured every dimension reports unavailable, which is the honest output |
| Closed-loop VUS validation | **Implemented** (`federated/vus_feedback.py`, `modules/local_validation.py`, `myconductor validate-vus`). A laboratory's validation of a VUS is ingested, linked to its isolate/lineage/site, and read back as site-local `Tier.PHENOTYPIC` evidence on future reports — reversibly, via a hash-chained ledger (`retract`). This is deliberately **not** the global catalogue: `federated/catalogue_update.py`'s confounding checks and expert-review queue are unchanged and unbypassed |
| Catalogue-version-aware governance | **Implemented** (`federated/catalogue_governance.py`, `myconductor governance-report`). Every catalogue-lane and local-validation piece of evidence now carries `catalogue_version` and `rule_id`; cross-site reconciliation classifies disagreement as catalogue-version drift versus a same-version pipeline/rule bug, and never resolves either by vote |
| Lineage-stratified continuous monitoring | **Implemented** in `mycobench` (`mycobench/monitoring.py`). Every `mycobench run` with paired phenotypes and lineage data folds into a persistent, per-drug/per-lineage accuracy state (`monitoring_state.json`), idempotent per cohort, with a `generalizable` flag (2+ independent lineages) distinct from the isolate-count power check. Rendered as a dashboard section in the existing offline HTML report — the closest thing this no-server tool has to a shared, cross-site view |
| Mechanism research queue | **Implemented** (`AnalysisReport.mechanism_queue`, populated from `modules/efflux.py`). Every efflux/regulatory `INDETERMINATE` now carries a named mechanistic hypothesis, the specific evidence gaps, and the assay(s) that would resolve it — generalising the VUS workbench's "research queue, not a dead end" idea beyond sequence variants |
| BAM → coverage-mask automation | **Implemented** (`io/bam_coverage.py`, `myconductor analyze --bam --bed`). Automates the README's own documented `mosdepth --by ... --thresholds ...` workflow; written from mosdepth's/samtools's documented output schema and **not run against real tool output** in this environment — same caveat as the engine adapters |
| FASTQ → VCF orchestration | **Implemented** (`io/read_calling.py`, `myconductor call-variants`). Thin orchestration over minimap2/bwa-mem2 + bcftools/GATK — not a new aligner or caller, and this project still does not pick these choices *for* a deployment: they are defaults, swappable per call, and **not run against real reads** in this environment. The resulting VCF is treated exactly like a VCF from anywhere else once it reaches `io/adapter.py` |
| A second organism profile | **Implemented**: an illustrative M. abscessus profile (`catalogue/organisms/mabscessus/`, `myconductor analyze --organism mabscessus`) covering three literature-grounded genes (erm(41) T28/C28 sequevar, rrl 2058/2059, rrs 1408 — see the profile's own `sources` field for citations). `ships_validated = False`, exactly like the bundled MTBC profile; five hand-picked entries, not a systematic review. Selecting an organism never changes another organism's behaviour — `modules/vus_workbench.py`'s gene→drug context is now scoped to genes the active profile actually declares, closing a gap where MTBC-specific gene knowledge would otherwise leak into any other organism |
| FHIR R4 output | **Structurally correct**, using real HL7 interpretation and data-absent-reason codes. **No LOINC or SNOMED codes**, because inventing them would be worse than omitting them — inject real ones via `TerminologyMap`. Not validated against any implementation guide. Deliberately does not carry the mechanism research queue — a laboratory work list, not a clinical observation |
| Read-level analysis (FASTQ/BAM) | Coverage and variant calling are now automated (see the two rows above). Species confirmation, contamination, mapping quality and mixed-infection detection are still `not_performed`: they are QC judgements about the alignment this project is not positioned to make correctly without its own validation data, unlike orchestrating an established aligner/caller with their own published accuracy |
| Trained VUS model | **Not implemented.** The previous hash-derived scorer has been removed; a synthetic annotator remains but cannot reach a report unless `demo_mode=True`, which stamps every page |
| Secure aggregation | `federated/transport.py` itself still sends a site's contribution to the coordinator in clear, exactly as before. A **reference implementation** of the pairwise-masking arithmetic (Bonawitz et al. 2017) now exists separately at `federated/secure_aggregation.py`, gated behind `acknowledged=True` — built at explicit user request to make the idea runnable and testable, **not reviewed by a cryptographer**, does not solve pairwise key agreement, and has no dropout tolerance. Nothing in the shipped transport uses it |
| Differential privacy, federated model training | **Not implemented**, and documented as absent in `federated/transport.py`. Cross-site *governance* (catalogue-version reconciliation, above) is implemented and is a different, non-cryptographic problem |
| Molecular docking | **No built-in scorer.** The previous fabricated values are gone; a real backend must be injected |
| Clinical evaluation of any component | **None** |

## Layout

```
myconductor/
  core/        models (the call/tier/identity vocabulary), router, pipeline
  io/          input adapter, callable-locus mask, QC, BAM->mask automation,
               FASTQ->VCF orchestration (read_calling)
  catalogue/   organism profile, drug->loci map, illustrative catalogue;
               organisms/ holds each additional bundled profile (mabscessus)
  modules/     catalogue, efflux (+ mechanism queue), heteroresistance (+
               calibrated posteriors), calibration, synthesis, vus_workbench,
               local_validation, discovery (cohort-level; NOT wired into the
               pipeline)
  adapters/    TB-Profiler, Mykrobe, AMRFinderPlus, WHO catalogue ingester
  reporting/   text report, FHIR R4 bundle
  federated/   isolate-level learning, governance, signed transport,
               closed-loop VUS validation (vus_feedback), catalogue-version
               governance (catalogue_governance), a secure-aggregation
               reference implementation NOT wired into the shipped
               transport (secure_aggregation)
docs/ARCHITECTURE.md   design, layer split, extension points
tests/                 unittest suite, including tests/test_safety_invariants.py
examples/run_demo.py   with and without coverage, plus the federated contract
```

## Validating it against real data

`mycobench` is the harness that measures this thing instead of asserting it. It
runs TB-Profiler and Myconductor over public mycobacterial isolates and reports
what the comparison supports — see [`docs/BENCHMARK.md`](docs/BENCHMARK.md).

```bash
mycobench validate-cohort --cohort nigeria-v1 --online --country Nigeria
mycobench estimate-download --cohort nigeria-v1
mycobench fetch-catalogue --out-dir data/who --write-profile
mycobench run --cohort nigeria-v1 --write-script run_v1.sh
```

Three shipped panels, every accession retrieved from NCBI rather than composed:

| Panel | Rows | Track | Download |
|---|---:|---|---:|
| `nigeria-v1` | 25 | concordance | ~10 GB |
| `nigeria-full` | 213 | concordance | ~40 GB |
| `ntm-controls` | 12 | species control | ~11 GB |

Two things about this are worth stating plainly. **The Nigerian panels cannot
measure accuracy** — their BioSamples carry no DST or MIC metadata, so they
measure agreement with TB-Profiler and nothing more; accuracy comes from the
CRyPTIC compendium via `mycobench fetch-phenotypes`. And **origin is confirmed
on the controlled `geo_loc_name` field, not free text**: an SRA search for
`Nigeria` also returns 37 runs deposited from the East China Sea.

`mycobench fetch-catalogue` also closes the null-coordinate gap, because WHO's
published files supply real genomic coordinates that the bundled illustrative
profile refuses to invent.

## Tests

```bash
python -m unittest discover -s tests -v
```

`tests/test_safety_invariants.py` is separate on purpose: each test there
corresponds to a specific defect in the earlier implementation, and a failure in
any of them is a clinical-safety regression rather than a functional bug. CI
also greps for the removed defects directly, so a reviewer does not have to
remember them.

Note what the suite does **not** do: it does not establish biological
correctness. It verifies that the software behaves as designed. Only comparison
against phenotypic data on real isolates can speak to correctness, and that work
has not been done.

## Scope and safety

Myconductor is a research and decision-support scaffold. It must not be used to
select treatment for a patient.

The report deliberately provides a **resistance evidence summary and guideline
eligibility assessment**, never a proposed regimen. Regimen construction depends
on treatment history, site of disease, age, pregnancy, comorbidity, drug
interactions, toxicity, drug availability, baseline ECG and laboratory findings,
confirmed or inferred cross-resistance, confidence in each call, and national
policy — none of which this tool sees. Clinical decision rules belong outside
it, versioned and independently governed.

Target discovery is a separate cohort-level research workflow. There is no code
path from a patient's report to a drug-discovery run, and `core/pipeline.py` does
not import `modules/discovery.py` — CI enforces this.

## Roadmap

Every part of phases 3–6 that was a software/architecture problem has now been
built; what remains gated is exactly what was always gated — real
measurements, external review and evaluation cohorts that only exist outside
this repository. Building a plausible-looking placeholder for any of those
would be worse than the gap, so none was added.

- **Phase 3** — done in software: coverage can now be derived directly from a
  BAM (`io/bam_coverage.py`), minority-allele calls carry a calibrated,
  explicit Bayesian posterior instead of a fixed rule
  (`modules/calibration.py`), and per-drug/per-lineage error rates are now
  continuously monitored across every cohort ingested, not measured once
  (`mycobench/monitoring.py`, see the HTML report's lineage-accuracy section).
  **Still gated on real data**: measured per-platform/per-drug/per-lineage
  limits of detection (the calibration table ships empty), the real WHO
  catalogue (`mycobench fetch-catalogue` ingests it; the bundled subset stays
  illustrative), and African-lineage isolates with paired phenotypes to
  actually populate the monitoring state and clear its pre-registered
  thresholds. FASTQ->VCF orchestration is now also implemented
  (`io/read_calling.py`, `myconductor call-variants`) as a swappable default
  over minimap2/bwa-mem2 + bcftools/GATK — still not validated against real
  reads in this environment, and still not a claim that these are the right
  choices for any given deployment.
- **Phase 4** — done in software: a validated VUS now feeds back as
  site-local, reversible evidence (`federated/vus_feedback.py`,
  `modules/local_validation.py`, `myconductor validate-vus`), closing the loop
  the workbench was missing. **Still gated**: real annotation sources
  (conservation, structure, population prevalence) before any VUS may be
  promoted to a *predictive* call across sites — this closed loop is
  explicitly local, not that promotion.
- **Phase 5** — done in software: multi-site catalogue-version governance
  (`federated/catalogue_governance.py`, `myconductor governance-report`),
  distinguishing catalogue-propagation gaps from pipeline/rule bugs; and a
  reference implementation of pairwise-masking secure-aggregation arithmetic
  (`federated/secure_aggregation.py`) exists to make the idea testable.
  **Still gated**: cryptographic review before secure aggregation can be
  *claimed* — the shipped transport (`federated/transport.py`) still sends
  every submission to the coordinator in clear and says so; the reference
  module does not solve pairwise key agreement or dropout tolerance, and is
  not wired into the transport. Also still needed: a data-sharing agreement
  with benefit-sharing agreed up front.
- **Phase 6** — done in software, to the same illustrative standard as the
  bundled MTBC profile: a second organism profile now ships
  (`catalogue/organisms/mabscessus/`, `myconductor analyze --organism
  mabscessus`), and the VUS workbench's gene->drug context is now scoped per
  profile so organisms cannot leak knowledge into each other. **Still
  gated**: this is five literature-grounded entries for interface
  demonstration, not the systematic review and published error rates a real
  M. abscessus deployment would need — `ships_validated = False`, same as
  every profile here. Each further organism still needs its own evaluation
  set.

## License

MIT — see `LICENSE`.
