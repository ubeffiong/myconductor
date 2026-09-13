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
| Engine adapters (TB-Profiler, Mykrobe, AMRFinderPlus) | **Written from published schemas, never run against real output.** They fail loudly on unrecognised schemas rather than mis-parsing. Add a golden-file test against your own pinned version before relying on one |
| WHO catalogue ingester | **Implemented**; the catalogue data is not redistributed here and must be obtained from WHO |
| VUS prioritisation | **Implemented** as ranking with explicit data gaps. With no annotator configured every dimension reports unavailable, which is the honest output |
| FHIR R4 output | **Structurally correct**, using real HL7 interpretation and data-absent-reason codes. **No LOINC or SNOMED codes**, because inventing them would be worse than omitting them — inject real ones via `TerminologyMap`. Not validated against any implementation guide |
| Read-level analysis (FASTQ/BAM) | **Not implemented.** Species confirmation, contamination, mapping quality, mixed infection and lineage are reported as `not_performed`, with what each would require |
| Trained VUS model | **Not implemented.** The previous hash-derived scorer has been removed; a synthetic annotator remains but cannot reach a report unless `demo_mode=True`, which stamps every page |
| Secure aggregation, differential privacy, federated model training | **Not implemented**, and documented as absent in `federated/transport.py` |
| Molecular docking | **No built-in scorer.** The previous fabricated values are gone; a real backend must be injected |
| Clinical evaluation of any component | **None** |

## Layout

```
myconductor/
  core/        models (the call/tier/identity vocabulary), router, pipeline
  io/          input adapter, callable-locus mask, QC
  catalogue/   organism profile, drug->loci map, illustrative catalogue
  modules/     catalogue, efflux, heteroresistance, synthesis, vus_workbench,
               discovery (cohort-level; NOT wired into the pipeline)
  adapters/    TB-Profiler, Mykrobe, AMRFinderPlus, WHO catalogue ingester
  reporting/   text report, FHIR R4 bundle
  federated/   isolate-level learning, governance, signed transport
docs/ARCHITECTURE.md   design, layer split, extension points
tests/                 unittest suite, including tests/test_safety_invariants.py
examples/run_demo.py   with and without coverage, plus the federated contract
```

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

Phases 0–5 of the rebuild are in place; phases 3–6 are gated on data and
external evaluation rather than on code:

- **Phase 3** needs FASTQ/BAM ingestion, the real catalogue, measured
  per-platform limits of detection, and per-drug error rates against
  pre-registered thresholds on independent isolates — with African lineages
  represented, where published tools are weakest.
- **Phase 4** needs real annotation sources before any VUS may be promoted to a
  predictive call.
- **Phase 5** needs cryptographic review before secure aggregation can be
  claimed, and a data-sharing agreement with benefit-sharing agreed up front.
- **Phase 6** adds one organism profile at a time, each with its own evaluation
  set and its own published error rates. Modules being swappable is an
  architectural property, not a biological claim.

## License

MIT — see `LICENSE`.
