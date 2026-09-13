# mycobench — validating Myconductor against real data

`mycobench` runs validated AMR engines over public mycobacterial isolates and
reports what the comparison actually supports. It exists because the honest
answer to "is Myconductor any good?" is currently *unknown*, and the only way
to change that is to measure it against real isolates with real phenotypes.

> Nothing here has been run end to end yet. The machinery is built and tested;
> the Docker daemon, the downloads and the compute are yours to start. What this
> document describes is runnable, not run.

## Two tracks, and why they must not be conflated

| Track | Needs | Measures | Panels |
|---|---|---|---|
| **Concordance** | nothing beyond reads | agreement with TB-Profiler | `nigeria-v1`, `nigeria-full` |
| **Accuracy** | paired DST or MIC | sensitivity, specificity, VME, ME | `cryptic-phenotyped` |

Concordance is cheap and available today. It is also **not a measure of
correctness**: two engines reading the same catalogue will agree with each other
whether or not either is right, and their errors are correlated rather than
independent. A concordance figure tells you the orchestration layer reproduces
a reference implementation. It tells you nothing about whether either matches
biology.

Accuracy requires a laboratory phenotype. The public Nigerian BioSamples do not
have one — a scan of 38 BioSamples across 12 BioProjects found no
susceptibility, resistance, DST or MIC attribute of any kind — so the accuracy
track is built from the CRyPTIC compendium, which pairs 12,289 sequenced
isolates with MICs for 13 drugs.

The report states at the top of the page which track produced its numbers,
because "96% agreement" reads like accuracy to anyone skimming.

## Prerequisites

```bash
# Docker must be running: the download and profile stages need it.
docker info >/dev/null && echo "daemon up"

# NCBI asks callers to identify themselves. A key raises the request
# allowance from 3/s to 10/s, which matters on a 213-row panel.
export NCBI_EMAIL="you@example.org"
export NCBI_API_KEY="..."     # https://www.ncbi.nlm.nih.gov/account/
```

Disk: budget roughly three times the download. `nigeria-v1` is ~10 GB of reads
and ~30 GB with alignments and tool output.

## The panels

See [`cohorts/README.md`](../cohorts/README.md) for composition and limitations.
In short:

- **`nigeria-v1`** — 25 isolates, 4 per BioProject, 16.8 Gbases. Frozen and
  NCBI-verified; its lock is committed. Measured download **9.6 GB**
  (7.8–10.6 GB), ~29 GB with working space.
- **`nigeria-full`** — all 213 geo-confirmed Nigerian isolates, 65.7 Gbases.
  Measured download **37.7 GB**, ~113 GB with working space. Severely
  study-clustered: 113 runs from one BioProject, giving an **effective *n* of
  2.7**. Capping at 4 per BioProject raises that to **7.5**.
- **`ntm-controls`** — 12 non-tuberculous mycobacteria that must be **refused**.

### Origin is confirmed, not assumed

An SRA free-text search for `Nigeria` returns 37 runs whose deposited
`geo_loc_name` is `China: East China Sea, Xiangshan Bay`. Discovery therefore
resolves each candidate's BioSample and requires the controlled `geo_loc_name`
field to confirm the country; the 134 BioSamples with no `geo_loc_name` are
excluded too, because unconfirmed is not confirmed.

```bash
mycobench discover-cohort --organism "Mycobacterium tuberculosis complex" \
  --country Nigeria --prefix ng --out-dir curation/nigeria
```

`rejected.tsv` keeps every excluded run with its reason and the evidence that
produced it, so the filter is auditable.

## Running it

### 1. Verify the panel against NCBI

```bash
mycobench validate-cohort --cohort nigeria-v1 --online --country Nigeria \
  --write-lock cohorts/nigeria-v1.lock.json
```

The lock records the evidence retrieved plus a checksum of the exact sheet
verified. Later, anyone can check the sheet has not drifted:

```bash
mycobench validate-cohort --cohort nigeria-v1 \
  --verify-lock cohorts/nigeria-v1.lock.json
```

### 2. See the download before starting it

```bash
mycobench estimate-download --cohort nigeria-v1
```

### 3. Get the real catalogue

The bundled catalogue is a 14-entry illustrative subset. The published WHO
catalogue is MIT-licensed and machine-readable, and it also supplies the
**genomic coordinates** Myconductor's bundled profile deliberately leaves null:

```bash
mycobench fetch-catalogue --out-dir data/who --write-profile
```

This pins WHO's repository to a commit, checksums each file, parses the
114-column master file by column *name*, and joins the coordinates file on
`variant` — keeping every coordinate spelling, because one variant has both an
MNV and an SNV form and matching only one defeats the purpose.

It also derives the per-drug tier 1 / tier 2 gene sets from WHO's own `tier`
column, replacing the illustrative gene lists.

### 4. Run the stages

Download and profile need their containers; interpret onward do not.

```bash
# Preview the exact commands instead of running them:
mycobench run --cohort nigeria-v1 --write-script run_v1.sh

# Or run stage by stage:
benchmark/docker_run.sh sra        mycobench run --cohort nigeria-v1 --stages download
benchmark/docker_run.sh tbprofiler mycobench run --cohort nigeria-v1 --stages profile
mycobench run --cohort nigeria-v1 --stages interpret validate report \
  --catalogue data/who/mtb_amr_catalogue.ingested.json
```

A failing isolate is recorded and skipped, never fatal: one transient SRA error
must not discard hours of successful downloads. Each stage writes
`<stage>_status.tsv` and aborts only when nothing succeeded.

### 5. The accuracy track

```bash
mycobench fetch-phenotypes --out-dir data/cryptic \
  --cohort-out cohorts/cryptic-phenotyped.tsv --limit 500

mycobench run --cohort cryptic-phenotyped \
  --phenotypes data/cryptic/cryptic_phenotypes.tsv
```

Three details decide whether the resulting numbers mean anything, and all three
are handled rather than glossed over:

- **Quality grade.** CRyPTIC grades each phenotype HIGH/MEDIUM/LOW. Only HIGH
  is admitted; measuring a predictor against a LOW-quality phenotype measures
  the phenotype.
- **Censored MICs.** Values arrive as `<=0.25` or `>4.0`. The censoring
  direction is retained rather than silently parsed into a measurement.
- **Withdrawn samples.** The release ships its own exclusion list, and anything
  on it is dropped before the join.

## How accuracy is scored

Myconductor can decline to call a drug. A naive confusion matrix has nowhere to
put that, and the usual workaround — folding "not assessed" into "susceptible" —
is exactly the defect the interpretation engine was rebuilt to remove. Putting
it back in the scoring would hide the failure it was designed to surface.

So **every accuracy figure is conditional on a call being made, and is reported
beside the call rate.** Neither means anything alone:

- a predictor that abstains on every hard isolate posts excellent accuracy on
  the easy remainder;
- a predictor that answers everything confidently posts a high call rate with
  bad errors.

The pre-registered targets require both, and a drug below its minimum evaluable
count is reported **underpowered** — a state distinct from failure, because not
having measured something is not the same as having measured it and found it
wanting.

### Error classes

| Name | Definition | Consequence |
|---|---|---|
| **VME** | predicted susceptible, phenotypically resistant | a patient is given a drug that will not work |
| **ME** | predicted resistant, phenotypically susceptible | a usable drug is withheld |

VME rate is *identically* `1 - sensitivity`, and ME rate is `1 - specificity`.
[`thresholds.py`](../mycobench/thresholds.py) therefore declares sensitivity and
specificity and **derives** the error ceilings. An earlier version declared both
independently and six of eleven drugs ended up with contradictory values, so one
bound silently overrode the other.

### Pre-registration

Targets are fixed in version control with a content hash the report prints, so a
reader can tell whether the bar was set before or after the measurement. They
are this project's declared targets, informed by WHO target product profile
expectations — **not WHO-endorsed thresholds**, and clearing one is not
regulatory acceptance.

### Continuous, lineage-stratified monitoring

A single run's accuracy figure says nothing about whether it holds across
lineages — a panel dominated by lineage 4 cannot speak to lineage 1. So every
`mycobench run` that supplies paired phenotypes *and* lineage data (TB-Profiler's
`main_lineage`/`sub_lineage`) folds its isolates into a persistent,
per-drug/per-lineage confusion-matrix state (`mycobench/monitoring.py`),
written to `results/monitoring_state.json`. Re-running the same cohort is a
no-op — it is keyed by cohort path, so isolates are never double-counted — but
a *new* cohort accumulates on top of everything ingested before it.

The state reports two things per drug: the usual `AccuracyResult` bounds
(sensitivity/specificity/VME/ME, Wilson intervals) **within each lineage**,
and whether the drug has been measured in enough independent lineages
(`MIN_LINEAGES_FOR_GENERALIZABILITY = 2`, mirroring the federated
catalogue-learning module's own confounding bar) to call the figure
`generalizable` at all. Both render as a dedicated section
(`lineage_accuracy.tsv` + the report's "Lineage-stratified accuracy" section)
in the same offline HTML report every run already produces — the closest
thing this no-server tool has to a dashboard every participating site can
regenerate and compare.

## Species controls

The NTM panel exists to be refused. TB-Profiler targets MTBC only and the
organism profile is written against H37Rv, so a drug verdict for
*M. abscessus* would mean the pipeline interpreted the wrong organism against
the wrong catalogue. Any `RESISTANT` or `SUSCEPTIBLE` call on those 12 isolates
is scored as a **failure**, not a result.

## What a run produces

```
results/
  run_manifest.json        versions, image digests, catalogue commit,
                           pre-registration hash, stage counts
  download_status.tsv      per isolate: ok / failed / skipped, with reasons
  profile_status.tsv
  interpret_status.tsv
  engine_calls.tsv         TB-Profiler's per-drug findings, and which parser read them
  myconductor_calls.tsv    Myconductor's per-drug verdict, tier and reason
  concordance.tsv          per-drug agreement
  accuracy.tsv             per-drug errors and verdict (only with phenotypes)
  lineage_accuracy.tsv     per-drug, per-lineage errors and generalizability
  monitoring_state.json    cumulative confusion counts across every cohort
                           ever ingested (see "Continuous, lineage-stratified
                           monitoring" above)
  species_control.tsv      NTM refusals
  mycobench_report.html    the report
  <sample>/                per-isolate TB-Profiler output and Myconductor report
```

Every figure in the report comes from those tables, which sit beside it so a
reader can recompute rather than trust.

## Parser validation status

| Parser | Status |
|---|---|
| `tb-profiler collate` TSV | **validated** against a real artifact from the TB-Profiler repository, whose exact header and `rpoB p.Ser450Leu (1.00)` value format are pinned in the tests |
| TB-Profiler results JSON | written from the documented schema, **not** validated against real output; marked as such on everything it emits, and fails loudly on an unrecognised shape |
| WHO master + coordinates | column lookup by name, tested against fixtures using the real column names |
| CRyPTIC reuse table | tested against the real header |

Which parser produced a call is recorded in `engine_calls.tsv`.

## The third track: conditional effect sizes for uncertain variants

The two tracks above measure Myconductor. This one measures *variants* — the
20,843 entries the WHO catalogue grades "Uncertain significance". They carry
coordinates, so they can be genotyped; nobody knows what they do, so they are
the gap worth attacking.

```bash
mycobench analyse \
  --catalogue data/who/mtb_amr_catalogue.ingested.json \
  --reuse-table data/cryptic/CRyPTIC_reuse_table_20240917.csv \
  --cache-dir /var/cache/mycobench \
  --out-dir results/analysis \
  --limit 400 --sample-seed 20260914 --jobs 8
```

For each candidate variant and drug it compares carriers' MICs against
non-carriers' **within each resistance background**, then pools those
comparisons weighted by how many pairs were actually decidable under censoring.
Holding the background fixed is what separates a determinant from a
hitchhiker: a variant that only ever travels with `rpoB S450L` has no variation
left to explain once that determinant is held constant.

Outcomes are five distinct states — `evidence-of-effect`, `no-evidence`,
`confounded`, `underpowered`, `not-estimable` — because collapsing the last
three into `no-evidence` is what lets a structurally unevaluable variant look
merely unproven.

**Nothing this track produces is a resistance call.** An effect size is a
discovery result. Promoting one to a clinical call requires curation and a
catalogue update, which is the federated ledger's job, not this one's.

### Sample deliberately, or the confounding checks do nothing

`--limit N` alone takes the *first* N rows, and the reuse table is ordered by
site and subject: the first forty isolates are all from `site.02`. Site then
never varies, so the site check has nothing to compare, and the sample is drawn
from one laboratory. `--sample-seed` draws at random instead and is
reproducible. On this release, thirty sampled isolates yielded 186 distinct
catalogued variants where forty sequential ones yielded 142.

Lineage is a harder problem: **the CRyPTIC reuse table has no lineage column**,
so every isolate loaded from it is untyped and the lineage check reports that
it could not be assessed. It does not report a diversity of 1.0 — that would be
indistinguishable from carriers genuinely confined to one lineage. Supply
lineages through `--metadata` to make that check live.

### What this costs, measured

The release's re-genotyped VCFs are **~20 MB compressed and ~178 MB
decompressed each**, about 1.27 million records, because every callable site in
the genome is present including the `0/0` reference calls. Those reference
calls are the point: a catalogued position called `0/0` is the positive
evidence that the locus was examined and the variant absent, which is what lets
a drug reach SUSCEPTIBLE rather than NOT_ASSESSED.

| | |
|---|---|
| Per isolate | ~20 MB download, ~4 s to parse |
| 400 isolates | ~7.4 GB, under an hour with `--jobs 8` |
| All 12,287 | **~248 GB** |

The full compendium is an infrastructure decision, not an incidental download.
Caching is by release-relative path, so a re-run costs nothing and a partial
run resumes.

### What 400 isolates can and cannot answer

Enough for common variants. Not enough for the drugs that matter most for new
regimens: CRyPTIC holds 71 bedaquiline-resistant isolates in 8,535, so a
400-isolate sample is expected to contain roughly three, and no stratified
comparison survives that. Rare-resistance questions need the full compendium
and the storage that implies.

Four drugs cannot support a location statistic at all on this release, because
more than half their observations sit below the lowest tested dilution —
amikacin 69%, delamanid 77%, clofazimine 56%, rifabutin 65% — and the pipeline
refuses them by name rather than reporting a median that is really the
censoring bound.

**Pretomanid has no reference standard on either side.** It is absent from the
WHO catalogue's 15 drugs and from CRyPTIC's 13 MIC columns: nothing supplies a
catalogued genotypic call for it, and nothing supplies a phenotype to score one
against. Its target in `thresholds.py` is carried so a report can say this
explicitly rather than omit the drug.

## What a run will not establish

- **Clinical validity.** Nothing in Myconductor has been clinically evaluated.
- **Correctness, from the concordance track alone.** Agreement between two
  implementations is not evidence either is right.
- **Independence.** Isolates from one BioProject share a sampling frame, a
  laboratory and often a transmission chain. The report prints an effective *n*
  so the clustering is visible rather than implied.
- **Transferability.** Lineage balance is not controlled, and lineage confounds
  both resistance prediction and catalogue transfer.
- **Provenance.** Geographic origin, collection date and isolation source are
  the submitter's deposited claims, not curator-verified facts.

## The gap worth closing

**Nigerian isolates with paired DST results.** They would convert the
concordance track into an accuracy track for African lineages, which is
precisely where published tools are least evaluated and where this project has
a genuine contribution to make. If you hold DST for any isolate in
`nigeria-full`, supplying it is worth more than any amount of further
engineering here.
