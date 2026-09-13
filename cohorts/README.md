# mycobench cohorts

Three panels ship here. Every accession in them was retrieved from NCBI, not
composed — regenerate any of them with `mycobench discover-cohort` and compare.

| Panel | Rows | Track | Download | Use it for |
|---|---:|---|---:|---|
| [`nigeria-v1.tsv`](nigeria-v1.tsv) | 25 | concordance | ~10 GB | reproducible headline concordance, per-BioProject capped |
| [`nigeria-full.tsv`](nigeria-full.tsv) | 213 | concordance | ~40 GB | full geo-confirmed Nigerian set |
| [`ntm-controls.tsv`](ntm-controls.tsv) | 12 | species control | ~11 GB | proving the pipeline refuses non-MTBC input |

Download sizes are estimates from NCBI base counts at ~0.615 bytes/base; run
`mycobench estimate-download` for the figure on your own disk.

`nigeria-v1` is **frozen**. Published numbers reference it by name and its lock
certifies that exact sheet. New isolates go into `nigeria-full`, a panel of your
own, or a `nigeria-v2` — never appended to v1.

## Two tracks, and why the distinction is load-bearing

**Concordance** compares Myconductor against TB-Profiler on the same isolates.
It needs no phenotype, and it measures agreement — not correctness. Two tools
sharing a catalogue agreeing tells you they read the same catalogue.

**Accuracy** compares a prediction against a laboratory phenotype. It needs
paired DST or MIC data, which the Nigerian BioSamples **do not have**: a scan of
38 BioSamples across 12 BioProjects found no susceptibility, resistance, DST or
MIC attribute of any kind. `phenotype_source` is therefore empty on every
Nigerian row, and no accuracy figure can be computed from these panels.

The accuracy track is built from CRyPTIC instead:

```bash
mycobench fetch-phenotypes --out cohorts/cryptic-phenotyped.tsv
```

That joins the CRyPTIC reuse table (12,289 isolates, 13 drugs, MIC and binary
phenotype with a per-drug quality grade) to public ENA runs. It is not committed
here because it is derived, versioned upstream, and large.

**Nigerian isolates with paired phenotypes would be the single most valuable
addition to this benchmark.** If you hold DST results for any isolate in
`nigeria-full`, supplying them is what converts the concordance track into an
accuracy track for African lineages — where published tools are least evaluated.

## Origin is confirmed on geo_loc_name, not free text

Searching SRA for `"Mycobacterium tuberculosis complex"[Organism] AND
Nigeria[All Fields]` returns 472 runs. That query is **not** a statement of
origin:

```
   155  geo_loc_name = 'Nigeria'
    37  geo_loc_name = 'China: East China Sea, Xiangshan Bay'
    22  geo_loc_name = 'Nigeria:BENUE'
    10  geo_loc_name = 'Nigeria:FCT'
     9  geo_loc_name = 'Nigeria:NIGER'
     5  geo_loc_name = 'Nigeria:PLATEAU'
     4  geo_loc_name = 'Nigeria:KOGI'
     3  geo_loc_name = 'Nigeria:KWARA'
     3  geo_loc_name = 'Nigeria:NASARAWA'
     2  geo_loc_name = 'Nigeria:Ibadan_OyoState'
    84  geo_loc_name absent
```

The 37 East China Sea runs match on free text alone — the word appears
elsewhere in their record. Accepting them would have put marine metagenomic
samples in a Nigerian clinical cohort and nothing downstream would have
noticed. `discover-cohort` therefore resolves every candidate's BioSample and
requires `geo_loc_name` to begin with the requested country; 134 BioSamples with
no `geo_loc_name` at all are also excluded, because absent is not confirmed.

After that filter: **213 runs, 9 BioProjects, 65.7 Gbases.**

## Known limitations, stated plainly

- **Severely study-clustered.** Of the 213 geo-confirmed runs, 113 (53%) come
  from `PRJNA1309160` and 56 from `PRJNA1116203`. Treating them as independent
  observations would badly overstate the effective sample size — measured, the
  full panel's **effective *n* is 2.7**, not 213. Capping at 4 per BioProject
  raises it to **7.5** across the same 9 studies, which is what `nigeria-v1`
  does. Before using `nigeria-full` for any interval or significance claim:

  ```bash
  mycobench review-candidates --candidates cohorts/nigeria-full.tsv \
    --out-dir review --max-per-bioproject 4
  ```

  Nothing is deleted — capped rows move to `held_back.tsv` with their reason.

- **No phenotypes.** See above. This is the binding limitation.
- **No publication reviewed.** Every row is tier B. Promoting one to tier A
  means reviewing its source publication, recording it in `source_study`, and
  re-running `validate-cohort --online --write-lock`.
- **Geography is deposited metadata**, not curator-verified provenance. A
  `geo_loc_name` of `Nigeria` is the submitter's claim.
- **Sub-national detail is uneven.** Some rows carry a state (BENUE, FCT,
  PLATEAU); most carry only `Nigeria`. No within-country stratification is
  possible from this metadata.
- **Lineage is not recorded here.** TB-Profiler assigns it during the run; until
  then these panels make no claim about lineage balance, which matters because
  lineage confounds both resistance prediction and catalogue transfer.
- **The NTM panel is not a benchmark.** Those 12 runs exist to be refused. A
  drug call on any of them is a pipeline failure, not a result.

## Verification

Each panel ships with a lock recording the NCBI evidence retrieved at
verification time. Check a sheet against its lock before trusting it:

```bash
mycobench validate-cohort --samples cohorts/nigeria-v1.tsv \
  --verify-lock cohorts/nigeria-v1.lock.json
```

A checksum mismatch means the sheet changed after verification. Regenerate the
lock only after an intentional edit:

```bash
mycobench validate-cohort --samples cohorts/nigeria-v1.tsv --online \
  --write-lock cohorts/nigeria-v1.lock.json
```

Set `NCBI_EMAIL` and `NCBI_API_KEY` first; without a key NCBI allows 3
requests/second and a 213-row panel takes noticeably longer.

## Running a panel

```bash
mycobench run --cohort nigeria-v1                      # concordance track
mycobench run --cohort ntm-controls                    # species controls
mycobench run --cohort nigeria-v1 --write-script run_ng_v1.sh   # inspect first
```

Reads are downloaded on demand into the configured data directory and are not
committed: they are large and already publicly hosted by NCBI.
