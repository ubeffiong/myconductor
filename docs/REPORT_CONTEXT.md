# HTML report context dataset

The isolate analysis is the source of calls, variants, evidence, coverage, and
quality-control findings. Cohort-level claims require a separate JSON document
with the schema identifier `myconductor.report-context.v1`. Pass it to the same
report command with `--report-context`; Myconductor does not create another
report.

```console
myconductor analyze sample.vcf --mask callable.tsv \
  --html sample-report.html --report-context programme-context.json
```

Start with [`examples/report_context.empty.json`](../examples/report_context.empty.json).
Collections may be omitted. An omitted or empty collection produces a labelled
empty state in the report.

| Collection | Grain | Required fields | Useful optional fields |
|---|---|---|---|
| `baseline_rows` | one drug benchmark | `drug`, `estimable` | coverage, error rate, sensitivity, specificity, predictive values, counts and confidence intervals |
| `lineage_rows` | one drug and lineage | `drug`, `lineage` | called count, diagnostic metrics, major and very-major error rates, statistical power |
| `prevalence_rows` | one period, drug, lineage and geography | `period`, `drug`, `resistant`, `total` | lineage, geography, source, supplied prevalence |
| `target_rows` | one molecular target | `target` or `gene` | essentiality, druggability, human homology, resistance liability, evidence, source |
| `watchlist_rows` | one privacy-gated watch-list aggregate | `drug` | variant, unresolved isolate count, contributing site count, evidence gaps, limitations, source |
| `measurability` | one drug | `drug`, `measurable` | missing side and consequence |
| `audit_entries` | one ledger event | timestamp and action | actor, target and entry hash |
| `federated_sites` | one privacy-gated site aggregate | site and submission count | anonymity threshold and status |
| `error_trend` | ordered cohort windows | numeric percentage | none |

`validation` is an object with optional `outcomes` and `timeline` arrays.
`reference_method` names the phenotypic method used by the benchmark.

Prevalence is calculated only when resistant and tested counts are supplied.
Target evidence is descriptive and cannot change a drug call. Candidate
mechanisms, epistasis notes, watch-list aggregates and variants of uncertain
significance remain unable to establish resistance. The loader rejects unknown
fields and incorrectly typed collections so misspellings cannot silently empty
a panel.

For a populated, deterministic inspection fixture, run
`python -m myconductor demo --full-report artifacts/full-synthetic-demo`.
The generated context is explicitly synthetic and is loaded back through this
validator before the HTML is rendered. The fixture renders a 12-isolate cohort
and writes both primary-sample JSON and `myconductor-full-synthetic-demo.cohort.json`,
so it is useful for visual review, report payload inspection and end-to-end
contract testing.
