# Myconductor architecture

Myconductor is not another variant profiler. It is the **orchestration-and-reasoning
layer** that sits above the existing tools and makes them behave as one adaptive
system. Today's tools (TB-Profiler, Mykrobe, MAST, tNGS pipelines) are static,
version-locked catalogue lookups: sequence in, per-drug S/R/indeterminate out, and
when a strain throws a variant of unknown significance (VUS) or defeats the whole
arsenal, the pipeline stops. Myconductor fills the gaps *between* those tools.

## Data flow

```
 adapter → known-variant → triage router →  ┬ VUS classifier  ┐
 (any input)  (catalogue)   (the brain)     └ efflux/regulatory┘
                                               │
                                     heteroresistance scan
                                               │
                                       regimen synthesis
                                               │
                                   effective regimen?
                                        │            │
                                    yes │            │ no
                                        ▼            ▼
                                clinical report   discovery loop
                                   + FHIR         (subtractive
                                        │          genomics + docking)
                                        └──────┬──────┘
                                               ▼
                                  federated living catalogue ↻
                                   (retrains the router)
```

## Modules

| Stage | File | What it does | Real-world backend to plug in |
|---|---|---|---|
| Adapter | `io/adapter.py` | Normalises VCF/TSV/JSON to canonical `Variant`s; coverage-gated QC | any variant caller |
| Known-variant | `modules/catalogue.py` | Catalogue lookup with WHO grades | full WHO catalogue / TB-Profiler output |
| **Triage router** | `core/router.py` | Dispatches each variant to the right lane | — (this is the novel core) |
| VUS classifier | `modules/vus_classifier.py` | Predicts resistance for unknown coding variants, with SHAP-style explanation | trained XGBoost / HANN |
| Feature engineering | `modules/features.py` | Turns a variant into biological features | SIFT/PolyPhen + AlphaFold |
| Efflux / regulatory | `modules/efflux.py` | Flags efflux overexpression and promoter effects | expression-inference model |
| Heteroresistance | `modules/heteroresistance.py` | Surfaces minority resistant subpopulations from VAF | low-frequency variant caller |
| Synthesis | `modules/synthesis.py` | Reconciles evidence, proposes a regimen | clinical guideline engine |
| Discovery | `modules/discovery.py` | Subtractive-genomics + docking when no regimen works | BLASTP + DEG + AutoDock Vina |
| Reporting | `reporting/` | Human-readable text + FHIR bundle | LIMS/EHR |
| Federated | `federated/catalogue_update.py` | Privacy-preserving learning + drift monitor | secure aggregation transport |

## Three design commitments

1. **Catalogued ≠ predicted.** The `Call` enum keeps graded catalogue calls and
   model predictions in separate states, marked `R`/`S` vs `R*`/`S*` everywhere,
   and FHIR observations from predictions are `preliminary`. Human-in-the-loop by
   construction.
2. **Closed diagnosis→discovery loop.** A failed regimen automatically seeds the
   discovery pipeline against *that strain's* targets. Diagnosis and discovery are
   one workflow.
3. **Learn without moving genomes.** Sites share only aggregated, de-identified
   variant→phenotype tallies. Enough concordant evidence promotes a variant into
   the catalogue; a drift monitor guards against silent model decay.

## Extending

Every lane implements the `VariantModule` interface (`modules/base.py`): accept a
`Variant`, return a `DrugEvidence` or `None`. Swapping the demo scorer for a
trained model, or the mini catalogue for the full WHO catalogue, touches nothing
else. The pipeline collaborators are all constructor-injected.

## Scope / honesty

This is a **research and decision-support scaffold**, not a clinical device. The
bundled catalogue is a tiny illustrative subset; the VUS scorer and docking scores
are transparent deterministic placeholders exercising the interfaces, not trained
or validated predictors. The hard part of turning this into a product is not the
code — it is training-data breadth (especially for bedaquiline/pretomanid),
clinical validation, and regulatory approval.
```
