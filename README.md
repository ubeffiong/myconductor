# Myconductor

**An orchestration-and-reasoning layer for tuberculosis antimicrobial resistance.**

*myco*bacterium + *conductor* — the layer that conducts the existing tools.

Today's TB-AMR tools (TB-Profiler, Mykrobe, MAST, tNGS pipelines) are all the
same kind of thing: a static, version-locked **catalogue lookup**. Sequence goes
in, each drug comes out susceptible / resistant / indeterminate, and when the
strain throws a variant of unknown significance (VUS) or defeats the whole drug
arsenal, the pipeline stops. Myconductor is **not another profiler** — it is the
layer that sits above them and makes them act as one adaptive system:

- a **triage router** that sends every variant to the right specialist lane
  (catalogue lookup, VUS ML classifier, or efflux/regulatory model);
- **heteroresistance detection** that surfaces low-frequency resistant
  subpopulations a consensus caller misses;
- **regimen synthesis** into BPaLM/BPaL/first-line;
- a **closed discovery loop** — when no regimen works, it automatically seeds a
  subtractive-genomics + docking pipeline against that strain's targets;
- a **federated living catalogue** that learns from new data without moving
  patient genomes off-site.

The core has **zero external dependencies** (Python 3.9+ standard library only),
so it runs offline — which is also the point for low-resource deployment.

## Install

```bash
pip install -e .          # or just run from the repo root with no install
```

## Quickstart

```bash
python -m myconductor demo                       # bundled end-to-end demo
python -m myconductor analyze path/to/input.vcf  # analyze your own variants
python -m myconductor analyze input.vcf --fhir out.json   # + FHIR bundle
```

Programmatic:

```python
from myconductor import Myconductor
from myconductor.reporting.render import render_text

report = Myconductor().analyze("myconductor/data/example_input.vcf")
print(render_text(report))
```

Input can be `.vcf`, `.tsv`, or `.json` (see `myconductor/data/` for examples and
`myconductor/io/adapter.py` for the field contract).

## What the demo shows

The bundled sample (`MTB-DEMO-001`) is a pre-XDR-style strain. Running it:

- routes 4 variants to the catalogue, 1 to the VUS classifier, 2 to the
  efflux/regulatory lane, and ignores a synonymous change;
- flags a **minority (15%) fluoroquinolone-resistant subpopulation** a consensus
  caller would miss;
- **drops a low-depth locus** at QC rather than mis-calling it susceptible;
- finds **no adequate regimen** and **triggers the discovery loop**, ranking
  novel targets by predicted docking affinity;
- separates catalogued calls (`R`/`S`) from model predictions (`R*`/`S*`)
  everywhere, including in the FHIR output.

## Tests

```bash
python -m unittest discover -s tests -v
```

## Layout

```
myconductor/
  core/       models, router (the brain), pipeline (the conductor)
  io/         input adapter + QC gating
  modules/    catalogue, vus_classifier, efflux, heteroresistance, synthesis, discovery
  reporting/  human-readable render + FHIR bundle
  federated/  privacy-preserving living-catalogue feedback + drift monitor
  catalogue/  illustrative WHO-style variant catalogue + drug regimens
  data/       example inputs
docs/ARCHITECTURE.md    design, module map, extension points
tests/                  unittest suite (no pytest required)
examples/run_demo.py    programmatic example
```

## Scope and safety

Myconductor is a **research and decision-support scaffold, not a clinical
device.** The bundled catalogue is a small illustrative subset (a real deployment
loads the full WHO catalogue of ~11,390 graded variants); the VUS scorer and
docking scores are transparent, deterministic placeholders that exercise the
model interfaces — they are **not trained or validated predictors**. Every report
is explicit that predicted calls are preliminary and that discovery-loop targets
are hypothesis-generating only. Turning this into a product is gated by
training-data breadth, clinical validation, and regulatory approval — not by the
code. See `docs/ARCHITECTURE.md`.

## License

MIT — see `LICENSE`.
