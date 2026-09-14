# Where machine learning stops, and why

This page exists because the same question keeps arriving: *why is there no
trained VUS model?* It has been asked as a gap list, as an implementation plan,
and as an audit finding. Three times is enough to conclude the reasoning was
not findable, so it is written down here in one place.

The short version: **Myconductor governs model output. It does not produce
it.** That boundary is a design decision with four reasons behind it, and none
of them is caution.

## What exists

| Capability | Where | Status |
|---|---|---|
| A tier for model evidence | `core/models.py::Tier.PREDICTED` | built |
| Model output as evidence | `core/models.py::InSilicoPrediction` | built |
| Feature attributions (SHAP and equivalents) | `InSilicoPrediction.feature_attributions` | built |
| Imported feature vectors | `modules/variant_features.py` | built |
| Model registry with governance | `federated/model_registry.py` | built |
| Approval gated on beating the incumbent | `model_registry.baseline_verdict` | built |
| The incumbent, measured | `mycobench/analysis/baseline.py` | built |
| Structural import from AlphaFold DB | `mycobench/alphafold.py` | built |
| Calibration against site phenotypes | `modules/prediction_calibration.py` | built |
| Model signals beside a discordance | `modules/discordance_context.py` | built |
| **A trained model** | — | **deliberately absent** |
| **Feature computation (FoldX, SIFT, folding)** | — | **deliberately absent** |
| **Classifier adapters (XGBoost, HANN, ensemble)** | — | **deliberately absent** |

Audits often report the first block as missing because they search for names
from an earlier proposal — `PredictedEvidence`, `CalibrationStatus`,
`VUSClassifierAdapter`, a `vus/features/` package. Those names were never
adopted. The capabilities they describe were, under the names above.

## Why no model runs here

**1. A model cannot change a single call, however good it is.**

Model output enters at `Tier.PREDICTED`, and `Tier.may_establish_resistance` is
false for that tier — enforced in `DrugEvidence.__post_init__`, not by
convention. The most capable classifier imaginable can therefore only ever emit
`INDETERMINATE`: it can withhold susceptibility and raise a hypothesis, which
is what the rule-based lanes already do. An in-repo XGBoost would add a large
dependency tree and buy **zero additional clinical calls**.

**2. The measured baseline says where a model could help, and it is narrow.**

A registered model must beat the catalogue on both coverage and error to be
approved. Measured on 400 CRyPTIC isolates:

| Drug | Catalogue coverage | Error rate | Room for a model |
|---|---:|---:|---|
| levofloxacin | 98% | 1.8% | almost none |
| moxifloxacin | 98% | 4.1% | almost none |
| ethambutol | 99% | 7.9% | none on coverage |
| kanamycin | 89% | 5.6% | thin |
| amikacin | 88% | 7.9% | thin |
| **isoniazid** | **46%** | 3.3% | **54 points of abstention** |
| **rifampicin** | **40%** | 15.4% | **60 points of abstention** |
| bedaquiline, linezolid, delamanid, clofazimine, rifabutin | — | not estimable | unmeasurable |

Two things follow. The only real headroom is rifampicin and isoniazid, and it
is **abstention, not error** — and those abstentions are recorded as "loci not
fully examined". That is a sequencing-coverage problem. More depth, or a panel
that covers the loci, closes it; a classifier cannot, because it has no
information about a locus nobody sequenced.

And for five drugs the baseline cannot be measured at all at this cohort size,
so a model for them could not be approved even if it existed — there is nothing
to compare it against.

**3. The zero-dependency core is CI-enforced**, and is why this installs
anywhere, including the offline and low-resource settings it is aimed at.
XGBoost, SHAP, PyTorch and AlphaFold are gigabytes plus GPU and database
infrastructure. FoldX and Rosetta are not free for commercial use; AlphaFold 3's
weights carry their own restrictions. None of that is visible in a plan that
lists them as "adapters".

**4. Governing every model is worth more than shipping one.** A repo that
trains one classifier becomes a competitor to every other classifier. The
registry, the falsification gate and the tier system work for any model anyone
brings — which is the part nobody else does.

## The one thing that must never be built

**Promotion of a prediction to `SUSCEPTIBLE`.**

This is the centrepiece of most integration proposals, usually with the
argument that a low calibrated probability asserts the *absence* of resistance
and is therefore the safe direction. In this system that is backwards.

`SUSCEPTIBLE` is the only call that admits a drug to a regimen
(`modules/synthesis.py`), so it is the **actionable** one:

- a wrong `RESISTANT` costs a usable drug — bad, but conservative;
- a wrong `SUSCEPTIBLE` puts a patient on a failing regimen, which is how
  amplified resistance and onward transmission of a resistant strain are made.

`INDETERMINATE` already says "unknown" and withholds the drug safely. Promotion
converts "we don't know" into "prescribe this".

There is a structural objection too. `SUSCEPTIBLE` requires positive evidence
that the drug's loci were *callable*, sourced from outside the variant list. A
model probability is not coverage evidence: if a locus was never callable, the
model was scoring features of a variant nobody could observe. Promotion would
route around the callable mask, which is the gate the whole design rests on.

`modules/prediction_calibration.py::refuse_promotion()` is a named, tested
refusal rather than an absent feature, so anyone looking for "how do I promote a
calibrated prediction" finds this reasoning instead of an unguarded gap.

## What calibration is allowed to do

Raise laboratory priority, and nothing else. A calibrated high probability of
resistance is a good reason to move a variant up the validation queue; the lab
confirms it, and the confirmation enters as `PHENOTYPIC` evidence, which *can*
establish resistance. That is the legitimate path from a model to a call, and
it runs through a laboratory.

A *low* probability returns no signal at all — not a reason to deprioritise. A
model being unexcited about a variant is not evidence the variant is harmless,
and letting it push work down the queue would be the model quietly deciding
what never gets tested.

Calibration also refuses on the count of **resistant** observations rather than
the total. The conventional "at least 30 paired observations" fails silently
here: bedaquiline resistance runs near 0.8% in CRyPTIC, so thirty pairs contain
a quarter of one resistant isolate. A refusal says how many more are needed.

## How to bring a model

1. Train it elsewhere, on data you can describe.
2. Measure the incumbent on your cohort: `mycobench baseline`. The manifest
   emits registry-ready baseline blocks.
3. Register the model with a performance row per drug and lineage, each naming
   the baseline it was measured against.
4. Submit, review, approve. Approval **refuses** unless the model answers at
   least as often *and* errs strictly less often than that baseline — a tie
   earns nothing, since it buys an unvalidated dependency for accuracy already
   on hand.
5. Supply predictions as `InSilicoPrediction`, with feature attributions and
   the feature vector the model saw.

The output will be `Tier.PREDICTED`, it will not establish resistance, and it
will appear beside discordances and in the VUS queue. That is the ceiling, and
it is the ceiling by design.

## What would change this

A specific, falsifiable finding — not a better model in the abstract:

- **A measurable baseline for the rare-resistance drugs.** Bedaquiline,
  linezolid and delamanid are where VUS work matters most and where the
  catalogue cannot currently be measured. That needs far more of the compendium
  than 400 isolates.
- **Evidence that abstention is closable by inference rather than coverage.**
  Today rifampicin and isoniazid abstain because loci were not fully examined.
  If someone shows a model recovering those cases without more sequencing, the
  argument in §2 weakens.
- **A query that the flat evidence model cannot answer.** That is when a
  knowledge graph earns its infrastructure; until then `reporting/semantic.py`
  is the interoperability that was actually missing.
