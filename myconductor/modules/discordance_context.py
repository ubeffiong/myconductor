"""Context beside a conflict, so a human can adjudicate it.

When two engines disagree about a drug, that disagreement is among the most
clinically useful things this platform reports, and it is deliberately left
unresolved: the call goes indeterminate and phenotypic testing is requested.
What has been missing is the material a clinician would want *while* reading
that conflict — what a model thought, what the structure shows — which today
sits elsewhere in the report, if it appears at all.

This module puts it next to the conflict. It does not resolve it.

Why this is display and nothing more
------------------------------------
The temptation with a tie-break signal is to let it break the tie. A model that
could convert a discordance into a call would be deciding clinical questions on
Tier.PREDICTED evidence through a side door, which is precisely the
black-box override the rest of this design refuses. So:

* ``Discordance.calls`` and ``Discordance.note`` are never modified;
* nothing here counts votes, weights sources or scores agreement;
* the context is phrased as observation, and every line names its source and
  says what it does not establish.

A reader can reasonably conclude "the model and the structure both point the
same way as TB-Profiler, so I will prioritise that phenotypic test". That is a
human making a judgement with more information, which is the goal. What they
cannot do is have the system make it for them.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Iterable, Optional, Sequence

from ..core.models import Discordance

#: Lines beyond this are dropped rather than flooding a report. Discordance is
#: meant to be read.
MAX_CONTEXT_LINES = 6

_NOT_ESTABLISHING = (
    "does not establish drug response; phenotypic testing still decides")


def _prediction_line(prediction, basis: Optional[dict] = None) -> str:
    confidence = ("" if prediction.confidence is None
                  else f" at confidence {prediction.confidence:.2f}")
    beaten = ""
    if basis and basis.get("baseline"):
        beaten = (f", approved against {basis['baseline']} on "
                  f"{basis.get('cohort', 'an evaluation cohort')}")
    top = sorted(prediction.feature_attributions.items(),
                 key=lambda kv: abs(kv[1]), reverse=True)[:3]
    drivers = ("; top features " + ", ".join(f"{n} {v:+.2f}" for n, v in top)
               if top else "")
    return (f"model {prediction.model_id} {prediction.model_version} predicts "
            f"{prediction.prediction}{confidence}{beaten}{drivers} — "
            f"{_NOT_ESTABLISHING}")


def _structural_line(annotation) -> str:
    distance = ""
    if annotation.ligand_distance is not None:
        distance = (f", {annotation.ligand_distance:g} "
                    f"{annotation.distance_unit} from "
                    f"{annotation.ligand_reference}")
    return (f"structure ({annotation.source} {annotation.source_version}) "
            f"reports {annotation.predicted_effect} at {annotation.location}"
            f"{distance} — {_NOT_ESTABLISHING}")


def context_for_drug(drug: str,
                     predictions: Sequence = (),
                     annotations: Sequence = (),
                     registry=None,
                     lineage: Optional[str] = None) -> tuple[str, ...]:
    """Observation lines relevant to one drug's conflict.

    A prediction is shown only if its model is currently approved. An
    unapproved model's output is not context, it is an unreviewed opinion, and
    putting it beside a clinical conflict would lend it standing it has not
    earned.
    """
    lines: list[str] = []

    for prediction in predictions:
        if prediction.drug != drug:
            continue
        basis = None
        if registry is not None:
            try:
                basis = registry.evidence_for(prediction, lineage)
            except ValueError:
                continue  # not approved for this scope: not shown at all
        lines.append(_prediction_line(prediction, basis))

    for annotation in annotations:
        lines.append(_structural_line(annotation))

    if len(lines) > MAX_CONTEXT_LINES:
        hidden = len(lines) - MAX_CONTEXT_LINES
        lines = lines[:MAX_CONTEXT_LINES]
        lines.append(f"({hidden} further signal(s) omitted; see the full "
                     f"findings)")
    return tuple(lines)


def annotate(discordances: Iterable[Discordance],
             predictions: Sequence = (),
             annotations_by_drug: Optional[dict] = None,
             registry=None,
             lineage: Optional[str] = None) -> list[Discordance]:
    """Attach context to each discordance, leaving the conflict itself intact.

    Returns new ``Discordance`` objects: the call tuple, the source tuple and
    the note are copied unchanged, and only ``context`` is populated.
    """
    annotations_by_drug = annotations_by_drug or {}
    annotated = []
    for discordance in discordances:
        context = context_for_drug(
            discordance.drug, predictions,
            annotations_by_drug.get(discordance.drug, ()), registry, lineage)
        annotated.append(replace(discordance, context=context)
                         if context else discordance)
    return annotated
