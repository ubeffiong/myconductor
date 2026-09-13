"""Import measured expression and reported interpretations; never infer drug response."""
import json
from dataclasses import asdict
from pathlib import Path
from ..core.models import Call, DrugEvidence, ExpressionEvidence, Lane, Tier


def load_expression(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("expression evidence must be a JSON list")
    return [ExpressionEvidence(**r) for r in data]


def reconcile_expression(results, expression_evidence, context):
    seen, findings = set(), []
    for observation in expression_evidence:
        if observation.observation_id in seen:
            raise ValueError("duplicate expression observation ID")
        seen.add(observation.observation_id)
        for key in ("sample_id", "isolate_id", "site_id", "organism"):
            if getattr(observation, key) != getattr(context, key):
                raise ValueError(f"expression {key} differs from current context")
        drugs = []
        for result in results:
            matched = [e for e in result.evidence if e.lane is Lane.EFFLUX_REGULATORY
                       and e.call is Call.INDETERMINATE and e.variant is not None
                       and (e.variant.gene == observation.gene or e.variant_key in observation.linked_variant_keys)]
            if not matched:
                continue
            drugs.append(result.drug)
            result.evidence.append(DrugEvidence(result.drug, Call.INDETERMINATE, Tier.PREDICTED, Lane.ENGINE,
                sample_id=context.sample_id, observation_id=observation.observation_id,
                rationale=f"Expression observation for {observation.gene}: reported {observation.reported_conclusion}. Drug response remains unestablished by expression.",
                metadata={"expression_observation": asdict(observation)},
                limitations=("Fold change alone does not establish overexpression, efflux activity or resistance.",)))
        findings.append(dict(asdict(observation), drugs=drugs, tier=Tier.PREDICTED.value,
                             matching_status="linked" if drugs else "unlinked", call_effect="none"))
    return findings
