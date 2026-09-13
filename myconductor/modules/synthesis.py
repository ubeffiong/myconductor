"""Regimen synthesis.

Reconciles all per-drug evidence into one verdict per drug, then proposes the
best available regimen. The reconciliation rule is conservative on purpose:
resistance evidence outweighs susceptibility, and a catalogued call always
beats a predicted one of equal strength.
"""
from __future__ import annotations

import json
from pathlib import Path

from ..core.models import Call, DrugEvidence, DrugResult, Regimen, Route

_DRUGS_PATH = Path(__file__).resolve().parent.parent / "catalogue" / "drugs.json"


def _rank(call: Call) -> int:
    # Higher = stronger evidence of resistance for reconciliation.
    return {
        Call.RESISTANT: 4,
        Call.PREDICTED_RESISTANT: 3,
        Call.INDETERMINATE: 2,
        Call.PREDICTED_SUSCEPTIBLE: 1,
        Call.SUSCEPTIBLE: 0,
    }[call]


def reconcile(evidence: list[DrugEvidence]) -> list[DrugResult]:
    by_drug: dict[str, list[DrugEvidence]] = {}
    for ev in evidence:
        by_drug.setdefault(ev.drug, []).append(ev)

    results: list[DrugResult] = []
    for drug, evs in sorted(by_drug.items()):
        # Prefer the strongest resistance signal; break ties toward catalogue.
        winner = max(
            evs,
            key=lambda e: (_rank(e.call), e.route == Route.CATALOGUE, e.confidence),
        )
        results.append(
            DrugResult(
                drug=drug, call=winner.call, confidence=winner.confidence, evidence=evs
            )
        )
    return results


class RegimenSynthesiser:
    def __init__(self, path: Path = _DRUGS_PATH):
        data = json.loads(Path(path).read_text())
        self.all_drugs: list[str] = data["drugs"]
        self.regimens: list[dict] = data["regimens"]

    def _effective(self, results_by_drug: dict[str, DrugResult], drug: str) -> bool:
        r = results_by_drug.get(drug)
        if r is None:
            return True  # no evidence of resistance -> presumed usable
        return not r.call.is_resistant and r.call != Call.INDETERMINATE

    def propose(self, results: list[DrugResult]) -> Regimen:
        by_drug = {r.drug: r for r in results}
        resistant = [r.drug for r in results if r.call.is_resistant]

        for regimen in self.regimens:
            drugs = regimen["drugs"]
            effective = [d for d in drugs if self._effective(by_drug, d)]
            if len(effective) >= regimen["min_effective"]:
                ineffective = [d for d in drugs if d not in effective]
                return Regimen(
                    proposed=effective,
                    effective_drugs=effective,
                    ineffective_drugs=ineffective,
                    adequate=True,
                    rationale=(
                        f"{regimen['name']}: {len(effective)}/{len(drugs)} companion "
                        f"drugs predicted effective (>= {regimen['min_effective']} required)."
                    ),
                )

        # Nothing adequate: report the largest salvage set we can muster.
        salvage = [d for d in self.all_drugs if self._effective(by_drug, d)]
        return Regimen(
            proposed=salvage,
            effective_drugs=salvage,
            ineffective_drugs=resistant,
            adequate=False,
            rationale=(
                "No standard regimen reaches its effective-drug floor. "
                f"Resistant to: {', '.join(resistant) or 'n/a'}. "
                "Escalating to the discovery loop."
            ),
        )
