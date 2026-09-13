"""Discovery loop.

Fires only when synthesis fails to find an adequate regimen -- the closed loop
where a diagnostic dead-end automatically seeds drug discovery. It runs an
in-silico subtractive-genomics funnel and a docking stub:

    proteome -> drop human homologs -> keep essential -> keep druggable -> dock

Each stage below is a labelled hook. In production:
    * homology       -> local BLASTP against the human proteome
    * essentiality   -> Database of Essential Genes (DEG) cross-reference
    * druggability   -> pocket detection over PDB / AlphaFold structures
    * docking        -> AutoDock Vina against a compound library

The demo ships a tiny illustrative proteome so the funnel produces real output
offline. The docking "scores" are placeholders, clearly flagged.
"""
from __future__ import annotations

from ..core.models import DiscoveryResult, DiscoveryTarget
from .features import _stable_unit

# Illustrative validated / candidate Mtb targets with annotation flags.
# essential/human_homolog/druggable would be computed, not hand-set, in prod.
_DEMO_PROTEOME = [
    {"protein": "DNA gyrase subunit B", "gene": "gyrB", "essential": True, "human_homolog": False, "druggable": True},
    {"protein": "Enoyl-ACP reductase", "gene": "inhA", "essential": True, "human_homolog": False, "druggable": True},
    {"protein": "Decaprenylphosphoryl-ribose oxidase", "gene": "dprE1", "essential": True, "human_homolog": False, "druggable": True},
    {"protein": "ATP synthase c-ring", "gene": "atpE", "essential": True, "human_homolog": False, "druggable": True},
    {"protein": "Peptidoglycan transpeptidase", "gene": "ponA1", "essential": True, "human_homolog": False, "druggable": True},
    {"protein": "Glutamate dehydrogenase", "gene": "gdh", "essential": False, "human_homolog": True, "druggable": True},
    {"protein": "Hypothetical conserved protein", "gene": "Rv1234", "essential": True, "human_homolog": False, "druggable": False},
]

_COMPOUNDS = ["repurpose:vancomycin", "repurpose:meropenem", "natural:berberine",
              "repurpose:fosfomycin", "natural:plumbagin"]


def _dock(gene: str) -> tuple[float, str]:
    """Placeholder docking. Returns (kcal/mol, best compound).

    Real implementation: AutoDock Vina against the target pocket. More negative
    is a stronger predicted binding affinity.
    """
    u = _stable_unit("dock:" + gene)
    score = round(-11.0 + 4.0 * u, 2)  # roughly -11 .. -7 kcal/mol
    compound = _COMPOUNDS[int(_stable_unit("cmpd:" + gene) * len(_COMPOUNDS)) % len(_COMPOUNDS)]
    return score, compound


def run(reason: str, triggered: bool) -> DiscoveryResult:
    if not triggered:
        return DiscoveryResult(triggered=False, reason=reason, targets=[])

    targets: list[DiscoveryTarget] = []
    for p in _DEMO_PROTEOME:
        # Subtractive funnel: must be essential, non-human, druggable.
        passes = p["essential"] and not p["human_homolog"] and p["druggable"]
        dock_score, compound = _dock(p["gene"]) if passes else (None, None)
        targets.append(
            DiscoveryTarget(
                protein=p["protein"],
                gene=p["gene"],
                essential=p["essential"],
                human_homolog=p["human_homolog"],
                druggable=p["druggable"],
                best_docking_kcal_mol=dock_score,
                top_compound=compound,
            )
        )
    # Rank surviving targets by predicted affinity.
    targets.sort(key=lambda t: (t.best_docking_kcal_mol is None,
                                t.best_docking_kcal_mol if t.best_docking_kcal_mol is not None else 0.0))
    return DiscoveryResult(triggered=True, reason=reason, targets=targets)
