"""Resistance-background stratification — separating causal from linked.

The problem this solves
-----------------------
A resistant isolate usually carries several variants: the determinant that
causes resistance, lineage markers it inherited, and hitchhikers in linkage
with the determinant. A marginal association — comparing carriers of a
candidate variant against everyone else — credits all of them equally, which
is how lineage markers end up in catalogues.

The fix is to compare like with like. For a candidate variant and a drug,
partition the isolates by their **resistance background**: the set of *already
graded* determinants for that drug which the isolate carries. Then ask whether
the candidate shifts MIC *within* a background, and pool those within-stratum
effects. A hitchhiker that only ever travels with ``rpoB S450L`` has no
variation left to explain once ``rpoB S450L`` is held fixed, and its pooled
effect collapses. A genuine determinant does not.

This is the logic behind WHO's solo-variant analysis, applied continuously: WHO
asks whether a variant appears alone and applies a threshold; this asks how
much it shifts MIC at a fixed background and reports an effect size with an
interval.

The stratification is only as good as the catalogue that defines "already
graded". With the real WHO catalogue that is 1,291 resistance-associated
entries across 15 drugs — which is why the background sets are small and most
isolates fall into the empty-background stratum, where the comparison is
cleanest.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional, Sequence

from . import stats

#: Catalogue calls that count as an established determinant for backgrounding.
#: 'indeterminate' deliberately does NOT: a variant of uncertain significance
#: is what we are trying to evaluate, so treating it as a known background
#: would assume the answer.
ESTABLISHED_CALLS = ("resistant",)


@dataclass
class DeterminantIndex:
    """Which variants are established determinants, per drug.

    Built from an ingested catalogue. Keyed on both the coordinate identity and
    the gene/label, because an isolate's genotype may be expressed either way
    and a background must not be missed because of a spelling difference.
    """

    by_drug_coordinate: dict[str, set[str]] = field(default_factory=dict)
    by_drug_label: dict[str, set[str]] = field(default_factory=dict)
    catalogue_version: str = ""
    n_entries: int = 0

    @classmethod
    def from_catalogue(cls, path: str | Path) -> "DeterminantIndex":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        index = cls(catalogue_version=data.get("catalogue_version", ""))
        for entry in data.get("variants", []):
            if entry.get("call") not in ESTABLISHED_CALLS:
                continue
            index.n_entries += 1
            label = f"{entry['gene']}_{entry['change']}"
            for drug in entry.get("drugs", []):
                index.by_drug_label.setdefault(drug, set()).add(label)
                for key in entry.get("coordinate_keys",
                                     [entry.get("coordinate_key")] or []):
                    if key:
                        index.by_drug_coordinate.setdefault(drug, set()).add(key)
        return index

    def determinants(self, drug: str) -> set[str]:
        return (self.by_drug_label.get(drug, set())
                | self.by_drug_coordinate.get(drug, set()))

    def background(self, drug: str,
                   genotype: Iterable[str]) -> frozenset[str]:
        """The established determinants for ``drug`` that this isolate carries."""
        known = self.determinants(drug)
        return frozenset(v for v in genotype if v in known)

    def drugs(self) -> list[str]:
        return sorted(set(self.by_drug_label) | set(self.by_drug_coordinate))


@dataclass
class Isolate:
    """One isolate's genotype, with optional lineage for confounding checks."""

    isolate_id: str
    genotype: frozenset[str] = frozenset()
    lineage: Optional[str] = None
    site: Optional[str] = None

    def carries(self, variant: str) -> bool:
        return variant in self.genotype


@dataclass
class Stratum:
    """One resistance background, and who does or does not carry the candidate."""

    background: frozenset[str]
    carriers: list[Isolate] = field(default_factory=list)
    non_carriers: list[Isolate] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.carriers) + len(self.non_carriers)

    @property
    def informative(self) -> bool:
        """A stratum with nobody on one side cannot contribute a comparison."""
        return bool(self.carriers) and bool(self.non_carriers)

    def label(self) -> str:
        if not self.background:
            return "(no established determinant)"
        return " + ".join(sorted(self.background))

    def lineage_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for isolate in self.carriers + self.non_carriers:
            key = isolate.lineage or "unknown"
            counts[key] = counts.get(key, 0) + 1
        return counts


@dataclass
class Stratification:
    """A candidate variant partitioned by resistance background for one drug."""

    variant: str
    drug: str
    strata: list[Stratum] = field(default_factory=list)
    n_carriers: int = 0
    n_non_carriers: int = 0

    @property
    def informative_strata(self) -> list[Stratum]:
        return [s for s in self.strata if s.informative]

    @property
    def n_informative(self) -> int:
        return len(self.informative_strata)

    @property
    def carriers_in_informative_strata(self) -> int:
        return sum(len(s.carriers) for s in self.informative_strata)

    @property
    def fully_confounded(self) -> bool:
        """Does the candidate never vary within any background?

        True means every carrier sits in a stratum with no non-carrier — the
        candidate is perfectly predicted by its background, so there is no
        within-background evidence about it at all. This is the hitchhiker
        signature, and it is a finding rather than a failure.
        """
        return self.n_carriers > 0 and self.n_informative == 0

    def lineage_diversity(self) -> float:
        """Effective number of lineages among carriers in informative strata."""
        counts: dict[str, int] = {}
        for stratum in self.informative_strata:
            for isolate in stratum.carriers:
                key = isolate.lineage or "unknown"
                counts[key] = counts.get(key, 0) + 1
        return stats.inverse_simpson(list(counts.values()))

    def describe(self) -> str:
        if self.fully_confounded:
            return (f"{self.variant}/{self.drug}: {self.n_carriers} carrier(s), "
                    f"none in a stratum containing a non-carrier — the variant "
                    f"never varies at fixed background")
        return (f"{self.variant}/{self.drug}: {self.n_carriers} carrier(s) "
                f"across {self.n_informative} informative stratum/strata "
                f"({self.carriers_in_informative_strata} usable)")


def stratify(variant: str, drug: str, isolates: Sequence[Isolate],
             index: DeterminantIndex,
             exclude_self: bool = True) -> Stratification:
    """Partition ``isolates`` by resistance background for one candidate.

    ``exclude_self`` removes the candidate from its own background, which
    matters when the candidate is itself already graded: without it, every
    carrier would land in a different stratum from every non-carrier by
    construction and nothing would be comparable.
    """
    buckets: dict[frozenset[str], Stratum] = {}
    result = Stratification(variant=variant, drug=drug)

    for isolate in isolates:
        background = index.background(drug, isolate.genotype)
        if exclude_self:
            background = frozenset(background - {variant})
        stratum = buckets.get(background)
        if stratum is None:
            stratum = Stratum(background=background)
            buckets[background] = stratum
        if isolate.carries(variant):
            stratum.carriers.append(isolate)
            result.n_carriers += 1
        else:
            stratum.non_carriers.append(isolate)
            result.n_non_carriers += 1

    result.strata = sorted(buckets.values(), key=lambda s: (-s.n, s.label()))
    return result


def cooccurrence(variant: str, drug: str, isolates: Sequence[Isolate],
                 index: DeterminantIndex) -> dict[str, float]:
    """Fraction of the candidate's carriers that also carry each determinant.

    A candidate co-occurring with an established determinant in nearly all its
    carriers cannot be credited with the phenotype, and this is the number that
    says so. Complements the stratification: ``fully_confounded`` is the
    categorical version of the same fact.
    """
    carriers = [i for i in isolates if i.carries(variant)]
    if not carriers:
        return {}
    known = index.determinants(drug)
    counts: dict[str, int] = {}
    for isolate in carriers:
        for other in isolate.genotype & known:
            if other == variant:
                continue
            counts[other] = counts.get(other, 0) + 1
    return {determinant: count / len(carriers)
            for determinant, count in sorted(counts.items(),
                                             key=lambda kv: -kv[1])}
