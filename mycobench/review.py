"""Study-dependence review: stop a clustered panel posing as independent data.

Of the 213 geo-confirmed Nigerian runs, 113 come from one BioProject and 56
from another. Treating those as 213 independent observations would badly
overstate the effective sample size behind any interval or significance claim —
isolates from one study share a sampling frame, a laboratory, a sequencing run
and often a transmission chain.

This module measures that clustering and produces a capped shortlist. It
deletes nothing: the full candidate table is always retained, because a cap is
a statement about independence, not about data quality.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Optional

from .cohort import read_rows


@dataclass
class StudyDependence:
    """How concentrated a candidate table is, per grouping field."""

    field_name: str
    counts: dict[str, int] = field(default_factory=dict)

    @property
    def n_groups(self) -> int:
        return len(self.counts)

    @property
    def n_rows(self) -> int:
        return sum(self.counts.values())

    @property
    def largest(self) -> tuple[str, int]:
        if not self.counts:
            return ("", 0)
        return max(self.counts.items(), key=lambda kv: kv[1])

    @property
    def largest_share(self) -> float:
        return (self.largest[1] / self.n_rows) if self.n_rows else 0.0

    @property
    def effective_n(self) -> float:
        """Inverse Simpson index over group sizes.

        A crude but honest summary: 213 rows spread evenly over 9 studies has an
        effective n near 9, not 213. It is a descriptive statistic, not a
        correction that licenses treating the rows as independent.
        """
        if not self.n_rows:
            return 0.0
        total = self.n_rows
        return 1.0 / sum((count / total) ** 2 for count in self.counts.values())

    def describe(self) -> str:
        group, size = self.largest
        return (f"{self.n_rows} row(s) across {self.n_groups} "
                f"{self.field_name}(s); largest is {group} with {size} "
                f"({self.largest_share:.0%}); effective n "
                f"{self.effective_n:.1f}")


def study_dependence(rows: Iterable[dict],
                     field_name: str = "bioproject") -> StudyDependence:
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        counts[(row.get(field_name) or "").strip() or "(unset)"] += 1
    return StudyDependence(field_name=field_name, counts=dict(counts))


@dataclass
class ReviewResult:
    shortlist: list[dict] = field(default_factory=list)
    held_back: list[dict] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)
    dependence: dict[str, StudyDependence] = field(default_factory=dict)

    @property
    def n_shortlist(self) -> int:
        return len(self.shortlist)

    def summary(self) -> str:
        lines = [f"shortlist {self.n_shortlist} row(s), "
                 f"{len(self.held_back)} held back"]
        for dependence in self.dependence.values():
            lines.append("  before: " + dependence.describe())
        after = study_dependence(self.shortlist)
        lines.append("  after : " + after.describe())
        return "\n".join(lines)


def _sort_key(row: dict) -> tuple:
    """Prefer deeper runs, then a stable id, so a cap is deterministic."""
    try:
        bases = int(row.get("bases") or 0)
    except ValueError:
        bases = 0
    return (-bases, (row.get("sample_id") or ""))


def review(rows: list[dict], max_per_bioproject: Optional[int] = None,
           max_per_organism: Optional[int] = None,
           max_total: Optional[int] = None) -> ReviewResult:
    """Cap a candidate table without deleting anything.

    Rows that exceed a cap move to ``held_back`` with the reason recorded. The
    caller writes both tables.
    """
    result = ReviewResult()
    result.dependence = {
        "bioproject": study_dependence(rows, "bioproject"),
        "organism": study_dependence(rows, "organism"),
    }

    per_project: dict[str, int] = defaultdict(int)
    per_organism: dict[str, int] = defaultdict(int)

    for row in sorted(rows, key=_sort_key):
        project = (row.get("bioproject") or "").strip()
        organism = (row.get("organism") or "").strip()

        if max_per_bioproject is not None and \
                per_project[project] >= max_per_bioproject:
            result.held_back.append(row)
            result.reasons.append(
                f"{row.get('sample_id')}: BioProject {project} already has "
                f"{max_per_bioproject} row(s) in the shortlist")
            continue
        if max_per_organism is not None and \
                per_organism[organism] >= max_per_organism:
            result.held_back.append(row)
            result.reasons.append(
                f"{row.get('sample_id')}: organism {organism} already has "
                f"{max_per_organism} row(s) in the shortlist")
            continue
        if max_total is not None and result.n_shortlist >= max_total:
            result.held_back.append(row)
            result.reasons.append(
                f"{row.get('sample_id')}: shortlist already has "
                f"{max_total} row(s)")
            continue

        result.shortlist.append(row)
        per_project[project] += 1
        per_organism[organism] += 1

    return result


def review_file(samples: str, max_per_bioproject: Optional[int] = None,
                max_per_organism: Optional[int] = None,
                max_total: Optional[int] = None) -> ReviewResult:
    rows, _ = read_rows(samples)
    return review(rows, max_per_bioproject, max_per_organism, max_total)


def dependence_rows(dependence: StudyDependence) -> list[dict]:
    """Long-form rows for the report and a written TSV."""
    total = dependence.n_rows or 1
    return [{"field": dependence.field_name, "group": group, "rows": count,
             "share": f"{count / total:.4f}"}
            for group, count in sorted(dependence.counts.items(),
                                       key=lambda kv: -kv[1])]


DEPENDENCE_COLUMNS = ("field", "group", "rows", "share")
