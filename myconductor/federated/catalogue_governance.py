"""Catalogue-version-aware multi-site reconciliation.

The problem
-----------
``federated/catalogue_update.py`` already governs how new evidence is learned
and promoted *within* a coordinated federation. It does not address a
different, more mundane problem: two sites running Myconductor against the
*same* variant can legitimately disagree because they are running different
WHO catalogue versions — v2 may grade a variant ``indeterminate`` while v3
grades the same coordinate ``resistant``. That disagreement is not a bug in
either site's pipeline, and it is a different problem from two sites running
the *same* catalogue version and still disagreeing (which usually is a
pipeline or rule bug). Conflating the two wastes the time of whoever has to
investigate, and "fragmented data-sharing frameworks" plus "inadequate
connectivity for timely database updates" are exactly the reviewed failure
modes this is meant to make legible rather than solve by fiat.

What this module does
----------------------
Treats catalogue version as a first-class coordinate. Every ``DrugEvidence``
already carries ``catalogue_version`` and ``rule_id`` (see
``core/models.py`` and ``modules/catalogue.py``/``modules/local_validation.py``,
which populate them). This module extracts one ``SiteCallRecord`` per
variant/drug pair a site reported on, and classifies cross-site disagreement
into exactly two categories:

``CATALOGUE_VERSION``
    the disagreeing sites ran different catalogue versions — the fix is
    catalogue propagation/update, not debugging.
``RULE``
    the disagreeing sites ran the *same* catalogue version (or neither
    reported one) and still produced different calls — the fix is in the
    pipeline or the rule, not the catalogue.

Like every other discordance in this codebase, this never resolves the
conflict by outvoting a minority site. It is a report, not an adjudication.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Optional

from ..core.models import AnalysisReport


@dataclass(frozen=True)
class SiteCallRecord:
    """One site's call on one variant/drug pair, with full provenance."""

    site_id: str
    isolate_id: Optional[str]
    variant_key: str
    gene: str
    drug: str
    catalogue_version: Optional[str]
    rule_id: Optional[str]
    call: str
    tier: str
    timestamp_utc: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "site_id": self.site_id, "isolate_id": self.isolate_id,
            "variant_key": self.variant_key, "gene": self.gene,
            "drug": self.drug, "catalogue_version": self.catalogue_version,
            "rule_id": self.rule_id, "call": self.call, "tier": self.tier,
            "timestamp_utc": self.timestamp_utc,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "SiteCallRecord":
        return cls(
            site_id=data["site_id"], isolate_id=data.get("isolate_id"),
            variant_key=data["variant_key"], gene=data.get("gene", ""),
            drug=data["drug"], catalogue_version=data.get("catalogue_version"),
            rule_id=data.get("rule_id"), call=data["call"],
            tier=data.get("tier", ""), timestamp_utc=data.get("timestamp_utc"),
        )


def site_call_records(report: AnalysisReport, site_id: str,
                      isolate_id: Optional[str] = None
                      ) -> list[SiteCallRecord]:
    """Extract one record per variant/drug pair with evidence in a report.

    A ``DrugResult`` with no variant-attached evidence (e.g. a clean,
    coverage-gated susceptible with no lane speaking) contributes nothing —
    there is no variant identity to reconcile across sites.
    """
    isolate_id = isolate_id or report.sample_id
    out: list[SiteCallRecord] = []
    seen: set[tuple[str, str]] = set()
    for result in report.drug_results:
        for ev in result.evidence:
            if ev.variant_key is None:
                continue
            key = (ev.variant_key, ev.drug)
            if key in seen:
                continue
            seen.add(key)
            out.append(SiteCallRecord(
                site_id=site_id, isolate_id=isolate_id,
                variant_key=ev.variant_key, gene=ev.variant.gene if ev.variant else "",
                drug=ev.drug, catalogue_version=ev.catalogue_version,
                rule_id=ev.rule_id, call=ev.call.value, tier=ev.tier.value,
                timestamp_utc=report.provenance.generated_utc,
            ))
    return out


class DiscordanceCategory(str, Enum):
    CATALOGUE_VERSION = "catalogue_version"
    RULE = "rule"


@dataclass
class CrossSiteDiscordance:
    variant_key: str
    drug: str
    category: DiscordanceCategory
    calls_by_site: dict[str, str]
    catalogue_versions_by_site: dict[str, Optional[str]]
    note: str


def reconcile_sites(records: Iterable[SiteCallRecord]
                    ) -> list[CrossSiteDiscordance]:
    """Classify every cross-site disagreement by its actual cause.

    Groups records by ``(variant_key, drug)``. Sites that agree are silent —
    this only reports where sites differ. Within a disagreeing group: if the
    calls that differ come from different catalogue versions, the category is
    ``CATALOGUE_VERSION``; if the differing calls share a catalogue version
    (or neither carries one), the category is ``RULE``. A group with both
    kinds of disagreement gets a note listing both rather than collapsing to
    one — the point is not to force a single label onto a mixed situation.
    """
    grouped: dict[tuple[str, str], list[SiteCallRecord]] = {}
    for record in records:
        grouped.setdefault((record.variant_key, record.drug), []).append(record)

    out: list[CrossSiteDiscordance] = []
    for (variant_key, drug), group in sorted(grouped.items()):
        calls_by_site = {r.site_id: r.call for r in group}
        if len(set(calls_by_site.values())) <= 1:
            continue  # every site agrees; nothing to report

        versions_by_site = {r.site_id: r.catalogue_version for r in group}
        distinct_calls = sorted(set(calls_by_site.values()))

        # For each pair of distinct calls, does a version difference explain
        # it, a rule difference, or both?
        version_explains = False
        rule_explains = False
        for i, call_a in enumerate(distinct_calls):
            for call_b in distinct_calls[i + 1:]:
                sites_a = [s for s, c in calls_by_site.items() if c == call_a]
                sites_b = [s for s, c in calls_by_site.items() if c == call_b]
                versions_a = {versions_by_site[s] for s in sites_a}
                versions_b = {versions_by_site[s] for s in sites_b}
                if versions_a.isdisjoint(versions_b):
                    version_explains = True
                else:
                    rule_explains = True

        if version_explains and not rule_explains:
            category = DiscordanceCategory.CATALOGUE_VERSION
            version_parts = [f"{s}={v or 'unknown'}"
                             for s, v in sorted(versions_by_site.items())]
            note = (
                f"disagreement tracks catalogue version "
                f"({', '.join(version_parts)}); resolve by propagating the "
                f"newer catalogue, not by debugging either site's pipeline"
            )
        elif rule_explains and not version_explains:
            category = DiscordanceCategory.RULE
            note = (
                f"sites report the same catalogue version(s) but different "
                f"calls ({', '.join(f'{s}={c}' for s, c in sorted(calls_by_site.items()))}); "
                f"this is a pipeline/rule discordance, not a catalogue-"
                f"propagation gap"
            )
        else:
            category = DiscordanceCategory.CATALOGUE_VERSION
            note = (
                "disagreement is not explained by catalogue version alone: "
                "some sites disagree despite sharing a version, and others "
                "differ by version — both a catalogue-propagation gap AND a "
                "rule discordance may be present; investigate both"
            )

        out.append(CrossSiteDiscordance(
            variant_key=variant_key, drug=drug, category=category,
            calls_by_site=calls_by_site,
            catalogue_versions_by_site=versions_by_site, note=note,
        ))
    return out
