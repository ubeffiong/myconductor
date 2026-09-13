"""Discover cohort candidates from NCBI, confirming origin rather than assuming it.

This is a discovery and evidence-producing tool, not a claim that every public
record is what it appears to be. It accepts a run only when:

* the platform is Illumina and the library is paired-end;
* the run carries a BioSample (without one, nothing can be confirmed);
* the organism matches the requested taxon; and
* where a country is requested, the BioSample's controlled ``geo_loc_name``
  field **confirms** it.

That last rule is the one worth spelling out. Searching SRA for
``"Mycobacterium tuberculosis complex"[Organism] AND Nigeria[All Fields]``
returns 472 runs, of which 37 have ``geo_loc_name = "China: East China Sea,
Xiangshan Bay"`` — the word matched somewhere else in the record. A further 134
BioSamples carry no ``geo_loc_name`` at all and are also rejected, because
unconfirmed is not confirmed. Accepting either group would have put marine
metagenomic samples into a Nigerian clinical cohort with nothing downstream
noticing.

Every rejection is written out with its reason and the evidence that produced
it, so the filter is auditable rather than a black box.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Iterable, Optional

from .cohort import is_mtbc
from .ncbi import (
    Credentials,
    NCBIError,
    bases_of,
    biosample_attributes,
    esearch,
    paired_illumina,
    runinfo,
)

ACCEPTED_COLUMNS = ("sample_id", "sra_run", "biosample", "bioproject",
                    "organism", "geo_loc_name", "collection_date",
                    "isolation_source", "platform", "layout", "bases",
                    "cohort_tier", "phenotype_source", "source_study",
                    "expected_outcome")
REJECTED_COLUMNS = ("sra_run", "biosample", "bioproject", "organism_query",
                    "country_query", "organism_evidence", "geo_evidence",
                    "reason")

PENDING_REVIEW = "NCBI_discovery_pending_publication_review"

#: Excluded from MTBC panels: an animal-adapted lineage is a different
#: epidemiological object, and mixing it into a human clinical cohort silently
#: changes what the cohort represents.
DEFAULT_EXCLUDE = ("variant bovis",)


@dataclass
class DiscoveryResult:
    accepted: list[dict] = field(default_factory=list)
    rejected: list[dict] = field(default_factory=list)
    searched: int = 0
    geo_distribution: dict[str, int] = field(default_factory=dict)

    @property
    def n_accepted(self) -> int:
        return len(self.accepted)

    def summary(self) -> str:
        lines = [f"searched {self.searched} run(s); accepted "
                 f"{self.n_accepted}, rejected {len(self.rejected)}"]
        reasons: dict[str, int] = {}
        for row in self.rejected:
            reasons[row["reason"]] = reasons.get(row["reason"], 0) + 1
        for reason, count in sorted(reasons.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {count:>5}  {reason}")
        if self.geo_distribution:
            lines.append("  geo_loc_name observed:")
            for value, count in sorted(self.geo_distribution.items(),
                                       key=lambda kv: -kv[1])[:10]:
                lines.append(f"  {count:>5}  {value or '(absent)'!r}")
        return "\n".join(lines)


def _sample_id(prefix: str, run: str) -> str:
    """Traceable by construction: the run accession is the identity."""
    return f"{prefix}_{run}" if prefix else run


def discover(organisms: Iterable[str], creds: Credentials,
             country: Optional[str] = None,
             max_runs: int = 600,
             prefix: str = "",
             exclude_organisms: tuple[str, ...] = DEFAULT_EXCLUDE,
             require_mtbc: bool = True,
             expected_outcome: str = "",
             min_bases: int = 50_000_000) -> DiscoveryResult:
    """Search SRA for candidate runs and confirm each one's deposited evidence."""
    result = DiscoveryResult()
    seen_runs: set[str] = set()

    for organism in organisms:
        term = f'"{organism}"[Organism]'
        if country:
            # The free-text clause only NARROWS the search; it never confirms
            # origin. Confirmation happens against geo_loc_name below.
            term += f" AND {country}[All Fields]"
        search = esearch("sra", term, creds, retmax=max_runs)
        time.sleep(creds.interval)
        if search["count"] == 0:
            continue
        rows = runinfo(search, creds, retmax=max_runs)
        time.sleep(creds.interval)
        result.searched += len(rows)

        candidates = paired_illumina(rows, exclude_organisms=exclude_organisms)
        rejected_here = {row["Run"] for row in rows} - {
            row["Run"] for row in candidates}
        for row in rows:
            if row["Run"] not in rejected_here:
                continue
            reason = "not paired-end Illumina with a BioSample"
            organism_name = row.get("ScientificName") or ""
            if any(token.lower() in organism_name.lower()
                   for token in exclude_organisms):
                reason = f"organism excluded by request ({organism_name})"
            result.rejected.append({
                "sra_run": row.get("Run", ""),
                "biosample": row.get("BioSample", ""),
                "bioproject": row.get("BioProject", ""),
                "organism_query": organism, "country_query": country or "",
                "organism_evidence": organism_name, "geo_evidence": "",
                "reason": reason})

        samples = biosample_attributes(
            [row["BioSample"] for row in candidates], creds)

        for row in candidates:
            run = row["Run"]
            if run in seen_runs:
                continue
            seen_runs.add(run)

            organism_name = row.get("ScientificName") or ""
            record = samples.get(row["BioSample"])
            geo = record.geo_loc_name if record else ""
            result.geo_distribution[geo] = result.geo_distribution.get(geo, 0) + 1

            def reject(reason: str) -> None:
                result.rejected.append({
                    "sra_run": run, "biosample": row.get("BioSample", ""),
                    "bioproject": row.get("BioProject", ""),
                    "organism_query": organism, "country_query": country or "",
                    "organism_evidence": organism_name, "geo_evidence": geo,
                    "reason": reason})

            if record is None:
                reject("BioSample record could not be retrieved")
                continue
            if require_mtbc and not is_mtbc(organism_name):
                reject(f"organism {organism_name!r} is not MTBC")
                continue
            if not require_mtbc and is_mtbc(organism_name):
                reject(f"organism {organism_name!r} is MTBC but a non-MTBC "
                       f"control panel was requested")
                continue
            if country and not record.confirms_country(country):
                reject(f"BioSample geo_loc_name {geo!r} does not confirm "
                       f"{country}" if geo else
                       f"BioSample carries no geo_loc_name, so {country} "
                       f"cannot be confirmed")
                continue
            if bases_of(row) < min_bases:
                reject(f"{bases_of(row):,} bases is below the {min_bases:,} "
                       f"minimum for a usable profile")
                continue

            phenotype = record.phenotype_attributes()
            result.accepted.append({
                "sample_id": _sample_id(prefix, run),
                "sra_run": run,
                "biosample": row.get("BioSample", ""),
                "bioproject": row.get("BioProject", ""),
                "organism": organism_name,
                "geo_loc_name": geo,
                "collection_date": record.collection_date,
                "isolation_source": record.isolation_source,
                "platform": row.get("Platform", ""),
                "layout": row.get("LibraryLayout", ""),
                "bases": row.get("bases", ""),
                # Discovery never reviews a publication, so every row is tier B
                # until a curator supplies a real study identifier.
                "cohort_tier": "B",
                "phenotype_source": (
                    "biosample:" + ",".join(sorted(phenotype)) if phenotype else ""),
                "source_study": PENDING_REVIEW,
                "expected_outcome": expected_outcome,
            })
    return result
