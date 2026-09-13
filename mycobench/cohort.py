"""Cohort schema, evidence verification, and the lock that certifies a sheet.

A cohort sheet is a claim about public data. The lock is what makes the claim
checkable: it records the NCBI evidence retrieved at verification time plus a
checksum of the exact sheet that was verified, so a later edit cannot silently
inherit the earlier verification.

Tier grades curation confidence, not pass/fail:

``A`` every deposited-evidence check passed **and** ``source_study`` names a
    reviewed publication — release-ready.
``B`` every check passed but no publication has been reviewed — complete, not
    curator-approved. Everything discovery produces is tier B.
``C`` not verified online at all; online verification always resolves to A or B.

Origin is verified against the controlled ``geo_loc_name`` field, never against
a free-text search. That distinction is not pedantry: an SRA free-text search
for ``Nigeria`` returns runs deposited from the East China Sea.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from .ncbi import (
    BIOPROJECT_ACCESSION,
    BIOSAMPLE_ACCESSION,
    RUN_ACCESSION,
    Credentials,
    NCBIError,
    biosample_attributes,
    runinfo_for_accessions,
)

LOCK_SCHEMA_VERSION = "1.0"

REQUIRED_COLUMNS = ("sample_id", "sra_run", "biosample", "bioproject",
                    "organism", "platform", "layout", "cohort_tier")
OPTIONAL_COLUMNS = ("geo_loc_name", "collection_date", "isolation_source",
                    "bases", "phenotype_source", "source_study",
                    "expected_outcome", "lineage")
ALL_COLUMNS = REQUIRED_COLUMNS + OPTIONAL_COLUMNS

#: A source_study naming one of these is a placeholder, not a citation.
UNREVIEWED_STUDY = re.compile(r"pending|needs?_|unverified|unknown|tbd|candidate",
                              re.IGNORECASE)

#: Control-panel contract values. A row declaring one is scored as a control,
#: not as a benchmark subject.
EXPECTED_OUTCOMES = ("species_mismatch_refused",)

MTBC_ORGANISMS = ("mycobacterium tuberculosis", "mycobacterium canettii",
                  "mycobacterium africanum", "mycobacterium bovis",
                  "mycobacterium tuberculosis complex")


class CohortError(ValueError):
    """The sheet, or its lock, is not usable."""


def read_rows(path: str | Path) -> tuple[list[dict], list[str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        stream = (line for line in handle
                  if line.strip() and not line.lstrip().startswith("#"))
        reader = csv.DictReader(stream, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_rows(path: str | Path, rows: Iterable[dict],
               columns: Iterable[str], header_comment: str = "") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(columns)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        if header_comment:
            for line in header_comment.strip().splitlines():
                handle.write(f"# {line}\n" if line.strip() else "#\n")
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    return path


def is_mtbc(organism: str) -> bool:
    value = (organism or "").strip().lower()
    return any(value.startswith(name) for name in MTBC_ORGANISMS)


def derive_tier(evidence_errors: list[str],
                source_study: Optional[str]) -> Optional[str]:
    if evidence_errors:
        return None
    study = (source_study or "").strip()
    return "B" if not study or UNREVIEWED_STUDY.search(study) else "A"


def schema_errors(rows: list[dict], fields: list[str]) -> list[str]:
    missing = [c for c in REQUIRED_COLUMNS if c not in fields]
    if missing:
        return ["cohort sheet must include: " + ", ".join(missing)]

    errors: list[str] = []
    seen_ids: set[str] = set()
    seen_runs: set[str] = set()
    for number, row in enumerate(rows, 2):
        sample_id = (row.get("sample_id") or "").strip()
        if not sample_id or sample_id in seen_ids:
            errors.append(f"row {number}: missing or duplicate sample_id")
        seen_ids.add(sample_id)

        run = (row.get("sra_run") or "").strip()
        if not RUN_ACCESSION.match(run):
            errors.append(f"row {number}: invalid SRA run accession {run!r}")
        elif run in seen_runs:
            errors.append(f"row {number}: duplicate sra_run {run}")
        seen_runs.add(run)

        if not BIOSAMPLE_ACCESSION.match((row.get("biosample") or "").strip()):
            errors.append(f"row {number}: invalid BioSample accession")
        if not BIOPROJECT_ACCESSION.match((row.get("bioproject") or "").strip()):
            errors.append(f"row {number}: invalid BioProject accession")
        if (row.get("platform") or "").strip().upper() != "ILLUMINA":
            errors.append(f"row {number}: platform must be ILLUMINA")
        if (row.get("layout") or "").strip().upper() != "PAIRED":
            errors.append(f"row {number}: layout must be PAIRED")
        if (row.get("cohort_tier") or "").strip() not in ("A", "B", "C"):
            errors.append(f"row {number}: cohort_tier must be A, B or C")

        bases = (row.get("bases") or "").strip()
        if bases:
            try:
                if int(bases) <= 0:
                    errors.append(f"row {number}: bases must be positive")
            except ValueError:
                errors.append(f"row {number}: bases must be an integer")

        expected = (row.get("expected_outcome") or "").strip()
        if expected and expected not in EXPECTED_OUTCOMES:
            errors.append(
                f"row {number}: expected_outcome {expected!r} is not one of "
                + ", ".join(EXPECTED_OUTCOMES))

        # A control row must not be MTBC, and a benchmark row must be: the
        # whole point of the control panel is that its organism is wrong for
        # the profile, and mixing the two silently would invalidate both.
        organism = row.get("organism") or ""
        if expected == "species_mismatch_refused" and is_mtbc(organism):
            errors.append(
                f"row {number}: {organism} is MTBC but declares "
                f"expected_outcome=species_mismatch_refused")
        if not expected and not is_mtbc(organism):
            errors.append(
                f"row {number}: {organism} is not MTBC; a non-MTBC row needs "
                f"expected_outcome=species_mismatch_refused or it will be "
                f"scored as a benchmark subject")
    return errors


@dataclass
class RowEvidence:
    sample_id: str
    run: dict = field(default_factory=dict)
    biosample: dict = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    derived_tier: Optional[str] = None


def verify_rows(rows: list[dict], creds: Credentials,
                country: Optional[str] = None) -> list[RowEvidence]:
    """Verify every row against deposited NCBI metadata, in two batched walks."""
    runs = [(row.get("sra_run") or "").strip() for row in rows]
    run_info = runinfo_for_accessions(runs, creds)
    samples = biosample_attributes(
        [(row.get("biosample") or "").strip() for row in rows], creds)

    evidence: list[RowEvidence] = []
    for row in rows:
        record = RowEvidence(sample_id=(row.get("sample_id") or "").strip())
        run = (row.get("sra_run") or "").strip()
        info = run_info.get(run)
        if info is None:
            record.errors.append(f"SRA run {run} not found at NCBI")
            evidence.append(record)
            continue

        record.run = {
            "run": run, "biosample": info.get("BioSample", ""),
            "bioproject": info.get("BioProject", ""),
            "sample_name": info.get("SampleName", ""),
            "platform": info.get("Platform", ""),
            "layout": info.get("LibraryLayout", ""),
            "bases": info.get("bases", ""),
            "organism": info.get("ScientificName", ""),
        }
        if info.get("BioSample") != (row.get("biosample") or "").strip():
            record.errors.append("SRA BioSample does not match the cohort row")
        if info.get("BioProject") != (row.get("bioproject") or "").strip():
            record.errors.append("SRA BioProject does not match the cohort row")
        if (info.get("Platform") or "").upper() != "ILLUMINA":
            record.errors.append("SRA platform is not ILLUMINA")
        if (info.get("LibraryLayout") or "").upper() != "PAIRED":
            record.errors.append("SRA library is not paired-end")

        declared_bases = (row.get("bases") or "").strip()
        if declared_bases and info.get("bases") and declared_bases != info["bases"]:
            record.errors.append(
                f"declared bases {declared_bases} does not match NCBI "
                f"{info['bases']}")

        sample = samples.get((row.get("biosample") or "").strip())
        if sample is None:
            record.errors.append("BioSample record could not be retrieved")
        else:
            record.biosample = {
                "accession": sample.accession,
                "geo_loc_name": sample.geo_loc_name,
                "collection_date": sample.collection_date,
                "host": sample.host,
                "isolation_source": sample.isolation_source,
                "phenotype_attributes": sample.phenotype_attributes(),
            }
            declared_geo = (row.get("geo_loc_name") or "").strip()
            if declared_geo and declared_geo != sample.geo_loc_name:
                record.errors.append(
                    f"declared geo_loc_name {declared_geo!r} does not match "
                    f"deposited {sample.geo_loc_name!r}")
            if country and not sample.confirms_country(country):
                record.errors.append(
                    f"BioSample geo_loc_name {sample.geo_loc_name!r} does not "
                    f"confirm {country}")

            # A row claiming a phenotype source must have one; a row claiming
            # none must not silently have had one available.
            declared_source = (row.get("phenotype_source") or "").strip()
            available = sample.phenotype_attributes()
            if not declared_source and available:
                record.errors.append(
                    "BioSample carries phenotype attribute(s) "
                    f"({', '.join(sorted(available))}) but phenotype_source is "
                    f"empty; record the source or explain the omission")

        record.derived_tier = derive_tier(record.errors, row.get("source_study"))
        declared_tier = (row.get("cohort_tier") or "").strip()
        if record.derived_tier and declared_tier != record.derived_tier:
            record.errors.append(
                f"declared cohort_tier {declared_tier} does not match the "
                f"evidence-derived tier {record.derived_tier} "
                + ("(source_study names a reviewed publication)"
                   if record.derived_tier == "A"
                   else "(source_study absent or pending review)"))
        evidence.append(record)
        time.sleep(creds.interval)
    return evidence


def write_lock(path: str | Path, samples_path: str | Path,
               evidence: list[RowEvidence], country: Optional[str] = None) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    source = Path(samples_path).read_bytes()
    payload = {
        "schema_version": LOCK_SCHEMA_VERSION,
        "sample_sheet": Path(samples_path).name,
        "sample_sheet_sha256": hashlib.sha256(source).hexdigest(),
        "country_confirmed": country or "",
        "evidence": [
            {"sample_id": record.sample_id, "run": record.run,
             "biosample": record.biosample, "derived_tier": record.derived_tier,
             "errors": record.errors}
            for record in evidence
        ],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def verify_lock(lock_path: str | Path, samples_path: str | Path) -> int:
    """Confirm a sheet is exactly the one its lock certifies."""
    try:
        lock = json.loads(Path(lock_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CohortError(f"cannot read verification lock: {exc}") from exc

    if str(lock.get("schema_version", "")) != LOCK_SCHEMA_VERSION:
        raise CohortError(
            f"lock schema is {lock.get('schema_version') or 'absent'}, "
            f"expected {LOCK_SCHEMA_VERSION}; re-run --online --write-lock")

    expected = lock.get("sample_sheet_sha256")
    if not isinstance(expected, str) or not expected:
        raise CohortError("lock has no sample_sheet_sha256")
    observed = hashlib.sha256(Path(samples_path).read_bytes()).hexdigest()
    if observed != expected:
        raise CohortError(
            "sample-sheet checksum differs from its lock; the sheet changed "
            "after verification. Re-run --online --write-lock")

    evidence = lock.get("evidence")
    if not isinstance(evidence, list) or not evidence:
        raise CohortError("lock has no evidence records")

    failed = [record.get("sample_id", "?") for record in evidence
              if record.get("errors")]
    if failed:
        raise CohortError(
            "lock records evidence errors for: " + ", ".join(failed[:5])
            + (" ..." if len(failed) > 5 else ""))
    return len(evidence)


def lock_bases(samples_path: str | Path) -> dict[str, int]:
    """Per-sample base counts from a sibling lock, for download estimation."""
    sheet = Path(samples_path)
    for candidate in (sheet.with_suffix(".lock.json"),
                      sheet.parent / (sheet.stem + ".lock.json")):
        if not candidate.is_file():
            continue
        try:
            data = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        out = {}
        for record in data.get("evidence", []):
            try:
                out[record.get("sample_id", "")] = int(
                    (record.get("run") or {}).get("bases") or 0)
            except (TypeError, ValueError):
                continue
        if out:
            return out
    return {}


def summarise(rows: list[dict]) -> dict:
    """Counts a report and a reviewer both need up front."""
    projects: dict[str, int] = {}
    organisms: dict[str, int] = {}
    tiers: dict[str, int] = {}
    total_bases = 0
    controls = 0
    phenotyped = 0
    for row in rows:
        projects[row.get("bioproject", "")] = projects.get(row.get("bioproject", ""), 0) + 1
        organisms[row.get("organism", "")] = organisms.get(row.get("organism", ""), 0) + 1
        tiers[row.get("cohort_tier", "")] = tiers.get(row.get("cohort_tier", ""), 0) + 1
        try:
            total_bases += int(row.get("bases") or 0)
        except ValueError:
            pass
        if (row.get("expected_outcome") or "").strip():
            controls += 1
        if (row.get("phenotype_source") or "").strip():
            phenotyped += 1
    largest = max(projects.values()) if projects else 0
    return {
        "n_rows": len(rows), "n_bioprojects": len(projects),
        "bioprojects": projects, "organisms": organisms, "tiers": tiers,
        "total_bases": total_bases, "n_controls": controls,
        "n_phenotyped": phenotyped,
        "largest_bioproject_share": (largest / len(rows)) if rows else 0.0,
    }
