"""CRyPTIC paired phenotypes — the only route here to an accuracy claim.

Public Nigerian BioSamples carry no DST or MIC metadata, so the Nigerian panels
can measure concordance and nothing more. The CRyPTIC compendium supplies what
is missing: 12,289 isolates with whole-genome sequencing and minimum inhibitory
concentrations for 13 antituberculars measured on one assay, released publicly
by the consortium.

Reuse terms
-----------
The consortium asks that work using this data cite the compendium paper. That
citation is recorded in the manifest this module writes, so it travels with the
derived cohort rather than living only in a README:

    The CRyPTIC Consortium, "A data compendium of Mycobacterium tuberculosis
    antibiotic resistance", https://doi.org/10.1101/2021.09.14.460274

Three details that decide whether a number means anything
---------------------------------------------------------
**Quality grade.** Every phenotype carries ``HIGH`` / ``MEDIUM`` / ``LOW``.
Measuring a predictor against a ``LOW``-quality phenotype measures the
phenotype, so only the grades in ``ACCEPTED_PHENOTYPE_QUALITY`` are admitted
and the rest are counted as not evaluable.

**Censored MICs.** Values arrive as ``<=0.25`` or ``>4.0`` — the true value is
outside the tested dilution range. Parsing those to ``0.25`` and ``4.0`` and
carrying on would silently convert a bound into a measurement, so the censoring
direction is retained alongside the number.

**Excluded samples.** The release ships its own withdrawal list. Anything on it
is dropped before the join, because a sample the consortium retracted is not
evidence.
"""
from __future__ import annotations

import csv
import hashlib
import io
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from . import USER_AGENT
from .thresholds import ACCEPTED_PHENOTYPE_QUALITY

FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reuse"
#: Pinned release. The directory holds several; naming one keeps an ingest
#: reproducible and makes an upgrade a deliberate act.
REUSE_TABLE = "CRyPTIC_reuse_table_20240917.csv"
EXCLUDED_TABLE = "CRyPTIC_excluded_samples_20220607.tsv"

CITATION = ('The CRyPTIC Consortium, "A data compendium of Mycobacterium '
            'tuberculosis antibiotic resistance", '
            'https://doi.org/10.1101/2021.09.14.460274')

#: CRyPTIC's three-letter codes mapped to the drug names Myconductor uses.
DRUG_CODES = {
    "AMI": "amikacin",
    "BDQ": "bedaquiline",
    "CFZ": "clofazimine",
    "DLM": "delamanid",
    "EMB": "ethambutol",
    "ETH": "ethionamide",
    "INH": "isoniazid",
    "KAN": "kanamycin",
    "LEV": "levofloxacin",
    "LZD": "linezolid",
    "MXF": "moxifloxacin",
    "RIF": "rifampicin",
    "RFB": "rifabutin",
}

REQUIRED_COLUMNS = ("ENA_RUN", "UNIQUEID", "ENA_SAMPLE")


class PhenotypeError(RuntimeError):
    """The phenotype release could not be fetched or did not parse."""


@dataclass(frozen=True)
class MIC:
    """One MIC, with its censoring preserved."""

    raw: str
    value: Optional[float]
    censoring: str = "none"   # none | left | right

    @property
    def is_censored(self) -> bool:
        return self.censoring != "none"

    def describe(self) -> str:
        return self.raw or "NA"


def parse_mic(raw: str) -> Optional[MIC]:
    """Parse an MIC, keeping ``<=`` / ``>`` as censoring rather than dropping it."""
    text = (raw or "").strip()
    if not text or text.upper() == "NA":
        return None
    censoring = "none"
    number = text
    for prefix, direction in (("<=", "left"), ("<", "left"),
                              (">=", "right"), (">", "right")):
        if text.startswith(prefix):
            censoring = direction
            number = text[len(prefix):]
            break
    try:
        return MIC(raw=text, value=float(number), censoring=censoring)
    except ValueError:
        return MIC(raw=text, value=None, censoring=censoring)


@dataclass
class DrugPhenotype:
    drug: str
    binary: Optional[str]          # R | S | None
    quality: Optional[str]         # HIGH | MEDIUM | LOW | None
    mic: Optional[MIC]

    @property
    def evaluable(self) -> bool:
        return (self.binary in ("R", "S")
                and (self.quality or "").upper() in ACCEPTED_PHENOTYPE_QUALITY)

    def reason_not_evaluable(self) -> Optional[str]:
        if self.binary not in ("R", "S"):
            return f"binary phenotype is {self.binary or 'absent'}"
        quality = (self.quality or "").upper()
        if quality not in ACCEPTED_PHENOTYPE_QUALITY:
            return (f"phenotype quality {quality or 'absent'} is below the "
                    f"accepted grade(s) "
                    f"{'/'.join(ACCEPTED_PHENOTYPE_QUALITY)}")
        return None


@dataclass
class IsolatePhenotypes:
    ena_run: str
    ena_sample: str
    unique_id: str
    drugs: dict[str, DrugPhenotype] = field(default_factory=dict)

    @property
    def evaluable_drugs(self) -> list[str]:
        return sorted(d for d, p in self.drugs.items() if p.evaluable)

    def binary_for(self, drug: str) -> Optional[str]:
        phenotype = self.drugs.get(drug)
        if phenotype is None or not phenotype.evaluable:
            return None
        return phenotype.binary


def _download(filename: str, dest_dir: str | Path,
              timeout: int = 600) -> tuple[Path, str, int]:
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    destination = dest_dir / filename
    request = urllib.request.Request(f"{FTP_BASE}/{filename}",
                                     headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except Exception as exc:  # noqa: BLE001 - reported with context
        raise PhenotypeError(
            f"could not fetch {filename} from the CRyPTIC release: {exc}") from exc
    destination.write_bytes(payload)
    return destination, hashlib.sha256(payload).hexdigest(), len(payload)


def fetch(dest_dir: str | Path) -> dict:
    """Download the pinned reuse table and the withdrawal list."""
    table, table_sha, table_size = _download(REUSE_TABLE, dest_dir)
    excluded, excluded_sha, excluded_size = _download(EXCLUDED_TABLE, dest_dir)
    return {
        "release": REUSE_TABLE,
        "citation": CITATION,
        "source": FTP_BASE,
        "files": {
            REUSE_TABLE: {"path": str(table), "sha256": table_sha,
                          "size": table_size},
            EXCLUDED_TABLE: {"path": str(excluded), "sha256": excluded_sha,
                             "size": excluded_size},
        },
    }


def read_excluded(path: str | Path) -> set[str]:
    """ENA sample accessions the consortium withdrew."""
    excluded: set[str] = set()
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("ENA_SAMPLE") or "").strip()
            if sample:
                excluded.add(sample)
    return excluded


def read_table(path: str | Path,
               excluded: Optional[set[str]] = None) -> list[IsolatePhenotypes]:
    """Parse the reuse table into per-isolate phenotype records."""
    excluded = excluded or set()
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLUMNS if c not in fields]
        if missing:
            raise PhenotypeError(
                f"{path}: reuse table missing column(s) {missing}; found "
                f"{fields[:6]}. This parser expects the "
                f"{REUSE_TABLE} layout.")

        isolates: list[IsolatePhenotypes] = []
        for row in reader:
            sample = (row.get("ENA_SAMPLE") or "").strip()
            if sample and sample in excluded:
                continue
            run = (row.get("ENA_RUN") or "").strip()
            if not run or run.upper() == "NONE":
                continue
            record = IsolatePhenotypes(
                ena_run=run, ena_sample=sample,
                unique_id=(row.get("UNIQUEID") or "").strip())
            for code, drug in DRUG_CODES.items():
                binary = (row.get(f"{code}_BINARY_PHENOTYPE") or "").strip().upper()
                quality = (row.get(f"{code}_PHENOTYPE_QUALITY") or "").strip().upper()
                record.drugs[drug] = DrugPhenotype(
                    drug=drug,
                    binary=binary if binary in ("R", "S") else None,
                    quality=quality or None,
                    mic=parse_mic(row.get(f"{code}_MIC") or ""))
            isolates.append(record)
    if not isolates:
        raise PhenotypeError(f"{path}: parsed no usable isolate rows")
    return isolates


def summarise(isolates: Iterable[IsolatePhenotypes]) -> dict:
    """Per-drug evaluable / resistant counts, before any prediction is made."""
    isolates = list(isolates)
    per_drug: dict[str, dict[str, int]] = {}
    for drug in DRUG_CODES.values():
        evaluable = resistant = susceptible = 0
        for isolate in isolates:
            phenotype = isolate.drugs.get(drug)
            if phenotype is None or not phenotype.evaluable:
                continue
            evaluable += 1
            if phenotype.binary == "R":
                resistant += 1
            else:
                susceptible += 1
        per_drug[drug] = {"evaluable": evaluable, "resistant": resistant,
                          "susceptible": susceptible}
    return {"n_isolates": len(isolates), "per_drug": per_drug}


COHORT_COLUMNS = ("sample_id", "sra_run", "biosample", "bioproject",
                  "organism", "geo_loc_name", "platform", "layout", "bases",
                  "cohort_tier", "phenotype_source", "source_study")


def build_cohort_rows(isolates: Iterable[IsolatePhenotypes],
                      run_info: dict[str, dict],
                      min_evaluable_drugs: int = 6) -> tuple[list[dict], list[str]]:
    """Turn phenotyped isolates into cohort rows, keeping only usable ones.

    ``run_info`` is NCBI run-info keyed by run accession, which is what supplies
    the BioSample, BioProject and base count. An isolate whose ENA run cannot be
    resolved at NCBI is rejected with that reason rather than carried with
    blank identifiers.
    """
    accepted: list[dict] = []
    rejected: list[str] = []
    for isolate in isolates:
        evaluable = isolate.evaluable_drugs
        if len(evaluable) < min_evaluable_drugs:
            rejected.append(
                f"{isolate.ena_run}: {len(evaluable)} drug(s) at accepted "
                f"phenotype quality, below the {min_evaluable_drugs} minimum")
            continue
        info = run_info.get(isolate.ena_run)
        if info is None:
            rejected.append(
                f"{isolate.ena_run}: run not resolvable at NCBI")
            continue
        if (info.get("Platform") or "").upper() != "ILLUMINA":
            rejected.append(f"{isolate.ena_run}: platform is not ILLUMINA")
            continue
        if (info.get("LibraryLayout") or "").upper() != "PAIRED":
            rejected.append(f"{isolate.ena_run}: library is not paired-end")
            continue
        accepted.append({
            "sample_id": f"cr_{isolate.ena_run}",
            "sra_run": isolate.ena_run,
            "biosample": info.get("BioSample", ""),
            "bioproject": info.get("BioProject", ""),
            "organism": info.get("ScientificName", "Mycobacterium tuberculosis"),
            "geo_loc_name": "",
            "platform": info.get("Platform", ""),
            "layout": info.get("LibraryLayout", ""),
            "bases": info.get("bases", ""),
            "cohort_tier": "B",
            "phenotype_source": f"cryptic:{REUSE_TABLE}",
            "source_study": "CRyPTIC_compendium_2022",
        })
    return accepted, rejected


def phenotype_table(isolates: Iterable[IsolatePhenotypes]) -> list[dict]:
    """Long-form phenotype rows for the validator and the report."""
    rows = []
    for isolate in isolates:
        for drug, phenotype in sorted(isolate.drugs.items()):
            rows.append({
                "sample_id": f"cr_{isolate.ena_run}",
                "sra_run": isolate.ena_run,
                "drug": drug,
                "phenotype": phenotype.binary or "",
                "quality": phenotype.quality or "",
                "mic": phenotype.mic.raw if phenotype.mic else "",
                "mic_value": (phenotype.mic.value
                              if phenotype.mic and phenotype.mic.value is not None
                              else ""),
                "mic_censoring": phenotype.mic.censoring if phenotype.mic else "",
                "evaluable": "yes" if phenotype.evaluable else "no",
                "not_evaluable_reason": phenotype.reason_not_evaluable() or "",
            })
    return rows


PHENOTYPE_COLUMNS = ("sample_id", "sra_run", "drug", "phenotype", "quality",
                     "mic", "mic_value", "mic_censoring", "evaluable",
                     "not_evaluable_reason")
