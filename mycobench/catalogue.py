"""Fetch and ingest the published WHO mutation catalogue.

The catalogue is not redistributed here. It is fetched from WHO's own
repository at a **pinned commit**, checksummed, and converted into the JSON
Myconductor's catalogue lane loads. Hand-transcribing a subset is how the
bundled illustrative file came to exist, and hand-transcription is exactly the
step that must never sit between a published catalogue and a clinical call.

Two files do the work
---------------------
``WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt``
    114-column TSV, one row per drug/variant pair. Carries the grading, the
    tier, the effect — and a ``genomic position`` column that is a *pointer*
    reading ``(see "Genomic_coordinates" sheet)`` rather than a coordinate.

``WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt``
    ``variant, chromosome, position, reference_nucleotide,
    alternative_nucleotide``. Joined on ``variant``, and one variant maps to
    **several** rows because the same change has both an MNV and an SNV
    representation. All of them are retained: a coordinate key that only
    matches one spelling defeats the purpose of having coordinates.

This is what closes Myconductor's null-coordinate gap. The bundled profile
declares no coordinates and refuses to invent any; an ingested catalogue
supplies real ones, so variants can be reconciled across engines on coordinate
identity instead of on a display label.

The tier column also yields the real per-drug gene sets, replacing the
illustrative ``drug_loci.json`` gene lists with WHO's own tier 1 / tier 2
assignment.
"""
from __future__ import annotations

import csv
import hashlib
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from . import USER_AGENT

#: WHO's own repository for the 2nd-edition catalogue, MIT licensed.
REPOSITORY = "GTB-tbsequencing/mutation-catalogue-2023"
#: Pinned so an ingest is reproducible. Bump deliberately, never silently.
PINNED_COMMIT = "0bb3914348c5a4c981859601447834c08f03ee3d"
RAW_BASE = f"https://raw.githubusercontent.com/{REPOSITORY}/{PINNED_COMMIT}"
DIRECTORY = "Final%20Result%20Files"

MASTER_FILE = "WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt"
COORDINATES_FILE = "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
LICENSE_FILE = "LICENSE.md"
WORKBOOK_FILE = "WHO-UCN-TB-2023.7-eng.xlsx"

CATALOGUE_VERSION = "WHO-UCN-TB-2023.7 (2nd edition)"
REFERENCE_ASSEMBLY = "NC_000962.3"

#: Master-file columns the ingester depends on, looked up by name rather than
#: index so a column insertion upstream cannot silently shift the meaning.
COL_DRUG = "drug"
COL_GENE = "gene"
COL_MUTATION = "mutation"
COL_VARIANT = "variant"
COL_TIER = "tier"
COL_EFFECT = "effect"
COL_FINAL_GRADING = "FINAL CONFIDENCE GRADING"
COL_INITIAL_GRADING = "INITIAL CONFIDENCE GRADING"
COL_SILENT = "Silent mutation"
COL_COMMENT = "Comment"

REQUIRED_MASTER_COLUMNS = (COL_DRUG, COL_GENE, COL_MUTATION, COL_VARIANT,
                           COL_TIER, COL_EFFECT, COL_FINAL_GRADING)
REQUIRED_COORDINATE_COLUMNS = ("variant", "chromosome", "position",
                               "reference_nucleotide", "alternative_nucleotide")

#: WHO's five grading groups, mapped onto Myconductor's call vocabulary.
#: Groups 1-2 associate with resistance; 3 is explicitly uncertain; 4-5 state
#: the variant is NOT associated, which is a susceptible-leaning statement
#: about the variant and still does not by itself make a drug usable.
GRADING_CALLS = {
    "1) assoc w r": "resistant",
    "2) assoc w r - interim": "resistant",
    "3) uncertain significance": "indeterminate",
    "4) not assoc w r - interim": "not_associated",
    "5) not assoc w r": "not_associated",
}

#: Ordinal weights for reconciliation only. The grading is not a probability
#: and this value is labelled as such wherever it surfaces.
GRADING_WEIGHT = {"resistant": 0.9, "indeterminate": 0.5,
                  "not_associated": 0.9}


class CatalogueError(RuntimeError):
    """The catalogue could not be fetched, or did not match the expected shape."""


@dataclass
class FetchedFile:
    name: str
    path: Path
    sha256: str
    size: int

    def describe(self) -> str:
        return f"{self.name} ({self.size:,} bytes, sha256:{self.sha256[:12]})"


def _url(filename: str) -> str:
    return f"{RAW_BASE}/{DIRECTORY}/{urllib.parse.quote(filename)}"


def fetch_file(filename: str, dest_dir: str | Path,
               timeout: int = 600) -> FetchedFile:
    """Download one catalogue file from the pinned commit and checksum it."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    destination = dest_dir / filename
    request = urllib.request.Request(_url(filename),
                                     headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except Exception as exc:  # noqa: BLE001 - reported with context below
        raise CatalogueError(
            f"could not fetch {filename} from {REPOSITORY}@"
            f"{PINNED_COMMIT[:8]}: {exc}") from exc
    destination.write_bytes(payload)
    return FetchedFile(filename, destination,
                       hashlib.sha256(payload).hexdigest(), len(payload))


def fetch_all(dest_dir: str | Path,
              include_workbook: bool = False) -> dict[str, FetchedFile]:
    """Fetch the files an ingest needs, plus the licence that governs reuse."""
    names = [MASTER_FILE, COORDINATES_FILE, LICENSE_FILE]
    if include_workbook:
        names.append(WORKBOOK_FILE)
    return {name: fetch_file(name, dest_dir) for name in names}


def read_coordinates(path: str | Path) -> dict[str, list[dict]]:
    """``variant -> [coordinate forms]``, all of them.

    One variant legitimately has several rows (an MNV spelling and an SNV
    spelling of the same change). Keeping only the first would make coordinate
    matching depend on which spelling an annotator emitted.

    Streamed rather than read whole: the master file this sits beside is 37 MB
    across 114 columns, and materialising every row as a dict costs far more
    memory than the file does on disk.
    """
    out: dict[str, list[dict]] = defaultdict(list)
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in REQUIRED_COORDINATE_COLUMNS
                   if c not in (reader.fieldnames or [])]
        if missing:
            raise CatalogueError(
                f"{path}: coordinates file missing column(s) {missing}; found "
                f"{(reader.fieldnames or [])[:8]}")
        for row in reader:
            variant = (row.get("variant") or "").strip()
            position = (row.get("position") or "").strip()
            if not variant or not position.isdigit():
                continue
            out[variant].append({
                "assembly": (row.get("chromosome")
                             or REFERENCE_ASSEMBLY).strip(),
                "position": int(position),
                "ref": (row.get("reference_nucleotide") or "").strip(),
                "alt": (row.get("alternative_nucleotide") or "").strip(),
            })
    return dict(out)


def coordinate_keys(forms: Iterable[dict]) -> list[str]:
    """Myconductor coordinate identity strings for every retained form."""
    keys = []
    for form in forms:
        if not (form.get("ref") and form.get("alt")):
            continue
        keys.append(f"{form['assembly']}:Chromosome:{form['position']}:"
                    f"{form['ref']}>{form['alt']}")
    return keys


@dataclass
class IngestResult:
    catalogue_version: str = CATALOGUE_VERSION
    entries: list[dict] = field(default_factory=list)
    drug_tiers: dict[str, dict[str, list[str]]] = field(default_factory=dict)
    skipped: list[str] = field(default_factory=list)
    drugs: set[str] = field(default_factory=set)
    provenance: dict = field(default_factory=dict)
    n_with_coordinates: int = 0

    @property
    def n_entries(self) -> int:
        return len(self.entries)

    def to_catalogue_json(self) -> dict:
        return {
            "catalogue_version": self.catalogue_version,
            "illustrative": False,
            "note": (f"Ingested from {REPOSITORY}@{PINNED_COMMIT[:8]}: "
                     f"{self.n_entries} variant entries across "
                     f"{len(self.drugs)} drug(s); "
                     f"{self.n_with_coordinates} carry genomic coordinates."),
            "reference_assembly": REFERENCE_ASSEMBLY,
            "provenance": self.provenance,
            "efflux_genes": ["Rv0678", "mmpR5", "pepQ", "Rv1979c"],
            "variants": self.entries,
        }

    def to_drug_loci_json(self, template: Optional[dict] = None) -> dict:
        """Per-drug tier 1 / tier 2 gene sets, taken from WHO's own tiering."""
        loci = sorted({gene for tiers in self.drug_tiers.values()
                       for genes in tiers.values() for gene in genes})
        payload = {
            "profile": "mtbc",
            "profile_version": f"who-2023.7+{PINNED_COMMIT[:8]}",
            "reference_assembly": REFERENCE_ASSEMBLY,
            "note": (f"Gene tiers derived from the WHO catalogue master file "
                     f"({REPOSITORY}@{PINNED_COMMIT[:8]}). Locus lengths are "
                     f"still null: the catalogue supplies variant coordinates, "
                     f"not gene spans, so a callable fraction must come from "
                     f"the input."),
            "loci": {gene: {"length_bp": None, "region": "coding"}
                     for gene in loci},
            "drugs": {drug: {"tier1": tiers.get("1", []),
                             "tier2": tiers.get("2", [])}
                      for drug, tiers in sorted(self.drug_tiers.items())},
            "efflux_regulators": (template or {}).get("efflux_regulators", {}),
            "promoter_loci": (template or {}).get("promoter_loci", {}),
        }
        return payload


def ingest(master_path: str | Path, coordinates_path: str | Path,
           provenance: Optional[dict] = None) -> IngestResult:
    """Convert the WHO master file into Myconductor's catalogue JSON."""
    coordinates = read_coordinates(coordinates_path)
    result = IngestResult(provenance=provenance or {})
    merged: dict[tuple[str, str], dict] = {}
    tiers: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))

    with open(master_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [c for c in REQUIRED_MASTER_COLUMNS if c not in fields]
        if missing:
            raise CatalogueError(
                f"{master_path}: master file missing column(s) {missing}. "
                f"Found {len(fields)} columns starting {fields[:6]}. The "
                f"ingester depends on these by name, so a renamed column is "
                f"reported rather than guessed at.")

        for row in reader:
            drug = (row.get(COL_DRUG) or "").strip().lower()
            gene = (row.get(COL_GENE) or "").strip()
            mutation = (row.get(COL_MUTATION) or "").strip()
            variant = (row.get(COL_VARIANT) or "").strip()
            grading = (row.get(COL_FINAL_GRADING) or "").strip()
            tier = (row.get(COL_TIER) or "").strip()

            if not (drug and gene and mutation):
                result.skipped.append(
                    f"incomplete row: {variant or '?'} / {drug or '?'}")
                continue
            if tier in ("1", "2"):
                tiers[drug][tier].add(gene)

            call = GRADING_CALLS.get(grading.lower())
            if call is None:
                result.skipped.append(
                    f"{variant or gene + '_' + mutation}/{drug}: grading "
                    f"{grading!r} is not one of WHO's five groups")
                continue

            result.drugs.add(drug)
            key = (gene, mutation)
            entry = merged.get(key)
            if entry is None:
                keys = coordinate_keys(coordinates.get(variant, []))
                entry = {
                    "gene": gene,
                    "change": mutation,
                    "variant": variant,
                    "drugs": [],
                    "call": call,
                    "who_grade": grading,
                    "confidence": GRADING_WEIGHT[call],
                    "confidence_note": ("ordinal weight derived from the WHO "
                                        "grading group; not a probability"),
                    "tier": tier,
                    "effect": (row.get(COL_EFFECT) or "").strip(),
                    "silent": (row.get(COL_SILENT) or "").strip().lower() in
                              ("1", "true", "yes"),
                }
                if keys:
                    entry["coordinate_key"] = keys[0]
                    entry["coordinate_keys"] = keys
                    result.n_with_coordinates += 1
                comment = (row.get(COL_COMMENT) or "").strip()
                if comment:
                    entry["comment"] = comment
                merged[key] = entry

            if drug not in entry["drugs"]:
                entry["drugs"].append(drug)
            # A variant graded resistance-associated for ANY drug keeps that
            # call: the entry is shared across its drugs and the safer call
            # must win.
            if call == "resistant" and entry["call"] != "resistant":
                entry["call"] = "resistant"
                entry["who_grade"] = grading
                entry["confidence"] = GRADING_WEIGHT["resistant"]

    if not merged:
        raise CatalogueError(
            f"{master_path}: no usable entries parsed "
            f"({len(result.skipped)} row(s) skipped; first: "
            f"{result.skipped[:2]})")

    result.entries = list(merged.values())
    result.drug_tiers = {drug: {tier: sorted(genes) for tier, genes in by_tier.items()}
                         for drug, by_tier in tiers.items()}
    return result


def ingest_from_dir(directory: str | Path,
                    fetched: Optional[dict[str, FetchedFile]] = None
                    ) -> IngestResult:
    directory = Path(directory)
    provenance = {
        "repository": REPOSITORY,
        "commit": PINNED_COMMIT,
        "catalogue_version": CATALOGUE_VERSION,
        "files": ({name: {"sha256": f.sha256, "size": f.size}
                   for name, f in fetched.items()} if fetched else {}),
        "licence": "MIT (see LICENSE.md in the WHO repository)",
    }
    return ingest(directory / MASTER_FILE, directory / COORDINATES_FILE,
                  provenance)


def write_outputs(result: IngestResult, catalogue_out: str | Path,
                  drug_loci_out: Optional[str | Path] = None,
                  drug_loci_template: Optional[str | Path] = None) -> list[Path]:
    written = []
    catalogue_out = Path(catalogue_out)
    catalogue_out.parent.mkdir(parents=True, exist_ok=True)
    catalogue_out.write_text(
        json.dumps(result.to_catalogue_json(), indent=2) + "\n",
        encoding="utf-8")
    written.append(catalogue_out)

    if drug_loci_out:
        template = None
        if drug_loci_template and Path(drug_loci_template).is_file():
            template = json.loads(
                Path(drug_loci_template).read_text(encoding="utf-8"))
        drug_loci_out = Path(drug_loci_out)
        drug_loci_out.write_text(
            json.dumps(result.to_drug_loci_json(template), indent=2) + "\n",
            encoding="utf-8")
        written.append(drug_loci_out)
    return written
