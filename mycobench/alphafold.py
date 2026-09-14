"""Precomputed AlphaFold structures, imported — not folded here.

The entire *M. tuberculosis* proteome is already in the AlphaFold Protein
Structure Database. Running AlphaFold locally would be recomputing a lookup, at
the cost of GPU infrastructure, multi-terabyte sequence databases and a
restrictive model licence. So this fetches what exists, caches it, and reads two
things out of it.

This lives in ``mycobench``, not ``myconductor``, because the core is offline by
construction: every network call in this project sits in the harness, and the
core consumes files. The output here is a structural-annotation file that
``myconductor.modules.structural_annotation`` loads like any other attributed
import.

What an AlphaFold model can and cannot tell you
-----------------------------------------------
It gives **pLDDT**, a per-residue confidence. That is genuinely useful and it is
also a guard: a variant sitting in a pLDDT 40 region has no reliable structural
interpretation at all, and reporting one anyway would be inventing detail. Low
confidence is emitted as an explicit unavailability, not as a low score.

It does **not** give ligand contacts. AlphaFold DB entries are *apo* — predicted
from sequence, with no drug, cofactor, nucleic acid or ion present. There is no
rifampicin in the rpoB model. A "distance to the rifampicin binding pocket"
derived from these coordinates would be measuring a distance to something that
is not in the file, and this module refuses to emit one rather than letting a
plausible number reach a report. Ligand geometry needs an experimental holo
structure; that is a different input, and the refusal says so.

What it will report is distance to a *named reference residue* supplied by the
caller — for instance, to a residue that published work identifies as lining a
pocket. That keeps the provenance where it belongs: on the caller's citation,
not on an inference this module invented.

Residue numbering does not line up, and the mismatch is silent
--------------------------------------------------------------
The catalogue's codon numbers are not UniProt's. Measured against the live
database: the WHO catalogue's rpoB Ser450 — the most consequential rifampicin
determinant there is — sits at residue **456** of UniProt P9WGY9. Position 450
of that model is a threonine. The offset is a consistent +6 across Leu430,
Ser431, Asp435, His445 and Ser450.

Nothing about that failure announces itself. Every residue number in range
resolves to *some* residue with *some* pLDDT, so a naive mapping produces a
confident, well-formatted, entirely wrong annotation. So ``annotation_rows``
requires the caller to state which residue it expects to find and refuses the
variant when the model disagrees. Supplying the expectation is a small cost;
the alternative is structural annotations silently attached to the wrong amino
acid.
"""
from __future__ import annotations

import json
import math
import os
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from . import USER_AGENT

#: AlphaFold DB serves one predicted model per UniProt accession.
ALPHAFOLD_BASE = "https://alphafold.ebi.ac.uk/files"
#: Metadata endpoint, used to discover which model version currently exists.
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"

#: Deliberately not a default. AlphaFold DB re-releases the whole database
#: periodically and retires the old file paths: a version hardcoded here was
#: v4, and by the time it was first run the database served v6 and every fetch
#: 404'd. So the version is either supplied by the caller — which is what
#: pinning a run means — or resolved from the API and then *recorded* in the
#: output, so a result still says which model it was computed from.
RESOLVE_LATEST = None

#: The database's own confidence bands. Below 70 the backbone placement is
#: unreliable; below 50 the region is usually disordered and carries no
#: interpretable geometry at all.
PLDDT_CONFIDENT = 70.0
PLDDT_UNRELIABLE = 50.0

#: AlphaFold entries are predicted without ligands, so nothing here can measure
#: a drug contact. Stated once, used wherever a caller asks for one.
NO_LIGAND_REASON = (
    "AlphaFold models are predicted apo: no drug, cofactor or ion is present "
    "in the coordinates, so no ligand distance can be measured from them. Use "
    "an experimental holo structure for that.")


class AlphaFoldError(RuntimeError):
    """A structure could not be fetched or read."""


@dataclass(frozen=True)
class Residue:
    """One residue's alpha-carbon position and model confidence."""

    number: int
    name: str
    x: float
    y: float
    z: float
    plddt: float

    @property
    def confident(self) -> bool:
        return self.plddt >= PLDDT_CONFIDENT

    @property
    def unreliable(self) -> bool:
        return self.plddt < PLDDT_UNRELIABLE

    def distance_to(self, other: "Residue") -> float:
        return math.dist((self.x, self.y, self.z), (other.x, other.y, other.z))


@dataclass
class Structure:
    """A predicted model, keyed by the accession it was predicted for."""

    accession: str
    version: str
    residues: dict[int, Residue] = field(default_factory=dict)
    source_url: str = ""

    @property
    def mean_plddt(self) -> Optional[float]:
        if not self.residues:
            return None
        return sum(r.plddt for r in self.residues.values()) / len(self.residues)

    def residue(self, number: int) -> Optional[Residue]:
        return self.residues.get(number)

    def describe(self) -> str:
        mean = self.mean_plddt
        return (f"AlphaFold {self.accession} ({self.version}): "
                f"{len(self.residues)} residue(s)"
                + ("" if mean is None else f", mean pLDDT {mean:.1f}"))


def latest_version(accession: str, timeout: int = 60) -> str:
    """Ask the database which model version it currently serves.

    Returned as ``vN`` so it can be recorded verbatim as the source version of
    everything derived from that model.
    """
    accession = accession.strip().upper()
    if not accession:
        raise AlphaFoldError("a UniProt accession is required")
    url = f"{ALPHAFOLD_API}/{accession}"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read())
    except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
        raise AlphaFoldError(
            f"could not resolve a model version for {accession}: {exc}") from exc
    if not payload:
        raise AlphaFoldError(f"AlphaFold DB has no entry for {accession}")
    version = payload[0].get("latestVersion")
    if version is None:
        raise AlphaFoldError(f"no model version reported for {accession}")
    return f"v{version}"


def model_url(accession: str, version: str) -> str:
    accession = accession.strip().upper()
    if not accession:
        raise AlphaFoldError("a UniProt accession is required")
    return f"{ALPHAFOLD_BASE}/AF-{accession}-F1-model_{version}.pdb"


def parse_pdb(text: str, accession: str, version: str,
              source_url: str = "") -> Structure:
    """Read alpha-carbon positions and pLDDT from a predicted model.

    AlphaFold writes the per-residue pLDDT into the B-factor column, which is
    why a plain PDB read is enough and no structure library is needed.
    """
    structure = Structure(accession=accession.upper(), version=version,
                          source_url=source_url)
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            number = int(line[22:26])
            residue = Residue(
                number=number, name=line[17:20].strip(),
                x=float(line[30:38]), y=float(line[38:46]),
                z=float(line[46:54]), plddt=float(line[60:66]))
        except ValueError as exc:
            raise AlphaFoldError(
                f"malformed ATOM record in {accession}: {line[:60]!r}") from exc
        structure.residues[number] = residue
    if not structure.residues:
        raise AlphaFoldError(
            f"no alpha-carbon records found for {accession}; the file may not "
            f"be an AlphaFold model")
    return structure


def fetch_structure(accession: str, cache_dir: str | Path,
                    version: Optional[str] = RESOLVE_LATEST,
                    timeout: int = 60,
                    cached_only: bool = False) -> Structure:
    """Download one predicted model into a cache, or reuse what is there.

    ``version`` pins the model. Left unset it is resolved from the database and
    recorded on the returned structure, so the result still identifies which
    model produced it.
    """
    cache_dir = Path(cache_dir)
    accession = accession.strip().upper()
    if version is None:
        if cached_only:
            raise AlphaFoldError(
                f"no version given for {accession} and cached_only is set; "
                f"resolving the latest version requires the network")
        version = latest_version(accession, timeout=timeout)
    destination = (cache_dir / f"AF-{accession}-F1-model_{version}.pdb").resolve()
    if not destination.is_relative_to(cache_dir.resolve()):
        raise AlphaFoldError("structure path escapes cache directory")

    url = model_url(accession, version)
    if destination.is_file() and destination.stat().st_size > 0:
        return parse_pdb(destination.read_text(encoding="utf-8"),
                         accession, version, url)
    if cached_only:
        raise AlphaFoldError(f"not cached: {accession}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise AlphaFoldError(f"could not fetch {url}: {exc}") from exc
    if not payload:
        raise AlphaFoldError(f"{url} returned no bytes")

    # Same atomic write as the VCF cache: a process killed mid-write must not
    # leave a truncated model that every later run accepts.
    temporary = destination.with_name(f"{destination.name}.{os.getpid()}.part")
    try:
        temporary.write_bytes(payload)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return parse_pdb(payload.decode("utf-8", errors="replace"), accession,
                     version, url)


def ligand_distance(*_args, **_kwargs):
    """Always refuses. Kept as a named refusal rather than an absent function.

    A contributor looking for "how do I get the binding-pocket distance" finds
    this and the reason, instead of an empty space they might fill with a
    number measured against a ligand that is not in the file.
    """
    raise AlphaFoldError(NO_LIGAND_REASON)


#: Three-letter to one-letter, so a caller may state either spelling.
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def _matches(model_name: str, expected: str) -> bool:
    """Does the model's residue match what the caller expected?"""
    model_name = model_name.strip().upper()
    expected = str(expected).strip().upper()
    if not expected:
        return False
    if len(expected) == 1:
        return _THREE_TO_ONE.get(model_name) == expected
    return model_name == expected


def annotation_rows(structure: Structure, variants: Iterable[dict],
                    organism: str, reference_assembly: str, timestamp: str,
                    reference_residues: Optional[dict[int, str]] = None
                    ) -> list[dict]:
    """Structural annotation rows for variants on this protein.

    ``variants`` are dicts carrying ``variant_key``, ``residue`` and
    ``expected_residue``. The expectation is mandatory and checked: catalogue
    codon numbering and UniProt numbering differ — by +6 for rpoB, measured —
    and every number in range resolves to some residue, so an unchecked mapping
    annotates the wrong amino acid without any sign that it has.

    A variant in an unreliable region is emitted with ``predicted_effect``
    saying the structure cannot support an interpretation, rather than being
    dropped silently or given a confident-looking location.

    ``reference_residues`` optionally maps a residue number to a name the
    *caller* cites as functionally important; distances are reported to those
    and attributed to the caller, never to an inference made here.
    """
    reference_residues = reference_residues or {}
    rows = []
    for variant in variants:
        number = variant.get("residue")
        key = variant.get("variant_key")
        expected = variant.get("expected_residue")
        if not key or number is None:
            raise AlphaFoldError(
                "each variant needs a variant_key and a residue number")
        if not expected:
            raise AlphaFoldError(
                f"{key}: expected_residue is required. Catalogue codon "
                f"numbering and UniProt numbering differ (rpoB by +6), and an "
                f"unchecked mapping silently annotates the wrong residue")
        residue = structure.residue(int(number))
        if residue is not None and not _matches(residue.name, expected):
            rows.append(_row(key, "unknown", organism, reference_assembly,
                             structure, timestamp,
                             f"NUMBERING MISMATCH: residue {number} of "
                             f"{structure.accession} is {residue.name}, not "
                             f"the expected {str(expected).upper()}. The "
                             f"catalogue's codon numbering does not align with "
                             f"this model's; no annotation is made. Supply the "
                             f"UniProt residue number for this variant"))
            continue
        if residue is None:
            rows.append(_row(key, "unknown", organism, reference_assembly,
                             structure, timestamp,
                             f"residue {number} is absent from the predicted "
                             f"model; no structural interpretation"))
            continue
        if residue.unreliable:
            rows.append(_row(key, "unknown", organism, reference_assembly,
                             structure, timestamp,
                             f"residue {number} sits in a pLDDT "
                             f"{residue.plddt:.0f} region; the model does not "
                             f"support a structural interpretation here",
                             confidence=residue.plddt / 100.0))
            continue

        nearest = None
        for ref_number, ref_name in sorted(reference_residues.items()):
            reference = structure.residue(ref_number)
            if reference is None:
                continue
            distance = residue.distance_to(reference)
            if nearest is None or distance < nearest[0]:
                nearest = (distance, ref_name)

        effect = (f"residue {number} modelled at pLDDT {residue.plddt:.0f}")
        if nearest:
            effect += (f"; {nearest[0]:.1f} A from {nearest[1]} as cited by "
                       f"the caller")
        rows.append(_row(key, "unknown", organism, reference_assembly,
                         structure, timestamp, effect,
                         confidence=residue.plddt / 100.0))
    return rows


def _row(variant_key: str, location: str, organism: str,
         reference_assembly: str, structure: Structure, timestamp: str,
         effect: str, confidence: Optional[float] = None) -> dict:
    return {
        "variant_key": variant_key,
        # Never "binding_pocket": an apo model cannot place a drug.
        "location": location,
        "predicted_effect": effect,
        "source": f"alphafold-db:{structure.accession}",
        "organism": organism,
        "reference_assembly": reference_assembly,
        "source_version": structure.version,
        "timestamp": timestamp,
        "confidence": confidence,
        # Deliberately absent; see NO_LIGAND_REASON.
        "ligand_distance": None,
        "distance_unit": None,
        "ligand_reference": None,
        "validation_status": "unknown",
    }


def write_annotations(rows: list[dict], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return path
