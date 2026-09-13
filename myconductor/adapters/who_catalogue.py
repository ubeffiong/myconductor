"""WHO mutation catalogue ingester.

The bundled ``mtb_amr_catalogue.json`` is a 14-entry illustrative subset. The
real catalogue (2nd edition) grades roughly 11,390 variants for 13 medicines
from more than 52,000 isolates. This module ingests it.

Why the data is not in this repository
--------------------------------------
The catalogue is published by WHO as a spreadsheet under its own terms. It is
not redistributed here, so it must be obtained from WHO and pointed at with
``--catalogue``. That is the honest arrangement: hand-transcribing a subset is
how the current file came to exist, and hand-transcription is precisely the
step that should never sit between a published catalogue and a clinical call.

``.xlsx`` requires the optional ``openpyxl`` dependency (``pip install
myconductor[catalogue]``). A CSV or TSV export needs no dependency at all.

Column names differ between the catalogue's sheets and editions, so the mapping
is explicit and overridable, and a missing required column is an error listing
what was found rather than a silent empty catalogue.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from ..core.models import Call
from .base import AdapterSchemaError

#: Accepted spellings for each field we need, lowercased.
DEFAULT_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "variant": ("variant", "mutation", "variant_name", "final_variant_name"),
    "gene": ("gene", "gene_name", "locus", "tier1_gene"),
    "drug": ("drug", "antibiotic", "medicine"),
    "grade": ("final confidence grading", "final_confidence_grading",
              "confidence", "grading", "who_confidence"),
    "position": ("genome position", "genome_position", "position", "pos"),
    "effect": ("effect", "consequence", "variant_type"),
    "tier": ("tier", "gene_tier"),
}

#: The catalogue's gradings, mapped onto our call vocabulary. Groups 1 and 2
#: are resistance-associated; 3 is uncertain; 4 and 5 are not associated.
GRADE_CALLS: dict[str, Call] = {
    "1) assoc w r": Call.RESISTANT,
    "2) assoc w r - interim": Call.RESISTANT,
    "assoc w r": Call.RESISTANT,
    "assoc w r - interim": Call.RESISTANT,
    "3) uncertain significance": Call.INDETERMINATE,
    "uncertain significance": Call.INDETERMINATE,
    "4) not assoc w r - interim": Call.SUSCEPTIBLE,
    "5) not assoc w r": Call.SUSCEPTIBLE,
    "not assoc w r": Call.SUSCEPTIBLE,
    "not assoc w r - interim": Call.SUSCEPTIBLE,
}

#: Confidence is not a probability in the catalogue; these are ordinal weights
#: for reconciliation only, and are labelled as such wherever they surface.
_GRADE_WEIGHT = {
    Call.RESISTANT: 0.9,
    Call.INDETERMINATE: 0.5,
    Call.SUSCEPTIBLE: 0.9,
}


@dataclass
class IngestResult:
    catalogue_version: str
    entries: list[dict] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    drugs: set[str] = field(default_factory=set)

    @property
    def n_entries(self) -> int:
        return len(self.entries)

    def to_bundled_json(self) -> dict:
        """Emit the shape ``modules.catalogue.CatalogueModule`` loads."""
        return {
            "catalogue_version": self.catalogue_version,
            "illustrative": False,
            "note": (f"Ingested from a published WHO catalogue export: "
                     f"{self.n_entries} entries across "
                     f"{len(self.drugs)} drug(s)."),
            "efflux_genes": ["Rv0678", "mmpR5", "pepQ", "Rv1979c"],
            "variants": self.entries,
        }

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_bundled_json(), indent=2))
        return path


def _resolve_columns(header: Iterable[str],
                     aliases: Optional[dict[str, tuple[str, ...]]] = None
                     ) -> dict[str, str]:
    aliases = aliases or DEFAULT_COLUMN_ALIASES
    lowered = {h.strip().lower(): h for h in header}
    resolved: dict[str, str] = {}
    for field_name, candidates in aliases.items():
        for candidate in candidates:
            if candidate in lowered:
                resolved[field_name] = lowered[candidate]
                break
    for required in ("variant", "drug", "grade"):
        if required not in resolved:
            raise AdapterSchemaError(
                f"catalogue export is missing a {required!r} column. "
                f"Tried: {aliases[required]}. "
                f"Columns found: {sorted(lowered)[:20]}. "
                f"Pass column_aliases= to override."
            )
    return resolved


def _split_variant(variant: str) -> tuple[Optional[str], Optional[str]]:
    """Split ``katG_p.Ser315Thr`` or ``rpoB_S450L`` into gene and change."""
    v = variant.strip()
    if "_" not in v:
        return None, None
    gene, change = v.split("_", 1)
    return (gene or None), (change or None)


def ingest_csv(path: str | Path, catalogue_version: Optional[str] = None,
               column_aliases: Optional[dict[str, tuple[str, ...]]] = None,
               delimiter: Optional[str] = None) -> IngestResult:
    """Ingest a CSV/TSV export of the WHO catalogue."""
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig")
    if delimiter is None:
        delimiter = "\t" if path.suffix.lower() in (".tsv", ".tab") else ","

    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    if not reader.fieldnames:
        raise AdapterSchemaError(f"{path}: no header row")
    cols = _resolve_columns(reader.fieldnames, column_aliases)

    result = IngestResult(
        catalogue_version=catalogue_version or f"ingested:{path.name}")
    merged: dict[tuple[str, str], dict] = {}

    for row in reader:
        variant = (row.get(cols["variant"]) or "").strip()
        drug = (row.get(cols["drug"]) or "").strip().lower()
        grade = (row.get(cols["grade"]) or "").strip()
        if not variant or not drug or not grade:
            result.skipped.append(f"incomplete row: {variant or '?'} / {drug or '?'}")
            continue

        call = GRADE_CALLS.get(grade.lower())
        if call is None:
            result.skipped.append(f"{variant}/{drug}: unrecognised grading {grade!r}")
            continue

        gene = (row.get(cols["gene"]) if "gene" in cols else None) or None
        gene_part, change = _split_variant(variant)
        gene = (gene or gene_part or "").strip()
        change = (change or variant).strip()
        if not gene:
            result.skipped.append(f"{variant}: cannot determine gene")
            continue

        key = (gene, change)
        entry = merged.get(key)
        if entry is None:
            entry = {
                "gene": gene,
                "change": change,
                "drugs": [],
                "call": call.value if call is not Call.SUSCEPTIBLE else "not_associated",
                "who_grade": grade,
                "confidence": _GRADE_WEIGHT.get(call, 0.5),
                "confidence_note": ("ordinal weight derived from the WHO "
                                    "grading; not a probability"),
            }
            position = row.get(cols["position"]) if "position" in cols else None
            if position and str(position).strip().isdigit():
                entry["genome_position"] = int(str(position).strip())
            if "tier" in cols and row.get(cols["tier"]):
                entry["tier"] = str(row[cols["tier"]]).strip()
            merged[key] = entry

        if drug not in entry["drugs"]:
            entry["drugs"].append(drug)
        result.drugs.add(drug)

        # A variant graded resistance-associated for any drug keeps that call.
        if call is Call.RESISTANT:
            entry["call"] = "resistant"
            entry["who_grade"] = grade
            entry["confidence"] = _GRADE_WEIGHT[Call.RESISTANT]

    if not merged:
        raise AdapterSchemaError(
            f"{path}: parsed no usable catalogue entries. "
            f"Skipped {len(result.skipped)} row(s); first few: "
            f"{result.skipped[:3]}"
        )

    result.entries = list(merged.values())
    return result


def ingest_xlsx(path: str | Path, sheet: Optional[str] = None,
                catalogue_version: Optional[str] = None,
                column_aliases: Optional[dict[str, tuple[str, ...]]] = None
                ) -> IngestResult:
    """Ingest the catalogue's published spreadsheet.

    Requires the optional ``openpyxl`` dependency; the error says so rather
    than failing obscurely.
    """
    try:
        from openpyxl import load_workbook  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise AdapterSchemaError(
            "reading .xlsx requires openpyxl: pip install 'myconductor[catalogue]'. "
            "Alternatively export the sheet to CSV/TSV and use ingest_csv, "
            "which needs no dependencies."
        ) from exc

    wb = load_workbook(Path(path), read_only=True, data_only=True)
    ws = wb[sheet] if sheet else wb[wb.sheetnames[0]]
    rows = ws.iter_rows(values_only=True)
    try:
        header = [str(c) if c is not None else "" for c in next(rows)]
    except StopIteration as exc:
        raise AdapterSchemaError(f"{path}: sheet is empty") from exc

    lines = [",".join(_csv_escape(h) for h in header)]
    for row in rows:
        lines.append(",".join(_csv_escape("" if c is None else str(c)) for c in row))

    import io
    reader = csv.DictReader(io.StringIO("\n".join(lines)))
    cols = _resolve_columns(reader.fieldnames or [], column_aliases)
    del cols  # validated above; ingest_csv re-resolves on the same header

    tmp = "\n".join(lines)
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile("w", suffix=".csv", delete=False,
                            encoding="utf-8") as fh:
        fh.write(tmp)
        tmp_path = fh.name
    try:
        return ingest_csv(
            tmp_path,
            catalogue_version=catalogue_version or f"ingested:{Path(path).name}",
            column_aliases=column_aliases,
            delimiter=",",
        )
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _csv_escape(value: str) -> str:
    if any(ch in value for ch in (',', '"', "\n")):
        return '"' + value.replace('"', '""') + '"'
    return value
