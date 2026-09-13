"""Parse TB-Profiler output into evidence Myconductor can reconcile.

Validation status, stated per parser
------------------------------------
``parse_collate`` is written against a **real TB-Profiler artifact**:
``tests/example_collate.txt`` in the TB-Profiler repository, whose header and
value format ("``rpoB p.Ser450Leu (1.00)``", "``-``" for no variant) are
reproduced in this module's tests. It is the parser to prefer.

``parse_results_json`` is written from the documented results schema and has
**not** been run against real output. It fails loudly on an unrecognised shape
rather than mis-parsing, and the schema-derived status is recorded on every
piece of evidence it emits.

The interpretation rule that matters
------------------------------------
TB-Profiler writes ``-`` for a drug with no resistance variant found. That is
**not** a susceptible call: it means no resistance mechanism was detected,
which is a different statement from "the loci were sequenced and are wild
type". Mapping ``-`` to susceptible is precisely the defect Myconductor was
rebuilt to remove, so this module maps it to ``NOT_ASSESSED`` and lets the
coverage gate decide. ``target_median_depth`` is carried through as the
coverage evidence that can license a susceptible call.
"""
from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

#: TB-Profiler drug columns, mapped to the names Myconductor's profile uses.
#: Drugs TB-Profiler reports that the mycobacterial profile does not cover are
#: mapped to None so they are recorded as uncompared rather than dropped.
COLLATE_DRUGS = {
    "rifampicin": "rifampicin",
    "rifapentine": None,
    "isoniazid": "isoniazid",
    "ethambutol": "ethambutol",
    "pyrazinamide": "pyrazinamide",
    "moxifloxacin": "moxifloxacin",
    "levofloxacin": None,
    "bedaquiline": "bedaquiline",
    "delamanid": None,
    "pretomanid": "pretomanid",
    "linezolid": "linezolid",
    "streptomycin": None,
    "amikacin": "amikacin",
    "kanamycin": None,
    "capreomycin": None,
    "clofazimine": "clofazimine",
    "ethionamide": None,
    "prothionamide": None,
}

REQUIRED_COLLATE_COLUMNS = ("sample", "main_lineage", "drtype")

#: "rpoB p.Ser450Leu (1.00)" — gene, change, allele fraction. Several variants
#: for one drug arrive comma-separated.
_VARIANT = re.compile(r"^\s*(?P<gene>\S+)\s+(?P<change>\S+)\s*"
                      r"(?:\((?P<freq>[0-9.]+)\))?\s*$")

NO_VARIANT = ("-", "", "NA", "None")


class EngineError(ValueError):
    """The engine output did not match the shape this parser understands."""


@dataclass
class EngineVariant:
    gene: str
    change: str
    frequency: Optional[float] = None

    def label(self) -> str:
        return f"{self.gene}_{self.change}"


@dataclass
class DrugFinding:
    drug: str
    variants: list[EngineVariant] = field(default_factory=list)
    supported: bool = True

    @property
    def resistance_detected(self) -> bool:
        return bool(self.variants)

    @property
    def minority(self) -> list[EngineVariant]:
        return [v for v in self.variants
                if v.frequency is not None and v.frequency < 0.90]


@dataclass
class EngineSample:
    """One isolate's TB-Profiler result, in this project's vocabulary."""

    sample: str
    main_lineage: str = ""
    sub_lineage: str = ""
    drtype: str = ""
    median_depth: Optional[float] = None
    pct_reads_mapped: Optional[float] = None
    num_dr_variants: Optional[int] = None
    num_other_variants: Optional[int] = None
    findings: dict[str, DrugFinding] = field(default_factory=dict)
    uncompared_drugs: dict[str, list[EngineVariant]] = field(default_factory=dict)
    parser: str = "collate"
    validated_parser: bool = True
    warnings: list[str] = field(default_factory=list)

    def resistance_calls(self) -> dict[str, str]:
        """``drug -> resistant | not_assessed``.

        Never ``susceptible``: TB-Profiler's absence of a resistance variant is
        not a coverage-backed susceptible call, and this module refuses to
        upgrade it into one.
        """
        out = {}
        for drug, finding in self.findings.items():
            out[drug] = "resistant" if finding.resistance_detected else "not_assessed"
        return out

    def all_variants(self) -> list[EngineVariant]:
        seen: dict[str, EngineVariant] = {}
        for finding in self.findings.values():
            for variant in finding.variants:
                seen[variant.label()] = variant
        for variants in self.uncompared_drugs.values():
            for variant in variants:
                seen[variant.label()] = variant
        return sorted(seen.values(), key=lambda v: v.label())

    def coverage_rows(self, loci: Iterable[str]) -> list[dict]:
        """Depth-table rows for Myconductor's callable mask.

        TB-Profiler collate reports one median depth across its whole target
        set, not per locus. That single figure is applied to every locus and
        the ``callable_fraction`` is left EMPTY on purpose: a genome-wide
        median says nothing about the breadth of any individual locus, and
        filling it in would manufacture the very evidence the mask exists to
        demand. Use per-gene coverage from the JSON output, or an independent
        depth table, to license susceptible calls.
        """
        if self.median_depth is None:
            return []
        return [{"locus": locus, "mean_depth": int(self.median_depth),
                 "callable_fraction": ""} for locus in sorted(set(loci))]


def _parse_variant_cell(cell: str) -> list[EngineVariant]:
    text = (cell or "").strip()
    if text in NO_VARIANT:
        return []
    variants = []
    for part in text.split(","):
        match = _VARIANT.match(part)
        if not match:
            continue
        frequency = match.group("freq")
        variants.append(EngineVariant(
            gene=match.group("gene"), change=match.group("change"),
            frequency=float(frequency) if frequency else None))
    return variants


def _optional_float(value: Optional[str]) -> Optional[float]:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _optional_int(value: Optional[str]) -> Optional[int]:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def parse_collate(path: str | Path) -> list[EngineSample]:
    """Parse ``tb-profiler collate`` TSV output. Validated against real output."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLLATE_COLUMNS if c not in fields]
        if missing:
            raise EngineError(
                f"{path}: collate output missing column(s) {missing}; found "
                f"{fields[:8]}. Expected `tb-profiler collate` TSV.")

        samples = []
        for row in reader:
            name = (row.get("sample") or "").strip()
            if not name:
                continue
            sample = EngineSample(
                sample=name,
                main_lineage=(row.get("main_lineage") or "").strip(),
                sub_lineage=(row.get("sub_lineage") or "").strip(),
                drtype=(row.get("drtype") or "").strip(),
                median_depth=_optional_float(row.get("target_median_depth")),
                pct_reads_mapped=_optional_float(row.get("pct_reads_mapped")),
                num_dr_variants=_optional_int(row.get("num_dr_variants")),
                num_other_variants=_optional_int(row.get("num_other_variants")),
                parser="collate", validated_parser=True)

            for column, mapped in COLLATE_DRUGS.items():
                if column not in fields:
                    continue
                variants = _parse_variant_cell(row.get(column, ""))
                if mapped is None:
                    if variants:
                        sample.uncompared_drugs[column] = variants
                    continue
                sample.findings[mapped] = DrugFinding(drug=mapped,
                                                      variants=variants)
            samples.append(sample)

    if not samples:
        raise EngineError(f"{path}: collate output contained no sample rows")
    return samples


def parse_results_json(path: str | Path) -> EngineSample:
    """Parse a single-sample TB-Profiler results JSON.

    Written from the documented schema and NOT validated against real output;
    ``validated_parser`` is False on what it returns, and the report says so.
    """
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise EngineError(f"{path}: not valid JSON") from exc
    if not isinstance(data, dict):
        raise EngineError(f"{path}: expected a JSON object at top level")

    variant_key = next((k for k in ("dr_variants", "drug_resistance_variants")
                        if isinstance(data.get(k), list)), None)
    if variant_key is None:
        raise EngineError(
            f"{path}: no 'dr_variants' or 'drug_resistance_variants' list. "
            f"Keys present: {sorted(data)[:10]}. This parser understands "
            f"TB-Profiler v4-v6 results JSON.")

    qc = data.get("qc") if isinstance(data.get("qc"), dict) else {}
    sample = EngineSample(
        sample=str(data.get("id") or data.get("sample_name") or Path(path).stem),
        main_lineage=str(data.get("main_lin") or data.get("main_lineage") or ""),
        sub_lineage=str(data.get("sublin") or data.get("sub_lineage") or ""),
        drtype=str(data.get("drtype") or ""),
        median_depth=_optional_float(qc.get("median_coverage")
                                     or qc.get("median_depth")),
        pct_reads_mapped=_optional_float(qc.get("percent_reads_mapped")
                                         or qc.get("pct_reads_mapped")),
        parser="results-json", validated_parser=False)
    sample.warnings.append(
        "parsed by a schema-derived JSON parser that has not been validated "
        "against real TB-Profiler output; prefer `tb-profiler collate` output")

    for raw in data[variant_key]:
        gene = raw.get("gene") or raw.get("gene_name")
        change = (raw.get("change") or raw.get("protein_change")
                  or raw.get("nucleotide_change"))
        if not (gene and change):
            sample.warnings.append(
                f"skipped a variant with no gene/change: {sorted(raw)[:6]}")
            continue
        frequency = raw.get("freq")
        variant = EngineVariant(
            gene=str(gene), change=str(change),
            frequency=float(frequency)
            if isinstance(frequency, (int, float)) else None)
        for annotation in (raw.get("drugs") or raw.get("annotation") or []):
            drug = (annotation.get("drug") or "").strip().lower()
            if not drug:
                continue
            mapped = COLLATE_DRUGS.get(drug, drug)
            if mapped is None:
                sample.uncompared_drugs.setdefault(drug, []).append(variant)
                continue
            finding = sample.findings.setdefault(
                mapped, DrugFinding(drug=mapped))
            if variant.label() not in {v.label() for v in finding.variants}:
                finding.variants.append(variant)

    # Drugs TB-Profiler assessed and found nothing for are absent from
    # dr_variants entirely; record them as not assessed rather than omitting
    # them, so a missing drug is never read as a clean result.
    for column, mapped in COLLATE_DRUGS.items():
        if mapped and mapped not in sample.findings:
            sample.findings[mapped] = DrugFinding(drug=mapped)
    return sample


def to_variant_tsv(sample: EngineSample) -> list[dict]:
    """Variants in the TSV shape Myconductor's input adapter reads."""
    rows = []
    for variant in sample.all_variants():
        rows.append({
            "gene": variant.gene,
            "change": variant.change,
            "vaf": f"{variant.frequency:.4f}" if variant.frequency is not None else "",
            "depth": int(sample.median_depth) if sample.median_depth else "",
            "region": "coding",
            "consequence": "",
            "platform": "illumina",
        })
    return rows


VARIANT_TSV_COLUMNS = ("gene", "change", "vaf", "depth", "region",
                       "consequence", "platform")
