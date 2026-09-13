"""TB-Profiler adapter.

Written from TB-Profiler's documented JSON results schema (v4-v6 key names,
which differ between versions). **Not validated against real output.**

Two deliberate choices:

*Resistance and uncertainty are imported; susceptibility is not inferred.*
TB-Profiler reports a per-drug verdict, but a "susceptible" verdict there has
the same coverage caveat it has everywhere. This adapter imports resistance and
uncertain calls as evidence and lets Myconductor's own callable-locus gate
decide susceptibility.

*Its coverage output becomes our mask.* Where TB-Profiler's QC block exposes
per-gene depth and callable fraction, the adapter converts it into a
``CallableMask``. That is orchestration doing real work: a validated engine's
coverage assessment is what licenses Myconductor's susceptible calls.
"""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Any, Optional

from ..core.models import (
    Call,
    Consequence,
    DrugEvidence,
    EngineRef,
    Lane,
    LocusCoverage,
    QCFinding,
    Tier,
    Variant,
    VariantIdentity,
)
from ..io.callable_mask import CallableMask
from .base import AdapterSchemaError, EngineAdapter, EngineReport

#: WHO-style confidence gradings, mapped to calls. Anything not listed here is
#: treated as uncertain rather than guessed at.
_GRADE_CALLS = {
    "assoc w r": Call.RESISTANT,
    "assoc w r - interim": Call.RESISTANT,
    "1) assoc w r": Call.RESISTANT,
    "2) assoc w r - interim": Call.RESISTANT,
    "uncertain significance": Call.INDETERMINATE,
    "3) uncertain significance": Call.INDETERMINATE,
    "not assoc w r": Call.SUSCEPTIBLE,
    "not assoc w r - interim": Call.SUSCEPTIBLE,
    "4) not assoc w r - interim": Call.SUSCEPTIBLE,
    "5) not assoc w r": Call.SUSCEPTIBLE,
}

_CONSEQUENCE = {
    "missense_variant": Consequence.MISSENSE,
    "synonymous_variant": Consequence.SYNONYMOUS,
    "stop_gained": Consequence.NONSENSE,
    "frameshift_variant": Consequence.FRAMESHIFT,
    "inframe_deletion": Consequence.INFRAME_INDEL,
    "inframe_insertion": Consequence.INFRAME_INDEL,
    "upstream_gene_variant": Consequence.UPSTREAM,
    "non_coding_transcript_exon_variant": Consequence.RRNA,
    "feature_ablation": Consequence.DELETION,
}


class TBProfilerAdapter(EngineAdapter):
    def __init__(self, version: str = "unknown",
                 catalogue_version: str = "unknown"):
        self.engine = EngineRef(
            name="tb-profiler", version=version,
            database="tbdb", database_version=catalogue_version,
        )

    def parse(self, path: str | Path) -> EngineReport:
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise AdapterSchemaError(
                f"{path}: not valid JSON; TB-Profiler results are expected as "
                f"the .results.json file"
            ) from exc

        if not isinstance(data, dict):
            raise AdapterSchemaError(f"{path}: expected a JSON object at top level")

        variant_key = _first_present(data, ("dr_variants", "drug_resistance_variants"))
        if variant_key is None:
            raise AdapterSchemaError(
                f"{path}: no 'dr_variants' or 'drug_resistance_variants' key. "
                f"Keys present: {sorted(data)[:12]}. This adapter understands "
                f"TB-Profiler v4-v6 results JSON."
            )

        engine = EngineRef(
            name=self.engine.name,
            version=str(data.get("tbprofiler_version",
                                 data.get("version", self.engine.version))),
            database=self.engine.database,
            database_version=str(data.get("db_version",
                                          {}).get("name", self.engine.database_version))
            if isinstance(data.get("db_version"), dict)
            else str(data.get("db_version", self.engine.database_version)),
        )

        report = EngineReport(
            engine=engine,
            source_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            sample_id=data.get("id") or data.get("sample_name"),
            lineage=_lineage_of(data),
        )

        for raw in data[variant_key]:
            report.evidence.extend(self._evidence_for(raw, engine, report))

        # Other (non-resistance) variants become candidates for the VUS
        # workbench rather than evidence.
        for raw in data.get("other_variants", []) or []:
            variant = _variant_of(raw)
            if variant is not None:
                report.variants.append(variant)

        report.mask = _mask_from_qc(data)
        report.qc.extend(_qc_of(data, report.mask))
        return report

    def _evidence_for(self, raw: dict, engine: EngineRef,
                      report: EngineReport) -> list[DrugEvidence]:
        variant = _variant_of(raw)
        if variant is None:
            report.warnings.append(
                f"skipped a dr_variant with no gene/change: {sorted(raw)[:8]}")
            return []
        report.variants.append(variant)
        if str(raw.get("filter", "pass")).lower() not in ("pass", "."):
            report.qc.append(QCFinding("engine_variant_filter", "fail", "TB-Profiler variant failed FILTER"))

        annotations = raw.get("drugs") or raw.get("annotation") or []
        if not annotations:
            report.warnings.append(
                f"{variant.label()}: no drug annotation on this variant")
            return []

        out: list[DrugEvidence] = []
        for ann in annotations:
            drug = (ann.get("drug") or "").strip().lower()
            if not drug:
                continue
            grade = (ann.get("confidence") or ann.get("who_confidence")
                     or ann.get("confers") or "")
            call = _GRADE_CALLS.get(str(grade).strip().lower())
            if call is None:
                call = Call.INDETERMINATE
                grade_note = (f"unrecognised confidence grading {grade!r}; "
                              f"treated as uncertain rather than guessed")
            else:
                grade_note = ""

            # A rule of this codebase: only graded evidence establishes
            # resistance. TB-Profiler's WHO gradings qualify; anything
            # ungraded drops to INDETERMINATE.
            tier = Tier.CATALOGUED if grade else Tier.NONE
            if call is Call.RESISTANT and not tier.may_establish_resistance:
                call = Call.INDETERMINATE

            limitations = [f"imported from {engine}"]
            if grade_note:
                limitations.append(grade_note)
            if variant.vaf is not None and variant.vaf < 0.9:
                limitations.append(
                    f"called at {variant.vaf:.0%} allele fraction; see the "
                    f"minority-allele assessment")

            out.append(DrugEvidence(
                drug=drug, call=call, tier=tier, lane=Lane.ENGINE,
                confidence=None, variant=variant.identity,
                who_grade=str(grade) or None, engine=engine,
                limitations=tuple(limitations), scope="variant",
                sample_id=report.sample_id, catalogue_version=engine.database_version,
                rationale=(f"TB-Profiler reports {variant.label()} for {drug} "
                           f"(grading: {grade or 'none'})."),
            ))
        return out


# -- helpers --------------------------------------------------------------
def _first_present(data: dict, keys: tuple[str, ...]) -> Optional[str]:
    for k in keys:
        if k in data and isinstance(data[k], list):
            return k
    return None


def _variant_of(raw: dict) -> Optional[Variant]:
    gene = raw.get("gene") or raw.get("gene_name") or raw.get("locus_tag")
    change = (raw.get("change") or raw.get("protein_change")
              or raw.get("nucleotide_change") or raw.get("hgvs_p"))
    if not gene or not change:
        return None

    pos = raw.get("genome_pos") or raw.get("pos")
    freq = raw.get("freq")
    depth = raw.get("depth") or raw.get("total_depth")
    alt_depth = raw.get("alt_depth") or raw.get("ad")
    if isinstance(alt_depth, list):
        alt_depth = alt_depth[-1] if alt_depth else None

    identity = VariantIdentity(
        gene=str(gene),
        chrom=raw.get("chrom") or "Chromosome",
        pos=int(pos) if isinstance(pos, (int, float)) else None,
        ref=raw.get("ref"),
        alt=raw.get("alt"),
        hgvs_p=str(change) if str(change).startswith(("p.", "P.")) else None,
        hgvs_c=str(change) if str(change).startswith(("c.", "n.", "r.")) else None,
        consequence=_CONSEQUENCE.get(str(raw.get("type", "")).lower(),
                                     Consequence.UNKNOWN),
    )
    if identity.hgvs_p is None and identity.hgvs_c is None:
        identity = VariantIdentity(
            gene=identity.gene, chrom=identity.chrom, pos=identity.pos,
            ref=identity.ref, alt=identity.alt, hgvs_p=str(change),
            consequence=identity.consequence,
        )

    return Variant(
        identity=identity,
        vaf=float(freq) if isinstance(freq, (int, float)) else None,
        depth=int(depth) if isinstance(depth, (int, float)) else None,
        alt_depth=int(alt_depth) if isinstance(alt_depth, (int, float)) else None,
    )


def _lineage_of(data: dict) -> Optional[str]:
    for key in ("sublin", "sub_lineage", "main_lin", "main_lineage"):
        value = data.get(key)
        if isinstance(value, str) and value:
            return value
    lineage = data.get("lineage")
    if isinstance(lineage, list) and lineage:
        last = lineage[-1]
        if isinstance(last, dict):
            return last.get("lineage") or last.get("lin")
    return None


def _mask_from_qc(data: dict) -> Optional[CallableMask]:
    """Convert TB-Profiler's per-gene coverage into a callable mask."""
    qc = data.get("qc")
    if not isinstance(qc, dict):
        return None
    rows = None
    for key in ("gene_coverage", "target_coverage", "gene_qc", "target_qc"):
        value = qc.get(key)
        if isinstance(value, list) and value:
            rows = value
            break
    if rows is None:
        return None

    loci: dict[str, LocusCoverage] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        locus = row.get("gene") or row.get("locus_tag") or row.get("target")
        if not locus:
            continue
        fraction = _as_fraction(
            next((row[k] for k in ("fraction", "percent_depth_pass", "fraction_covered")
                  if row.get(k) is not None), None)
        )
        depth = row.get("median_depth") or row.get("depth")
        loci[str(locus)] = LocusCoverage(
            locus=str(locus),
            mean_depth=int(depth) if isinstance(depth, (int, float)) else None,
            callable_fraction=fraction,
            source="tb-profiler qc",
        )
    return CallableMask(loci, source="tb-profiler qc") if loci else None


def _as_fraction(value: Any) -> Optional[float]:
    if not isinstance(value, (int, float)):
        return None
    v = float(value)
    return v / 100.0 if v > 1.0 else v


def _qc_of(data: dict, mask: Optional[CallableMask]) -> list[QCFinding]:
    findings: list[QCFinding] = []
    qc = data.get("qc") if isinstance(data.get("qc"), dict) else {}
    depth = qc.get("median_coverage") or qc.get("median_depth")
    if isinstance(depth, (int, float)):
        findings.append(QCFinding(
            "engine_median_depth", "pass" if depth >= 20 else "warn",
            f"TB-Profiler median depth {depth}x"))
    mapped = qc.get("percent_reads_mapped") or qc.get("pct_reads_mapped")
    if isinstance(mapped, (int, float)):
        findings.append(QCFinding(
            "engine_reads_mapped", "pass" if mapped >= 90 else "warn",
            f"TB-Profiler mapped {mapped}% of reads"))
    findings.append(QCFinding(
        "engine_gene_coverage",
        "pass" if mask else "warn",
        "per-gene coverage imported as a callable mask" if mask else
        "TB-Profiler output carried no per-gene coverage; susceptibility will "
        "need a mask from elsewhere"))
    return findings
