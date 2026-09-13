"""Input adapter: normalise heterogeneous inputs to canonical variants.

What this layer is responsible for
----------------------------------
Variant identity. Two engines can only be reconciled if the same biological
change resolves to the same key, which means the adapter has to do real
normalisation rather than pass labels through:

* **multiallelic records** are split into one variant per alternate allele
  (a single VCF line with ``ALT=A,G`` is two variants, not one);
* **indels are left-aligned and trimmed** so ``CTT -> C`` and ``TT -> .``
  representations of the same deletion agree;
* **FILTER and QUAL** are carried and honoured — a record failing the caller's
  own filters is not evidence;
* **consequence** is normalised to an enum, and synonymy is derived rather
  than trusted from a free-text flag;
* **multi-sample VCFs** require the sample to be named, instead of silently
  using column 10.

What it refuses to do
---------------------
Guess. Structural variants, symbolic alternate alleles and copy-number records
cannot be represented by the current model, so they are **rejected with a
recorded reason** rather than coerced into a SNV-shaped object that would then
be interpreted as one. Likewise the adapter never infers coverage from the
variant list; susceptibility requires a callable mask (``io.callable_mask``).

Not implemented, and reported as such: reference-genome liftover between
assemblies, and HGVS re-normalisation across annotator dialects.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..core.models import (
    Consequence,
    MTB_ASSEMBLY,
    QCFinding,
    Region,
    Variant,
    VariantIdentity,
)

_REGION_ALIASES = {
    "coding": Region.CODING, "cds": Region.CODING, "exon": Region.CODING,
    "promoter": Region.PROMOTER, "upstream": Region.PROMOTER,
    "intergenic": Region.INTERGENIC, "noncoding": Region.INTERGENIC,
    "rrna": Region.RRNA, "rna": Region.RRNA,
}

_CONSEQUENCE_ALIASES = {
    "missense": Consequence.MISSENSE, "missense_variant": Consequence.MISSENSE,
    "nonsynonymous": Consequence.MISSENSE, "non_synonymous": Consequence.MISSENSE,
    "synonymous": Consequence.SYNONYMOUS, "silent": Consequence.SYNONYMOUS,
    "synonymous_variant": Consequence.SYNONYMOUS,
    "nonsense": Consequence.NONSENSE, "stop_gained": Consequence.NONSENSE,
    "stop": Consequence.NONSENSE,
    "frameshift": Consequence.FRAMESHIFT,
    "frameshift_variant": Consequence.FRAMESHIFT,
    "inframe_deletion": Consequence.INFRAME_INDEL,
    "inframe_insertion": Consequence.INFRAME_INDEL,
    "inframe": Consequence.INFRAME_INDEL,
    "upstream": Consequence.UPSTREAM, "upstream_gene_variant": Consequence.UPSTREAM,
    "rrna": Consequence.RRNA, "rrna_variant": Consequence.RRNA,
    "feature_ablation": Consequence.DELETION, "gene_deletion": Consequence.DELETION,
}

#: Symbolic and structural ALT forms the current variant model cannot express.
_SYMBOLIC_ALT_PREFIXES = ("<", "[", "]", ".")


class AdapterError(ValueError):
    """The input could not be parsed at all."""


@dataclass
class RejectedRecord:
    raw: str
    reason: str


@dataclass
class AdapterResult:
    sample_id: str
    variants: list[Variant] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    qc: list[QCFinding] = field(default_factory=list)
    rejected: list[RejectedRecord] = field(default_factory=list)
    platform: Optional[str] = None
    assembly: str = MTB_ASSEMBLY
    unresolved_coordinates: int = 0


# -- normalisation --------------------------------------------------------
def normalise_alleles(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Trim and left-align a simple indel or substitution.

    Removes the shared suffix, then the shared prefix (advancing ``pos``),
    always leaving at least one base on each side. This is what makes
    ``(100, "CTT", "C")`` and ``(101, "TT", "")``-style representations of the
    same event converge on one key.
    """
    ref, alt = ref.upper(), alt.upper()
    # Shared suffix.
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # Shared prefix.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def _consequence_from(effect: Optional[str], change: str,
                      region: Region) -> Consequence:
    """Normalise a consequence, deriving it where the input omits one."""
    if effect:
        mapped = _CONSEQUENCE_ALIASES.get(effect.strip().lower())
        if mapped:
            return mapped

    c = (change or "").strip()
    lowered = c.lower()
    if "fs" in lowered or "frameshift" in lowered:
        return Consequence.FRAMESHIFT
    if "stop" in lowered or c.endswith("*"):
        return Consequence.NONSENSE
    if "del" in lowered or "ins" in lowered or "dup" in lowered:
        return Consequence.INFRAME_INDEL
    if region is Region.RRNA:
        return Consequence.RRNA
    if region in (Region.PROMOTER, Region.INTERGENIC) or c.startswith(("c.-", "-")):
        return Consequence.UPSTREAM
    # Derive synonymy: "A542A" has the same reference and alternate residue.
    if _is_synonymous_protein_change(c):
        return Consequence.SYNONYMOUS
    if _is_missense_protein_change(c):
        return Consequence.MISSENSE
    return Consequence.UNKNOWN


def _protein_parts(change: str) -> Optional[tuple[str, str, str]]:
    c = change.strip()
    if c.startswith("p."):
        c = c[2:]
    if len(c) < 3 or not c[0].isalpha():
        return None
    digits = "".join(ch for ch in c if ch.isdigit())
    if not digits:
        return None
    head = c[: c.index(digits[0])] if digits[0] in c else ""
    tail = c[c.index(digits) + len(digits):] if digits in c else ""
    if not head:
        return None
    return head, digits, tail


def _is_synonymous_protein_change(change: str) -> bool:
    parts = _protein_parts(change)
    if not parts:
        return False
    head, _, tail = parts
    return bool(tail) and head.upper() == tail.upper()


def _is_missense_protein_change(change: str) -> bool:
    parts = _protein_parts(change)
    if not parts:
        return False
    head, _, tail = parts
    return bool(tail) and tail.isalpha() and head.upper() != tail.upper()


def _region_of(value: Optional[str], change: str) -> Region:
    if value:
        mapped = _REGION_ALIASES.get(value.strip().lower())
        if mapped:
            return mapped
    if (change or "").strip().startswith(("c.-", "-")):
        return Region.PROMOTER
    return Region.CODING


def _opt_float(v: Optional[str]) -> Optional[float]:
    if v in (None, "", "."):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _opt_int(v: Optional[str]) -> Optional[int]:
    if v in (None, "", "."):
        return None
    try:
        return int(float(v))
    except ValueError:
        return None


def _info_dict(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            out[k.strip().upper()] = v.strip()
        elif token.strip():
            out[token.strip().upper()] = "true"
    return out


# -- VCF ------------------------------------------------------------------
def _parse_vcf(text: str, sample: Optional[str] = None) -> AdapterResult:
    result = AdapterResult(sample_id="sample")
    header: list[str] = []
    sample_cols: list[str] = []

    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line.strip():
            continue

        if line.startswith("##"):
            if line.startswith("##sampleID="):
                result.sample_id = line.split("=", 1)[1].strip()
            elif line.lower().startswith("##platform="):
                result.platform = line.split("=", 1)[1].strip()
            elif line.lower().startswith("##reference="):
                result.assembly = line.split("=", 1)[1].strip()
            continue

        if line.startswith("#CHROM"):
            header = line.lstrip("#").split("\t")
            sample_cols = header[9:] if len(header) > 9 else []
            continue

        f = line.split("\t")
        if len(f) < 8:
            result.rejected.append(RejectedRecord(
                line, f"fewer than 8 VCF columns ({len(f)})"))
            continue

        chrom, pos_s, _id, ref, alt_field, qual_s, filt, info = f[:8]
        meta = _info_dict(info)

        # Structural records. Previously rejected wholesale, which meant the
        # pipeline could not see insertional inactivation of mmpR5 — a leading
        # route to bedaquiline resistance. A gene-disrupting SV in a named
        # gene is now represented as a truncating change, because that is what
        # it is; only records with no interpretable gene are still refused.
        svtype = meta.get("SVTYPE", "").upper()
        if svtype or ("END" in meta and not alt_field.strip()):
            gene_name = meta.get("GENE", "")
            disrupting = svtype in ("DEL", "INS", "DUP", "INV", "BND", "CNV")
            if not gene_name or not disrupting:
                result.rejected.append(RejectedRecord(
                    line,
                    f"structural record (SVTYPE={svtype or 'absent'}) with no "
                    f"interpretable gene; not representable"))
                continue
            consequence = (Consequence.INSERTION if svtype in ("INS", "DUP")
                           else Consequence.DELETION)
            label = meta.get("AACHANGE") or f"{svtype.lower()}_disruption"
            identity = VariantIdentity(
                gene=gene_name, assembly=result.assembly, chrom=chrom,
                pos=_opt_int(pos_s), hgvs_p=label, consequence=consequence)
            result.variants.append(Variant(
                identity=identity,
                depth=_opt_int(meta.get("DP")),
                region=_region_of(meta.get("REGION"), label),
                filters=tuple(f for f in filt.split(";") if f) if filt else (),
                qual=_opt_float(qual_s),
                platform=result.platform))
            result.warnings.append(
                f"{gene_name}: structural record SVTYPE={svtype} represented as "
                f"a gene-disrupting {consequence.value}. Breakpoint-level "
                f"evidence is NOT assessed, so the disruption is asserted from "
                f"the caller's own claim.")
            continue

        # Sample column selection for multi-sample files.
        fmt_map: dict[str, str] = {}
        if len(f) >= 10 and sample_cols:
            if len(sample_cols) > 1 and sample is None:
                raise AdapterError(
                    f"multi-sample VCF with {len(sample_cols)} samples "
                    f"({', '.join(sample_cols)}); pass sample=<name> to choose "
                    f"one. Refusing to guess."
                )
            target = sample or sample_cols[0]
            if target not in sample_cols:
                raise AdapterError(
                    f"sample {target!r} not in VCF; available: "
                    f"{', '.join(sample_cols)}"
                )
            col = 9 + sample_cols.index(target)
            if col < len(f):
                fmt_map = dict(zip(f[8].split(":"), f[col].split(":")))
            if sample:
                result.sample_id = target
        elif result.sample_id == "sample" and sample_cols:
            result.sample_id = sample_cols[0]

        filters = tuple(x for x in filt.split(";") if x) if filt else ()
        qual = _opt_float(qual_s)
        gene = meta.get("GENE", chrom)
        region = _region_of(meta.get("REGION"), meta.get("AACHANGE", ""))
        pos = _opt_int(pos_s)

        # Multiallelic split: one variant per alternate allele.
        alts = [a.strip() for a in alt_field.split(",") if a.strip()]
        if not alts:
            result.rejected.append(RejectedRecord(line, "no alternate allele"))
            continue

        af_values = _split_numeric(fmt_map.get("AF") or meta.get("AF"))
        ad_values = _split_numeric(fmt_map.get("AD"))
        depth = _opt_int(fmt_map.get("DP") or meta.get("DP"))
        gq = _opt_int(fmt_map.get("GQ"))

        for i, alt in enumerate(alts):
            if alt.startswith(_SYMBOLIC_ALT_PREFIXES):
                result.rejected.append(RejectedRecord(
                    line, f"symbolic alternate allele {alt!r}; not "
                          f"representable by the current variant model"))
                continue

            n_pos, n_ref, n_alt = (pos, ref, alt)
            if pos is not None:
                n_pos, n_ref, n_alt = normalise_alleles(pos, ref, alt)

            change = meta.get("AACHANGE") or meta.get("HGVS_P") or ""
            hgvs_c = meta.get("HGVS_C")
            if not change and not hgvs_c and n_pos is not None:
                change = f"{n_ref}{n_pos}{n_alt}"

            consequence = _consequence_from(meta.get("EFFECT"), change, region)

            vaf = af_values[i] if i < len(af_values) else (
                af_values[0] if len(af_values) == 1 else None)
            alt_depth = None
            if len(ad_values) >= len(alts) + 1:
                alt_depth = int(ad_values[i + 1])
                if depth is None:
                    depth = int(sum(ad_values))
            if vaf is None and alt_depth is not None and depth:
                vaf = alt_depth / depth

            identity = VariantIdentity(
                gene=gene,
                assembly=result.assembly,
                chrom=chrom,
                pos=n_pos,
                ref=n_ref if n_pos is not None else None,
                alt=n_alt if n_pos is not None else None,
                hgvs_c=hgvs_c or (change if change.startswith(("c.", "-")) else None),
                hgvs_p=(change if change and not change.startswith(("c.", "-")) else None),
                consequence=consequence,
            )
            if not identity.is_coordinate_resolved:
                result.unresolved_coordinates += 1

            result.variants.append(Variant(
                identity=identity,
                vaf=vaf,
                depth=depth,
                alt_depth=alt_depth,
                region=region,
                filters=filters,
                qual=qual,
                genotype_quality=gq,
                platform=result.platform,
            ))

    if not header:
        raise AdapterError("no #CHROM header line found; is this a VCF?")
    return result


def _split_numeric(value: Optional[str]) -> list[float]:
    if not value:
        return []
    out: list[float] = []
    for part in value.split(","):
        f = _opt_float(part)
        if f is not None:
            out.append(f)
    return out


# -- TSV / JSON -----------------------------------------------------------
def _parse_tsv(text: str) -> AdapterResult:
    rows = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not rows:
        raise AdapterError("empty TSV")
    header = [h.strip().lower() for h in rows[0].split("\t")]
    for required in ("gene", "change"):
        if required not in header:
            raise AdapterError(f"TSV missing required column {required!r}")

    result = AdapterResult(sample_id="sample")
    for ln in rows[1:]:
        cols = dict(zip(header, [c.strip() for c in ln.split("\t")]))
        result.platform = result.platform or (cols.get("platform") or None)
        change = cols["change"]
        region = _region_of(cols.get("region"), change)
        variant = Variant.of(
            cols["gene"], change,
            consequence=_consequence_from(cols.get("consequence") or cols.get("effect"),
                                          change, region),
            region=region,
            vaf=_opt_float(cols.get("vaf")),
            depth=_opt_int(cols.get("depth")),
            alt_depth=_opt_int(cols.get("alt_depth")),
            platform=cols.get("platform") or None,
        )
        if not variant.identity.is_coordinate_resolved:
            result.unresolved_coordinates += 1
        result.variants.append(variant)
    return result


def _parse_json(text: str) -> AdapterResult:
    data = json.loads(text)
    if "variants" not in data:
        raise AdapterError("JSON input has no 'variants' key")
    result = AdapterResult(
        sample_id=data.get("sample_id", "sample"),
        platform=data.get("platform"),
        assembly=data.get("assembly", MTB_ASSEMBLY),
    )
    for v in data["variants"]:
        change = v.get("change", "")
        region = _region_of(v.get("region"), change)
        variant = Variant.of(
            v["gene"], change,
            assembly=result.assembly,
            chrom=v.get("chrom"),
            pos=v.get("pos"),
            ref=v.get("ref"),
            alt=v.get("alt"),
            consequence=_consequence_from(v.get("consequence") or v.get("effect"),
                                          change, region),
            region=region,
            vaf=v.get("vaf"),
            depth=v.get("depth"),
            alt_depth=v.get("alt_depth"),
            platform=v.get("platform") or result.platform,
        )
        if not variant.identity.is_coordinate_resolved:
            result.unresolved_coordinates += 1
        result.variants.append(variant)
    return result


# -- entry point ----------------------------------------------------------
def load(path: str | Path, depth_floor: int = 10,
         sample: Optional[str] = None,
         platform: Optional[str] = None) -> AdapterResult:
    """Parse an input file and apply record-level QC gating.

    Low-depth and filter-failing records are removed from the evidence path and
    recorded. Removing them does not make a drug susceptible: with the variant
    gone, the drug falls through to the callable-locus check, which reports
    NOT_ASSESSED unless a mask says otherwise.
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    suffix = path.suffix.lower()

    if suffix == ".vcf":
        result = _parse_vcf(text, sample=sample)
    elif suffix in (".tsv", ".txt"):
        result = _parse_tsv(text)
    elif suffix == ".json":
        result = _parse_json(text)
    else:
        raise AdapterError(
            f"unsupported input format {suffix!r}; expected .vcf, .tsv or .json"
        )

    if platform:
        result.platform = platform
        for v in result.variants:
            v.platform = platform

    kept: list[Variant] = []
    dropped_depth, dropped_filter = [], []
    for v in result.variants:
        if not v.passed_filters:
            dropped_filter.append(v)
            result.warnings.append(
                f"{v.label()}: FILTER={','.join(v.filters)}; record excluded "
                f"from evidence (the caller did not trust it)."
            )
            continue
        if v.depth is not None and v.depth < depth_floor:
            dropped_depth.append(v)
            result.warnings.append(
                f"{v.label()}: depth {v.depth} below floor {depth_floor}x; "
                f"record excluded from evidence. This does NOT make the "
                f"corresponding drug susceptible — it becomes NOT_ASSESSED "
                f"unless a callable mask shows the locus was covered."
            )
            continue
        kept.append(v)
    result.variants = kept

    result.qc.extend(_qc_findings(result, depth_floor, dropped_depth,
                                  dropped_filter))
    return result


def _qc_findings(result: AdapterResult, depth_floor: int,
                 dropped_depth: list[Variant],
                 dropped_filter: list[Variant]) -> list[QCFinding]:
    qc = [
        QCFinding("input_parsed", "pass",
                  f"{len(result.variants)} variant record(s) accepted"),
    ]
    if dropped_depth:
        qc.append(QCFinding(
            "record_depth", "warn",
            f"{len(dropped_depth)} record(s) below {depth_floor}x excluded: "
            + ", ".join(v.label() for v in dropped_depth)))
    if dropped_filter:
        qc.append(QCFinding(
            "record_filters", "warn",
            f"{len(dropped_filter)} record(s) failed caller FILTER: "
            + ", ".join(v.label() for v in dropped_filter)))
    if result.rejected:
        qc.append(QCFinding(
            "unrepresentable_records", "warn",
            f"{len(result.rejected)} record(s) rejected: "
            + "; ".join(sorted({r.reason for r in result.rejected}))))
    if result.unresolved_coordinates:
        qc.append(QCFinding(
            "variant_coordinates", "warn",
            f"{result.unresolved_coordinates} variant(s) carry no genomic "
            f"coordinates, so they can only be matched by gene/label. "
            f"Cross-engine reconciliation on these is unreliable."))
    qc.append(QCFinding(
        "sequencing_platform",
        "pass" if result.platform else "warn",
        f"platform: {result.platform}" if result.platform else
        "no platform declared (##platform= in VCF, or platform= in TSV/JSON); "
        "minority-allele assessment cannot apply a limit of detection"))
    return qc
