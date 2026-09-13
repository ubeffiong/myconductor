"""Adapter layer: normalise heterogeneous inputs to canonical Variants.

Accepts a minimal VCF, a TSV, or JSON so the demo doesn't depend on a real
variant caller. Applies coverage-gated QC: a locus sequenced below the depth
floor is dropped and recorded as a warning, mirroring the "only report a drug
susceptible when its loci were actually sequenced" rule that good pipelines use.

Supported inputs
----------------
VCF   : uses INFO fields GENE, AACHANGE, REGION plus FORMAT AF/DP when present.
TSV   : columns gene, change, [vaf], [depth], [region], [silent]
JSON  : {"sample_id": "...", "variants": [ {gene, change, ...}, ... ]}
"""
from __future__ import annotations

import json
from pathlib import Path

from ..core.models import Region, Variant

_REGION_ALIASES = {
    "coding": Region.CODING,
    "cds": Region.CODING,
    "promoter": Region.PROMOTER,
    "intergenic": Region.INTERGENIC,
    "rrna": Region.RRNA,
}


def _region(value: str | None) -> Region:
    if not value:
        return Region.CODING
    return _REGION_ALIASES.get(value.strip().lower(), Region.CODING)


class AdapterResult:
    def __init__(self, sample_id: str, variants: list[Variant], warnings: list[str]):
        self.sample_id = sample_id
        self.variants = variants
        self.warnings = warnings


def _gate(variants: list[Variant], depth_floor: int) -> tuple[list[Variant], list[str]]:
    kept, warnings = [], []
    for v in variants:
        if v.depth and v.depth < depth_floor:
            warnings.append(
                f"{v.key()}: depth {v.depth} < floor {depth_floor}; locus dropped "
                f"(call would be indeterminate, not susceptible)."
            )
            continue
        kept.append(v)
    return kept, warnings


def _parse_tsv(text: str) -> tuple[str, list[Variant]]:
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    header = [h.strip().lower() for h in lines[0].split("\t")]
    variants = []
    for line in lines[1:]:
        cols = dict(zip(header, [c.strip() for c in line.split("\t")]))
        variants.append(
            Variant(
                gene=cols["gene"],
                change=cols["change"],
                vaf=float(cols.get("vaf", "1.0")),
                depth=int(cols.get("depth", "0")),
                region=_region(cols.get("region")),
                silent=cols.get("silent", "").lower() in ("1", "true", "yes"),
            )
        )
    return "sample", variants


def _parse_json(text: str) -> tuple[str, list[Variant]]:
    data = json.loads(text)
    sid = data.get("sample_id", "sample")
    variants = [
        Variant(
            gene=v["gene"],
            change=v["change"],
            vaf=float(v.get("vaf", 1.0)),
            depth=int(v.get("depth", 0)),
            region=_region(v.get("region")),
            silent=bool(v.get("silent", False)),
        )
        for v in data["variants"]
    ]
    return sid, variants


def _info(info: str) -> dict[str, str]:
    out = {}
    for token in info.split(";"):
        if "=" in token:
            k, val = token.split("=", 1)
            out[k.strip().upper()] = val.strip()
    return out


def _parse_vcf(text: str) -> tuple[str, list[Variant]]:
    sid = "sample"
    variants = []
    for line in text.splitlines():
        if line.startswith("##"):
            if line.startswith("##sampleID="):
                sid = line.split("=", 1)[1].strip()
            continue
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        chrom, pos, _id, ref, alt, _qual, _filt, info = f[:8]
        meta = _info(info)
        vaf, depth = 1.0, 0
        if len(f) >= 10:  # FORMAT + sample
            fmt = f[8].split(":")
            vals = f[9].split(":")
            d = dict(zip(fmt, vals))
            if "AF" in d:
                try:
                    vaf = float(d["AF"])
                except ValueError:
                    pass
            if "DP" in d:
                try:
                    depth = int(d["DP"])
                except ValueError:
                    pass
        variants.append(
            Variant(
                gene=meta.get("GENE", chrom),
                change=meta.get("AACHANGE", f"{ref}{pos}{alt}"),
                genomic_pos=int(pos) if pos.isdigit() else None,
                ref=ref,
                alt=alt,
                vaf=vaf,
                depth=depth,
                region=_region(meta.get("REGION")),
                silent=meta.get("EFFECT", "").lower() in ("synonymous", "silent"),
            )
        )
    return sid, variants


def load(path: str | Path, depth_floor: int = 10) -> AdapterResult:
    path = Path(path)
    text = path.read_text()
    suffix = path.suffix.lower()
    if suffix == ".vcf":
        sid, variants = _parse_vcf(text)
    elif suffix in (".tsv", ".txt"):
        sid, variants = _parse_tsv(text)
    elif suffix == ".json":
        sid, variants = _parse_json(text)
    else:
        raise ValueError(f"Unsupported input format: {suffix}")
    kept, warnings = _gate(variants, depth_floor)
    return AdapterResult(sid, kept, warnings)
