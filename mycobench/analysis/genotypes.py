"""Per-isolate genotypes from CRyPTIC VCFs — the link phenotypes need.

The reuse table pairs an isolate with its MICs and with a *path* to its VCF; it
does not carry variant calls. Without those calls there is nothing to stratify,
so this module is what makes the effect-size analysis runnable rather than
theoretical.

These files are large, and planning around a wrong figure wastes a lot of time
-----------------------------------------------------------------------------
The release's re-genotyped VCFs are **~20 MB compressed and ~178 MB
decompressed each**, about 1.27 million records, because every callable site in
the genome is present including the ``0/0`` reference calls. Measured, not
assumed. The full 12,287-isolate compendium is therefore roughly **248 GB**,
not the fraction of a gigabyte a small-file assumption suggests, and caching
all of it is an infrastructure decision rather than an incidental download.

Those reference calls are the reason the files are worth their size: a
catalogued position called ``0/0`` is positive evidence that the locus was
examined and the variant was absent, which is what ``assessed_variants``
carries and what lets a drug reach SUSCEPTIBLE instead of NOT_ASSESSED.

Genotyping is restricted to catalogued positions, deliberately
-------------------------------------------------------------
These VCFs carry coordinates, not gene annotations. Rather than bolt on a
variant annotator — another unvalidated component — this module inverts the
WHO catalogue's own coordinate file into a ``(position, ref, alt) -> variant``
index and reports which *catalogued* variants an isolate carries.

That is a real limitation and it is worth stating plainly: a variant absent
from the catalogue is invisible here, so this cannot discover entirely novel
determinants. It is nonetheless exactly what the immediate task needs — the
20,843 entries the real catalogue grades as "Uncertain significance" all have
coordinates, and evaluating those is the VUS gap.

Matching ignores chromosome naming
----------------------------------
The catalogue writes ``NC_000962.3`` while a VCF may write ``Chromosome``,
``NC_000962.3`` or ``MTB_anc``. There is one chromosome in *M. tuberculosis*,
so the index keys on position and alleles and records the observed contig names
for the report rather than failing to match on a naming convention.
"""
from __future__ import annotations

import gzip
import json
import os
import random
import urllib.error
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional, Sequence

from .. import USER_AGENT
from ..phenotypes import FTP_BASE
from .strata import Isolate

#: The reuse table's VCF paths are relative to the release root, one level
#: above the reuse directory.
RELEASE_BASE = FTP_BASE.rsplit("/", 1)[0]

#: FILTER values that still count as a call. CRyPTIC's masked VCFs use PASS
#: and a set of mask reasons; anything else is treated as no call.
ACCEPTED_FILTERS = frozenset({"PASS", ".", ""})


class GenotypeError(RuntimeError):
    """A VCF could not be fetched or parsed."""


def _trim(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Left-align and trim, so the two sources agree on one spelling."""
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


@dataclass
class CoordinateIndex:
    """``(position, ref, alt) -> catalogued variant label(s)``.

    Built from the ingested catalogue, which retains every coordinate spelling
    WHO publishes for a variant — so an isolate's SNV representation and the
    catalogue's MNV representation both resolve to the same label.
    """

    by_allele: dict[tuple[int, str, str], set[str]] = field(default_factory=dict)
    #: Every position any catalogued key occupies. Used to skip the ~98% of
    #: records in a re-genotyped VCF that cannot match anything, without
    #: parsing them. See ``parse_vcf`` for why this is exact rather than
    #: approximate.
    positions: set[int] = field(default_factory=set)
    catalogue_version: str = ""
    n_variants: int = 0
    n_keys: int = 0

    def __post_init__(self) -> None:
        """Derive ``positions`` from ``by_allele`` when it was not supplied.

        ``parse_vcf`` skips records whose position is absent from this set, so
        an index carrying allele keys but an empty position set would match
        **nothing at all** — silently, with no error and an empty genotype.
        Callers that assemble ``by_allele`` directly rather than through
        ``from_catalogue`` must not have to know that a second structure exists
        to keep in step. ``from_catalogue`` fills ``by_allele`` after
        construction and maintains both as it goes, so nothing is derived here.
        """
        if self.by_allele and not self.positions:
            self.positions = {position for (position, _ref, _alt)
                              in self.by_allele}

    @classmethod
    def from_catalogue(cls, path: str | Path) -> "CoordinateIndex":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        index = cls(catalogue_version=data.get("catalogue_version", ""))
        for entry in data.get("variants", []):
            keys = entry.get("coordinate_keys")
            if not keys:
                key = entry.get("coordinate_key")
                keys = [key] if key else []
            if not keys:
                continue
            label = f"{entry['gene']}_{entry['change']}"
            index.n_variants += 1
            for raw in keys:
                parsed = cls._parse_key(raw)
                if parsed is None:
                    continue
                index.by_allele.setdefault(parsed, set()).add(label)
                index.positions.add(parsed[0])
                # Also index the trimmed spelling, so an MNV in the catalogue
                # still matches an isolate's minimal representation.
                trimmed = _trim(*parsed)
                if trimmed != parsed:
                    index.by_allele.setdefault(trimmed, set()).add(label)
                    index.positions.add(trimmed[0])
        index.n_keys = len(index.by_allele)
        return index

    @staticmethod
    def _parse_key(key: str) -> Optional[tuple[int, str, str]]:
        """``assembly:contig:pos:ref>alt`` -> ``(pos, ref, alt)``."""
        try:
            _, _, position, alleles = key.split(":", 3)
            ref, alt = alleles.split(">", 1)
            return int(position), ref.upper(), alt.upper()
        except (ValueError, AttributeError):
            return None

    def lookup(self, position: int, ref: str, alt: str) -> set[str]:
        """Catalogued labels for one allele, under either spelling.

        Called once per alternate allele, and a re-genotyped record carries a
        median of thirteen, so this is the hottest path in the whole load. Two
        shortcuts keep it honest and cheap:

        ``_trim`` cannot change anything unless **both** sides are multi-base —
        each of its loops requires ``len(ref) > 1 and len(alt) > 1`` — so an
        ordinary SNV skips it outright rather than paying for a call that is
        guaranteed to return its argument. And a miss returns without building
        a set, because the overwhelming majority of alleles match nothing.
        """
        ref, alt = ref.upper(), alt.upper()
        found = self.by_allele.get((position, ref, alt))
        if len(ref) > 1 and len(alt) > 1:
            trimmed = _trim(position, ref, alt)
            if trimmed != (position, ref, alt):
                shifted = self.by_allele.get(trimmed)
                if shifted:
                    return set(shifted) | found if found else set(shifted)
        return set(found) if found else set()

    def describe(self) -> str:
        return (f"{self.n_variants} catalogued variant(s) with coordinates, "
                f"{self.n_keys} allele key(s) indexed "
                f"(catalogue {self.catalogue_version})")


# -- fetching -------------------------------------------------------------
def vcf_url(relative_path: str) -> str:
    """Absolute URL for a reuse-table VCF path."""
    cleaned = relative_path.strip()
    while cleaned.startswith("../"):
        cleaned = cleaned[3:]
    return f"{RELEASE_BASE}/{cleaned}"


def fetch_vcf(relative_path: str, cache_dir: str | Path,
              timeout: int = 120, cached_only: bool = False) -> Path:
    """Download one VCF into a local cache, or reuse what is already there.

    Caching is by the release-relative path, so a re-run of the analysis costs
    nothing and a partial run resumes.
    """
    cache_dir = Path(cache_dir)
    cleaned = relative_path.strip()
    while cleaned.startswith("../"):
        cleaned = cleaned[3:]
    destination = (cache_dir / cleaned).resolve()
    if not destination.is_relative_to(cache_dir.resolve()):
        raise GenotypeError("VCF path escapes cache directory")
    if destination.is_file() and destination.stat().st_size > 0:
        return destination
    if cached_only:
        raise GenotypeError(f"not cached: {cleaned}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(vcf_url(relative_path),
                                     headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise GenotypeError(
            f"could not fetch {vcf_url(relative_path)}: {exc}") from exc
    if not payload:
        raise GenotypeError(f"{vcf_url(relative_path)} returned no bytes")
    # Write to a temporary file and rename, rather than writing the final path
    # directly. These are ~20 MB, so the write is not instantaneous, and the
    # cache is trusted on nothing more than "the file exists and is non-empty".
    # A process killed mid-write would otherwise leave a truncated VCF at the
    # canonical path, which every later run would accept and fail to parse --
    # permanently, since nothing invalidates it. os.replace is atomic on POSIX
    # and Windows, so a reader sees either no file or the whole file. It also
    # makes concurrent fetches of one path safe, which matters now that
    # prefetch_vcfs runs several workers at once.
    temporary = destination.with_name(f"{destination.name}.{os.getpid()}.part")
    try:
        temporary.write_bytes(payload)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


# -- parsing --------------------------------------------------------------
@dataclass
class VCFGenotype:
    variants: frozenset[str] = frozenset()
    #: Alternate alleles carried **at coordinates the catalogue could match**,
    #: not genome-wide. The genome-wide count is in ``n_sites``. This is the
    #: denominator that makes ``match_rate`` diagnostic: a low rate here means
    #: the right positions carried alleles the catalogue does not list, which is
    #: what a coordinate-system or reference mismatch looks like. A genome-wide
    #: denominator would bury that under the trivial fact that most variation
    #: falls outside resistance loci.
    n_records: int = 0
    n_matched: int = 0
    n_filtered_out: int = 0
    #: Every data record in the file, including reference calls.
    n_sites: int = 0
    contigs: frozenset[str] = frozenset()
    assessed_variants: frozenset[str] = frozenset()

    @property
    def match_rate(self) -> float:
        return self.n_matched / self.n_records if self.n_records else 0.0


def parse_vcf(path: str | Path, index: CoordinateIndex) -> VCFGenotype:
    """Read a (gzipped) VCF and report which catalogued variants it carries.

    These are re-genotyped VCFs: every callable site in the genome appears,
    including the ``0/0`` reference calls, so one file is ~1.27 million records
    and ~178 MB decompressed. Those reference calls are the point — a catalogued
    position called ``0/0`` is positive evidence that the locus was examined and
    the variant was absent, which is what ``assessed_variants`` carries and what
    lets a drug reach SUSCEPTIBLE at all. They are not skipped.

    What is skipped is the ~98% of records that cannot match the catalogue
    under any spelling. The test is exact, not heuristic: ``_trim`` only moves a
    record's position while ``len(ref) > 1``, so a **single-base REF cannot
    shift**, and for those an exact position match is provably sufficient. Any
    record with a multi-base REF is parsed in full regardless of position,
    because trimming could carry it onto a catalogued coordinate. Dropping the
    position test entirely would be correct too — just twelve times slower.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    variants: set[str] = set()
    assessed: set[str] = set()
    contigs: set[str] = set()
    sites = records = matched = filtered = 0
    positions = index.positions
    # Contig is read by comparing a prefix rather than slicing every line: a
    # VCF is sorted by contig, so this allocates once per contig instead of
    # 1.27 million times, while still noticing a second contig if one appears.
    seen_prefix = "\x00"

    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line or line[0] == "#":
                    continue
                # Locate the first four tabs without allocating substrings.
                # Only records that survive the position test below are worth
                # the cost of splitting into fields.
                tab1 = line.find("\t")
                if tab1 < 0:
                    continue
                if not line.startswith(seen_prefix):
                    chrom = line[:tab1]
                    seen_prefix = chrom + "\t"
                    contigs.add(chrom)
                tab2 = line.find("\t", tab1 + 1)
                tab3 = line.find("\t", tab2 + 1)
                tab4 = line.find("\t", tab3 + 1)
                if tab2 < 0 or tab3 < 0 or tab4 < 0:
                    continue
                sites += 1

                # Reject what cannot possibly match, before paying for a split
                # and a median of thirteen per-allele lookups.
                #
                # A single-base REF cannot be shifted by _trim at all, so an
                # exact position match decides it. A multi-base REF can be
                # left-trimmed forward by at most len(REF) - 1, so it can only
                # reach positions in [P, P + len(REF) - 1]; if the catalogue
                # holds none of those, no spelling of this record matches.
                # Both tests are exact: neither can discard a real match.
                position = line[tab1 + 1:tab2]
                if not position.isdigit():
                    continue
                start = int(position)
                if tab4 - tab3 == 2:
                    if start not in positions:
                        continue
                elif not any(p in positions
                             for p in range(start, start + tab4 - tab3 - 1)):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                _chrom, position, _id, ref, alt_field = fields[:5]
                if not position.isdigit():
                    continue
                filter_value = fields[6] if len(fields) > 6 else ""
                filters = {f for f in filter_value.split(";") if f}
                if filters and not filters <= ACCEPTED_FILTERS:
                    filtered += 1
                    continue

                alleles = None
                if len(fields) > 10:
                    raise GenotypeError("cohort VCF must contain exactly one sample")
                if len(fields) == 10 and "GT" in fields[8].split(":"):
                    fmt = dict(zip(fields[8].split(":"), fields[9].split(":")))
                    gt = fmt.get("GT", ".").replace("|", "/").split("/")
                    if not gt or any(not a.isdigit() for a in gt):
                        filtered += 1
                        continue
                    alleles = {int(a) for a in gt}
                    if max(alleles) > len(alt_field.split(",")):
                        raise GenotypeError("VCF genotype allele index exceeds ALT count")
                for allele_index, alt in enumerate(alt_field.split(","), 1):
                    alt = alt.strip()
                    # A symbolic or absent ALT carries no allele to match.
                    if not alt or alt in (".", "<NON_REF>", "*") \
                            or alt.startswith("<"):
                        continue
                    found = index.lookup(int(position), ref, alt)
                    assessed |= found
                    if alleles is not None and allele_index not in alleles:
                        continue
                    records += 1
                    if found:
                        matched += 1
                        variants |= found
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise GenotypeError(f"could not parse {path}: {exc}") from exc

    return VCFGenotype(variants=frozenset(variants), n_records=records,
                       n_matched=matched, n_filtered_out=filtered,
                       n_sites=sites,
                       contigs=frozenset(contigs), assessed_variants=frozenset(assessed))


# -- driver ---------------------------------------------------------------
@dataclass
class GenotypeLoad:
    isolates: dict[str, Isolate] = field(default_factory=dict)
    failures: list[str] = field(default_factory=list)
    n_requested: int = 0
    n_records: int = 0
    n_matched: int = 0
    n_sites: int = 0
    contigs: set[str] = field(default_factory=set)

    @property
    def n_loaded(self) -> int:
        return len(self.isolates)

    @property
    def match_rate(self) -> float:
        return self.n_matched / self.n_records if self.n_records else 0.0

    @property
    def carriers_per_variant(self) -> dict[str, int]:
        counts: dict[str, int] = defaultdict(int)
        for isolate in self.isolates.values():
            for variant in isolate.genotype:
                counts[variant] += 1
        return dict(counts)

    def describe(self) -> str:
        lines = [
            f"{self.n_loaded}/{self.n_requested} isolate(s) genotyped, "
            f"{len(self.failures)} failure(s)",
            f"{self.n_matched}/{self.n_records} alternate allele(s) at "
            f"catalogued coordinates matched a catalogued allele "
            f"({self.match_rate:.1%}); {self.n_sites} site(s) read",
        ]
        counts = self.carriers_per_variant
        if counts:
            lines.append(f"{len(counts)} distinct catalogued variant(s) observed")
        if self.contigs:
            lines.append(f"contig name(s) seen: {', '.join(sorted(self.contigs))}")
        return "\n".join(lines)


def prefetch_vcfs(paths: Sequence[str], cache_dir: str | Path, jobs: int = 4,
                  progress_every: int = 50) -> int:
    """Populate the cache concurrently. Returns how many were downloaded.

    Fetching is network-bound and one file takes several seconds, so a serial
    load spends most of its wall clock waiting. Parsing stays serial — it is
    CPU-bound and threads would not help it — so this only moves the waiting.

    Failures are **not** raised. A path that cannot be fetched here is simply
    left out of the cache, and the serial pass that follows reports it as that
    isolate's failure with the same message it would have produced anyway. That
    keeps one transient FTP error from discarding a partial load of thousands.
    """
    cache_dir = Path(cache_dir)
    missing = []
    for path in paths:
        try:
            fetch_vcf(path, cache_dir, cached_only=True)
        except GenotypeError:
            missing.append(path)
    if not missing:
        return 0

    if progress_every:
        print(f"[genotypes] fetching {len(missing)} VCF(s) with "
              f"{jobs} worker(s)", flush=True)
    done = 0
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(fetch_vcf, path, cache_dir): path
                   for path in missing}
        for future in as_completed(futures):
            done += 1
            try:
                future.result()
            except GenotypeError:
                pass  # reported by the serial pass, against its isolate
            if progress_every and done % progress_every == 0:
                print(f"[genotypes] fetched {done}/{len(missing)}", flush=True)
    return done


def load_genotypes(rows: Iterable[dict], index: CoordinateIndex,
                   cache_dir: str | Path,
                   limit: Optional[int] = None,
                   prefer_regenotyped: bool = True,
                   progress_every: int = 250, cached_only: bool = False,
                   jobs: int = 1,
                   sample_seed: Optional[int] = None) -> GenotypeLoad:
    """Genotype isolates from reuse-table rows.

    ``rows`` are dicts from the CRyPTIC reuse table, needing ``ENA_RUN`` and a
    VCF path. A row whose VCF cannot be fetched or parsed is **recorded and
    skipped**, never fatal: one transient FTP error must not discard a partial
    load of thousands.

    ``prefer_regenotyped`` uses the regenotyped VCF where present, which is the
    consortium's own reconciled call set.

    ``jobs`` above one downloads the selected VCFs concurrently before parsing
    any of them. It changes nothing about the result, only how long it waits.

    ``sample_seed`` draws the ``limit`` isolates at random instead of taking
    the first ones. The reuse table is ordered by site and subject, so a plain
    prefix is drawn from one or two laboratories: every isolate then shares a
    site, lineage diversity collapses, and the site and lineage confounding
    checks silently have nothing to compare. Sampling is seeded so a run stays
    reproducible.
    """
    result = GenotypeLoad()
    eligible: list[tuple[str, str, dict]] = []
    for row in rows:
        if limit is not None and sample_seed is None and len(eligible) >= limit:
            break
        run = (row.get("ENA_RUN") or "").strip()
        if not run or run.upper() == "NONE":
            continue
        path = ((row.get("REGENOTYPED_VCF") if prefer_regenotyped else "")
                or row.get("VCF") or "").strip()
        if not path:
            result.failures.append(f"{run}: no VCF path in the reuse table")
            continue
        eligible.append((run, path, row))

    if sample_seed is not None and limit is not None and len(eligible) > limit:
        # Sort first: the draw must not depend on the order rows arrived in.
        eligible.sort(key=lambda item: item[0])
        eligible = random.Random(sample_seed).sample(eligible, limit)
        if progress_every:
            print(f"[genotypes] sampled {limit} isolate(s) with seed "
                  f"{sample_seed}", flush=True)
    selected = eligible

    if jobs > 1 and not cached_only and selected:
        prefetch_vcfs([path for _run, path, _row in selected], cache_dir,
                      jobs=jobs)

    for run, path, row in selected:
        result.n_requested += 1
        try:
            local = fetch_vcf(path, cache_dir, cached_only=cached_only)
            genotype = parse_vcf(local, index)
        except GenotypeError as exc:
            result.failures.append(str(exc))
            continue

        result.n_records += genotype.n_records
        result.n_matched += genotype.n_matched
        result.n_sites += genotype.n_sites
        result.contigs |= set(genotype.contigs)
        result.isolates[run] = Isolate(
            isolate_id=f"cr_{run}", genotype=genotype.variants,
            assessed_variants=genotype.assessed_variants,
            lineage=(row.get("LINEAGE") or None),
            site=(row.get("UNIQUEID") or "").split(".")[1]
                 if (row.get("UNIQUEID") or "").startswith("site.") else None)
        if progress_every and result.n_requested % progress_every == 0:
            print(f"[genotypes] {result.n_requested} isolate(s): "
                  f"{result.match_rate:.1%} of alleles catalogued", flush=True)
    return result
