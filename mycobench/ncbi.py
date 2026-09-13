"""NCBI E-utilities client: paced, retried, and honest about what it proves.

Two lessons are baked in.

**Pace the failure path too.** NCBI allows 3 requests/second without an API key
and 10 with one. A cohort walk costs several requests per row, so callers must
pace themselves rather than rely on 429 retries to absorb the overrun — and an
unpaced *error* path is exactly what turns transient rate limiting into a
cascade.

**A free-text match is not a statement of origin.** Searching SRA for
``Nigeria[All Fields]`` returns runs whose deposited ``geo_loc_name`` is
``China: East China Sea, Xiangshan Bay`` — the word appears elsewhere in the
record. ``biosample_attributes`` exists so discovery can confirm origin against
the controlled field instead, and an absent ``geo_loc_name`` is treated as
unconfirmed rather than as a pass.
"""
from __future__ import annotations

import csv
import http.client
import json
import os
import random
import re
import time
from dataclasses import dataclass, field
from typing import Iterable, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from . import USER_AGENT

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

RETRYABLE_STATUS = (429, 500, 502, 503, 504)
#: A dropped keep-alive surfaces as http.client.RemoteDisconnected, which is
#: not a URLError; without it here one flaky connection discards a whole walk.
TRANSIENT = (URLError, TimeoutError, ConnectionError, http.client.HTTPException,
             json.JSONDecodeError, UnicodeDecodeError)

RUN_ACCESSION = re.compile(r"^(SRR|ERR|DRR)\d+$")
BIOSAMPLE_ACCESSION = re.compile(r"^SAM[NED]A?\d+$")
BIOPROJECT_ACCESSION = re.compile(r"^PRJ[NED][AB]\d+$")


class NCBIError(RuntimeError):
    """A request failed after its retries, or returned something unusable."""


def request_interval(api_key: Optional[str] = None) -> float:
    """Seconds to pause between per-record request bursts."""
    return 0.11 if api_key else 0.34


@dataclass
class Credentials:
    email: Optional[str] = None
    api_key: Optional[str] = None

    @classmethod
    def from_env(cls, email: Optional[str] = None,
                 api_key: Optional[str] = None) -> "Credentials":
        return cls(email=email or os.environ.get("NCBI_EMAIL"),
                   api_key=api_key or os.environ.get("NCBI_API_KEY"))

    @property
    def interval(self) -> float:
        return request_interval(self.api_key)

    def params(self) -> dict:
        out = {}
        if self.email:
            out["email"] = self.email
        if self.api_key:
            out["api_key"] = self.api_key
        return out

    def describe(self) -> str:
        if self.api_key:
            return "NCBI API key present (10 requests/second allowance)"
        return ("no NCBI API key (3 requests/second); set NCBI_API_KEY to walk "
                "large panels faster")


def fetch(request: Request, timeout: int = 60, retries: int = 4,
          label: str = "NCBI request", parse=None):
    """Read a URL with bounded retry across HTTP and connection failures."""
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=timeout) as response:
                payload = response.read()
            return parse(payload) if parse else payload
        except (HTTPError, *TRANSIENT) as exc:
            fatal = isinstance(exc, HTTPError) and exc.code not in RETRYABLE_STATUS
            if attempt == retries - 1 or fatal:
                raise NCBIError(
                    f"{label} failed after {attempt + 1} attempt(s): {exc}") from exc
            retry_after = (exc.headers.get("Retry-After")
                           if isinstance(exc, HTTPError) else None)
            time.sleep(float(retry_after) if retry_after and retry_after.isdigit()
                       else min(30, 2 ** attempt + random.random()))


def _get(endpoint: str, params: dict, creds: Credentials, timeout: int = 60,
         label: str = "NCBI request", parse=None):
    request = Request(f"{endpoint}?{urlencode({**params, **creds.params()})}",
                      headers={"User-Agent": USER_AGENT})
    return fetch(request, timeout, 4, label, parse=parse)


def esearch(db: str, term: str, creds: Credentials, retmax: int = 500) -> dict:
    """Search one database, returning the history handle for a later efetch."""
    result = _get(f"{EUTILS}/esearch.fcgi",
                  {"db": db, "term": term, "retmax": retmax,
                   "usehistory": "y", "retmode": "json"},
                  creds, label=f"esearch {db}", parse=json.loads)
    payload = result.get("esearchresult") or {}
    if "webenv" not in payload:
        raise NCBIError(f"esearch on {db} returned no history handle for {term!r}")
    return {"count": int(payload.get("count", 0)),
            "webenv": payload["webenv"],
            "query_key": payload.get("querykey") or payload.get("query_key"),
            "ids": payload.get("idlist", [])}


def runinfo(db_search: dict, creds: Credentials,
            retmax: int = 500) -> list[dict]:
    """Fetch SRA run-info rows for a previous esearch, as parsed CSV."""
    text = _get(f"{EUTILS}/efetch.fcgi",
                {"db": "sra", "WebEnv": db_search["webenv"],
                 "query_key": db_search["query_key"],
                 "rettype": "runinfo", "retmode": "text", "retmax": retmax},
                creds, timeout=300, label="efetch runinfo",
                parse=lambda payload: payload.decode("utf-8", "replace"))
    rows = [row for row in csv.DictReader(text.splitlines()) if row.get("Run")]
    if not rows:
        raise NCBIError("run-info returned no rows")
    return rows


def runinfo_for_accessions(accessions: Iterable[str],
                           creds: Credentials) -> dict[str, dict]:
    """Run-info keyed by run accession, for an explicit list of runs."""
    accessions = [a for a in accessions if a]
    out: dict[str, dict] = {}
    for start in range(0, len(accessions), 100):
        batch = accessions[start:start + 100]
        term = " OR ".join(f"{a}[Accession]" for a in batch)
        search = esearch("sra", term, creds, retmax=len(batch) * 2)
        time.sleep(creds.interval)
        for row in runinfo(search, creds, retmax=len(batch) * 2):
            out[row["Run"]] = row
        time.sleep(creds.interval)
    return out


@dataclass
class BioSampleRecord:
    accession: str
    attributes: dict[str, str] = field(default_factory=dict)

    @property
    def geo_loc_name(self) -> str:
        return self.attributes.get("geo_loc_name", "")

    @property
    def collection_date(self) -> str:
        return self.attributes.get("collection_date", "")

    @property
    def host(self) -> str:
        return self.attributes.get("host", "")

    @property
    def isolation_source(self) -> str:
        return self.attributes.get("isolation_source", "")

    def confirms_country(self, country: str) -> bool:
        """Does the controlled geo_loc_name field confirm this country?

        Requires the value to *begin* with the country, so ``Nigeria`` and
        ``Nigeria:BENUE`` both pass while a country mentioned incidentally
        elsewhere in the record does not. An absent value never passes:
        unconfirmed is not confirmed.
        """
        value = self.geo_loc_name.strip()
        if not value:
            return False
        return value.lower().startswith(country.strip().lower())

    def phenotype_attributes(self) -> dict[str, str]:
        """Any deposited attribute that looks like a DST or MIC result.

        Public mycobacterial BioSamples very rarely carry one; this exists so a
        cohort can record that it looked, rather than assuming.
        """
        pattern = re.compile(
            r"(?i)suscept|resistan|\bdst\b|\bmic\b|phenotyp|"
            r"rifamp|isoniaz|ethambut|pyrazinam|fluoroquinolon")
        return {k: v for k, v in self.attributes.items() if pattern.search(k)}


_BIOSAMPLE_BLOCK = re.compile(r'(?s)<BioSample[^>]*?accession="([^"]+)".*?</BioSample>')
_ATTRIBUTE = re.compile(r'attribute_name="([^"]+)"[^>]*>([^<]*)')


def biosample_attributes(accessions: Iterable[str], creds: Credentials,
                         batch_size: int = 100) -> dict[str, BioSampleRecord]:
    """Resolve deposited BioSample attributes, batched and paced.

    Accessions are resolved to UIDs first. ``efetch`` on ``db=biosample``
    wants UIDs: passing an accession happens to work for ``SAMN`` records but
    returns an **empty document** — not an error — for ``SAMEA`` and ``SAMD``
    ones deposited through ENA and DDBJ. Skipping the esearch step therefore
    produces a silently empty result that reads as "this sample has no
    metadata" when it means "the lookup was malformed".

    Parsed with a regex rather than an XML parser on purpose: the payload is a
    concatenation of independent BioSample documents, and a malformed record in
    one must not discard the whole batch.
    """
    accessions = sorted({a for a in accessions if a})
    out: dict[str, BioSampleRecord] = {}
    for start in range(0, len(accessions), batch_size):
        batch = accessions[start:start + batch_size]
        search = esearch("biosample",
                         " OR ".join(f"{a}[Accession]" for a in batch),
                         creds, retmax=len(batch) * 2)
        time.sleep(creds.interval)
        uids = search.get("ids") or []
        if not uids:
            continue
        text = _get(f"{EUTILS}/efetch.fcgi",
                    {"db": "biosample", "id": ",".join(uids), "retmode": "xml"},
                    creds, timeout=300, label="efetch biosample",
                    parse=lambda payload: payload.decode("utf-8", "replace"))
        for match in _BIOSAMPLE_BLOCK.finditer(text):
            accession = match.group(1)
            attributes = {name: value.strip()
                          for name, value in _ATTRIBUTE.findall(match.group(0))}
            out[accession] = BioSampleRecord(accession, attributes)
        time.sleep(creds.interval)
    return out


def paired_illumina(rows: Iterable[dict],
                    exclude_organisms: tuple[str, ...] = ()) -> list[dict]:
    """Keep only paired-end Illumina runs that carry a BioSample.

    A run with no BioSample cannot have its origin confirmed, so it is dropped
    here rather than carried forward unverifiable.
    """
    kept = []
    for row in rows:
        if (row.get("Platform") or "").upper() != "ILLUMINA":
            continue
        if (row.get("LibraryLayout") or "").upper() != "PAIRED":
            continue
        if not row.get("BioSample"):
            continue
        organism = (row.get("ScientificName") or "")
        if any(token.lower() in organism.lower() for token in exclude_organisms):
            continue
        kept.append(row)
    return kept


def bases_of(row: dict) -> int:
    try:
        return int(row.get("bases") or 0)
    except (TypeError, ValueError):
        return 0
