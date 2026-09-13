"""Site authentication and disclosure control for federated submissions.

Read this before deploying
--------------------------
What this module provides: authenticated, integrity-protected, replay-resistant
submission of aggregated observations between a site and a coordinator, plus a
k-anonymity gate that refuses to release cells too small to be safe.

What it does **not** provide, and must not be described as providing:

* **Secure aggregation.** The coordinator sees each site's contribution in
  clear here. A real secure-aggregation protocol (pairwise masking, threshold
  homomorphic encryption) is a cryptographic design task requiring review by a
  cryptographer and a vetted library — not something to hand-roll and call
  production-ready. A **reference implementation** of the pairwise-masking
  arithmetic now exists at ``federated/secure_aggregation.py``, built at
  explicit user request to make the idea runnable and testable — its own
  module docstring is equally explicit that it is unreviewed, does not solve
  pairwise key agreement, and has no dropout tolerance. Nothing in this
  module (``transport.py``) uses it; a submission's payload still travels to
  the coordinator in clear, as documented throughout this file.
* **Differential privacy.** No noise is added and no privacy budget is
  tracked. Adding calibrated noise to counts this small would destroy the
  signal; the correct control at this scale is the k-anonymity gate below plus
  a data-sharing agreement.
* **Federated model training.** Nothing here trains or averages model
  parameters. The federation shares aggregated *evidence*, which is a
  deliberately weaker and more auditable contract.

The residual risk that remains
------------------------------
Removing patient identifiers does not make a record anonymous. A rare variant,
combined with a site and a phenotype, can single out one patient — and rarity
is exactly what makes a variant scientifically interesting, so the risk
concentrates on the records the federation most wants. ``k_anonymity_violations``
is the mitigation, and it is a floor, not a guarantee.

HMAC-SHA256 with a per-site shared secret is used for authentication. That
requires secure secret distribution and rotation, which is an operational
matter this module cannot solve.
"""
from __future__ import annotations

import hashlib
import hmac
import json
import secrets
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import Iterable, Optional

from .catalogue_update import IsolateObservation


class TransportError(ValueError):
    """A submission could not be authenticated or accepted."""


# -- sites ----------------------------------------------------------------
@dataclass(frozen=True)
class Site:
    site_id: str
    name: str
    country: Optional[str] = None
    #: Shared secret for HMAC. Distribution and rotation are out of scope.
    secret: bytes = b""
    accredited: bool = False


class SiteRegistry:
    def __init__(self) -> None:
        self._sites: dict[str, Site] = {}

    def register(self, site: Site) -> Site:
        if not site.secret:
            raise TransportError(
                f"site {site.site_id!r} has no shared secret; unauthenticated "
                f"submissions are not accepted")
        self._sites[site.site_id] = site
        return site

    def get(self, site_id: str) -> Site:
        site = self._sites.get(site_id)
        if site is None:
            raise TransportError(f"unknown site {site_id!r}")
        return site

    @property
    def site_ids(self) -> set[str]:
        return set(self._sites)

    @staticmethod
    def new_secret() -> bytes:
        return secrets.token_bytes(32)


# -- signing --------------------------------------------------------------
def _canonical(payload: list[dict]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _mac_body(site_id: str, nonce: str, timestamp_utc: str,
              payload: list[dict]) -> bytes:
    return _canonical([
        {"site_id": site_id, "nonce": nonce, "timestamp_utc": timestamp_utc},
        *payload,
    ]).encode()


@dataclass
class Submission:
    site_id: str
    nonce: str
    timestamp_utc: str
    payload: list[dict] = field(default_factory=list)
    signature: str = ""

    @property
    def observations(self) -> list[IsolateObservation]:
        return [IsolateObservation(
            isolate_id=p["isolate_id"], site_id=p["site_id"], drug=p["drug"],
            phenotype=p.get("phenotype"),
            variant_keys=tuple(p.get("variant_keys", ())),
            lineage=p.get("lineage"), dst_method=p.get("dst_method"),
            critical_concentration=p.get("critical_concentration"),
            mic=p.get("mic"), mic_unit=p.get("mic_unit", "mg/L"),
            lab_quality=p.get("lab_quality", "unknown"),
        ) for p in self.payload]


def build_submission(site: Site,
                     observations: Iterable[IsolateObservation]) -> Submission:
    """Package and sign a site's observations."""
    payload = []
    for obs in observations:
        if obs.site_id != site.site_id:
            raise TransportError(
                f"observation for isolate {obs.isolate_id!r} carries site "
                f"{obs.site_id!r} but is being submitted by {site.site_id!r}")
        payload.append({
            "isolate_id": obs.isolate_id, "site_id": obs.site_id,
            "drug": obs.drug, "phenotype": obs.phenotype,
            "variant_keys": list(obs.variant_keys), "lineage": obs.lineage,
            "dst_method": obs.dst_method,
            "critical_concentration": obs.critical_concentration,
            "mic": obs.mic, "mic_unit": obs.mic_unit,
            "lab_quality": obs.lab_quality,
        })

    nonce = secrets.token_hex(16)
    timestamp = datetime.now(timezone.utc).isoformat()
    signature = hmac.new(
        site.secret, _mac_body(site.site_id, nonce, timestamp, payload),
        hashlib.sha256,
    ).hexdigest()
    return Submission(site.site_id, nonce, timestamp, payload, signature)


class SubmissionVerifier:
    """Verifies signature, freshness and nonce novelty."""

    def __init__(self, registry: SiteRegistry, max_age: timedelta = timedelta(hours=24)):
        self.registry = registry
        self.max_age = max_age
        self._seen_nonces: set[str] = set()

    def verify(self, submission: Submission) -> list[IsolateObservation]:
        site = self.registry.get(submission.site_id)

        expected = hmac.new(
            site.secret,
            _mac_body(submission.site_id, submission.nonce,
                      submission.timestamp_utc, submission.payload),
            hashlib.sha256,
        ).hexdigest()
        if not hmac.compare_digest(expected, submission.signature):
            raise TransportError(
                f"signature mismatch for site {submission.site_id!r}; "
                f"submission rejected")

        if submission.nonce in self._seen_nonces:
            raise TransportError(
                f"nonce {submission.nonce!r} already seen; replay rejected")

        try:
            sent = datetime.fromisoformat(submission.timestamp_utc)
        except ValueError as exc:
            raise TransportError(
                f"unparseable timestamp {submission.timestamp_utc!r}") from exc
        if sent.tzinfo is None:
            sent = sent.replace(tzinfo=timezone.utc)
        age = datetime.now(timezone.utc) - sent
        if age > self.max_age:
            raise TransportError(
                f"submission is {age} old, beyond the {self.max_age} window")
        if age < timedelta(hours=-1):
            raise TransportError("submission timestamp is in the future")

        self._seen_nonces.add(submission.nonce)
        return submission.observations


# -- disclosure control ---------------------------------------------------
@dataclass
class DisclosureThresholds:
    min_isolates_per_cell: int = 5
    min_sites_per_cell: int = 2


def k_anonymity_violations(
    observations: Iterable[IsolateObservation],
    thresholds: Optional[DisclosureThresholds] = None,
) -> list[str]:
    """Cells too small to release.

    A cell is a variant–drug pair. One below the floor should be withheld or
    suppressed before aggregation leaves the coordinator, because a rare
    variant plus a site plus a phenotype can identify a patient.
    """
    t = thresholds or DisclosureThresholds()
    cells: dict[tuple[str, str], tuple[set[str], set[str]]] = {}
    for obs in observations:
        if not obs.usable:
            continue
        for variant_key in set(obs.variant_keys):
            isolates, sites = cells.setdefault((variant_key, obs.drug),
                                               (set(), set()))
            isolates.add(obs.isolate_id)
            sites.add(obs.site_id)

    violations = []
    for (variant_key, drug), (isolates, sites) in sorted(cells.items()):
        if len(isolates) < t.min_isolates_per_cell:
            violations.append(
                f"{variant_key}/{drug}: {len(isolates)} isolate(s), below the "
                f"{t.min_isolates_per_cell}-isolate disclosure floor")
        if len(sites) < t.min_sites_per_cell:
            violations.append(
                f"{variant_key}/{drug}: {len(sites)} site(s), below the "
                f"{t.min_sites_per_cell}-site disclosure floor")
    return violations


def suppress_small_cells(
    observations: Iterable[IsolateObservation],
    thresholds: Optional[DisclosureThresholds] = None,
) -> tuple[list[IsolateObservation], list[str]]:
    """Drop observations whose only variant cells are below the floor."""
    t = thresholds or DisclosureThresholds()
    observations = list(observations)
    counts: dict[tuple[str, str], set[str]] = {}
    for obs in observations:
        if not obs.usable:
            continue
        for variant_key in set(obs.variant_keys):
            counts.setdefault((variant_key, obs.drug), set()).add(obs.isolate_id)

    kept, suppressed = [], []
    for obs in observations:
        if not obs.usable:
            continue
        safe = tuple(
            v for v in obs.variant_keys
            if len(counts.get((v, obs.drug), set())) >= t.min_isolates_per_cell
        )
        if not safe:
            suppressed.append(
                f"{obs.isolate_id}/{obs.drug}: all variant cells below the "
                f"disclosure floor")
            continue
        if len(safe) != len(obs.variant_keys):
            suppressed.append(
                f"{obs.isolate_id}/{obs.drug}: "
                f"{len(obs.variant_keys) - len(safe)} rare variant(s) withheld")
        kept.append(IsolateObservation(
            isolate_id=obs.isolate_id, site_id=obs.site_id, drug=obs.drug,
            phenotype=obs.phenotype, variant_keys=safe, lineage=obs.lineage,
            dst_method=obs.dst_method,
            critical_concentration=obs.critical_concentration,
            mic=obs.mic, mic_unit=obs.mic_unit, lab_quality=obs.lab_quality,
        ))
    return kept, suppressed
