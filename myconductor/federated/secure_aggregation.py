"""Reference implementation of pairwise-masking secure aggregation.

READ THIS BEFORE USE
---------------------
This is a **reference implementation of the masking arithmetic** at the core
of Bonawitz et al. 2017, "Practical Secure Aggregation for Privacy-Preserving
Machine Learning" (ACM CCS) — written, at explicit user request, to make the
idea concrete and testable in this codebase. It is **not** a reviewed,
production-grade secure-aggregation library. `federated/transport.py`
already states that real secure aggregation "is a cryptographic design task
requiring review by a cryptographer and a vetted library — not something to
hand-roll here." This module hand-rolls the arithmetic anyway, gated behind
``acknowledged=True`` exactly the way ``modules/features.py::DemoAnnotator``
gates its own synthetic output — its presence here means "the idea made
runnable and testable," not "safe to deploy."

Three things this module does **not** solve, and must not be read as solving:

1. **Pairwise key agreement.** Correct masking requires every pair of sites
   to share a secret the *coordinator does not know*. In the published
   protocol that comes from authenticated Diffie-Hellman (e.g. X25519)
   directly between the two sites. This module does not perform that
   exchange — it takes already-established pairwise secrets as input.
   :func:`derive_test_secrets` below produces secrets the coordinator
   *could* reconstruct (it derives them from a seed the coordinator holds),
   which is fine for exercising the arithmetic in a test and unsafe for
   anything else. A real deployment needs a vetted DH/KEM implementation —
   this project has not added one, and doing so is exactly the kind of
   choice this module's own existence should not be read as pre-empting.

2. **Per-site secret confidentiality.** Each site must hold only the
   pairwise secrets it is actually a party to (see
   :func:`secrets_for_site`), never a table of every pair in the
   federation — a party holding all pairwise secrets could unmask
   everyone's individual contribution, which defeats the entire point.

3. **Dropout tolerance.** The published protocol's practical value is a
   Shamir-secret-sharing recovery path so the sum is still computable when
   clients disappear mid-round, without ever unmasking an individual
   contribution. That machinery is not implemented here.
   :func:`unmask_sum` only produces a correct result when every site that
   contributed a mask against another site also submitted its own masked
   vector in the same round; there is no way to detect a silent dropout
   from this module alone.

What it does get right
-----------------------
The masking arithmetic: pairwise additive masks over a fixed modulus,
expanded from a shared secret via HMAC-based PRG (a standard, reviewable
construction, not invented here), which cancel out exactly in the sum
regardless of contribution order — so a coordinator that only ever sees
masked per-site vectors can recover the true aggregate and nothing else,
*given* points 1-3 above are actually handled by the deployment.
"""
from __future__ import annotations

import hashlib
import hmac
from dataclasses import dataclass, field
from typing import Iterable, Optional

#: Counts in this codebase's federated context (isolate tallies) are small;
#: this modulus gives enormous headroom without needing bignum tuning.
DEFAULT_MODULUS = 1 << 64


class SecureAggregationError(ValueError):
    """A masking/unmasking precondition was violated."""


@dataclass(frozen=True)
class PairwiseSecret:
    """A secret shared between exactly two sites. See point 1 above for how
    this must actually be established in a real deployment (not here)."""

    site_a: str
    site_b: str
    secret: bytes

    @property
    def pair(self) -> tuple[str, str]:
        return tuple(sorted((self.site_a, self.site_b)))


def derive_test_secrets(site_ids: Iterable[str],
                        coordinator_seed: bytes) -> list[PairwiseSecret]:
    """FOR TESTING THE ARITHMETIC ONLY.

    These secrets are derived from ``coordinator_seed`` alone, so anyone
    holding that seed — including, by construction, whoever runs this
    function — can reconstruct every pairwise secret and unmask every site's
    individual contribution. That is the opposite of what secure aggregation
    is for. Use this only to exercise :class:`PairwiseMaskingAggregator` and
    :func:`unmask_sum` end to end in a test; never in anything resembling a
    real multi-site deployment.
    """
    ids = sorted(set(site_ids))
    out = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            seed = hmac.new(coordinator_seed, f"{a}|{b}".encode("utf-8"),
                           hashlib.sha256).digest()
            out.append(PairwiseSecret(a, b, seed))
    return out


def secrets_for_site(site_id: str,
                     pairwise: Iterable[PairwiseSecret]) -> dict[str, bytes]:
    """The subset of pairwise secrets ``site_id`` is actually a party to.

    A real site must never hold, receive, or be constructed with any other
    site's pairwise secrets — only this. This helper exists to build that
    correctly-scoped view from a full set (e.g. from
    :func:`derive_test_secrets` in a test that simulates every site in one
    process); it is not a substitute for keeping secrets separated by
    process/site boundary in a real deployment.
    """
    out: dict[str, bytes] = {}
    for ps in pairwise:
        if ps.site_a == site_id:
            out[ps.site_b] = ps.secret
        elif ps.site_b == site_id:
            out[ps.site_a] = ps.secret
    return out


def _prg(seed: bytes, length: int, modulus: int) -> list[int]:
    """Expand ``seed`` into ``length`` pseudorandom integers mod ``modulus``,
    via HMAC-SHA256 in counter mode — a standard PRG construction, not
    invented here."""
    out: list[int] = []
    counter = 0
    while len(out) < length:
        block = hmac.new(seed, counter.to_bytes(8, "big"),
                         hashlib.sha256).digest()
        for i in range(0, len(block), 8):
            if len(out) >= length:
                break
            out.append(int.from_bytes(block[i:i + 8], "big") % modulus)
        counter += 1
    return out


@dataclass
class PairwiseMaskingAggregator:
    """One site's role in one round of pairwise-masking secure aggregation.

    ``pairwise_secrets`` must contain exactly this site's own secrets (see
    :func:`secrets_for_site`) — one entry per other participating site,
    keyed by that site's id. Passing every pair's secret to every site
    reconstructs the exact failure mode this protocol exists to prevent.
    """

    site_id: str
    other_site_ids: list[str]
    pairwise_secrets: dict[str, bytes]
    modulus: int = DEFAULT_MODULUS
    acknowledged: bool = False

    def __post_init__(self) -> None:
        if not self.acknowledged:
            raise SecureAggregationError(
                "PairwiseMaskingAggregator is an unreviewed reference "
                "implementation of a published protocol, not a vetted "
                "secure-aggregation library (see this module's docstring). "
                "Pass acknowledged=True to confirm you understand that "
                "before using it."
            )
        missing = [s for s in self.other_site_ids
                  if s not in self.pairwise_secrets]
        if missing:
            raise SecureAggregationError(
                f"{self.site_id} is missing a pairwise secret for: "
                f"{missing}. Every pair of participating sites needs one "
                f"established before masking (see point 1 in this module's "
                f"docstring)."
            )

    def mask(self, values: list[int]) -> list[int]:
        """Mask this site's integer vector for one round.

        Adds (or subtracts, by lexicographic site-id order) a PRG-expanded
        pairwise mask for every other participating site, so the masks
        cancel exactly once every site's masked vector is summed
        (:func:`unmask_sum`) — provided every site that starts the round
        also finishes it (see point 3 in the module docstring).
        """
        masked = list(values)
        for other in self.other_site_ids:
            secret = self.pairwise_secrets[other]
            mask_vec = _prg(secret, len(values), self.modulus)
            sign = 1 if self.site_id < other else -1
            for i in range(len(values)):
                masked[i] = (masked[i] + sign * mask_vec[i]) % self.modulus
        return masked


def unmask_sum(masked_vectors: Iterable[list[int]],
              modulus: int = DEFAULT_MODULUS) -> list[int]:
    """The coordinator's side: sum every site's masked vector.

    Pairwise masks cancel to zero **only if** every site that masked against
    another site also contributed its own masked vector to this same round
    (see point 3 in the module docstring — there is no dropout tolerance
    here). If a site drops out after masking began, the result is wrong
    silently, not flagged: that gap is exactly what the published protocol's
    Shamir-secret-sharing recovery path exists to close, and it is not
    implemented here.
    """
    vectors = list(masked_vectors)
    if not vectors:
        raise SecureAggregationError("no vectors to aggregate")
    length = len(vectors[0])
    if any(len(v) != length for v in vectors):
        raise SecureAggregationError(
            "every site's vector must be the same length")
    total = [0] * length
    for vector in vectors:
        for i, value in enumerate(vector):
            total[i] = (total[i] + value) % modulus
    return total
