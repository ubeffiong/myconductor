"""Packaging and provenance invariants.

Two things here that no other test covers, both of which failed silently.

**One distribution, one version.** ``myconductor`` and ``mycobench`` ship from a
single ``pyproject.toml``. When their ``__version__`` strings drift apart, the
wrong one is stamped into the analysis manifest, into run reports and into the
User-Agent sent to NCBI and EBI — and a result whose recorded version does not
identify the code that produced it cannot be reproduced from its own output.
``mycobench`` sat at 0.1.0 while the distribution was 0.2.0.

**A manifest must be enough to run the analysis again.** The codebase pins
bootstrap and permutation seeds, a catalogue commit, and fixture hashes so a
figure can be regenerated. The manifest recorded none of its own inputs, which
broke that chain at the last step: without the sampling seed nothing records
*which* isolates were drawn, so a 400-isolate result could not be reproduced
even in principle.
"""
from __future__ import annotations

import re
import unittest
from pathlib import Path

import mycobench
import myconductor
from mycobench.analysis import stats
from mycobench.analysis.pipeline import AnalysisInputs, AnalysisResult

ROOT = Path(__file__).resolve().parent.parent


def pyproject_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert match, "pyproject.toml has no version"
    return match.group(1)


class VersionTests(unittest.TestCase):
    def test_both_packages_report_the_distribution_version(self):
        expected = pyproject_version()
        self.assertEqual(myconductor.__version__, expected)
        self.assertEqual(mycobench.__version__, expected,
                         "mycobench.__version__ drifted from pyproject.toml; "
                         "it is stamped into manifests and the NCBI User-Agent")

    def test_user_agent_carries_the_real_version(self):
        self.assertIn(mycobench.__version__, mycobench.USER_AGENT)


class ProvenanceTests(unittest.TestCase):
    """The manifest must record what it would take to run this again."""

    REQUIRED = (
        "mycobench_version", "generated_at", "catalogue", "phenotypes",
        "reuse_table", "limit", "sample_seed", "min_carriers",
        "bootstrap_seed", "permutation_seed", "accepted_phenotype_quality",
    )

    def test_result_carries_a_provenance_slot(self):
        self.assertEqual(AnalysisResult().provenance, {})

    def test_manifest_records_every_reproducibility_input(self):
        from mycobench.analysis import pipeline
        source = Path(pipeline.__file__).read_text(encoding="utf-8")
        for key in self.REQUIRED:
            self.assertIn(f'"{key}"', source,
                          f"{key} is not recorded in the manifest provenance")

    def test_seeds_are_fixed_constants_not_clock_reads(self):
        # A seed that moves between runs makes every interval unreproducible.
        self.assertIsInstance(stats.BOOTSTRAP_SEED, int)
        self.assertIsInstance(stats.PERMUTATION_SEED, int)
        self.assertNotEqual(stats.BOOTSTRAP_SEED, stats.PERMUTATION_SEED)

    def test_inputs_expose_everything_provenance_claims_to_record(self):
        """Guards against a field being recorded that no longer exists."""
        fields = AnalysisInputs.__dataclass_fields__
        for name in ("catalogue", "phenotypes", "reuse_table", "limit",
                     "sample_seed", "min_carriers", "drugs", "partition",
                     "independent_clusters"):
            self.assertIn(name, fields)


if __name__ == "__main__":
    unittest.main()
