"""mycobench — a real-data validation harness for Myconductor.

Runs validated AMR engines over public mycobacterial isolates and reports what
the comparison actually supports. Two tracks, kept separate because they answer
different questions:

**Concordance** — Myconductor against TB-Profiler on the same isolates. Needs
no phenotype. Measures agreement, not correctness: two tools reading the same
catalogue agreeing tells you they read the same catalogue.

**Accuracy** — a prediction against a laboratory phenotype. Needs paired DST or
MIC data, which public Nigerian BioSamples do not carry, so this track is built
from the CRyPTIC compendium instead.

A third panel of non-tuberculous mycobacteria exists to be **refused**: a drug
call on a non-MTBC genome is a pipeline failure, not a result.

Nothing here fabricates data. Accessions come from NCBI, the catalogue comes
from WHO's published files, and phenotypes come from CRyPTIC — each pinned to a
version and recorded in a lock.
"""
from __future__ import annotations

#: Must equal myconductor.__version__ and the version in pyproject.toml. Both
#: packages ship from one distribution, so a separate number here is not a
#: second version, it is a wrong one: it is stamped into the analysis manifest,
#: into run reports, and into the User-Agent sent to NCBI and EBI. A run whose
#: recorded version does not identify the code that produced it cannot be
#: reproduced from its own outputs. tests/test_packaging.py holds the three in
#: step.
__version__ = "0.2.0"

#: Contact details NCBI asks callers to send. Set NCBI_EMAIL/NCBI_API_KEY in
#: the environment; an API key raises the request allowance from 3/s to 10/s.
USER_AGENT = f"mycobench/{__version__} (+https://github.com/ubeffiong/myconductor)"

__all__ = ["__version__", "USER_AGENT"]
