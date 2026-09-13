"""Minimal programmatic use of Myconductor.

Run from the repo root:  python examples/run_demo.py
"""
import sys
from pathlib import Path

# Allow running without `pip install` from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from myconductor import Myconductor
from myconductor.federated.catalogue_update import extract_shareable
from myconductor.reporting.fhir import to_fhir_json
from myconductor.reporting.render import render_text

DATA = Path(__file__).resolve().parent.parent / "myconductor" / "data"


def main() -> None:
    conductor = Myconductor()
    report = conductor.analyze(DATA / "example_input.vcf")

    print(render_text(report))

    # Emit a FHIR bundle for LIMS/EHR integration.
    Path("demo_report.fhir.json").write_text(to_fhir_json(report))
    print("\nWrote demo_report.fhir.json")

    # Show the privacy-preserving federated contract: only aggregated,
    # de-identified variant->phenotype tallies would ever leave the site.
    confirmed = {"pncA_D12A": True}  # e.g. later phenotypic DST confirmed R
    shareable = extract_shareable(report, confirmed)
    print(f"\nShareable observations (no genomes leave the node): {shareable}")


if __name__ == "__main__":
    main()
