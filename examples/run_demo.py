"""Programmatic use of Myconductor, including the parts that refuse to answer.

Run from the repo root:  python examples/run_demo.py

Demonstrates, in order:
  1. the same input analysed WITHOUT coverage evidence -- nothing susceptible;
  2. the same input WITH a callable mask -- specific drugs become susceptible;
  3. the federated contract at isolate level, including a promotion that is
     correctly refused because the association only holds in one lineage.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from myconductor import Myconductor
from myconductor.core.models import Call
from myconductor.federated.catalogue_update import (
    IsolateObservation,
    aggregate,
    evaluate_all,
)
from myconductor.federated.transport import (
    Site,
    SiteRegistry,
    SubmissionVerifier,
    build_submission,
    k_anonymity_violations,
)
from myconductor.reporting.fhir import to_fhir_json
from myconductor.reporting.render import render_text

DATA = Path(__file__).resolve().parent.parent / "myconductor" / "data"
VCF = DATA / "example_input.vcf"
MASK = DATA / "example_callable.tsv"


def without_coverage() -> None:
    print("=" * 74)
    print(" 1. No coverage evidence supplied")
    print("=" * 74)
    report = Myconductor(platform="illumina").analyze(VCF)
    susceptible = [r.drug for r in report.drug_results if r.permits_use]
    not_assessed = [r.drug for r in report.drug_results
                    if r.call is Call.NOT_ASSESSED]
    print(f"  susceptible   : {susceptible or 'none — as it should be'}")
    print(f"  not assessed  : {', '.join(not_assessed)}")
    print("  Nothing is susceptible because nothing showed the loci were")
    print("  sequenced. The old code called all of these 'effective'.\n")


def with_coverage() -> None:
    print("=" * 74)
    print(" 2. With a callable mask")
    print("=" * 74)
    report = Myconductor(platform="illumina").analyze(VCF, mask_path=MASK)
    print(render_text(report))

    out = Path("demo_report.fhir.json")
    out.write_text(to_fhir_json(report))
    print(f"\nWrote {out} "
          f"(note the dataAbsentReason entries for drugs with no verdict)\n")


def federated() -> None:
    print("=" * 74)
    print(" 3. Federated evidence at isolate level")
    print("=" * 74)

    # A causal variant, seen across two lineages and three sites.
    observations = []
    for i in range(80):
        observations.append(IsolateObservation(
            isolate_id=f"ISO-R{i}", site_id=f"SITE-{i % 3}",
            drug="rifampicin", phenotype=True,
            variant_keys=("rpoB_S450L",),
            lineage="lineage2" if i % 2 else "lineage4",
            dst_method="MGIT", lab_quality="accredited",
        ))
    for i in range(20):
        observations.append(IsolateObservation(
            isolate_id=f"ISO-S{i}", site_id=f"SITE-{i % 3}",
            drug="rifampicin", phenotype=False,
            variant_keys=("rpoB_V170F",),
            lineage="lineage2" if i % 2 else "lineage4",
            dst_method="MGIT", lab_quality="accredited",
        ))
    # A lineage marker: only ever seen in one lineage.
    for i in range(60):
        observations.append(IsolateObservation(
            isolate_id=f"ISO-L{i}", site_id=f"SITE-{i % 3}",
            drug="rifampicin", phenotype=True,
            variant_keys=("Rv1234_A1B",), lineage="lineage2",
            dst_method="MGIT", lab_quality="accredited",
        ))

    # Signed transport between a site and the coordinator.
    registry = SiteRegistry()
    site = registry.register(Site(
        site_id="SITE-0", name="Demo Laboratory 0",
        secret=SiteRegistry.new_secret(), accredited=True,
    ))
    site_0 = [o for o in observations if o.site_id == "SITE-0"]
    submission = build_submission(site, site_0)
    accepted = SubmissionVerifier(registry).verify(submission)
    print(f"  signed submission from {site.site_id}: "
          f"{len(accepted)} observation(s) accepted")

    violations = k_anonymity_violations(observations)
    print(f"  disclosure-floor violations: {len(violations)}")

    aggregated = aggregate([observations])
    for candidate in evaluate_all(aggregated, known_resistance_variants={"rpoB_S450L"}):
        verdict = "ELIGIBLE" if candidate.eligible else "BLOCKED"
        print(f"\n  [{verdict}] {candidate.variant_key} / {candidate.drug}")
        print(f"      {candidate.n_resistant}R / {candidate.n_susceptible}S "
              f"across {candidate.n_sites} site(s), "
              f"lineages supporting: {candidate.supporting_lineages or 'none'}")
        for reason in candidate.blocked_reasons:
            print(f"      - {reason}")
    print("\n  Note Rv1234_A1B: 60 resistant isolates, 100% association, three")
    print("  sites — and refused, on the lineage check alone. It is only ever")
    print("  seen in lineage2, so it cannot be separated from a marker of a")
    print("  resistant-enriched lineage. The old variant-level counter would")
    print("  have promoted it into the catalogue and shipped it to every site.")
    print("  rpoB_S450L, the causal variant, clears the bar — and still needs")
    print("  expert curation before release.\n")


def main() -> None:
    without_coverage()
    with_coverage()
    federated()


if __name__ == "__main__":
    main()
