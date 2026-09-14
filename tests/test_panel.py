"""Targeted panels: a locus the assay never amplified is never callable.

The gate that lets a drug reach SUSCEPTIBLE is the callable mask, and on
whole-genome data that is the right gate — every locus was at least attempted.
On a targeted panel it is not sufficient alone, because a mask can be wrong in
the one direction that matters. A whole-genome BED lists loci the panel never
amplified. A generic depth table can report coverage built from off-target
reads or barcode bleed. An engine adapter can supply QC for its full gene list
rather than for the panel actually run.

Any of those lets a drug the assay never examined pass the callable gate and
reach SUSCEPTIBLE — the original defect of this codebase, arriving through a
file instead of through a ``return True``. So the panel intersects the mask and
the narrower of the two always wins.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from myconductor.core.models import LocusCoverage
from myconductor.io.callable_mask import CallableMask
from myconductor.io.panel import (
    Panel, PanelError, restrict, unsupported_drugs,
)

PANEL_DIR = Path(__file__).resolve().parent.parent / "myconductor/catalogue/panels"


def panel(**kwargs) -> Panel:
    base = dict(name="Test panel", version="1", assay="targeted-amplicon",
                source="unit test", loci=frozenset({"rpoB", "katG", "inhA"}))
    base.update(kwargs)
    return Panel(**base)


def mask(loci, fraction=0.99, depth=120) -> CallableMask:
    return CallableMask(
        {locus: LocusCoverage(locus=locus, mean_depth=depth,
                              callable_fraction=fraction, source="depth-table")
         for locus in loci}, source="depth-table")


class RestrictionTests(unittest.TestCase):
    """The case that motivates the whole module."""

    def test_coverage_claimed_off_panel_is_discarded(self):
        # A whole-genome BED says pncA was covered. The panel never amplified
        # it, so pyrazinamide must not be callable on that basis.
        whole_genome = mask({"rpoB", "katG", "inhA", "pncA", "gyrA"})
        self.assertTrue(whole_genome.coverage_for("pncA").callable_fraction)

        restricted, report = restrict(whole_genome, panel())
        self.assertEqual(restricted.coverage_for("pncA").source, "absent")
        self.assertIn("pncA", report.discarded)
        self.assertTrue(report.restricted)

    def test_the_restricted_mask_fails_the_callable_gate_off_panel(self):
        restricted, _ = restrict(mask({"rpoB", "pncA"}), panel())
        callable_, _coverages, reasons = restricted.assess(
            ["pncA"], depth_floor=20, fraction_floor=0.95)
        self.assertFalse(callable_)
        self.assertTrue(any("absent from the callable mask" in r
                            for r in reasons))

    def test_on_panel_coverage_survives_intact(self):
        restricted, report = restrict(mask({"rpoB", "katG", "inhA"}), panel())
        callable_, _coverages, reasons = restricted.assess(
            ["rpoB"], depth_floor=20, fraction_floor=0.95)
        self.assertTrue(callable_, reasons)
        self.assertFalse(report.restricted)
        self.assertEqual(sorted(report.kept), ["inhA", "katG", "rpoB"])

    def test_being_on_the_panel_is_not_evidence_of_coverage(self):
        """A panel says a locus *could* have been sequenced, not that it was."""
        restricted, _ = restrict(mask({"rpoB"}), panel())
        callable_, _coverages, reasons = restricted.assess(
            ["katG"], depth_floor=20, fraction_floor=0.95)
        self.assertFalse(callable_)
        self.assertTrue(reasons)

    def test_the_restricted_source_records_what_narrowed_it(self):
        restricted, _ = restrict(mask({"rpoB", "pncA"}), panel())
        self.assertIn("restricted to Test panel 1", restricted.source)

    def test_the_report_explains_itself(self):
        _restricted, report = restrict(mask({"rpoB", "pncA", "gyrA"}), panel())
        self.assertIn("never amplified", report.describe())
        self.assertIn("2 locus/loci", report.describe())

    def test_an_unrestricted_report_says_so(self):
        _restricted, report = restrict(mask({"rpoB"}), panel())
        self.assertIn("every mask locus is on the panel", report.describe())


class DrugSupportTests(unittest.TestCase):
    def test_a_drug_whose_loci_are_all_on_panel_is_supported(self):
        self.assertTrue(panel().supports_drug(["rpoB"]))
        self.assertTrue(panel().supports_drug(["katG", "inhA"]))

    def test_a_partially_covered_drug_is_not_supported(self):
        """Resistance could be at the locus that was not amplified."""
        self.assertFalse(panel().supports_drug(["katG", "fabG1", "inhA"]))

    def test_a_drug_with_no_declared_loci_is_not_supported(self):
        self.assertFalse(panel().supports_drug([]))

    def test_unsupported_drugs_name_the_missing_loci(self):
        missing = unsupported_drugs(panel(), {
            "rifampicin": ["rpoB"],
            "pyrazinamide": ["pncA"],
            "isoniazid": ["katG", "fabG1", "inhA"],
        })
        self.assertNotIn("rifampicin", missing)
        self.assertEqual(missing["pyrazinamide"], ["pncA"])
        self.assertEqual(missing["isoniazid"], ["fabG1"])

    def test_a_drug_with_no_loci_at_all_is_listed_as_unsupported(self):
        self.assertIn("mystery", unsupported_drugs(panel(), {"mystery": []}))


class DeclarationTests(unittest.TestCase):
    def test_an_empty_panel_is_refused(self):
        """Safe but useless, and far more likely to be a mistake."""
        with self.assertRaises(PanelError) as caught:
            panel(loci=frozenset())
        self.assertIn("declares no loci", str(caught.exception))

    def test_an_unknown_assay_kind_is_refused(self):
        with self.assertRaises(PanelError):
            panel(assay="nanopore-ish")

    def test_provenance_is_mandatory(self):
        for attribute in ("name", "version", "assay", "source"):
            with self.assertRaises(PanelError):
                panel(**{attribute: "   "})

    def test_a_misspelled_field_is_refused_not_ignored(self):
        """A typo'd 'loci' key would otherwise declare an empty panel."""
        with self.assertRaises(PanelError) as caught:
            Panel.from_dict({"name": "p", "version": "1",
                             "assay": "targeted-amplicon", "source": "s",
                             "locii": ["rpoB"]})
        self.assertIn("unknown panel field", str(caught.exception))

    def test_loci_must_be_a_list(self):
        with self.assertRaises(PanelError):
            Panel.from_dict({"name": "p", "version": "1",
                             "assay": "targeted-amplicon", "source": "s",
                             "loci": "rpoB"})

    def test_round_trips_through_a_file(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "panel.json"
            path.write_text(json.dumps({
                "name": "p", "version": "1", "assay": "targeted-amplicon",
                "source": "s", "loci": ["rpoB", "katG"]}), encoding="utf-8")
            loaded = Panel.load(path)
        self.assertEqual(loaded.loci, frozenset({"rpoB", "katG"}))

    def test_an_unreadable_panel_is_refused_with_its_path(self):
        with self.assertRaises(PanelError):
            Panel.load("no/such/panel.json")


class VerificationTests(unittest.TestCase):
    """Over-declaring is the dangerous direction, so say which list this is."""

    def test_panels_are_unverified_unless_stated(self):
        self.assertFalse(panel().verified)
        self.assertIn("UNVERIFIED locus list", panel().describe())

    def test_a_verified_panel_says_so(self):
        self.assertIn("verified", panel(verified=True).describe())

    def test_the_bundled_panels_ship_unverified(self):
        """They are reconstructed from published descriptions, not designs."""
        for name in ("deeplex-myctb", "first-line-core"):
            loaded = Panel.load(PANEL_DIR / f"{name}.json")
            self.assertFalse(loaded.verified, name)
            self.assertIn("ILLUSTRATIVE", loaded.note)

    def test_the_bundled_panels_name_only_loci_the_profile_knows(self):
        from myconductor.catalogue.profile import load_profile
        profile = load_profile()
        for name in ("deeplex-myctb", "first-line-core"):
            loaded = Panel.load(PANEL_DIR / f"{name}.json")
            unknown = sorted(loaded.loci - set(profile.loci))
            self.assertEqual(unknown, [], f"{name} names unknown loci")


class NarrowPanelTests(unittest.TestCase):
    """End to end against the bundled minimal panel and the real profile."""

    def test_a_first_line_panel_cannot_support_second_line_drugs(self):
        from myconductor.catalogue.profile import load_profile
        profile = load_profile()
        loaded = Panel.load(PANEL_DIR / "first-line-core.json")
        required = {drug: profile.required_loci(drug) for drug in profile.drugs}
        missing = unsupported_drugs(loaded, required)

        # Covered by this narrow panel.
        self.assertNotIn("rifampicin", missing)
        self.assertNotIn("isoniazid", missing)
        self.assertNotIn("ethambutol", missing)
        # Not covered, and therefore not susceptible-able.
        for drug in ("bedaquiline", "linezolid", "moxifloxacin",
                     "amikacin", "pyrazinamide", "clofazimine"):
            self.assertIn(drug, missing, drug)

    def test_a_whole_genome_mask_cannot_rescue_an_off_panel_drug(self):
        """The headline case, with the real profile and a real panel."""
        from myconductor.catalogue.profile import load_profile
        profile = load_profile()
        loaded = Panel.load(PANEL_DIR / "first-line-core.json")
        # A mask generated against the whole genome claims everything.
        everything = mask(set(profile.loci))
        restricted, report = restrict(everything, loaded)

        callable_, _coverages, reasons = restricted.assess(
            profile.required_loci("bedaquiline"),
            depth_floor=20, fraction_floor=0.95)
        self.assertFalse(callable_)
        self.assertTrue(reasons)
        self.assertIn("atpE", report.discarded)


class EndToEndTests(unittest.TestCase):
    """Through the real pipeline, which is where it has to hold."""

    def _vcf(self, root: Path) -> Path:
        path = root / "sample.vcf"
        path.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "NC_000962.3\t761155\t.\tC\tT\t60\tPASS\tGENE=rpoB\tGT:DP\t1:100\n",
            encoding="utf-8")
        return path

    def _depth_table(self, root: Path, loci) -> Path:
        path = root / "depth.tsv"
        rows = ["locus\tmean_depth\tcallable_fraction"]
        rows += [f"{locus}\t120\t0.99" for locus in sorted(loci)]
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        return path

    def test_off_panel_drugs_become_not_assessed_through_the_pipeline(self):
        """The whole point, exercised through Myconductor.analyze.

        The depth table claims every locus in the profile was callable at 120x
        — the shape a whole-genome pipeline produces against a targeted run.
        Unrestricted, every drug is merely INDETERMINATE, which reads as "we
        looked and could not conclude". With a first-line-only panel the assay
        never amplified atpE, Rv0678, pncA, gyrA, rrs, rrl or rplC, so those
        drugs must say NOT_ASSESSED: we never looked.

        Note the distinction being protected. Neither state admits a drug to a
        regimen, so nothing unsafe happens either way here — but a clinician
        reading INDETERMINATE may reasonably order a repeat or a deeper run,
        where NOT_ASSESSED tells them no amount of sequencing on *this assay*
        will ever answer. The panel is what makes that difference reportable.
        """
        from myconductor.catalogue.profile import load_profile
        from myconductor.core.models import Call
        from myconductor.core.pipeline import Myconductor

        profile = load_profile()
        conductor = Myconductor()
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            vcf = self._vcf(root)
            depth = self._depth_table(root, profile.loci)
            unrestricted = {r.drug: r.call for r
                            in conductor.analyze(vcf, mask_path=depth).drug_results}
            restricted = {r.drug: r.call for r in conductor.analyze(
                vcf, mask_path=depth,
                panel_path=PANEL_DIR / "first-line-core.json").drug_results}

        off_panel = ("bedaquiline", "clofazimine", "linezolid", "moxifloxacin",
                     "amikacin", "pyrazinamide", "pretomanid")
        for drug in off_panel:
            self.assertIsNot(restricted[drug], Call.SUSCEPTIBLE, drug)
            self.assertIs(restricted[drug], Call.NOT_ASSESSED, drug)
            self.assertIsNot(unrestricted[drug], Call.NOT_ASSESSED,
                             f"{drug}: precondition, the mask alone did not "
                             f"already withhold assessment")

    def test_a_panel_never_makes_a_call_less_cautious(self):
        """Restriction may only move calls toward withholding, never away.

        A panel that *granted* something the unrestricted mask withheld would
        mean the narrowing had added information, which is impossible: it only
        ever removes coverage claims.
        """
        from myconductor.catalogue.profile import load_profile
        from myconductor.core.models import Call
        from myconductor.core.pipeline import Myconductor

        profile = load_profile()
        conductor = Myconductor()
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            vcf = self._vcf(root)
            depth = self._depth_table(root, profile.loci)
            unrestricted = {r.drug: r.call for r
                            in conductor.analyze(vcf, mask_path=depth).drug_results}
            restricted = {r.drug: r.call for r in conductor.analyze(
                vcf, mask_path=depth,
                panel_path=PANEL_DIR / "first-line-core.json").drug_results}

        for drug, call in restricted.items():
            if call is Call.SUSCEPTIBLE:
                self.assertIs(unrestricted[drug], Call.SUSCEPTIBLE, drug)

    def test_the_report_says_which_drugs_the_assay_cannot_cover(self):
        """"This assay does not cover pyrazinamide" beats a bare NOT_ASSESSED."""
        from myconductor.catalogue.profile import load_profile
        from myconductor.core.pipeline import Myconductor

        profile = load_profile()
        conductor = Myconductor()
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            report = conductor.analyze(
                self._vcf(root),
                mask_path=self._depth_table(root, profile.loci),
                panel_path=PANEL_DIR / "first-line-core.json")

        self.assertEqual(report.panel["name"], "First-line core")
        self.assertFalse(report.panel["verified"])
        self.assertIn("pyrazinamide", report.panel["unsupported_drugs"])
        self.assertIn("pncA", report.panel["unsupported_drugs"]["pyrazinamide"])
        self.assertIn("atpE", report.panel["discarded"])

    def test_an_on_panel_drug_still_reaches_a_verdict(self):
        """Restriction must not break the drugs the assay does cover."""
        from myconductor.catalogue.profile import load_profile
        from myconductor.core.models import Call
        from myconductor.core.pipeline import Myconductor

        profile = load_profile()
        conductor = Myconductor()
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            report = conductor.analyze(
                self._vcf(root),
                mask_path=self._depth_table(root, profile.loci),
                panel_path=PANEL_DIR / "first-line-core.json")
        results = {r.drug: r for r in report.drug_results}
        self.assertIn(results["ethambutol"].call,
                      (Call.SUSCEPTIBLE, Call.INDETERMINATE))

    def test_no_panel_leaves_the_report_unchanged(self):
        from myconductor.catalogue.profile import load_profile
        from myconductor.core.pipeline import Myconductor

        profile = load_profile()
        conductor = Myconductor()
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            report = conductor.analyze(
                self._vcf(root),
                mask_path=self._depth_table(root, profile.loci))
        self.assertEqual(report.panel, {})


if __name__ == "__main__":
    unittest.main()
