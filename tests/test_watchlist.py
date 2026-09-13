import tempfile
import unittest
from pathlib import Path
from myconductor.federated.watchlist import WatchlistStore, DiscordanceTicket, aggregate
from myconductor.federated.transport import DisclosureThresholds
from dataclasses import asdict


class WatchlistTests(unittest.TestCase):
    def ticket(self, key="1", site="site1", isolate="i"):
        return DiscordanceTicket(key, "s", isolate, site, "test", "test-drug", "resistant",
                                  "susceptible", "report-hash", "2026-09-01")

    def test_resolution_reopen_and_persistence(self):
        store = WatchlistStore()
        store.append("opened", asdict(self.ticket()))
        store.transition("1", "under_investigation", "reviewer", "started")
        store.transition("1", "resolved", "reviewer", "reviewed", "assay discrepancy documented")
        store.transition("1", "under_investigation", "reviewer", "new evidence")
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "watchlist.json"
            store.save(path)
            self.assertEqual(store.state(), WatchlistStore.load(path).state())

    def test_small_cells_suppressed_and_duplicates_do_not_inflate(self):
        ticket = asdict(self.ticket())
        result = aggregate({"site1": [ticket] * 20}, DisclosureThresholds(2, 2))
        self.assertEqual(result["released"], [])
        self.assertTrue(result["has_suppressed_cells"])

    def test_released_output_contains_only_counts(self):
        result = aggregate({"site1":[asdict(self.ticket())],
                            "site2":[asdict(self.ticket("2", "site2", "i"))]}, DisclosureThresholds(2, 2))
        self.assertEqual(result["released"][0]["unresolved_isolates"], 2)
        self.assertNotIn("isolate_id", str(result))
