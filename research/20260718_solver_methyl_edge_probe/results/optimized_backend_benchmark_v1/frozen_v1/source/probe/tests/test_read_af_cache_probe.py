#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import sys
import unittest


HERE = pathlib.Path(__file__).resolve().parent
PROBE_ROOT = HERE.parent
CURRENT_ROOT = PROBE_ROOT.parent / "20260716_read_linked_hypercube_exact_likelihood_validation"
sys.path.insert(0, str(PROBE_ROOT / "scripts"))
sys.path.insert(0, str(CURRENT_ROOT / "scripts"))

from read_af_cache_probe import instrument_case, load_module  # noqa: E402


RANKER = CURRENT_ROOT / "scripts" / "build_m2_patterns_and_rank.py"


class ReadAfCacheProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module(str(RANKER), "ranker_cache_probe_test")

    def test_same_vertex_set_two_edges_fits_once(self):
        row = instrument_case(
            self.module,
            "same_v",
            [("AR", 3), ("RA", 3), ("AA", 3)],
            2,
        )
        self.assertEqual(2, row["raw_tree_candidates_T"])
        self.assertEqual(1, row["distinct_vertex_sets_V"])
        self.assertEqual(2, row["top_edge_variants"])
        self.assertEqual(1, row["primary_fit_calls"])
        self.assertTrue(row["check_one_primary_fit_per_vertex_set"])

    def test_two_vertex_sets_fit_twice(self):
        row = instrument_case(self.module, "two_v", [("AA", 3)], 2)
        self.assertEqual(2, row["distinct_vertex_sets_V"])
        self.assertEqual(2, row["primary_fit_calls"])
        self.assertTrue(row["check_one_primary_fit_per_vertex_set"])


if __name__ == "__main__":
    unittest.main()
