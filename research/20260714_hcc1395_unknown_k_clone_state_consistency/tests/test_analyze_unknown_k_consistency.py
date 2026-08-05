#!/usr/bin/env python3
"""Unit tests for retained-pattern state bounds."""

from __future__ import annotations

import importlib.util
import itertools
import random
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_unknown_k_consistency.py"
SPEC = importlib.util.spec_from_file_location("unknown_k", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class StateBoundTests(unittest.TestCase):
    def test_conflict(self) -> None:
        self.assertTrue(MODULE.pattern_conflict("AXX", "RXX"))
        self.assertFalse(MODULE.pattern_conflict("AXX", "XAX"))
        self.assertFalse(MODULE.pattern_conflict("XRX", "XXA"))

    def test_exact_minimum_cover(self) -> None:
        self.assertEqual(MODULE.minimum_compatible_state_count([]), 0)
        self.assertEqual(MODULE.minimum_compatible_state_count(["AXX", "XAX"]), 1)
        self.assertEqual(MODULE.minimum_compatible_state_count(["AXX", "RXX"]), 2)
        self.assertEqual(MODULE.minimum_compatible_state_count(["AA", "AR", "RA"]), 3)
        self.assertEqual(
            MODULE.minimum_compatible_state_count(["AXX", "RAX", "RRA", "AAA"]), 3
        )

    def test_maximum_support_respects_read_count(self) -> None:
        self.assertEqual(MODULE.maximum_supported_state_count({"AXX": 2}, 3), 2)
        self.assertEqual(MODULE.maximum_supported_state_count({"AXX": 10}, 3), 4)
        self.assertEqual(
            MODULE.maximum_supported_state_count({"XXX": 20}, 3, mutation_only=True), 7
        )

    def test_full_catalog_is_identified_only_without_x(self) -> None:
        full = {
            "family": "1",
            "obs_subreads": {"RR": 5, "AR": 3},
            "n_reads": 8,
            "trees": [],
        }
        partial = {
            "family": "1",
            "obs_subreads": {"RX": 5, "AR": 3},
            "n_reads": 8,
            "trees": [],
        }
        full_result = MODULE.analyze_lineage(full, 2)
        partial_result = MODULE.analyze_lineage(partial, 2)
        self.assertTrue(full_result["catalog_identified_under_retained_table"])
        self.assertFalse(partial_result["catalog_identified_under_retained_table"])
        self.assertEqual(full_result["k_state_min_retained"], 2)

    def test_exact_minimum_matches_bruteforce_random_small_cases(self) -> None:
        rng = random.Random(20260714)
        for k in range(1, 5):
            all_patterns = ["".join(chars) for chars in itertools.product("RAX", repeat=k)]
            states = MODULE.full_states(k)
            for _ in range(40):
                patterns = rng.sample(all_patterns, rng.randint(1, min(8, len(all_patterns))))
                expected = None
                for size in range(1, len(states) + 1):
                    if any(
                        all(any(MODULE.compatible(pattern, state) for state in catalog) for pattern in patterns)
                        for catalog in itertools.combinations(states, size)
                    ):
                        expected = size
                        break
                self.assertEqual(MODULE.minimum_compatible_state_count(patterns), expected)


if __name__ == "__main__":
    unittest.main()
