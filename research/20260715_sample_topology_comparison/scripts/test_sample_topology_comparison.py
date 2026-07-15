#!/usr/bin/env python3
"""Unit tests for sample topology comparison metrics."""

from __future__ import annotations

import importlib.util
import unittest
from collections import Counter
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("analyze_sample_topology_comparison.py")
SPEC = importlib.util.spec_from_file_location("sample_topology_comparison", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MetricContractTest(unittest.TestCase):
    def test_total_variation_boundaries(self) -> None:
        left = {"a": 1.0, "b": 0.0}
        self.assertEqual(MODULE.total_variation(left, left), 0.0)
        self.assertEqual(MODULE.total_variation(left, {"a": 0.0, "b": 1.0}), 1.0)

    def test_distances_are_symmetric(self) -> None:
        left = {"a": 0.75, "b": 0.25}
        right = {"a": 0.50, "b": 0.50}
        self.assertAlmostEqual(MODULE.total_variation(left, right), MODULE.total_variation(right, left))
        self.assertAlmostEqual(MODULE.jensen_shannon_distance(left, right), MODULE.jensen_shannon_distance(right, left))
        self.assertAlmostEqual(MODULE.hellinger_distance(left, right), MODULE.hellinger_distance(right, left))

    def test_confusion_perfect_agreement(self) -> None:
        confusion = Counter({("a", "a"): 4, ("b", "b"): 6})
        result = MODULE.confusion_statistics(confusion, ("a", "b"))
        self.assertEqual(result["n"], 10)
        self.assertEqual(result["raw_agreement"], 1.0)
        self.assertEqual(result["cohen_kappa"], 1.0)
        self.assertEqual(result["matched_marginal_tvd"], 0.0)

    def test_confusion_kappa_adjusts_prevalence(self) -> None:
        confusion = Counter({("a", "a"): 80, ("a", "b"): 10, ("b", "a"): 10})
        result = MODULE.confusion_statistics(confusion, ("a", "b"))
        self.assertAlmostEqual(result["raw_agreement"], 0.8)
        self.assertLess(result["cohen_kappa"], result["raw_agreement"])

    def test_mean_profiles_first_averages_members(self) -> None:
        result = MODULE.mean_profiles(
            [{"a": 0.8, "b": 0.2}, {"a": 0.4, "b": 0.6}],
            ("a", "b"),
        )
        self.assertAlmostEqual(result["a"], 0.6)
        self.assertAlmostEqual(result["b"], 0.4)


if __name__ == "__main__":
    unittest.main()
