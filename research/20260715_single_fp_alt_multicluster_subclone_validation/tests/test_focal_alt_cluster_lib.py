#!/usr/bin/env python3
"""Deterministic tests for focal ALT null-validated clustering."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "focal_alt_cluster_lib.py"
SPEC = importlib.util.spec_from_file_location("focal_alt_cluster_lib", SCRIPT)
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


class FocalAltClusterTest(unittest.TestCase):
    def test_clear_two_group_pattern_is_detected(self) -> None:
        first = np.tile(np.asarray([0.02, 0.03, 0.01, 0.04, 0.02, 0.03]), (4, 1))
        second = np.tile(np.asarray([0.98, 0.97, 0.99, 0.96, 0.98, 0.97]), (4, 1))
        methylation = np.vstack([first, second])
        distance = M.bernoulli_distance(methylation)
        result = M.analyze_phylo(distance, methylation, seeds=3, rnull=8)
        self.assertGreaterEqual(result["coarse_ng"], 2)
        self.assertGreaterEqual(result["modal_fraction"], M.MODAL_CONFIDENCE)
        self.assertIn("modal_assignment_ari_median", result)
        self.assertGreaterEqual(result["modal_assignment_ari_median"], 0.0)

    def test_single_pattern_is_not_forced_by_null_method(self) -> None:
        rng = np.random.default_rng(7)
        methylation = np.clip(rng.normal(0.1, 0.01, size=(8, 12)), 0.01, 0.99)
        distance = M.bernoulli_distance(methylation)
        result = M.analyze_phylo(distance, methylation, seeds=3, rnull=8)
        self.assertEqual(result["coarse_ng"], 1)

    def test_peel_removes_incomplete_pair(self) -> None:
        distance = np.asarray(
            [
                [0.0, 0.1, 0.2],
                [0.1, 0.0, -1.0],
                [0.2, -1.0, 0.0],
            ]
        )
        kept = M.peel_complete(distance)
        self.assertLess(len(kept), 3)

    def test_hp_family_keeps_extended_tags(self) -> None:
        self.assertEqual(M.hp_family("1-2"), "HP1-side")
        self.assertEqual(M.hp_family("2-2"), "HP2-side")
        self.assertEqual(M.hp_family("3"), "HP3-ambiguous")
        self.assertEqual(M.hp_family("."), "untagged")

    def test_row_circular_null_preserves_each_read_values_and_missingness(self) -> None:
        methylation = np.asarray(
            [
                [0.1, np.nan, 0.8, 0.9, 0.2],
                [np.nan, 0.3, 0.4, 0.7, 0.6],
            ]
        )
        permuted = M.permute_methylation(
            methylation, np.random.default_rng(19), mode="row_circular"
        )
        for original, shifted in zip(methylation, permuted):
            original_values = sorted(original[np.isfinite(original)].tolist())
            shifted_values = sorted(shifted[np.isfinite(shifted)].tolist())
            self.assertEqual(original_values, shifted_values)
            self.assertEqual(np.isnan(original).sum(), np.isnan(shifted).sum())

    def test_row_circular_null_samples_identity_shift(self) -> None:
        class IdentityRng:
            @staticmethod
            def integers(low: int, high: int) -> int:
                self.assertEqual((low, high), (0, 4))
                return 0

        methylation = np.asarray([[0.1, np.nan, 0.8, 0.9]])
        permuted = M.permute_methylation(methylation, IdentityRng(), mode="row_circular")
        np.testing.assert_equal(permuted, methylation)

    def test_unknown_null_mode_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            M.permute_methylation(
                np.asarray([[0.1, 0.9], [0.2, 0.8]]),
                np.random.default_rng(3),
                mode="unknown",
            )

    def test_empty_null_distribution_fails_closed(self) -> None:
        methylation = np.full((6, 4), np.nan)
        distance = np.full((6, 6), 0.1, dtype=float)
        np.fill_diagonal(distance, 0.0)
        distance[:3, 3:] = 0.9
        distance[3:, :3] = 0.9
        result = M.analyze_phylo(distance, methylation, seeds=2, rnull=4)
        self.assertEqual(result["coarse_ng"], 1)
        self.assertTrue(
            any(
                trace.get("failure") == "insufficient_valid_null"
                for trace in result["coarse_split_trace"]
            )
        )


if __name__ == "__main__":
    unittest.main()
