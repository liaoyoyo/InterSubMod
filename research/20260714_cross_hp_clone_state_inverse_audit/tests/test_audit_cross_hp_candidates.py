#!/usr/bin/env python3
"""Contract tests for the cross-HP candidate audit."""

from __future__ import annotations

import importlib.util
import math
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_cross_hp_candidates.py"
SPEC = importlib.util.spec_from_file_location("audit_cross_hp_candidates", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class CrossHpAuditTest(unittest.TestCase):
    def test_updated_toy_is_direct_plus_sister_and_recovers_mixture(self) -> None:
        direct = {
            "n_trees": 1,
            "trees": [{"edges": [["ROOT", "AR"], ["AR", "AA"]]}],
            # pi=(0.10,0.20,0.25,0.15,0.30) scaled to 100 full reads.
            "obs_populations": {"RR": 25, "AR": 50, "AA": 25},
        }
        sister = {
            "n_trees": 1,
            "trees": [{"edges": [["ROOT", "AR"], ["ROOT", "RA"]]}],
            "obs_populations": {"RR": 55, "AR": 15, "RA": 30},
        }
        self.assertEqual(MODULE.classify_tree(direct["trees"][0], 2)["shape_class"], "direct_only")
        self.assertEqual(MODULE.classify_tree(sister["trees"][0], 2)["shape_class"], "sister_only")
        result = MODULE.infer_two_site_catalog(direct, sister, minread=3)
        self.assertTrue(result["evaluable"])
        self.assertTrue(result["plugin_nonnegative"])
        expected = {
            "background_pi0": 0.10,
            "clone1_pi1": 0.20,
            "clone2_pi2": 0.25,
            "clone3_pi3": 0.15,
            "clone4_pi4": 0.30,
        }
        for key, value in expected.items():
            self.assertTrue(math.isclose(result["plugin_estimates"][key], value, abs_tol=1e-12))

    def test_unexpected_full_state_fails_fixed_catalog(self) -> None:
        direct = {
            "n_trees": 1,
            "trees": [{"edges": [["ROOT", "AR"], ["AR", "AA"]]}],
            "obs_populations": {"RR": 25, "AR": 45, "RA": 5, "AA": 25},
        }
        sister = {
            "n_trees": 1,
            "trees": [{"edges": [["ROOT", "AR"], ["ROOT", "RA"]]}],
            "obs_populations": {"RR": 55, "AR": 15, "RA": 30},
        }
        result = MODULE.infer_two_site_catalog(direct, sister, minread=3)
        self.assertFalse(result["evaluable"])
        self.assertEqual(result["unexpected_direct"], ["RA"])

    def test_tied_direct_orientation_is_not_tree_unique_inference(self) -> None:
        direct = {
            "n_trees": 2,
            "trees": [
                {"edges": [["ROOT", "AR"], ["AR", "AA"]]},
                {"edges": [["ROOT", "RA"], ["RA", "AA"]]},
            ],
            "obs_populations": {"RR": 20, "AR": 40, "RA": 10, "AA": 30},
        }
        sister = {
            "n_trees": 1,
            "trees": [{"edges": [["ROOT", "AR"], ["ROOT", "RA"]]}],
            "obs_populations": {"RR": 50, "AR": 20, "RA": 30},
        }
        result = MODULE.infer_two_site_catalog(direct, sister, minread=3)
        self.assertFalse(result["evaluable"])
        self.assertEqual(result["reason"], "tree_not_unique")


if __name__ == "__main__":
    unittest.main()
