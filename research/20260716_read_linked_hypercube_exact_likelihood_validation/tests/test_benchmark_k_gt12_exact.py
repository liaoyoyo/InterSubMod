#!/usr/bin/env python3
"""Unit tests for the resource-bounded k>12 exact-efficiency pilot."""

from __future__ import annotations

import importlib.util
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "benchmark_k_gt12_exact.py"
SPEC = importlib.util.spec_from_file_location("benchmark_k_gt12_exact", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class BenchmarkCaseTests(unittest.TestCase):
    def test_sparse_case_has_eight_active_bits(self):
        case = MODULE.build_case(16, "sparse_active", 7)
        self.assertEqual(case["active_k"], 8)
        self.assertEqual(case["expected_variables"], 256)
        self.assertTrue(all(len(pattern) == 16 for pattern in case["full_patterns"]))
        self.assertTrue(all(len(pattern) == 16 for pattern in case["partial_patterns"]))

    def test_dense_case_has_all_bits_active(self):
        case = MODULE.build_case(13, "dense_active", 11)
        self.assertEqual(case["active_k"], 13)
        self.assertEqual(case["evidence_active_k"], 13)
        self.assertEqual(case["expected_variables"], 8192)

    def test_case_generation_is_deterministic(self):
        self.assertEqual(
            MODULE.build_case(14, "dense_active", 1234),
            MODULE.build_case(14, "dense_active", 1234),
        )

    def test_validate_row_requires_certificate_for_complete_enumeration(self):
        row = {
            "expected_variables": 8,
            "max_sets": 4,
            "solve": {
                "status": "OPTIMAL",
                "n_variables": 8,
                "objective": 2.0,
                "objective_bound": 2.0,
                "mip_gap": 0.0,
            },
            "enumeration": {
                "status": "MAX_SETS_REACHED",
                "complete": True,
                "n_vertex_sets": 4,
                "max_sets": 4,
            },
        }
        self.assertIn("enumeration completeness lacks infeasibility certificate", MODULE.validate_row(row))


if __name__ == "__main__":
    unittest.main()
