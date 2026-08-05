#!/usr/bin/env python3
"""Focused unit tests for the compact HCC C++ parity gate."""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import unittest
from collections import Counter

MODULE_PATH = Path(__file__).with_name("validate_hcc1395_cpp_topology_af_gate.py")
SPEC = importlib.util.spec_from_file_location(
    "validate_hcc1395_cpp_topology_af_gate", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


class HccCppTopologyAfGateTest(unittest.TestCase):
    def test_old_family_reconstruction_matches_pinned_cpp_digest(self) -> None:
        group = {
            "region_id": "chr1|PS=100|HP=1|U_SYN:B0001",
            "hp_family": "1",
            "positions": [101, 102],
            "populations_by_hp": {},
            "subread_groups_by_hp": {"1": {"AX": 5, "XA": 5}},
        }
        unit = {
            "region": group["region_id"],
            "n_hidden": 2,
            "n_trees": 3,
            "trees": [
                {"edges": [["ROOT", "H_AR"], ["ROOT", "H_RA"]]},
                {"edges": [["ROOT", "H_AR"], ["H_AR", "H_AA"]]},
                {"edges": [["ROOT", "H_RA"], ["H_RA", "H_AA"]]},
            ],
        }

        observed = gate.old_complete_contract(unit, group)

        self.assertEqual(observed["active_bit_count"], 2)
        self.assertEqual(observed["objective_h"], 2)
        self.assertEqual(observed["minimum_vertex_set_count"], 3)
        self.assertEqual(observed["total_tree_count"], 3)
        self.assertEqual(
            observed["minimum_family_sha256"],
            "b0d4bb34ab635f53c33b1acc858cb12f7e22dd94c53f1e65d047fe4a0e1ae6a5",
        )

    def test_capped_policy_accepts_complete_or_clean_abstain_only(self) -> None:
        complete = {
            "objective_status": "OBJECTIVE_CERTIFIED",
            "family_status": "FAMILY_COMPLETE",
            "objective_h": 5,
            "minimum_vertex_set_count": 4,
            "minimum_family_sha256": "a" * 64,
            "total_tree_count": "9",
        }
        abstain = {
            "objective_status": "OBJECTIVE_CERTIFIED",
            "family_status": "FAMILY_INCOMPLETE_CAP",
            "unit_status": "family_incomplete",
            "read_af_status": "not_evaluable_family_incomplete",
            "minimum_vertex_set_count": None,
            "minimum_family_sha256": None,
            "total_tree_count": None,
            "best_score_fraction": None,
            "best_tree_tie_count": None,
            "best_tree_unique": None,
        }

        self.assertEqual(gate.validate_capped_row(complete), ("cpp_complete", []))
        self.assertEqual(gate.validate_capped_row(abstain), ("cpp_abstain", []))

        leaking = dict(abstain)
        leaking["best_score_fraction"] = "0/1"
        outcome, issues = gate.validate_capped_row(leaking)
        self.assertEqual(outcome, "cpp_abstain")
        self.assertIn("abstain_leaks_best_score_fraction", issues)

    def test_read_af_parity_detects_exact_score_drift(self) -> None:
        oracle = [
            {
                "region": "ranked",
                "af_status": "usable",
                "top_score_fraction": "0/1",
                "n_top_exact": 3,
                "read_af_fractions": ["1/2", "1/2"],
            },
            {
                "region": "zero",
                "af_status": "zero_denominator",
            },
        ]
        cpp = {
            "ranked": {
                "read_af_status": "ranked_complete",
                "best_score_fraction": "0/1",
                "best_tree_tie_count": "3",
                "best_tree_unique": False,
                "af_coverage": [
                    {"fraction": "1/2"},
                    {"fraction": "1/2"},
                ],
            },
            "zero": {
                "read_af_status": "zero_denominator",
                "best_score_fraction": None,
                "best_tree_tie_count": None,
                "best_tree_unique": None,
            },
        }
        counts: Counter[str] = Counter()
        examples = []

        summary, _ = gate.compare_read_af(
            oracle, cpp, counts, examples, limit=10
        )

        self.assertEqual(summary["mismatch_fields"], 0)
        self.assertFalse(counts)

        drifted = copy.deepcopy(cpp)
        drifted["ranked"]["best_score_fraction"] = "1/1"
        drift_counts: Counter[str] = Counter()
        drift_examples = []
        drift_summary, _ = gate.compare_read_af(
            oracle, drifted, drift_counts, drift_examples, limit=10
        )
        self.assertEqual(drift_summary["mismatch_fields"], 1)
        self.assertEqual(
            drift_counts["read_af:best_score_fraction"],
            1,
        )


if __name__ == "__main__":
    unittest.main()
