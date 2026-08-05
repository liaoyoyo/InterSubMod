#!/usr/bin/env python3
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_m0_likelihood_census import (  # noqa: E402
    _selection_status,
    edge_set_digest,
    node_to_vertex,
    score_lineage,
    tree_vertex_set,
)


class CandidateParsingTests(unittest.TestCase):
    def test_node_and_vertex_set_parsing(self):
        tree = {"edges": [["ROOT", "H_AR"], ["H_AR", "AA"]], "recurrence": []}
        self.assertEqual(node_to_vertex("ROOT", 2), 0)
        self.assertEqual(node_to_vertex("H_AR", 2), 1)
        self.assertEqual(node_to_vertex("AA", 2), 3)
        self.assertEqual(tree_vertex_set(tree, 2), (0, 1, 3))

    def test_edge_digest_is_order_invariant(self):
        left = {"edges": [["ROOT", "H_AR"], ["H_AR", "AA"]], "recurrence": [1]}
        right = {"edges": [["H_AR", "AA"], ["ROOT", "H_AR"]], "recurrence": [1]}
        self.assertEqual(edge_set_digest(left), edge_set_digest(right))


class M0ScoringTests(unittest.TestCase):
    def _region(self):
        return {"region": "chr1:10-20", "chrom": "chr1", "start": 10, "end": 20, "n_sSNV": 2}

    def test_same_vertex_edges_collapse_and_remain_unresolved(self):
        lineage = {
            "family": "1",
            "display_trees_complete": True,
            "n_reads": 100,
            "n_trees": 2,
            "n_trees_stored": 2,
            "obs_subreads": {"RR": 60, "AR": 30, "AA": 10},
            "trees": [
                {"edges": [["ROOT", "H_AR"], ["H_AR", "AA"]], "recurrence": []},
                {"edges": [["ROOT", "H_RA"], ["H_RA", "AA"]], "recurrence": []},
            ],
        }
        row = score_lineage("S", self._region(), lineage, error_rate=0.01, tie_tolerance=1e-6, frozen_solver=None)
        self.assertEqual(row["raw_tree_candidates_T"], 2)
        self.assertEqual(row["distinct_vertex_sets_V"], 2)
        self.assertEqual(row["best_vertex_sets"], 1)
        self.assertEqual(row["selection_status"], "LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE")

    def test_identical_vertex_set_edge_variants_score_once(self):
        lineage = {
            "family": "2",
            "display_trees_complete": True,
            "n_reads": 100,
            "n_trees": 2,
            "n_trees_stored": 2,
            "obs_subreads": {"RRR": 50, "ARR": 25, "AAR": 15, "AAA": 10},
            "trees": [
                {"edges": [["ROOT", "H_ARR"], ["H_ARR", "H_AAR"], ["H_AAR", "AAA"]], "recurrence": []},
                {"edges": [["ROOT", "H_ARR"], ["H_ARR", "H_AAR"], ["H_AAR", "AAA"]], "recurrence": [2]},
            ],
        }
        row = score_lineage("S", {**self._region(), "n_sSNV": 3}, lineage, error_rate=0.01, tie_tolerance=1e-6, frozen_solver=None)
        self.assertEqual(row["distinct_vertex_sets_V"], 1)
        self.assertEqual(row["top_edge_variants"], 2)
        self.assertEqual(row["selection_status"], "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED")

    def test_selection_status_contract(self):
        self.assertEqual(_selection_status(1, 1, 1, 1), "T1_CANDIDATE_UNIQUE")
        self.assertEqual(
            _selection_status(4, 2, 2, 4),
            "LIKELIHOOD_TIED_VERTEX_SETS",
        )


if __name__ == "__main__":
    unittest.main()
