#!/usr/bin/env python3

import math
import pathlib
import sys
import types
import unittest
from unittest import mock

import numpy as np


BASE = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = BASE / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import compressed_vaf_rank_probe as probe  # noqa: E402
import run_compressed_vaf_rank_probe as runner  # noqa: E402
from compressed_vaf_rank_probe import (  # noqa: E402
    PerfectTracebackDag,
    brute_force_perfect_vertex_sets,
    exhaustive_current_rank,
    rank_perfect_family_lazy,
)


def asymmetric_quality_groups(k):
    groups = []
    for alt_bit in range(k):
        pattern = "".join(
            "A" if bit == alt_bit else "R" for bit in range(k)
        )
        groups.append((pattern, tuple(35 for _ in range(k)), k - alt_bit))
    groups.append(("A" * k, tuple(35 for _ in range(k)), 1))
    return tuple(groups)


class PerfectTracebackDagOracleTest(unittest.TestCase):
    def assert_dag_matches_independent_oracle(
        self,
        *,
        full_patterns,
        partial_patterns,
        k,
        structural_alt_universe_mask,
        require_branching=False,
    ):
        dag = PerfectTracebackDag(
            full_patterns=full_patterns,
            partial_patterns=partial_patterns,
            k=k,
            structural_alt_universe_mask=structural_alt_universe_mask,
        )
        branches = dag.branches()
        seen = {}
        observed_branching = False
        for branch in branches:
            branch_count = 0
            possible = set(branch.possible_vertices)
            for parents in dag.iter_branch(branch):
                branch_count += 1
                vertices = dag.vertex_set(parents)
                self.assertTrue(set(vertices) <= possible)
                self.assertNotIn(vertices, seen)
                seen[vertices] = branch.identity()
                child_counts = {}
                for parent in parents:
                    child_counts[parent] = child_counts.get(parent, 0) + 1
                observed_branching |= any(
                    count > 1 for count in child_counts.values()
                )
            self.assertEqual(branch_count, branch.completion_count)
        brute = brute_force_perfect_vertex_sets(
            full_patterns=full_patterns,
            partial_patterns=partial_patterns,
            k=k,
            structural_alt_universe_mask=structural_alt_universe_mask,
        )
        self.assertEqual(tuple(sorted(seen)), brute)
        self.assertEqual(
            sum(branch.completion_count for branch in branches),
            dag.forest_count(dag.full_mask),
        )
        self.assertEqual(len(seen), dag.forest_count(dag.full_mask))
        if require_branching:
            self.assertTrue(observed_branching)

    def test_m1_to_m5_traceback_family_equals_independent_parent_oracle(self):
        for m in range(1, 6):
            with self.subTest(m=m):
                full = []
                partial = ["A" * m]
                dag = PerfectTracebackDag(
                    full_patterns=full,
                    partial_patterns=partial,
                    k=m,
                    structural_alt_universe_mask=(1 << m) - 1,
                )
                traceback = tuple(sorted(dag.iter_all_vertex_sets()))
                brute = brute_force_perfect_vertex_sets(
                    full_patterns=full,
                    partial_patterns=partial,
                    k=m,
                    structural_alt_universe_mask=(1 << m) - 1,
                )
                self.assertEqual(traceback, brute)
                self.assertEqual(len(traceback), dag.forest_count(dag.full_mask))

    def test_structural_fixture_matrix_branch_partition_and_containment(self):
        fixtures = (
            {
                "name": "branching",
                "full_patterns": [],
                "partial_patterns": ["AXX", "XAX", "XXA"],
                "k": 3,
                "structural_alt_universe_mask": 0b111,
                "require_branching": True,
            },
            {
                "name": "mandatory_full",
                "full_patterns": ["AARR"],
                "partial_patterns": ["XXAA"],
                "k": 4,
                "structural_alt_universe_mask": 0b1111,
            },
            {
                "name": "mixed_rax",
                "full_patterns": [],
                "partial_patterns": ["ARX", "XAA", "RXA"],
                "k": 3,
                "structural_alt_universe_mask": 0b111,
            },
            {
                "name": "gapped_active",
                "full_patterns": [],
                "partial_patterns": ["AXXXX", "XXAXX", "XXXXA"],
                "k": 5,
                "structural_alt_universe_mask": 0b10101,
                "require_branching": True,
            },
        )
        for fixture in fixtures:
            with self.subTest(fixture=fixture["name"]):
                kwargs = dict(fixture)
                kwargs.pop("name")
                self.assert_dag_matches_independent_oracle(**kwargs)

    def test_recurrence_required_family_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "recurrence-aware"):
            PerfectTracebackDag(
                full_patterns=[],
                partial_patterns=["AA", "AR", "RA"],
                k=2,
                structural_alt_universe_mask=3,
            )

    def test_top_level_branches_conserve_exact_count(self):
        dag = PerfectTracebackDag(
            full_patterns=[],
            partial_patterns=["AXX", "XAX", "XXA"],
            k=3,
            structural_alt_universe_mask=7,
        )
        self.assertEqual(
            sum(branch.completion_count for branch in dag.branches()),
            dag.forest_count(dag.full_mask),
        )


class CurrentLikelihoodEquivalenceTest(unittest.TestCase):
    def test_m1_to_m5_best_score_tree_and_complete_ties_match(self):
        for m in range(1, 6):
            with self.subTest(m=m):
                partial = ["A" * m]
                quality = asymmetric_quality_groups(m)
                brute = brute_force_perfect_vertex_sets(
                    full_patterns=[],
                    partial_patterns=partial,
                    k=m,
                    structural_alt_universe_mask=(1 << m) - 1,
                )
                oracle = exhaustive_current_rank(
                    vertex_sets=brute,
                    quality_groups=quality,
                    k=m,
                    tie_tolerance=1e-6,
                )
                lazy = rank_perfect_family_lazy(
                    full_patterns=[],
                    partial_patterns=partial,
                    k=m,
                    quality_groups=quality,
                    structural_alt_universe_mask=(1 << m) - 1,
                    tie_tolerance=1e-6,
                    max_candidates=10_000,
                    deadline_seconds=30,
                    max_tie_ids=10_000,
                    enable_numerical_pruning=False,
                )
                self.assertTrue(lazy["ranking_complete"])
                self.assertTrue(lazy["machine_exact_by_full_enumeration"])
                self.assertFalse(lazy["numerical_bound_certified"])
                self.assertEqual(lazy["pruned_candidate_count"], 0)
                self.assertEqual(
                    lazy["logical_family_count"], oracle["candidate_count"]
                )
                self.assertAlmostEqual(
                    lazy["best_log_likelihood"],
                    oracle["best_log_likelihood"],
                    places=10,
                )
                self.assertEqual(
                    lazy["winner_vertex_set_ids"],
                    oracle["winner_vertex_set_ids"],
                )
                self.assertEqual(lazy["tie_count"], oracle["tie_count"])

    def test_symmetric_two_chain_case_returns_complete_two_way_tie(self):
        quality = (
            ("AR", (40, 40), 3),
            ("RA", (40, 40), 3),
            ("AA", (40, 40), 1),
        )
        result = rank_perfect_family_lazy(
            full_patterns=[],
            partial_patterns=["AA"],
            k=2,
            quality_groups=quality,
            structural_alt_universe_mask=3,
            max_candidates=10,
            deadline_seconds=10,
        )
        self.assertTrue(result["ranking_complete"])
        self.assertEqual(result["logical_family_count"], 2)
        self.assertEqual(result["tie_count"], 2)

    def test_float_pruning_is_diagnostic_and_never_machine_exact(self):
        result = rank_perfect_family_lazy(
            full_patterns=[],
            partial_patterns=["AAAAA"],
            k=5,
            quality_groups=asymmetric_quality_groups(5),
            structural_alt_universe_mask=31,
            max_candidates=1000,
            deadline_seconds=10,
            bound_mode="union_mixture_float",
            enable_numerical_pruning=True,
        )
        self.assertGreater(result["pruned_candidate_count"], 0)
        self.assertTrue(result["traceback_accounting_complete"])
        self.assertFalse(result["numerical_bound_certified"])
        self.assertFalse(result["search_complete"])
        self.assertFalse(result["tie_class_complete"])
        self.assertFalse(result["ranking_complete"])
        self.assertEqual(
            result["status"],
            "INCOMPLETE_UNCERTIFIED_FLOAT_BOUND_PRUNING",
        )
        self.assertIsNone(result["best_log_likelihood"])
        self.assertEqual(result["winner_vertex_set_ids"], [])
        self.assertFalse(result["diagnostic_incumbent"]["authoritative"])

    def test_candidate_cap_is_fail_closed_and_publishes_no_winner(self):
        result = rank_perfect_family_lazy(
            full_patterns=[],
            partial_patterns=["AAAAA"],
            k=5,
            quality_groups=asymmetric_quality_groups(5),
            structural_alt_universe_mask=31,
            max_candidates=2,
            deadline_seconds=10,
            enable_numerical_pruning=False,
        )
        self.assertEqual(result["status"], "INCOMPLETE_MAX_CANDIDATES")
        self.assertFalse(result["search_complete"])
        self.assertFalse(result["ranking_complete"])
        self.assertIsNone(result["best_log_likelihood"])
        self.assertEqual(result["winner_vertex_set_ids"], [])
        self.assertIsNotNone(result["diagnostic_incumbent"])
        self.assertFalse(result["diagnostic_incumbent"]["authoritative"])

    def test_tie_output_cap_is_fail_closed(self):
        quality = (
            ("AR", (40, 40), 3),
            ("RA", (40, 40), 3),
            ("AA", (40, 40), 1),
        )
        result = rank_perfect_family_lazy(
            full_patterns=[],
            partial_patterns=["AA"],
            k=2,
            quality_groups=quality,
            structural_alt_universe_mask=3,
            max_candidates=10,
            deadline_seconds=10,
            max_tie_ids=1,
            enable_numerical_pruning=False,
        )
        self.assertEqual(result["status"], "INCOMPLETE_TIE_CLASS_OUTPUT_CAP")
        self.assertTrue(result["search_complete"])
        self.assertFalse(result["tie_class_complete"])
        self.assertFalse(result["ranking_complete"])
        self.assertEqual(result["winner_vertex_set_ids"], [])

    def test_rowwise_float_bound_is_explicitly_not_certified(self):
        result = rank_perfect_family_lazy(
            full_patterns=[],
            partial_patterns=["AAAA"],
            k=4,
            quality_groups=asymmetric_quality_groups(4),
            structural_alt_universe_mask=15,
            max_candidates=100,
            deadline_seconds=10,
            bound_mode="rowwise",
            enable_numerical_pruning=True,
        )
        self.assertFalse(result["numerical_bound_certified"])
        if result["pruned_candidate_count"]:
            self.assertFalse(result["ranking_complete"])
            self.assertEqual(result["winner_vertex_set_ids"], [])
        else:
            self.assertTrue(result["machine_exact_by_full_enumeration"])

    def test_control_types_are_strict(self):
        base = {
            "full_patterns": [],
            "partial_patterns": ["AA"],
            "k": 2,
            "quality_groups": asymmetric_quality_groups(2),
            "structural_alt_universe_mask": 3,
        }
        invalid = (
            {"max_candidates": True},
            {"max_candidates": 2.0},
            {"max_tie_ids": True},
            {"max_tie_ids": 2.0},
            {"deadline_seconds": True},
            {"deadline_seconds": np.float64(1.0)},
            {"tie_tolerance": True},
            {"tie_tolerance": np.float64(1e-6)},
            {"enable_numerical_pruning": 1},
        )
        for override in invalid:
            with self.subTest(override=override):
                with self.assertRaises(ValueError):
                    rank_perfect_family_lazy(**base, **override)

    def test_candidate_likelihood_and_gap_must_be_finite(self):
        quality = (("A", (35,), 1),)
        invalid = (
            (math.nan, 0.0),
            (-1.0, math.nan),
            (-1.0, math.inf),
            (-1.0, -1e-9),
        )
        for log_likelihood, gap in invalid:
            with self.subTest(log_likelihood=log_likelihood, gap=gap):
                fake = types.SimpleNamespace(
                    converged=True,
                    monotone=True,
                    log_likelihood=log_likelihood,
                    global_log_likelihood_gap_bound=gap,
                    optimizer_status="TEST",
                )
                with mock.patch.object(
                    probe,
                    "fit_quality_aware_mixture",
                    return_value=fake,
                ):
                    with self.assertRaises(RuntimeError):
                        probe._score_candidate(
                            quality,
                            (0, 1),
                            k=1,
                            minimum_error_rate=1e-6,
                            maximum_error_rate=0.25,
                        )


class RunnerIntegrityTest(unittest.TestCase):
    def test_source_snapshot_includes_hypercube_and_detects_drift(self):
        paths = runner.source_paths()
        self.assertIn("hypercube_exact", paths)
        before = {
            "source": {"path": "/tmp/source.py", "sha256": "a" * 64}
        }
        after = {
            "source": {"path": "/tmp/source.py", "sha256": "b" * 64}
        }
        with self.assertRaisesRegex(RuntimeError, "source"):
            runner.assert_snapshot_unchanged(
                before,
                after,
                label="test sources",
            )

    def test_v2_default_and_numpy_runtime_binding_exist(self):
        self.assertIn(
            "compressed_vaf_rank_probe_v2",
            runner.DEFAULT_OUTPUT.parts,
        )
        self.assertRegex(np.__version__, r"^\d+\.\d+")


if __name__ == "__main__":
    unittest.main()
