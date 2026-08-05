#!/usr/bin/env python3
"""Regression and oracle tests for the optimized bounded research backend."""

from __future__ import annotations

import itertools
import math
import pathlib
import random
import sys
import unittest


HERE = pathlib.Path(__file__).resolve().parent
PROBE_ROOT = HERE.parent
sys.path.insert(0, str(PROBE_ROOT / "scripts"))

from optimized_hypercube_backend import (  # noqa: E402
    BitsetObligationBnb,
    evaluate_evolutionary_mode,
    fixed_node_parent_mapping,
    optimal_family_digest,
    prepare_bitset_problem,
    project_unary_hidden_chains,
    solve_group_terminal_subset_dp,
    solve_with_certificates,
)
from solver_probe import brute_force_optimal, build_instance  # noqa: E402


def family_digest(instance, vertex_sets) -> str:
    return optimal_family_digest(instance.k, vertex_sets)


def strict_eligibility() -> dict:
    return {
        "allele_specific_cn_known": True,
        "copy_neutral_diploid_both_homologs_retained": True,
        "total_cn": 2,
        "major_cn": 1,
        "minor_cn": 1,
        "loh_absent": True,
        "clonal_loh_absent": True,
        "subclonal_loh_absent": True,
        "clonal_deletion_absent": True,
        "subclonal_deletion_absent": True,
        "clonal_amplification_absent": True,
        "subclonal_amplification_absent": True,
        "clonal_cna_absent": True,
        "subclonal_cna_absent": True,
        "wgd_boundary_absent": True,
        "mutation_loss_absent": True,
        "phasing_qc_pass": True,
        "mapping_qc_pass": True,
        "base_quality_qc_pass": True,
        "strand_qc_pass": True,
        "read_independence_qc_pass": True,
        "allele_specific_duplicated_copies_absent": True,
        "somatic_variant_qc_pass": True,
    }


class OptimizedStructuralSolverTest(unittest.TestCase):
    def assert_matches_brute(self, full, partial, k):
        instance = build_instance(full, partial, k)
        brute = brute_force_optimal(instance, combination_cap=2_000_000)
        self.assertTrue(brute["complete"])
        problem = prepare_bitset_problem(instance)
        dp = solve_group_terminal_subset_dp(problem, max_q=8)
        self.assertTrue(dp.objective_certified)
        self.assertEqual(brute["objective"], dp.objective)
        bnb = BitsetObligationBnb(
            problem,
            time_limit_seconds=10.0,
            certified_target_objective=dp.objective,
        ).run()
        self.assertTrue(bnb.objective_certified)
        self.assertTrue(bnb.family_complete)
        self.assertEqual("CANDIDATE_SET_COMPLETE", bnb.status)
        self.assertEqual(brute["objective"], bnb.objective)
        self.assertEqual(
            family_digest(instance, brute["vertex_sets"]),
            bnb.family_digest,
        )

    def test_declared_oracles_match_objective_and_full_family(self):
        cases = [
            ([], [], 2),
            (["AA"], [], 2),
            ([], ["AX", "XA"], 2),
            (["AAA"], [], 3),
            (["RRA", "AAA"], ["RAX"], 3),
            ([], ["AAA"], 3),
            (["AR", "RA", "AA"], [], 2),
            (["AAAA"], [], 4),
        ]
        for full, partial, k in cases:
            with self.subTest(full=full, partial=partial, k=k):
                self.assert_matches_brute(full, partial, k)

    def test_seeded_k_le_4_matches_independent_bruteforce(self):
        rng = random.Random(20260718)
        checked = 0
        for _ in range(180):
            k = rng.randint(1, 4)
            full_pool = [
                "".join(symbols)
                for symbols in itertools.product("RA", repeat=k)
            ]
            partial_pool = [
                "".join(symbols)
                for symbols in itertools.product("RAX", repeat=k)
                if "X" in symbols
            ]
            full = rng.sample(full_pool, rng.randint(0, min(3, len(full_pool))))
            partial = rng.sample(
                partial_pool,
                rng.randint(0, min(4, len(partial_pool))),
            )
            instance = build_instance(full, partial, k)
            brute = brute_force_optimal(instance, combination_cap=500_000)
            if not brute["complete"]:
                continue
            problem = prepare_bitset_problem(instance)
            dp = solve_group_terminal_subset_dp(problem, max_q=8)
            self.assertEqual(brute["objective"], dp.objective, (full, partial, k))
            bnb = BitsetObligationBnb(
                problem,
                time_limit_seconds=10.0,
                certified_target_objective=dp.objective,
            ).run()
            self.assertTrue(bnb.family_complete, (full, partial, k, bnb.status))
            self.assertEqual(
                family_digest(instance, brute["vertex_sets"]),
                bnb.family_digest,
                (full, partial, k),
            )
            checked += 1
        self.assertGreaterEqual(checked, 150)

    def test_duplicate_dominance_required_hit_and_singleton_reductions(self):
        instance = build_instance(
            ["ARR"],
            ["AXX", "AXX", "AAX", "ARA"],
            3,
        )
        brute = brute_force_optimal(instance)
        problem = prepare_bitset_problem(instance)
        self.assertGreaterEqual(problem.reduction.duplicate_groups_removed, 1)
        self.assertGreaterEqual(problem.reduction.required_hit_groups_removed, 1)
        self.assertGreaterEqual(problem.reduction.singleton_groups_forced, 1)
        dp = solve_group_terminal_subset_dp(problem)
        bnb = BitsetObligationBnb(
            problem,
            certified_target_objective=dp.objective,
        ).run()
        self.assertEqual(brute["objective"], dp.objective)
        self.assertEqual(
            family_digest(instance, brute["vertex_sets"]),
            bnb.family_digest,
        )

    def test_small_q_route_is_objective_only(self):
        instance = build_instance([], ["AX", "XA"], 2)
        result = solve_group_terminal_subset_dp(prepare_bitset_problem(instance))
        self.assertEqual("OPTIMAL_VALUE_CERTIFIED", result.status)
        self.assertTrue(result.objective_certified)
        self.assertFalse(result.family_complete)
        self.assertEqual(2, result.objective)

    def test_q_limit_routes_away_without_false_certificate(self):
        instance = build_instance(
            [],
            ["AXXX", "XAXX", "XXAX", "XXXA"],
            4,
        )
        result = solve_group_terminal_subset_dp(
            prepare_bitset_problem(instance),
            max_q=2,
        )
        self.assertEqual("ROUTE_NOT_ELIGIBLE_Q_GT_MAX", result.status)
        self.assertFalse(result.objective_certified)
        self.assertIsNone(result.objective)

    def test_dp_resource_gate_routes_away_before_allocation(self):
        instance = build_instance(["AAAA"], [], 4)
        result = solve_group_terminal_subset_dp(
            prepare_bitset_problem(instance),
            max_table_cells=1,
        )
        self.assertEqual("ROUTE_NOT_ELIGIBLE_RESOURCE_GATE", result.status)
        self.assertFalse(result.objective_certified)
        self.assertGreater(result.stats.table_cells, result.stats.max_table_cells)

    def test_max_set_cap_preserves_objective_but_blocks_ranking(self):
        instance = build_instance(["AAAA"], [], 4)
        receipt = solve_with_certificates(
            instance,
            max_sets=1,
            time_limit_seconds=10.0,
        )
        objective = receipt["objective_certificate"]
        family = receipt["family_certificate"]
        self.assertTrue(objective["objective_certified"])
        self.assertEqual(3, objective["objective"])
        self.assertEqual("CANDIDATE_SET_INCOMPLETE_MAX_SETS", family["status"])
        self.assertFalse(family["family_complete"])
        self.assertFalse(receipt["ranking_allowed"])
        self.assertEqual("ABSTAIN_INCOMPLETE_FAMILY", receipt["failure_policy"])

    def test_legacy_cap_counterexample_is_resolved_by_objective_first_route(self):
        # The legacy one-pass probe can encounter a cost-3 set before the true
        # h*=2 family for this case.  The small-q DP proves h* first; the true
        # optimum family has one set and can therefore complete at max_sets=1.
        instance = build_instance(["RAA"], ["ARX", "AXA"], 3)
        brute = brute_force_optimal(instance)
        self.assertEqual(2, brute["objective"])
        receipt = solve_with_certificates(
            instance,
            max_sets=1,
            time_limit_seconds=10.0,
        )
        objective = receipt["objective_certificate"]
        family = receipt["family_certificate"]
        self.assertEqual(2, objective["objective"])
        self.assertEqual(2, family["objective"])
        self.assertNotEqual(3, family["objective"])
        self.assertEqual("CANDIDATE_SET_COMPLETE", family["status"])
        self.assertTrue(receipt["ranking_allowed"])

    def test_total_deadline_before_objective_blocks_ranking(self):
        instance = build_instance(["AAAA"], [], 4)
        receipt = solve_with_certificates(
            instance,
            max_sets=None,
            time_limit_seconds=1e-12,
        )
        objective = receipt["objective_certificate"]
        family = receipt["family_certificate"]
        self.assertEqual("OBJECTIVE_INCOMPLETE_DEADLINE", objective["status"])
        self.assertFalse(objective["objective_certified"])
        self.assertEqual("NO_FEASIBLE_CERTIFICATE_DEADLINE", family["status"])
        self.assertFalse(family["family_complete"])
        self.assertFalse(receipt["ranking_allowed"])
        self.assertEqual("prepare+subset_dp+bitset_bnb", receipt["deadline_scope"])

    def test_bnb_deadline_after_dp_preserves_objective_but_blocks_ranking(self):
        problem = prepare_bitset_problem(build_instance(["AAAA"], [], 4))
        objective = solve_group_terminal_subset_dp(problem)
        family = BitsetObligationBnb(
            problem,
            time_limit_seconds=1e-12,
            certified_target_objective=objective.objective,
        ).run()
        self.assertEqual("CANDIDATE_SET_INCOMPLETE_DEADLINE", family.status)
        self.assertTrue(family.objective_certified)
        self.assertEqual(3, family.objective)
        self.assertFalse(family.family_complete)

    def test_direct_bnb_without_incumbent_does_not_claim_feasible(self):
        problem = prepare_bitset_problem(build_instance(["AAAA"], [], 4))
        capped = BitsetObligationBnb(
            problem,
            max_sets=0,
            time_limit_seconds=10.0,
        ).run()
        deadline = BitsetObligationBnb(
            problem,
            max_sets=None,
            time_limit_seconds=1e-12,
        ).run()
        self.assertEqual(
            "NO_FEASIBLE_CERTIFICATE_MAX_SETS",
            capped.status,
        )
        self.assertEqual(
            "NO_FEASIBLE_CERTIFICATE_DEADLINE",
            deadline.status,
        )
        for result in (capped, deadline):
            self.assertFalse(result.objective_certified)
            self.assertIsNone(result.objective)
            self.assertIsNone(result.incumbent_objective)
            self.assertFalse(result.family_complete)

    def test_nonfinite_deadlines_are_rejected_in_every_public_route(self):
        problem = prepare_bitset_problem(build_instance(["AA"], [], 2))
        for invalid in (math.nan, math.inf, -math.inf):
            with self.subTest(route="bnb", invalid=invalid):
                with self.assertRaises(ValueError):
                    BitsetObligationBnb(
                        problem,
                        time_limit_seconds=invalid,
                    )
            with self.subTest(route="dp", invalid=invalid):
                with self.assertRaises(ValueError):
                    solve_group_terminal_subset_dp(
                        problem,
                        time_limit_seconds=invalid,
                    )
            with self.subTest(route="wrapper", invalid=invalid):
                with self.assertRaises(ValueError):
                    solve_with_certificates(
                        build_instance(["AA"], [], 2),
                        time_limit_seconds=invalid,
                    )

    def test_direct_dp_cannot_certify_after_tiny_deadline(self):
        problem = prepare_bitset_problem(build_instance(["AAAA"], [], 4))
        result = solve_group_terminal_subset_dp(
            problem,
            time_limit_seconds=1e-12,
        )
        self.assertEqual("OBJECTIVE_INCOMPLETE_DEADLINE", result.status)
        self.assertFalse(result.objective_certified)
        self.assertIsNone(result.objective)


class ParentMappingTest(unittest.TestCase):
    def test_fixed_node_additive_score_without_cartesian_materialization(self):
        selected = {0, 1, 2, 3}
        scores = {
            (0, 1): 0.5,
            (0, 2): 0.4,
            (1, 3): 2.0,
            (2, 3): 1.0,
        }
        result = fixed_node_parent_mapping(selected, edge_scores=scores)
        self.assertEqual(2, result["total_parent_tree_count"])
        self.assertEqual(1, result["best_parent_tree_count"])
        self.assertAlmostEqual(2.9, result["best_total_score"])
        self.assertEqual(1, result["canonical_best_parent"]["3"])
        self.assertFalse(
            result["complexity"]["parent_cartesian_product_materialized"]
        )

    def test_fixed_node_tie_count_and_softmax(self):
        selected = {0, 1, 2, 3}
        scores = {
            (0, 1): 0.0,
            (0, 2): 0.0,
            (1, 3): 1.0,
            (2, 3): 1.0,
        }
        result = fixed_node_parent_mapping(selected, edge_scores=scores)
        self.assertEqual(2, result["best_parent_tree_count"])
        self.assertFalse(result["unique_best_tree"])
        posterior = result["parent_posterior"]["3"]
        self.assertAlmostEqual(0.5, posterior["1"])
        self.assertAlmostEqual(0.5, posterior["2"])

    def test_missing_edge_score_fails_closed(self):
        with self.assertRaises(ValueError):
            fixed_node_parent_mapping(
                {0, 1, 2, 3},
                edge_scores={(0, 1): 0.0, (0, 2): 0.0, (1, 3): 1.0},
            )

    def test_invalid_tolerance_or_scaled_score_fails_closed(self):
        scores = {
            (0, 1): 0.0,
            (0, 2): 0.0,
            (1, 3): 1.0,
            (2, 3): 1.0,
        }
        for tolerance in (-1.0, math.inf, math.nan):
            with self.subTest(tolerance=tolerance):
                with self.assertRaises(ValueError):
                    fixed_node_parent_mapping(
                        {0, 1, 2, 3},
                        edge_scores=scores,
                        tie_tolerance=tolerance,
                    )
        overflow_scores = dict(scores)
        overflow_scores[(1, 3)] = 1e308
        with self.assertRaises(ValueError):
            fixed_node_parent_mapping(
                {0, 1, 2, 3},
                edge_scores=overflow_scores,
                beta=2.0,
            )
        with self.assertRaises(ValueError):
            fixed_node_parent_mapping(
                {0, 1, 3},
                edge_scores={(0, 1): 1e308, (1, 3): 1e308},
                beta=1.0,
            )

    def test_exhaustive_parent_best_matches_factorized_result(self):
        rng = random.Random(1729)
        for k in range(1, 5):
            vertices = range(1 << k)
            for _ in range(40):
                selected = {0}
                for vertex in vertices:
                    if vertex and rng.random() < 0.45:
                        selected.add(vertex)
                if any(
                    not any(
                        (vertex ^ (1 << bit)) in selected
                        for bit in range(k)
                        if vertex & (1 << bit)
                    )
                    for vertex in selected - {0}
                ):
                    continue
                parent_lists = {
                    child: [
                        child ^ (1 << bit)
                        for bit in range(k)
                        if child & (1 << bit)
                        and (child ^ (1 << bit)) in selected
                    ]
                    for child in selected - {0}
                }
                scores = {
                    (parent, child): rng.randint(-3, 3) / 2.0
                    for child, parents in parent_lists.items()
                    for parent in parents
                }
                explicit = [
                    sum(scores[(parent, child)] for child, parent in zip(parent_lists, choices))
                    for choices in itertools.product(*parent_lists.values())
                ]
                result = fixed_node_parent_mapping(selected, edge_scores=scores)
                self.assertAlmostEqual(max(explicit), result["best_total_score"])
                self.assertEqual(
                    sum(abs(value - max(explicit)) <= 1e-12 for value in explicit),
                    result["best_parent_tree_count"],
                )


class ProjectionAndBiologyGateTest(unittest.TestCase):
    def test_unary_hidden_chain_becomes_unordered_multi_mutation_edge(self):
        projection = project_unary_hidden_chains(
            {0, 1, 3, 7},
            {1: 0, 3: 1, 7: 3},
            {0, 7},
        )
        self.assertEqual([1, 3], projection["collapsed_vertices"])
        self.assertEqual(1, len(projection["projected_edges"]))
        edge = projection["projected_edges"][0]
        self.assertEqual("MULTI_MUTATION_EDGE_EQUIVALENCE", edge["edge_type"])
        self.assertEqual("UNRESOLVED_NO_READ_EVIDENCE", edge["mutation_order"])
        self.assertEqual([0, 1, 2], edge["mutation_bits"])
        self.assertEqual(6, edge["unresolved_order_count"])
        self.assertFalse(projection["objective_changed"])
        self.assertFalse(
            projection["projected_model_objective_equivalence_claimed"]
        )

    def test_branchpoint_is_not_collapsed(self):
        projection = project_unary_hidden_chains(
            {0, 1, 3, 5},
            {1: 0, 3: 1, 5: 1},
            {0, 3, 5},
        )
        self.assertIn(1, projection["kept_vertices"])
        self.assertNotIn(1, projection["collapsed_vertices"])

    def test_strict_mode_requires_cn_loh_gate(self):
        result = evaluate_evolutionary_mode(
            {0, 1, 2, 3},
            k=2,
            mode="M1_STRICT_INFINITE_SITES",
            eligibility={},
        )
        self.assertEqual("ABSTAIN_CN_LOH_GATE", result["status"])
        self.assertFalse(result["constraint_applied"])

    def test_strict_diamond_is_infeasible_only_after_gate(self):
        gates = strict_eligibility()
        diamond = evaluate_evolutionary_mode(
            {0, 1, 2, 3},
            k=2,
            mode="M1_STRICT_INFINITE_SITES",
            eligibility=gates,
        )
        chain = evaluate_evolutionary_mode(
            {0, 1, 3},
            k=2,
            mode="M1_STRICT_INFINITE_SITES",
            eligibility=gates,
        )
        self.assertEqual("STRICT_INFEASIBLE", diamond["status"])
        self.assertEqual("STRICT_COMPATIBLE", chain["status"])
        self.assertFalse(diamond["publication_allowed"])
        self.assertTrue(chain["publication_allowed"])

        disconnected = evaluate_evolutionary_mode(
            {0, 3},
            k=2,
            mode="M1_STRICT_INFINITE_SITES",
            eligibility=gates,
        )
        self.assertEqual("STRICT_INFEASIBLE", disconnected["status"])
        self.assertFalse(disconnected["publication_allowed"])
        self.assertEqual(
            [3],
            disconnected["orphan_states_without_selected_predecessor"],
        )

    def test_strict_mode_rejects_invalid_states_and_adverse_cn_evidence(self):
        with self.assertRaises(ValueError):
            evaluate_evolutionary_mode(
                {1},
                k=1,
                mode="M1_STRICT_INFINITE_SITES",
                eligibility=strict_eligibility(),
            )
        with self.assertRaises(ValueError):
            evaluate_evolutionary_mode(
                {0, 4},
                k=2,
                mode="M1_STRICT_INFINITE_SITES",
                eligibility=strict_eligibility(),
            )
        adverse = strict_eligibility()
        adverse.update({"minor_cn": 0, "clonal_loh": True, "deletion": True})
        result = evaluate_evolutionary_mode(
            {0, 1},
            k=1,
            mode="M1_STRICT_INFINITE_SITES",
            eligibility=adverse,
        )
        self.assertEqual("ABSTAIN_CN_LOH_GATE", result["status"])
        self.assertIn("minor_cn", result["invalid_or_missing_copy_state"])
        self.assertIn("clonal_loh", result["contradictory_adverse_evidence"])
        self.assertFalse(result["publication_allowed"])

    def test_dollo_mode_remains_explicitly_unresolved(self):
        result = evaluate_evolutionary_mode(
            {0, 1},
            k=1,
            mode="M2_LOSS_SUPPORTED_DOLLO",
        )
        self.assertEqual(
            "UNRESOLVED_LOSS_SUPPORTED_MODEL_NOT_IMPLEMENTED",
            result["status"],
        )
        self.assertFalse(result["publication_allowed"])


if __name__ == "__main__":
    unittest.main()
