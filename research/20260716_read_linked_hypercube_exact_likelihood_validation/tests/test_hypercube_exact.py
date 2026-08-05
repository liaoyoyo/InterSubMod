#!/usr/bin/env python3
import hashlib
import itertools
import json
import math
import pathlib
import random
import sys
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import hypercube_exact as hypercube_module  # noqa: E402
from hypercube_exact import (  # noqa: E402
    SymbolicPattern,
    enumerate_optimal_vertex_sets,
    fit_emission_mixture_certified,
    fit_vertex_mixture,
    fit_vertex_mixture_slsqp,
    mixture_kkt_certificate,
    observed_alt_universe,
    parse_full,
    solve_min_hidden,
    submasks,
)


def exhaustive_minimum_vertex_sets(full_patterns, partial_patterns, k, universe_mode):
    """Independent small-k oracle; deliberately does not use MILP reductions."""
    universe_mask = (
        observed_alt_universe(full_patterns, partial_patterns)
        if universe_mode == "observed_alt"
        else (1 << k) - 1
    )
    vertices = submasks(universe_mask)
    mandatory = {0, *(parse_full(pattern) for pattern in full_patterns)}
    optional = tuple(vertex for vertex in vertices if vertex not in mandatory)
    groups = tuple(SymbolicPattern.from_string(pattern) for pattern in partial_patterns)
    optimum = None
    solutions = set()
    for selected_bits in range(1 << len(optional)):
        selected = mandatory | {
            vertex for offset, vertex in enumerate(optional) if selected_bits & (1 << offset)
        }
        objective = len(selected - mandatory)
        if optimum is not None and objective > optimum:
            continue
        connected = all(
            vertex == 0
            or any(
                vertex & (1 << bit) and (vertex ^ (1 << bit)) in selected
                for bit in range(k)
            )
            for vertex in selected
        )
        covered = all(
            any(group.compatible(vertex, universe_mask) for vertex in selected)
            for group in groups
        )
        if not connected or not covered:
            continue
        if optimum is None or objective < optimum:
            optimum = objective
            solutions = set()
        solutions.add(tuple(sorted(selected)))
    return optimum, solutions


def explicit_joint_completion_worlds(full_patterns, partial_patterns, k, universe_mode):
    """Independent completion-world oracle used only to prove symbolic equivalence.

    A partial completion is *not* promoted to a free full-read terminal: every
    selected state outside the original full-pattern mandatory set keeps its
    extra-vertex cost.  Only worlds attaining the global minimum are pooled.
    """
    universe_mask = (
        observed_alt_universe(full_patterns, partial_patterns)
        if universe_mode == "observed_alt"
        else (1 << k) - 1
    )
    vertices = submasks(universe_mask)
    mandatory = {0, *(parse_full(pattern) for pattern in full_patterns)}
    optional = tuple(vertex for vertex in vertices if vertex not in mandatory)
    groups = tuple(SymbolicPattern.from_string(pattern) for pattern in partial_patterns)
    completions = tuple(group.enumerate_completions(universe_mask) for group in groups)
    assignments = itertools.product(*completions) if completions else [tuple()]
    world_results = []
    for assignment in assignments:
        world_optimum = None
        world_solutions = set()
        required = set(assignment)
        for selected_bits in range(1 << len(optional)):
            selected = mandatory | {
                vertex for offset, vertex in enumerate(optional)
                if selected_bits & (1 << offset)
            }
            if not required <= selected:
                continue
            connected = all(
                vertex == 0
                or any(
                    vertex & (1 << bit) and (vertex ^ (1 << bit)) in selected
                    for bit in range(k)
                )
                for vertex in selected
            )
            if not connected:
                continue
            objective = len(selected - mandatory)
            if world_optimum is None or objective < world_optimum:
                world_optimum = objective
                world_solutions = set()
            if objective == world_optimum:
                world_solutions.add(tuple(sorted(selected)))
        world_results.append((tuple(assignment), world_optimum, world_solutions))
    global_optimum = min(result[1] for result in world_results)
    global_solutions = {
        selected
        for _, objective, solutions in world_results
        if objective == global_optimum
        for selected in solutions
    }
    return global_optimum, global_solutions, world_results


class SymbolicPatternTests(unittest.TestCase):
    def test_symbolic_matches_explicit_for_all_patterns_through_k6(self):
        for k in range(1, 7):
            universe = (1 << k) - 1
            for symbols in itertools.product("RAX", repeat=k):
                pattern = SymbolicPattern.from_string("".join(symbols))
                symbolic = {v for v in submasks(universe) if pattern.compatible(v, universe)}
                explicit = set(pattern.enumerate_completions(universe))
                self.assertEqual(symbolic, explicit)
                self.assertEqual(len(symbolic), 1 << pattern.n_free)

    def test_effective_completion_count_respects_solver_universe(self):
        pattern = SymbolicPattern.from_string("RAX")
        self.assertEqual(pattern.n_completions(), 2)
        self.assertEqual(pattern.n_completions(0b010), 1)


class MilpTests(unittest.TestCase):
    @staticmethod
    def _fake_solution(status, objective, selected=()):
        payload = {
            "status": status,
            "objective": objective,
            "selected_vertices": tuple(selected),
        }
        return SimpleNamespace(
            status=status,
            objective=objective,
            selected_vertices=tuple(selected),
            to_dict=lambda: dict(payload),
        )

    def test_rra_aaa_rax_unique_minimum(self):
        result = enumerate_optimal_vertex_sets(["RRA", "AAA"], ["RAX"], 3, max_sets=20)
        self.assertTrue(result["complete"])
        self.assertEqual(result["objective"], 1)
        self.assertEqual({tuple(v) for v in result["vertex_sets"]}, {(0, 4, 6, 7)})

    def test_rra_aaa_rax_symbolic_matches_global_completion_world_minimum(self):
        objective, expected_sets, worlds = explicit_joint_completion_worlds(
            ["RRA", "AAA"], ["RAX"], 3, "all_loci"
        )
        by_completion = {
            assignment[0]: (world_objective, len(world_sets))
            for assignment, world_objective, world_sets in worlds
        }
        self.assertEqual(by_completion[parse_full("RAR")], (2, 3))
        self.assertEqual(by_completion[parse_full("RAA")], (1, 1))
        actual = enumerate_optimal_vertex_sets(
            ["RRA", "AAA"], ["RAX"], 3,
            universe_mode="all_loci", max_sets=20,
        )
        self.assertTrue(actual["complete"])
        self.assertEqual(actual["objective"], objective)
        self.assertEqual({tuple(values) for values in actual["vertex_sets"]}, expected_sets)

    def test_ax_xa_requires_joint_coverage_and_deduplicated_global_worlds(self):
        objective, expected_sets, worlds = explicit_joint_completion_worlds(
            [], ["AX", "XA"], 2, "all_loci"
        )
        self.assertEqual(len(worlds), 4)
        self.assertEqual(sum(len(world_sets) for _, _, world_sets in worlds), 5)
        self.assertEqual(objective, 2)
        self.assertEqual(len(expected_sets), 3)
        actual = enumerate_optimal_vertex_sets(
            [], ["AX", "XA"], 2,
            universe_mode="all_loci", max_sets=20,
        )
        self.assertTrue(actual["complete"])
        self.assertEqual(actual["objective"], objective)
        self.assertEqual({tuple(values) for values in actual["vertex_sets"]}, expected_sets)
        groups = [SymbolicPattern.from_string(pattern) for pattern in ("AX", "XA")]
        per_pattern_first_success = ({0, parse_full("AR")}, {0, parse_full("RA")})
        self.assertTrue(all(
            not all(any(group.compatible(vertex) for vertex in selected) for group in groups)
            for selected in per_pattern_first_success
        ))

    def test_aaa_has_two_hidden_nodes(self):
        result = solve_min_hidden(["AAA"], [], 3)
        self.assertEqual(result.status, "OPTIMAL")
        self.assertEqual(result.objective, 2.0)
        self.assertEqual(result.mip_gap, 0.0)

    def test_observed_alt_mode_reproduces_x_only_dimension_collapse(self):
        result = solve_min_hidden([], ["RAX"], 3, universe_mode="observed_alt")
        self.assertEqual(result.universe_mask, 0b010)
        self.assertEqual(result.selected_vertices, (0, 2))

    def test_exact_group_reductions_have_auditable_counts(self):
        result = solve_min_hidden(
            [],
            ["AAXX", "AAXX", "AXXX", "XXAX"],
            4,
            universe_mode="all_loci",
        )
        self.assertEqual(result.n_partial_groups_input, 4)
        self.assertEqual(result.n_partial_groups_duplicate_removed, 1)
        self.assertEqual(result.n_partial_groups_dominated_removed, 1)
        self.assertEqual(result.n_partial_groups_required_hit_removed, 0)
        self.assertEqual(result.n_partial_groups_forced_removed, 0)
        self.assertEqual(result.n_partial_groups_active, 2)

    def test_singleton_group_is_forced_but_keeps_its_objective_cost(self):
        result = solve_min_hidden([], ["AAX", "AAA"], 3, universe_mode="all_loci")
        self.assertEqual(result.n_forced_vertices, 1)
        self.assertEqual(result.n_partial_groups_forced_removed, 1)
        self.assertEqual(result.n_partial_groups_required_hit_removed, 1)
        self.assertEqual(result.n_partial_groups_active, 0)
        self.assertEqual(result.objective, 3.0)

    def test_mandatory_states_remove_already_satisfied_groups(self):
        result = solve_min_hidden(["AARR"], ["AAXX", "XXXX"], 4, universe_mode="all_loci")
        self.assertEqual(result.n_partial_groups_required_hit_removed, 2)
        self.assertEqual(result.n_partial_groups_active, 0)

    def test_fixed_objective_uses_sparse_exact_no_good(self):
        first = solve_min_hidden(["AAAA"], [], 4, universe_mode="all_loci")
        second = solve_min_hidden(
            ["AAAA"],
            [],
            4,
            universe_mode="all_loci",
            objective_value=3,
            excluded_vertex_sets=[first.selected_vertices],
        )
        self.assertEqual(second.status, "OPTIMAL")
        self.assertEqual(second.objective, 3.0)
        self.assertEqual(second.n_no_good_nonzeros, 3)
        self.assertNotEqual(second.selected_vertices, first.selected_vertices)

    def test_seeded_small_k_matches_independent_exhaustive_oracle(self):
        rng = random.Random(20260716)
        for k in range(1, 4):
            all_full = ["".join(symbols) for symbols in itertools.product("RA", repeat=k)]
            all_partial = [
                "".join(symbols)
                for symbols in itertools.product("RAX", repeat=k)
                if "X" in symbols
            ]
            for case_index in range(6):
                full = rng.sample(all_full, rng.randint(0, min(2, len(all_full))))
                partial = [rng.choice(all_partial) for _ in range(rng.randint(0, 4))]
                if partial and case_index % 2 == 0:
                    partial.append(partial[0])
                mode = "observed_alt" if case_index % 2 == 0 else "all_loci"
                expected_objective, expected_sets = exhaustive_minimum_vertex_sets(
                    full, partial, k, mode
                )
                actual = enumerate_optimal_vertex_sets(
                    full,
                    partial,
                    k,
                    universe_mode=mode,
                    max_sets=512,
                    time_limit_seconds=5,
                )
                self.assertTrue(actual["complete"], (k, case_index, full, partial, mode))
                self.assertEqual(actual["objective"], expected_objective)
                self.assertEqual(
                    {tuple(values) for values in actual["vertex_sets"]},
                    expected_sets,
                    (k, case_index, full, partial, mode),
                )

    def test_complete_chain_enumeration_keeps_all_24_optima(self):
        with mock.patch.object(
            hypercube_module,
            "_build_problem",
            wraps=hypercube_module._build_problem,
        ) as build_problem, mock.patch.object(
            hypercube_module,
            "_solve_selection",
            wraps=hypercube_module._solve_selection,
        ) as solve_selection:
            result = enumerate_optimal_vertex_sets(
                ["AAAA"], [], 4, universe_mode="all_loci", max_sets=32
            )
        self.assertTrue(result["complete"])
        self.assertEqual(result["objective"], 3)
        self.assertEqual(len(result["vertex_sets"]), math.factorial(4))
        self.assertEqual(build_problem.call_count, 1)
        self.assertEqual(solve_selection.call_count, math.factorial(4) + 1)
        ordered_digest = hashlib.sha256(
            json.dumps(result["vertex_set_ids"], separators=(",", ":")).encode()
        ).hexdigest()
        self.assertEqual(
            ordered_digest,
            "a1ec0da9ba2928269e1154af269891c7f171edc764770baf82148b60c2c93a86",
        )

    def test_zero_nonzero_objective_and_no_good_rows_are_not_dropped(self):
        result = solve_min_hidden(
            [],
            [],
            1,
            universe_mode="observed_alt",
            objective_value=0,
            excluded_vertex_sets=[(0,)],
        )
        self.assertEqual(result.status, "INFEASIBLE")
        self.assertEqual(result.n_constraints, 2)
        self.assertEqual(result.n_no_good_nonzeros, 0)

    def test_duplicate_exclusions_preserve_input_row_count_and_sparse_nnz(self):
        first = solve_min_hidden(["AAAA"], [], 4, universe_mode="all_loci")
        result = solve_min_hidden(
            ["AAAA"],
            [],
            4,
            universe_mode="all_loci",
            objective_value=3,
            excluded_vertex_sets=[
                first.selected_vertices,
                first.selected_vertices,
            ],
        )
        self.assertEqual(result.status, "OPTIMAL")
        self.assertEqual(result.n_no_good_nonzeros, 6)
        self.assertEqual(result.n_constraints, first.n_constraints + 3)

    def test_enumeration_limit_and_completion_statuses_remain_fail_closed(self):
        first_optimal = self._fake_solution("OPTIMAL", 1.0, (0, 1, 3))
        limit = self._fake_solution("LIMIT_REACHED", 1.0, (0, 2, 3))
        infeasible = self._fake_solution("INFEASIBLE", None)
        cases = (
            (
                "first-limit",
                [limit],
                32,
                False,
                0,
                None,
            ),
            (
                "next-limit",
                [first_optimal, limit],
                32,
                False,
                1,
                "LIMIT_REACHED",
            ),
            (
                "next-infeasible",
                [first_optimal, infeasible],
                32,
                True,
                1,
                None,
            ),
            (
                "max-sets-zero",
                [first_optimal],
                0,
                False,
                1,
                "MAX_SETS_REACHED",
            ),
            (
                "max-sets-one",
                [first_optimal],
                1,
                False,
                1,
                "MAX_SETS_REACHED",
            ),
        )
        for (
            label,
            side_effect,
            max_sets,
            expected_complete,
            expected_sets,
            expected_stop,
        ) in cases:
            with self.subTest(label=label), mock.patch.object(
                hypercube_module,
                "_solve_selection",
                side_effect=side_effect,
            ) as solve_selection:
                result = enumerate_optimal_vertex_sets(
                    ["AA"],
                    [],
                    2,
                    universe_mode="all_loci",
                    max_sets=max_sets,
                    time_limit_seconds=7.25,
                )
            self.assertIs(result["complete"], expected_complete)
            self.assertEqual(len(result["vertex_sets"]), expected_sets)
            stop_status = result.get("stop_status")
            observed_stop = (
                stop_status.get("status")
                if isinstance(stop_status, dict)
                else stop_status
            )
            self.assertEqual(observed_stop, expected_stop)
            self.assertTrue(
                all(
                    call.kwargs["time_limit_seconds"] == 7.25
                    for call in solve_selection.call_args_list
                )
            )

    def test_zero_objective_is_complete_even_when_max_sets_is_zero(self):
        result = enumerate_optimal_vertex_sets(
            [],
            [],
            1,
            universe_mode="observed_alt",
            max_sets=0,
            time_limit_seconds=5,
        )
        self.assertTrue(result["complete"])
        self.assertEqual(result["objective"], 0)
        self.assertEqual(result["vertex_sets"], [(0,)])


class LikelihoodTests(unittest.TestCase):
    def test_same_vertex_set_is_edge_invariant(self):
        patterns = [("RR", 40), ("AR", 20), ("RA", 10), ("AA", 5), ("AX", 7)]
        fit_a = fit_vertex_mixture(patterns, [0, 1, 2, 3])
        fit_b = fit_vertex_mixture(patterns, [3, 1, 0, 2])
        weights_b = {v: w for v, w in zip(fit_b.vertices, fit_b.weights)}
        self.assertAlmostEqual(fit_a.log_likelihood, fit_b.log_likelihood, places=10)
        for vertex, weight in zip(fit_a.vertices, fit_a.weights):
            self.assertAlmostEqual(weight, weights_b[vertex], places=10)

    def test_split_merge_count_invariance(self):
        vertices = [0, 1, 2, 3]
        merged = fit_vertex_mixture([("RR", 20), ("AX", 10)], vertices)
        split = fit_vertex_mixture([("RR", 7), ("RR", 13), ("AX", 4), ("AX", 6)], vertices)
        self.assertAlmostEqual(merged.log_likelihood, split.log_likelihood, places=10)
        for left, right in zip(merged.weights, split.weights):
            self.assertAlmostEqual(left, right, places=10)

    def test_all_x_does_not_change_mle(self):
        vertices = [0, 1]
        base = fit_vertex_mixture([("R", 30), ("A", 10)], vertices)
        with_x = fit_vertex_mixture([("R", 30), ("A", 10), ("X", 100)], vertices)
        self.assertAlmostEqual(base.log_likelihood, with_x.log_likelihood, places=10)
        for left, right in zip(base.weights, with_x.weights):
            self.assertAlmostEqual(left, right, places=10)

    def test_em_is_monotone_and_simplex_valid(self):
        fit = fit_vertex_mixture([("RRR", 50), ("ARR", 20), ("AAX", 11), ("AAA", 4)], [0, 1, 3, 7])
        self.assertTrue(fit.converged)
        self.assertTrue(fit.monotone)
        self.assertAlmostEqual(sum(fit.weights), 1.0, places=12)
        self.assertTrue(all(weight >= 0 for weight in fit.weights))

    def test_slsqp_matches_em_log_likelihood(self):
        patterns = [("RRR", 50), ("ARR", 20), ("AAX", 11), ("AAA", 4)]
        vertices = [0, 1, 3, 7]
        em = fit_vertex_mixture(patterns, vertices)
        fast = fit_vertex_mixture_slsqp(patterns, vertices)
        self.assertTrue(fast.converged)
        self.assertTrue(fast.monotone)
        self.assertAlmostEqual(em.log_likelihood, fast.log_likelihood, places=7)
        self.assertAlmostEqual(sum(fast.weights), 1.0, places=12)
        self.assertLessEqual(fast.global_log_likelihood_gap_bound, 1e-8)

    def test_boundary_optimum_is_simplex_and_globally_certified(self):
        fit = fit_vertex_mixture_slsqp([("R", 100)], [0, 1])
        self.assertTrue(fit.converged)
        self.assertTrue(fit.monotone)
        self.assertAlmostEqual(sum(fit.weights), 1.0, places=12)
        self.assertLessEqual(fit.global_log_likelihood_gap_bound, 1e-8)
        self.assertLessEqual(fit.simplex_residual, 1e-12)
        self.assertLess(fit.weights[1], 1e-8)

    def test_rank_deficient_weights_are_not_claimed_identifiable(self):
        emission = np.asarray([[0.9, 0.9], [0.1, 0.1]], dtype=float)
        counts = np.asarray([9.0, 1.0], dtype=float)
        fit = fit_emission_mixture_certified(emission, counts, [0, 1])
        self.assertTrue(fit.converged)
        self.assertFalse(fit.mixture_weights_identifiable)
        self.assertEqual(fit.augmented_emission_rank, 1)
        self.assertLessEqual(fit.global_log_likelihood_gap_bound, 1e-8)

    def test_slsqp_status_failure_can_be_rescued_by_kkt_refinement(self):
        emission = np.asarray([[0.99, 0.01], [0.01, 0.99]], dtype=float)
        counts = np.asarray([80.0, 20.0], dtype=float)
        fake = SimpleNamespace(
            x=np.asarray([0.5, 0.5], dtype=float),
            success=False,
            status=8,
            nit=1,
        )
        with mock.patch("hypercube_exact.minimize", return_value=fake):
            fit = fit_emission_mixture_certified(emission, counts, [0, 1])
        self.assertFalse(fit.slsqp_success)
        self.assertEqual(fit.slsqp_status, 8)
        self.assertTrue(fit.converged)
        self.assertEqual(fit.optimizer_status, "CERTIFIED_PAIRWISE_REFINEMENT")
        self.assertLessEqual(fit.global_log_likelihood_gap_bound, 1e-8)
        self.assertGreater(fit.weights[0], fit.weights[1])

    def test_high_vertex_rank_deficient_case_is_deterministic(self):
        # 125 mixture states but only 8 observable pattern rows deliberately
        # stress the non-identifiable high-V regime.
        rng = np.random.default_rng(20260716)
        emission = rng.uniform(0.05, 1.0, size=(8, 125))
        counts = np.asarray([21, 13, 8, 5, 3, 2, 1, 1], dtype=float)
        first = fit_emission_mixture_certified(emission, counts, range(125))
        second = fit_emission_mixture_certified(emission, counts, range(125))
        self.assertTrue(first.converged)
        self.assertFalse(first.mixture_weights_identifiable)
        self.assertLessEqual(first.global_log_likelihood_gap_bound, 1e-8)
        self.assertAlmostEqual(first.log_likelihood, second.log_likelihood, places=12)
        np.testing.assert_allclose(first.weights, second.weights, atol=1e-12, rtol=0.0)

    def test_kkt_gap_is_a_global_log_likelihood_bound(self):
        emission = np.asarray([[0.9, 0.2], [0.1, 0.8]], dtype=float)
        counts = np.asarray([7.0, 3.0], dtype=float)
        arbitrary = np.asarray([0.5, 0.5], dtype=float)
        optimum = fit_emission_mixture_certified(emission, counts, [0, 1])
        certificate = mixture_kkt_certificate(emission, counts, arbitrary)
        arbitrary_ll = float(np.dot(counts, np.log(emission @ arbitrary)))
        self.assertLessEqual(
            optimum.log_likelihood - arbitrary_ll,
            certificate["global_log_likelihood_gap_bound"] + 1e-10,
        )


if __name__ == "__main__":
    unittest.main()
