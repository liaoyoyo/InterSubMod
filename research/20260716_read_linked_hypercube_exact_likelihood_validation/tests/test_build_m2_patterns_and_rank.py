#!/usr/bin/env python3
import csv
import gzip
import json
import pathlib
import sys
import tempfile
import unittest
from collections import Counter
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_m2_patterns_and_rank import (  # noqa: E402
    Component,
    METHOD_CONTRACT,
    StructuralEnumerationCache,
    UnitEvidence,
    _SemanticListHasher,
    _pattern_rows,
    _quality_emission_matrix,
    conditional_candidate_ranking_bootstrap,
    build_evidence,
    fit_quality_aware_mixture,
    rank_unit,
    run,
    semantic_digest,
    sha256_path,
    summarize_runtime_values,
    summarize_units,
    topology_classes_for_vertex_set,
    validate_output_directory_contract,
)
from hypercube_exact import MixtureFit  # noqa: E402
from run_full_m2_ranking import (  # noqa: E402
    aggregate_primary_runtime_diagnostics,
    aggregate_rank_receipts,
)


def make_unit(pattern_counts, k, *, basis="hp1", family="1"):
    component = Component("D", "chr1", 3, basis, "C", 10, 10 + k - 1, k, tuple(range(k)))
    unit = UnitEvidence(component, family)
    for pattern, count in pattern_counts:
        codes = {index: symbol for index, symbol in enumerate(pattern) if symbol in {"R", "A"}}
        qualities = {index: 30 for index in codes}
        for _ in range(count):
            unit.add_projection(codes, qualities)
    return unit


class SymbolicDownstreamContractTests(unittest.TestCase):
    @staticmethod
    def _cache_kwargs(**overrides):
        values = {
            "full_patterns": ["AA"],
            "partial_patterns": [],
            "k": 2,
            "structural_alt_universe_mask": 0b11,
            "exact_k_max": 12,
            "max_vertex_sets": 32,
            "solver_time_limit_seconds": 5.0,
            "universe_mode": "observed_alt",
            "context": {
                "minread": 3,
                "threshold": 3,
                "unit_identity": ("D", "chr1", "PS_HP1", "100", "C", "1"),
            },
        }
        values.update(overrides)
        return values

    def test_structural_cache_complete_hit_is_deepcopy_isolated(self):
        solver = mock.Mock(wraps=__import__("hypercube_exact").enumerate_optimal_vertex_sets)
        cache = StructuralEnumerationCache(
            solver_source_sha256="solver-v1",
            enumerate_function=solver,
        )
        first, first_meta = cache.get_or_enumerate(**self._cache_kwargs())
        expected_n_sets = len(first["vertex_sets"])
        first["first"]["status"] = "MUTATED"
        first["vertex_sets"].append((999,))
        second, second_meta = cache.get_or_enumerate(**self._cache_kwargs())
        self.assertEqual(solver.call_count, 1)
        self.assertEqual(first_meta["structural_cache_outcome"], "MISS_STORED_COMPLETE")
        self.assertEqual(second_meta["structural_cache_outcome"], "HIT_COMPLETE")
        self.assertEqual(second["first"]["status"], "OPTIMAL")
        self.assertEqual(len(second["vertex_sets"]), expected_n_sets)
        self.assertEqual(cache.diagnostics()["lookups"], 2)
        self.assertEqual(cache.diagnostics()["hits"], 1)

    def test_structural_cache_never_negative_caches_incomplete(self):
        solver = mock.Mock(
            return_value={
                "first": {"status": "LIMIT_REACHED", "objective": None},
                "complete": False,
                "vertex_sets": [],
            }
        )
        cache = StructuralEnumerationCache(
            solver_source_sha256="solver-v1",
            enumerate_function=solver,
        )
        for _ in range(2):
            result, metadata = cache.get_or_enumerate(**self._cache_kwargs())
            self.assertFalse(result["complete"])
            self.assertEqual(
                metadata["structural_cache_outcome"],
                "MISS_INCOMPLETE_NOT_STORED",
            )
        diagnostics = cache.diagnostics()
        self.assertEqual(solver.call_count, 2)
        self.assertEqual(diagnostics["stores_complete"], 0)
        self.assertEqual(diagnostics["rejected_incomplete"], 2)
        self.assertEqual(diagnostics["entries_final"], 0)

    def test_structural_cache_digest_collision_cannot_create_false_hit(self):
        solver = mock.Mock(wraps=__import__("hypercube_exact").enumerate_optimal_vertex_sets)
        cache = StructuralEnumerationCache(
            solver_source_sha256="solver-v1",
            diagnostic_digest=lambda _payload: "forced-collision",
            enumerate_function=solver,
        )
        first, _ = cache.get_or_enumerate(**self._cache_kwargs())
        second, _ = cache.get_or_enumerate(
            **self._cache_kwargs(
                full_patterns=["A"],
                k=1,
                structural_alt_universe_mask=0b1,
            )
        )
        self.assertNotEqual(first["vertex_set_ids"], second["vertex_set_ids"])
        self.assertEqual(solver.call_count, 2)
        self.assertEqual(cache.diagnostics()["digest_collision_guard_events"], 1)

    def test_structural_cache_solver_contract_parameters_are_keyed(self):
        solver = mock.Mock(wraps=__import__("hypercube_exact").enumerate_optimal_vertex_sets)
        cache = StructuralEnumerationCache(
            solver_source_sha256="solver-v1",
            enumerate_function=solver,
        )
        cache.get_or_enumerate(**self._cache_kwargs())
        cache.get_or_enumerate(
            **self._cache_kwargs(max_vertex_sets=33)
        )
        cache.get_or_enumerate(
            **self._cache_kwargs(solver_time_limit_seconds=5.5)
        )
        cache.get_or_enumerate(
            **self._cache_kwargs(full_patterns=["AR", "AA"])
        )
        self.assertEqual(solver.call_count, 4)
        self.assertEqual(cache.diagnostics()["hits"], 0)
        self.assertEqual(cache.diagnostics()["misses"], 4)

    def test_same_structure_different_counts_reuses_only_enumeration(self):
        cache = StructuralEnumerationCache(solver_source_sha256="solver-v1")
        first = rank_unit(
            make_unit([("AA", 3)], 2),
            fixed_error_grid=(),
            solver_time_limit_seconds=5,
            structural_enumeration_cache=cache,
        )
        second = rank_unit(
            make_unit([("AA", 6)], 2),
            fixed_error_grid=(),
            solver_time_limit_seconds=5,
            structural_enumeration_cache=cache,
        )
        diagnostics = cache.diagnostics()
        self.assertEqual(diagnostics["solver_invocations"], 1)
        self.assertEqual(diagnostics["hits"], 1)
        self.assertNotEqual(first["best_log_likelihood"], second["best_log_likelihood"])
        self.assertNotEqual(first["unit_semantic_sha256"], second["unit_semantic_sha256"])

    def test_monotonic_runtime_segments_are_separate_and_exactly_summarized(self):
        # V=1 and no fixed-error grid gives exactly six instrumentation clock reads:
        # unit start; candidate start/end; primary likelihood start/end; unit end.
        with mock.patch(
            "build_m2_patterns_and_rank._monotonic_ns",
            side_effect=(0, 1_000_000_000, 3_000_000_000, 4_000_000_000,
                         7_000_000_000, 10_000_000_000),
        ):
            row = rank_unit(
                make_unit([("A", 3)], 1),
                fixed_error_grid=(),
                solver_time_limit_seconds=5,
            )
        timing = row["_runtime_diagnostics"]
        self.assertEqual(timing["candidate_generation_elapsed_seconds"], 2.0)
        self.assertEqual(timing["likelihood_fit_elapsed_seconds"], 3.0)
        self.assertEqual(timing["unit_total_elapsed_seconds"], 10.0)
        with mock.patch(
            "build_m2_patterns_and_rank._monotonic_ns",
            side_effect=(0, 2_000_000_000, 5_000_000_000, 8_000_000_000,
                         12_000_000_000, 20_000_000_000),
        ):
            differently_timed = rank_unit(
                make_unit([("A", 3)], 1),
                fixed_error_grid=(),
                solver_time_limit_seconds=5,
            )
        self.assertNotEqual(timing, differently_timed["_runtime_diagnostics"])
        self.assertEqual(row["unit_semantic_sha256"], differently_timed["unit_semantic_sha256"])
        summary = summarize_runtime_values([4.0, 1.0, 5.0, 2.0, 3.0])
        self.assertEqual(
            summary,
            {"n": 5, "sum": 15.0, "max": 5.0, "p50": 3.0, "p95": 5.0, "p99": 5.0},
        )

    def test_streaming_semantic_list_hash_matches_materialized_list_hash(self):
        rows = [
            {"text": "中文", "number": 1, "flag": True},
            {"text": "second", "number": 2.5, "missing": None},
        ]
        streaming = _SemanticListHasher()
        for row in rows:
            streaming.update(row)
        self.assertEqual(streaming.n_rows, len(rows))
        self.assertEqual(streaming.hexdigest(), semantic_digest(rows))

    def test_rank_unit_detail_sinks_preserve_rows_and_public_result(self):
        retained = rank_unit(make_unit([("AA", 3)], 2), solver_time_limit_seconds=5)
        streamed_candidates = []
        streamed_responsibilities = []
        streamed = rank_unit(
            make_unit([("AA", 3)], 2),
            solver_time_limit_seconds=5,
            candidate_record_sink=streamed_candidates.append,
            responsibility_record_sink=streamed_responsibilities.append,
            retain_detail_records=False,
        )
        self.assertEqual(
            {key: value for key, value in retained.items() if not key.startswith("_")},
            {key: value for key, value in streamed.items() if not key.startswith("_")},
        )
        self.assertEqual(streamed_candidates, retained["_candidate_records"])
        self.assertEqual(streamed_responsibilities, retained["_responsibility_records"])
        self.assertEqual(streamed["_candidate_records"], [])
        self.assertEqual(streamed["_responsibility_records"], [])

    def test_rax_is_one_symbolic_group_with_two_conceptual_completions(self):
        unit = make_unit([("RRA", 3), ("AAA", 3), ("RAX", 3)], 3)
        rows = _pattern_rows(unit, minread=3)
        rax = [row for row in rows if row["pattern_rax"] == "RAX"]
        self.assertEqual(sum(row["n_molecules"] for row in rax), 3)
        self.assertEqual({row["n_conceptual_completions"] for row in rax}, {2})
        result = rank_unit(unit, solver_time_limit_seconds=5)
        self.assertEqual(result["minimum_hidden_nodes"], 1)
        self.assertEqual(result["raw_tree_candidates_T"], 1)
        self.assertEqual(result["distinct_vertex_sets_V"], 1)
        self.assertEqual(result["selection_status"], "T1_CANDIDATE_UNIQUE")
        self.assertEqual(result["partial_group_coverage_denominator"], 1)
        self.assertEqual(result["partial_groups_covered"], 1)
        self.assertEqual(result["partial_groups_unsatisfied"], 0)

    def test_effective_partial_completions_can_be_less_than_two_power_u(self):
        unit = make_unit([("RAX", 3)], 3)
        row = _pattern_rows(unit, minread=3)[0]
        self.assertEqual(row["n_conceptual_completions"], 2)
        self.assertEqual(row["n_effective_free_in_structural_alt_universe"], 0)
        self.assertEqual(row["n_effective_completions"], 1)

    def test_all_x_is_excluded_but_conserved(self):
        unit = make_unit([("XXX", 5)], 3)
        result = rank_unit(unit)
        self.assertEqual(result["molecule_component_projections"], 5)
        self.assertEqual(result["informative_scoring_molecules"], 0)
        self.assertEqual(result["all_x_excluded_molecules"], 5)
        self.assertEqual(result["selection_status"], "NO_INFORMATIVE_SCORING_MOLECULE")

    def test_same_vertex_set_different_edges_are_not_scored_separately(self):
        # {RR, AR, RA, AA}: AA has two possible parents, but there is one vertex set.
        unit = make_unit([("AR", 3), ("RA", 3), ("AA", 3)], 2)
        result = rank_unit(unit, solver_time_limit_seconds=5)
        self.assertEqual(result["raw_tree_candidates_T"], 2)
        self.assertEqual(result["distinct_vertex_sets_V"], 1)
        self.assertEqual(result["top_edge_variants"], 2)
        self.assertEqual(result["selection_status"], "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED")

    def test_multiple_vertex_sets_are_fitted_before_convergence_is_inspected(self):
        # AA alone permits the two minimum paths RR->AR->AA and RR->RA->AA.
        # Both candidates must be fitted and retained as a likelihood tie; this
        # exercises the real V>1 branch before all_converged is assigned.
        unit = make_unit([("AA", 3)], 2)
        result = rank_unit(unit, solver_time_limit_seconds=5)
        self.assertEqual(result["raw_tree_candidates_T"], 2)
        self.assertEqual(result["distinct_vertex_sets_V"], 2)
        self.assertEqual(result["best_vertex_sets"], 2)
        self.assertEqual(result["selection_status"], "LIKELIHOOD_TIED_VERTEX_SETS")
        self.assertTrue(result["all_fits_converged"])
        self.assertTrue(result["all_fits_monotone"])

    def test_primary_nonconvergence_skips_bootstrap(self):
        unit = make_unit([("AA", 3)], 2)

        def nonconverged_fit(_patterns, vertices, **_kwargs):
            vertices = tuple(vertices)
            return MixtureFit(
                vertices,
                tuple(1.0 / len(vertices) for _ in vertices),
                0.0,
                False,
                0,
                1,
                True,
            )

        with mock.patch(
            "build_m2_patterns_and_rank.fit_quality_aware_mixture",
            side_effect=nonconverged_fit,
        ), mock.patch(
            "build_m2_patterns_and_rank.conditional_candidate_ranking_bootstrap"
        ) as bootstrap:
            result = rank_unit(
                unit,
                conditional_candidate_ranking_bootstrap_replicates=3,
                solver_time_limit_seconds=5,
            )
        self.assertEqual(result["selection_status"], "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE")
        self.assertEqual(
            result["conditional_candidate_ranking_bootstrap_status"],
            "NOT_RUN_PRIMARY_NONCONVERGENCE",
        )
        self.assertEqual(result["conditional_candidate_ranking_bootstrap_replicates"], 0)
        bootstrap.assert_not_called()

    def test_minread_applies_only_to_structure_not_scoring(self):
        unit = make_unit([("AR", 3), ("RA", 2)], 2)
        result = rank_unit(unit, minread=3, solver_time_limit_seconds=5)
        self.assertEqual(result["informative_scoring_molecules"], 5)
        self.assertEqual(result["structural_retained_molecules"], 3)
        self.assertEqual(result["below_minread_scoring_molecules"], 2)
        self.assertEqual(result["n_structural_pattern_groups"], 1)

    def test_k_greater_than_twelve_is_explicit_abstain(self):
        unit = make_unit([("A" * 13, 3)], 13)
        result = rank_unit(unit, exact_k_max=12)
        self.assertEqual(result["selection_status"], "LOCAL_ONLY_OR_ABSTAIN_K_GT_EXACT_LIMIT")
        self.assertFalse(result["candidate_generation_complete"])
        self.assertIsNone(result["raw_tree_candidates_T"])

    def test_raw_k_over_limit_but_effective_alt_k_small_is_solved(self):
        unit = make_unit([("A" + "R" * 12, 3)], 13)
        result = rank_unit(unit, exact_k_max=1, solver_time_limit_seconds=5)
        self.assertTrue(result["candidate_generation_complete"])
        self.assertEqual(result["k_component_sites"], 13)
        self.assertEqual(result["k_observed_alt_active"], 1)
        self.assertEqual(result["k_scoring_alt_observed"], 1)
        self.assertEqual(result["n_not_structural_alt_active_sites"], 12)

    def test_exact_patterns_are_not_compatibility_pooled_for_minread(self):
        unit = make_unit([("AXR", 1), ("ARX", 1), ("AXX", 1)], 3)
        result = rank_unit(unit, minread=3)
        self.assertEqual(result["selection_status"], "STRUCTURE_ABSTAIN_NO_MINREAD_ALT_PATTERN")
        self.assertEqual(result["structural_retained_molecules"], 0)

    def test_r_only_dimensions_do_not_change_candidate_ranking(self):
        compact = make_unit([("AA", 3)], 2)
        expanded = make_unit([("AARR", 3)], 4)
        compact_result = rank_unit(compact, solver_time_limit_seconds=5)
        expanded_result = rank_unit(expanded, solver_time_limit_seconds=5)
        self.assertEqual(compact_result["selection_status"], expanded_result["selection_status"])
        self.assertEqual(compact_result["distinct_vertex_sets_V"], expanded_result["distinct_vertex_sets_V"])
        self.assertEqual(compact_result["k_observed_alt_active"], expanded_result["k_observed_alt_active"])

    def test_quality_emission_uses_per_call_phred_and_marginalizes_x(self):
        matrix, counts = _quality_emission_matrix(
            [("AX", (40, -1), 2)],
            [0, 1, 3],
            minimum_error_rate=1e-6,
            maximum_error_rate=0.25,
        )
        self.assertEqual(tuple(counts), (2.0,))
        denominator = 1.0 - 2.0 * 1e-4 / 3.0
        self.assertAlmostEqual(matrix[0, 0], (1e-4 / 3.0) / denominator, places=10)
        self.assertAlmostEqual(matrix[0, 1], (1.0 - 1e-4) / denominator, places=10)
        # The unknown second site must not distinguish vertices 1 and 3.
        self.assertAlmostEqual(matrix[0, 1], matrix[0, 2], places=14)

    def test_quality_mixture_has_global_kkt_certificate(self):
        fit = fit_quality_aware_mixture(
            [("R", (35,), 80), ("A", (35,), 20)],
            [0, 1],
        )
        self.assertTrue(fit.converged)
        self.assertTrue(fit.monotone)
        self.assertLessEqual(fit.global_log_likelihood_gap_bound, 1e-8)
        self.assertLessEqual(fit.simplex_residual, 1e-12)

    def test_projection_and_cell_mass_conservation(self):
        unit = make_unit([("RAX", 3), ("XXX", 2)], 3)
        result = rank_unit(unit, minread=3, solver_time_limit_seconds=5)
        self.assertEqual(
            result["molecule_component_projections"],
            result["informative_scoring_molecules"] + result["all_x_excluded_molecules"],
        )
        self.assertEqual(
            result["projected_cells"],
            result["fixed_ra_cells"] + result["non_ra_cells_marginalized"],
        )
        self.assertEqual(result["bq_scoring_molecules"], result["informative_scoring_molecules"])

    def test_topology_classes_are_analytical_over_parent_choices(self):
        self.assertEqual(topology_classes_for_vertex_set((0, 1)), ("single-only",))
        self.assertEqual(topology_classes_for_vertex_set((0, 1, 2)), ("sister-only",))
        self.assertEqual(topology_classes_for_vertex_set((0, 1, 3)), ("direct-only",))
        self.assertEqual(topology_classes_for_vertex_set((0, 1, 2, 3)), ("sister+direct",))

    def test_seeded_molecule_bootstrap_is_deterministic_and_separate_from_error_grid(self):
        quality_counts = [
            ("RRA", (30, 30, 30), 5),
            ("RAX", (30, 30, -1), 5),
            ("RAA", (30, 30, 30), 2),
        ]
        kwargs = {
            "k": 3,
            "primary_top_vertices": [(0, 4, 6)],
            "replicates": 3,
            "seed": 17,
            "minimum_error_rate": 1e-6,
            "maximum_error_rate": 0.25,
            "tie_tolerance": 1e-6,
        }
        first = conditional_candidate_ranking_bootstrap(
            quality_counts, [(0, 4, 6), (0, 2, 4)], **kwargs
        )
        second = conditional_candidate_ranking_bootstrap(
            quality_counts, [(0, 4, 6), (0, 2, 4)], **kwargs
        )
        self.assertEqual(first, second)
        self.assertEqual(first["replicates"], 3)
        self.assertEqual(first["seed"], 17)

    def test_nonconverged_candidates_do_not_leak_into_topology_aggregate(self):
        unit = make_unit([("AA", 3)], 2)

        def nonconverged_fit(_patterns, vertices, **_kwargs):
            vertices = tuple(vertices)
            return MixtureFit(vertices, tuple(1 / len(vertices) for _ in vertices), 0.0, False, 0, 1, True)

        with mock.patch(
            "build_m2_patterns_and_rank.fit_quality_aware_mixture", side_effect=nonconverged_fit
        ):
            row = rank_unit(unit, solver_time_limit_seconds=5)
        summary = summarize_units([row])
        self.assertEqual(summary["topology_evaluated_units"], 0)
        self.assertEqual(summary["topology_class_inclusion_counts"], {})
        self.assertEqual(summary["coarse_topology_unique_class_counts"], {})


class FileContractTests(unittest.TestCase):
    def test_preflight_created_ranking_directory_contract_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            missing = root / "missing"
            with self.assertRaisesRegex(FileExistsError, "unavailable"):
                validate_output_directory_contract(missing, require_existing_empty=True)

            empty = root / "empty"
            empty.mkdir()
            validate_output_directory_contract(empty, require_existing_empty=True)
            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                validate_output_directory_contract(empty)

            (empty / "prior-output").write_text("x", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "not empty"):
                validate_output_directory_contract(empty, require_existing_empty=True)

            target = root / "target"
            target.mkdir()
            symlink = root / "symlink"
            symlink.symlink_to(target, target_is_directory=True)
            with self.assertRaisesRegex(FileExistsError, "not a real directory"):
                validate_output_directory_contract(symlink, require_existing_empty=True)

    @staticmethod
    def _write_gzip_tsv(path, fields, rows):
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)

    def _fixture(self, root):
        calls = root / "D.chr1.molecule_sparse_calls.tsv.gz"
        sites = root / "D.chr1.site_catalog.tsv.gz"
        components = root / "D.chr1.components.tsv.gz"
        membership = root / "D.chr1.site_component_membership.tsv.gz"
        self._write_gzip_tsv(
            sites,
            ("site_index", "chrom", "pos1", "ref", "alt"),
            [
                {"site_index": index, "chrom": "chr1", "pos1": 10 + index, "ref": "C", "alt": "T"}
                for index in range(3)
            ],
        )
        self._write_gzip_tsv(
            components,
            (
                "dataset", "chrom", "component_basis", "phase_set", "phase_set_status",
                "inference_role", "threshold", "component_id", "start1", "end1", "k",
            ),
            [{
                "dataset": "D", "chrom": "chr1", "component_basis": "PS_HP1", "phase_set": "100",
                "phase_set_status": "KNOWN_PS_PRIMARY", "inference_role": "PRIMARY_PS_AWARE", "threshold": 3,
                "component_id": "C", "start1": 10, "end1": 12, "k": 3,
            }],
        )
        self._write_gzip_tsv(
            membership,
            (
                "dataset", "chrom", "component_basis", "phase_set", "phase_set_status",
                "inference_role", "threshold", "site_index", "pos1", "component_id",
            ),
            [
                {
                    "dataset": "D", "chrom": "chr1", "component_basis": "PS_HP1", "phase_set": "100",
                    "phase_set_status": "KNOWN_PS_PRIMARY", "inference_role": "PRIMARY_PS_AWARE", "threshold": 3,
                    "site_index": index, "pos1": 10 + index, "component_id": "C",
                }
                for index in range(3)
            ],
        )
        call_rows = []
        for pattern in ("RRA", "AAA"):
            for repeat in range(3):
                call_rows.append({
                    "molecule_id": f"{pattern}{repeat}", "hp_family": "1", "phase_set": "100", "site_indices": "0,1,2",
                    "call_codes": pattern, "base_qualities": "30,30,30",
                })
        for repeat in range(3):
            call_rows.append({
                "molecule_id": f"RAX{repeat}", "hp_family": "1", "phase_set": "100", "site_indices": "0,1",
                "call_codes": "RA", "base_qualities": "30,30",
            })
        call_rows.append({
            "molecule_id": "missing-reasons", "hp_family": "1", "phase_set": "100", "site_indices": "0,1,2",
            "call_codes": "ODL", "base_qualities": ",,",
        })
        self._write_gzip_tsv(
            calls,
            ("molecule_id", "hp_family", "phase_set", "site_indices", "call_codes", "base_qualities"),
            call_rows,
        )
        outputs = {
            path.name: {"size_bytes": path.stat().st_size, "sha256": sha256_path(path)}
            for path in (calls, sites, components, membership)
        }
        receipt_path = root / "receipt.json"
        receipt_path.write_text(
            json.dumps({
                "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
                "schema_version": "1.2.0",
                "scope": {"dataset": "D", "chrom": "chr1"},
                "parameters": {"component_linkage_bases": ["PS_HP1", "PS_HP2"]},
                "all_pass": True,
                "outputs": outputs,
                "receipt_integrity": {
                    "scheme": "external_sha256_sidecar_v1",
                    "sidecar_name": "receipt.json.sha256",
                },
            }),
            encoding="utf-8",
        )
        (root / "receipt.json.sha256").write_text(
            f"{sha256_path(receipt_path)}  receipt.json\n", encoding="ascii"
        )

    def test_end_to_end_receipt_has_funnels_quality_policy_and_hashes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            self._fixture(input_dir)
            receipt = run(
                input_dir,
                output_dir,
                families=("1",),
                structural_exact_pattern_minreads=(1, 2, 3, 5),
                primary_structural_exact_pattern_minread=3,
                solver_time_limit_seconds=5,
            )
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(receipt["parameters"]["method_contract"], METHOD_CONTRACT)
            self.assertTrue(
                receipt["parameters"]["method_contract"]
                ["active_compatible_vertex_indices_materialized_for_sparse_rows"]
            )
            self.assertFalse(
                receipt["parameters"]["method_contract"]
                ["completion_wise_tree_worlds_materialized"]
            )
            self.assertNotIn(
                "no_" + "partial_completions_materialized", receipt["checks"]
            )
            self.assertEqual(
                receipt["scope"]["component_basis_mode"],
                "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
            )
            self.assertEqual(receipt["aggregate"]["molecule_component_projections"], 10)
            self.assertEqual(receipt["aggregate"]["informative_scoring_molecules"], 9)
            self.assertEqual(receipt["aggregate"]["all_x_excluded_molecules"], 1)
            self.assertEqual(receipt["aggregate"]["bq_scoring_molecules"], 9)
            optimizer = receipt["candidate_evidence_counts"]["optimizer_diagnostics"]
            self.assertEqual(
                optimizer["candidate_fits"],
                optimizer["globally_kkt_certified_candidate_fits"],
            )
            self.assertLessEqual(optimizer["max_global_log_likelihood_gap_bound"], 1e-8)
            self.assertTrue(
                receipt["checks"][
                    "all_converged_candidate_fits_have_global_kkt_certificate"
                ]
            )
            self.assertEqual(receipt["aggregate"]["projected_call_class_counts"]["O"], 1)
            self.assertEqual(receipt["aggregate"]["projected_call_class_counts"]["D"], 1)
            self.assertEqual(receipt["aggregate"]["projected_call_class_counts"]["L"], 1)
            self.assertEqual(len(receipt["outputs"]["m2_component_hp_rank_units.tsv.gz"]["semantic_sha256"]), 64)
            self.assertIn("m2_compressed_vertex_set_candidates.tsv.gz", receipt["outputs"])
            self.assertIn("m2_winning_pattern_state_responsibilities.tsv.gz", receipt["outputs"])
            self.assertIn("m2_unit_runtime_diagnostics.tsv.gz", receipt["outputs"])
            runtime = receipt["runtime_diagnostics"]
            self.assertEqual(runtime["clock"], "time.monotonic_ns")
            self.assertEqual(runtime["unit"], "seconds")
            structural_cache = runtime["structural_enumeration_cache"]
            self.assertEqual(structural_cache["lookups"], 3)
            self.assertEqual(structural_cache["misses"], 1)
            self.assertEqual(structural_cache["hits"], 2)
            self.assertEqual(structural_cache["solver_invocations"], 1)
            self.assertEqual(structural_cache["stores_complete"], 1)
            self.assertEqual(structural_cache["rejected_incomplete"], 0)
            self.assertEqual(structural_cache["same_unit_cross_minread_hits"], 2)
            self.assertEqual(
                runtime["scopes"]["primary_unit_evaluations"]["n_unit_evaluations"], 1
            )
            self.assertEqual(
                runtime["scopes"]["all_structural_minread_unit_evaluations"]
                ["n_unit_evaluations"],
                4,
            )
            for scope in runtime["scopes"].values():
                for metric in (
                    "candidate_generation_elapsed_seconds",
                    "likelihood_fit_elapsed_seconds",
                    "unit_total_elapsed_seconds",
                ):
                    self.assertEqual(scope[metric]["n"], scope["n_unit_evaluations"])
                    self.assertGreaterEqual(scope[metric]["p95"], 0.0)
            with gzip.open(
                output_dir / "m2_unit_runtime_diagnostics.tsv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                timing_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(timing_rows), 4)
            self.assertEqual(sum(row["structural_minread_role"] == "PRIMARY" for row in timing_rows), 1)
            timing_aggregate = aggregate_primary_runtime_diagnostics([{
                "dataset": "D", "chrom": "chr1", "status": "PASS", "receipt": receipt,
            }])
            self.assertEqual(timing_aggregate["n_child_runtime_files"], 1)
            self.assertEqual(timing_aggregate["n_unit_evaluations"], 1)
            self.assertEqual(receipt["partial_pattern_funnel"]["partial_groups_unsatisfied"], 0)
            self.assertEqual(receipt["schema_version"], "2.0.0")
            self.assertEqual(
                set(receipt["sensitivity_by_structural_exact_pattern_minread"]),
                {"1", "2", "3", "5"},
            )
            for output_name, count_key in (
                ("m2_compressed_vertex_set_candidates.tsv.gz", "compressed_vertex_set_candidate_rows"),
                ("m2_winning_pattern_state_responsibilities.tsv.gz", "posterior_responsibility_rows"),
            ):
                output_path = output_dir / output_name
                with gzip.open(output_path, "rt", encoding="utf-8") as handle:
                    self.assertEqual(
                        sum(1 for _row in csv.DictReader(handle, delimiter="\t")),
                        receipt["candidate_evidence_counts"][count_key],
                    )
                self.assertEqual(receipt["outputs"][output_name]["sha256"], sha256_path(output_path))
            self.assertEqual(list(root.glob(".output.staging.*")), [])

    def test_cache_enabled_and_disabled_scientific_outputs_are_semantically_equal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            input_dir = root / "input"
            input_dir.mkdir()
            self._fixture(input_dir)
            cached = run(
                input_dir,
                root / "cached",
                families=("1",),
                structural_exact_pattern_minreads=(1, 2, 3, 5),
                primary_structural_exact_pattern_minread=3,
                solver_time_limit_seconds=5,
                enable_structural_enumeration_cache=True,
            )
            uncached = run(
                input_dir,
                root / "uncached",
                families=("1",),
                structural_exact_pattern_minreads=(1, 2, 3, 5),
                primary_structural_exact_pattern_minread=3,
                solver_time_limit_seconds=5,
                enable_structural_enumeration_cache=False,
            )
            self.assertTrue(cached["all_pass"])
            self.assertTrue(uncached["all_pass"])
            for output_name in (
                "m2_symbolic_pattern_counts.tsv.gz",
                "m2_component_hp_rank_units.tsv.gz",
                "m2_compressed_vertex_set_candidates.tsv.gz",
                "m2_winning_pattern_state_responsibilities.tsv.gz",
            ):
                self.assertEqual(
                    cached["outputs"][output_name]["semantic_sha256"],
                    uncached["outputs"][output_name]["semantic_sha256"],
                    output_name,
                )
            self.assertEqual(cached["aggregate"], uncached["aggregate"])
            self.assertGreater(
                cached["runtime_diagnostics"]["structural_enumeration_cache"]["hits"],
                0,
            )
            self.assertFalse(
                uncached["runtime_diagnostics"]["structural_enumeration_cache"]["enabled"]
            )

    def test_ranking_preserves_preflight_output_directory_inode(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            inode_before = output_dir.stat().st_ino
            self._fixture(input_dir)
            receipt = run(
                input_dir,
                output_dir,
                families=("1",),
                structural_exact_pattern_minreads=(1, 2, 3, 5),
                primary_structural_exact_pattern_minread=3,
                solver_time_limit_seconds=5,
                require_existing_empty_output_dir=True,
            )
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(output_dir.stat().st_ino, inode_before)
            self.assertTrue((output_dir / "receipt.json").is_file())
            self.assertEqual(list(root.glob(".output.staging.*")), [])

    def test_hp_basis_only_creates_matching_family_units(self):
        with tempfile.TemporaryDirectory() as temporary:
            calls = pathlib.Path(temporary) / "calls.tsv.gz"
            self._write_gzip_tsv(
                calls,
                ("molecule_id", "hp_family", "phase_set", "site_indices", "call_codes", "base_qualities"),
                [
                    {"molecule_id": "m1", "hp_family": "1", "phase_set": "100", "site_indices": "0", "call_codes": "A", "base_qualities": "30"},
                    {"molecule_id": "m2", "hp_family": "2", "phase_set": "200", "site_indices": "0", "call_codes": "A", "base_qualities": "30"},
                ],
            )
            hp1 = Component("D", "chr1", 3, "PS_HP1", "C1", 10, 10, 1, (0,), "100", "KNOWN_PS_PRIMARY", "PRIMARY_PS_AWARE")
            hp2 = Component("D", "chr1", 3, "PS_HP2", "C2", 10, 10, 1, (0,), "200", "KNOWN_PS_PRIMARY", "PRIMARY_PS_AWARE")
            units, counts = build_evidence(
                calls,
                {("PS_HP1", "100", 3, "C1"): hp1, ("PS_HP2", "200", 3, "C2"): hp2},
                {("PS_HP1", "100", 3): {0: "C1"}, ("PS_HP2", "200", 3): {0: "C2"}},
                ("1", "2"),
                {0},
            )
            self.assertEqual(
                set(units),
                {("PS_HP1", "100", 3, "C1", "1"), ("PS_HP2", "200", 3, "C2", "2")},
            )
            self.assertEqual(counts["molecule_component_family_projections"], 2)

    def test_primary_evidence_is_strictly_family_and_phase_set_aware(self):
        with tempfile.TemporaryDirectory() as temporary:
            calls = pathlib.Path(temporary) / "calls.tsv.gz"
            self._write_gzip_tsv(
                calls,
                ("molecule_id", "hp_family", "phase_set", "site_indices", "call_codes", "base_qualities"),
                [
                    {"molecule_id": "ps100", "hp_family": "1", "phase_set": "100", "site_indices": "0", "call_codes": "A", "base_qualities": "30"},
                    {"molecule_id": "ps200", "hp_family": "1", "phase_set": "200", "site_indices": "0", "call_codes": "R", "base_qualities": "30"},
                    {"molecule_id": "missing", "hp_family": "1", "phase_set": "", "site_indices": "0", "call_codes": "A", "base_qualities": "30"},
                    {"molecule_id": "wrong-family", "hp_family": "2", "phase_set": "100", "site_indices": "0", "call_codes": "A", "base_qualities": "30"},
                ],
            )
            c100 = Component("D", "chr1", 3, "PS_HP1", "C100", 10, 10, 1, (0,), "100", "KNOWN_PS_PRIMARY", "PRIMARY_PS_AWARE")
            c200 = Component("D", "chr1", 3, "PS_HP1", "C200", 10, 10, 1, (0,), "200", "KNOWN_PS_PRIMARY", "PRIMARY_PS_AWARE")
            units, counts = build_evidence(
                calls,
                {("PS_HP1", "100", 3, "C100"): c100, ("PS_HP1", "200", 3, "C200"): c200},
                {("PS_HP1", "100", 3): {0: "C100"}, ("PS_HP1", "200", 3): {0: "C200"}},
                ("1", "2"),
                {0},
            )
            self.assertEqual(units[("PS_HP1", "100", 3, "C100", "1")].pattern_counts, Counter({"A": 1}))
            self.assertEqual(units[("PS_HP1", "200", 3, "C200", "1")].pattern_counts, Counter({"R": 1}))
            self.assertEqual(counts["molecule_component_family_projections"], 2)
            self.assertEqual(counts["sparse_molecule_rows_excluded_by_component_or_phase_set_contract"], 2)

    def test_evidence_route_index_visits_only_matching_exact_ps(self):
        with tempfile.TemporaryDirectory() as temporary:
            calls = pathlib.Path(temporary) / "calls.tsv.gz"
            self._write_gzip_tsv(
                calls,
                ("molecule_id", "hp_family", "phase_set", "site_indices", "call_codes", "base_qualities"),
                [{
                    "molecule_id": "match-42", "hp_family": "1", "phase_set": "42",
                    "site_indices": "0", "call_codes": "A", "base_qualities": "30",
                }],
            )
            components = {}
            mappings = {}
            for index in range(100):
                ps = str(index)
                component_id = f"C{index}"
                component = Component(
                    "D", "chr1", 3, "PS_HP1", component_id, 10, 10, 1, (0,),
                    ps, "KNOWN_PS_PRIMARY", "PRIMARY_PS_AWARE",
                )
                components[("PS_HP1", ps, 3, component_id)] = component
                mappings[("PS_HP1", ps, 3)] = {0: component_id}
            units, counts = build_evidence(calls, components, mappings, ("1",), {0})
            self.assertEqual(
                units[("PS_HP1", "42", 3, "C42", "1")].pattern_counts,
                Counter({"A": 1}),
            )
            self.assertTrue(all(
                not unit.pattern_counts
                for key, unit in units.items()
                if key != ("PS_HP1", "42", 3, "C42", "1")
            ))
            self.assertEqual(counts["evidence_component_mapping_keys_total"], 100)
            self.assertEqual(counts["evidence_mapping_routes_visited"], 1)
            self.assertEqual(counts["evidence_mapping_keys_naive_scan_reference"], 100)
            self.assertEqual(counts["molecule_component_family_projections"], 1)

    def test_full_aggregator_stratifies_dataset_basis_and_threshold(self):
        summary = {
            "n_component_hp_units": 2,
            "n_components": 2,
            "molecule_component_projections": 10,
            "informative_scoring_molecules": 9,
            "all_x_excluded_molecules": 1,
            "selection_status_counts": {"T1_CANDIDATE_UNIQUE": 2},
            "candidate_generation_status_counts": {"EXACT_OPTIMAL_VERTEX_SETS_COMPLETE": 2},
            "k_route_counts": {"EXACT_K_LE12": 2},
            "projected_call_class_counts": {"R_or_A": 12, "O": 1},
            "bootstrap_status_counts": {"COMPLETE": 1, "NOT_APPLICABLE_V1": 1},
            "topology_class_inclusion_counts": {"direct-only": 2},
            "topology_derivation_status_counts": {
                "ANALYTICAL_COMPLETE_OVER_ALL_WINNING_VERTEX_SETS_AND_PARENT_CHOICES": 2
            },
        }
        results = [
            {
                "dataset": dataset,
                "chrom": "chr1",
                "status": "REUSED_PASS",
                "receipt": {
                    "input_counts": {
                        "sparse_molecule_rows_total": 10,
                        "selected_sparse_call_code_counts": {"R": 6, "A": 6, "O": 1},
                    },
                    "aggregate_by_linkage_basis_threshold": {"HP1": {"3": summary}},
                },
            }
            for dataset in ("D1", "D2")
        ]
        aggregate = aggregate_rank_receipts(results)
        group = aggregate["by_linkage_basis_threshold"]["HP1"]["3"]
        self.assertEqual(group["n_component_hp_units"], 4)
        self.assertEqual(group["informative_scoring_molecules"], 18)
        self.assertEqual(group["topology_class_inclusion_counts"]["direct-only"], 4)
        self.assertEqual(set(aggregate["by_dataset"]), {"D1", "D2"})

    def test_summary_quality_primary_partition_is_conservative(self):
        rows = [
            rank_unit(make_unit([("AR", 3)], 2), solver_time_limit_seconds=5),
            rank_unit(
                make_unit([("AR", 3), ("RA", 3), ("AA", 3)], 2),
                solver_time_limit_seconds=5,
            ),
            rank_unit(make_unit([("AA", 3)], 2), solver_time_limit_seconds=5),
            rank_unit(make_unit([("XXX", 3)], 3), solver_time_limit_seconds=5),
        ]
        summary = summarize_units(rows)
        self.assertEqual(summary["quality_primary_unique_vertex_units"], 2)
        self.assertEqual(summary["quality_primary_tied_vertex_units"], 1)
        self.assertEqual(summary["rank_abstain_units"], 1)
        self.assertEqual(
            summary["quality_primary_unique_vertex_units"]
            + summary["quality_primary_tied_vertex_units"]
            + summary["rank_abstain_units"],
            summary["n_component_hp_units"],
        )


if __name__ == "__main__":
    unittest.main()
