#!/usr/bin/env python3
"""Tests for the isolated solver stress-panel exporter and runner."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import tempfile
import types
import unittest
from unittest import mock


REPO = pathlib.Path("/big7_disk/liaoyoyo2001/InterSubMod")
SCRIPT_DIR = (
    REPO / "research/20260718_solver_methyl_edge_probe/scripts"
)


def load(path: pathlib.Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


exporter = load(
    SCRIPT_DIR / "export_solver_stress_panel.py",
    "export_solver_stress_panel_test",
)
runner = load(
    SCRIPT_DIR / "run_solver_stress_panel.py",
    "run_solver_stress_panel_test",
)


def source_manifest():
    census = exporter.DEFAULT_CENSUS
    receipt = __import__("json").loads(census.read_text(encoding="utf-8"))
    return {
        "integrity": {"manifest_content_sha256": "unit-test-manifest"},
        "bindings": {
            "source_files": {
                "solver_probe": {
                    "path": str(exporter.DEFAULT_PROBE_SOURCE),
                },
                "optimized_backend": {
                    "path": str(exporter.DEFAULT_OPTIMIZED_SOURCE),
                },
                "frozen_current_solver": {
                    "path": receipt["source_files"]["hypercube_source"]["path"],
                },
            }
        },
    }


def case_for(full, partial, k, factorial=None):
    probe = load(exporter.DEFAULT_PROBE_SOURCE, "solver_probe_case_builder")
    instance = probe.build_instance(full, partial, k)
    payload = exporter.structural_payload(
        k=k,
        universe_mask=instance.universe_mask,
        full_patterns=full,
        partial_patterns=partial,
    )
    digest = exporter.semantic_sha256(payload)
    return {
        "case_id": f"test_{digest[:12]}",
        "structural_key_sha256": digest,
        "structural_input": payload,
        "factorial_oracle": factorial,
    }


class ExporterTests(unittest.TestCase):
    def test_structural_payload_is_order_invariant(self):
        left = exporter.structural_payload(
            k=3,
            universe_mask=7,
            full_patterns=["ARA", "AAR", "ARA"],
            partial_patterns=["AXA", "AAX"],
        )
        right = exporter.structural_payload(
            k=3,
            universe_mask=7,
            full_patterns=["AAR", "ARA"],
            partial_patterns=["AAX", "AXA"],
        )
        self.assertEqual(left, right)
        self.assertEqual(
            exporter.semantic_sha256(left),
            exporter.semantic_sha256(right),
        )

    def test_manifest_content_hash_excludes_only_its_own_value(self):
        document = {
            "x": 1,
            "integrity": {"manifest_content_sha256": "old"},
        }
        first = exporter.manifest_content_sha256(document)
        document["integrity"]["manifest_content_sha256"] = first
        self.assertEqual(exporter.manifest_content_sha256(document), first)
        document["x"] = 2
        self.assertNotEqual(exporter.manifest_content_sha256(document), first)

    def test_selected_record_block_audit_proves_nonterminal_contiguous_blocks(self):
        a = ("a",)
        b = ("b",)
        eof = ("eof",)
        audit = exporter.audit_selected_record_blocks(
            [a, a, b, b, eof],
            selected_keys={a, b},
            required_keys={a, b},
            stream_name="fixture",
        )
        self.assertTrue(audit["all_pass"])
        self.assertTrue(
            audit["checks"][
                "selected_blocks_end_strictly_before_truncation_boundary"
            ]
        )

    def test_selected_record_block_audit_rejects_reopened_or_terminal_block(self):
        a = ("a",)
        b = ("b",)
        audit = exporter.audit_selected_record_blocks(
            [a, b, a],
            selected_keys={a},
            required_keys={a},
            stream_name="fixture",
        )
        self.assertFalse(audit["all_pass"])
        self.assertFalse(
            audit["checks"]["selected_records_form_single_contiguous_block"]
        )
        self.assertFalse(
            audit["checks"]["selected_unit_is_not_terminal_eof_unit"]
        )


class RunnerTests(unittest.TestCase):
    def test_unbounded_current_and_optimized_match_factorial_family(self):
        manifest = source_manifest()
        case = case_for(
            ["AAAA"],
            [],
            4,
            {
                "kind": "single_terminal_permutation_family",
                "terminal_weight": 4,
                "expected_optimal_vertex_sets": 24,
            },
        )
        rows = {}
        for backend in runner.BACKENDS:
            rows[backend] = runner.run_backend_case(
                manifest,
                case,
                backend=backend,
                deadline_seconds=10.0,
                max_sets=None,
                q_max=8,
            )
            self.assertTrue(rows[backend]["family_complete"])
            self.assertTrue(rows[backend]["ranking_allowed"])
            self.assertFalse(rows[backend]["parent_trees_materialized"])
            self.assertEqual(rows[backend]["optimal_family_count"], 24)
            self.assertTrue(rows[backend]["factorial_oracle_pass"])
        self.assertIn("visited_states", rows["optimized"]["family_certificate_stats"])
        self.assertIn("table_cells", rows["optimized"]["objective_dp_stats"])
        self.assertEqual(
            rows["current"]["objective"],
            rows["optimized"]["objective"],
        )
        self.assertEqual(
            rows["current"]["optimal_family_digest"],
            rows["optimized"]["optimal_family_digest"],
        )

    def test_numeric_cap_is_fail_closed_for_both_backends(self):
        manifest = source_manifest()
        case = case_for(["AAAA"], [], 4)
        for backend in runner.BACKENDS:
            row = runner.run_backend_case(
                manifest,
                case,
                backend=backend,
                deadline_seconds=10.0,
                max_sets=1,
                q_max=8,
            )
            self.assertFalse(row["family_complete"])
            self.assertFalse(row["ranking_allowed"])
            self.assertFalse(row["incomplete_ranked"])
            self.assertIn("MAX_SETS", row["status"])

    def test_family_objective_is_authoritative_when_dp_routes_away(self):
        manifest = source_manifest()
        case = case_for(["AAAA"], [], 4)
        row = runner.run_backend_case(
            manifest,
            case,
            backend="optimized",
            deadline_seconds=10.0,
            max_sets=None,
            q_max=0,
        )
        self.assertFalse(row["dp_objective_certified"])
        self.assertIsNone(row["dp_objective"])
        self.assertTrue(row["family_complete"])
        self.assertTrue(row["objective_certified"])
        self.assertEqual(row["objective"], 3)
        self.assertTrue(row["ranking_allowed"])

    def test_acceptance_summary_requires_exact_fast_complete_pair(self):
        manifest = {
            "integrity": {"manifest_content_sha256": "manifest"},
            "cases": [
                {
                    "case_id": "control",
                    "selection_reasons": [
                        "COMPLETE_WITH_VERTEX_FAMILY_GE_100"
                    ],
                    "factorial_oracle": {
                        "expected_optimal_vertex_sets": 24,
                    },
                }
            ],
        }
        common = {
            "case_id": "control",
            "repeat": 1,
            "family_complete": True,
            "objective": 3,
            "optimal_family_digest": "same",
            "worker_completed": True,
            "incomplete_ranked": False,
            "factorial_oracle": {
                "expected_optimal_vertex_sets": 24,
            },
            "factorial_oracle_pass": True,
            "optimal_family_count": 24,
        }
        rows = [
            {
                **common,
                "backend": "current",
                "solver_elapsed_seconds": 10.0,
                "ru_maxrss_kib": 1000,
            },
            {
                **common,
                "backend": "optimized",
                "solver_elapsed_seconds": 1.0,
                "ru_maxrss_kib": 1200,
            },
        ]
        receipt = runner.build_summary(
            manifest,
            rows,
            selected_case_ids=["control"],
            requested_backends=list(runner.BACKENDS),
            deadline_seconds=30.0,
            max_sets=None,
        )
        self.assertTrue(receipt["all_pass"])
        self.assertEqual(
            receipt["performance"][
                "total_completion_time_speedup_current_over_optimized"
            ],
            10.0,
        )
        self.assertEqual(
            receipt["performance"]["reported_ratio_kind"],
            "BOTH_BACKENDS_COMPLETE_TOTAL_WALL_SPEEDUP",
        )
        self.assertEqual(
            receipt["exactness"][
                "control_objective_or_digest_mismatches"
            ],
            [],
        )

    def test_incomplete_current_ratio_is_labeled_conservative_lower_bound(self):
        manifest = {
            "integrity": {"manifest_content_sha256": "manifest"},
            "cases": [
                {
                    "case_id": "tail",
                    "selection_reasons": [
                        "33_INCOMPLETE_CANDIDATE_GENERATION"
                    ],
                    "factorial_oracle": {
                        "expected_optimal_vertex_sets": 24,
                    },
                }
            ],
        }
        rows = [
            {
                "case_id": "tail",
                "repeat": 1,
                "backend": "current",
                "family_complete": False,
                "status": "CANDIDATE_SET_INCOMPLETE_DEADLINE",
                "worker_completed": True,
                "incomplete_ranked": False,
                "solver_elapsed_seconds": 30.0,
                "subprocess_elapsed_seconds": 30.0,
                "factorial_oracle": {
                    "expected_optimal_vertex_sets": 24,
                },
            },
            {
                "case_id": "tail",
                "repeat": 1,
                "backend": "optimized",
                "family_complete": True,
                "status": "CANDIDATE_SET_COMPLETE",
                "worker_completed": True,
                "incomplete_ranked": False,
                "solver_elapsed_seconds": 1.0,
                "subprocess_elapsed_seconds": 1.0,
                "factorial_oracle": {
                    "expected_optimal_vertex_sets": 24,
                },
                "factorial_oracle_pass": True,
                "optimal_family_count": 24,
            },
        ]
        receipt = runner.build_summary(
            manifest,
            rows,
            selected_case_ids=["tail"],
            requested_backends=list(runner.BACKENDS),
            deadline_seconds=30.0,
            max_sets=None,
            expected_repeats=1,
        )
        performance = receipt["performance"]
        self.assertIsNone(
            performance[
                "total_completion_time_speedup_current_over_optimized"
            ]
        )
        self.assertEqual(
            performance[
                "conservative_lower_bound_speedup_current_over_optimized"
            ],
            30.0,
        )
        self.assertTrue(performance["current_deadline_censored"])
        self.assertEqual(
            performance["reported_ratio_kind"],
            "CONSERVATIVE_LOWER_BOUND_CURRENT_OVER_OPTIMIZED",
        )

    def test_factorial_gate_does_not_pass_empty_selection_vacuously(self):
        manifest = {
            "integrity": {"manifest_content_sha256": "manifest"},
            "cases": [
                {
                    "case_id": "plain",
                    "selection_reasons": [],
                    "factorial_oracle": None,
                }
            ],
        }
        receipt = runner.build_summary(
            manifest,
            [],
            selected_case_ids=["plain"],
            requested_backends=list(runner.BACKENDS),
            deadline_seconds=30.0,
            max_sets=None,
            expected_repeats=1,
        )
        self.assertFalse(
            receipt["checks"][
                "factorial_oracles_expected_count_positive_and_all_pass"
            ]
        )

    def test_repeated_schedule_interleaves_backend_order(self):
        schedule = runner.build_run_schedule(
            ["a"],
            list(runner.BACKENDS),
            repeats=2,
        )
        orders = {
            row["repeat"]: tuple(row["backend_order_for_case_repeat"])
            for row in schedule
        }
        self.assertEqual(orders[1], ("current", "optimized"))
        self.assertEqual(orders[2], ("optimized", "current"))

    def test_authority_pointer_binds_pointer_manifest_and_runner_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest_dir = root / "authoritative_r3"
            manifest_dir.mkdir()
            manifest_path = manifest_dir / "manifest.json"
            runner_sha = runner.sha256_file(runner.SCRIPT_PATH)
            document = {
                "schema_name": "intersubmod.solver_stress_panel_manifest",
                "checks": {"fixture": True},
                "authority": {"status": "AUTHORITATIVE_R3"},
                "bindings": {
                    "input_files": {},
                    "source_files": {
                        "stress_runner": {
                            "path": str(runner.SCRIPT_PATH),
                            "sha256": runner_sha,
                        }
                    },
                },
                "integrity": {},
            }
            document["integrity"]["manifest_content_sha256"] = (
                runner.manifest_content_sha256(document)
            )
            runner.write_json_with_sha256(manifest_path, document)
            pointer_path = root / runner.AUTHORITY_POINTER_NAME
            pointer = {
                "schema": (
                    "intersubmod.solver_stress_panel.authority_pointer.v1"
                ),
                "status": "AUTHORITATIVE_R3",
                "authoritative_manifest": "authoritative_r3/manifest.json",
                "authoritative_manifest_file_sha256": (
                    runner.sha256_file(manifest_path)
                ),
                "authoritative_manifest_content_sha256": document[
                    "integrity"
                ]["manifest_content_sha256"],
                "authoritative_runner_sha256": runner_sha,
            }
            runner.write_json_with_sha256(pointer_path, pointer)
            resolved, observed, authority = runner.resolve_authority_pointer(
                pointer_path,
                verify_inputs=False,
            )
            self.assertEqual(resolved, manifest_path.resolve())
            self.assertEqual(observed["authority"]["status"], "AUTHORITATIVE_R3")
            self.assertEqual(authority["runner_sha256"], runner_sha)

    def test_controller_interruption_writes_partial_marker_and_sidecar(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = pathlib.Path(tmp) / "run"
            output.mkdir()
            args = types.SimpleNamespace(
                output_dir=output,
                authority_pointer=None,
                manifest=pathlib.Path("/tmp/manifest.json"),
                authority_status="AUTHORITATIVE_R3",
                panel="primary",
                selection="all",
                backends=list(runner.BACKENDS),
                repeats=2,
                deadline=30.0,
                max_sets=None,
                q_max=8,
            )
            with mock.patch.object(
                runner,
                "_controller_main_impl",
                side_effect=KeyboardInterrupt(),
            ):
                with self.assertRaises(KeyboardInterrupt):
                    runner.controller_main(args)
            marker = output / "PARTIAL_RUN.json"
            self.assertTrue(marker.is_file())
            self.assertTrue(marker.with_name(marker.name + ".sha256").is_file())
            self.assertEqual(
                json.loads(marker.read_text())["reason"],
                "KeyboardInterrupt",
            )

    def test_max_sets_parser_supports_none_and_positive_integer(self):
        self.assertIsNone(runner.parse_max_sets("none"))
        self.assertIsNone(runner.parse_max_sets("unbounded"))
        self.assertEqual(runner.parse_max_sets("256"), 256)


if __name__ == "__main__":
    unittest.main()
