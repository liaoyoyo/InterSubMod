#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


TOPIC = Path(__file__).resolve().parents[1]
SCRIPT_PATH = TOPIC / "scripts" / "build_final_numeric_summary.py"
FIXTURE_PATH = Path(__file__).resolve().parent / "fixtures" / "final_numeric_summary_bundle.json"
SPEC = importlib.util.spec_from_file_location("build_final_numeric_summary", SCRIPT_PATH)
assert SPEC and SPEC.loader
SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def write_json_receipt(path: Path, payload: dict) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    digest = sha256_path(path)
    path.with_name(f"{path.name}.sha256").write_text(
        f"{digest}  {path.name}\n", encoding="ascii"
    )
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": digest}


def file_identity(path: Path) -> dict:
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha256_path(path)}


def partial_funnel() -> dict:
    return {
        "definitions": {
            "conceptual": "2^u conceptual completions are not independent tree worlds",
            "effective": "one sparse compatible-state hit row per retained exact pattern",
            "universe_source": "synthetic fixture",
        },
        "unique_rax_pattern_groups": {
            "denominator": 1,
            "full": 0,
            "partial": 1,
            "u_number_of_X_distribution": {"1": 1},
            "conceptual_completions_2_pow_u_distribution": {"2": 1},
            "conceptual_completions_weighted_total": 2,
            "observed_alt_effective_completions_distribution": {"2": 1},
            "observed_alt_effective_completions_weighted_total": 2,
            "observed_alt_effective_zero_due_to_fixed_alt_outside_structural_universe": 0,
        },
    }


def rank_cell(outcome: str) -> dict:
    tied = outcome == "tied"
    raw_t = 4 if tied else 1
    vertex_sets = 2 if tied else 1
    values = {field: 0 for field in SUMMARY.SUM_FIELDS}
    values.update(
        {
            "n_component_hp_units": 1,
            "n_components": 1,
            "molecule_component_projections": 4,
            "informative_scoring_molecules": 3,
            "all_x_excluded_molecules": 1,
            "structural_retained_molecules": 2,
            "below_minread_scoring_molecules": 1,
            "bq_scoring_molecules": 3,
            "raw_tree_candidates_T_complete_units": raw_t,
            "distinct_vertex_sets_V_complete_units": vertex_sets,
            "solver_complete_units": 1,
            "quality_primary_unique_vertex_units": 0 if tied else 1,
            "quality_primary_tied_vertex_units": 1 if tied else 0,
            "fixed_error_grid_stable_units": 1,
            "fixed_error_grid_evaluated_units": 1,
            "coarse_topology_class_unique_units": 1,
            "parent_edge_assignment_unique_units": 0 if tied else 1,
            "exact_topology_proven_unique_units": 0 if tied else 1,
            "topology_evaluated_units": 1,
            "topology_class_inclusion_counts_denominator": 1,
            "k_component_sites_total": 2,
            "k_observed_alt_active_total": 2,
            "k_scoring_alt_observed_total": 2,
            "structural_partial_pattern_groups": 1,
            "partial_group_coverage_denominator": 1,
            "partial_groups_covered": 1,
        }
    )
    counters = {field: {} for field in SUMMARY.COUNTER_FIELDS}
    counters.update(
        {
            "selection_status_counts": {"LIKELIHOOD_TIED" if tied else "LIKELIHOOD_UNIQUE": 1},
            "candidate_generation_status_counts": {"EXACT_OPTIMAL_VERTEX_SETS_COMPLETE": 1},
            "k_route_counts": {"EXACT_K_LE12": 1},
            "projected_call_class_counts": {"R_or_A": 4},
            "coarse_topology_unique_class_counts": {"direct-only": 1},
            "topology_class_inclusion_counts": {"direct-only": 1},
            "topology_derivation_status_counts": {"ANALYTICAL_COMPLETE": 1},
            "exact_topology_uniqueness_status_counts": {
                "EDGE_NONIDENTIFIABLE" if tied else "PROVEN_UNIQUE": 1
            },
        }
    )
    return {**values, **counters, "partial_pattern_funnel": partial_funnel()}


def candidate_row(
    *, unit: str, family: str, candidate_id: str, vertex_id: str, roles: dict,
    winner_status: str, topology: list[str], parent_choice_count: int = 1,
) -> dict[str, str]:
    states = {
        key: "".join("A" if int(key) & (1 << bit) else "R" for bit in range(2))
        for key in roles
    }
    return {
        "unit_key": unit,
        "dataset": "SYNTH",
        "chrom": "chr1",
        "component_id": f"C{family}",
        "threshold": "3",
        "hp_family": family,
        "ps": "10",
        "candidate_id": candidate_id,
        "vertex_set_id": vertex_id,
        "vertex_states": json.dumps(states, sort_keys=True, separators=(",", ":")),
        "vertex_roles": json.dumps(roles, sort_keys=True, separators=(",", ":")),
        "parent_choice_count": str(parent_choice_count),
        "profile_log_likelihood": "-1",
        "relative_log_likelihood": "0",
        "mixture_weights_pi": json.dumps({key: 1.0 / len(roles) for key in roles}, separators=(",", ":")),
        "winner_status": winner_status,
        "tie_group": "TOP_TIE" if winner_status == "TIED_WINNER" else "",
        "coarse_topology_class": json.dumps(topology, separators=(",", ":")),
        "candidate_set_complete": "true",
    }


def build_bundle(root: Path) -> dict[str, Path]:
    fixture = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
    dataset = fixture["dataset"]
    chrom = fixture["chromosome"]
    threshold = str(fixture["bridge_threshold"])
    extraction_root = root / "extraction"
    ranking_root = root / "ranking"

    component = fixture["component_cell"]
    extraction_child = {
        "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
        "schema_version": "1.2.0",
        "scope": {"dataset": dataset, "chrom": chrom},
        "counts": fixture["extraction_counts"],
        "component_summary_by_linkage_basis": {
            basis: {threshold: component} for basis in SUMMARY.PRIMARY_BASES
        },
        "checks": {"synthetic_fixture_conserves": True},
        "all_pass": True,
    }
    extraction_child_path = extraction_root / "samples" / dataset / chrom / "receipt.json"
    extraction_child_identity = write_json_receipt(extraction_child_path, extraction_child)
    component_aggregate = {
        basis: {threshold: component} for basis in SUMMARY.PRIMARY_BASES
    }
    extraction_aggregate = {
        "counts": fixture["extraction_counts"],
        "component_summary_by_linkage_basis": component_aggregate,
        "component_summary_by_threshold": {},
        "phase_set_contract_totals": {},
        "legacy_cross_phase_set_aggregation_audit": {},
        "by_dataset": {
            dataset: {
                "task_status_counts": {"PASS": 1},
                "counts": fixture["extraction_counts"],
            }
        },
    }
    extraction_full = {
        "schema_name": "intersubmod.m2_full_extraction_receipt",
        "schema_version": "1.2.0",
        "scope": {"datasets": [dataset], "chromosomes": [chrom], "expected_tasks": 1},
        "aggregate": extraction_aggregate,
        "results": [{"dataset": dataset, "chrom": chrom, "status": "PASS", "receipt": extraction_child}],
        "checks": {"synthetic_full_extraction": True},
        "all_pass": True,
    }
    extraction_full_path = extraction_root / "full_extraction_receipt.json"
    extraction_full_identity = write_json_receipt(extraction_full_path, extraction_full)

    runtime_path = ranking_root / "samples" / dataset / chrom / "m2_unit_runtime_diagnostics.tsv.gz"
    runtime_path.parent.mkdir(parents=True, exist_ok=True)
    runtime_rows = []
    for index, item in enumerate(fixture["runtime_rows"], start=1):
        runtime_rows.append(
            {
                "dataset": dataset,
                "chrom": chrom,
                "threshold": threshold,
                "component_basis": item["basis"],
                "phase_set": "10",
                "component_id": f"C{item['family']}",
                "family": item["family"],
                "structural_exact_pattern_minread": str(
                    fixture["primary_structural_exact_pattern_minread"]
                ),
                "structural_minread_role": "PRIMARY",
                "candidate_generation_invoked": "true",
                "likelihood_fit_invoked": "true",
                **{
                    metric: str(item[metric]) for metric in SUMMARY.RUNTIME_METRICS
                },
            }
        )
    with gzip.open(runtime_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY.RUNTIME_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(runtime_rows)
    runtime_identity = file_identity(runtime_path)
    runtime_values = {
        metric: [float(row[metric]) for row in runtime_rows] for metric in SUMMARY.RUNTIME_METRICS
    }
    runtime_metrics = {
        metric: SUMMARY.runtime_summary(values) for metric, values in runtime_values.items()
    }
    runtime_invoked = {
        metric: runtime_metrics[metric]
        for metric in (
            "candidate_generation_elapsed_seconds",
            "likelihood_fit_elapsed_seconds",
        )
    }
    child_runtime = {
        "schema_name": "intersubmod.m2_unit_runtime_diagnostics",
        "schema_version": "1.0.0",
        "clock": "time.monotonic_ns",
        "unit": "seconds",
        "per_unit_output": runtime_path.name,
        "scopes": {
            "primary_unit_evaluations": {
                "n_unit_evaluations": 2,
                **runtime_metrics,
            },
            "all_structural_minread_unit_evaluations": {
                "n_unit_evaluations": 2,
                **runtime_metrics,
            },
        },
        "primary_invoked_segment_scopes": runtime_invoked,
    }
    hp1 = rank_cell("unique")
    hp2 = rank_cell("tied")
    input_counts = {
        "sparse_molecule_rows_total": 4,
        "sparse_molecule_rows_known_ps": 4,
        "sparse_molecule_rows_missing_ps": 0,
        "sparse_molecule_rows_included_in_at_least_one_selected_unit": 4,
        "selected_sparse_call_code_counts": {"A": 4, "R": 4},
        "hp_family_rows": {"1": 2, "2": 2},
    }
    ranking_child = {
        "schema_name": "intersubmod.m2_symbolic_patterns_vertex_rank_receipt",
        "schema_version": "2.0.0",
        "scope": {"dataset": [dataset], "chrom": [chrom]},
        "input_counts": input_counts,
        "aggregate_by_linkage_basis_threshold": {
            "PS_HP1": {threshold: hp1},
            "PS_HP2": {threshold: hp2},
        },
        "partial_pattern_funnel_by_linkage_basis_threshold": {
            "PS_HP1": {threshold: partial_funnel()},
            "PS_HP2": {threshold: partial_funnel()},
        },
        "runtime_diagnostics": child_runtime,
        "outputs": {runtime_path.name: runtime_identity},
        "checks": {"synthetic_rank_child": True},
        "all_pass": True,
    }
    ranking_child_path = ranking_root / "samples" / dataset / chrom / "receipt.json"
    ranking_child_identity = write_json_receipt(ranking_child_path, ranking_child)

    hp1_unit = "SYNTH|chr01|PS_HP1|PS=10|B003|C1|HP1|M3"
    hp2_unit = "SYNTH|chr01|PS_HP2|PS=10|B003|C2|HP2|M3"
    candidate_rows = [
        candidate_row(
            unit=hp1_unit,
            family="1",
            candidate_id="C000001",
            vertex_id="V1",
            roles={"0": ["root"], "3": ["full-observed"]},
            winner_status="UNIQUE_WINNER",
            topology=["direct-only"],
        ),
        candidate_row(
            unit=hp2_unit,
            family="2",
            candidate_id="C000001",
            vertex_id="V2A",
            roles={"0": ["root"], "1": ["partial-compatible"], "3": ["full-observed"]},
            winner_status="TIED_WINNER",
            topology=["direct-only"],
            parent_choice_count=2,
        ),
        candidate_row(
            unit=hp2_unit,
            family="2",
            candidate_id="C000002",
            vertex_id="V2B",
            roles={"0": ["root"], "2": ["connector"], "3": ["full-observed"]},
            winner_status="TIED_WINNER",
            topology=["direct-only"],
            parent_choice_count=2,
        ),
    ]
    candidate_path = ranking_root / "m2_ps_aware_candidate_table.tsv.gz"
    semantic = hashlib.sha256()
    with gzip.open(candidate_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY.CANDIDATE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in candidate_rows:
            writer.writerow(row)
            semantic.update(
                json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
                + b"\n"
            )
    candidate_identity = file_identity(candidate_path)
    candidate_metadata = {
        "schema_version": "2.0.0",
        "format": "tsv.gz",
        "columns": list(SUMMARY.CANDIDATE_COLUMNS),
        "sort_order": "unit_key,candidate_id",
        **candidate_identity,
        "semantic_sha256": semantic.hexdigest(),
        "n_rows": 3,
        "n_units": 2,
    }
    input_flat = {key: value for key, value in input_counts.items() if isinstance(value, int)}
    ranking_aggregate = {
        "input_call_funnel": input_flat,
        "input_sparse_call_code_counts": input_counts["selected_sparse_call_code_counts"],
        "by_linkage_basis_threshold": {
            "PS_HP1": {threshold: hp1},
            "PS_HP2": {threshold: hp2},
        },
        "by_dataset": {
            dataset: {
                "input_call_funnel": input_flat,
                "input_sparse_call_code_counts": input_counts["selected_sparse_call_code_counts"],
                "input_hp_family_rows": input_counts["hp_family_rows"],
                "by_linkage_basis_threshold": {
                    "PS_HP1": {threshold: hp1},
                    "PS_HP2": {threshold: hp2},
                },
            }
        },
    }
    full_runtime = {
        "schema_name": "intersubmod.m2_full_primary_runtime_diagnostics",
        "schema_version": "1.0.0",
        "clock": "time.monotonic_ns",
        "unit": "seconds",
        "n_child_runtime_files": 1,
        "n_unit_evaluations": 2,
        "metrics": runtime_metrics,
        "metrics_when_invoked": runtime_invoked,
    }
    ranking_full = {
        "schema_name": "intersubmod.m2_full_ranking_receipt",
        "schema_version": "2.0.0",
        "scope": {"datasets": [dataset], "chromosomes": [chrom], "expected_tasks": 1},
        "run_contract": {
            "thresholds": [3],
            "parameters": {"primary_structural_exact_pattern_minread": 3},
        },
        "aggregate": ranking_aggregate,
        "candidate_table": candidate_metadata,
        "runtime_diagnostics": full_runtime,
        "results": [
            {
                "dataset": dataset,
                "chrom": chrom,
                "status": "PASS",
                "rank_receipt": {
                    "path": ranking_child_identity["path"],
                    "sha256": ranking_child_identity["sha256"],
                },
            }
        ],
        "checks": {"synthetic_full_ranking": True},
        "all_pass": True,
    }
    ranking_full_path = ranking_root / "full_ranking_receipt.json"
    ranking_full_identity = write_json_receipt(ranking_full_path, ranking_full)

    final = {
        "schema_name": "intersubmod.m2_full_independent_verification",
        "schema_version": "1.0.0",
        "scope": {"datasets": [dataset], "chromosomes": [chrom], "expected_tasks": 1},
        "release_binding": {
            "validation_evidence_eligible": True,
            "deep_release_verification": {"all_pass": True},
        },
        "extraction": {
            "receipt_path": extraction_full_identity["path"],
            "receipt_sha256": extraction_full_identity["sha256"],
            "task_index": [
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "child_receipt_sha256": extraction_child_identity["sha256"],
                }
            ],
            "recomputed_aggregate": extraction_aggregate,
        },
        "ranking": {
            "receipt_path": ranking_full_identity["path"],
            "receipt_sha256": ranking_full_identity["sha256"],
            "task_index": [
                {
                    "dataset": dataset,
                    "chrom": chrom,
                    "rank_receipt_sha256": ranking_child_identity["sha256"],
                }
            ],
            "recomputed": {
                "input_call_funnel": input_flat,
                "input_sparse_call_code_counts": input_counts["selected_sparse_call_code_counts"],
                "by_linkage_basis_threshold": {
                    "PS_HP1": {threshold: hp1},
                    "PS_HP2": {threshold: hp2},
                },
            },
            "candidate_table": {
                "path": candidate_identity["path"],
                "size_bytes": candidate_identity["size_bytes"],
                "sha256": candidate_identity["sha256"],
                "semantic_sha256": semantic.hexdigest(),
                "n_rows": 3,
                "n_units": 2,
                "all_rows_match_independent_child_reconstruction": True,
                "winner_partitions_conserved": True,
            },
            "runtime_diagnostics": {
                "n_child_runtime_files": 1,
                "n_unit_evaluations": 2,
                "metrics": runtime_metrics,
                "metrics_when_invoked": runtime_invoked,
                "all_child_and_full_runtime_summaries_independently_recomputed": True,
            },
        },
        "checks": {"synthetic_independent_verifier": True},
        "all_pass": True,
    }
    final_path = root / "final_verification.json"
    write_json_receipt(final_path, final)
    return {
        "extraction_root": extraction_root,
        "ranking_root": ranking_root,
        "final": final_path,
        "candidate": candidate_path,
        "extraction_child": extraction_child_path,
        "ranking_child": ranking_child_path,
        "runtime": runtime_path,
        "terminal_extraction": extraction_full_path,
        "terminal_ranking": ranking_full_path,
    }


class FinalNumericSummaryTest(unittest.TestCase):
    def test_authenticated_bundle_derives_component_hstar_and_tied_topology(self):
        expected = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))["expected"]
        with tempfile.TemporaryDirectory() as tmp:
            paths = build_bundle(Path(tmp))
            result = SUMMARY.build_summary(
                paths["extraction_root"],
                paths["ranking_root"],
                paths["final"],
                require_canonical_scope=False,
            )
        component = result["extraction"]["by_dataset"]["SYNTH"][
            "component_by_linkage_basis_threshold"
        ]["PS_HP1"]["3"]
        self.assertEqual(component["k_greater_than_1"], 1)
        hp1 = result["ranking"]["by_dataset"]["SYNTH"][
            "by_HP_basis_and_bridge_threshold"
        ]["PS_HP1"]["3"]
        hp2 = result["ranking"]["by_dataset"]["SYNTH"][
            "by_HP_basis_and_bridge_threshold"
        ]["PS_HP2"]["3"]
        self.assertEqual(hp1["ranking_outcome"]["unique_first"], expected["HP1_unique"])
        self.assertEqual(hp2["ranking_outcome"]["tied_first"], expected["HP2_tied"])
        self.assertEqual(
            hp2["candidate_structure"]["candidate_table"]["tied_by_coarse_topology"]["consistent"],
            expected["HP2_tied_topology_consistent"],
        )
        self.assertEqual(
            hp2["candidate_structure"]["candidate_table"]["n_parent_edge_trees_T"],
            4,
        )
        self.assertEqual(
            hp1["candidate_structure"]["candidate_table"]["topology"][
                "exact_topology_proven_unique_units"
            ],
            1,
        )
        self.assertEqual(
            result["ranking"]["overall_candidate_table"]["h_star_distribution"],
            expected["h_star_distribution"],
        )
        partition = result["ranking"]["overall_candidate_table"]["tree_vertex_partition"]
        self.assertEqual(partition["counts"], expected["tree_vertex_partition"]["counts"])
        self.assertEqual(
            partition["denominator"], expected["tree_vertex_partition"]["denominator"]
        )
        self.assertEqual(
            partition["shares"]["T_EQ_1_V_EQ_1"]["value"], 0.5
        )
        self.assertEqual(
            hp1["candidate_structure"]["candidate_table"]["tree_vertex_partition"][
                "counts"
            ],
            {
                "T_EQ_1_V_EQ_1": 1,
                "T_GT_1_V_EQ_1": 0,
                "T_GT_1_V_GT_1": 0,
            },
        )
        overall_hp2 = result["ranking"]["overall_by_HP_basis_and_bridge_threshold"][
            "PS_HP2"
        ]["3"]
        self.assertEqual(
            overall_hp2["candidate_structure"]["candidate_table"][
                "tree_vertex_partition"
            ]["counts"],
            {
                "T_EQ_1_V_EQ_1": 0,
                "T_GT_1_V_EQ_1": 0,
                "T_GT_1_V_GT_1": 1,
            },
        )
        self.assertIsNone(
            result["unsupported_or_nonidentifiable"][
                "exact_parent_edge_topology_for_tied_vertex_sets"
            ]["value"]
        )
        self.assertTrue(result["all_pass"])

    def test_candidate_byte_tamper_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = build_bundle(Path(tmp))
            with paths["candidate"].open("ab") as handle:
                handle.write(b"tamper")
            with self.assertRaisesRegex(SUMMARY.SummaryError, "candidate table: recorded SHA-256 mismatch"):
                SUMMARY.build_summary(
                    paths["extraction_root"], paths["ranking_root"], paths["final"],
                    require_canonical_scope=False,
                )

    def test_extraction_child_byte_tamper_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = build_bundle(Path(tmp))
            with paths["extraction_child"].open("ab") as handle:
                handle.write(b"tamper")
            with self.assertRaisesRegex(SUMMARY.SummaryError, "extraction child .* sidecar"):
                SUMMARY.build_summary(
                    paths["extraction_root"], paths["ranking_root"], paths["final"],
                    require_canonical_scope=False,
                )

    def test_ranking_child_byte_tamper_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = build_bundle(Path(tmp))
            with paths["ranking_child"].open("ab") as handle:
                handle.write(b"tamper")
            with self.assertRaisesRegex(SUMMARY.SummaryError, "ranking child .* sidecar"):
                SUMMARY.build_summary(
                    paths["extraction_root"], paths["ranking_root"], paths["final"],
                    require_canonical_scope=False,
                )

    def test_runtime_byte_tamper_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = build_bundle(Path(tmp))
            with paths["runtime"].open("ab") as handle:
                handle.write(b"tamper")
            with self.assertRaisesRegex(SUMMARY.SummaryError, "runtime .* recorded SHA-256 mismatch"):
                SUMMARY.build_summary(
                    paths["extraction_root"], paths["ranking_root"], paths["final"],
                    require_canonical_scope=False,
                )

    def test_tied_winners_with_different_coarse_classes_are_inconsistent(self):
        unit = "SYNTH|chr01|PS_HP2|PS=10|B003|C2|HP2|M3"
        rows = [
            candidate_row(
                unit=unit,
                family="2",
                candidate_id="C000001",
                vertex_id="A",
                roles={"0": ["root"], "1": ["connector"]},
                winner_status="TIED_WINNER",
                topology=["direct-only"],
            ),
            candidate_row(
                unit=unit,
                family="2",
                candidate_id="C000002",
                vertex_id="B",
                roles={"0": ["root"], "2": ["connector"]},
                winner_status="TIED_WINNER",
                topology=["sister-only"],
            ),
        ]
        strata = {}
        overall = SUMMARY.new_candidate_accumulator()
        SUMMARY.finish_candidate_unit(rows, ("SYNTH",), ("chr1",), strata, overall)
        frozen = SUMMARY.freeze_candidate_accumulator(overall)
        self.assertEqual(frozen["tied_by_coarse_topology"]["consistent"], 0)
        self.assertEqual(frozen["tied_by_coarse_topology"]["inconsistent"], 1)
        self.assertEqual(frozen["topology"]["coarse_class_multiple_units"], 1)
        self.assertEqual(
            frozen["topology"]["coarse_ambiguous_class_set_counts"],
            {"sister-only|direct-only": 1},
        )

    def test_parent_choice_count_must_be_positive(self):
        unit = "SYNTH|chr01|PS_HP1|PS=10|B003|C1|HP1|M3"
        row = candidate_row(
            unit=unit,
            family="1",
            candidate_id="C000001",
            vertex_id="V1",
            roles={"0": ["root"], "3": ["full-observed"]},
            winner_status="UNIQUE_WINNER",
            topology=["direct-only"],
            parent_choice_count=0,
        )
        with self.assertRaisesRegex(SUMMARY.SummaryError, "must be positive"):
            SUMMARY.finish_candidate_unit(
                [row], ("SYNTH",), ("chr1",), {}, SUMMARY.new_candidate_accumulator()
            )

    def test_tree_vertex_partition_covers_all_three_mutually_exclusive_states(self):
        strata = {}
        overall = SUMMARY.new_candidate_accumulator()
        rows_by_unit = (
            [
                candidate_row(
                    unit="SYNTH|chr01|PS_HP1|PS=10|B003|C1|HP1|M3",
                    family="1",
                    candidate_id="C000001",
                    vertex_id="V1",
                    roles={"0": ["root"], "3": ["full-observed"]},
                    winner_status="UNIQUE_WINNER",
                    topology=["direct-only"],
                    parent_choice_count=1,
                )
            ],
            [
                candidate_row(
                    unit="SYNTH|chr01|PS_HP1|PS=10|B003|C2|HP1|M3",
                    family="1",
                    candidate_id="C000001",
                    vertex_id="V2",
                    roles={"0": ["root"], "3": ["full-observed"]},
                    winner_status="UNIQUE_WINNER",
                    topology=["direct-only"],
                    parent_choice_count=2,
                )
            ],
            [
                candidate_row(
                    unit="SYNTH|chr01|PS_HP2|PS=10|B003|C3|HP2|M3",
                    family="2",
                    candidate_id="C000001",
                    vertex_id="V3A",
                    roles={"0": ["root"], "1": ["connector"]},
                    winner_status="TIED_WINNER",
                    topology=["direct-only"],
                    parent_choice_count=1,
                ),
                candidate_row(
                    unit="SYNTH|chr01|PS_HP2|PS=10|B003|C3|HP2|M3",
                    family="2",
                    candidate_id="C000002",
                    vertex_id="V3B",
                    roles={"0": ["root"], "2": ["connector"]},
                    winner_status="TIED_WINNER",
                    topology=["direct-only"],
                    parent_choice_count=1,
                ),
            ],
        )
        for rows in rows_by_unit:
            SUMMARY.finish_candidate_unit(
                rows, ("SYNTH",), ("chr1",), strata, overall
            )
        frozen = SUMMARY.freeze_candidate_accumulator(overall)
        partition = frozen["tree_vertex_partition"]
        self.assertEqual(
            partition["counts"],
            {
                "T_EQ_1_V_EQ_1": 1,
                "T_GT_1_V_EQ_1": 1,
                "T_GT_1_V_GT_1": 1,
            },
        )
        self.assertEqual(partition["denominator"], 3)
        self.assertEqual(
            {key: value["value"] for key, value in partition["shares"].items()},
            {key: 1 / 3 for key in SUMMARY.TREE_VERTEX_BUCKETS},
        )
        self.assertEqual(frozen["n_candidate_vertex_sets_V"], 4)
        self.assertEqual(frozen["n_parent_edge_trees_T"], 5)

    def test_tree_vertex_partition_tampering_fails_closed(self):
        target = SUMMARY.new_candidate_accumulator()
        target["n_units"] = 1
        target["n_candidate_vertex_sets_V"] = 1
        target["n_parent_edge_trees_T"] = 1
        target["unique_first"] = 1
        target["topology_evaluated_units"] = 1
        target["coarse_topology_class_unique_units"] = 1
        target["parent_edge_assignment_unique_units"] = 1
        target["exact_topology_proven_unique_units"] = 1
        target["h_star_distribution"][0] = 1
        target["tree_vertex_partition_counts"]["T_EQ_1_V_EQ_1"] = 1
        valid = SUMMARY.freeze_candidate_accumulator(target)["tree_vertex_partition"]

        mutations = {
            "partition": lambda payload: payload["counts"].update({"UNKNOWN": 0}),
            "denominator": lambda payload: payload.update({"denominator": 2}),
            "share": lambda payload: payload["shares"]["T_EQ_1_V_EQ_1"].update(
                {"numerator": 99}
            ),
            "nonconservation": lambda payload: payload["counts"].update(
                {"T_GT_1_V_EQ_1": 1}
            ),
        }
        expected_errors = {
            "partition": "count bucket keys mismatch",
            "denominator": "denominator mismatch",
            "share": "share mismatch",
            "nonconservation": "do not conserve denominator",
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                tampered = json.loads(json.dumps(valid))
                mutate(tampered)
                with self.assertRaisesRegex(SUMMARY.SummaryError, expected_errors[name]):
                    SUMMARY.verify_tree_vertex_partition(tampered, 1, "tampered")

    def test_candidate_aggregate_T_less_than_V_fails_closed(self):
        target = SUMMARY.new_candidate_accumulator()
        target["n_units"] = 1
        target["n_candidate_vertex_sets_V"] = 2
        target["n_parent_edge_trees_T"] = 1
        with self.assertRaisesRegex(SUMMARY.SummaryError, "T>=V"):
            SUMMARY.freeze_candidate_accumulator(target)

    def test_h_star_role_vocabulary_fails_closed(self):
        row = candidate_row(
            unit="SYNTH|chr01|PS_HP1|PS=10|B003|C1|HP1|M3",
            family="1",
            candidate_id="C000001",
            vertex_id="V1",
            roles={"0": ["root"], "1": ["mystery-role"]},
            winner_status="UNIQUE_WINNER",
            topology=["direct-only"],
        )
        with self.assertRaisesRegex(SUMMARY.SummaryError, "unknown vertex role"):
            SUMMARY.minimum_extra_states(row, "bad-role")

    def test_component_distribution_conservation_rejects_mismatch(self):
        bad = {
            "n_components": 2,
            "n_singletons_k1": 0,
            "n_multisite_k_gt1": 1,
            "max_k_component_sites": 2,
            "k_component_sites_distribution": {"2": 1},
        }
        with self.assertRaisesRegex(SUMMARY.SummaryError, "k distribution sum mismatch"):
            SUMMARY.verify_component_cell(bad, "bad")

    def test_formal_mode_rejects_noncanonical_scope_before_presentation(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = build_bundle(Path(tmp))
            with self.assertRaisesRegex(SUMMARY.SummaryError, "seven-dataset scope"):
                SUMMARY.build_summary(
                    paths["extraction_root"], paths["ranking_root"], paths["final"],
                    require_canonical_scope=True,
                )

    def test_output_is_exclusive_and_sidecar_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "summary.json"
            identity = SUMMARY.write_json_and_sidecar_exclusive(
                output,
                {"schema_name": "x", "checks": {"x": True}, "all_pass": True},
            )
            self.assertEqual(identity["sha256"], sha256_path(output))
            self.assertEqual(output.stat().st_mode & 0o777, 0o444)
            with self.assertRaisesRegex(SUMMARY.SummaryError, "already exists"):
                SUMMARY.write_json_and_sidecar_exclusive(output, {"all_pass": True})


if __name__ == "__main__":
    unittest.main()
