#!/usr/bin/env python3
"""Build a reproducible receipt for adaptive hypercube solver routing.

The receipt deliberately separates:
1. mathematical operation/memory proxies,
2. bounded real-fixture solver measurements,
3. canonical-v5 topology counts,
4. historical layered-v2 VAF-ranking counts.

It does not change or execute the production solver.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List


REPO_ROOT = Path(__file__).resolve().parents[3]
PROBE_ROOT = REPO_ROOT / "research" / "20260718_solver_methyl_edge_probe"

DEFAULT_INPUTS = {
    "solver_probe": PROBE_ROOT / "results" / "solver_probe_receipt.json",
    "current_scipy": PROBE_ROOT / "results" / "current_scipy_real_benchmark.json",
    "subcube_oracle": PROBE_ROOT / "results" / "parent_mapping_subcube_plan_oracle_receipt.json",
    "canonical_v5": (
        REPO_ROOT
        / "research"
        / "20260710_layered_reconstruction_v2"
        / "current_layered_topology_v3_raw_all_v1.json"
    ),
    "historical_vaf": (
        REPO_ROOT
        / "research"
        / "20260711_read_group_C_tree_T_topology_report"
        / "data"
        / "vaf_top_tie_census.json"
    ),
    "historical_final_shape": (
        REPO_ROOT
        / "research"
        / "20260712_vaf_selected_shape_four_class_census"
        / "data"
        / "20260712_vaf_final_single_shape_four_class_census.json"
    ),
    "m2_no_go_record": (
        REPO_ROOT
        / "research"
        / "20260716_read_linked_hypercube_exact_likelihood_validation"
        / "20260717_M2_frozen_v2正式pilot_NO_GO驗證紀錄_01.md"
    ),
}


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def proxy_ops(m: int, q: int) -> int:
    """Terminal-subset DP operation-count proxy."""
    vertices = 1 << m
    edges = m * (1 << (m - 1))
    return (3**q) * vertices + (1 << q) * edges


def dense_bytes(m: int, q: int, bytes_per_cell: int) -> int:
    return (1 << q) * (1 << m) * bytes_per_cell


def maximum_m_for_gate(
    q: int,
    max_ops: int,
    max_bytes: int,
    bytes_per_cell: int,
    search_max_m: int = 30,
) -> int:
    allowed = [
        m
        for m in range(1, search_max_m + 1)
        if proxy_ops(m, q) <= max_ops
        and dense_bytes(m, q, bytes_per_cell) <= max_bytes
    ]
    return max(allowed) if allowed else 0


def full_cube_parent_log10(m: int) -> float:
    return sum(math.comb(m, rank) * math.log10(rank) for rank in range(1, m + 1))


def parse_number(text: str) -> float:
    return float(text.replace(",", ""))


def parse_m2_tail_record(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    completed = re.search(
        r"Candidate＋likelihood皆invoked\s*\|\s*([\d,]+)\s*\|\s*([\d,.]+) s\s*\|\s*([\d,.]+) s",
        text,
    )
    tails = re.search(
        r"Candidate invoked、likelihood未invoked\s*\|\s*([\d,]+)\s*\|\s*([\d,.]+) s\s*\|\s*([\d,.]+) s",
        text,
    )
    not_invoked = re.search(
        r"Solver未invoked\s*\|\s*([\d,]+)\s*\|\s*([\d,.]+) s\s*\|\s*([\d,.]+) s",
        text,
    )
    share = re.search(r"占minread=1 candidate-generation時間的`([\d.]+)%`", text)
    if not all((completed, tails, not_invoked, share)):
        raise RuntimeError(f"Could not parse audited M2 tail table from {path}")

    completed_units = int(parse_number(completed.group(1)))
    completed_candidate_seconds = parse_number(completed.group(2))
    completed_likelihood_seconds = parse_number(completed.group(3))
    tail_units = int(parse_number(tails.group(1)))
    tail_candidate_seconds = parse_number(tails.group(2))
    not_invoked_units = int(parse_number(not_invoked.group(1)))
    total_units = completed_units + tail_units + not_invoked_units
    total_candidate_seconds = completed_candidate_seconds + tail_candidate_seconds
    reported_share = parse_number(share.group(1)) / 100.0
    computed_share = tail_candidate_seconds / total_candidate_seconds
    if abs(reported_share - computed_share) > 1e-6:
        raise RuntimeError(
            f"Tail-time share mismatch: reported={reported_share}, computed={computed_share}"
        )

    return {
        "grain": "component basis × bridge threshold × component, structural-minread=1",
        "diagnostic_only": True,
        "total_unit_instances": total_units,
        "completed_candidate_and_likelihood_units": completed_units,
        "tail_candidate_only_units": tail_units,
        "solver_not_invoked_units": not_invoked_units,
        "tail_unit_fraction": tail_units / total_units,
        "completed_candidate_seconds": completed_candidate_seconds,
        "tail_candidate_seconds": tail_candidate_seconds,
        "total_candidate_seconds": total_candidate_seconds,
        "completed_likelihood_seconds": completed_likelihood_seconds,
        "tail_candidate_time_fraction": computed_share,
        "candidate_time_upper_bound_reduction_if_tail_cost_zero": (
            total_candidate_seconds / completed_candidate_seconds
        ),
        "candidate_plus_likelihood_upper_bound_reduction_if_tail_cost_zero": (
            (total_candidate_seconds + completed_likelihood_seconds)
            / (completed_candidate_seconds + completed_likelihood_seconds)
        ),
        "claim_limit": (
            "Truncated frozen-v2 run without ranking receipt; valid only for "
            "performance diagnosis, not topology or biological proportions."
        ),
    }


def index_by(items: Iterable[Dict[str, Any]], key: str) -> Dict[str, Dict[str, Any]]:
    return {str(item[key]): item for item in items}


def build_receipt(max_ops: int, max_mib: int, bytes_per_cell: int) -> Dict[str, Any]:
    inputs = {name: path.resolve() for name, path in DEFAULT_INPUTS.items()}
    for path in inputs.values():
        if not path.exists():
            raise FileNotFoundError(path)

    solver_probe = load_json(inputs["solver_probe"])
    current_scipy = load_json(inputs["current_scipy"])
    subcube_oracle = load_json(inputs["subcube_oracle"])
    canonical_v5 = load_json(inputs["canonical_v5"])
    historical_vaf = load_json(inputs["historical_vaf"])
    historical_shape = load_json(inputs["historical_final_shape"])
    m2_tail = parse_m2_tail_record(inputs["m2_no_go_record"])

    solver_cases = index_by(solver_probe["real"], "case_id")
    scipy_cases = index_by(current_scipy["cases"], "case_id")
    subcube_cases = index_by(
        subcube_oracle["real_fixture_census"]["cases"], "case_id"
    )

    real_cases: List[Dict[str, Any]] = []
    for case_id in sorted(set(solver_cases) & set(scipy_cases) & set(subcube_cases)):
        solver_case = solver_cases[case_id]
        scipy_case = scipy_cases[case_id]
        subcube_case = subcube_cases[case_id]
        m = int(solver_case["effective_m"])
        raw_k = int(subcube_case["raw_k"])
        hidden_optimum = int(solver_case["bnb"]["objective"])
        free_vertices = (1 << m) - 1
        naive_layers = sum(
            math.comb(free_vertices, size)
            for size in range(0, hidden_optimum + 1)
        )
        visited = int(solver_case["bnb"]["stats"]["visited_states"])
        bnb_seconds = float(solver_case["bnb"]["elapsed_seconds"])
        scipy_seconds = float(scipy_case["elapsed_seconds"])
        real_cases.append(
            {
                "case_id": case_id,
                "raw_k": raw_k,
                "effective_m": m,
                "raw_vertices": 1 << raw_k,
                "active_vertices": 1 << m,
                "active_vertex_reduction_fraction": 1.0 - (1 << m) / (1 << raw_k),
                "raw_edges": raw_k * (1 << (raw_k - 1)),
                "active_edges": m * (1 << (m - 1)),
                "active_edge_reduction_fraction": (
                    1.0
                    - (m * (1 << (m - 1)))
                    / (raw_k * (1 << (raw_k - 1)))
                ),
                "input_groups": int(subcube_case["input_groups"]),
                "reduced_groups": int(subcube_case["reduced_groups"]),
                "group_reduction_fraction": (
                    1.0
                    - int(subcube_case["reduced_groups"])
                    / int(subcube_case["input_groups"])
                ),
                "minimum_extra_states": hidden_optimum,
                "naive_subset_layers_through_h_star": naive_layers,
                "bnb_visited_states": visited,
                "bnb_search_state_reduction_fraction": 1.0 - visited / naive_layers,
                "bnb_search_state_ratio": naive_layers / visited,
                "optimal_vertex_sets": int(solver_case["bnb"]["n_vertex_sets"]),
                "current_scipy_seconds": scipy_seconds,
                "bnb_seconds": bnb_seconds,
                "single_run_wall_ratio_current_over_bnb": scipy_seconds / bnb_seconds,
                "complete_digest_match": bool(
                    scipy_case["verification_rerun"]["matches_bnb_and_persistent"]
                ),
            }
        )

    max_bytes = max_mib * 1024 * 1024
    adaptive_rows = []
    for q in range(1, 11):
        max_m = maximum_m_for_gate(
            q=q,
            max_ops=max_ops,
            max_bytes=max_bytes,
            bytes_per_cell=bytes_per_cell,
        )
        adaptive_rows.append(
            {
                "q": q,
                "maximum_m": max_m,
                "proxy_ops_at_maximum_m": proxy_ops(max_m, q) if max_m else 0,
                "dense_bytes_at_maximum_m": (
                    dense_bytes(max_m, q, bytes_per_cell) if max_m else 0
                ),
                "dense_mib_at_maximum_m": (
                    dense_bytes(max_m, q, bytes_per_cell) / 1024 / 1024
                    if max_m
                    else 0.0
                ),
            }
        )

    full_parent_rows = []
    for m in [1, 2, 3, 4, 5, 6, 8, 10, 12]:
        log10_value = full_cube_parent_log10(m)
        exact_value = None
        if log10_value < 35:
            exact_value = math.prod(
                rank ** math.comb(m, rank) for rank in range(1, m + 1)
            )
        full_parent_rows.append(
            {
                "m": m,
                "vertices": 1 << m,
                "edges": m * (1 << (m - 1)),
                "full_cube_parent_tree_count": exact_value,
                "full_cube_parent_tree_log10": log10_value,
            }
        )

    canonical = canonical_v5["canonical"]
    aggregate = canonical["aggregate"]
    sample_rows = []
    for sample in canonical["samples"]:
        sample_rows.append(
            {
                "sample": sample["sample"],
                "W_primary": int(sample["W_primary"]),
                "complete_regions": int(sample["complete_regions"]),
                "incomplete_regions": int(sample["incomplete_regions"]),
                "incomplete_fraction": (
                    int(sample["incomplete_regions"]) / int(sample["W_primary"])
                ),
                "share_of_all_incomplete": (
                    int(sample["incomplete_regions"])
                    / int(aggregate["incomplete_regions"])
                ),
            }
        )
    sample_rows.sort(key=lambda row: row["incomplete_regions"], reverse=True)

    vaf_aggregate = historical_vaf["aggregate"]
    top_classes = vaf_aggregate["exact_top_classes_all_evaluable"]
    historical_shape_aggregate = historical_shape["aggregate"]
    one_topology = (
        int(top_classes["unique_first"])
        + int(top_classes["tied_first_same_topology"])
    )

    checks = {
        "solver_probe_real_cases_complete": all(
            case["bnb"]["complete"] for case in solver_probe["real"]
        ),
        "solver_fixture_digests_match": all(
            case["complete_digest_match"] for case in real_cases
        ),
        "subcube_oracle_pass": (
            subcube_oracle["real_fixture_census"]["status"] == "PASS"
            and subcube_oracle["downward_closure"]["status"] == "PASS"
            and subcube_oracle["parent_factorization"]["status"] == "PASS"
        ),
        "canonical_v5_pass": bool(canonical_v5["all_pass"]),
        "canonical_conservation": (
            int(aggregate["complete_regions"])
            + int(aggregate["incomplete_regions"])
            == int(aggregate["W_primary"])
        ),
        "historical_vaf_conservation": (
            sum(int(value) for value in top_classes.values())
            == int(vaf_aggregate["evaluable_ambiguous_regions"])
        ),
        "historical_final_shape_conservation": (
            int(historical_shape_aggregate["final_single_shape_regions"])
            + int(historical_shape_aggregate["unresolved_regions"])
            == int(historical_shape_aggregate["complete_regions"])
        ),
        "adaptive_rows_within_gate": all(
            row["proxy_ops_at_maximum_m"] <= max_ops
            and row["dense_bytes_at_maximum_m"] <= max_bytes
            for row in adaptive_rows
        ),
    }
    checks["all_pass"] = all(checks.values())

    return {
        "schema": "intersubmod.adaptive_solver_limit_assessment.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "Task Type A exploratory assessment; formulas plus two bounded real "
            "fixtures and historical diagnostics; not production validation"
        ),
        "goals": ["G3", "G4", "G5"],
        "inputs": {
            name: {
                "path": str(path.relative_to(REPO_ROOT)),
                "sha256": sha256(path),
            }
            for name, path in inputs.items()
        },
        "definitions": {
            "raw_k": "sSNV count before active-bit compression",
            "effective_m": "observed-ALT active dimensions used by the exact solver",
            "q": "reduced partial groups plus mandatory terminal states in the proposed subset-DP route",
            "V": "number of distinct optimal vertex sets; full materialization costs at least Omega(V)",
            "proxy_warning": (
                "DP operation count and dense 16-byte cells are planning proxies, "
                "not measured runtime or Python memory."
            ),
        },
        "theory": {
            "vertices_formula": "P=2^m",
            "edges_formula": "E=m*2^(m-1)",
            "subset_dp_time_proxy": "3^q*2^m + 2^q*m*2^(m-1)",
            "subset_dp_space_cells": "2^q*2^m",
            "output_lower_bound": "Omega(V)",
            "adaptive_gate": {
                "max_proxy_ops": max_ops,
                "max_dense_mib": max_mib,
                "bytes_per_dense_cell": bytes_per_cell,
                "rows": adaptive_rows,
                "status": "PROPOSED_PILOT_GATE_NOT_RUNTIME_GUARANTEE",
            },
            "full_cube_parent_assignments": full_parent_rows,
            "unobserved_unary_chain_equivalence": [
                {
                    "k": k,
                    "ordered_hamming_chains": math.factorial(k),
                    "unordered_multi_mutation_edge_equivalence_classes": 1,
                    "model_change": True,
                }
                for k in [4, 8, 12]
            ],
        },
        "bounded_real_fixture_evidence": {
            "cases": real_cases,
            "limitations": current_scipy["limitations"],
            "downward_closure_real_fixture_reduction": "0% in both fixtures",
        },
        "m2_frozen_v2_tail_diagnostic": m2_tail,
        "canonical_v5_context": {
            "version": "LongPhase-S recalibrated FILTER=PASS v5",
            "scope": canonical_v5["scope"],
            "W_primary": int(aggregate["W_primary"]),
            "complete_regions": int(aggregate["complete_regions"]),
            "incomplete_regions": int(aggregate["incomplete_regions"]),
            "incomplete_fraction": (
                int(aggregate["incomplete_regions"]) / int(aggregate["W_primary"])
            ),
            "primary_units": int(aggregate["primary_units"]),
            "partial_only_units": int(aggregate["primary_units_partial_only"]),
            "partial_only_fraction": (
                int(aggregate["primary_units_partial_only"])
                / int(aggregate["primary_units"])
            ),
            "topology_multiple_exact_multiple": int(
                aggregate["topology_classes"]["topology_multiple_exact_multiple"]
            ),
            "topology_multiple_fraction_of_complete": (
                int(aggregate["topology_classes"]["topology_multiple_exact_multiple"])
                / int(aggregate["complete_regions"])
            ),
            "sample_incomplete_rows": sample_rows,
            "top_two_incomplete_samples": [row["sample"] for row in sample_rows[:2]],
            "top_two_share_of_all_incomplete": (
                sum(row["incomplete_regions"] for row in sample_rows[:2])
                / int(aggregate["incomplete_regions"])
            ),
        },
        "historical_layered_v2_vaf_context": {
            "version": "historical layered-v2; not current canonical v5",
            "complete_regions": int(vaf_aggregate["complete_regions"]),
            "ambiguous_complete_regions": int(
                vaf_aggregate["ambiguous_complete_regions"]
            ),
            "evaluable_ambiguous_regions": int(
                vaf_aggregate["evaluable_ambiguous_regions"]
            ),
            "first_rank_one_topology": one_topology,
            "first_rank_one_topology_fraction": (
                one_topology / int(vaf_aggregate["evaluable_ambiguous_regions"])
            ),
            "cross_topology_ties": int(
                top_classes["tied_first_different_topology"]
            ),
            "historical_final_single_shape_regions": int(
                historical_shape_aggregate["final_single_shape_regions"]
            ),
            "historical_unresolved_regions": int(
                historical_shape_aggregate["unresolved_regions"]
            ),
            "claim_limit": (
                "Useful for ranking-endpoint design only; must not be combined "
                "with canonical-v5 denominators as a current success rate."
            ),
        },
        "recommended_router": [
            {
                "priority": 1,
                "route": "closed-form/trivial bypass",
                "condition": "k=1 or already forced unique state",
                "exactness": "unchanged",
            },
            {
                "priority": 2,
                "route": "objective subset-DP then exact B&B/ZDD enumeration",
                "condition": "small reduced q and predicted operation/memory gate passes",
                "exactness": "requires small-m exhaustive equivalence oracle",
            },
            {
                "priority": 3,
                "route": "obligation-driven exact B&B with bitsets",
                "condition": "partial constraints strong and h* expected small",
                "exactness": "bounded fixtures exact; full tail panel pending",
            },
            {
                "priority": 4,
                "route": "current MILP fallback",
                "condition": "within current effective-m contract and other exact routes decline",
                "exactness": "current reference",
            },
            {
                "priority": 5,
                "route": "compressed family/count-only or abstain",
                "condition": "V exceeds materialization budget or deadline",
                "exactness": (
                    "exact only if representation is uncapped and complete; "
                    "otherwise complete=false"
                ),
            },
        ],
        "recommended_limit_policy": {
            "production_default": "keep effective m<=12 until dual-pilot validation",
            "first_exploratory_relaxation": (
                "pilot m<=15 when q<=6 and predicted gate passes"
            ),
            "second_exploratory_relaxation": (
                "consider m<=17 when q<=4 after packed-memory benchmark"
            ),
            "theoretical_small_q_ceiling_under_proxy": (
                "q<=2 permits m<=19 under the planning proxy; not approved for production"
            ),
            "raw_k_policy": (
                "do not cap raw k when exact active-bit compression yields an "
                "effective m that passes the route gate"
            ),
            "candidate_cap_policy": (
                "do not simply raise max_sets=256; split objective, uniqueness, "
                "count/compressed-family, and full-enumeration modes"
            ),
        },
        "checks": checks,
        "verdict": (
            "PROBE_ADAPTIVE_RELAXATION"
            if checks["all_pass"]
            else "NO_GO_RECEIPT_CHECK_FAILED"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=PROBE_ROOT / "results" / "adaptive_solver_limit_assessment_receipt.json",
    )
    parser.add_argument("--max-proxy-ops", type=int, default=50_000_000)
    parser.add_argument("--max-dense-mib", type=int, default=512)
    parser.add_argument("--bytes-per-dense-cell", type=int, default=16)
    args = parser.parse_args()

    receipt = build_receipt(
        max_ops=args.max_proxy_ops,
        max_mib=args.max_dense_mib,
        bytes_per_cell=args.bytes_per_dense_cell,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    print(
        "PASS"
        f" verdict={receipt['verdict']}"
        f" q_rows={len(receipt['theory']['adaptive_gate']['rows'])}"
        f" real_cases={len(receipt['bounded_real_fixture_evidence']['cases'])}"
        f" output={args.output.resolve()}"
    )
    return 0 if receipt["checks"]["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
