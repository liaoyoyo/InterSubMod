#!/usr/bin/env python3
"""Bounded comparison of recurrence-allowed and strict perfect-phylogeny models.

This script is an exploratory pilot. It reads canonical-v5 JSON artifacts but
never mutates them. The strict model adds the rooted three-gamete condition to
the existing minimum-extra-vertex MILP contract.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import importlib.util
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack


REPO_ROOT = Path(__file__).resolve().parents[3]
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANONICAL_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5"
)
EXACT_SOLVER_PATH = (
    REPO_ROOT
    / "research/20260716_read_linked_hypercube_exact_likelihood_validation"
    / "scripts/hypercube_exact.py"
)
TARGET_SAMPLES = ("H2009", "COLO829")
TARGET_M = (31, 63, 127, 255)


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("hypercube_exact_current", path)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise RuntimeError(f"Cannot load module: {path}")
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H = load_module(EXACT_SOLVER_PATH)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_or_none(value: Any) -> float | None:
    if value is None:
        return None
    numeric = float(value)
    return numeric if math.isfinite(numeric) else None


def popcount(value: int) -> int:
    return bin(int(value)).count("1")


def parse_cap_pool(cap_reason: str) -> int | None:
    match = re.search(r"C\((\d+),[34]\)", cap_reason)
    if match:
        return int(match.group(1))
    match = re.search(r"extra_cap=\d+\(pool=(\d+)\)", cap_reason)
    return int(match.group(1)) if match else None


def select_real_cases(canonical_root: Path) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for sample in TARGET_SAMPLES:
        layered_path = (
            canonical_root
            / "samples"
            / sample
            / f"layered_reconstruction_{sample}.json"
        )
        with layered_path.open(encoding="utf-8") as handle:
            layered = json.load(handle)
        buckets: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for detail in layered["detail"]:
            if not detail.get("is_primary_lineage") or not detail.get("capped"):
                continue
            pool = parse_cap_pool(str(detail.get("cap_reason") or ""))
            if pool in TARGET_M:
                buckets[int(pool)].append(detail)
        for pool in TARGET_M:
            ordered = sorted(
                buckets[pool],
                key=lambda detail: (
                    int(str(detail["chrom"])[3:]),
                    int(detail["start"]),
                    int(detail["end"]),
                    str(detail["family"]),
                ),
            )
            if not ordered:
                raise RuntimeError(f"Missing real pilot stratum: {sample}, M={pool}")
            selected.append({"sample": sample, "expected_pool": pool, "detail": ordered[0]})
    return selected


def load_group_map(canonical_root: Path) -> dict[tuple[str, str, int, int], dict[str, Any]]:
    group_map: dict[tuple[str, str, int, int], dict[str, Any]] = {}
    for sample in TARGET_SAMPLES:
        pattern = str(canonical_root / "samples" / sample / "mlhp_part_*.json")
        for filename in sorted(glob.glob(pattern)):
            with Path(filename).open(encoding="utf-8") as handle:
                document = json.load(handle)
            for group in document["groups"]:
                key = (sample, group["chrom"], int(group["start"]), int(group["end"]))
                group_map[key] = group
    return group_map


def make_real_instance(
    selected: dict[str, Any],
    group_map: dict[tuple[str, str, int, int], dict[str, Any]],
) -> dict[str, Any]:
    sample = selected["sample"]
    detail = selected["detail"]
    family = str(detail["family"])
    key = (sample, detail["chrom"], int(detail["start"]), int(detail["end"]))
    group = group_map[key]
    full_counts = dict((group.get("populations_by_hp") or {}).get(family, {}) or {})
    partial_counts = dict((group.get("subread_groups_by_hp") or {}).get(family, {}) or {})
    patterns = list(full_counts) + list(partial_counts)
    if not patterns:
        raise RuntimeError(f"No structural patterns for {key}, HP{family}")
    k = len(patterns[0])
    if any(len(pattern) != k for pattern in patterns):
        raise RuntimeError(f"Inconsistent pattern width for {key}, HP{family}")
    instance = {
        "case_key": (
            f"{sample}|{detail['chrom']}:{detail['start']}-{detail['end']}|HP{family}"
        ),
        "sample": sample,
        "chrom": detail["chrom"],
        "start": int(detail["start"]),
        "end": int(detail["end"]),
        "family": family,
        "k": k,
        "full": sorted(full_counts),
        "partial": sorted(partial_counts),
        "full_counts": full_counts,
        "partial_counts": partial_counts,
        "legacy_cap_reason": detail.get("cap_reason"),
        "legacy_greedy_hidden": detail.get("n_hidden"),
        "expected_pool": int(selected["expected_pool"]),
    }
    built = build_base_problem(instance)
    observed_pool = len(built["vertices"]) - len(built["mandatory"])
    if observed_pool != instance["expected_pool"]:
        raise AssertionError(
            f"Pool mismatch for {instance['case_key']}: "
            f"{observed_pool} != {instance['expected_pool']}"
        )
    return instance


def make_synthetic_instance(
    case_key: str,
    full: Sequence[str],
    partial: Sequence[str],
    k: int,
) -> dict[str, Any]:
    return {
        "case_key": case_key,
        "sample": "SYNTHETIC",
        "chrom": None,
        "start": None,
        "end": None,
        "family": None,
        "k": int(k),
        "full": list(full),
        "partial": list(partial),
        "full_counts": {pattern: 1 for pattern in full},
        "partial_counts": {pattern: 1 for pattern in partial},
        "legacy_cap_reason": None,
        "legacy_greedy_hidden": None,
        "expected_pool": None,
    }


def build_base_problem(
    instance: dict[str, Any],
    *,
    objective_value: int | None = None,
    excluded_vertex_sets: Sequence[Sequence[int]] = (),
) -> dict[str, Any]:
    groups = tuple(H.SymbolicPattern.from_string(pattern) for pattern in instance["partial"])
    (
        vertices,
        objective,
        bounds,
        constraints,
        universe_mask,
        reduction,
        n_no_good_nonzeros,
    ) = H._build_problem(
        tuple(instance["full"]),
        groups,
        int(instance["k"]),
        "observed_alt",
        objective_value=objective_value,
        excluded_vertex_sets=excluded_vertex_sets,
    )
    mandatory = {0, *(H.parse_full(pattern) for pattern in instance["full"])}
    return {
        "vertices": tuple(int(value) for value in vertices),
        "objective": np.asarray(objective, dtype=float),
        "bounds": bounds,
        "constraints": constraints,
        "universe_mask": int(universe_mask),
        "groups": groups,
        "mandatory": mandatory,
        "reduction": reduction,
        "n_no_good_nonzeros": int(n_no_good_nonzeros),
    }


def strict_perfect_checker(
    instance: dict[str, Any],
    selected_vertices: Iterable[int],
    universe_mask: int,
) -> dict[str, Any]:
    selected = set(int(value) for value in selected_vertices)
    mandatory = {0, *(H.parse_full(pattern) for pattern in instance["full"])}
    missing_mandatory = sorted(mandatory - selected)
    groups = [H.SymbolicPattern.from_string(pattern) for pattern in instance["partial"]]
    uncovered_groups = [
        group.pattern
        for group in groups
        if not any(group.compatible(vertex, universe_mask) for vertex in selected)
    ]
    closure_failures = []
    for vertex in sorted(selected):
        if vertex == 0:
            continue
        predecessors = [
            vertex ^ (1 << bit)
            for bit in range(int(instance["k"]))
            if vertex & (1 << bit)
        ]
        if not any(predecessor in selected for predecessor in predecessors):
            closure_failures.append(vertex)
    active_bits = [
        bit for bit in range(int(instance["k"])) if universe_mask & (1 << bit)
    ]
    three_gamete_failures = []
    for offset, left in enumerate(active_bits):
        for right in active_bits[offset + 1 :]:
            observed = {
                (
                    1 if vertex & (1 << left) else 0,
                    1 if vertex & (1 << right) else 0,
                )
                for vertex in selected
            }
            if {(1, 0), (0, 1), (1, 1)} <= observed:
                three_gamete_failures.append([left, right])
    return {
        "pass": not (
            missing_mandatory
            or uncovered_groups
            or closure_failures
            or three_gamete_failures
        ),
        "missing_mandatory": missing_mandatory,
        "uncovered_groups": uncovered_groups,
        "closure_failures": closure_failures,
        "rooted_three_gamete_failures": three_gamete_failures,
    }


def augment_with_perfect_constraints(base: dict[str, Any]):
    vertices = base["vertices"]
    n_z = len(vertices)
    active_bits = [
        bit
        for bit in range(max(1, int(base["universe_mask"]).bit_length()))
        if base["universe_mask"] & (1 << bit)
    ]
    pairs = [
        (left, right)
        for offset, left in enumerate(active_bits)
        for right in active_bits[offset + 1 :]
    ]
    categories = ((1, 0), (0, 1), (1, 1))
    q_index: dict[tuple[int, int, tuple[int, int]], int] = {}
    for left, right in pairs:
        for category in categories:
            q_index[(left, right, category)] = n_z + len(q_index)
    n_q = len(q_index)

    base_matrix = hstack(
        [
            csr_matrix(base["constraints"].A),
            csr_matrix((base["constraints"].A.shape[0], n_q)),
        ],
        format="csr",
    )
    row_idx: list[int] = []
    col_idx: list[int] = []
    data: list[float] = []
    row_lb: list[float] = []
    row_ub: list[float] = []
    row = 0
    for left, right in pairs:
        for category in categories:
            member_indices = [
                index
                for index, vertex in enumerate(vertices)
                if (
                    1 if vertex & (1 << left) else 0,
                    1 if vertex & (1 << right) else 0,
                )
                == category
            ]
            for index in member_indices:
                row_idx.append(row)
                col_idx.append(index)
                data.append(1.0)
            row_idx.append(row)
            col_idx.append(q_index[(left, right, category)])
            data.append(-float(max(1, len(member_indices))))
            row_lb.append(-np.inf)
            row_ub.append(0.0)
            row += 1
        for category in categories:
            row_idx.append(row)
            col_idx.append(q_index[(left, right, category)])
            data.append(1.0)
        row_lb.append(-np.inf)
        row_ub.append(2.0)
        row += 1
    perfect_matrix = coo_matrix(
        (data, (row_idx, col_idx)), shape=(row, n_z + n_q)
    ).tocsr()
    matrix = vstack([base_matrix, perfect_matrix], format="csr")
    lower = np.concatenate(
        [np.asarray(base["constraints"].lb), np.asarray(row_lb, dtype=float)]
    )
    upper = np.concatenate(
        [np.asarray(base["constraints"].ub), np.asarray(row_ub, dtype=float)]
    )
    objective = np.concatenate([base["objective"], np.zeros(n_q, dtype=float)])
    bounds = Bounds(
        np.concatenate([np.asarray(base["bounds"].lb), np.zeros(n_q, dtype=float)]),
        np.concatenate([np.asarray(base["bounds"].ub), np.ones(n_q, dtype=float)]),
    )
    return (
        objective,
        bounds,
        LinearConstraint(matrix, lower, upper),
        n_z,
        n_q,
        len(perfect_matrix.indptr) - 1,
    )


def solve_strict_perfect(
    instance: dict[str, Any],
    *,
    time_limit_seconds: float,
    objective_value: int | None = None,
    excluded_vertex_sets: Sequence[Sequence[int]] = (),
) -> dict[str, Any]:
    base = build_base_problem(
        instance,
        objective_value=objective_value,
        excluded_vertex_sets=excluded_vertex_sets,
    )
    objective, bounds, constraints, n_z, n_q, n_perfect_constraints = (
        augment_with_perfect_constraints(base)
    )
    started = time.perf_counter()
    result = milp(
        c=objective,
        integrality=np.ones(len(objective), dtype=np.uint8),
        bounds=bounds,
        constraints=constraints,
        options={
            "time_limit": float(time_limit_seconds),
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )
    runtime = time.perf_counter() - started
    selected = tuple(
        base["vertices"][index]
        for index, value in enumerate(
            (result.x[:n_z] if result.x is not None else [])
        )
        if value > 0.5
    )
    status = H._status_name(int(result.status), bool(result.success))
    checker = (
        strict_perfect_checker(
            instance,
            selected,
            int(base["universe_mask"]),
        )
        if selected
        else None
    )
    return {
        "status": status,
        "status_code": int(result.status),
        "message": str(result.message),
        "objective": finite_or_none(result.fun),
        "objective_bound": finite_or_none(
            getattr(result, "mip_dual_bound", None)
        ),
        "mip_gap": finite_or_none(getattr(result, "mip_gap", None)),
        "mip_node_count": (
            int(result.mip_node_count)
            if getattr(result, "mip_node_count", None) is not None
            else None
        ),
        "runtime_seconds": runtime,
        "selected_vertices": list(selected),
        "vertex_set_id": (
            H.vertex_set_digest(int(instance["k"]), selected) if selected else None
        ),
        "parent_choice_count": (
            int(H.parent_choice_count(selected)) if selected else None
        ),
        "checker": checker,
        "n_z_variables": n_z,
        "n_perfect_indicator_variables": n_q,
        "n_constraints": int(constraints.A.shape[0]),
        "n_perfect_constraints": int(n_perfect_constraints),
        "universe_mask": int(base["universe_mask"]),
        "effective_universe_u": popcount(base["universe_mask"]),
        "n_partial_groups_input": int(base["reduction"].n_input),
        "n_partial_groups_active": len(base["reduction"].groups),
        "n_partial_groups_duplicate_removed": int(
            base["reduction"].n_duplicate_removed
        ),
        "n_partial_groups_dominated_removed": int(
            base["reduction"].n_dominated_removed
        ),
        "n_partial_groups_required_hit_removed": int(
            base["reduction"].n_required_hit_removed
        ),
        "n_partial_groups_forced_removed": int(
            base["reduction"].n_forced_removed
        ),
    }


def solve_recurrence_allowed(
    instance: dict[str, Any], time_limit_seconds: float
) -> dict[str, Any]:
    started = time.perf_counter()
    solution = H.solve_min_hidden(
        instance["full"],
        instance["partial"],
        int(instance["k"]),
        universe_mode="observed_alt",
        time_limit_seconds=float(time_limit_seconds),
    )
    runtime = time.perf_counter() - started
    checker = None
    if solution.selected_vertices:
        checker = strict_perfect_checker(
            instance,
            solution.selected_vertices,
            int(solution.universe_mask),
        )
    return {
        "status": solution.status,
        "status_code": solution.status_code,
        "message": solution.message,
        "objective": solution.objective,
        "objective_bound": solution.objective_bound,
        "mip_gap": solution.mip_gap,
        "mip_node_count": solution.mip_node_count,
        "runtime_seconds": runtime,
        "selected_vertices": list(solution.selected_vertices),
        "vertex_set_id": solution.vertex_set_id,
        "parent_choice_count": (
            int(H.parent_choice_count(solution.selected_vertices))
            if solution.selected_vertices
            else None
        ),
        "structural_checker_pass": bool(
            checker
            and not checker["missing_mandatory"]
            and not checker["uncovered_groups"]
            and not checker["closure_failures"]
        ),
        "strict_perfect_checker": checker,
        "n_z_variables": solution.n_variables,
        "n_constraints": solution.n_constraints,
        "universe_mask": solution.universe_mask,
        "effective_universe_u": popcount(solution.universe_mask),
        "n_partial_groups_input": solution.n_partial_groups_input,
        "n_partial_groups_active": solution.n_partial_groups_active,
        "n_partial_groups_duplicate_removed": (
            solution.n_partial_groups_duplicate_removed
        ),
        "n_partial_groups_dominated_removed": (
            solution.n_partial_groups_dominated_removed
        ),
        "n_partial_groups_required_hit_removed": (
            solution.n_partial_groups_required_hit_removed
        ),
        "n_partial_groups_forced_removed": (
            solution.n_partial_groups_forced_removed
        ),
    }


def enumerate_strict_optimal_sets(
    instance: dict[str, Any],
    h_star: int,
    *,
    total_time_limit_seconds: float,
    per_solve_time_limit_seconds: float,
    max_sets: int,
) -> dict[str, Any]:
    if h_star == 0:
        return {
            "status": "COMPLETE_ZERO_OBJECTIVE",
            "complete": True,
            "h_star": 0,
            "n_vertex_sets": 1,
            "total_parent_trees": 1,
            "runtime_seconds": 0.0,
        }
    started = time.perf_counter()
    found: list[tuple[int, ...]] = []
    parent_counts: list[int] = []
    stop_status: Any = None
    while len(found) < max_sets:
        remaining = total_time_limit_seconds - (time.perf_counter() - started)
        if remaining <= 0:
            stop_status = "GLOBAL_TIME_LIMIT"
            break
        solution = solve_strict_perfect(
            instance,
            time_limit_seconds=min(per_solve_time_limit_seconds, remaining),
            objective_value=h_star,
            excluded_vertex_sets=found,
        )
        if solution["status"] == "INFEASIBLE":
            return {
                "status": "COMPLETE_INFEASIBLE_AFTER_NO_GOODS",
                "complete": True,
                "h_star": h_star,
                "n_vertex_sets": len(found),
                "total_parent_trees": sum(parent_counts),
                "all_parent_counts_one": all(value == 1 for value in parent_counts),
                "runtime_seconds": time.perf_counter() - started,
            }
        if solution["status"] != "OPTIMAL":
            stop_status = {
                "status": solution["status"],
                "objective_bound": solution["objective_bound"],
                "mip_gap": solution["mip_gap"],
            }
            break
        selected = tuple(sorted(int(value) for value in solution["selected_vertices"]))
        if selected in found:
            raise AssertionError("No-good cut failed to exclude a prior strict solution")
        found.append(selected)
        parent_counts.append(int(solution["parent_choice_count"]))
    if len(found) >= max_sets:
        stop_status = "MAX_SETS_REACHED"
    return {
        "status": "INCOMPLETE",
        "complete": False,
        "stop_status": stop_status,
        "h_star": h_star,
        "n_vertex_sets": len(found),
        "parent_tree_count_lower_bound": sum(parent_counts),
        "all_parent_counts_one": all(value == 1 for value in parent_counts),
        "runtime_seconds": time.perf_counter() - started,
    }


def summarize_instance(
    instance: dict[str, Any],
    recurrence_allowed: dict[str, Any],
    strict_perfect: dict[str, Any],
) -> dict[str, Any]:
    status_differs = recurrence_allowed["status"] != strict_perfect["status"]
    model_a_optimum_is_strict = bool(
        recurrence_allowed["status"] == "OPTIMAL"
        and recurrence_allowed["selected_vertices"]
        and recurrence_allowed["strict_perfect_checker"]
        and recurrence_allowed["strict_perfect_checker"]["pass"]
    )
    proven_change = (
        recurrence_allowed["status"] == "OPTIMAL"
        and strict_perfect["status"] == "INFEASIBLE"
    ) or (
        recurrence_allowed["status"] == "OPTIMAL"
        and strict_perfect["status"] == "OPTIMAL"
        and round(float(recurrence_allowed["objective"]))
        != round(float(strict_perfect["objective"]))
    )
    return {
        "case_key": instance["case_key"],
        "sample": instance["sample"],
        "chrom": instance["chrom"],
        "start": instance["start"],
        "end": instance["end"],
        "family": instance["family"],
        "k": instance["k"],
        "full_patterns": instance["full"],
        "partial_patterns": instance["partial"],
        "n_full_patterns": len(instance["full"]),
        "n_partial_patterns": len(instance["partial"]),
        "legacy_cap_reason": instance["legacy_cap_reason"],
        "legacy_greedy_hidden": instance["legacy_greedy_hidden"],
        "expected_pool": instance["expected_pool"],
        "recurrence_allowed": recurrence_allowed,
        "strict_perfect": strict_perfect,
        "solver_status_differs": status_differs,
        "strict_proven_changes_feasibility_or_objective": proven_change,
        "strict_optimum_cross_certified_by_model_a": model_a_optimum_is_strict,
        "strict_optimum_cross_certified_objective": (
            int(round(float(recurrence_allowed["objective"])))
            if model_a_optimum_is_strict
            else None
        ),
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    wall_started = time.perf_counter()
    synthetic_instances = [
        make_synthetic_instance("SYN_sisters_10_01", ["AR", "RA"], [], 2),
        make_synthetic_instance(
            "SYN_recurrence_10_01_11", ["AR", "RA", "AA"], [], 2
        ),
        make_synthetic_instance("SYN_partial_AX_XA", [], ["AX", "XA"], 2),
    ]
    real_selected = select_real_cases(args.canonical_root)
    group_map = load_group_map(args.canonical_root)
    real_instances = [
        make_real_instance(selected, group_map) for selected in real_selected
    ]

    synthetic_results = []
    for instance in synthetic_instances:
        model_a = solve_recurrence_allowed(instance, args.phase_a_time_limit)
        model_b = solve_strict_perfect(
            instance, time_limit_seconds=args.phase_a_time_limit
        )
        synthetic_results.append(summarize_instance(instance, model_a, model_b))
        print(
            "SYNTHETIC",
            instance["case_key"],
            f"A={model_a['status']}:{model_a['objective']}",
            f"B={model_b['status']}:{model_b['objective']}",
            flush=True,
        )

    real_results = []
    for instance in real_instances:
        model_a = solve_recurrence_allowed(instance, args.phase_a_time_limit)
        model_b = solve_strict_perfect(
            instance, time_limit_seconds=args.phase_a_time_limit
        )
        real_results.append(summarize_instance(instance, model_a, model_b))
        print(
            "REAL",
            instance["case_key"],
            f"M={instance['expected_pool']}",
            f"A={model_a['status']}:{model_a['objective']}",
            f"B={model_b['status']}:{model_b['objective']}",
            flush=True,
        )

    strict_phase_b = []
    by_key = {instance["case_key"]: instance for instance in real_instances}
    for result in real_results:
        if (
            result["expected_pool"] == 31
            and result["strict_perfect"]["status"] == "OPTIMAL"
        ):
            h_star = int(round(float(result["strict_perfect"]["objective"])))
            enumeration = enumerate_strict_optimal_sets(
                by_key[result["case_key"]],
                h_star,
                total_time_limit_seconds=args.phase_b_total_time_limit,
                per_solve_time_limit_seconds=args.phase_b_per_solve_time_limit,
                max_sets=args.max_sets,
            )
            strict_phase_b.append(
                {"case_key": result["case_key"], **enumeration}
            )
            print(
                "STRICT_PHASE_B",
                result["case_key"],
                enumeration["status"],
                f"sets={enumeration['n_vertex_sets']}",
                flush=True,
            )

    summary = {
        "synthetic_cases": len(synthetic_results),
        "real_cases": len(real_results),
        "real_model_a_optimal": sum(
            row["recurrence_allowed"]["status"] == "OPTIMAL"
            for row in real_results
        ),
        "real_model_b_optimal": sum(
            row["strict_perfect"]["status"] == "OPTIMAL"
            for row in real_results
        ),
        "real_model_a_limit": sum(
            row["recurrence_allowed"]["status"] == "LIMIT_REACHED"
            for row in real_results
        ),
        "real_model_b_limit": sum(
            row["strict_perfect"]["status"] == "LIMIT_REACHED"
            for row in real_results
        ),
        "real_model_b_infeasible": sum(
            row["strict_perfect"]["status"] == "INFEASIBLE"
            for row in real_results
        ),
        "real_solver_status_differs": sum(
            row["solver_status_differs"]
            for row in real_results
        ),
        "real_strict_proven_changes_feasibility_or_objective": sum(
            row["strict_proven_changes_feasibility_or_objective"]
            for row in real_results
        ),
        "real_both_models_certified_optimal": sum(
            row["recurrence_allowed"]["status"] == "OPTIMAL"
            and row["strict_perfect"]["status"] == "OPTIMAL"
            for row in real_results
        ),
        "real_both_certified_same_objective": sum(
            row["recurrence_allowed"]["status"] == "OPTIMAL"
            and row["strict_perfect"]["status"] == "OPTIMAL"
            and round(float(row["recurrence_allowed"]["objective"]))
            == round(float(row["strict_perfect"]["objective"]))
            for row in real_results
        ),
        "real_strict_optimum_cross_certified_by_model_a": sum(
            row["strict_optimum_cross_certified_by_model_a"]
            for row in real_results
        ),
        "all_model_a_incumbents_structurally_valid": all(
            row["recurrence_allowed"]["structural_checker_pass"]
            for row in real_results
            if row["recurrence_allowed"]["selected_vertices"]
        ),
        "all_model_b_incumbents_strictly_valid": all(
            row["strict_perfect"]["checker"]["pass"]
            for row in real_results
            if row["strict_perfect"]["selected_vertices"]
        ),
    }
    return {
        "schema": "intersubmod.perfect_phylogeny_constraint_pilot.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "A_EXPLORATORY_PILOT",
        "scope_flag": "PARTIAL_NOT_COMPREHENSIVE_VALIDATION",
        "claim_ceiling": (
            "Algorithmic feasibility and bounded runtime evidence only; "
            "not biological truth and not a clone/subclone prevalence estimate."
        ),
        "inputs": {
            "canonical_root": str(args.canonical_root),
            "exact_solver_path": str(EXACT_SOLVER_PATH),
            "exact_solver_sha256": sha256(EXACT_SOLVER_PATH),
            "samples": list(TARGET_SAMPLES),
            "target_extra_pools": list(TARGET_M),
        },
        "parameters": {
            "phase_a_time_limit_seconds_per_solve": args.phase_a_time_limit,
            "phase_b_total_time_limit_seconds_per_case": (
                args.phase_b_total_time_limit
            ),
            "phase_b_time_limit_seconds_per_solve": (
                args.phase_b_per_solve_time_limit
            ),
            "max_sets": args.max_sets,
            "universe_mode": "observed_alt",
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
        "models": {
            "A_recurrence_allowed": (
                "ROOT/full mandatory + partial-group coverage + selected-node "
                "Hamming-1 predecessor closure; a mutation may arise on multiple branches."
            ),
            "B_strict_perfect": (
                "Model A plus rooted three-gamete exclusion for every active "
                "mutation pair; equivalent here to no global recurrent acquisition."
            ),
        },
        "environment": {
            "scipy_version": scipy.__version__,
            "numpy_version": np.__version__,
        },
        "synthetic_results": synthetic_results,
        "real_results": real_results,
        "strict_phase_b": strict_phase_b,
        "summary": summary,
        "wall_runtime_seconds": time.perf_counter() - wall_started,
        "limitations": [
            "Eight deterministic real units are a pilot, not a population estimate.",
            "LIMIT_REACHED is an incumbent and bound, not a certified optimum.",
            "Phase B is complete only after an INFEASIBLE no-more-solutions certificate.",
            "Strict perfect phylogeny is a changed biological/model contract, not safe preprocessing.",
            "CNA, LOH, mutation loss, recurrence, and allele-specific amplification can make strict perfect phylogeny biologically false.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--canonical-root",
        type=Path,
        default=DEFAULT_CANONICAL_ROOT,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "results" / "pilot_receipt.json",
    )
    parser.add_argument("--phase-a-time-limit", type=float, default=8.0)
    parser.add_argument("--phase-b-total-time-limit", type=float, default=15.0)
    parser.add_argument("--phase-b-per-solve-time-limit", type=float, default=3.0)
    parser.add_argument("--max-sets", type=int, default=512)
    args = parser.parse_args()
    if (
        args.phase_a_time_limit <= 0
        or args.phase_b_total_time_limit <= 0
        or args.phase_b_per_solve_time_limit <= 0
        or args.max_sets < 1
    ):
        raise ValueError("All time limits and max-sets must be positive")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    receipt = run(args)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print("OUTPUT", args.output, flush=True)
    print("SHA256", sha256(args.output), flush=True)
    print("SUMMARY", json.dumps(receipt["summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
