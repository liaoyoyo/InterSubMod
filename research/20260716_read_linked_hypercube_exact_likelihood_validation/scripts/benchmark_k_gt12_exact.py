#!/usr/bin/env python3
"""Resource-bounded k>12 pilot for the exact symbolic hypercube MILP.

This script is deliberately independent from the production ranker/extractor.  The
default run creates two deterministic seeded cases for each k=13..16, executes every
MILP stage in a disposable subprocess, and records explicit timeout/abstention rather
than silently treating an incumbent as an optimum.

The output is an exploratory performance receipt, not biological validation.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import math
import os
import pathlib
import random
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from typing import Any


HERE = pathlib.Path(__file__).resolve().parent
CORE_PATH = HERE / "hypercube_exact.py"
DEFAULT_PRIOR_PILOT = HERE.parent / "results" / "pilot" / "pilot_receipt.json"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def popcount(value: int) -> int:
    return bin(int(value)).count("1")


def state_to_pattern_local(state: int, k: int) -> str:
    return "".join("A" if state & (1 << bit) else "R" for bit in range(k))


def observed_alt_mask_local(full_patterns: list[str], partial_patterns: list[str]) -> int:
    mask = 0
    for pattern in full_patterns + partial_patterns:
        for bit, symbol in enumerate(pattern):
            if symbol == "A":
                mask |= 1 << bit
    return mask


def vertex_set_digest_local(k: int, vertices: list[int] | tuple[int, ...] | set[int]) -> str:
    payload = {"k": int(k), "vertices": sorted(int(vertex) for vertex in vertices)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def build_case(k: int, case_kind: str, seed: int) -> dict[str, Any]:
    """Build a deterministic feasible stress case without invoking SciPy."""
    if k < 2:
        raise ValueError("k must be at least 2")
    if case_kind not in {"sparse_active", "dense_active"}:
        raise ValueError(f"unsupported case_kind: {case_kind}")
    rng = random.Random(seed)
    active_m = min(8, k) if case_kind == "sparse_active" else k
    active_bits = sorted(rng.sample(range(k), active_m))
    active_mask = sum(1 << bit for bit in active_bits)

    # The all-active terminal makes the intended active-bit count auditable.  Two
    # additional terminals create branching pressure without requiring a known
    # objective value in advance.
    states = {active_mask}
    while len(states) < min(3, 1 << active_m):
        state = 0
        for bit in active_bits:
            if rng.random() < 0.45:
                state |= 1 << bit
        if state:
            states.add(state)
    full_patterns = [state_to_pattern_local(state, k) for state in sorted(states)]

    partial_patterns: list[str] = []
    for _ in range(3):
        symbols = ["X"] * k
        n_fixed = max(2, active_m // 3)
        for bit in rng.sample(active_bits, n_fixed):
            symbols[bit] = "A" if rng.random() < 0.5 else "R"
        partial_patterns.append("".join(symbols))

    observed_mask = observed_alt_mask_local(full_patterns, partial_patterns)
    return {
        "case_kind": case_kind,
        "seed": int(seed),
        "k": int(k),
        "universe_mode": "observed_alt",
        "full_patterns": full_patterns,
        "partial_patterns": partial_patterns,
        "intended_active_mask": int(active_mask),
        "observed_alt_mask": int(observed_mask),
        "evidence_active_k": popcount(observed_mask),
        "active_k": popcount(observed_mask),
        "expected_variables": 1 << popcount(observed_mask),
    }


def _load_core():
    sys.path.insert(0, str(HERE))
    from hypercube_exact import solve_min_hidden  # pylint: disable=import-outside-toplevel

    return solve_min_hidden


def _worker_solve(spec: dict[str, Any]) -> dict[str, Any]:
    solve_min_hidden = _load_core()
    started = time.perf_counter()
    result = solve_min_hidden(
        spec["full_patterns"],
        spec["partial_patterns"],
        int(spec["k"]),
        universe_mode=spec["universe_mode"],
        time_limit_seconds=float(spec["solver_time_limit_seconds"]),
    )
    payload = result.to_dict()
    payload["runtime_seconds"] = time.perf_counter() - started
    return payload


def _worker_enumerate(spec: dict[str, Any]) -> dict[str, Any]:
    """Enumerate optimum vertex sets under one global worker deadline."""
    solve_min_hidden = _load_core()
    started = time.perf_counter()
    deadline = started + float(spec["enumeration_time_limit_seconds"])
    max_sets = int(spec["max_sets"])

    remaining = max(0.05, deadline - time.perf_counter())
    first = solve_min_hidden(
        spec["full_patterns"],
        spec["partial_patterns"],
        int(spec["k"]),
        universe_mode=spec["universe_mode"],
        time_limit_seconds=remaining,
    )
    if first.status != "OPTIMAL" or first.objective is None:
        return {
            "status": f"FIRST_{first.status}",
            "complete": False,
            "n_vertex_sets": 0,
            "max_sets": max_sets,
            "runtime_seconds": time.perf_counter() - started,
            "first": first.to_dict(),
        }

    optimum = int(round(first.objective))
    sets = [first.selected_vertices]
    stop_solution: dict[str, Any] | None = None
    status = "MAX_SETS_REACHED"
    complete = False
    while len(sets) < max_sets:
        remaining = deadline - time.perf_counter()
        if remaining <= 0.05:
            status = "GLOBAL_TIME_LIMIT_REACHED"
            break
        candidate = solve_min_hidden(
            spec["full_patterns"],
            spec["partial_patterns"],
            int(spec["k"]),
            universe_mode=spec["universe_mode"],
            time_limit_seconds=remaining,
            objective_value=optimum,
            excluded_vertex_sets=sets,
        )
        stop_solution = candidate.to_dict()
        if candidate.status == "INFEASIBLE":
            status = "COMPLETE_INFEASIBLE_CERTIFICATE"
            complete = True
            break
        if candidate.status != "OPTIMAL":
            status = f"STOP_{candidate.status}"
            break
        sets.append(candidate.selected_vertices)

    return {
        "status": status,
        "complete": complete,
        "objective": optimum,
        "n_vertex_sets": len(sets),
        "max_sets": max_sets,
        "vertex_sets": [list(vertices) for vertices in sets],
        "vertex_set_ids": [vertex_set_digest_local(int(spec["k"]), vertices) for vertices in sets],
        "runtime_seconds": time.perf_counter() - started,
        "first": first.to_dict(),
        "stop_solution": stop_solution,
    }


def worker_main(encoded_spec: str, stage: str) -> int:
    # Nice only this disposable benchmark child.  Thread limits are also set by the
    # parent before NumPy/SciPy import.
    try:
        os.nice(10)
    except OSError:
        pass
    spec = json.loads(base64.urlsafe_b64decode(encoded_spec.encode()).decode())
    payload = _worker_solve(spec) if stage == "solve" else _worker_enumerate(spec)
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    return 0


def run_stage(
    spec: dict[str, Any],
    stage: str,
    hard_timeout_seconds: float,
) -> dict[str, Any]:
    encoded = base64.urlsafe_b64encode(json.dumps(spec, sort_keys=True).encode()).decode()
    command = [sys.executable, str(pathlib.Path(__file__).resolve()), "--worker", stage, "--spec", encoded]
    env = dict(os.environ)
    env.update(
        {
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=float(hard_timeout_seconds),
            env=env,
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "HARD_TIMEOUT",
            "runtime_seconds": time.perf_counter() - started,
            "hard_timeout_seconds": float(hard_timeout_seconds),
        }
    if completed.returncode != 0:
        return {
            "status": "WORKER_ERROR",
            "returncode": completed.returncode,
            "runtime_seconds": time.perf_counter() - started,
            "stderr": completed.stderr[-4000:],
        }
    try:
        payload = json.loads(completed.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        return {
            "status": "WORKER_OUTPUT_ERROR",
            "runtime_seconds": time.perf_counter() - started,
            "error": str(exc),
            "stdout": completed.stdout[-4000:],
            "stderr": completed.stderr[-4000:],
        }
    payload["supervisor_runtime_seconds"] = time.perf_counter() - started
    return payload


def validate_row(row: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    solve = row["solve"]
    if solve.get("n_variables") is not None and solve["n_variables"] != row["expected_variables"]:
        failures.append("n_variables does not equal 2^active_k")
    if solve.get("status") == "OPTIMAL":
        objective = solve.get("objective")
        bound = solve.get("objective_bound")
        gap = solve.get("mip_gap")
        if objective is None or bound is None or not math.isclose(float(objective), float(bound), abs_tol=1e-8):
            failures.append("OPTIMAL without matching primal/dual objective")
        if gap is not None and float(gap) > 1e-9:
            failures.append("OPTIMAL with nonzero MIP gap")
    enumeration = row.get("enumeration", {})
    if enumeration.get("complete") and enumeration.get("status") != "COMPLETE_INFEASIBLE_CERTIFICATE":
        failures.append("enumeration completeness lacks infeasibility certificate")
    if enumeration.get("n_vertex_sets", 0) > enumeration.get("max_sets", row.get("max_sets", 0)):
        failures.append("enumeration exceeded max_sets")
    return failures


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=pathlib.Path)
    parser.add_argument("--k-min", type=int, default=13)
    parser.add_argument("--k-max", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260716)
    parser.add_argument("--stage-timeout-seconds", type=float, default=12.0)
    parser.add_argument("--solver-time-limit-seconds", type=float, default=10.0)
    parser.add_argument("--enumeration-time-limit-seconds", type=float, default=10.0)
    parser.add_argument("--max-sets", type=int, default=8)
    parser.add_argument(
        "--all-loci-sensitivity-k",
        type=int,
        default=13,
        help="Repeat the sparse case at this k in the full all-loci universe; use 0 to disable.",
    )
    parser.add_argument("--total-wall-limit-seconds", type=float, default=240.0)
    parser.add_argument("--worker", choices=("solve", "enumerate"))
    parser.add_argument("--spec")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.worker:
        if not args.spec:
            raise SystemExit("--spec is required with --worker")
        return worker_main(args.spec, args.worker)
    if args.output_dir is None:
        raise SystemExit("--output-dir is required")
    if args.k_min < 2 or args.k_max < args.k_min:
        raise SystemExit("invalid k range")
    if args.max_sets < 1:
        raise SystemExit("--max-sets must be positive")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    load_before = os.getloadavg()
    rows: list[dict[str, Any]] = []
    for k in range(args.k_min, args.k_max + 1):
        case_kinds = ["sparse_active", "dense_active"]
        if k == args.all_loci_sensitivity_k:
            case_kinds.append("sparse_all_loci_sensitivity")
        for kind_idx, case_kind in enumerate(case_kinds):
            base_kind = "sparse_active" if case_kind == "sparse_all_loci_sensitivity" else case_kind
            case_seed = int(args.seed + 1000 * k + (0 if case_kind == "sparse_all_loci_sensitivity" else kind_idx))
            case = build_case(k, base_kind, case_seed)
            if case_kind == "sparse_all_loci_sensitivity":
                case.update(
                    {
                        "case_kind": case_kind,
                        "source_case_kind": "sparse_active",
                        "universe_mode": "all_loci",
                        "active_k": k,
                        "expected_variables": 1 << k,
                    }
                )
            elapsed = time.perf_counter() - started
            remaining_total = args.total_wall_limit_seconds - elapsed
            if remaining_total <= 1.0:
                case["solve"] = {"status": "TOTAL_WALL_BUDGET_ABSTAIN", "runtime_seconds": 0.0}
                case["enumeration"] = {
                    "status": "NOT_RUN_SOLVE_NOT_OPTIMAL",
                    "complete": False,
                    "n_vertex_sets": 0,
                    "max_sets": args.max_sets,
                    "runtime_seconds": 0.0,
                }
            else:
                worker_spec = dict(case)
                worker_spec.update(
                    {
                        "solver_time_limit_seconds": min(args.solver_time_limit_seconds, max(0.1, remaining_total - 0.5)),
                        "enumeration_time_limit_seconds": args.enumeration_time_limit_seconds,
                        "max_sets": args.max_sets,
                    }
                )
                solve_hard = min(args.stage_timeout_seconds, max(0.1, remaining_total - 0.25))
                case["solve"] = run_stage(worker_spec, "solve", solve_hard)
                if case["solve"].get("status") == "OPTIMAL":
                    remaining_total = args.total_wall_limit_seconds - (time.perf_counter() - started)
                    if remaining_total <= 1.0:
                        case["enumeration"] = {
                            "status": "TOTAL_WALL_BUDGET_ABSTAIN",
                            "complete": False,
                            "n_vertex_sets": 0,
                            "max_sets": args.max_sets,
                            "runtime_seconds": 0.0,
                        }
                    else:
                        worker_spec["enumeration_time_limit_seconds"] = min(
                            args.enumeration_time_limit_seconds, max(0.1, remaining_total - 0.5)
                        )
                        enum_hard = min(args.stage_timeout_seconds, max(0.1, remaining_total - 0.25))
                        case["enumeration"] = run_stage(worker_spec, "enumerate", enum_hard)
                else:
                    case["enumeration"] = {
                        "status": "NOT_RUN_SOLVE_NOT_OPTIMAL",
                        "complete": False,
                        "n_vertex_sets": 0,
                        "max_sets": args.max_sets,
                        "runtime_seconds": 0.0,
                    }
            case["validation_failures"] = validate_row(case)
            case["row_valid"] = not case["validation_failures"]
            rows.append(case)

    cases_path = output_dir / "cases.jsonl"
    with cases_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")

    optimal_rows = [row for row in rows if row["solve"].get("status") == "OPTIMAL"]
    certified_rows = [
        row
        for row in optimal_rows
        if row["solve"].get("objective") is not None
        and row["solve"].get("objective_bound") is not None
        and math.isclose(
            float(row["solve"]["objective"]), float(row["solve"]["objective_bound"]), abs_tol=1e-8
        )
        and (row["solve"].get("mip_gap") is None or float(row["solve"]["mip_gap"]) <= 1e-9)
    ]
    sensitivity: dict[str, Any] = {"status": "NOT_REQUESTED_OR_OUTSIDE_K_RANGE"}
    observed_row = next(
        (
            row
            for row in rows
            if row["k"] == args.all_loci_sensitivity_k and row["case_kind"] == "sparse_active"
        ),
        None,
    )
    all_loci_row = next(
        (row for row in rows if row["case_kind"] == "sparse_all_loci_sensitivity"),
        None,
    )
    if observed_row is not None and all_loci_row is not None:
        observed_objective = observed_row["solve"].get("objective")
        all_loci_objective = all_loci_row["solve"].get("objective")
        sensitivity = {
            "status": (
                "EVALUABLE"
                if observed_row["solve"].get("status") == "OPTIMAL"
                and all_loci_row["solve"].get("status") == "OPTIMAL"
                else "ABSTAIN_NONOPTIMAL_STAGE"
            ),
            "k": args.all_loci_sensitivity_k,
            "observed_alt_active_k": observed_row["evidence_active_k"],
            "observed_alt_variables": observed_row["expected_variables"],
            "all_loci_variables": all_loci_row["expected_variables"],
            "observed_alt_objective": observed_objective,
            "all_loci_objective": all_loci_objective,
            "objective_equal": (
                None
                if observed_objective is None or all_loci_objective is None
                else math.isclose(float(observed_objective), float(all_loci_objective), abs_tol=1e-8)
            ),
            "identity_comparison_complete": bool(observed_row["enumeration"].get("complete"))
            and bool(all_loci_row["enumeration"].get("complete")),
        }
        observed_sets = observed_row["enumeration"].get("vertex_sets", [])
        all_loci_sets = all_loci_row["enumeration"].get("vertex_sets", [])
        if observed_sets and all_loci_sets:
            projection_mask = int(observed_row["observed_alt_mask"])
            observed_ids = {
                vertex_set_digest_local(int(observed_row["k"]), set(vertices)) for vertices in observed_sets
            }
            projected_ids = {
                vertex_set_digest_local(
                    int(all_loci_row["k"]), {int(vertex) & projection_mask for vertex in vertices}
                )
                for vertices in all_loci_sets
            }
            sensitivity.update(
                {
                    "observed_enumerated_n": len(observed_ids),
                    "all_loci_enumerated_n": len(all_loci_sets),
                    "all_loci_projected_unique_n": len(projected_ids),
                    "projected_digest_intersection_n": len(observed_ids & projected_ids),
                    "identity_result_is_lower_bound_if_incomplete": not sensitivity[
                        "identity_comparison_complete"
                    ],
                }
            )
    receipt = {
        "schema": "intersubmod.k_gt12_exact_efficiency_pilot_receipt",
        "schema_version": "1.0.0",
        "scope": "EXPLORATORY_PILOT_NOT_FINAL_VALIDATION",
        "partial_flag": True,
        "serves_goals": ["G3", "G4", "G5"],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(shlex.quote(item) for item in sys.argv),
        "inputs": {
            "benchmark_script": str(pathlib.Path(__file__).resolve()),
            "benchmark_script_sha256": sha256(pathlib.Path(__file__).resolve()),
            "hypercube_exact": str(CORE_PATH.resolve()),
            "hypercube_exact_sha256": sha256(CORE_PATH),
            "prior_pilot_receipt": str(DEFAULT_PRIOR_PILOT.resolve()),
            "prior_pilot_receipt_sha256": sha256(DEFAULT_PRIOR_PILOT),
        },
        "configuration": {
            "k_min": args.k_min,
            "k_max": args.k_max,
            "case_kinds": ["sparse_active", "dense_active"],
            "all_loci_sensitivity_k": args.all_loci_sensitivity_k,
            "base_seed": args.seed,
            "stage_timeout_seconds": args.stage_timeout_seconds,
            "solver_time_limit_seconds": args.solver_time_limit_seconds,
            "enumeration_time_limit_seconds": args.enumeration_time_limit_seconds,
            "max_sets": args.max_sets,
            "total_wall_limit_seconds": args.total_wall_limit_seconds,
            "blas_threads": 1,
        },
        "resource_snapshot": {
            "logical_cpus": os.cpu_count(),
            "load_average_before": list(load_before),
            "load_average_after": list(os.getloadavg()),
        },
        "outputs": {
            "cases_jsonl": str(cases_path),
            "cases_jsonl_sha256": sha256(cases_path),
        },
        "summary": {
            "n_cases": len(rows),
            "n_row_contract_pass": sum(row["row_valid"] for row in rows),
            "n_solve_optimal": len(optimal_rows),
            "n_solve_certified_gap_zero": len(certified_rows),
            "n_solve_abstain_or_limit": len(rows) - len(optimal_rows),
            "n_enumeration_complete": sum(bool(row["enumeration"].get("complete")) for row in rows),
            "n_enumeration_max_sets": sum(
                row["enumeration"].get("status") == "MAX_SETS_REACHED" for row in rows
            ),
            "max_active_k": max(row["active_k"] for row in rows),
            "max_variables": max(row["expected_variables"] for row in rows),
            "wall_seconds": time.perf_counter() - started,
            "within_total_wall_limit": (time.perf_counter() - started) <= args.total_wall_limit_seconds + 1.0,
            "pilot_contract_pass": all(row["row_valid"] for row in rows),
            "scientific_validation": False,
        },
        "all_loci_sensitivity": sensitivity,
    }
    receipt_path = output_dir / "receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "receipt.sha256").write_text(f"{sha256(receipt_path)}  receipt.json\n", encoding="utf-8")
    print(json.dumps(receipt["summary"], sort_keys=True))
    print(f"receipt={receipt_path}")
    return 0 if receipt["summary"]["pilot_contract_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
