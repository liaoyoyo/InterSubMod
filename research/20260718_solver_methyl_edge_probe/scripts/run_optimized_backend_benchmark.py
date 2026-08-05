#!/usr/bin/env python3
"""Reproducible current-vs-optimized bounded benchmark and receipt builder."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.util
import json
import os
import pathlib
import platform
import resource
import statistics
import subprocess
import sys
import time
from typing import Any


REPO = pathlib.Path("/big7_disk/liaoyoyo2001/InterSubMod")
PROBE_PATH = REPO / "research/20260718_solver_methyl_edge_probe/scripts/solver_probe.py"
OPTIMIZED_PATH = (
    REPO
    / "research/20260718_solver_methyl_edge_probe/scripts/optimized_hypercube_backend.py"
)
CURRENT_SOLVER_PATH = (
    REPO
    / "research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/hypercube_exact.py"
)
FIXTURE_PATH = (
    REPO
    / "research/20260718_solver_methyl_edge_probe/tests/fixtures/real_units.json"
)
PLAN_PATH = (
    REPO
    / "research/20260718_solver_methyl_edge_probe/20260718_Hypercube邊與subcube改良研究計畫_01.md"
)
TEST_PATH = (
    REPO
    / "research/20260718_solver_methyl_edge_probe/tests/test_optimized_hypercube_backend.py"
)
ORACLE_SCRIPT_PATH = (
    REPO
    / "research/20260718_solver_methyl_edge_probe/scripts/verify_optimized_backend_oracles.py"
)
ORACLE_RECEIPT_PATH = (
    REPO
    / "research/20260718_solver_methyl_edge_probe/results/optimized_backend_benchmark_v1/oracle_receipt.json"
)
CASES = ("AAAA", "H2009_M31", "COLO829_M31")
BACKENDS = ("current_scipy_milp", "optimized_dp_bitset_bnb")
THREAD_ENV = {
    "PYTHONHASHSEED": "0",
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(path: pathlib.Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def canonical_family_digest(k: int, vertex_sets, probe) -> tuple[str, list[str]]:
    identifiers = sorted(probe.vertex_set_digest(k, values) for values in vertex_sets)
    payload = {
        "schema": "intersubmod.optimal_vertex_family_digest.v1",
        "k": int(k),
        "vertex_set_ids": identifiers,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest(), identifiers


def load_case(case_id: str, probe):
    if case_id == "AAAA":
        full = ["AAAA"]
        partial: list[str] = []
        k = 4
        expected_h = 3
        max_sets = 256
    else:
        document = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
        fixture = next(row for row in document["units"] if row["case_id"] == case_id)
        full = list(fixture["full_counts"])
        partial = list(fixture["partial_counts"])
        k = int(fixture["k"])
        expected_h = int(fixture["expected_h"])
        max_sets = int(fixture["max_sets"])
    return (
        probe.build_instance(full, partial, k),
        full,
        partial,
        k,
        expected_h,
        max_sets,
    )


def worker(args: argparse.Namespace) -> int:
    probe = load_module(PROBE_PATH, "solver_probe")
    instance, full, partial, k, expected_h, max_sets = load_case(
        args.worker_case, probe
    )
    harness_started_wall = time.perf_counter()
    harness_started_cpu = time.process_time()

    if args.worker_backend == "current_scipy_milp":
        result = probe.current_scipy_enumerate(
            str(CURRENT_SOLVER_PATH),
            full,
            partial,
            k,
            max_sets=max_sets,
            time_limit_seconds=args.time_limit,
        )
        status = str(result["status"])
        complete = bool(result["complete"])
        objective_certified = complete
        family_complete = complete
        h_star = result.get("objective")
        vertex_sets = [tuple(values) for values in result.get("vertex_sets", [])]
        solver_elapsed = float(result["elapsed_seconds"])
        diagnostics = {
            "model_builds": result.get("model_builds"),
            "solve_calls": result.get("solve_calls"),
        }
        q = None
        prepare_elapsed = None
        objective_dp_elapsed = None
    else:
        optimized = load_module(OPTIMIZED_PATH, "optimized_hypercube_backend")
        solver_started = time.perf_counter()
        result = optimized.solve_with_certificates(
            instance,
            q_max=args.q_max,
            time_limit_seconds=args.time_limit,
            max_sets=max_sets,
        )
        solver_elapsed = time.perf_counter() - solver_started
        objective = result["objective_certificate"]
        family = result["family_certificate"]
        status = str(family["status"])
        complete = bool(family["family_complete"])
        objective_certified = bool(objective["objective_certified"])
        family_complete = bool(family["family_complete"])
        h_star = objective["objective"]
        vertex_sets = [tuple(values) for values in family["vertex_sets"]]
        diagnostics = dict(family["stats"])
        diagnostics.update(
            {
                "model_builds": 0,
                "solve_calls": 0,
                "downward_closure_vertices": result["reduction"][
                    "downward_closure_vertices"
                ],
                "input_vertices": result["reduction"]["input_vertices"],
                "input_groups": result["reduction"]["input_groups"],
                "active_groups": result["reduction"]["active_groups"],
            }
        )
        q = result["q"]
        prepare_elapsed = result["prepare_elapsed_seconds"]
        objective_dp_elapsed = objective["elapsed_seconds"]

    harness_cpu = time.process_time() - harness_started_cpu
    harness_wall = time.perf_counter() - harness_started_wall
    digest, identifiers = canonical_family_digest(k, vertex_sets, probe)
    all_valid = all(check["pass"] for check in (
        probe.check_selected(instance, values) for values in vertex_sets
    ))

    oracle_digest = None
    if args.worker_case == "AAAA":
        oracle = probe.brute_force_optimal(instance)
        if not oracle["complete"]:
            raise RuntimeError("AAAA brute-force oracle did not complete")
        oracle_digest, _ = canonical_family_digest(k, oracle["vertex_sets"], probe)

    if not objective_certified or h_star != expected_h:
        raise RuntimeError(
            f"objective certificate failed: certified={objective_certified} "
            f"observed={h_star} expected={expected_h}"
        )
    if not complete or not family_complete:
        raise RuntimeError(
            f"incomplete family refused by benchmark: status={status}"
        )
    if not all_valid:
        raise RuntimeError("backend emitted an independently invalid vertex set")
    if oracle_digest is not None and digest != oracle_digest:
        raise RuntimeError("AAAA family differs from brute-force oracle")

    receipt = {
        "schema": "intersubmod.optimized_backend_benchmark.worker.v1",
        "case_id": args.worker_case,
        "backend": args.worker_backend,
        "status": status,
        "complete": complete,
        "objective_certified": objective_certified,
        "family_complete": family_complete,
        "ranking_allowed": family_complete,
        "capped": not family_complete,
        "h_star": h_star,
        "optimal_family_count": len(identifiers),
        "canonical_family_digest": digest,
        "canonical_digest_schema": "intersubmod.optimal_vertex_family_digest.v1",
        "oracle_digest": oracle_digest,
        "all_sets_independently_valid": all_valid,
        "k": k,
        "effective_m": probe.popcount(instance.universe_mask),
        "n_vertices": len(instance.vertices),
        "q": q,
        "max_sets": max_sets,
        "total_unit_deadline_seconds": args.time_limit,
        "solver_elapsed_seconds": solver_elapsed,
        "prepare_elapsed_seconds": prepare_elapsed,
        "objective_dp_elapsed_seconds": objective_dp_elapsed,
        "harness_wall_seconds": harness_wall,
        "harness_cpu_seconds": harness_cpu,
        "ru_maxrss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "diagnostics": diagnostics,
        "source_sha256": {
            "probe": sha256_file(PROBE_PATH),
            "optimized_backend": sha256_file(OPTIMIZED_PATH),
            "current_solver": sha256_file(CURRENT_SOLVER_PATH),
            "fixtures": sha256_file(FIXTURE_PATH),
        },
    }
    args.worker_output.parent.mkdir(parents=True, exist_ok=True)
    args.worker_output.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))
    return 0


def percentile(values: list[float], quantile: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * quantile
    lower = int(math_floor(position))
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def math_floor(value: float) -> int:
    return int(value // 1)


def summarize_metric(rows: list[dict[str, Any]], field: str) -> dict[str, float]:
    values = [float(row[field]) for row in rows]
    return {
        "n": len(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
        "p95": percentile(values, 0.95),
        "mean": statistics.mean(values),
        "cv": (
            statistics.pstdev(values) / statistics.mean(values)
            if len(values) > 1 and statistics.mean(values)
            else 0.0
        ),
    }


def read_imported_runs(
    directory: pathlib.Path,
    repeats: int,
    expected_backend: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case_id in CASES:
        for repeat in range(1, repeats + 1):
            path = directory / f"{case_id}.r{repeat}.json"
            if not path.exists():
                raise FileNotFoundError(f"missing imported baseline repeat: {path}")
            row = json.loads(path.read_text(encoding="utf-8"))
            if row["backend"] != expected_backend:
                raise ValueError(
                    f"imported backend mismatch in {path}: "
                    f"{row['backend']} != {expected_backend}"
                )
            row["imported_baseline_path"] = str(path)
            rows.append(row)
    return rows


def run_repeats(
    script: pathlib.Path,
    output_dir: pathlib.Path,
    backend: str,
    repeats: int,
    warmups: int,
    time_limit: float,
    q_max: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    env = os.environ.copy()
    env.update(THREAD_ENV)
    for case_id in CASES:
        for index in range(warmups + repeats):
            label = f"warmup{index + 1}" if index < warmups else f"r{index - warmups + 1}"
            output = output_dir / backend / f"{case_id}.{label}.json"
            stdout_path = output.with_suffix(".stdout")
            stderr_path = output.with_suffix(".stderr")
            command = [
                sys.executable,
                str(script),
                "--worker-backend",
                backend,
                "--worker-case",
                case_id,
                "--worker-output",
                str(output),
                "--time-limit",
                str(time_limit),
                "--q-max",
                str(q_max),
            ]
            process_started = time.perf_counter()
            completed = subprocess.run(
                command,
                cwd=str(REPO),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            process_elapsed = time.perf_counter() - process_started
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_path.write_text(completed.stdout, encoding="utf-8")
            stderr_path.write_text(completed.stderr, encoding="utf-8")
            if completed.returncode != 0:
                raise RuntimeError(
                    f"benchmark worker failed ({backend}/{case_id}/{label}): "
                    f"{completed.stderr[-2000:]}"
                )
            row = json.loads(output.read_text(encoding="utf-8"))
            row["repeat_label"] = label
            row["process_elapsed_seconds"] = process_elapsed
            row["stdout_path"] = str(stdout_path)
            row["stderr_path"] = str(stderr_path)
            if index >= warmups:
                rows.append(row)
    return rows


def package_environment() -> dict[str, Any]:
    versions: dict[str, str | None] = {}
    for name in ("numpy", "scipy"):
        try:
            module = importlib.import_module(name)
            versions[name] = str(module.__version__)
        except ModuleNotFoundError:
            versions[name] = None
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "packages": versions,
        "thread_environment": THREAD_ENV,
    }


def suite(args: argparse.Namespace) -> int:
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    script = pathlib.Path(__file__).resolve()

    if args.import_current_dir:
        current_rows = read_imported_runs(
            args.import_current_dir.resolve(),
            args.repeats,
            "current_scipy_milp",
        )
        current_provenance = {
            "mode": "IMPORTED_INDEPENDENT_BASELINE",
            "path": str(args.import_current_dir.resolve()),
        }
    else:
        current_rows = run_repeats(
            script,
            output_dir,
            "current_scipy_milp",
            args.repeats,
            args.warmups,
            args.time_limit,
            args.q_max,
        )
        current_provenance = {"mode": "RUN_BY_THIS_HARNESS"}
    if args.import_optimized_dir:
        optimized_rows = read_imported_runs(
            args.import_optimized_dir.resolve(),
            args.repeats,
            "optimized_dp_bitset_bnb",
        )
        optimized_provenance = {
            "mode": "IMPORTED_THIS_HARNESS_RAW_RUNS",
            "path": str(args.import_optimized_dir.resolve()),
        }
    else:
        optimized_rows = run_repeats(
            script,
            output_dir,
            "optimized_dp_bitset_bnb",
            args.repeats,
            args.warmups,
            args.time_limit,
            args.q_max,
        )
        optimized_provenance = {"mode": "RUN_BY_THIS_HARNESS"}

    rows_by_backend_case: dict[str, dict[str, list[dict[str, Any]]]] = {
        backend: {case_id: [] for case_id in CASES}
        for backend in BACKENDS
    }
    for row in current_rows + optimized_rows:
        rows_by_backend_case[row["backend"]][row["case_id"]].append(row)

    comparisons: list[dict[str, Any]] = []
    all_exact = True
    all_repeat_stable = True
    for case_id in CASES:
        current_case = rows_by_backend_case["current_scipy_milp"][case_id]
        optimized_case = rows_by_backend_case["optimized_dp_bitset_bnb"][case_id]
        current_digests = {row["canonical_family_digest"] for row in current_case}
        optimized_digests = {row["canonical_family_digest"] for row in optimized_case}
        current_h = {row["h_star"] for row in current_case}
        optimized_h = {row["h_star"] for row in optimized_case}
        current_counts = {row["optimal_family_count"] for row in current_case}
        optimized_counts = {row["optimal_family_count"] for row in optimized_case}
        repeat_stable = all(
            len(values) == 1
            for values in (
                current_digests,
                optimized_digests,
                current_h,
                optimized_h,
                current_counts,
                optimized_counts,
            )
        )
        exact = (
            repeat_stable
            and current_digests == optimized_digests
            and current_h == optimized_h
            and current_counts == optimized_counts
            and all(row["complete"] for row in current_case + optimized_case)
        )
        all_repeat_stable &= repeat_stable
        all_exact &= exact
        current_solver = summarize_metric(current_case, "solver_elapsed_seconds")
        optimized_solver = summarize_metric(optimized_case, "solver_elapsed_seconds")
        current_wall = summarize_metric(current_case, "harness_wall_seconds")
        optimized_wall = summarize_metric(optimized_case, "harness_wall_seconds")
        current_rss = summarize_metric(current_case, "ru_maxrss_kib")
        optimized_rss = summarize_metric(optimized_case, "ru_maxrss_kib")
        comparisons.append(
            {
                "case_id": case_id,
                "exact_match": exact,
                "repeat_stable": repeat_stable,
                "h_star": next(iter(current_h)) if len(current_h) == 1 else None,
                "optimal_family_count": (
                    next(iter(current_counts)) if len(current_counts) == 1 else None
                ),
                "canonical_family_digest": (
                    next(iter(current_digests)) if len(current_digests) == 1 else None
                ),
                "current": {
                    "solver_elapsed_seconds": current_solver,
                    "harness_wall_seconds": current_wall,
                    "ru_maxrss_kib": current_rss,
                },
                "optimized": {
                    "solver_elapsed_seconds": optimized_solver,
                    "harness_wall_seconds": optimized_wall,
                    "ru_maxrss_kib": optimized_rss,
                },
                "median_speedup_solver": (
                    current_solver["median"] / optimized_solver["median"]
                ),
                "median_speedup_harness": (
                    current_wall["median"] / optimized_wall["median"]
                ),
                "median_rss_ratio_optimized_over_current": (
                    optimized_rss["median"] / current_rss["median"]
                ),
            }
        )

    receipt = {
        "schema": "intersubmod.optimized_backend_benchmark.suite.v1",
        "created_at_epoch_seconds": time.time(),
        "task_type": "A_EXPLORATORY_PILOT",
        "scope": {
            "partial": True,
            "tier": "BOUNDED_PROBE",
            "production_claim_allowed": False,
            "cases": list(CASES),
            "repeats": args.repeats,
            "warmups_for_new_runs": args.warmups,
            "claim_ceiling": (
                "Exact equivalence and local timing for AAAA plus two M31 fixtures; "
                "not 33-tail, full-sample, canonical, or production validation."
            ),
        },
        "source_sha256": {
            "probe": sha256_file(PROBE_PATH),
            "optimized_backend": sha256_file(OPTIMIZED_PATH),
            "current_solver": sha256_file(CURRENT_SOLVER_PATH),
            "fixtures": sha256_file(FIXTURE_PATH),
            "benchmark_harness": sha256_file(script),
            "optimized_tests": sha256_file(TEST_PATH),
            "oracle_script": sha256_file(ORACLE_SCRIPT_PATH),
            "oracle_receipt": sha256_file(ORACLE_RECEIPT_PATH),
        },
        "planning_context_snapshot": {
            "path": str(PLAN_PATH),
            "sha256_at_receipt_build": sha256_file(PLAN_PATH),
            "mutable_context_not_executable_provenance": True,
        },
        "environment": package_environment(),
        "comparison_contract": {
            "same_cases": True,
            "same_fixture_bytes": True,
            "same_objective": True,
            "same_max_sets": True,
            "same_thread_environment": True,
            "declared_time_limit_seconds": args.time_limit,
            "time_limit_semantics_same": False,
            "time_limit_semantics": {
                "current_scipy_milp": "per MILP solve",
                "optimized_dp_bitset_bnb": "per unit across prepare+subset DP+B&B",
            },
            "completed_before_any_limit": all(
                row["complete"] for row in current_rows + optimized_rows
            ),
            "timing_comparable_with_caveat": all_exact,
            "current_provenance": current_provenance,
            "optimized_provenance": optimized_provenance,
        },
        "comparisons": comparisons,
        "certificates": {
            "objective": {
                "status": "PROVEN_OPTIMAL" if all_exact else "FAIL_MISMATCH",
                "mismatch_count": sum(not row["exact_match"] for row in comparisons),
            },
            "candidate_family": {
                "status": "COMPLETE" if all_exact else "FAIL_MISMATCH",
                "all_repeats_stable": all_repeat_stable,
                "incomplete_ranked_count": 0,
            },
            "edge": {
                "status": "TESTED_SEPARATELY_BY_PARENT_ORACLE",
                "real_edge_claim_allowed": False,
            },
            "biology": {
                "status": "NOT_TESTED",
                "ancestry_claim_allowed": False,
            },
        },
        "fail_closed": {
            "incomplete_ranked_count": 0,
            "winner_emitted_on_incomplete": False,
            "cap_and_deadline_negative_tests_required": True,
        },
        "promotion_gate": {
            "overall": "PASS_FOR_BOUNDED_DUAL_PILOT" if all_exact else "FAIL_CLOSED",
            "canonical_or_production_promotion_allowed": False,
            "authorized_next_action": (
                "Run the frozen 33-tail plus H2009/H1437 stress panel."
                if all_exact
                else "Diagnose exact mismatch before any expansion."
            ),
            "prohibited_actions": [
                "Replace the canonical solver",
                "Start a production/full run from this bounded receipt",
                "Claim biological ancestry improvement",
            ],
        },
        "raw_runs": current_rows + optimized_rows,
        "limitations": [
            "Only three deterministic cases were timed.",
            "Imported current baseline and optimized measurements were sequential but not interleaved.",
            "The current source files are untracked; byte hashes, not git HEAD, bind identity.",
            "Subset DP certifies the objective only; bitset B&B certifies the complete family.",
        ],
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "verdict": receipt["promotion_gate"]["overall"],
                "all_exact": all_exact,
                "comparisons": [
                    {
                        "case_id": row["case_id"],
                        "speedup": row["median_speedup_solver"],
                        "rss_ratio": row[
                            "median_rss_ratio_optimized_over_current"
                        ],
                    }
                    for row in comparisons
                ],
            },
            sort_keys=True,
        )
    )
    return 0 if all_exact else 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-backend", choices=BACKENDS)
    parser.add_argument("--worker-case", choices=CASES)
    parser.add_argument("--worker-output", type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path)
    parser.add_argument("--receipt", type=pathlib.Path)
    parser.add_argument("--import-current-dir", type=pathlib.Path)
    parser.add_argument("--import-optimized-dir", type=pathlib.Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--time-limit", type=float, default=30.0)
    parser.add_argument("--q-max", type=int, default=8)
    args = parser.parse_args()
    if args.worker_backend:
        if not args.worker_case or not args.worker_output:
            parser.error("worker mode requires --worker-case and --worker-output")
    elif not args.output_dir or not args.receipt:
        parser.error("suite mode requires --output-dir and --receipt")
    if args.repeats < 1 or args.warmups < 0:
        parser.error("repeats must be >=1 and warmups must be >=0")
    return args


def main() -> int:
    args = parse_args()
    return worker(args) if args.worker_backend else suite(args)


if __name__ == "__main__":
    raise SystemExit(main())
