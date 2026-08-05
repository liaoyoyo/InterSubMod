#!/usr/bin/env python3
"""Run a fail-closed current-vs-optimized structural solver stress panel.

Each structural key is executed in an isolated process.  A unit-level deadline
covers the full enumeration loop, incomplete families are never rankable, and
only a compressed optimal-family digest is persisted.  Parent trees are never
materialized.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import resource
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Mapping, Sequence


SCRIPT_PATH = pathlib.Path(__file__).resolve()
THREAD_ENV = {
    "PYTHONHASHSEED": "0",
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}
BACKENDS = ("current", "optimized")
AUTHORITY_POINTER_NAME = "AUTHORITATIVE_MANIFEST.json"


class StressPanelError(RuntimeError):
    """Raised when a source, manifest, solver, or acceptance gate fails."""


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def semantic_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def manifest_content_sha256(document: Mapping[str, Any]) -> str:
    payload = copy.deepcopy(dict(document))
    payload.setdefault("integrity", {}).pop("manifest_content_sha256", None)
    return semantic_sha256(payload)


def verify_byte_sidecar(path: pathlib.Path) -> str:
    sidecar = path.with_name(path.name + ".sha256")
    if not sidecar.is_file():
        raise StressPanelError(f"SHA-256 sidecar is missing: {sidecar}")
    fields = sidecar.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2 or fields[1] != path.name:
        raise StressPanelError(f"malformed SHA-256 sidecar: {sidecar}")
    observed = sha256_file(path)
    if fields[0] != observed:
        raise StressPanelError(
            f"byte hash mismatch for {path}: {observed} != {fields[0]}"
        )
    return observed


def load_module(path: pathlib.Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise StressPanelError(f"cannot import source: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def verify_manifest(
    manifest_path: pathlib.Path,
    *,
    require_sidecar: bool = True,
    verify_inputs: bool = True,
    expected_authority_status: str | None = None,
    expected_runner_sha256: str | None = None,
) -> Dict[str, Any]:
    if not manifest_path.is_file():
        raise StressPanelError(f"manifest is missing: {manifest_path}")
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    if document.get("schema_name") != "intersubmod.solver_stress_panel_manifest":
        raise StressPanelError("unexpected stress-panel manifest schema")
    expected_content = document.get("integrity", {}).get(
        "manifest_content_sha256"
    )
    observed_content = manifest_content_sha256(document)
    if expected_content != observed_content:
        raise StressPanelError(
            f"manifest content hash mismatch: {observed_content} != {expected_content}"
        )
    if not all(document.get("checks", {}).values()):
        raise StressPanelError("manifest checks are not all pass")

    if require_sidecar:
        verify_byte_sidecar(manifest_path)

    authority_status = document.get("authority", {}).get("status")
    if expected_authority_status is not None:
        if authority_status != expected_authority_status:
            raise StressPanelError(
                "manifest authority status mismatch: "
                f"{authority_status!r} != {expected_authority_status!r}"
            )
        if not str(authority_status).startswith("AUTHORITATIVE_"):
            raise StressPanelError(
                f"manifest is not authoritative: {authority_status!r}"
            )

    categories = ["source_files"]
    if verify_inputs:
        categories.insert(0, "input_files")
    for category in categories:
        for name, record in document["bindings"][category].items():
            path = pathlib.Path(record["path"])
            if not path.is_file():
                raise StressPanelError(f"bound {category}/{name} missing: {path}")
            observed = sha256_file(path)
            if observed != record["sha256"]:
                raise StressPanelError(
                    f"bound {category}/{name} hash mismatch: "
                    f"{observed} != {record['sha256']}"
                )
    bound_runner = document["bindings"]["source_files"].get("stress_runner")
    if bound_runner is None:
        raise StressPanelError("manifest does not bind stress_runner")
    current_runner_sha = sha256_file(SCRIPT_PATH)
    if bound_runner["sha256"] != current_runner_sha:
        raise StressPanelError(
            "manifest stress_runner binding does not match executing runner: "
            f"{bound_runner['sha256']} != {current_runner_sha}"
        )
    if (
        expected_runner_sha256 is not None
        and current_runner_sha != expected_runner_sha256
    ):
        raise StressPanelError(
            "authority pointer runner hash does not match executing runner: "
            f"{expected_runner_sha256} != {current_runner_sha}"
        )
    return document


def resolve_authority_pointer(
    pointer_path: pathlib.Path,
    *,
    verify_inputs: bool = True,
) -> tuple[pathlib.Path, Dict[str, Any], Dict[str, Any]]:
    """Resolve and cryptographically bind pointer, manifest, and runner."""
    if not pointer_path.is_file():
        raise StressPanelError(f"authority pointer is missing: {pointer_path}")
    pointer_file_sha = verify_byte_sidecar(pointer_path)
    pointer = json.loads(pointer_path.read_text(encoding="utf-8"))
    if pointer.get("schema") != "intersubmod.solver_stress_panel.authority_pointer.v1":
        raise StressPanelError("unexpected authority-pointer schema")
    status = str(pointer.get("status", ""))
    if not status.startswith("AUTHORITATIVE_"):
        raise StressPanelError(f"authority pointer is not active: {status!r}")
    manifest_path = (
        pointer_path.parent / str(pointer["authoritative_manifest"])
    ).resolve()
    expected_manifest_file_sha = pointer.get(
        "authoritative_manifest_file_sha256"
    )
    if expected_manifest_file_sha is None:
        raise StressPanelError(
            "authority pointer lacks authoritative_manifest_file_sha256"
        )
    observed_manifest_file_sha = sha256_file(manifest_path)
    if observed_manifest_file_sha != expected_manifest_file_sha:
        raise StressPanelError(
            "authority manifest byte hash mismatch: "
            f"{observed_manifest_file_sha} != {expected_manifest_file_sha}"
        )
    runner_sha = str(pointer.get("authoritative_runner_sha256", ""))
    manifest = verify_manifest(
        manifest_path,
        verify_inputs=verify_inputs,
        expected_authority_status=status,
        expected_runner_sha256=runner_sha,
    )
    if (
        manifest["integrity"]["manifest_content_sha256"]
        != pointer.get("authoritative_manifest_content_sha256")
    ):
        raise StressPanelError("authority pointer manifest content hash mismatch")
    return manifest_path, manifest, {
        "pointer_path": str(pointer_path.resolve()),
        "pointer_file_sha256": pointer_file_sha,
        "status": status,
        "manifest_file_sha256": observed_manifest_file_sha,
        "runner_sha256": runner_sha,
    }


def resolve_direct_manifest(
    manifest_path: pathlib.Path,
    *,
    authority_status: str | None,
    verify_inputs: bool = True,
) -> tuple[pathlib.Path, Dict[str, Any], Dict[str, Any]]:
    """Allow direct manifests only with explicit, current authority status."""
    if authority_status is None:
        raise StressPanelError(
            "direct --manifest use requires --authority-status; "
            "--authority-pointer is preferred"
        )
    pointer_path = manifest_path.parent.parent / AUTHORITY_POINTER_NAME
    if pointer_path.is_file():
        authoritative_path, manifest, authority = resolve_authority_pointer(
            pointer_path,
            verify_inputs=verify_inputs,
        )
        if authoritative_path != manifest_path.resolve():
            raise StressPanelError(
                f"direct manifest is superseded by {authoritative_path}"
            )
        if authority["status"] != authority_status:
            raise StressPanelError(
                "explicit authority status does not match current pointer"
            )
        return authoritative_path, manifest, authority
    manifest = verify_manifest(
        manifest_path,
        verify_inputs=verify_inputs,
        expected_authority_status=authority_status,
    )
    return manifest_path.resolve(), manifest, {
        "pointer_path": None,
        "pointer_file_sha256": None,
        "status": authority_status,
        "manifest_file_sha256": sha256_file(manifest_path),
        "runner_sha256": sha256_file(SCRIPT_PATH),
    }


def family_digest(probe, k: int, vertex_sets: Iterable[Iterable[int]]) -> str:
    identifiers = sorted(
        probe.vertex_set_digest(k, values) for values in vertex_sets
    )
    return semantic_sha256(
        {
            "schema": "intersubmod.optimal_vertex_family_digest.v1",
            "k": int(k),
            "vertex_set_ids": identifiers,
        }
    )


def _remaining(deadline: float) -> float:
    return max(0.0, deadline - time.perf_counter())


def current_enumerate_with_total_deadline(
    current,
    *,
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
    k: int,
    time_limit_seconds: float,
    max_sets: int | None,
) -> Dict[str, Any]:
    """Enumerate with one wall-clock budget across all rebuilt MILP models."""
    if not math.isfinite(time_limit_seconds) or time_limit_seconds <= 0:
        raise ValueError("time_limit_seconds must be finite and positive")
    if max_sets is not None and max_sets < 1:
        raise ValueError("max_sets must be positive or None")
    started = time.perf_counter()
    deadline = started + float(time_limit_seconds)
    model_builds = 0
    sets: list[tuple[int, ...]] = []
    objective = None
    objective_certified = False
    status = "NO_RESULT"

    def solve(*, objective_value=None, excluded=()):
        nonlocal model_builds
        remaining = _remaining(deadline)
        if remaining <= 0:
            return None
        model_builds += 1
        return current.solve_min_hidden(
            full_patterns,
            partial_patterns,
            int(k),
            universe_mode="observed_alt",
            time_limit_seconds=max(1e-9, remaining),
            objective_value=objective_value,
            excluded_vertex_sets=excluded,
        )

    first = solve()
    if first is None:
        status = "CANDIDATE_SET_INCOMPLETE_DEADLINE"
    elif _remaining(deadline) <= 0:
        status = "CANDIDATE_SET_INCOMPLETE_DEADLINE"
    elif first.status != "OPTIMAL" or first.objective is None:
        status = f"FIRST_{first.status}"
    else:
        objective = int(round(float(first.objective)))
        objective_certified = True
        sets = [tuple(sorted(first.selected_vertices))]
        if objective == 0:
            status = "CANDIDATE_SET_COMPLETE"
        else:
            while True:
                if max_sets is not None and len(sets) >= max_sets:
                    status = "CANDIDATE_SET_INCOMPLETE_MAX_SETS"
                    break
                next_solution = solve(
                    objective_value=objective,
                    excluded=sets,
                )
                if next_solution is None or _remaining(deadline) <= 0:
                    status = "CANDIDATE_SET_INCOMPLETE_DEADLINE"
                    break
                if next_solution.status == "INFEASIBLE":
                    status = "CANDIDATE_SET_COMPLETE"
                    break
                if next_solution.status != "OPTIMAL":
                    status = f"CANDIDATE_SET_INCOMPLETE_{next_solution.status}"
                    break
                candidate = tuple(sorted(next_solution.selected_vertices))
                if candidate in sets:
                    raise StressPanelError(
                        "current solver repeated an excluded vertex set"
                    )
                sets.append(candidate)

    elapsed = time.perf_counter() - started
    complete = status == "CANDIDATE_SET_COMPLETE"
    return {
        "backend_detail": "frozen_current_scipy_rebuild_loop",
        "status": status,
        "complete": complete,
        "objective_certified": objective_certified,
        "family_complete": complete,
        "objective": objective,
        "vertex_sets": sets,
        "solver_elapsed_seconds": elapsed,
        "model_builds": model_builds,
        "total_unit_deadline_exceeded": elapsed >= time_limit_seconds,
    }


def _load_case_instance(
    manifest: Mapping[str, Any],
    case: Mapping[str, Any],
):
    source = manifest["bindings"]["source_files"]
    probe = load_module(pathlib.Path(source["solver_probe"]["path"]), "solver_probe")
    structural = case["structural_input"]
    observed_key = semantic_sha256(structural)
    if observed_key != case["structural_key_sha256"]:
        raise StressPanelError(
            f"structural key mismatch for {case['case_id']}: "
            f"{observed_key} != {case['structural_key_sha256']}"
        )
    instance = probe.build_instance(
        list(structural["full_patterns"]),
        list(structural["partial_patterns"]),
        int(structural["k"]),
        universe_mode="observed_alt",
    )
    if instance.universe_mask != int(
        structural["structural_alt_universe_mask"]
    ):
        raise StressPanelError(f"universe-mask mismatch for {case['case_id']}")
    return probe, instance


def run_backend_case(
    manifest: Mapping[str, Any],
    case: Mapping[str, Any],
    *,
    backend: str,
    deadline_seconds: float,
    max_sets: int | None,
    q_max: int,
) -> Dict[str, Any]:
    if backend not in BACKENDS:
        raise ValueError(f"unsupported backend: {backend}")
    probe, instance = _load_case_instance(manifest, case)
    source = manifest["bindings"]["source_files"]

    if backend == "current":
        current = load_module(
            pathlib.Path(source["frozen_current_solver"]["path"]),
            "_stress_panel_frozen_current",
        )
        result = current_enumerate_with_total_deadline(
            current,
            full_patterns=list(instance.full_patterns),
            partial_patterns=list(instance.partial_patterns),
            k=instance.k,
            time_limit_seconds=deadline_seconds,
            max_sets=max_sets,
        )
    else:
        optimized = load_module(
            pathlib.Path(source["optimized_backend"]["path"]),
            "optimized_hypercube_backend",
        )
        result_raw = optimized.solve_with_certificates(
            instance,
            q_max=q_max,
            time_limit_seconds=deadline_seconds,
            max_sets=max_sets,
        )
        family = result_raw["family_certificate"]
        objective = result_raw["objective_certificate"]
        result = {
            "backend_detail": family["backend"],
            "status": family["status"],
            "complete": bool(family["family_complete"]),
            # The family certificate is authoritative.  It carries either the
            # subset-DP target or the objective proved by a complete B&B
            # traversal.  The DP record remains diagnostic because it may
            # deliberately route away on q/resource gates.
            "objective_certified": bool(family["objective_certified"]),
            "family_complete": bool(family["family_complete"]),
            "objective": family["objective"],
            "vertex_sets": [
                tuple(values) for values in family.get("vertex_sets", [])
            ],
            "solver_elapsed_seconds": float(
                result_raw["total_elapsed_seconds"]
            ),
            "model_builds": 0,
            "total_unit_deadline_exceeded": bool(
                result_raw["total_unit_deadline_exceeded"]
            ),
            "q": result_raw["q"],
            "reduction": result_raw["reduction"],
            "objective_status": family["status"],
            "dp_objective_status": objective["status"],
            "dp_objective_certified": bool(
                objective["objective_certified"]
            ),
            "dp_objective": objective["objective"],
            "family_certificate_stats": copy.deepcopy(
                family.get("stats", {})
            ),
            "objective_dp_stats": copy.deepcopy(
                objective.get("stats", {})
            ),
        }

    vertex_sets = [tuple(values) for values in result.pop("vertex_sets")]
    independent_checks = [
        probe.check_selected(instance, values) for values in vertex_sets
    ]
    if not all(check["pass"] for check in independent_checks):
        raise StressPanelError(
            f"backend emitted invalid vertex set for {case['case_id']}"
        )
    complete = bool(result["family_complete"])
    ranking_allowed = complete
    incomplete_ranked = (not complete) and ranking_allowed
    if incomplete_ranked:
        raise StressPanelError("incomplete family was marked rankable")
    digest = (
        family_digest(probe, instance.k, vertex_sets) if vertex_sets else None
    )
    factorial = case.get("factorial_oracle")
    factorial_pass = None
    if factorial is not None and complete:
        factorial_pass = len(vertex_sets) == int(
            factorial["expected_optimal_vertex_sets"]
        )

    return {
        "schema": "intersubmod.solver_stress_panel.worker.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "case_id": case["case_id"],
        "structural_key_sha256": case["structural_key_sha256"],
        "backend": backend,
        **result,
        "optimal_family_count": len(vertex_sets),
        "optimal_family_digest": digest,
        "ranking_allowed": ranking_allowed,
        "incomplete_ranked": incomplete_ranked,
        "all_sets_independently_valid": all(
            check["pass"] for check in independent_checks
        ),
        "parent_trees_materialized": False,
        "max_sets": max_sets,
        "total_unit_deadline_seconds": float(deadline_seconds),
        "k": instance.k,
        "effective_m": probe.popcount(instance.universe_mask),
        "n_vertices": len(instance.vertices),
        "factorial_oracle": factorial,
        "factorial_oracle_pass": factorial_pass,
        "ru_maxrss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "manifest_content_sha256": manifest["integrity"][
            "manifest_content_sha256"
        ],
    }


def percentile(values: Sequence[float], quantile: float) -> float | None:
    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * quantile
    lower = int(math.floor(position))
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def _case_index(manifest: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {case["case_id"]: case for case in manifest["cases"]}


def select_case_ids(
    manifest: Mapping[str, Any],
    *,
    panel: str,
    selection: str,
) -> list[str]:
    if panel == "both":
        selected = set(manifest["panels"]["primary"]["case_ids"])
        selected.update(manifest["panels"]["sensitivity"]["case_ids"])
    else:
        selected = set(manifest["panels"][panel]["case_ids"])
    cases = _case_index(manifest)

    def eligible(case: Mapping[str, Any]) -> bool:
        reasons = set(case["selection_reasons"])
        if selection == "all":
            return True
        if selection == "controls":
            return "COMPLETE_WITH_VERTEX_FAMILY_GE_100" in reasons
        if selection == "tails":
            return "33_INCOMPLETE_CANDIDATE_GENERATION" in reasons
        if selection == "factorial":
            return case.get("factorial_oracle") is not None
        raise ValueError(selection)

    return sorted(case_id for case_id in selected if eligible(cases[case_id]))


def build_summary(
    manifest: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
    *,
    selected_case_ids: Sequence[str],
    requested_backends: Sequence[str],
    deadline_seconds: float,
    max_sets: int | None,
    expected_repeats: int | None = None,
) -> Dict[str, Any]:
    case_index = _case_index(manifest)
    control_cases = {
        case_id
        for case_id in selected_case_ids
        if "COMPLETE_WITH_VERTEX_FAMILY_GE_100"
        in case_index[case_id]["selection_reasons"]
    }
    by_key = {
        (row["case_id"], int(row["repeat"]), row["backend"]): row
        for row in rows
    }
    observed_repeat_values = sorted({int(row["repeat"]) for row in rows})
    if expected_repeats is None:
        expected_repeats = max(observed_repeat_values, default=0)
    repeat_values = list(range(1, int(expected_repeats) + 1))
    control_digest_mismatches: list[Dict[str, Any]] = []
    control_incomplete: list[Dict[str, Any]] = []
    paired_controls: list[tuple[Mapping[str, Any], Mapping[str, Any]]] = []
    if set(BACKENDS) <= set(requested_backends):
        for case_id in sorted(control_cases):
            for repeat in repeat_values:
                current = by_key.get((case_id, repeat, "current"))
                optimized = by_key.get((case_id, repeat, "optimized"))
                if current is None or optimized is None:
                    control_incomplete.append(
                        {
                            "case_id": case_id,
                            "repeat": repeat,
                            "reason": "MISSING_BACKEND_RESULT",
                        }
                    )
                    continue
                if not current.get("family_complete") or not optimized.get(
                    "family_complete"
                ):
                    control_incomplete.append(
                        {
                            "case_id": case_id,
                            "repeat": repeat,
                            "current_status": current.get("status"),
                            "optimized_status": optimized.get("status"),
                        }
                    )
                    continue
                paired_controls.append((current, optimized))
                if (
                    current.get("objective") != optimized.get("objective")
                    or current.get("optimal_family_digest")
                    != optimized.get("optimal_family_digest")
                ):
                    control_digest_mismatches.append(
                        {
                            "case_id": case_id,
                            "repeat": repeat,
                            "current_objective": current.get("objective"),
                            "optimized_objective": optimized.get("objective"),
                            "current_digest": current.get(
                                "optimal_family_digest"
                            ),
                            "optimized_digest": optimized.get(
                                "optimal_family_digest"
                            ),
                        }
                    )

    effective_wall: Dict[str, list[float]] = {name: [] for name in BACKENDS}
    solver_wall: Dict[str, list[float]] = {name: [] for name in BACKENDS}
    for row in rows:
        backend = str(row["backend"])
        value = float(
            row.get(
                "subprocess_elapsed_seconds",
                row.get("solver_elapsed_seconds", deadline_seconds),
            )
        )
        if not row.get("worker_completed"):
            value = max(value, min(deadline_seconds, value))
        effective_wall[backend].append(value)
        solver_wall[backend].append(
            float(row.get("solver_elapsed_seconds", deadline_seconds))
        )
    observed_wall_totals = {
        backend: math.fsum(effective_wall[backend])
        for backend in requested_backends
    }
    backend_incomplete = {
        backend: [
            row
            for row in rows
            if row["backend"] == backend and not row.get("family_complete")
        ]
        for backend in requested_backends
    }
    backend_all_complete = {
        backend: not backend_incomplete[backend]
        for backend in requested_backends
    }
    completion_wall_totals = {
        backend: (
            observed_wall_totals[backend]
            if backend_all_complete[backend]
            else None
        )
        for backend in requested_backends
    }
    completion_time_speedup = None
    if (
        completion_wall_totals.get("current") is not None
        and completion_wall_totals.get("optimized") is not None
        and float(completion_wall_totals["optimized"]) > 0
    ):
        completion_time_speedup = float(
            completion_wall_totals["current"]
        ) / float(completion_wall_totals["optimized"])
    conservative_lower_bound_speedup = None
    if (
        backend_incomplete.get("current")
        and backend_all_complete.get("optimized")
        and observed_wall_totals.get("optimized", 0) > 0
    ):
        conservative_lower_bound_speedup = (
            observed_wall_totals["current"]
            / observed_wall_totals["optimized"]
        )
    current_deadline_censored = [
        {
            "case_id": row["case_id"],
            "repeat": row["repeat"],
            "status": row.get("status"),
        }
        for row in backend_incomplete.get("current", [])
        if "DEADLINE" in str(row.get("status", ""))
    ]
    performance_ratio = (
        completion_time_speedup
        if completion_time_speedup is not None
        else conservative_lower_bound_speedup
    )
    performance_ratio_kind = (
        "BOTH_BACKENDS_COMPLETE_TOTAL_WALL_SPEEDUP"
        if completion_time_speedup is not None
        else (
            "CONSERVATIVE_LOWER_BOUND_CURRENT_OVER_OPTIMIZED"
            if conservative_lower_bound_speedup is not None
            else "NOT_COMPARABLE"
        )
    )
    p95 = {
        backend: percentile(effective_wall[backend], 0.95)
        for backend in requested_backends
    }
    p95_not_worse = (
        p95.get("current") is not None
        and p95.get("optimized") is not None
        and float(p95["optimized"]) <= float(p95["current"])
    )
    current_rss = [
        float(current["ru_maxrss_kib"])
        for current, _ in paired_controls
        if current.get("ru_maxrss_kib") is not None
    ]
    optimized_rss = [
        float(optimized["ru_maxrss_kib"])
        for _, optimized in paired_controls
        if optimized.get("ru_maxrss_kib") is not None
    ]
    current_rss_p95 = percentile(current_rss, 0.95)
    optimized_rss_p95 = percentile(optimized_rss, 0.95)
    rss_ratio = (
        optimized_rss_p95 / current_rss_p95
        if current_rss_p95 and optimized_rss_p95 is not None
        else None
    )
    incomplete_ranked = [
        {
            "case_id": row["case_id"],
            "repeat": row["repeat"],
            "backend": row["backend"],
        }
        for row in rows
        if row.get("incomplete_ranked")
    ]
    optimized_incomplete = [
        {
            "case_id": row["case_id"],
            "repeat": row["repeat"],
            "status": row.get("status"),
        }
        for row in rows
        if row["backend"] == "optimized" and not row.get("family_complete")
    ]
    factorial_checks = [
        {
            "case_id": row["case_id"],
            "repeat": row["repeat"],
            "expected": row["factorial_oracle"][
                "expected_optimal_vertex_sets"
            ],
            "observed": row.get("optimal_family_count"),
            "pass": row.get("factorial_oracle_pass"),
        }
        for row in rows
        if row["backend"] == "optimized"
        and row.get("factorial_oracle") is not None
    ]
    selected_factorial_cases = {
        case_id
        for case_id in selected_case_ids
        if case_index[case_id].get("factorial_oracle") is not None
    }
    expected_factorial_check_count = (
        len(selected_factorial_cases) * int(expected_repeats)
        if "optimized" in requested_backends
        else 0
    )
    factorial_gate_pass = (
        expected_factorial_check_count > 0
        and len(factorial_checks) == expected_factorial_check_count
        and all(check["pass"] is True for check in factorial_checks)
    )

    comparable = set(BACKENDS) <= set(requested_backends)
    expected_worker_rows = (
        len(selected_case_ids)
        * len(requested_backends)
        * int(expected_repeats)
    )
    checks = {
        "both_backends_requested": comparable,
        "all_expected_worker_rows_present": len(rows) == expected_worker_rows,
        "controls_all_paired_and_complete": (
            comparable
            and not control_incomplete
            and len(paired_controls)
            == len(control_cases) * int(expected_repeats)
        ),
        "controls_zero_objective_or_digest_mismatch": (
            comparable and not control_digest_mismatches
        ),
        "optimized_all_selected_families_complete": not optimized_incomplete,
        "performance_ratio_ge_5x": (
            performance_ratio is not None and performance_ratio >= 5.0
        ),
        "p95_wall_not_worse": p95_not_worse,
        "control_p95_rss_ratio_le_1_5x": (
            rss_ratio is not None and rss_ratio <= 1.5
        ),
        "zero_incomplete_ranked": not incomplete_ranked,
        "factorial_oracles_expected_count_positive_and_all_pass": (
            factorial_gate_pass
        ),
    }
    return {
        "schema_name": "intersubmod.solver_stress_panel_receipt",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": (
            "PASS_STRESS_PANEL" if all(checks.values()) else "FAIL_CLOSED"
        ),
        "manifest_content_sha256": manifest["integrity"][
            "manifest_content_sha256"
        ],
        "configuration": {
            "selected_case_ids": list(selected_case_ids),
            "requested_backends": list(requested_backends),
            "deadline_seconds": deadline_seconds,
            "max_sets": max_sets,
            "parent_trees_materialized": False,
        },
        "coverage": {
            "selected_distinct_structural_keys": len(selected_case_ids),
            "control_distinct_structural_keys": len(control_cases),
            "worker_rows": len(rows),
            "expected_worker_rows": expected_worker_rows,
            "worker_failures": sum(
                not row.get("worker_completed") for row in rows
            ),
            "selected_factorial_distinct_structural_keys": len(
                selected_factorial_cases
            ),
            "expected_factorial_checks": expected_factorial_check_count,
            "observed_factorial_checks": len(factorial_checks),
        },
        "exactness": {
            "control_incomplete": control_incomplete,
            "control_objective_or_digest_mismatches": (
                control_digest_mismatches
            ),
            "optimized_incomplete": optimized_incomplete,
            "incomplete_ranked": incomplete_ranked,
            "factorial_oracle_checks": factorial_checks,
            "current_incomplete": [
                {
                    "case_id": row["case_id"],
                    "repeat": row["repeat"],
                    "status": row.get("status"),
                }
                for row in backend_incomplete.get("current", [])
            ],
            "current_deadline_censored": current_deadline_censored,
        },
        "performance": {
            "acceptance_wall_basis": (
                "end-to-end isolated worker subprocess wall, including "
                "startup/import/source verification"
            ),
            "observed_wall_seconds_including_censored_runs": (
                observed_wall_totals
            ),
            "completion_wall_seconds": completion_wall_totals,
            "total_completion_time_speedup_current_over_optimized": (
                completion_time_speedup
            ),
            "conservative_lower_bound_speedup_current_over_optimized": (
                conservative_lower_bound_speedup
            ),
            "reported_ratio_kind": performance_ratio_kind,
            "reported_current_over_optimized_ratio": performance_ratio,
            "current_has_incomplete_runs": bool(
                backend_incomplete.get("current")
            ),
            "current_deadline_censored": bool(current_deadline_censored),
            "current_deadline_censored_run_count": len(
                current_deadline_censored
            ),
            "censoring_interpretation": (
                "When current is incomplete, the reported ratio is a "
                "conservative lower bound, not a ratio of two completion "
                "times."
            ),
            "p95_effective_wall_seconds": p95,
            "p95_not_worse": p95_not_worse,
            "total_solver_unit_wall_seconds": {
                backend: math.fsum(solver_wall[backend])
                for backend in requested_backends
            },
            "p95_solver_unit_wall_seconds": {
                backend: percentile(solver_wall[backend], 0.95)
                for backend in requested_backends
            },
            "control_ru_maxrss_kib_p95": {
                "current": current_rss_p95,
                "optimized": optimized_rss_p95,
            },
            "control_rss_ratio_optimized_over_current": rss_ratio,
        },
        "checks": checks,
        "all_pass": all(checks.values()),
        "runs": list(rows),
    }


def parse_max_sets(value: str) -> int | None:
    normalized = value.strip().lower()
    if normalized in {"none", "unbounded", "null"}:
        return None
    parsed = int(normalized)
    if parsed < 1:
        raise argparse.ArgumentTypeError("max-sets must be positive or none")
    return parsed


def resolve_manifest_arguments(
    *,
    authority_pointer: pathlib.Path | None,
    manifest_path: pathlib.Path | None,
    authority_status: str | None,
    verify_inputs: bool,
) -> tuple[pathlib.Path, Dict[str, Any], Dict[str, Any]]:
    if authority_pointer is not None:
        if manifest_path is not None:
            raise StressPanelError(
                "specify --authority-pointer or --manifest, not both"
            )
        return resolve_authority_pointer(
            authority_pointer,
            verify_inputs=verify_inputs,
        )
    if manifest_path is None:
        raise StressPanelError("an authority pointer or manifest is required")
    return resolve_direct_manifest(
        manifest_path,
        authority_status=authority_status,
        verify_inputs=verify_inputs,
    )


def write_json_with_sha256(
    path: pathlib.Path,
    document: Mapping[str, Any],
) -> pathlib.Path:
    encoded = (
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    with path.open("xb") as handle:
        handle.write(encoded)
    sidecar = path.with_name(path.name + ".sha256")
    with sidecar.open("x", encoding="utf-8") as handle:
        handle.write(f"{hashlib.sha256(encoded).hexdigest()}  {path.name}\n")
    return sidecar


def build_run_schedule(
    selected_case_ids: Sequence[str],
    requested_backends: Sequence[str],
    repeats: int,
) -> list[Dict[str, Any]]:
    """Balance first/second backend position when repeated timing is requested."""
    schedule: list[Dict[str, Any]] = []
    for case_index, case_id in enumerate(selected_case_ids):
        for repeat in range(1, repeats + 1):
            order = list(requested_backends)
            if (
                repeats > 1
                and set(BACKENDS) <= set(order)
                and (case_index + repeat) % 2 == 0
            ):
                order.reverse()
            for backend_position, backend in enumerate(order, start=1):
                schedule.append(
                    {
                        "schedule_index": len(schedule) + 1,
                        "case_id": case_id,
                        "repeat": repeat,
                        "backend": backend,
                        "backend_position_within_repeat": backend_position,
                        "backend_order_for_case_repeat": list(order),
                    }
                )
    return schedule


def worker_main(args: argparse.Namespace) -> int:
    # The controller binds the large truncated inputs once.  Per-case workers
    # re-verify the immutable manifest and executable sources, avoiding an
    # irrelevant repeated scan of the diagnostic gzip files.
    _, manifest, authority = resolve_manifest_arguments(
        authority_pointer=args.worker_authority_pointer,
        manifest_path=args.worker_manifest,
        authority_status=args.worker_authority_status,
        verify_inputs=False,
    )
    case = _case_index(manifest).get(args.worker_case_id)
    if case is None:
        raise StressPanelError(f"unknown case id: {args.worker_case_id}")
    result = run_backend_case(
        manifest,
        case,
        backend=args.worker_backend,
        deadline_seconds=args.deadline,
        max_sets=args.max_sets,
        q_max=args.q_max,
    )
    result["authority"] = authority
    if args.worker_output.exists():
        raise StressPanelError(f"worker output exists: {args.worker_output}")
    args.worker_output.parent.mkdir(parents=True, exist_ok=True)
    with args.worker_output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


def _controller_main_impl(args: argparse.Namespace) -> int:
    manifest_path, manifest, authority = resolve_manifest_arguments(
        authority_pointer=args.authority_pointer,
        manifest_path=args.manifest,
        authority_status=args.authority_status,
        verify_inputs=True,
    )
    selected_case_ids = select_case_ids(
        manifest,
        panel=args.panel,
        selection=args.selection,
    )
    if not selected_case_ids:
        raise StressPanelError("case selection is empty")
    if args.output_dir.exists():
        raise StressPanelError(
            f"write-once run output already exists: {args.output_dir}"
        )
    args.output_dir.mkdir(parents=True, exist_ok=False)
    raw_dir = args.output_dir / "raw"
    raw_dir.mkdir()
    cases = _case_index(manifest)
    rows: list[Dict[str, Any]] = []
    env = os.environ.copy()
    env.update(THREAD_ENV)
    requested_backends = tuple(args.backends)
    if len(set(requested_backends)) != len(requested_backends):
        raise StressPanelError("duplicate backends are not allowed")
    schedule = build_run_schedule(
        selected_case_ids,
        requested_backends,
        args.repeats,
    )

    for job in schedule:
        case_id = str(job["case_id"])
        repeat = int(job["repeat"])
        backend = str(job["backend"])
        stem = f"{case_id}.r{repeat}.{backend}"
        worker_output = raw_dir / f"{stem}.json"
        stdout_path = raw_dir / f"{stem}.stdout"
        stderr_path = raw_dir / f"{stem}.stderr"
        command = [
            sys.executable,
            str(SCRIPT_PATH),
            "--worker-case-id",
            case_id,
            "--worker-backend",
            backend,
            "--worker-output",
            str(worker_output),
            "--deadline",
            str(args.deadline),
            "--max-sets",
            "none" if args.max_sets is None else str(args.max_sets),
            "--q-max",
            str(args.q_max),
        ]
        if args.authority_pointer is not None:
            command[2:2] = [
                "--worker-authority-pointer",
                str(args.authority_pointer),
            ]
        else:
            command[2:2] = [
                "--worker-manifest",
                str(manifest_path),
                "--worker-authority-status",
                str(args.authority_status),
            ]
        started = time.perf_counter()
        try:
            completed = subprocess.run(
                command,
                cwd=str(pathlib.Path.cwd()),
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=args.deadline + args.subprocess_grace,
                check=False,
            )
            timed_out = False
            stdout = completed.stdout
            stderr = completed.stderr
            returncode = completed.returncode
        except subprocess.TimeoutExpired as error:
            timed_out = True
            stdout = (
                error.stdout.decode()
                if isinstance(error.stdout, bytes)
                else (error.stdout or "")
            )
            stderr = (
                error.stderr.decode()
                if isinstance(error.stderr, bytes)
                else (error.stderr or "")
            )
            returncode = 124
        subprocess_elapsed = time.perf_counter() - started
        stdout_path.write_text(stdout, encoding="utf-8")
        stderr_path.write_text(stderr, encoding="utf-8")
        if returncode == 0 and worker_output.is_file():
            row = json.loads(worker_output.read_text(encoding="utf-8"))
            row["worker_completed"] = True
        else:
            row = {
                "schema": "intersubmod.solver_stress_panel.worker.v1",
                "case_id": case_id,
                "structural_key_sha256": cases[case_id][
                    "structural_key_sha256"
                ],
                "backend": backend,
                "status": (
                    "WORKER_HARD_DEADLINE"
                    if timed_out
                    else "WORKER_ERROR"
                ),
                "complete": False,
                "family_complete": False,
                "objective_certified": False,
                "ranking_allowed": False,
                "incomplete_ranked": False,
                "parent_trees_materialized": False,
                "optimal_family_count": 0,
                "optimal_family_digest": None,
                "solver_elapsed_seconds": min(
                    subprocess_elapsed,
                    args.deadline,
                ),
                "worker_completed": False,
                "worker_returncode": returncode,
                "worker_error_tail": stderr[-2000:],
            }
        row["repeat"] = repeat
        row["subprocess_elapsed_seconds"] = subprocess_elapsed
        row["worker_returncode"] = returncode
        row["worker_stdout_path"] = str(stdout_path.resolve())
        row["worker_stderr_path"] = str(stderr_path.resolve())
        row["schedule_index"] = int(job["schedule_index"])
        row["backend_position_within_repeat"] = int(
            job["backend_position_within_repeat"]
        )
        row["backend_order_for_case_repeat"] = list(
            job["backend_order_for_case_repeat"]
        )
        rows.append(row)

    receipt = build_summary(
        manifest,
        rows,
        selected_case_ids=selected_case_ids,
        requested_backends=requested_backends,
        deadline_seconds=args.deadline,
        max_sets=args.max_sets,
        expected_repeats=args.repeats,
    )
    receipt["configuration"].update(
        {
            "panel": args.panel,
            "selection": args.selection,
            "repeats": args.repeats,
            "q_max": args.q_max,
            "subprocess_grace_seconds": args.subprocess_grace,
            "authority": authority,
            "execution_schedule": schedule,
            "backend_order_policy": (
                "ALTERNATE_FIRST_BACKEND_BY_CASE_AND_REPEAT"
                if args.repeats > 1
                else "REQUESTED_ORDER_SINGLE_REPEAT"
            ),
        }
    )
    receipt_path = args.output_dir / "receipt.json"
    receipt_sidecar = write_json_with_sha256(receipt_path, receipt)
    print(
        json.dumps(
            {
                "status": receipt["verdict"],
                "output": str(receipt_path.resolve()),
                "output_sha256_sidecar": str(receipt_sidecar.resolve()),
                "selected_keys": len(selected_case_ids),
                "worker_rows": len(rows),
                "control_digest_mismatches": len(
                    receipt["exactness"][
                        "control_objective_or_digest_mismatches"
                    ]
                ),
                "optimized_incomplete": len(
                    receipt["exactness"]["optimized_incomplete"]
                ),
                "reported_ratio": receipt["performance"][
                    "reported_current_over_optimized_ratio"
                ],
                "reported_ratio_kind": receipt["performance"][
                    "reported_ratio_kind"
                ],
                "p95_not_worse": receipt["performance"]["p95_not_worse"],
                "rss_ratio": receipt["performance"][
                    "control_rss_ratio_optimized_over_current"
                ],
                "incomplete_ranked": len(
                    receipt["exactness"]["incomplete_ranked"]
                ),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if receipt["all_pass"] else 3


def _write_partial_run_marker(
    args: argparse.Namespace,
    error: BaseException,
) -> None:
    if args.output_dir is None or not args.output_dir.is_dir():
        return
    if (args.output_dir / "receipt.json").exists():
        return
    marker = args.output_dir / "PARTIAL_RUN.json"
    if marker.exists():
        return
    raw_dir = args.output_dir / "raw"
    completed_worker_outputs = (
        len(list(raw_dir.glob("*.json"))) if raw_dir.is_dir() else 0
    )
    document = {
        "schema": "intersubmod.solver_stress_panel.partial_run.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PARTIAL_RUN_NOT_ACCEPTANCE_EVIDENCE",
        "reason": type(error).__name__,
        "error": str(error),
        "completed_worker_output_files": completed_worker_outputs,
        "configuration": {
            "authority_pointer": (
                str(args.authority_pointer.resolve())
                if args.authority_pointer is not None
                else None
            ),
            "manifest": (
                str(args.manifest.resolve())
                if args.manifest is not None
                else None
            ),
            "authority_status": args.authority_status,
            "panel": args.panel,
            "selection": args.selection,
            "backends": list(args.backends),
            "repeats": args.repeats,
            "deadline_seconds": args.deadline,
            "max_sets": args.max_sets,
            "q_max": args.q_max,
        },
    }
    write_json_with_sha256(marker, document)


def controller_main(args: argparse.Namespace) -> int:
    try:
        return _controller_main_impl(args)
    except BaseException as error:
        _write_partial_run_marker(args, error)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path)
    parser.add_argument("--authority-pointer", type=pathlib.Path)
    parser.add_argument("--authority-status")
    parser.add_argument("--output-dir", type=pathlib.Path)
    parser.add_argument(
        "--panel",
        choices=("primary", "sensitivity", "both"),
        default="primary",
    )
    parser.add_argument(
        "--selection",
        choices=("all", "controls", "tails", "factorial"),
        default="all",
    )
    parser.add_argument(
        "--backends",
        nargs="+",
        choices=BACKENDS,
        default=list(BACKENDS),
    )
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--deadline", type=float, default=30.0)
    parser.add_argument("--max-sets", type=parse_max_sets, default=256)
    parser.add_argument("--q-max", type=int, default=8)
    parser.add_argument("--subprocess-grace", type=float, default=5.0)

    parser.add_argument("--worker-manifest", type=pathlib.Path)
    parser.add_argument("--worker-authority-pointer", type=pathlib.Path)
    parser.add_argument("--worker-authority-status")
    parser.add_argument("--worker-case-id")
    parser.add_argument("--worker-backend", choices=BACKENDS)
    parser.add_argument("--worker-output", type=pathlib.Path)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if (
            args.worker_manifest is not None
            or args.worker_authority_pointer is not None
        ):
            required = (
                args.worker_case_id,
                args.worker_backend,
                args.worker_output,
            )
            if any(value is None for value in required):
                parser.error("worker mode requires case/backend/output")
            return worker_main(args)
        if (
            args.manifest is None
            and args.authority_pointer is None
        ) or args.output_dir is None:
            parser.error(
                "controller mode requires an authority pointer/manifest "
                "and --output-dir"
            )
        if args.repeats < 1:
            parser.error("--repeats must be positive")
        if not math.isfinite(args.deadline) or args.deadline <= 0:
            parser.error("--deadline must be finite and positive")
        if not math.isfinite(args.subprocess_grace) or args.subprocess_grace < 0:
            parser.error("--subprocess-grace must be finite and nonnegative")
        return controller_main(args)
    except KeyboardInterrupt:
        print(
            json.dumps(
                {"status": "PARTIAL_RUN", "error": "KeyboardInterrupt"},
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 130
    except Exception as error:
        print(
            json.dumps(
                {"status": "FAIL_CLOSED", "error": str(error)},
                ensure_ascii=False,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
