#!/usr/bin/env python3
"""Run the isolated compressed primary-likelihood feasibility probe."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from compressed_vaf_rank_probe import (
    brute_force_perfect_vertex_sets,
    exhaustive_current_rank,
    rank_perfect_family_lazy,
)
from run_perfect_family_count_validation import (
    MANIFEST_SCHEMA,
    POINTER_SCHEMA,
    REQUIRED_AUTHORITY_STATUS,
    ValidationError as AuthorityValidationError,
    _index_manifest_cases,
    _require_sha256,
    _resolve_authority_pointer,
)


BASE = Path(__file__).resolve().parents[1]
DEFAULT_POINTER = (
    BASE
    / "results"
    / "solver_stress_panel_v1"
    / "AUTHORITATIVE_MANIFEST.json"
)
DEFAULT_OUTPUT = (
    BASE
    / "results"
    / "compressed_vaf_rank_probe_v2"
    / "receipt.json"
)
HARD_CASE_ID = "k11_m11_09b1da787e58efed"
RECEIPT_SCHEMA = "intersubmod.compressed_vaf_rank_probe.receipt.v2"
SCRIPT_PATH = Path(__file__).resolve()
PROTOTYPE_PATH = SCRIPT_PATH.with_name("compressed_vaf_rank_probe.py")
HISTORICAL_RANKER_PATH = (
    BASE.parents[1]
    / "docs"
    / "methodology"
    / "_assets"
    / "20260627_subclone_4axis_teaching"
    / "scripts"
    / "read_af_tree_ordering_multisample.py"
)
M2_RANKER_PATH = (
    BASE.parents[1]
    / "research"
    / "20260716_read_linked_hypercube_exact_likelihood_validation"
    / "scripts"
    / "build_m2_patterns_and_rank.py"
)
HYPERCUBE_EXACT_PATH = M2_RANKER_PATH.with_name("hypercube_exact.py")
PERFECT_COUNTER_PATH = SCRIPT_PATH.with_name("perfect_family_count.py")
AUTHORITY_VALIDATOR_PATH = SCRIPT_PATH.with_name(
    "run_perfect_family_count_validation.py"
)


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def source_paths() -> dict[str, Path]:
    return {
        "runner": SCRIPT_PATH,
        "prototype": PROTOTYPE_PATH,
        "historical_read_af_heuristic": HISTORICAL_RANKER_PATH,
        "current_m2_ranker": M2_RANKER_PATH,
        "hypercube_exact": HYPERCUBE_EXACT_PATH,
        "perfect_family_counter": PERFECT_COUNTER_PATH,
        "r3_authority_validator": AUTHORITY_VALIDATOR_PATH,
    }


def snapshot_files(paths: dict[str, Path]) -> dict[str, dict[str, str]]:
    snapshot: dict[str, dict[str, str]] = {}
    for name, raw_path in sorted(paths.items()):
        path = raw_path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"bound source is missing: {path}")
        snapshot[name] = {
            "path": str(path),
            "sha256": sha256_path(path),
        }
    return snapshot


def assert_snapshot_unchanged(
    before: dict[str, dict[str, str]],
    after: dict[str, dict[str, str]],
    *,
    label: str,
) -> None:
    if before != after:
        changed = sorted(
            name
            for name in set(before) | set(after)
            if before.get(name) != after.get(name)
        )
        raise RuntimeError(
            f"{label} changed during execution: {', '.join(changed)}"
        )


def _quality_groups(k: int):
    groups = []
    for alt_bit in range(k):
        pattern = "".join(
            "A" if bit == alt_bit else "R" for bit in range(k)
        )
        groups.append((pattern, tuple(35 for _ in range(k)), k - alt_bit))
    groups.append(("A" * k, tuple(35 for _ in range(k)), 1))
    return tuple(groups)


def run_small_oracles() -> dict[str, Any]:
    started = time.perf_counter()
    cases = []
    for m in range(1, 6):
        partial = ["A" * m]
        quality = _quality_groups(m)
        family = brute_force_perfect_vertex_sets(
            full_patterns=[],
            partial_patterns=partial,
            k=m,
            structural_alt_universe_mask=(1 << m) - 1,
        )
        exhaustive = exhaustive_current_rank(
            vertex_sets=family,
            quality_groups=quality,
            k=m,
            tie_tolerance=1e-6,
        )
        full_evaluation = rank_perfect_family_lazy(
            full_patterns=[],
            partial_patterns=partial,
            k=m,
            quality_groups=quality,
            structural_alt_universe_mask=(1 << m) - 1,
            tie_tolerance=1e-6,
            max_candidates=10_000,
            deadline_seconds=30,
            max_tie_ids=10_000,
            enable_numerical_pruning=False,
        )
        float_pruned_diagnostic = rank_perfect_family_lazy(
            full_patterns=[],
            partial_patterns=partial,
            k=m,
            quality_groups=quality,
            structural_alt_universe_mask=(1 << m) - 1,
            tie_tolerance=1e-6,
            max_candidates=10_000,
            deadline_seconds=30,
            max_tie_ids=10_000,
            bound_mode="union_mixture_float",
            enable_numerical_pruning=True,
        )
        best_delta = abs(
            float(full_evaluation["best_log_likelihood"])
            - float(exhaustive["best_log_likelihood"])
        )
        oracle_confirmed = (
            full_evaluation["ranking_complete"]
            and full_evaluation["machine_exact_by_full_enumeration"]
            and full_evaluation["logical_family_count"] == len(family)
            and best_delta <= 1e-10
            and full_evaluation["winner_vertex_set_ids"]
            == exhaustive["winner_vertex_set_ids"]
            and full_evaluation["tie_count"] == exhaustive["tie_count"]
        )
        diagnostic_pruned = (
            float_pruned_diagnostic["pruned_candidate_count"] > 0
        )
        diagnostic_fail_closed = (
            (
                not float_pruned_diagnostic["ranking_complete"]
                and float_pruned_diagnostic["best_log_likelihood"] is None
                and float_pruned_diagnostic["winner_vertex_set_ids"] == []
                and not float_pruned_diagnostic[
                    "numerical_bound_certified"
                ]
            )
            if diagnostic_pruned
            else float_pruned_diagnostic["ranking_complete"]
        )
        record = {
            "m": m,
            "family_count": len(family),
            "traceback_count": full_evaluation["logical_family_count"],
            "evaluated_candidate_count": full_evaluation[
                "evaluated_candidate_count"
            ],
            "pruned_candidate_count": full_evaluation[
                "pruned_candidate_count"
            ],
            "exhaustive_best_log_likelihood": exhaustive[
                "best_log_likelihood"
            ],
            "full_evaluation_best_log_likelihood": full_evaluation[
                "best_log_likelihood"
            ],
            "best_score_abs_delta": best_delta,
            "exhaustive_winner_ids": exhaustive["winner_vertex_set_ids"],
            "full_evaluation_winner_ids": full_evaluation[
                "winner_vertex_set_ids"
            ],
            "exhaustive_tie_count": exhaustive["tie_count"],
            "full_evaluation_tie_count": full_evaluation["tie_count"],
            "oracle_confirmed": oracle_confirmed,
            "float_pruned_diagnostic": {
                "evaluated_candidate_count": float_pruned_diagnostic[
                    "evaluated_candidate_count"
                ],
                "pruned_candidate_count": float_pruned_diagnostic[
                    "pruned_candidate_count"
                ],
                "numerical_bound_certified": float_pruned_diagnostic[
                    "numerical_bound_certified"
                ],
                "ranking_complete": float_pruned_diagnostic[
                    "ranking_complete"
                ],
                "status": float_pruned_diagnostic["status"],
                "fail_closed": diagnostic_fail_closed,
            },
            "pass": oracle_confirmed and diagnostic_fail_closed,
            "full_evaluation_elapsed_seconds": full_evaluation[
                "elapsed_seconds"
            ],
        }
        cases.append(record)

    tie_quality = (
        ("AR", (40, 40), 3),
        ("RA", (40, 40), 3),
        ("AA", (40, 40), 1),
    )
    tie_case = rank_perfect_family_lazy(
        full_patterns=[],
        partial_patterns=["AA"],
        k=2,
        quality_groups=tie_quality,
        structural_alt_universe_mask=3,
        max_candidates=10,
        deadline_seconds=10,
        enable_numerical_pruning=False,
    )
    tie_pass = (
        tie_case["ranking_complete"]
        and tie_case["logical_family_count"] == 2
        and tie_case["tie_count"] == 2
    )
    return {
        "scope": "m=1..5 independent parent-vector exhaustive oracle",
        "oracle_confirmation_contract": (
            "independent parent-vector family enumeration plus exhaustive "
            "current-endpoint candidate scoring; does not certify general "
            "float-pruned core search"
        ),
        "cases": cases,
        "symmetric_tie_control": {
            "logical_family_count": tie_case["logical_family_count"],
            "tie_count": tie_case["tie_count"],
            "ranking_complete": tie_case["ranking_complete"],
            "pass": tie_pass,
        },
        "n_cases": len(cases) + 1,
        "n_pass": sum(bool(case["pass"]) for case in cases)
        + int(tie_pass),
        "all_pass": all(bool(case["pass"]) for case in cases)
        and tie_pass,
        "elapsed_seconds": time.perf_counter() - started,
    }


def _parse_quality_vector(value: str, *, k: int) -> tuple[int, ...]:
    tokens = value.split(",")
    if len(tokens) != k:
        raise ValueError("frozen BQ vector length differs from k")
    return tuple(-1 if token == "" else int(token) for token in tokens)


def load_hard_case(
    pointer_path: Path,
) -> tuple[
    dict[str, Any],
    tuple[tuple[str, tuple[int, ...], int], ...],
    dict[str, Any],
]:
    """Load one hard unit only through the fail-closed R3 authority chain."""
    try:
        manifest_path, manifest, pointer_binding = _resolve_authority_pointer(
            pointer_path
        )
        cases = _index_manifest_cases(manifest)
    except (AuthorityValidationError, OSError, json.JSONDecodeError) as error:
        raise ValueError(f"R3 authority verification failed: {error}") from error
    if HARD_CASE_ID not in cases:
        raise ValueError(f"R3 manifest lacks hard case {HARD_CASE_ID}")
    case = dict(cases[HARD_CASE_ID])
    units = sorted(
        (
            unit
            for unit in manifest["units"]
            if unit["case_id"] == HARD_CASE_ID
        ),
        key=lambda unit: unit["unit_id"],
    )
    if not units:
        raise ValueError(f"R3 manifest has no unit for {HARD_CASE_ID}")
    declared_unit_ids = case.get("unit_ids")
    observed_unit_ids = [unit["unit_id"] for unit in units]
    if (
        not isinstance(declared_unit_ids, list)
        or not declared_unit_ids
        or sorted(declared_unit_ids) != observed_unit_ids
        or len(set(observed_unit_ids)) != len(observed_unit_ids)
    ):
        raise ValueError("R3 hard-case unit IDs are missing, duplicated, or inconsistent")
    unit = units[0]
    identity = unit["identity"]
    try:
        pattern_binding = manifest["bindings"]["input_files"]["patterns"]
        pattern_path_value = pattern_binding["path"]
        expected_pattern_sha = _require_sha256(
            pattern_binding["sha256"],
            "R3 patterns binding SHA-256",
        )
    except (KeyError, TypeError, AuthorityValidationError) as error:
        raise ValueError("R3 patterns binding is malformed") from error
    if not isinstance(pattern_path_value, str) or not pattern_path_value:
        raise ValueError("R3 patterns binding path is invalid")
    pattern_path = Path(pattern_path_value).resolve()
    if not pattern_path.is_file():
        raise ValueError(f"R3 frozen patterns source is missing: {pattern_path}")
    if sha256_path(pattern_path) != expected_pattern_sha:
        raise ValueError("frozen pattern source SHA mismatch")
    pointer_path = pointer_path.resolve()
    pointer_sidecar = pointer_path.with_name(pointer_path.name + ".sha256")
    manifest_sidecar = manifest_path.with_name(manifest_path.name + ".sha256")
    hard_input_paths = {
        "authority_pointer": pointer_path,
        "authority_pointer_sidecar": pointer_sidecar,
        "manifest": manifest_path,
        "manifest_sidecar": manifest_sidecar,
        "patterns": pattern_path,
    }
    input_snapshot_before = snapshot_files(hard_input_paths)
    if (
        input_snapshot_before["authority_pointer"]["sha256"]
        != pointer_binding["sha256"]
        or input_snapshot_before["manifest"]["sha256"]
        != _require_sha256(
            json.loads(pointer_path.read_text(encoding="utf-8"))[
                "authoritative_manifest_file_sha256"
            ],
            "pointer manifest file SHA-256",
        )
        or input_snapshot_before["patterns"]["sha256"]
        != expected_pattern_sha
    ):
        raise ValueError("R3 hard input snapshot differs from verified bindings")

    groups = []
    matched_started = False
    structural_k = int(case["structural_input"]["k"])
    with gzip.open(pattern_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            matches = (
                row["dataset"] == identity["dataset"]
                and row["chrom"] == identity["chrom"]
                and row["threshold"] == identity["threshold"]
                and row["component_basis"] == identity["component_basis"]
                and row["phase_set"] == identity["phase_set"]
                and row["component_id"] == identity["component_id"]
                and row["family"] == identity["family"]
                and row["structural_exact_pattern_minread"]
                == identity["structural_exact_pattern_minread"]
            )
            if matched_started and not matches:
                break
            if not matches:
                continue
            matched_started = True
            if row["scoring_eligible"].lower() != "true":
                continue
            k = int(row["k"])
            if k != structural_k:
                raise ValueError("frozen BQ row k differs from R3 structural case")
            groups.append(
                (
                    row["pattern_rax"],
                    _parse_quality_vector(
                        row["fixed_base_qualities"],
                        k=k,
                    ),
                    int(row["n_molecules"]),
                )
            )
    if not groups:
        raise ValueError("no frozen hard-case scoring quality groups")
    input_snapshot_after = snapshot_files(hard_input_paths)
    assert_snapshot_unchanged(
        input_snapshot_before,
        input_snapshot_after,
        label="R3 hard inputs",
    )
    return (
        case,
        tuple(groups),
        {
            "authority_pointer": pointer_binding["path"],
            "authority_pointer_sha256": pointer_binding["sha256"],
            "authority_pointer_sidecar": str(pointer_sidecar),
            "authority_pointer_sidecar_sha256": sha256_path(pointer_sidecar),
            "authority_pointer_schema": POINTER_SCHEMA,
            "authority_pointer_status": REQUIRED_AUTHORITY_STATUS,
            "manifest": str(manifest_path),
            "manifest_sha256": sha256_path(manifest_path),
            "manifest_sidecar": str(manifest_sidecar),
            "manifest_sidecar_sha256": sha256_path(manifest_sidecar),
            "manifest_schema": MANIFEST_SCHEMA,
            "manifest_authority_status": manifest["authority"]["status"],
            "manifest_content_sha256": manifest["integrity"][
                "manifest_content_sha256"
            ],
            "manifest_checks_all_true": str(
                bool(manifest["checks"])
                and all(value is True for value in manifest["checks"].values())
            ).lower(),
            "patterns": str(pattern_path),
            "patterns_sha256": expected_pattern_sha,
            "unit_id": unit["unit_id"],
            "integrity_file_snapshot": input_snapshot_before,
            "input_hashes_pre_equals_post": True,
        },
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--authority-pointer",
        type=Path,
        default=DEFAULT_POINTER,
    )
    parser.add_argument(
        "--run-hard",
        action="store_true",
        help="Run the bounded hard-case diagnostic; disabled by default.",
    )
    parser.add_argument("--hard-max-candidates", type=int, default=20)
    parser.add_argument("--hard-deadline", type=float, default=5.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if "compressed_vaf_rank_probe_v1" in args.output.resolve().parts:
        raise ValueError("v1 output is immutable; write a new v2 output")
    started = time.perf_counter()
    source_snapshot_before = snapshot_files(source_paths())
    oracle = run_small_oracles()
    hard: dict[str, Any]
    hard_bindings: dict[str, Any] = {}
    if args.run_hard:
        hard = {
            "status": "NOT_RUN_NO_PROCESS_WATCHDOG",
            "reason": (
                "the in-process deadline is checked only after DP and all "
                "top-level branch/possible-vertex construction; without an "
                "outer process watchdog it is not a hard wall-time bound"
            ),
            "requested": True,
            "requested_max_candidates": args.hard_max_candidates,
            "requested_deadline_seconds": args.hard_deadline,
            "deadline_is_hard_wall_bound": False,
            "ranking_complete": False,
        }
    else:
        hard = {
            "status": "NOT_RUN_P1_REVIEW_GATE",
            "reason": (
                "09b remains disabled until an outer process watchdog and P1 "
                "review are complete"
            ),
            "requested": False,
            "deadline_is_hard_wall_bound": False,
            "ranking_complete": False,
        }
    source_snapshot_after = snapshot_files(source_paths())
    assert_snapshot_unchanged(
        source_snapshot_before,
        source_snapshot_after,
        label="bound sources",
    )
    source_snapshot_prewrite = snapshot_files(source_paths())
    assert_snapshot_unchanged(
        source_snapshot_before,
        source_snapshot_prewrite,
        label="bound sources before receipt write",
    )
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "created_at_unix": time.time(),
        "task_type": "A_EXPLORATORY_PILOT_PARTIAL",
        "scope": (
            "isolated recurrence-free perfect families; current primary BQ "
            "profile mixture likelihood only"
        ),
        "verdict": (
            "PASS_SMALL_ORACLE_HARD_NOT_RUN_P1_GATE"
            if oracle["all_pass"]
            else "FAIL_SMALL_ORACLE"
        ),
        "endpoint_contract": {
            "same": (
                "primary BQ-aware candidate vertex-set profile "
                "log-likelihood and tie_tolerance membership"
            ),
            "not_same": (
                "fixed-error sensitivity, bootstrap, responsibilities, "
                "all-candidate TSV, and full M2 release receipt"
            ),
        },
        "complexity": {
            "traceback_count_dp_time": "O(3^m poly(m))",
            "traceback_count_dp_space": "O(2^m + m^2)",
            "top_level_branch_and_possible_vertex_materialization": (
                "worst-case approximately O(m * 3^m); all top-level branches "
                "and each possible_vertices tuple are materialized before "
                "enumeration"
            ),
            "recursive_branch_and_bound": False,
            "machine_complete_ranking_worst_case": (
                "Theta(F * Fit(m,n_patterns)); without directed-rounding "
                "certificates every candidate must be evaluated"
            ),
            "float_pruning_role": (
                "diagnostic only; ordinary float relaxation may reduce "
                "evaluations but cannot publish authoritative best/ties"
            ),
            "full_tie_class_lower_bound": (
                "Omega(number_of_tied_vertex_sets) unless an additional exact "
                "tie-class compression theorem is supplied"
            ),
        },
        "small_oracle": oracle,
        "hard_case": hard,
        "hard_bindings": hard_bindings,
        "numerical_certificate_policy": {
            "ordinary_float_upper_bound_machine_certified": False,
            "authoritative_core_result_requires": (
                "evaluated_candidate_count == logical_family_count, "
                "pruned_candidate_count == 0, search_complete, and "
                "tie_class_complete"
            ),
            "small_oracle_scope": (
                "oracle_confirmed is an external m<=5 exhaustive comparison; "
                "it does not certify general float-pruned searches"
            ),
        },
        "deadline_contract": {
            "core_deadline_is_hard_wall_bound": False,
            "uncovered_work": (
                "integer subset DP plus construction/materialization of all "
                "top-level branches and possible_vertices; float UB fitting "
                "also occurs before loop deadline checks"
            ),
            "hard_policy": "NOT_RUN until an outer process watchdog exists",
        },
        "hard_authority_contract": {
            "pointer_schema": POINTER_SCHEMA,
            "required_status": REQUIRED_AUTHORITY_STATUS,
            "manifest_schema": MANIFEST_SCHEMA,
            "verification": (
                "pointer byte sidecar + schema/status + relative-path "
                "containment + manifest byte digest/sidecar + semantic digest "
                "+ all frozen checks + case/unit structural validation"
            ),
        },
        "source_bindings": source_snapshot_before,
        "runtime_bindings": {
            "numpy_version": np.__version__,
            "python_version": platform.python_version(),
            "python_executable": sys.executable,
        },
        "integrity": {
            "source_hashes_pre_equals_post": True,
            "source_hashes_pre_equals_prewrite": True,
            "hard_input_hashes_pre_equals_post": (
                hard_bindings.get(
                    "input_hashes_pre_equals_post",
                    "NOT_APPLICABLE_HARD_NOT_RUN",
                )
            ),
        },
        "claim_ceiling": (
            "Core best/ties are machine-complete only after full candidate "
            "evaluation with no float-bound pruning; m<=5 oracle_confirmed is "
            "a separate external exhaustive result."
        ),
        "canonical_or_production_promotion_allowed": False,
        "elapsed_seconds": time.perf_counter() - started,
    }
    if not math.isfinite(receipt["elapsed_seconds"]):
        raise AssertionError("receipt elapsed time is non-finite")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(args.output, flags, 0o444)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump(
            receipt,
            handle,
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
        handle.write("\n")
    digest = sha256_path(args.output)
    sidecar = args.output.with_name(args.output.name + ".sha256")
    descriptor = os.open(sidecar, flags, 0o444)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        handle.write(f"{digest}  {args.output.name}\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "verdict": receipt["verdict"],
                "small_oracle_pass": oracle["all_pass"],
                "hard_status": hard["status"],
                "elapsed_seconds": receipt["elapsed_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
