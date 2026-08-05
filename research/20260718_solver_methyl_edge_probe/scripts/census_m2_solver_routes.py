#!/usr/bin/env python3
"""Reproduce diagnostic solver-route statistics from the truncated M2 pilot.

This script is deliberately read-only with respect to the frozen pilot.  It
uses the frozen hypercube implementation to reconstruct observed-ALT active
dimensions and exact-preserving partial-group reductions.  The input gzip
streams are known to be truncated, so only complete TSV records are consumed
and the output is explicitly marked diagnostic-only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Mapping, Sequence, Tuple


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FROZEN_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260716_read_linked_hypercube_exact_likelihood_validation/"
    "m2_frozen_release_v2"
)
DEFAULT_RANKING_ROOT = (
    DEFAULT_FROZEN_ROOT
    / "pilots"
    / "HCC1395_DORADO_chr6"
    / "ranking_bootstrap0"
)
DEFAULT_HYPERCUBE_SOURCE = (
    DEFAULT_FROZEN_ROOT
    / "release_contract"
    / "code_snapshot"
    / "research"
    / "20260716_read_linked_hypercube_exact_likelihood_validation"
    / "scripts"
    / "hypercube_exact.py"
)
DEFAULT_OUTPUT = (
    REPO_ROOT
    / "research"
    / "20260718_solver_methyl_edge_probe"
    / "results"
    / "m2_solver_route_census_receipt.json"
)

UNIT_FIELDS = (
    "dataset",
    "chrom",
    "threshold",
    "component_basis",
    "phase_set",
    "component_id",
    "family",
    "structural_exact_pattern_minread",
)
IDENTITY_FIELDS = UNIT_FIELDS[:-1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"Expected boolean TSV field, received {value!r}")
    return normalized == "true"


def unit_key(row: Mapping[str, str]) -> Tuple[str, ...]:
    return tuple(row[field] for field in UNIT_FIELDS)


def identity_key(key: Sequence[str]) -> Tuple[str, ...]:
    return tuple(key[:-1])


def iter_truncated_tsv_gzip(
    path: Path,
) -> Tuple[Iterator[Dict[str, str]], Dict[str, Any]]:
    """Yield complete records and retain an audit status for a truncated gzip."""
    status: Dict[str, Any] = {
        "path": str(path.resolve()),
        "complete_records": 0,
        "malformed_records_skipped": 0,
        "truncation_detected": False,
        "terminal_error": "",
    }

    def records() -> Iterator[Dict[str, str]]:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            try:
                header_line = handle.readline()
            except (EOFError, gzip.BadGzipFile) as error:
                status["truncation_detected"] = True
                status["terminal_error"] = type(error).__name__
                return
            if not header_line:
                return
            header = next(csv.reader([header_line], delimiter="\t"))
            while True:
                try:
                    line = handle.readline()
                except (EOFError, gzip.BadGzipFile) as error:
                    status["truncation_detected"] = True
                    status["terminal_error"] = type(error).__name__
                    break
                if not line:
                    break
                values = next(csv.reader([line], delimiter="\t"))
                if len(values) != len(header):
                    status["malformed_records_skipped"] += 1
                    continue
                status["complete_records"] += 1
                yield dict(zip(header, values))

    return records(), status


def load_frozen_hypercube(path: Path):
    module_name = "_intersubmod_frozen_m2_hypercube_for_census"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import frozen hypercube source: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def proxy_ops(m: int, q: int) -> int:
    if m < 1:
        return 0
    return (3**q) * (1 << m) + (1 << q) * m * (1 << (m - 1))


def dense_bytes(m: int, q: int, bytes_per_cell: int) -> int:
    return (1 << q) * (1 << m) * bytes_per_cell


def median_int(values: Iterable[int]) -> float | None:
    materialized = list(values)
    return statistics.median(materialized) if materialized else None


def structural_key(
    *,
    k: int,
    universe_mask: int,
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
) -> Tuple[Any, ...]:
    """Equality authority when solver source and run parameters are fixed."""
    return (
        int(k),
        int(universe_mask),
        tuple(sorted(full_patterns)),
        tuple(sorted(partial_patterns)),
    )


def classify_runtime(row: Mapping[str, str]) -> str:
    candidate = parse_bool(row["candidate_generation_invoked"])
    likelihood = parse_bool(row["likelihood_fit_invoked"])
    if candidate and likelihood:
        return "complete"
    if candidate and not likelihood:
        return "tail"
    if not candidate and not likelihood:
        return "no_solver"
    return "invalid_likelihood_without_candidate"


def build_receipt(
    ranking_root: Path,
    hypercube_source: Path,
    *,
    exact_k_max: int,
    max_vertex_sets: int,
    solver_time_limit_seconds: float,
    max_proxy_ops: int,
    max_dense_mib: int,
    bytes_per_dense_cell: int,
) -> Dict[str, Any]:
    paths = {
        "runtime": ranking_root / "m2_unit_runtime_diagnostics.tsv.gz",
        "patterns": ranking_root / "m2_symbolic_pattern_counts.tsv.gz",
        "candidates": ranking_root / "m2_compressed_vertex_set_candidates.tsv.gz",
        "hypercube_source": hypercube_source,
    }
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)

    frozen = load_frozen_hypercube(hypercube_source)

    patterns_by_unit: Dict[Tuple[str, ...], Dict[str, Any]] = {}
    pattern_rows, pattern_status = iter_truncated_tsv_gzip(paths["patterns"])
    for row in pattern_rows:
        key = unit_key(row)
        entry = patterns_by_unit.setdefault(
            key,
            {"k": int(row["k"]), "retained_patterns": set()},
        )
        if entry["k"] != int(row["k"]):
            raise AssertionError(f"Inconsistent k for {key}")
        if parse_bool(row["structural_retained"]):
            entry["retained_patterns"].add(row["pattern_rax"])

    candidates_by_unit: Dict[Tuple[str, ...], set[str]] = defaultdict(set)
    candidate_rows, candidate_status = iter_truncated_tsv_gzip(paths["candidates"])
    for row in candidate_rows:
        candidates_by_unit[unit_key(row)].add(row["vertex_set_id"])

    runtime_rows: list[Dict[str, Any]] = []
    runtime_stream, runtime_status = iter_truncated_tsv_gzip(paths["runtime"])
    for order, row in enumerate(runtime_stream):
        key = unit_key(row)
        if key not in patterns_by_unit:
            raise AssertionError(f"Runtime unit lacks symbolic-pattern rows: {key}")
        pattern_entry = patterns_by_unit[key]
        retained = sorted(pattern_entry["retained_patterns"])
        full = tuple(pattern for pattern in retained if "X" not in pattern)
        partial = tuple(pattern for pattern in retained if "X" in pattern)
        k = int(pattern_entry["k"])
        universe_mask = int(frozen.observed_alt_universe(full, partial))
        mandatory = {0, *(frozen.parse_full(pattern) for pattern in full)}
        groups = tuple(frozen.SymbolicPattern.from_string(pattern) for pattern in partial)
        reduction = frozen._reduce_partial_groups(groups, universe_mask, mandatory)
        required = mandatory | set(reduction.forced_vertices)
        nonroot_required = sorted(required - {0})
        single_terminal_weight = None
        if len(nonroot_required) == 1 and not reduction.groups:
            single_terminal_weight = bin(int(nonroot_required[0])).count("1")

        runtime_rows.append(
            {
                "order": order,
                "key": key,
                "identity": identity_key(key),
                "minread": int(key[-1]),
                "classification": classify_runtime(row),
                "candidate_seconds": float(
                    row["candidate_generation_elapsed_seconds"]
                ),
                "likelihood_seconds": float(row["likelihood_fit_elapsed_seconds"]),
                "unit_total_seconds": float(row["unit_total_elapsed_seconds"]),
                "k": k,
                "m": bin(int(universe_mask)).count("1"),
                "universe_mask": universe_mask,
                "full_patterns": full,
                "partial_patterns": partial,
                "input_partial_groups": len(partial),
                "reduced_partial_groups": len(reduction.groups),
                "forced_vertices": len(reduction.forced_vertices),
                "total_obligations": (
                    len(mandatory - {0})
                    + len(reduction.forced_vertices)
                    + len(reduction.groups)
                ),
                "structural_key": structural_key(
                    k=k,
                    universe_mask=universe_mask,
                    full_patterns=full,
                    partial_patterns=partial,
                ),
                "candidate_vertex_sets": len(candidates_by_unit.get(key, set())),
                "single_terminal_weight": single_terminal_weight,
            }
        )

    def summarize_minread(minread: int) -> Dict[str, Any]:
        rows = [row for row in runtime_rows if row["minread"] == minread]
        by_class = Counter(row["classification"] for row in rows)
        complete = [row for row in rows if row["classification"] == "complete"]
        tails = [row for row in rows if row["classification"] == "tail"]
        no_solver = [row for row in rows if row["classification"] == "no_solver"]
        invalid = [
            row
            for row in rows
            if row["classification"] == "invalid_likelihood_without_candidate"
        ]
        no_solver_reasons = Counter(
            "m_zero"
            if row["m"] == 0
            else "m_gt_exact_limit"
            if row["m"] > exact_k_max
            else "other"
            for row in no_solver
        )

        total_candidate = math.fsum(row["candidate_seconds"] for row in rows)
        complete_candidate = math.fsum(
            row["candidate_seconds"] for row in complete
        )
        tail_candidate = math.fsum(row["candidate_seconds"] for row in tails)
        total_likelihood = math.fsum(row["likelihood_seconds"] for row in rows)
        total_unit = math.fsum(row["unit_total_seconds"] for row in rows)
        complete_unit = math.fsum(row["unit_total_seconds"] for row in complete)
        tail_unit = math.fsum(row["unit_total_seconds"] for row in tails)
        no_solver_unit = math.fsum(row["unit_total_seconds"] for row in no_solver)

        v_counts = Counter(row["candidate_vertex_sets"] for row in complete)
        high_v = [
            row for row in complete if row["candidate_vertex_sets"] >= 100
        ]
        hard_tail = tails + high_v
        hard_tail_seconds = math.fsum(
            row["candidate_seconds"] for row in hard_tail
        )

        factorial_complete = [
            row for row in complete if row["single_terminal_weight"] is not None
        ]
        factorial_mismatches = [
            {
                "unit_key": list(row["key"]),
                "weight": row["single_terminal_weight"],
                "expected": math.factorial(row["single_terminal_weight"]),
                "observed": row["candidate_vertex_sets"],
            }
            for row in factorial_complete
            if math.factorial(row["single_terminal_weight"])
            != row["candidate_vertex_sets"]
        ]
        factorial_tails = [
            row for row in tails if row["single_terminal_weight"] is not None
        ]

        def cache_scope(scope_rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
            seen: set[Tuple[Any, ...]] = set()
            saved_seconds = 0.0
            key_groups: Dict[
                Tuple[Any, ...], list[Dict[str, Any]]
            ] = defaultdict(list)
            for scope_row in scope_rows:
                key_groups[scope_row["structural_key"]].append(scope_row)
                if scope_row["structural_key"] in seen:
                    saved_seconds += scope_row["candidate_seconds"]
                else:
                    seen.add(scope_row["structural_key"])
            largest_key_rows = max(
                key_groups.values(),
                key=lambda group: (
                    len(group),
                    math.fsum(item["candidate_seconds"] for item in group),
                ),
                default=[],
            )
            total_seconds = math.fsum(
                scope_row["candidate_seconds"] for scope_row in scope_rows
            )
            return {
                "invoked_units": len(scope_rows),
                "distinct_exact_keys": len(key_groups),
                "call_reduction_fraction": (
                    1.0 - len(key_groups) / len(scope_rows)
                    if scope_rows
                    else 0.0
                ),
                "sequential_repeat_units": len(scope_rows) - len(key_groups),
                "sequential_repeat_candidate_seconds": saved_seconds,
                "sequential_repeat_candidate_time_fraction": (
                    saved_seconds / total_seconds if total_seconds else 0.0
                ),
                "largest_key_repetitions": len(largest_key_rows),
                "largest_key_candidate_seconds": math.fsum(
                    item["candidate_seconds"] for item in largest_key_rows
                ),
                "largest_key_sha256": (
                    semantic_sha256(largest_key_rows[0]["structural_key"])
                    if largest_key_rows
                    else ""
                ),
            }

        max_dense_bytes = max_dense_mib * 1024 * 1024
        over_cap = [row for row in no_solver if row["m"] > exact_k_max]
        adaptive_eligible = [
            row
            for row in over_cap
            if proxy_ops(row["m"], row["total_obligations"])
            <= max_proxy_ops
            and dense_bytes(
                row["m"],
                row["total_obligations"],
                bytes_per_dense_cell,
            )
            <= max_dense_bytes
        ]
        candidate_units_for_minread = {
            key for key in candidates_by_unit if int(key[-1]) == minread
        }
        complete_units_for_minread = {row["key"] for row in complete}

        return {
            "minread": minread,
            "runtime_units": len(rows),
            "classification": dict(sorted(by_class.items())),
            "invalid_likelihood_without_candidate": len(invalid),
            "no_solver_reasons": dict(sorted(no_solver_reasons.items())),
            "candidate_seconds": {
                "all": total_candidate,
                "complete": complete_candidate,
                "tail": tail_candidate,
                "tail_fraction": (
                    tail_candidate / total_candidate if total_candidate else 0.0
                ),
            },
            "likelihood_seconds_all": total_likelihood,
            "unit_total_seconds": {
                "all": total_unit,
                "complete": complete_unit,
                "tail": tail_unit,
                "no_solver": no_solver_unit,
            },
            "group_reduction": {
                "all_input_partial_groups": sum(
                    row["input_partial_groups"] for row in rows
                ),
                "all_reduced_partial_groups": sum(
                    row["reduced_partial_groups"] for row in rows
                ),
                "tail_input_partial_groups": sum(
                    row["input_partial_groups"] for row in tails
                ),
                "tail_reduced_partial_groups": sum(
                    row["reduced_partial_groups"] for row in tails
                ),
            },
            "tail_structure": {
                "units": len(tails),
                "effective_m_min": min((row["m"] for row in tails), default=None),
                "effective_m_median": median_int(row["m"] for row in tails),
                "effective_m_max": max((row["m"] for row in tails), default=None),
                "reduced_q_median": median_int(
                    row["reduced_partial_groups"] for row in tails
                ),
                "total_obligations_median": median_int(
                    row["total_obligations"] for row in tails
                ),
            },
            "candidate_vertex_sets": {
                "candidate_rows": sum(
                    row["candidate_vertex_sets"] for row in complete
                ),
                "v_equal_1_units": v_counts.get(1, 0),
                "v_gt_1_units": len(complete) - v_counts.get(1, 0),
                "v_ge_100_units": len(high_v),
                "maximum_v": max(v_counts, default=0),
                "candidate_unit_keys_match_complete_runtime": (
                    candidate_units_for_minread == complete_units_for_minread
                ),
                "candidate_unit_key_symmetric_difference": len(
                    candidate_units_for_minread ^ complete_units_for_minread
                ),
                "v_ge_100_candidate_time_fraction_of_complete": (
                    math.fsum(row["candidate_seconds"] for row in high_v)
                    / complete_candidate
                    if complete_candidate
                    else 0.0
                ),
            },
            "combined_hard_tail": {
                "definition": "candidate-incomplete tails plus complete units with V>=100",
                "units": len(hard_tail),
                "unit_fraction": len(hard_tail) / len(rows) if rows else 0.0,
                "candidate_seconds": hard_tail_seconds,
                "candidate_time_fraction": (
                    hard_tail_seconds / total_candidate
                    if total_candidate
                    else 0.0
                ),
            },
            "single_terminal_factorial": {
                "complete_units": len(factorial_complete),
                "weight_distribution": dict(
                    sorted(
                        Counter(
                            row["single_terminal_weight"]
                            for row in factorial_complete
                        ).items()
                    )
                ),
                "mismatches": factorial_mismatches,
                "tail_units": [
                    {
                        "unit_key": list(row["key"]),
                        "weight": row["single_terminal_weight"],
                        "analytic_vertex_sets": math.factorial(
                            row["single_terminal_weight"]
                        ),
                        "exceeds_max_vertex_sets": (
                            math.factorial(row["single_terminal_weight"])
                            > max_vertex_sets
                        ),
                    }
                    for row in factorial_tails
                ],
            },
            "structural_cache_opportunity": {
                "complete": cache_scope(complete),
                "tail": cache_scope(tails),
                "combined": cache_scope(complete + tails),
                "incomplete_semantics": (
                    "Tail keys cannot be stored as complete results; use "
                    "checkpoint/resume, same-key single scheduling, or an "
                    "explicitly incomplete memo."
                ),
            },
            "adaptive_m_relaxation": {
                "m_gt_exact_limit_units": len(over_cap),
                "eligible_units": len(adaptive_eligible),
                "gate": {
                    "max_proxy_ops": max_proxy_ops,
                    "max_dense_mib": max_dense_mib,
                    "bytes_per_dense_cell": bytes_per_dense_cell,
                },
                "eligible_cases": [
                    {
                        "unit_key": list(row["key"]),
                        "m": row["m"],
                        "q_reduced_partial": row["reduced_partial_groups"],
                        "q_total_obligations": row["total_obligations"],
                        "proxy_ops": proxy_ops(
                            row["m"], row["total_obligations"]
                        ),
                        "dense_bytes": dense_bytes(
                            row["m"],
                            row["total_obligations"],
                            bytes_per_dense_cell,
                        ),
                    }
                    for row in adaptive_eligible
                ],
            },
        }

    minread1 = summarize_minread(1)
    minread2 = summarize_minread(2)
    tails_by_minread = {
        minread: {
            row["identity"]
            for row in runtime_rows
            if row["minread"] == minread and row["classification"] == "tail"
        }
        for minread in (1, 2)
    }
    overlap = tails_by_minread[1] & tails_by_minread[2]
    union = tails_by_minread[1] | tails_by_minread[2]
    tail_rows_by_minread = {
        minread: [
            row
            for row in runtime_rows
            if row["minread"] == minread and row["classification"] == "tail"
        ]
        for minread in (1, 2)
    }
    tail_keys_by_minread = {
        minread: {row["structural_key"] for row in rows}
        for minread, rows in tail_rows_by_minread.items()
    }
    tail_key_by_identity = {
        minread: {row["identity"]: row["structural_key"] for row in rows}
        for minread, rows in tail_rows_by_minread.items()
    }

    checks = {
        "runtime_minread1_7058": minread1["runtime_units"] == 7058,
        "runtime_minread2_7058": minread2["runtime_units"] == 7058,
        "minread1_classification_expected": minread1["classification"]
        == {"complete": 4852, "no_solver": 2173, "tail": 33},
        "minread2_classification_expected": minread2["classification"]
        == {"complete": 4435, "no_solver": 2590, "tail": 33},
        "minread1_group_reduction_expected": minread1["group_reduction"]
        ["all_input_partial_groups"]
        == 12561
        and minread1["group_reduction"]["all_reduced_partial_groups"] == 1173,
        "minread1_tail_group_reduction_expected": minread1["group_reduction"]
        ["tail_input_partial_groups"]
        == 1039
        and minread1["group_reduction"]["tail_reduced_partial_groups"] == 127,
        "minread1_factorial_complete_zero_mismatch": not minread1[
            "single_terminal_factorial"
        ]["mismatches"],
        "minread1_factorial_complete_expected": minread1[
            "single_terminal_factorial"
        ]["complete_units"]
        == 4380,
        "minread1_factorial_tail_expected": len(
            minread1["single_terminal_factorial"]["tail_units"]
        )
        == 3,
        "minread1_factorial_tails_exceed_cap": all(
            row["exceeds_max_vertex_sets"]
            for row in minread1["single_terminal_factorial"]["tail_units"]
        ),
        "minread1_overcap_expected": minread1["adaptive_m_relaxation"]
        ["m_gt_exact_limit_units"]
        == 31,
        "minread1_adaptive_eligible_expected": minread1[
            "adaptive_m_relaxation"
        ]["eligible_units"]
        == 2,
        "minread1_candidate_units_match": minread1["candidate_vertex_sets"]
        ["candidate_unit_keys_match_complete_runtime"],
        "minread2_candidate_units_match": minread2["candidate_vertex_sets"]
        ["candidate_unit_keys_match_complete_runtime"],
        "tail_overlap_expected": len(overlap) == 20 and len(union) == 46,
        "tail_cross_minread_keys_are_distinct": not (
            tail_keys_by_minread[1] & tail_keys_by_minread[2]
        )
        and all(
            tail_key_by_identity[1][identity]
            != tail_key_by_identity[2][identity]
            for identity in overlap
        ),
    }

    source_files = {
        name: {
            "path": str(path.resolve()),
            "sha256": sha256(path),
        }
        for name, path in paths.items()
    }
    source_files["runtime"]["parse_status"] = runtime_status
    source_files["patterns"]["parse_status"] = pattern_status
    source_files["candidates"]["parse_status"] = candidate_status

    return {
        "schema_name": "intersubmod.m2_solver_route_census_receipt",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "A_EXPLORATORY_PILOT",
        "goals": ["G3", "G4", "G5"],
        "verdict": (
            "PASS_DIAGNOSTIC_CENSUS"
            if all(checks.values())
            else "FAIL_DIAGNOSTIC_CENSUS"
        ),
        "scope": {
            "dataset": "HCC1395_DORADO",
            "chrom": "chr6",
            "bootstrap": 0,
            "eligibility": "DIAGNOSTIC_ONLY_TRUNCATED_FROZEN_V2",
            "not_claimed": [
                "complete ranking result",
                "cross-sample performance",
                "topology or clone/subclone proportion",
            ],
        },
        "solver_contract": {
            "universe_mode": "observed_alt",
            "exact_k_max": exact_k_max,
            "max_vertex_sets": max_vertex_sets,
            "solver_time_limit_seconds": solver_time_limit_seconds,
            "structural_key_authority_when_fixed": [
                "k",
                "structural_alt_universe_mask",
                "sorted full_patterns",
                "sorted partial_patterns",
            ],
            "full_cache_key_must_also_bind": [
                "cache schema",
                "solver source SHA-256",
                "universe mode",
                "exact_k_max",
                "max_vertex_sets",
                "time limit",
            ],
        },
        "source_files": source_files,
        "minread1": minread1,
        "minread2": minread2,
        "tail_stability": {
            "minread1_tail_units": len(tails_by_minread[1]),
            "minread2_tail_units": len(tails_by_minread[2]),
            "identity_intersection": len(overlap),
            "identity_union": len(union),
            "overlap_fraction_of_each_tail_set": (
                len(overlap) / len(tails_by_minread[1])
                if tails_by_minread[1]
                else 0.0
            ),
            "jaccard": len(overlap) / len(union) if union else 0.0,
            "minread1_unique_structural_keys": len(tail_keys_by_minread[1]),
            "minread2_unique_structural_keys": len(tail_keys_by_minread[2]),
            "structural_key_intersection": len(
                tail_keys_by_minread[1] & tail_keys_by_minread[2]
            ),
            "same_identity_and_same_structural_key": sum(
                tail_key_by_identity[1][identity]
                == tail_key_by_identity[2][identity]
                for identity in overlap
            ),
            "combined_candidate_seconds": (
                minread1["candidate_seconds"]["all"]
                + minread2["candidate_seconds"]["all"]
            ),
            "combined_tail_candidate_seconds": (
                minread1["candidate_seconds"]["tail"]
                + minread2["candidate_seconds"]["tail"]
            ),
            "combined_tail_candidate_time_fraction": (
                minread1["candidate_seconds"]["tail"]
                + minread2["candidate_seconds"]["tail"]
            )
            / (
                minread1["candidate_seconds"]["all"]
                + minread2["candidate_seconds"]["all"]
            ),
            "independence_warning": (
                "minread=1 and minread=2 are sensitivity settings on the same "
                "frozen chr6 pilot, not independent samples. Recurrent unit "
                "identity does not imply the same structural solver input."
            ),
        },
        "checks": checks,
        "all_pass": all(checks.values()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking-root", type=Path, default=DEFAULT_RANKING_ROOT)
    parser.add_argument(
        "--hypercube-source", type=Path, default=DEFAULT_HYPERCUBE_SOURCE
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--exact-k-max", type=int, default=12)
    parser.add_argument("--max-vertex-sets", type=int, default=256)
    parser.add_argument("--solver-time-limit-seconds", type=float, default=30.0)
    parser.add_argument("--max-proxy-ops", type=int, default=50_000_000)
    parser.add_argument("--max-dense-mib", type=int, default=512)
    parser.add_argument("--bytes-per-dense-cell", type=int, default=16)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = build_receipt(
        args.ranking_root.resolve(),
        args.hypercube_source.resolve(),
        exact_k_max=args.exact_k_max,
        max_vertex_sets=args.max_vertex_sets,
        solver_time_limit_seconds=args.solver_time_limit_seconds,
        max_proxy_ops=args.max_proxy_ops,
        max_dense_mib=args.max_dense_mib,
        bytes_per_dense_cell=args.bytes_per_dense_cell,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        "PASS"
        if receipt["all_pass"]
        else "FAIL",
        f"minread1_units={receipt['minread1']['runtime_units']}",
        f"minread1_tail={receipt['minread1']['classification'].get('tail', 0)}",
        f"hard_tail={receipt['minread1']['combined_hard_tail']['units']}",
        f"adaptive_eligible={receipt['minread1']['adaptive_m_relaxation']['eligible_units']}",
        f"output={args.output.resolve()}",
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
