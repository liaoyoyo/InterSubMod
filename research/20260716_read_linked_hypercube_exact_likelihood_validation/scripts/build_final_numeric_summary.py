#!/usr/bin/env python3
"""Build a fail-closed presentation-stage numeric summary for full M2 evidence.

This script does not alter or replace the frozen exact-11 release verifier.  It
accepts only an already passing independent-verification receipt, authenticates
the terminal extraction/ranking receipts and every child receipt again, then
derives presentation counts from the authenticated child aggregates and the
canonical candidate/runtime tables.

Important semantic boundary:
  * HP1/HP2 and bridge thresholds are separate unit-evaluation strata.  Their
    sums are not deduplicated biological-region counts.
  * h* is the minimum extra-state objective reconstructed from vertex roles;
    it is not a hidden-clone count.
  * tied x topology is a coarse-class comparison over tied winning vertex
    sets.  Exact parent-edge topology is not recoverable from the canonical
    candidate table and is therefore emitted as null with a reason.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import stat
from array import array
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping, Sequence


DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
AUTOSOMES = tuple(f"chr{value}" for value in range(1, 23))
PRIMARY_BASES = ("PS_HP1", "PS_HP2")
RUNTIME_METRICS = (
    "candidate_generation_elapsed_seconds",
    "likelihood_fit_elapsed_seconds",
    "unit_total_elapsed_seconds",
)
RUNTIME_COLUMNS = (
    "dataset",
    "chrom",
    "threshold",
    "component_basis",
    "phase_set",
    "component_id",
    "family",
    "structural_exact_pattern_minread",
    "structural_minread_role",
    "candidate_generation_invoked",
    "likelihood_fit_invoked",
    *RUNTIME_METRICS,
)
CANDIDATE_COLUMNS = (
    "unit_key",
    "dataset",
    "chrom",
    "component_id",
    "threshold",
    "hp_family",
    "ps",
    "candidate_id",
    "vertex_set_id",
    "vertex_states",
    "vertex_roles",
    "parent_choice_count",
    "profile_log_likelihood",
    "relative_log_likelihood",
    "mixture_weights_pi",
    "winner_status",
    "tie_group",
    "coarse_topology_class",
    "candidate_set_complete",
)
SUM_FIELDS = (
    "n_component_hp_units",
    "n_components",
    "molecule_component_projections",
    "informative_scoring_molecules",
    "all_x_excluded_molecules",
    "structural_retained_molecules",
    "below_minread_scoring_molecules",
    "bq_scoring_molecules",
    "non_ra_cells_marginalized",
    "raw_tree_candidates_T_complete_units",
    "distinct_vertex_sets_V_complete_units",
    "solver_complete_units",
    "solver_incomplete_or_not_run_units",
    "quality_primary_unique_vertex_units",
    "quality_primary_tied_vertex_units",
    "rank_abstain_units",
    "fixed_error_grid_stable_units",
    "fixed_error_grid_evaluated_units",
    "conditional_candidate_ranking_bootstrap_complete_units",
    "conditional_candidate_ranking_bootstrap_not_run_units",
    "conditional_candidate_ranking_bootstrap_nonconverged_units",
    "coarse_topology_class_unique_units",
    "coarse_topology_multiple_class_units",
    "parent_edge_assignment_unique_units",
    "exact_topology_proven_unique_units",
    "topology_evaluated_units",
    "topology_class_inclusion_counts_denominator",
    "k_component_sites_total",
    "k_observed_alt_active_total",
    "k_scoring_alt_observed_total",
    "not_structural_alt_active_sites_total",
    "structural_partial_pattern_groups",
    "partial_group_coverage_denominator",
    "partial_groups_covered",
    "partial_groups_unsatisfied",
)
COUNTER_FIELDS = (
    "selection_status_counts",
    "candidate_generation_status_counts",
    "k_route_counts",
    "projected_call_class_counts",
    "conditional_candidate_ranking_bootstrap_status_counts",
    "topology_class_inclusion_counts",
    "coarse_topology_unique_class_counts",
    "coarse_topology_ambiguous_class_set_counts",
    "topology_derivation_status_counts",
    "exact_topology_uniqueness_status_counts",
)
TOPOLOGY_CLASS_ORDER = {
    "single-only": 0,
    "sister-only": 1,
    "direct-only": 2,
    "sister+direct": 3,
}
VERTEX_ROLE_VOCABULARY = {
    "root",
    "full-observed",
    "partial-compatible",
    "connector",
}
TREE_VERTEX_BUCKETS = (
    "T_EQ_1_V_EQ_1",
    "T_GT_1_V_EQ_1",
    "T_GT_1_V_GT_1",
)
TREE_VERTEX_BUCKET_DEFINITIONS = {
    "T_EQ_1_V_EQ_1": (
        "exactly one distinct globally minimum-extra-state vertex set (V=1) and "
        "exactly one feasible parent-edge tree assignment (T=1)"
    ),
    "T_GT_1_V_EQ_1": (
        "exactly one distinct globally minimum-extra-state vertex set (V=1), but "
        "more than one feasible parent-edge tree assignment (T>1)"
    ),
    "T_GT_1_V_GT_1": (
        "more than one distinct globally minimum-extra-state vertex set (V>1); "
        "because every vertex set has at least one parent-edge assignment, T>=V>1"
    ),
}
TREE_VERTEX_PARTITION_DEFINITION = (
    "Mutually exclusive per solver-complete candidate-table unit. V is the number of "
    "candidate rows (distinct optimal vertex sets), and T is the sum of "
    "parent_choice_count over those rows. The invariant T>=V>=1 leaves exactly these "
    "three buckets. This is a structural candidate partition, not a clone-count claim."
)


class SummaryError(RuntimeError):
    """Raised when persisted evidence cannot support a presentation number."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SummaryError(message)


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _stable_regular_file(path: Path, label: str) -> os.stat_result:
    try:
        before = path.lstat()
    except OSError as exc:
        raise SummaryError(f"{label}: unavailable file {path}: {exc}") from exc
    require(stat.S_ISREG(before.st_mode), f"{label}: not a regular file: {path}")
    require(not path.is_symlink(), f"{label}: symlink is forbidden: {path}")
    return before


def _assert_unchanged(path: Path, before: os.stat_result, label: str) -> None:
    after = path.lstat()
    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    require(identity_before == identity_after, f"{label}: file changed while reading: {path}")


def verify_sha256_sidecar(path: Path, label: str) -> str:
    before = _stable_regular_file(path, label)
    digest = sha256_path(path)
    _assert_unchanged(path, before, label)
    sidecar = path.with_name(f"{path.name}.sha256")
    side_before = _stable_regular_file(sidecar, f"{label} sidecar")
    try:
        fields = sidecar.read_text(encoding="ascii", errors="strict").strip().split()
    except OSError as exc:
        raise SummaryError(f"{label}: cannot read sidecar: {exc}") from exc
    _assert_unchanged(sidecar, side_before, f"{label} sidecar")
    require(
        fields == [digest, path.name],
        f"{label}: SHA-256 sidecar does not exactly bind {path.name}",
    )
    return digest


def load_authenticated_json(
    path: Path, schema_name: str, schema_version: str, label: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    path = path.resolve()
    digest = verify_sha256_sidecar(path, label)
    before = path.lstat()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeError) as exc:
        raise SummaryError(f"{label}: invalid JSON: {exc}") from exc
    _assert_unchanged(path, before, label)
    require(isinstance(payload, dict), f"{label}: JSON root is not an object")
    require(payload.get("schema_name") == schema_name, f"{label}: schema_name mismatch")
    require(payload.get("schema_version") == schema_version, f"{label}: schema_version mismatch")
    require(payload.get("all_pass") is True, f"{label}: all_pass is not true")
    checks = payload.get("checks") or {}
    require(
        isinstance(checks, dict) and checks and all(value is True for value in checks.values()),
        f"{label}: one or more declared checks are not true",
    )
    return payload, {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": digest,
    }


def resolve_recorded_path(raw: Any, base: Path, label: str) -> Path:
    require(isinstance(raw, str) and raw, f"{label}: missing recorded path")
    path = Path(raw)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def verify_recorded_identity(identity: Mapping[str, Any], base: Path, label: str) -> dict[str, Any]:
    path = resolve_recorded_path(identity.get("path"), base, label)
    before = _stable_regular_file(path, label)
    digest = sha256_path(path)
    _assert_unchanged(path, before, label)
    require(identity.get("sha256") == digest, f"{label}: recorded SHA-256 mismatch")
    require(int(identity.get("size_bytes", -1)) == before.st_size, f"{label}: recorded size mismatch")
    return {"path": str(path), "size_bytes": before.st_size, "sha256": digest}


def expected_pairs(datasets: Sequence[str], chromosomes: Sequence[str]) -> list[tuple[str, str]]:
    return [(dataset, chrom) for dataset in datasets for chrom in chromosomes]


def result_index(
    rows: Sequence[Mapping[str, Any]], datasets: Sequence[str], chromosomes: Sequence[str], label: str
) -> dict[tuple[str, str], Mapping[str, Any]]:
    expected = set(expected_pairs(datasets, chromosomes))
    index: dict[tuple[str, str], Mapping[str, Any]] = {}
    for row in rows:
        key = (str(row.get("dataset")), str(row.get("chrom")))
        require(key not in index, f"{label}: duplicate task {key}")
        index[key] = row
    require(set(index) == expected, f"{label}: task Cartesian scope mismatch")
    return index


def integer_mapping(payload: Mapping[str, Any]) -> dict[str, int]:
    output = {
        str(key): int(value)
        for key, value in payload.items()
        if isinstance(value, int) and not isinstance(value, bool)
    }
    require(all(value >= 0 for value in output.values()), "negative count in integer mapping")
    return output


def verify_extraction_count_conservation(counts: Mapping[str, Any], label: str) -> dict[str, bool]:
    def value(key: str) -> int:
        return int(counts.get(key, 0))

    checks = {
        "raw_alignment_class": value("raw_overlapping_alignments")
        == value("alignment_class_primary")
        + value("alignment_class_secondary")
        + value("alignment_class_supplementary")
        + value("alignment_class_unmapped"),
        "raw_filter_funnel": value("raw_overlapping_alignments")
        == value("excluded_by_flag")
        + value("mapq_rejected_after_flag")
        + value("canonical_eligible_alignments"),
        "eligible_rows": value("canonical_eligible_alignments")
        == value("molecule_sparse_rows_written")
        == value("unique_molecule_ids"),
        "fixed_calls_within_sparse_calls": value("fixed_ra_calls")
        <= value("site_call_rows_sparse"),
        "alt_calls_within_fixed_calls": value("alt_calls") <= value("fixed_ra_calls"),
        "sidecar_exact_matches": value("sidecar_exact_matches")
        == value("canonical_eligible_alignments"),
        "sidecar_missing_zero": value("sidecar_missing") == 0,
    }
    failed = sorted(key for key, passed in checks.items() if not passed)
    require(not failed, f"{label}: extraction count conservation failed: {failed}")
    return checks


def ratio(numerator: int | float, denominator: int | float, denominator_label: str) -> dict[str, Any]:
    require(numerator >= 0, f"negative numerator: {denominator_label}")
    require(denominator >= 0, f"negative denominator: {denominator_label}")
    if denominator == 0:
        return {
            "numerator": numerator,
            "denominator": denominator,
            "value": None,
            "percent": None,
            "denominator_label": denominator_label,
            "reason": "denominator_is_zero",
        }
    value = numerator / denominator
    return {
        "numerator": numerator,
        "denominator": denominator,
        "value": value,
        "percent": value * 100.0,
        "denominator_label": denominator_label,
        "reason": None,
    }


def null_metric(reason: str) -> dict[str, Any]:
    return {"value": None, "reason": reason}


def merge_numeric_tree(target: MutableMapping[str, Any], source: Mapping[str, Any], path: str = "") -> None:
    for key, value in source.items():
        current = f"{path}/{key}" if path else str(key)
        if key == "definitions":
            if key in target:
                require(target[key] == value, f"definition drift at {current}")
            else:
                target[key] = dict(value)
        elif isinstance(value, Mapping):
            child = target.setdefault(key, {})
            require(isinstance(child, dict), f"numeric-tree type drift at {current}")
            merge_numeric_tree(child, value, current)
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            target[key] = target.get(key, 0) + value
        elif key in target:
            require(target[key] == value, f"scalar drift at {current}")
        else:
            target[key] = value


def verify_component_cell(cell: Mapping[str, Any], label: str) -> None:
    distribution = {
        int(key): int(value)
        for key, value in (cell.get("k_component_sites_distribution") or {}).items()
    }
    n_components = int(cell.get("n_components", 0))
    require(
        all(key >= 1 and value >= 0 for key, value in distribution.items()),
        f"{label}: invalid k distribution domain",
    )
    require(
        all(
            value >= 0
            for value in cell.values()
            if isinstance(value, int) and not isinstance(value, bool)
        ),
        f"{label}: negative component count",
    )
    require(sum(distribution.values()) == n_components, f"{label}: k distribution sum mismatch")
    require(int(cell.get("n_singletons_k1", 0)) == distribution.get(1, 0), f"{label}: k=1 mismatch")
    require(
        int(cell.get("n_multisite_k_gt1", 0))
        == sum(value for key, value in distribution.items() if key > 1),
        f"{label}: k>1 mismatch",
    )
    require(
        int(cell.get("max_k_component_sites", 0)) == (max(distribution) if distribution else 0),
        f"{label}: max k mismatch",
    )
    require(
        int(cell.get("n_active_site_memberships", -1))
        == sum(key * value for key, value in distribution.items()),
        f"{label}: active-site membership mass mismatch",
    )
    require(
        int(cell.get("n_effective_k_route_deferred_to_ranker", -1)) == n_components,
        f"{label}: deferred effective-k route count mismatch",
    )


def new_component_accumulator() -> dict[str, Any]:
    return {"sums": Counter(), "max": Counter(), "distribution": Counter()}


def add_component_cell(target: MutableMapping[str, Any], cell: Mapping[str, Any], label: str) -> None:
    verify_component_cell(cell, label)
    for key, value in cell.items():
        if isinstance(value, int) and not isinstance(value, bool):
            if key in {"max_k_component_sites", "max_k"}:
                target["max"][key] = max(target["max"].get(key, 0), int(value))
            else:
                target["sums"][key] += int(value)
    target["distribution"].update(
        {
            int(key): int(value)
            for key, value in (cell.get("k_component_sites_distribution") or {}).items()
        }
    )


def freeze_component_cell(target: Mapping[str, Any]) -> dict[str, Any]:
    distribution = {str(key): value for key, value in sorted(target["distribution"].items())}
    return {
        **dict(target["sums"]),
        **dict(target["max"]),
        "k_component_sites_distribution": distribution,
        "k_distribution": distribution,
    }


def presentation_component_cell(cell: Mapping[str, Any]) -> dict[str, Any]:
    n_components = int(cell.get("n_components", 0))
    k1 = int(cell.get("n_singletons_k1", 0))
    k_gt1 = int(cell.get("n_multisite_k_gt1", 0))
    return {
        "n_components": n_components,
        "k_equals_1": k1,
        "k_greater_than_1": k_gt1,
        "k_distribution": dict(cell.get("k_component_sites_distribution") or {}),
        "max_k_component_sites": int(cell.get("max_k_component_sites", 0)),
        "active_site_membership_mass": int(cell.get("n_active_site_memberships", 0)),
        "k_equals_1_share": ratio(k1, n_components, "components_in_same_dataset_basis_threshold"),
        "k_greater_than_1_share": ratio(
            k_gt1, n_components, "components_in_same_dataset_basis_threshold"
        ),
    }


def new_rank_accumulator() -> dict[str, Any]:
    return {
        "sums": Counter(),
        "counters": {field: Counter() for field in COUNTER_FIELDS},
        "partial_pattern_funnel": {},
    }


def add_rank_cell(
    target: MutableMapping[str, Any], summary: Mapping[str, Any], partial: Mapping[str, Any] | None
) -> None:
    values = {field: summary.get(field, 0) for field in SUM_FIELDS}
    require(
        all(
            isinstance(value, int) and not isinstance(value, bool) and value >= 0
            for value in values.values()
        ),
        "rank cell contains a non-integer or negative count",
    )
    target["sums"].update({field: int(value) for field, value in values.items()})
    for field in COUNTER_FIELDS:
        target["counters"][field].update(integer_mapping(summary.get(field) or {}))
    if partial is not None:
        merge_numeric_tree(target["partial_pattern_funnel"], partial)


def freeze_rank_cell(target: Mapping[str, Any]) -> dict[str, Any]:
    return {
        **dict(target["sums"]),
        **{
            field: dict(sorted(target["counters"][field].items()))
            for field in COUNTER_FIELDS
        },
        "partial_pattern_funnel": target["partial_pattern_funnel"],
    }


def verify_rank_conservation(cell: Mapping[str, Any], label: str) -> dict[str, bool]:
    units = int(cell.get("n_component_hp_units", 0))
    checks = {
        "selection_partition": sum((cell.get("selection_status_counts") or {}).values()) == units,
        "unique_tied_abstain_partition": (
            int(cell.get("quality_primary_unique_vertex_units", 0))
            + int(cell.get("quality_primary_tied_vertex_units", 0))
            + int(cell.get("rank_abstain_units", 0))
            == units
        ),
        "solver_partition": (
            int(cell.get("solver_complete_units", 0))
            + int(cell.get("solver_incomplete_or_not_run_units", 0))
            == units
        ),
        "projection_funnel": int(cell.get("molecule_component_projections", 0))
        == int(cell.get("informative_scoring_molecules", 0))
        + int(cell.get("all_x_excluded_molecules", 0)),
        "scoring_funnel": int(cell.get("informative_scoring_molecules", 0))
        == int(cell.get("structural_retained_molecules", 0))
        + int(cell.get("below_minread_scoring_molecules", 0)),
        "bq_equals_informative": int(cell.get("bq_scoring_molecules", 0))
        == int(cell.get("informative_scoring_molecules", 0)),
        "T_not_less_than_V": int(cell.get("raw_tree_candidates_T_complete_units", 0))
        >= int(cell.get("distinct_vertex_sets_V_complete_units", 0)),
        "topology_partition": int(cell.get("coarse_topology_class_unique_units", 0))
        + int(cell.get("coarse_topology_multiple_class_units", 0))
        == int(cell.get("topology_evaluated_units", 0)),
        "partial_coverage": int(cell.get("partial_groups_covered", 0))
        + int(cell.get("partial_groups_unsatisfied", 0))
        == int(cell.get("partial_group_coverage_denominator", 0)),
        "partial_unsatisfied_zero": int(cell.get("partial_groups_unsatisfied", 0)) == 0,
        "k_route_partition": sum((cell.get("k_route_counts") or {}).values()) == units,
        "effective_k_mass": int(cell.get("k_component_sites_total", 0))
        == int(cell.get("k_observed_alt_active_total", 0))
        + int(cell.get("not_structural_alt_active_sites_total", 0)),
        "scoring_and_structural_k_order": int(cell.get("k_observed_alt_active_total", 0))
        <= int(cell.get("k_scoring_alt_observed_total", 0))
        <= int(cell.get("k_component_sites_total", 0)),
        "topology_denominator": int(cell.get("topology_class_inclusion_counts_denominator", 0))
        == int(cell.get("topology_evaluated_units", 0)),
        "coarse_unique_class_counts": sum(
            (cell.get("coarse_topology_unique_class_counts") or {}).values()
        )
        == int(cell.get("coarse_topology_class_unique_units", 0)),
        "coarse_ambiguous_class_sets": sum(
            (cell.get("coarse_topology_ambiguous_class_set_counts") or {}).values()
        )
        == int(cell.get("coarse_topology_multiple_class_units", 0)),
        "exact_parent_edge_uniqueness": int(cell.get("exact_topology_proven_unique_units", 0))
        == int(cell.get("parent_edge_assignment_unique_units", 0))
        <= int(cell.get("topology_evaluated_units", 0)),
    }
    failed = sorted(key for key, value in checks.items() if not value)
    require(not failed, f"{label}: conservation failed: {failed}")
    return checks


def compare_rank_raw(actual: Mapping[str, Any], expected: Mapping[str, Any], label: str) -> None:
    for field in SUM_FIELDS:
        require(
            int(actual.get(field, 0)) == int(expected.get(field, 0)),
            f"{label}/{field}: aggregate mismatch",
        )
    for field in COUNTER_FIELDS:
        require(
            (actual.get(field) or {}) == (expected.get(field) or {}),
            f"{label}/{field}: counter mismatch",
        )
    require(
        (actual.get("partial_pattern_funnel") or {})
        == (expected.get("partial_pattern_funnel") or {}),
        f"{label}/partial_pattern_funnel: mismatch",
    )


def runtime_summary(values: Iterable[float]) -> dict[str, Any]:
    data = [float(value) for value in values]
    n = len(data)
    if n == 0:
        return {"n": 0, "sum": 0.0, "max": None, "p50": None, "p95": None, "p99": None}
    require(all(math.isfinite(value) and value >= 0.0 for value in data), "invalid runtime value")
    total = math.fsum(data)
    data.sort()

    def nearest(probability: float) -> float:
        return data[min(n - 1, math.ceil(probability * n) - 1)]

    return {
        "n": n,
        "sum": total,
        "max": data[-1],
        "p50": nearest(0.50),
        "p95": nearest(0.95),
        "p99": nearest(0.99),
    }


def new_runtime_accumulator() -> dict[str, Any]:
    return {
        "all": {metric: array("d") for metric in RUNTIME_METRICS},
        "invoked": {
            metric: array("d")
            for metric in (
                "candidate_generation_elapsed_seconds",
                "likelihood_fit_elapsed_seconds",
            )
        },
    }


def add_runtime_value(
    target: MutableMapping[str, Any], parsed: Mapping[str, float], invocation: Mapping[str, bool]
) -> None:
    for metric in RUNTIME_METRICS:
        target["all"][metric].append(parsed[metric])
    for metric, invoked in invocation.items():
        if invoked:
            target["invoked"][metric].append(parsed[metric])


def freeze_runtime(target: Mapping[str, Any]) -> dict[str, Any]:
    n = len(target["all"][RUNTIME_METRICS[0]])
    require(all(len(target["all"][metric]) == n for metric in RUNTIME_METRICS), "runtime row mismatch")
    return {
        "n_unit_evaluations": n,
        "metrics": {
            metric: runtime_summary(target["all"][metric]) for metric in RUNTIME_METRICS
        },
        "metrics_when_invoked": {
            metric: runtime_summary(values) for metric, values in target["invoked"].items()
        },
        "interpretation": (
            "process-local monotonic wall-clock diagnostic; environment/load dependent, "
            "not a biological result"
        ),
        "memory_model": (
            "runtime TSV rows are streamed; exact nearest-rank quantiles retain compact "
            "float64 value arrays (O(N) values) but do not retain full TSV rows"
        ),
    }


def parse_topology_classes(raw: str, label: str) -> tuple[str, ...]:
    try:
        values = json.loads(raw)
    except (TypeError, json.JSONDecodeError) as exc:
        raise SummaryError(f"{label}: invalid coarse topology JSON") from exc
    require(
        isinstance(values, list)
        and values
        and all(isinstance(value, str) and value for value in values),
        f"{label}: coarse topology class list is empty or malformed",
    )
    require(
        set(values).issubset(TOPOLOGY_CLASS_ORDER),
        f"{label}: unknown coarse topology class",
    )
    return tuple(sorted(set(values), key=TOPOLOGY_CLASS_ORDER.__getitem__))


def minimum_extra_states(row: Mapping[str, str], label: str) -> int:
    try:
        states = json.loads(row["vertex_states"])
        roles = json.loads(row["vertex_roles"])
    except (TypeError, json.JSONDecodeError) as exc:
        raise SummaryError(f"{label}: malformed vertex state/role JSON") from exc
    require(
        isinstance(states, dict)
        and isinstance(roles, dict)
        and states
        and set(states) == set(roles),
        f"{label}: vertex state/role key mismatch",
    )
    state_lengths: set[int] = set()
    for key, state in states.items():
        require(
            isinstance(key, str)
            and key.isdigit()
            and isinstance(state, str)
            and state
            and set(state).issubset({"R", "A"}),
            f"{label}: malformed vertex state",
        )
        require(
            int(key) == sum((symbol == "A") << bit for bit, symbol in enumerate(state)),
            f"{label}: vertex bitmask/state mismatch",
        )
        state_lengths.add(len(state))
    require(len(state_lengths) == 1, f"{label}: inconsistent vertex state lengths")
    for key, role_values in roles.items():
        require(
            isinstance(role_values, list)
            and role_values
            and all(isinstance(value, str) and value for value in role_values),
            f"{label}: malformed vertex role list",
        )
        role_set = set(role_values)
        require(
            len(role_set) == len(role_values)
            and role_set.issubset(VERTEX_ROLE_VOCABULARY),
            f"{label}: duplicate or unknown vertex role",
        )
        require(("root" in role_set) == (int(key) == 0), f"{label}: root role mismatch")
        if "connector" in role_set:
            require(role_set == {"connector"}, f"{label}: connector role overlap")
        elif "root" not in role_set and "full-observed" not in role_set:
            require(
                "partial-compatible" in role_set,
                f"{label}: extra state has neither partial-compatible nor connector role",
            )
    return sum(
        "root" not in role_values and "full-observed" not in role_values
        for role_values in roles.values()
    )


def new_candidate_accumulator() -> dict[str, Any]:
    return {
        "n_units": 0,
        "n_candidate_vertex_sets_V": 0,
        "n_parent_edge_trees_T": 0,
        "unique_first": 0,
        "tied_first": 0,
        "solver_complete_optimizer_abstain": 0,
        "tied_topology_consistent": 0,
        "tied_topology_inconsistent": 0,
        "topology_evaluated_units": 0,
        "coarse_topology_class_unique_units": 0,
        "coarse_topology_multiple_class_units": 0,
        "parent_edge_assignment_unique_units": 0,
        "exact_topology_proven_unique_units": 0,
        "tree_vertex_partition_counts": Counter(),
        "h_star_distribution": Counter(),
        "coarse_topology_class_inclusion_counts": Counter(),
        "coarse_topology_unique_class_counts": Counter(),
        "coarse_topology_ambiguous_class_set_counts": Counter(),
    }


def finish_candidate_unit(
    rows: Sequence[Mapping[str, str]],
    datasets: Sequence[str],
    chromosomes: Sequence[str],
    by_stratum: MutableMapping[tuple[str, str, str], dict[str, Any]],
    all_candidates: MutableMapping[str, Any],
) -> None:
    if not rows:
        return
    first = rows[0]
    unit_key = first["unit_key"]
    dataset = first["dataset"]
    chrom = first["chrom"]
    family = first["hp_family"]
    basis = f"PS_HP{family}"
    threshold = str(int(first["threshold"]))
    require(dataset in datasets and chrom in chromosomes, f"candidate {unit_key}: scope mismatch")
    require(basis in PRIMARY_BASES, f"candidate {unit_key}: unsupported HP family")
    require(first["ps"] not in {"", ".", "NA", "None", "UNKNOWN", "unknown"}, f"candidate {unit_key}: missing PS")
    require(
        all(
            row["unit_key"] == unit_key
            and row["dataset"] == dataset
            and row["chrom"] == chrom
            and row["component_id"] == first["component_id"]
            and row["hp_family"] == family
            and row["ps"] == first["ps"]
            and str(int(row["threshold"])) == threshold
            for row in rows
        ),
        f"candidate {unit_key}: within-unit scope drift",
    )
    require(
        len({row["candidate_id"] for row in rows}) == len(rows)
        and len({row["vertex_set_id"] for row in rows}) == len(rows),
        f"candidate {unit_key}: duplicate candidate or vertex-set identifier",
    )
    try:
        parent_choice_counts = [int(row["parent_choice_count"]) for row in rows]
    except (KeyError, TypeError, ValueError) as exc:
        raise SummaryError(f"candidate {unit_key}: invalid parent_choice_count") from exc
    require(
        all(value >= 1 for value in parent_choice_counts),
        f"candidate {unit_key}: parent_choice_count must be positive",
    )
    n_vertex_sets_V = len(rows)
    n_parent_edge_trees_T = sum(parent_choice_counts)
    require(
        n_parent_edge_trees_T >= n_vertex_sets_V >= 1,
        f"candidate {unit_key}: T>=V>=1 invariant failed",
    )
    if n_parent_edge_trees_T == 1 and n_vertex_sets_V == 1:
        tree_vertex_bucket = "T_EQ_1_V_EQ_1"
    elif n_parent_edge_trees_T > 1 and n_vertex_sets_V == 1:
        tree_vertex_bucket = "T_GT_1_V_EQ_1"
    elif n_parent_edge_trees_T > 1 and n_vertex_sets_V > 1:
        tree_vertex_bucket = "T_GT_1_V_GT_1"
    else:  # The explicit branch makes the three-state contract fail closed.
        raise SummaryError(f"candidate {unit_key}: unsupported T/V partition state")
    require(
        all(row["candidate_set_complete"].lower() == "true" for row in rows),
        f"candidate {unit_key}: incomplete candidate set in canonical table",
    )
    h_values = {minimum_extra_states(row, unit_key) for row in rows}
    require(len(h_values) == 1, f"candidate {unit_key}: h* differs across candidates")
    h_star = next(iter(h_values))
    statuses = [row["winner_status"].upper() for row in rows]
    unique_rows = [row for row in rows if row["winner_status"].upper() == "UNIQUE_WINNER"]
    tied_rows = [row for row in rows if row["winner_status"].upper() == "TIED_WINNER"]
    abstain_rows = [row for row in rows if row["winner_status"].upper() == "ABSTAIN_UNIT_OPTIMIZER"]
    allowed = {"UNIQUE_WINNER", "TIED_WINNER", "NON_WINNER", "ABSTAIN_UNIT_OPTIMIZER"}
    require(set(statuses).issubset(allowed), f"candidate {unit_key}: unknown winner status")
    if abstain_rows:
        require(len(abstain_rows) == len(rows), f"candidate {unit_key}: mixed abstain partition")
        outcome = "abstain"
    elif unique_rows:
        require(len(unique_rows) == 1 and not tied_rows, f"candidate {unit_key}: invalid unique partition")
        outcome = "unique"
    else:
        require(len(tied_rows) >= 2, f"candidate {unit_key}: invalid tied partition")
        outcome = "tied"

    targets = (by_stratum.setdefault((dataset, basis, threshold), new_candidate_accumulator()), all_candidates)
    for target in targets:
        target["n_units"] += 1
        target["n_candidate_vertex_sets_V"] += n_vertex_sets_V
        target["n_parent_edge_trees_T"] += n_parent_edge_trees_T
        target["tree_vertex_partition_counts"][tree_vertex_bucket] += 1
        target["h_star_distribution"][h_star] += 1
        if outcome == "unique":
            target["unique_first"] += 1
        elif outcome == "tied":
            target["tied_first"] += 1
        else:
            target["solver_complete_optimizer_abstain"] += 1

    winner_rows = unique_rows or tied_rows
    if winner_rows:
        union: set[str] = set()
        winner_parent_choices = 0
        for row in winner_rows:
            classes = parse_topology_classes(row["coarse_topology_class"], unit_key)
            union.update(classes)
            winner_parent_choices += int(row["parent_choice_count"])
        ordered_union = tuple(sorted(union, key=TOPOLOGY_CLASS_ORDER.__getitem__))
        for target in targets:
            target["topology_evaluated_units"] += 1
            target["coarse_topology_class_inclusion_counts"].update(ordered_union)
            if len(ordered_union) == 1:
                target["coarse_topology_class_unique_units"] += 1
                target["coarse_topology_unique_class_counts"][ordered_union[0]] += 1
            else:
                target["coarse_topology_multiple_class_units"] += 1
                target["coarse_topology_ambiguous_class_set_counts"][
                    "|".join(ordered_union)
                ] += 1
            if winner_parent_choices == 1:
                target["parent_edge_assignment_unique_units"] += 1
                target["exact_topology_proven_unique_units"] += 1
        if outcome == "tied":
            key = "tied_topology_consistent" if len(union) == 1 else "tied_topology_inconsistent"
            for target in targets:
                target[key] += 1


def verify_tree_vertex_partition(
    payload: Mapping[str, Any], expected_denominator: int, label: str
) -> dict[str, bool]:
    require(isinstance(payload, Mapping), f"{label}: partition is not an object")
    counts = payload.get("counts")
    shares = payload.get("shares")
    require(isinstance(counts, Mapping), f"{label}: counts is not an object")
    require(isinstance(shares, Mapping), f"{label}: shares is not an object")
    require(set(counts) == set(TREE_VERTEX_BUCKETS), f"{label}: count bucket keys mismatch")
    require(set(shares) == set(TREE_VERTEX_BUCKETS), f"{label}: share bucket keys mismatch")
    require(
        all(isinstance(counts[key], int) and not isinstance(counts[key], bool) for key in TREE_VERTEX_BUCKETS),
        f"{label}: bucket count must be an integer",
    )
    require(all(int(counts[key]) >= 0 for key in TREE_VERTEX_BUCKETS), f"{label}: negative bucket count")
    require(
        isinstance(payload.get("denominator"), int)
        and not isinstance(payload.get("denominator"), bool)
        and int(payload["denominator"]) == expected_denominator,
        f"{label}: denominator mismatch",
    )
    require(
        sum(int(counts[key]) for key in TREE_VERTEX_BUCKETS) == expected_denominator,
        f"{label}: bucket counts do not conserve denominator",
    )
    require(
        payload.get("definition") == TREE_VERTEX_PARTITION_DEFINITION
        and payload.get("bucket_definitions") == TREE_VERTEX_BUCKET_DEFINITIONS,
        f"{label}: definition mismatch",
    )
    for key in TREE_VERTEX_BUCKETS:
        require(
            shares[key] == ratio(int(counts[key]), expected_denominator, "solver_complete_candidate_units"),
            f"{label}: share mismatch for {key}",
        )
    return {
        "exact_three_bucket_domain": True,
        "denominator_matches_solver_complete_units": True,
        "bucket_counts_conserve": True,
        "shares_recompute_exactly": True,
        "definitions_are_canonical": True,
    }


def build_tree_vertex_partition(target: Mapping[str, Any], n_units: int) -> dict[str, Any]:
    raw_counts = target.get("tree_vertex_partition_counts") or {}
    require(
        isinstance(raw_counts, Mapping) and set(raw_counts).issubset(TREE_VERTEX_BUCKETS),
        "candidate tree/vertex partition contains an unknown bucket",
    )
    counts = {key: int(raw_counts.get(key, 0)) for key in TREE_VERTEX_BUCKETS}
    payload = {
        "counts": counts,
        "denominator": n_units,
        "shares": {
            key: ratio(counts[key], n_units, "solver_complete_candidate_units")
            for key in TREE_VERTEX_BUCKETS
        },
        "definition": TREE_VERTEX_PARTITION_DEFINITION,
        "bucket_definitions": dict(TREE_VERTEX_BUCKET_DEFINITIONS),
    }
    verify_tree_vertex_partition(payload, n_units, "candidate tree/vertex partition")
    return payload


def freeze_candidate_accumulator(target: Mapping[str, Any]) -> dict[str, Any]:
    n_units = int(target["n_units"])
    n_vertex_sets_V = int(target["n_candidate_vertex_sets_V"])
    n_parent_edge_trees_T = int(target["n_parent_edge_trees_T"])
    unique = int(target["unique_first"])
    tied = int(target["tied_first"])
    optimizer_abstain = int(target["solver_complete_optimizer_abstain"])
    topology_evaluated = int(target["topology_evaluated_units"])
    coarse_unique = int(target["coarse_topology_class_unique_units"])
    coarse_multiple = int(target["coarse_topology_multiple_class_units"])
    if n_units == 0:
        require(
            n_vertex_sets_V == 0 and n_parent_edge_trees_T == 0,
            "empty candidate accumulator has nonzero T or V",
        )
    else:
        require(
            n_parent_edge_trees_T >= n_vertex_sets_V >= n_units,
            "candidate aggregate T>=V>=number of solver-complete units invariant failed",
        )
    require(unique + tied + optimizer_abstain == n_units, "candidate outcome partition mismatch")
    require(
        int(target["tied_topology_consistent"])
        + int(target["tied_topology_inconsistent"])
        == tied,
        "tied coarse-topology partition mismatch",
    )
    require(
        coarse_unique + coarse_multiple == topology_evaluated,
        "candidate coarse-topology partition mismatch",
    )
    require(
        sum(target["h_star_distribution"].values()) == n_units,
        "candidate h* distribution mismatch",
    )
    tree_vertex_partition = build_tree_vertex_partition(target, n_units)
    return {
        "n_solver_complete_candidate_units": n_units,
        "n_candidate_vertex_sets_V": n_vertex_sets_V,
        "n_parent_edge_trees_T": n_parent_edge_trees_T,
        "tree_vertex_partition": tree_vertex_partition,
        "unique_first": unique,
        "tied_first": tied,
        "solver_complete_optimizer_abstain": optimizer_abstain,
        "unique_first_share": ratio(unique, n_units, "solver_complete_candidate_units"),
        "tied_first_share": ratio(tied, n_units, "solver_complete_candidate_units"),
        "optimizer_abstain_share": ratio(
            optimizer_abstain, n_units, "solver_complete_candidate_units"
        ),
        "tied_by_coarse_topology": {
            "consistent": int(target["tied_topology_consistent"]),
            "inconsistent": int(target["tied_topology_inconsistent"]),
            "denominator": tied,
            "consistent_share": ratio(
                int(target["tied_topology_consistent"]), tied, "tied_first_solver_complete_units"
            ),
            "inconsistent_share": ratio(
                int(target["tied_topology_inconsistent"]), tied, "tied_first_solver_complete_units"
            ),
            "definition": (
                "consistent iff the union of coarse topology classes across all tied winning "
                "vertex sets has cardinality one; this does not claim exact parent-edge uniqueness"
            ),
        },
        "topology": {
            "evaluated_units": topology_evaluated,
            "coarse_class_unique_units": coarse_unique,
            "coarse_class_multiple_units": coarse_multiple,
            "coarse_unique_class_counts": dict(
                sorted(target["coarse_topology_unique_class_counts"].items())
            ),
            "coarse_ambiguous_class_set_counts": dict(
                sorted(target["coarse_topology_ambiguous_class_set_counts"].items())
            ),
            "parent_edge_assignment_unique_units": int(
                target["parent_edge_assignment_unique_units"]
            ),
            "exact_topology_proven_unique_units": int(
                target["exact_topology_proven_unique_units"]
            ),
        },
        "h_star_distribution": {
            str(key): value for key, value in sorted(target["h_star_distribution"].items())
        },
        "coarse_topology_class_inclusion_counts": dict(
            sorted(target["coarse_topology_class_inclusion_counts"].items())
        ),
    }


def parse_candidate_table(
    metadata: Mapping[str, Any], ranking_path: Path, datasets: Sequence[str], chromosomes: Sequence[str]
) -> tuple[dict[tuple[str, str, str], dict[str, Any]], dict[str, Any], dict[str, Any]]:
    require(metadata.get("schema_version") == "2.0.0", "candidate metadata schema mismatch")
    require(metadata.get("format") == "tsv.gz", "candidate metadata format mismatch")
    require(tuple(metadata.get("columns") or ()) == CANDIDATE_COLUMNS, "candidate metadata columns mismatch")
    path = resolve_recorded_path(metadata.get("path"), ranking_path.parent, "candidate table")
    identity = verify_recorded_identity(metadata, ranking_path.parent, "candidate table")
    semantic = hashlib.sha256()
    by_stratum: dict[tuple[str, str, str], dict[str, Any]] = {}
    overall = new_candidate_accumulator()
    n_rows = 0
    n_units = 0
    previous_unit: str | None = None
    current: list[dict[str, str]] = []
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require(tuple(reader.fieldnames or ()) == CANDIDATE_COLUMNS, "candidate table header mismatch")
            for row in reader:
                n_rows += 1
                semantic.update(
                    json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
                    + b"\n"
                )
                unit_key = row["unit_key"]
                if previous_unit is None:
                    previous_unit = unit_key
                elif unit_key != previous_unit:
                    require(unit_key > previous_unit, "candidate unit sort order regression")
                    finish_candidate_unit(current, datasets, chromosomes, by_stratum, overall)
                    n_units += 1
                    current = []
                    previous_unit = unit_key
                current.append(dict(row))
            if current:
                finish_candidate_unit(current, datasets, chromosomes, by_stratum, overall)
                n_units += 1
    except (OSError, UnicodeError, csv.Error, gzip.BadGzipFile) as exc:
        raise SummaryError(f"candidate table read failed: {exc}") from exc
    require(n_rows == int(metadata.get("n_rows", -1)), "candidate n_rows mismatch")
    require(n_units == int(metadata.get("n_units", -1)), "candidate n_units mismatch")
    require(semantic.hexdigest() == metadata.get("semantic_sha256"), "candidate semantic SHA mismatch")
    identity["semantic_sha256"] = semantic.hexdigest()
    identity["n_rows"] = n_rows
    identity["n_units"] = n_units
    return by_stratum, overall, identity


def presentation_rank_cell(
    cell: Mapping[str, Any], candidate: Mapping[str, Any] | None
) -> dict[str, Any]:
    units = int(cell.get("n_component_hp_units", 0))
    complete = int(cell.get("solver_complete_units", 0))
    incomplete = int(cell.get("solver_incomplete_or_not_run_units", 0))
    unique = int(cell.get("quality_primary_unique_vertex_units", 0))
    tied = int(cell.get("quality_primary_tied_vertex_units", 0))
    abstain = int(cell.get("rank_abstain_units", 0))
    projections = int(cell.get("molecule_component_projections", 0))
    informative = int(cell.get("informative_scoring_molecules", 0))
    structural = int(cell.get("structural_retained_molecules", 0))
    below = int(cell.get("below_minread_scoring_molecules", 0))
    candidate_accumulator = candidate or new_candidate_accumulator()
    require(
        int(candidate_accumulator["n_units"]) == complete,
        "presentation candidate-table solver-complete unit binding mismatch",
    )
    require(
        int(candidate_accumulator["n_candidate_vertex_sets_V"])
        == int(cell.get("distinct_vertex_sets_V_complete_units", 0)),
        "presentation candidate-table V binding mismatch",
    )
    require(
        int(candidate_accumulator["n_parent_edge_trees_T"])
        == int(cell.get("raw_tree_candidates_T_complete_units", 0)),
        "presentation candidate-table T binding mismatch",
    )
    candidate_payload = freeze_candidate_accumulator(candidate_accumulator)
    return {
        "units": {
            "n_component_hp_unit_evaluations": units,
            "solver_complete": complete,
            "solver_incomplete_or_not_run": incomplete,
            "solver_complete_share": ratio(complete, units, "component_hp_unit_evaluations"),
        },
        "molecule_funnel": {
            "component_projections": projections,
            "informative_scoring_molecules": informative,
            "all_X_excluded_molecules": int(cell.get("all_x_excluded_molecules", 0)),
            "structural_retained_molecules": structural,
            "below_structural_minread_but_scored_molecules": below,
            "informative_share_of_projections": ratio(
                informative, projections, "molecule_component_projections"
            ),
            "structural_share_of_informative": ratio(
                structural, informative, "informative_scoring_molecules"
            ),
            "below_minread_share_of_informative": ratio(
                below, informative, "informative_scoring_molecules"
            ),
        },
        "partial_read_funnel": {
            "structural_partial_pattern_groups": int(
                cell.get("structural_partial_pattern_groups", 0)
            ),
            "coverage_denominator": int(cell.get("partial_group_coverage_denominator", 0)),
            "covered": int(cell.get("partial_groups_covered", 0)),
            "unsatisfied": int(cell.get("partial_groups_unsatisfied", 0)),
            "covered_share": ratio(
                int(cell.get("partial_groups_covered", 0)),
                int(cell.get("partial_group_coverage_denominator", 0)),
                "structural_partial_group_constraints",
            ),
            "full_detail": cell.get("partial_pattern_funnel") or {},
        },
        "candidate_structure": {
            "raw_parent_edge_trees_T_complete_units": int(
                cell.get("raw_tree_candidates_T_complete_units", 0)
            ),
            "distinct_optimal_vertex_sets_V_complete_units": int(
                cell.get("distinct_vertex_sets_V_complete_units", 0)
            ),
            "mean_T_per_solver_complete_unit": ratio(
                int(cell.get("raw_tree_candidates_T_complete_units", 0)),
                complete,
                "solver_complete_units",
            ),
            "mean_V_per_solver_complete_unit": ratio(
                int(cell.get("distinct_vertex_sets_V_complete_units", 0)),
                complete,
                "solver_complete_units",
            ),
            "candidate_table": candidate_payload,
        },
        "ranking_outcome": {
            "unique_first": unique,
            "tied_first": tied,
            "abstain_all_causes": abstain,
            "unique_share": ratio(unique, units, "component_hp_unit_evaluations"),
            "tied_share": ratio(tied, units, "component_hp_unit_evaluations"),
            "abstain_share": ratio(abstain, units, "component_hp_unit_evaluations"),
            "selection_status_counts": dict(cell.get("selection_status_counts") or {}),
        },
        "topology": {
            "evaluated_units": int(cell.get("topology_evaluated_units", 0)),
            "coarse_class_unique_units": int(
                cell.get("coarse_topology_class_unique_units", 0)
            ),
            "coarse_class_multiple_units": int(
                cell.get("coarse_topology_multiple_class_units", 0)
            ),
            "coarse_unique_class_counts": dict(
                cell.get("coarse_topology_unique_class_counts") or {}
            ),
            "coarse_ambiguous_class_set_counts": dict(
                cell.get("coarse_topology_ambiguous_class_set_counts") or {}
            ),
            "parent_edge_assignment_unique_units": int(
                cell.get("parent_edge_assignment_unique_units", 0)
            ),
            "exact_topology_proven_unique_units": int(
                cell.get("exact_topology_proven_unique_units", 0)
            ),
            "coarse_class_unique_share": ratio(
                int(cell.get("coarse_topology_class_unique_units", 0)),
                int(cell.get("topology_evaluated_units", 0)),
                "topology_evaluated_units",
            ),
            "coarse_class_multiple_share": ratio(
                int(cell.get("coarse_topology_multiple_class_units", 0)),
                int(cell.get("topology_evaluated_units", 0)),
                "topology_evaluated_units",
            ),
            "exact_topology_proven_unique_share": ratio(
                int(cell.get("exact_topology_proven_unique_units", 0)),
                int(cell.get("topology_evaluated_units", 0)),
                "topology_evaluated_units",
            ),
        },
        "effective_k": {
            "component_site_mass": int(cell.get("k_component_sites_total", 0)),
            "observed_ALT_active_mass": int(cell.get("k_observed_alt_active_total", 0)),
            "not_structural_ALT_active_mass": int(
                cell.get("not_structural_alt_active_sites_total", 0)
            ),
            "k_route_counts": dict(cell.get("k_route_counts") or {}),
            "observed_ALT_active_share_of_component_site_mass": ratio(
                int(cell.get("k_observed_alt_active_total", 0)),
                int(cell.get("k_component_sites_total", 0)),
                "component_site_mass",
            ),
            "not_structural_ALT_active_share_of_component_site_mass": ratio(
                int(cell.get("not_structural_alt_active_sites_total", 0)),
                int(cell.get("k_component_sites_total", 0)),
                "component_site_mass",
            ),
            "route_shares_of_unit_evaluations": {
                str(key): ratio(int(value), units, "component_hp_unit_evaluations")
                for key, value in sorted((cell.get("k_route_counts") or {}).items())
            },
        },
    }


def _nested_component(
    cells: Mapping[tuple[str, str], Mapping[str, Any]]
) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = defaultdict(dict)
    for (basis, threshold), value in sorted(cells.items()):
        output[basis][threshold] = freeze_component_cell(value)
    return dict(output)


def _nested_rank(cells: Mapping[tuple[str, str], Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = defaultdict(dict)
    for (basis, threshold), value in sorted(cells.items()):
        output[basis][threshold] = freeze_rank_cell(value)
    return dict(output)


def _lookup_nested(payload: Mapping[str, Any], basis: str, threshold: str, label: str) -> Mapping[str, Any]:
    try:
        value = payload[basis][threshold]
    except (KeyError, TypeError) as exc:
        raise SummaryError(f"{label}: missing {basis}/{threshold}") from exc
    require(isinstance(value, Mapping), f"{label}: invalid cell {basis}/{threshold}")
    return value


def build_summary(
    extraction_root: Path,
    ranking_root: Path,
    final_verification_path: Path,
    *,
    require_canonical_scope: bool = True,
) -> dict[str, Any]:
    extraction_root = extraction_root.resolve()
    ranking_root = ranking_root.resolve()
    final_verification_path = final_verification_path.resolve()
    final, final_identity = load_authenticated_json(
        final_verification_path,
        "intersubmod.m2_full_independent_verification",
        "1.0.0",
        "final independent verification",
    )
    scope = final.get("scope") or {}
    datasets = tuple(scope.get("datasets") or ())
    chromosomes = tuple(scope.get("chromosomes") or ())
    require(datasets and chromosomes, "final verification scope is empty")
    require(
        int(scope.get("expected_tasks", -1)) == len(datasets) * len(chromosomes),
        "final verification expected_tasks mismatch",
    )
    if require_canonical_scope:
        require(datasets == DATASETS, "final presentation requires exact seven-dataset scope")
        require(chromosomes == AUTOSOMES, "final presentation requires chr1-chr22 scope")
    release = final.get("release_binding") or {}
    require(release.get("validation_evidence_eligible") is True, "release is not validation-evidence eligible")
    require(
        ((release.get("deep_release_verification") or {}).get("all_pass")) is True,
        "release binding deep verification is not passing",
    )

    extraction_path = extraction_root / "full_extraction_receipt.json"
    ranking_path = ranking_root / "full_ranking_receipt.json"
    extraction, extraction_identity = load_authenticated_json(
        extraction_path,
        "intersubmod.m2_full_extraction_receipt",
        "1.2.0",
        "terminal full extraction",
    )
    ranking, ranking_identity = load_authenticated_json(
        ranking_path,
        "intersubmod.m2_full_ranking_receipt",
        "2.0.0",
        "terminal full ranking",
    )
    for payload, label in ((extraction, "extraction"), (ranking, "ranking")):
        payload_scope = payload.get("scope") or {}
        require(tuple(payload_scope.get("datasets") or ()) == datasets, f"{label} dataset scope mismatch")
        require(tuple(payload_scope.get("chromosomes") or ()) == chromosomes, f"{label} chromosome scope mismatch")
        require(int(payload_scope.get("expected_tasks", -1)) == len(datasets) * len(chromosomes), f"{label} task count mismatch")

    final_extraction = final.get("extraction") or {}
    final_ranking = final.get("ranking") or {}
    require(
        resolve_recorded_path(final_extraction.get("receipt_path"), final_verification_path.parent, "final extraction binding")
        == extraction_path
        and final_extraction.get("receipt_sha256") == extraction_identity["sha256"],
        "final verifier does not bind the supplied terminal extraction receipt",
    )
    require(
        resolve_recorded_path(final_ranking.get("receipt_path"), final_verification_path.parent, "final ranking binding")
        == ranking_path
        and final_ranking.get("receipt_sha256") == ranking_identity["sha256"],
        "final verifier does not bind the supplied terminal ranking receipt",
    )
    require(
        final_extraction.get("recomputed_aggregate") == extraction.get("aggregate"),
        "final verifier extraction recomputation differs from terminal aggregate",
    )

    expected = expected_pairs(datasets, chromosomes)
    extraction_results = result_index(extraction.get("results") or [], datasets, chromosomes, "extraction results")
    ranking_results = result_index(ranking.get("results") or [], datasets, chromosomes, "ranking results")
    extraction_final_tasks = result_index(final_extraction.get("task_index") or [], datasets, chromosomes, "final extraction task index")
    ranking_final_tasks = result_index(final_ranking.get("task_index") or [], datasets, chromosomes, "final ranking task index")

    extraction_counts_global: Counter[str] = Counter()
    extraction_counts_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    extraction_global_cells: dict[tuple[str, str], dict[str, Any]] = {}
    extraction_dataset_cells: dict[str, dict[tuple[str, str], dict[str, Any]]] = defaultdict(dict)
    extraction_child_sources: list[dict[str, Any]] = []

    for dataset, chrom in expected:
        result = extraction_results[(dataset, chrom)]
        require(result.get("status") in {"PASS", "REUSED_PASS"}, f"extraction task not passing: {dataset}/{chrom}")
        child_path = extraction_root / "samples" / dataset / chrom / "receipt.json"
        child, identity = load_authenticated_json(
            child_path,
            "intersubmod.lossless_read_linkage_chromosome_receipt",
            "1.2.0",
            f"extraction child {dataset}/{chrom}",
        )
        require(result.get("receipt") == child, f"extraction embedded child drift: {dataset}/{chrom}")
        child_scope = child.get("scope") or {}
        require(child_scope.get("dataset") == dataset and child_scope.get("chrom") == chrom, f"extraction child scope mismatch: {dataset}/{chrom}")
        require(extraction_final_tasks[(dataset, chrom)].get("child_receipt_sha256") == identity["sha256"], f"final extraction child SHA mismatch: {dataset}/{chrom}")
        extraction_child_sources.append({"dataset": dataset, "chrom": chrom, **identity})
        counts = integer_mapping(child.get("counts") or {})
        verify_extraction_count_conservation(counts, f"extraction child {dataset}/{chrom}")
        extraction_counts_global.update(counts)
        extraction_counts_dataset[dataset].update(counts)
        summaries = child.get("component_summary_by_linkage_basis") or {}
        require(summaries, f"extraction child has no component summary: {dataset}/{chrom}")
        for basis, thresholds in summaries.items():
            require(isinstance(thresholds, Mapping) and thresholds, f"invalid component thresholds: {dataset}/{chrom}/{basis}")
            for threshold, cell in thresholds.items():
                key = (str(basis), str(threshold))
                add_component_cell(
                    extraction_global_cells.setdefault(key, new_component_accumulator()),
                    cell,
                    f"{dataset}/{chrom}/{basis}/{threshold}",
                )
                add_component_cell(
                    extraction_dataset_cells[dataset].setdefault(key, new_component_accumulator()),
                    cell,
                    f"{dataset}/{chrom}/{basis}/{threshold}",
                )

    reconstructed_extraction_cells = _nested_component(extraction_global_cells)
    extraction_aggregate = extraction.get("aggregate") or {}
    require(
        reconstructed_extraction_cells
        == (extraction_aggregate.get("component_summary_by_linkage_basis") or {}),
        "component aggregate does not equal sum of authenticated extraction children",
    )
    require(dict(extraction_counts_global) == (extraction_aggregate.get("counts") or {}), "extraction count aggregate mismatch")
    for dataset in datasets:
        require(
            dict(extraction_counts_dataset[dataset])
            == (((extraction_aggregate.get("by_dataset") or {}).get(dataset) or {}).get("counts") or {}),
            f"extraction dataset count mismatch: {dataset}",
        )

    rank_global_cells: dict[tuple[str, str], dict[str, Any]] = {}
    rank_dataset_cells: dict[str, dict[tuple[str, str], dict[str, Any]]] = defaultdict(dict)
    input_global: Counter[str] = Counter()
    input_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    call_global: Counter[str] = Counter()
    call_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    hp_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    ranking_child_sources: list[dict[str, Any]] = []
    runtime_sources: list[dict[str, Any]] = []
    runtime_global = new_runtime_accumulator()
    runtime_dataset = {dataset: new_runtime_accumulator() for dataset in datasets}
    runtime_stratum: dict[tuple[str, str, str], dict[str, Any]] = {}

    for dataset, chrom in expected:
        result = ranking_results[(dataset, chrom)]
        require(result.get("status") in {"PASS", "REUSED_PASS"}, f"ranking task not passing: {dataset}/{chrom}")
        declared_child = result.get("rank_receipt") or {}
        child_path = ranking_root / "samples" / dataset / chrom / "receipt.json"
        child, identity = load_authenticated_json(
            child_path,
            "intersubmod.m2_symbolic_patterns_vertex_rank_receipt",
            "2.0.0",
            f"ranking child {dataset}/{chrom}",
        )
        require(resolve_recorded_path(declared_child.get("path"), ranking_path.parent, f"ranking child path {dataset}/{chrom}") == child_path, f"ranking child path binding mismatch: {dataset}/{chrom}")
        require(declared_child.get("sha256") == identity["sha256"], f"ranking child result SHA mismatch: {dataset}/{chrom}")
        require(ranking_final_tasks[(dataset, chrom)].get("rank_receipt_sha256") == identity["sha256"], f"final ranking child SHA mismatch: {dataset}/{chrom}")
        child_scope = child.get("scope") or {}
        require(child_scope.get("dataset") == [dataset] and child_scope.get("chrom") == [chrom], f"ranking child scope mismatch: {dataset}/{chrom}")
        ranking_child_sources.append({"dataset": dataset, "chrom": chrom, **identity})

        input_counts = child.get("input_counts") or {}
        input_global.update(integer_mapping(input_counts))
        input_dataset[dataset].update(integer_mapping(input_counts))
        calls = integer_mapping(input_counts.get("selected_sparse_call_code_counts") or {})
        call_global.update(calls)
        call_dataset[dataset].update(calls)
        hp_dataset[dataset].update(integer_mapping(input_counts.get("hp_family_rows") or {}))
        partial = child.get("partial_pattern_funnel_by_linkage_basis_threshold") or {}
        summaries = child.get("aggregate_by_linkage_basis_threshold") or {}
        require(summaries, f"ranking child has no primary aggregate: {dataset}/{chrom}")
        for basis, thresholds in summaries.items():
            for threshold, cell in thresholds.items():
                key = (str(basis), str(threshold))
                funnel = (partial.get(str(basis)) or {}).get(str(threshold))
                add_rank_cell(rank_global_cells.setdefault(key, new_rank_accumulator()), cell, funnel)
                add_rank_cell(
                    rank_dataset_cells[dataset].setdefault(key, new_rank_accumulator()), cell, funnel
                )

        runtime = child.get("runtime_diagnostics") or {}
        require(runtime.get("schema_name") == "intersubmod.m2_unit_runtime_diagnostics", f"runtime schema mismatch: {dataset}/{chrom}")
        require(runtime.get("schema_version") == "1.0.0", f"runtime schema version mismatch: {dataset}/{chrom}")
        runtime_name = runtime.get("per_unit_output")
        runtime_identity_declared = (child.get("outputs") or {}).get(runtime_name) or {}
        verified_runtime_identity = verify_recorded_identity(
            runtime_identity_declared, child_path.parent, f"runtime {dataset}/{chrom}"
        )
        runtime_sources.append({"dataset": dataset, "chrom": chrom, **verified_runtime_identity})
        child_runtime = new_runtime_accumulator()
        all_scope_runtime = new_runtime_accumulator()
        runtime_path = Path(verified_runtime_identity["path"])
        try:
            with gzip.open(runtime_path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                require(tuple(reader.fieldnames or ()) == RUNTIME_COLUMNS, f"runtime columns mismatch: {dataset}/{chrom}")
                for row in reader:
                    require(row["dataset"] == dataset and row["chrom"] == chrom, f"runtime row scope mismatch: {dataset}/{chrom}")
                    parsed = {metric: float(row[metric]) for metric in RUNTIME_METRICS}
                    require(
                        all(math.isfinite(value) and value >= 0.0 for value in parsed.values())
                        and parsed["candidate_generation_elapsed_seconds"]
                        + parsed["likelihood_fit_elapsed_seconds"]
                        <= parsed["unit_total_elapsed_seconds"] + 1e-9,
                        f"runtime segment conservation failure: {dataset}/{chrom}",
                    )
                    flags = {
                        "candidate_generation_elapsed_seconds": row["candidate_generation_invoked"],
                        "likelihood_fit_elapsed_seconds": row["likelihood_fit_invoked"],
                    }
                    require(set(flags.values()).issubset({"true", "false"}), f"runtime invocation flag invalid: {dataset}/{chrom}")
                    invocation = {key: value == "true" for key, value in flags.items()}
                    add_runtime_value(all_scope_runtime, parsed, invocation)
                    if row["structural_minread_role"] == "PRIMARY":
                        basis = row["component_basis"]
                        threshold = str(int(row["threshold"]))
                        require(basis in PRIMARY_BASES, f"primary runtime basis mismatch: {dataset}/{chrom}")
                        add_runtime_value(child_runtime, parsed, invocation)
                        add_runtime_value(runtime_dataset[dataset], parsed, invocation)
                        add_runtime_value(runtime_global, parsed, invocation)
                        add_runtime_value(
                            runtime_stratum.setdefault(
                                (dataset, basis, threshold), new_runtime_accumulator()
                            ),
                            parsed,
                            invocation,
                        )
        except (OSError, UnicodeError, csv.Error, gzip.BadGzipFile, ValueError) as exc:
            raise SummaryError(f"runtime read failed {dataset}/{chrom}: {exc}") from exc
        child_runtime_frozen = freeze_runtime(child_runtime)
        declared_primary = ((runtime.get("scopes") or {}).get("primary_unit_evaluations") or {})
        require(
            {
                "n_unit_evaluations": child_runtime_frozen["n_unit_evaluations"],
                **child_runtime_frozen["metrics"],
            }
            == declared_primary,
            f"runtime primary summary mismatch: {dataset}/{chrom}",
        )
        require(
            child_runtime_frozen["metrics_when_invoked"]
            == (runtime.get("primary_invoked_segment_scopes") or {}),
            f"runtime invoked summary mismatch: {dataset}/{chrom}",
        )
        all_frozen = freeze_runtime(all_scope_runtime)
        require(
            {"n_unit_evaluations": all_frozen["n_unit_evaluations"], **all_frozen["metrics"]}
            == (((runtime.get("scopes") or {}).get("all_structural_minread_unit_evaluations")) or {}),
            f"runtime all-minread summary mismatch: {dataset}/{chrom}",
        )

    recomputed_rank = _nested_rank(rank_global_cells)
    rank_aggregate = ranking.get("aggregate") or {}
    declared_rank = rank_aggregate.get("by_linkage_basis_threshold") or {}
    require(set(recomputed_rank) == set(declared_rank), "ranking global basis set mismatch")
    conservation_checks: list[dict[str, Any]] = []
    for basis, thresholds in recomputed_rank.items():
        require(set(thresholds) == set(declared_rank[basis]), f"ranking threshold set mismatch: {basis}")
        for threshold, cell in thresholds.items():
            compare_rank_raw(declared_rank[basis][threshold], cell, f"ranking/global/{basis}/{threshold}")
            conservation_checks.append(
                {
                    "scope": f"global/{basis}/{threshold}",
                    "checks": verify_rank_conservation(cell, f"global/{basis}/{threshold}"),
                }
            )
    final_recomputed_rank = ((final_ranking.get("recomputed") or {}).get("by_linkage_basis_threshold") or {})
    require(final_recomputed_rank == recomputed_rank, "final verifier ranking recomputation differs from authenticated children")
    require(dict(input_global) == (rank_aggregate.get("input_call_funnel") or {}), "ranking input funnel mismatch")
    require(dict(sorted(call_global.items())) == (rank_aggregate.get("input_sparse_call_code_counts") or {}), "ranking call-code aggregate mismatch")

    recomputed_rank_dataset: dict[str, dict[str, dict[str, Any]]] = {}
    for dataset in datasets:
        nested = _nested_rank(rank_dataset_cells[dataset])
        recomputed_rank_dataset[dataset] = nested
        declared_dataset = ((rank_aggregate.get("by_dataset") or {}).get(dataset) or {})
        declared_cells = declared_dataset.get("by_linkage_basis_threshold") or {}
        require(set(nested) == set(declared_cells), f"ranking dataset basis set mismatch: {dataset}")
        for basis, thresholds in nested.items():
            require(set(thresholds) == set(declared_cells[basis]), f"ranking dataset threshold set mismatch: {dataset}/{basis}")
            for threshold, cell in thresholds.items():
                compare_rank_raw(declared_cells[basis][threshold], cell, f"ranking/{dataset}/{basis}/{threshold}")
                conservation_checks.append(
                    {
                        "scope": f"{dataset}/{basis}/{threshold}",
                        "checks": verify_rank_conservation(cell, f"{dataset}/{basis}/{threshold}"),
                    }
                )
        require(dict(input_dataset[dataset]) == (declared_dataset.get("input_call_funnel") or {}), f"ranking dataset input funnel mismatch: {dataset}")
        require(dict(sorted(call_dataset[dataset].items())) == (declared_dataset.get("input_sparse_call_code_counts") or {}), f"ranking dataset call-code mismatch: {dataset}")
        require(dict(sorted(hp_dataset[dataset].items())) == (declared_dataset.get("input_hp_family_rows") or {}), f"ranking dataset HP-row mismatch: {dataset}")

    candidate_metadata = ranking.get("candidate_table") or {}
    final_candidate = final_ranking.get("candidate_table") or {}
    for key in ("path", "size_bytes", "sha256", "semantic_sha256", "n_rows", "n_units"):
        require(
            final_candidate.get(key) == candidate_metadata.get(key),
            f"final candidate-table binding differs from ranking terminal: {key}",
        )
    candidate_by_stratum, candidate_overall, candidate_identity = parse_candidate_table(
        candidate_metadata, ranking_path, datasets, chromosomes
    )
    require(
        candidate_identity["sha256"] == candidate_metadata.get("sha256")
        and candidate_identity["semantic_sha256"] == candidate_metadata.get("semantic_sha256"),
        "candidate table hash binding failed",
    )
    for dataset in datasets:
        for basis, thresholds in recomputed_rank_dataset[dataset].items():
            for threshold, cell in thresholds.items():
                candidate = candidate_by_stratum.get((dataset, basis, threshold), new_candidate_accumulator())
                complete = int(cell.get("solver_complete_units", 0))
                incomplete = int(cell.get("solver_incomplete_or_not_run_units", 0))
                require(int(candidate["n_units"]) == complete, f"candidate unit coverage mismatch: {dataset}/{basis}/{threshold}")
                require(int(candidate["n_candidate_vertex_sets_V"]) == int(cell.get("distinct_vertex_sets_V_complete_units", 0)), f"candidate V mismatch: {dataset}/{basis}/{threshold}")
                require(int(candidate["n_parent_edge_trees_T"]) == int(cell.get("raw_tree_candidates_T_complete_units", 0)), f"candidate T mismatch: {dataset}/{basis}/{threshold}")
                require(int(candidate["unique_first"]) == int(cell.get("quality_primary_unique_vertex_units", 0)), f"candidate unique mismatch: {dataset}/{basis}/{threshold}")
                require(int(candidate["tied_first"]) == int(cell.get("quality_primary_tied_vertex_units", 0)), f"candidate tied mismatch: {dataset}/{basis}/{threshold}")
                require(
                    int(candidate["solver_complete_optimizer_abstain"])
                    == int(cell.get("rank_abstain_units", 0)) - incomplete,
                    f"candidate optimizer-abstain mismatch: {dataset}/{basis}/{threshold}",
                )
                require(
                    int(candidate["tied_topology_consistent"])
                    + int(candidate["tied_topology_inconsistent"])
                    == int(candidate["tied_first"]),
                    f"tied topology partition mismatch: {dataset}/{basis}/{threshold}",
                )
                for candidate_key, cell_key in (
                    ("topology_evaluated_units", "topology_evaluated_units"),
                    ("coarse_topology_class_unique_units", "coarse_topology_class_unique_units"),
                    ("coarse_topology_multiple_class_units", "coarse_topology_multiple_class_units"),
                    ("parent_edge_assignment_unique_units", "parent_edge_assignment_unique_units"),
                    ("exact_topology_proven_unique_units", "exact_topology_proven_unique_units"),
                ):
                    require(
                        int(candidate[candidate_key]) == int(cell.get(cell_key, 0)),
                        f"candidate {candidate_key} mismatch: {dataset}/{basis}/{threshold}",
                    )
                for candidate_key, cell_key in (
                    ("coarse_topology_class_inclusion_counts", "topology_class_inclusion_counts"),
                    ("coarse_topology_unique_class_counts", "coarse_topology_unique_class_counts"),
                    (
                        "coarse_topology_ambiguous_class_set_counts",
                        "coarse_topology_ambiguous_class_set_counts",
                    ),
                ):
                    require(
                        dict(sorted(candidate[candidate_key].items()))
                        == dict(sorted((cell.get(cell_key) or {}).items())),
                        f"candidate {candidate_key} mismatch: {dataset}/{basis}/{threshold}",
                    )
    candidate_global_expected = {
        "n_units": sum(
            int(cell.get("solver_complete_units", 0))
            for thresholds in recomputed_rank.values()
            for cell in thresholds.values()
        ),
        "n_candidate_vertex_sets_V": sum(
            int(cell.get("distinct_vertex_sets_V_complete_units", 0))
            for thresholds in recomputed_rank.values()
            for cell in thresholds.values()
        ),
        "n_parent_edge_trees_T": sum(
            int(cell.get("raw_tree_candidates_T_complete_units", 0))
            for thresholds in recomputed_rank.values()
            for cell in thresholds.values()
        ),
    }
    for key, expected_value in candidate_global_expected.items():
        require(
            int(candidate_overall[key]) == expected_value,
            f"candidate table global {key} binding mismatch",
        )

    frozen_runtime_global = freeze_runtime(runtime_global)
    final_runtime = final_ranking.get("runtime_diagnostics") or {}
    require(
        frozen_runtime_global["n_unit_evaluations"] == int(final_runtime.get("n_unit_evaluations", -1))
        and frozen_runtime_global["metrics"] == (final_runtime.get("metrics") or {})
        and frozen_runtime_global["metrics_when_invoked"]
        == (final_runtime.get("metrics_when_invoked") or {}),
        "runtime global summary differs from final independent verifier",
    )
    declared_full_runtime = ranking.get("runtime_diagnostics") or {}
    require(
        frozen_runtime_global["n_unit_evaluations"]
        == int(declared_full_runtime.get("n_unit_evaluations", -1))
        and frozen_runtime_global["metrics"] == (declared_full_runtime.get("metrics") or {})
        and frozen_runtime_global["metrics_when_invoked"]
        == (declared_full_runtime.get("metrics_when_invoked") or {}),
        "runtime global summary differs from ranking terminal",
    )
    for dataset in datasets:
        for basis, thresholds in recomputed_rank_dataset[dataset].items():
            for threshold, cell in thresholds.items():
                count = len(
                    runtime_stratum.get((dataset, basis, threshold), new_runtime_accumulator())["all"]
                    [RUNTIME_METRICS[0]]
                )
                require(count == int(cell.get("n_component_hp_units", 0)), f"runtime/unit count mismatch: {dataset}/{basis}/{threshold}")

    extraction_presentation_by_dataset: dict[str, Any] = {}
    for dataset in datasets:
        nested = _nested_component(extraction_dataset_cells[dataset])
        extraction_presentation_by_dataset[dataset] = {
            "counts": dict(extraction_counts_dataset[dataset]),
            "component_by_linkage_basis_threshold": {
                basis: {
                    threshold: presentation_component_cell(cell)
                    for threshold, cell in thresholds.items()
                }
                for basis, thresholds in nested.items()
            },
        }

    ranking_presentation_by_dataset: dict[str, Any] = {}
    for dataset in datasets:
        ranking_presentation_by_dataset[dataset] = {
            "input_call_funnel": dict(input_dataset[dataset]),
            "input_sparse_call_code_counts": dict(sorted(call_dataset[dataset].items())),
            "input_hp_family_rows": dict(sorted(hp_dataset[dataset].items())),
            "by_HP_basis_and_bridge_threshold": {
                basis: {
                    threshold: presentation_rank_cell(
                        cell,
                        candidate_by_stratum.get((dataset, basis, threshold)),
                    )
                    for threshold, cell in thresholds.items()
                }
                for basis, thresholds in recomputed_rank_dataset[dataset].items()
            },
            "runtime_all_primary_unit_evaluations": freeze_runtime(runtime_dataset[dataset]),
            "runtime_by_HP_basis_and_bridge_threshold": {
                basis: {
                    threshold: freeze_runtime(
                        runtime_stratum.get(
                            (dataset, basis, threshold), new_runtime_accumulator()
                        )
                    )
                    for threshold in thresholds
                }
                for basis, thresholds in recomputed_rank_dataset[dataset].items()
            },
        }

    overall_rank_presentation = {
        basis: {
            threshold: presentation_rank_cell(
                cell,
                # Merge already accumulated per-dataset candidate strata for this cell.
                _merge_candidate_strata(candidate_by_stratum, datasets, basis, threshold),
            )
            for threshold, cell in thresholds.items()
        }
        for basis, thresholds in recomputed_rank.items()
    }

    producer_path = Path(__file__).resolve()
    producer_identity = {
        "path": str(producer_path),
        "size_bytes": producer_path.stat().st_size,
        "sha256": sha256_path(producer_path),
    }
    return {
        "schema_name": "intersubmod.m2_final_numeric_summary",
        "schema_version": "1.0.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "scope": {
            "datasets": list(datasets),
            "chromosomes": list(chromosomes),
            "expected_tasks_per_stage": len(expected),
            "n_technical_datasets": len(datasets),
            "n_biological_samples": 6 if datasets == DATASETS else None,
            "scope_is_canonical_full_7_dataset_chr1_22": datasets == DATASETS
            and chromosomes == AUTOSOMES,
        },
        "definitions": {
            "component": "read-linked component within one declared linkage basis and bridge threshold",
            "k": "number of component sites; ranking effective-k may be smaller after structural ALT routing",
            "unit": "HP family x known PS x read-linked component x bridge threshold",
            "T": "sum of feasible parent-edge tree assignments over solver-complete units",
            "V": "sum of distinct globally minimum-extra-state vertex sets over solver-complete units",
            "h_star": (
                "minimum extra-state objective: vertices whose roles contain neither root nor "
                "full-observed; includes partial-compatible representatives/connectors and is not "
                "a hidden-clone count"
            ),
            "partial_read": (
                "a distinct R/A/X pattern contributes one joint subcube group constraint when it "
                "meets structural minread; completion worlds are not enumerated separately"
            ),
            "threshold_warning": (
                "same genomic evidence is re-evaluated at bridge thresholds; summing thresholds "
                "counts evaluations, not deduplicated regions"
            ),
            "HP_warning": (
                "PS_HP1 and PS_HP2 are separate haplotype unit evaluations; their sum is not a "
                "deduplicated cellular clone or biological-region count"
            ),
        },
        "primary_parameters": {
            "bridge_thresholds": list(
                (ranking.get("run_contract") or {}).get("thresholds")
                or []
            ),
            "primary_structural_exact_pattern_minread": int(
                (((ranking.get("run_contract") or {}).get("parameters") or {}).get(
                    "primary_structural_exact_pattern_minread", -1
                ))
            ),
            "HP_bases": list(PRIMARY_BASES),
        },
        "extraction": {
            "overall_counts": dict(extraction_counts_global),
            "overall_component_by_linkage_basis_threshold": {
                basis: {
                    threshold: presentation_component_cell(cell)
                    for threshold, cell in thresholds.items()
                }
                for basis, thresholds in reconstructed_extraction_cells.items()
            },
            "by_dataset": extraction_presentation_by_dataset,
        },
        "ranking": {
            "overall_input_call_funnel": dict(input_global),
            "overall_input_sparse_call_code_counts": dict(sorted(call_global.items())),
            "overall_by_HP_basis_and_bridge_threshold": overall_rank_presentation,
            "overall_candidate_table": freeze_candidate_accumulator(candidate_overall),
            "overall_runtime_all_primary_unit_evaluations": frozen_runtime_global,
            "by_dataset": ranking_presentation_by_dataset,
        },
        "unsupported_or_nonidentifiable": {
            "deduplicated_biological_regions_across_thresholds": null_metric(
                "threshold strata reuse evidence and no cross-threshold region identity table is authenticated"
            ),
            "cellular_HP1_HP2_clone_pairing": null_metric(
                "HP-tagged reads do not establish cell-level pairing between homologous haplotypes"
            ),
            "exact_parent_edge_topology_for_tied_vertex_sets": null_metric(
                "canonical candidate table stores parent_choice_count but not every parent-edge assignment"
            ),
            "global_most_likely_structure_across_h_star_plus_1_or_more": null_metric(
                "likelihood ranking is conditional on the globally minimum h* candidate set"
            ),
            "independent_VAF_confirmation_term": null_metric(
                "same read-pattern evidence already enters the conditional likelihood; adding same-read VAF would double count"
            ),
            "formal_full_run_peak_RSS_CPU_and_disk_envelope": null_metric(
                "scientific receipts authenticate per-unit monotonic runtime but do not contain a formal full-run process/resource telemetry series"
            ),
        },
        "source_ledger": {
            "producer": producer_identity,
            "final_independent_verification": final_identity,
            "terminal_extraction": extraction_identity,
            "terminal_ranking": ranking_identity,
            "candidate_table": candidate_identity,
            "extraction_children": extraction_child_sources,
            "ranking_children": ranking_child_sources,
            "runtime_files": runtime_sources,
        },
        "checks": {
            "canonical_full_scope": (
                datasets == DATASETS and chromosomes == AUTOSOMES
            ) if require_canonical_scope else True,
            "final_verifier_all_checks_true_and_release_eligible": True,
            "terminal_receipts_authenticated_and_bound_to_final_verifier": True,
            "all_extraction_child_receipts_authenticated_and_reaggregated": True,
            "all_ranking_child_receipts_authenticated_and_reaggregated": True,
            "all_component_k_distributions_conserve": True,
            "all_rank_cells_conserve": all(
                all(value is True for value in row["checks"].values())
                for row in conservation_checks
            ),
            "candidate_table_physical_and_semantic_hashes_match": True,
            "candidate_table_T_V_outcome_and_h_star_bindings_conserve": True,
            "tree_vertex_partition_conserves_and_binds_authenticated_T_V": True,
            "tied_coarse_topology_partition_conserves": True,
            "all_runtime_files_authenticated_and_summaries_recomputed": True,
            "unsupported_claims_are_explicit_null_not_zero": True,
        },
        "all_pass": True,
    }


def _merge_candidate_strata(
    strata: Mapping[tuple[str, str, str], Mapping[str, Any]],
    datasets: Sequence[str],
    basis: str,
    threshold: str,
) -> dict[str, Any]:
    merged = new_candidate_accumulator()
    for dataset in datasets:
        source = strata.get((dataset, basis, threshold))
        if source is None:
            continue
        for key in (
            "n_units",
            "n_candidate_vertex_sets_V",
            "n_parent_edge_trees_T",
            "unique_first",
            "tied_first",
            "solver_complete_optimizer_abstain",
            "tied_topology_consistent",
            "tied_topology_inconsistent",
            "topology_evaluated_units",
            "coarse_topology_class_unique_units",
            "coarse_topology_multiple_class_units",
            "parent_edge_assignment_unique_units",
            "exact_topology_proven_unique_units",
        ):
            merged[key] += int(source[key])
        merged["h_star_distribution"].update(source["h_star_distribution"])
        merged["tree_vertex_partition_counts"].update(
            source["tree_vertex_partition_counts"]
        )
        merged["coarse_topology_class_inclusion_counts"].update(
            source["coarse_topology_class_inclusion_counts"]
        )
        merged["coarse_topology_unique_class_counts"].update(
            source["coarse_topology_unique_class_counts"]
        )
        merged["coarse_topology_ambiguous_class_set_counts"].update(
            source["coarse_topology_ambiguous_class_set_counts"]
        )
    return merged


def write_json_and_sidecar_exclusive(path: Path, payload: Mapping[str, Any]) -> dict[str, Any]:
    path = path.resolve()
    sidecar = path.with_name(f"{path.name}.sha256")
    require(not path.exists() and not sidecar.exists(), "output or output sidecar already exists")
    path.parent.mkdir(parents=True, exist_ok=True)
    final_payload = dict(payload)
    final_payload["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": sidecar.name,
    }
    serialized = json.dumps(
        final_payload,
        ensure_ascii=False,
        indent=2,
        allow_nan=False,
    ) + "\n"
    with path.open("x", encoding="utf-8") as handle:
        handle.write(serialized)
        handle.flush()
        os.fsync(handle.fileno())
    digest = sha256_path(path)
    with sidecar.open("x", encoding="ascii") as handle:
        handle.write(f"{digest}  {path.name}\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.chmod(path, 0o444)
    os.chmod(sidecar, 0o444)
    return {"path": str(path), "sha256": digest, "sidecar": str(sidecar)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extraction-root", required=True, type=Path)
    parser.add_argument("--ranking-root", required=True, type=Path)
    parser.add_argument("--final-verification-receipt", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = build_summary(
        args.extraction_root,
        args.ranking_root,
        args.final_verification_receipt,
        require_canonical_scope=True,
    )
    require(summary.get("all_pass") is True, "summary did not reach passing state")
    identity = write_json_and_sidecar_exclusive(args.output, summary)
    print(json.dumps(identity, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
