#!/usr/bin/env python3
"""Build the fail-closed professor-facing Hypercube/likelihood HTML report.

The final report is deliberately gated on all Task-B evidence.  Missing or
invalid full-M2 evidence never becomes a zero and never produces the final
filename.  ``--allow-partial`` is an explicit preview escape hatch; it adds a
prominent PARTIAL ribbon and redirects a final-looking filename to a
``.partial-preview.html`` filename.

Only Python's standard library is used.  Charts are accessible inline SVG and
the resulting HTML has no remote script, font, image, or stylesheet dependency.
"""

from __future__ import annotations

import argparse
import ctypes
import csv
import datetime as dt
import errno
import gzip
import hashlib
import html
import json
import math
import os
import re
import stat
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from urllib.parse import quote
from zoneinfo import ZoneInfo


DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
BIOLOGICAL_IDS = {
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}
FINAL_REPORT_NAME = "20260716_read_linked_hypercube_exact_likelihood全量驗證報告_01.html"
DISPLAY_BRIDGE_THRESHOLD = "3"

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
TOPOLOGY_CLASS_ORDER = {
    "single-only": 0,
    "sister-only": 1,
    "direct-only": 2,
    "sister+direct": 3,
}

RUNTIME_METRICS = (
    "candidate_generation_elapsed_seconds",
    "likelihood_fit_elapsed_seconds",
    "unit_total_elapsed_seconds",
)

EXPECTED_PARTIAL_FUNNEL_DEFINITIONS = {
    "conceptual": (
        "a pattern with u X has 2^u full-cube conceptual completions; completion-wise tree worlds "
        "are not materialized, while active compatible vertex indices are materialized per MILP "
        "construction for each reduced sparse group-hit row"
    ),
    "effective": (
        "primary observed-alt universe: 2^popcount(free_mask & structural_alt_universe); "
        "zero when a scoring pattern fixes ALT outside that minread-specific universe"
    ),
    "universe_source": "exact R/A/X structural pattern count >= 3",
}

EXPECTED_METHOD_CONTRACT = {
    "schema_name": "intersubmod.m2_partial_read_likelihood_method_contract",
    "schema_version": "1.0.0",
    "structural_group_scope": "DISTINCT_EXACT_RAX_COUNT_GE_THRESHOLD",
    "active_compatible_vertex_indices_materialized_for_sparse_rows": True,
    "completion_wise_tree_worlds_materialized": False,
    "cross_read_cartesian_products_materialized": False,
    "likelihood_primary_term": (
        "BQ_SYMMETRIC_SUBSTITUTION_CONDITIONAL_RA_READ_PATTERN_MIXTURE_PROFILE_LIKELIHOOD"
    ),
    "same_molecule_vaf_added_as_separate_term": False,
    "parent_edges_or_transitions_scored": False,
    "missing_calls_marginalized": True,
}

REQUIRED_RANKING_CHECKS = (
    "all_tasks_pass",
    "all_datasets_present",
    "all_autosomes_present",
    "all_child_receipts_schema_2_0",
    "all_upstream_output_hashes_verified",
    "parameters_match_extraction",
    "inputs_match_extraction_receipt",
    "aggregate_conserved",
    "all_154_child_method_contracts_identical_and_source_bound",
    "k_gt12_never_claimed_global_optimal",
    "same_read_vaf_not_double_counted",
    "no_cross_ps_pattern_pooling",
    "known_ps_never_mixed",
    "missing_ps_separate_diagnostic",
    "conditional_ra_model_not_claimed_full_generative",
    "bootstrap_is_fixed_candidate_set_ranking_stability",
    "topology_inclusion_counts_not_used_as_composition",
    "candidate_table_hash_verified",
    "candidate_table_row_schema_complete",
    "candidate_table_covers_all_rankable_units",
    "all_child_runtime_diagnostics_contracts_valid",
    "full_runtime_aggregate_covers_all_154_children",
    "full_runtime_aggregate_covers_all_primary_unit_evaluations",
    "full_runtime_aggregate_recomputed_from_streamed_unit_rows",
)

CANDIDATE_TABLE_COLUMNS = (
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
REQUIRED_RANKING_CELL_KEYS = (
    "n_components",
    "n_component_hp_units",
    "molecule_component_projections",
    "informative_scoring_molecules",
    "all_x_excluded_molecules",
    "structural_retained_molecules",
    "below_minread_scoring_molecules",
    "bq_scoring_molecules",
    "projected_call_class_counts",
    "selection_status_counts",
    "candidate_generation_status_counts",
    "solver_complete_incomplete_counts",
    "k_route_counts",
    "solver_complete_units",
    "solver_incomplete_or_not_run_units",
    "quality_primary_unique_vertex_units",
    "quality_primary_tied_vertex_units",
    "rank_abstain_units",
    "raw_tree_candidates_T_complete_units",
    "distinct_vertex_sets_V_complete_units",
    "likelihood_stability_counts",
    "topology_evaluated_units",
    "coarse_topology_unique_units",
    "coarse_topology_unique_class_counts",
    "coarse_topology_ambiguous_units",
    "ambiguous_topology_class_set_counts",
    "topology_derivation_status_counts",
    "exact_topology_uniqueness_status_counts",
    "parent_edge_assignment_unique_units",
    "exact_topology_proven_unique_units",
    "topology_class_inclusion_counts_denominator",
    "k_component_sites_total",
    "k_observed_alt_active_total",
    "k_scoring_alt_observed_total",
    "not_structural_alt_active_sites_total",
    "structural_partial_pattern_groups",
    "partial_group_coverage_denominator",
    "partial_groups_covered",
    "partial_groups_unsatisfied",
    "partial_pattern_denominators",
    "partial_u_distribution",
    "denominator_map",
    "conservation_checks",
    "all_conserved",
    "partial_pattern_funnel",
)

PARTIAL_FUNNEL_GRAINS = (
    "unique_rax_pattern_groups",
    "bq_quality_pattern_groups",
    "molecule_projections",
    "structural_unique_rax_pattern_groups",
)

EXPECTED_DENOMINATOR_MAP_KEYS = {
    "unit_denominator",
    "molecule_projection_denominator",
    "informative_molecule_denominator",
    "solver_complete_unit_denominator",
    "topology_evaluated_unit_denominator",
    "partial_group_coverage_denominator",
    "topology_inclusion_denominator",
}

EXPECTED_CONSERVATION_CHECK_KEYS = {
    "selection_status_sum_equals_units",
    "unique_tied_abstain_partition_units",
    "solver_complete_plus_notrun_equals_units",
    "projection_funnel_conserved",
    "structural_scoring_funnel_conserved",
    "bq_molecules_equal_informative",
    "raw_T_not_less_than_distinct_V",
    "topology_unique_plus_multiple_equals_evaluated",
    "partial_coverage_conserved_and_zero_unsatisfied",
    "k_route_partition_equals_units",
}

NUMERIC_SUMMARY_REQUIRED_CHECKS = {
    "canonical_full_scope",
    "final_verifier_all_checks_true_and_release_eligible",
    "terminal_receipts_authenticated_and_bound_to_final_verifier",
    "all_extraction_child_receipts_authenticated_and_reaggregated",
    "all_ranking_child_receipts_authenticated_and_reaggregated",
    "all_component_k_distributions_conserve",
    "all_rank_cells_conserve",
    "candidate_table_physical_and_semantic_hashes_match",
    "candidate_table_T_V_outcome_and_h_star_bindings_conserve",
    "tied_coarse_topology_partition_conserves",
    "all_runtime_files_authenticated_and_summaries_recomputed",
    "unsupported_claims_are_explicit_null_not_zero",
}

NUMERIC_SUMMARY_UNSUPPORTED_KEYS = {
    "deduplicated_biological_regions_across_thresholds",
    "cellular_HP1_HP2_clone_pairing",
    "exact_parent_edge_topology_for_tied_vertex_sets",
    "global_most_likely_structure_across_h_star_plus_1_or_more",
    "independent_VAF_confirmation_term",
    "formal_full_run_peak_RSS_CPU_and_disk_envelope",
}

NUMERIC_SUMMARY_RATIO_KEYS = {
    "numerator",
    "denominator",
    "value",
    "percent",
    "denominator_label",
    "reason",
}


class ReportGateError(RuntimeError):
    """Raised when evidence cannot support the requested report status."""


@dataclass(frozen=True)
class Source:
    source_id: str
    label: str
    path: Path
    sha256: str
    scope: str
    payload: Any | None = None


@dataclass
class Assessment:
    sources: dict[str, Source] = field(default_factory=dict)
    hard_issues: list[str] = field(default_factory=list)
    final_issues: list[str] = field(default_factory=list)

    @property
    def final_ready(self) -> bool:
        return not self.hard_issues and not self.final_issues


@dataclass(frozen=True)
class BuildResult:
    output_path: Path
    final_ready: bool
    hard_issues: tuple[str, ...]
    final_issues: tuple[str, ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_nonnegative_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _is_finite_nonnegative_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and float(value) >= 0.0
    )


def _require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def _source(
    source_id: str,
    label: str,
    path: Path,
    scope: str,
    *,
    json_payload: bool = True,
) -> Source:
    if not path.is_file():
        raise FileNotFoundError(path)
    return Source(
        source_id=source_id,
        label=label,
        path=path.resolve(),
        sha256=_sha256(path),
        scope=scope,
        payload=_load_json(path) if json_payload else None,
    )


def _sidecar_authenticated_source(
    source_id: str, label: str, path: Path, scope: str
) -> tuple[Source | None, list[str]]:
    """Read a JSON source only when its adjacent SHA-256 sidecar is exact."""

    errors: list[str] = []
    try:
        source = _source(source_id, label, path, scope)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return None, [f"{label} 無法讀取：{exc}"]
    sidecar = source.path.with_name(f"{source.path.name}.sha256")
    if not sidecar.is_file():
        return None, [f"{label} 缺 SHA-256 sidecar：{sidecar}"]
    try:
        fields = sidecar.read_text(encoding="ascii", errors="strict").strip().split()
    except (OSError, UnicodeError) as exc:
        return None, [f"{label} SHA-256 sidecar 無法讀取：{exc}"]
    _require(
        fields == [source.sha256, source.path.name],
        f"{label} SHA-256 sidecar 未精確綁定檔名與內容",
        errors,
    )
    return source, errors


def _new_candidate_aggregate() -> dict[str, Any]:
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


def _parse_candidate_topology(raw: str, label: str) -> tuple[str, ...]:
    try:
        values = json.loads(raw)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ReportGateError(f"M2 candidate topology JSON 錯誤：{label}") from exc
    if (
        not isinstance(values, list)
        or not values
        or any(not isinstance(value, str) or not value for value in values)
        or not set(values).issubset(TOPOLOGY_CLASS_ORDER)
    ):
        raise ReportGateError(f"M2 candidate topology class 無效：{label}")
    return tuple(sorted(set(values), key=TOPOLOGY_CLASS_ORDER.__getitem__))


def _freeze_candidate_aggregate(target: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: int(target[key])
        for key in (
            "n_units", "n_candidate_vertex_sets_V", "n_parent_edge_trees_T",
            "unique_first", "tied_first", "solver_complete_optimizer_abstain",
            "tied_topology_consistent", "tied_topology_inconsistent",
            "topology_evaluated_units", "coarse_topology_class_unique_units",
            "coarse_topology_multiple_class_units",
            "parent_edge_assignment_unique_units",
            "exact_topology_proven_unique_units",
        )
    } | {
        "tree_vertex_partition_counts": {
            key: int(target["tree_vertex_partition_counts"].get(key, 0))
            for key in TREE_VERTEX_BUCKETS
        },
        "h_star_distribution": {
            str(key): int(value)
            for key, value in sorted(target["h_star_distribution"].items())
        },
        "coarse_topology_class_inclusion_counts": dict(
            sorted(target["coarse_topology_class_inclusion_counts"].items())
        ),
        "coarse_topology_unique_class_counts": dict(
            sorted(target["coarse_topology_unique_class_counts"].items())
        ),
        "coarse_topology_ambiguous_class_set_counts": dict(
            sorted(target["coarse_topology_ambiguous_class_set_counts"].items())
        ),
    }


def _candidate_table_source(ranking_source: Source) -> tuple[Source | None, list[str]]:
    """Verify and independently aggregate the canonical candidate table.

    This stream is deliberately independent of S16.  In addition to bounded
    display examples it reconstructs every dataset×HP×threshold candidate
    count, T/V partition, winner state, h*, and topology uniqueness cell.
    """
    errors: list[str] = []
    metadata = (ranking_source.payload or {}).get("candidate_table") or {}
    _require(metadata.get("schema_version") == "2.0.0", "M2 candidate table schema_version 不是 2.0.0", errors)
    _require(metadata.get("format") == "tsv.gz", "M2 candidate table format 不是 tsv.gz", errors)
    _require(tuple(metadata.get("columns") or ()) == CANDIDATE_TABLE_COLUMNS, "M2 candidate table columns/order 不符", errors)
    _require(metadata.get("sort_order") == "unit_key,candidate_id", "M2 candidate table sort_order 不符", errors)
    _require(_is_nonnegative_int(metadata.get("n_rows")) and metadata.get("n_rows", 0) > 0, "M2 candidate table n_rows 缺失", errors)
    _require(_is_nonnegative_int(metadata.get("n_units")) and metadata.get("n_units", 0) > 0, "M2 candidate table n_units 缺失", errors)
    _require(isinstance(metadata.get("sha256"), str) and len(metadata.get("sha256", "")) == 64, "M2 candidate table sha256 缺失", errors)
    _require(isinstance(metadata.get("semantic_sha256"), str) and len(metadata.get("semantic_sha256", "")) == 64, "M2 candidate table semantic_sha256 缺失", errors)
    raw_path = metadata.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        errors.append("M2 candidate table path 缺失")
        return None, errors
    path = Path(raw_path)
    if not path.is_absolute():
        path = ranking_source.path.parent / path
    path = path.resolve()
    if not path.is_file():
        errors.append(f"M2 candidate table 不存在：{path}")
        return None, errors
    _require(path.stat().st_size == metadata.get("size_bytes"), "M2 candidate table size_bytes 不符", errors)
    actual_sha = _sha256(path)
    _require(actual_sha == metadata.get("sha256"), "M2 candidate table SHA-256 不符", errors)
    if errors:
        return None, errors

    row_count = 0
    unit_count = 0
    previous_unit: str | None = None
    current_rows: list[dict[str, str]] = []
    example_groups: dict[str, list[dict[str, str]]] = {}
    minimum_extra_state_distribution: Counter[int] = Counter()
    by_stratum: dict[tuple[str, str, str], dict[str, Any]] = {}
    overall = _new_candidate_aggregate()

    def finish_unit(rows: list[dict[str, str]]) -> None:
        if not rows:
            return
        first = rows[0]
        unit_key = first["unit_key"]
        dataset = first["dataset"]
        basis = f'PS_HP{first["hp_family"]}'
        try:
            threshold = str(int(first["threshold"]))
        except (TypeError, ValueError) as exc:
            errors.append(f"M2 candidate threshold 無效：{unit_key}")
            return
        if dataset not in DATASETS or basis not in {"PS_HP1", "PS_HP2"} or threshold not in {"1", "2", "3", "5"}:
            errors.append(f"M2 candidate stratum 無效：{unit_key}")
            return
        if any(
            row["unit_key"] != unit_key
            or row["dataset"] != dataset
            or row["chrom"] != first["chrom"]
            or row["component_id"] != first["component_id"]
            or row["hp_family"] != first["hp_family"]
            or row["ps"] != first["ps"]
            or str(int(row["threshold"])) != threshold
            for row in rows
        ):
            errors.append(f"M2 candidate 同一 unit scope 漂移：{unit_key}")
            return
        if len({row["candidate_id"] for row in rows}) != len(rows) or len({row["vertex_set_id"] for row in rows}) != len(rows):
            errors.append(f"M2 candidate unit 有重複 candidate/vertex set：{unit_key}")
            return
        unit_extra_counts: set[int] = set()
        parent_choice_counts: list[int] = []
        for row in rows:
            try:
                states = json.loads(row["vertex_states"])
                roles = json.loads(row["vertex_roles"])
            except (TypeError, json.JSONDecodeError):
                errors.append(f"M2 candidate row vertex/role JSON 錯誤：{row['unit_key']}")
                return
            if (
                not isinstance(states, dict)
                or not isinstance(roles, dict)
                or set(states) != set(roles)
                or not roles
                or any(
                    not isinstance(role_values, list)
                    or not role_values
                    or any(not isinstance(role, str) or not role for role in role_values)
                    for role_values in roles.values()
                )
            ):
                errors.append(f"M2 candidate row vertex/role production schema 錯誤：{row['unit_key']}")
                return
            # This exactly mirrors the solver objective: root and fully observed
            # mandatory states cost zero; a partial representative or connector
            # is objective-bearing.  It is an extra-state count, not a clone count.
            unit_extra_counts.add(
                sum(
                    "root" not in role_values and "full-observed" not in role_values
                    for role_values in roles.values()
                )
            )
            if row["candidate_set_complete"].lower() != "true":
                errors.append(f"M2 canonical candidate table 含 incomplete candidate：{row['unit_key']}")
                return
            try:
                parent_choice_counts.append(int(row["parent_choice_count"]))
            except (TypeError, ValueError):
                errors.append(f"M2 candidate parent_choice_count 無效：{unit_key}")
                return
        if len(unit_extra_counts) != 1:
            errors.append(f"M2 同一 unit 的 minimum-extra-state objective 不一致：{unit_key}")
            return
        if any(value < 1 for value in parent_choice_counts):
            errors.append(f"M2 candidate parent_choice_count 必須 >=1：{unit_key}")
            return
        h_star = next(iter(unit_extra_counts))
        minimum_extra_state_distribution[h_star] += 1
        statuses = [row["winner_status"].upper() for row in rows]
        if not set(statuses).issubset({"UNIQUE_WINNER", "TIED_WINNER", "NON_WINNER", "ABSTAIN_UNIT_OPTIMIZER"}):
            errors.append(f"M2 candidate winner status 無效：{unit_key}")
            return
        unique_rows = [row for row in rows if row["winner_status"].upper() == "UNIQUE_WINNER"]
        tied_rows = [row for row in rows if row["winner_status"].upper() == "TIED_WINNER"]
        abstain_rows = [row for row in rows if row["winner_status"].upper() == "ABSTAIN_UNIT_OPTIMIZER"]
        if abstain_rows:
            if len(abstain_rows) != len(rows):
                errors.append(f"M2 candidate abstain partition 混合其他 status：{unit_key}")
                return
            category = "ABSTAIN"
        elif unique_rows:
            if len(unique_rows) != 1 or tied_rows:
                errors.append(f"M2 candidate unique winner partition 無效：{unit_key}")
                return
            category = "UNIQUE"
        else:
            if len(tied_rows) < 2:
                errors.append(f"M2 candidate tied winner partition 無效：{unit_key}")
                return
            category = "TIE"

        n_v = len(rows)
        n_t = sum(parent_choice_counts)
        if not n_t >= n_v >= 1:
            errors.append(f"M2 candidate T>=V>=1 不成立：{unit_key}")
            return
        if n_t == 1 and n_v == 1:
            tree_vertex_bucket = "T_EQ_1_V_EQ_1"
        elif n_t > 1 and n_v == 1:
            tree_vertex_bucket = "T_GT_1_V_EQ_1"
        elif n_t > 1 and n_v > 1:
            tree_vertex_bucket = "T_GT_1_V_GT_1"
        else:
            errors.append(f"M2 candidate T/V 落入未定義 partition：{unit_key}")
            return

        targets = (
            by_stratum.setdefault((dataset, basis, threshold), _new_candidate_aggregate()),
            overall,
        )
        for target in targets:
            target["n_units"] += 1
            target["n_candidate_vertex_sets_V"] += n_v
            target["n_parent_edge_trees_T"] += n_t
            target["tree_vertex_partition_counts"][tree_vertex_bucket] += 1
            target["h_star_distribution"][h_star] += 1
            if category == "UNIQUE":
                target["unique_first"] += 1
            elif category == "TIE":
                target["tied_first"] += 1
            else:
                target["solver_complete_optimizer_abstain"] += 1

        winning_rows = unique_rows or tied_rows
        if winning_rows:
            union: set[str] = set()
            try:
                for row in winning_rows:
                    union.update(_parse_candidate_topology(row["coarse_topology_class"], unit_key))
            except ReportGateError as exc:
                errors.append(str(exc))
                return
            ordered_union = tuple(sorted(union, key=TOPOLOGY_CLASS_ORDER.__getitem__))
            winning_parent_choices = sum(int(row["parent_choice_count"]) for row in winning_rows)
            for target in targets:
                target["topology_evaluated_units"] += 1
                target["coarse_topology_class_inclusion_counts"].update(ordered_union)
                if len(ordered_union) == 1:
                    target["coarse_topology_class_unique_units"] += 1
                    target["coarse_topology_unique_class_counts"][ordered_union[0]] += 1
                else:
                    target["coarse_topology_multiple_class_units"] += 1
                    target["coarse_topology_ambiguous_class_set_counts"]["|".join(ordered_union)] += 1
                if winning_parent_choices == 1:
                    target["parent_edge_assignment_unique_units"] += 1
                    target["exact_topology_proven_unique_units"] += 1
            if category == "TIE":
                tie_key = "tied_topology_consistent" if len(ordered_union) == 1 else "tied_topology_inconsistent"
                for target in targets:
                    target[tie_key] += 1
        if category not in example_groups:
            example_groups[category] = [dict(row) for row in rows[:12]]

    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != CANDIDATE_TABLE_COLUMNS:
                errors.append("M2 candidate table 實際 header/order 不符")
                return None, errors
            for row in reader:
                row_count += 1
                unit_key = row["unit_key"]
                if previous_unit is None or unit_key != previous_unit:
                    if previous_unit is not None:
                        finish_unit(current_rows)
                        current_rows = []
                        if unit_key < previous_unit:
                            errors.append("M2 candidate table 未依 unit_key 排序")
                    unit_count += 1
                    previous_unit = unit_key
                current_rows.append(row)
                if row["dataset"] not in DATASETS or row["chrom"] not in AUTOSOMES:
                    errors.append(f"M2 candidate row scope 不符：{unit_key}")
                    break
                if row["ps"] in {"", ".", "NA", "None", "unknown", "UNKNOWN"}:
                    errors.append(f"M2 primary candidate table 含 missing PS：{unit_key}")
                    break
                try:
                    parent_choices = int(row["parent_choice_count"])
                    profile_ll = float(row["profile_log_likelihood"])
                    relative_ll = float(row["relative_log_likelihood"])
                    weights = json.loads(row["mixture_weights_pi"])
                except (ValueError, TypeError, json.JSONDecodeError):
                    errors.append(f"M2 candidate row 數值/π 格式錯誤：{unit_key}")
                    break
                weight_values = list(weights.values()) if isinstance(weights, dict) else weights if isinstance(weights, list) else []
                if (
                    parent_choices < 1
                    or not math.isfinite(profile_ll)
                    or not math.isfinite(relative_ll)
                    or not weight_values
                    or any(not isinstance(value, (int, float)) or value < 0 or value > 1 for value in weight_values)
                    or not math.isclose(sum(float(value) for value in weight_values), 1.0, abs_tol=1e-5)
                ):
                    errors.append(f"M2 candidate row parent/LL/π constraint 錯誤：{unit_key}")
                    break
                if not row["vertex_states"] or not row["vertex_roles"] or not row["coarse_topology_class"]:
                    errors.append(f"M2 candidate row vertex/topology 欄位缺失：{unit_key}")
                    break
            finish_unit(current_rows)
    except (OSError, UnicodeError, csv.Error) as exc:
        errors.append(f"M2 candidate table 無法讀取：{exc}")
    _require(row_count == metadata.get("n_rows"), "M2 candidate table 實際 n_rows 不符", errors)
    _require(unit_count == metadata.get("n_units"), "M2 candidate table 實際 n_units 不符", errors)
    _require(
        sum(minimum_extra_state_distribution.values()) == metadata.get("n_units"),
        "M2 candidate table minimum-extra-state h* 分布未覆蓋所有 units",
        errors,
    )
    _require(bool(example_groups), "M2 candidate table 無可展示的 unit example", errors)
    if errors:
        return None, errors
    return (
        Source(
            source_id="S9",
            label="M2 PS-aware candidate table",
            path=path,
            sha256=actual_sha,
            scope=(
                "unit key、minimum-extra-state h*、vertex/role、parent-choice count、"
                "LL/relative score、π、winner/tie、coarse class"
            ),
            payload={
                "metadata": metadata,
                "example_groups": example_groups,
                "minimum_extra_state_distribution": {
                    str(key): value
                    for key, value in sorted(minimum_extra_state_distribution.items())
                },
                "by_stratum": {
                    key: _freeze_candidate_aggregate(value)
                    for key, value in sorted(by_stratum.items())
                },
                "overall": _freeze_candidate_aggregate(overall),
            },
        ),
        [],
    )


def _validate_canonical(payload: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.current_layered_topology_summary", "canonical schema_name 不符", errors)
    _require(payload.get("all_pass") is True, "canonical all_pass 不是 true", errors)
    _require(payload.get("task_type") == "B_comprehensive_validation", "canonical 不是 Task-B comprehensive validation", errors)
    canonical = payload.get("canonical") or {}
    aggregate = canonical.get("aggregate") or {}
    samples = canonical.get("samples") or []
    _require(aggregate.get("dataset_count") == 7, "canonical dataset_count 不是 7", errors)
    _require(isinstance(samples, list) and len(samples) == 7, "canonical samples 不是 7 datasets", errors)
    if isinstance(samples, list):
        _require({row.get("sample") for row in samples} == set(DATASETS), "canonical dataset 名單不完整", errors)
    required_counts = (
        "tree_input_records",
        "autosomal_biallelic_sSNV",
        "retained_sSNV",
        "W_tree",
        "W_primary",
        "no_primary_lineage",
        "primary_units",
        "complete_regions",
        "incomplete_regions",
    )
    for key in required_counts:
        _require(_is_nonnegative_int(aggregate.get(key)), f"canonical aggregate.{key} 缺失或型別錯誤", errors)
    if all(_is_nonnegative_int(aggregate.get(key)) for key in ("W_primary", "complete_regions", "incomplete_regions")):
        _require(
            aggregate["complete_regions"] + aggregate["incomplete_regions"] == aggregate["W_primary"],
            "canonical complete + incomplete 不等於 W_primary",
            errors,
        )
    if all(_is_nonnegative_int(aggregate.get(key)) for key in ("W_tree", "W_primary", "no_primary_lineage")):
        _require(
            aggregate["W_primary"] + aggregate["no_primary_lineage"] == aggregate["W_tree"],
            "canonical W_primary + no_primary_lineage 不等於 W_tree",
            errors,
        )
    topology = aggregate.get("topology_classes") or {}
    if _is_nonnegative_int(aggregate.get("W_primary")):
        topo_keys = (
            "exact_and_topology_unique",
            "topology_unique_exact_multiple",
            "topology_multiple_exact_multiple",
            "incomplete",
        )
        _require(all(_is_nonnegative_int(topology.get(key)) for key in topo_keys), "canonical topology_classes 不完整", errors)
        if all(_is_nonnegative_int(topology.get(key)) for key in topo_keys):
            _require(sum(topology[key] for key in topo_keys) == aggregate["W_primary"], "canonical topology class 總和不等於 W_primary", errors)
    if isinstance(samples, list) and len(samples) == 7:
        for key in required_counts:
            if all(_is_nonnegative_int(row.get(key)) for row in samples) and _is_nonnegative_int(aggregate.get(key)):
                _require(sum(row[key] for row in samples) == aggregate[key], f"canonical per-dataset {key} 不守恆", errors)
        dataset_topology_keys = (
            "exact_and_topology_unique",
            "topology_unique_exact_multiple",
            "topology_multiple_exact_multiple",
            "incomplete",
            "impossible_exact_unique_topology_multiple",
        )
        for row in samples:
            dataset = row.get("sample", "UNKNOWN")
            dataset_counts = (
                "W_tree",
                "W_primary",
                "no_primary_lineage",
                "complete_regions",
                "incomplete_regions",
            )
            counts_valid = all(_is_nonnegative_int(row.get(key)) for key in dataset_counts)
            _require(counts_valid, f"canonical {dataset} internal region counts 不完整", errors)
            if counts_valid:
                _require(
                    row["complete_regions"] + row["incomplete_regions"] == row["W_primary"],
                    f"canonical {dataset} complete + incomplete 不等於 W_primary",
                    errors,
                )
                _require(
                    row["W_primary"] + row["no_primary_lineage"] == row["W_tree"],
                    f"canonical {dataset} W_primary + no_primary_lineage 不等於 W_tree",
                    errors,
                )
            dataset_topology = row.get("topology_classes") or {}
            topology_valid = all(_is_nonnegative_int(dataset_topology.get(key)) for key in dataset_topology_keys)
            _require(topology_valid, f"canonical {dataset} topology_classes 不完整", errors)
            if topology_valid and _is_nonnegative_int(row.get("W_primary")):
                _require(
                    sum(dataset_topology[key] for key in dataset_topology_keys) == row["W_primary"],
                    f"canonical {dataset} topology 互斥分類總和不等於 W_primary",
                    errors,
                )
        _require(all(row.get("all_invariants_pass") is True for row in samples), "canonical 有 dataset invariant 未通過", errors)
    return errors


def _validate_funnel(
    payload: Mapping[str, Any],
    canonical_payload: Mapping[str, Any],
    canonical_sha256: str,
) -> list[str]:
    """Validate the source-backed raw-ClairS-to-retained funnel receipt."""
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod_current_ssnv_funnel_receipt", "funnel schema_name 不符", errors)
    _require(payload.get("schema_version") == "1.0.0", "funnel schema_version 不是 1.0.0", errors)
    _require(payload.get("task_type") == "B_comprehensive_validation", "funnel 不是 Task-B", errors)
    _require(payload.get("all_pass") is True, "funnel all_pass 不是 true", errors)
    scope = payload.get("scope") or {}
    _require(tuple(scope.get("datasets") or ()) == DATASETS, "funnel datasets/order 不符", errors)
    _require(scope.get("dataset_count") == 7, "funnel dataset_count 不是 7", errors)
    _require(scope.get("biological_sample_count") == 6, "funnel biological_sample_count 不是 6", errors)
    checks = payload.get("checks") or {}
    _require(bool(checks) and all(value is True for value in checks.values()), "funnel checks 未全數 true", errors)

    inputs = payload.get("inputs") or {}
    canonical_input = inputs.get("canonical_json") or {}
    _require(canonical_input.get("sha256") == canonical_sha256, "funnel canonical input SHA-256 不符", errors)
    ledger_sources = inputs.get("site_ledger_summaries") or []
    _require(isinstance(ledger_sources, list) and len(ledger_sources) == 7, "funnel site ledger sources 不是 7 份", errors)
    if isinstance(ledger_sources, list):
        _require(tuple(row.get("dataset") for row in ledger_sources) == DATASETS, "funnel site ledger dataset/order 不符", errors)
        for row in ledger_sources:
            path_value = row.get("path")
            expected_sha = row.get("sha256")
            if not isinstance(path_value, str) or not path_value:
                errors.append(f"funnel {row.get('dataset')} ledger path 缺失")
                continue
            path = Path(path_value)
            if not path.is_file():
                errors.append(f"funnel source ledger 不存在：{path}")
                continue
            _require(_sha256(path) == expected_sha, f"funnel source ledger SHA-256 不符：{path}", errors)

    aggregate = payload.get("aggregate") or {}
    count_keys = (
        "raw_clairs_records",
        "longphase_s_recalibrated_pass",
        "autosomal_biallelic_sSNV",
        "retained_sSNV",
    )
    _require(all(_is_nonnegative_int(aggregate.get(key)) for key in count_keys), "funnel aggregate counts 不完整", errors)
    if all(_is_nonnegative_int(aggregate.get(key)) for key in count_keys):
        raw, tree, autosomal, retained = (int(aggregate[key]) for key in count_keys)
        _require(raw >= tree >= autosomal >= retained, "funnel aggregate 不是 raw≥PASS≥autosomal≥retained", errors)
        canonical_aggregate = (canonical_payload.get("canonical") or {}).get("aggregate") or {}
        _require(tree == canonical_aggregate.get("tree_input_records"), "funnel PASS 與 canonical tree_input 不符", errors)
        _require(autosomal == canonical_aggregate.get("autosomal_biallelic_sSNV"), "funnel autosomal 與 canonical 不符", errors)
        _require(retained == canonical_aggregate.get("retained_sSNV"), "funnel retained 與 canonical 不符", errors)
        branches = aggregate.get("branch_counts") or {}
        branch_keys = (
            "excluded_by_longphase_filter",
            "out_of_scope_non_autosomal",
            "max_snv_excluded",
            "positional_singleton",
            "retained",
        )
        _require(all(_is_nonnegative_int(branches.get(key)) for key in branch_keys), "funnel aggregate branch_counts 不完整", errors)
        if all(_is_nonnegative_int(branches.get(key)) for key in branch_keys):
            _require(sum(branches[key] for key in branch_keys) == raw, "funnel aggregate raw branches 不守恆", errors)
            _require(tree == raw - branches["excluded_by_longphase_filter"], "funnel raw→PASS branch 不守恆", errors)
            _require(autosomal == tree - branches["out_of_scope_non_autosomal"], "funnel PASS→autosomal branch 不守恆", errors)
            _require(
                autosomal == branches["max_snv_excluded"] + branches["positional_singleton"] + branches["retained"],
                "funnel autosomal retained/singleton/MAX_SNV branches 不守恆",
                errors,
            )
        ratios = aggregate.get("relative_ratios") or {}
        totals = aggregate.get("total_ratios_over_raw") or {}
        expected_ratios = {
            "longphase_pass_over_raw": tree / raw if raw else None,
            "autosomal_over_longphase_pass": autosomal / tree if tree else None,
            "retained_over_autosomal": retained / autosomal if autosomal else None,
        }
        expected_totals = {
            "autosomal": autosomal / raw if raw else None,
            "retained": retained / raw if raw else None,
        }
        _require(
            all(isinstance(ratios.get(key), (int, float)) and math.isclose(ratios[key], value, abs_tol=1e-12) for key, value in expected_ratios.items()),
            "funnel relative ratios 不是由 counts 重算",
            errors,
        )
        _require(
            all(isinstance(totals.get(key), (int, float)) and math.isclose(totals[key], value, abs_tol=1e-12) for key, value in expected_totals.items()),
            "funnel total ratios 不是由 counts 重算",
            errors,
        )

    dataset_rows = payload.get("datasets") or []
    canonical_rows = {
        row.get("sample"): row for row in ((canonical_payload.get("canonical") or {}).get("samples") or [])
    }
    _require(isinstance(dataset_rows, list) and len(dataset_rows) == 7, "funnel dataset rows 不是 7", errors)
    if isinstance(dataset_rows, list):
        _require(tuple(row.get("dataset") for row in dataset_rows) == DATASETS, "funnel dataset rows/order 不符", errors)
        for row in dataset_rows:
            dataset = row.get("dataset")
            values = [row.get(key) for key in count_keys]
            _require(all(_is_nonnegative_int(value) for value in values), f"funnel {dataset} counts 不完整", errors)
            if all(_is_nonnegative_int(value) for value in values):
                raw, tree, autosomal, retained = (int(value) for value in values)
                _require(raw >= tree >= autosomal >= retained, f"funnel {dataset} 順序不守恆", errors)
                canonical_row = canonical_rows.get(dataset) or {}
                _require(tree == canonical_row.get("tree_input_records"), f"funnel {dataset} PASS 與 canonical 不符", errors)
                _require(autosomal == canonical_row.get("autosomal_biallelic_sSNV"), f"funnel {dataset} autosomal 與 canonical 不符", errors)
                _require(retained == canonical_row.get("retained_sSNV"), f"funnel {dataset} retained 與 canonical 不符", errors)
        if all(isinstance(row, Mapping) for row in dataset_rows):
            for key in count_keys:
                if all(_is_nonnegative_int(row.get(key)) for row in dataset_rows) and _is_nonnegative_int(aggregate.get(key)):
                    _require(sum(row[key] for row in dataset_rows) == aggregate[key], f"funnel per-dataset {key} 加總不符", errors)
    return errors


def _validate_funnel_verification(
    payload: Mapping[str, Any],
    verification_source: Source,
    canonical_source: Source,
    funnel_source: Source,
) -> list[str]:
    """Validate the independent 322-check current-funnel recomputation receipt."""

    errors: list[str] = []
    _require(
        payload.get("schema_name") == "intersubmod_current_ssnv_funnel_independent_verification",
        "funnel verifier schema_name 不符",
        errors,
    )
    _require(payload.get("schema_version") == "1.0.0", "funnel verifier schema_version 不是 1.0.0", errors)
    _require(payload.get("task_type") == "B_comprehensive_validation", "funnel verifier 不是 Task-B", errors)
    _require(payload.get("all_pass") is True, "funnel verifier all_pass 不是 true", errors)

    scope = payload.get("scope") or {}
    _require(tuple(scope.get("datasets") or ()) == DATASETS, "funnel verifier datasets/order 不符", errors)
    _require(scope.get("dataset_count") == 7, "funnel verifier dataset_count 不是 7", errors)
    _require(scope.get("biological_sample_count") == 6, "funnel verifier biological_sample_count 不是 6", errors)
    _require(
        scope.get("chromosomes") == "chr1-chr22 for autosomal_biallelic_sSNV and downstream",
        "funnel verifier chromosome scope 不符",
        errors,
    )

    checks = payload.get("checks") or {}
    _require(isinstance(checks, dict) and len(checks) == 322, "funnel verifier checks 不是 322 項", errors)
    if isinstance(checks, dict):
        _require(all(value is True for value in checks.values()), "funnel verifier checks 未全數 true", errors)
    _require(payload.get("n_checks") == 322, "funnel verifier n_checks 不是 322", errors)
    _require(payload.get("n_pass") == 322, "funnel verifier n_pass 不是 322", errors)
    _require(payload.get("n_fail") == 0, "funnel verifier n_fail 不是 0", errors)
    _require(payload.get("failures") == [], "funnel verifier failures 不是空清單", errors)

    inputs = payload.get("inputs") or {}
    bindings = (
        ("canonical_json", canonical_source, "canonical"),
        ("receipt_under_test", funnel_source, "producer funnel"),
    )
    for input_key, expected_source, label in bindings:
        binding = inputs.get(input_key) or {}
        raw_path = binding.get("path")
        _require(isinstance(raw_path, str) and bool(raw_path), f"funnel verifier {label} path 缺失", errors)
        if isinstance(raw_path, str) and raw_path:
            recorded_path = _resolve_recorded_path(raw_path, verification_source.path)
            _require(recorded_path == expected_source.path, f"funnel verifier {label} path binding 不符", errors)
        _require(
            binding.get("sha256") == expected_source.sha256,
            f"funnel verifier {label} SHA-256 binding 不符",
            errors,
        )

    recomputed_aggregate = ((payload.get("recomputed") or {}).get("aggregate") or {})
    producer_aggregate = (funnel_source.payload or {}).get("aggregate") or {}
    count_keys = (
        "raw_clairs_records",
        "longphase_s_recalibrated_pass",
        "autosomal_biallelic_sSNV",
        "retained_sSNV",
    )
    _require(
        all(recomputed_aggregate.get(key) == producer_aggregate.get(key) for key in count_keys),
        "funnel verifier recomputed aggregate counts 與 producer funnel 不符",
        errors,
    )
    _require(
        recomputed_aggregate.get("branch_counts") == producer_aggregate.get("branch_counts"),
        "funnel verifier recomputed branch_counts 與 producer funnel 不符",
        errors,
    )
    _require(
        recomputed_aggregate.get("relative_ratios") == producer_aggregate.get("relative_ratios"),
        "funnel verifier recomputed relative_ratios 與 producer funnel 不符",
        errors,
    )
    _require(
        recomputed_aggregate.get("total_ratios_over_raw") == producer_aggregate.get("total_ratios_over_raw"),
        "funnel verifier recomputed total_ratios 與 producer funnel 不符",
        errors,
    )
    return errors


def _validate_pilot(payload: Mapping[str, Any], pilot_source: Source) -> list[str]:
    errors: list[str] = []
    _require(payload.get("all_pass") is True, "symbolic/MILP pilot all_pass 不是 true", errors)
    _require(payload.get("scope") == "PILOT_NOT_FINAL_VALIDATION", "pilot scope 未明示 PILOT_NOT_FINAL_VALIDATION", errors)
    symbolic = payload.get("symbolic_exhaustive") or {}
    cross = payload.get("legacy_milp_crosscheck") or {}
    k9 = payload.get("k9_k12") or {}
    likelihood = payload.get("likelihood") or {}
    _require(symbolic.get("mismatches") == 0 and symbolic.get("status") == "PASS", "symbolic exhaustive control 未通過", errors)
    _require(cross.get("n_vertex_set_mismatch") == 0 and cross.get("n_fail") == 0, "frozen-vs-MILP cross-check 未通過", errors)
    _require(k9.get("n_fail") == 0 and k9.get("n_pass") == k9.get("n_cases"), "k9-k12 exact pilot 未全數通過", errors)
    _require(likelihood.get("status") == "PASS", "likelihood control 未通過", errors)
    implementation = payload.get("implementation") or {}
    for path_key, sha_key, label in (
        ("runner", "runner_sha256", "runner"),
        ("research_solver", "research_solver_sha256", "research solver"),
    ):
        raw_path = implementation.get(path_key)
        expected_sha = implementation.get(sha_key)
        _require(isinstance(raw_path, str) and bool(raw_path), f"pilot implementation {label} path 缺失", errors)
        _require(
            isinstance(expected_sha, str) and len(expected_sha) == 64,
            f"pilot implementation {label} SHA-256 缺失",
            errors,
        )
        if isinstance(raw_path, str) and raw_path:
            path = _resolve_recorded_path(raw_path, pilot_source.path)
            if not path.is_file():
                errors.append(f"pilot implementation {label} live file 不存在：{path}")
            elif isinstance(expected_sha, str):
                _require(
                    _sha256(path) == expected_sha,
                    f"pilot implementation {label} live SHA-256 不符",
                    errors,
                )
    return errors


def _validate_m0(payload: Mapping[str, Any]) -> tuple[list[str], list[str]]:
    hard: list[str] = []
    final: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.hypercube_m0_likelihood_census_receipt", "M0 schema_name 不符", hard)
    _require(payload.get("all_pass") is True, "M0 all_pass 不是 true", hard)
    aggregate = payload.get("aggregate") or {}
    units = aggregate.get("n_hp_lineage_units")
    statuses = aggregate.get("selection_status_counts") or {}
    _require(_is_nonnegative_int(units), "M0 n_hp_lineage_units 缺失", hard)
    _require(isinstance(statuses, dict) and all(_is_nonnegative_int(v) for v in statuses.values()), "M0 status counts 不完整", hard)
    if _is_nonnegative_int(units) and isinstance(statuses, dict) and all(_is_nonnegative_int(v) for v in statuses.values()):
        _require(sum(statuses.values()) == units, "M0 status counts 不等於 HP lineage units", hard)
    selected = payload.get("selected_datasets") or []
    population = payload.get("population") or {}
    excluded = population.get("excluded_hp_lineage_unit_counts") or {}
    _require(_is_nonnegative_int(population.get("n_primary_mutation_regions")), "M0 primary mutation regions 缺失", hard)
    _require(_is_nonnegative_int(population.get("n_fully_m0_eligible_regions")), "M0 fully eligible regions 缺失", hard)
    _require(_is_nonnegative_int(population.get("n_regions_with_any_eligible_hp_lineage")), "M0 any-eligible regions 缺失", hard)
    _require(_is_nonnegative_int(excluded.get("CAPPED")), "M0 CAPPED HP-unit count 缺失", hard)
    if payload.get("full_task_b_scope") is not True:
        final.append("M0 receipt 不是 full_task_b_scope=true")
    if set(selected) != set(DATASETS) or len(selected) != 7:
        final.append(f"M0 只有 {len(selected)} / 7 datasets")
    by_dataset = aggregate.get("by_dataset") or {}
    if payload.get("full_task_b_scope") is True:
        _require(set(by_dataset) == set(DATASETS), "M0 full receipt 的 by_dataset 不完整", hard)
        by_dataset_statuses: Counter[str] = Counter()
        by_dataset_units = 0
        for dataset in DATASETS:
            cell = by_dataset.get(dataset) or {}
            cell_units = cell.get("n_hp_lineage_units")
            cell_statuses = cell.get("selection_status_counts") or {}
            _require(_is_nonnegative_int(cell_units), f"M0 by_dataset.{dataset} eligible units 缺失", hard)
            _require(
                isinstance(cell_statuses, dict) and all(_is_nonnegative_int(value) for value in cell_statuses.values()),
                f"M0 by_dataset.{dataset} status counts 不完整",
                hard,
            )
            if _is_nonnegative_int(cell_units) and isinstance(cell_statuses, dict) and all(
                _is_nonnegative_int(value) for value in cell_statuses.values()
            ):
                _require(sum(cell_statuses.values()) == cell_units, f"M0 by_dataset.{dataset} status 不守恆", hard)
                by_dataset_units += int(cell_units)
                by_dataset_statuses.update(cell_statuses)
        if _is_nonnegative_int(units):
            _require(by_dataset_units == units, "M0 by_dataset eligible units 加總不符", hard)
        if isinstance(statuses, dict):
            _require(dict(by_dataset_statuses) == dict(statuses), "M0 by_dataset status 加總不符", hard)
    outputs = payload.get("outputs") or {}
    _require(isinstance(outputs.get("rows"), str) and bool(outputs.get("rows")), "M0 rows path 缺失", hard)
    _require(_is_nonnegative_int(outputs.get("rows_size_bytes")) and outputs.get("rows_size_bytes", 0) > 0, "M0 rows_size_bytes 缺失", hard)
    _require(isinstance(outputs.get("rows_sha256"), str) and len(outputs.get("rows_sha256", "")) == 64, "M0 rows_sha256 缺失", hard)
    return hard, final


def _resolve_recorded_path(raw_path: str, receipt_path: Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path.resolve()
    cwd_candidate = path.resolve()
    if cwd_candidate.exists():
        return cwd_candidate
    return (receipt_path.parent / path).resolve()


def _semantic_json_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _stable_regular_identity(path: Path, label: str) -> tuple[bytes, os.stat_result]:
    """Read one regular non-symlink file and reject an in-read identity change."""

    try:
        before = path.lstat()
    except OSError as exc:
        raise ReportGateError(f"{label} 不存在或無法 stat：{path}：{exc}") from exc
    if not stat.S_ISREG(before.st_mode) or path.is_symlink():
        raise ReportGateError(f"{label} 不是 regular non-symlink file：{path}")
    try:
        payload = path.read_bytes()
        after = path.lstat()
    except OSError as exc:
        raise ReportGateError(f"{label} 無法穩定讀取：{path}：{exc}") from exc
    before_id = (
        before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
        before.st_ctime_ns,
    )
    after_id = (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if before_id != after_id or len(payload) != before.st_size:
        raise ReportGateError(f"{label} 在讀取期間發生 identity/size 變動：{path}")
    return payload, before


def _require_immutable_release_file(
    observed: os.stat_result,
    label: str,
    errors: list[str],
) -> None:
    """Require the physical single-link 0444 boundary promised by S13."""

    _require(
        stat.S_IMODE(observed.st_mode) == 0o444,
        f"{label} physical mode 不是 0444",
        errors,
    )
    _require(observed.st_nlink == 1, f"{label} physical st_nlink 不是 1", errors)


def _require_immutable_release_parents(
    contract_root: Path,
    parent: Path,
    label: str,
    errors: list[str],
    checked: set[Path],
) -> None:
    """Authenticate every lexical parent from one release file to its root."""

    root = contract_root.absolute()
    current = parent.absolute()
    try:
        current.relative_to(root)
    except ValueError:
        errors.append(f"{label} parent 逃逸 release contract root")
        return
    while True:
        if current not in checked:
            try:
                observed = current.lstat()
            except OSError as exc:
                errors.append(f"{label} parent 無法 stat：{current}：{exc}")
                return
            _require(
                stat.S_ISDIR(observed.st_mode)
                and not stat.S_ISLNK(observed.st_mode),
                f"{label} parent 不是 regular non-symlink directory：{current}",
                errors,
            )
            _require(
                stat.S_IMODE(observed.st_mode) == 0o555,
                f"{label} parent contract mode 不是 0555：{current}",
                errors,
            )
            checked.add(current)
        if current == root:
            break
        current = current.parent


def _release_manifest_source(
    verification_source: Source,
) -> tuple[Source | None, list[str]]:
    """Physically authenticate the release manifest persisted through S13.

    S13's booleans are necessary but not sufficient for presentation.  This
    reader follows the persisted path, authenticates the adjacent sidecar,
    checks canonical authority/scope, and re-hashes every immutable snapshot
    copy using the size/SHA identities stored in the manifest.
    """

    errors: list[str] = []
    binding = (verification_source.payload or {}).get("release_binding") or {}
    identity = binding.get("release_manifest") or {}
    raw_path = identity.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        return None, ["M2 S13 release_manifest path 缺失"]
    path = _resolve_recorded_path(raw_path, verification_source.path)
    contract_root = path.parent
    checked_release_dirs: set[Path] = set()
    try:
        raw, observed = _stable_regular_identity(path, "M2 frozen release manifest")
        document = json.loads(raw.decode("utf-8"))
    except ReportGateError as exc:
        return None, [str(exc)]
    except (UnicodeError, json.JSONDecodeError) as exc:
        return None, [f"M2 frozen release manifest JSON 無法解析：{exc}"]
    _require_immutable_release_file(observed, "M2 frozen release manifest", errors)
    _require_immutable_release_parents(
        contract_root, path.parent, "M2 frozen release manifest", errors,
        checked_release_dirs,
    )
    digest = hashlib.sha256(raw).hexdigest()
    _require(digest == identity.get("sha256"), "M2 frozen release manifest physical SHA-256 與 S13 不符", errors)
    _require(
        _semantic_json_sha256(document) == identity.get("semantic_sha256"),
        "M2 frozen release manifest semantic SHA-256 與 S13 不符",
        errors,
    )

    side_identity = identity.get("sidecar") or {}
    sidecar = path.with_name(f"{path.name}.sha256")
    _require(
        isinstance(side_identity.get("path"), str)
        and _resolve_recorded_path(side_identity["path"], verification_source.path) == sidecar.resolve(),
        "M2 frozen release sidecar path 不是 manifest 相鄰固定檔名",
        errors,
    )
    try:
        side_raw, side_stat = _stable_regular_identity(sidecar, "M2 frozen release manifest sidecar")
        side_fields = side_raw.decode("ascii", errors="strict").strip().split()
    except (ReportGateError, UnicodeError) as exc:
        errors.append(str(exc))
        side_raw = b""
        side_fields = []
    else:
        _require_immutable_release_file(
            side_stat, "M2 frozen release manifest sidecar", errors
        )
        _require_immutable_release_parents(
            contract_root, sidecar.parent, "M2 frozen release manifest sidecar",
            errors, checked_release_dirs,
        )
    _require(side_fields == [digest, path.name], "M2 frozen release sidecar 未精確綁定 manifest 檔名與內容", errors)
    _require(
        hashlib.sha256(side_raw).hexdigest() == side_identity.get("sha256"),
        "M2 frozen release sidecar physical SHA-256 與 S13 不符",
        errors,
    )

    expected_scope = {
        "technical_datasets": 7,
        "biological_samples": 6,
        "chromosomes": 22,
        "chromosome_names": list(AUTOSOMES),
        "tasks": 154,
        "datasets": list(DATASETS),
    }
    _require(isinstance(document, dict), "M2 frozen release manifest root 不是 object", errors)
    if isinstance(document, dict):
        _require(
            document.get("schema_name") == "intersubmod.m2_release_run_manifest"
            and document.get("schema_version") == "1.0.0"
            and document.get("task_type") == "B_COMPREHENSIVE_VALIDATION",
            "M2 frozen release manifest schema/task type 不符",
            errors,
        )
        _require(
            document.get("authority_mode") == binding.get("authority_mode") == "CANONICAL_V5_FROZEN"
            and document.get("validation_evidence_eligible") is True
            and binding.get("validation_evidence_eligible") is True,
            "M2 frozen release manifest 不是 canonical validation authority",
            errors,
        )
        _require(document.get("all_pass") is True, "M2 frozen release manifest all_pass 不是 true", errors)
        _require(document.get("scope") == expected_scope, "M2 frozen release manifest exact 7×chr1-22 scope 不符", errors)
        manifest_checks = document.get("checks") or {}
        _require(bool(manifest_checks) and all(value is True for value in manifest_checks.values()), "M2 frozen release manifest checks 未全數 true", errors)
        _require(
            document.get("receipt_integrity") == {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": sidecar.name,
                "covers": path.name,
            },
            "M2 frozen release manifest receipt_integrity 不符",
            errors,
        )

        snapshot = document.get("source_snapshot") or {}
        entries = snapshot.get("entries") or []
        bound_sources = binding.get("snapshot_sources") or {}
        _require(
            snapshot.get("n_files") == 11
            and isinstance(entries, list)
            and len(entries) == 11
            and isinstance(bound_sources, dict)
            and len(bound_sources) == 11,
            "M2 frozen release snapshot 不是 exact-11",
            errors,
        )
        seen_roles: set[str] = set()
        if isinstance(entries, list) and isinstance(bound_sources, dict):
            for entry in entries:
                if not isinstance(entry, dict):
                    errors.append("M2 frozen release snapshot entry 不是 object")
                    continue
                role = entry.get("role")
                copy_record = entry.get("snapshot") or {}
                relative = copy_record.get("path")
                if not isinstance(role, str) or not role or role in seen_roles:
                    errors.append("M2 frozen release snapshot role 缺失或重複")
                    continue
                seen_roles.add(role)
                if not isinstance(relative, str) or not relative:
                    errors.append(f"M2 frozen release snapshot/{role} path 缺失")
                    continue
                relative_path = Path(relative)
                if relative_path.is_absolute() or ".." in relative_path.parts:
                    errors.append(f"M2 frozen release snapshot/{role} path 逃逸 contract root")
                    continue
                physical = path.parent / relative_path
                try:
                    physical.resolve(strict=True).relative_to(contract_root.resolve(strict=True))
                except (OSError, ValueError):
                    errors.append(f"M2 frozen release snapshot/{role} physical path 逃逸 contract root")
                    continue
                try:
                    copy_raw, copy_stat = _stable_regular_identity(
                        physical, f"M2 frozen release snapshot/{role}"
                    )
                except ReportGateError as exc:
                    errors.append(str(exc))
                    continue
                _require(
                    copy_record.get("mode_octal") == "0444"
                    and copy_record.get("st_nlink") == 1,
                    f"M2 frozen release snapshot/{role} manifest 未宣告 mode 0444/st_nlink 1",
                    errors,
                )
                _require_immutable_release_file(
                    copy_stat, f"M2 frozen release snapshot/{role}", errors
                )
                _require_immutable_release_parents(
                    contract_root, physical.parent,
                    f"M2 frozen release snapshot/{role}", errors,
                    checked_release_dirs,
                )
                copy_sha = hashlib.sha256(copy_raw).hexdigest()
                bound = bound_sources.get(role) or {}
                _require(copy_stat.st_size == copy_record.get("size_bytes"), f"M2 frozen release snapshot/{role} physical size 不符 manifest", errors)
                _require(copy_sha == copy_record.get("sha256"), f"M2 frozen release snapshot/{role} physical SHA-256 不符 manifest", errors)
                _require(
                    isinstance(bound.get("path"), str)
                    and Path(bound["path"]).resolve() == physical.resolve()
                    and bound.get("sha256") == copy_sha,
                    f"M2 frozen release snapshot/{role} identity 未由 S13 精確綁定",
                    errors,
                )
        _require(seen_roles == set(bound_sources), "M2 frozen release snapshot role set 與 S13 不符", errors)

        canonical = document.get("canonical_manifest") or {}
        canonical_copy = canonical.get("immutable_copy") or {}
        canonical_relative = canonical_copy.get("path")
        if isinstance(canonical_relative, str) and canonical_relative:
            canonical_relative_path = Path(canonical_relative)
            if canonical_relative_path.is_absolute() or ".." in canonical_relative_path.parts:
                errors.append("M2 frozen canonical manifest immutable copy path 逃逸 contract root")
                canonical_path = None
            else:
                canonical_path = path.parent / canonical_relative_path
            if canonical_path is not None:
                try:
                    canonical_path.resolve(strict=True).relative_to(
                        contract_root.resolve(strict=True)
                    )
                except (OSError, ValueError):
                    errors.append("M2 frozen canonical manifest physical path 逃逸 contract root")
                    canonical_path = None
        else:
            canonical_path = None
            errors.append("M2 frozen canonical manifest immutable copy identity 缺失")
        if canonical_path is not None:
            try:
                canonical_raw, canonical_stat = _stable_regular_identity(
                    canonical_path, "M2 frozen canonical manifest copy"
                )
            except ReportGateError as exc:
                errors.append(str(exc))
            else:
                _require(
                    canonical_copy.get("mode_octal") == "0444"
                    and canonical_copy.get("st_nlink") == 1,
                    "M2 frozen canonical manifest copy manifest 未宣告 mode 0444/st_nlink 1",
                    errors,
                )
                _require_immutable_release_file(
                    canonical_stat, "M2 frozen canonical manifest copy", errors
                )
                _require_immutable_release_parents(
                    contract_root, canonical_path.parent,
                    "M2 frozen canonical manifest copy", errors,
                    checked_release_dirs,
                )
                canonical_sha = hashlib.sha256(canonical_raw).hexdigest()
                canonical_binding = binding.get("canonical_input_manifest") or {}
                _require(canonical_stat.st_size == canonical_copy.get("size_bytes"), "M2 frozen canonical manifest copy size 不符", errors)
                _require(canonical_sha == canonical_copy.get("sha256") == canonical.get("expected_sha256"), "M2 frozen canonical manifest copy SHA/authority 不符", errors)
                _require(
                    isinstance(canonical_binding.get("path"), str)
                    and Path(canonical_binding["path"]).resolve() == canonical_path.resolve()
                    and canonical_binding.get("sha256") == canonical_sha,
                    "M2 S13 canonical_input_manifest identity 不符 physical copy",
                    errors,
                )

    deep = binding.get("deep_release_verification") or {}
    _require(
        deep.get("mode") == "FROZEN_FREEZER_VERIFY_RELEASE_CONTRACT"
        and deep.get("release_manifest_sha256") == digest
        and deep.get("verified_snapshot_files") == 11
        and deep.get("all_pass") is True,
        "M2 S13 deep_release_verification identity/scope 不完整",
        errors,
    )
    freezer_path_raw = deep.get("freezer_path")
    if isinstance(freezer_path_raw, str) and freezer_path_raw and _is_sha256(deep.get("freezer_sha256")):
        freezer_path = Path(freezer_path_raw)
        if not freezer_path.is_absolute():
            freezer_path = verification_source.path.parent / freezer_path
        try:
            freezer_raw, freezer_stat = _stable_regular_identity(
                freezer_path, "M2 frozen deep verifier/freezer"
            )
        except ReportGateError as exc:
            errors.append(str(exc))
        else:
            _require_immutable_release_file(
                freezer_stat, "M2 frozen deep verifier/freezer", errors
            )
            _require_immutable_release_parents(
                contract_root, freezer_path.parent,
                "M2 frozen deep verifier/freezer", errors,
                checked_release_dirs,
            )
            _require(hashlib.sha256(freezer_raw).hexdigest() == deep["freezer_sha256"], "M2 frozen deep verifier/freezer physical SHA-256 不符", errors)
    else:
        errors.append("M2 S13 deep_release_verification freezer identity 缺失")

    if errors:
        return None, errors
    return (
        Source(
            source_id="S17",
            label="Physically authenticated canonical M2 frozen release manifest",
            path=path,
            sha256=digest,
            scope=(
                f"CANONICAL_V5_FROZEN；7 datasets/6 biological samples×chr1-22；"
                f"exact-11 snapshots；manifest_size={observed.st_size} bytes"
            ),
            payload=document,
        ),
        [],
    )


def _resource_session_source(
    verification_source: Source,
) -> tuple[Source | None, list[str]]:
    """Load the persisted extraction-session resource-gate attestation.

    The independent verifier already authenticates every session/batch gate.
    The presentation layer re-reads the session receipt by the path+SHA bound in
    that verifier so the professor-facing zero-conflict values are not prose-only.
    """

    errors: list[str] = []
    extraction = (verification_source.payload or {}).get("extraction") or {}
    orchestration = extraction.get("orchestration") or {}
    identity = orchestration.get("session_receipt") or {}
    raw_path = identity.get("path")
    expected_sha = identity.get("sha256")
    if not isinstance(raw_path, str) or not raw_path or not _is_sha256(expected_sha):
        return None, ["M2 verifier extraction orchestration session identity 不完整"]
    path = _resolve_recorded_path(raw_path, verification_source.path)
    if not path.is_file():
        return None, [f"M2 extraction orchestration session receipt 不存在：{path}"]
    try:
        source = _source(
            "S15",
            "M2 extraction zero-conflict resource-gate session attestation",
            path,
            "formal extraction session 啟動前；zero conflicting processes/roots；非 RSS/CPU/disk telemetry",
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return None, [f"M2 extraction resource session receipt 無法讀取：{exc}"]
    _require(source.sha256 == expected_sha, "M2 extraction resource session receipt SHA-256 與 verifier 綁定不符", errors)
    session = source.payload or {}
    expected_session_keys = {
        "schema_name", "schema_version", "stage", "session_id", "session_nonce",
        "created_at_utc", "created_monotonic_ns", "host_boot_id",
        "release_manifest", "release_binding_semantic_sha256",
        "run_contract_semantic_sha256", "scope", "output_root",
        "producer_sources", "scheduler_policy", "parent_extraction",
        "resource_gate", "receipt_integrity",
    }
    _require(set(session) == expected_session_keys, "M2 extraction resource session schema keys 漂移", errors)
    _require(
        session.get("schema_name") == "intersubmod.m2_orchestration_session"
        and session.get("schema_version") == "1.0.0"
        and session.get("stage") == "extraction",
        "M2 extraction resource session schema/stage 不符",
        errors,
    )
    _require(session.get("session_id") == orchestration.get("session_id"), "M2 extraction resource session_id 與 verifier 不符", errors)
    _require(session.get("parent_extraction") is None, "M2 extraction resource session 不應有 parent extraction", errors)
    scope = session.get("scope") or {}
    _require(
        scope == {
            "datasets": list(DATASETS),
            "chromosomes": list(AUTOSOMES),
            "expected_tasks": 154,
        },
        "M2 extraction resource session scope 不是 7×chr1-22",
        errors,
    )
    gate = session.get("resource_gate") or {}
    expected_gate_keys = {
        "ignore_resource_gate", "process_count", "root_count", "representatives",
        "observed_at_utc", "observed_monotonic_ns", "host_boot_id",
        "process_snapshot_sha256",
    }
    _require(set(gate) == expected_gate_keys, "M2 extraction resource_gate keys 漂移", errors)
    _require(gate.get("ignore_resource_gate") is False, "M2 extraction resource gate 被 bypass", errors)
    _require(
        gate.get("process_count") == 0
        and gate.get("root_count") == 0
        and gate.get("representatives") == [],
        "M2 extraction resource gate 不是 zero-conflict attestation",
        errors,
    )
    snapshot = {
        "process_count": gate.get("process_count"),
        "root_count": gate.get("root_count"),
        "representatives": gate.get("representatives"),
    }
    _require(
        gate.get("process_snapshot_sha256") == _semantic_json_sha256(snapshot),
        "M2 extraction resource gate process_snapshot_sha256 不符",
        errors,
    )
    _require(
        isinstance(gate.get("observed_at_utc"), str)
        and bool(gate.get("observed_at_utc"))
        and _is_nonnegative_int(gate.get("observed_monotonic_ns")),
        "M2 extraction resource gate observation time 無效",
        errors,
    )
    _require(
        gate.get("host_boot_id") == session.get("host_boot_id")
        and _is_nonnegative_int(session.get("created_monotonic_ns"))
        and _is_nonnegative_int(gate.get("observed_monotonic_ns"))
        and session.get("created_monotonic_ns") >= gate.get("observed_monotonic_ns"),
        "M2 extraction resource gate 與 session 的 host/time 綁定不符",
        errors,
    )
    return source, errors


def _m0_rows_source(m0_source: Source) -> tuple[Source | None, list[str]]:
    errors: list[str] = []
    outputs = (m0_source.payload or {}).get("outputs") or {}
    raw_path = outputs.get("rows")
    if not isinstance(raw_path, str) or not raw_path:
        return None, ["M0 rows path 缺失"]
    path = _resolve_recorded_path(raw_path, m0_source.path)
    if not path.is_file():
        return None, [f"M0 rows TSV 不存在：{path}"]
    _require(path.stat().st_size == outputs.get("rows_size_bytes"), "M0 rows TSV size_bytes 不符", errors)
    actual_sha = _sha256(path)
    _require(actual_sha == outputs.get("rows_sha256"), "M0 rows TSV SHA-256 不符", errors)
    if errors:
        return None, errors
    return (
        Source(
            source_id="S12",
            label="M0 HP-unit census rows",
            path=path,
            sha256=actual_sha,
            scope="64,973 eligible HP-lineage-unit row-level census",
            payload=None,
        ),
        [],
    )


def _validate_m0_verification(payload: Mapping[str, Any], m0_source: Source, rows_source: Source) -> list[str]:
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.hypercube_m0_independent_verification", "M0 verifier schema_name 不符", errors)
    _require(payload.get("schema_version") == "1.1.0", "M0 verifier schema_version 不是 1.1.0", errors)
    _require(payload.get("verdict") == "PASS", "M0 verifier verdict 不是 PASS", errors)
    _require(payload.get("candidate_mode") == "deep", "M0 verifier 不是 deep candidate mode", errors)
    scope = payload.get("scope") or {}
    _require(scope.get("full_task_b_scope") is True, "M0 verifier 不是 full Task-B", errors)
    _require(tuple(scope.get("selected_datasets") or ()) == DATASETS, "M0 verifier datasets/order 不符", errors)
    _require(scope.get("row_schema_matches") is True, "M0 verifier row schema 未通過", errors)
    inputs = payload.get("inputs") or {}
    _require(inputs.get("receipt_sha256") == m0_source.sha256, "M0 verifier receipt SHA-256 不符", errors)
    _require(inputs.get("receipt_size_bytes") == m0_source.path.stat().st_size, "M0 verifier receipt size 不符", errors)
    _require(inputs.get("rows_sha256") == rows_source.sha256, "M0 verifier rows SHA-256 不符", errors)
    _require(inputs.get("rows_size_bytes") == rows_source.path.stat().st_size, "M0 verifier rows size 不符", errors)
    checks = payload.get("checks") or {}
    _require(bool(checks) and checks.get("n_errors") == 0, "M0 verifier n_errors 非 0", errors)
    _require(
        bool(checks) and all(value is True for key, value in checks.items() if key != "n_errors"),
        "M0 verifier checks 未全數 true",
        errors,
    )
    recomputed = payload.get("independently_recomputed_aggregate") or {}
    m0_aggregate = (m0_source.payload or {}).get("aggregate") or {}
    _require(recomputed.get("n_hp_lineage_units") == m0_aggregate.get("n_hp_lineage_units"), "M0 verifier eligible units 與 receipt 不符", errors)
    _require(recomputed.get("selection_status_counts") == m0_aggregate.get("selection_status_counts"), "M0 verifier status counts 與 receipt 不符", errors)
    _require(recomputed.get("by_dataset") == m0_aggregate.get("by_dataset"), "M0 verifier by_dataset 與 receipt 不符", errors)
    categorical = payload.get("categorical_conservation") or {}
    _require(categorical.get("partition_conserves") is True, "M0 verifier selection partition 未守恆", errors)
    _require(categorical.get("T_V_partition_conserves") is True, "M0 verifier T/V partition 未守恆", errors)
    reconciliation = payload.get("canonical_reconciliation") or {}
    _require(reconciliation.get("eligible_plus_excluded_equals_primary_hp_units") is True, "M0 verifier eligible+capped 與 canonical 不守恆", errors)
    _require(reconciliation.get("missing_tsv_units") == 0 and reconciliation.get("extra_tsv_units") == 0, "M0 verifier TSV unit 集合不符", errors)
    _require(reconciliation.get("n_candidate_units_deep_checked") == m0_aggregate.get("n_hp_lineage_units"), "M0 verifier 未 deep-check 全部 eligible units", errors)
    return errors


def _validate_m2_extraction(payload: Mapping[str, Any], *, require_ps_aware: bool = True) -> list[str]:
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.m2_full_extraction_receipt", "M2 extraction schema_name 不符", errors)
    if require_ps_aware:
        _require(payload.get("schema_version") == "1.2.0", "M2 extraction full schema_version 不是 PS-aware 1.2.0", errors)
    _require(payload.get("all_pass") is True, "M2 extraction all_pass 不是 true", errors)
    scope = payload.get("scope") or {}
    _require(tuple(scope.get("datasets") or ()) == DATASETS, "M2 extraction datasets/sort 不符", errors)
    _require(tuple(scope.get("chromosomes") or ()) == AUTOSOMES, "M2 extraction chromosomes 不是 chr1-22", errors)
    _require(scope.get("expected_tasks") == 154, "M2 extraction expected_tasks 不是 154", errors)
    if require_ps_aware:
        _require(scope.get("child_schema_version") == "1.2.0", "M2 extraction child receipts 不是 schema 1.2.0", errors)
    status = payload.get("task_status_counts") or {}
    passing = int(status.get("PASS", 0)) + int(status.get("REUSED_PASS", 0))
    _require(set(status) <= {"PASS", "REUSED_PASS"}, "M2 extraction status 含 PASS/REUSED_PASS 以外類別", errors)
    _require(passing == 154 and sum(int(value) for value in status.values()) == 154, "M2 extraction 不是恰好 154/154 PASS", errors)
    _require(payload.get("n_results") == 154, "M2 extraction n_results 不是 154", errors)
    task_index = payload.get("results") or []
    _require(isinstance(task_index, list) and len(task_index) == 154, "M2 extraction results/task index 不是 154 rows", errors)
    if isinstance(task_index, list) and len(task_index) == 154:
        observed = [(row.get("dataset"), row.get("chrom")) for row in task_index]
        expected = [(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES]
        _require(observed == expected, "M2 extraction results/task index scope/order 不符", errors)
        _require(len(set(observed)) == 154, "M2 extraction results/task index 有重複 task", errors)
        for row in task_index:
            _require(row.get("status") in {"PASS", "REUSED_PASS"}, "M2 extraction task index 含非 passing status", errors)
            receipt = row.get("receipt") or {}
            if require_ps_aware:
                _require(receipt.get("schema_version") == "1.2.0", "M2 extraction task index child schema 非 1.2.0", errors)
    checks = payload.get("checks") or {}
    if require_ps_aware:
        _require(checks.get("all_child_receipts_schema_1_2") is True, "M2 extraction 未證明 154 child receipts 皆為 schema 1.2.0", errors)
    _require(bool(checks) and all(value is True for value in checks.values()), "M2 extraction checks 未全數 true", errors)
    aggregate = payload.get("aggregate") or {}
    by_dataset = aggregate.get("by_dataset") or {}
    _require(set(by_dataset) == set(DATASETS), "M2 extraction by_dataset 不完整", errors)
    summaries = aggregate.get("component_summary_by_linkage_basis") or {}
    if require_ps_aware:
        _require(all(basis in summaries for basis in ("PS_HP1", "PS_HP2")), "M2 extraction 缺 PS_HP1/PS_HP2 primary component summaries", errors)
        for basis in ("PS_HP1", "PS_HP2"):
            cells = summaries.get(basis) or {}
            _require(set(cells) == {"1", "2", "3", "5"}, f"M2 extraction {basis} threshold grid 不完整", errors)
            for threshold, cell in cells.items():
                distribution = cell.get("k_component_sites_distribution") or {}
                required_component_counts = ("n_components", "n_singletons_k1", "n_multisite_k_gt1", "max_k_component_sites")
                _require(all(_is_nonnegative_int(cell.get(key)) for key in required_component_counts), f"M2 extraction {basis}/{threshold} component metrics 缺失", errors)
                _require(isinstance(distribution, dict) and all(str(key).isdigit() and _is_nonnegative_int(value) for key, value in distribution.items()), f"M2 extraction {basis}/{threshold} k distribution 無效", errors)
                if all(_is_nonnegative_int(cell.get(key)) for key in required_component_counts) and isinstance(distribution, dict) and all(
                    str(key).isdigit() and _is_nonnegative_int(value) for key, value in distribution.items()
                ):
                    _require(sum(distribution.values()) == cell["n_components"], f"M2 extraction {basis}/{threshold} k distribution 不守恆", errors)
                    _require(cell["n_singletons_k1"] + cell["n_multisite_k_gt1"] == cell["n_components"], f"M2 extraction {basis}/{threshold} k=1/k>1 不守恆", errors)
                    expected_max = max((int(key) for key, value in distribution.items() if value), default=0)
                    _require(cell["max_k_component_sites"] == expected_max, f"M2 extraction {basis}/{threshold} max_k 不是 distribution maximum", errors)
    counts = aggregate.get("counts") or {}
    required_counts = (
        "raw_overlapping_alignments",
        "canonical_eligible_alignments",
        "unique_molecule_ids",
        "sidecar_exact_matches",
        "sidecar_missing",
        "fixed_ra_calls",
        "alt_calls",
    )
    _require(all(_is_nonnegative_int(counts.get(key)) for key in required_counts), "M2 extraction funnel counts 不完整", errors)
    if all(_is_nonnegative_int(counts.get(key)) for key in required_counts):
        _require(counts["raw_overlapping_alignments"] >= counts["canonical_eligible_alignments"], "M2 extraction raw < eligible", errors)
        _require(counts["canonical_eligible_alignments"] == counts["unique_molecule_ids"], "M2 extraction eligible alignment != unique molecule", errors)
        _require(counts["sidecar_exact_matches"] == counts["canonical_eligible_alignments"], "M2 extraction sidecar exact != eligible", errors)
        _require(counts["sidecar_missing"] == 0, "M2 extraction sidecar_missing 非 0", errors)
        _require(counts["fixed_ra_calls"] >= counts["alt_calls"], "M2 extraction ALT calls > fixed R/A calls", errors)
    if set(by_dataset) == set(DATASETS):
        for dataset in DATASETS:
            cell = by_dataset.get(dataset) or {}
            task_counts = cell.get("task_status_counts") or {}
            dataset_counts = cell.get("counts") or {}
            _require(set(task_counts) <= {"PASS", "REUSED_PASS"}, f"M2 extraction by_dataset.{dataset} 含非 passing task status", errors)
            _require(sum(int(value) for value in task_counts.values()) == 22, f"M2 extraction by_dataset.{dataset} 不是 22 chromosomes", errors)
            _require(all(_is_nonnegative_int(dataset_counts.get(key)) for key in required_counts), f"M2 extraction by_dataset.{dataset} counts 不完整", errors)
        for key in required_counts:
            if all(_is_nonnegative_int((by_dataset[dataset].get("counts") or {}).get(key)) for dataset in DATASETS):
                _require(
                    sum(by_dataset[dataset]["counts"][key] for dataset in DATASETS) == counts.get(key),
                    f"M2 extraction by_dataset {key} 加總不符 aggregate",
                    errors,
                )
    by_basis = aggregate.get("component_summary_by_linkage_basis") or {}
    _require("pooled" in by_basis and isinstance(by_basis.get("pooled"), dict), "M2 extraction 缺 pooled component summaries", errors)
    if isinstance(by_basis.get("pooled"), dict):
        _require(set(by_basis["pooled"]) >= {"1", "2", "3", "5"}, "M2 extraction 缺 bridge thresholds 1/2/3/5", errors)
    return errors


def _count_mapping(value: Any) -> bool:
    return isinstance(value, dict) and all(
        isinstance(key, str) and _is_nonnegative_int(count)
        for key, count in value.items()
    )


def _is_power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


def _validate_partial_completion_grain(
    grain: Any, label: str, errors: list[str]
) -> None:
    required = {
        "denominator", "full", "partial", "u_number_of_X_distribution",
        "conceptual_completions_2_pow_u_distribution",
        "conceptual_completions_weighted_total",
        "observed_alt_effective_completions_distribution",
        "observed_alt_effective_completions_weighted_total",
        "observed_alt_effective_zero_due_to_fixed_alt_outside_structural_universe",
    }
    _require(isinstance(grain, dict) and set(grain) == required, f"{label} keys 不完整", errors)
    if not isinstance(grain, dict) or set(grain) != required:
        return
    scalar_keys = (
        "denominator", "full", "partial",
        "conceptual_completions_weighted_total",
        "observed_alt_effective_completions_weighted_total",
        "observed_alt_effective_zero_due_to_fixed_alt_outside_structural_universe",
    )
    _require(all(_is_nonnegative_int(grain.get(key)) for key in scalar_keys), f"{label} scalar count 型別錯誤", errors)
    u_distribution = grain["u_number_of_X_distribution"]
    conceptual_distribution = grain["conceptual_completions_2_pow_u_distribution"]
    effective_distribution = grain["observed_alt_effective_completions_distribution"]
    _require(
        _count_mapping(u_distribution)
        and all(key.isdigit() for key in u_distribution),
        f"{label} u distribution 型別錯誤",
        errors,
    )
    _require(
        _count_mapping(conceptual_distribution)
        and all(key.isdigit() and _is_power_of_two(int(key)) for key in conceptual_distribution),
        f"{label} conceptual completion distribution 不是 2^u counts",
        errors,
    )
    _require(
        _count_mapping(effective_distribution)
        and all(
            key.isdigit() and (int(key) == 0 or _is_power_of_two(int(key)))
            for key in effective_distribution
        ),
        f"{label} effective completion distribution 含非 0/2^m key",
        errors,
    )
    if not (
        all(_is_nonnegative_int(grain.get(key)) for key in scalar_keys)
        and _count_mapping(u_distribution)
        and _count_mapping(conceptual_distribution)
        and _count_mapping(effective_distribution)
    ):
        return
    denominator = grain["denominator"]
    _require(grain["full"] + grain["partial"] == denominator, f"{label} full + partial != denominator", errors)
    _require(sum(u_distribution.values()) == denominator, f"{label} u distribution 不守恆", errors)
    _require(sum(conceptual_distribution.values()) == denominator, f"{label} conceptual distribution 不守恆", errors)
    _require(sum(effective_distribution.values()) == denominator, f"{label} effective distribution 不守恆", errors)
    expected_conceptual = {
        str(1 << int(u)): count for u, count in u_distribution.items()
    }
    _require(
        conceptual_distribution == expected_conceptual,
        f"{label} conceptual distribution 與 u→2^u 不符",
        errors,
    )
    _require(grain["full"] == u_distribution.get("0", 0), f"{label} full != u=0", errors)
    conceptual_weighted = sum(int(key) * count for key, count in conceptual_distribution.items())
    effective_weighted = sum(int(key) * count for key, count in effective_distribution.items())
    _require(
        grain["conceptual_completions_weighted_total"] == conceptual_weighted,
        f"{label} conceptual weighted total 無法由 distribution 重算",
        errors,
    )
    _require(
        grain["observed_alt_effective_completions_weighted_total"] == effective_weighted,
        f"{label} effective weighted total 無法由 distribution 重算",
        errors,
    )
    _require(effective_weighted <= conceptual_weighted, f"{label} effective completions > conceptual completions", errors)
    _require(
        grain["observed_alt_effective_zero_due_to_fixed_alt_outside_structural_universe"]
        == effective_distribution.get("0", 0),
        f"{label} effective-zero count 不符 distribution",
        errors,
    )


def _validate_partial_pattern_funnel(
    funnel: Any, cell: Mapping[str, Any], label: str, errors: list[str]
) -> None:
    required = {
        "definitions", *PARTIAL_FUNNEL_GRAINS, "units_denominator",
        "units_with_partial_structural_groups", "partial_group_coverage_denominator",
        "partial_groups_covered", "partial_groups_unsatisfied",
    }
    _require(isinstance(funnel, dict) and set(funnel) == required, f"M2 ranking {label} partial_pattern_funnel keys 不完整", errors)
    if not isinstance(funnel, dict) or set(funnel) != required:
        return
    _require(
        funnel.get("definitions") == EXPECTED_PARTIAL_FUNNEL_DEFINITIONS,
        f"M2 ranking {label} partial completion definitions 漂移",
        errors,
    )
    for grain_name in PARTIAL_FUNNEL_GRAINS:
        _validate_partial_completion_grain(
            funnel[grain_name], f"M2 ranking {label} partial/{grain_name}", errors
        )
    scalar_keys = (
        "units_denominator", "units_with_partial_structural_groups",
        "partial_group_coverage_denominator", "partial_groups_covered",
        "partial_groups_unsatisfied",
    )
    _require(all(_is_nonnegative_int(funnel.get(key)) for key in scalar_keys), f"M2 ranking {label} partial funnel scalar 型別錯誤", errors)
    if not all(_is_nonnegative_int(funnel.get(key)) for key in scalar_keys):
        return
    _require(funnel["units_denominator"] == cell.get("n_component_hp_units"), f"M2 ranking {label} partial units denominator 不符", errors)
    _require(funnel["units_with_partial_structural_groups"] <= funnel["units_denominator"], f"M2 ranking {label} partial-unit count > units", errors)
    _require(
        funnel["partial_groups_covered"] + funnel["partial_groups_unsatisfied"]
        == funnel["partial_group_coverage_denominator"],
        f"M2 ranking {label} partial group coverage 不守恆",
        errors,
    )
    _require(funnel["partial_groups_unsatisfied"] == 0, f"M2 ranking {label} partial groups 有 unsatisfied", errors)
    denominator_links = {
        "unique_patterns": "unique_rax_pattern_groups",
        "quality_groups": "bq_quality_pattern_groups",
        "molecule_projections": "molecule_projections",
    }
    old_denominators = cell.get("partial_pattern_denominators") or {}
    old_u = cell.get("partial_u_distribution") or {}
    for old_name, grain_name in denominator_links.items():
        grain = funnel[grain_name]
        _require(old_denominators.get(old_name) == grain.get("denominator"), f"M2 ranking {label} partial denominator alias {old_name} 不符", errors)
        _require(old_u.get(old_name) == grain.get("u_number_of_X_distribution"), f"M2 ranking {label} partial u alias {old_name} 不符", errors)
    _require(
        funnel["molecule_projections"].get("denominator") == cell.get("informative_scoring_molecules"),
        f"M2 ranking {label} partial molecule denominator != informative scoring molecules",
        errors,
    )
    _require(
        funnel["structural_unique_rax_pattern_groups"].get("partial")
        == cell.get("structural_partial_pattern_groups"),
        f"M2 ranking {label} structural partial groups 不符",
        errors,
    )
    for key in ("partial_group_coverage_denominator", "partial_groups_covered", "partial_groups_unsatisfied"):
        _require(funnel[key] == cell.get(key), f"M2 ranking {label} partial funnel/{key} 不符 top-level", errors)


def _validate_ranking_cell(cell: Mapping[str, Any], label: str, errors: list[str]) -> None:
    for key in REQUIRED_RANKING_CELL_KEYS:
        _require(key in cell, f"M2 ranking {label} 缺 {key}", errors)
    numeric_keys = (
        "n_components", "n_component_hp_units", "molecule_component_projections",
        "informative_scoring_molecules", "all_x_excluded_molecules",
        "structural_retained_molecules", "below_minread_scoring_molecules",
        "bq_scoring_molecules", "solver_complete_units",
        "solver_incomplete_or_not_run_units", "quality_primary_unique_vertex_units",
        "quality_primary_tied_vertex_units", "rank_abstain_units",
        "raw_tree_candidates_T_complete_units", "distinct_vertex_sets_V_complete_units",
        "topology_evaluated_units", "coarse_topology_unique_units",
        "coarse_topology_ambiguous_units", "parent_edge_assignment_unique_units",
        "exact_topology_proven_unique_units", "topology_class_inclusion_counts_denominator",
        "k_component_sites_total", "k_observed_alt_active_total",
        "k_scoring_alt_observed_total", "not_structural_alt_active_sites_total",
        "structural_partial_pattern_groups", "partial_group_coverage_denominator",
        "partial_groups_covered", "partial_groups_unsatisfied",
    )
    flat_mapping_keys = (
        "projected_call_class_counts", "selection_status_counts",
        "candidate_generation_status_counts", "solver_complete_incomplete_counts",
        "k_route_counts", "likelihood_stability_counts",
        "coarse_topology_unique_class_counts", "ambiguous_topology_class_set_counts",
        "topology_derivation_status_counts", "exact_topology_uniqueness_status_counts",
        "partial_pattern_denominators",
    )
    numeric_valid = all(_is_nonnegative_int(cell.get(key)) for key in numeric_keys)
    _require(numeric_valid, f"M2 ranking {label} count 型別錯誤", errors)
    mappings_valid = all(_count_mapping(cell.get(key)) for key in flat_mapping_keys)
    _require(mappings_valid, f"M2 ranking {label} 分類 counts 型別錯誤", errors)
    if numeric_valid:
        n_units = cell["n_component_hp_units"]
        _require(cell["molecule_component_projections"] == cell["informative_scoring_molecules"] + cell["all_x_excluded_molecules"], f"M2 ranking {label} projection funnel 不守恆", errors)
        _require(cell["informative_scoring_molecules"] == cell["structural_retained_molecules"] + cell["below_minread_scoring_molecules"], f"M2 ranking {label} structural/scoring funnel 不守恆", errors)
        _require(cell["bq_scoring_molecules"] == cell["informative_scoring_molecules"], f"M2 ranking {label} BQ scoring 不守恆", errors)
        _require(cell["raw_tree_candidates_T_complete_units"] >= cell["distinct_vertex_sets_V_complete_units"], f"M2 ranking {label} T < V", errors)
        _require(cell["solver_complete_units"] + cell["solver_incomplete_or_not_run_units"] == n_units, f"M2 ranking {label} solver unit partition 不守恆", errors)
        _require(cell["quality_primary_unique_vertex_units"] + cell["quality_primary_tied_vertex_units"] + cell["rank_abstain_units"] == n_units, f"M2 ranking {label} unique/tied/abstain 不守恆", errors)
        _require(cell["topology_evaluated_units"] == cell["coarse_topology_unique_units"] + cell["coarse_topology_ambiguous_units"], f"M2 ranking {label} topology evaluated != unique + ambiguous", errors)
        _require(cell["parent_edge_assignment_unique_units"] == cell["exact_topology_proven_unique_units"], f"M2 ranking {label} parent-edge unique != exact-topology proven unique", errors)
        _require(cell["parent_edge_assignment_unique_units"] <= cell["topology_evaluated_units"], f"M2 ranking {label} parent-edge unique > topology evaluated", errors)
        _require(cell["topology_class_inclusion_counts_denominator"] == cell["topology_evaluated_units"], f"M2 ranking {label} topology inclusion denominator 不符", errors)
        _require(cell["k_component_sites_total"] == cell["k_observed_alt_active_total"] + cell["not_structural_alt_active_sites_total"], f"M2 ranking {label} effective-k site mass 不守恆", errors)
        _require(cell["k_scoring_alt_observed_total"] >= cell["k_observed_alt_active_total"], f"M2 ranking {label} scoring ALT k < structural effective k", errors)
        _require(cell["partial_groups_covered"] + cell["partial_groups_unsatisfied"] == cell["partial_group_coverage_denominator"], f"M2 ranking {label} partial coverage top-level 不守恆", errors)
        _require(cell["partial_groups_unsatisfied"] == 0, f"M2 ranking {label} partial unsatisfied 非 0", errors)
    if mappings_valid and numeric_valid:
        n_units = cell["n_component_hp_units"]
        for field in ("selection_status_counts", "candidate_generation_status_counts", "likelihood_stability_counts", "topology_derivation_status_counts", "exact_topology_uniqueness_status_counts"):
            _require(sum(cell[field].values()) == n_units, f"M2 ranking {label} {field} != HP units", errors)
        solver = cell["solver_complete_incomplete_counts"]
        _require(
            solver == {
                "COMPLETE": cell["solver_complete_units"],
                "INCOMPLETE_OR_ABSTAIN": cell["solver_incomplete_or_not_run_units"],
            },
            f"M2 ranking {label} solver complete/incomplete alias 不符",
            errors,
        )
        expected_k_routes = {"EXACT_COMPLETE", "EXACT_INCOMPLETE", "NOT_RUN_NO_STRUCTURE", "GT_EXACT_LIMIT"}
        routes = cell["k_route_counts"]
        _require(set(routes) == expected_k_routes and sum(routes.values()) == n_units, f"M2 ranking {label} effective-k route partition 不完整", errors)
        _require(routes.get("EXACT_COMPLETE") == cell["solver_complete_units"], f"M2 ranking {label} EXACT_COMPLETE != solver complete", errors)
        _require(sum(value for key, value in routes.items() if key != "EXACT_COMPLETE") == cell["solver_incomplete_or_not_run_units"], f"M2 ranking {label} non-complete k routes != solver incomplete", errors)
        unique_classes = cell["coarse_topology_unique_class_counts"]
        _require(sum(unique_classes.values()) == cell["coarse_topology_unique_units"], f"M2 ranking {label} coarse unique classes 非互斥守恆", errors)
        ambiguous_sets = cell["ambiguous_topology_class_set_counts"]
        _require(sum(ambiguous_sets.values()) == cell["coarse_topology_ambiguous_units"], f"M2 ranking {label} ambiguous class-set counts 不守恆", errors)
        exact_status = cell["exact_topology_uniqueness_status_counts"]
        allowed_exact_statuses = {
            "PROVEN_UNIQUE_BY_SINGLE_WINNING_PARENT_EDGE_ASSIGNMENT",
            "NOT_EVALUATED_CANONICAL_SHAPE_ISOMORPHISM_FOR_MULTIPLE_EDGE_ASSIGNMENTS",
            "ABSTAIN_PRIMARY_LIKELIHOOD_NONCONVERGENCE",
            "NOT_EVALUATED_WITHOUT_COMPLETE_WINNING_VERTEX_SET",
        }
        _require(set(exact_status) <= allowed_exact_statuses, f"M2 ranking {label} exact-topology status 有未定義類別", errors)
        _require(exact_status.get("PROVEN_UNIQUE_BY_SINGLE_WINNING_PARENT_EDGE_ASSIGNMENT", 0) == cell["exact_topology_proven_unique_units"], f"M2 ranking {label} exact-topology proven count 不符", errors)
        _require(exact_status.get("NOT_EVALUATED_CANONICAL_SHAPE_ISOMORPHISM_FOR_MULTIPLE_EDGE_ASSIGNMENTS", 0) == cell["topology_evaluated_units"] - cell["parent_edge_assignment_unique_units"], f"M2 ranking {label} multi-parent exact-topology not-evaluated count 不符", errors)
    partial_denominators = cell.get("partial_pattern_denominators") or {}
    expected_denominators = {"unique_patterns", "quality_groups", "molecule_projections"}
    _require(_count_mapping(partial_denominators) and set(partial_denominators) == expected_denominators, f"M2 ranking {label} partial 三種 denominator 不完整", errors)
    partial_distributions = cell.get("partial_u_distribution") or {}
    _require(isinstance(partial_distributions, dict) and set(partial_distributions) == expected_denominators, f"M2 ranking {label} partial u distributions 不完整", errors)
    if isinstance(partial_distributions, dict) and set(partial_distributions) == expected_denominators and _count_mapping(partial_denominators):
        for denominator_name in sorted(expected_denominators):
            distribution = partial_distributions[denominator_name]
            _require(_count_mapping(distribution) and bool(distribution) and all(u.isdigit() for u in distribution), f"M2 ranking {label} partial u distribution {denominator_name} 型別錯誤", errors)
            if _count_mapping(distribution):
                _require(sum(distribution.values()) == partial_denominators[denominator_name], f"M2 ranking {label} partial u distribution {denominator_name} 不守恆", errors)
    denominator_map = cell.get("denominator_map") or {}
    _require(_count_mapping(denominator_map) and set(denominator_map) == EXPECTED_DENOMINATOR_MAP_KEYS, f"M2 ranking {label} denominator_map 不完整", errors)
    if numeric_valid and _count_mapping(denominator_map) and set(denominator_map) == EXPECTED_DENOMINATOR_MAP_KEYS:
        expected_map = {
            "unit_denominator": cell["n_component_hp_units"],
            "molecule_projection_denominator": cell["molecule_component_projections"],
            "informative_molecule_denominator": cell["informative_scoring_molecules"],
            "solver_complete_unit_denominator": cell["solver_complete_units"],
            "topology_evaluated_unit_denominator": cell["topology_evaluated_units"],
            "partial_group_coverage_denominator": cell["partial_group_coverage_denominator"],
            "topology_inclusion_denominator": cell["topology_class_inclusion_counts_denominator"],
        }
        _require(denominator_map == expected_map, f"M2 ranking {label} denominator_map values 不符 counts", errors)
    checks = cell.get("conservation_checks") or {}
    _require(isinstance(checks, dict) and set(checks) == EXPECTED_CONSERVATION_CHECK_KEYS and all(value is True for value in checks.values()), f"M2 ranking {label} conservation_checks 不完整或未全通過", errors)
    _require(cell.get("all_conserved") is True, f"M2 ranking {label} all_conserved 不是 true", errors)
    _validate_partial_pattern_funnel(cell.get("partial_pattern_funnel"), cell, label, errors)


def _validate_runtime_summary(
    summary: Any, label: str, errors: list[str], *, expected_n: int | None = None
) -> None:
    keys = {"n", "sum", "max", "p50", "p95", "p99"}
    _require(isinstance(summary, dict) and set(summary) == keys, f"{label} runtime summary keys 不完整", errors)
    if not isinstance(summary, dict) or set(summary) != keys:
        return
    n = summary["n"]
    total = summary["sum"]
    _require(_is_nonnegative_int(n) and _is_finite_nonnegative_number(total), f"{label} runtime n/sum 無效", errors)
    if not (_is_nonnegative_int(n) and _is_finite_nonnegative_number(total)):
        return
    if expected_n is not None:
        _require(n == expected_n, f"{label} runtime n 不符 expected units", errors)
    if n == 0:
        _require(total == 0 and all(summary[key] is None for key in ("max", "p50", "p95", "p99")), f"{label} empty runtime summary 不合約", errors)
        return
    ordered = [summary[key] for key in ("p50", "p95", "p99", "max")]
    _require(all(_is_finite_nonnegative_number(value) for value in ordered), f"{label} runtime quantile 無效", errors)
    if all(_is_finite_nonnegative_number(value) for value in ordered):
        _require(ordered == sorted(ordered), f"{label} runtime p50/p95/p99/max 次序錯誤", errors)
        _require(total + 1e-12 >= summary["max"], f"{label} runtime sum < max", errors)


def _validate_runtime_diagnostics(
    runtime: Any, expected_units: int, label: str, errors: list[str]
) -> None:
    required = {
        "schema_name", "schema_version", "clock", "unit", "scope",
        "quantile_definition", "interpretation", "aggregation_memory_contract",
        "n_child_runtime_files", "n_unit_evaluations", "metrics",
        "metrics_when_invoked", "all_child_summaries_recomputed_from_per_unit_rows",
    }
    _require(isinstance(runtime, dict) and set(runtime) == required, f"{label} runtime_diagnostics keys 不完整", errors)
    if not isinstance(runtime, dict) or set(runtime) != required:
        return
    _require(runtime["schema_name"] == "intersubmod.m2_full_primary_runtime_diagnostics" and runtime["schema_version"] == "1.0.0", f"{label} runtime schema 不符", errors)
    _require(runtime["clock"] == "time.monotonic_ns" and runtime["unit"] == "seconds", f"{label} runtime clock/unit 不符", errors)
    _require(runtime["n_child_runtime_files"] == 154, f"{label} runtime 未覆蓋 154 child files", errors)
    _require(runtime["n_unit_evaluations"] == expected_units, f"{label} runtime unit evaluations 不符 aggregate units", errors)
    _require(runtime["all_child_summaries_recomputed_from_per_unit_rows"] is True, f"{label} runtime 未由 per-unit rows 重算", errors)
    _require(all(isinstance(runtime[key], str) and runtime[key] for key in ("scope", "quantile_definition", "interpretation", "aggregation_memory_contract")), f"{label} runtime metadata 缺失", errors)
    metrics = runtime["metrics"]
    invoked = runtime["metrics_when_invoked"]
    _require(isinstance(metrics, dict) and set(metrics) == set(RUNTIME_METRICS), f"{label} runtime metrics 不完整", errors)
    _require(isinstance(invoked, dict) and set(invoked) == set(RUNTIME_METRICS[:2]), f"{label} invoked runtime metrics 不完整", errors)
    if isinstance(metrics, dict):
        for metric in RUNTIME_METRICS:
            _validate_runtime_summary(metrics.get(metric), f"{label}/{metric}", errors, expected_n=expected_units)
    if isinstance(invoked, dict):
        for metric in RUNTIME_METRICS[:2]:
            _validate_runtime_summary(invoked.get(metric), f"{label}/when_invoked/{metric}", errors)
            summary = invoked.get(metric)
            if isinstance(summary, dict) and _is_nonnegative_int(summary.get("n")):
                _require(summary["n"] <= expected_units, f"{label}/when_invoked/{metric} n > all units", errors)


def _validate_m2_ranking(payload: Mapping[str, Any], extraction_sha256: str) -> list[str]:
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.m2_full_ranking_receipt", "M2 ranking schema_name 不符", errors)
    _require(payload.get("schema_version") == "2.0.0", "M2 ranking full schema_version 不是 PS-aware 2.0.0", errors)
    _require(payload.get("all_pass") is True, "M2 ranking all_pass 不是 true", errors)
    scope = payload.get("scope") or {}
    _require(tuple(scope.get("datasets") or ()) == DATASETS, "M2 ranking datasets/sort 不符", errors)
    _require(tuple(scope.get("chromosomes") or ()) == AUTOSOMES, "M2 ranking chromosomes 不是 chr1-22", errors)
    _require(scope.get("expected_tasks") == 154 and scope.get("n_results") == 154, "M2 ranking task scope 不是 154/154", errors)
    _require(scope.get("child_schema_version") == "2.0.0", "M2 ranking child receipts 不是 schema 2.0.0", errors)
    _require(
        scope.get("primary_unit") == "HP_family×known_PS×read_linked_component×threshold",
        "M2 ranking primary_unit 不是 HP_family×known_PS×read_linked_component×threshold",
        errors,
    )
    _require(scope.get("missing_ps_policy") == "SEPARATE_DIAGNOSTIC_NOT_PRIMARY", "M2 ranking missing PS 未分離為 diagnostic", errors)
    _require(
        scope.get("primary_likelihood") == "BQ_SYMMETRIC_SUBSTITUTION_CONDITIONAL_RA",
        "M2 ranking 未明示 BQ symmetric-substitution conditional R/A model",
        errors,
    )
    upstream = payload.get("upstream_extraction_receipt") or {}
    _require(upstream.get("sha256") == extraction_sha256, "M2 ranking 綁定的 extraction receipt SHA-256 不符", errors)
    run_contract = payload.get("run_contract") or {}
    _require(
        run_contract.get("method_contract") == EXPECTED_METHOD_CONTRACT,
        "M2 ranking method contract 不符（partial-read/VAF/edge/missing semantics 漂移）",
        errors,
    )
    checks = payload.get("checks") or {}
    for key in REQUIRED_RANKING_CHECKS:
        _require(checks.get(key) is True, f"M2 ranking checks.{key} 不是 true", errors)
    _require(bool(checks) and all(value is True for value in checks.values()), "M2 ranking 有額外 check 未通過", errors)
    task_index = payload.get("task_index") or []
    _require(isinstance(task_index, list) and len(task_index) == 154, "M2 ranking task_index 不是 154 rows", errors)
    if isinstance(task_index, list) and len(task_index) == 154:
        observed = {(row.get("dataset"), row.get("chrom")) for row in task_index}
        expected = {(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES}
        _require(observed == expected, "M2 ranking task_index 沒有覆蓋完整 dataset×chr", errors)
        for row in task_index:
            _require(row.get("all_pass") is True, "M2 ranking per-chrom receipt 有未通過", errors)
            _require(row.get("schema_version") == "2.0.0", "M2 ranking per-chrom receipt 非 schema 2.0.0", errors)
            _require(row.get("parameters_match_extraction") is True, "M2 ranking per-chrom parameters 不一致", errors)
            _require(row.get("input_hashes_match_extraction") is True, "M2 ranking per-chrom input hashes 不一致", errors)
            _require(row.get("upstream_outputs_verified") is True, "M2 ranking per-chrom upstream hashes 未驗證", errors)
            _require(row.get("no_cross_ps_pattern_pooling") is True, "M2 ranking per-chrom 混合不同 known PS", errors)
            _require(row.get("known_ps_never_mixed") is True, "M2 ranking per-chrom known PS 未隔離", errors)
            _require(row.get("missing_ps_separate_diagnostic") is True, "M2 ranking per-chrom missing PS 未分離", errors)
            _require(row.get("runtime_diagnostics_contract_valid") is True, "M2 ranking per-chrom runtime diagnostics contract 無效", errors)
            _require(row.get("method_contract_matches") is True, "M2 ranking per-chrom method contract 漂移", errors)
            _require(row.get("ranker_source_bound") is True, "M2 ranking per-chrom 未綁定 ranker source", errors)
    aggregate = payload.get("aggregate") or {}
    cells = aggregate.get("by_linkage_basis_and_threshold") or {}
    _require(isinstance(cells, dict) and bool(cells), "M2 ranking aggregate 缺 linkage_basis×threshold", errors)
    if isinstance(cells, dict):
        _require(all(basis in cells for basis in ("PS_HP1", "PS_HP2")), "M2 ranking aggregate 缺 PS_HP1/PS_HP2 primary strata", errors)
        for basis in ("PS_HP1", "PS_HP2"):
            if basis in cells:
                _require(set(cells[basis]) == {"1", "2", "3", "5"}, f"M2 ranking {basis} threshold grid 不是 1/2/3/5", errors)
        for basis, thresholds in cells.items():
            _require(isinstance(thresholds, dict) and bool(thresholds), f"M2 ranking basis {basis} 無 threshold", errors)
            if isinstance(thresholds, dict):
                for threshold, cell in thresholds.items():
                    _require(isinstance(cell, dict), f"M2 ranking {basis}/{threshold} 不是 object", errors)
                    if isinstance(cell, dict):
                        _validate_ranking_cell(cell, f"{basis}/{threshold}", errors)
    by_dataset = payload.get("by_dataset") or {}
    _require(set(by_dataset) == set(DATASETS), "M2 ranking by_dataset 不完整", errors)
    for dataset, value in by_dataset.items():
        nested = value.get("by_linkage_basis_and_threshold") if isinstance(value, dict) else None
        _require(isinstance(nested, dict) and bool(nested), f"M2 ranking by_dataset.{dataset} 缺分層結果", errors)
        if isinstance(nested, dict):
            _require(all(basis in nested for basis in ("PS_HP1", "PS_HP2")), f"M2 ranking by_dataset.{dataset} 缺 PS_HP1/PS_HP2", errors)
            for basis in ("PS_HP1", "PS_HP2"):
                thresholds = nested.get(basis)
                _require(isinstance(thresholds, dict) and set(thresholds) == {"1", "2", "3", "5"}, f"M2 ranking by_dataset.{dataset}.{basis} threshold grid 不完整", errors)
                if isinstance(thresholds, dict):
                    for threshold, cell in thresholds.items():
                        _require(isinstance(cell, dict) and bool(cell), f"M2 ranking by_dataset.{dataset}.{basis}/{threshold} cell 為空", errors)
                        if isinstance(cell, dict) and cell:
                            _validate_ranking_cell(cell, f"by_dataset.{dataset}.{basis}/{threshold}", errors)
        funnel = value.get("sample_funnel") if isinstance(value, dict) else None
        _require(isinstance(funnel, dict), f"M2 ranking by_dataset.{dataset} 缺 sample_funnel", errors)
        if isinstance(funnel, dict):
            required_counts = (
                "raw_sparse_molecules",
                "ps_known_molecules",
                "ps_missing_molecules",
                "hp_included_molecules",
                "hp_excluded_molecules",
            )
            _require(all(_is_nonnegative_int(funnel.get(key)) for key in required_counts), f"M2 ranking by_dataset.{dataset} molecule funnel 型別錯誤", errors)
            if all(_is_nonnegative_int(funnel.get(key)) for key in required_counts):
                _require(funnel["raw_sparse_molecules"] == funnel["ps_known_molecules"] + funnel["ps_missing_molecules"], f"M2 ranking by_dataset.{dataset} PS funnel 不守恆", errors)
                _require(funnel["raw_sparse_molecules"] == funnel["hp_included_molecules"] + funnel["hp_excluded_molecules"], f"M2 ranking by_dataset.{dataset} HP funnel 不守恆", errors)
            call_classes = funnel.get("call_class_counts") or {}
            _require(set(call_classes) == {"R", "A", "O", "D", "S", "L", "X"}, f"M2 ranking by_dataset.{dataset} R/A/O/D/S/L/X 不完整", errors)
            _require(all(_is_nonnegative_int(count) for count in call_classes.values()), f"M2 ranking by_dataset.{dataset} call-class count 型別錯誤", errors)
            _require(funnel.get("biological_id") == BIOLOGICAL_IDS.get(dataset), f"M2 ranking by_dataset.{dataset} biological_id 不符", errors)
    if set(by_dataset) == set(DATASETS) and isinstance(cells, dict):
        for basis in ("PS_HP1", "PS_HP2"):
            for threshold in ("1", "2", "3", "5"):
                aggregate_cell = (cells.get(basis) or {}).get(threshold)
                dataset_cells = [
                    (((by_dataset.get(dataset) or {}).get("by_linkage_basis_and_threshold") or {}).get(basis) or {}).get(threshold)
                    for dataset in DATASETS
                ]
                if isinstance(aggregate_cell, dict) and all(isinstance(cell, dict) and bool(cell) for cell in dataset_cells):
                    try:
                        recomputed = _sum_count_trees(dataset_cells)
                    except ReportGateError as exc:
                        errors.append(f"M2 ranking per-dataset aggregate 無法重算：{exc}")
                    else:
                        _require(recomputed == aggregate_cell, f"M2 ranking per-dataset {basis}/{threshold} 加總不符 aggregate", errors)
    primary_unit_evaluations = 0
    if isinstance(cells, dict):
        for thresholds in cells.values():
            if isinstance(thresholds, dict):
                for cell in thresholds.values():
                    if isinstance(cell, dict) and _is_nonnegative_int(cell.get("n_component_hp_units")):
                        primary_unit_evaluations += cell["n_component_hp_units"]
    _validate_runtime_diagnostics(
        payload.get("runtime_diagnostics"),
        primary_unit_evaluations,
        "M2 ranking full",
        errors,
    )
    return errors


def _validate_m2_pilot(payload: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.lossless_read_linkage_chromosome_receipt", "M2 pilot extraction schema_name 不符", errors)
    _require(payload.get("schema_version") == "1.0.0", "Legacy M2 pilot schema_version 不是 1.0.0", errors)
    _require(payload.get("all_pass") is True, "M2 pilot extraction all_pass 不是 true", errors)
    scope = payload.get("scope") or {}
    _require(scope.get("dataset") == "HCC1954" and scope.get("chrom") == "chr22", "M2 pilot extraction 不是 HCC1954 chr22", errors)
    counts = payload.get("counts") or {}
    raw_hp = counts.get("raw_HP_counts")
    _require(isinstance(raw_hp, dict) and all(_is_nonnegative_int(value) for value in raw_hp.values()), "M2 pilot raw_HP_counts 缺失", errors)
    if isinstance(raw_hp, dict) and all(_is_nonnegative_int(value) for value in raw_hp.values()):
        _require("." in raw_hp, "M2 pilot raw_HP_counts 缺 missing HP='.' 類別", errors)
        _require(
            _is_nonnegative_int(counts.get("canonical_eligible_alignments")),
            "M2 pilot canonical_eligible_alignments 缺失",
            errors,
        )
        if _is_nonnegative_int(counts.get("canonical_eligible_alignments")):
            _require(
                sum(raw_hp.values()) == counts["canonical_eligible_alignments"],
                "M2 pilot raw HP values（含 '.'）總和不等於 canonical eligible alignments",
                errors,
            )
    return errors


def _validate_m2_verification(
    payload: Mapping[str, Any],
    extraction_source: Source,
    ranking_source: Source,
    candidate_source: Source,
) -> list[str]:
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.m2_full_independent_verification", "M2 verifier schema_name 不符", errors)
    _require(payload.get("schema_version") == "1.0.0", "M2 verifier schema_version 不是 1.0.0", errors)
    _require(payload.get("task_type") == "B_COMPREHENSIVE_VALIDATION", "M2 verifier 不是 Task-B", errors)
    _require(payload.get("all_pass") is True, "M2 verifier all_pass 不是 true", errors)
    release_binding = payload.get("release_binding") or {}
    _require(
        release_binding.get("schema_name") == "intersubmod.m2_release_binding"
        and release_binding.get("schema_version") == "1.0.0"
        and release_binding.get("authority_mode") == "CANONICAL_V5_FROZEN"
        and release_binding.get("validation_evidence_eligible") is True,
        "M2 verifier 未綁定可用於 validation 的 canonical release contract",
        errors,
    )
    deep_release = release_binding.get("deep_release_verification") or {}
    _require(
        deep_release.get("all_pass") is True
        and deep_release.get("verified_snapshot_files") == 11,
        "M2 verifier 未由 frozen freezer 深度重驗 11 個 snapshot files",
        errors,
    )
    post_identity = payload.get("post_input_identity") or {}
    _require(
        post_identity.get("exact_snapshot_equal") is True
        and post_identity.get("n_artifact_roles") == 42
        and _is_sha256(post_identity.get("input_identity_snapshot_sha256")),
        "M2 verifier 未證明 POST 與 frozen PRE 的 42-role snapshot 完全一致",
        errors,
    )
    for boundary_name in ("post_receipt", "immutable_pre_receipt"):
        boundary = post_identity.get(boundary_name) or {}
        _require(
            _is_sha256(boundary.get("sha256"))
            and _is_sha256(boundary.get("semantic_sha256")),
            f"M2 verifier {boundary_name} 身分不完整",
            errors,
        )
    _require(bool(post_identity.get("claim_boundary")), "M2 verifier 缺 PRE/POST claim boundary", errors)
    scope = payload.get("scope") or {}
    _require(tuple(scope.get("datasets") or ()) == DATASETS, "M2 verifier datasets/order 不符", errors)
    _require(tuple(scope.get("chromosomes") or ()) == AUTOSOMES, "M2 verifier chromosomes/order 不符", errors)
    _require(scope.get("expected_tasks") == 154, "M2 verifier expected_tasks 不是 154", errors)
    independence = payload.get("verification_independence") or {}
    _require(independence.get("imports_production_aggregator") is False, "M2 verifier 匯入 production aggregator", errors)
    _require(independence.get("imports_production_ranker") is False, "M2 verifier 匯入 production ranker", errors)
    _require(independence.get("reads_bam") is False, "M2 verifier 不應重讀 BAM", errors)
    _require(independence.get("recomputes_from_child_receipts_and_tables") is True, "M2 verifier 未由 child receipts/tables 重算", errors)
    _require(independence.get("candidate_table_compared_rowwise_to_independent_child_reconstruction") is True, "M2 verifier 未逐列重建 candidate table", errors)
    _require(
        independence.get("method_contract_verification_mode")
        == "STATIC_EXACT_CONTRACT_AND_ACTUAL_SOURCE_SHA_BINDING",
        "M2 verifier 未明示 method-contract static/source-SHA 驗證模式",
        errors,
    )
    _require(
        independence.get("profile_likelihood_recomputed_from_patterns_states_pi") is True,
        "M2 verifier 未獨立由 patterns/states/pi 重算 profile likelihood",
        errors,
    )
    _require(
        bool(independence.get("profile_likelihood_claim_boundary")),
        "M2 verifier 缺 profile-likelihood claim boundary",
        errors,
    )
    extraction = payload.get("extraction") or {}
    ranking = payload.get("ranking") or {}
    _require(extraction.get("receipt_sha256") == extraction_source.sha256, "M2 verifier extraction receipt SHA-256 不符", errors)
    _require(ranking.get("receipt_sha256") == ranking_source.sha256, "M2 verifier ranking receipt SHA-256 不符", errors)
    _require(extraction.get("n_tasks") == 154, "M2 verifier extraction n_tasks 不是 154", errors)
    _require(ranking.get("n_tasks") == 154, "M2 verifier ranking n_tasks 不是 154", errors)
    _require(ranking.get("all_aggregate_cells_conserved") is True, "M2 verifier ranking aggregate cells 未守恆", errors)
    expected_chain = [8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 154]
    for stage_name, stage_payload in (("extraction", extraction), ("ranking", ranking)):
        orchestration = stage_payload.get("orchestration") or {}
        _require(orchestration.get("stage") == stage_name, f"M2 verifier {stage_name} orchestration stage 不符", errors)
        _require(
            orchestration.get("legal_cumulative_counts") == expected_chain
            and orchestration.get("n_batches") == len(expected_chain)
            and orchestration.get("n_attested_children") == 154,
            f"M2 verifier {stage_name} orchestration 8→16 batch/154-child chain 不完整",
            errors,
        )
        for flag in (
            "all_child_receipts_and_outputs_rehashed",
            "all_grants_and_completions_session_bound",
            "all_checkpoints_and_terminal_chain_bound",
            "no_open_orphan_or_preseeded_child_accepted",
        ):
            _require(
                orchestration.get(flag) is True,
                f"M2 verifier {stage_name} orchestration.{flag} 不是 true",
                errors,
            )
        _require(
            _is_sha256(orchestration.get("session_id"))
            and _is_sha256((orchestration.get("session_receipt") or {}).get("sha256"))
            and _is_sha256((orchestration.get("terminal_receipt") or {}).get("sha256"))
            and bool(orchestration.get("claim_boundary")),
            f"M2 verifier {stage_name} orchestration identity/claim boundary 不完整",
            errors,
        )
    method_verification = ranking.get("method_contract_verification") or {}
    _require(
        method_verification.get("contract") == EXPECTED_METHOD_CONTRACT,
        "M2 verifier method contract exact comparison 不符",
        errors,
    )
    _require(
        method_verification.get("n_children_exactly_matched_and_source_bound") == 154
        and method_verification.get("all_children_exactly_matched_and_source_bound") is True,
        "M2 verifier 未驗證 154 child contracts/source binding",
        errors,
    )
    candidate = ranking.get("candidate_table") or {}
    _require(candidate.get("sha256") == candidate_source.sha256, "M2 verifier candidate table SHA-256 不符", errors)
    _require(candidate.get("all_rows_match_independent_child_reconstruction") is True, "M2 verifier candidate rows 未逐列吻合", errors)
    _require(candidate.get("winner_partitions_conserved") is True, "M2 verifier candidate winner partition 未守恆", errors)
    profile = ranking.get("profile_likelihood_independent_recomputation") or {}
    _require(profile.get("n_children") == 154, "M2 verifier profile likelihood 未覆蓋 154 child receipts", errors)
    _require(
        profile.get("n_units") == candidate.get("n_units")
        and profile.get("n_candidates") == candidate.get("n_rows"),
        "M2 verifier profile likelihood unit/candidate coverage 不等於 canonical candidate table",
        errors,
    )
    for key in (
        "n_units", "n_candidates", "n_pattern_rows", "n_scoring_molecules",
        "peak_pattern_rows_per_unit", "peak_candidates_per_unit", "peak_states_per_candidate",
        "peak_emission_cells_per_candidate",
    ):
        _require(_is_nonnegative_int(profile.get(key)), f"M2 verifier profile likelihood {key} 型別/範圍錯誤", errors)
    for key in (
        "max_abs_ll_delta", "max_abs_relative_weight_delta", "max_abs_gap_delta",
        "max_abs_simplex_residual_delta",
    ):
        _require(
            _is_finite_nonnegative_number(profile.get(key)),
            f"M2 verifier profile likelihood {key} 非有限非負數",
            errors,
        )
    for key in (
        "all_profile_likelihoods_and_certificates_match",
        "all_relative_weights_match",
        "all_winner_tie_partitions_match",
    ):
        _require(profile.get(key) is True, f"M2 verifier profile likelihood {key} 不是 true", errors)
    _require(
        isinstance(profile.get("child_summaries"), list)
        and len(profile.get("child_summaries") or ()) == 154,
        "M2 verifier profile likelihood child summaries 不是 154 筆",
        errors,
    )
    _require(
        isinstance(profile.get("count_unit_definitions"), dict)
        and isinstance(profile.get("numeric_contract"), dict)
        and bool(profile.get("streaming_memory_contract")),
        "M2 verifier profile likelihood 缺 count/numeric/streaming contract",
        errors,
    )
    verifier_runtime = ranking.get("runtime_diagnostics") or {}
    ranking_runtime = (ranking_source.payload or {}).get("runtime_diagnostics") or {}
    expected_runtime_keys = {
        "n_child_runtime_files", "n_unit_evaluations", "metrics",
        "metrics_when_invoked",
        "all_child_and_full_runtime_summaries_independently_recomputed",
    }
    _require(
        isinstance(verifier_runtime, dict) and set(verifier_runtime) == expected_runtime_keys,
        "M2 verifier runtime_diagnostics keys 不完整",
        errors,
    )
    if isinstance(verifier_runtime, dict) and set(verifier_runtime) == expected_runtime_keys:
        _require(verifier_runtime["n_child_runtime_files"] == 154, "M2 verifier runtime 未覆蓋 154 children", errors)
        _require(
            verifier_runtime["n_unit_evaluations"] == ranking_runtime.get("n_unit_evaluations"),
            "M2 verifier runtime unit count 與 producer 不符",
            errors,
        )
        _require(verifier_runtime["metrics"] == ranking_runtime.get("metrics"), "M2 verifier runtime metrics 與獨立重算不符", errors)
        _require(verifier_runtime["metrics_when_invoked"] == ranking_runtime.get("metrics_when_invoked"), "M2 verifier invoked runtime metrics 與獨立重算不符", errors)
        _require(verifier_runtime["all_child_and_full_runtime_summaries_independently_recomputed"] is True, "M2 verifier runtime 未獨立重算", errors)
    checks = payload.get("checks") or {}
    for key in (
        "runtime_diagnostics_independently_recomputed",
        "profile_likelihood_recomputed_independently",
        "profile_likelihood_certificates_match",
        "profile_relative_weights_match",
        "profile_winner_tie_partitions_match",
        "release_contract_authenticated_and_eligible",
        "release_contract_all_snapshot_sources_rehashed",
        "extraction_and_ranking_bound_to_same_release",
        "full_runner_dependency_paths_and_shas_match_release",
        "frozen_scientific_and_scheduler_parameters_match",
        "post_input_identity_authenticated_and_exactly_equals_frozen_pre",
        "extraction_orchestration_session_batch_grant_completion_chain_verified",
        "ranking_orchestration_session_batch_grant_completion_chain_verified",
        "ranking_orchestration_bound_to_verified_extraction_session",
    ):
        _require(checks.get(key) is True, f"M2 verifier checks.{key} 不是 true", errors)
    _require(bool(checks) and all(value is True for value in checks.values()), "M2 verifier checks 未全數 true", errors)
    return errors


def _validate_numeric_ratio(
    payload: Any,
    numerator: int | float,
    denominator: int | float,
    label: str,
    errors: list[str],
) -> None:
    _require(
        isinstance(payload, dict) and set(payload) == NUMERIC_SUMMARY_RATIO_KEYS,
        f"M2 numeric summary {label} ratio schema 不完整",
        errors,
    )
    if not isinstance(payload, dict) or set(payload) != NUMERIC_SUMMARY_RATIO_KEYS:
        return
    _require(
        payload.get("numerator") == numerator
        and payload.get("denominator") == denominator,
        f"M2 numeric summary {label} ratio numerator/denominator 不符",
        errors,
    )
    _require(
        isinstance(payload.get("denominator_label"), str)
        and bool(payload.get("denominator_label")),
        f"M2 numeric summary {label} ratio denominator_label 缺失",
        errors,
    )
    if denominator == 0:
        _require(
            payload.get("value") is None
            and payload.get("percent") is None
            and payload.get("reason") == "denominator_is_zero",
            f"M2 numeric summary {label} zero-denominator ratio 未明示 NA",
            errors,
        )
        return
    expected_value = float(numerator) / float(denominator)
    expected_percent = expected_value * 100.0
    _require(
        _is_finite_nonnegative_number(payload.get("value"))
        and math.isclose(float(payload["value"]), expected_value, rel_tol=1e-12, abs_tol=1e-12),
        f"M2 numeric summary {label} ratio value 不可重算",
        errors,
    )
    _require(
        _is_finite_nonnegative_number(payload.get("percent"))
        and math.isclose(float(payload["percent"]), expected_percent, rel_tol=1e-12, abs_tol=1e-10),
        f"M2 numeric summary {label} ratio percent 不可重算",
        errors,
    )
    _require(payload.get("reason") is None, f"M2 numeric summary {label} 非零分母卻有 NA reason", errors)


def _validate_numeric_component_cell(
    cell: Any, label: str, errors: list[str]
) -> None:
    required = {
        "n_components", "k_equals_1", "k_greater_than_1", "k_distribution",
        "max_k_component_sites", "active_site_membership_mass",
        "k_equals_1_share", "k_greater_than_1_share",
    }
    _require(isinstance(cell, dict) and set(cell) == required, f"M2 numeric summary {label} component keys 不完整", errors)
    if not isinstance(cell, dict) or set(cell) != required:
        return
    count_keys = (
        "n_components", "k_equals_1", "k_greater_than_1",
        "max_k_component_sites", "active_site_membership_mass",
    )
    _require(all(_is_nonnegative_int(cell.get(key)) for key in count_keys), f"M2 numeric summary {label} component count 無效", errors)
    distribution = cell.get("k_distribution")
    distribution_valid = (
        _count_mapping(distribution)
        and all(key.isdigit() and int(key) >= 1 for key in distribution)
    )
    _require(distribution_valid, f"M2 numeric summary {label} k distribution 無效", errors)
    if not all(_is_nonnegative_int(cell.get(key)) for key in count_keys) or not distribution_valid:
        return
    n_components = cell["n_components"]
    _require(sum(distribution.values()) == n_components, f"M2 numeric summary {label} k distribution 不守恆", errors)
    _require(cell["k_equals_1"] == distribution.get("1", 0), f"M2 numeric summary {label} k=1 不符 distribution", errors)
    _require(
        cell["k_greater_than_1"] == sum(value for key, value in distribution.items() if int(key) > 1),
        f"M2 numeric summary {label} k>1 不符 distribution",
        errors,
    )
    _require(cell["k_equals_1"] + cell["k_greater_than_1"] == n_components, f"M2 numeric summary {label} k=1+k>1 != components", errors)
    _require(cell["max_k_component_sites"] == max((int(key) for key in distribution), default=0), f"M2 numeric summary {label} max k 不符", errors)
    _require(
        cell["active_site_membership_mass"]
        == sum(int(key) * value for key, value in distribution.items()),
        f"M2 numeric summary {label} active-site mass 不符",
        errors,
    )
    _validate_numeric_ratio(cell["k_equals_1_share"], cell["k_equals_1"], n_components, f"{label}/k=1", errors)
    _validate_numeric_ratio(cell["k_greater_than_1_share"], cell["k_greater_than_1"], n_components, f"{label}/k>1", errors)


def _validate_numeric_rank_cell(cell: Any, label: str, errors: list[str]) -> None:
    required = {
        "units", "molecule_funnel", "partial_read_funnel", "candidate_structure",
        "ranking_outcome", "topology", "effective_k",
    }
    _require(isinstance(cell, dict) and set(cell) == required, f"M2 numeric summary {label} rank cell keys 不完整", errors)
    if not isinstance(cell, dict) or set(cell) != required:
        return
    units = cell.get("units") or {}
    outcome = cell.get("ranking_outcome") or {}
    effective = cell.get("effective_k") or {}
    structure = cell.get("candidate_structure") or {}
    candidate = structure.get("candidate_table") or {}
    exact_nested_keys = {
        "units": {
            "n_component_hp_unit_evaluations", "solver_complete",
            "solver_incomplete_or_not_run", "solver_complete_share",
        },
        "molecule_funnel": {
            "component_projections", "informative_scoring_molecules",
            "all_X_excluded_molecules", "structural_retained_molecules",
            "below_structural_minread_but_scored_molecules",
            "informative_share_of_projections", "structural_share_of_informative",
            "below_minread_share_of_informative",
        },
        "partial_read_funnel": {
            "structural_partial_pattern_groups", "coverage_denominator", "covered",
            "unsatisfied", "covered_share", "full_detail",
        },
        "candidate_structure": {
            "raw_parent_edge_trees_T_complete_units",
            "distinct_optimal_vertex_sets_V_complete_units",
            "mean_T_per_solver_complete_unit", "mean_V_per_solver_complete_unit",
            "candidate_table",
        },
        "ranking_outcome": {
            "unique_first", "tied_first", "abstain_all_causes", "unique_share",
            "tied_share", "abstain_share", "selection_status_counts",
        },
        "topology": {
            "evaluated_units", "coarse_class_unique_units",
            "coarse_class_multiple_units", "coarse_unique_class_counts",
            "coarse_ambiguous_class_set_counts",
            "parent_edge_assignment_unique_units",
            "exact_topology_proven_unique_units", "coarse_class_unique_share",
            "coarse_class_multiple_share", "exact_topology_proven_unique_share",
        },
        "effective_k": {
            "component_site_mass", "observed_ALT_active_mass",
            "not_structural_ALT_active_mass", "k_route_counts",
            "observed_ALT_active_share_of_component_site_mass",
            "not_structural_ALT_active_share_of_component_site_mass",
            "route_shares_of_unit_evaluations",
        },
    }
    for group, required_keys in exact_nested_keys.items():
        value = cell.get(group)
        _require(
            isinstance(value, dict) and set(value) == required_keys,
            f"M2 numeric summary {label}/{group} keys 不完整",
            errors,
        )
    unit_count = units.get("n_component_hp_unit_evaluations")
    complete = units.get("solver_complete")
    incomplete = units.get("solver_incomplete_or_not_run")
    unique = outcome.get("unique_first")
    tied = outcome.get("tied_first")
    abstain = outcome.get("abstain_all_causes")
    scalar_values = (unit_count, complete, incomplete, unique, tied, abstain)
    _require(all(_is_nonnegative_int(value) for value in scalar_values), f"M2 numeric summary {label} unit/outcome count 無效", errors)
    if not all(_is_nonnegative_int(value) for value in scalar_values):
        return
    _require(complete + incomplete == unit_count, f"M2 numeric summary {label} solver partition 不守恆", errors)
    _require(unique + tied + abstain == unit_count, f"M2 numeric summary {label} ranking partition 不守恆", errors)
    _validate_numeric_ratio(units.get("solver_complete_share"), complete, unit_count, f"{label}/solver complete", errors)
    for key, value in (("unique_share", unique), ("tied_share", tied), ("abstain_share", abstain)):
        _validate_numeric_ratio(outcome.get(key), value, unit_count, f"{label}/{key}", errors)

    molecule = cell.get("molecule_funnel") or {}
    molecule_counts = (
        molecule.get("component_projections"),
        molecule.get("informative_scoring_molecules"),
        molecule.get("all_X_excluded_molecules"),
        molecule.get("structural_retained_molecules"),
        molecule.get("below_structural_minread_but_scored_molecules"),
    )
    _require(all(_is_nonnegative_int(value) for value in molecule_counts), f"M2 numeric summary {label} molecule funnel count 無效", errors)
    if all(_is_nonnegative_int(value) for value in molecule_counts):
        projections, informative, all_x, structural, below = molecule_counts
        _require(projections == informative + all_x, f"M2 numeric summary {label} projection funnel 不守恆", errors)
        _require(informative == structural + below, f"M2 numeric summary {label} structural/scoring funnel 不守恆", errors)
        _validate_numeric_ratio(molecule.get("informative_share_of_projections"), informative, projections, f"{label}/informative projection", errors)
        _validate_numeric_ratio(molecule.get("structural_share_of_informative"), structural, informative, f"{label}/structural informative", errors)
        _validate_numeric_ratio(molecule.get("below_minread_share_of_informative"), below, informative, f"{label}/below-minread informative", errors)

    partial = cell.get("partial_read_funnel") or {}
    partial_counts = (
        partial.get("structural_partial_pattern_groups"),
        partial.get("coverage_denominator"), partial.get("covered"),
        partial.get("unsatisfied"),
    )
    _require(all(_is_nonnegative_int(value) for value in partial_counts), f"M2 numeric summary {label} partial funnel count 無效", errors)
    if all(_is_nonnegative_int(value) for value in partial_counts):
        _, partial_denominator, covered, unsatisfied = partial_counts
        _require(covered + unsatisfied == partial_denominator, f"M2 numeric summary {label} partial coverage 不守恆", errors)
        _validate_numeric_ratio(partial.get("covered_share"), covered, partial_denominator, f"{label}/partial covered", errors)
    _require(isinstance(partial.get("full_detail"), dict), f"M2 numeric summary {label} partial full_detail 缺失", errors)

    mass_keys = (
        "component_site_mass", "observed_ALT_active_mass",
        "not_structural_ALT_active_mass",
    )
    _require(all(_is_nonnegative_int(effective.get(key)) for key in mass_keys), f"M2 numeric summary {label} effective-k mass 無效", errors)
    routes = effective.get("k_route_counts") or {}
    expected_routes = {"EXACT_COMPLETE", "EXACT_INCOMPLETE", "NOT_RUN_NO_STRUCTURE", "GT_EXACT_LIMIT"}
    _require(_count_mapping(routes) and set(routes) == expected_routes, f"M2 numeric summary {label} effective-k routes 不完整", errors)
    if all(_is_nonnegative_int(effective.get(key)) for key in mass_keys):
        _require(
            effective["component_site_mass"]
            == effective["observed_ALT_active_mass"] + effective["not_structural_ALT_active_mass"],
            f"M2 numeric summary {label} effective-k mass 不守恆",
            errors,
        )
        _validate_numeric_ratio(
            effective.get("observed_ALT_active_share_of_component_site_mass"),
            effective["observed_ALT_active_mass"], effective["component_site_mass"],
            f"{label}/effective-k active mass", errors,
        )
        _validate_numeric_ratio(
            effective.get("not_structural_ALT_active_share_of_component_site_mass"),
            effective["not_structural_ALT_active_mass"], effective["component_site_mass"],
            f"{label}/effective-k nonstructural mass", errors,
        )
    if _count_mapping(routes):
        _require(sum(routes.values()) == unit_count, f"M2 numeric summary {label} effective-k route partition 不守恆", errors)
        _require(routes.get("EXACT_COMPLETE") == complete, f"M2 numeric summary {label} EXACT_COMPLETE != solver complete", errors)
        _require(
            sum(count for route, count in routes.items() if route != "EXACT_COMPLETE") == incomplete,
            f"M2 numeric summary {label} non-complete routes != solver incomplete",
            errors,
        )
        route_ratios = effective.get("route_shares_of_unit_evaluations") or {}
        _require(isinstance(route_ratios, dict) and set(route_ratios) == set(routes), f"M2 numeric summary {label} route ratios 不完整", errors)
        if isinstance(route_ratios, dict):
            for route, count in routes.items():
                _validate_numeric_ratio(route_ratios.get(route), count, unit_count, f"{label}/route/{route}", errors)

    outer_t = structure.get("raw_parent_edge_trees_T_complete_units")
    outer_v = structure.get("distinct_optimal_vertex_sets_V_complete_units")
    candidate_count_keys = (
        "n_solver_complete_candidate_units", "n_candidate_vertex_sets_V",
        "n_parent_edge_trees_T", "unique_first", "tied_first",
        "solver_complete_optimizer_abstain",
    )
    _require(
        _is_nonnegative_int(outer_t)
        and _is_nonnegative_int(outer_v)
        and all(_is_nonnegative_int(candidate.get(key)) for key in candidate_count_keys),
        f"M2 numeric summary {label} candidate count 無效",
        errors,
    )
    if not (
        _is_nonnegative_int(outer_t)
        and _is_nonnegative_int(outer_v)
        and all(_is_nonnegative_int(candidate.get(key)) for key in candidate_count_keys)
    ):
        return
    optimizer_abstain = candidate["solver_complete_optimizer_abstain"]
    _require(outer_t >= outer_v, f"M2 numeric summary {label} T < V", errors)
    _require(candidate["n_solver_complete_candidate_units"] == complete, f"M2 numeric summary {label} candidate units != solver complete", errors)
    _require(candidate["n_candidate_vertex_sets_V"] == outer_v and candidate["n_parent_edge_trees_T"] == outer_t, f"M2 numeric summary {label} candidate T/V 不符 rank cell", errors)
    _require(candidate["unique_first"] == unique and candidate["tied_first"] == tied, f"M2 numeric summary {label} candidate unique/tied 不符 rank outcome", errors)
    _require(optimizer_abstain == abstain - incomplete, f"M2 numeric summary {label} optimizer abstain 不符", errors)
    _require(unique + tied + optimizer_abstain == complete, f"M2 numeric summary {label} candidate outcome partition 不守恆", errors)
    for key, count in (("unique_first_share", unique), ("tied_first_share", tied), ("optimizer_abstain_share", optimizer_abstain)):
        _validate_numeric_ratio(candidate.get(key), count, complete, f"{label}/candidate/{key}", errors)
    _validate_numeric_ratio(structure.get("mean_T_per_solver_complete_unit"), outer_t, complete, f"{label}/mean T", errors)
    _validate_numeric_ratio(structure.get("mean_V_per_solver_complete_unit"), outer_v, complete, f"{label}/mean V", errors)

    required_candidate_keys = {
        "n_solver_complete_candidate_units", "n_candidate_vertex_sets_V",
        "n_parent_edge_trees_T", "tree_vertex_partition", "unique_first",
        "tied_first", "solver_complete_optimizer_abstain", "unique_first_share",
        "tied_first_share", "optimizer_abstain_share", "tied_by_coarse_topology",
        "topology", "h_star_distribution", "coarse_topology_class_inclusion_counts",
    }
    _require(set(candidate) == required_candidate_keys, f"M2 numeric summary {label} candidate table keys 不完整", errors)
    tree_partition = candidate.get("tree_vertex_partition") or {}
    required_partition_keys = {
        "counts", "denominator", "shares", "definition", "bucket_definitions",
    }
    _require(set(tree_partition) == required_partition_keys, f"M2 numeric summary {label} T/V partition keys 不完整", errors)
    partition_counts = tree_partition.get("counts") or {}
    partition_shares = tree_partition.get("shares") or {}
    _require(_count_mapping(partition_counts) and set(partition_counts) == set(TREE_VERTEX_BUCKETS), f"M2 numeric summary {label} T/V partition counts 不完整", errors)
    _require(isinstance(partition_shares, dict) and set(partition_shares) == set(TREE_VERTEX_BUCKETS), f"M2 numeric summary {label} T/V partition shares 不完整", errors)
    _require(tree_partition.get("denominator") == complete, f"M2 numeric summary {label} T/V partition denominator 不符", errors)
    if _count_mapping(partition_counts) and set(partition_counts) == set(TREE_VERTEX_BUCKETS):
        _require(sum(partition_counts.values()) == complete, f"M2 numeric summary {label} T/V partition 不守恆", errors)
        for bucket in TREE_VERTEX_BUCKETS:
            _validate_numeric_ratio(partition_shares.get(bucket), partition_counts[bucket], complete, f"{label}/T-V/{bucket}", errors)
    _require(
        tree_partition.get("definition") == TREE_VERTEX_PARTITION_DEFINITION
        and tree_partition.get("bucket_definitions") == TREE_VERTEX_BUCKET_DEFINITIONS,
        f"M2 numeric summary {label} T/V partition definition 漂移",
        errors,
    )
    h_distribution = candidate.get("h_star_distribution") or {}
    _require(_count_mapping(h_distribution) and all(key.isdigit() for key in h_distribution), f"M2 numeric summary {label} h* distribution 無效", errors)
    if _count_mapping(h_distribution):
        _require(sum(h_distribution.values()) == complete, f"M2 numeric summary {label} h* distribution 不守恆", errors)
    tied_topology = candidate.get("tied_by_coarse_topology") or {}
    consistent = tied_topology.get("consistent")
    inconsistent = tied_topology.get("inconsistent")
    tied_denominator = tied_topology.get("denominator")
    _require(
        _is_nonnegative_int(consistent)
        and _is_nonnegative_int(inconsistent)
        and tied_denominator == tied
        and consistent + inconsistent == tied,
        f"M2 numeric summary {label} tied×coarse-Topo partition 不守恆",
        errors,
    )
    if _is_nonnegative_int(consistent) and _is_nonnegative_int(inconsistent) and tied_denominator == tied:
        _validate_numeric_ratio(tied_topology.get("consistent_share"), consistent, tied, f"{label}/tied Topo consistent", errors)
        _validate_numeric_ratio(tied_topology.get("inconsistent_share"), inconsistent, tied, f"{label}/tied Topo inconsistent", errors)
    _require(isinstance(tied_topology.get("definition"), str) and bool(tied_topology.get("definition")), f"M2 numeric summary {label} tied×Topo definition 缺失", errors)

    topology = cell.get("topology") or {}
    topology_candidate = candidate.get("topology") or {}
    required_candidate_topology_keys = {
        "evaluated_units", "coarse_class_unique_units", "coarse_class_multiple_units",
        "coarse_unique_class_counts", "coarse_ambiguous_class_set_counts",
        "parent_edge_assignment_unique_units", "exact_topology_proven_unique_units",
    }
    _require(set(topology_candidate) == required_candidate_topology_keys, f"M2 numeric summary {label} candidate topology keys 不完整", errors)
    evaluated = topology.get("evaluated_units")
    coarse_unique = topology.get("coarse_class_unique_units")
    coarse_multiple = topology.get("coarse_class_multiple_units")
    parent_unique = topology.get("parent_edge_assignment_unique_units")
    exact_unique = topology.get("exact_topology_proven_unique_units")
    _require(all(_is_nonnegative_int(value) for value in (evaluated, coarse_unique, coarse_multiple, parent_unique, exact_unique)), f"M2 numeric summary {label} topology count 無效", errors)
    unique_classes = topology.get("coarse_unique_class_counts") or {}
    ambiguous_sets = topology.get("coarse_ambiguous_class_set_counts") or {}
    _require(_count_mapping(unique_classes), f"M2 numeric summary {label} coarse unique class counts 無效", errors)
    _require(_count_mapping(ambiguous_sets), f"M2 numeric summary {label} ambiguous class-set counts 無效", errors)
    if all(_is_nonnegative_int(value) for value in (evaluated, coarse_unique, coarse_multiple, parent_unique, exact_unique)):
        _require(coarse_unique + coarse_multiple == evaluated, f"M2 numeric summary {label} topology partition 不守恆", errors)
        _require(parent_unique == exact_unique <= evaluated, f"M2 numeric summary {label} parent/exact uniqueness 不守恆", errors)
        if _count_mapping(unique_classes):
            _require(sum(unique_classes.values()) == coarse_unique, f"M2 numeric summary {label} unique topology classes 不守恆", errors)
        if _count_mapping(ambiguous_sets):
            _require(sum(ambiguous_sets.values()) == coarse_multiple, f"M2 numeric summary {label} ambiguous topology sets 不守恆", errors)
        _validate_numeric_ratio(topology.get("coarse_class_unique_share"), coarse_unique, evaluated, f"{label}/Topo unique", errors)
        _validate_numeric_ratio(topology.get("coarse_class_multiple_share"), coarse_multiple, evaluated, f"{label}/Topo multiple", errors)
        _validate_numeric_ratio(topology.get("exact_topology_proven_unique_share"), exact_unique, evaluated, f"{label}/exact Topo unique", errors)
    for key in required_candidate_topology_keys:
        _require(topology_candidate.get(key) == topology.get(key), f"M2 numeric summary {label} candidate/topology {key} 不符", errors)
    inclusion = candidate.get("coarse_topology_class_inclusion_counts") or {}
    _require(_count_mapping(inclusion), f"M2 numeric summary {label} topology inclusion counts 無效", errors)


def _validate_component_against_raw(
    cell: Mapping[str, Any], raw: Mapping[str, Any], label: str, errors: list[str]
) -> None:
    expected = {
        "n_components": raw.get("n_components"),
        "k_equals_1": raw.get("n_singletons_k1"),
        "k_greater_than_1": raw.get("n_multisite_k_gt1"),
        "k_distribution": raw.get("k_component_sites_distribution") or {},
        "max_k_component_sites": raw.get("max_k_component_sites"),
        "active_site_membership_mass": raw.get("n_active_site_memberships"),
    }
    for key, value in expected.items():
        _require(cell.get(key) == value, f"M2 numeric summary {label}/{key} 未精確 reconcile extraction raw cell", errors)


def _validate_rank_against_raw(
    cell: Mapping[str, Any], raw: Mapping[str, Any], label: str, errors: list[str]
) -> None:
    scalar_bindings = (
        ("units", "n_component_hp_unit_evaluations", "n_component_hp_units"),
        ("units", "solver_complete", "solver_complete_units"),
        ("units", "solver_incomplete_or_not_run", "solver_incomplete_or_not_run_units"),
        ("molecule_funnel", "component_projections", "molecule_component_projections"),
        ("molecule_funnel", "informative_scoring_molecules", "informative_scoring_molecules"),
        ("molecule_funnel", "all_X_excluded_molecules", "all_x_excluded_molecules"),
        ("molecule_funnel", "structural_retained_molecules", "structural_retained_molecules"),
        ("molecule_funnel", "below_structural_minread_but_scored_molecules", "below_minread_scoring_molecules"),
        ("partial_read_funnel", "structural_partial_pattern_groups", "structural_partial_pattern_groups"),
        ("partial_read_funnel", "coverage_denominator", "partial_group_coverage_denominator"),
        ("partial_read_funnel", "covered", "partial_groups_covered"),
        ("partial_read_funnel", "unsatisfied", "partial_groups_unsatisfied"),
        ("candidate_structure", "raw_parent_edge_trees_T_complete_units", "raw_tree_candidates_T_complete_units"),
        ("candidate_structure", "distinct_optimal_vertex_sets_V_complete_units", "distinct_vertex_sets_V_complete_units"),
        ("ranking_outcome", "unique_first", "quality_primary_unique_vertex_units"),
        ("ranking_outcome", "tied_first", "quality_primary_tied_vertex_units"),
        ("ranking_outcome", "abstain_all_causes", "rank_abstain_units"),
        ("topology", "evaluated_units", "topology_evaluated_units"),
        ("topology", "coarse_class_unique_units", "coarse_topology_class_unique_units"),
        ("topology", "coarse_class_multiple_units", "coarse_topology_multiple_class_units"),
        ("topology", "parent_edge_assignment_unique_units", "parent_edge_assignment_unique_units"),
        ("topology", "exact_topology_proven_unique_units", "exact_topology_proven_unique_units"),
        ("effective_k", "component_site_mass", "k_component_sites_total"),
        ("effective_k", "observed_ALT_active_mass", "k_observed_alt_active_total"),
        ("effective_k", "not_structural_ALT_active_mass", "not_structural_alt_active_sites_total"),
    )
    for group, summary_key, raw_key in scalar_bindings:
        _require(
            (cell.get(group) or {}).get(summary_key) == raw.get(raw_key),
            f"M2 numeric summary {label}/{group}.{summary_key} 未精確 reconcile S7 raw cell",
            errors,
        )
    mapping_bindings = (
        ("ranking_outcome", "selection_status_counts", "selection_status_counts"),
        ("topology", "coarse_unique_class_counts", "coarse_topology_unique_class_counts"),
        ("topology", "coarse_ambiguous_class_set_counts", "coarse_topology_ambiguous_class_set_counts"),
        ("effective_k", "k_route_counts", "k_route_counts"),
    )
    for group, summary_key, raw_key in mapping_bindings:
        _require(
            (cell.get(group) or {}).get(summary_key) == (raw.get(raw_key) or {}),
            f"M2 numeric summary {label}/{group}.{summary_key} 未精確 reconcile S7 raw cell",
            errors,
        )
    _require(
        (cell.get("partial_read_funnel") or {}).get("full_detail")
        == (raw.get("partial_pattern_funnel") or {}),
        f"M2 numeric summary {label}/partial full_detail 未精確 reconcile S7 raw cell",
        errors,
    )


def _validate_candidate_against_stream(
    candidate: Mapping[str, Any], streamed: Mapping[str, Any], label: str,
    errors: list[str],
) -> None:
    scalar_bindings = {
        "n_solver_complete_candidate_units": "n_units",
        "n_candidate_vertex_sets_V": "n_candidate_vertex_sets_V",
        "n_parent_edge_trees_T": "n_parent_edge_trees_T",
        "unique_first": "unique_first",
        "tied_first": "tied_first",
        "solver_complete_optimizer_abstain": "solver_complete_optimizer_abstain",
    }
    for summary_key, stream_key in scalar_bindings.items():
        _require(
            candidate.get(summary_key) == streamed.get(stream_key),
            f"M2 numeric summary {label}/candidate.{summary_key} 未精確 reconcile S9 stream",
            errors,
        )
    partition = candidate.get("tree_vertex_partition") or {}
    _require(
        partition.get("counts") == streamed.get("tree_vertex_partition_counts"),
        f"M2 numeric summary {label}/candidate tree_vertex_partition 未精確 reconcile S9 stream",
        errors,
    )
    _require(
        candidate.get("h_star_distribution") == streamed.get("h_star_distribution"),
        f"M2 numeric summary {label}/candidate h* 未精確 reconcile S9 stream",
        errors,
    )
    tied = candidate.get("tied_by_coarse_topology") or {}
    _require(
        tied.get("consistent") == streamed.get("tied_topology_consistent")
        and tied.get("inconsistent") == streamed.get("tied_topology_inconsistent"),
        f"M2 numeric summary {label}/candidate tied×Topo 未精確 reconcile S9 stream",
        errors,
    )
    topology = candidate.get("topology") or {}
    topology_bindings = {
        "evaluated_units": "topology_evaluated_units",
        "coarse_class_unique_units": "coarse_topology_class_unique_units",
        "coarse_class_multiple_units": "coarse_topology_multiple_class_units",
        "coarse_unique_class_counts": "coarse_topology_unique_class_counts",
        "coarse_ambiguous_class_set_counts": "coarse_topology_ambiguous_class_set_counts",
        "parent_edge_assignment_unique_units": "parent_edge_assignment_unique_units",
        "exact_topology_proven_unique_units": "exact_topology_proven_unique_units",
    }
    for summary_key, stream_key in topology_bindings.items():
        _require(
            topology.get(summary_key) == streamed.get(stream_key),
            f"M2 numeric summary {label}/candidate topology.{summary_key} 未精確 reconcile S9 stream",
            errors,
        )
    _require(
        candidate.get("coarse_topology_class_inclusion_counts")
        == streamed.get("coarse_topology_class_inclusion_counts"),
        f"M2 numeric summary {label}/candidate topology inclusion 未精確 reconcile S9 stream",
        errors,
    )


def _authenticated_json_from_ledger(
    identity: Any, base: Path, label: str, errors: list[str]
) -> tuple[dict[str, Any] | None, Path | None]:
    if not isinstance(identity, dict):
        errors.append(f"{label} identity 不是 object")
        return None, None
    raw_path = identity.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        errors.append(f"{label} path 缺失")
        return None, None
    path = _resolve_recorded_path(raw_path, base)
    try:
        raw, observed = _stable_regular_identity(path, label)
    except ReportGateError as exc:
        errors.append(str(exc))
        return None, None
    digest = hashlib.sha256(raw).hexdigest()
    _require(observed.st_size == identity.get("size_bytes"), f"{label} size_bytes 不符", errors)
    _require(digest == identity.get("sha256"), f"{label} SHA-256 不符", errors)
    sidecar = path.with_name(f"{path.name}.sha256")
    try:
        side_raw, _ = _stable_regular_identity(sidecar, f"{label} sidecar")
        fields = side_raw.decode("ascii", errors="strict").strip().split()
    except (ReportGateError, UnicodeError) as exc:
        errors.append(str(exc))
        return None, path
    _require(fields == [digest, path.name], f"{label} sidecar 未精確綁定", errors)
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        errors.append(f"{label} JSON 無法解析：{exc}")
        return None, path
    _require(isinstance(payload, dict), f"{label} JSON root 不是 object", errors)
    return payload if isinstance(payload, dict) else None, path


def _reaggregate_extraction_components_from_children(
    summary_payload: Mapping[str, Any], summary_source: Source,
    extraction_source: Source, errors: list[str],
) -> dict[tuple[str, str, str], dict[str, Any]]:
    ledger_rows = (summary_payload.get("source_ledger") or {}).get("extraction_children")
    _require(isinstance(ledger_rows, list) and len(ledger_rows) == 154, "M2 numeric summary extraction_children ledger 不是 154 筆", errors)
    if not isinstance(ledger_rows, list):
        return {}
    result_rows = (extraction_source.payload or {}).get("results") or []
    result_index = {
        (row.get("dataset"), row.get("chrom")): row
        for row in result_rows if isinstance(row, dict)
    }
    expected = {(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES}
    observed: set[tuple[str, str]] = set()
    accumulators: dict[tuple[str, str, str], dict[str, Any]] = {}
    for identity in ledger_rows:
        dataset = identity.get("dataset") if isinstance(identity, dict) else None
        chrom = identity.get("chrom") if isinstance(identity, dict) else None
        key = (dataset, chrom)
        if key in observed:
            errors.append(f"M2 numeric summary extraction child ledger 重複：{dataset}/{chrom}")
            continue
        observed.add(key)
        child, _ = _authenticated_json_from_ledger(
            identity, summary_source.path, f"M2 extraction child {dataset}/{chrom}", errors
        )
        if child is None:
            continue
        _require(
            child.get("schema_name") == "intersubmod.lossless_read_linkage_chromosome_receipt"
            and child.get("schema_version") == "1.2.0"
            and child.get("all_pass") is True,
            f"M2 extraction child {dataset}/{chrom} schema/all_pass 不符",
            errors,
        )
        scope = child.get("scope") or {}
        _require(scope.get("dataset") == dataset and scope.get("chrom") == chrom, f"M2 extraction child {dataset}/{chrom} scope 不符", errors)
        terminal_row = result_index.get(key) or {}
        _require(terminal_row.get("receipt") == child, f"M2 extraction child {dataset}/{chrom} 未精確綁定 S6 embedded receipt", errors)
        summaries = child.get("component_summary_by_linkage_basis") or {}
        for basis in ("PS_HP1", "PS_HP2"):
            for threshold in ("1", "2", "3", "5"):
                raw_cell = (summaries.get(basis) or {}).get(threshold)
                if not isinstance(raw_cell, dict):
                    errors.append(f"M2 extraction child {dataset}/{chrom}/{basis}/{threshold} component cell 缺失")
                    continue
                target = accumulators.setdefault(
                    (str(dataset), basis, threshold),
                    {
                        "n_components": 0, "n_singletons_k1": 0,
                        "n_multisite_k_gt1": 0, "max_k_component_sites": 0,
                        "n_active_site_memberships": 0,
                        "k_component_sites_distribution": Counter(),
                    },
                )
                for count_key in (
                    "n_components", "n_singletons_k1", "n_multisite_k_gt1",
                    "n_active_site_memberships",
                ):
                    value = raw_cell.get(count_key)
                    _require(_is_nonnegative_int(value), f"M2 extraction child {dataset}/{chrom}/{basis}/{threshold}.{count_key} 無效", errors)
                    if _is_nonnegative_int(value):
                        target[count_key] += value
                max_k = raw_cell.get("max_k_component_sites")
                _require(_is_nonnegative_int(max_k), f"M2 extraction child {dataset}/{chrom}/{basis}/{threshold}.max_k 無效", errors)
                if _is_nonnegative_int(max_k):
                    target["max_k_component_sites"] = max(target["max_k_component_sites"], max_k)
                distribution = raw_cell.get("k_component_sites_distribution") or {}
                _require(_count_mapping(distribution), f"M2 extraction child {dataset}/{chrom}/{basis}/{threshold}.k distribution 無效", errors)
                if _count_mapping(distribution):
                    target["k_component_sites_distribution"].update(distribution)
    _require(observed == expected, "M2 numeric summary extraction child ledger scope 不是完整 7×chr1-22", errors)
    return {
        key: {
            **{name: value for name, value in target.items() if name != "k_component_sites_distribution"},
            "k_component_sites_distribution": dict(sorted(target["k_component_sites_distribution"].items(), key=lambda item: int(item[0]))),
        }
        for key, target in accumulators.items()
    }


def _validate_summary_source_binding(
    identity: Any,
    expected: Source,
    summary_source: Source,
    label: str,
    errors: list[str],
) -> None:
    _require(isinstance(identity, dict), f"M2 numeric summary source_ledger.{label} 缺失", errors)
    if not isinstance(identity, dict):
        return
    raw_path = identity.get("path")
    _require(isinstance(raw_path, str) and bool(raw_path), f"M2 numeric summary {label} path 缺失", errors)
    if isinstance(raw_path, str) and raw_path:
        _require(
            _resolve_recorded_path(raw_path, summary_source.path) == expected.path,
            f"M2 numeric summary {label} path 未綁定教授版 input",
            errors,
        )
    _require(identity.get("sha256") == expected.sha256, f"M2 numeric summary {label} SHA-256 未精確綁定教授版 input", errors)
    _require(identity.get("size_bytes") == expected.path.stat().st_size, f"M2 numeric summary {label} size_bytes 不符教授版 input", errors)


def _validate_numeric_summary_producer_binding(
    identity: Any,
    summary_source: Source,
    errors: list[str],
) -> None:
    """Bind S16 to the numeric producer frozen beside this HTML builder."""

    _require(
        isinstance(identity, dict),
        "M2 numeric summary source_ledger.producer 缺失",
        errors,
    )
    if not isinstance(identity, dict):
        return
    frozen_sibling = Path(__file__).resolve().with_name(
        "build_final_numeric_summary.py"
    )
    try:
        frozen_payload, frozen_observed = _stable_regular_identity(
            frozen_sibling, "M2 numeric summary frozen sibling producer"
        )
    except ReportGateError as exc:
        errors.append(str(exc))
        return
    raw_path = identity.get("path")
    _require(
        isinstance(raw_path, str) and bool(raw_path),
        "M2 numeric summary producer path 缺失",
        errors,
    )
    if not isinstance(raw_path, str) or not raw_path:
        return
    recorded = Path(raw_path)
    if not recorded.is_absolute():
        recorded = summary_source.path.parent / recorded
    try:
        recorded_payload, recorded_observed = _stable_regular_identity(
            recorded, "M2 numeric summary recorded live producer"
        )
    except ReportGateError as exc:
        errors.append(str(exc))
        return
    recorded_sha = hashlib.sha256(recorded_payload).hexdigest()
    _require(
        identity.get("size_bytes") == recorded_observed.st_size,
        "M2 numeric summary producer size_bytes 與 recorded live producer 不符",
        errors,
    )
    _require(
        identity.get("sha256") == recorded_sha,
        "M2 numeric summary producer SHA-256 與 recorded live producer 不符",
        errors,
    )
    _require(
        recorded_observed.st_size == frozen_observed.st_size
        and recorded_sha == hashlib.sha256(frozen_payload).hexdigest(),
        "M2 numeric summary recorded live producer bytes 與 HTML builder 同目錄 frozen sibling 不符",
        errors,
    )


def _validate_m2_numeric_summary(
    payload: Mapping[str, Any],
    source: Source,
    extraction_source: Source,
    ranking_source: Source,
    verification_source: Source,
    candidate_source: Source,
) -> list[str]:
    errors: list[str] = []
    _require(payload.get("schema_name") == "intersubmod.m2_final_numeric_summary", "M2 numeric summary schema_name 不符", errors)
    _require(payload.get("schema_version") == "1.0.0", "M2 numeric summary schema_version 不符", errors)
    _require(payload.get("task_type") == "B_COMPREHENSIVE_VALIDATION", "M2 numeric summary 不是 Task-B", errors)
    _require(payload.get("all_pass") is True, "M2 numeric summary all_pass 不是 true", errors)
    scope = payload.get("scope") or {}
    _require(tuple(scope.get("datasets") or ()) == DATASETS, "M2 numeric summary datasets/order 不符", errors)
    _require(tuple(scope.get("chromosomes") or ()) == AUTOSOMES, "M2 numeric summary chromosomes/order 不是 chr1-22", errors)
    _require(
        scope.get("expected_tasks_per_stage") == 154
        and scope.get("n_technical_datasets") == 7
        and scope.get("n_biological_samples") == 6
        and scope.get("scope_is_canonical_full_7_dataset_chr1_22") is True,
        "M2 numeric summary canonical 7-dataset×chr1-22 scope 不完整",
        errors,
    )
    checks = payload.get("checks") or {}
    _require(
        isinstance(checks, dict)
        and NUMERIC_SUMMARY_REQUIRED_CHECKS <= set(checks)
        and all(value is True for value in checks.values()),
        "M2 numeric summary checks 不完整或未全數 true",
        errors,
    )
    parameters = payload.get("primary_parameters") or {}
    _require(parameters.get("bridge_thresholds") == [1, 2, 3, 5], "M2 numeric summary bridge thresholds 不是 1/2/3/5", errors)
    _require(parameters.get("HP_bases") == ["PS_HP1", "PS_HP2"], "M2 numeric summary HP bases 不符", errors)
    _require(parameters.get("primary_structural_exact_pattern_minread") == 3, "M2 numeric summary structural minread 不是 3", errors)
    integrity = payload.get("receipt_integrity") or {}
    _require(
        integrity == {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{source.path.name}.sha256",
        },
        "M2 numeric summary receipt_integrity/sidecar_name 不符",
        errors,
    )
    ledger = payload.get("source_ledger") or {}
    _validate_numeric_summary_producer_binding(ledger.get("producer"), source, errors)
    _validate_summary_source_binding(ledger.get("terminal_extraction"), extraction_source, source, "terminal_extraction", errors)
    _validate_summary_source_binding(ledger.get("terminal_ranking"), ranking_source, source, "terminal_ranking", errors)
    _validate_summary_source_binding(ledger.get("final_independent_verification"), verification_source, source, "final_independent_verification", errors)
    _validate_summary_source_binding(ledger.get("candidate_table"), candidate_source, source, "candidate_table", errors)

    unsupported = payload.get("unsupported_or_nonidentifiable") or {}
    _require(isinstance(unsupported, dict) and set(unsupported) == NUMERIC_SUMMARY_UNSUPPORTED_KEYS, "M2 numeric summary unsupported/null metrics 不完整", errors)
    if isinstance(unsupported, dict):
        for key in NUMERIC_SUMMARY_UNSUPPORTED_KEYS:
            metric = unsupported.get(key)
            _require(
                isinstance(metric, dict)
                and metric.get("value") is None
                and isinstance(metric.get("reason"), str)
                and bool(metric.get("reason")),
                f"M2 numeric summary unsupported.{key} 未以 null+reason 表示",
                errors,
            )

    extraction = payload.get("extraction") or {}
    extraction_overall = extraction.get("overall_component_by_linkage_basis_threshold") or {}
    extraction_by_dataset = extraction.get("by_dataset") or {}
    ranking = payload.get("ranking") or {}
    ranking_overall = ranking.get("overall_by_HP_basis_and_bridge_threshold") or {}
    ranking_by_dataset = ranking.get("by_dataset") or {}
    extraction_child_cells = _reaggregate_extraction_components_from_children(
        payload, source, extraction_source, errors
    )
    s6_component_cells = ((extraction_source.payload or {}).get("aggregate") or {}).get("component_summary_by_linkage_basis") or {}
    s7_aggregate_cells = ((ranking_source.payload or {}).get("aggregate") or {}).get("by_linkage_basis_and_threshold") or {}
    s7_by_dataset = (ranking_source.payload or {}).get("by_dataset") or {}
    streamed_by_stratum = (candidate_source.payload or {}).get("by_stratum") or {}
    _require(set(extraction_by_dataset) == set(DATASETS), "M2 numeric summary extraction.by_dataset 不完整", errors)
    _require(set(ranking_by_dataset) == set(DATASETS), "M2 numeric summary ranking.by_dataset 不完整", errors)
    threshold_set = {"1", "2", "3", "5"}
    for basis in ("PS_HP1", "PS_HP2"):
        _require(isinstance(extraction_overall.get(basis), dict) and set(extraction_overall.get(basis, {})) == threshold_set, f"M2 numeric summary extraction overall {basis} thresholds 不完整", errors)
        _require(isinstance(ranking_overall.get(basis), dict) and set(ranking_overall.get(basis, {})) == threshold_set, f"M2 numeric summary ranking overall {basis} thresholds 不完整", errors)
        for threshold in sorted(threshold_set, key=int):
            overall_component = (extraction_overall.get(basis) or {}).get(threshold)
            overall_rank_cell = (ranking_overall.get(basis) or {}).get(threshold)
            _validate_numeric_component_cell(overall_component, f"extraction/overall/{basis}/{threshold}", errors)
            _validate_numeric_rank_cell(overall_rank_cell, f"ranking/overall/{basis}/{threshold}", errors)
            raw_component = (s6_component_cells.get(basis) or {}).get(threshold)
            raw_rank_overall = (s7_aggregate_cells.get(basis) or {}).get(threshold)
            _require(isinstance(raw_component, dict), f"S6 extraction raw component cell 缺失：{basis}/{threshold}", errors)
            _require(isinstance(raw_rank_overall, dict), f"S7 ranking raw aggregate cell 缺失：{basis}/{threshold}", errors)
            if isinstance(overall_component, dict) and isinstance(raw_component, dict):
                _validate_component_against_raw(overall_component, raw_component, f"extraction/overall/{basis}/{threshold}", errors)
            if isinstance(overall_rank_cell, dict) and isinstance(raw_rank_overall, dict):
                _validate_rank_against_raw(overall_rank_cell, raw_rank_overall, f"ranking/overall/{basis}/{threshold}", errors)
            dataset_components: list[Mapping[str, Any]] = []
            dataset_rank_cells: list[Mapping[str, Any]] = []
            for dataset in DATASETS:
                component_nested = ((extraction_by_dataset.get(dataset) or {}).get("component_by_linkage_basis_threshold") or {})
                rank_nested = ((ranking_by_dataset.get(dataset) or {}).get("by_HP_basis_and_bridge_threshold") or {})
                _require(isinstance(component_nested.get(basis), dict) and set(component_nested.get(basis, {})) == threshold_set, f"M2 numeric summary extraction/{dataset}/{basis} thresholds 不完整", errors)
                _require(isinstance(rank_nested.get(basis), dict) and set(rank_nested.get(basis, {})) == threshold_set, f"M2 numeric summary ranking/{dataset}/{basis} thresholds 不完整", errors)
                component_cell = (component_nested.get(basis) or {}).get(threshold)
                rank_cell = (rank_nested.get(basis) or {}).get(threshold)
                _validate_numeric_component_cell(component_cell, f"extraction/{dataset}/{basis}/{threshold}", errors)
                _validate_numeric_rank_cell(rank_cell, f"ranking/{dataset}/{basis}/{threshold}", errors)
                raw_component_dataset = extraction_child_cells.get((dataset, basis, threshold))
                raw_rank_dataset = ((((s7_by_dataset.get(dataset) or {}).get("by_linkage_basis_and_threshold") or {}).get(basis) or {}).get(threshold))
                _require(isinstance(raw_component_dataset, dict), f"S16/S6-child extraction raw cell 缺失：{dataset}/{basis}/{threshold}", errors)
                _require(isinstance(raw_rank_dataset, dict), f"S7 ranking raw dataset cell 缺失：{dataset}/{basis}/{threshold}", errors)
                if isinstance(component_cell, dict) and isinstance(raw_component_dataset, dict):
                    _validate_component_against_raw(component_cell, raw_component_dataset, f"extraction/{dataset}/{basis}/{threshold}", errors)
                if isinstance(rank_cell, dict) and isinstance(raw_rank_dataset, dict):
                    _validate_rank_against_raw(rank_cell, raw_rank_dataset, f"ranking/{dataset}/{basis}/{threshold}", errors)
                    streamed = streamed_by_stratum.get((dataset, basis, threshold), _freeze_candidate_aggregate(_new_candidate_aggregate()))
                    _validate_candidate_against_stream(
                        ((rank_cell.get("candidate_structure") or {}).get("candidate_table") or {}),
                        streamed,
                        f"ranking/{dataset}/{basis}/{threshold}",
                        errors,
                    )
                if isinstance(component_cell, dict):
                    dataset_components.append(component_cell)
                if isinstance(rank_cell, dict):
                    dataset_rank_cells.append(rank_cell)
            if isinstance(overall_component, dict) and len(dataset_components) == len(DATASETS):
                for key in ("n_components", "k_equals_1", "k_greater_than_1", "active_site_membership_mass"):
                    _require(overall_component.get(key) == sum(cell.get(key, 0) for cell in dataset_components), f"M2 numeric summary extraction {basis}/{threshold} by-dataset {key} 加總不符 overall", errors)
                merged_distribution: Counter[str] = Counter()
                for cell in dataset_components:
                    merged_distribution.update(cell.get("k_distribution") or {})
                _require(dict(merged_distribution) == overall_component.get("k_distribution"), f"M2 numeric summary extraction {basis}/{threshold} by-dataset k distribution 加總不符 overall", errors)
                _require(overall_component.get("max_k_component_sites") == max((cell.get("max_k_component_sites", 0) for cell in dataset_components), default=0), f"M2 numeric summary extraction {basis}/{threshold} by-dataset max k 不符 overall", errors)
            if isinstance(overall_rank_cell, dict) and len(dataset_rank_cells) == len(DATASETS):
                scalar_paths = (
                    ("units", "n_component_hp_unit_evaluations"),
                    ("units", "solver_complete"),
                    ("units", "solver_incomplete_or_not_run"),
                    ("ranking_outcome", "unique_first"),
                    ("ranking_outcome", "tied_first"),
                    ("ranking_outcome", "abstain_all_causes"),
                    ("effective_k", "component_site_mass"),
                    ("effective_k", "observed_ALT_active_mass"),
                    ("effective_k", "not_structural_ALT_active_mass"),
                    ("candidate_structure", "raw_parent_edge_trees_T_complete_units"),
                    ("candidate_structure", "distinct_optimal_vertex_sets_V_complete_units"),
                )
                for group, key in scalar_paths:
                    _require(
                        overall_rank_cell.get(group, {}).get(key)
                        == sum(cell.get(group, {}).get(key, 0) for cell in dataset_rank_cells),
                        f"M2 numeric summary ranking {basis}/{threshold} by-dataset {group}.{key} 加總不符 overall",
                        errors,
                    )
                overall_candidate = overall_rank_cell.get("candidate_structure", {}).get("candidate_table", {})
                dataset_candidates = [cell.get("candidate_structure", {}).get("candidate_table", {}) for cell in dataset_rank_cells]
                for key in ("n_solver_complete_candidate_units", "unique_first", "tied_first", "solver_complete_optimizer_abstain"):
                    _require(overall_candidate.get(key) == sum(candidate.get(key, 0) for candidate in dataset_candidates), f"M2 numeric summary ranking {basis}/{threshold} candidate {key} 加總不符 overall", errors)
                merged_h: Counter[str] = Counter()
                for candidate in dataset_candidates:
                    merged_h.update(candidate.get("h_star_distribution") or {})
                _require(dict(merged_h) == overall_candidate.get("h_star_distribution"), f"M2 numeric summary ranking {basis}/{threshold} h* 加總不符 overall", errors)
                for key in ("consistent", "inconsistent", "denominator"):
                    _require(
                        (overall_candidate.get("tied_by_coarse_topology") or {}).get(key)
                        == sum((candidate.get("tied_by_coarse_topology") or {}).get(key, 0) for candidate in dataset_candidates),
                        f"M2 numeric summary ranking {basis}/{threshold} tied×Topo {key} 加總不符 overall",
                        errors,
                    )
                merged_stream = _new_candidate_aggregate()
                for dataset in DATASETS:
                    stream_cell = streamed_by_stratum.get((dataset, basis, threshold))
                    if stream_cell is None:
                        continue
                    for key in (
                        "n_units", "n_candidate_vertex_sets_V", "n_parent_edge_trees_T",
                        "unique_first", "tied_first", "solver_complete_optimizer_abstain",
                        "tied_topology_consistent", "tied_topology_inconsistent",
                        "topology_evaluated_units", "coarse_topology_class_unique_units",
                        "coarse_topology_multiple_class_units",
                        "parent_edge_assignment_unique_units",
                        "exact_topology_proven_unique_units",
                    ):
                        merged_stream[key] += stream_cell[key]
                    for key in TREE_VERTEX_BUCKETS:
                        merged_stream["tree_vertex_partition_counts"][key] += stream_cell["tree_vertex_partition_counts"].get(key, 0)
                    for key in (
                        "h_star_distribution", "coarse_topology_class_inclusion_counts",
                        "coarse_topology_unique_class_counts",
                        "coarse_topology_ambiguous_class_set_counts",
                    ):
                        merged_stream[key].update(stream_cell[key])
                _validate_candidate_against_stream(
                    overall_candidate,
                    _freeze_candidate_aggregate(merged_stream),
                    f"ranking/overall/{basis}/{threshold}",
                    errors,
                )
    return errors


def assess_sources(
    *,
    canonical_json: Path,
    funnel_receipt: Path | None,
    funnel_verification_receipt: Path | None,
    m0_receipt: Path | None,
    m0_verification_receipt: Path | None,
    pilot_receipt: Path,
    method_audit: Path,
    literature_audit: Path,
    m2_extraction_receipt: Path | None,
    m2_ranking_receipt: Path | None,
    m2_verification_receipt: Path | None,
    m2_numeric_summary: Path | None,
    m2_pilot_extraction_receipt: Path | None,
) -> Assessment:
    result = Assessment()
    required_specs = (
        ("S1", "Current canonical layered topology", canonical_json, "chr1–22；7 datasets / 6 biological samples", True),
        ("S3", "Symbolic/MILP/likelihood pilot receipt", pilot_receipt, "演算法控制；PILOT_NOT_FINAL_VALIDATION", True),
        ("S4", "Partial-read method audit", method_audit, "方法語意、反例與 solver/likelihood contract", False),
        ("S5", "Primary-literature and claim-boundary audit", literature_audit, "方法文獻與可宣稱邊界；不替代 quantitative receipts", False),
    )
    for source_id, label, path, scope, is_json in required_specs:
        try:
            result.sources[source_id] = _source(source_id, label, path, scope, json_payload=is_json)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"{label} 無法讀取：{exc}")
    if "S1" in result.sources:
        result.hard_issues.extend(_validate_canonical(result.sources["S1"].payload))
    if "S3" in result.sources:
        result.hard_issues.extend(_validate_pilot(result.sources["S3"].payload, result.sources["S3"]))

    if funnel_receipt is None:
        result.final_issues.append("缺少 source-backed current sSNV funnel receipt")
    elif "S1" not in result.sources:
        result.hard_issues.append("有 current funnel receipt，但沒有可驗證的 canonical JSON")
    else:
        try:
            result.sources["S10"] = _source(
                "S10",
                "Source-backed current ClairS-to-tree funnel receipt",
                funnel_receipt,
                "7 datasets；raw ClairS → LongPhase-S PASS → chr1-22 biallelic → retained",
            )
            result.hard_issues.extend(
                _validate_funnel(
                    result.sources["S10"].payload,
                    result.sources["S1"].payload,
                    result.sources["S1"].sha256,
                )
            )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"Current funnel receipt 無法讀取：{exc}")

    if funnel_verification_receipt is None:
        result.final_issues.append("缺少 current sSNV funnel independent verification receipt")
    elif "S1" not in result.sources or "S10" not in result.sources:
        result.final_issues.append("Current funnel independent verifier 無法套用：canonical/funnel source 尚未通過")
    else:
        try:
            result.sources["S14"] = _source(
                "S14",
                "Current sSNV funnel independent verification receipt",
                funnel_verification_receipt,
                "7 dataset ledgers；322 checks；canonical/producer path+SHA binding；branch conservation",
            )
            result.hard_issues.extend(
                _validate_funnel_verification(
                    result.sources["S14"].payload,
                    result.sources["S14"],
                    result.sources["S1"],
                    result.sources["S10"],
                )
            )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"Current funnel independent verifier receipt 無法讀取：{exc}")

    if m0_receipt is None:
        result.final_issues.append("缺少 M0 full Task-B receipt")
    else:
        try:
            result.sources["S2"] = _source("S2", "M0 ambiguity census receipt", m0_receipt, "舊 50 kb / capped baseline 的工程歧義 census")
            hard, final = _validate_m0(result.sources["S2"].payload)
            result.hard_issues.extend(hard)
            result.final_issues.extend(final)
            if "S1" in result.sources:
                canonical_aggregate = result.sources["S1"].payload["canonical"]["aggregate"]
                m0_payload = result.sources["S2"].payload
                m0_population = m0_payload.get("population") or {}
                capped_units = (m0_population.get("excluded_hp_lineage_unit_counts") or {}).get("CAPPED")
                if _is_nonnegative_int(capped_units) and _is_nonnegative_int((m0_payload.get("aggregate") or {}).get("n_hp_lineage_units")):
                    _require(
                        m0_payload["aggregate"]["n_hp_lineage_units"] + capped_units == canonical_aggregate.get("primary_units"),
                        "M0 eligible + capped HP units 不等於 canonical primary_units",
                        result.hard_issues,
                    )
                _require(
                    m0_population.get("n_primary_mutation_regions") == canonical_aggregate.get("W_primary"),
                    "M0 primary regions 不等於 canonical W_primary",
                    result.hard_issues,
                )
                _require(
                    m0_population.get("n_fully_m0_eligible_regions") == canonical_aggregate.get("complete_regions"),
                    "M0 fully eligible regions 不等於 canonical complete_regions",
                    result.hard_issues,
                )
            rows_source, rows_errors = _m0_rows_source(result.sources["S2"])
            result.hard_issues.extend(rows_errors)
            if rows_source is not None:
                result.sources["S12"] = rows_source
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"M0 receipt 無法讀取：{exc}")

    if m0_verification_receipt is None:
        result.final_issues.append("缺少 M0 independent deep verification receipt")
    elif "S2" not in result.sources or "S12" not in result.sources:
        result.hard_issues.append("有 M0 verifier receipt，但缺少可驗證的 M0 receipt/rows")
    else:
        try:
            result.sources["S11"] = _source(
                "S11",
                "M0 independent deep verification receipt",
                m0_verification_receipt,
                "row-level aggregate、T/V partition、canonical reconciliation、64,973 candidate deep checks",
            )
            result.hard_issues.extend(
                _validate_m0_verification(
                    result.sources["S11"].payload,
                    result.sources["S2"],
                    result.sources["S12"],
                )
            )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"M0 independent verifier receipt 無法讀取：{exc}")

    if m2_extraction_receipt is None:
        result.final_issues.append("缺少 M2 full extraction receipt")
    else:
        try:
            result.sources["S6"] = _source("S6", "M2 full read-linked extraction receipt", m2_extraction_receipt, "chr1–22 × 7 datasets；154 chromosome tasks")
            extraction_payload = result.sources["S6"].payload
            extraction_version = extraction_payload.get("schema_version")
            if extraction_version == "1.2.0":
                result.hard_issues.extend(_validate_m2_extraction(extraction_payload, require_ps_aware=True))
            elif extraction_version == "1.1.0":
                result.hard_issues.extend(_validate_m2_extraction(extraction_payload, require_ps_aware=False))
                result.final_issues.append("M2 extraction schema 1.1.0 僅可 diagnostic；final 要求 child/full 1.2.0")
            else:
                result.hard_issues.append(f"M2 extraction 不支援的 schema_version：{extraction_version}")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"M2 extraction receipt 無法讀取：{exc}")

    if m2_ranking_receipt is None:
        result.final_issues.append("缺少 M2 full ranking receipt")
    elif "S6" not in result.sources:
        result.hard_issues.append("有 M2 ranking receipt，但沒有可驗證的 extraction receipt")
    else:
        try:
            result.sources["S7"] = _source("S7", "M2 full symbolic/likelihood ranking receipt", m2_ranking_receipt, "linkage basis × bridge threshold × dataset；154 tasks")
            ranking_payload = result.sources["S7"].payload
            ranking_version = ranking_payload.get("schema_version")
            if ranking_version == "2.0.0":
                result.hard_issues.extend(_validate_m2_ranking(ranking_payload, result.sources["S6"].sha256))
                candidate_source, candidate_errors = _candidate_table_source(result.sources["S7"])
                result.hard_issues.extend(candidate_errors)
                if candidate_source is not None:
                    result.sources["S9"] = candidate_source
            elif ranking_version == "1.0.0" and ranking_payload.get("schema_name") == "intersubmod.m2_full_ranking_receipt" and ranking_payload.get("all_pass") is True:
                result.final_issues.append("M2 ranking schema 1.0.0 混合 PS 風險未關閉；只可 diagnostic，final 要求 child/full 2.0.0")
            else:
                result.hard_issues.append(f"M2 ranking 不支援或未通過的 schema_version：{ranking_version}")
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"M2 ranking receipt 無法讀取：{exc}")

    if m2_verification_receipt is None:
        result.final_issues.append("缺少 M2 full independent verification receipt")
    elif not all(source_id in result.sources for source_id in ("S6", "S7", "S9")):
        result.final_issues.append("M2 independent verifier 無法套用：full extraction/ranking/candidate sources 尚未全部通過")
    else:
        try:
            result.sources["S13"] = _source(
                "S13",
                "M2 full independent verification receipt",
                m2_verification_receipt,
                "154 extraction + 154 ranking child receipts；aggregate reconstruction；candidate row-by-row reconstruction",
            )
            result.hard_issues.extend(
                _validate_m2_verification(
                    result.sources["S13"].payload,
                    result.sources["S6"],
                    result.sources["S7"],
                    result.sources["S9"],
                )
            )
            resource_source, resource_errors = _resource_session_source(
                result.sources["S13"]
            )
            result.hard_issues.extend(resource_errors)
            if resource_source is not None:
                result.sources["S15"] = resource_source
            release_source, release_errors = _release_manifest_source(
                result.sources["S13"]
            )
            result.hard_issues.extend(release_errors)
            if release_source is not None:
                result.sources["S17"] = release_source
                if resource_source is not None:
                    session_release = (resource_source.payload or {}).get("release_manifest") or {}
                    s13_release = (result.sources["S13"].payload.get("release_binding") or {}).get("release_manifest") or {}
                    _require(
                        session_release == s13_release,
                        "M2 extraction session release_manifest 未精確綁定 S13 physical release",
                        result.hard_issues,
                    )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"M2 independent verifier receipt 無法讀取：{exc}")

    if m2_numeric_summary is None:
        result.final_issues.append("缺少 M2 authenticated final numeric summary")
    elif not all(source_id in result.sources for source_id in ("S6", "S7", "S9", "S13")):
        result.final_issues.append(
            "M2 final numeric summary 無法套用：terminal extraction/ranking/final verifier 尚未全部通過"
        )
    else:
        numeric_source, numeric_errors = _sidecar_authenticated_source(
            "S16",
            "M2 authenticated final numeric summary",
            m2_numeric_summary,
            "presentation-stage全量數字；per-dataset×HP×bridge component/k、effective-k、h*、unique/tied×coarse-Topo",
        )
        result.hard_issues.extend(numeric_errors)
        if numeric_source is not None:
            result.sources["S16"] = numeric_source
            result.hard_issues.extend(
                _validate_m2_numeric_summary(
                    numeric_source.payload,
                    numeric_source,
                    result.sources["S6"],
                    result.sources["S7"],
                    result.sources["S13"],
                    result.sources["S9"],
                )
            )

    if m2_pilot_extraction_receipt is not None:
        try:
            result.sources["S8"] = _source("S8", "LEGACY schema-1.0 M2 HCC1954 chr22 extraction pilot", m2_pilot_extraction_receipt, "單一 dataset × 單一 chromosome；legacy diagnostic context only；不可替代 schema 1.2 resource pilot")
            result.hard_issues.extend(_validate_m2_pilot(result.sources["S8"].payload))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            result.hard_issues.append(f"M2 pilot extraction receipt 無法讀取：{exc}")
    return result


def _h(value: Any) -> str:
    return html.escape(str(value), quote=True)


def _n(value: int | float | None) -> str:
    if value is None:
        return "未評估"
    if isinstance(value, float) and not value.is_integer():
        return f"{value:,.4g}"
    return f"{int(value):,}"


def _pct(numerator: int | float | None, denominator: int | float | None, digits: int = 2) -> str:
    if numerator is None or denominator in (None, 0):
        return "不適用"
    return f"{100 * float(numerator) / float(denominator):.{digits}f}%"


def _metric_cell(
    value: int | float | None,
    unit: str,
    total_value: int | float | None,
    total_label: str,
    relative_value: int | float | None,
    relative_label: str,
) -> str:
    return (
        f'<span class="metric-main">{_h(_n(value))} <small>{_h(unit)}</small></span>'
        f'<span class="denom">總分母：{_h(total_label)} {_h(_n(total_value))}（{_h(_pct(value, total_value))}）</span>'
        f'<span class="denom">相對分母：{_h(relative_label)} {_h(_n(relative_value))}（{_h(_pct(value, relative_value))}）</span>'
    )


def _count_cell(value: int | float, unit: str, note: str) -> str:
    """Render an absolute count where a percentage denominator is meaningless."""
    return (
        f'<span class="metric-main">{_h(_n(value))} <small>{_h(unit)}</small></span>'
        f'<span class="denom">{_h(note)}</span>'
    )


def _source_ref(source_id: str) -> str:
    return f'<sup class="source-ref"><a href="#source-{_h(source_id)}" aria-label="查看來源 {_h(source_id)}">{_h(source_id)}</a></sup>'


def _bar_svg(
    rows: Sequence[tuple[str, int]],
    *,
    title: str,
    description: str,
    unit: str,
    color: str = "#b33424",
    width: int = 980,
) -> str:
    if not rows:
        return ""
    max_value = max(value for _, value in rows) or 1
    left, right, top, row_h = 205, 105, 48, 38
    inner = width - left - right
    height = top + row_h * len(rows) + 28
    svg_id = hashlib.sha1(title.encode("utf-8")).hexdigest()[:10]
    elements = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-labelledby="bar-title-{svg_id} bar-desc-{svg_id}">',
        f'<title id="bar-title-{svg_id}">{_h(title)}</title>',
        f'<desc id="bar-desc-{svg_id}">{_h(description)}</desc>',
    ]
    for index, (label, value) in enumerate(rows):
        y = top + index * row_h
        bar_w = inner * value / max_value
        elements.extend(
            (
                f'<text x="0" y="{y + 18}" class="svg-label">{_h(label)}</text>',
                f'<line x1="{left}" x2="{left + inner}" y1="{y + 12}" y2="{y + 12}" class="svg-track"/>',
                f'<rect x="{left}" y="{y + 3}" width="{bar_w:.2f}" height="18" fill="{color}"/>',
                f'<text x="{left + bar_w + 8:.2f}" y="{y + 18}" class="svg-value">{_h(_n(value))} {_h(unit)}</text>',
            )
        )
    elements.append("</svg>")
    return "".join(elements)


def _stacked_svg(
    rows: Sequence[tuple[str, int, str]],
    *,
    title: str,
    description: str,
    total: int,
    width: int = 980,
) -> str:
    x0, y0, bar_w, bar_h = 18, 52, width - 36, 42
    svg_id = hashlib.sha1(title.encode("utf-8")).hexdigest()[:10]
    elements = [
        f'<svg class="chart" viewBox="0 0 {width} {142 + len(rows) * 24}" role="img" aria-labelledby="stack-title-{svg_id} stack-desc-{svg_id}">',
        f'<title id="stack-title-{svg_id}">{_h(title)}</title>',
        f'<desc id="stack-desc-{svg_id}">{_h(description)}</desc>',
    ]
    cursor = x0
    for label, value, color in rows:
        segment = bar_w * value / total if total else 0
        elements.append(f'<rect x="{cursor:.2f}" y="{y0}" width="{segment:.2f}" height="{bar_h}" fill="{color}" stroke="#171511" stroke-width="1"/>')
        cursor += segment
    for index, (label, value, color) in enumerate(rows):
        y = 126 + index * 24
        elements.extend(
            (
                f'<rect x="18" y="{y - 12}" width="13" height="13" fill="{color}" stroke="#171511"/>',
                f'<text x="42" y="{y}" class="svg-label">{_h(label)} · {_h(_n(value))}（{_h(_pct(value, total))}）</text>',
            )
        )
    elements.append("</svg>")
    return "".join(elements)


def _partial_u_heatmap_svg(
    distributions: Mapping[str, Mapping[str, int]],
    denominators: Mapping[str, int],
) -> str:
    """Render u / 2^u composition without mixing the three denominators."""
    denominator_order = tuple(
        name for name in (
            "unique_patterns", "quality_groups", "molecule_projections",
            "structural_patterns",
        )
        if name in distributions and name in denominators
    )
    row_labels = {
        "unique_patterns": "Unique patterns",
        "quality_groups": "BQ-quality groups",
        "molecule_projections": "Molecule projections",
        "structural_patterns": "Structural unique patterns",
    }
    u_values = sorted({int(u) for name in denominator_order for u in distributions[name]})
    left, top, cell_w, cell_h = 190, 76, 68, 62
    width = max(980, left + cell_w * len(u_values) + 20)
    height = top + cell_h * len(denominator_order) + 42
    svg_id = hashlib.sha1(json.dumps(distributions, sort_keys=True).encode("utf-8")).hexdigest()[:10]
    elements = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-labelledby="partial-u-title-{svg_id} partial-u-desc-{svg_id}">',
        f'<title id="partial-u-title-{svg_id}">Partial-pattern u and 2^u distribution by denominator</title>',
        f'<desc id="partial-u-desc-{svg_id}">Three separate rows show unique patterns, quality groups, and molecule projections. Each cell reports count and within-row share for u unknown positions and 2 to the u conceptual completions.</desc>',
    ]
    for column, u_value in enumerate(u_values):
        x = left + column * cell_w + cell_w / 2
        elements.append(f'<text x="{x:.1f}" y="30" text-anchor="middle" class="svg-label">u={u_value}</text>')
        elements.append(f'<text x="{x:.1f}" y="50" text-anchor="middle" class="svg-note">2ᵘ={2 ** u_value:,}</text>')
    for row, denominator_name in enumerate(denominator_order):
        y = top + row * cell_h
        denominator = denominators[denominator_name]
        elements.append(f'<text x="0" y="{y + 24}" class="svg-label">{_h(row_labels[denominator_name])}</text>')
        elements.append(f'<text x="0" y="{y + 44}" class="svg-note">denom. {_h(_n(denominator))}</text>')
        for column, u_value in enumerate(u_values):
            value = int(distributions[denominator_name].get(str(u_value), 0))
            fraction = value / denominator if denominator else 0.0
            opacity = 0.10 + 0.72 * math.sqrt(fraction)
            x = left + column * cell_w
            elements.append(f'<rect x="{x}" y="{y}" width="{cell_w - 4}" height="{cell_h - 5}" fill="#b33424" fill-opacity="{opacity:.3f}" stroke="#171511"/>')
            elements.append(f'<text x="{x + (cell_w - 4)/2:.1f}" y="{y + 24}" text-anchor="middle" class="svg-value">{_h(_n(value))}</text>')
            elements.append(f'<text x="{x + (cell_w - 4)/2:.1f}" y="{y + 44}" text-anchor="middle" class="svg-note">{_h(_pct(value, denominator, 1))}</text>')
    elements.append("</svg>")
    return "".join(elements)


def _workflow_svg() -> str:
    stages = (
        ("LongPhase-S PASS", "sSNV + tagged evidence"),
        ("read-linked component", "unique molecule bridge"),
        ("R/A/O/D/S/L/X", "lossless sparse calls"),
        ("symbolic group", "N ∩ G(p) ≠ ∅"),
        ("minimal hidden", "optimal sets N（count = V）"),
        ("likelihood", "unique / tie / abstain"),
    )
    width, box_w, gap, height = 1120, 160, 27, 196
    elements = [
        f'<svg class="process-svg" viewBox="0 0 {width} {height}" role="img" aria-labelledby="flow-title flow-desc">',
        '<title id="flow-title">Read-linked Hypercube reconstruction flow</title>',
        '<desc id="flow-desc">From LongPhase-S PASS sites through molecule-defined components, symbolic group-Steiner candidate generation, and likelihood ranking.</desc>',
    ]
    for index, (title, subtitle) in enumerate(stages):
        x = 8 + index * (box_w + gap)
        if index:
            elements.append(f'<path d="M {x - gap + 4} 91 H {x - 8}" class="flow-arrow" marker-end="url(#arrow)"/>')
        elements.extend(
            (
                f'<rect x="{x}" y="45" width="{box_w}" height="92" class="flow-box"/>',
                f'<text x="{x + 12}" y="76" class="flow-title">{_h(title)}</text>',
                f'<text x="{x + 12}" y="104" class="flow-sub">{_h(subtitle)}</text>',
            )
        )
    elements.insert(3, '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L7,3 z" fill="#b33424"/></marker></defs>')
    elements.append("</svg>")
    return "".join(elements)


def _partial_svg() -> str:
    return """
<svg class="process-svg" viewBox="0 0 1060 330" role="img" aria-labelledby="partial-title partial-desc">
  <title id="partial-title">Partial read 的 symbolic group 與 likelihood marginalization</title>
  <desc id="partial-desc">RAX conceptually matches RAR and RAA. The solver uses one joint group constraint whose sparse row lists compatible active vertices, then likelihood marginalizes over the selected vertex set.</desc>
  <defs><marker id="p-arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L7,3 z" fill="#b33424"/></marker></defs>
  <text x="24" y="38" class="svg-kicker">一條 partial pattern</text>
  <rect x="24" y="58" width="160" height="72" class="flow-box accent"/>
  <text x="75" y="102" class="cube-code">RAX</text>
  <path d="M184 94 H252" class="flow-arrow" marker-end="url(#p-arrow)"/>
  <text x="264" y="38" class="svg-kicker">概念上的 2ᵘ 個 states</text>
  <rect x="264" y="58" width="130" height="72" class="flow-box"/>
  <rect x="414" y="58" width="130" height="72" class="flow-box"/>
  <text x="300" y="102" class="cube-code">RAR</text>
  <text x="450" y="102" class="cube-code">RAA</text>
  <text x="264" y="157" class="svg-note">不逐一另跑 completion 世界，也不跨 reads 做笛卡兒積</text>
  <path d="M544 94 H612" class="flow-arrow" marker-end="url(#p-arrow)"/>
  <text x="624" y="38" class="svg-kicker">一條 joint group constraint</text>
  <rect x="624" y="58" width="398" height="72" class="flow-box dark"/>
  <text x="654" y="102" class="cube-formula">N ∩ G(RAX) ≠ ∅</text>
  <path d="M824 130 V194" class="flow-arrow" marker-end="url(#p-arrow)"/>
  <text x="24" y="222" class="svg-kicker">全域目標</text>
  <rect x="24" y="242" width="326" height="64" class="flow-box"/>
  <text x="52" y="281" class="cube-formula">min hidden vertices</text>
  <path d="M350 274 H418" class="flow-arrow" marker-end="url(#p-arrow)"/>
  <rect x="430" y="242" width="258" height="64" class="flow-box"/>
  <text x="454" y="271" class="cube-formula">全部 optimal sets N</text>
  <text x="454" y="294" class="svg-note">僅限 complete = true</text>
  <path d="M688 274 H756" class="flow-arrow" marker-end="url(#p-arrow)"/>
  <rect x="768" y="222" width="254" height="84" class="flow-box accent"/>
  <text x="793" y="248" class="cube-formula">P(p|N,π) =</text>
  <text x="793" y="272" class="cube-formula">Σᵥ∈N πᵥP(p|v)</text>
  <text x="793" y="298" class="svg-note">X marginalize；read 只計一次</text>
</svg>
"""


def _relative_href(source_path: Path, output_path: Path) -> str:
    relative = os.path.relpath(source_path, start=output_path.parent)
    return quote(relative.replace(os.sep, "/"), safe="/._-~")


def _sum_count_trees(values: Sequence[Any]) -> Any:
    """Merge production ranking cells without treating metadata as counts.

    Ranking cells are mostly sparse count trees, but the production receipt also
    contains invariant strings (for example partial-pattern definitions) and
    boolean conservation contracts.  Counts are additive, invariant scalars
    must be present and exactly equal, and booleans are combined by logical AND.
    Missing sparse numeric categories are interpreted as zero; missing metadata
    is rejected instead of being silently invented.
    """

    def merge(nodes: Sequence[Any], expected: int, path: tuple[str, ...]) -> Any:
        if not nodes:
            raise ReportGateError(
                f"ranking cell 聚合缺少值：{'/'.join(path) or '<root>'}"
            )
        if all(isinstance(node, Mapping) for node in nodes):
            keys = sorted({key for node in nodes for key in node})
            merged: dict[str, Any] = {}
            for key in keys:
                present = [node[key] for node in nodes if key in node]
                merged[key] = merge(present, expected, (*path, str(key)))
            return merged

        missing = expected - len(nodes)
        if all(isinstance(node, bool) for node in nodes):
            if missing:
                raise ReportGateError(
                    f"ranking cell boolean contract 缺漏：{'/'.join(path)}"
                )
            return all(nodes)
        if any(isinstance(node, bool) for node in nodes):
            raise ReportGateError(
                f"ranking cell boolean/count 型別漂移：{'/'.join(path)}"
            )
        if all(_is_nonnegative_int(node) for node in nodes):
            return sum(int(node) for node in nodes)
        if all(isinstance(node, str) for node in nodes):
            if missing or any(node != nodes[0] for node in nodes[1:]):
                raise ReportGateError(
                    f"ranking cell invariant definition 不一致：{'/'.join(path)}"
                )
            return nodes[0]
        raise ReportGateError(
            f"ranking cells 含無法以 count/AND/exact-equal 合併的欄位："
            f"{'/'.join(path) or '<root>'}"
        )

    return merge(tuple(values), len(values), ())


def _primary_ranking_cell(payload: Mapping[str, Any]) -> tuple[str, str, Mapping[str, Any]]:
    """Return the professor-display cell, not a frozen primary bridge choice."""
    cells = payload["aggregate"]["by_linkage_basis_and_threshold"]
    bases = ("PS_HP1", "PS_HP2")
    threshold = DISPLAY_BRIDGE_THRESHOLD
    combined = _sum_count_trees([cells[basis][threshold] for basis in bases])
    return "PS_HP1 + PS_HP2", threshold, combined


def _canonical_sections(
    canonical: Mapping[str, Any],
    funnel_source: Source | None,
    funnel_verification_source: Source | None,
) -> str:
    aggregate = canonical["canonical"]["aggregate"]
    samples = canonical["canonical"]["samples"]
    funnel = funnel_source.payload if funnel_source is not None else None
    funnel_aggregate = (funnel or {}).get("aggregate") or {}
    funnel_by_dataset = {row["dataset"]: row for row in ((funnel or {}).get("datasets") or [])}
    topology = aggregate["topology_classes"]
    retained_chart = _bar_svg(
        [(row["sample"], row["retained_sSNV"]) for row in samples],
        title="各 dataset retained sSNV",
        description="Seven dataset counts from the current canonical chr1-22 summary; HCC1395 and HCC1395_DORADO are pipeline datasets for one biological sample.",
        unit="sSNV",
    )
    topo_rows = (
        ("T=1 且 Topo=1", topology["exact_and_topology_unique"], "#171511"),
        ("T>1 但 Topo=1", topology["topology_unique_exact_multiple"], "#b33424"),
        ("T>1 且 Topo>1", topology["topology_multiple_exact_multiple"], "#d98b78"),
        ("Incomplete／capped", topology["incomplete"], "#d9d1c2"),
    )
    topology_chart = _stacked_svg(
        topo_rows,
        title="Current canonical 50 kb baseline 的候選與拓撲狀態",
        description="Composition of W_primary regions; complete topology shares and incomplete regions are shown in one total denominator.",
        total=aggregate["W_primary"],
    )
    rows = []
    for row in samples:
        role = "同一 HCC1395 生物樣本的 Dorado pipeline dataset" if row["sample"] == "HCC1395_DORADO" else "dataset"
        funnel_row = funnel_by_dataset.get(row["sample"])
        upstream_cells = ""
        if funnel_row is not None:
            upstream_cells = (
                f'<td>{_metric_cell(funnel_row["raw_clairs_records"], "VCF record", funnel_aggregate["raw_clairs_records"], "7-dataset raw ClairS records", funnel_row["raw_clairs_records"], "該 dataset raw records")}</td>'
                f'<td>{_metric_cell(funnel_row["longphase_s_recalibrated_pass"], "sSNV record", funnel_aggregate["longphase_s_recalibrated_pass"], "7-dataset LongPhase-S PASS", funnel_row["raw_clairs_records"], "該 dataset raw records")}</td>'
            )
        rows.append(
            "<tr>"
            f'<th scope="row">{_h(row["sample"])}<span class="denom">Biological ID：{_h(BIOLOGICAL_IDS[row["sample"]])}<br>{_h(role)}</span></th>'
            f'{upstream_cells}'
            f'<td>{_metric_cell(row["autosomal_biallelic_sSNV"], "sSNV", aggregate["autosomal_biallelic_sSNV"], "7-dataset autosomal sSNV", row["tree_input_records"], "該 dataset input records")}</td>'
            f'<td>{_metric_cell(row["retained_sSNV"], "sSNV", aggregate["retained_sSNV"], "7-dataset retained sSNV", row["autosomal_biallelic_sSNV"], "該 dataset autosomal sSNV")}</td>'
            f'<td>{_metric_cell(row["W_tree"], "region", aggregate["W_tree"], "7-dataset legacy W_tree", row["W_tree"], "該 dataset W_tree")}</td>'
            f'<td>{_metric_cell(row["complete_regions"], "region", aggregate["complete_regions"], "7-dataset complete regions", row["W_primary"], "該 dataset W_primary")}</td>'
            "</tr>"
        )
    if funnel:
        raw = funnel_aggregate["raw_clairs_records"]
        tree = funnel_aggregate["longphase_s_recalibrated_pass"]
        autosomal = funnel_aggregate["autosomal_biallelic_sSNV"]
        retained = funnel_aggregate["retained_sSNV"]
        branches = funnel_aggregate["branch_counts"]
        funnel_metrics = f"""
  <div class="metric-ledger" aria-label="Source-backed site funnel headline metrics">
    <div>{_metric_cell(raw, "raw VCF record", raw, "全部 raw ClairS records", raw, "漏斗起點")}</div>
    <div>{_metric_cell(tree, "PASS sSNV record", raw, "全部 raw ClairS records", raw, "上一層 raw records")}</div>
    <div>{_metric_cell(autosomal, "chr1–22 biallelic sSNV", raw, "全部 raw ClairS records", tree, "上一層 LongPhase-S PASS")}</div>
    <div>{_metric_cell(retained, "retained sSNV", raw, "全部 raw ClairS records", autosomal, "上一層 chr1–22 biallelic sSNV")}</div>
  </div>
  <details open><summary>完整 site-level 漏斗：每個排除分支的單位、數量與分母</summary><div class="table-wrap"><table><thead><tr><th>分支</th><th>數量、總比例、相對比例</th><th>精確意義</th></tr></thead><tbody>
    <tr><th scope="row">LongPhase-S 非 PASS</th><td>{_metric_cell(branches["excluded_by_longphase_filter"], "sSNV record", raw, "全部 raw ClairS records", raw, "raw ClairS records")}</td><td>recalibration 後未進入 tree-input PASS 的 records。</td></tr>
    <tr><th scope="row">非 chr1–22 biallelic SNV</th><td>{_metric_cell(branches["out_of_scope_non_autosomal"], "sSNV record", raw, "全部 raw ClairS records", tree, "LongPhase-S PASS records")}</td><td>chrX／其他 contig／非 biallelic SNV；是分析範圍排除，不是 ONT 讀不到。</td></tr>
    <tr><th scope="row">MAX_SNV branch excluded</th><td>{_metric_cell(branches["max_snv_excluded"], "sSNV record", raw, "全部 raw ClairS records", autosomal, "chr1–22 biallelic sSNV")}</td><td>位於超過舊 MAX_SNV=8 計算容量分支的 sSNV records；不是 capped region 數。</td></tr>
    <tr><th scope="row">Positional singleton</th><td>{_metric_cell(branches["positional_singleton"], "sSNV record", raw, "全部 raw ClairS records", autosomal, "chr1–22 biallelic sSNV")}</td><td>舊 positional grouping 中沒有多位點建樹條件的單點 records。</td></tr>
    <tr><th scope="row">Retained</th><td>{_metric_cell(branches["retained"], "sSNV record", raw, "全部 raw ClairS records", autosomal, "chr1–22 biallelic sSNV")}</td><td>進入 current 50 kb candidate reconstruction 的 sSNV records。</td></tr>
  </tbody></table></div></details>
  {f'<p class="visual-note"><strong>獨立重算：</strong>{_h(funnel_verification_source.payload["n_pass"])} / {_h(funnel_verification_source.payload["n_checks"])} checks PASS，0 failures；canonical 與 producer funnel 的 path/SHA binding 均通過。{_source_ref("S14")}</p>' if funnel_verification_source is not None else ''}
"""
        upstream_headers = "<th>Raw ClairS</th><th>LongPhase-S PASS</th>"
        funnel_reference = _source_ref("S10") + (
            _source_ref("S14") if funnel_verification_source is not None else ""
        )
    else:
        funnel_metrics = '<p class="unavailable">Source-backed raw ClairS funnel receipt 未提供；raw 起點與排除分支保持未評估。</p>'
        upstream_headers = ""
        funnel_reference = ""
    return f"""
<section id="baseline" class="report-section">
  <p class="section-index">01 / CURRENT BASELINE</p>
  <h2>Current KPI 是 50 kb positional baseline；不是 read-defined clone region</h2>
  <p class="takeaway"><strong>這一層回答「現有 canonical 候選有多少」；不能回答「molecule 真正連成哪些區域」。</strong>
  sSNV 的篩選比例與 region 的完整性比例分開計算，避免把不同單位串成假漏斗。{_source_ref("S1")}{funnel_reference}</p>
  {funnel_metrics}
  <h3>單位由 site 轉為 region 後，重新設定分母</h3>
  <div class="table-wrap"><table><thead><tr><th>Region funnel</th><th>數量、總比例、相對比例</th><th>守恆關係</th></tr></thead><tbody>
    <tr><th scope="row">Legacy W_tree</th><td>{_metric_cell(aggregate["W_tree"], "region", aggregate["W_tree"], "全部 legacy W_tree regions", aggregate["W_tree"], "region 漏斗起點")}</td><td>site→region，不能沿用 sSNV 分母。</td></tr>
    <tr><th scope="row">W_primary</th><td>{_metric_cell(aggregate["W_primary"], "region", aggregate["W_tree"], "全部 legacy W_tree regions", aggregate["W_tree"], "legacy W_tree regions")}</td><td>W_primary + no-primary = W_tree。</td></tr>
    <tr><th scope="row">No primary lineage</th><td>{_metric_cell(aggregate["no_primary_lineage"], "region", aggregate["W_tree"], "全部 legacy W_tree regions", aggregate["W_tree"], "legacy W_tree regions")}</td><td>{aggregate["W_primary"]:,} + {aggregate["no_primary_lineage"]:,} = {aggregate["W_tree"]:,}。</td></tr>
    <tr><th scope="row">Complete candidate enumeration</th><td>{_metric_cell(aggregate["complete_regions"], "region", aggregate["W_primary"], "全部 W_primary regions", aggregate["W_primary"], "W_primary regions")}</td><td>complete + incomplete = W_primary。</td></tr>
    <tr><th scope="row">Incomplete / capped-limited candidate enumeration</th><td>{_metric_cell(aggregate["incomplete_regions"], "region", aggregate["W_primary"], "全部 W_primary regions", aggregate["W_primary"], "W_primary regions")}</td><td>此 {aggregate["incomplete_regions"]:,} regions 不等於 MAX_SNV branch 的 sSNV record 數，也不可直接稱唯一拓撲。</td></tr>
  </tbody></table></div>
  <h3>各 dataset retained sSNV</h3>
  <p class="visual-note">長條比較相同單位的 sSNV 數量；精確分母與樣本內保留率見下表。HCC1395 與 HCC1395_DORADO 是兩個資料流程、但只是一個 biological sample。</p>
  {retained_chart}
  <details open><summary>7 datasets / 6 biological samples：精確數量與兩種分母</summary>
    <div class="table-wrap"><table><thead><tr><th>Dataset / biological ID</th>{upstream_headers}<th>chr1–22 biallelic sSNV</th><th>Retained sSNV</th><th>Legacy W_tree</th><th>Complete regions</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
  </details>
  <h3>T 與 Topo 必須分開讀</h3>
  <p class="visual-note"><strong>Legacy T</strong> 是 parent-edge 樹候選數；<strong>Topo</strong> 是去除節點標籤後的形狀類別數。T&gt;1 但 Topo=1 代表邊或標籤候選仍多、形狀已一致；Topo&gt;1 才是形狀尚未唯一。Topo&gt;1 為 <strong>{topology["topology_multiple_exact_multiple"]:,}</strong> regions：占 W_primary {_pct(topology["topology_multiple_exact_multiple"], aggregate["W_primary"])}，或占 complete {_pct(topology["topology_multiple_exact_multiple"], aggregate["complete_regions"])}。兩個分母都必須明示。</p>
  {topology_chart}
  <details><summary>50 kb 與舊 k≤8 cap 到底代表什麼</summary><ul>
    <li><strong>50 kb 是 positional heuristic，不是 ONT 物理讀長上限。</strong> Current W 先按相鄰 sSNV gap&gt;50 kb 切開，當下完全不看 read；同一 W 的總 span 反而可以因鏈式相鄰而超過50 kb。</li>
    <li><strong>舊 MAX_SNV=8 是計算 cap，不是已證明的生物邊界。</strong> 把較大區域取8點不等於把它切成互相獨立且可拼回的 exact 子樹；所以 capped/incomplete 不可併入唯一拓撲比例。</li>
    <li><strong>M2 不先把 component 截成8點。</strong> 它保留整個 chromosome site catalog，以實際 unique-molecule bridge 定 component；solver 再依 observed-ALT effective k 路由。超出 exact gate 時標 local/abstain或使用有certificate的分解，不把任意切片冒充全域最佳。</li>
  </ul></details>
</section>
"""


def _m0_section(source: Source | None, canonical: Mapping[str, Any]) -> str:
    if source is None:
        return """
<section id="m0" class="report-section unavailable"><p class="section-index">02 / M0 AMBIGUITY CENSUS</p><h2>M0 全量 receipt 尚未提供</h2><p>缺值保持「未評估」，不以 0 取代。M0 只用舊 thresholded alignment-exposure patterns，最終 M2 molecule likelihood 不可由此代替。</p></section>
"""
    payload = source.payload
    aggregate = payload["aggregate"]
    population = payload["population"]
    canonical_aggregate = canonical["canonical"]["aggregate"]
    capped_hp_units = int((population.get("excluded_hp_lineage_unit_counts") or {}).get("CAPPED", 0))
    statuses = aggregate["selection_status_counts"]
    selected = payload.get("selected_datasets") or []
    status_rows = [(key, value) for key, value in sorted(statuses.items(), key=lambda item: (-item[1], item[0]))]
    chart = _bar_svg(status_rows, title="M0 candidate-selection status", description="Engineering ambiguity census over eligible HP lineage units; this is not lossless molecule likelihood.", unit="HP unit", color="#171511")
    scope_word = "全 7 datasets" if payload.get("full_task_b_scope") else f"subset：{len(selected)} / 7 datasets（{', '.join(selected)}）"
    rows = "".join(
        f'<tr><th scope="row">{_h(key)}</th><td>{_metric_cell(value, "HP unit", aggregate["n_hp_lineage_units"], "M0 eligible HP units", aggregate["n_hp_lineage_units"], "selection-status 母體")}</td></tr>'
        for key, value in status_rows
    )
    t1 = int(statuses.get("T1_CANDIDATE_UNIQUE", 0))
    edge_equivalent = int(statuses.get("T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED", 0))
    likelihood_unique = int(statuses.get("LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE", 0)) + int(
        statuses.get("LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED", 0)
    )
    likelihood_tied = int(statuses.get("LIKELIHOOD_TIED_VERTEX_SETS", 0))
    optimizer_abstain = int(statuses.get("RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE", 0))
    v_gt1 = likelihood_unique + likelihood_tied + optimizer_abstain
    t_gt1 = edge_equivalent + v_gt1
    categorical_rows = "".join(
        (
            '<tr><th scope="row">T=1、V=1</th><td>'
            + _metric_cell(t1, "HP unit", aggregate["n_hp_lineage_units"], "M0 eligible HP units", aggregate["n_hp_lineage_units"], "M0 eligible HP units")
            + '</td></tr>',
            '<tr><th scope="row">T&gt;1、V=1（只有 edge 不唯一）</th><td>'
            + _metric_cell(edge_equivalent, "HP unit", aggregate["n_hp_lineage_units"], "M0 eligible HP units", t_gt1, "T>1 units")
            + '</td></tr>',
            '<tr><th scope="row">T&gt;1、V&gt;1：likelihood 唯一第一 vertex set</th><td>'
            + _metric_cell(likelihood_unique, "HP unit", aggregate["n_hp_lineage_units"], "M0 eligible HP units", v_gt1, "T>1/V>1 units")
            + '</td></tr>',
            '<tr><th scope="row">T&gt;1、V&gt;1：likelihood 並列第一 vertex sets</th><td>'
            + _metric_cell(likelihood_tied, "HP unit", aggregate["n_hp_lineage_units"], "M0 eligible HP units", v_gt1, "T>1/V>1 units")
            + '</td></tr>',
            '<tr><th scope="row">T&gt;1、V&gt;1：optimizer abstain</th><td>'
            + _metric_cell(optimizer_abstain, "HP unit", aggregate["n_hp_lineage_units"], "M0 eligible HP units", v_gt1, "T>1/V>1 units")
            + '</td></tr>',
        )
    )
    dataset_rows = []
    for dataset in DATASETS:
        cell = (aggregate.get("by_dataset") or {}).get(dataset) or {}
        if not _is_nonnegative_int(cell.get("n_hp_lineage_units")):
            continue
        denominator = int(cell["n_hp_lineage_units"])
        ds_status = cell.get("selection_status_counts") or {}
        ds_t1 = int(ds_status.get("T1_CANDIDATE_UNIQUE", 0))
        ds_edge = int(ds_status.get("T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED", 0))
        ds_unique = int(ds_status.get("LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE", 0)) + int(
            ds_status.get("LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED", 0)
        )
        ds_tied = int(ds_status.get("LIKELIHOOD_TIED_VERTEX_SETS", 0))
        ds_abstain = int(ds_status.get("RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE", 0))
        ds_t_gt1 = denominator - ds_t1
        ds_v_gt1 = ds_unique + ds_tied + ds_abstain
        role = "HCC1395 的 Dorado pipeline dataset；非獨立生物樣本" if dataset == "HCC1395_DORADO" else "pipeline dataset"
        dataset_rows.append(
            '<tr>'
            f'<th scope="row">{_h(dataset)}<span class="denom">Biological ID：{_h(BIOLOGICAL_IDS[dataset])}<br>{_h(role)}</span></th>'
            f'<td>{_metric_cell(denominator, "HP unit", aggregate["n_hp_lineage_units"], "7-dataset eligible HP units", denominator, "該 dataset eligible HP units")}</td>'
            f'<td>{_metric_cell(ds_t1, "HP unit", aggregate["n_hp_lineage_units"], "7-dataset eligible HP units", denominator, "該 dataset eligible HP units")}</td>'
            f'<td>{_metric_cell(ds_edge, "HP unit", aggregate["n_hp_lineage_units"], "7-dataset eligible HP units", ds_t_gt1, "該 dataset T>1 units")}</td>'
            f'<td>{_metric_cell(ds_unique, "HP unit", aggregate["n_hp_lineage_units"], "7-dataset eligible HP units", ds_v_gt1, "該 dataset V>1 units")}</td>'
            f'<td>{_metric_cell(ds_tied, "HP unit", aggregate["n_hp_lineage_units"], "7-dataset eligible HP units", ds_v_gt1, "該 dataset V>1 units")}</td>'
            f'<td>{_metric_cell(ds_abstain, "HP unit", aggregate["n_hp_lineage_units"], "7-dataset eligible HP units", ds_v_gt1, "該 dataset V>1 units")}</td>'
            '</tr>'
        )
    dataset_table = ""
    if dataset_rows:
        dataset_table = (
            '<details><summary>各 dataset 的 T/V、likelihood unique/tie 與 abstain</summary>'
            '<div class="table-wrap"><table><thead><tr><th>Dataset / biological ID</th><th>Eligible HP units</th><th>T=1/V=1</th><th>T&gt;1/V=1</th><th>V&gt;1 unique</th><th>V&gt;1 tied</th><th>Optimizer abstain</th></tr></thead>'
            f'<tbody>{"".join(dataset_rows)}</tbody></table></div></details>'
        )
    return f"""
<section id="m0" class="report-section">
  <p class="section-index">02 / M0 AMBIGUITY CENSUS</p>
  <h2>M0 先分清「vertex ambiguity」與「edge-only ambiguity」</h2>
  <p class="takeaway"><strong>Scope：{_h(scope_word)}。</strong> M0 仍是舊 50 kb、MINREAD thresholded exposure 的工程 census；它可驗證計算與歧義分類，但不能升級為 lossless molecule likelihood 或真實 clone。數字同時由 row-level TSV 與 independent deep verifier 重算。{_source_ref("S2")}{_source_ref("S11")}{_source_ref("S12")}</p>
  <details open><summary>Primary HP-unit 漏斗：先說清楚 72,994 與 64,973 的差別</summary><div class="table-wrap"><table><thead><tr><th>HP-lineage-unit 層級</th><th>數量、總比例、相對比例</th><th>意義</th></tr></thead><tbody>
    <tr><th scope="row">Canonical primary units</th><td>{_metric_cell(canonical_aggregate["primary_units"], "HP unit", canonical_aggregate["primary_units"], "全部 primary HP units", canonical_aggregate["primary_units"], "HP-unit 漏斗起點")}</td><td>所有 W_primary 內的 primary HP lineage units；不是 region 或 clone。</td></tr>
    <tr><th scope="row">M0 eligible</th><td>{_metric_cell(aggregate["n_hp_lineage_units"], "HP unit", canonical_aggregate["primary_units"], "全部 primary HP units", canonical_aggregate["primary_units"], "全部 primary HP units")}</td><td>候選列舉完整、可進入 M0 likelihood census。</td></tr>
    <tr><th scope="row">Capped / excluded</th><td>{_metric_cell(capped_hp_units, "HP unit", canonical_aggregate["primary_units"], "全部 primary HP units", canonical_aggregate["primary_units"], "全部 primary HP units")}</td><td>{aggregate["n_hp_lineage_units"]:,} + {capped_hp_units:,} = {canonical_aggregate["primary_units"]:,}；不得把 excluded 當 0 或加入唯一率。</td></tr>
  </tbody></table></div></details>
  <p class="visual-note">狀態分母是 eligible HP lineage units，不是 region，也不是 cell clone。不同 HP unit 可屬於同一 region。</p>
  {chart}
  <details open><summary>把 T、V 與 likelihood 排名拆成互斥層級</summary><div class="table-wrap"><table><thead><tr><th>判讀層級</th><th>數量、總分母、相對分母</th></tr></thead><tbody>{categorical_rows}</tbody></table></div></details>
  <details open><summary>M0 selection status 精確數值</summary><div class="table-wrap"><table><thead><tr><th>Status</th><th>數量、總分母、相對分母</th></tr></thead><tbody>{rows}</tbody></table></div></details>
  {dataset_table}
</section>
"""


def _m2_extraction_section(source: Source | None) -> str:
    if source is None:
        return """
<section id="m2-extraction" class="report-section unavailable"><p class="section-index">03 / M2 READ-LINKED EXTRACTION</p><h2>M2 154-task extraction 尚未完成或未提供</h2><p>缺值保持「未評估」，不以 0 取代。因此目前不能把 canonical 50 kb region 改稱 read-linked component，也不能填寫 &gt;50 kb bridge、component k distribution 或 molecule funnel 的全量數字。</p></section>
"""
    payload = source.payload
    legacy_note = ""
    if payload.get("schema_version") != "1.2.0":
        legacy_note = '<p class="takeaway"><strong>DIAGNOSTIC ONLY：</strong>此 extraction receipt 是 schema 1.1.0；尚未證明 154 child receipts 全部符合 PS-aware 1.2.0，不能解除 PARTIAL。</p>'
    aggregate = payload["aggregate"]
    counts = aggregate["counts"]
    summaries = aggregate["component_summary_by_linkage_basis"]
    primary_available = all(basis in summaries for basis in ("PS_HP1", "PS_HP2"))
    if primary_available:
        component_cells: list[tuple[str, str, dict[str, Any], int]] = []
        thresholds = sorted(set(summaries["PS_HP1"]) | set(summaries["PS_HP2"]), key=int)
        for threshold in thresholds:
            hp_cells = [summaries[basis].get(threshold, {}) for basis in ("PS_HP1", "PS_HP2")]
            distribution: Counter[int] = Counter()
            for cell in hp_cells:
                distribution.update(
                    {int(key): int(value) for key, value in (cell.get("k_component_sites_distribution") or cell.get("k_distribution") or {}).items()}
                )
            combined_cell = {
                "n_components": sum(int(cell.get("n_components", 0)) for cell in hp_cells),
                "n_singletons_k1": sum(int(cell.get("n_singletons_k1", 0)) for cell in hp_cells),
                "n_multisite_k_gt1": sum(int(cell.get("n_multisite_k_gt1", 0)) for cell in hp_cells),
                "n_component_site_k_gt12": sum(value for key, value in distribution.items() if key > 12),
                "max_k_component_sites": max((int(cell.get("max_k_component_sites", 0)) for cell in hp_cells), default=0),
            }
            combined_n = combined_cell["n_components"]
            for basis, cell in zip(("PS_HP1", "PS_HP2"), hp_cells):
                component_cells.append((basis, threshold, dict(cell), combined_n))
            component_cells.append(("HP1+HP2 display sum", threshold, combined_cell, combined_n))
        component_scope = "PS-aware primary 分開顯示 PS_HP1 與 PS_HP2；combined 只是教授版 display sum。site membership 可在不同 PS 重複，故不是 unique genomic region 數"
    else:
        component_cells = [
            ("legacy pooled", threshold, dict(cell), int(cell.get("n_components", 0)))
            for threshold, cell in summaries["pooled"].items()
        ]
        component_scope = "legacy pooled diagnostic；缺少 PS_HP1/PS_HP2 primary summaries"
    funnel = (
        ("raw target-overlapping alignments", counts["raw_overlapping_alignments"], counts["raw_overlapping_alignments"], "raw alignments"),
        ("canonical eligible molecules", counts["unique_molecule_ids"], counts["raw_overlapping_alignments"], "raw alignments"),
        ("sidecar exact matches", counts["sidecar_exact_matches"], counts["unique_molecule_ids"], "eligible molecules"),
        ("fixed R/A calls", counts["fixed_ra_calls"], counts["fixed_ra_calls"], "fixed R/A call rows"),
        ("ALT calls", counts["alt_calls"], counts["fixed_ra_calls"], "fixed R/A calls"),
    )
    funnel_rows = "".join(
        f'<tr><th scope="row">{_h(label)}</th><td>{_metric_cell(value, "alignment/molecule/call（依列名）", counts["raw_overlapping_alignments"] if index < 3 else counts["fixed_ra_calls"], "該漏斗同單位總母體", relative, relative_label)}</td></tr>'
        for index, (label, value, relative, relative_label) in enumerate(funnel)
    )
    component_rows = []
    for basis_label, threshold, cell, threshold_total in component_cells:
        component_site_k_gt12 = cell.get("n_component_site_k_gt12")
        if component_site_k_gt12 is None:
            distribution = cell.get("k_component_sites_distribution") or cell.get("k_distribution") or {}
            component_site_k_gt12 = sum(int(value) for key, value in distribution.items() if int(key) > 12)
        component_rows.append(
            "<tr>"
            f'<th scope="row">{_h(basis_label)}<span class="denom">bridge ≥{_h(threshold)} unique molecule(s)</span></th>'
            f'<td>{_metric_cell(cell.get("n_components"), "component", threshold_total, "該 threshold HP1+HP2 display sum", cell.get("n_components"), "該 HP stratum components")}</td>'
            f'<td>{_metric_cell(cell.get("n_singletons_k1"), "component", cell.get("n_components"), "該 threshold components", cell.get("n_components"), "該 threshold components")}</td>'
            f'<td>{_metric_cell(cell.get("n_multisite_k_gt1"), "component", cell.get("n_components"), "該 threshold components", cell.get("n_components"), "該 threshold components")}</td>'
            f'<td>{_metric_cell(component_site_k_gt12, "component", cell.get("n_components"), "該 threshold components", cell.get("n_multisite_k_gt1"), "該 threshold multi-site components")}</td>'
            f'<td><span class="metric-main">{_h(_n(cell.get("max_k_component_sites")))} <small>site k</small></span><span class="denom">此 threshold 的最大 component-site k；不等於 effective solver k</span></td>'
            "</tr>"
        )
    return f"""
<section id="m2-extraction" class="report-section">
  <p class="section-index">03 / M2 READ-LINKED EXTRACTION</p>
  <h2>M2 以 unique molecule bridge 重建 component；50 kb 只保留作舊基準</h2>
  <p class="takeaway"><strong>154 / 154 chromosome tasks 通過後，區域單位才正式由 positional W 轉為 read-linked component。</strong> bridge thresholds 1/2/3/5 全保留，不能只挑一個結果講故事。{_source_ref("S6")}</p>
  {legacy_note}
  <details open><summary>Alignment → molecule → call 漏斗</summary><div class="table-wrap"><table><thead><tr><th>Stage</th><th>數量、總分母、相對分母</th></tr></thead><tbody>{funnel_rows}</tbody></table></div></details>
  <h3>PS-aware bridge threshold sensitivity</h3>
  <p class="visual-note">{_h(component_scope)}。threshold 是跨相鄰 active-site cut 的 unique molecule 支持下限；數值越高，component 通常越碎。這裡的 k 是 component-site k；真正 exact/abstain 路由依 ranking 的 observed-ALT effective k 判斷。</p>
  <div class="table-wrap"><table><thead><tr><th>HP stratum / bridge rule</th><th>Components</th><th>k=1</th><th>k&gt;1</th><th>Component-site k&gt;12</th><th>Max component-site k</th></tr></thead><tbody>{''.join(component_rows)}</tbody></table></div>
</section>
"""


def _candidate_examples_html(candidate_source: Source | None) -> str:
    if candidate_source is None:
        return ""
    metadata = candidate_source.payload["metadata"]
    groups = candidate_source.payload["example_groups"]
    minimum_extra_state_distribution = candidate_source.payload[
        "minimum_extra_state_distribution"
    ]
    rows: list[str] = []
    for category in ("UNIQUE", "TIE", "ABSTAIN", "OTHER"):
        for row in groups.get(category, []):
            parent_choices = int(row["parent_choice_count"])
            rows.append(
                "<tr>"
                f'<th scope="row"><code>{_h(row["unit_key"])}</code><span class="denom">{_h(row["dataset"])} · {_h(row["chrom"])} · HP {_h(row["hp_family"])} · known PS {_h(row["ps"])} · threshold {_h(row["threshold"])}</span></th>'
                f'<td><code>{_h(row["vertex_states"])}</code><span class="denom">roles：{_h(row["vertex_roles"])}</span></td>'
                f'<td>{_count_cell(parent_choices, "parent choices", "此 candidate/vertex set 的可行 parent assignments；絕對數")}</td>'
                f'<td><span class="metric-main">{_h(row["profile_log_likelihood"])} <small>log-likelihood</small></span><span class="denom">相對同 unit candidates：{_h(row["relative_log_likelihood"])} ΔlogL</span></td>'
                f'<td><code>{_h(row["mixture_weights_pi"])}</code><span class="denom">latent expected proportions；不是 cell clone assignment</span></td>'
                f'<td>{_h(row["winner_status"])}<span class="denom">tie group：{_h(row["tie_group"] or "—")}<br>coarse class：{_h(row["coarse_topology_class"])}</span></td>'
                "</tr>"
            )
    mean_rows_per_unit = metadata["n_rows"] / metadata["n_units"]
    mean_rows_note = (
        f'{metadata["n_rows"]} rows ÷ {metadata["n_units"]} units；算術比值，不是百分比'
    )
    minimum_extra_rows = "".join(
        "<tr>"
        f'<th scope="row">h* = {_h(extra_count)}</th>'
        f'<td>{_metric_cell(count, "solver-complete unit", metadata["n_units"], "candidate-table units", metadata["n_units"], "solver-complete candidate-bearing units")}</td>'
        "</tr>"
        for extra_count, count in sorted(
            minimum_extra_state_distribution.items(), key=lambda item: int(item[0])
        )
    )
    return f"""
  <details open><summary>實際區域候選：hashed candidate table 的 bounded examples</summary>
    <p><strong>π weights 只能稱 latent expected proportions。</strong> 它們是在固定候選 vertex set 與 conditional R/A emission 下的混合權重，不是把每個 molecule 指派到真實 cell clone。<strong>Winner 只是同一 unit、同一 minimum-extra-state h* 候選集內 likelihood 最高</strong>；沒有與 h*+1、h*+2 等所有可行結構比較，不能稱為「跨所有複雜度的全域最可能真實拓撲」。完整表由 S9 的 compressed-file SHA-256 鎖定；此處每種 outcome 只展示第一個 deterministic unit，避免把 examples 當總體比例。{_source_ref("S9")}</p>
    <div class="metric-ledger" aria-label="Candidate-table absolute counts and arithmetic ratio">
      <div>{_count_cell(metadata["n_rows"], "candidate rows", "hashed candidate-table rows；絕對數")}</div>
      <div>{_count_cell(metadata["n_units"], "PS-aware units", "candidate-table units；絕對數")}</div>
      <div>{_count_cell(mean_rows_per_unit, "rows / unit", mean_rows_note)}</div>
    </div>
    <h3>Minimum-extra-state count h*</h3>
    <p class="visual-note"><code>h*</code> 是 solver objective：候選 N 中除了 root 與 full-observed mandatory states 以外的 state 數。因此它同時包含 <strong>partial-compatible representative</strong> 與 <strong>connector</strong>；這是 minimum-extra-state count，<strong>不是 hidden clone 數</strong>。分母只是 hashed candidate table 覆蓋的 solver-complete units。</p>
    <div class="table-wrap"><table><thead><tr><th>Minimum-extra-state objective</th><th>Count / denominator</th></tr></thead><tbody>{minimum_extra_rows}</tbody></table></div>
    <p class="visual-note"><strong>Parent edges 此處只保留 <code>parent_choice_count</code> 計數，並未在 canonical candidate table 展開每一組 parent-edge assignment。</strong> 所以 parent choices &gt;1 只能報告 edge non-identifiability，不能假裝已選出唯一邊。</p>
    <div class="table-wrap"><table><thead><tr><th>Unit key / scope</th><th>Vertex states / roles</th><th>Parent choices</th><th>LL / relative score</th><th>π weights</th><th>Winner / tie / coarse class</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
  </details>
"""


def _sample_funnel_html(payload: Mapping[str, Any]) -> str:
    by_dataset = payload["by_dataset"]
    funnels = {dataset: by_dataset[dataset]["sample_funnel"] for dataset in DATASETS}
    aggregate = {
        key: sum(funnels[dataset][key] for dataset in DATASETS)
        for key in (
            "raw_sparse_molecules",
            "ps_known_molecules",
            "ps_missing_molecules",
            "hp_included_molecules",
            "hp_excluded_molecules",
        )
    }
    global_call_total = sum(sum(funnels[dataset]["call_class_counts"].values()) for dataset in DATASETS)
    chart = _bar_svg(
        [(dataset, funnels[dataset]["raw_sparse_molecules"]) for dataset in DATASETS],
        title="Raw sparse molecules · 7 datasets / 6 biological samples",
        description="Seven pipeline datasets are shown. HCC1395 and HCC1395_DORADO share one biological ID and are not independent biological samples.",
        unit="molecule",
    )
    rows = []
    for dataset in DATASETS:
        funnel = funnels[dataset]
        sample_call_total = sum(funnel["call_class_counts"].values())
        call_lines = "".join(
            f'<span class="denom"><strong>{_h(call_class)}</strong> {_h(_n(value))} call cell · relative {_h(_pct(value, sample_call_total))} / sample call cells · total {_h(_pct(value, global_call_total))} / 7-dataset call cells</span>'
            for call_class, value in funnel["call_class_counts"].items()
        )
        role = "HCC1395 的 Dorado pipeline dataset；非獨立生物樣本" if dataset == "HCC1395_DORADO" else "pipeline dataset"
        rows.append(
            "<tr>"
            f'<th scope="row">{_h(dataset)}<span class="denom">Biological ID：{_h(funnel["biological_id"])}<br>{_h(role)}</span></th>'
            f'<td>{_metric_cell(funnel["raw_sparse_molecules"], "molecule", aggregate["raw_sparse_molecules"], "7-dataset raw sparse molecules", funnel["raw_sparse_molecules"], "該 dataset raw sparse molecules")}</td>'
            f'<td>{_metric_cell(funnel["ps_known_molecules"], "molecule", aggregate["ps_known_molecules"], "7-dataset known-PS molecules", funnel["raw_sparse_molecules"], "該 dataset raw sparse molecules")}{_metric_cell(funnel["ps_missing_molecules"], "molecule", aggregate["ps_missing_molecules"], "7-dataset missing-PS molecules", funnel["raw_sparse_molecules"], "該 dataset raw sparse molecules")}</td>'
            f'<td>{_metric_cell(funnel["hp_included_molecules"], "molecule", aggregate["hp_included_molecules"], "7-dataset molecules entering ≥1 PS-aware primary unit", funnel["raw_sparse_molecules"], "該 dataset raw sparse molecules")}{_metric_cell(funnel["hp_excluded_molecules"], "molecule", aggregate["hp_excluded_molecules"], "7-dataset molecules entering no PS-aware primary unit", funnel["raw_sparse_molecules"], "該 dataset raw sparse molecules")}</td>'
            f'<td><span class="metric-main">{_h(_n(sample_call_total))} <small>call cell</small></span><span class="denom">總分母：7-dataset call cells {_h(_n(global_call_total))}（{_h(_pct(sample_call_total, global_call_total))}）</span>{call_lines}</td>'
            "</tr>"
        )
    return f"""
  <h3>7 datasets / 6 biological samples：raw molecule、PS、primary-unit inclusion 與 call-class 漏斗</h3>
  <p class="visual-note">這是 7 個 pipeline datasets 的工程比較；HCC1395 與 HCC1395_DORADO 映射到同一 biological ID。PS known/missing 指 raw sparse row 是否帶任意非缺失 PS tag；「進入 primary unit」則還同時要求 HP1/2、known PS 且觸及至少一個所選 PS-aware component。兩者是平行守恆切面，不可串成單一路徑百分比。JSON 為相容舊介面仍使用 <code>hp_included_molecules</code> 欄名。</p>
  {chart}
  <details open><summary>各 dataset 的 raw sparse molecules、PS、primary-unit inclusion、R/A/O/D/S/L/X</summary><div class="table-wrap"><table><thead><tr><th>Dataset / biological ID</th><th>Raw sparse molecules</th><th>Any PS tag known / missing</th><th>進入 ≥1 PS-aware primary unit / 未進入</th><th>R/A/O/D/S/L/X call cells</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></details>
"""


def _distribution_lines(
    distribution: Mapping[str, int], denominator: int, *, prefix: str = ""
) -> str:
    return "".join(
        f'<span class="denom">{_h(prefix)}{_h(key)}：{_h(_n(value))} · '
        f'{_h(_pct(value, denominator))} / 該 grain</span>'
        for key, value in sorted(distribution.items(), key=lambda item: int(item[0]))
    )


def _ranking_hp_summary_html(
    hp_cells: Mapping[str, Mapping[str, Any]], combined: Mapping[str, Any], threshold: str
) -> str:
    combined_units = combined["n_component_hp_units"]
    rows: list[str] = []
    for label, cell in (*hp_cells.items(), ("HP1+HP2 display sum", combined)):
        units = cell["n_component_hp_units"]
        rows.append(
            "<tr>"
            f'<th scope="row">{_h(label)}<span class="denom">bridge ≥{_h(threshold)} display convention</span></th>'
            f'<td>{_metric_cell(units, "HP/PS component unit", combined_units, "HP1+HP2 display-sum units", units, "該 HP stratum units")}</td>'
            f'<td>{_metric_cell(cell["solver_complete_units"], "unit", combined_units, "HP1+HP2 display-sum units", units, "該 HP stratum units")}{_metric_cell(cell["solver_incomplete_or_not_run_units"], "unit", combined_units, "HP1+HP2 display-sum units", units, "該 HP stratum units")}</td>'
            f'<td>{_metric_cell(cell["quality_primary_unique_vertex_units"], "unit", combined_units, "HP1+HP2 display-sum units", units, "該 HP stratum units")}{_metric_cell(cell["quality_primary_tied_vertex_units"], "unit", combined_units, "HP1+HP2 display-sum units", units, "該 HP stratum units")}{_metric_cell(cell["rank_abstain_units"], "unit", combined_units, "HP1+HP2 display-sum units", units, "該 HP stratum units")}</td>'
            f'<td>{_count_cell(cell["raw_tree_candidates_T_complete_units"], "T", "solver-complete units 的 parent-edge candidate sum；不是 unit 分母")}{_count_cell(cell["distinct_vertex_sets_V_complete_units"], "V", "solver-complete units 的 distinct minimum-h* vertex-set sum；不是 T 的比例")}</td>'
            f'<td>{_metric_cell(cell["topology_evaluated_units"], "unit", combined_units, "HP1+HP2 display-sum units", units, "該 HP stratum units")}{_metric_cell(cell["coarse_topology_unique_units"], "unit", combined["topology_evaluated_units"], "HP1+HP2 topology-evaluated units", cell["topology_evaluated_units"], "該 HP topology-evaluated units")}{_metric_cell(cell["coarse_topology_ambiguous_units"], "unit", combined["topology_evaluated_units"], "HP1+HP2 topology-evaluated units", cell["topology_evaluated_units"], "該 HP topology-evaluated units")}</td>'
            "</tr>"
        )
    return (
        '<details open><summary>HP1 與 HP2 分開：solver、likelihood、T/V 與 topology 分母</summary>'
        '<div class="table-wrap"><table><thead><tr><th>HP stratum</th><th>Primary units</th>'
        '<th>Solver complete / incomplete</th><th>Likelihood unique / tied / abstain</th>'
        '<th>T / V（絕對數）</th><th>Topology evaluated / coarse unique / ambiguous</th></tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table></div></details>'
    )


def _structural_scoring_funnel_html(
    hp_cells: Mapping[str, Mapping[str, Any]], combined: Mapping[str, Any], threshold: str
) -> str:
    rows: list[str] = []
    for label, cell in (*hp_cells.items(), ("HP1+HP2 display sum", combined)):
        projections = cell["molecule_component_projections"]
        informative = cell["informative_scoring_molecules"]
        rows.append(
            "<tr>"
            f'<th scope="row">{_h(label)}<span class="denom">bridge ≥{_h(threshold)}；單位是 molecule projection</span></th>'
            f'<td>{_metric_cell(projections, "projection", combined["molecule_component_projections"], "HP1+HP2 projection sum", projections, "該 HP projection 起點")}</td>'
            f'<td>{_metric_cell(cell["all_x_excluded_molecules"], "projection", combined["molecule_component_projections"], "HP1+HP2 projection sum", projections, "該 HP 全部 projections")}</td>'
            f'<td>{_metric_cell(informative, "projection", combined["informative_scoring_molecules"], "HP1+HP2 informative projections", projections, "該 HP 全部 projections")}</td>'
            f'<td>{_metric_cell(cell["structural_retained_molecules"], "projection", combined["informative_scoring_molecules"], "HP1+HP2 informative projections", informative, "該 HP informative projections")}</td>'
            f'<td>{_metric_cell(cell["below_minread_scoring_molecules"], "projection", combined["informative_scoring_molecules"], "HP1+HP2 informative projections", informative, "該 HP informative projections")}</td>'
            f'<td>{_metric_cell(cell["bq_scoring_molecules"], "projection", combined["bq_scoring_molecules"], "HP1+HP2 BQ-scored projections", informative, "該 HP informative projections")}</td>'
            "</tr>"
        )
    return f"""
  <h3>Structural evidence 與 scoring evidence 分開守恆</h3>
  <p class="visual-note"><strong>同一 molecule 可投影到不同 threshold/component，所以這些不是全基因組 unique read 數。</strong> All-X 無 R/A 資訊而排除；informative 再互斥分成「達 structural minread」與「低於門檻但仍進 likelihood scoring」。BQ scoring 應等於全部 informative projections。</p>
  <div class="table-wrap"><table><thead><tr><th>HP stratum</th><th>All projections</th><th>All-X excluded</th><th>Informative scoring</th><th>Structural retained</th><th>Below-minread but scored</th><th>BQ-scored</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
"""


def _solver_effective_k_html(
    hp_cells: Mapping[str, Mapping[str, Any]], combined: Mapping[str, Any], threshold: str
) -> str:
    rows: list[str] = []
    for label, cell in (*hp_cells.items(), ("HP1+HP2 display sum", combined)):
        units = cell["n_component_hp_units"]
        component_k = cell["k_component_sites_total"]
        effective_k = cell["k_observed_alt_active_total"]
        scoring_k = cell["k_scoring_alt_observed_total"]
        route_lines = "".join(
            f'<span class="denom"><code>{_h(route)}</code> {_h(_n(value))} · {_h(_pct(value, units))} / 該 HP units</span>'
            for route, value in cell["k_route_counts"].items()
        )
        rows.append(
            "<tr>"
            f'<th scope="row">{_h(label)}<span class="denom">bridge ≥{_h(threshold)}</span></th>'
            f'<td>{_metric_cell(cell["solver_complete_units"], "unit", combined["n_component_hp_units"], "HP1+HP2 units", units, "該 HP units")}{_metric_cell(cell["solver_incomplete_or_not_run_units"], "unit", combined["n_component_hp_units"], "HP1+HP2 units", units, "該 HP units")}</td>'
            f'<td>{route_lines}</td>'
            f'<td>{_count_cell(component_k, "site-coordinate projection", f"sum k_component_sites；mean {component_k / units if units else 0:.3f} / unit")}</td>'
            f'<td>{_count_cell(effective_k, "active ALT-coordinate projection", f"sum effective k；mean {effective_k / units if units else 0:.3f} / unit；{_pct(effective_k, component_k)} / component-site mass")}</td>'
            f'<td>{_count_cell(scoring_k, "observed-ALT coordinate projection", f"sum scoring ALT k；mean {scoring_k / units if units else 0:.3f} / unit")}</td>'
            f'<td>{_metric_cell(cell["not_structural_alt_active_sites_total"], "site-coordinate projection", combined["k_component_sites_total"], "HP1+HP2 component-site mass", component_k, "該 HP component-site mass")}</td>'
            "</tr>"
        )
    return f"""
  <h3>Solver complete/incomplete 與 effective-k 路由</h3>
  <p class="visual-note"><code>k_component_sites</code> 是 component 全部座標；<code>k_observed_alt_active</code> 才是 structural minread 下進 exact solver 的 effective k。下表的 k totals 是跨 unit 的 coordinate projections 加總，不是 unique genomic sites。<code>GT_EXACT_LIMIT</code> 才表示 effective k&gt;12；component-site k&gt;12 不等於 solver 必然 abstain。</p>
  <div class="table-wrap"><table><thead><tr><th>HP stratum</th><th>Complete / incomplete</th><th>Effective-k route</th><th>Σ component k</th><th>Σ structural effective k</th><th>Σ scoring ALT k</th><th>Σ non-structural coordinates</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
"""


def _partial_completion_funnel_html(
    hp_cells: Mapping[str, Mapping[str, Any]], combined: Mapping[str, Any]
) -> str:
    grain_labels = {
        "unique_rax_pattern_groups": "Unique R/A/X patterns",
        "bq_quality_pattern_groups": "BQ-quality pattern groups",
        "molecule_projections": "Molecule projections",
        "structural_unique_rax_pattern_groups": "Structural unique R/A/X patterns",
    }
    rows: list[str] = []
    combined_funnel = combined["partial_pattern_funnel"]
    for hp_label, cell in (*hp_cells.items(), ("HP1+HP2 display sum", combined)):
        funnel = cell["partial_pattern_funnel"]
        for grain_name in PARTIAL_FUNNEL_GRAINS:
            grain = funnel[grain_name]
            total_grain = combined_funnel[grain_name]["denominator"]
            denominator = grain["denominator"]
            conceptual_total = grain["conceptual_completions_weighted_total"]
            effective_total = grain["observed_alt_effective_completions_weighted_total"]
            rows.append(
                "<tr>"
                f'<th scope="row">{_h(hp_label)}<span class="denom">{_h(grain_labels[grain_name])}</span></th>'
                f'<td>{_metric_cell(denominator, "observation", total_grain, "HP1+HP2 該 grain observations", denominator, "該 HP/grain observations")}</td>'
                f'<td>{_metric_cell(grain["full"], "observation", total_grain, "HP1+HP2 該 grain observations", denominator, "該 HP/grain observations")}{_metric_cell(grain["partial"], "observation", total_grain, "HP1+HP2 該 grain observations", denominator, "該 HP/grain observations")}</td>'
                f'<td>{_distribution_lines(grain["u_number_of_X_distribution"], denominator, prefix="u=")}</td>'
                f'<td>{_distribution_lines(grain["conceptual_completions_2_pow_u_distribution"], denominator, prefix="2^u=")}<span class="metric-main">{_h(_n(conceptual_total))} <small>conceptual completion mass</small></span></td>'
                f'<td>{_distribution_lines(grain["observed_alt_effective_completions_distribution"], denominator, prefix="effective=")}<span class="metric-main">{_h(_n(effective_total))} <small>effective completion mass</small></span><span class="denom">effective / conceptual weighted mass：{_h(_pct(effective_total, conceptual_total))}</span></td>'
                f'<td>{_metric_cell(grain["observed_alt_effective_zero_due_to_fixed_alt_outside_structural_universe"], "observation", total_grain, "HP1+HP2 該 grain observations", denominator, "該 HP/grain observations")}</td>'
                "</tr>"
            )
    coverage = combined_funnel
    return f"""
  <h3>Partial-read conceptual 與 effective completions：四種 grain 全數列出</h3>
  <p class="visual-note"><strong>Conceptual</strong> 是單一 pattern 的 full cube <code>2^u</code>；<strong>effective</strong> 是投影到 minread-specific observed-ALT universe 後實際可相容的 active states。Weighted mass 是「每筆 observation 的 completion 數」加總，不是 tree 數、clone 數或獨立機率分母。</p>
  <p class="visual-note">Combined partial structural coverage：{_h(_n(coverage["partial_groups_covered"]))} / {_h(_n(coverage["partial_group_coverage_denominator"]))} groups covered（{_h(_pct(coverage["partial_groups_covered"], coverage["partial_group_coverage_denominator"]))}）；unsatisfied {_h(_n(coverage["partial_groups_unsatisfied"]))}。Units with partial structural groups：{_h(_n(coverage["units_with_partial_structural_groups"]))} / {_h(_n(coverage["units_denominator"]))}（{_h(_pct(coverage["units_with_partial_structural_groups"], coverage["units_denominator"]))}）。</p>
  <details open><summary>HP1 / HP2 × 四種 observation grain 的完整 completion 數值</summary><div class="table-wrap"><table><thead><tr><th>HP / observation grain</th><th>Denominator</th><th>Full / partial</th><th>u distribution</th><th>Conceptual 2^u distribution / weighted total</th><th>Effective distribution / weighted total</th><th>Effective zero</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></details>
"""


def _topology_identifiability_html(
    hp_cells: Mapping[str, Mapping[str, Any]], combined: Mapping[str, Any]
) -> str:
    rows: list[str] = []
    for label, cell in (*hp_cells.items(), ("HP1+HP2 display sum", combined)):
        evaluated = cell["topology_evaluated_units"]
        parent_unique = cell["parent_edge_assignment_unique_units"]
        multi_parent = evaluated - parent_unique
        exact_lines = "".join(
            f'<span class="denom"><code>{_h(status)}</code> {_h(_n(value))} · {_h(_pct(value, cell["n_component_hp_units"]))} / 該 HP units</span>'
            for status, value in cell["exact_topology_uniqueness_status_counts"].items()
        )
        rows.append(
            "<tr>"
            f'<th scope="row">{_h(label)}</th>'
            f'<td>{_metric_cell(evaluated, "unit", combined["n_component_hp_units"], "HP1+HP2 units", cell["n_component_hp_units"], "該 HP units")}</td>'
            f'<td>{_metric_cell(parent_unique, "unit", combined["topology_evaluated_units"], "HP1+HP2 topology-evaluated units", evaluated, "該 HP topology-evaluated units")}</td>'
            f'<td>{_metric_cell(multi_parent, "unit", combined["topology_evaluated_units"], "HP1+HP2 topology-evaluated units", evaluated, "該 HP topology-evaluated units")}</td>'
            f'<td>{_metric_cell(cell["exact_topology_proven_unique_units"], "unit", combined["topology_evaluated_units"], "HP1+HP2 topology-evaluated units", evaluated, "該 HP topology-evaluated units")}</td>'
            f'<td>{exact_lines}</td>'
            "</tr>"
        )
    return f"""
  <h3>Exact topology 與 parent-choice uniqueness 的可識別邊界</h3>
  <p class="visual-note"><strong>Exact-topology proven unique 目前只在 winning candidates 合計剛好一組 labeled parent-edge assignment 時成立。</strong>若有多組 parent choices，receipt 標記的是 <code>NOT_EVALUATED_CANONICAL_SHAPE_ISOMORPHISM...</code>，這不等於已證明 exact topology 不唯一。Coarse class 唯一也不可替代 exact labeled topology 唯一。</p>
  <div class="table-wrap"><table><thead><tr><th>HP stratum</th><th>Topology evaluated</th><th>Single parent assignment</th><th>Multiple parent assignments</th><th>Exact topology proven unique</th><th>Exact-topology status partition</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
"""


def _runtime_and_resource_html(
    ranking_source: Source,
    verification_source: Source | None,
    resource_source: Source | None,
) -> str:
    runtime = ranking_source.payload["runtime_diagnostics"]
    runtime_rows: list[str] = []
    metric_labels = {
        "candidate_generation_elapsed_seconds": "Candidate generation",
        "likelihood_fit_elapsed_seconds": "Likelihood fit",
        "unit_total_elapsed_seconds": "Unit total",
    }
    for metric in RUNTIME_METRICS:
        summary = runtime["metrics"][metric]
        runtime_rows.append(
            "<tr>"
            f'<th scope="row">{_h(metric_labels[metric])}<span class="denom">all primary-unit evaluations</span></th>'
            f'<td>{_h(_n(summary["n"]))}</td><td>{_h(_n(summary["sum"]))}</td>'
            f'<td>{_h(_n(summary["p50"]))}</td><td>{_h(_n(summary["p95"]))}</td>'
            f'<td>{_h(_n(summary["p99"]))}</td><td>{_h(_n(summary["max"]))}</td>'
            "</tr>"
        )
        if metric in runtime["metrics_when_invoked"]:
            invoked = runtime["metrics_when_invoked"][metric]
            runtime_rows.append(
                "<tr>"
                f'<th scope="row">{_h(metric_labels[metric])}<span class="denom">only when invoked</span></th>'
                f'<td>{_h(_n(invoked["n"]))}</td><td>{_h(_n(invoked["sum"]))}</td>'
                f'<td>{_h(_n(invoked["p50"]))}</td><td>{_h(_n(invoked["p95"]))}</td>'
                f'<td>{_h(_n(invoked["p99"]))}</td><td>{_h(_n(invoked["max"]))}</td>'
                "</tr>"
            )
    independent_note = (
        f'獨立 verifier 已由 154 child per-unit TSV 重算並精確吻合。{_source_ref("S13")}'
        if verification_source is not None else "未提供獨立 runtime verifier。"
    )
    resource_html = '<p class="visual-note">無可用 resource-conflict session attestation。</p>'
    if verification_source is not None and resource_source is not None:
        gate = resource_source.payload["resource_gate"]
        orchestration = verification_source.payload["extraction"]["orchestration"]
        n_batches = orchestration["n_batches"]
        resource_html = f"""
  <div class="metric-ledger" aria-label="Resource conflict gate attestations">
    <div>{_count_cell(gate["process_count"], "conflicting process", "persisted extraction-session gate observation")}</div>
    <div>{_count_cell(gate["root_count"], "conflicting process root", "persisted extraction-session gate observation")}</div>
    <div>{_count_cell(1 + n_batches, "zero-conflict launch attestation", f"1 session + {n_batches} extraction batches；各批由 independent verifier 驗證")}</div>
    <div>{_count_cell(0, "resource-gate bypass", "ignore_resource_gate=false")}</div>
  </div>
  <p class="visual-note">Session observation time：{_h(gate["observed_at_utc"])}；process snapshot SHA-256：<code>{_h(gate["process_snapshot_sha256"])}</code>。Ranking 沒有另啟一個 resource gate，而是綁定已驗證的 extraction parent session。{_source_ref("S15")}{_source_ref("S13")}</p>
"""
    return f"""
  <h3>Ranking runtime 與 resource-conflict attestation</h3>
  <p class="visual-note">Runtime 是 <code>time.monotonic_ns</code> 的 process-local wall time，使用 exact empirical nearest-rank p50/p95/p99；{independent_note}這是當次環境/負載下的效能診斷，不是科學結果或跨機器 reproducibility claim。</p>
  <div class="table-wrap"><table><thead><tr><th>Runtime segment</th><th>n evaluations</th><th>Σ seconds</th><th>p50 s</th><th>p95 s</th><th>p99 s</th><th>max s</th></tr></thead><tbody>{''.join(runtime_rows)}</tbody></table></div>
  {resource_html}
  <p class="visual-note"><strong>Resource-conflict gate 只證明啟動觀察時沒有已知衝突作業；不是 peak RSS、CPU-hours、I/O 或 disk-footprint telemetry。</strong>本報告不會把「0 conflicts」誤寫為「0 resource usage」。</p>
"""


def _numeric_summary_per_dataset_html(source: Source | None) -> str:
    if source is None:
        return '<p class="visual-note">尚未提供 sidecar-authenticated final numeric summary；per-dataset all-threshold ledger 保持未評估。</p>'
    payload = source.payload
    extraction = payload["extraction"]
    ranking = payload["ranking"]
    overall_components = extraction["overall_component_by_linkage_basis_threshold"]
    extraction_by_dataset = extraction["by_dataset"]
    overall_ranking = ranking["overall_by_HP_basis_and_bridge_threshold"]
    ranking_by_dataset = ranking["by_dataset"]
    component_rows: list[str] = []
    structure_rows: list[str] = []
    route_order = ("EXACT_COMPLETE", "EXACT_INCOMPLETE", "NOT_RUN_NO_STRUCTURE", "GT_EXACT_LIMIT")

    def scope_header(dataset: str | None) -> str:
        if dataset is None:
            return '<th scope="row"><strong>OVERALL · 7 datasets</strong><span class="denom">6 biological samples；HCC1395 pipelines 未去重</span></th>'
        return (
            f'<th scope="row">{_h(dataset)}'
            f'<span class="denom">Biological ID：{_h(BIOLOGICAL_IDS[dataset])}</span></th>'
        )

    def topology_mapping_lines(
        mapping: Mapping[str, int], total_mapping: Mapping[str, int],
        relative_denominator: int, total_denominator: int,
        *, include_all_classes: bool,
    ) -> str:
        keys = (
            tuple(TOPOLOGY_CLASS_ORDER)
            if include_all_classes
            else tuple(sorted(set(mapping) | set(total_mapping)))
        )
        if not keys:
            return _count_cell(0, "ambiguous class-set unit", "此 stratum 無 ambiguous topology class set")
        return "".join(
            _metric_cell(
                int(mapping.get(key, 0)),
                f"{key} unit",
                total_denominator,
                "同 HP×threshold 的 7-dataset topology denominator",
                relative_denominator,
                "該 stratum topology denominator",
            )
            for key in keys
        )

    def t_or_v_cell(
        value: int, total: int, unit_denominator: int, symbol: str,
    ) -> str:
        total_pct = _pct(value, total)
        mean = "不適用" if unit_denominator == 0 else f"{value / unit_denominator:.4g}"
        return (
            f'<span class="metric-main">{_n(value)} <small>Σ{_h(symbol)}</small></span>'
            f'<span class="denom">總比例：{_h(total_pct)}（分母：同 HP×threshold 的 7-dataset Σ{_h(symbol)}={_n(total)}）</span>'
            f'<span class="denom">相對量：{_h(mean)} {_h(symbol)}/solver-complete unit（分母={_n(unit_denominator)}）</span>'
        )

    # Eight aggregate rows are first-class rows, not a separately hand-copied KPI.
    scopes: tuple[str | None, ...] = (None, *DATASETS)
    for dataset in scopes:
        for basis in ("PS_HP1", "PS_HP2"):
            for threshold in ("1", "2", "3", "5"):
                component = (
                    overall_components[basis][threshold]
                    if dataset is None
                    else extraction_by_dataset[dataset]["component_by_linkage_basis_threshold"][basis][threshold]
                )
                component_total = overall_components[basis][threshold]
                n_components = component["n_components"]
                total_components = component_total["n_components"]
                component_rows.append(
                    "<tr>"
                    f'{scope_header(dataset)}'
                    f'<td><code>{_h(basis)}</code></td><td>≥{_h(threshold)}</td>'
                    f'<td>{_metric_cell(n_components, "component", total_components, f"7-dataset {basis} bridge≥{threshold} components", n_components, "該 stratum components")}</td>'
                    f'<td>{_metric_cell(component["k_equals_1"], "component", total_components, f"7-dataset {basis} bridge≥{threshold} components", n_components, "該 stratum components")}</td>'
                    f'<td>{_metric_cell(component["k_greater_than_1"], "component", total_components, f"7-dataset {basis} bridge≥{threshold} components", n_components, "該 stratum components")}</td>'
                    f'<td>{_count_cell(component["max_k_component_sites"], "site k", "該 stratum 最大 component-site k；不是 effective k")}</td>'
                    "</tr>"
                )

                rank_cell = (
                    overall_ranking[basis][threshold]
                    if dataset is None
                    else ranking_by_dataset[dataset]["by_HP_basis_and_bridge_threshold"][basis][threshold]
                )
                rank_total = overall_ranking[basis][threshold]
                units = rank_cell["units"]["n_component_hp_unit_evaluations"]
                total_units = rank_total["units"]["n_component_hp_unit_evaluations"]
                route_cells = "".join(
                    _metric_cell(
                        rank_cell["effective_k"]["k_route_counts"][route],
                        f"{route} unit",
                        total_units,
                        f"7-dataset {basis} bridge≥{threshold} HP/PS units",
                        units,
                        "該 stratum units",
                    )
                    for route in route_order
                )
                candidate = rank_cell["candidate_structure"]["candidate_table"]
                candidate_total = rank_total["candidate_structure"]["candidate_table"]
                candidate_units = candidate["n_solver_complete_candidate_units"]
                total_candidate_units = candidate_total["n_solver_complete_candidate_units"]
                tied_denominator = candidate["tied_by_coarse_topology"]["denominator"]
                total_tied_denominator = candidate_total["tied_by_coarse_topology"]["denominator"]
                h_lines = "".join(
                    _metric_cell(
                        count,
                        f"h*={h_star} unit",
                        total_candidate_units,
                        f"7-dataset {basis} bridge≥{threshold} solver-complete candidate units",
                        candidate_units,
                        "該 stratum solver-complete candidate units",
                    )
                    for h_star, count in sorted(candidate["h_star_distribution"].items(), key=lambda item: int(item[0]))
                )
                tree_partition = candidate["tree_vertex_partition"]
                total_tree_partition = candidate_total["tree_vertex_partition"]
                partition_lines = "".join(
                    _metric_cell(
                        tree_partition["counts"][bucket],
                        f"{bucket} unit",
                        total_tree_partition["denominator"],
                        f"7-dataset {basis} bridge≥{threshold} solver-complete units",
                        tree_partition["denominator"],
                        "該 stratum solver-complete units",
                    )
                    for bucket in TREE_VERTEX_BUCKETS
                )
                topology = candidate["topology"]
                topology_total = candidate_total["topology"]
                unique_class_lines = topology_mapping_lines(
                    topology["coarse_unique_class_counts"],
                    topology_total["coarse_unique_class_counts"],
                    topology["coarse_class_unique_units"],
                    topology_total["coarse_class_unique_units"],
                    include_all_classes=True,
                )
                ambiguous_lines = topology_mapping_lines(
                    topology["coarse_ambiguous_class_set_counts"],
                    topology_total["coarse_ambiguous_class_set_counts"],
                    topology["coarse_class_multiple_units"],
                    topology_total["coarse_class_multiple_units"],
                    include_all_classes=False,
                )
                structure_rows.append(
                    "<tr>"
                    f'{scope_header(dataset)}'
                    f'<td><code>{_h(basis)}</code></td><td>≥{_h(threshold)}</td>'
                    f'<td>{_metric_cell(units, "HP/PS component unit", total_units, f"7-dataset {basis} bridge≥{threshold} units", units, "該 stratum units")}</td>'
                    f'<td>{route_cells}</td>'
                    f'<td>{t_or_v_cell(candidate["n_parent_edge_trees_T"], candidate_total["n_parent_edge_trees_T"], candidate_units, "T")}</td>'
                    f'<td>{t_or_v_cell(candidate["n_candidate_vertex_sets_V"], candidate_total["n_candidate_vertex_sets_V"], candidate_units, "V")}</td>'
                    f'<td>{partition_lines}</td>'
                    f'<td>{h_lines}</td>'
                    f'<td>{_metric_cell(candidate["unique_first"], "candidate unit", total_candidate_units, f"7-dataset {basis} bridge≥{threshold} candidate units", candidate_units, "該 stratum candidate units")}</td>'
                    f'<td>{_metric_cell(candidate["tied_first"], "candidate unit", total_candidate_units, f"7-dataset {basis} bridge≥{threshold} candidate units", candidate_units, "該 stratum candidate units")}</td>'
                    f'<td>{_metric_cell(rank_cell["ranking_outcome"]["abstain_all_causes"], "all-cause abstain unit", total_units, f"7-dataset {basis} bridge≥{threshold} units", units, "該 stratum units")}{_metric_cell(candidate["solver_complete_optimizer_abstain"], "optimizer abstain unit", total_candidate_units, f"7-dataset {basis} bridge≥{threshold} candidate units", candidate_units, "該 stratum solver-complete units")}</td>'
                    f'<td>{unique_class_lines}</td><td>{ambiguous_lines}</td>'
                    f'<td>{_metric_cell(topology["parent_edge_assignment_unique_units"], "parent-edge unique unit", topology_total["evaluated_units"], f"7-dataset {basis} bridge≥{threshold} topology-evaluated units", topology["evaluated_units"], "該 stratum topology-evaluated units")}{_metric_cell(topology["exact_topology_proven_unique_units"], "exact-topology proven unique unit", topology_total["evaluated_units"], f"7-dataset {basis} bridge≥{threshold} topology-evaluated units", topology["evaluated_units"], "該 stratum topology-evaluated units")}</td>'
                    f'<td>{_metric_cell(candidate["tied_by_coarse_topology"]["consistent"], "tied Topo-consistent unit", total_tied_denominator, f"7-dataset {basis} bridge≥{threshold} tied-first units", tied_denominator, "該 stratum tied-first units")}{_metric_cell(candidate["tied_by_coarse_topology"]["inconsistent"], "tied Topo-inconsistent unit", total_tied_denominator, f"7-dataset {basis} bridge≥{threshold} tied-first units", tied_denominator, "該 stratum tied-first units")}</td>'
                    "</tr>"
                )
    unsupported_labels = {
        "deduplicated_biological_regions_across_thresholds": "跨 thresholds 去重 biological regions",
        "cellular_HP1_HP2_clone_pairing": "Cell-level HP1/HP2 clone pairing",
        "exact_parent_edge_topology_for_tied_vertex_sets": "Tied vertex sets 的 exact parent-edge topology",
        "global_most_likely_structure_across_h_star_plus_1_or_more": "跨 h*+1 以上複雜度的 global-most-likely structure",
        "independent_VAF_confirmation_term": "同 read 的獨立 VAF confirmation term",
        "formal_full_run_peak_RSS_CPU_and_disk_envelope": "Formal full-run peak RSS/CPU/disk envelope",
    }
    unsupported_rows = "".join(
        f'<tr><th scope="row">{_h(unsupported_labels[key])}</th><td><strong>未估計（null）</strong></td><td>{_h(payload["unsupported_or_nonidentifiable"][key]["reason"])}</td></tr>'
        for key in sorted(NUMERIC_SUMMARY_UNSUPPORTED_KEYS)
    )
    return f"""
  <h3>Authenticated final numeric summary：每 dataset × HP × bridge 的完整數字</h3>
  <p class="takeaway"><strong>以下不是把 thresholds 或 HP1/HP2 加總成 biological regions。</strong>8 個 OVERALL rows 是 2 HP×4 thresholds 的 7-dataset 技術母體；其後 56 rows 是 dataset×HP×threshold。所有值都由 S16 逐格對 S6/S7 raw cells 與 S9 canonical candidate stream 重算後才顯示。「總分母」固定為同 HP×threshold 的 7-dataset 母體，「相對分母」固定為該 stratum 的同類單位。{_source_ref("S16")}{_source_ref("S9")}</p>
  <details open><summary>Extraction component 與 k=1/k&gt;1（8 overall + 56 dataset strata）</summary><div class="table-wrap"><table><thead><tr><th>Scope</th><th>HP basis</th><th>Bridge</th><th>Components</th><th>k=1</th><th>k&gt;1</th><th>Max k</th></tr></thead><tbody>{''.join(component_rows)}</tbody></table></div></details>
  <details open><summary>完整 T/V、effective route、h*、ranking 與 topology ledger（8 overall + 56 dataset strata）</summary><div class="table-wrap"><table><thead><tr><th>Scope</th><th>HP basis</th><th>Bridge</th><th>Units</th><th>Effective-k route（每值有 inline label）</th><th>ΣT</th><th>ΣV</th><th>T/V 三分桶</th><th>h* distribution</th><th>Unique first</th><th>Tied first</th><th>Abstain</th><th>Coarse unique class counts（完整四類）</th><th>Ambiguous class-set counts（完整 observed sets）</th><th>Parent-edge / exact unique</th><th>Tied coarse-Topo</th></tr></thead><tbody>{''.join(structure_rows)}</tbody></table></div></details>
  <details open><summary>目前不可由正式證據辨識的量：null + reason</summary><div class="table-wrap"><table><thead><tr><th>Metric</th><th>Value</th><th>Reason</th></tr></thead><tbody>{unsupported_rows}</tbody></table></div></details>
"""


def _m2_ranking_section(
    source: Source | None,
    candidate_source: Source | None = None,
    verification_source: Source | None = None,
    resource_source: Source | None = None,
    numeric_source: Source | None = None,
) -> str:
    if source is None:
        return """
<section id="m2-ranking" class="report-section unavailable"><p class="section-index">04 / M2 SYMBOLIC + LIKELIHOOD</p><h2>M2 PS-aware full ranking receipt 尚未完成或未提供</h2><p>尚不能填寫 V/T、likelihood unique/tie/abstain、solver complete/incomplete 或固定候選集條件下的 BQ/error-grid ranking stability；任何此類空缺均不是 0。</p></section>
"""
    payload = source.payload
    if payload.get("schema_version") != "2.0.0":
        return """
<section id="m2-ranking" class="report-section unavailable"><p class="section-index">04 / M2 SYMBOLIC + LIKELIHOOD</p><h2>舊 ranking schema 僅作 diagnostic；數字刻意不升級為 current result</h2><p>schema 1.0.0 可能把不同 PS block 的 patterns 混在同一 unit。終版要求 child/full schema 2.0.0、primary unit＝<code>HP_family×known_PS×read_linked_component×threshold</code>，missing PS 另列 diagnostic；在此之前不展示 likelihood unique/tie/abstain 聚合數字。</p></section>
"""
    basis, threshold, cell = _primary_ranking_cell(payload)
    hp_cells = {
        family: payload["aggregate"]["by_linkage_basis_and_threshold"][family][threshold]
        for family in ("PS_HP1", "PS_HP2")
    }
    statuses = sorted(cell["selection_status_counts"].items(), key=lambda item: (-item[1], item[0]))
    chart = _bar_svg(statuses, title=f"M2 likelihood status · {basis} display sum · bridge ≥{threshold}", description="Selection status over the HP1+HP2 display sum at the professor-view bridge convention; HP-specific denominators are reported separately in the table.", unit="component/HP unit")
    topology_unique_chart = _bar_svg(
        sorted(cell["coarse_topology_unique_class_counts"].items(), key=lambda item: (-item[1], item[0])),
        title=f"Mutually exclusive coarse topology classes · {basis} · bridge ≥{threshold}",
        description="Only units with one unique coarse topology class are included. The denominator is coarse_topology_unique_units, not all topology-evaluated units.",
        unit="coarse-unique unit",
        color="#171511",
    )
    combined_partial = cell["partial_pattern_funnel"]
    partial_chart_distributions = {
        "unique_patterns": combined_partial["unique_rax_pattern_groups"]["u_number_of_X_distribution"],
        "quality_groups": combined_partial["bq_quality_pattern_groups"]["u_number_of_X_distribution"],
        "molecule_projections": combined_partial["molecule_projections"]["u_number_of_X_distribution"],
        "structural_patterns": combined_partial["structural_unique_rax_pattern_groups"]["u_number_of_X_distribution"],
    }
    partial_chart_denominators = {
        "unique_patterns": combined_partial["unique_rax_pattern_groups"]["denominator"],
        "quality_groups": combined_partial["bq_quality_pattern_groups"]["denominator"],
        "molecule_projections": combined_partial["molecule_projections"]["denominator"],
        "structural_patterns": combined_partial["structural_unique_rax_pattern_groups"]["denominator"],
    }
    partial_chart = _partial_u_heatmap_svg(partial_chart_distributions, partial_chart_denominators)
    status_rows = "".join(
        f'<tr><th scope="row">{_h(key)}</th><td>{_metric_cell(value, "component/HP unit", cell["n_component_hp_units"], "primary cell HP units", cell["n_component_hp_units"], "selection-status 母體")}</td></tr>'
        for key, value in statuses
    )
    strata_rows = []
    for stratum_basis, thresholds in payload["aggregate"]["by_linkage_basis_and_threshold"].items():
        for stratum_threshold, value in sorted(thresholds.items(), key=lambda item: int(item[0])):
            strata_rows.append(
                "<tr>"
                f'<th scope="row">{_h(stratum_basis)} / ≥{_h(stratum_threshold)}</th>'
                f'<td>{_count_cell(value["n_components"], "component", "絕對數；component 不是 percentage category")}</td>'
                f'<td>{_count_cell(value["n_component_hp_units"], "HP unit", "絕對數；同一 component 可有多個 HP/PS units")}</td>'
                f'<td>{_count_cell(value["raw_tree_candidates_T_complete_units"], "tree candidate T", "只在 complete units 加總；T 不是分母")}</td>'
                f'<td>{_count_cell(value["distinct_vertex_sets_V_complete_units"], "vertex set V", "只在 complete units 加總；V 不是 T 的百分比")}</td>'
                "</tr>"
            )
    topology_unique_rows = "".join(
        f'<tr><th scope="row">{_h(key)}</th><td>{_metric_cell(value, "coarse-unique unit", cell["coarse_topology_unique_units"], "coarse-topology unique units", cell["topology_evaluated_units"], "topology-evaluated units")}</td></tr>'
        for key, value in sorted(cell["coarse_topology_unique_class_counts"].items(), key=lambda item: (-item[1], item[0]))
    )
    topology_ambiguous_rows = "".join(
        f'<tr><th scope="row">{_h(key)}</th><td>{_metric_cell(value, "ambiguous unit", cell["topology_evaluated_units"], "topology-evaluated units", cell["coarse_topology_ambiguous_units"], "coarse-topology ambiguous units")}</td></tr>'
        for key, value in sorted(cell["ambiguous_topology_class_set_counts"].items(), key=lambda item: (-item[1], item[0]))
    )
    partial_denominator_rows = "".join(
        f'<tr><th scope="row">{_h(key)}</th><td>{_metric_cell(value, key.replace("_", " "), value, "該觀察粒度總母體", value, "u 分布加總")}</td></tr>'
        for key, value in cell["partial_pattern_denominators"].items()
    )
    dataset_result_rows = []
    for dataset in DATASETS:
        nested = payload["by_dataset"][dataset]["by_linkage_basis_and_threshold"]
        for family in ("PS_HP1", "PS_HP2"):
            dataset_cell = nested[family][threshold]
            unit_denominator = dataset_cell["n_component_hp_units"]
            family_total_units = hp_cells[family]["n_component_hp_units"]
            status_lines = "".join(
                f'<span class="denom"><code>{_h(key)}</code> {_h(_n(value))} · {_h(_pct(value, unit_denominator))} / 該 dataset×HP units</span>'
                for key, value in sorted(dataset_cell["selection_status_counts"].items(), key=lambda item: (-item[1], item[0]))
            )
            topology_denominator = dataset_cell["topology_evaluated_units"]
            unique_denominator = dataset_cell["coarse_topology_unique_units"]
            topology_lines = (
                f'<span class="denom">Topo evaluated {_h(_n(topology_denominator))} · {_h(_pct(topology_denominator, unit_denominator))} / units</span>'
                f'<span class="denom">Coarse unique {_h(_n(unique_denominator))} · {_h(_pct(unique_denominator, topology_denominator))} / evaluated</span>'
                f'<span class="denom">Coarse ambiguous {_h(_n(dataset_cell["coarse_topology_ambiguous_units"]))} · {_h(_pct(dataset_cell["coarse_topology_ambiguous_units"], topology_denominator))} / evaluated</span>'
                f'<span class="denom">Single parent assignment {_h(_n(dataset_cell["parent_edge_assignment_unique_units"]))} · {_h(_pct(dataset_cell["parent_edge_assignment_unique_units"], topology_denominator))} / evaluated</span>'
                + "".join(
                    f'<span class="denom"><code>{_h(key)}</code> {_h(_n(value))} · {_h(_pct(value, unique_denominator))} / coarse-unique</span>'
                    for key, value in sorted(dataset_cell["coarse_topology_unique_class_counts"].items(), key=lambda item: (-item[1], item[0]))
                )
            )
            role = "同一 HCC1395 biological sample 的 Dorado dataset" if dataset == "HCC1395_DORADO" else "pipeline dataset"
            dataset_result_rows.append(
                "<tr>"
                f'<th scope="row">{_h(dataset)}<span class="denom">Biological ID：{_h(BIOLOGICAL_IDS[dataset])}<br>{_h(role)}</span></th>'
                f'<td><code>{_h(family)}</code></td>'
                f'<td>{_metric_cell(unit_denominator, "HP/PS component unit", family_total_units, f"7-dataset {family} units", unit_denominator, "該 dataset×HP units")}</td>'
                f'<td>{_metric_cell(dataset_cell["solver_complete_units"], "unit", family_total_units, f"7-dataset {family} units", unit_denominator, "該 dataset×HP units")}{_metric_cell(dataset_cell["solver_incomplete_or_not_run_units"], "unit", family_total_units, f"7-dataset {family} units", unit_denominator, "該 dataset×HP units")}</td>'
                f'<td>{status_lines}</td>'
                f'<td>{_count_cell(dataset_cell["raw_tree_candidates_T_complete_units"], "T", "solver-complete units 的 parent-edge candidate sum")}{_count_cell(dataset_cell["distinct_vertex_sets_V_complete_units"], "V", "solver-complete units 的 distinct minimum-h* vertex-set sum")}</td>'
                f'<td>{topology_lines}</td>'
                "</tr>"
            )
    verification_html = ""
    if verification_source is not None:
        profile = (
            (verification_source.payload.get("ranking") or {}).get(
                "profile_likelihood_independent_recomputation"
            )
            or {}
        )
        required_profile_keys = {
            "n_children", "n_units", "n_candidates", "n_pattern_rows",
            "n_scoring_molecules", "max_abs_ll_delta",
            "max_abs_relative_weight_delta", "max_abs_gap_delta",
            "max_abs_simplex_residual_delta", "peak_pattern_rows_per_unit",
            "peak_candidates_per_unit", "peak_states_per_candidate",
            "peak_emission_cells_per_candidate",
        }
        if required_profile_keys <= set(profile):
            verification_html = f"""
  <h3>獨立 verifier 已從 pattern、BQ、state 與 π 重算 likelihood</h3>
  <p class="takeaway"><strong>這不是只比對 producer receipt。</strong> Verifier 不 import production ranker，逐一串流 {_h(_n(profile["n_children"]))} 個 dataset×chromosome child，重新計算每個 candidate 的 conditional BQ profile likelihood、simplex residual、Frank–Wolfe/KKT gap、relative weight 與 winner/tie partition。{_source_ref("S13")}</p>
  <div class="metric-ledger" aria-label="Independent profile likelihood recomputation coverage">
    <div>{_count_cell(profile["n_units"], "candidate-bearing primary unit", "HP family × known PS × read-linked component × threshold；與 candidate table units 完全相等")}</div>
    <div>{_count_cell(profile["n_candidates"], "candidate vertex-set fit", "與 canonical candidate table rows 完全相等")}</div>
    <div>{_count_cell(profile["n_pattern_rows"], "R/A/X + fixed-BQ group", "只加總實際 join 到 candidate-bearing units 的 pattern rows")}</div>
    <div>{_count_cell(profile["n_scoring_molecules"], "molecule projection", "跨 unit/threshold 的投影加總；不是全域 unique molecules")}</div>
  </div>
  <div class="table-wrap"><table><thead><tr><th>獨立重算差值／streaming peak</th><th>觀察值</th><th>判讀</th></tr></thead><tbody>
    <tr><th scope="row">max |Δ profile LL|</th><td>{_h(_n(profile["max_abs_ll_delta"]))}</td><td>stored π 下獨立重算與 persisted LL 的最大絕對差</td></tr>
    <tr><th scope="row">max |Δ relative weight|</th><td>{_h(_n(profile["max_abs_relative_weight_delta"]))}</td><td>candidate 正規化後 relative likelihood weight 差</td></tr>
    <tr><th scope="row">max |Δ KKT gap|</th><td>{_h(_n(profile["max_abs_gap_delta"]))}</td><td>同一 concave mixture objective 的 global-gap certificate 差</td></tr>
    <tr><th scope="row">max |Δ simplex residual|</th><td>{_h(_n(profile["max_abs_simplex_residual_delta"]))}</td><td>π 非負且總和為 1 的數值殘差差</td></tr>
    <tr><th scope="row">Peak per unit</th><td>{_h(_n(profile["peak_pattern_rows_per_unit"]))} patterns / {_h(_n(profile["peak_candidates_per_unit"]))} candidates</td><td>只保留一個 primary unit；非全表 materialization</td></tr>
    <tr><th scope="row">Peak per candidate</th><td>{_h(_n(profile["peak_states_per_candidate"]))} states / {_h(_n(profile["peak_emission_cells_per_candidate"]))} emission cells</td><td>streaming memory contract 的最大局部工作集</td></tr>
  </tbody></table></div>
  <p class="visual-note">此驗證證明在保存的 likelihood contract 下，LL、KKT certificate、relative weights 與唯一／並列分類可被另一份實作重算；它仍不把 mutation-state snapshot 升格為唯一 parent edge、真實 clone 數或 cell-level HP1/HP2 配對。</p>
"""
    terminal_dataset_detail = "" if numeric_source is not None else f'<details open><summary>各 dataset × PS_HP1/PS_HP2 的 solver、likelihood、T/V 與 topology</summary><div class="table-wrap"><table><thead><tr><th>Dataset / biological ID</th><th>HP stratum</th><th>Primary units</th><th>Solver complete / incomplete</th><th>Likelihood selection status</th><th>T / V（絕對數）</th><th>Topology / parent choice</th></tr></thead><tbody>{"".join(dataset_result_rows)}</tbody></table></div></details>'
    return f"""
<section id="m2-ranking" class="report-section">
  <p class="section-index">04 / M2 SYMBOLIC + LIKELIHOOD</p>
  <h2>先合併相同 vertex set N，再判斷 likelihood 唯一、並列或 abstain</h2>
  <p class="takeaway"><strong>教授版 display convention：bridge ≥{_h(threshold)} unique molecules，並且 PS_HP1 與 PS_HP2 分開顯示。</strong> 這個 ≥3 不是 frozen primary bridge threshold；formal release 保留 1/2/3/5 全部 sensitivity strata。它也不同於 frozen <code>primary_structural_exact_pattern_minread=3</code>。Primary unit 是 <code>HP_family×known_PS×read_linked_component×threshold</code>；missing PS 分離為 diagnostic。T 是 parent-edge tree 數；N 是一個 mutation-state vertex set；V 是不同 optimal N 的數量。Likelihood winner 只在已證明同一 minimum-extra-state h* 的候選 N 中比較，沒有與 h*+1 等更複雜的所有可行結構比較。相同 N 的不同 edges 對 snapshot reads 必須同分；canonical candidate table 只保留 parent-choice 計數，並未展開所有 edges，因此不能用同 reads 算出的 VAF 再加權一次或硬選唯一邊。{_source_ref("S7")}</p>
  <div class="metric-ledger" aria-label="Combined PS_HP1 and PS_HP2 primary ranking counts">
    <div>{_count_cell(cell["n_component_hp_units"], "HP/PS component unit", "PS_HP1+PS_HP2；bridge threshold primary")}</div>
    <div>{_count_cell(cell["raw_tree_candidates_T_complete_units"], "tree candidate T", "complete units 的 parent-edge candidate 絕對總數；不是比例")}</div>
    <div>{_count_cell(cell["distinct_vertex_sets_V_complete_units"], "vertex set V", "complete units 的 distinct optimal state-set 絕對總數；不是 T 的比例")}</div>
    <div>{_count_cell(cell["topology_evaluated_units"], "topology-evaluated unit", "只有可完整評估 topology 的 primary units")}</div>
  </div>
  <p class="visual-note">BQ emission 是 <strong>symmetric-substitution conditional R/A model</strong>：只在已固定的 R/A calls 上條件化計分；它不是包含 indel、other allele、mapping 與 missingness mechanism 的完整 Phred generative model。</p>
  {_ranking_hp_summary_html(hp_cells, cell, threshold)}
  <p class="visual-note">下圖只是 HP1+HP2 display-sum 的絕對數概覽；正式分母要看上表各 HP stratum。optimizer、k&gt;12 或候選未完整時必須留在 abstain 類別。</p>
  {chart}
  {terminal_dataset_detail}
  <details open><summary>Primary likelihood status 精確數值</summary><div class="table-wrap"><table><thead><tr><th>Status</th><th>數量、總分母、相對分母</th></tr></thead><tbody>{status_rows}</tbody></table></div></details>
  {_structural_scoring_funnel_html(hp_cells, cell, threshold)}
  {_solver_effective_k_html(hp_cells, cell, threshold)}
  {_sample_funnel_html(payload)}
  <details><summary>全部 linkage basis × bridge threshold 的 T / V 稽核表</summary><div class="table-wrap"><table><thead><tr><th>Stratum</th><th>Components</th><th>HP units</th><th>Raw T</th><th>Distinct V</th></tr></thead><tbody>{''.join(strata_rows)}</tbody></table></div></details>
  <h3>Coarse topology composition 只使用互斥 unique classes</h3>
  <p class="visual-note"><code>topology_class_inclusion_counts</code> 可重疊，禁止拿來畫比例。下圖只用 <code>coarse_topology_unique_class_counts</code>，圖內分母是 {_h(_n(cell["coarse_topology_unique_units"]))} coarse-unique units；相對全部 {_h(_n(cell["topology_evaluated_units"]))} topology-evaluated units 的比例另見表格。</p>
  {topology_unique_chart}
  <div class="claim-grid"><div><h3>互斥 unique class</h3><div class="table-wrap"><table><thead><tr><th>Coarse class</th><th>Count / denominators</th></tr></thead><tbody>{topology_unique_rows}</tbody></table></div></div><div><h3>Ambiguous class-set（另列）</h3><div class="table-wrap"><table><thead><tr><th>Possible class set</th><th>Count / denominators</th></tr></thead><tbody>{topology_ambiguous_rows}</tbody></table></div></div></div>
  {_topology_identifiability_html(hp_cells, cell)}
  <h3>Partial-read u / 2ᵘ 分布使用四個不同分母</h3>
  <p class="visual-note">熱圖每一列獨立正規化：unique patterns 回答 pattern 多樣性、quality groups 保留同 pattern 的 BQ emission group、molecule projections 回答實際 observation mass、structural patterns 只保留達 minread 者。四者不可互換，色深與百分比只在同一列比較。</p>
  {partial_chart}
  <details><summary>Partial-read 三種 legacy alias denominator（與四種完整 grain 對照）</summary><div class="table-wrap"><table><thead><tr><th>Observation grain</th><th>Count / denominator</th></tr></thead><tbody>{partial_denominator_rows}</tbody></table></div></details>
  {_partial_completion_funnel_html(hp_cells, cell)}
  {_numeric_summary_per_dataset_html(numeric_source)}
  {_candidate_examples_html(candidate_source)}
  {verification_html}
  {_runtime_and_resource_html(source, verification_source, resource_source)}
</section>
"""


def _pilot_and_hp_section(pilot: Mapping[str, Any], m2_pilot: Source | None) -> str:
    symbolic = pilot["symbolic_exhaustive"]
    cross = pilot["legacy_milp_crosscheck"]
    k9 = pilot["k9_k12"]
    raw_hp_table = ""
    if m2_pilot is not None:
        counts = m2_pilot.payload["counts"]
        raw_hp = counts["raw_HP_counts"]
        eligible_by_raw_hp = counts["canonical_eligible_alignments"]
        rows = "".join(
            f'<tr><th scope="row">HP={_h(tag)}</th><td>{_metric_cell(value, "eligible alignment", eligible_by_raw_hp, "HCC1954 chr22 eligible alignments by raw HP value（含 HP=. missing）", eligible_by_raw_hp, "本 pilot 全部 eligible alignments by raw HP value")}</td></tr>'
            for tag, value in sorted(raw_hp.items())
        )
        raw_hp_table = f"""
<details open><summary>LEGACY schema 1.0：HCC1954 chr22 raw HP tags（單染色體 context example）</summary>
  <p><strong>此 S8 只能說明舊 pilot 的 raw HP value 組成，不是目前 PS-aware schema 1.2/2.0 的資源或結果驗證。</strong> 分母是 {_h(_n(eligible_by_raw_hp))} 個 canonical eligible alignments，等於 <code>raw_HP_counts</code> 所有 raw values 的總和，<strong>明確包含 missing <code>HP=.</code></strong>。Primary family 1 合併 <code>{{1, 1-1, 1-2}}</code>；family 2 合併 <code>{{2, 2-1, 2-2}}</code>。<code>1-1/2-1</code> 是 LongPhase-S 由 somatic 資訊產生的 <strong>somatic sub-haplotype tag</strong>，不是第三條染色體、也不是一個獨立 clone；因此不可再把它當作正交驗證。{_source_ref("S8")}</p>
  <div class="table-wrap"><table><thead><tr><th>Raw tag</th><th>數量、總分母、相對分母</th></tr></thead><tbody>{rows}</tbody></table></div>
</details>"""
    return f"""
<section id="method" class="report-section">
  <p class="section-index">05 / METHOD CONTRACT</p>
  <h2>Partial read 有 2ᵘ 個概念 completions；實作使用一條 joint group constraint</h2>
  <p class="takeaway"><strong>使用者提出的「每種 completion 各跑一次，任一成功就通過」不等於全域最小解。</strong> 只有 count 達 <code>structural_exact_pattern_minread</code> 的 distinct exact R/A/X partial pattern 才形成 structural group constraint；所有 groups 一次聯合求最少 hidden。所有 informative molecules（包含低於結構門檻者）仍各自一次進 likelihood。只有 solver 證明列舉完整時，才保留全部同分 global-optimal vertex sets N，否則整個 unit abstain。likelihood 對 X 做 marginalization。{_source_ref("S4")}</p>
  {_partial_svg()}
  <div class="equation-block"><p><strong>Structural constraint</strong></p><p><code>G(p) = {{v : ((v XOR alt_mask) AND fixed_mask) = 0}}</code></p><p><code>N ∩ G(p) ≠ ∅</code></p><p><strong>Likelihood</strong></p><p><code>P(p | N, π) = Σᵥ∈N πᵥ P(p | v, error)</code></p></div>
  <p class="visual-note">正式 M2 v2 的 pattern pool 以 <code>HP_family×known_PS×read_linked_component×threshold</code> 為 primary unit；missing PS 不與 known PS 混合。BQ 只參數化 symmetric-substitution conditional R/A emission，並不代表完整 ONT/Phred generative error model。</p>
  <details open><summary>實際 partial-read 處理與較高效做法</summary>
    <ul><li>單一 pattern 有 u 個 X 時，full-cube 概念 states 是 2ᵘ；current solver 的有效數量是 <code>2^popcount(free_mask ∩ structural_ALT_universe)</code>。每次 MILP construction 會為每個 reduced group 列出 compatible active-vertex indices 作為一條 sparse hit row，但<strong>不會把每個 completion 當成一個獨立 tree world 各跑一次</strong>。</li><li>多條 reads 的 naive joint completion 是 <code>2^(Σuᵢ)</code>；正確模型是每個達結構 minread 的 distinct exact R/A/X pattern 一個 <code>Σ compatible xᵥ ≥ 1</code>，所有 constraints 聯合求解，避免跨 reads 笛卡兒積；低於結構門檻的 informative molecules 仍保留於 scoring。</li><li>目前已做 exact-preserving duplicate/dominance、required-hit/singleton forcing、active-bit predecessor、fixed-hidden sparse no-good 與 h=0 early-complete；這些只刪除等價或被支配約束，不改 global optimum。</li><li>固定最小 hidden objective 後，以 no-good cuts 列出 distinct optimal vertex sets N；若遇 candidate cap、time limit、effective k&gt;12 或未完成列舉，<code>complete=false</code>，不得把已找到的部分集合稱為「全部合理結果」。</li><li>相同 N、不同 parent edges 對 snapshot read pattern 不可辨，必須標 <code>EDGE_NONIDENTIFIABLE</code>，不能硬選。k 更大時可考慮 BDD/ZDD、separator DP 或 branch-and-price，但必須保留 completeness/lower-bound certificate。</li></ul>
  </details>
  <h3>已通過的演算法控制不是全量生物驗證</h3>
  <div class="table-wrap"><table><thead><tr><th>Control</th><th>結果與分母</th></tr></thead><tbody>
    <tr><th scope="row">Symbolic vs explicit membership</th><td>{_metric_cell(symbolic["n_patterns"] - symbolic["mismatches"], "matching pattern", symbolic["n_patterns"], "exhaustive patterns", symbolic["n_patterns"], "tested patterns")}</td></tr>
    <tr><th scope="row">State checks</th><td>{_metric_cell(symbolic["n_state_checks"], "state check", symbolic["n_state_checks"], "executed state checks", symbolic["n_state_checks"], "executed state checks")}</td></tr>
    <tr><th scope="row">Frozen solver vs MILP</th><td>{_metric_cell(cross["n_vertex_set_checks"] - cross["n_vertex_set_mismatch"], "matching vertex-set check", cross["n_vertex_set_checks"], "comparable vertex-set checks", cross["n_vertex_set_checks"], "tested checks")}</td></tr>
    <tr><th scope="row">k=9–12 exact pilot</th><td>{_metric_cell(k9["n_pass"], "PASS case", k9["n_cases"], "k9–12 cases", k9["n_cases"], "tested cases")}</td></tr>
  </tbody></table></div>
  {raw_hp_table}
  <details><summary>HP family、PS 與 circularity 定義</summary><ul>
    <li><strong>Primary family 1</strong>：合併 raw tags <code>1, 1-1, 1-2</code>；<strong>family 2</strong>：合併 <code>2, 2-1, 2-2</code>。</li>
    <li><code>1-1/2-1</code> 等 sub-haplotype tag 是 somatic-aware phasing 的產物；它可用來分層 read evidence，但不能再當成獨立 somatic truth。</li>
    <li>HP3 是輔助 tag family，不是第三條 homolog；PS 是 phase-block QC context，不是 topology edge 或 lineage label。</li>
    <li>HP1/HP2 在區域內的 state 組合仍不能自動完成 cell-level 配對；沒有單細胞或其他正交證據時，不可宣稱完整 clone composition。</li>
  </ul></details>
</section>
"""


def _source_ledger(sources: Mapping[str, Source], output_path: Path) -> str:
    rows = []
    for source_id in sorted(sources, key=lambda value: int(value[1:])):
        source = sources[source_id]
        href = _relative_href(source.path, output_path)
        rows.append(
            f'<tr id="source-{_h(source_id)}"><th scope="row">{_h(source_id)}</th><td><a href="{_h(href)}">{_h(source.label)}</a><span class="denom">{_h(source.scope)}</span></td><td><code>{_h(source.sha256)}</code></td></tr>'
        )
    return f"""
<section id="sources" class="report-section source-ledger">
  <p class="section-index">08 / EVIDENCE LEDGER</p>
  <h2>每個 quantitative claim 只回指本地 JSON / receipt</h2>
  <p>方法文獻只支持模型與 claim boundary；不替代本地數據。SHA-256 是本次 HTML 生成時的內容身分。</p>
  <div class="table-wrap"><table><thead><tr><th>ID</th><th>來源與 scope</th><th>SHA-256</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
</section>
"""


def _definitions_section() -> str:
    rows = (
        ("Dataset / biological sample", "Dataset 是一套 pipeline 輸入與輸出；biological sample 是實際生物材料。HCC1395 與 HCC1395_DORADO 是 2 datasets，但只有 1 個 biological sample。"),
        ("S / sSNV", "Somatic single-nucleotide variant 位點數；單位是 site，不是 region、read 或 clone。"),
        ("chr1–22 / chrX scope", "本輪 primary 固定 chr1–22，讓 HP1/HP2 與兩套同源染色體的工程口徑一致。chrX 受樣本性別、非PAR hemizygosity、PAR、X-inactivation與癌症CN/LOH影響，需獨立定義copy-state/HP預期後再分析；分開處理不是因為ONT不能讀chrX。癌症autosome也可能有CN/LOH，因此chr1–22仍需CN sensitivity，不能直接假設處處二倍體。"),
        ("W", "Current canonical 以相鄰位置 gap≤50 kb 形成的 legacy positional region。M2 不把 W 當 read-defined region。"),
        ("Read-linked component", "在同一 dataset×chromosome×HP family×known PS 內，由 unique molecule 跨 cut 支持連成的分析區域；threshold 1/2/3/5 分別保存。"),
        ("k / effective k", "component-site k 是該 component 內 site 數；effective k 是達 structural minread 且曾觀測 ALT 的 active 維度。Exact solver route 依 effective k，不依 raw component-site k。"),
        ("HP family / PS", "HP1/HP2 是 haplotype family；PS 是 phase-set block。相同 HP 數字跨不同 PS 不保證方向一致，因此 primary tree 不跨 known PS 混合。"),
        ("R / A / O / D / S / L / X", "REF、ALT、other base、deletion、reference skip、low BQ、uncovered/unknown。Structural與likelihood固定資訊只用R/A；其他類別保留漏斗並以missing處理。"),
        ("C / structural pattern groups", "若沿用 C，僅指達門檻的不同 read R/A/X pattern／group terminals 數；partial groups 可重疊，所以 C 不是 clone 數。正式輸出改用 n_structural_pattern_groups 避免歧義。"),
        ("Legacy C_region（禁止與新 C 混用）", "舊 layered reconstruction 曾以 region 內各 HP tree-candidate 數的 joint product 表示候選組合量；它不是 read pattern group 數。教授版與論文應寫 legacy_joint_tree_combinations，避免再用 C_region。"),
        ("G(p) / u / 2ᵘ", "p 是一條 partial pattern；G(p) 是與它相容的 state group。u 是 X 數；full cube 概念 completion 為2ᵘ，observed-ALT structural universe 的 effective completion 可更少。"),
        ("N / hidden-extra vertex", "N 是一個候選 mutation-state vertex set，必含 root與達 structural minread 的完整 observations，且每個保留的partial structural group至少命中一個state。非root、非完整observed的額外state統稱hidden/extra；其中也可能是partial-supported，不應一律解讀為未觀察clone。"),
        ("h* / minimum-extra-state count", "h* 是全域最小 objective 下，N 中除 root 與 full-observed mandatory states 外的 state 數；包含 partial-compatible representative 與 connector。它不是 hidden clone 數。Likelihood winner 只在同一 h* 候選集內排名，不代表跨 h*+1 等所有複雜度的全域最可能真實拓撲。"),
        ("V", "全域最少 hidden/extra objective 下，不同 optimal mutation-state vertex sets N 的數量。V>1 才需要問 read likelihood 能否排序不同 state sets。"),
        ("T", "在所有 optimal vertex sets 上展開每個節點 parent choice 後的原始 parent-edge tree 候選總數。T>1、V=1 表示只有 edge/順序仍不唯一。"),
        ("Topo", "去除部分標籤後的樹形類別。Coarse class、exact labeled topology與shape-isomorphism是不同層級，報告不可互換。"),
        ("π / likelihood", "π 是固定候選 vertex set 下擬合的 latent expected mixture proportion；不是每條read的硬指派、cell fraction或已確認clone proportion。Likelihood唯一只代表在已列舉candidate set與宣告error model下有唯一第一。"),
        ("Unique / tied / abstain", "Unique：唯一第一vertex set；tied：多個vertex sets在容許差內並列；abstain：solver enumeration、effective-k或optimizer證據不足，保留未判定而不硬選。"),
    )
    body = "".join(f'<tr><th scope="row">{_h(term)}</th><td>{_h(definition)}</td></tr>' for term, definition in rows)
    return f"""
<section id="definitions" class="report-section">
  <p class="section-index">00 / DEFINITIONS</p>
  <h2>先固定單位與名詞，避免把 site、region、HP unit、state 與 clone 混為一談</h2>
  <p class="takeaway"><strong>本報告的核心單位會隨層級改變。</strong> 每個數字旁都必須同時讀取 unit 與 denominator；名稱相似不代表可以直接相除。</p>
  <details open><summary>名詞與符號完整定義</summary><div class="table-wrap"><table><thead><tr><th>名詞 / 符號</th><th>本報告中的精確意思</th></tr></thead><tbody>{body}</tbody></table></div></details>
</section>
"""


def _css() -> str:
    return r"""
:root{--paper:#f4efe4;--paper-2:#ebe3d3;--ink:#171511;--muted:#625e56;--rule:#9a9182;--accent:#b33424;--accent-soft:#d98b78;--white:#fffdf7;--focus:#005fcc}
*{box-sizing:border-box}html{scroll-behavior:smooth;background:var(--paper);color:var(--ink)}body{margin:0;font-family:"Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif;line-height:1.72;background:linear-gradient(90deg,var(--paper-2) 0 18px,var(--paper) 18px 100%)}a{color:var(--accent);text-underline-offset:.18em}a:focus-visible,summary:focus-visible{outline:3px solid var(--focus);outline-offset:3px}.shell{max-width:1180px;margin:0 auto;padding:44px 54px 90px}.masthead{border-top:10px solid var(--ink);border-bottom:1px solid var(--ink);padding:25px 0 30px;position:relative}.masthead:after{content:"";position:absolute;right:0;top:25px;width:96px;height:96px;border:1px solid var(--accent);border-radius:50%;background:radial-gradient(circle at center,var(--accent) 0 3px,transparent 4px)}.eyebrow,.section-index{font-family:ui-monospace,"SFMono-Regular",Consolas,monospace;letter-spacing:.14em;text-transform:uppercase;font-size:.76rem;color:var(--accent);font-weight:700}.masthead h1{font-family:"Noto Serif TC","Songti TC","PMingLiU",serif;font-size:clamp(2.25rem,5vw,4.8rem);line-height:1.02;max-width:920px;margin:.28em 0 .24em;letter-spacing:-.04em}.dek{font-family:"Noto Serif TC","Songti TC","PMingLiU",serif;font-size:1.18rem;max-width:820px;margin:0}.meta-line{display:flex;gap:22px;flex-wrap:wrap;margin-top:22px;color:var(--muted);font-size:.9rem}.status-ribbon{position:sticky;z-index:20;top:0;background:var(--accent);color:#fff;padding:11px 18px;font-weight:800;letter-spacing:.08em;text-align:center;box-shadow:0 4px 0 rgba(23,21,17,.18)}.status-ribbon.final{background:var(--ink)}.toc{display:flex;gap:8px 20px;flex-wrap:wrap;padding:18px 0;border-bottom:1px solid var(--rule);font-size:.88rem}.toc a{color:var(--ink);text-decoration:none;border-bottom:2px solid transparent}.toc a:hover{border-color:var(--accent)}.technical-summary{padding:34px 0 28px;border-bottom:4px double var(--ink)}.technical-summary h2,.report-section h2{font-family:"Noto Serif TC","Songti TC","PMingLiU",serif;line-height:1.16;letter-spacing:-.025em}.technical-summary h2{font-size:2rem;margin:.2em 0}.summary-grid{display:grid;grid-template-columns:1.25fr .75fr;gap:36px}.verdict{font-size:1.12rem}.gate-list{margin:0;padding:18px 22px;background:var(--ink);color:var(--paper)}.gate-list li+li{margin-top:8px}.report-section{padding:58px 0 48px;border-bottom:1px solid var(--ink)}.report-section h2{font-size:clamp(1.8rem,3vw,2.85rem);max-width:900px;margin:.16em 0 .48em}.report-section h3{font-family:"Noto Serif TC","Songti TC","PMingLiU",serif;font-size:1.45rem;margin-top:38px}.takeaway{font-size:1.07rem;max-width:940px;border-left:5px solid var(--accent);padding-left:20px}.source-ref a{display:inline-block;font:700 .68rem ui-monospace,monospace;text-decoration:none;border:1px solid currentColor;padding:0 .26em;margin-left:.3em;vertical-align:super}.metric-ledger{display:grid;grid-template-columns:repeat(4,1fr);border-top:1px solid var(--ink);border-bottom:1px solid var(--ink);margin:30px 0}.metric-ledger>div{padding:20px 18px 20px 0}.metric-ledger>div+div{border-left:1px solid var(--rule);padding-left:18px}.metric-main{display:block;font-family:ui-monospace,"SFMono-Regular",Consolas,monospace;font-size:1.45rem;font-weight:800;line-height:1.2}.metric-main small{font:600 .7rem "Noto Sans TC",sans-serif;color:var(--muted)}.denom{display:block;font-size:.72rem;line-height:1.42;color:var(--muted);margin-top:4px}.visual-note{font-size:.92rem;color:var(--muted);max-width:880px}.chart,.process-svg{display:block;width:100%;height:auto;margin:18px 0 30px;background:var(--white);border-top:1px solid var(--ink);border-bottom:1px solid var(--ink)}.svg-label{font:600 14px "Noto Sans TC",sans-serif;fill:var(--ink)}.svg-value{font:700 13px ui-monospace,monospace;fill:var(--ink)}.svg-track{stroke:#d9d1c2;stroke-width:18}.flow-box{fill:var(--white);stroke:var(--ink);stroke-width:1.5}.flow-box.accent{fill:#f1d8cf;stroke:var(--accent)}.flow-box.dark{fill:var(--ink);stroke:var(--ink)}.flow-box.dark+text.cube-formula{fill:var(--paper)}.flow-arrow{stroke:var(--accent);stroke-width:2;fill:none}.flow-title{font:700 16px "Noto Sans TC",sans-serif;fill:var(--ink)}.flow-sub{font:12px ui-monospace,monospace;fill:var(--muted)}.svg-kicker{font:700 14px "Noto Sans TC",sans-serif;fill:var(--accent)}.cube-code{font:800 24px ui-monospace,monospace;fill:var(--ink)}.cube-formula{font:700 18px ui-monospace,monospace;fill:var(--ink)}.svg-note{font:13px "Noto Sans TC",sans-serif;fill:var(--muted)}.equation-block{display:grid;grid-template-columns:auto 1fr;gap:0 24px;border-top:1px solid var(--rule);border-bottom:1px solid var(--rule);padding:14px 0;margin:24px 0}.equation-block p{margin:4px 0}.equation-block code{font-size:1rem}.table-wrap{overflow-x:auto;margin:20px 0}table{border-collapse:collapse;width:100%;font-size:.86rem;background:rgba(255,253,247,.45)}th,td{text-align:left;vertical-align:top;padding:14px 12px;border-bottom:1px solid var(--rule)}thead th{font-family:ui-monospace,monospace;font-size:.72rem;letter-spacing:.04em;text-transform:uppercase;border-bottom:2px solid var(--ink)}tbody th{min-width:165px}.source-ledger code{font-size:.69rem;word-break:break-all}details{border-top:1px solid var(--rule);margin-top:26px}summary{cursor:pointer;padding:14px 0;font-weight:750;color:var(--ink)}details[open] summary{color:var(--accent)}.unavailable{background:repeating-linear-gradient(-45deg,transparent 0 15px,rgba(179,52,36,.045) 15px 30px);padding-left:24px;padding-right:24px}.claim-grid{display:grid;grid-template-columns:1fr 1fr;gap:30px}.claim-grid>div{border-top:4px solid var(--ink);padding-top:13px}.claim-grid .ceiling{border-color:var(--accent)}.footer{padding-top:34px;color:var(--muted);font-size:.82rem}.print-only{display:none}
@media(max-width:800px){body{background:var(--paper)}.shell{padding:25px 20px 60px}.masthead:after{display:none}.summary-grid,.claim-grid{grid-template-columns:1fr}.metric-ledger{grid-template-columns:1fr 1fr}.metric-ledger>div:nth-child(3){border-left:0}.metric-ledger>div:nth-child(n+3){border-top:1px solid var(--rule)}.chart{min-width:760px}.process-svg{min-width:840px}.report-section{overflow-x:auto}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
@media print{@page{size:A4;margin:14mm 13mm}body{background:#fff;color:#000}.status-ribbon{position:static;box-shadow:none}.shell{max-width:none;padding:0}.toc{display:none}.masthead{break-after:avoid}.report-section{break-before:page;padding-top:10mm}.report-section:first-of-type{break-before:auto}.chart,.process-svg{max-width:100%;break-inside:avoid;background:#fff}.table-wrap{overflow:visible}table{font-size:7.5pt}th,td{padding:5pt}.metric-ledger{break-inside:avoid}.source-ledger code{font-size:5.8pt}details>summary{display:none}details>*{display:block!important}.print-only{display:block}.unavailable{background:#fff;border-left:3px solid #000}}
"""


def _render_html(assessment: Assessment, output_path: Path) -> str:
    canonical = assessment.sources["S1"].payload
    pilot = assessment.sources["S3"].payload
    final_ready = assessment.final_ready
    now = dt.datetime.now(ZoneInfo("Asia/Taipei"))
    status = "FINAL · TASK-B EVIDENCE GATE PASS" if final_ready else "PARTIAL PREVIEW · NOT VALIDATION EVIDENCE"
    gate_issues = assessment.final_issues if assessment.final_issues else ["所有 final gate 已通過"]
    issue_items = "".join(f"<li>{_h(issue)}</li>" for issue in gate_issues)
    m0 = assessment.sources.get("S2")
    extraction = assessment.sources.get("S6")
    ranking = assessment.sources.get("S7")
    m2_pilot = assessment.sources.get("S8")
    scope = canonical.get("scope", "7 datasets / 6 biological samples / chr1-22")
    quantitative_source_ids = [
        source_id
        for source_id in ("S1", "S10", "S14", "S2", "S11", "S12", "S3", "S6", "S7", "S13", "S15", "S16", "S17", "S8", "S9")
        if source_id in assessment.sources
    ]
    quantitative_source_text = " / ".join(quantitative_source_ids)
    title = "Read-linked Hypercube Steiner 與 likelihood 全量驗證"
    return f"""<!doctype html>
<html lang="zh-Hant">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="light">
  <meta name="generator" content="InterSubMod build_validated_html_report.py">
  <title>{_h(title)}</title>
  <style>{_css()}</style>
</head>
<body>
<div class="status-ribbon {'final' if final_ready else ''}" role="status">{_h(status)}</div>
<main class="shell">
  <header class="masthead">
    <p class="eyebrow">InterSubMod · Professor methods brief · G3 / G4 / G5</p>
    <h1>{_h(title)}</h1>
    <p class="dek">從 current 50 kb candidate baseline，轉向以 ONT unique molecules 定義區域；在 Boolean hypercube 上求最少 hidden mutation states，再只對 PS-aware、資料可辨識的 vertex sets 做 conditional R/A likelihood 排序。</p>
    <div class="meta-line"><span>Evidence cutoff：{_h(now.strftime('%Y-%m-%d %H:%M %Z'))}</span><span>Scope：{_h(scope)}</span><span>7 datasets ≠ 7 biological samples</span></div>
  </header>
  <nav class="toc" aria-label="報告目錄"><a href="#summary">技術摘要</a><a href="#definitions">定義</a><a href="#baseline">Current baseline</a><a href="#m0">M0</a><a href="#m2-extraction">M2 extraction</a><a href="#m2-ranking">M2 ranking</a><a href="#method">Partial read</a><a href="#claims">主張邊界</a><a href="#sources">來源</a></nav>
  <section id="summary" class="technical-summary">
    <p class="section-index">TECHNICAL SUMMARY</p>
    <h2>{'全量資料鏈已閉合；結果可用於教授版方法與數量報告' if final_ready else '方法與 baseline 可解釋；full-M2 數據鏈尚未閉合'}</h2>
    <div class="summary-grid"><div class="verdict">
      <p><strong>核心判斷：</strong>partial read 不是拆成多棵獨立樹再「任一成功」；達 structural minread 的 distinct pattern 形成一個 subcube group constraint。solver 必須聯合滿足所有 groups、先最小化 hidden vertices、再列完全部同分 optimal vertex sets N；V 只是這些 N 的數量。ranking unit 必須隔離 known PS；conditional R/A likelihood 對 partial dimension 做 marginalization，不把同一 read 的 VAF 再算第二次。</p>
      <p><strong>區域判斷：</strong>canonical 50 kb W 只作歷史比較。M2 終版的區域必須由 unique molecule bridge 定義，並把 bridge threshold、k&gt;12 local/abstain 與所有漏斗分母一併報告。</p>
    </div><ul class="gate-list"><li><strong>Artifact status</strong>：{_h(status)}</li>{issue_items}</ul></div>
  </section>
  <section class="report-section" aria-label="整體流程"><p class="section-index">METHOD AT A GLANCE</p><h2>六步驟把「可連接」與「可排序」分開</h2><p class="visual-note">流程前半定義 read-linked evidence 與可行候選；後半才問候選是否可由 molecule patterns 區分。</p>{_workflow_svg()}</section>
  {_definitions_section()}
  {_canonical_sections(canonical, assessment.sources.get("S10"), assessment.sources.get("S14"))}
  {_m0_section(m0, canonical)}
  {_m2_extraction_section(extraction)}
  {_m2_ranking_section(
      ranking,
      assessment.sources.get("S9"),
      assessment.sources.get("S13"),
      assessment.sources.get("S15"),
      assessment.sources.get("S16"),
  )}
  {_pilot_and_hp_section(pilot, m2_pilot)}
  <section id="claims" class="report-section">
    <p class="section-index">06 / CLAIM CEILING</p><h2>最終可稱 solver-certified regional mutation-state candidates；不可直接稱真實 clone tree</h2>
    <div class="claim-grid"><div><h3>可以合理宣稱</h3><ul><li>在宣告的 molecule、HP、PS、error、bridge-threshold 與 solver contract 下，可行且最少 hidden/extra 的 regional mutation-state vertex sets。</li><li>哪些 vertex sets 在 read-pattern likelihood 下呈唯一第一、並列第一或 fail-closed abstain。</li><li>7 technical datasets 均以同一份可重算工程契約處理；HCC1395 與 Dorado 必須分開列示並映射到同一 biological ID。</li></ul></div><div class="ceiling"><h3>目前不可宣稱</h3><ul><li>真實 cellular clone 數或完整 HP1/HP2 cell-level 配對。</li><li>相同 vertex set 中唯一 parent edge；snapshot patterns 對此不可辨。</li><li>甲基已「確認」拓撲；除非有預先定義且獨立的正交驗證。</li><li>k&gt;12 local blocks 等同 global optimum。</li><li>只看 HCC1395 / Dorado aggregate counts 就宣稱 technical concordance；需另做 genomic matched-unit 分析。</li></ul></div></div>
    <p class="visual-note">方法與原始文獻的對齊、以及何者是工程結果／生物主張，見 {_source_ref('S5')}。本報告可見的 quantitative results 只引用 {_h(quantitative_source_text)}。</p>
  </section>
  <section class="report-section"><p class="section-index">07 / NEXT DECISION</p><h2>{'固定 PS-aware primary stratum，評估固定候選集條件下的 ranking stability' if final_ready else 'M0 legacy baseline 已完成；補齊 M2 extraction 1.2 與 ranking 2.0 全量 receipts 後，才能解除 PARTIAL'}</h2><ol><li>保持 1/2/3/5 bridge thresholds 全部輸出；預先指定 PS-aware primary stratum，其他作 sensitivity。</li><li>保存 k route、solver complete/incomplete、likelihood tie/abstain；bootstrap／BQ／error-grid 只稱 <strong>conditional-on-fixed-candidate-set ranking stability</strong>，不可稱整條 tree stability。</li><li>若需確認 parent edge 或真實 clone pairing，另設獨立 lineage／single-cell／multi-region／methylation validation；不可從同 reads 重複計分。</li></ol></section>
  {_source_ledger(assessment.sources, output_path)}
  <footer class="footer"><p>Standalone HTML · no external JavaScript / font / image dependency · A4 print stylesheet · generated from SHA-identified local evidence. {'FINAL artifact 已通過本報告宣告的 Task-B receipt gates；生物主張仍受 claim ceiling 限制。' if final_ready else 'PARTIAL preview 不可作論文、口試或 validation evidence。'}</p></footer>
</main>
</body>
</html>
"""


def _partial_path(output_path: Path) -> Path:
    name = output_path.name
    lowered = name.lower()
    if "partial" in lowered or "preview" in lowered:
        return output_path
    return output_path.with_name(f"{output_path.stem}.partial-preview{output_path.suffix or '.html'}")


def _rename_noreplace(source: Path, destination: Path) -> None:
    """Atomically publish without ever replacing an existing path."""

    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise RuntimeError("renameat2(RENAME_NOREPLACE) is unavailable")
    renameat2.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    renameat2.restype = ctypes.c_int
    result = renameat2(
        -100,
        os.fsencode(source),
        -100,
        os.fsencode(destination),
        1,
    )
    if result == 0:
        return
    observed_errno = ctypes.get_errno()
    if observed_errno in {errno.EEXIST, errno.ENOTEMPTY}:
        raise FileExistsError(f"refusing to overwrite existing report: {destination}")
    raise OSError(observed_errno, os.strerror(observed_errno), str(destination))


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _archive_existing_report(path: Path) -> Path:
    """Preserve an explicitly superseded report; never delete or truncate it."""

    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    archived = path.with_name(f"{path.name}.superseded.{stamp}.{os.getpid()}")
    _rename_noreplace(path, archived)
    return archived


def build_report(
    *,
    canonical_json: Path,
    funnel_receipt: Path | None,
    funnel_verification_receipt: Path | None,
    m0_receipt: Path | None,
    m0_verification_receipt: Path | None,
    pilot_receipt: Path,
    method_audit: Path,
    literature_audit: Path,
    m2_extraction_receipt: Path | None,
    m2_ranking_receipt: Path | None,
    m2_verification_receipt: Path | None,
    m2_numeric_summary: Path | None,
    m2_pilot_extraction_receipt: Path | None,
    output_path: Path,
    allow_partial: bool = False,
    overwrite: bool = False,
) -> BuildResult:
    assessment = assess_sources(
        canonical_json=canonical_json,
        funnel_receipt=funnel_receipt,
        funnel_verification_receipt=funnel_verification_receipt,
        m0_receipt=m0_receipt,
        m0_verification_receipt=m0_verification_receipt,
        pilot_receipt=pilot_receipt,
        method_audit=method_audit,
        literature_audit=literature_audit,
        m2_extraction_receipt=m2_extraction_receipt,
        m2_ranking_receipt=m2_ranking_receipt,
        m2_verification_receipt=m2_verification_receipt,
        m2_numeric_summary=m2_numeric_summary,
        m2_pilot_extraction_receipt=m2_pilot_extraction_receipt,
    )
    if assessment.hard_issues:
        raise ReportGateError("hard evidence validation failed:\n- " + "\n- ".join(assessment.hard_issues))
    if not assessment.final_ready and not allow_partial:
        raise ReportGateError("final report gate failed:\n- " + "\n- ".join(assessment.final_issues))
    actual_output = output_path if assessment.final_ready else _partial_path(output_path)
    actual_output.parent.mkdir(parents=True, exist_ok=True)
    parent_state = os.lstat(actual_output.parent)
    if not stat.S_ISDIR(parent_state.st_mode) or stat.S_ISLNK(parent_state.st_mode):
        raise RuntimeError(f"report output parent is not a real directory: {actual_output.parent}")
    if os.path.lexists(actual_output):
        if not overwrite:
            raise FileExistsError(f"refusing to overwrite existing report: {actual_output}")
        _archive_existing_report(actual_output)
    document = _render_html(assessment, actual_output)
    if "None" in document or re.search(r"(?<![A-Za-z])nan(?![A-Za-z])", document, flags=re.IGNORECASE):
        raise RuntimeError("rendered report contains an unhandled missing numeric value")
    if "http://" in document or "https://" in document or "<script" in document.lower():
        raise RuntimeError("standalone report unexpectedly contains remote URL or script")
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=actual_output.parent, prefix=f".{actual_output.name}.", suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(document)
        handle.flush()
        os.fsync(handle.fileno())
        os.fchmod(handle.fileno(), 0o444)
    _rename_noreplace(temporary, actual_output)
    _fsync_directory(actual_output.parent)
    return BuildResult(
        output_path=actual_output,
        final_ready=assessment.final_ready,
        hard_issues=tuple(assessment.hard_issues),
        final_issues=tuple(assessment.final_issues),
    )


def _args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-json", required=True, type=Path)
    parser.add_argument("--funnel-receipt", type=Path)
    parser.add_argument("--funnel-verification-receipt", type=Path)
    parser.add_argument("--m0-receipt", type=Path)
    parser.add_argument("--m0-verification-receipt", type=Path)
    parser.add_argument("--pilot-receipt", required=True, type=Path)
    parser.add_argument("--method-audit", required=True, type=Path)
    parser.add_argument("--literature-audit", required=True, type=Path)
    parser.add_argument("--m2-extraction-receipt", type=Path)
    parser.add_argument("--m2-ranking-receipt", type=Path)
    parser.add_argument("--m2-verification-receipt", type=Path)
    parser.add_argument("--m2-numeric-summary", required=True, type=Path)
    parser.add_argument("--m2-pilot-extraction-receipt", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = _args()
    try:
        result = build_report(
            canonical_json=args.canonical_json,
            funnel_receipt=args.funnel_receipt,
            funnel_verification_receipt=args.funnel_verification_receipt,
            m0_receipt=args.m0_receipt,
            m0_verification_receipt=args.m0_verification_receipt,
            pilot_receipt=args.pilot_receipt,
            method_audit=args.method_audit,
            literature_audit=args.literature_audit,
            m2_extraction_receipt=args.m2_extraction_receipt,
            m2_ranking_receipt=args.m2_ranking_receipt,
            m2_verification_receipt=args.m2_verification_receipt,
            m2_numeric_summary=args.m2_numeric_summary,
            m2_pilot_extraction_receipt=args.m2_pilot_extraction_receipt,
            output_path=args.output,
            allow_partial=args.allow_partial,
            overwrite=args.overwrite,
        )
    except (ReportGateError, FileExistsError, FileNotFoundError, ValueError) as exc:
        print(json.dumps({"all_pass": False, "error": str(exc)}, ensure_ascii=False, indent=2))
        raise SystemExit(2)
    print(
        json.dumps(
            {
                # A successfully rendered preview is not an evidence-gate pass.
                # Keep these fields separate so downstream automation cannot
                # accidentally promote a partial artifact to a final result.
                "generation_pass": True,
                "all_pass": result.final_ready,
                "final_ready": result.final_ready,
                "artifact_status": "FINAL_VALIDATION_EVIDENCE" if result.final_ready else "PARTIAL_PREVIEW_NOT_VALIDATION_EVIDENCE",
                "output": str(result.output_path),
                "output_size_bytes": result.output_path.stat().st_size,
                "output_sha256": _sha256(result.output_path),
                "final_issues": list(result.final_issues),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
