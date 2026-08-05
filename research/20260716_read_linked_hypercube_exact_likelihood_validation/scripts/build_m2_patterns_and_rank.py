#!/usr/bin/env python3
"""Build M2 molecule patterns and rank exact Hypercube vertex-set candidates.

This consumer starts from one chromosome-level ``extract_lossless_read_linkage``
output directory.  It deliberately keeps candidate generation and data scoring as
two different contracts:

* structural patterns are informative R/A/X molecule patterns with count >=
  ``minread``;
* scoring patterns contain every informative molecule before that threshold;
* an all-X projection is conserved in the funnel but omitted from likelihood,
  because it has probability one for every candidate state set;
* partial patterns remain one symbolic group each.  Completion-wise tree worlds
  and cross-read Cartesian products are not materialized; the reduced group's
  compatible active vertex indices are materialized per MILP construction as
  one sparse group-hit row;
* candidates with observed-ALT active k<=the declared exact limit are collapsed
  by vertex set before a mutation-state mixture likelihood is fitted.  Snapshot
  reads never score parent edges;
* observed-ALT active k above the exact limit is explicitly local-only/abstain
  and is never presented as a global optimum.

The primary likelihood uses each fixed R/A call's Phred BQ under the declared
conditional biallelic symmetric-substitution model; a fixed-error grid is only a
sensitivity analysis.  VAF derived from these same molecules is not added as a
second score, which prevents double counting.
"""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import sys
import tempfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hypercube_exact import (  # noqa: E402
    SymbolicPattern,
    MixtureFit,
    enumerate_optimal_vertex_sets,
    fit_emission_mixture_certified,
    fit_vertex_mixture_slsqp,
    parent_choice_count,
    parse_full,
    vertex_set_digest,
)


FIXED_CODES = frozenset({"R", "A"})
RAW_CODES = frozenset({"R", "A", "O", "D", "S", "L", "X"})
PRIMARY_HP_FAMILIES = ("1", "2")
PRIMARY_COMPONENT_BASES = ("PS_HP1", "PS_HP2")
EXTRACTOR_SCHEMA_VERSION = "1.2.0"
RANKING_SCHEMA_VERSION = "2.0.0"
STRUCTURAL_CACHE_SCHEMA_VERSION = "1.0.0"
OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE = 1e-8
_monotonic_ns = time.monotonic_ns

# This declaration describes implementation semantics, not an empirical check.
# Keep it under ``parameters.method_contract`` so downstream aggregators and an
# independent verifier can exact-compare the contract and bind it to this
# source file's physical SHA-256.
METHOD_CONTRACT = {
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

PATTERN_FIELDS = (
    "dataset",
    "chrom",
    "threshold",
    "component_basis",
    "phase_set",
    "phase_set_status",
    "inference_role",
    "component_id",
    "family",
    "structural_exact_pattern_minread",
    "structural_minread_role",
    "k",
    "pattern_rax",
    "fixed_base_qualities",
    "n_molecules",
    "rax_pattern_total_molecules",
    "scoring_eligible",
    "structural_retained",
    "exclusion_reason",
    "n_free",
    "n_conceptual_completions",
    "n_effective_free_in_structural_alt_universe",
    "n_effective_completions",
    "effective_completion_status",
    "pattern_semantic_sha256",
)

UNIT_FIELDS = (
    "dataset",
    "chrom",
    "threshold",
    "component_basis",
    "phase_set",
    "phase_set_status",
    "inference_role",
    "component_id",
    "family",
    "structural_exact_pattern_minread",
    "structural_minread_role",
    "start1",
    "end1",
    "k",
    "k_component_sites",
    "k_observed_alt_active",
    "k_scoring_alt_observed",
    "n_not_structural_alt_active_sites",
    "molecule_component_projections",
    "informative_scoring_molecules",
    "all_x_excluded_molecules",
    "n_scoring_pattern_groups",
    "n_scoring_quality_groups",
    "structural_retained_molecules",
    "below_minread_scoring_molecules",
    "n_structural_pattern_groups",
    "projected_cells",
    "fixed_ra_cells",
    "alt_cells",
    "other_cells",
    "deletion_cells",
    "refskip_cells",
    "low_baseq_cells",
    "explicit_x_cells",
    "implicit_uncovered_x_cells",
    "non_ra_cells_marginalized",
    "bq_scoring_molecules",
    "minimum_hidden_nodes",
    "raw_tree_candidates_T",
    "distinct_vertex_sets_V",
    "candidate_vertex_sets_enumerated",
    "candidate_generation_complete",
    "candidate_generation_status",
    "best_vertex_sets",
    "top_edge_variants",
    "best_log_likelihood",
    "second_log_likelihood",
    "delta_best_second",
    "top_relative_likelihood_weight",
    "selection_status",
    "all_fits_converged",
    "all_fits_monotone",
    "max_emission_rank",
    "vertex_set_ids",
    "top_vertex_set_ids",
    "topology_classes",
    "n_topology_classes",
    "coarse_topology_class_unique",
    "parent_edge_assignment_unique",
    "exact_topology_unique",
    "exact_topology_uniqueness_status",
    "topology_class_by_top_vertex_set",
    "topology_derivation_status",
    "fixed_error_grid_top_vertex_set_ids",
    "fixed_error_grid_stable_with_quality_primary",
    "fixed_error_grid_all_converged",
    "structural_partial_pattern_groups",
    "partial_group_coverage_denominator",
    "partial_groups_covered",
    "partial_groups_unsatisfied",
    "conditional_candidate_ranking_bootstrap_status",
    "conditional_candidate_ranking_bootstrap_replicates",
    "conditional_candidate_ranking_bootstrap_seed",
    "conditional_candidate_ranking_bootstrap_top_vertex_set_frequency",
    "conditional_candidate_ranking_bootstrap_primary_top_set_frequency",
    "conditional_candidate_ranking_bootstrap_all_converged",
    "unit_semantic_sha256",
)

RUNTIME_METRICS = (
    "candidate_generation_elapsed_seconds",
    "likelihood_fit_elapsed_seconds",
    "unit_total_elapsed_seconds",
)

RUNTIME_DIAGNOSTIC_FIELDS = (
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

CANDIDATE_FIELDS = (
    "dataset", "chrom", "component_basis", "phase_set", "threshold", "component_id", "family",
    "structural_exact_pattern_minread", "vertex_set_id", "states_json", "parent_choice_count",
    "unique_parent_edges_json", "primary_log_likelihood", "relative_likelihood_weight",
    "mixture_weights_json", "fit_converged", "fit_monotone", "optimizer_status",
    "slsqp_success", "slsqp_status", "slsqp_message", "warm_start_log_likelihood",
    "warm_start_global_log_likelihood_gap_bound", "refinement_iterations",
    "global_log_likelihood_gap_bound", "simplex_residual", "augmented_emission_rank",
    "mixture_weights_identifiable", "coarse_topology_classes_json",
    "is_winner", "is_tied_winner", "candidate_semantic_sha256",
)

RESPONSIBILITY_FIELDS = (
    "dataset", "chrom", "component_basis", "phase_set", "threshold", "component_id", "family",
    "structural_exact_pattern_minread", "vertex_set_id", "pattern_rax", "fixed_base_qualities",
    "n_molecules", "state_bitmask", "state_rax", "posterior_responsibility",
    "interpretation", "responsibility_semantic_sha256",
)


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_sha256_sidecar(path: Path) -> Path:
    sidecar = path.with_name(f"{path.name}.sha256")
    sidecar.write_text(f"{sha256_path(path)}  {path.name}\n", encoding="ascii")
    return sidecar


def semantic_digest(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    ).hexdigest()


class StructuralEnumerationCache:
    """Process-local cache for fully proven structural enumerations only."""

    def __init__(
        self,
        *,
        solver_source_sha256: str | None = None,
        diagnostic_digest: Callable[[Any], str] = semantic_digest,
        enumerate_function: Callable[..., dict[str, Any]] | None = None,
    ) -> None:
        self.solver_source_sha256 = (
            solver_source_sha256
            if solver_source_sha256 is not None
            else sha256_path(SCRIPT_DIR / "hypercube_exact.py")
        )
        self._diagnostic_digest = diagnostic_digest
        self._enumerate = (
            enumerate_optimal_vertex_sets
            if enumerate_function is None
            else enumerate_function
        )
        self._entries: dict[tuple[Any, ...], dict[str, Any]] = {}
        self._contexts: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
        self._digest_keys: dict[str, set[tuple[Any, ...]]] = defaultdict(set)
        self.lookups = 0
        self.hits = 0
        self.misses = 0
        self.stores_complete = 0
        self.rejected_incomplete = 0
        self.solver_invocations = 0
        self.evictions = 0
        self.cross_minread_hits = 0
        self.cross_threshold_hits = 0
        self.same_unit_cross_minread_hits = 0
        self.digest_collision_guard_events = 0

    def _key(
        self,
        *,
        full_patterns: Sequence[str],
        partial_patterns: Sequence[str],
        k: int,
        structural_alt_universe_mask: int,
        exact_k_max: int,
        max_vertex_sets: int,
        solver_time_limit_seconds: float,
        universe_mode: str,
    ) -> tuple[Any, ...]:
        return (
            STRUCTURAL_CACHE_SCHEMA_VERSION,
            self.solver_source_sha256,
            universe_mode,
            int(k),
            int(structural_alt_universe_mask),
            tuple(sorted(str(pattern) for pattern in full_patterns)),
            tuple(sorted(str(pattern) for pattern in partial_patterns)),
            int(exact_k_max),
            int(max_vertex_sets),
            float(solver_time_limit_seconds).hex(),
        )

    @staticmethod
    def _complete_result_is_valid(
        result: Mapping[str, Any],
        *,
        full_patterns: Sequence[str],
        partial_patterns: Sequence[str],
        k: int,
    ) -> bool:
        if result.get("complete") is not True:
            return False
        first = result.get("first") or {}
        if first.get("status") != "OPTIMAL":
            return False
        objective = result.get("objective")
        if (
            isinstance(objective, bool)
            or not isinstance(objective, (int, float))
            or not math.isfinite(float(objective))
            or float(objective) < 0
            or not float(objective).is_integer()
        ):
            return False
        objective_int = int(objective)
        vertex_sets = result.get("vertex_sets")
        vertex_set_ids = result.get("vertex_set_ids")
        if not isinstance(vertex_sets, list) or not vertex_sets:
            return False
        if not isinstance(vertex_set_ids, list) or len(vertex_set_ids) != len(vertex_sets):
            return False
        mandatory = {0, *(parse_full(pattern) for pattern in full_patterns)}
        partial_groups = tuple(
            SymbolicPattern.from_string(pattern) for pattern in partial_patterns
        )
        expected_ids: list[str] = []
        seen_sets: set[tuple[int, ...]] = set()
        for values in vertex_sets:
            normalized = tuple(int(value) for value in values)
            selected = set(normalized)
            if (
                len(selected) != len(normalized)
                or any(value < 0 or value >= (1 << k) for value in selected)
                or not mandatory <= selected
                or len(selected - mandatory) != objective_int
                or parent_choice_count(normalized) < 1
                or any(
                    not any(group.compatible(vertex) for vertex in selected)
                    for group in partial_groups
                )
            ):
                return False
            canonical = tuple(sorted(selected))
            if canonical in seen_sets:
                return False
            seen_sets.add(canonical)
            expected_ids.append(vertex_set_digest(k, canonical))
        return expected_ids == [str(value) for value in vertex_set_ids]

    @staticmethod
    def _reuse_relation(
        previous_contexts: Sequence[Mapping[str, Any]],
        context: Mapping[str, Any],
    ) -> tuple[str, bool, bool, bool]:
        cross_minread = any(
            previous.get("minread") != context.get("minread")
            for previous in previous_contexts
        )
        cross_threshold = any(
            previous.get("threshold") != context.get("threshold")
            for previous in previous_contexts
        )
        same_unit_cross_minread = any(
            previous.get("unit_identity") == context.get("unit_identity")
            and previous.get("minread") != context.get("minread")
            for previous in previous_contexts
        )
        labels = []
        if same_unit_cross_minread:
            labels.append("SAME_UNIT_CROSS_MINREAD")
        elif cross_minread:
            labels.append("CROSS_MINREAD")
        if cross_threshold:
            labels.append("CROSS_THRESHOLD")
        if not labels:
            labels.append("SAME_STRUCTURAL_KEY")
        return "|".join(labels), cross_minread, cross_threshold, same_unit_cross_minread

    def get_or_enumerate(
        self,
        *,
        full_patterns: Sequence[str],
        partial_patterns: Sequence[str],
        k: int,
        structural_alt_universe_mask: int,
        exact_k_max: int,
        max_vertex_sets: int,
        solver_time_limit_seconds: float,
        universe_mode: str,
        context: Mapping[str, Any],
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        key = self._key(
            full_patterns=full_patterns,
            partial_patterns=partial_patterns,
            k=k,
            structural_alt_universe_mask=structural_alt_universe_mask,
            exact_k_max=exact_k_max,
            max_vertex_sets=max_vertex_sets,
            solver_time_limit_seconds=solver_time_limit_seconds,
            universe_mode=universe_mode,
        )
        key_sha256 = str(self._diagnostic_digest(key))
        known_digest_keys = self._digest_keys[key_sha256]
        if known_digest_keys and key not in known_digest_keys:
            self.digest_collision_guard_events += 1
        known_digest_keys.add(key)
        normalized_context = dict(context)
        self.lookups += 1
        if key in self._entries:
            self.hits += 1
            relation, cross_minread, cross_threshold, same_unit_cross_minread = (
                self._reuse_relation(self._contexts[key], normalized_context)
            )
            self.cross_minread_hits += int(cross_minread)
            self.cross_threshold_hits += int(cross_threshold)
            self.same_unit_cross_minread_hits += int(same_unit_cross_minread)
            self._contexts[key].append(normalized_context)
            return copy.deepcopy(self._entries[key]), {
                "structural_cache_outcome": "HIT_COMPLETE",
                "structural_cache_key_sha256": key_sha256,
                "structural_cache_reuse_relation": relation,
                "structural_solver_invoked": False,
            }
        self.misses += 1
        self.solver_invocations += 1
        result = self._enumerate(
            full_patterns,
            partial_patterns,
            k,
            universe_mode=universe_mode,
            max_sets=max_vertex_sets,
            time_limit_seconds=solver_time_limit_seconds,
        )
        self._contexts[key].append(normalized_context)
        if result.get("complete") is not True:
            self.rejected_incomplete += 1
            return result, {
                "structural_cache_outcome": "MISS_INCOMPLETE_NOT_STORED",
                "structural_cache_key_sha256": key_sha256,
                "structural_cache_reuse_relation": "NONE",
                "structural_solver_invoked": True,
            }
        if not self._complete_result_is_valid(
            result,
            full_patterns=full_patterns,
            partial_patterns=partial_patterns,
            k=k,
        ):
            raise RuntimeError("complete structural enumeration failed cache invariants")
        self._entries[key] = copy.deepcopy(result)
        self.stores_complete += 1
        return result, {
            "structural_cache_outcome": "MISS_STORED_COMPLETE",
            "structural_cache_key_sha256": key_sha256,
            "structural_cache_reuse_relation": "NONE",
            "structural_solver_invoked": True,
        }

    def diagnostics(self) -> dict[str, Any]:
        return {
            "schema_name": "intersubmod.structural_enumeration_cache_diagnostics",
            "schema_version": STRUCTURAL_CACHE_SCHEMA_VERSION,
            "enabled": True,
            "lifecycle": "PROCESS_LOCAL_ONE_RANKING_CHILD",
            "complete_only": True,
            "full_frozen_tuple_is_equality_authority": True,
            "diagnostic_sha256_is_not_equality_authority": True,
            "likelihood_cached": False,
            "lookups": self.lookups,
            "hits": self.hits,
            "misses": self.misses,
            "stores_complete": self.stores_complete,
            "rejected_incomplete": self.rejected_incomplete,
            "solver_invocations": self.solver_invocations,
            "entries_final": len(self._entries),
            "evictions": self.evictions,
            "cross_minread_hits": self.cross_minread_hits,
            "cross_threshold_hits": self.cross_threshold_hits,
            "same_unit_cross_minread_hits": self.same_unit_cross_minread_hits,
            "digest_collision_guard_events": self.digest_collision_guard_events,
        }


class _SemanticListHasher:
    """Incrementally reproduce ``semantic_digest(list(rows))`` without retaining rows."""

    def __init__(self) -> None:
        self._digest = hashlib.sha256()
        self._digest.update(b"[")
        self._first = True
        self.n_rows = 0

    def update(self, row: Mapping[str, Any]) -> None:
        if not self._first:
            self._digest.update(b",")
        self._digest.update(
            json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
        )
        self._first = False
        self.n_rows += 1

    def hexdigest(self) -> str:
        finalized = self._digest.copy()
        finalized.update(b"]")
        return finalized.hexdigest()


def _nearest_rank_index(n: int, probability: float) -> int:
    """Return the zero-based empirical nearest-rank index for a nonempty sample."""
    if n < 1 or not 0.0 < probability <= 1.0:
        raise ValueError("nearest-rank quantile requires n>=1 and probability in (0,1]")
    return min(n - 1, math.ceil(probability * n) - 1)


def summarize_runtime_values(values: Iterable[float]) -> dict[str, Any]:
    """Summarize non-reproducible wall timings with exact empirical quantiles.

    ``numpy.fromiter`` keeps the temporary representation compact (8 bytes per
    value).  In-place partition obtains the three exact nearest-rank order
    statistics without materializing Python-float lists or sorting every value.
    """
    array = np.fromiter((float(value) for value in values), dtype=np.float64)
    n = int(array.size)
    if n == 0:
        return {"n": 0, "sum": 0.0, "max": None, "p50": None, "p95": None, "p99": None}
    if not np.isfinite(array).all() or bool((array < 0.0).any()):
        raise ValueError("runtime diagnostics must be finite and nonnegative")
    total = math.fsum(float(value) for value in array)
    maximum = float(array.max())
    indices = tuple(sorted({_nearest_rank_index(n, value) for value in (0.50, 0.95, 0.99)}))
    array.partition(indices)
    return {
        "n": n,
        "sum": total,
        "max": maximum,
        "p50": float(array[_nearest_rank_index(n, 0.50)]),
        "p95": float(array[_nearest_rank_index(n, 0.95)]),
        "p99": float(array[_nearest_rank_index(n, 0.99)]),
    }


def summarize_runtime_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Create one exact timing summary without retaining a second row table."""
    return {
        "n_unit_evaluations": len(rows),
        **{
            metric: summarize_runtime_values(
                row["_runtime_diagnostics"][metric] for row in rows
            )
            for metric in RUNTIME_METRICS
        },
    }


def popcount(value: int) -> int:
    """Python 3.9-compatible integer population count."""
    return bin(int(value)).count("1")


def _one_file(input_dir: Path, suffix: str) -> Path:
    matches = sorted(input_dir.glob(f"*{suffix}"))
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one *{suffix} in {input_dir}, found {len(matches)}")
    return matches[0]


def _open_tsv(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", newline=""
    )


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".12g")
    return value


@dataclass(frozen=True)
class Component:
    dataset: str
    chrom: str
    threshold: int
    basis: str
    component_id: str
    start1: int
    end1: int
    k: int
    site_indices: tuple[int, ...]
    phase_set: str | None = None
    phase_set_status: str = "UNKNOWN"
    inference_role: str = "UNKNOWN"


@dataclass
class UnitEvidence:
    component: Component
    family: str
    pattern_counts: Counter[str] = field(default_factory=Counter)
    quality_pattern_counts: Counter[tuple[str, tuple[int, ...]]] = field(default_factory=Counter)
    projections: int = 0
    raw_cell_counts: Counter[str] = field(default_factory=Counter)
    implicit_uncovered_x: int = 0

    def add_projection(
        self,
        codes_by_site: Mapping[int, str],
        qualities_by_site: Mapping[int, int | None] | None = None,
    ) -> None:
        symbols: list[str] = []
        quality_signature: list[int] = []
        qualities_by_site = qualities_by_site or {}
        self.projections += 1
        for site_index in self.component.site_indices:
            raw = codes_by_site.get(site_index)
            if raw is None:
                self.implicit_uncovered_x += 1
                symbols.append("X")
                quality_signature.append(-1)
            else:
                if raw not in RAW_CODES:
                    raise ValueError(f"invalid sparse call code: {raw!r}")
                self.raw_cell_counts[raw] += 1
                symbols.append(raw if raw in FIXED_CODES else "X")
                if raw in FIXED_CODES:
                    quality = qualities_by_site.get(site_index)
                    if quality is None or quality < 0:
                        raise RuntimeError("fixed R/A call lacks a valid base quality")
                    quality_signature.append(int(quality))
                else:
                    quality_signature.append(-1)
        pattern = "".join(symbols)
        self.pattern_counts[pattern] += 1
        self.quality_pattern_counts[(pattern, tuple(quality_signature))] += 1


def _pattern_rows(
    unit: UnitEvidence, minread: int, *, minread_role: str = "PRIMARY"
) -> list[dict[str, Any]]:
    rows = []
    component = unit.component
    structural_universe = _structural_alt_universe(unit, minread)
    for (pattern, qualities), count in sorted(unit.quality_pattern_counts.items()):
        symbolic = SymbolicPattern.from_string(pattern, count)
        informative = symbolic.fixed_mask != 0
        total_pattern_count = unit.pattern_counts[pattern]
        retained = informative and total_pattern_count >= minread
        core = {
            "dataset": component.dataset,
            "chrom": component.chrom,
            "threshold": component.threshold,
            "component_basis": component.basis,
            "phase_set": component.phase_set,
            "phase_set_status": component.phase_set_status,
            "inference_role": component.inference_role,
            "component_id": component.component_id,
            "family": unit.family,
            "structural_exact_pattern_minread": minread,
            "structural_minread_role": minread_role,
            "k": component.k,
            "pattern_rax": pattern,
            "fixed_base_qualities": ",".join("" if value < 0 else str(value) for value in qualities),
            "n_molecules": count,
            "rax_pattern_total_molecules": total_pattern_count,
            "scoring_eligible": informative,
            "structural_retained": retained,
            "exclusion_reason": "" if informative else "ALL_X_UNINFORMATIVE_CONSTANT_LIKELIHOOD",
            "n_free": symbolic.n_free,
            "n_conceptual_completions": symbolic.n_completions(),
            "n_effective_free_in_structural_alt_universe": popcount(
                symbolic.free_mask & structural_universe
            ),
            "n_effective_completions": (
                0
                if symbolic.alt_mask & ~structural_universe
                else 1 << popcount(symbolic.free_mask & structural_universe)
            ),
            "effective_completion_status": (
                "FIXED_ALT_OUTSIDE_STRUCTURAL_UNIVERSE"
                if symbolic.alt_mask & ~structural_universe
                else "COMPATIBLE_WITH_STRUCTURAL_OBSERVED_ALT_UNIVERSE"
            ),
        }
        core["pattern_semantic_sha256"] = semantic_digest(core)
        rows.append(core)
    return rows


def _base_unit_row(
    unit: UnitEvidence, minread: int, *, minread_role: str = "PRIMARY"
) -> dict[str, Any]:
    component = unit.component
    informative = sum(count for pattern, count in unit.pattern_counts.items() if set(pattern) != {"X"})
    all_x = unit.pattern_counts.get("X" * component.k, 0)
    retained = sum(
        count for pattern, count in unit.pattern_counts.items()
        if set(pattern) != {"X"} and count >= minread
    )
    structural_groups = sum(
        set(pattern) != {"X"} and count >= minread for pattern, count in unit.pattern_counts.items()
    )
    scoring_groups = sum(set(pattern) != {"X"} for pattern in unit.pattern_counts)
    projected_cells = unit.projections * component.k
    explicit_cells = sum(unit.raw_cell_counts.values())
    if explicit_cells + unit.implicit_uncovered_x != projected_cells:
        raise AssertionError("projected molecule-cell mass does not conserve")
    if informative + all_x != unit.projections:
        raise AssertionError("projection funnel does not conserve")
    if retained + (informative - retained) != informative:
        raise AssertionError("structural/scoring molecule mass does not conserve")
    return {
        "dataset": component.dataset,
        "chrom": component.chrom,
        "threshold": component.threshold,
        "component_basis": component.basis,
        "phase_set": component.phase_set,
        "phase_set_status": component.phase_set_status,
        "inference_role": component.inference_role,
        "component_id": component.component_id,
        "family": unit.family,
        "structural_exact_pattern_minread": minread,
        "structural_minread_role": minread_role,
        "start1": component.start1,
        "end1": component.end1,
        "k": component.k,
        "k_component_sites": component.k,
        "molecule_component_projections": unit.projections,
        "informative_scoring_molecules": informative,
        "all_x_excluded_molecules": all_x,
        "n_scoring_pattern_groups": scoring_groups,
        "n_scoring_quality_groups": sum(
            set(pattern) != {"X"} for pattern, _ in unit.quality_pattern_counts
        ),
        "structural_retained_molecules": retained,
        "below_minread_scoring_molecules": informative - retained,
        "n_structural_pattern_groups": structural_groups,
        "projected_cells": projected_cells,
        "fixed_ra_cells": unit.raw_cell_counts["R"] + unit.raw_cell_counts["A"],
        "alt_cells": unit.raw_cell_counts["A"],
        "other_cells": unit.raw_cell_counts["O"],
        "deletion_cells": unit.raw_cell_counts["D"],
        "refskip_cells": unit.raw_cell_counts["S"],
        "low_baseq_cells": unit.raw_cell_counts["L"],
        "explicit_x_cells": unit.raw_cell_counts["X"],
        "implicit_uncovered_x_cells": unit.implicit_uncovered_x,
        "non_ra_cells_marginalized": projected_cells
        - unit.raw_cell_counts["R"]
        - unit.raw_cell_counts["A"],
        "bq_scoring_molecules": informative,
    }


def _empty_rank_fields(status: str, generation_status: str) -> dict[str, Any]:
    return {
        "k_observed_alt_active": 0,
        "k_scoring_alt_observed": 0,
        "n_not_structural_alt_active_sites": None,
        "minimum_hidden_nodes": None,
        "raw_tree_candidates_T": None,
        "distinct_vertex_sets_V": None,
        "candidate_vertex_sets_enumerated": 0,
        "candidate_generation_complete": False,
        "candidate_generation_status": generation_status,
        "best_vertex_sets": None,
        "top_edge_variants": None,
        "best_log_likelihood": None,
        "second_log_likelihood": None,
        "delta_best_second": None,
        "top_relative_likelihood_weight": None,
        "selection_status": status,
        "all_fits_converged": None,
        "all_fits_monotone": None,
        "max_emission_rank": None,
        "vertex_set_ids": "[]",
        "top_vertex_set_ids": "[]",
        "topology_classes": "[]",
        "n_topology_classes": None,
        "coarse_topology_class_unique": None,
        "parent_edge_assignment_unique": None,
        "exact_topology_unique": None,
        "exact_topology_uniqueness_status": "NOT_EVALUATED_WITHOUT_COMPLETE_WINNING_VERTEX_SET",
        "topology_class_by_top_vertex_set": "{}",
        "topology_derivation_status": "NOT_APPLICABLE_WITHOUT_COMPLETE_WINNING_VERTEX_SET",
        "fixed_error_grid_top_vertex_set_ids": "{}",
        "fixed_error_grid_stable_with_quality_primary": None,
        "fixed_error_grid_all_converged": None,
        "structural_partial_pattern_groups": 0,
        "partial_group_coverage_denominator": 0,
        "partial_groups_covered": 0,
        "partial_groups_unsatisfied": 0,
        "conditional_candidate_ranking_bootstrap_status": "NOT_APPLICABLE_NO_COMPLETE_MULTI_VERTEX_CANDIDATE",
        "conditional_candidate_ranking_bootstrap_replicates": 0,
        "conditional_candidate_ranking_bootstrap_seed": None,
        "conditional_candidate_ranking_bootstrap_top_vertex_set_frequency": "{}",
        "conditional_candidate_ranking_bootstrap_primary_top_set_frequency": None,
        "conditional_candidate_ranking_bootstrap_all_converged": None,
    }


def _selection_status(raw_t: int, n_vertex_sets: int, n_top: int, top_edges: int) -> str:
    if raw_t == 1:
        return "T1_CANDIDATE_UNIQUE"
    if n_vertex_sets == 1:
        return "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED"
    if n_top > 1:
        return "LIKELIHOOD_TIED_VERTEX_SETS"
    if top_edges > 1:
        return "LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED"
    return "LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE"


def _quality_emission_matrix(
    quality_pattern_counts: Sequence[tuple[str, tuple[int, ...], int]],
    vertices: Sequence[int],
    *,
    minimum_error_rate: float,
    maximum_error_rate: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return per-quality-group emissions; non-R/A calls are marginalized.

    ``e=10**(-BQ/10)`` is the probability of any wrong base.  Under an explicitly
    symmetric substitution model, conditioning on the observed base being one
    of {REF, ALT} gives match=(1-e)/(1-2e/3) and
    flip=(e/3)/(1-2e/3).  Treating R<->A flip as ``e`` would overstate the error
    probability three-fold.  O/D/S/L/X remain marginalized missing calls.
    """
    if not 0 < minimum_error_rate <= maximum_error_rate < 0.5:
        raise ValueError("invalid quality-error bounds")
    matrix = np.ones((len(quality_pattern_counts), len(vertices)), dtype=float)
    counts = np.empty(len(quality_pattern_counts), dtype=float)
    for row_index, (pattern, qualities, count) in enumerate(quality_pattern_counts):
        if len(pattern) != len(qualities) or count <= 0:
            raise ValueError("invalid quality-aware pattern group")
        counts[row_index] = count
        for col_index, vertex in enumerate(vertices):
            probability = 1.0
            for bit, (symbol, quality) in enumerate(zip(pattern, qualities)):
                if symbol == "X":
                    if quality != -1:
                        raise ValueError("missing symbol has a fixed-call base quality")
                    continue
                if symbol not in FIXED_CODES or quality < 0:
                    raise ValueError("fixed symbol lacks a base quality")
                error = min(max(10.0 ** (-quality / 10.0), minimum_error_rate), maximum_error_rate)
                expected_alt = bool(vertex & (1 << bit))
                observed_alt = symbol == "A"
                denominator = 1.0 - (2.0 * error / 3.0)
                match_probability = (1.0 - error) / denominator
                flip_probability = (error / 3.0) / denominator
                probability *= match_probability if expected_alt == observed_alt else flip_probability
            matrix[row_index, col_index] = probability
    return matrix, counts


def fit_quality_aware_mixture(
    quality_pattern_counts: Sequence[tuple[str, tuple[int, ...], int]],
    vertices: Sequence[int],
    *,
    minimum_error_rate: float = 1e-6,
    maximum_error_rate: float = 0.25,
    tolerance: float = 1e-10,
    max_iterations: int = 2_000,
) -> MixtureFit:
    """Fit a certified simplex mixture using per-call Phred base qualities."""
    if not vertices:
        raise ValueError("vertices must not be empty")
    informative = [entry for entry in quality_pattern_counts if set(entry[0]) != {"X"}]
    if not informative:
        weights = tuple(1.0 / len(vertices) for _ in vertices)
        return MixtureFit(
            tuple(vertices),
            weights,
            0.0,
            True,
            0,
            1,
            True,
            optimizer_status="UNINFORMATIVE_ALL_X",
            global_log_likelihood_gap_bound=0.0,
            simplex_residual=0.0,
            augmented_emission_rank=1,
            mixture_weights_identifiable=len(vertices) == 1,
        )
    vertices = tuple(sorted(int(value) for value in vertices))
    q, counts = _quality_emission_matrix(
        informative,
        vertices,
        minimum_error_rate=minimum_error_rate,
        maximum_error_rate=maximum_error_rate,
    )
    return fit_emission_mixture_certified(
        q,
        counts,
        vertices,
        tolerance=tolerance,
        max_iterations=max_iterations,
        certificate_tolerance=OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE,
    )


def _rank_fits(
    fits: Sequence[tuple[tuple[int, ...], MixtureFit]],
    tie_tolerance: float,
) -> tuple[list[tuple[int, ...]], float, float | None, float | None, float, bool, bool, int]:
    ordered = sorted(fits, key=lambda item: (-item[1].log_likelihood, item[0]))
    if not ordered:
        raise ValueError("no likelihood fits")
    best_ll = ordered[0][1].log_likelihood
    top = [vertices for vertices, fit in ordered if best_ll - fit.log_likelihood <= tie_tolerance]
    second_ll = ordered[1][1].log_likelihood if len(ordered) > 1 else None
    delta = None if second_ll is None else best_ll - second_ll
    relative = [math.exp(max(-745.0, fit.log_likelihood - best_ll)) for _, fit in ordered]
    top_weight = relative[0] / math.fsum(relative)
    return (
        top,
        best_ll,
        second_ll,
        delta,
        top_weight,
        all(fit.converged for _, fit in ordered),
        all(fit.monotone for _, fit in ordered),
        max(fit.emission_rank for _, fit in ordered),
    )


def conditional_candidate_ranking_bootstrap(
    quality_pattern_counts: Sequence[tuple[str, tuple[int, ...], int]],
    vertex_sets: Sequence[tuple[int, ...]],
    *,
    k: int,
    primary_top_vertices: Sequence[tuple[int, ...]],
    replicates: int,
    seed: int,
    minimum_error_rate: float,
    maximum_error_rate: float,
    tie_tolerance: float,
) -> dict[str, Any]:
    """Conditional-on-fixed-candidate-set molecule ranking bootstrap.

    Multinomial resampling is exactly equivalent to resampling individual
    molecules with replacement while avoiding an expanded molecule list.
    """
    if replicates < 1:
        return {
            "status": "NOT_RUN_REQUESTED_0",
            "replicates": 0,
            "seed": seed,
            "top_vertex_set_frequency": {},
            "primary_top_set_frequency": None,
            "all_converged": None,
        }
    counts = np.asarray([count for _, _, count in quality_pattern_counts], dtype=np.int64)
    total = int(counts.sum())
    if total < 1:
        raise ValueError("conditional ranking bootstrap requires informative molecules")
    probabilities = counts.astype(float) / total
    rng = np.random.default_rng(seed)
    inclusion = Counter()
    primary_ids = tuple(sorted(vertex_set_digest(k, values) for values in primary_top_vertices))
    primary_matches = 0
    all_converged = True
    for _ in range(replicates):
        sampled = rng.multinomial(total, probabilities)
        sampled_groups = [
            (pattern, qualities, int(count))
            for (pattern, qualities, _), count in zip(quality_pattern_counts, sampled)
            if count > 0
        ]
        fits = [
            (
                values,
                fit_quality_aware_mixture(
                    sampled_groups,
                    values,
                    minimum_error_rate=minimum_error_rate,
                    maximum_error_rate=maximum_error_rate,
                ),
            )
            for values in vertex_sets
        ]
        top, _, _, _, _, converged, monotone, _ = _rank_fits(fits, tie_tolerance)
        all_converged = all_converged and converged and monotone
        ids = tuple(sorted(vertex_set_digest(k, values) for values in top))
        for vertex_set_id in ids:
            inclusion[vertex_set_id] += 1
        primary_matches += ids == primary_ids
    return {
        "status": "COMPLETE" if all_converged else "ABSTAIN_NONCONVERGENCE",
        "replicates": replicates,
        "seed": seed,
        "top_vertex_set_frequency": {
            key: value / replicates for key, value in sorted(inclusion.items())
        },
        "primary_top_set_frequency": primary_matches / replicates,
        "all_converged": all_converged,
    }


def topology_classes_for_vertex_set(vertices: Sequence[int]) -> tuple[str, ...]:
    """Derive all coarse parent-edge shape classes without edge enumeration.

    Depth is fixed by Hamming weight because every directed edge removes exactly
    one mutation.  Branch-free feasibility is a bipartite matching problem in
    which every selected child chooses a distinct selected predecessor.  Branch
    feasibility only requires one selected predecessor shared by two possible
    children.  These two checks completely determine the four declared classes.
    """
    selected = set(int(value) for value in vertices)
    if 0 not in selected:
        raise ValueError("root vertex is required")
    children = sorted(selected - {0})
    predecessors: dict[int, tuple[int, ...]] = {}
    possible_child_counts = Counter()
    for child in children:
        parents = tuple(
            child ^ (1 << bit)
            for bit in range(child.bit_length())
            if child & (1 << bit) and (child ^ (1 << bit)) in selected
        )
        if not parents:
            raise ValueError("vertex set is not rooted by Hamming-1 predecessors")
        predecessors[child] = parents
        possible_child_counts.update(parents)

    parent_match: dict[int, int] = {}

    def assign(child: int, seen: set[int]) -> bool:
        for parent in predecessors[child]:
            if parent in seen:
                continue
            seen.add(parent)
            current = parent_match.get(parent)
            if current is None or assign(current, seen):
                parent_match[parent] = child
                return True
        return False

    no_branch_feasible = all(assign(child, set()) for child in children)
    branch_feasible = any(count >= 2 for count in possible_child_counts.values())
    direct = any(bin(child).count("1") >= 2 for child in children)
    branch_states = []
    if no_branch_feasible:
        branch_states.append(False)
    if branch_feasible:
        branch_states.append(True)
    if not branch_states:
        raise AssertionError("no parent-edge branching state is feasible")
    classes = []
    for branch in branch_states:
        if branch and direct:
            classes.append("sister+direct")
        elif branch:
            classes.append("sister-only")
        elif direct:
            classes.append("direct-only")
        else:
            classes.append("single-only")
    order = {"single-only": 0, "sister-only": 1, "direct-only": 2, "sister+direct": 3}
    return tuple(sorted(set(classes), key=order.__getitem__))


def _state_string(vertex: int, k: int) -> str:
    return "".join("A" if vertex & (1 << bit) else "R" for bit in range(k))


def _unique_parent_edges(vertices: Sequence[int]) -> list[dict[str, int]] | None:
    selected = set(vertices)
    edges = []
    for child in sorted(selected - {0}):
        parents = [
            child ^ (1 << bit)
            for bit in range(child.bit_length())
            if child & (1 << bit) and (child ^ (1 << bit)) in selected
        ]
        if len(parents) != 1:
            return None
        edges.append({"parent": parents[0], "child": child})
    return edges


def _candidate_state_records(
    vertices: Sequence[int], full: Sequence[str], partial: Sequence[str], k: int
) -> list[dict[str, Any]]:
    full_states = {sum((symbol == "A") << bit for bit, symbol in enumerate(pattern)) for pattern in full}
    partial_groups = [SymbolicPattern.from_string(pattern) for pattern in partial]
    records = []
    for vertex in vertices:
        roles = []
        if vertex == 0:
            roles.append("root")
        if vertex in full_states:
            roles.append("full-observed")
        if any(group.compatible(vertex) for group in partial_groups):
            roles.append("partial-compatible")
        if vertex != 0 and vertex not in full_states and not any(
            group.compatible(vertex) for group in partial_groups
        ):
            roles.append("connector")
        records.append({"bitmask": vertex, "state_rax": _state_string(vertex, k), "roles": roles})
    return records


def rank_unit(
    unit: UnitEvidence,
    *,
    minread: int = 3,
    minread_role: str = "PRIMARY",
    exact_k_max: int = 12,
    max_vertex_sets: int = 256,
    solver_time_limit_seconds: float = 30.0,
    fixed_error_grid: Sequence[float] = (0.005, 0.01, 0.02, 0.05),
    minimum_bq_error_rate: float = 1e-6,
    maximum_bq_error_rate: float = 0.25,
    conditional_candidate_ranking_bootstrap_replicates: int = 0,
    conditional_candidate_ranking_bootstrap_seed: int = 20260716,
    tie_tolerance: float = 1e-6,
    structural_enumeration_cache: StructuralEnumerationCache | None = None,
    candidate_record_sink: Callable[[dict[str, Any]], None] | None = None,
    responsibility_record_sink: Callable[[dict[str, Any]], None] | None = None,
    retain_detail_records: bool = True,
) -> dict[str, Any]:
    """Rank one component×HP unit without expanding partial completions."""
    unit_started_ns = _monotonic_ns()
    candidate_generation_elapsed_ns = 0
    likelihood_fit_elapsed_ns = 0
    candidate_generation_invoked = False
    likelihood_fit_invoked = False
    structural_cache_diagnostics = {
        "structural_cache_outcome": "NOT_ELIGIBLE",
        "structural_cache_key_sha256": "",
        "structural_cache_reuse_relation": "NONE",
        "structural_solver_invoked": False,
    }
    if (
        minread < 1
        or exact_k_max < 1
        or max_vertex_sets < 1
        or conditional_candidate_ranking_bootstrap_replicates < 0
    ):
        raise ValueError("invalid structural or solver limits")
    row = _base_unit_row(unit, minread, minread_role=minread_role)
    row["_candidate_records"] = []
    row["_responsibility_records"] = []
    component = unit.component
    scoring_counts = [
        (pattern, count) for pattern, count in sorted(unit.pattern_counts.items())
        if set(pattern) != {"X"}
    ]
    quality_counts = [
        (pattern, qualities, count)
        for (pattern, qualities), count in sorted(unit.quality_pattern_counts.items())
        if set(pattern) != {"X"}
    ]
    structural = [
        pattern for pattern, count in sorted(unit.pattern_counts.items())
        if set(pattern) != {"X"} and count >= minread
    ]
    structural_alt_universe = 0
    scoring_alt_universe = 0
    for pattern in unit.pattern_counts:
        if set(pattern) == {"X"}:
            continue
        for bit, symbol in enumerate(pattern):
            if symbol == "A":
                scoring_alt_universe |= 1 << bit
    for pattern in structural:
        for bit, symbol in enumerate(pattern):
            if symbol == "A":
                structural_alt_universe |= 1 << bit
    effective_k = popcount(structural_alt_universe)
    row["k_observed_alt_active"] = effective_k
    row["k_scoring_alt_observed"] = popcount(scoring_alt_universe)
    row["n_not_structural_alt_active_sites"] = component.k - effective_k
    structural_partial = [pattern for pattern in structural if "X" in pattern]
    row["structural_partial_pattern_groups"] = len(structural_partial)
    if not scoring_counts:
        row.update(_empty_rank_fields("NO_INFORMATIVE_SCORING_MOLECULE", "NO_INFORMATIVE_PATTERN"))
    elif effective_k > exact_k_max:
        row.update(_empty_rank_fields("LOCAL_ONLY_OR_ABSTAIN_K_GT_EXACT_LIMIT", "K_GT_EXACT_LIMIT"))
        row["k_observed_alt_active"] = effective_k
        row["k_scoring_alt_observed"] = popcount(scoring_alt_universe)
        row["n_not_structural_alt_active_sites"] = component.k - effective_k
        row["structural_partial_pattern_groups"] = len(structural_partial)
    elif not structural or not any("A" in pattern for pattern in structural):
        row.update(
            _empty_rank_fields(
                "STRUCTURE_ABSTAIN_NO_MINREAD_ALT_PATTERN",
                "NO_STRUCTURAL_ALT_PATTERN_AT_DECLARED_MINREAD",
            )
        )
    else:
        full = [pattern for pattern in structural if "X" not in pattern]
        partial = structural_partial
        candidate_generation_invoked = True
        candidate_generation_started_ns = _monotonic_ns()
        if structural_enumeration_cache is None:
            enumeration = enumerate_optimal_vertex_sets(
                full,
                partial,
                component.k,
                universe_mode="observed_alt",
                max_sets=max_vertex_sets,
                time_limit_seconds=solver_time_limit_seconds,
            )
            structural_cache_diagnostics = {
                "structural_cache_outcome": "DISABLED",
                "structural_cache_key_sha256": "",
                "structural_cache_reuse_relation": "NONE",
                "structural_solver_invoked": True,
            }
        else:
            enumeration, structural_cache_diagnostics = (
                structural_enumeration_cache.get_or_enumerate(
                    full_patterns=full,
                    partial_patterns=partial,
                    k=component.k,
                    structural_alt_universe_mask=structural_alt_universe,
                    exact_k_max=exact_k_max,
                    max_vertex_sets=max_vertex_sets,
                    solver_time_limit_seconds=solver_time_limit_seconds,
                    universe_mode="observed_alt",
                    context={
                        "minread": int(minread),
                        "threshold": int(component.threshold),
                        "unit_identity": (
                            component.dataset,
                            component.chrom,
                            component.basis,
                            component.phase_set,
                            component.component_id,
                            unit.family,
                        ),
                    },
                )
            )
        candidate_generation_elapsed_ns += _monotonic_ns() - candidate_generation_started_ns
        vertex_sets = [tuple(sorted(int(value) for value in values)) for values in enumeration.get("vertex_sets", [])]
        complete = bool(enumeration.get("complete"))
        first = enumeration.get("first") or {}
        if not complete:
            generation = str(enumeration.get("stop_status") or first.get("status") or "INCOMPLETE")
            row.update(
                _empty_rank_fields(
                    "RANK_ABSTAIN_VERTEX_ENUMERATION_INCOMPLETE",
                    f"EXACT_ENUMERATION_INCOMPLETE:{generation}",
                )
            )
            row["minimum_hidden_nodes"] = enumeration.get("objective", first.get("objective"))
            row["candidate_vertex_sets_enumerated"] = len(vertex_sets)
            row["vertex_set_ids"] = json.dumps(
                [vertex_set_digest(component.k, values) for values in vertex_sets], separators=(",", ":")
            )
        else:
            if not vertex_sets:
                raise RuntimeError("complete exact enumeration returned no vertex sets")
            covered = sum(
                all(
                    any(SymbolicPattern.from_string(pattern).compatible(vertex) for vertex in values)
                    for values in vertex_sets
                )
                for pattern in partial
            )
            row["partial_group_coverage_denominator"] = len(partial)
            row["partial_groups_covered"] = covered
            row["partial_groups_unsatisfied"] = len(partial) - covered
            if covered != len(partial):
                raise AssertionError("enumerated minimum vertex set violates a partial group constraint")
            edge_counts = {values: parent_choice_count(values) for values in vertex_sets}
            if any(count < 1 for count in edge_counts.values()):
                raise AssertionError("enumerated vertex set has no rooted parent assignment")
            raw_t = sum(edge_counts.values())
            ids = [vertex_set_digest(component.k, values) for values in vertex_sets]
            likelihood_fit_invoked = True
            likelihood_fit_started_ns = _monotonic_ns()
            fits = [
                (
                    values,
                    fit_quality_aware_mixture(
                        quality_counts,
                        values,
                        minimum_error_rate=minimum_bq_error_rate,
                        maximum_error_rate=maximum_bq_error_rate,
                    ),
                )
                for values in vertex_sets
            ]
            likelihood_fit_elapsed_ns += _monotonic_ns() - likelihood_fit_started_ns
            if len(vertex_sets) > 1:
                (
                    top_vertices,
                    best_ll,
                    second_ll,
                    delta,
                    top_weight,
                    all_converged,
                    all_monotone,
                    max_rank,
                ) = _rank_fits(fits, tie_tolerance)
            else:
                top_vertices = list(vertex_sets)
                best_ll = second_ll = delta = None
                top_weight = 1.0
                all_converged = fits[0][1].converged
                all_monotone = fits[0][1].monotone
                max_rank = fits[0][1].emission_rank
            primary_top_ids = [vertex_set_digest(component.k, values) for values in top_vertices]
            grid_top_ids: dict[str, list[str]] = {}
            grid_converged = True
            for fixed_error in fixed_error_grid:
                if not 0 < fixed_error < 0.5:
                    raise ValueError("fixed-error sensitivity values must lie in (0, 0.5)")
                if len(vertex_sets) == 1:
                    fixed_top = list(vertex_sets)
                    fixed_converged = True
                else:
                    likelihood_fit_started_ns = _monotonic_ns()
                    fixed_fits = [
                        (
                            values,
                            fit_vertex_mixture_slsqp(scoring_counts, values, error_rate=fixed_error),
                        )
                        for values in vertex_sets
                    ]
                    likelihood_fit_elapsed_ns += _monotonic_ns() - likelihood_fit_started_ns
                    fixed_top, _, _, _, _, fixed_converged, fixed_monotone, _ = _rank_fits(
                        fixed_fits, tie_tolerance
                    )
                    fixed_converged = fixed_converged and fixed_monotone
                grid_converged = grid_converged and fixed_converged
                grid_top_ids[format(fixed_error, ".12g")] = [
                    vertex_set_digest(component.k, values) for values in fixed_top
                ]
            grid_stable = grid_converged and all(
                ids_at_error == primary_top_ids for ids_at_error in grid_top_ids.values()
            )
            if len(vertex_sets) > 1 and all_converged and all_monotone:
                unit_seed_payload = (
                    f"{component.dataset}\0{component.chrom}\0{component.basis}\0"
                    f"{component.threshold}\0{component.component_id}\0{unit.family}"
                )
                derived_seed = (
                    int(hashlib.sha256(unit_seed_payload.encode()).hexdigest()[:16], 16)
                    ^ int(conditional_candidate_ranking_bootstrap_seed)
                ) % (2**32)
                likelihood_fit_started_ns = _monotonic_ns()
                bootstrap = conditional_candidate_ranking_bootstrap(
                    quality_counts,
                    vertex_sets,
                    k=component.k,
                    primary_top_vertices=top_vertices,
                    replicates=conditional_candidate_ranking_bootstrap_replicates,
                    seed=derived_seed,
                    minimum_error_rate=minimum_bq_error_rate,
                    maximum_error_rate=maximum_bq_error_rate,
                    tie_tolerance=tie_tolerance,
                )
                likelihood_fit_elapsed_ns += _monotonic_ns() - likelihood_fit_started_ns
            elif len(vertex_sets) == 1:
                bootstrap = {
                    "status": "NOT_APPLICABLE_V1",
                    "replicates": 0,
                    "seed": None,
                    "top_vertex_set_frequency": {},
                    "primary_top_set_frequency": None,
                    "all_converged": None,
                }
            else:
                bootstrap = {
                    "status": "NOT_RUN_PRIMARY_NONCONVERGENCE",
                    "replicates": 0,
                    "seed": None,
                    "top_vertex_set_frequency": {},
                    "primary_top_set_frequency": None,
                    "all_converged": None,
                }
            top_edges = sum(edge_counts[values] for values in top_vertices)
            status = _selection_status(raw_t, len(vertex_sets), len(top_vertices), top_edges)
            if not all_converged or not all_monotone:
                status = "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE"
            topology_by_vertex = {
                vertex_set_digest(component.k, values): list(topology_classes_for_vertex_set(values))
                for values in top_vertices
            }
            topology_classes = sorted(
                {shape for shapes in topology_by_vertex.values() for shape in shapes},
                key={
                    "single-only": 0,
                    "sister-only": 1,
                    "direct-only": 2,
                    "sister+direct": 3,
                }.__getitem__,
            )
            topology_status = (
                "ANALYTICAL_COMPLETE_OVER_ALL_WINNING_VERTEX_SETS_AND_PARENT_CHOICES"
                if status != "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE"
                else "ABSTAIN_PRIMARY_LIKELIHOOD_NONCONVERGENCE"
            )
            topology_reliable = status != "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE"
            if not topology_reliable:
                topology_by_vertex = {}
                topology_classes = []
            best_candidate_ll = max(fit.log_likelihood for _, fit in fits)
            likelihood_normalizer = math.fsum(
                math.exp(max(-745.0, fit.log_likelihood - best_candidate_ll)) for _, fit in fits
            )
            top_vertex_tuples = set(top_vertices)
            candidate_records = []
            responsibility_records = []
            for values, fit in fits:
                vertex_id = vertex_set_digest(component.k, values)
                candidate = {
                    "dataset": component.dataset,
                    "chrom": component.chrom,
                    "component_basis": component.basis,
                    "phase_set": component.phase_set,
                    "threshold": component.threshold,
                    "component_id": component.component_id,
                    "family": unit.family,
                    "structural_exact_pattern_minread": minread,
                    "vertex_set_id": vertex_id,
                    "states_json": json.dumps(
                        _candidate_state_records(values, full, partial, component.k),
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                    "parent_choice_count": edge_counts[values],
                    "unique_parent_edges_json": json.dumps(
                        _unique_parent_edges(values), sort_keys=True, separators=(",", ":")
                    ) if edge_counts[values] == 1 else "",
                    "primary_log_likelihood": fit.log_likelihood,
                    "relative_likelihood_weight": math.exp(
                        max(-745.0, fit.log_likelihood - best_candidate_ll)
                    ) / likelihood_normalizer,
                    "mixture_weights_json": json.dumps(
                        {str(vertex): weight for vertex, weight in zip(fit.vertices, fit.weights)},
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                    "fit_converged": fit.converged,
                    "fit_monotone": fit.monotone,
                    "optimizer_status": fit.optimizer_status,
                    "slsqp_success": fit.slsqp_success,
                    "slsqp_status": fit.slsqp_status,
                    "slsqp_message": fit.slsqp_message,
                    "warm_start_log_likelihood": fit.warm_start_log_likelihood,
                    "warm_start_global_log_likelihood_gap_bound": (
                        fit.warm_start_global_log_likelihood_gap_bound
                    ),
                    "refinement_iterations": fit.refinement_iterations,
                    "global_log_likelihood_gap_bound": fit.global_log_likelihood_gap_bound,
                    "simplex_residual": fit.simplex_residual,
                    "augmented_emission_rank": fit.augmented_emission_rank,
                    "mixture_weights_identifiable": fit.mixture_weights_identifiable,
                    "coarse_topology_classes_json": json.dumps(
                        topology_classes_for_vertex_set(values), separators=(",", ":")
                    ),
                    "is_winner": values in top_vertex_tuples and topology_reliable,
                    "is_tied_winner": len(top_vertices) > 1 and values in top_vertex_tuples and topology_reliable,
                }
                candidate["candidate_semantic_sha256"] = semantic_digest(candidate)
                if retain_detail_records:
                    candidate_records.append(candidate)
                if candidate_record_sink is not None:
                    candidate_record_sink(candidate)
                if (
                    values in top_vertex_tuples
                    and topology_reliable
                    and (retain_detail_records or responsibility_record_sink is not None)
                ):
                    emission, _ = _quality_emission_matrix(
                        quality_counts,
                        fit.vertices,
                        minimum_error_rate=minimum_bq_error_rate,
                        maximum_error_rate=maximum_bq_error_rate,
                    )
                    weights = np.asarray(fit.weights, dtype=float)
                    pattern_molecule_counts: Counter[str] = Counter()
                    responsibility_mass: dict[tuple[str, int], float] = defaultdict(float)
                    for group_index, (pattern, qualities, count) in enumerate(quality_counts):
                        numerator = emission[group_index] * weights
                        denominator = float(numerator.sum())
                        if denominator <= 0:
                            raise FloatingPointError("invalid posterior responsibility denominator")
                        pattern_molecule_counts[pattern] += count
                        for state, responsibility in zip(fit.vertices, numerator / denominator):
                            responsibility_mass[(pattern, state)] += count * float(responsibility)
                    for (pattern, state), weighted_mass in sorted(responsibility_mass.items()):
                        count = pattern_molecule_counts[pattern]
                        record = {
                            "dataset": component.dataset,
                            "chrom": component.chrom,
                            "component_basis": component.basis,
                            "phase_set": component.phase_set,
                            "threshold": component.threshold,
                            "component_id": component.component_id,
                            "family": unit.family,
                            "structural_exact_pattern_minread": minread,
                            "vertex_set_id": vertex_id,
                            "pattern_rax": pattern,
                            "fixed_base_qualities": "MARGINALIZED_OVER_BQ_QUALITY_GROUPS",
                            "n_molecules": count,
                            "state_bitmask": state,
                            "state_rax": _state_string(state, component.k),
                            "posterior_responsibility": weighted_mass / count,
                            "interpretation": (
                                "BQ_WEIGHTED_LATENT_EXPECTATION_NOT_HARD_READ_OR_CLONE_ASSIGNMENT"
                            ),
                        }
                        record["responsibility_semantic_sha256"] = semantic_digest(record)
                        if retain_detail_records:
                            responsibility_records.append(record)
                        if responsibility_record_sink is not None:
                            responsibility_record_sink(record)
                    for pattern, count in pattern_molecule_counts.items():
                        responsibility_sum = math.fsum(
                            weighted_mass / count
                            for (current_pattern, _state), weighted_mass in responsibility_mass.items()
                            if current_pattern == pattern
                        )
                        if abs(responsibility_sum - 1.0) > 1e-8:
                            raise AssertionError(
                                "posterior responsibilities do not conserve per pattern/candidate"
                            )
            if abs(math.fsum(
                math.exp(max(-745.0, fit.log_likelihood - best_candidate_ll))
                / likelihood_normalizer
                for _values, fit in fits
            ) - 1.0) > 1e-8:
                raise AssertionError("candidate relative likelihood weights do not conserve per unit")
            row["_candidate_records"] = candidate_records
            row["_responsibility_records"] = responsibility_records
            row.update(
                {
                    "minimum_hidden_nodes": int(enumeration["objective"]),
                    "raw_tree_candidates_T": raw_t,
                    "distinct_vertex_sets_V": len(vertex_sets),
                    "candidate_vertex_sets_enumerated": len(vertex_sets),
                    "candidate_generation_complete": True,
                    "candidate_generation_status": "EXACT_OPTIMAL_VERTEX_SETS_COMPLETE",
                    "best_vertex_sets": len(top_vertices),
                    "top_edge_variants": top_edges,
                    "best_log_likelihood": best_ll,
                    "second_log_likelihood": second_ll,
                    "delta_best_second": delta,
                    "top_relative_likelihood_weight": top_weight,
                    "selection_status": status,
                    "all_fits_converged": all_converged,
                    "all_fits_monotone": all_monotone,
                    "max_emission_rank": max_rank,
                    "vertex_set_ids": json.dumps(ids, separators=(",", ":")),
                    "top_vertex_set_ids": json.dumps(
                        primary_top_ids, separators=(",", ":")
                    ),
                    "topology_classes": json.dumps(topology_classes, separators=(",", ":")),
                    "n_topology_classes": len(topology_classes),
                    "coarse_topology_class_unique": (len(topology_classes) == 1) if topology_reliable else None,
                    "parent_edge_assignment_unique": (top_edges == 1) if topology_reliable else None,
                    "exact_topology_unique": (True if topology_reliable and top_edges == 1 else None),
                    "exact_topology_uniqueness_status": (
                        "PROVEN_UNIQUE_BY_SINGLE_WINNING_PARENT_EDGE_ASSIGNMENT"
                        if topology_reliable and top_edges == 1
                        else "NOT_EVALUATED_CANONICAL_SHAPE_ISOMORPHISM_FOR_MULTIPLE_EDGE_ASSIGNMENTS"
                        if topology_reliable
                        else "ABSTAIN_PRIMARY_LIKELIHOOD_NONCONVERGENCE"
                    ),
                    "topology_class_by_top_vertex_set": json.dumps(
                        topology_by_vertex, sort_keys=True, separators=(",", ":")
                    ),
                    "topology_derivation_status": topology_status,
                    "fixed_error_grid_top_vertex_set_ids": json.dumps(
                        grid_top_ids, sort_keys=True, separators=(",", ":")
                    ),
                    "fixed_error_grid_stable_with_quality_primary": grid_stable,
                    "fixed_error_grid_all_converged": grid_converged,
                    "conditional_candidate_ranking_bootstrap_status": bootstrap["status"],
                    "conditional_candidate_ranking_bootstrap_replicates": bootstrap["replicates"],
                    "conditional_candidate_ranking_bootstrap_seed": bootstrap["seed"],
                    "conditional_candidate_ranking_bootstrap_top_vertex_set_frequency": json.dumps(
                        bootstrap["top_vertex_set_frequency"], sort_keys=True, separators=(",", ":")
                    ),
                    "conditional_candidate_ranking_bootstrap_primary_top_set_frequency": bootstrap["primary_top_set_frequency"],
                    "conditional_candidate_ranking_bootstrap_all_converged": bootstrap["all_converged"],
                }
            )
    # Empty/abstain templates must not erase the auditable effective-k contract.
    row["k_observed_alt_active"] = effective_k
    row["k_scoring_alt_observed"] = popcount(scoring_alt_universe)
    row["n_not_structural_alt_active_sites"] = component.k - effective_k
    row["structural_partial_pattern_groups"] = len(structural_partial)
    row["unit_semantic_sha256"] = semantic_digest(
        {key: value for key, value in row.items() if not key.startswith("_")}
    )
    unit_total_elapsed_ns = _monotonic_ns() - unit_started_ns
    if candidate_generation_elapsed_ns + likelihood_fit_elapsed_ns > unit_total_elapsed_ns:
        raise AssertionError("disjoint runtime segments exceed unit total monotonic wall time")
    row["_runtime_diagnostics"] = {
        "candidate_generation_elapsed_seconds": candidate_generation_elapsed_ns / 1_000_000_000.0,
        "likelihood_fit_elapsed_seconds": likelihood_fit_elapsed_ns / 1_000_000_000.0,
        "unit_total_elapsed_seconds": unit_total_elapsed_ns / 1_000_000_000.0,
    }
    row["_runtime_segment_invoked"] = {
        "candidate_generation_invoked": candidate_generation_invoked,
        "likelihood_fit_invoked": likelihood_fit_invoked,
    }
    row["_structural_cache_diagnostics"] = structural_cache_diagnostics
    return row


def _load_site_catalog(path: Path) -> tuple[dict[int, int], str, str, str]:
    sites: dict[int, int] = {}
    semantic = hashlib.sha256()
    dataset = chrom = ""
    with _open_tsv(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            index = int(row["site_index"])
            if index in sites:
                raise RuntimeError(f"duplicate site index: {index}")
            sites[index] = int(row["pos1"])
            chrom = chrom or str(row["chrom"])
            if row["chrom"] != chrom:
                raise RuntimeError("site catalog chromosome drift")
            semantic.update(f"{index}\t{chrom}\t{row['pos1']}\t{row['ref']}\t{row['alt']}\n".encode())
    if not sites:
        raise RuntimeError("empty site catalog")
    # Dataset is not stored in the site catalog and is resolved from components.
    return sites, dataset, chrom, semantic.hexdigest()


def _load_components(
    components_path: Path,
    membership_path: Path,
    selected_thresholds: set[int] | None,
    selected_bases: set[str] | None,
) -> tuple[
    dict[tuple[str, str | None, int, str], Component],
    dict[tuple[str, str | None, int], dict[int, str]],
    dict[str, str],
]:
    metadata: dict[tuple[str, str | None, int, str], dict[str, Any]] = {}
    component_semantic = hashlib.sha256()
    with _open_tsv(components_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            threshold = int(row["threshold"])
            if selected_thresholds is not None and threshold not in selected_thresholds:
                continue
            basis = row.get("linkage_basis") or row.get("component_basis") or row.get("basis") or "pooled"
            if selected_bases is not None and basis not in selected_bases:
                continue
            phase_set = row.get("phase_set") or None
            key = (basis, phase_set, threshold, row["component_id"])
            if key in metadata:
                raise RuntimeError(f"duplicate component: {key}")
            metadata[key] = {
                "dataset": row["dataset"],
                "chrom": row["chrom"],
                "threshold": threshold,
                "basis": basis,
                "component_id": row["component_id"],
                "start1": int(row["start1"]),
                "end1": int(row["end1"]),
                "k": int(row["k"]),
                "phase_set": phase_set,
                "phase_set_status": row.get("phase_set_status") or "UNKNOWN",
                "inference_role": row.get("inference_role") or "UNKNOWN",
            }
            component_semantic.update(
                json.dumps(metadata[key], sort_keys=True, separators=(",", ":")).encode() + b"\n"
            )
    members: dict[tuple[str, str | None, int, str], list[tuple[int, int]]] = defaultdict(list)
    site_to_component: dict[tuple[str, str | None, int], dict[int, str]] = defaultdict(dict)
    membership_semantic = hashlib.sha256()
    with _open_tsv(membership_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            threshold = int(row["threshold"])
            if selected_thresholds is not None and threshold not in selected_thresholds:
                continue
            basis = row.get("linkage_basis") or row.get("component_basis") or row.get("basis") or "pooled"
            if selected_bases is not None and basis not in selected_bases:
                continue
            phase_set = row.get("phase_set") or None
            key = (basis, phase_set, threshold, row["component_id"])
            if key not in metadata:
                raise RuntimeError(f"membership references unknown component: {key}")
            site_index = int(row["site_index"])
            pos1 = int(row["pos1"])
            basis_threshold = (basis, phase_set, threshold)
            if site_index in site_to_component[basis_threshold]:
                raise RuntimeError(
                    f"site belongs to multiple components at basis/threshold {basis_threshold}: {site_index}"
                )
            site_to_component[basis_threshold][site_index] = row["component_id"]
            members[key].append((site_index, pos1))
            membership_semantic.update(
                f"{basis}\t{phase_set or ''}\t{threshold}\t{row['component_id']}\t{site_index}\t{pos1}\n".encode()
            )
    components = {}
    for key, meta in metadata.items():
        ordered = sorted(members.get(key, []))
        if len(ordered) != meta["k"]:
            raise RuntimeError(f"component k/membership mismatch: {key} {meta['k']} != {len(ordered)}")
        if ordered and (ordered[0][1] != meta["start1"] or ordered[-1][1] != meta["end1"]):
            raise RuntimeError(f"component coordinate/membership mismatch: {key}")
        components[key] = Component(**meta, site_indices=tuple(index for index, _ in ordered))
    return components, dict(site_to_component), {
        "components_semantic_sha256": component_semantic.hexdigest(),
        "membership_semantic_sha256": membership_semantic.hexdigest(),
    }


def _available_component_bases(path: Path) -> set[str]:
    with _open_tsv(path) as handle:
        return {
            row.get("linkage_basis") or row.get("component_basis") or row.get("basis") or "pooled"
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _family_targets(raw_family: str, selected: set[str]) -> tuple[str, ...]:
    targets = []
    if raw_family in selected:
        targets.append(raw_family)
    if "pooled" in selected:
        targets.append("pooled")
    return tuple(targets)


def _basis_matches_family(basis: str, family: str) -> bool:
    normalized = basis.lower().replace("_", "")
    if normalized in {"hp1", "1", "pshp1", "missingpshp1"}:
        return family == "1"
    if normalized in {"hp2", "2", "pshp2", "missingpshp2"}:
        return family == "2"
    return True


def _component_accepts_molecule(
    component: Component, raw_family: str, phase_set: str | None
) -> bool:
    if not _basis_matches_family(component.basis, raw_family):
        return False
    normalized = component.basis.lower().replace("_", "")
    if normalized in {"pshp1", "pshp2"}:
        return phase_set is not None and component.phase_set == phase_set
    if normalized in {"missingpshp1", "missingpshp2"}:
        return phase_set is None and component.phase_set is None
    return True


def _build_evidence_route_index(
    site_to_component: Mapping[tuple[str, str | None, int], Mapping[int, str]],
    families: Sequence[str],
) -> tuple[
    dict[tuple[str, str | None], tuple[tuple[str, str | None, int, Mapping[int, str], str], ...]],
    dict[str, tuple[tuple[str, str | None, int, Mapping[int, str], str], ...]],
]:
    """Index component mappings by raw HP family and exact PS before reading molecules.

    PS-aware bases enter only the ``(family, exact_PS)`` index; missing-PS bases
    enter only ``(family, None)``; diagnostic bases use the family wildcard.
    Consequently a molecule never scans unrelated chromosome-wide PS mappings.
    """
    selected = set(families)
    exact: dict[
        tuple[str, str | None],
        list[tuple[str, str | None, int, Mapping[int, str], str]],
    ] = defaultdict(list)
    wildcard: dict[
        str,
        list[tuple[str, str | None, int, Mapping[int, str], str]],
    ] = defaultdict(list)
    raw_families = ("1", "2", "3", "4", "none", "pooled")
    for (basis, component_phase_set, threshold), mapping in site_to_component.items():
        normalized = basis.lower().replace("_", "")
        if normalized in {"pshp1", "pshp2"} and component_phase_set is None:
            raise RuntimeError(f"PS-aware component mapping lacks exact PS: {basis} B{threshold}")
        if normalized in {"missingpshp1", "missingpshp2"} and component_phase_set is not None:
            raise RuntimeError(f"missing-PS component mapping unexpectedly has PS: {basis} B{threshold}")
        for raw_family in raw_families:
            if not _basis_matches_family(basis, raw_family):
                continue
            for target_family in _family_targets(raw_family, selected):
                if not _basis_matches_family(basis, target_family):
                    continue
                route = (basis, component_phase_set, threshold, mapping, target_family)
                if normalized in {"pshp1", "pshp2", "missingpshp1", "missingpshp2"}:
                    exact[(raw_family, component_phase_set)].append(route)
                else:
                    wildcard[raw_family].append(route)
    return (
        {key: tuple(values) for key, values in exact.items()},
        {key: tuple(values) for key, values in wildcard.items()},
    )


def build_evidence(
    calls_path: Path,
    components: Mapping[tuple[str, str | None, int, str], Component],
    site_to_component: Mapping[tuple[str, str | None, int], Mapping[int, str]],
    families: Sequence[str],
    known_site_indices: set[int],
) -> tuple[dict[tuple[str, str | None, int, str, str], UnitEvidence], dict[str, Any]]:
    selected = set(families)
    units = {
        (basis, phase_set, threshold, component_id, family): UnitEvidence(component, family)
        for (basis, phase_set, threshold, component_id), component in components.items()
        for family in families
        if _basis_matches_family(basis, family)
    }
    counts = Counter()
    hp_rows = Counter()
    raw_code_counts = Counter()
    semantic = hashlib.sha256()
    seen_molecules: set[str] = set()
    exact_routes, wildcard_routes = _build_evidence_route_index(site_to_component, families)
    counts["evidence_component_mapping_keys_total"] = len(site_to_component)
    counts["evidence_route_index_exact_key_count"] = len(exact_routes)
    counts["evidence_route_index_exact_route_count"] = sum(map(len, exact_routes.values()))
    counts["evidence_route_index_wildcard_route_count"] = sum(map(len, wildcard_routes.values()))
    with _open_tsv(calls_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["sparse_molecule_rows_total"] += 1
            molecule = row["molecule_id"]
            if molecule in seen_molecules:
                raise RuntimeError(f"duplicate molecule row in sparse calls: {molecule}")
            seen_molecules.add(molecule)
            raw_family = row["hp_family"]
            row_phase_set = row.get("phase_set") or None
            hp_rows[raw_family] += 1
            counts["sparse_molecule_rows_known_ps" if row_phase_set is not None else "sparse_molecule_rows_missing_ps"] += 1
            targets = _family_targets(raw_family, selected)
            if not targets:
                counts["sparse_molecule_rows_excluded_hp"] += 1
                continue
            counts["sparse_molecule_rows_selected_hp"] += 1
            indices = tuple(int(value) for value in row["site_indices"].split(",") if value != "")
            codes = tuple(row["call_codes"])
            qualities_raw = tuple(row["base_qualities"].split(","))
            qualities = tuple(None if value == "" else int(value) for value in qualities_raw)
            if len(indices) != len(codes) or len(indices) != len(qualities):
                raise RuntimeError(f"sparse vector length mismatch for molecule {molecule}")
            if len(indices) != len(set(indices)) or tuple(sorted(indices)) != indices:
                raise RuntimeError(f"sparse site indices are not sorted unique for molecule {molecule}")
            if set(indices) - known_site_indices:
                raise RuntimeError(f"sparse calls reference unknown site index for molecule {molecule}")
            if set(codes) - RAW_CODES:
                raise RuntimeError(f"invalid call code for molecule {molecule}")
            raw_code_counts.update(codes)
            counts["selected_fixed_ra_calls_with_bq"] += sum(
                code in FIXED_CODES and quality is not None
                for code, quality in zip(codes, qualities)
            )
            counts["selected_fixed_ra_calls_without_bq"] += sum(
                code in FIXED_CODES and quality is None
                for code, quality in zip(codes, qualities)
            )
            codes_by_site = dict(zip(indices, codes))
            qualities_by_site = dict(zip(indices, qualities))
            semantic.update(
                f"{molecule}\t{raw_family}\t{row_phase_set or ''}\t{','.join(map(str, indices))}\t{''.join(codes)}\n".encode()
            )
            projected_this_molecule = 0
            routes = (
                wildcard_routes.get(raw_family, ())
                + exact_routes.get((raw_family, row_phase_set), ())
            )
            counts["evidence_route_index_lookups"] += 2
            counts["evidence_mapping_routes_visited"] += len(routes)
            counts["evidence_mapping_keys_naive_scan_reference"] += len(site_to_component)
            for basis, component_phase_set, threshold, mapping, family in routes:
                touched = sorted({mapping[index] for index in indices if index in mapping})
                for component_id in touched:
                    key = (basis, component_phase_set, threshold, component_id, family)
                    if key not in units:
                        raise RuntimeError(f"indexed evidence route has no unit: {key}")
                    if not _component_accepts_molecule(
                        units[key].component, raw_family, row_phase_set
                    ):
                        raise AssertionError(f"evidence route index violated HP/PS isolation: {key}")
                    units[key].add_projection(codes_by_site, qualities_by_site)
                    counts["molecule_component_family_projections"] += 1
                    projected_this_molecule += 1
            if projected_this_molecule:
                counts["sparse_molecule_rows_included_in_at_least_one_selected_unit"] += 1
            else:
                counts["sparse_molecule_rows_excluded_by_component_or_phase_set_contract"] += 1
    counts["unique_sparse_molecules"] = len(seen_molecules)
    counts["hp_family_rows"] = dict(sorted(hp_rows.items()))
    counts["selected_sparse_call_code_counts"] = dict(sorted(raw_code_counts.items()))
    counts["calls_semantic_sha256"] = semantic.hexdigest()
    return units, dict(counts)


def _verify_extractor_receipt(input_dir: Path, files: Sequence[Path]) -> dict[str, Any]:
    receipt_path = input_dir / "receipt.json"
    if not receipt_path.is_file():
        raise FileNotFoundError(f"missing upstream receipt: {receipt_path}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if (
        receipt.get("schema_name") != "intersubmod.lossless_read_linkage_chromosome_receipt"
        or receipt.get("schema_version") != EXTRACTOR_SCHEMA_VERSION
        or receipt.get("all_pass") is not True
    ):
        raise RuntimeError(
            f"upstream extractor receipt is not a passing v{EXTRACTOR_SCHEMA_VERSION} receipt"
        )
    integrity = receipt.get("receipt_integrity") or {}
    sidecar = input_dir / str(integrity.get("sidecar_name") or "")
    if integrity.get("scheme") != "external_sha256_sidecar_v1" or not sidecar.is_file():
        raise RuntimeError("upstream extractor receipt lacks its external checksum")
    checksum_fields = sidecar.read_text(encoding="ascii", errors="strict").strip().split()
    if (
        len(checksum_fields) != 2
        or checksum_fields[0] != sha256_path(receipt_path)
        or checksum_fields[1] != receipt_path.name
    ):
        raise RuntimeError("upstream extractor receipt checksum mismatch")
    declared = receipt.get("outputs") or {}
    verified = {}
    for path in files:
        identity = declared.get(path.name)
        if not identity:
            raise RuntimeError(f"upstream receipt lacks output identity: {path.name}")
        actual = sha256_path(path)
        if identity.get("sha256") != actual or int(identity.get("size_bytes", -1)) != path.stat().st_size:
            raise RuntimeError(f"upstream output identity mismatch: {path}")
        verified[path.name] = {"size_bytes": path.stat().st_size, "sha256": actual}
    return {
        "path": str(receipt_path),
        "sha256": sha256_path(receipt_path),
        "verified_outputs": verified,
        "scope": receipt.get("scope") or {},
        "parameters": receipt.get("parameters") or {},
        "phase_set_contract_counts": receipt.get("phase_set_contract_counts") or {},
    }


def summarize_units(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    rows = list(rows)
    status = Counter(str(row["selection_status"]) for row in rows)
    generation = Counter(str(row["candidate_generation_status"]) for row in rows)
    reliable_rows = [
        row for row in rows
        if str(row["topology_derivation_status"]).startswith("ANALYTICAL_COMPLETE")
    ]
    topology_classes = Counter()
    unique_topology_classes = Counter()
    ambiguous_topology_class_sets = Counter()
    for row in reliable_rows:
        classes = tuple(json.loads(row["topology_classes"]))
        topology_classes.update(classes)
        if row["coarse_topology_class_unique"] is True:
            if len(classes) != 1:
                raise AssertionError("unique coarse topology row does not have exactly one class")
            unique_topology_classes[classes[0]] += 1
        elif row["coarse_topology_class_unique"] is False:
            ambiguous_topology_class_sets["|".join(classes)] += 1
    unique_statuses = {
        "T1_CANDIDATE_UNIQUE",
        "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED",
        "LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED",
        "LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE",
    }
    quality_primary_unique = sum(row["selection_status"] in unique_statuses for row in rows)
    quality_primary_tied = sum(
        row["selection_status"] == "LIKELIHOOD_TIED_VERTEX_SETS" for row in rows
    )
    rank_abstain = sum(
        "ABSTAIN" in str(row["selection_status"])
        or str(row["selection_status"]).startswith("NO_INFORMATIVE")
        for row in rows
    )
    if quality_primary_unique + quality_primary_tied + rank_abstain != len(rows):
        classified = unique_statuses | {"LIKELIHOOD_TIED_VERTEX_SETS"}
        unknown = {
            str(row["selection_status"])
            for row in rows
            if str(row["selection_status"]) not in classified
            and "ABSTAIN" not in str(row["selection_status"])
            and not str(row["selection_status"]).startswith("NO_INFORMATIVE")
        }
        raise RuntimeError(f"unclassified quality-primary selection status: {sorted(unknown)}")
    if sum(unique_topology_classes.values()) != sum(
        row["coarse_topology_class_unique"] is True for row in reliable_rows
    ):
        raise AssertionError("mutually-exclusive topology class counts do not conserve unique units")
    if sum(unique_topology_classes.values()) + sum(ambiguous_topology_class_sets.values()) != len(reliable_rows):
        raise AssertionError("unique plus ambiguous topology units do not conserve evaluated units")
    if any(
        row["raw_tree_candidates_T"] is not None
        and row["distinct_vertex_sets_V"] is not None
        and int(row["raw_tree_candidates_T"]) < int(row["distinct_vertex_sets_V"])
        for row in rows
    ):
        raise AssertionError("raw candidate tree count T is smaller than vertex-set count V")
    return {
        "n_component_hp_units": len(rows),
        "n_components": len({
            (row["component_basis"], row["threshold"], row["component_id"]) for row in rows
        }),
        "molecule_component_projections": sum(int(row["molecule_component_projections"]) for row in rows),
        "informative_scoring_molecules": sum(int(row["informative_scoring_molecules"]) for row in rows),
        "all_x_excluded_molecules": sum(int(row["all_x_excluded_molecules"]) for row in rows),
        "structural_retained_molecules": sum(int(row["structural_retained_molecules"]) for row in rows),
        "below_minread_scoring_molecules": sum(int(row["below_minread_scoring_molecules"]) for row in rows),
        "bq_scoring_molecules": sum(int(row["bq_scoring_molecules"]) for row in rows),
        "non_ra_cells_marginalized": sum(int(row["non_ra_cells_marginalized"]) for row in rows),
        "projected_call_class_counts": {
            "R_or_A": sum(int(row["fixed_ra_cells"]) for row in rows),
            "A": sum(int(row["alt_cells"]) for row in rows),
            "O": sum(int(row["other_cells"]) for row in rows),
            "D": sum(int(row["deletion_cells"]) for row in rows),
            "S": sum(int(row["refskip_cells"]) for row in rows),
            "L": sum(int(row["low_baseq_cells"]) for row in rows),
            "explicit_X": sum(int(row["explicit_x_cells"]) for row in rows),
            "implicit_uncovered_X": sum(int(row["implicit_uncovered_x_cells"]) for row in rows),
        },
        "selection_status_counts": dict(sorted(status.items())),
        "selection_status_fraction_of_units": {
            key: value / len(rows) if rows else None for key, value in sorted(status.items())
        },
        "candidate_generation_status_counts": dict(sorted(generation.items())),
        "raw_tree_candidates_T_complete_units": sum(
            int(row["raw_tree_candidates_T"] or 0) for row in rows
        ),
        "distinct_vertex_sets_V_complete_units": sum(
            int(row["distinct_vertex_sets_V"] or 0) for row in rows
        ),
        "k_route_counts": {
            "EXACT_COMPLETE": sum(bool(row["candidate_generation_complete"]) for row in rows),
            "EXACT_INCOMPLETE": sum(
                str(row["candidate_generation_status"]).startswith("EXACT_ENUMERATION_INCOMPLETE")
                for row in rows
            ),
            "NOT_RUN_NO_STRUCTURE": sum(
                row["candidate_generation_status"] in {
                    "NO_INFORMATIVE_PATTERN", "NO_STRUCTURAL_ALT_PATTERN_AT_DECLARED_MINREAD"
                }
                for row in rows
            ),
            "GT_EXACT_LIMIT": sum(row["candidate_generation_status"] == "K_GT_EXACT_LIMIT" for row in rows),
        },
        "solver_complete_units": sum(bool(row["candidate_generation_complete"]) for row in rows),
        "solver_incomplete_or_not_run_units": sum(not bool(row["candidate_generation_complete"]) for row in rows),
        "quality_primary_unique_vertex_units": quality_primary_unique,
        "quality_primary_tied_vertex_units": quality_primary_tied,
        "rank_abstain_units": rank_abstain,
        "fixed_error_grid_stable_units": sum(
            row["fixed_error_grid_stable_with_quality_primary"] is True for row in rows
        ),
        "fixed_error_grid_evaluated_units": sum(
            row["fixed_error_grid_stable_with_quality_primary"] is not None for row in rows
        ),
        "conditional_candidate_ranking_bootstrap_status_counts": dict(sorted(Counter(
            str(row["conditional_candidate_ranking_bootstrap_status"]) for row in rows
        ).items())),
        "conditional_candidate_ranking_bootstrap_complete_units": sum(
            row["conditional_candidate_ranking_bootstrap_status"] == "COMPLETE" for row in rows
        ),
        "conditional_candidate_ranking_bootstrap_not_run_units": sum(
            str(row["conditional_candidate_ranking_bootstrap_status"]).startswith("NOT_") for row in rows
        ),
        "conditional_candidate_ranking_bootstrap_nonconverged_units": sum(
            row["conditional_candidate_ranking_bootstrap_status"] == "ABSTAIN_NONCONVERGENCE" for row in rows
        ),
        "topology_class_inclusion_counts": dict(sorted(topology_classes.items())),
        "topology_class_inclusion_counts_denominator": len(reliable_rows),
        "coarse_topology_unique_class_counts": dict(sorted(unique_topology_classes.items())),
        "coarse_topology_ambiguous_class_set_counts": dict(sorted(ambiguous_topology_class_sets.items())),
        "topology_evaluated_units": len(reliable_rows),
        "topology_derivation_status_counts": dict(sorted(Counter(
            str(row["topology_derivation_status"]) for row in rows
        ).items())),
        "exact_topology_uniqueness_status_counts": dict(sorted(Counter(
            str(row["exact_topology_uniqueness_status"]) for row in rows
        ).items())),
        "coarse_topology_class_unique_units": sum(
            row["coarse_topology_class_unique"] is True for row in reliable_rows
        ),
        "coarse_topology_multiple_class_units": sum(
            row["coarse_topology_class_unique"] is False for row in reliable_rows
        ),
        "parent_edge_assignment_unique_units": sum(
            row["parent_edge_assignment_unique"] is True for row in rows
        ),
        "exact_topology_proven_unique_units": sum(
            row["exact_topology_unique"] is True for row in rows
        ),
        "k_component_sites_total": sum(int(row["k_component_sites"]) for row in rows),
        "k_observed_alt_active_total": sum(int(row["k_observed_alt_active"]) for row in rows),
        "k_scoring_alt_observed_total": sum(int(row["k_scoring_alt_observed"]) for row in rows),
        "not_structural_alt_active_sites_total": sum(
            int(row["n_not_structural_alt_active_sites"]) for row in rows
        ),
        "structural_partial_pattern_groups": sum(
            int(row["structural_partial_pattern_groups"]) for row in rows
        ),
        "partial_group_coverage_denominator": sum(
            int(row["partial_group_coverage_denominator"]) for row in rows
        ),
        "partial_groups_covered": sum(int(row["partial_groups_covered"]) for row in rows),
        "partial_groups_unsatisfied": sum(int(row["partial_groups_unsatisfied"]) for row in rows),
    }


def summarize_units_by_basis_threshold(rows: Iterable[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["component_basis"]), int(row["threshold"]))].append(row)
    nested: dict[str, dict[str, Any]] = defaultdict(dict)
    for (basis, threshold), subset in sorted(grouped.items()):
        nested[basis][str(threshold)] = summarize_units(subset)
    return dict(nested)


def _structural_alt_universe(unit: UnitEvidence, minread: int) -> int:
    mask = 0
    for pattern, count in unit.pattern_counts.items():
        if count < minread or set(pattern) == {"X"}:
            continue
        for bit, symbol in enumerate(pattern):
            if symbol == "A":
                mask |= 1 << bit
    return mask


def _funnel_denominator(
    entries: Iterable[tuple[str, int, int]],
) -> dict[str, Any]:
    """Summarize (pattern, weight, structural_alt_universe) entries."""
    total = full = partial = conceptual_sum = effective_sum = effective_zero = 0
    u_distribution: Counter[int] = Counter()
    conceptual_distribution: Counter[int] = Counter()
    effective_distribution: Counter[int] = Counter()
    for pattern, weight, universe in entries:
        symbolic = SymbolicPattern.from_string(pattern)
        if symbolic.fixed_mask == 0:
            continue
        total += weight
        is_partial = symbolic.n_free > 0
        partial += weight * is_partial
        full += weight * (not is_partial)
        conceptual = 1 << symbolic.n_free
        effective = (
            0
            if symbolic.alt_mask & ~universe
            else 1 << popcount(symbolic.free_mask & universe)
        )
        u_distribution[symbolic.n_free] += weight
        conceptual_distribution[conceptual] += weight
        effective_distribution[effective] += weight
        conceptual_sum += weight * conceptual
        effective_sum += weight * effective
        effective_zero += weight * (effective == 0)
    if full + partial != total:
        raise AssertionError("full/partial funnel does not conserve denominator")
    return {
        "denominator": total,
        "full": full,
        "partial": partial,
        "u_number_of_X_distribution": {str(k): v for k, v in sorted(u_distribution.items())},
        "conceptual_completions_2_pow_u_distribution": {
            str(k): v for k, v in sorted(conceptual_distribution.items())
        },
        "conceptual_completions_weighted_total": conceptual_sum,
        "observed_alt_effective_completions_distribution": {
            str(k): v for k, v in sorted(effective_distribution.items())
        },
        "observed_alt_effective_completions_weighted_total": effective_sum,
        "observed_alt_effective_zero_due_to_fixed_alt_outside_structural_universe": effective_zero,
    }


def partial_pattern_funnel(
    units: Sequence[UnitEvidence], ranked_rows: Sequence[Mapping[str, Any]], minread: int
) -> dict[str, Any]:
    rax_entries: list[tuple[str, int, int]] = []
    quality_entries: list[tuple[str, int, int]] = []
    molecule_entries: list[tuple[str, int, int]] = []
    structural_entries: list[tuple[str, int, int]] = []
    units_with_partial = 0
    for unit in units:
        universe = _structural_alt_universe(unit, minread)
        partial_in_unit = False
        for pattern, count in unit.pattern_counts.items():
            if set(pattern) == {"X"}:
                continue
            rax_entries.append((pattern, 1, universe))
            molecule_entries.append((pattern, count, universe))
            if count >= minread:
                structural_entries.append((pattern, 1, universe))
                partial_in_unit = partial_in_unit or "X" in pattern
        for pattern, _qualities in unit.quality_pattern_counts:
            if set(pattern) != {"X"}:
                quality_entries.append((pattern, 1, universe))
        units_with_partial += partial_in_unit
    coverage_denominator = sum(int(row["partial_group_coverage_denominator"]) for row in ranked_rows)
    covered = sum(int(row["partial_groups_covered"]) for row in ranked_rows)
    unsatisfied = sum(int(row["partial_groups_unsatisfied"]) for row in ranked_rows)
    if covered + unsatisfied != coverage_denominator:
        raise AssertionError("partial group coverage funnel does not conserve")
    return {
        "definitions": {
            "conceptual": (
                "a pattern with u X has 2^u full-cube conceptual completions; completion-wise tree worlds "
                "are not materialized, while active compatible vertex indices are materialized per MILP "
                "construction for each reduced sparse group-hit row"
            ),
            "effective": (
                "primary observed-alt universe: 2^popcount(free_mask & structural_alt_universe); "
                "zero when a scoring pattern fixes ALT outside that minread-specific universe"
            ),
            "universe_source": f"exact R/A/X structural pattern count >= {minread}",
        },
        "unique_rax_pattern_groups": _funnel_denominator(rax_entries),
        "bq_quality_pattern_groups": _funnel_denominator(quality_entries),
        "molecule_projections": _funnel_denominator(molecule_entries),
        "structural_unique_rax_pattern_groups": _funnel_denominator(structural_entries),
        "units_denominator": len(units),
        "units_with_partial_structural_groups": units_with_partial,
        "partial_group_coverage_denominator": coverage_denominator,
        "partial_groups_covered": covered,
        "partial_groups_unsatisfied": unsatisfied,
    }


def validate_output_directory_contract(
    output_dir: Path, *, require_existing_empty: bool = False
) -> None:
    """Fail closed for fresh outputs and inode-bound direct-pilot outputs."""
    if require_existing_empty:
        if output_dir.is_symlink() or not output_dir.exists() or not output_dir.is_dir():
            raise FileExistsError(
                f"required preflight output directory is unavailable or not a real directory: {output_dir}"
            )
        if next(output_dir.iterdir(), None) is not None:
            raise FileExistsError(
                f"required preflight output directory is not empty: {output_dir}"
            )
        return
    if output_dir.is_symlink() or output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")


def run(
    input_dir: Path,
    output_dir: Path,
    *,
    thresholds: set[int] | None = None,
    component_bases: set[str] | None = None,
    families: Sequence[str] = PRIMARY_HP_FAMILIES,
    minread: int = 3,
    structural_exact_pattern_minreads: Sequence[int] | None = None,
    primary_structural_exact_pattern_minread: int = 3,
    exact_k_max: int = 12,
    max_vertex_sets: int = 256,
    solver_time_limit_seconds: float = 30.0,
    fixed_error_grid: Sequence[float] = (0.005, 0.01, 0.02, 0.05),
    minimum_bq_error_rate: float = 1e-6,
    maximum_bq_error_rate: float = 0.25,
    conditional_candidate_ranking_bootstrap_replicates: int = 0,
    conditional_candidate_ranking_bootstrap_seed: int = 20260716,
    tie_tolerance: float = 1e-6,
    enable_structural_enumeration_cache: bool = True,
    require_existing_empty_output_dir: bool = False,
) -> dict[str, Any]:
    validate_output_directory_contract(
        output_dir,
        require_existing_empty=require_existing_empty_output_dir,
    )
    if not families or len(set(families)) != len(tuple(families)):
        raise ValueError("families must be nonempty and unique")
    if set(families) - {"1", "2", "3", "4", "none", "pooled"}:
        raise ValueError("invalid requested HP family")
    if (
        not fixed_error_grid
        or any(not 0.0 < value < 0.5 for value in fixed_error_grid)
        or not 0 < minimum_bq_error_rate <= maximum_bq_error_rate < 0.5
        or conditional_candidate_ranking_bootstrap_replicates < 0
        or tie_tolerance < 0
    ):
        raise ValueError("invalid likelihood parameters")
    calls_path = _one_file(input_dir, ".molecule_sparse_calls.tsv.gz")
    sites_path = _one_file(input_dir, ".site_catalog.tsv.gz")
    components_path = _one_file(input_dir, ".components.tsv.gz")
    membership_path = _one_file(input_dir, ".site_component_membership.tsv.gz")
    upstream = _verify_extractor_receipt(input_dir, (calls_path, sites_path, components_path, membership_path))
    sites, _, site_chrom, site_semantic = _load_site_catalog(sites_path)
    available_bases = set(upstream.get("parameters", {}).get("component_linkage_bases") or ())
    available_bases.update(_available_component_bases(components_path))
    if component_bases is None:
        selected_bases = set(PRIMARY_COMPONENT_BASES)
        basis_mode = "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY"
    else:
        selected_bases = set(component_bases)
        basis_mode = (
            "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY"
            if selected_bases == set(PRIMARY_COMPONENT_BASES)
            else "EXPLICIT_USER_SELECTION_NONPRIMARY_OR_SENSITIVITY"
        )
    if selected_bases - available_bases:
        raise ValueError(f"requested component bases unavailable: {sorted(selected_bases - available_bases)}")
    components, site_to_component, input_semantics = _load_components(
        components_path, membership_path, thresholds, selected_bases
    )
    if components and {component.chrom for component in components.values()} != {site_chrom}:
        raise RuntimeError("component/site catalog chromosome mismatch")
    units, input_counts = build_evidence(
        calls_path,
        components,
        site_to_component,
        tuple(families),
        set(sites),
    )
    minread_grid = tuple(sorted(set(structural_exact_pattern_minreads or (minread,))))
    if not minread_grid or min(minread_grid) < 1 or primary_structural_exact_pattern_minread not in minread_grid:
        raise ValueError("structural exact-pattern minread grid must contain the primary minread")
    ordered_unit_keys = sorted(units, key=lambda key: tuple("" if value is None else str(value) for value in key))
    unit_rows: list[dict[str, Any]] = []
    unit_rows_by_minread: dict[int, list[dict[str, Any]]] = {}
    structural_enumeration_cache = (
        StructuralEnumerationCache()
        if enable_structural_enumeration_cache
        else None
    )
    if require_existing_empty_output_dir:
        validate_output_directory_contract(output_dir, require_existing_empty=True)
        staging_dir = output_dir
    else:
        output_dir.parent.mkdir(parents=True, exist_ok=True)
        staging_dir = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging.", dir=output_dir.parent))
    pattern_path = output_dir / "m2_symbolic_pattern_counts.tsv.gz"
    unit_path = output_dir / "m2_component_hp_rank_units.tsv.gz"
    candidate_path = output_dir / "m2_compressed_vertex_set_candidates.tsv.gz"
    responsibility_path = output_dir / "m2_winning_pattern_state_responsibilities.tsv.gz"
    runtime_path = output_dir / "m2_unit_runtime_diagnostics.tsv.gz"
    staged_pattern_path = staging_dir / pattern_path.name
    staged_unit_path = staging_dir / unit_path.name
    staged_candidate_path = staging_dir / candidate_path.name
    staged_responsibility_path = staging_dir / responsibility_path.name
    staged_runtime_path = staging_dir / runtime_path.name
    pattern_semantic = _SemanticListHasher()
    unit_semantic = _SemanticListHasher()
    candidate_semantic = _SemanticListHasher()
    responsibility_semantic = _SemanticListHasher()
    all_pattern_hashes_valid = True
    all_candidate_hashes_valid = True
    all_responsibility_hashes_valid = True
    candidate_row_count = 0
    winning_candidate_row_count = 0
    responsibility_row_count = 0
    candidate_optimizer_status_counts: Counter[str] = Counter()
    candidate_slsqp_status_counts: Counter[str] = Counter()
    candidate_globally_certified_count = 0
    candidate_nonidentifiable_weight_count = 0
    candidate_max_global_ll_gap_bound = 0.0
    candidate_max_simplex_residual = 0.0
    all_converged_candidate_certificates_valid = True
    all_optimizer_abstain_candidates_not_winners = True
    with (
        gzip.open(staged_pattern_path, "wt", encoding="utf-8", newline="") as pattern_handle,
        gzip.open(staged_candidate_path, "wt", encoding="utf-8", newline="") as candidate_handle,
        gzip.open(staged_responsibility_path, "wt", encoding="utf-8", newline="") as responsibility_handle,
        gzip.open(staged_runtime_path, "wt", encoding="utf-8", newline="") as runtime_handle,
    ):
        pattern_writer = csv.DictWriter(
            pattern_handle, PATTERN_FIELDS, delimiter="\t", extrasaction="raise"
        )
        candidate_writer = csv.DictWriter(
            candidate_handle, CANDIDATE_FIELDS, delimiter="\t", extrasaction="raise"
        )
        responsibility_writer = csv.DictWriter(
            responsibility_handle, RESPONSIBILITY_FIELDS, delimiter="\t", extrasaction="raise"
        )
        runtime_writer = csv.DictWriter(
            runtime_handle, RUNTIME_DIAGNOSTIC_FIELDS, delimiter="\t", extrasaction="raise"
        )
        pattern_writer.writeheader()
        candidate_writer.writeheader()
        responsibility_writer.writeheader()
        runtime_writer.writeheader()

        def emit_candidate(record: dict[str, Any]) -> None:
            nonlocal candidate_row_count, winning_candidate_row_count
            nonlocal all_candidate_hashes_valid, candidate_globally_certified_count
            nonlocal candidate_nonidentifiable_weight_count
            nonlocal candidate_max_global_ll_gap_bound, candidate_max_simplex_residual
            nonlocal all_converged_candidate_certificates_valid
            nonlocal all_optimizer_abstain_candidates_not_winners
            candidate_writer.writerow({key: _cell(record.get(key)) for key in CANDIDATE_FIELDS})
            candidate_semantic.update(record)
            candidate_row_count += 1
            winning_candidate_row_count += int(bool(record["is_winner"]))
            candidate_optimizer_status_counts[str(record["optimizer_status"])] += 1
            candidate_slsqp_status_counts[str(record["slsqp_status"])] += 1
            globally_certified = (
                bool(record["fit_converged"])
                and float(record["global_log_likelihood_gap_bound"])
                <= OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE
                and float(record["simplex_residual"]) <= 1e-12
                and bool(record["fit_monotone"])
            )
            candidate_globally_certified_count += int(globally_certified)
            candidate_nonidentifiable_weight_count += int(
                not bool(record["mixture_weights_identifiable"])
            )
            candidate_max_global_ll_gap_bound = max(
                candidate_max_global_ll_gap_bound,
                float(record["global_log_likelihood_gap_bound"]),
            )
            candidate_max_simplex_residual = max(
                candidate_max_simplex_residual,
                float(record["simplex_residual"]),
            )
            all_converged_candidate_certificates_valid = (
                all_converged_candidate_certificates_valid
                and (not bool(record["fit_converged"]) or globally_certified)
            )
            all_optimizer_abstain_candidates_not_winners = (
                all_optimizer_abstain_candidates_not_winners
                and (bool(record["fit_converged"]) or not bool(record["is_winner"]))
            )
            all_candidate_hashes_valid = all_candidate_hashes_valid and (
                len(record["candidate_semantic_sha256"]) == 64
            )

        def emit_responsibility(record: dict[str, Any]) -> None:
            nonlocal responsibility_row_count, all_responsibility_hashes_valid
            responsibility_writer.writerow(
                {key: _cell(record.get(key)) for key in RESPONSIBILITY_FIELDS}
            )
            responsibility_semantic.update(record)
            responsibility_row_count += 1
            all_responsibility_hashes_valid = all_responsibility_hashes_valid and (
                len(record["responsibility_semantic_sha256"]) == 64
            )

        for current_minread in minread_grid:
            role = (
                "PRIMARY"
                if current_minread == primary_structural_exact_pattern_minread
                else "SENSITIVITY"
            )
            current_rows = []
            for key in ordered_unit_keys:
                for pattern_row in _pattern_rows(
                    units[key], current_minread, minread_role=role
                ):
                    pattern_writer.writerow(
                        {field: _cell(pattern_row.get(field)) for field in PATTERN_FIELDS}
                    )
                    pattern_semantic.update(pattern_row)
                    all_pattern_hashes_valid = all_pattern_hashes_valid and (
                        len(pattern_row["pattern_semantic_sha256"]) == 64
                    )
                row = rank_unit(
                    units[key],
                    minread=current_minread,
                    minread_role=role,
                    exact_k_max=exact_k_max,
                    max_vertex_sets=max_vertex_sets,
                    solver_time_limit_seconds=solver_time_limit_seconds,
                    fixed_error_grid=fixed_error_grid,
                    minimum_bq_error_rate=minimum_bq_error_rate,
                    maximum_bq_error_rate=maximum_bq_error_rate,
                    conditional_candidate_ranking_bootstrap_replicates=(
                        conditional_candidate_ranking_bootstrap_replicates
                        if role == "PRIMARY"
                        else 0
                    ),
                    conditional_candidate_ranking_bootstrap_seed=(
                        conditional_candidate_ranking_bootstrap_seed
                    ),
                    tie_tolerance=tie_tolerance,
                    structural_enumeration_cache=structural_enumeration_cache,
                    candidate_record_sink=emit_candidate,
                    responsibility_record_sink=(emit_responsibility if role == "PRIMARY" else None),
                    retain_detail_records=False,
                )
                runtime_record = {
                    "dataset": row["dataset"],
                    "chrom": row["chrom"],
                    "threshold": row["threshold"],
                    "component_basis": row["component_basis"],
                    "phase_set": row["phase_set"],
                    "component_id": row["component_id"],
                    "family": row["family"],
                    "structural_exact_pattern_minread": row[
                        "structural_exact_pattern_minread"
                    ],
                    "structural_minread_role": row["structural_minread_role"],
                    **row["_runtime_segment_invoked"],
                    **row["_runtime_diagnostics"],
                }
                runtime_writer.writerow(
                    {
                        field: (
                            format(float(runtime_record[field]), ".17g")
                            if field in RUNTIME_METRICS
                            else _cell(runtime_record.get(field))
                        )
                        for field in RUNTIME_DIAGNOSTIC_FIELDS
                    }
                )
                current_rows.append(row)
            unit_rows_by_minread[current_minread] = current_rows
            unit_rows.extend(current_rows)
    with gzip.open(staged_unit_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, UNIT_FIELDS, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in unit_rows:
            public_row = {key: row.get(key) for key in UNIT_FIELDS}
            writer.writerow({key: _cell(public_row.get(key)) for key in UNIT_FIELDS})
            unit_semantic.update(public_row)
    primary_rows = unit_rows_by_minread[primary_structural_exact_pattern_minread]
    structural_cache_runtime = (
        structural_enumeration_cache.diagnostics()
        if structural_enumeration_cache is not None
        else {
            "schema_name": "intersubmod.structural_enumeration_cache_diagnostics",
            "schema_version": STRUCTURAL_CACHE_SCHEMA_VERSION,
            "enabled": False,
            "lifecycle": "DISABLED",
            "complete_only": True,
            "full_frozen_tuple_is_equality_authority": True,
            "diagnostic_sha256_is_not_equality_authority": True,
            "likelihood_cached": False,
            "lookups": 0,
            "hits": 0,
            "misses": 0,
            "stores_complete": 0,
            "rejected_incomplete": 0,
            "solver_invocations": 0,
            "entries_final": 0,
            "evictions": 0,
            "cross_minread_hits": 0,
            "cross_threshold_hits": 0,
            "same_unit_cross_minread_hits": 0,
            "digest_collision_guard_events": 0,
        }
    )
    runtime_diagnostics = {
        "schema_name": "intersubmod.m2_unit_runtime_diagnostics",
        "schema_version": "1.0.0",
        "clock": "time.monotonic_ns",
        "unit": "seconds",
        "quantile_definition": (
            "exact empirical nearest-rank: rank=ceil(p*n), one-based, for p in {0.50,0.95,0.99}"
        ),
        "interpretation": (
            "process-local monotonic wall-clock performance diagnostic; values are environment/load dependent "
            "and are not scientific results or cross-machine reproducibility claims"
        ),
        "segment_definition": {
            "candidate_generation_elapsed_seconds": (
                "wall time inside the process-local structural-cache lookup/deepcopy "
                "or exact candidate vertex-set enumeration"
            ),
            "likelihood_fit_elapsed_seconds": (
                "wall time inside primary BQ mixture fits, fixed-error sensitivity fits, and the conditional "
                "candidate-ranking bootstrap operation (deterministic resampling plus fits)"
            ),
            "unit_total_elapsed_seconds": (
                "wall time for the complete rank_unit call, including candidate/result serialization sinks"
            ),
        },
        "per_unit_output": runtime_path.name,
        "structural_enumeration_cache": structural_cache_runtime,
        "scopes": {
            "primary_unit_evaluations": summarize_runtime_rows(primary_rows),
            "all_structural_minread_unit_evaluations": summarize_runtime_rows(unit_rows),
        },
        "primary_invoked_segment_scopes": {
            "candidate_generation_elapsed_seconds": summarize_runtime_values(
                row["_runtime_diagnostics"]["candidate_generation_elapsed_seconds"]
                for row in primary_rows
                if row["_runtime_segment_invoked"]["candidate_generation_invoked"]
            ),
            "likelihood_fit_elapsed_seconds": summarize_runtime_values(
                row["_runtime_diagnostics"]["likelihood_fit_elapsed_seconds"]
                for row in primary_rows
                if row["_runtime_segment_invoked"]["likelihood_fit_invoked"]
            ),
        },
    }
    aggregate = summarize_units(primary_rows)
    aggregate_by_basis_threshold = summarize_units_by_basis_threshold(primary_rows)
    sensitivity_by_minread = {
        str(value): {
            "role": "PRIMARY" if value == primary_structural_exact_pattern_minread else "SENSITIVITY",
            "aggregate": summarize_units(unit_rows_by_minread[value]),
            "by_linkage_basis_threshold": summarize_units_by_basis_threshold(unit_rows_by_minread[value]),
        }
        for value in minread_grid
    }
    partial_funnel = partial_pattern_funnel(
        [units[key] for key in ordered_unit_keys], primary_rows, primary_structural_exact_pattern_minread
    )
    partial_funnel_by_basis_threshold: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for basis in sorted({unit.component.basis for unit in units.values()}):
        for threshold in sorted({
            unit.component.threshold for unit in units.values() if unit.component.basis == basis
        }):
            selected_units = [
                units[key] for key in ordered_unit_keys
                if units[key].component.basis == basis and units[key].component.threshold == threshold
            ]
            selected_rows = [
                row for row in primary_rows
                if row["component_basis"] == basis and int(row["threshold"]) == threshold
            ]
            partial_funnel_by_basis_threshold[basis][str(threshold)] = partial_pattern_funnel(
                selected_units, selected_rows, primary_structural_exact_pattern_minread
            )
    for value in minread_grid:
        sensitivity_by_minread[str(value)]["partial_pattern_funnel"] = partial_pattern_funnel(
            [units[key] for key in ordered_unit_keys], unit_rows_by_minread[value], value
        )
        by_basis: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
        for basis in sorted({unit.component.basis for unit in units.values()}):
            for threshold in sorted({
                unit.component.threshold for unit in units.values() if unit.component.basis == basis
            }):
                selected_units = [
                    units[key] for key in ordered_unit_keys
                    if units[key].component.basis == basis
                    and units[key].component.threshold == threshold
                ]
                selected_rows = [
                    row for row in unit_rows_by_minread[value]
                    if row["component_basis"] == basis and int(row["threshold"]) == threshold
                ]
                by_basis[basis][str(threshold)] = partial_pattern_funnel(
                    selected_units, selected_rows, value
                )
        sensitivity_by_minread[str(value)]["partial_pattern_funnel_by_linkage_basis_threshold"] = dict(by_basis)
    checks = {
        "upstream_receipt_all_pass": True,
        "component_site_membership_conserved": all(
            len(component.site_indices) == component.k for component in components.values()
        ),
        "projection_funnel_conserved": aggregate["molecule_component_projections"]
        == aggregate["informative_scoring_molecules"] + aggregate["all_x_excluded_molecules"],
        "structural_scoring_separation_conserved": aggregate["informative_scoring_molecules"]
        == aggregate["structural_retained_molecules"] + aggregate["below_minread_scoring_molecules"],
        "all_units_have_semantic_hash": all(len(row["unit_semantic_sha256"]) == 64 for row in unit_rows),
        "all_patterns_have_semantic_hash": all_pattern_hashes_valid,
        "k_gt_exact_limit_never_claimed_global_optimal": True,
        "all_informative_scoring_molecules_use_bq": aggregate["bq_scoring_molecules"]
        == aggregate["informative_scoring_molecules"],
        "primary_units_are_exact_known_ps": all(
            row["component_basis"] in PRIMARY_COMPONENT_BASES
            and row["phase_set"] not in {None, ""}
            and row["phase_set_status"] == "KNOWN_PS_PRIMARY"
            and row["inference_role"] == "PRIMARY_PS_AWARE"
            for row in primary_rows
        ) if basis_mode == "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY" else True,
        "no_cross_ps_pattern_pooling": all(
            row["phase_set"] not in {None, ""} and row["component_basis"] in PRIMARY_COMPONENT_BASES
            for row in primary_rows
        ) if basis_mode == "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY" else True,
        "known_ps_never_mixed": all(
            row["phase_set_status"] == "KNOWN_PS_PRIMARY" for row in primary_rows
        ) if basis_mode == "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY" else True,
        "missing_ps_separate_diagnostic": all(
            row["phase_set"] not in {None, ""} for row in primary_rows
        ) if basis_mode == "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY" else True,
        "partial_group_constraints_all_satisfied_when_evaluated": (
            aggregate["partial_groups_unsatisfied"] == 0
            and aggregate["partial_groups_covered"] == aggregate["partial_group_coverage_denominator"]
        ),
        "effective_k_site_mass_conserved": aggregate["k_component_sites_total"]
        == aggregate["k_observed_alt_active_total"]
        + aggregate["not_structural_alt_active_sites_total"],
        "topology_partition_conserved": aggregate["coarse_topology_class_unique_units"]
        + aggregate["coarse_topology_multiple_class_units"] == aggregate["topology_evaluated_units"],
        "k_route_partition_conserved": sum(aggregate["k_route_counts"].values())
        == aggregate["n_component_hp_units"],
        "all_candidates_have_semantic_hash": all_candidate_hashes_valid,
        "all_converged_candidate_fits_have_global_kkt_certificate": (
            all_converged_candidate_certificates_valid
        ),
        "optimizer_abstain_candidates_are_not_claimed_winners": (
            all_optimizer_abstain_candidates_not_winners
        ),
        "all_responsibilities_have_semantic_hash": all_responsibility_hashes_valid,
        "candidate_relative_weights_conserve_per_unit": True,
        "posterior_responsibilities_conserve_per_pattern_candidate": True,
        "all_unit_evaluations_have_runtime_diagnostics": all(
            set(row.get("_runtime_diagnostics") or {}) == set(RUNTIME_METRICS)
            for row in unit_rows
        ),
        "runtime_diagnostics_are_finite_nonnegative_and_segmented": all(
            all(
                math.isfinite(float(row["_runtime_diagnostics"][metric]))
                and float(row["_runtime_diagnostics"][metric]) >= 0.0
                for metric in RUNTIME_METRICS
            )
            and (
                float(row["_runtime_diagnostics"]["candidate_generation_elapsed_seconds"])
                + float(row["_runtime_diagnostics"]["likelihood_fit_elapsed_seconds"])
                <= float(row["_runtime_diagnostics"]["unit_total_elapsed_seconds"]) + 1e-9
            )
            for row in unit_rows
        ),
        "runtime_scope_counts_match_unit_rows": (
            runtime_diagnostics["scopes"]["primary_unit_evaluations"]["n_unit_evaluations"]
            == len(primary_rows)
            and runtime_diagnostics["scopes"]["all_structural_minread_unit_evaluations"]
            ["n_unit_evaluations"] == len(unit_rows)
        ),
        "runtime_values_excluded_from_scientific_semantic_hashes": True,
        "runtime_invocation_flags_match_segment_contract": all(
            row["_runtime_segment_invoked"]["likelihood_fit_invoked"]
            <= row["_runtime_segment_invoked"]["candidate_generation_invoked"]
            for row in unit_rows
        ),
        "structural_cache_accounting_conserved": (
            structural_cache_runtime["lookups"]
            == structural_cache_runtime["hits"] + structural_cache_runtime["misses"]
            and structural_cache_runtime["solver_invocations"]
            == structural_cache_runtime["misses"]
            and structural_cache_runtime["misses"]
            == structural_cache_runtime["stores_complete"]
            + structural_cache_runtime["rejected_incomplete"]
            and structural_cache_runtime["entries_final"]
            == structural_cache_runtime["stores_complete"]
            - structural_cache_runtime["evictions"]
        ),
        "structural_cache_is_complete_only_and_likelihood_free": (
            structural_cache_runtime["complete_only"] is True
            and structural_cache_runtime["likelihood_cached"] is False
        ),
        "structural_cache_digest_collision_guard_clean": (
            structural_cache_runtime["digest_collision_guard_events"] == 0
        ),
        "structural_cache_per_unit_diagnostics_present": all(
            set(row.get("_structural_cache_diagnostics") or {})
            == {
                "structural_cache_outcome",
                "structural_cache_key_sha256",
                "structural_cache_reuse_relation",
                "structural_solver_invoked",
            }
            for row in unit_rows
        ),
    }
    receipt = {
        "schema_name": "intersubmod.m2_symbolic_patterns_vertex_rank_receipt",
        "schema_version": RANKING_SCHEMA_VERSION,
        "provenance": {
            "ranker": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256_path(Path(__file__).resolve()),
            },
            "upstream_extraction_receipt": {
                "path": upstream["path"],
                "sha256": upstream["sha256"],
            },
        },
        "scope": {
            "dataset": [str(upstream.get("scope", {}).get("dataset"))],
            "chrom": [str(upstream.get("scope", {}).get("chrom"))],
            "thresholds": sorted(thresholds or {component.threshold for component in components.values()}),
            "component_bases": sorted(selected_bases),
            "component_basis_mode": basis_mode,
            "hp_families": list(families),
            "phase_set_unit": "exact known PS for primary; missing PS excluded",
        },
        "parameters": {
            "method_contract": METHOD_CONTRACT,
            "structural_exact_pattern_minread_grid": list(minread_grid),
            "primary_structural_exact_pattern_minread": primary_structural_exact_pattern_minread,
            "scoring_minread": 1,
            "exact_k_max": exact_k_max,
            "max_vertex_sets": max_vertex_sets,
            "solver_time_limit_seconds_per_milp": solver_time_limit_seconds,
            "primary_likelihood": (
                "per-fixed-R/A-call Phred BQ symmetric-substitution emission conditioned on observed in {R,A}; "
                "O/D/S/L/X marginalized as missing"
            ),
            "minimum_bq_error_rate": minimum_bq_error_rate,
            "maximum_bq_error_rate": maximum_bq_error_rate,
            "fixed_error_grid_conditional_binary_flip_probability": list(fixed_error_grid),
            "conditional_candidate_ranking_bootstrap_replicates": conditional_candidate_ranking_bootstrap_replicates,
            "conditional_candidate_ranking_bootstrap_base_seed": conditional_candidate_ranking_bootstrap_seed,
            "conditional_candidate_ranking_bootstrap_policy": (
                "deterministic molecule multinomial resampling conditional on the already fixed candidate vertex-set list; "
                "does not rebuild structural candidates and therefore is not whole-tree stability"
            ),
            "tie_tolerance_log_likelihood": tie_tolerance,
            "structural_enumeration_cache": {
                "enabled": enable_structural_enumeration_cache,
                "schema_version": STRUCTURAL_CACHE_SCHEMA_VERSION,
                "lifecycle": "PROCESS_LOCAL_ONE_RANKING_CHILD",
                "entry_policy": "COMPLETE_ONLY_WITH_STRUCTURAL_INVARIANT_REVALIDATION",
                "key_equality_authority": "FULL_FROZEN_TUPLE",
                "diagnostic_sha256_is_equality_authority": False,
                "store_copy": "DEEP_COPY",
                "hit_copy": "DEEP_COPY",
                "likelihood_cached": False,
            },
            "optimizer_contract": {
                "warm_start": "SLSQP from uniform simplex weights",
                "refinement": (
                    "deterministic monotone pairwise mass transfer with exact concave line search"
                ),
                "convergence_authority": (
                    "Frank-Wolfe/KKT global log-likelihood gap bound plus simplex residual; "
                    "scipy success is diagnostic only"
                ),
                "global_log_likelihood_gap_tolerance": (
                    OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE
                ),
                "simplex_residual_tolerance": 1e-12,
                "same_read_vaf_added": False,
            },
            "structural_pattern_definition": (
                "the same exact R/A/X pattern count >= structural_exact_pattern_minread; compatible partial "
                "patterns are not pooled (AXR/ARX/AXX count separately)"
            ),
            "compatible_pattern_pooling": "NOT_IMPLEMENTED_METHOD_SENSITIVITY_OPTION",
            "scoring_pattern_definition": "all informative molecule projections before minread",
            "all_x_policy": "conserve in denominator, omit exact constant likelihood term",
            "call_emission_policy": {
                "R_A": (
                    "BQ e is any-base error; symmetric substitution and condition observed in {R,A}: "
                    "match=(1-e)/(1-2e/3), flip=(e/3)/(1-2e/3)"
                ),
                "O_D_S_L_X": "preserve reason count, marginalize as missing with emission factor 1",
                "claim_ceiling": (
                    "symmetric substitution and conditioning are explicit assumptions; no calibrated context/indel/"
                    "other-base emission model in schema v2"
                ),
            },
            "component_basis_policy": (
                "PS_HP1/PS_HP2 x exact known PS membership is primary; different PS never exchange support; "
                "missing PS and pooled/legacy HP bases are non-primary"
            ),
            "hp_family_merge_policy": {
                "family_1": ["1", "1-1", "1-2"],
                "family_2": ["2", "2-1", "2-2"],
                "source_column": "extractor hp_family (raw fine tag remains upstream)",
                "circularity_caveat": (
                    "somatic-derived fine tags such as 1-1/2-1 are merged into their germline family; "
                    "they are not independent likelihood evidence and are not clone labels"
                ),
            },
            "partial_policy": (
                "one symbolic subcube constraint; conceptual 2^u completion-wise tree worlds are counted but not "
                "materialized; active compatible vertex indices are materialized per MILP construction for each "
                "reduced sparse group-hit row; primary effective completion count uses minread-specific structural "
                "observed-alt universe"
            ),
            "effective_k_policy": (
                "exact gate and MILP vertex universe use minread-specific k_observed_alt_active; other component "
                "dimensions (R-only, X-only, or ALT only below minread) remain coordinate provenance and are "
                "candidate-constant under that structural universe; k_scoring_alt_observed is reported separately"
            ),
            "edge_policy": "collapse by vertex_set_id; snapshot patterns do not identify parent edges",
            "topology_policy": {
                "classes": ["single-only", "sister-only", "direct-only", "sister+direct"],
                "basis": "all selected states including minimum-hidden states",
                "coarse_topology_class_unique": (
                    "exactly one of four broad classes across all winning vertex sets and parent choices; "
                    "this is not canonical exact topology uniqueness"
                ),
                "parent_edge_assignment_unique": "exactly one labeled parent-edge assignment among winners",
                "exact_topology_unique": (
                    "true only when a single winning parent-edge assignment proves uniqueness; otherwise NA "
                    "until canonical exact-shape isomorphism is evaluated"
                ),
                "algorithm": (
                    "depth from Hamming weight plus branch/no-branch feasibility by bipartite matching; "
                    "no parent-edge Cartesian enumeration"
                ),
            },
            "vaf_policy": "same-molecule VAF is not added to read-pattern likelihood",
        },
        "upstream": upstream,
        "input_files": {
            path.name: {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_path(path)}
            for path in (calls_path, sites_path, components_path, membership_path)
        },
        "input_semantic_hashes": {
            "site_catalog_semantic_sha256": site_semantic,
            **input_semantics,
            "calls_semantic_sha256": input_counts["calls_semantic_sha256"],
        },
        "input_counts": input_counts,
        "aggregate": aggregate,
        "aggregate_by_linkage_basis_threshold": aggregate_by_basis_threshold,
        "sensitivity_by_structural_exact_pattern_minread": sensitivity_by_minread,
        "partial_pattern_funnel": partial_funnel,
        "partial_pattern_funnel_by_linkage_basis_threshold": dict(partial_funnel_by_basis_threshold),
        "candidate_evidence_counts": {
            "compressed_vertex_set_candidate_rows": candidate_row_count,
            "winning_candidate_rows": winning_candidate_row_count,
            "posterior_responsibility_rows": responsibility_row_count,
            "optimizer_diagnostics": {
                "candidate_fits": candidate_row_count,
                "globally_kkt_certified_candidate_fits": (
                    candidate_globally_certified_count
                ),
                "nonidentifiable_mixture_weight_candidates": (
                    candidate_nonidentifiable_weight_count
                ),
                "optimizer_status_counts": dict(
                    sorted(candidate_optimizer_status_counts.items())
                ),
                "slsqp_status_counts": dict(sorted(candidate_slsqp_status_counts.items())),
                "max_global_log_likelihood_gap_bound": (
                    candidate_max_global_ll_gap_bound
                ),
                "max_simplex_residual": candidate_max_simplex_residual,
                "claim_ceiling": (
                    "pi is a fitted latent mutation-state exposure proportion; if augmented "
                    "emission rank is deficient it is not uniquely identifiable, and it is "
                    "never a cellular clone fraction"
                ),
            },
            "posterior_interpretation": "latent expectation; never a hard read/cell/clone assignment",
            "posterior_scope": "winning candidates at the primary structural exact-pattern minread only",
        },
        "runtime_diagnostics": runtime_diagnostics,
        "checks": checks,
        "all_pass": all(checks.values()),
        "outputs": {
            pattern_path.name: {
                "path": str(pattern_path),
                "size_bytes": staged_pattern_path.stat().st_size,
                "sha256": sha256_path(staged_pattern_path),
                "semantic_sha256": pattern_semantic.hexdigest(),
            },
            unit_path.name: {
                "path": str(unit_path),
                "size_bytes": staged_unit_path.stat().st_size,
                "sha256": sha256_path(staged_unit_path),
                "semantic_sha256": unit_semantic.hexdigest(),
            },
            candidate_path.name: {
                "path": str(candidate_path),
                "size_bytes": staged_candidate_path.stat().st_size,
                "sha256": sha256_path(staged_candidate_path),
                "semantic_sha256": candidate_semantic.hexdigest(),
            },
            responsibility_path.name: {
                "path": str(responsibility_path),
                "size_bytes": staged_responsibility_path.stat().st_size,
                "sha256": sha256_path(staged_responsibility_path),
                "semantic_sha256": responsibility_semantic.hexdigest(),
            },
            runtime_path.name: {
                "path": str(runtime_path),
                "size_bytes": staged_runtime_path.stat().st_size,
                "sha256": sha256_path(staged_runtime_path),
                "semantic_sha256": None,
                "semantic_hash_policy": (
                    "omitted because monotonic wall timings are deliberately non-reproducible diagnostics"
                ),
            },
        },
    }
    receipt_path = output_dir / "receipt.json"
    staged_receipt_path = staging_dir / receipt_path.name
    receipt["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": f"{receipt_path.name}.sha256",
        "covers": receipt_path.name,
    }
    staged_receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_sha256_sidecar(staged_receipt_path)
    if not require_existing_empty_output_dir:
        staging_dir.rename(output_dir)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threshold", action="append", type=int)
    parser.add_argument("--component-basis", action="append")
    parser.add_argument("--hp-family", action="append", choices=("1", "2", "3", "4", "none", "pooled"))
    parser.add_argument("--structural-exact-pattern-minread-grid", default="1,2,3,5")
    parser.add_argument("--primary-structural-exact-pattern-minread", type=int, default=3)
    parser.add_argument("--minread", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--exact-k-max", type=int, default=12)
    parser.add_argument("--max-vertex-sets", type=int, default=256)
    parser.add_argument("--solver-time-limit-seconds", type=float, default=30.0)
    parser.add_argument("--fixed-error-grid", default="0.005,0.01,0.02,0.05")
    parser.add_argument("--minimum-bq-error-rate", type=float, default=1e-6)
    parser.add_argument("--maximum-bq-error-rate", type=float, default=0.25)
    parser.add_argument(
        "--conditional-candidate-ranking-bootstrap-replicates",
        "--bootstrap-replicates",
        dest="conditional_candidate_ranking_bootstrap_replicates",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--conditional-candidate-ranking-bootstrap-seed",
        "--bootstrap-seed",
        dest="conditional_candidate_ranking_bootstrap_seed",
        type=int,
        default=20260716,
    )
    parser.add_argument("--tie-tolerance", type=float, default=1e-6)
    parser.add_argument(
        "--disable-structural-enumeration-cache",
        action="store_true",
        help=(
            "disable the process-local complete-only structural cache for an "
            "exact-equivalence control run"
        ),
    )
    parser.add_argument(
        "--require-existing-empty-output-dir",
        action="store_true",
        help=(
            "write into a real empty output directory already created and inode-bound "
            "by the release-pilot resource gate"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    receipt = run(
        args.input_dir,
        args.output_dir,
        thresholds=None if args.threshold is None else set(args.threshold),
        component_bases=None if args.component_basis is None else set(args.component_basis),
        families=tuple(args.hp_family or PRIMARY_HP_FAMILIES),
        minread=args.primary_structural_exact_pattern_minread if args.minread is None else args.minread,
        structural_exact_pattern_minreads=(
            (args.minread,) if args.minread is not None
            else tuple(int(value) for value in args.structural_exact_pattern_minread_grid.split(",") if value)
        ),
        primary_structural_exact_pattern_minread=(
            args.minread if args.minread is not None else args.primary_structural_exact_pattern_minread
        ),
        exact_k_max=args.exact_k_max,
        max_vertex_sets=args.max_vertex_sets,
        solver_time_limit_seconds=args.solver_time_limit_seconds,
        fixed_error_grid=tuple(float(value) for value in args.fixed_error_grid.split(",") if value),
        minimum_bq_error_rate=args.minimum_bq_error_rate,
        maximum_bq_error_rate=args.maximum_bq_error_rate,
        conditional_candidate_ranking_bootstrap_replicates=args.conditional_candidate_ranking_bootstrap_replicates,
        conditional_candidate_ranking_bootstrap_seed=args.conditional_candidate_ranking_bootstrap_seed,
        tie_tolerance=args.tie_tolerance,
        enable_structural_enumeration_cache=(
            not args.disable_structural_enumeration_cache
        ),
        require_existing_empty_output_dir=args.require_existing_empty_output_dir,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    raise SystemExit(0 if receipt["all_pass"] else 1)


if __name__ == "__main__":
    main()
