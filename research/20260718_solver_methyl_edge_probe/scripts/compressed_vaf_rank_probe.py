#!/usr/bin/env python3
"""Isolated exact/lazy probe for the current BQ-aware vertex-set likelihood.

This module does not modify the canonical M2 router.  It applies only to the
recurrence-free perfect-family layer certified by ``perfect_family_count``.
The structural family is represented by the same subset-DP traceback DAG and
is traversed lazily.  A branch may be pruned only with a mathematically safe
likelihood upper bound.

The scored endpoint is exactly the primary endpoint used by
``build_m2_patterns_and_rank.fit_quality_aware_mixture``:

    max_pi sum_i n_i log(sum_{v in V} pi_v Q_iv).

It does not reproduce the fixed-error grid, bootstrap, responsibilities, or
one-output-row-per-candidate release contract.
"""

from __future__ import annotations

import dataclasses
import hashlib
import itertools
import json
import math
import sys
import time
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

import numpy as np


BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parents[1]
M2_SCRIPTS = (
    REPO
    / "research"
    / "20260716_read_linked_hypercube_exact_likelihood_validation"
    / "scripts"
)
if str(M2_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(M2_SCRIPTS))

from build_m2_patterns_and_rank import (  # noqa: E402
    OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE,
    _quality_emission_matrix,
    fit_quality_aware_mixture,
)
from hypercube_exact import vertex_set_digest  # noqa: E402


SCHEMA = "intersubmod.compressed_vaf_rank_probe.v1"
SCOPE = "ISOLATED_PERFECT_FAMILY_PRIMARY_BQ_LIKELIHOOD_NOT_PRODUCTION"
QualityGroup = tuple[str, tuple[int, ...], int]


def _popcount(value: int) -> int:
    return int(value).bit_count()


def _validate_pattern(pattern: str, *, k: int, allow_x: bool) -> None:
    if len(pattern) != k:
        raise ValueError("pattern length differs from k")
    allowed = {"R", "A", "X"} if allow_x else {"R", "A"}
    invalid = set(pattern) - allowed
    if invalid:
        raise ValueError(f"invalid pattern symbols: {sorted(invalid)}")


def _induced_components(mask: int, adjacency: Sequence[int]) -> tuple[int, ...]:
    components: list[int] = []
    remaining = int(mask)
    while remaining:
        component = 0
        frontier = remaining & -remaining
        while frontier:
            component |= frontier
            neighbours = 0
            cursor = frontier
            while cursor:
                bit = cursor & -cursor
                neighbours |= adjacency[bit.bit_length() - 1]
                cursor ^= bit
            frontier = neighbours & mask & ~component
        components.append(component)
        remaining &= ~component
    return tuple(components)


def _merge_parent_vectors(
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> tuple[int, ...]:
    merged = []
    for left_value, right_value in zip(left, right):
        if left_value != -2 and right_value != -2:
            raise AssertionError("traceback subforests overlap")
        merged.append(left_value if left_value != -2 else right_value)
    return tuple(merged)


@dataclasses.dataclass(frozen=True)
class TraceBranch:
    """One disjoint top-level traceback-DAG branch."""

    first_block: int
    first_root: int
    remaining_mask: int
    completion_count: int
    possible_vertices: tuple[int, ...]

    def identity(self) -> str:
        payload = dataclasses.asdict(self)
        return hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()


class PerfectTracebackDag:
    """Count and lazily enumerate valid recurrence-free mutation forests."""

    def __init__(
        self,
        *,
        full_patterns: Sequence[str],
        partial_patterns: Sequence[str],
        k: int,
        structural_alt_universe_mask: int | None = None,
        max_m: int = 20,
    ) -> None:
        if type(k) is not int or k < 1:
            raise ValueError("k must be a positive integer")
        if type(max_m) is not int or max_m < 1:
            raise ValueError("max_m must be a positive integer")
        self.k = k
        for pattern in full_patterns:
            _validate_pattern(pattern, k=self.k, allow_x=False)
        for pattern in partial_patterns:
            _validate_pattern(pattern, k=self.k, allow_x=True)
        observed_alt_mask = 0
        for pattern in tuple(full_patterns) + tuple(partial_patterns):
            for bit, symbol in enumerate(pattern):
                if symbol == "A":
                    observed_alt_mask |= 1 << bit
        if structural_alt_universe_mask is None:
            universe_mask = observed_alt_mask
        else:
            if type(structural_alt_universe_mask) is not int:
                raise ValueError("structural ALT universe must be an integer")
            universe_mask = structural_alt_universe_mask
            if universe_mask != observed_alt_mask:
                raise ValueError(
                    "structural ALT universe differs from observed ALT union"
                )
        self.universe_mask = universe_mask
        self.active_bits = tuple(
            bit for bit in range(self.k) if universe_mask & (1 << bit)
        )
        self.m = len(self.active_bits)
        if self.m < 1:
            raise ValueError("perfect-family ranking requires active mutations")
        if self.m > max_m:
            raise ValueError(f"effective_m={self.m} exceeds gate {max_m}")
        compact_index = {
            original: compact for compact, original in enumerate(self.active_bits)
        }
        adjacency = [0 for _ in range(self.m)]
        forbidden = [0 for _ in range(self.m)]
        for pattern in tuple(full_patterns) + tuple(partial_patterns):
            alt = [
                compact_index[bit]
                for bit, symbol in enumerate(pattern)
                if symbol == "A"
            ]
            ref = [
                compact_index[bit]
                for bit, symbol in enumerate(pattern)
                if symbol == "R" and bit in compact_index
            ]
            for left_index, left in enumerate(alt):
                for right in alt[left_index + 1 :]:
                    adjacency[left] |= 1 << right
                    adjacency[right] |= 1 << left
            for ancestor in ref:
                for descendant in alt:
                    forbidden[ancestor] |= 1 << descendant
        self.adjacency = tuple(adjacency)
        self.forbidden = tuple(forbidden)
        self.full_mask = (1 << self.m) - 1
        if self.forest_count(self.full_mask) < 1:
            raise ValueError(
                "perfect traceback family is empty; recurrence-aware layer required"
            )

    def _component_blocks(self, mask: int) -> Iterator[int]:
        components = _induced_components(mask, self.adjacency)
        first = components[0]
        others = components[1:]
        for choice in range(1 << len(others)):
            block = first
            for index, component in enumerate(others):
                if choice & (1 << index):
                    block |= component
            yield block

    @lru_cache(maxsize=None)
    def tree_count(self, mask: int) -> int:
        if not mask:
            return 0
        total = 0
        roots = mask
        while roots:
            root_bit = roots & -roots
            root = root_bit.bit_length() - 1
            descendants = mask ^ root_bit
            if not (self.forbidden[root] & descendants):
                total += self.forest_count(descendants)
            roots ^= root_bit
        return total

    @lru_cache(maxsize=None)
    def forest_count(self, mask: int) -> int:
        if not mask:
            return 1
        total = 0
        for block in self._component_blocks(mask):
            total += self.tree_count(block) * self.forest_count(mask ^ block)
        return total

    def _empty_parent_vector(self) -> tuple[int, ...]:
        return tuple(-2 for _ in range(self.m))

    def _iter_tree(self, mask: int) -> Iterator[tuple[int, ...]]:
        roots = mask
        while roots:
            root_bit = roots & -roots
            root = root_bit.bit_length() - 1
            descendants = mask ^ root_bit
            roots ^= root_bit
            if self.forbidden[root] & descendants:
                continue
            for child_forest in self._iter_forest(descendants):
                values = list(child_forest)
                values[root] = -1
                for mutation in range(self.m):
                    if descendants & (1 << mutation) and values[mutation] == -1:
                        values[mutation] = root
                yield tuple(values)

    def _iter_forest(self, mask: int) -> Iterator[tuple[int, ...]]:
        if not mask:
            yield self._empty_parent_vector()
            return
        for block in self._component_blocks(mask):
            for tree in self._iter_tree(block):
                for remainder in self._iter_forest(mask ^ block):
                    yield _merge_parent_vectors(tree, remainder)

    def _expand_state(self, compact_state: int) -> int:
        expanded = 0
        for compact_bit, original_bit in enumerate(self.active_bits):
            if compact_state & (1 << compact_bit):
                expanded |= 1 << original_bit
        return expanded

    def vertex_set(
        self,
        parents: Sequence[int],
    ) -> tuple[int, ...]:
        if len(parents) != self.m:
            raise ValueError("parent vector length differs from effective m")
        memo: dict[int, int] = {}
        visiting: set[int] = set()

        def state(mutation: int) -> int:
            if mutation in memo:
                return memo[mutation]
            if mutation in visiting:
                raise ValueError("parent vector contains a cycle")
            visiting.add(mutation)
            parent = int(parents[mutation])
            if parent == -1:
                value = 1 << mutation
            elif 0 <= parent < self.m and parent != mutation:
                value = state(parent) | (1 << mutation)
            else:
                raise ValueError("invalid selected parent")
            visiting.remove(mutation)
            memo[mutation] = value
            return value

        states = {0}
        states.update(self._expand_state(state(mutation)) for mutation in range(self.m))
        if len(states) != self.m + 1:
            raise AssertionError("perfect forest did not produce m+1 distinct states")
        return tuple(sorted(states))

    def branches(self) -> tuple[TraceBranch, ...]:
        branches: list[TraceBranch] = []
        for block in self._component_blocks(self.full_mask):
            remainder = self.full_mask ^ block
            roots = block
            while roots:
                root_bit = roots & -roots
                root = root_bit.bit_length() - 1
                roots ^= root_bit
                child_mask = block ^ root_bit
                if self.forbidden[root] & child_mask:
                    continue
                completion_count = (
                    self.forest_count(child_mask)
                    * self.forest_count(remainder)
                )
                if not completion_count:
                    continue
                possible_compact = {0}
                subset = child_mask
                while True:
                    possible_compact.add(root_bit | subset)
                    if subset == 0:
                        break
                    subset = (subset - 1) & child_mask
                subset = remainder
                while subset:
                    possible_compact.add(subset)
                    subset = (subset - 1) & remainder
                branches.append(
                    TraceBranch(
                        first_block=block,
                        first_root=root,
                        remaining_mask=remainder,
                        completion_count=completion_count,
                        possible_vertices=tuple(
                            sorted(self._expand_state(value) for value in possible_compact)
                        ),
                    )
                )
        if sum(branch.completion_count for branch in branches) != self.forest_count(
            self.full_mask
        ):
            raise AssertionError("top-level traceback branches do not conserve count")
        return tuple(branches)

    def iter_branch(self, branch: TraceBranch) -> Iterator[tuple[int, ...]]:
        root_bit = 1 << branch.first_root
        child_mask = branch.first_block ^ root_bit
        produced = 0
        for child_forest in self._iter_forest(child_mask):
            tree_values = list(child_forest)
            tree_values[branch.first_root] = -1
            for mutation in range(self.m):
                if child_mask & (1 << mutation) and tree_values[mutation] == -1:
                    tree_values[mutation] = branch.first_root
            tree = tuple(tree_values)
            for remainder in self._iter_forest(branch.remaining_mask):
                produced += 1
                yield _merge_parent_vectors(tree, remainder)
        if produced != branch.completion_count:
            raise AssertionError(
                "traceback branch enumeration differs from exact DP count"
            )

    def iter_all_vertex_sets(self) -> Iterator[tuple[int, ...]]:
        for branch in self.branches():
            for parents in self.iter_branch(branch):
                yield self.vertex_set(parents)


@dataclasses.dataclass(frozen=True)
class CandidateScore:
    vertex_set_id: str
    vertices: tuple[int, ...]
    log_likelihood: float
    global_gap_bound: float
    optimizer_status: str


def _validate_quality_groups(
    quality_groups: Sequence[QualityGroup],
    *,
    k: int,
) -> tuple[QualityGroup, ...]:
    if not quality_groups:
        raise ValueError("quality groups must not be empty")
    normalized = []
    for pattern, qualities, count in quality_groups:
        _validate_pattern(pattern, k=k, allow_x=True)
        qualities = tuple(qualities)
        if len(qualities) != k:
            raise ValueError("quality vector length differs from k")
        if any(type(value) is not int for value in qualities):
            raise ValueError("quality values must be built-in integers")
        if type(count) is not int or count < 1:
            raise ValueError("quality-group count must be positive")
        for symbol, quality in zip(pattern, qualities):
            if (symbol == "X") != (quality == -1):
                raise ValueError("X/quality missingness contract violated")
        normalized.append((pattern, qualities, count))
    return tuple(normalized)


def branch_upper_bound(
    quality_groups: Sequence[QualityGroup],
    vertices: Sequence[int],
    *,
    minimum_error_rate: float,
    maximum_error_rate: float,
    mode: str,
) -> tuple[float, str, bool]:
    """Return a float relaxation, not a directed-rounding certificate.

    The formulas are upper bounds over exact real arithmetic.  NumPy/SciPy
    ordinary floating-point evaluation does not guarantee outward rounding,
    so the returned numerical certificate flag is deliberately false.
    """
    emission, counts = _quality_emission_matrix(
        quality_groups,
        vertices,
        minimum_error_rate=minimum_error_rate,
        maximum_error_rate=maximum_error_rate,
    )
    if mode == "rowwise":
        upper = float(np.dot(counts, np.log(np.max(emission, axis=1))))
        if not math.isfinite(upper):
            raise FloatingPointError("rowwise likelihood bound is non-finite")
        return (
            upper,
            "ROWWISE_MAX_EMISSION_FLOAT_RELAXATION_NOT_OUTWARD_ROUNDED",
            False,
        )
    if mode != "union_mixture_float":
        raise ValueError("unknown upper-bound mode")
    fit = fit_quality_aware_mixture(
        quality_groups,
        vertices,
        minimum_error_rate=minimum_error_rate,
        maximum_error_rate=maximum_error_rate,
    )
    if (
        type(fit.converged) is not bool
        or not fit.converged
        or type(fit.monotone) is not bool
        or not fit.monotone
        or not math.isfinite(fit.log_likelihood)
        or not math.isfinite(fit.global_log_likelihood_gap_bound)
        or fit.global_log_likelihood_gap_bound < 0
    ):
        raise RuntimeError("union-mixture float relaxation is not numerically valid")
    upper = float(fit.log_likelihood + fit.global_log_likelihood_gap_bound)
    if not math.isfinite(upper):
        raise FloatingPointError("union-mixture float relaxation is non-finite")
    return (
        upper,
        "UNION_MIXTURE_CONVEX_HULL_FLOAT_RELAXATION_NOT_OUTWARD_ROUNDED",
        False,
    )


def _score_candidate(
    quality_groups: Sequence[QualityGroup],
    vertices: tuple[int, ...],
    *,
    k: int,
    minimum_error_rate: float,
    maximum_error_rate: float,
) -> CandidateScore:
    fit = fit_quality_aware_mixture(
        quality_groups,
        vertices,
        minimum_error_rate=minimum_error_rate,
        maximum_error_rate=maximum_error_rate,
    )
    if (
        type(fit.converged) is not bool
        or not fit.converged
        or type(fit.monotone) is not bool
        or not fit.monotone
        or not math.isfinite(fit.log_likelihood)
        or not math.isfinite(fit.global_log_likelihood_gap_bound)
        or fit.global_log_likelihood_gap_bound < 0
        or fit.global_log_likelihood_gap_bound
        > OPTIMIZER_GLOBAL_LL_GAP_TOLERANCE
    ):
        raise RuntimeError("candidate optimizer result is not finite/KKT-bounded")
    return CandidateScore(
        vertex_set_id=vertex_set_digest(k, vertices),
        vertices=vertices,
        log_likelihood=float(fit.log_likelihood),
        global_gap_bound=float(fit.global_log_likelihood_gap_bound),
        optimizer_status=str(fit.optimizer_status),
    )


def rank_perfect_family_lazy(
    *,
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
    k: int,
    quality_groups: Sequence[QualityGroup],
    structural_alt_universe_mask: int | None = None,
    tie_tolerance: float = 1e-6,
    max_candidates: int = 10_000,
    deadline_seconds: float = 30.0,
    max_tie_ids: int = 10_000,
    bound_mode: str = "union_mixture_float",
    enable_numerical_pruning: bool = False,
    minimum_error_rate: float = 1e-6,
    maximum_error_rate: float = 0.25,
) -> dict[str, Any]:
    """Rank one perfect family and publish only machine-complete enumeration."""
    started = time.perf_counter()
    if (
        type(deadline_seconds) not in (int, float)
        or not math.isfinite(deadline_seconds)
        or deadline_seconds <= 0
        or type(max_candidates) is not int
        or max_candidates < 1
        or type(max_tie_ids) is not int
        or max_tie_ids < 1
        or type(tie_tolerance) not in (int, float)
        or not math.isfinite(tie_tolerance)
        or tie_tolerance < 0
    ):
        raise ValueError("invalid bounded-ranking controls")
    if type(enable_numerical_pruning) is not bool:
        raise ValueError("enable_numerical_pruning must be a bool")
    if type(bound_mode) is not str:
        raise ValueError("bound_mode must be a string")
    if (
        type(minimum_error_rate) not in (int, float)
        or type(maximum_error_rate) not in (int, float)
        or not math.isfinite(minimum_error_rate)
        or not math.isfinite(maximum_error_rate)
        or minimum_error_rate <= 0
        or maximum_error_rate <= minimum_error_rate
        or maximum_error_rate >= 0.5
    ):
        raise ValueError("invalid error-rate bounds")
    quality_groups = _validate_quality_groups(quality_groups, k=k)
    dag = PerfectTracebackDag(
        full_patterns=full_patterns,
        partial_patterns=partial_patterns,
        k=k,
        structural_alt_universe_mask=structural_alt_universe_mask,
    )
    total = dag.forest_count(dag.full_mask)
    branches = []
    for branch in dag.branches():
        if enable_numerical_pruning:
            upper, bound_status, bound_certified = branch_upper_bound(
                quality_groups,
                branch.possible_vertices,
                minimum_error_rate=minimum_error_rate,
                maximum_error_rate=maximum_error_rate,
                mode=bound_mode,
            )
        else:
            upper = math.inf
            bound_status = "NOT_COMPUTED_FULL_ENUMERATION_MODE"
            bound_certified = False
        branches.append(
            (
                upper,
                branch.identity(),
                branch,
                bound_status,
                bound_certified,
            )
        )
    if enable_numerical_pruning:
        branches.sort(key=lambda item: (-item[0], item[1]))
    else:
        branches.sort(key=lambda item: item[1])
    precomputation_elapsed = time.perf_counter() - started

    scores: list[CandidateScore] = []
    evaluated = 0
    pruned = 0
    processed_branches = 0
    stop_status: str | None = None
    incumbent = -math.inf
    seen_ids: set[str] = set()
    bound_epsilon = 1e-12
    numerical_bound_certified = bool(branches) and all(
        certified for _ub, _identity, _branch, _status, certified in branches
    )

    for upper, _identity, branch, _bound_status, _bound_certified in branches:
        if time.perf_counter() - started >= deadline_seconds:
            stop_status = "INCOMPLETE_DEADLINE"
            break
        if (
            enable_numerical_pruning
            and incumbent > -math.inf
            and upper
            < (
            incumbent - tie_tolerance - bound_epsilon
            )
        ):
            pruned += branch.completion_count
            processed_branches += 1
            continue
        branch_evaluated = 0
        for parents in dag.iter_branch(branch):
            if time.perf_counter() - started >= deadline_seconds:
                stop_status = "INCOMPLETE_DEADLINE"
                break
            if evaluated >= max_candidates:
                stop_status = "INCOMPLETE_MAX_CANDIDATES"
                break
            if (
                enable_numerical_pruning
                and incumbent > -math.inf
                and upper
                < (
                incumbent - tie_tolerance - bound_epsilon
                )
            ):
                pruned += branch.completion_count - branch_evaluated
                break
            vertices = dag.vertex_set(parents)
            score = _score_candidate(
                quality_groups,
                vertices,
                k=k,
                minimum_error_rate=minimum_error_rate,
                maximum_error_rate=maximum_error_rate,
            )
            if score.vertex_set_id in seen_ids:
                raise AssertionError("perfect traceback produced duplicate vertex set")
            seen_ids.add(score.vertex_set_id)
            scores.append(score)
            evaluated += 1
            branch_evaluated += 1
            incumbent = max(incumbent, score.log_likelihood)
        if stop_status:
            break
        processed_branches += 1

    traceback_accounting_complete = (
        stop_status is None
        and processed_branches == len(branches)
        and evaluated + pruned == total
    )
    relied_on_numerical_pruning = pruned > 0
    machine_enumeration_complete = (
        traceback_accounting_complete
        and (
            not relied_on_numerical_pruning
            or numerical_bound_certified
        )
    )
    search_complete = machine_enumeration_complete
    if stop_status is not None:
        status = stop_status or "INCOMPLETE_TRACEBACK_ACCOUNTING"
        best_score = None
        winner_ids: list[str] = []
        tie_count = None
        tie_complete = False
    elif not traceback_accounting_complete:
        status = "INCOMPLETE_TRACEBACK_ACCOUNTING"
        best_score = None
        winner_ids = []
        tie_count = None
        tie_complete = False
    elif relied_on_numerical_pruning and not numerical_bound_certified:
        status = "INCOMPLETE_UNCERTIFIED_FLOAT_BOUND_PRUNING"
        best_score = None
        winner_ids = []
        tie_count = None
        tie_complete = False
    elif not scores:
        raise AssertionError("complete search has no evaluated candidate")
    else:
        best_score = max(score.log_likelihood for score in scores)
        ties = sorted(
            score.vertex_set_id
            for score in scores
            if best_score - score.log_likelihood <= tie_tolerance
        )
        tie_count = len(ties)
        tie_complete = tie_count <= max_tie_ids
        winner_ids = ties if tie_complete else []
        status = (
            "EXACT_BEST_AND_COMPLETE_TIE_CLASS"
            if tie_complete
            else "INCOMPLETE_TIE_CLASS_OUTPUT_CAP"
        )

    ranking_complete = search_complete and tie_complete
    diagnostic_incumbent = None
    if scores:
        diagnostic_best = max(score.log_likelihood for score in scores)
        diagnostic_incumbent = {
            "log_likelihood": diagnostic_best,
            "vertex_set_ids": sorted(
                score.vertex_set_id
                for score in scores
                if diagnostic_best - score.log_likelihood <= tie_tolerance
            )[:max_tie_ids],
            "authoritative": False,
        }
    return {
        "schema": SCHEMA,
        "scope": SCOPE,
        "status": status,
        "endpoint_same_as_current_primary_bq_likelihood": True,
        "full_m2_release_endpoint_same": False,
        "effective_m": dag.m,
        "logical_family_count": total,
        "traceback_branch_count": len(branches),
        "traceback_count_conserved": sum(
            branch.completion_count
            for _ub, _identity, branch, _status, _certified in branches
        )
        == total,
        "bound_mode": bound_mode,
        "bound_contract": (
            "ordinary float relaxation; no outward rounding certificate"
        ),
        "numerical_bound_certified": numerical_bound_certified,
        "numerical_pruning_requested": enable_numerical_pruning,
        "numerical_pruning_used": relied_on_numerical_pruning,
        "top_level_branch_bounds_only": True,
        "recursive_branch_and_bound": False,
        "possible_vertices_materialized_for_all_top_level_branches": True,
        "top_level_materialization_worst_case": "O(m * 3^m)",
        "evaluated_candidate_count": evaluated,
        "pruned_candidate_count": pruned,
        "processed_branch_count": processed_branches,
        "traceback_accounting_complete": traceback_accounting_complete,
        "machine_enumeration_complete": machine_enumeration_complete,
        "search_complete": search_complete,
        "tie_class_complete": tie_complete,
        "ranking_complete": ranking_complete,
        "exact_publishable": ranking_complete,
        "best_log_likelihood": best_score if ranking_complete else None,
        "winner_vertex_set_ids": winner_ids if ranking_complete else [],
        "tie_count": tie_count if ranking_complete else None,
        "diagnostic_incumbent": (
            None if ranking_complete else diagnostic_incumbent
        ),
        "candidate_materialization_complete": evaluated == total,
        "machine_exact_by_full_enumeration": (
            ranking_complete and evaluated == total and pruned == 0
        ),
        "logical_family_materialized": False,
        "deadline_seconds": deadline_seconds,
        "deadline_is_hard_wall_bound": False,
        "deadline_scope": (
            "checked in top-level/traceback enumeration loops after DP and "
            "all top-level branch/possible-vertex construction"
        ),
        "precomputation_elapsed_seconds": precomputation_elapsed,
        "max_candidates": max_candidates,
        "max_tie_ids": max_tie_ids,
        "tie_tolerance": tie_tolerance,
        "elapsed_seconds": time.perf_counter() - started,
    }


def brute_force_perfect_vertex_sets(
    *,
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
    k: int,
    structural_alt_universe_mask: int | None = None,
) -> tuple[tuple[int, ...], ...]:
    """Independent small-m exhaustive oracle over all parent vectors."""
    observed = 0
    for pattern in tuple(full_patterns) + tuple(partial_patterns):
        for bit, symbol in enumerate(pattern):
            if symbol == "A":
                observed |= 1 << bit
    universe = observed if structural_alt_universe_mask is None else int(
        structural_alt_universe_mask
    )
    if universe != observed:
        raise ValueError("structural ALT universe differs from observed ALT union")
    active = tuple(bit for bit in range(k) if universe & (1 << bit))
    m = len(active)
    if m > 5:
        raise ValueError("brute-force oracle is gated at m<=5")

    def compatible(pattern: str, state: int) -> bool:
        return all(
            symbol == "X"
            or (symbol == "A") == bool(state & (1 << bit))
            for bit, symbol in enumerate(pattern)
        )

    family: set[tuple[int, ...]] = set()
    choices = tuple(range(-1, m))
    for parents in itertools.product(choices, repeat=m):
        if any(parent == mutation for mutation, parent in enumerate(parents)):
            continue
        memo: dict[int, int] = {}
        valid = True

        def compact_state(mutation: int, stack: set[int]) -> int:
            nonlocal valid
            if mutation in memo:
                return memo[mutation]
            if mutation in stack:
                valid = False
                return 0
            parent = parents[mutation]
            if parent == -1:
                value = 1 << mutation
            elif 0 <= parent < m:
                value = compact_state(parent, stack | {mutation}) | (
                    1 << mutation
                )
            else:
                valid = False
                return 0
            memo[mutation] = value
            return value

        compact_states = [
            compact_state(mutation, set()) for mutation in range(m)
        ]
        if not valid:
            continue
        states = {0}
        for compact in compact_states:
            expanded = 0
            for compact_bit, original_bit in enumerate(active):
                if compact & (1 << compact_bit):
                    expanded |= 1 << original_bit
            states.add(expanded)
        if len(states) != m + 1:
            continue
        if all(
            any(compatible(pattern, state) for state in states)
            for pattern in tuple(full_patterns) + tuple(partial_patterns)
        ):
            family.add(tuple(sorted(states)))
    return tuple(sorted(family))


def exhaustive_current_rank(
    *,
    vertex_sets: Iterable[tuple[int, ...]],
    quality_groups: Sequence[QualityGroup],
    k: int,
    tie_tolerance: float,
    minimum_error_rate: float = 1e-6,
    maximum_error_rate: float = 0.25,
) -> dict[str, Any]:
    """Materialized oracle using the current primary candidate fit."""
    scores = [
        _score_candidate(
            quality_groups,
            tuple(vertices),
            k=k,
            minimum_error_rate=minimum_error_rate,
            maximum_error_rate=maximum_error_rate,
        )
        for vertices in vertex_sets
    ]
    if not scores:
        raise ValueError("oracle family is empty")
    best = max(score.log_likelihood for score in scores)
    winners = sorted(
        score.vertex_set_id
        for score in scores
        if best - score.log_likelihood <= tie_tolerance
    )
    return {
        "best_log_likelihood": best,
        "winner_vertex_set_ids": winners,
        "tie_count": len(winners),
        "candidate_count": len(scores),
    }
