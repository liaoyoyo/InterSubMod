#!/usr/bin/env python3
"""Exact pilot utilities for rooted group-Steiner trees on a Boolean hypercube.

This research module is deliberately independent from the canonical v5 solver.
It provides:

* symbolic R/A/X subcubes that materialize only the active compatible vertex
  indices needed by each reduced sparse MILP row, not completion-wise tree
  worlds or a cross-read Cartesian product;
* a SciPy/HiGHS MILP for the minimum-extra-vertex objective;
* optional enumeration of all optimal vertex sets;
* a mask-conditioned read-pattern mixture likelihood.

The likelihood observes mutation states, not evolutionary transitions. Candidates
with the same vertex set must therefore receive the same score regardless of edges.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp, minimize
from scipy.sparse import coo_matrix, csr_matrix, vstack


def _popcount(value: int) -> int:
    """Python 3.9-compatible population count."""
    return bin(int(value)).count("1")


@dataclass(frozen=True)
class SymbolicPattern:
    """A subcube represented by observed and free bit masks."""

    k: int
    fixed_mask: int
    alt_mask: int
    free_mask: int
    pattern: str
    count: int = 1

    @classmethod
    def from_string(cls, pattern: str, count: int = 1) -> "SymbolicPattern":
        if not pattern:
            raise ValueError("pattern must not be empty")
        if count <= 0:
            raise ValueError("count must be positive")
        invalid = sorted(set(pattern) - {"R", "A", "X"})
        if invalid:
            raise ValueError(f"invalid pattern symbols: {invalid}")
        fixed_mask = 0
        alt_mask = 0
        for bit, symbol in enumerate(pattern):
            if symbol in {"R", "A"}:
                fixed_mask |= 1 << bit
            if symbol == "A":
                alt_mask |= 1 << bit
        all_mask = (1 << len(pattern)) - 1
        return cls(
            k=len(pattern),
            fixed_mask=fixed_mask,
            alt_mask=alt_mask,
            free_mask=all_mask ^ fixed_mask,
            pattern=pattern,
            count=int(count),
        )

    @property
    def n_free(self) -> int:
        return _popcount(self.free_mask)

    def compatible(self, vertex: int, universe_mask: int | None = None) -> bool:
        if vertex < 0 or vertex >= (1 << self.k):
            return False
        if universe_mask is not None and vertex & ~universe_mask:
            return False
        return ((vertex ^ self.alt_mask) & self.fixed_mask) == 0

    def n_completions(self, universe_mask: int | None = None) -> int:
        free = self.free_mask if universe_mask is None else self.free_mask & universe_mask
        return 1 << _popcount(free)

    def enumerate_completions(self, universe_mask: int | None = None) -> tuple[int, ...]:
        """Testing/reference helper; production constraints use ``compatible``."""
        active = (1 << self.k) - 1 if universe_mask is None else universe_mask
        return tuple(v for v in submasks(active) if self.compatible(v, active))


def parse_full(pattern: str) -> int:
    if not pattern or set(pattern) - {"R", "A"}:
        raise ValueError(f"full pattern must contain only R/A: {pattern!r}")
    state = 0
    for bit, symbol in enumerate(pattern):
        if symbol == "A":
            state |= 1 << bit
    return state


def state_to_pattern(state: int, k: int) -> str:
    return "".join("A" if state & (1 << bit) else "R" for bit in range(k))


def submasks(mask: int) -> tuple[int, ...]:
    values = []
    current = mask
    while True:
        values.append(current)
        if current == 0:
            break
        current = (current - 1) & mask
    return tuple(sorted(values))


def observed_alt_universe(full_patterns: Iterable[str], partial_patterns: Iterable[str]) -> int:
    mask = 0
    for pattern in list(full_patterns) + list(partial_patterns):
        for bit, symbol in enumerate(pattern):
            if symbol == "A":
                mask |= 1 << bit
    return mask


def vertex_set_digest(k: int, vertices: Iterable[int]) -> str:
    payload = {"k": int(k), "vertices": sorted(int(v) for v in vertices)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@dataclass
class MilpSolution:
    status: str
    status_code: int
    message: str
    objective: float | None
    objective_bound: float | None
    mip_gap: float | None
    mip_node_count: int | None
    selected_vertices: tuple[int, ...]
    vertex_set_id: str | None
    n_variables: int
    n_constraints: int
    universe_mask: int
    universe_mode: str
    n_partial_groups_input: int
    n_partial_groups_active: int
    n_partial_groups_duplicate_removed: int
    n_partial_groups_dominated_removed: int
    n_partial_groups_required_hit_removed: int
    n_partial_groups_forced_removed: int
    n_forced_vertices: int
    n_no_good_nonzeros: int

    def to_dict(self) -> dict:
        return asdict(self)


def _status_name(status_code: int, success: bool) -> str:
    if success and status_code == 0:
        return "OPTIMAL"
    return {
        1: "LIMIT_REACHED",
        2: "INFEASIBLE",
        3: "UNBOUNDED",
        4: "NUMERICAL_OR_OTHER",
    }.get(int(status_code), "UNKNOWN")


@dataclass(frozen=True)
class _PartialGroupReduction:
    """Exact presolve bookkeeping for structural partial-read constraints."""

    groups: tuple[SymbolicPattern, ...]
    forced_vertices: tuple[int, ...]
    n_input: int
    n_duplicate_removed: int
    n_dominated_removed: int
    n_required_hit_removed: int
    n_forced_removed: int


def _effective_group_signature(group: SymbolicPattern, universe_mask: int) -> tuple[int, int]:
    """Return the fixed/value masks after inactive dimensions collapse to REF."""
    if group.alt_mask & ~universe_mask:
        raise ValueError(f"partial group has ALT outside solver universe: {group.pattern}")
    fixed_mask = group.fixed_mask & universe_mask
    return fixed_mask, group.alt_mask & fixed_mask


def _signature_contains_vertex(signature: tuple[int, int], vertex: int) -> bool:
    fixed_mask, alt_mask = signature
    return ((vertex ^ alt_mask) & fixed_mask) == 0


def _signature_is_subset(
    left: tuple[int, int],
    right: tuple[int, int],
) -> bool:
    """Return whether subcube ``left`` is contained in subcube ``right``."""
    left_fixed, left_alt = left
    right_fixed, right_alt = right
    return (right_fixed & ~left_fixed) == 0 and ((left_alt ^ right_alt) & right_fixed) == 0


def _reduce_partial_groups(
    partial_patterns: Sequence[SymbolicPattern],
    universe_mask: int,
    mandatory: set[int],
) -> _PartialGroupReduction:
    """Remove constraints that are algebraically redundant without changing feasibility.

    Counts never enter the structural objective, so identical effective subcubes are
    duplicate constraints.  If G1 is contained in G2, satisfying G1 also satisfies
    G2.  A singleton subcube is represented exactly by fixing its only vertex to one;
    it remains objective-bearing unless it was already a full-read mandatory state.
    """
    representatives: dict[tuple[int, int], SymbolicPattern] = {}
    for group in partial_patterns:
        signature = _effective_group_signature(group, universe_mask)
        representatives.setdefault(signature, group)
    n_duplicate_removed = len(partial_patterns) - len(representatives)

    remaining: list[tuple[tuple[int, int], SymbolicPattern]] = []
    forced_vertices: set[int] = set()
    n_required_hit_removed = 0
    n_forced_removed = 0
    for signature, group in representatives.items():
        if any(_signature_contains_vertex(signature, vertex) for vertex in mandatory):
            n_required_hit_removed += 1
            continue
        fixed_mask, alt_mask = signature
        if (universe_mask & ~fixed_mask) == 0:
            forced_vertices.add(alt_mask)
            n_forced_removed += 1
            continue
        remaining.append((signature, group))

    if forced_vertices:
        still_uncovered = []
        for signature, group in remaining:
            if any(_signature_contains_vertex(signature, vertex) for vertex in forced_vertices):
                n_required_hit_removed += 1
            else:
                still_uncovered.append((signature, group))
        remaining = still_uncovered

    # A strict subset has fewer completions, so ascending cardinality ensures a
    # possible dominator is considered before the subcube it makes redundant.
    remaining.sort(
        key=lambda item: (
            1 << _popcount(universe_mask & ~item[0][0]),
            item[0][0],
            item[0][1],
        )
    )
    kept: list[tuple[tuple[int, int], SymbolicPattern]] = []
    n_dominated_removed = 0
    for signature, group in remaining:
        if any(_signature_is_subset(kept_signature, signature) for kept_signature, _ in kept):
            n_dominated_removed += 1
            continue
        kept.append((signature, group))

    return _PartialGroupReduction(
        groups=tuple(group for _, group in kept),
        forced_vertices=tuple(sorted(forced_vertices)),
        n_input=len(partial_patterns),
        n_duplicate_removed=n_duplicate_removed,
        n_dominated_removed=n_dominated_removed,
        n_required_hit_removed=n_required_hit_removed,
        n_forced_removed=n_forced_removed,
    )


def _compatible_vertices(group: SymbolicPattern, universe_mask: int) -> tuple[int, ...]:
    """Materialize only one group's symbolic subcube, not the entire universe scan."""
    fixed_mask, alt_mask = _effective_group_signature(group, universe_mask)
    free_mask = universe_mask & ~fixed_mask
    return tuple(sorted(alt_mask | free_state for free_state in submasks(free_mask)))


@dataclass(frozen=True)
class _PreparedSelectionProblem:
    """Fixed MILP payload shared by every solve in one exact enumeration."""

    k: int
    vertices: tuple[int, ...]
    index: dict[int, int]
    vertex_universe: frozenset[int]
    mandatory: frozenset[int]
    objective: np.ndarray
    integrality: np.ndarray
    bounds: Bounds
    base_matrix: csr_matrix
    base_lb: np.ndarray
    base_ub: np.ndarray
    universe_mask: int
    universe_mode: str
    reduction: _PartialGroupReduction


@dataclass
class _SelectionConstraintState:
    """Per-enumeration constraint rows appended to an immutable prepared base."""

    matrix: csr_matrix
    lb: np.ndarray
    ub: np.ndarray
    objective_value: int | None
    n_no_good_nonzeros: int = 0


def _sparse_row(row: dict[int, float], n_columns: int) -> csr_matrix:
    """Build one sparse row, preserving an explicit all-zero constraint row."""
    columns = [column for column, value in row.items() if value]
    data = [float(row[column]) for column in columns]
    return coo_matrix(
        (data, ([0] * len(columns), columns)),
        shape=(1, n_columns),
        dtype=float,
    ).tocsr()


def _append_constraint_row(
    state: _SelectionConstraintState,
    row: dict[int, float],
    lower_bound: float,
    upper_bound: float,
) -> None:
    state.matrix = vstack(
        (state.matrix, _sparse_row(row, state.matrix.shape[1])),
        format="csr",
    )
    state.lb = np.append(state.lb, float(lower_bound))
    state.ub = np.append(state.ub, float(upper_bound))


def _build_problem(
    full_patterns: Sequence[str],
    partial_patterns: Sequence[SymbolicPattern],
    k: int,
    universe_mode: str,
) -> _PreparedSelectionProblem:
    if universe_mode not in {"observed_alt", "all_loci"}:
        raise ValueError("universe_mode must be observed_alt or all_loci")
    for pattern in full_patterns:
        if len(pattern) != k:
            raise ValueError("full pattern length differs from k")
    for group in partial_patterns:
        if group.k != k:
            raise ValueError("partial pattern length differs from k")

    universe_mask = (
        observed_alt_universe(full_patterns, (g.pattern for g in partial_patterns))
        if universe_mode == "observed_alt"
        else (1 << k) - 1
    )
    vertices = submasks(universe_mask)
    index = {vertex: idx for idx, vertex in enumerate(vertices)}
    vertex_universe = frozenset(vertices)
    mandatory = frozenset({0, *(parse_full(pattern) for pattern in full_patterns)})
    if any(vertex not in index for vertex in mandatory):
        raise AssertionError("mandatory state not represented in solver universe")

    reduction = _reduce_partial_groups(partial_patterns, universe_mask, mandatory)
    forced_vertices = set(reduction.forced_vertices)
    required_vertices = set(mandatory) | forced_vertices

    n = len(vertices)
    lower_bounds = np.zeros(n, dtype=float)
    upper_bounds = np.ones(n, dtype=float)
    for vertex in required_vertices:
        lower_bounds[index[vertex]] = 1.0
        upper_bounds[index[vertex]] = 1.0

    objective = np.array([0.0 if vertex in mandatory else 1.0 for vertex in vertices], dtype=float)
    rows: list[dict[int, float]] = []
    row_lb: list[float] = []
    row_ub: list[float] = []

    # Every selected non-root state must have a selected Hamming-1 predecessor.
    # The DAG strictly decreases in popcount, so this local condition implies
    # root connectivity without a separate multi-commodity flow.
    active_bits = tuple(bit for bit in range(k) if universe_mask & (1 << bit))
    for vertex in vertices:
        if vertex == 0:
            continue
        row = {index[vertex]: 1.0}
        for bit in active_bits:
            if vertex & (1 << bit):
                pred = vertex ^ (1 << bit)
                if pred in index:
                    row[index[pred]] = row.get(index[pred], 0.0) - 1.0
        rows.append(row)
        row_lb.append(-np.inf)
        row_ub.append(0.0)

    # Each partial read is a group terminal: select at least one compatible state.
    for group in reduction.groups:
        compatible_indices = [index[v] for v in _compatible_vertices(group, universe_mask)]
        if not compatible_indices:
            raise ValueError(f"partial group has no completion in universe: {group.pattern}")
        rows.append({idx: 1.0 for idx in compatible_indices})
        row_lb.append(1.0)
        row_ub.append(np.inf)

    data: list[float] = []
    row_idx: list[int] = []
    col_idx: list[int] = []
    for ridx, row in enumerate(rows):
        for cidx, value in row.items():
            if value:
                row_idx.append(ridx)
                col_idx.append(cidx)
                data.append(value)
    matrix = coo_matrix((data, (row_idx, col_idx)), shape=(len(rows), n)).tocsr()
    return _PreparedSelectionProblem(
        k=k,
        vertices=vertices,
        index=index,
        vertex_universe=vertex_universe,
        mandatory=mandatory,
        objective=objective,
        integrality=np.ones(len(vertices), dtype=np.uint8),
        bounds=Bounds(lower_bounds, upper_bounds),
        base_matrix=matrix,
        base_lb=np.asarray(row_lb, dtype=float),
        base_ub=np.asarray(row_ub, dtype=float),
        universe_mask=universe_mask,
        universe_mode=universe_mode,
        reduction=reduction,
    )


def _new_constraint_state(
    problem: _PreparedSelectionProblem,
    objective_value: int | None = None,
) -> _SelectionConstraintState:
    state = _SelectionConstraintState(
        matrix=problem.base_matrix,
        lb=problem.base_lb,
        ub=problem.base_ub,
        objective_value=objective_value,
    )
    if objective_value is not None:
        _append_constraint_row(
            state,
            {
                idx: value
                for idx, value in enumerate(problem.objective)
                if value
            },
            float(objective_value),
            float(objective_value),
        )
    return state


def _append_no_good(
    problem: _PreparedSelectionProblem,
    state: _SelectionConstraintState,
    selected: Sequence[int],
) -> None:
    """Append one exact exclusion in discovery order without rebuilding base rows."""
    selected_set = set(selected)
    if not selected_set <= problem.vertex_universe:
        raise ValueError("excluded vertex set contains a state outside the solver universe")
    selected_extra = selected_set - problem.mandatory
    if (
        state.objective_value is not None
        and problem.mandatory <= selected_set
        and len(selected_extra) == int(state.objective_value)
    ):
        # With sum(extra x)=h fixed, sum(x over the prior h extra
        # vertices)<=h-1 excludes exactly that solution.  The row can be empty
        # at h=0 and must still be appended because its upper bound is -1.
        row = {problem.index[vertex]: 1.0 for vertex in selected_extra}
        lower_bound = -np.inf
        upper_bound = float(state.objective_value) - 1.0
    else:
        # General exact exclusion when a cardinality equality is unavailable.
        row = {
            idx: (-1.0 if vertex in selected_set else 1.0)
            for idx, vertex in enumerate(problem.vertices)
        }
        lower_bound = 1.0 - len(selected_set)
        upper_bound = np.inf
    _append_constraint_row(state, row, lower_bound, upper_bound)
    state.n_no_good_nonzeros += len(row)


def _solve_selection(
    problem: _PreparedSelectionProblem,
    state: _SelectionConstraintState,
    *,
    time_limit_seconds: float,
) -> MilpSolution:
    constraints = LinearConstraint(state.matrix, state.lb, state.ub)
    result = milp(
        c=problem.objective,
        integrality=problem.integrality,
        bounds=problem.bounds,
        constraints=constraints,
        options={"time_limit": float(time_limit_seconds), "mip_rel_gap": 0.0, "presolve": True},
    )
    selected = tuple(
        problem.vertices[idx]
        for idx, value in enumerate(result.x if result.x is not None else [])
        if value > 0.5
    )
    status = _status_name(int(result.status), bool(result.success))
    objective_result = None if result.fun is None else float(result.fun)
    reduction = problem.reduction
    return MilpSolution(
        status=status,
        status_code=int(result.status),
        message=str(result.message),
        objective=objective_result,
        objective_bound=(None if getattr(result, "mip_dual_bound", None) is None else float(result.mip_dual_bound)),
        mip_gap=(None if getattr(result, "mip_gap", None) is None else float(result.mip_gap)),
        mip_node_count=(None if getattr(result, "mip_node_count", None) is None else int(result.mip_node_count)),
        selected_vertices=selected,
        vertex_set_id=(vertex_set_digest(problem.k, selected) if selected else None),
        n_variables=len(problem.vertices),
        n_constraints=int(state.matrix.shape[0]),
        universe_mask=problem.universe_mask,
        universe_mode=problem.universe_mode,
        n_partial_groups_input=reduction.n_input,
        n_partial_groups_active=len(reduction.groups),
        n_partial_groups_duplicate_removed=reduction.n_duplicate_removed,
        n_partial_groups_dominated_removed=reduction.n_dominated_removed,
        n_partial_groups_required_hit_removed=reduction.n_required_hit_removed,
        n_partial_groups_forced_removed=reduction.n_forced_removed,
        n_forced_vertices=len(reduction.forced_vertices),
        n_no_good_nonzeros=state.n_no_good_nonzeros,
    )


def solve_min_hidden(
    full_patterns: Sequence[str],
    partial_pattern_strings: Sequence[str],
    k: int,
    *,
    universe_mode: str = "observed_alt",
    time_limit_seconds: float = 30.0,
    objective_value: int | None = None,
    excluded_vertex_sets: Sequence[Sequence[int]] = (),
) -> MilpSolution:
    groups = tuple(SymbolicPattern.from_string(pattern) for pattern in partial_pattern_strings)
    problem = _build_problem(
        tuple(full_patterns),
        groups,
        k,
        universe_mode,
    )
    state = _new_constraint_state(problem, objective_value=objective_value)
    for selected in excluded_vertex_sets:
        _append_no_good(problem, state, selected)
    return _solve_selection(
        problem,
        state,
        time_limit_seconds=time_limit_seconds,
    )


def enumerate_optimal_vertex_sets(
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
    k: int,
    *,
    universe_mode: str = "observed_alt",
    max_sets: int = 256,
    time_limit_seconds: float = 30.0,
) -> dict:
    groups = tuple(SymbolicPattern.from_string(pattern) for pattern in partial_patterns)
    problem = _build_problem(
        tuple(full_patterns),
        groups,
        k,
        universe_mode,
    )
    first = _solve_selection(
        problem,
        _new_constraint_state(problem),
        time_limit_seconds=time_limit_seconds,
    )
    if first.status != "OPTIMAL" or first.objective is None:
        return {"first": first.to_dict(), "complete": False, "vertex_sets": []}
    optimum = int(round(first.objective))
    sets = [first.selected_vertices]
    if optimum == 0:
        return {
            "first": first.to_dict(),
            "objective": optimum,
            "complete": True,
            "vertex_sets": [tuple(v) for v in sets],
            "vertex_set_ids": [vertex_set_digest(k, v) for v in sets],
        }
    fixed_state = _new_constraint_state(problem, objective_value=optimum)
    _append_no_good(problem, fixed_state, first.selected_vertices)
    while len(sets) < max_sets:
        next_solution = _solve_selection(
            problem,
            fixed_state,
            time_limit_seconds=time_limit_seconds,
        )
        if next_solution.status == "INFEASIBLE":
            return {
                "first": first.to_dict(),
                "objective": optimum,
                "complete": True,
                "vertex_sets": [tuple(v) for v in sets],
                "vertex_set_ids": [vertex_set_digest(k, v) for v in sets],
            }
        if next_solution.status != "OPTIMAL":
            return {
                "first": first.to_dict(),
                "objective": optimum,
                "complete": False,
                "stop_status": next_solution.to_dict(),
                "vertex_sets": [tuple(v) for v in sets],
                "vertex_set_ids": [vertex_set_digest(k, v) for v in sets],
            }
        sets.append(next_solution.selected_vertices)
        if len(sets) < max_sets:
            _append_no_good(problem, fixed_state, next_solution.selected_vertices)
    return {
        "first": first.to_dict(),
        "objective": optimum,
        "complete": False,
        "stop_status": "MAX_SETS_REACHED",
        "vertex_sets": [tuple(v) for v in sets],
        "vertex_set_ids": [vertex_set_digest(k, v) for v in sets],
    }


def parent_choice_count(vertices: Iterable[int]) -> int:
    selected = set(vertices)
    count = 1
    for vertex in selected:
        if vertex == 0:
            continue
        n_parents = sum(
            1
            for bit in range(vertex.bit_length())
            if vertex & (1 << bit) and (vertex ^ (1 << bit)) in selected
        )
        if n_parents == 0:
            return 0
        count *= n_parents
    return count


@dataclass
class MixtureFit:
    vertices: tuple[int, ...]
    weights: tuple[float, ...]
    log_likelihood: float
    converged: bool
    iterations: int
    emission_rank: int
    monotone: bool
    optimizer_status: str = "UNSPECIFIED"
    slsqp_success: bool | None = None
    slsqp_status: int | None = None
    slsqp_message: str | None = None
    warm_start_log_likelihood: float = math.nan
    warm_start_global_log_likelihood_gap_bound: float = math.inf
    refinement_iterations: int = 0
    global_log_likelihood_gap_bound: float = math.inf
    simplex_residual: float = math.inf
    augmented_emission_rank: int = 0
    mixture_weights_identifiable: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


def mixture_kkt_certificate(
    emission: np.ndarray,
    counts: np.ndarray,
    weights: np.ndarray,
) -> dict[str, float | int | bool]:
    """Return a global first-order certificate for a simplex mixture MLE.

    For ``ell(pi)=sum_i n_i log((Q pi)_i)``, concavity gives

    ``ell(pi*) - ell(pi) <= max_j grad_j ell(pi) - pi' grad ell(pi)``.

    The right-hand side is the Frank-Wolfe/KKT gap.  It is zero exactly at a
    KKT point (up to floating-point tolerance), and therefore certifies a
    *global* optimum without treating an optimizer's status flag as proof.
    Under an exactly normalized simplex, ``pi' grad ell(pi)=sum_i n_i``.
    """
    emission = np.asarray(emission, dtype=float)
    counts = np.asarray(counts, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if emission.ndim != 2 or counts.shape != (emission.shape[0],):
        raise ValueError("emission/count shape mismatch")
    if weights.shape != (emission.shape[1],):
        raise ValueError("emission/weight shape mismatch")
    if np.any(emission < 0) or not np.all(np.isfinite(emission)):
        raise ValueError("emission matrix must be finite and nonnegative")
    if np.any(counts <= 0) or not np.all(np.isfinite(counts)):
        raise ValueError("counts must be finite and positive")
    denominator = emission @ weights
    if np.any(denominator <= 0) or not np.all(np.isfinite(denominator)):
        raise FloatingPointError("non-positive or non-finite likelihood denominator")
    gradient = emission.T @ (counts / denominator)
    weighted_gradient = float(np.dot(weights, gradient))
    gap = max(0.0, float(np.max(gradient) - weighted_gradient))
    simplex_residual = max(
        abs(float(weights.sum()) - 1.0),
        max(0.0, -float(np.min(weights))),
    )
    augmented_rank = int(
        np.linalg.matrix_rank(
            np.vstack((emission, np.ones((1, emission.shape[1]), dtype=float)))
        )
    )
    return {
        "global_log_likelihood_gap_bound": gap,
        "simplex_residual": simplex_residual,
        "emission_rank": int(np.linalg.matrix_rank(emission)),
        "augmented_emission_rank": augmented_rank,
        "mixture_weights_identifiable": augmented_rank == emission.shape[1],
    }


def _pairwise_simplex_refinement(
    emission: np.ndarray,
    counts: np.ndarray,
    initial: np.ndarray,
    *,
    certificate_tolerance: float,
    max_iterations: int,
) -> tuple[np.ndarray, int, bool]:
    """Monotonically transfer mass between two states until KKT-certified.

    At each iteration, mass is moved from the active state with the smallest
    likelihood gradient to the state with the largest gradient.  An exact
    one-dimensional concave line search is solved by bisection.  Unlike plain
    EM, this can activate a state whose current weight is zero and tends to
    remove boundary-state false failures quickly.
    """
    weights = np.asarray(initial, dtype=float).copy()
    weights = np.clip(weights, 0.0, 1.0)
    total = float(weights.sum())
    if total <= 0 or not math.isfinite(total):
        weights = np.full(emission.shape[1], 1.0 / emission.shape[1], dtype=float)
    else:
        weights /= total
    monotone = True
    previous_ll = float(np.dot(counts, np.log(emission @ weights)))
    for iteration in range(max_iterations + 1):
        certificate = mixture_kkt_certificate(emission, counts, weights)
        if certificate["global_log_likelihood_gap_bound"] <= certificate_tolerance:
            return weights, iteration, monotone
        if iteration == max_iterations:
            return weights, iteration, monotone

        denominator = emission @ weights
        gradient = emission.T @ (counts / denominator)
        destination = int(np.argmax(gradient))
        active = np.flatnonzero(weights > 0.0)
        if not len(active):
            return weights, iteration, False
        source = int(active[np.argmin(gradient[active])])
        if source == destination:
            return weights, iteration, False

        direction_emission = emission[:, destination] - emission[:, source]
        upper = float(weights[source])
        derivative_upper = float(
            np.dot(counts, direction_emission / (denominator + upper * direction_emission))
        )
        if derivative_upper >= 0.0:
            step = upper
        else:
            lower = 0.0
            # The derivative is monotone decreasing along this concave line.
            for _ in range(80):
                midpoint = (lower + upper) / 2.0
                derivative = float(
                    np.dot(
                        counts,
                        direction_emission / (denominator + midpoint * direction_emission),
                    )
                )
                if derivative > 0.0:
                    lower = midpoint
                else:
                    upper = midpoint
            step = (lower + upper) / 2.0

        updated = weights.copy()
        updated[source] -= step
        updated[destination] += step
        updated[np.abs(updated) < 1e-18] = 0.0
        updated /= float(updated.sum())
        current_ll = float(np.dot(counts, np.log(emission @ updated)))
        if current_ll < previous_ll - 1e-10 * max(1.0, abs(previous_ll)):
            monotone = False
            return weights, iteration, monotone
        weights = updated
        previous_ll = max(previous_ll, current_ll)
    raise AssertionError("unreachable pairwise refinement state")


def fit_emission_mixture_certified(
    emission: np.ndarray,
    counts: np.ndarray,
    vertices: Sequence[int],
    *,
    tolerance: float = 1e-10,
    max_iterations: int = 2_000,
    certificate_tolerance: float = 1e-8,
    refinement_max_iterations: int = 10_000,
) -> MixtureFit:
    """Fit a concave simplex mixture and certify its global LL error bound.

    SLSQP is only a warm start.  Its status flag is recorded but is not the
    convergence criterion.  A deterministic, monotone pairwise refinement
    then supplies a simplex/KKT (Frank-Wolfe) gap certificate.  This preserves
    the likelihood objective and does not add VAF or any second use of reads.
    """
    emission = np.asarray(emission, dtype=float)
    counts = np.asarray(counts, dtype=float)
    vertices = tuple(int(vertex) for vertex in vertices)
    if emission.ndim != 2 or emission.shape[1] != len(vertices):
        raise ValueError("emission/vertex shape mismatch")
    if counts.shape != (emission.shape[0],):
        raise ValueError("emission/count shape mismatch")
    if not vertices or not len(counts):
        raise ValueError("vertices and counts must not be empty")
    if np.any(emission <= 0) or not np.all(np.isfinite(emission)):
        raise ValueError("certified mixture requires finite positive emissions")
    if np.any(counts <= 0) or not np.all(np.isfinite(counts)):
        raise ValueError("counts must be finite and positive")
    if tolerance <= 0 or certificate_tolerance <= 0:
        raise ValueError("optimizer tolerances must be positive")

    initial = np.full(len(vertices), 1.0 / len(vertices), dtype=float)

    def objective(weights: np.ndarray) -> float:
        denominator = emission @ weights
        if np.any(denominator <= 0) or not np.all(np.isfinite(denominator)):
            return np.inf
        return -float(np.dot(counts, np.log(denominator)))

    def gradient(weights: np.ndarray) -> np.ndarray:
        denominator = emission @ weights
        if np.any(denominator <= 0) or not np.all(np.isfinite(denominator)):
            raise FloatingPointError("non-positive or non-finite likelihood denominator")
        return -(emission.T @ (counts / denominator))

    initial_objective = objective(initial)
    result = minimize(
        objective,
        initial,
        method="SLSQP",
        jac=gradient,
        bounds=[(0.0, 1.0)] * len(vertices),
        constraints={
            "type": "eq",
            "fun": lambda value: float(value.sum() - 1.0),
            "jac": lambda value: np.ones_like(value),
        },
        options={"ftol": tolerance, "maxiter": int(max_iterations), "disp": False},
    )
    candidate = np.asarray(getattr(result, "x", initial), dtype=float)
    if candidate.shape != initial.shape or not np.all(np.isfinite(candidate)):
        candidate = initial.copy()
    candidate = np.clip(candidate, 0.0, 1.0)
    candidate_total = float(candidate.sum())
    candidate = candidate / candidate_total if candidate_total > 0 else initial.copy()
    pre_certificate = mixture_kkt_certificate(emission, counts, candidate)

    refined, refinement_iterations, refinement_monotone = _pairwise_simplex_refinement(
        emission,
        counts,
        candidate,
        certificate_tolerance=certificate_tolerance,
        max_iterations=refinement_max_iterations,
    )
    final_objective = objective(refined)
    certificate = mixture_kkt_certificate(emission, counts, refined)
    monotone = (
        refinement_monotone
        and final_objective <= initial_objective + 1e-8
        and final_objective <= objective(candidate) + 1e-8
    )
    certified = (
        certificate["global_log_likelihood_gap_bound"] <= certificate_tolerance
        and certificate["simplex_residual"] <= 1e-12
        and monotone
    )
    if certified and pre_certificate["global_log_likelihood_gap_bound"] <= certificate_tolerance:
        optimizer_status = "CERTIFIED_SLSQP_WARM_START"
    elif certified:
        optimizer_status = "CERTIFIED_PAIRWISE_REFINEMENT"
    else:
        optimizer_status = "ABSTAIN_KKT_GAP_NOT_CERTIFIED"
    return MixtureFit(
        vertices=vertices,
        weights=tuple(float(value) for value in refined),
        log_likelihood=-final_objective,
        converged=certified,
        iterations=int(getattr(result, "nit", 0)) + refinement_iterations,
        emission_rank=int(certificate["emission_rank"]),
        monotone=monotone,
        optimizer_status=optimizer_status,
        slsqp_success=bool(getattr(result, "success", False)),
        slsqp_status=int(getattr(result, "status", -1)),
        slsqp_message=str(getattr(result, "message", "")),
        warm_start_log_likelihood=-objective(candidate),
        warm_start_global_log_likelihood_gap_bound=float(
            pre_certificate["global_log_likelihood_gap_bound"]
        ),
        refinement_iterations=refinement_iterations,
        global_log_likelihood_gap_bound=float(
            certificate["global_log_likelihood_gap_bound"]
        ),
        simplex_residual=float(certificate["simplex_residual"]),
        augmented_emission_rank=int(certificate["augmented_emission_rank"]),
        mixture_weights_identifiable=bool(certificate["mixture_weights_identifiable"]),
    )


def emission_matrix(
    patterns: Sequence[SymbolicPattern],
    vertices: Sequence[int],
    error_rate: float,
) -> np.ndarray:
    if not 0.0 < error_rate < 0.5:
        raise ValueError("error_rate must be between 0 and 0.5")
    matrix = np.empty((len(patterns), len(vertices)), dtype=float)
    for row, pattern in enumerate(patterns):
        n_fixed = _popcount(pattern.fixed_mask)
        for col, vertex in enumerate(vertices):
            mismatches = _popcount((vertex ^ pattern.alt_mask) & pattern.fixed_mask)
            matches = n_fixed - mismatches
            matrix[row, col] = ((1.0 - error_rate) ** matches) * (error_rate ** mismatches)
    return matrix


def fit_vertex_mixture(
    pattern_counts: Sequence[tuple[str, int]],
    vertices: Sequence[int],
    *,
    error_rate: float = 0.01,
    tolerance: float = 1e-11,
    max_iterations: int = 20_000,
) -> MixtureFit:
    if not vertices:
        raise ValueError("vertices must not be empty")
    all_patterns = tuple(SymbolicPattern.from_string(pattern, count) for pattern, count in pattern_counts)
    if not all_patterns:
        raise ValueError("pattern_counts must not be empty")
    k_values = {pattern.k for pattern in all_patterns}
    if len(k_values) != 1:
        raise ValueError("all patterns must have the same length")
    if any(vertex >= (1 << next(iter(k_values))) for vertex in vertices):
        raise ValueError("vertex exceeds pattern dimension")

    # An all-X observation has P(Y|v)=1 for every state and contributes
    # log(sum(pi))=0. Removing it is exact and avoids tolerance-dependent EM
    # self-reinforcement from an uninformative latent assignment.
    patterns = tuple(pattern for pattern in all_patterns if pattern.fixed_mask != 0)
    if not patterns:
        uniform = tuple(1.0 / len(vertices) for _ in vertices)
        return MixtureFit(
            vertices=tuple(int(v) for v in vertices),
            weights=uniform,
            log_likelihood=0.0,
            converged=True,
            iterations=0,
            emission_rank=1,
            monotone=True,
            optimizer_status="UNINFORMATIVE_ALL_X",
            global_log_likelihood_gap_bound=0.0,
            simplex_residual=0.0,
            augmented_emission_rank=1,
            mixture_weights_identifiable=len(vertices) == 1,
        )

    q = emission_matrix(patterns, vertices, error_rate)
    counts = np.asarray([pattern.count for pattern in patterns], dtype=float)
    pi = np.full(len(vertices), 1.0 / len(vertices), dtype=float)
    previous = -np.inf
    monotone = True
    converged = False
    ll = -np.inf
    for iteration in range(1, max_iterations + 1):
        weighted = q * pi[np.newaxis, :]
        denominator = weighted.sum(axis=1)
        if np.any(denominator <= 0) or not np.all(np.isfinite(denominator)):
            raise FloatingPointError("non-positive or non-finite likelihood denominator")
        responsibility = weighted / denominator[:, np.newaxis]
        expected = (counts[:, np.newaxis] * responsibility).sum(axis=0)
        pi_new = expected / counts.sum()
        denominator_new = (q * pi_new[np.newaxis, :]).sum(axis=1)
        ll = float(np.dot(counts, np.log(denominator_new)))
        if ll < previous - 1e-9:
            monotone = False
        if math.isfinite(previous) and abs(ll - previous) <= tolerance:
            pi = pi_new
            converged = True
            break
        pi = pi_new
        previous = ll
    certificate = mixture_kkt_certificate(q, counts, pi)
    return MixtureFit(
        vertices=tuple(int(v) for v in vertices),
        weights=tuple(float(v) for v in pi),
        log_likelihood=ll,
        converged=converged,
        iterations=iteration,
        emission_rank=int(np.linalg.matrix_rank(q)),
        monotone=monotone,
        optimizer_status="EM_LEGACY_LL_INCREMENT_STOP" if converged else "EM_LEGACY_MAX_ITERATIONS",
        global_log_likelihood_gap_bound=float(
            certificate["global_log_likelihood_gap_bound"]
        ),
        simplex_residual=float(certificate["simplex_residual"]),
        augmented_emission_rank=int(certificate["augmented_emission_rank"]),
        mixture_weights_identifiable=bool(certificate["mixture_weights_identifiable"]),
    )


def fit_vertex_mixture_slsqp(
    pattern_counts: Sequence[tuple[str, int]],
    vertices: Sequence[int],
    *,
    error_rate: float = 0.01,
    tolerance: float = 1e-10,
    max_iterations: int = 2_000,
) -> MixtureFit:
    """Fit the same likelihood with SLSQP plus a global KKT certificate."""
    if not vertices:
        raise ValueError("vertices must not be empty")
    all_patterns = tuple(SymbolicPattern.from_string(pattern, count) for pattern, count in pattern_counts)
    if not all_patterns:
        raise ValueError("pattern_counts must not be empty")
    k_values = {pattern.k for pattern in all_patterns}
    if len(k_values) != 1:
        raise ValueError("all patterns must have the same length")
    k = next(iter(k_values))
    vertices = tuple(sorted(int(vertex) for vertex in vertices))
    if any(vertex < 0 or vertex >= (1 << k) for vertex in vertices):
        raise ValueError("vertex exceeds pattern dimension")
    patterns = tuple(pattern for pattern in all_patterns if pattern.fixed_mask != 0)
    if not patterns:
        uniform = tuple(1.0 / len(vertices) for _ in vertices)
        return MixtureFit(
            vertices,
            uniform,
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

    q = emission_matrix(patterns, vertices, error_rate)
    counts = np.asarray([pattern.count for pattern in patterns], dtype=float)
    return fit_emission_mixture_certified(
        q,
        counts,
        vertices,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
