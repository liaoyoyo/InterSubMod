#!/usr/bin/env python3
"""Exact-preserving optimized backend for the bounded 2026-07-18 probe.

This module remains isolated from the canonical production router.  It adds:

* integer-bitset structural domains and dynamic obligation antichains;
* a one-traversal exact branch-and-bound for the full optimal vertex-set family;
* a q-parameterized directed-Steiner subset DP for an objective-only certificate;
* fixed-node-set additive parent mapping without parent-tree materialization;
* an explicitly projection-only representation of unsupported unary chains; and
* fail-closed evolutionary-model gates.

The structural estimand is unchanged: minimize the number of selected mutation
states outside the root/full-read mandatory set, while every selected non-root
state has a selected Hamming-1 predecessor and every partial group is hit.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping, Sequence

from solver_probe import ExactInstance, check_selected, popcount, vertex_set_digest


INF = 10**12


def _positive_finite_seconds(value: float, *, name: str) -> float:
    """Normalize a finite positive duration and reject fail-open NaN/Inf."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be finite and positive")
    try:
        normalized = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite and positive") from error
    if normalized <= 0 or not math.isfinite(normalized):
        raise ValueError(f"{name} must be finite and positive")
    return normalized


def mask_count(mask: int) -> int:
    """Return the population count without requiring Python 3.10 int.bit_count."""
    return bin(int(mask)).count("1")


def iter_mask_indices(mask: int):
    """Yield set-bit indices from least to greatest."""
    value = int(mask)
    while value:
        lowest = value & -value
        yield lowest.bit_length() - 1
        value ^= lowest


def optimal_family_digest(k: int, vertex_sets: Iterable[Iterable[int]]) -> str:
    """Hash the unordered family of canonical per-vertex-set digests."""
    identifiers = sorted(vertex_set_digest(k, values) for values in vertex_sets)
    payload = {
        "schema": "intersubmod.optimal_vertex_family_digest.v1",
        "k": int(k),
        "vertex_set_ids": identifiers,
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class ReductionStats:
    input_vertices: int
    downward_closure_vertices: int
    input_groups: int
    active_groups: int
    duplicate_groups_removed: int
    required_hit_groups_removed: int
    dominated_groups_removed: int
    singleton_groups_forced: int
    forced_vertices: int


@dataclass(frozen=True)
class PreparedBitsetProblem:
    """Bitset form of one exact structural instance."""

    instance: ExactInstance
    vertices: tuple[int, ...]
    vertex_to_index: dict[int, int]
    active_mask: int
    root_index: int
    original_mandatory_mask: int
    required_mask: int
    cost_mask: int
    group_masks: tuple[int, ...]
    predecessor_masks: tuple[int, ...]
    children_masks: tuple[int, ...]
    node_costs: tuple[int, ...]
    terminal_domains: tuple[int, ...]
    ranks_descending: tuple[int, ...]
    reduction: ReductionStats

    @property
    def q(self) -> int:
        return len(self.terminal_domains)

    def selected_vertices(self, selected_mask: int) -> tuple[int, ...]:
        return tuple(self.vertices[index] for index in iter_mask_indices(selected_mask))

    def selected_cost(self, selected_mask: int) -> int:
        return mask_count(selected_mask & self.cost_mask)


def _group_antichain(group_masks: Sequence[int]) -> tuple[list[int], int]:
    """Keep inclusion-minimal nonzero domains; hitting them hits every superset."""
    unique = sorted(set(int(mask) for mask in group_masks if mask), key=lambda value: (mask_count(value), value))
    kept: list[int] = []
    removed = len(group_masks) - len(unique)
    for domain in unique:
        if any((smaller & domain) == smaller for smaller in kept):
            removed += 1
        else:
            kept.append(domain)
    return kept, removed


def prepare_bitset_problem(instance: ExactInstance) -> PreparedBitsetProblem:
    """Apply exact structural reductions and construct fixed integer bitsets."""
    full_vertices = tuple(instance.vertices)
    full_index = {vertex: index for index, vertex in enumerate(full_vertices)}

    original_mandatory_full = 0
    for vertex in instance.mandatory:
        original_mandatory_full |= 1 << full_index[vertex]
    required_full = original_mandatory_full

    raw_groups: list[int] = []
    for group in instance.groups:
        mask = 0
        for vertex in group:
            mask |= 1 << full_index[vertex]
        if not mask:
            raise ValueError("partial group has no vertex in the exact universe")
        raw_groups.append(mask)

    input_group_count = len(raw_groups)
    unique_groups = sorted(set(raw_groups))
    duplicate_removed = input_group_count - len(unique_groups)
    groups = list(unique_groups)
    required_hit_removed = 0
    dominated_removed = 0
    singleton_groups_forced = 0
    forced_full = 0

    while True:
        uncovered: list[int] = []
        for group in groups:
            if group & required_full:
                required_hit_removed += 1
            else:
                uncovered.append(group)
        groups, removed = _group_antichain(uncovered)
        dominated_removed += removed

        singleton_union = 0
        singleton_count = 0
        for group in groups:
            if mask_count(group) == 1:
                singleton_union |= group
                singleton_count += 1
        new_forced = singleton_union & ~required_full
        if not new_forced:
            break
        singleton_groups_forced += singleton_count
        forced_full |= new_forced
        required_full |= new_forced

    # Remove the singleton groups hit by the last forcing pass.
    final_groups: list[int] = []
    for group in groups:
        if group & required_full:
            required_hit_removed += 1
        else:
            final_groups.append(group)
    groups, removed = _group_antichain(final_groups)
    dominated_removed += removed

    # Every minimum solution can be represented inside the ancestor closure of
    # mandatory/forced vertices and the still-active group domains.
    seed_values: set[int] = {
        full_vertices[index] for index in iter_mask_indices(required_full)
    }
    for group in groups:
        seed_values.update(full_vertices[index] for index in iter_mask_indices(group))
    closure_values: set[int] = {0}
    for seed in seed_values:
        current = seed
        while True:
            closure_values.add(current)
            if current == 0:
                break
            current = (current - 1) & seed

    vertices = tuple(vertex for vertex in full_vertices if vertex in closure_values)
    vertex_to_index = {vertex: index for index, vertex in enumerate(vertices)}
    active_mask = (1 << len(vertices)) - 1
    root_index = vertex_to_index[0]

    original_mandatory_mask = 0
    required_mask = 0
    for vertex in instance.mandatory:
        original_mandatory_mask |= 1 << vertex_to_index[vertex]
    for full_idx in iter_mask_indices(required_full):
        required_mask |= 1 << vertex_to_index[full_vertices[full_idx]]
    cost_mask = active_mask & ~original_mandatory_mask

    group_masks: list[int] = []
    for group in groups:
        active_group = 0
        for full_idx in iter_mask_indices(group):
            vertex = full_vertices[full_idx]
            if vertex in vertex_to_index:
                active_group |= 1 << vertex_to_index[vertex]
        if not active_group:
            raise AssertionError("downward closure removed an active group")
        group_masks.append(active_group)

    predecessor_masks: list[int] = []
    children_masks = [0 for _ in vertices]
    for child_index, vertex in enumerate(vertices):
        predecessor_mask = 0
        for predecessor in instance.predecessors(vertex):
            if predecessor in vertex_to_index:
                parent_index = vertex_to_index[predecessor]
                predecessor_mask |= 1 << parent_index
                children_masks[parent_index] |= 1 << child_index
        predecessor_masks.append(predecessor_mask)

    node_costs = tuple(
        0 if vertex in instance.mandatory else 1 for vertex in vertices
    )
    terminal_domains: list[int] = [
        1 << index
        for index in iter_mask_indices(required_mask)
        if index != root_index
    ]
    terminal_domains.extend(group_masks)

    ranks_descending = tuple(
        sorted(
            range(len(vertices)),
            key=lambda index: (popcount(vertices[index]), vertices[index]),
            reverse=True,
        )
    )
    return PreparedBitsetProblem(
        instance=instance,
        vertices=vertices,
        vertex_to_index=vertex_to_index,
        active_mask=active_mask,
        root_index=root_index,
        original_mandatory_mask=original_mandatory_mask,
        required_mask=required_mask,
        cost_mask=cost_mask,
        group_masks=tuple(group_masks),
        predecessor_masks=tuple(predecessor_masks),
        children_masks=tuple(children_masks),
        node_costs=node_costs,
        terminal_domains=tuple(terminal_domains),
        ranks_descending=ranks_descending,
        reduction=ReductionStats(
            input_vertices=len(full_vertices),
            downward_closure_vertices=len(vertices),
            input_groups=input_group_count,
            active_groups=len(group_masks),
            duplicate_groups_removed=duplicate_removed,
            required_hit_groups_removed=required_hit_removed,
            dominated_groups_removed=dominated_removed,
            singleton_groups_forced=singleton_groups_forced,
            forced_vertices=mask_count(forced_full),
        ),
    )


@dataclass
class BitsetBnbStats:
    visited_states: int = 0
    memo_hits: int = 0
    propagated_singletons: int = 0
    raw_obligations: int = 0
    antichain_obligations: int = 0
    dynamic_dominated_removed: int = 0
    pruned_infeasible: int = 0
    pruned_cost: int = 0
    pruned_bound: int = 0
    feasible_leaves: int = 0
    max_depth: int = 0
    deadline_checks: int = 0


@dataclass
class BitsetBnbResult:
    backend: str
    k: int
    status: str
    complete: bool
    objective_certified: bool
    family_complete: bool
    objective: int | None
    incumbent_objective: int | None
    vertex_sets: list[tuple[int, ...]]
    elapsed_seconds: float
    lower_bound_name: str
    stats: BitsetBnbStats
    reduction: ReductionStats
    family_digest: str | None

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["vertex_sets"] = [list(values) for values in self.vertex_sets]
        result["vertex_set_ids"] = [
            vertex_set_digest(self.k, values) for values in self.vertex_sets
        ]
        return result


class BitsetObligationBnb:
    """Exact obligation B&B using integer bitsets and dynamic antichains."""

    def __init__(
        self,
        problem: PreparedBitsetProblem,
        *,
        time_limit_seconds: float = 30.0,
        max_sets: int | None = None,
        certified_target_objective: int | None = None,
    ) -> None:
        validated_time_limit = _positive_finite_seconds(
            time_limit_seconds,
            name="time_limit_seconds",
        )
        if max_sets is not None and max_sets < 0:
            raise ValueError("max_sets must be nonnegative or None")
        if certified_target_objective is not None and certified_target_objective < 0:
            raise ValueError("certified_target_objective must be nonnegative")
        self.problem = problem
        self.time_limit_seconds = validated_time_limit
        self.max_sets = max_sets
        self.certified_target = certified_target_objective
        self.best = (
            math.inf
            if certified_target_objective is None
            else int(certified_target_objective)
        )
        self.solutions: set[tuple[int, ...]] = set()
        self.memo: set[tuple[int, int]] = set()
        self.stats = BitsetBnbStats()
        self.timed_out = False
        self.set_cap_reached = False
        self.deadline = math.inf

    def _deadline_reached(self) -> bool:
        self.stats.deadline_checks += 1
        if time.perf_counter() >= self.deadline:
            self.timed_out = True
            return True
        return False

    def _root_reachable(
        self,
        index: int,
        excluded: int,
        memo: dict[int, bool],
    ) -> bool:
        if index in memo:
            return memo[index]
        bit = 1 << index
        if excluded & bit:
            value = False
        elif index == self.problem.root_index:
            value = True
        else:
            value = any(
                self._root_reachable(parent, excluded, memo)
                for parent in iter_mask_indices(self.problem.predecessor_masks[index])
            )
        memo[index] = value
        return value

    def _obligations(
        self,
        selected: int,
        excluded: int,
    ) -> tuple[list[tuple[str, int, int]], bool]:
        raw: list[tuple[str, int, int]] = []
        reachable_memo: dict[int, bool] = {}
        for group_index, group in enumerate(self.problem.group_masks):
            if not (selected & group):
                domain = group & ~excluded
                domain = sum(
                    1 << index
                    for index in iter_mask_indices(domain)
                    if self._root_reachable(index, excluded, reachable_memo)
                )
                if not domain:
                    return [], False
                raw.append(("group", group_index, domain))
        root_bit = 1 << self.problem.root_index
        for child_index in iter_mask_indices(selected & ~root_bit):
            predecessors = self.problem.predecessor_masks[child_index]
            if not (selected & predecessors):
                domain = predecessors & ~excluded
                domain = sum(
                    1 << index
                    for index in iter_mask_indices(domain)
                    if self._root_reachable(index, excluded, reachable_memo)
                )
                if not domain:
                    return [], False
                raw.append(("parent", child_index, domain))

        self.stats.raw_obligations += len(raw)
        ordered = sorted(
            raw,
            key=lambda item: (
                mask_count(item[2]),
                item[2],
                0 if item[0] == "group" else 1,
                item[1],
            ),
        )
        kept: list[tuple[str, int, int]] = []
        for obligation in ordered:
            domain = obligation[2]
            if any((smaller[2] & domain) == smaller[2] for smaller in kept):
                self.stats.dynamic_dominated_removed += 1
            else:
                kept.append(obligation)
        self.stats.antichain_obligations += len(kept)
        return kept, True

    def _propagate(
        self,
        selected: int,
        excluded: int,
    ) -> tuple[int, int, list[tuple[str, int, int]]] | None:
        if selected & excluded:
            return None
        current = selected
        while True:
            reachable_memo: dict[int, bool] = {}
            if any(
                not self._root_reachable(index, excluded, reachable_memo)
                for index in iter_mask_indices(current)
            ):
                return None
            obligations, possible = self._obligations(current, excluded)
            if not possible:
                return None
            singleton_union = 0
            for _, _, domain in obligations:
                if mask_count(domain) == 1:
                    singleton_union |= domain
            new_values = singleton_union & ~current
            if not new_values:
                return current, excluded, obligations
            if new_values & excluded:
                return None
            current |= new_values
            self.stats.propagated_singletons += mask_count(new_values)

    def _connection_cost(
        self,
        index: int,
        selected: int,
        excluded: int,
        memo: dict[int, float],
    ) -> float:
        if index in memo:
            return memo[index]
        bit = 1 << index
        if excluded & bit:
            value = math.inf
        elif index == self.problem.root_index:
            value = 0.0
        else:
            parent_cost = min(
                (
                    self._connection_cost(parent, selected, excluded, memo)
                    for parent in iter_mask_indices(self.problem.predecessor_masks[index])
                ),
                default=math.inf,
            )
            value = parent_cost + (
                0.0 if selected & bit else float(self.problem.node_costs[index])
            )
        memo[index] = value
        return value

    def _lower_bound(
        self,
        selected: int,
        excluded: int,
        obligations: Sequence[tuple[str, int, int]],
    ) -> int:
        current_cost = self.problem.selected_cost(selected)
        if not obligations:
            return current_cost

        domains = [domain for _, _, domain in obligations]
        used = 0
        packing = 0
        for domain in sorted(domains, key=lambda value: (mask_count(value), value)):
            if not (used & domain):
                used |= domain
                packing += 1

        coverage = [0 for _ in self.problem.vertices]
        for domain in domains:
            for index in iter_mask_indices(domain):
                coverage[index] += 1
        maximum_coverage = max(coverage, default=1)
        cover_bound = math.ceil(len(domains) / max(1, maximum_coverage))

        connection_memo: dict[int, float] = {}
        selected_connection = max(
            (
                self._connection_cost(index, selected, excluded, connection_memo)
                for index in iter_mask_indices(selected)
            ),
            default=0.0,
        )
        group_connection = max(
            (
                min(
                    self._connection_cost(index, selected, excluded, connection_memo)
                    for index in iter_mask_indices(domain)
                )
                for kind, _, domain in obligations
                if kind == "group"
            ),
            default=0.0,
        )
        additional = max(
            float(packing),
            float(cover_bound),
            selected_connection,
            group_connection,
        )
        if not math.isfinite(additional):
            return INF
        return current_cost + int(math.ceil(additional - 1e-12))

    def _candidate_order(
        self,
        domain: int,
        selected: int,
        excluded: int,
        obligations: Sequence[tuple[str, int, int]],
    ) -> list[int]:
        connection_memo: dict[int, float] = {}
        return sorted(
            iter_mask_indices(domain),
            key=lambda index: (
                self._connection_cost(index, selected, excluded, connection_memo),
                -sum(bool(other & (1 << index)) for _, _, other in obligations),
                popcount(self.problem.vertices[index]),
                self.problem.vertices[index],
            ),
        )

    def _search(self, selected: int, excluded: int, depth: int) -> None:
        if self._deadline_reached():
            return
        if self.max_sets is not None and len(self.solutions) >= self.max_sets:
            self.set_cap_reached = True
            return

        propagated = self._propagate(selected, excluded)
        if propagated is None:
            self.stats.pruned_infeasible += 1
            return
        if self._deadline_reached():
            return
        selected, excluded, obligations = propagated
        current_cost = self.problem.selected_cost(selected)
        if current_cost > self.best:
            self.stats.pruned_cost += 1
            return
        key = (selected, excluded)
        if key in self.memo:
            self.stats.memo_hits += 1
            return
        self.memo.add(key)
        self.stats.visited_states += 1
        self.stats.max_depth = max(self.stats.max_depth, depth)

        lower_bound = self._lower_bound(selected, excluded, obligations)
        if self._deadline_reached():
            return
        if lower_bound > self.best:
            self.stats.pruned_bound += 1
            return
        if not obligations:
            if self._deadline_reached():
                return
            self.stats.feasible_leaves += 1
            vertices = self.problem.selected_vertices(selected)
            validation = check_selected(self.problem.instance, vertices)
            if not validation["pass"]:
                raise AssertionError(f"bitset B&B emitted an invalid set: {validation}")
            cost = self.problem.selected_cost(selected)
            if self.certified_target is not None and cost < self.certified_target:
                raise AssertionError("subset-DP target exceeds a feasible structural objective")
            if cost < self.best:
                self.best = cost
                self.solutions = {vertices}
            elif cost == self.best:
                self.solutions.add(vertices)
            return

        chosen = min(
            obligations,
            key=lambda item: (
                mask_count(item[2]),
                0 if item[0] == "group" else 1,
                item[1],
            ),
        )
        earlier = 0
        for candidate_index in self._candidate_order(
            chosen[2], selected, excluded, obligations
        ):
            if self.timed_out or self.set_cap_reached:
                return
            candidate = 1 << candidate_index
            self._search(selected | candidate, excluded | earlier, depth + 1)
            earlier |= candidate

    def run(self) -> BitsetBnbResult:
        started = time.perf_counter()
        self.deadline = started + self.time_limit_seconds
        self._search(self.problem.required_mask, 0, 0)
        # A unique final leaf or the final DP-like scan can finish after the
        # last recursive entry check.  Never issue a complete certificate then.
        if time.perf_counter() >= self.deadline:
            self.timed_out = True
        elapsed = time.perf_counter() - started
        traversal_complete = not self.timed_out and not self.set_cap_reached
        finite_incumbent = math.isfinite(self.best) and bool(self.solutions)
        incumbent = int(self.best) if finite_incumbent else None
        objective_certified = self.certified_target is not None or (
            traversal_complete and finite_incumbent
        )
        objective = (
            int(self.certified_target)
            if self.certified_target is not None
            else (incumbent if objective_certified else None)
        )
        family_complete = traversal_complete and finite_incumbent

        if self.timed_out:
            if objective_certified:
                status = "CANDIDATE_SET_INCOMPLETE_DEADLINE"
            elif finite_incumbent:
                status = "FEASIBLE_UNPROVEN_DEADLINE"
            else:
                status = "NO_FEASIBLE_CERTIFICATE_DEADLINE"
        elif self.set_cap_reached:
            if objective_certified:
                status = "CANDIDATE_SET_INCOMPLETE_MAX_SETS"
            elif finite_incumbent:
                status = "FEASIBLE_UNPROVEN_MAX_SETS"
            else:
                status = "NO_FEASIBLE_CERTIFICATE_MAX_SETS"
        elif family_complete:
            status = "CANDIDATE_SET_COMPLETE"
        else:
            status = "INFEASIBLE"

        if traversal_complete and self.certified_target is not None and not self.solutions:
            raise AssertionError("certified target produced no feasible structural family")
        vertex_sets = sorted(self.solutions)
        return BitsetBnbResult(
            backend="bitset_obligation_bnb_v1",
            k=self.problem.instance.k,
            status=status,
            complete=family_complete,
            objective_certified=objective_certified,
            family_complete=family_complete,
            objective=objective,
            incumbent_objective=incumbent,
            vertex_sets=vertex_sets,
            elapsed_seconds=elapsed,
            lower_bound_name="max(disjoint_packing,coverage,root_connection)",
            stats=self.stats,
            reduction=self.problem.reduction,
            family_digest=(
                optimal_family_digest(self.problem.instance.k, vertex_sets)
                if vertex_sets
                else None
            ),
        )


@dataclass
class SubsetDpStats:
    q: int
    vertex_count: int
    table_cells: int
    directed_edge_count: int
    estimated_dense_bytes: int
    proxy_operations: int
    max_table_cells: int
    max_dense_bytes: int
    max_proxy_operations: int
    merge_attempts: int = 0
    edge_relaxations: int = 0
    finite_cells: int = 0
    deadline_checks: int = 0


@dataclass
class SubsetDpResult:
    backend: str
    status: str
    objective_certified: bool
    family_complete: bool
    objective: int | None
    elapsed_seconds: float
    stats: SubsetDpStats
    reduction: ReductionStats

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def solve_group_terminal_subset_dp(
    problem: PreparedBitsetProblem,
    *,
    max_q: int = 8,
    max_table_cells: int = 1_048_576,
    max_dense_bytes: int = 512 * 1024 * 1024,
    bytes_per_dense_cell: int = 16,
    max_proxy_operations: int = 50_000_000,
    time_limit_seconds: float | None = None,
) -> SubsetDpResult:
    """Prove h* using a node-weighted directed-Steiner subset DP.

    ``dp[S][v]`` excludes the cost of ``v`` itself.  Moving from a parent to a
    child pays the child's cost, so merging two terminal subsets at ``v`` does
    not double count ``v``.
    """
    started = time.perf_counter()
    q = problem.q
    directed_edge_count = sum(mask_count(mask) for mask in problem.children_masks)
    table_cells = (1 << q) * len(problem.vertices)
    estimated_dense_bytes = table_cells * int(bytes_per_dense_cell)
    proxy_operations = (
        (3**q) * len(problem.vertices)
        + (1 << q) * directed_edge_count
    )
    stats = SubsetDpStats(
        q=q,
        vertex_count=len(problem.vertices),
        table_cells=table_cells,
        directed_edge_count=directed_edge_count,
        estimated_dense_bytes=estimated_dense_bytes,
        proxy_operations=proxy_operations,
        max_table_cells=int(max_table_cells),
        max_dense_bytes=int(max_dense_bytes),
        max_proxy_operations=int(max_proxy_operations),
    )
    if (
        max_q < 0
        or max_table_cells < 1
        or max_dense_bytes < 1
        or bytes_per_dense_cell < 1
        or max_proxy_operations < 1
    ):
        raise ValueError("subset-DP route limits must be positive")
    validated_time_limit = (
        None
        if time_limit_seconds is None
        else _positive_finite_seconds(
            time_limit_seconds,
            name="time_limit_seconds",
        )
    )
    deadline = (
        math.inf
        if validated_time_limit is None
        else started + validated_time_limit
    )

    def deadline_reached() -> bool:
        stats.deadline_checks += 1
        return time.perf_counter() >= deadline

    def deadline_result() -> SubsetDpResult:
        return SubsetDpResult(
            backend="group_terminal_subset_dp_v1",
            status="OBJECTIVE_INCOMPLETE_DEADLINE",
            objective_certified=False,
            family_complete=False,
            objective=None,
            elapsed_seconds=time.perf_counter() - started,
            stats=stats,
            reduction=problem.reduction,
        )

    if q > max_q:
        return SubsetDpResult(
            backend="group_terminal_subset_dp_v1",
            status="ROUTE_NOT_ELIGIBLE_Q_GT_MAX",
            objective_certified=False,
            family_complete=False,
            objective=None,
            elapsed_seconds=time.perf_counter() - started,
            stats=stats,
            reduction=problem.reduction,
        )
    if (
        table_cells > max_table_cells
        or estimated_dense_bytes > max_dense_bytes
        or proxy_operations > max_proxy_operations
    ):
        return SubsetDpResult(
            backend="group_terminal_subset_dp_v1",
            status="ROUTE_NOT_ELIGIBLE_RESOURCE_GATE",
            objective_certified=False,
            family_complete=False,
            objective=None,
            elapsed_seconds=time.perf_counter() - started,
            stats=stats,
            reduction=problem.reduction,
        )
    if deadline_reached():
        return deadline_result()
    if q == 0:
        return SubsetDpResult(
            backend="group_terminal_subset_dp_v1",
            status="OPTIMAL_VALUE_CERTIFIED",
            objective_certified=True,
            family_complete=False,
            objective=0,
            elapsed_seconds=time.perf_counter() - started,
            stats=stats,
            reduction=problem.reduction,
        )

    n_vertices = len(problem.vertices)
    table = [[INF for _ in range(n_vertices)] for _ in range(1 << q)]
    if deadline_reached():
        return deadline_result()
    for terminal_index, domain in enumerate(problem.terminal_domains):
        row = table[1 << terminal_index]
        for vertex_index in iter_mask_indices(domain):
            row[vertex_index] = 0
        if deadline_reached():
            return deadline_result()

    masks = sorted(range(1, 1 << q), key=lambda value: (mask_count(value), value))
    for terminal_mask in masks:
        if deadline_reached():
            return deadline_result()
        row = table[terminal_mask]
        submask = (terminal_mask - 1) & terminal_mask
        while submask:
            other = terminal_mask ^ submask
            if submask < other:
                left = table[submask]
                right = table[other]
                for vertex_index in range(n_vertices):
                    stats.merge_attempts += 1
                    if (
                        stats.merge_attempts & 1023 == 0
                        and deadline_reached()
                    ):
                        return deadline_result()
                    candidate = left[vertex_index] + right[vertex_index]
                    if candidate < row[vertex_index]:
                        row[vertex_index] = candidate
            submask = (submask - 1) & terminal_mask

        # The hypercube is a DAG.  Descending rank finalizes every child before
        # its predecessors, giving an exact nonnegative shortest-path closure.
        for vertex_index in problem.ranks_descending:
            for child_index in iter_mask_indices(problem.children_masks[vertex_index]):
                stats.edge_relaxations += 1
                if (
                    stats.edge_relaxations & 1023 == 0
                    and deadline_reached()
                ):
                    return deadline_result()
                if row[child_index] >= INF:
                    continue
                candidate = problem.node_costs[child_index] + row[child_index]
                if candidate < row[vertex_index]:
                    row[vertex_index] = candidate
        if deadline_reached():
            return deadline_result()

    objective = table[(1 << q) - 1][problem.root_index]
    stats.finite_cells = sum(
        value < INF for row in table[1:] for value in row
    )
    if deadline_reached():
        return deadline_result()
    if objective >= INF:
        status = "INFEASIBLE_CERTIFIED"
        certified_objective = None
    else:
        status = "OPTIMAL_VALUE_CERTIFIED"
        certified_objective = int(objective)
    return SubsetDpResult(
        backend="group_terminal_subset_dp_v1",
        status=status,
        objective_certified=True,
        family_complete=False,
        objective=certified_objective,
        elapsed_seconds=time.perf_counter() - started,
        stats=stats,
        reduction=problem.reduction,
    )


def solve_with_certificates(
    instance: ExactInstance,
    *,
    q_max: int = 8,
    time_limit_seconds: float = 30.0,
    max_sets: int | None = None,
    max_dp_table_cells: int = 1_048_576,
    max_dp_dense_bytes: int = 512 * 1024 * 1024,
    max_dp_proxy_operations: int = 50_000_000,
) -> dict[str, Any]:
    """Route objective proof through small-q DP and enumerate via bitset B&B."""
    validated_time_limit = _positive_finite_seconds(
        time_limit_seconds,
        name="time_limit_seconds",
    )
    total_started = time.perf_counter()
    prepare_started = time.perf_counter()
    problem = prepare_bitset_problem(instance)
    prepare_elapsed = time.perf_counter() - prepare_started
    remaining_for_dp = max(
        1e-12,
        validated_time_limit - (time.perf_counter() - total_started),
    )
    objective = solve_group_terminal_subset_dp(
        problem,
        max_q=q_max,
        max_table_cells=max_dp_table_cells,
        max_dense_bytes=max_dp_dense_bytes,
        max_proxy_operations=max_dp_proxy_operations,
        time_limit_seconds=remaining_for_dp,
    )
    target = objective.objective if objective.objective_certified else None
    remaining_for_bnb = max(
        1e-12,
        validated_time_limit - (time.perf_counter() - total_started),
    )
    family = BitsetObligationBnb(
        problem,
        time_limit_seconds=remaining_for_bnb,
        max_sets=max_sets,
        certified_target_objective=target,
    ).run()
    if objective.objective_certified and family.incumbent_objective is not None:
        if objective.objective != family.incumbent_objective:
            raise AssertionError("subset-DP and bitset B&B objective mismatch")
    total_elapsed = time.perf_counter() - total_started
    total_deadline_exceeded = total_elapsed >= validated_time_limit
    if total_deadline_exceeded and family.family_complete:
        family.complete = False
        family.family_complete = False
        if family.objective_certified:
            family.status = "CANDIDATE_SET_INCOMPLETE_DEADLINE"
        elif family.incumbent_objective is not None:
            family.status = "FEASIBLE_UNPROVEN_DEADLINE"
        else:
            family.status = "NO_FEASIBLE_CERTIFICATE_DEADLINE"
    return {
        "schema": "intersubmod.optimized_hypercube_backend.solve.v1",
        "scope": "BOUNDED_RESEARCH_BACKEND_NOT_PRODUCTION",
        "prepare_elapsed_seconds": prepare_elapsed,
        "total_elapsed_seconds": total_elapsed,
        "total_unit_deadline_seconds": validated_time_limit,
        "total_unit_deadline_exceeded": total_deadline_exceeded,
        "deadline_scope": "prepare+subset_dp+bitset_bnb",
        "reduction": asdict(problem.reduction),
        "q": problem.q,
        "objective_certificate": objective.to_dict(),
        "family_certificate": family.to_dict(),
        "ranking_allowed": bool(family.family_complete),
        "failure_policy": (
            "ALLOW_COMPLETE_FAMILY_RANKING"
            if family.family_complete
            else "ABSTAIN_INCOMPLETE_FAMILY"
        ),
    }


def fixed_node_parent_mapping(
    selected_vertices: Iterable[int],
    *,
    edge_scores: Mapping[tuple[int, int], float] | None = None,
    beta: float = 1.0,
    tie_tolerance: float = 1e-12,
) -> dict[str, Any]:
    """Factor an additive fixed-N parent score child by child in O(E_N)."""
    if beta <= 0 or not math.isfinite(beta):
        raise ValueError("beta must be finite and positive")
    if tie_tolerance < 0 or not math.isfinite(tie_tolerance):
        raise ValueError("tie_tolerance must be finite and nonnegative")
    selected = tuple(sorted(set(int(value) for value in selected_vertices)))
    selected_set = set(selected)
    if 0 not in selected_set:
        raise ValueError("fixed node set must include root 0")

    legal_parents: dict[int, tuple[int, ...]] = {}
    best_parents: dict[int, tuple[int, ...]] = {}
    posterior: dict[int, dict[int, float]] = {}
    total_tree_count = 1
    best_tree_count = 1
    best_total_score = 0.0
    log_partition = 0.0

    for child in selected:
        if child == 0:
            continue
        parents = tuple(
            sorted(
                child ^ (1 << bit)
                for bit in range(child.bit_length())
                if child & (1 << bit)
                and (child ^ (1 << bit)) in selected_set
            )
        )
        if not parents:
            raise ValueError(f"selected child {child} has no selected predecessor")
        legal_parents[child] = parents
        scores: list[float] = []
        for parent in parents:
            if edge_scores is None:
                score = 0.0
            else:
                key = (parent, child)
                if key not in edge_scores:
                    raise ValueError(f"missing local edge score for {key}")
                score = float(edge_scores[key])
                if not math.isfinite(score):
                    raise ValueError(f"non-finite local edge score for {key}")
            if not math.isfinite(beta * score):
                raise ValueError(
                    f"scaled local edge score overflows for {(parent, child)}"
                )
            scores.append(score)
        maximum = max(scores)
        ties = tuple(
            parent
            for parent, score in zip(parents, scores)
            if abs(score - maximum) <= tie_tolerance
        )
        best_parents[child] = ties
        total_tree_count *= len(parents)
        best_tree_count *= len(ties)
        best_total_score += maximum
        if not math.isfinite(best_total_score):
            raise ValueError("best parent score accumulation overflows")

        shifted = [beta * (score - maximum) for score in scores]
        normalizer = sum(math.exp(value) for value in shifted)
        log_partition += beta * maximum + math.log(normalizer)
        if not math.isfinite(log_partition):
            raise ValueError("parent log-partition accumulation overflows")
        posterior[child] = {
            parent: math.exp(value) / normalizer
            for parent, value in zip(parents, shifted)
        }

    return {
        "schema": "intersubmod.fixed_node_parent_mapping.v1",
        "selected_vertices": list(selected),
        "legal_parents": {
            str(child): list(parents) for child, parents in legal_parents.items()
        },
        "best_parents": {
            str(child): list(parents) for child, parents in best_parents.items()
        },
        "canonical_best_parent": {
            str(child): min(parents) for child, parents in best_parents.items()
        },
        "total_parent_tree_count": total_tree_count,
        "best_parent_tree_count": best_tree_count,
        "best_total_score": best_total_score,
        "unique_best_tree": best_tree_count == 1,
        "log_partition": log_partition,
        "parent_posterior": {
            str(child): {str(parent): probability for parent, probability in values.items()}
            for child, values in posterior.items()
        },
        "complexity": {
            "legal_edges_scored_once": sum(len(values) for values in legal_parents.values()),
            "parent_cartesian_product_materialized": False,
        },
    }


def project_unary_hidden_chains(
    selected_vertices: Iterable[int],
    parent_of: Mapping[int, int],
    evidence_vertices: Iterable[int],
) -> dict[str, Any]:
    """Collapse unsupported unary connectors as an unresolved projection only."""
    selected = set(int(value) for value in selected_vertices)
    evidence = set(int(value) for value in evidence_vertices) | {0}
    if 0 not in selected:
        raise ValueError("selected vertices must include root")
    if set(parent_of) != selected - {0}:
        raise ValueError("parent_of must map every selected non-root vertex exactly once")

    children: dict[int, list[int]] = {vertex: [] for vertex in selected}
    for child, parent in parent_of.items():
        child = int(child)
        parent = int(parent)
        if child not in selected or parent not in selected:
            raise ValueError("parent mapping references a vertex outside the selected set")
        difference = child ^ parent
        if parent & ~child or mask_count(difference) != 1:
            raise ValueError("parent mapping must use monotone Hamming-1 edges")
        children[parent].append(child)

    kept = {
        vertex
        for vertex in selected
        if vertex in evidence or len(children[vertex]) != 1
    }
    kept.add(0)
    collapsed = selected - kept
    projected_edges: list[dict[str, Any]] = []
    for child in sorted(kept - {0}):
        cursor = child
        mutations: list[int] = []
        traversed: list[int] = []
        while cursor not in kept or cursor == child:
            parent = int(parent_of[cursor])
            difference = cursor ^ parent
            mutations.append(difference.bit_length() - 1)
            if cursor != child:
                traversed.append(cursor)
            cursor = parent
            if cursor in kept:
                break
        projected_edges.append(
            {
                "parent": cursor,
                "child": child,
                "mutation_bits": sorted(mutations),
                "n_mutations": len(mutations),
                "unresolved_order_count": math.factorial(len(mutations)),
                "edge_type": (
                    "UNIT_MUTATION_EDGE"
                    if len(mutations) == 1
                    else "MULTI_MUTATION_EDGE_EQUIVALENCE"
                ),
                "mutation_order": (
                    "IDENTIFIED_BY_UNIT_EDGE"
                    if len(mutations) == 1
                    else "UNRESOLVED_NO_READ_EVIDENCE"
                ),
                "collapsed_internal_vertices": sorted(traversed),
            }
        )
    return {
        "schema": "intersubmod.unary_hidden_chain_projection.v1",
        "status": "PROJECTION_ONLY_NOT_STRUCTURAL_SOLVER_REDUCTION",
        "original_selected_vertices": sorted(selected),
        "evidence_vertices": sorted(evidence & selected),
        "kept_vertices": sorted(kept),
        "collapsed_vertices": sorted(collapsed),
        "projected_edges": projected_edges,
        "objective_changed": False,
        "objective_field_semantics": (
            "The original structural candidate and its objective are retained."
        ),
        "projected_model_objective_equivalence_claimed": False,
        "equivalence_scope": (
            "One fixed parent tree with no evidence on collapsed unary vertices; "
            "not a cross-candidate family compressor."
        ),
        "original_candidate_retained_for_audit": True,
    }


def _validated_selected_state_set(
    selected_vertices: Iterable[int],
    *,
    k: int,
) -> tuple[int, ...]:
    if not isinstance(k, int) or isinstance(k, bool) or k < 1:
        raise ValueError("k must be a positive integer")
    values = tuple(selected_vertices)
    if any(not isinstance(value, int) or isinstance(value, bool) for value in values):
        raise ValueError("selected mutation states must be integers")
    selected = tuple(sorted(set(values)))
    if 0 not in selected:
        raise ValueError("selected mutation states must include root 0")
    upper = 1 << k
    outside = [value for value in selected if value < 0 or value >= upper]
    if outside:
        raise ValueError(f"selected mutation states outside k-bit universe: {outside}")
    return selected


def strict_gain_once_compatible(
    selected_vertices: Iterable[int],
    *,
    k: int,
) -> tuple[bool, list[tuple[int, int]]]:
    """Apply the rooted binary three-gamete condition to a fixed node set."""
    selected = _validated_selected_state_set(selected_vertices, k=k)
    violations: list[tuple[int, int]] = []
    for left in range(k):
        for right in range(left + 1, k):
            categories = set()
            for vertex in selected:
                pair = (
                    1 if vertex & (1 << left) else 0,
                    1 if vertex & (1 << right) else 0,
                )
                if pair != (0, 0):
                    categories.add(pair)
            if {(1, 0), (0, 1), (1, 1)} <= categories:
                violations.append((left, right))
    return not violations, violations


def evaluate_evolutionary_mode(
    selected_vertices: Iterable[int],
    *,
    k: int,
    mode: str,
    eligibility: Mapping[str, bool] | None = None,
) -> dict[str, Any]:
    """Fail closed for strict/Dollo biology modes outside their evidence gate."""
    allowed_modes = {
        "M0_RECURRENCE_ALLOWED",
        "M1_STRICT_INFINITE_SITES",
        "M2_LOSS_SUPPORTED_DOLLO",
    }
    if mode not in allowed_modes:
        raise ValueError(f"unknown evolutionary mode: {mode}")
    selected = _validated_selected_state_set(selected_vertices, k=k)
    if mode == "M0_RECURRENCE_ALLOWED":
        return {
            "mode": mode,
            "status": "ELIGIBLE_BASELINE",
            "constraint_applied": True,
            "publication_label": "rooted_monotone_molecule_state_topology",
        }
    if mode == "M2_LOSS_SUPPORTED_DOLLO":
        return {
            "mode": mode,
            "status": "UNRESOLVED_LOSS_SUPPORTED_MODEL_NOT_IMPLEMENTED",
            "constraint_applied": False,
            "publication_allowed": False,
            "reason": "CNA/LOH regions require explicit loss and copy-state evidence",
        }

    required_true_gates = (
        "allele_specific_cn_known",
        "copy_neutral_diploid_both_homologs_retained",
        "clonal_loh_absent",
        "loh_absent",
        "subclonal_loh_absent",
        "clonal_deletion_absent",
        "subclonal_deletion_absent",
        "clonal_amplification_absent",
        "clonal_cna_absent",
        "subclonal_cna_absent",
        "subclonal_amplification_absent",
        "wgd_boundary_absent",
        "mutation_loss_absent",
        "phasing_qc_pass",
        "mapping_qc_pass",
        "base_quality_qc_pass",
        "strand_qc_pass",
        "read_independence_qc_pass",
        "allele_specific_duplicated_copies_absent",
        "somatic_variant_qc_pass",
    )
    eligibility = dict(eligibility or {})
    failed = [
        name for name in required_true_gates
        if eligibility.get(name) is not True
    ]
    expected_cn = {"total_cn": 2, "major_cn": 1, "minor_cn": 1}
    invalid_cn = [
        name
        for name, expected in expected_cn.items()
        if (
            isinstance(eligibility.get(name), bool)
            or eligibility.get(name) != expected
        )
    ]
    adverse_evidence_fields = (
        "clonal_loh",
        "clonal_loh_present",
        "subclonal_loh",
        "subclonal_loh_present",
        "deletion",
        "deletion_present",
        "amplification",
        "amplification_present",
        "clonal_cna",
        "subclonal_cna",
        "wgd",
        "wgd_boundary",
        "wgd_boundary_present",
        "duplicated_copies",
        "allele_specific_duplicated_copies",
        "allele_specific_duplicated_copies_present",
        "mutation_loss",
        "mutation_loss_present",
    )
    contradictory = [
        name for name in adverse_evidence_fields
        if eligibility.get(name) not in (None, False, 0)
    ]
    if failed or invalid_cn or contradictory:
        return {
            "mode": mode,
            "status": "ABSTAIN_CN_LOH_GATE",
            "constraint_applied": False,
            "publication_allowed": False,
            "failed_or_missing_gates": failed,
            "invalid_or_missing_copy_state": invalid_cn,
            "contradictory_adverse_evidence": contradictory,
        }
    compatible, violations = strict_gain_once_compatible(selected, k=k)
    orphan_states = [
        vertex
        for vertex in selected
        if vertex != 0
        and not any(
            (vertex ^ (1 << bit)) in selected
            for bit in range(k)
            if vertex & (1 << bit)
        )
    ]
    compatible = compatible and not orphan_states
    return {
        "mode": mode,
        "status": "STRICT_COMPATIBLE" if compatible else "STRICT_INFEASIBLE",
        "constraint_applied": True,
        "publication_allowed": compatible,
        "violating_mutation_pairs": [list(values) for values in violations],
        "orphan_states_without_selected_predecessor": orphan_states,
        "publication_label": "strict-perfect-compatible molecule-state topology",
    }
