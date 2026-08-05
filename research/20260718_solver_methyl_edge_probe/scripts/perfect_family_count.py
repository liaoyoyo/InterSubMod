#!/usr/bin/env python3
"""Exact family counter for the recurrence-free ``h*=m-d`` fast path.

This module is deliberately isolated from the canonical M2 router.  It counts
optimal *vertex-set* families without materializing them, but only after the
following contract is proved by construction:

1. the active mutation universe is the union of ALT coordinates in the full
   and partial read patterns; and
2. at least one recurrence-free rooted mutation forest satisfies every read
   and contains every mandatory full-read state.

Let ``m`` be the number of active mutations.  Every feasible rooted state tree
must acquire each active mutation at least once, so it has at least ``m``
non-root vertices.  A recurrence-free forest constructed by this module has
exactly ``m`` mutation edges.  If there are ``d`` distinct mandatory non-root
states, the solver objective (non-mandatory vertices) is at least ``m-d``.
Therefore a nonzero count certifies ``h*=m-d``; root-only cases have ``d=0``
and hence ``h*=m``.  At that objective every mutation edge label occurs
exactly once, so mutation forests and selected Boolean-state vertex sets are
in bijection.

For a partial pattern with ALT set ``A`` and REF set ``R``, a recurrence-free
mutation forest hits the induced subcube exactly when:

* every pair of mutations in ``A`` is ancestor-comparable; and
* no mutation in ``R`` is an ancestor of a mutation in ``A``.

The first condition forms an undirected comparability graph ``C``.  The second
forms a directed forbidden-ancestor relation ``F``.  The subset dynamic
program uses:

``T[S]``
    Number of valid mutation trees spanning exactly ``S``.

``W[S]``
    Number of valid unordered mutation forests under one virtual parent.

For a tree root ``x``, ``x`` is legal only if it is not forbidden to be an
ancestor of any member of ``S - {x}``, hence

``T[S] = sum_x W[S - {x}]`` over legal roots.

For ``W[S]``, connected components of ``C[S]`` are indivisible: splitting one
between two child trees would make a required-comparable pair incomparable.
The recurrence chooses the child block containing the first component and
then recurses on the remainder, counting every unordered partition once.

The worst-case time is ``O(3^m poly(m))`` and space is ``O(2^m + m^2)``.
When the count is zero, the module does not enumerate recurrence-allowed
``h*>m`` families and returns ``ABSTAIN_UNSUPPORTED_RECURRENCE_COUNT``.
"""

from __future__ import annotations

import dataclasses
import math
import time
from typing import Any, Mapping, Sequence


EXACT_STATUS = "EXACT_H_EQ_M_PERFECT_FAMILY_COUNT"
EXACT_MANDATORY_STATUS = (
    "EXACT_H_EQ_M_MINUS_MANDATORY_PERFECT_FAMILY_COUNT"
)
EXACT_STATUSES = frozenset({EXACT_STATUS, EXACT_MANDATORY_STATUS})
ABSTAIN_STATUS = "ABSTAIN_UNSUPPORTED_RECURRENCE_COUNT"
SCOPE = "ISOLATED_PERFECT_EVENT_H_EQ_M_MINUS_D_FAST_PATH_NOT_PRODUCTION"


def _popcount(value: int) -> int:
    return bin(int(value)).count("1")


def _validate_pattern(pattern: str, *, k: int, allow_x: bool) -> None:
    if len(pattern) != k:
        raise ValueError("pattern length differs from k")
    allowed = {"R", "A", "X"} if allow_x else {"R", "A"}
    invalid = set(pattern) - allowed
    if invalid:
        raise ValueError(f"invalid pattern symbols: {sorted(invalid)}")


@dataclasses.dataclass(frozen=True)
class PerfectFamilyCountResult:
    scope: str
    status: str
    exact: bool
    ranking_allowed: bool
    k: int
    effective_m: int
    mandatory_nonroot_state_count: int
    objective: int | None
    perfect_family_count: int
    comparable_pair_count: int
    forbidden_ancestor_pair_count: int
    subset_count: int
    elapsed_seconds: float
    proof_basis: str

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)


def _induced_components(mask: int, adjacency: Sequence[int]) -> tuple[int, ...]:
    """Return connected-component masks of the graph induced by ``mask``."""
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


def count_perfect_families(
    *,
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
    k: int,
    structural_alt_universe_mask: int | None = None,
    max_m: int = 20,
) -> PerfectFamilyCountResult:
    """Count all recurrence-free optimal vertex-set families.

    A nonzero count is both a construction with ``m`` mutation edges and an
    objective certificate because every active mutation must occur on at least
    one edge.  With ``d`` mandatory non-root full-read states, the solver
    objective is therefore ``h*=m-d``.  Zero means that the recurrence-free
    layer is empty; the function abstains rather than making a claim about the
    larger recurrence-allowed family.
    """
    started = time.perf_counter()
    if type(k) is not int or k < 1:
        raise ValueError("k must be a positive integer")
    if type(max_m) is not int or max_m < 1:
        raise ValueError("max_m must be a positive integer")

    for pattern in full_patterns:
        _validate_pattern(pattern, k=k, allow_x=False)
    for pattern in partial_patterns:
        _validate_pattern(pattern, k=k, allow_x=True)
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
        if universe_mask < 0 or universe_mask >= (1 << k):
            raise ValueError("structural ALT universe lies outside k")
        if universe_mask != observed_alt_mask:
            raise ValueError(
                "structural ALT universe differs from observed ALT union"
            )

    active_bits = tuple(
        bit for bit in range(k) if universe_mask & (1 << bit)
    )
    effective_m = len(active_bits)
    if effective_m > max_m:
        raise ValueError(
            f"effective_m={effective_m} exceeds exact counter gate {max_m}"
        )
    compact_index = {
        original_bit: compact_bit
        for compact_bit, original_bit in enumerate(active_bits)
    }
    adjacency = [0 for _ in range(effective_m)]
    forbidden = [0 for _ in range(effective_m)]
    mandatory_nonroot_states: set[int] = set()

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
        if pattern in full_patterns:
            mandatory_mask = 0
            for mutation in alt:
                mandatory_mask |= 1 << mutation
            if mandatory_mask:
                mandatory_nonroot_states.add(mandatory_mask)
        for left_index, left in enumerate(alt):
            for right in alt[left_index + 1 :]:
                adjacency[left] |= 1 << right
                adjacency[right] |= 1 << left
        for ancestor in ref:
            for descendant in alt:
                forbidden[ancestor] |= 1 << descendant

    state_count = 1 << effective_m
    tree_counts = [0 for _ in range(state_count)]
    forest_counts = [0 for _ in range(state_count)]
    forest_counts[0] = 1

    # Bottom-up order guarantees that every forest term used by T[S] is
    # smaller than S, while T[S] is available before W[S] uses its one-block
    # partition.
    for size in range(1, effective_m + 1):
        for mask in range(1, state_count):
            if _popcount(mask) != size:
                continue

            roots = mask
            while roots:
                root_bit = roots & -roots
                root = root_bit.bit_length() - 1
                descendants = mask ^ root_bit
                if not (forbidden[root] & descendants):
                    tree_counts[mask] += forest_counts[descendants]
                roots ^= root_bit

            components = _induced_components(mask, adjacency)
            first_component = components[0]
            other_components = components[1:]
            for choice in range(1 << len(other_components)):
                block = first_component
                for index, component in enumerate(other_components):
                    if choice & (1 << index):
                        block |= component
                forest_counts[mask] += (
                    tree_counts[block] * forest_counts[mask ^ block]
                )

    family_count = forest_counts[state_count - 1]
    exact = family_count > 0
    mandatory_count = len(mandatory_nonroot_states)
    if exact:
        status = (
            EXACT_STATUS
            if mandatory_count == 0
            else EXACT_MANDATORY_STATUS
        )
    else:
        status = ABSTAIN_STATUS
    objective = effective_m - mandatory_count if exact else None
    return PerfectFamilyCountResult(
        scope=SCOPE,
        status=status,
        exact=exact,
        ranking_allowed=False,
        k=k,
        effective_m=effective_m,
        mandatory_nonroot_state_count=mandatory_count,
        objective=objective,
        perfect_family_count=family_count,
        comparable_pair_count=sum(_popcount(row) for row in adjacency) // 2,
        forbidden_ancestor_pair_count=sum(
            _popcount(row) for row in forbidden
        ),
        subset_count=state_count,
        elapsed_seconds=time.perf_counter() - started,
        proof_basis=(
            "nonzero recurrence-free construction gives h<=m-d; every active "
            "ALT requires an acquisition edge and d mandatory non-root states "
            "are objective-free, so h>=m-d"
            if exact
            else "no recurrence-free h=m family; recurrence-allowed counting "
            "is outside this fast path"
        ),
    )


def count_manifest_case(
    case: Mapping[str, Any],
    *,
    max_m: int = 20,
) -> PerfectFamilyCountResult:
    """Count one frozen stress-panel manifest case."""
    structural = case.get("structural_input")
    if not isinstance(structural, Mapping):
        raise ValueError("manifest case lacks structural_input")
    return count_perfect_families(
        full_patterns=tuple(structural.get("full_patterns", ())),
        partial_patterns=tuple(structural.get("partial_patterns", ())),
        k=int(structural["k"]),
        structural_alt_universe_mask=int(
            structural["structural_alt_universe_mask"]
        ),
        max_m=max_m,
    )


def assert_finite_elapsed(result: PerfectFamilyCountResult) -> None:
    """Fail closed if a timing field cannot be serialized as finite evidence."""
    if result.elapsed_seconds < 0 or not math.isfinite(result.elapsed_seconds):
        raise AssertionError("counter elapsed time is not finite and nonnegative")
