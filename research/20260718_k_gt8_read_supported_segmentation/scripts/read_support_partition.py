#!/usr/bin/env python3
"""Deterministic bounded partitioning for an ordered read-support hypergraph.

The ordered vertices are genomic sites.  A hyperedge is one read-pattern
constraint over one or more sites.  Blocks are contiguous, cover every site
exactly once, and contain at most ``max_block_size`` sites.

The dynamic program maximizes the weight of hyperedges that remain wholly
inside one block.  Consequently, a hyperedge crossing multiple cut boundaries
is still charged only once as a lost constraint.  Ties are resolved in this
fixed order:

1. greater retained weight;
2. greater retained pattern count;
3. fewer blocks;
4. greater sum of genomic gaps selected as cuts;
5. lexicographically smaller zero-based cut-index tuple.

An edge whose ordered interval spans more than ``max_block_size`` sites cannot
fit in any valid disjoint block and is reported as
``unavoidable_span_gt_max_block_size``.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
from decimal import Context, Decimal, InvalidOperation
from typing import Iterable, Optional, Sequence, Tuple


RETAINED = "retained"
CUT = "cut"
UNAVOIDABLE = "unavoidable_span_gt_max_block_size"


def _decimal_weight(value: object) -> Decimal:
    """Convert a numeric support weight to an exact, deterministic Decimal."""
    if isinstance(value, bool):
        raise ValueError("constraint weight must be numeric, not bool")
    try:
        converted = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, TypeError, ValueError) as exc:
        raise ValueError(f"invalid constraint weight: {value!r}") from exc
    if not converted.is_finite() or converted < 0:
        raise ValueError("constraint weight must be finite and non-negative")
    return converted


def _exact_weight_context(weights: Sequence[Decimal]) -> Context:
    """Return a local context that can exactly add every non-negative weight.

    Input weights are finite Decimals, but the process-wide Decimal context may
    still round their sum.  The DP and the final disposition audit aggregate
    constraints in different orders, so such rounding can make mathematically
    identical sums differ in their last digit.  Because all weights are
    non-negative, a precision spanning the smallest input exponent through the
    largest possible summed adjusted exponent is sufficient for every subset
    sum used by the algorithm.
    """

    nonzero = tuple(weight for weight in weights if weight)
    if not nonzero:
        return Context(prec=1)
    min_exponent = min(weight.as_tuple().exponent for weight in nonzero)
    max_adjusted = max(weight.adjusted() for weight in nonzero)
    carry_digits = len(str(len(nonzero)))
    precision = max_adjusted - min_exponent + carry_digits + 2
    return Context(prec=max(1, precision))


def _exact_weight_sum(
    weights: Iterable[Decimal], context: Context
) -> Decimal:
    """Sum finite weights exactly under a caller-sized local context."""

    total = Decimal(0)
    for weight in weights:
        total = context.add(total, weight)
    return total


@dataclass(frozen=True)
class Hyperedge:
    """One uniquely identified read-pattern constraint.

    ``sites`` contains genomic coordinates from the ordered vertex catalog.
    ``weight`` usually represents unique-molecule support.  ``pattern_count``
    is a separate secondary objective, useful when one row aggregates a
    read-pattern group.
    """

    constraint_id: str
    sites: Tuple[int, ...]
    weight: Decimal = Decimal(1)
    pattern_count: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.constraint_id, str) or not self.constraint_id:
            raise ValueError("constraint_id must be a non-empty string")
        raw_sites = tuple(self.sites)
        if not raw_sites:
            raise ValueError("hyperedge sites must not be empty")
        if any(isinstance(site, bool) or not isinstance(site, int) for site in raw_sites):
            raise ValueError("hyperedge sites must be integer coordinates")
        if len(raw_sites) != len(set(raw_sites)):
            raise ValueError("hyperedge sites must be unique")
        if isinstance(self.pattern_count, bool) or not isinstance(self.pattern_count, int):
            raise ValueError("pattern_count must be an integer")
        if self.pattern_count < 1:
            raise ValueError("pattern_count must be positive")
        object.__setattr__(self, "sites", tuple(sorted(raw_sites)))
        object.__setattr__(self, "weight", _decimal_weight(self.weight))


@dataclass(frozen=True)
class ConstraintDisposition:
    """Final, mutually exclusive disposition of one hyperedge."""

    constraint_id: str
    status: str
    span_sites: int
    crossed_cut_count: int
    block_index: Optional[int]
    weight: Decimal
    pattern_count: int


@dataclass(frozen=True)
class PartitionBlock:
    """One contiguous block in the final partition."""

    block_index: int
    start_index: int
    end_index_exclusive: int
    positions: Tuple[int, ...]
    retained_constraint_ids: Tuple[str, ...]
    retained_weight: Decimal
    retained_pattern_count: int

    @property
    def k(self) -> int:
        return len(self.positions)


@dataclass(frozen=True)
class PartitionResult:
    """Auditable result of the ordered-hypergraph dynamic program."""

    positions: Tuple[int, ...]
    max_block_size: int
    blocks: Tuple[PartitionBlock, ...]
    cut_indices: Tuple[int, ...]
    cut_gaps: Tuple[int, ...]
    retained_constraint_ids: Tuple[str, ...]
    lost_constraint_ids: Tuple[str, ...]
    unavoidable_constraint_ids: Tuple[str, ...]
    total_weight: Decimal
    retained_weight: Decimal
    lost_weight: Decimal
    total_pattern_count: int
    retained_pattern_count: int
    lost_pattern_count: int
    dispositions: Tuple[ConstraintDisposition, ...]

    @property
    def retention_ratio(self) -> Decimal:
        if self.total_weight == 0:
            return Decimal(1)
        return self.retained_weight / self.total_weight


@dataclass(frozen=True)
class _IndexedHyperedge:
    edge: Hyperedge
    lo: int
    hi: int

    @property
    def span_sites(self) -> int:
        return self.hi - self.lo + 1


@dataclass(frozen=True)
class _State:
    retained_weight: Decimal
    retained_pattern_count: int
    block_count: int
    cut_gap_sum: int
    cuts: Tuple[int, ...]


def _validate_positions(positions: Sequence[int]) -> Tuple[int, ...]:
    normalized = tuple(positions)
    if not normalized:
        raise ValueError("positions must not be empty")
    if any(isinstance(position, bool) or not isinstance(position, int) for position in normalized):
        raise ValueError("positions must be integer coordinates")
    if any(right <= left for left, right in zip(normalized, normalized[1:])):
        raise ValueError("positions must be strictly increasing and unique")
    return normalized


def _state_is_better(candidate: _State, incumbent: Optional[_State]) -> bool:
    if incumbent is None:
        return True
    candidate_key = (
        candidate.retained_weight,
        candidate.retained_pattern_count,
        -candidate.block_count,
        candidate.cut_gap_sum,
    )
    incumbent_key = (
        incumbent.retained_weight,
        incumbent.retained_pattern_count,
        -incumbent.block_count,
        incumbent.cut_gap_sum,
    )
    if candidate_key != incumbent_key:
        return candidate_key > incumbent_key
    return candidate.cuts < incumbent.cuts


def partition_ordered_hypergraph(
    positions: Sequence[int],
    constraints: Iterable[Hyperedge],
    *,
    max_block_size: int = 8,
) -> PartitionResult:
    """Partition ordered sites into optimal contiguous bounded blocks.

    A constraint is retained exactly when all of its sites occur in one final
    block.  The objective sums each retained constraint once; lost constraints
    are also classified once, independently of how many cuts they cross.
    """

    ordered_positions = _validate_positions(positions)
    if isinstance(max_block_size, bool) or not isinstance(max_block_size, int):
        raise ValueError("max_block_size must be an integer")
    if max_block_size < 1:
        raise ValueError("max_block_size must be positive")

    raw_edges = tuple(constraints)
    if any(not isinstance(edge, Hyperedge) for edge in raw_edges):
        raise ValueError("constraints must contain Hyperedge instances")
    edges = tuple(sorted(raw_edges, key=lambda edge: edge.constraint_id))
    ids = tuple(edge.constraint_id for edge in edges)
    if len(ids) != len(set(ids)):
        raise ValueError("constraint_id values must be unique")
    weight_context = _exact_weight_context(
        tuple(edge.weight for edge in edges)
    )

    position_index = {position: index for index, position in enumerate(ordered_positions)}
    indexed_edges = []
    interval_metrics = {}
    for edge in edges:
        unknown = tuple(site for site in edge.sites if site not in position_index)
        if unknown:
            raise ValueError(
                f"constraint {edge.constraint_id!r} references sites outside positions: {unknown}"
            )
        indices = tuple(position_index[site] for site in edge.sites)
        indexed = _IndexedHyperedge(edge=edge, lo=min(indices), hi=max(indices))
        indexed_edges.append(indexed)
        if indexed.span_sites <= max_block_size:
            key = (indexed.lo, indexed.hi)
            prior_weight, prior_patterns = interval_metrics.get(
                key, (Decimal(0), 0)
            )
            interval_metrics[key] = (
                weight_context.add(prior_weight, edge.weight),
                prior_patterns + edge.pattern_count,
            )

    # Only O(n * K^2) gain cells are needed because blocks cannot exceed K.
    block_gain = {}
    n_sites = len(ordered_positions)
    for start in range(n_sites):
        retained_weight = Decimal(0)
        retained_patterns = 0
        for end_inclusive in range(start, min(n_sites, start + max_block_size)):
            for lo in range(start, end_inclusive + 1):
                weight, patterns = interval_metrics.get(
                    (lo, end_inclusive), (Decimal(0), 0)
                )
                retained_weight = weight_context.add(
                    retained_weight, weight
                )
                retained_patterns += patterns
            block_gain[(start, end_inclusive + 1)] = (
                retained_weight,
                retained_patterns,
            )

    states = [None] * (n_sites + 1)
    states[0] = _State(Decimal(0), 0, 0, 0, ())
    for end_exclusive in range(1, n_sites + 1):
        best = None
        first_start = max(0, end_exclusive - max_block_size)
        for start in range(first_start, end_exclusive):
            previous = states[start]
            if previous is None:
                raise AssertionError("dynamic-program prefix state is missing")
            gain_weight, gain_patterns = block_gain[(start, end_exclusive)]
            if start:
                cuts = previous.cuts + (start,)
                gap_sum = (
                    previous.cut_gap_sum
                    + ordered_positions[start]
                    - ordered_positions[start - 1]
                )
            else:
                cuts = previous.cuts
                gap_sum = previous.cut_gap_sum
            candidate = _State(
                retained_weight=weight_context.add(
                    previous.retained_weight, gain_weight
                ),
                retained_pattern_count=previous.retained_pattern_count
                + gain_patterns,
                block_count=previous.block_count + 1,
                cut_gap_sum=gap_sum,
                cuts=cuts,
            )
            if _state_is_better(candidate, best):
                best = candidate
        states[end_exclusive] = best

    final_state = states[n_sites]
    if final_state is None:
        raise AssertionError("dynamic program produced no terminal state")
    cut_indices = final_state.cuts
    boundaries = (0,) + cut_indices + (n_sites,)
    raw_blocks = tuple(zip(boundaries, boundaries[1:]))

    retained_by_block = [[] for _ in raw_blocks]
    dispositions = []
    retained_ids = []
    lost_ids = []
    unavoidable_ids = []
    retained_weight_check = Decimal(0)
    retained_patterns_check = 0

    for indexed in indexed_edges:
        crossed_cut_count = (
            bisect_right(cut_indices, indexed.hi)
            - bisect_right(cut_indices, indexed.lo)
        )
        candidate_block = bisect_right(cut_indices, indexed.lo)
        block_start, block_end_exclusive = raw_blocks[candidate_block]
        block_index = (
            candidate_block
            if block_start <= indexed.lo and indexed.hi < block_end_exclusive
            else None
        )
        if block_index is not None:
            status = RETAINED
            retained_by_block[block_index].append(indexed.edge.constraint_id)
            retained_ids.append(indexed.edge.constraint_id)
            retained_weight_check = weight_context.add(
                retained_weight_check, indexed.edge.weight
            )
            retained_patterns_check += indexed.edge.pattern_count
        else:
            if crossed_cut_count == 0:
                raise AssertionError("lost constraint does not cross a cut")
            lost_ids.append(indexed.edge.constraint_id)
            if indexed.span_sites > max_block_size:
                status = UNAVOIDABLE
                unavoidable_ids.append(indexed.edge.constraint_id)
            else:
                status = CUT
        dispositions.append(
            ConstraintDisposition(
                constraint_id=indexed.edge.constraint_id,
                status=status,
                span_sites=indexed.span_sites,
                crossed_cut_count=crossed_cut_count,
                block_index=block_index,
                weight=indexed.edge.weight,
                pattern_count=indexed.edge.pattern_count,
            )
        )

    if retained_weight_check != final_state.retained_weight:
        raise AssertionError("retained-weight objective does not match dispositions")
    if retained_patterns_check != final_state.retained_pattern_count:
        raise AssertionError("retained-pattern objective does not match dispositions")

    edge_by_id = {edge.constraint_id: edge for edge in edges}
    blocks = []
    for block_index, (start, end_exclusive) in enumerate(raw_blocks):
        block_ids = tuple(sorted(retained_by_block[block_index]))
        blocks.append(
            PartitionBlock(
                block_index=block_index,
                start_index=start,
                end_index_exclusive=end_exclusive,
                positions=ordered_positions[start:end_exclusive],
                retained_constraint_ids=block_ids,
                retained_weight=_exact_weight_sum(
                    (edge_by_id[constraint_id].weight for constraint_id in block_ids),
                    weight_context,
                ),
                retained_pattern_count=sum(
                    edge_by_id[constraint_id].pattern_count
                    for constraint_id in block_ids
                ),
            )
        )

    if tuple(position for block in blocks for position in block.positions) != ordered_positions:
        raise AssertionError("site partition does not conserve ordered positions")
    if any(block.k > max_block_size for block in blocks):
        raise AssertionError("block-size contract was violated")

    total_weight = _exact_weight_sum(
        (edge.weight for edge in edges), weight_context
    )
    total_patterns = sum(edge.pattern_count for edge in edges)
    retained_ids_tuple = tuple(sorted(retained_ids))
    lost_ids_tuple = tuple(sorted(lost_ids))
    if set(retained_ids_tuple) & set(lost_ids_tuple):
        raise AssertionError("constraint appears in both retained and lost sets")
    if set(retained_ids_tuple) | set(lost_ids_tuple) != set(ids):
        raise AssertionError("constraint disposition does not conserve IDs")

    return PartitionResult(
        positions=ordered_positions,
        max_block_size=max_block_size,
        blocks=tuple(blocks),
        cut_indices=cut_indices,
        cut_gaps=tuple(
            ordered_positions[cut] - ordered_positions[cut - 1]
            for cut in cut_indices
        ),
        retained_constraint_ids=retained_ids_tuple,
        lost_constraint_ids=lost_ids_tuple,
        unavoidable_constraint_ids=tuple(sorted(unavoidable_ids)),
        total_weight=total_weight,
        retained_weight=retained_weight_check,
        lost_weight=weight_context.subtract(
            total_weight, retained_weight_check
        ),
        total_pattern_count=total_patterns,
        retained_pattern_count=retained_patterns_check,
        lost_pattern_count=total_patterns - retained_patterns_check,
        dispositions=tuple(dispositions),
    )
