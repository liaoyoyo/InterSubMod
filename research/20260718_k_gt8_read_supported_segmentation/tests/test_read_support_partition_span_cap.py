from decimal import Decimal
from itertools import combinations
from pathlib import Path
import random
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import read_support_partition as k_only  # noqa: E402
import read_support_partition_span_cap as span_cap  # noqa: E402


def _disposition_map(result):
    return {row.constraint_id: row for row in result.dispositions}


def _brute_force_cut_indices(
    positions,
    constraints,
    *,
    max_block_size,
    max_block_span_bp,
):
    """Return the exact optimum by enumerating every possible cut subset."""

    position_index = {
        position: index for index, position in enumerate(positions)
    }
    best_score = None
    best_cuts = None
    for n_cuts in range(len(positions)):
        for cuts in combinations(range(1, len(positions)), n_cuts):
            boundaries = (0,) + cuts + (len(positions),)
            blocks = tuple(zip(boundaries, boundaries[1:]))
            if any(
                end - start > max_block_size for start, end in blocks
            ):
                continue
            if max_block_span_bp is not None and any(
                positions[end - 1] - positions[start] > max_block_span_bp
                for start, end in blocks
            ):
                continue

            retained_weight = Decimal(0)
            retained_patterns = 0
            for edge in constraints:
                indices = tuple(position_index[site] for site in edge.sites)
                lo = min(indices)
                hi = max(indices)
                if any(start <= lo and hi < end for start, end in blocks):
                    retained_weight += edge.weight
                    retained_patterns += edge.pattern_count
            cut_gap_sum = sum(
                positions[cut] - positions[cut - 1] for cut in cuts
            )
            score = (
                retained_weight,
                retained_patterns,
                -len(blocks),
                cut_gap_sum,
            )
            if (
                best_score is None
                or score > best_score
                or (score == best_score and cuts < best_cuts)
            ):
                best_score = score
                best_cuts = cuts
    assert best_cuts is not None
    return best_cuts


def test_large_gap_forces_cut_and_charges_bridge_once():
    positions = (100, 110, 120, 10_000, 10_010)
    constraints = (
        span_cap.Hyperedge("left", (100, 120), weight=11),
        span_cap.Hyperedge("bridge", (120, 10_000), weight=5),
        span_cap.Hyperedge("right", (10_000, 10_010), weight=13),
    )

    result = span_cap.partition_ordered_hypergraph(
        positions,
        constraints,
        max_block_size=8,
        max_block_span_bp=30,
    )

    assert result.cut_indices == (3,)
    assert [block.positions for block in result.blocks] == [
        (100, 110, 120),
        (10_000, 10_010),
    ]
    assert all(block.span_bp <= 30 for block in result.blocks)
    assert result.retained_weight == Decimal(24)
    assert result.lost_weight == Decimal(5)
    bridge = _disposition_map(result)["bridge"]
    assert bridge.status == span_cap.UNAVOIDABLE_SPAN_CAP
    assert bridge.crossed_cut_count == 1
    assert result.unavoidable_constraint_ids == ("bridge",)


def test_singletons_remain_feasible_at_zero_bp_cap():
    positions = (100, 1_000, 2_000)
    constraints = (
        span_cap.Hyperedge("one", (100,), weight=2),
        span_cap.Hyperedge("two", (1_000,), weight=3),
        span_cap.Hyperedge("long", (100, 2_000), weight=7),
    )

    result = span_cap.partition_ordered_hypergraph(
        positions,
        constraints,
        max_block_size=8,
        max_block_span_bp=0,
    )

    assert result.cut_indices == (1, 2)
    assert [block.k for block in result.blocks] == [1, 1, 1]
    assert [block.span_bp for block in result.blocks] == [0, 0, 0]
    assert result.retained_constraint_ids == ("one", "two")
    assert result.lost_weight == Decimal(7)
    assert (
        _disposition_map(result)["long"].status
        == span_cap.UNAVOIDABLE_SPAN_CAP
    )


def test_span_cap_is_inclusive():
    result = span_cap.partition_ordered_hypergraph(
        (100, 200),
        (span_cap.Hyperedge("pair", (100, 200), weight=4),),
        max_block_size=8,
        max_block_span_bp=100,
    )

    assert result.cut_indices == ()
    assert result.retained_constraint_ids == ("pair",)


def test_span_cap_matches_bruteforce_oracle():
    generator = random.Random(20260718)
    for case_index in range(250):
        n_sites = generator.randint(1, 7)
        positions = [generator.randint(1, 10)]
        for _ in range(n_sites - 1):
            positions.append(positions[-1] + generator.randint(1, 25))
        positions = tuple(positions)
        max_block_size = generator.randint(1, min(4, n_sites))
        max_block_span_bp = generator.choice(
            (None, 0, 5, 10, 20, 40, 80)
        )

        constraints = []
        for edge_index in range(generator.randint(0, 9)):
            edge_size = generator.randint(1, n_sites)
            edge_sites = tuple(
                sorted(generator.sample(positions, edge_size))
            )
            constraints.append(
                span_cap.Hyperedge(
                    f"case_{case_index}_edge_{edge_index}",
                    edge_sites,
                    weight=generator.randint(0, 7),
                    pattern_count=generator.randint(1, 4),
                )
            )
        constraints = tuple(constraints)

        expected_cuts = _brute_force_cut_indices(
            positions,
            constraints,
            max_block_size=max_block_size,
            max_block_span_bp=max_block_span_bp,
        )
        observed = span_cap.partition_ordered_hypergraph(
            positions,
            constraints,
            max_block_size=max_block_size,
            max_block_span_bp=max_block_span_bp,
        )

        assert observed.cut_indices == expected_cuts
        assert all(block.k <= max_block_size for block in observed.blocks)
        if max_block_span_bp is not None:
            assert all(
                block.span_bp <= max_block_span_bp
                for block in observed.blocks
            )
        assert (
            observed.retained_weight + observed.lost_weight
            == observed.total_weight
        )
        assert (
            set(observed.retained_constraint_ids)
            | set(observed.lost_constraint_ids)
            == {edge.constraint_id for edge in constraints}
        )


def _normalized_result(result):
    return {
        "positions": result.positions,
        "max_block_size": result.max_block_size,
        "blocks": tuple(
            (
                block.block_index,
                block.start_index,
                block.end_index_exclusive,
                block.positions,
                block.retained_constraint_ids,
                block.retained_weight,
                block.retained_pattern_count,
            )
            for block in result.blocks
        ),
        "cut_indices": result.cut_indices,
        "cut_gaps": result.cut_gaps,
        "retained_constraint_ids": result.retained_constraint_ids,
        "lost_constraint_ids": result.lost_constraint_ids,
        "unavoidable_constraint_ids": result.unavoidable_constraint_ids,
        "total_weight": result.total_weight,
        "retained_weight": result.retained_weight,
        "lost_weight": result.lost_weight,
        "total_pattern_count": result.total_pattern_count,
        "retained_pattern_count": result.retained_pattern_count,
        "lost_pattern_count": result.lost_pattern_count,
        "dispositions": tuple(
            (
                row.constraint_id,
                row.status,
                row.span_sites,
                row.crossed_cut_count,
                row.block_index,
                row.weight,
                row.pattern_count,
            )
            for row in result.dispositions
        ),
    }


def test_none_cap_has_semantic_parity_with_current_k_only_dp():
    generator = random.Random(1870)
    for case_index in range(200):
        n_sites = generator.randint(1, 9)
        positions = [generator.randint(1, 10)]
        for _ in range(n_sites - 1):
            positions.append(positions[-1] + generator.randint(1, 100))
        positions = tuple(positions)
        max_block_size = generator.randint(1, min(8, n_sites))
        edge_specs = []
        for edge_index in range(generator.randint(0, 12)):
            edge_size = generator.randint(1, n_sites)
            edge_specs.append(
                (
                    f"case_{case_index}_edge_{edge_index}",
                    tuple(sorted(generator.sample(positions, edge_size))),
                    generator.randint(0, 9),
                    generator.randint(1, 4),
                )
            )

        old_edges = tuple(
            k_only.Hyperedge(
                edge_id,
                sites,
                weight=weight,
                pattern_count=pattern_count,
            )
            for edge_id, sites, weight, pattern_count in edge_specs
        )
        new_edges = tuple(
            span_cap.Hyperedge(
                edge_id,
                sites,
                weight=weight,
                pattern_count=pattern_count,
            )
            for edge_id, sites, weight, pattern_count in edge_specs
        )
        old_result = k_only.partition_ordered_hypergraph(
            positions,
            old_edges,
            max_block_size=max_block_size,
        )
        new_result = span_cap.partition_ordered_hypergraph(
            positions,
            new_edges,
            max_block_size=max_block_size,
            max_block_span_bp=None,
        )

        assert new_result.max_block_span_bp is None
        assert _normalized_result(new_result) == _normalized_result(old_result)


@pytest.mark.parametrize("cap", (-1, True, 1.5, "100"))
def test_invalid_span_cap_fails_closed(cap):
    with pytest.raises(
        ValueError,
        match="max_block_span_bp",
    ):
        span_cap.partition_ordered_hypergraph(
            (100,),
            (),
            max_block_span_bp=cap,
        )
