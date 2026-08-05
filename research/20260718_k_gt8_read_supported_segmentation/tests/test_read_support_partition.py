from decimal import Decimal, localcontext
from pathlib import Path
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from read_support_partition import (  # noqa: E402
    CUT,
    RETAINED,
    UNAVOIDABLE,
    Hyperedge,
    partition_ordered_hypergraph,
)


def disposition_map(result):
    return {row.constraint_id: row for row in result.dispositions}


def test_weak_bridge_selects_six_plus_six_and_retains_strong_groups():
    positions = tuple(range(1, 13))
    constraints = (
        Hyperedge("left", tuple(range(1, 7)), weight=20),
        Hyperedge("weak_bridge", (6, 7), weight=3),
        Hyperedge("right", tuple(range(7, 13)), weight=18),
    )

    result = partition_ordered_hypergraph(
        positions, constraints, max_block_size=8
    )

    assert result.cut_indices == (6,)
    assert [block.positions for block in result.blocks] == [
        tuple(range(1, 7)),
        tuple(range(7, 13)),
    ]
    assert result.retained_weight == Decimal(38)
    assert result.lost_weight == Decimal(3)
    assert result.retention_ratio == Decimal(38) / Decimal(41)
    assert result.retained_constraint_ids == ("left", "right")
    assert result.lost_constraint_ids == ("weak_bridge",)
    assert disposition_map(result)["weak_bridge"].status == CUT


def test_one_long_constraint_is_lost_once_even_when_crossing_multiple_cuts():
    positions = tuple(range(1, 25))
    constraints = (Hyperedge("long_read", (1, 24), weight=5),)

    result = partition_ordered_hypergraph(
        positions, constraints, max_block_size=8
    )

    assert result.cut_indices == (8, 16)
    assert [block.k for block in result.blocks] == [8, 8, 8]
    assert result.total_weight == Decimal(5)
    assert result.retained_weight == Decimal(0)
    assert result.lost_weight == Decimal(5)
    assert result.lost_constraint_ids == ("long_read",)
    assert result.unavoidable_constraint_ids == ("long_read",)
    row = disposition_map(result)["long_read"]
    assert row.status == UNAVOIDABLE
    assert row.crossed_cut_count == 2
    assert row.span_sites == 24


def test_single_hyperedge_spanning_nine_ordered_sites_is_unavoidable():
    positions = tuple(range(10, 19))
    result = partition_ordered_hypergraph(
        positions,
        (Hyperedge("nine_site_span", (10, 18), weight=7),),
        max_block_size=8,
    )

    assert len(result.blocks) == 2
    assert result.lost_weight == Decimal(7)
    assert result.unavoidable_constraint_ids == ("nine_site_span",)
    assert disposition_map(result)["nine_site_span"].status == UNAVOIDABLE


def test_k_at_most_limit_is_identity_partition_and_retains_everything():
    positions = (100, 200, 300, 400, 500, 600, 700, 800)
    constraints = (
        Hyperedge("ends", (100, 800), weight=9, pattern_count=2),
        Hyperedge("middle", (300, 400, 500), weight=4, pattern_count=3),
    )

    result = partition_ordered_hypergraph(positions, constraints)

    assert len(result.blocks) == 1
    assert result.blocks[0].positions == positions
    assert result.cut_indices == ()
    assert result.retained_weight == Decimal(13)
    assert result.retained_pattern_count == 5
    assert result.lost_constraint_ids == ()
    assert all(row.status == RETAINED for row in result.dispositions)


def test_retained_weight_precedes_retained_pattern_count():
    positions = (1, 2, 3)
    constraints = (
        Hyperedge("heavier", (1, 2), weight=6, pattern_count=1),
        Hyperedge("more_patterns", (2, 3), weight=5, pattern_count=100),
    )

    result = partition_ordered_hypergraph(
        positions, constraints, max_block_size=2
    )

    assert result.cut_indices == (2,)
    assert result.retained_constraint_ids == ("heavier",)


def test_pattern_count_breaks_equal_weight_tie():
    positions = (1, 2, 3)
    constraints = (
        Hyperedge("one_pattern", (1, 2), weight=5, pattern_count=1),
        Hyperedge("two_patterns", (2, 3), weight=5, pattern_count=2),
    )

    result = partition_ordered_hypergraph(
        positions, constraints, max_block_size=2
    )

    assert result.cut_indices == (1,)
    assert result.retained_constraint_ids == ("two_patterns",)
    assert result.retained_pattern_count == 2


def test_fewer_blocks_then_larger_cut_gap_then_lexicographic_cut():
    gap_positions = (1, 2, 100, 101, 102, 103)
    gap_result = partition_ordered_hypergraph(
        gap_positions, (), max_block_size=4
    )

    assert len(gap_result.blocks) == 2
    assert gap_result.cut_indices == (2,)
    assert gap_result.cut_gaps == (98,)

    equal_gap_positions = tuple(range(1, 7))
    lex_result = partition_ordered_hypergraph(
        equal_gap_positions, (), max_block_size=4
    )

    assert len(lex_result.blocks) == 2
    assert lex_result.cut_indices == (2,)
    assert lex_result.cut_gaps == (1,)


def test_zero_loss_multi_block_partition_conserves_every_constraint():
    positions = tuple(range(1, 13))
    constraints = (
        Hyperedge("cluster_a", (1, 2, 3, 4), weight=10),
        Hyperedge("cluster_b", (5, 6, 7, 8), weight=11),
        Hyperedge("cluster_c", (9, 10, 11, 12), weight=12),
    )

    result = partition_ordered_hypergraph(
        positions, constraints, max_block_size=4
    )

    assert result.cut_indices == (4, 8)
    assert result.retained_weight == result.total_weight == Decimal(33)
    assert result.lost_weight == Decimal(0)
    assert result.lost_constraint_ids == ()
    assert all(block.k <= 4 for block in result.blocks)


def test_constraint_order_does_not_change_partition_or_dispositions():
    positions = tuple(range(1, 13))
    constraints = (
        Hyperedge("left", tuple(range(1, 7)), weight=20),
        Hyperedge("weak", (6, 7), weight=3),
        Hyperedge("right", tuple(range(7, 13)), weight=18),
    )

    forward = partition_ordered_hypergraph(
        positions, constraints, max_block_size=8
    )
    reverse = partition_ordered_hypergraph(
        positions, reversed(constraints), max_block_size=8
    )

    assert reverse.cut_indices == forward.cut_indices
    assert reverse.retained_constraint_ids == forward.retained_constraint_ids
    assert reverse.lost_constraint_ids == forward.lost_constraint_ids
    assert reverse.dispositions == forward.dispositions


def test_finite_log_weights_are_summed_exactly_despite_aggregation_order():
    with localcontext() as context:
        context.prec = 28
        ln_2 = Decimal(2).ln()
        ln_6 = Decimal(6).ln()

    # Constraint-ID order is ln(2), ln(2), ln(6), while block aggregation
    # groups ln(2) + ln(6) at position 1 before adding position 2.  At the
    # default precision those two mathematically equivalent orders differ by
    # 1E-27 and previously tripped the disposition audit.
    constraints = (
        Hyperedge("a_ln2_pos1", (1,), weight=ln_2),
        Hyperedge("b_ln2_pos2", (2,), weight=ln_2),
        Hyperedge("c_ln6_pos1", (1,), weight=ln_6),
    )

    result = partition_ordered_hypergraph(
        (1, 2), constraints, max_block_size=2
    )

    expected = Decimal("3.1780538303479456196469416010")
    assert result.cut_indices == ()
    assert result.total_weight == expected
    assert result.retained_weight == expected
    assert result.blocks[0].retained_weight == expected
    assert result.lost_weight == Decimal(0)
    assert all(row.status == RETAINED for row in result.dispositions)


@pytest.mark.parametrize(
    "positions,constraints,error",
    [
        ((2, 1), (), "strictly increasing"),
        ((1, 1), (), "strictly increasing"),
        ((1, 2), (Hyperedge("outside", (1, 3)),), "outside positions"),
        (
            (1, 2),
            (Hyperedge("duplicate", (1,)), Hyperedge("duplicate", (2,))),
            "must be unique",
        ),
    ],
)
def test_invalid_catalog_or_constraint_fails_closed(
    positions, constraints, error
):
    with pytest.raises(ValueError, match=error):
        partition_ordered_hypergraph(positions, constraints)


def test_hyperedge_input_contract_fails_closed():
    with pytest.raises(ValueError, match="must not be empty"):
        Hyperedge("empty", ())
    with pytest.raises(ValueError, match="must be unique"):
        Hyperedge("duplicate_sites", (1, 1))
    with pytest.raises(ValueError, match="non-negative"):
        Hyperedge("negative", (1,), weight=-1)
    with pytest.raises(ValueError, match="positive"):
        Hyperedge("no_pattern", (1,), pattern_count=0)
    with pytest.raises(ValueError, match="Hyperedge instances"):
        partition_ordered_hypergraph((1, 2), ("not-a-hyperedge",))
