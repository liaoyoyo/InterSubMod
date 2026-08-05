from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "tools" / "strict_endpoint_graph.py"
SPEC = importlib.util.spec_from_file_location("test_strict_endpoint_graph_core", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
GRAPH = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = GRAPH
SPEC.loader.exec_module(GRAPH)


def molecule(molecule_id: str, *calls: tuple[int, str]):
    return GRAPH.MoleculeCalls(molecule_id=molecule_id, calls=calls)


def support_by_pair(rows):
    return {(row.site_i, row.site_j): row for row in rows}


def test_ab_and_bc_transitively_connect_one_component() -> None:
    molecules = [molecule(f"ab-{index}", (0, "A"), (1, "R")) for index in range(3)]
    molecules += [molecule(f"bc-{index}", (1, "A"), (2, "R")) for index in range(3)]

    support = GRAPH.accumulate_pair_support(molecules)

    assert GRAPH.connected_components((0, 1, 2), support, threshold=3) == ((0, 1, 2),)
    assert (0, 2) not in support_by_pair(support)


def test_ac_endpoint_evidence_does_not_passively_connect_b() -> None:
    support = GRAPH.accumulate_pair_support(
        [molecule(f"ac-{index}", (10, "R"), (30, "A")) for index in range(3)]
    )

    assert GRAPH.connected_components((10, 20, 30), support, threshold=3) == (
        (10, 30),
        (20,),
    )


def test_threshold_two_edge_splits_at_threshold_three_and_partitions_are_nested() -> None:
    support = GRAPH.accumulate_pair_support(
        [
            molecule("m1", (0, "R"), (1, "R")),
            molecule("m2", (0, "A"), (1, "A")),
        ]
    )

    partitions = GRAPH.components_for_thresholds(
        (0, 1),
        support,
        thresholds=(3, 1, 2),
    )

    assert [partition.threshold for partition in partitions] == [1, 2, 3]
    assert partitions[0].components == ((0, 1),)
    assert partitions[1].components == ((0, 1),)
    assert partitions[2].components == ((0,), (1,))


def test_duplicate_molecule_id_fails_closed() -> None:
    repeated = molecule("same-molecule", (0, "R"), (1, "A"))

    with pytest.raises(ValueError, match="duplicate molecule_id"):
        GRAPH.accumulate_pair_support((repeated, repeated))


def test_more_than_50kb_distance_does_not_split_supported_endpoints() -> None:
    far_left = 100
    far_right = 100_101
    support = GRAPH.accumulate_pair_support(
        [molecule(f"long-{index}", (far_left, "R"), (far_right, "A")) for index in range(3)]
    )

    assert far_right - far_left > 50_000
    assert GRAPH.connected_components((far_left, far_right), support, threshold=3) == (
        (far_left, far_right),
    )


def test_interleaving_non_contiguous_components_are_preserved() -> None:
    molecules = [molecule(f"even-{index}", (0, "R"), (2, "A")) for index in range(3)]
    molecules += [molecule(f"odd-{index}", (1, "A"), (3, "R")) for index in range(3)]
    support = GRAPH.accumulate_pair_support(molecules)

    assert GRAPH.connected_components((0, 1, 2, 3), support, threshold=3) == (
        (0, 2),
        (1, 3),
    )


def test_pair_state_counts_conserve_total_and_use_sorted_endpoint_orientation() -> None:
    support = GRAPH.accumulate_pair_support(
        [
            molecule("rr", (4, "R"), (9, "R")),
            molecule("ra", (4, "R"), (9, "A")),
            molecule("ar", (4, "A"), (9, "R")),
            molecule("aa", (4, "A"), (9, "A")),
            molecule("aa-2", (4, "A"), (9, "A")),
        ]
    )

    assert len(support) == 1
    row = support[0]
    assert row.as_row() == {
        "site_i": 4,
        "site_j": 9,
        "total": 5,
        "RR": 1,
        "RA": 1,
        "AR": 1,
        "AA": 2,
    }
    assert row.total == sum(row.state_counts.values())


def test_digest_is_deterministic_under_molecule_and_site_input_order() -> None:
    molecules = (
        molecule("m1", (1, "R"), (3, "A")),
        molecule("m2", (1, "A"), (2, "A"), (3, "R")),
    )
    forward = GRAPH.accumulate_pair_support(molecules)
    reverse = GRAPH.accumulate_pair_support(reversed(molecules))

    assert forward == reverse
    assert GRAPH.strict_graph_digest((1, 2, 3), forward) == GRAPH.strict_graph_digest(
        (3, 1, 2),
        reversed(reverse),
    )


@pytest.mark.parametrize(
    ("calls", "message"),
    [
        (((1, "R"), (1, "A")), "strictly increasing"),
        (((2, "R"), (1, "A")), "strictly increasing"),
        (((1, "X"),), "code must be R or A"),
    ],
)
def test_molecule_calls_reject_non_contractual_calls(calls, message) -> None:
    with pytest.raises(ValueError, match=message):
        GRAPH.MoleculeCalls("invalid", calls)


def test_pair_support_constructor_enforces_state_conservation() -> None:
    with pytest.raises(ValueError, match="state conservation failed"):
        GRAPH.PairSupport(0, 1, total=3, rr=1, ra=1, ar=0, aa=0)


def test_streaming_accumulator_matches_batch_and_reports_mass() -> None:
    molecules = (
        molecule("m1", (1, "R"), (3, "A")),
        molecule("m2", (1, "A"), (2, "A"), (3, "R")),
    )
    accumulator = GRAPH.EndpointPairAccumulator()
    for row in molecules:
        accumulator.add(row)

    assert accumulator.molecule_count == 2
    assert accumulator.fixed_call_count == 5
    assert accumulator.site_indices == (1, 2, 3)
    assert accumulator.pair_support() == GRAPH.accumulate_pair_support(molecules)
