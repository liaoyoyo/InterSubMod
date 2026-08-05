#!/usr/bin/env python3
"""Exhaustive small-k oracle for the rooted three-gamete condition."""

from __future__ import annotations

import itertools
from collections import Counter


def predecessors(vertex: int, k: int, selected: set[int]) -> list[int]:
    return [
        vertex ^ (1 << bit)
        for bit in range(k)
        if vertex & (1 << bit) and vertex ^ (1 << bit) in selected
    ]


def is_root_connected(selected: set[int], k: int) -> bool:
    return 0 in selected and all(
        vertex == 0 or predecessors(vertex, k, selected)
        for vertex in selected
    )


def rooted_three_gamete_pass(selected: set[int], k: int) -> bool:
    for left in range(k):
        for right in range(left + 1, k):
            categories = {
                (
                    1 if vertex & (1 << left) else 0,
                    1 if vertex & (1 << right) else 0,
                )
                for vertex in selected
            }
            if {(1, 0), (0, 1), (1, 1)} <= categories:
                return False
    return True


def has_global_unique_acquisition_tree(selected: set[int], k: int) -> bool:
    ordered = sorted(selected - {0})
    choices = [predecessors(vertex, k, selected) for vertex in ordered]
    for parents in itertools.product(*choices):
        acquired_bits = []
        for vertex, parent in zip(ordered, parents):
            edge = vertex ^ parent
            acquired_bits.append(edge.bit_length() - 1)
        if all(count <= 1 for count in Counter(acquired_bits).values()):
            return True
    return False


def main() -> None:
    subsets_checked = 0
    connected_checked = 0
    for k in range(1, 4):
        nonroot = list(range(1, 1 << k))
        for mask in range(1 << len(nonroot)):
            selected = {0}
            selected.update(
                vertex
                for index, vertex in enumerate(nonroot)
                if mask & (1 << index)
            )
            subsets_checked += 1
            if not is_root_connected(selected, k):
                continue
            connected_checked += 1
            three_gamete = rooted_three_gamete_pass(selected, k)
            tree_exists = has_global_unique_acquisition_tree(selected, k)
            if three_gamete != tree_exists:
                raise AssertionError(
                    f"Oracle mismatch at k={k}, selected={sorted(selected)}, "
                    f"three_gamete={three_gamete}, tree_exists={tree_exists}"
                )
    print(
        f"PASS subsets={subsets_checked} root_connected={connected_checked} "
        "mismatches=0 k_max=3"
    )


if __name__ == "__main__":
    main()
