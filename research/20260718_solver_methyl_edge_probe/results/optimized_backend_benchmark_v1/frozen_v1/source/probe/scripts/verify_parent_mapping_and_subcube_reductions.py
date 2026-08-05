#!/usr/bin/env python3
"""Verify parent-map factorization and a safe downward-closure reduction."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path
from typing import Iterable, Sequence


REPO = Path(__file__).resolve().parents[3]
CURRENT_SOLVER = (
    REPO
    / "research/20260716_read_linked_hypercube_exact_likelihood_validation"
    / "scripts/hypercube_exact.py"
)
REAL_FIXTURES = (
    REPO
    / "research/20260718_solver_methyl_edge_probe"
    / "tests/fixtures/real_units.json"
)


def popcount(value: int) -> int:
    return bin(value).count("1")


def predecessors(vertex: int, k: int, selected: set[int]) -> list[int]:
    return [
        vertex ^ (1 << bit)
        for bit in range(k)
        if vertex & (1 << bit) and (vertex ^ (1 << bit)) in selected
    ]


def rooted(selected: set[int], k: int) -> bool:
    return all(predecessors(vertex, k, selected) for vertex in selected - {0})


def edge_score(parent: int, child: int) -> int:
    """Return a deterministic local score with intentional ties."""
    return ((parent * 17 + child * 13) % 7) - 3


def verify_parent_factorization(k_max: int = 3) -> dict[str, object]:
    checked = 0
    for k in range(1, k_max + 1):
        universe = list(range(1 << k))
        for mask in range(1 << (len(universe) - 1)):
            selected = {0} | {
                universe[index + 1]
                for index in range(len(universe) - 1)
                if mask >> index & 1
            }
            if not rooted(selected, k):
                continue
            children = sorted(selected - {0})
            domains = [predecessors(child, k, selected) for child in children]
            explicit_best = max(
                sum(edge_score(parent, child) for parent, child in zip(choice, children))
                for choice in itertools.product(*domains)
            )
            analytic_best = sum(
                max(edge_score(parent, child) for parent in domain)
                for child, domain in zip(children, domains)
            )
            if explicit_best != analytic_best:
                raise AssertionError(
                    f"parent factorization mismatch k={k} selected={sorted(selected)}"
                )
            checked += 1
    return {
        "status": "PASS",
        "k_max": k_max,
        "root_connected_sets": checked,
        "explicit_best_equals_parent_map_best": True,
    }


def verify_strict_diamond() -> dict[str, object]:
    selected = {0, 1, 2, 3}
    assignments = []
    for parent_of_three in predecessors(3, 2, selected):
        edges = [(0, 1), (0, 2), (parent_of_three, 3)]
        gain_counts = {
            str(bit): sum(
                1
                for parent, child in edges
                if not (parent >> bit) & 1 and (child >> bit) & 1
            )
            for bit in range(2)
        }
        assignments.append(
            {
                "parent_of_11": parent_of_three,
                "gain_counts": gain_counts,
                "strict": all(count <= 1 for count in gain_counts.values()),
            }
        )
    if any(item["strict"] for item in assignments):
        raise AssertionError("diamond unexpectedly has a strict mutation-once parent map")
    return {
        "status": "PASS",
        "selected": [0, 1, 2, 3],
        "parent_assignments": assignments,
        "strict_feasible": False,
    }


def downward(vertex: int, universe: Iterable[int]) -> set[int]:
    return {candidate for candidate in universe if candidate & ~vertex == 0}


def pattern_members(pattern: Sequence[str], k: int) -> set[int]:
    return {
        vertex
        for vertex in range(1 << k)
        if all(
            symbol == "X"
            or ((vertex >> index) & 1) == (symbol == "A")
            for index, symbol in enumerate(pattern)
        )
    }


def verify_downward_closure(k_max: int = 3) -> dict[str, object]:
    cases = 0
    optimal_sets = 0
    for k in range(1, k_max + 1):
        universe = set(range(1 << k))
        partial_groups = [
            pattern_members(pattern, k)
            for pattern in itertools.product("RAX", repeat=k)
            if "X" in pattern
        ]
        mandatory_families = [
            set(values)
            for size in range(3)
            for values in itertools.combinations(sorted(universe - {0}), size)
        ]
        group_families = (
            [tuple()]
            + [(group,) for group in partial_groups]
            + list(itertools.combinations_with_replacement(partial_groups, 2))
        )
        for full_states in mandatory_families:
            mandatory = {0} | full_states
            optional = sorted(universe - mandatory)
            for groups in group_families:
                best_cost = None
                best_sets: list[set[int]] = []
                for mask in range(1 << len(optional)):
                    selected = mandatory | {
                        vertex for index, vertex in enumerate(optional) if mask >> index & 1
                    }
                    if not rooted(selected, k):
                        continue
                    if any(not (selected & group) for group in groups):
                        continue
                    cost = len(selected - mandatory)
                    if best_cost is None or cost < best_cost:
                        best_cost = cost
                        best_sets = [selected]
                    elif cost == best_cost:
                        best_sets.append(selected)

                terminal_union = set(full_states)
                for group in groups:
                    terminal_union |= group
                closure = {0}
                for terminal in terminal_union:
                    closure |= downward(terminal, universe)
                for selected in best_sets:
                    if not selected <= closure:
                        raise AssertionError(
                            f"downward closure mismatch k={k} selected={sorted(selected)}"
                        )
                cases += 1
                optimal_sets += len(best_sets)
    return {
        "status": "PASS",
        "k_max": k_max,
        "cases": cases,
        "optimal_sets": optimal_sets,
        "all_minimum_sets_within_downward_closure": True,
    }


def current_solver_real_fixture_census() -> dict[str, object]:
    sys.path.insert(0, str(CURRENT_SOLVER.parent))
    import hypercube_exact  # pylint: disable=import-error,import-outside-toplevel

    units = json.loads(REAL_FIXTURES.read_text())["units"]
    rows = []
    for unit in units:
        patterns = list(unit.get("full_counts", {})) + list(unit.get("partial_counts", {}))
        active_mask = 0
        for pattern in patterns:
            for index, symbol in enumerate(pattern):
                if symbol == "A":
                    active_mask |= 1 << index
        groups = [
            hypercube_exact.SymbolicPattern.from_string(pattern)
            for pattern in unit.get("partial_counts", {})
        ]
        mandatory = {0}
        for pattern in unit.get("full_counts", {}):
            mandatory.add(
                hypercube_exact.SymbolicPattern.from_string(pattern).alt_mask & active_mask
            )
        reduced = hypercube_exact._reduce_partial_groups(  # pylint: disable=protected-access
            groups, active_mask, mandatory
        )

        active_bits = [
            index for index in range(unit["k"]) if active_mask & (1 << index)
        ]
        m = len(active_bits)
        universe = set(range(1 << m))

        def projected_members(pattern: str) -> set[int]:
            return {
                vertex
                for vertex in universe
                if all(
                    pattern[raw_index] == "X"
                    or ((vertex >> projected_index) & 1)
                    == (pattern[raw_index] == "A")
                    for projected_index, raw_index in enumerate(active_bits)
                )
            }

        terminal_union: set[int] = set()
        for pattern in patterns:
            terminal_union |= projected_members(pattern)
        closure = {0}
        for terminal in terminal_union:
            closure |= downward(terminal, universe)
        rows.append(
            {
                "case_id": unit["case_id"],
                "raw_k": unit["k"],
                "effective_m": m,
                "active_vertices": 1 << m,
                "input_groups": len(groups),
                "reduced_groups": len(reduced.groups),
                "forced_vertices": len(reduced.forced_vertices),
                "downward_closure_vertices": len(closure),
                "downward_closure_reduction_fraction": 1.0 - len(closure) / (1 << m),
            }
        )
    return {
        "status": "PASS",
        "current_solver_sha256": hashlib.sha256(CURRENT_SOLVER.read_bytes()).hexdigest(),
        "fixture_path": str(REAL_FIXTURES),
        "cases": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = {
        "schema": "intersubmod.parent_mapping_subcube_plan_oracle.v1",
        "scope": "bounded_exploratory_not_production_validation",
        "inputs": {
            "current_solver": str(CURRENT_SOLVER),
            "real_fixtures": str(REAL_FIXTURES),
        },
        "parent_factorization": verify_parent_factorization(),
        "strict_diamond": verify_strict_diamond(),
        "downward_closure": verify_downward_closure(),
        "real_fixture_census": current_solver_real_fixture_census(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        "PASS"
        f" parent_sets={receipt['parent_factorization']['root_connected_sets']}"
        f" closure_cases={receipt['downward_closure']['cases']}"
        f" closure_optimal_sets={receipt['downward_closure']['optimal_sets']}"
        f" real_cases={len(receipt['real_fixture_census']['cases'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
