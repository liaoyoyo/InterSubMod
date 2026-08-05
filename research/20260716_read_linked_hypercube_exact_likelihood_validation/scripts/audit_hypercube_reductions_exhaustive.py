#!/usr/bin/env python3
"""Deterministic exhaustive audit for Hypercube group presolve and no-good cuts.

The reducer under test is imported from ``hypercube_exact``.  The original
group predicate is reconstructed directly from the R/A/X characters so a bug
in ``SymbolicPattern.compatible`` cannot cancel a reducer bug.  Runtime is a
diagnostic only; deterministic counts and zero mismatches are the pass gates.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
import time
from pathlib import Path
from typing import Any

from hypercube_exact import SymbolicPattern, _reduce_partial_groups, submasks


EXPECTED = {
    "presolve_cases": 61_340,
    "selected_set_predicate_checks": 1_979_356,
    "sparse_no_good_pairs": 23_909,
    "dense_no_good_pairs": 21_844,
    "removal_events": {
        "duplicate": 19_918,
        "dominated": 1_554,
        "mandatory_or_forced_hit": 112_111,
        "singleton_forced": 4_848,
    },
}


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def direct_pattern_matches(pattern: str, vertex: int, universe_mask: int) -> bool:
    """Independent character-level group membership predicate."""
    if vertex < 0 or vertex & ~universe_mask:
        return False
    for bit, symbol in enumerate(pattern):
        state_alt = bool(vertex & (1 << bit))
        if symbol == "R" and state_alt:
            return False
        if symbol == "A" and not state_alt:
            return False
        if symbol not in {"R", "A", "X"}:
            raise ValueError(f"invalid pattern symbol: {symbol!r}")
    return True


def audit_presolve() -> dict[str, Any]:
    case_count = 0
    predicate_checks = 0
    mismatches = 0
    removals = {
        "duplicate": 0,
        "dominated": 0,
        "mandatory_or_forced_hit": 0,
        "singleton_forced": 0,
    }
    for k in range(1, 4):
        all_partial = [
            SymbolicPattern.from_string("".join(symbols))
            for symbols in itertools.product("RAX", repeat=k)
            if "X" in symbols
        ]
        for universe_mask in range(1 << k):
            vertices = submasks(universe_mask)
            index = {vertex: offset for offset, vertex in enumerate(vertices)}
            valid_groups = [
                group for group in all_partial
                if not (group.alt_mask & ~universe_mask)
            ]
            nonroot = [vertex for vertex in vertices if vertex]
            group_multisets = [
                groups
                for n_groups in range(4)
                for groups in itertools.combinations_with_replacement(
                    valid_groups, n_groups
                )
            ]
            direct_group_masks = {
                group.pattern: sum(
                    1 << index[vertex]
                    for vertex in vertices
                    if direct_pattern_matches(group.pattern, vertex, universe_mask)
                )
                for group in valid_groups
            }
            for n_extra_mandatory in range(min(2, len(nonroot)) + 1):
                for extra_mandatory in itertools.combinations(
                    nonroot, n_extra_mandatory
                ):
                    mandatory = {0, *extra_mandatory}
                    mandatory_bits = sum(
                        1 << index[vertex] for vertex in mandatory
                    )
                    optional = [
                        vertex for vertex in vertices if vertex not in mandatory
                    ]
                    selected_bitsets = [
                        mandatory_bits
                        | sum(
                            1 << index[vertex]
                            for offset, vertex in enumerate(optional)
                            if selected_mask & (1 << offset)
                        )
                        for selected_mask in range(1 << len(optional))
                    ]
                    for groups in group_multisets:
                        reduced = _reduce_partial_groups(
                            groups, universe_mask, mandatory
                        )
                        case_count += 1
                        removals["duplicate"] += reduced.n_duplicate_removed
                        removals["dominated"] += reduced.n_dominated_removed
                        removals["mandatory_or_forced_hit"] += (
                            reduced.n_required_hit_removed
                        )
                        removals["singleton_forced"] += reduced.n_forced_removed
                        original_masks = [
                            direct_group_masks[group.pattern] for group in groups
                        ]
                        reduced_masks = [
                            direct_group_masks[group.pattern]
                            for group in reduced.groups
                        ]
                        forced_bits = sum(
                            1 << index[vertex]
                            for vertex in reduced.forced_vertices
                        )
                        for selected in selected_bitsets:
                            original_feasible = all(
                                selected & group_mask
                                for group_mask in original_masks
                            )
                            reduced_feasible = (
                                (selected & forced_bits) == forced_bits
                                and all(
                                    selected & group_mask
                                    for group_mask in reduced_masks
                                )
                            )
                            predicate_checks += 1
                            mismatches += original_feasible != reduced_feasible
    return {
        "presolve_cases": case_count,
        "selected_set_predicate_checks": predicate_checks,
        "mismatches": mismatches,
        "removal_events": removals,
    }


def audit_sparse_no_good() -> dict[str, int]:
    pairs = 0
    mismatches = 0
    for n_variables in range(1, 9):
        for n_mandatory in range(min(3, n_variables) + 1):
            mandatory = set(range(n_mandatory))
            optional = tuple(range(n_mandatory, n_variables))
            for hidden in range(len(optional) + 1):
                sets = [
                    mandatory | set(values)
                    for values in itertools.combinations(optional, hidden)
                ]
                for excluded in sets:
                    for candidate in sets:
                        cut_accepts = (
                            sum(
                                vertex in candidate
                                for vertex in excluded - mandatory
                            )
                            <= hidden - 1
                        )
                        expected_accepts = candidate != excluded
                        pairs += 1
                        mismatches += cut_accepts != expected_accepts
    return {"pairs": pairs, "mismatches": mismatches}


def audit_dense_no_good() -> dict[str, int]:
    pairs = 0
    mismatches = 0
    for n_variables in range(1, 8):
        sets = [
            {
                offset
                for offset in range(n_variables)
                if selected_mask & (1 << offset)
            }
            for selected_mask in range(1 << n_variables)
        ]
        for excluded in sets:
            for candidate in sets:
                lhs = sum(
                    offset in candidate
                    for offset in range(n_variables)
                    if offset not in excluded
                ) - sum(offset in candidate for offset in excluded)
                cut_accepts = lhs >= 1 - len(excluded)
                expected_accepts = candidate != excluded
                pairs += 1
                mismatches += cut_accepts != expected_accepts
    return {"pairs": pairs, "mismatches": mismatches}


def build_receipt() -> dict[str, Any]:
    started = time.monotonic()
    presolve = audit_presolve()
    sparse = audit_sparse_no_good()
    dense = audit_dense_no_good()
    script_path = Path(__file__).resolve()
    solver_path = script_path.parent / "hypercube_exact.py"
    checks = {
        "presolve_case_count_matches": (
            presolve["presolve_cases"] == EXPECTED["presolve_cases"]
        ),
        "predicate_count_matches": (
            presolve["selected_set_predicate_checks"]
            == EXPECTED["selected_set_predicate_checks"]
        ),
        "presolve_zero_mismatches": presolve["mismatches"] == 0,
        "removal_event_counts_match": (
            presolve["removal_events"] == EXPECTED["removal_events"]
        ),
        "sparse_no_good_pair_count_matches": (
            sparse["pairs"] == EXPECTED["sparse_no_good_pairs"]
        ),
        "sparse_no_good_zero_mismatches": sparse["mismatches"] == 0,
        "dense_no_good_pair_count_matches": (
            dense["pairs"] == EXPECTED["dense_no_good_pairs"]
        ),
        "dense_no_good_zero_mismatches": dense["mismatches"] == 0,
    }
    return {
        "schema_name": "intersubmod.hypercube_reduction_exhaustive_audit",
        "schema_version": "1.0.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION_METHOD_AUDIT",
        "scope": {
            "presolve_k": [1, 2, 3],
            "all_universe_masks": True,
            "partial_group_multiset_size_max": 3,
            "extra_mandatory_states_max": 2,
            "sparse_no_good_variable_count": [1, 8],
            "dense_no_good_variable_count": [1, 7],
        },
        "independence_boundary": {
            "reducer_under_test": "hypercube_exact._reduce_partial_groups",
            "original_membership_oracle": "direct R/A/X character comparison",
            "symbolic_membership_separately_tested_through_k": 6,
        },
        "sources": {
            "audit_script": {
                "path": str(script_path),
                "sha256": sha256_path(script_path),
            },
            "solver": {
                "path": str(solver_path.resolve()),
                "sha256": sha256_path(solver_path),
            },
        },
        "runtime": {
            "python": platform.python_version(),
            "wall_seconds_diagnostic_only": time.monotonic() - started,
        },
        "presolve": presolve,
        "sparse_fixed_cardinality_no_good": sparse,
        "general_dense_no_good": dense,
        "checks": checks,
        "all_pass": bool(checks) and all(value is True for value in checks.values()),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    receipt = build_receipt()
    if args.json_output is not None:
        output = args.json_output.resolve()
        sidecar = output.with_name(output.name + ".sha256")
        if not args.overwrite:
            for target in (output, sidecar):
                if target.exists():
                    raise FileExistsError(f"refusing to overwrite: {target}")
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(receipt, ensure_ascii=False, indent=2) + "\n"
        output.write_text(payload, encoding="utf-8")
        sidecar.write_text(
            f"{sha256_path(output)}  {output.name}\n", encoding="utf-8"
        )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    if not receipt["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
