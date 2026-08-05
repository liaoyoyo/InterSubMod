#!/usr/bin/env python3
"""Run the pre-decision exact/symbolic/likelihood pilot.

This is a marked PILOT. It reads canonical v5 artifacts but never mutates them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import pathlib
import random
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from hypercube_exact import (  # noqa: E402
    SymbolicPattern,
    enumerate_optimal_vertex_sets,
    fit_vertex_mixture,
    solve_min_hidden,
    state_to_pattern,
    submasks,
    vertex_set_digest,
)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_legacy_solver(path: pathlib.Path):
    spec = importlib.util.spec_from_file_location("canonical_v5_tree_solver", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import legacy solver: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def legacy_vertex_sets(result: dict, k: int) -> set[str]:
    digests = set()
    for node_set in result.get("_feasible_N", []):
        vertices = []
        for alt_set in node_set:
            state = sum(1 << bit for bit in alt_set)
            vertices.append(state)
        digests.add(vertex_set_digest(k, vertices))
    return digests


def compare_case(legacy, full: list[str], partial: list[str], k: int, label: str, max_sets: int = 128) -> dict:
    started = time.perf_counter()
    old = legacy.enumerate_min_trees(full, partial, k, tree_cap=1)
    row = {
        "case_type": label,
        "k": k,
        "n_full": len(full),
        "n_partial": len(partial),
        "legacy_capped": bool(old.get("capped")),
        "legacy_n_hidden": old.get("n_hidden"),
        "legacy_n_vertex_sets": old.get("n_feasible_N"),
        "objective_match": None,
        "vertex_sets_checked": False,
        "vertex_sets_match": None,
    }
    if old.get("capped"):
        row.update({"status": "SKIP_LEGACY_CAPPED", "runtime_seconds": time.perf_counter() - started})
        return row
    new = solve_min_hidden(full, partial, k, universe_mode="observed_alt", time_limit_seconds=10)
    row.update(
        {
            "milp_status": new.status,
            "milp_objective": new.objective,
            "milp_bound": new.objective_bound,
            "milp_gap": new.mip_gap,
            "objective_match": new.status == "OPTIMAL"
            and int(round(new.objective if new.objective is not None else -999)) == int(old["n_hidden"]),
        }
    )
    n_sets = int(old.get("n_feasible_N", 0))
    if new.status == "OPTIMAL" and 0 < n_sets <= max_sets:
        enumeration = enumerate_optimal_vertex_sets(
            full,
            partial,
            k,
            universe_mode="observed_alt",
            max_sets=n_sets + 1,
            time_limit_seconds=10,
        )
        row["vertex_sets_checked"] = True
        row["milp_vertex_sets_complete"] = bool(enumeration.get("complete"))
        new_digests = set(enumeration.get("vertex_set_ids", []))
        old_digests = legacy_vertex_sets(old, k)
        row["vertex_sets_match"] = bool(enumeration.get("complete")) and new_digests == old_digests
        row["milp_n_vertex_sets"] = len(new_digests)
    row["status"] = (
        "PASS"
        if row["objective_match"] and (not row["vertex_sets_checked"] or row["vertex_sets_match"])
        else "FAIL"
    )
    row["runtime_seconds"] = time.perf_counter() - started
    return row


def run_symbolic_exhaustive() -> dict:
    n_patterns = 0
    n_state_checks = 0
    mismatches = 0
    started = time.perf_counter()
    for k in range(1, 9):
        universe = (1 << k) - 1
        states = submasks(universe)
        for symbols in itertools.product("RAX", repeat=k):
            pattern = SymbolicPattern.from_string("".join(symbols))
            symbolic = tuple(v for v in states if pattern.compatible(v, universe))
            explicit = pattern.enumerate_completions(universe)
            n_patterns += 1
            n_state_checks += len(states)
            if symbolic != explicit or len(symbolic) != pattern.n_completions(universe):
                mismatches += 1
    return {
        "k_min": 1,
        "k_max": 8,
        "n_patterns": n_patterns,
        "n_state_checks": n_state_checks,
        "mismatches": mismatches,
        "status": "PASS" if mismatches == 0 else "FAIL",
        "runtime_seconds": time.perf_counter() - started,
    }


def run_random_crosscheck(legacy, seed: int = 20260716, per_k: int = 5) -> list[dict]:
    rng = random.Random(seed)
    rows = []
    for k in range(2, 9):
        for case_idx in range(per_k):
            n_full = rng.randint(1, min(3, (1 << k) - 1))
            states = rng.sample(range(1, 1 << k), n_full)
            full = [state_to_pattern(state, k) for state in states]
            partial = []
            for _ in range(rng.randint(0, 3)):
                symbols = [rng.choice("RAX") for _ in range(k)]
                partial.append("".join(symbols))
            rows.append(compare_case(legacy, full, partial, k, f"seeded_k{k}_{case_idx}"))
    return rows


def run_canonical_crosscheck(legacy, canonical_root: pathlib.Path, per_dataset: int = 5) -> list[dict]:
    rows: list[dict] = []
    samples_root = canonical_root / "samples"
    for sample_dir in sorted(path for path in samples_root.iterdir() if path.is_dir()):
        accepted = 0
        for part_path in sorted(sample_dir.glob("mlhp_part_*.json")):
            payload = json.loads(part_path.read_text())
            for group in payload.get("groups", []):
                k = int(group.get("n_sSNV", 0))
                if not 2 <= k <= 8 or int(group.get("n_sSNV_cap_excluded", 0)) != 0:
                    continue
                for hp in ("1", "2"):
                    full_dict = group.get("populations_by_hp", {}).get(hp, {})
                    partial_dict = group.get("subread_groups_by_hp", {}).get(hp, {})
                    if not full_dict and not partial_dict:
                        continue
                    row = compare_case(
                        legacy,
                        list(full_dict.keys()),
                        list(partial_dict.keys()),
                        k,
                        "canonical_v5",
                    )
                    row.update(
                        {
                            "sample": payload.get("sample", sample_dir.name),
                            "chrom": group.get("chrom"),
                            "start": group.get("start"),
                            "end": group.get("end"),
                            "hp": hp,
                            "source": str(part_path),
                        }
                    )
                    rows.append(row)
                    if row["status"] == "PASS":
                        accepted += 1
                    if accepted >= per_dataset:
                        break
                if accepted >= per_dataset:
                    break
            if accepted >= per_dataset:
                break
    return rows


def run_k9_k12() -> list[dict]:
    rows = []
    for k in range(9, 13):
        chain_states = [(1 << rank) - 1 for rank in range(1, k + 1)]
        even_states = [(1 << rank) - 1 for rank in range(2, k + 1, 2)]
        if k % 2:
            even_states.append((1 << k) - 1)
        singleton_states = [1 << bit for bit in range(k)]
        endpoint_face = "A" + ("X" * (k - 2)) + "A"
        cases = [
            ("full_prefix_chain", [state_to_pattern(v, k) for v in chain_states], [], 0),
            ("even_prefix_missing_connectors", [state_to_pattern(v, k) for v in even_states], [], k // 2),
            ("singleton_full_plus_partial_face", [state_to_pattern(v, k) for v in singleton_states], [endpoint_face], 1),
        ]
        for case_name, full, partial, expected in cases:
            started = time.perf_counter()
            result = solve_min_hidden(full, partial, k, universe_mode="all_loci", time_limit_seconds=20)
            objective = None if result.objective is None else int(round(result.objective))
            rows.append(
                {
                    "case_type": "k9_k12",
                    "case_name": case_name,
                    "k": k,
                    "n_full": len(full),
                    "n_partial": len(partial),
                    "expected_objective": expected,
                    "milp_status": result.status,
                    "milp_objective": objective,
                    "milp_bound": result.objective_bound,
                    "milp_gap": result.mip_gap,
                    "mip_node_count": result.mip_node_count,
                    "n_variables": result.n_variables,
                    "n_constraints": result.n_constraints,
                    "runtime_seconds": time.perf_counter() - started,
                    "status": "PASS" if result.status == "OPTIMAL" and objective == expected and (result.mip_gap or 0) <= 1e-9 else "FAIL",
                }
            )
    return rows


def run_likelihood_checks() -> dict:
    same_a = fit_vertex_mixture([("RR", 40), ("AR", 20), ("RA", 10), ("AA", 5), ("AX", 7)], [0, 1, 2, 3])
    same_b = fit_vertex_mixture([("RR", 40), ("AR", 20), ("RA", 10), ("AA", 5), ("AX", 7)], [3, 1, 0, 2])
    partial_a = fit_vertex_mixture([("RX", 20), ("AX", 20)], [0, 1])
    partial_b = fit_vertex_mixture([("RX", 20), ("AX", 20)], [0, 3])
    true_set = fit_vertex_mixture([("RR", 50), ("AR", 30), ("AA", 20)], [0, 1, 3])
    wrong_set = fit_vertex_mixture([("RR", 50), ("AR", 30), ("AA", 20)], [0, 2, 3])
    checks = {
        "same_vertex_order_delta": abs(same_a.log_likelihood - same_b.log_likelihood),
        "partial_missing_dimension_delta": abs(partial_a.log_likelihood - partial_b.log_likelihood),
        "true_minus_wrong_log_likelihood": true_set.log_likelihood - wrong_set.log_likelihood,
        "all_converged": all(x.converged for x in (same_a, same_b, partial_a, partial_b, true_set, wrong_set)),
        "all_monotone": all(x.monotone for x in (same_a, same_b, partial_a, partial_b, true_set, wrong_set)),
    }
    checks["status"] = (
        "PASS"
        if checks["same_vertex_order_delta"] <= 1e-10
        and checks["partial_missing_dimension_delta"] <= 1e-10
        and checks["true_minus_wrong_log_likelihood"] > 0
        and checks["all_converged"]
        and checks["all_monotone"]
        else "FAIL"
    )
    return checks


def write_tsv(path: pathlib.Path, rows: list[dict]) -> None:
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-root", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--canonical-per-dataset", type=int, default=5)
    args = parser.parse_args()

    canonical_root = args.canonical_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    runner_path = pathlib.Path(__file__).resolve()
    research_solver_path = HERE / "hypercube_exact.py"
    solver_path = canonical_root / "source_bundle/files/imported/003_tree_enumeration_solver.py"
    manifest_path = canonical_root / "input_manifest.snapshot.json"
    legacy = load_legacy_solver(solver_path)

    symbolic = run_symbolic_exhaustive()
    golden = [
        compare_case(legacy, ["RRA", "AAA"], ["RAX"], 3, "golden_rra_aaa_rax"),
        compare_case(legacy, ["AAA"], [], 3, "golden_aaa"),
        compare_case(legacy, ["RRA"], [], 3, "golden_rra"),
    ]
    random_rows = run_random_crosscheck(legacy)
    canonical_rows = run_canonical_crosscheck(legacy, canonical_root, args.canonical_per_dataset)
    large_rows = run_k9_k12()
    likelihood = run_likelihood_checks()
    all_compare = golden + random_rows + canonical_rows

    summary = {
        "schema_version": "pilot-receipt-v1",
        "scope": "PILOT_NOT_FINAL_VALIDATION",
        "canonical_root": str(canonical_root),
        "canonical_manifest": str(manifest_path),
        "canonical_manifest_sha256": sha256(manifest_path),
        "canonical_solver": str(solver_path),
        "canonical_solver_sha256": sha256(solver_path),
        "implementation": {
            "runner": str(runner_path),
            "runner_sha256": sha256(runner_path),
            "research_solver": str(research_solver_path),
            "research_solver_sha256": sha256(research_solver_path),
        },
        "symbolic_exhaustive": symbolic,
        "legacy_milp_crosscheck": {
            "n_cases": len(all_compare),
            "n_pass": sum(row["status"] == "PASS" for row in all_compare),
            "n_fail": sum(row["status"] == "FAIL" for row in all_compare),
            "n_skip_legacy_capped": sum(row["status"] == "SKIP_LEGACY_CAPPED" for row in all_compare),
            "n_vertex_set_checks": sum(bool(row.get("vertex_sets_checked")) for row in all_compare),
            "n_vertex_set_mismatch": sum(row.get("vertex_sets_match") is False for row in all_compare),
            "datasets": sorted({row.get("sample") for row in canonical_rows if row.get("sample")}),
            "canonical_n_cases": len(canonical_rows),
            "canonical_n_pass": sum(row["status"] == "PASS" for row in canonical_rows),
        },
        "k9_k12": {
            "n_cases": len(large_rows),
            "n_pass": sum(row["status"] == "PASS" for row in large_rows),
            "n_fail": sum(row["status"] == "FAIL" for row in large_rows),
            "max_runtime_seconds": max(row["runtime_seconds"] for row in large_rows),
            "max_variables": max(row["n_variables"] for row in large_rows),
        },
        "likelihood": likelihood,
    }
    summary["all_pass"] = (
        symbolic["status"] == "PASS"
        and summary["legacy_milp_crosscheck"]["n_fail"] == 0
        and summary["k9_k12"]["n_fail"] == 0
        and likelihood["status"] == "PASS"
        and len(summary["legacy_milp_crosscheck"]["datasets"]) == 7
    )

    write_tsv(output_dir / "legacy_milp_crosscheck.tsv", all_compare)
    write_tsv(output_dir / "k9_k12_exact_cases.tsv", large_rows)
    (output_dir / "pilot_receipt.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
