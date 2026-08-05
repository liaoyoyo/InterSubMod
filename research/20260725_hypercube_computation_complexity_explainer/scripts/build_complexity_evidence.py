#!/usr/bin/env python3
"""Recompute the evidence used by the hypercube computation explainer.

This is a read-only analysis of the frozen 2026-07-24 exact-PS×HP MLHP inputs,
the recurrence-allowed C++ topology outputs, and the exact best-tree signature
census.  It does not run the solver and does not alter canonical evidence.
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort-root", required=True, type=Path)
    parser.add_argument("--signature-summary", required=True, type=Path)
    parser.add_argument("--signature-receipt", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def nearest_rank(values: Iterable[int], probability: float) -> int:
    ordered = sorted(values)
    if not ordered:
        return 0
    index = max(0, math.ceil(probability * len(ordered)) - 1)
    return ordered[index]


def pct(numerator: int, denominator: int) -> float:
    return 0.0 if denominator == 0 else 100.0 * numerator / denominator


def full_cube_tree_count(q: int) -> int:
    result = 1
    for depth in range(1, q + 1):
        result *= depth ** math.comb(q, depth)
    return result


def formula_table() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for q in range(13):
        states = 1 << q
        edges = 0 if q == 0 else q * (1 << (q - 1))
        trees = full_cube_tree_count(q)
        tree_text = str(trees)
        log10_trees = 0.0
        if trees > 0:
            log10_trees = sum(
                math.comb(q, depth) * math.log10(depth)
                for depth in range(1, q + 1)
            )
        rows.append(
            {
                "q": q,
                "hypercube_states": states,
                "directed_unit_flip_edges": edges,
                "maximum_one_partial_group_states": states,
                "connection_lb_submask_loose_bound": 3**q,
                "full_cube_parent_tree_count": tree_text,
                "full_cube_parent_tree_digits": len(tree_text),
                "full_cube_parent_tree_log10": log10_trees,
            }
        )
    return rows


def active_indices(
    k: int, full_patterns: dict[str, int], partial_patterns: dict[str, int]
) -> list[int]:
    active = [False] * k
    for pattern in (*full_patterns.keys(), *partial_patterns.keys()):
        if len(pattern) != k:
            raise ValueError(f"pattern length {len(pattern)} != k={k}")
        for index, symbol in enumerate(pattern):
            if symbol == "A":
                active[index] = True
    return [index for index, value in enumerate(active) if value]


def new_q_bucket() -> dict[str, Any]:
    return {
        "units": 0,
        "family_complete": 0,
        "resource_abstain": 0,
        "ranked": 0,
        "partial_units": 0,
        "full_patterns": 0,
        "partial_patterns": 0,
        "partial_group_state_entries": 0,
        "partial_group_state_entries_values": [],
        "search_nodes": [],
        "solver_microseconds": [],
        "minimum_vertex_set_counts": [],
        "total_tree_counts": [],
        "best_tree_tie_counts": [],
    }


def compact_q_bucket(q: int, bucket: dict[str, Any]) -> dict[str, Any]:
    units = bucket["units"]
    nodes = bucket.pop("search_nodes")
    microseconds = bucket.pop("solver_microseconds")
    partial_values = bucket.pop("partial_group_state_entries_values")
    vertex_counts = bucket.pop("minimum_vertex_set_counts")
    total_trees = bucket.pop("total_tree_counts")
    best_ties = bucket.pop("best_tree_tie_counts")
    bucket.update(
        {
            "q": q,
            "family_complete_pct": pct(bucket["family_complete"], units),
            "resource_abstain_pct": pct(bucket["resource_abstain"], units),
            "ranked_pct": pct(bucket["ranked"], units),
            "search_nodes_p50": nearest_rank(nodes, 0.50),
            "search_nodes_p95": nearest_rank(nodes, 0.95),
            "search_nodes_p99": nearest_rank(nodes, 0.99),
            "search_nodes_max": max(nodes, default=0),
            "solver_ms_p50": nearest_rank(microseconds, 0.50) / 1000.0,
            "solver_ms_p95": nearest_rank(microseconds, 0.95) / 1000.0,
            "solver_ms_p99": nearest_rank(microseconds, 0.99) / 1000.0,
            "solver_ms_max": max(microseconds, default=0) / 1000.0,
            "partial_group_entries_p95": nearest_rank(partial_values, 0.95),
            "partial_group_entries_max": max(partial_values, default=0),
            "minimum_vertex_set_count_max": max(vertex_counts, default=0),
            "total_tree_count_max": str(max(total_trees, default=0)),
            "best_tree_tie_count_max": str(max(best_ties, default=0)),
        }
    )
    return bucket


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} is not a JSON object")
    return value


def analyze(args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, bool]]:
    cohort_root = args.cohort_root.resolve(strict=True)
    summary_path = cohort_root / "summary" / "all7_summary.json"
    cohort_receipt_path = cohort_root / "cohort_receipt.json"
    summary = load_json(summary_path)
    cohort_receipt = load_json(cohort_receipt_path)
    signature_summary = load_json(args.signature_summary)
    signature_receipt = load_json(args.signature_receipt)

    sources = [
        identity(summary_path),
        identity(cohort_receipt_path),
        identity(args.signature_summary),
        identity(args.signature_receipt),
    ]
    q_buckets: dict[int, dict[str, Any]] = collections.defaultdict(new_q_bucket)
    aggregate = collections.Counter()
    extreme_partial: dict[str, Any] | None = None
    maxima: dict[str, dict[str, Any]] = {}

    for sample in SAMPLES:
        sample_dir = cohort_root / "samples" / sample
        mlhp_path = sample_dir / f"{sample}.exact_ps_mlhp.json"
        topology_path = sample_dir / f"{sample}.topology.jsonl"
        sources.extend((identity(mlhp_path), identity(topology_path)))
        mlhp = load_json(mlhp_path)
        groups = mlhp.get("groups")
        if not isinstance(groups, list):
            raise ValueError(f"{mlhp_path} has no groups array")

        with topology_path.open(encoding="utf-8") as topology_handle:
            topology_lines = iter(topology_handle)
            for group_index, group in enumerate(groups):
                try:
                    topology = json.loads(next(topology_lines))
                except StopIteration as error:
                    raise ValueError(f"{topology_path} ended before group {group_index}") from error
                if topology.get("group_index") != group_index:
                    raise ValueError(
                        f"{topology_path} group_index mismatch at {group_index}: "
                        f"{topology.get('group_index')}"
                    )

                hp = group["hp_family"]
                k = len(group["positions"])
                full = (group.get("populations_by_hp") or {}).get(hp) or {}
                partial = (group.get("subread_groups_by_hp") or {}).get(hp) or {}
                active = active_indices(k, full, partial)
                q = len(active)
                if topology.get("original_bit_count") != k:
                    raise ValueError(f"{sample} group {group_index}: original k mismatch")
                if topology.get("active_bit_count") != q:
                    raise ValueError(f"{sample} group {group_index}: active q mismatch")

                unknown_counts = [
                    sum(pattern[index] == "X" for index in active)
                    for pattern in partial
                ]
                partial_state_entries = sum(1 << unknown for unknown in unknown_counts)
                joint_log2 = sum(unknown_counts)
                bucket = q_buckets[q]
                bucket["units"] += 1
                bucket["family_complete"] += topology.get("family_status") == "FAMILY_COMPLETE"
                bucket["resource_abstain"] += (
                    topology.get("family_status") == "ABSTAIN_RESOURCE_LIMIT"
                )
                bucket["ranked"] += topology.get("read_af_status") == "ranked_complete"
                bucket["partial_units"] += bool(partial)
                bucket["full_patterns"] += len(full)
                bucket["partial_patterns"] += len(partial)
                bucket["partial_group_state_entries"] += partial_state_entries
                bucket["partial_group_state_entries_values"].append(partial_state_entries)
                bucket["search_nodes"].append(int(topology.get("search_nodes", 0)))
                bucket["solver_microseconds"].append(
                    int(topology.get("solver_elapsed_microseconds", 0))
                )
                if isinstance(topology.get("minimum_vertex_set_count"), int):
                    bucket["minimum_vertex_set_counts"].append(
                        topology["minimum_vertex_set_count"]
                    )
                if topology.get("total_tree_count") is not None:
                    bucket["total_tree_counts"].append(int(topology["total_tree_count"]))
                if topology.get("best_tree_tie_count") is not None:
                    bucket["best_tree_tie_counts"].append(
                        int(topology["best_tree_tie_count"])
                    )

                aggregate["groups"] += 1
                aggregate["groups_with_partial"] += bool(partial)
                aggregate["full_patterns"] += len(full)
                aggregate["partial_patterns"] += len(partial)
                aggregate["partial_group_state_entries"] += partial_state_entries
                aggregate["maximum_unknowns_in_one_partial"] = max(
                    aggregate["maximum_unknowns_in_one_partial"],
                    max(unknown_counts, default=0),
                )
                aggregate["maximum_partial_groups_in_one_unit"] = max(
                    aggregate["maximum_partial_groups_in_one_unit"], len(partial)
                )
                aggregate["maximum_full_patterns_in_one_unit"] = max(
                    aggregate["maximum_full_patterns_in_one_unit"], len(full)
                )
                aggregate["maximum_joint_completion_log2"] = max(
                    aggregate["maximum_joint_completion_log2"], joint_log2
                )
                aggregate["maximum_partial_group_state_entries_per_unit"] = max(
                    aggregate["maximum_partial_group_state_entries_per_unit"],
                    partial_state_entries,
                )

                candidate_extreme = {
                    "sample": sample,
                    "group_index": group_index,
                    "unit_id": group["unit_id"],
                    "original_k": k,
                    "active_q": q,
                    "full_pattern_count": len(full),
                    "partial_group_count": len(partial),
                    "sum_partial_group_state_entries": partial_state_entries,
                    "naive_joint_completion_log2": joint_log2,
                    "search_nodes": topology.get("search_nodes"),
                    "solver_elapsed_microseconds": topology.get(
                        "solver_elapsed_microseconds"
                    ),
                    "objective_status": topology.get("objective_status"),
                    "family_status": topology.get("family_status"),
                    "read_af_status": topology.get("read_af_status"),
                }
                if (
                    extreme_partial is None
                    or partial_state_entries
                    > extreme_partial["sum_partial_group_state_entries"]
                ):
                    extreme_partial = candidate_extreme

                for field in (
                    "minimum_vertex_set_count",
                    "total_tree_count",
                    "best_tree_tie_count",
                    "search_nodes",
                    "solver_elapsed_microseconds",
                ):
                    raw_value = topology.get(field)
                    if raw_value is None:
                        continue
                    numeric_value = int(raw_value)
                    if field not in maxima or numeric_value > maxima[field]["value"]:
                        maxima[field] = {
                            "value": numeric_value,
                            "sample": sample,
                            "group_index": group_index,
                            "unit_id": group["unit_id"],
                            "active_q": q,
                            "family_status": topology.get("family_status"),
                            "read_af_status": topology.get("read_af_status"),
                        }

            try:
                extra = next(topology_lines)
            except StopIteration:
                extra = None
            if extra is not None:
                raise ValueError(f"{topology_path} has extra topology rows")

    compact_by_q = [
        compact_q_bucket(q, q_buckets[q]) for q in sorted(q_buckets)
    ]
    expected_q = {
        int(key): int(value)
        for key, value in summary["totals"]["active_k_distribution"].items()
    }
    observed_q = {row["q"]: row["units"] for row in compact_by_q}
    cohort_totals = summary["totals"]
    signature_cohort = signature_summary["cohort"]

    checks = {
        "sample_set_exact": tuple(summary["scope"]["samples"]) == SAMPLES,
        "group_total_matches_summary": (
            aggregate["groups"] == cohort_totals["groups_total"] == 98_955
        ),
        "active_q_distribution_matches_summary": observed_q == expected_q,
        "family_complete_matches_summary": (
            sum(row["family_complete"] for row in compact_by_q)
            == cohort_totals["objective_certified_units"]
        ),
        "resource_abstain_matches_summary": (
            sum(row["resource_abstain"] for row in compact_by_q)
            == cohort_totals["objective_abstain_units"]
            == 10_717
        ),
        "ranked_matches_summary": (
            sum(row["ranked"] for row in compact_by_q)
            == cohort_totals["ranked_units"]
            == signature_cohort["ranked_units"]
            == 71_955
        ),
        "signature_receipt_pass": (
            signature_receipt["checks"]["all_pass"] is True
            and signature_receipt["scope"]["global_best_trees_enumerated"]
            == signature_cohort["best_tree_total"]
        ),
        "technical_run_pass_but_not_validation_eligible": (
            cohort_receipt["technical_all_pass"] is True
            and cohort_receipt["validation_evidence_eligible"] is False
            and cohort_receipt["all_mutation_bearing_families_complete"] is False
        ),
        "formula_q12_boundary": (
            formula_table()[-1]["hypercube_states"] == 4096
            and formula_table()[-1]["directed_unit_flip_edges"] == 24576
        ),
    }

    evidence = {
        "schema_name": "intersubmod.hypercube_computation_explainer.evidence",
        "schema_version": "1.0.0",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "task_type": "F_DEMO_ILLUSTRATION",
        "service_goals": ["G4", "G5"],
        "scope": {
            "datasets": list(SAMPLES),
            "biological_ids": 6,
            "chromosomes": "chr1-22",
            "analysis_unit": (
                "dataset × chromosome × exact nonmissing PS × HP1/HP2 × "
                "read-linked bounded block"
            ),
            "solver_model": (
                "recurrence-allowed exact obligation B&B; active q<=12; "
                "max_search_nodes=1000; max_family_size=100000"
            ),
            "status": (
                "technical all7 run; not topology-complete and not production "
                "strict directed-topology validation"
            ),
        },
        "sources": sources,
        "formula_table": formula_table(),
        "actual": {
            "cohort_totals": cohort_totals,
            "cohort_wall_seconds": cohort_receipt["wall_seconds"],
            "partial_group_representation": dict(aggregate),
            "extreme_partial_unit": extreme_partial,
            "by_active_q": compact_by_q,
            "observed_maxima": maxima,
            "signature_census": {
                "ranked_units": signature_cohort["ranked_units"],
                "global_best_trees_enumerated": signature_cohort["best_tree_total"],
                "maximum_best_tree_tie_count": signature_cohort[
                    "maximum_best_tree_tie_count"
                ],
                "resolution": signature_cohort["resolution"],
            },
        },
        "checks": checks,
        "all_pass": all(checks.values()),
        "claim_boundary": (
            "Computation counts describe mutation-state candidate ambiguity and "
            "resource use. They are not clone counts, cellular prevalence, or a "
            "validated strict perfect-phylogeny result."
        ),
    }
    return evidence, checks


def write_json(path: Path, value: dict[str, Any]) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    args = parse_args()
    evidence, checks = analyze(args)
    write_json(args.output, evidence)
    receipt = {
        "schema_name": "intersubmod.hypercube_computation_explainer.receipt",
        "schema_version": "1.0.0",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output": identity(args.output),
        "checks": checks,
        "all_pass": all(checks.values()),
    }
    write_json(args.receipt, receipt)
    print(
        json.dumps(
            {
                "all_pass": receipt["all_pass"],
                "groups": evidence["actual"]["cohort_totals"]["groups_total"],
                "partial_patterns": evidence["actual"][
                    "partial_group_representation"
                ]["partial_patterns"],
                "partial_group_state_entries": evidence["actual"][
                    "partial_group_representation"
                ]["partial_group_state_entries"],
                "ranked": evidence["actual"]["cohort_totals"]["ranked_units"],
                "resource_abstain": evidence["actual"]["cohort_totals"][
                    "objective_abstain_units"
                ],
                "output": str(args.output.resolve()),
                "receipt": str(args.receipt.resolve()),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
