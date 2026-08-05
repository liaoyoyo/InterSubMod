#!/usr/bin/env python3
"""Read-only audit from strict two-site edges to solver-visible patterns."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
from typing import Iterable


AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
BlockKey = tuple[str, str]
RouteKey = tuple[str, str, str, int]
Pair = tuple[int, int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--min-read", default=3, type=int)
    return parser.parse_args()


def read_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def parse_positions(value: str) -> tuple[int, ...]:
    positions = tuple(int(token) for token in value.split(",") if token)
    if positions != tuple(sorted(set(positions))):
        raise ValueError("positions must be strictly increasing and unique")
    return positions


def read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected JSON object")
    return value


def connected_components(vertices: tuple[int, ...], edges: set[Pair]) -> tuple[tuple[int, ...], ...]:
    adjacency = {vertex: set() for vertex in vertices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(vertices)
    components: list[tuple[int, ...]] = []
    while remaining:
        seed = min(remaining)
        stack = [seed]
        remaining.remove(seed)
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return tuple(sorted(components))


def pct(numerator: int, denominator: int) -> float | None:
    return round(100.0 * numerator / denominator, 6) if denominator else None


def main() -> int:
    args = parse_args()
    if args.min_read < 1:
        raise ValueError("min-read must be positive")
    root = args.partition_root.resolve(strict=True)

    blocks: dict[BlockKey, dict[str, object]] = {}
    route: dict[RouteKey, BlockKey] = {}
    patterns: dict[BlockKey, Counter[str]] = defaultdict(Counter)
    source_sites: dict[str, set[int]] = defaultdict(set)
    receipt_checks = Counter()

    for chrom in AUTOSOMES:
        chrom_dir = root / "chromosomes" / chrom
        partition_receipt = read_json(chrom_dir / "python_partition" / "receipt.json")
        strict_receipt = read_json(chrom_dir / "strict_regions" / "receipt.json")
        extraction_receipt = read_json(chrom_dir / "extraction" / "receipt.json")
        comparison = read_json(chrom_dir / "comparison.json")
        receipt_checks["chromosomes"] += 1
        receipt_checks["partition_pass"] += partition_receipt.get("all_pass") is True
        receipt_checks["strict_pass"] += strict_receipt.get("all_pass") is True
        receipt_checks["extraction_pass"] += extraction_receipt.get("all_pass") is True
        receipt_checks["comparison_pass"] += (
            comparison.get("all_pass") is True
            and int(comparison.get("mismatch_count", -1)) == 0
        )

        block_path = chrom_dir / "python_partition" / "blocks.tsv.gz"
        for row in read_tsv_gz(block_path):
            block_key = (chrom, row["block_id"])
            positions = parse_positions(row["positions1"])
            if len(positions) != int(row["k"]):
                raise ValueError(f"{block_path}: k mismatch for {row['block_id']}")
            if block_key in blocks:
                raise ValueError(f"duplicate block key: {block_key}")
            blocks[block_key] = {
                "chrom": chrom,
                "block_id": row["block_id"],
                "unit_id": row["unit_id"],
                "component_id": row["component_id"],
                "linkage_basis": row["linkage_basis"],
                "hp_family": row["hp_family"],
                "phase_set": row["phase_set"],
                "positions": positions,
                "k": len(positions),
            }
            source_sites[row["component_id"]].update(positions)
            for position in positions:
                route_key = (chrom, row["hp_family"], row["phase_set"], position)
                if route_key in route:
                    raise ValueError(f"duplicate route: {route_key}")
                route[route_key] = block_key

    internal_strict_edges: dict[BlockKey, set[Pair]] = defaultdict(set)
    source_edge_counts: dict[str, Counter[str]] = defaultdict(Counter)
    strict_edge_counts = Counter()
    cross_block_examples: list[dict[str, object]] = []
    for chrom in AUTOSOMES:
        chrom_dir = root / "chromosomes" / chrom
        edge_candidates = tuple(
            (chrom_dir / "strict_regions").glob("*.endpoint_edges.tsv.gz")
        )
        if len(edge_candidates) != 1:
            raise ValueError(f"{chrom}: expected one strict edge file")
        for row in read_tsv_gz(edge_candidates[0]):
            if row["passes_primary_threshold"] != "true":
                continue
            strict_edge_counts["primary_edges"] += 1
            left = int(row["pos_i1"])
            right = int(row["pos_j1"])
            pair = (left, right)
            left_key = route.get((chrom, row["hp_family"], row["phase_set"], left))
            right_key = route.get((chrom, row["hp_family"], row["phase_set"], right))
            if left_key is None or right_key is None:
                strict_edge_counts["endpoint_route_missing"] += 1
                continue
            left_block = blocks[left_key]
            right_block = blocks[right_key]
            if left_block["component_id"] != right_block["component_id"]:
                raise ValueError("strict edge endpoints route to different source components")
            component_id = str(left_block["component_id"])
            source_edge_counts[component_id]["total"] += 1
            if left_key == right_key:
                strict_edge_counts["within_block"] += 1
                source_edge_counts[component_id]["within_block"] += 1
                internal_strict_edges[left_key].add(pair)
            else:
                strict_edge_counts["cross_block"] += 1
                source_edge_counts[component_id]["cross_block"] += 1
                if len(cross_block_examples) < 12:
                    cross_block_examples.append(
                        {
                            "chrom": chrom,
                            "component_id": component_id,
                            "left": left,
                            "right": right,
                            "left_block": left_block["block_id"],
                            "right_block": right_block["block_id"],
                            "support_total": int(row["support_total"]),
                        }
                    )

    for chrom in AUTOSOMES:
        chrom_dir = root / "chromosomes" / chrom
        molecule_candidates = tuple(
            (chrom_dir / "extraction").glob("*.molecule_sparse_calls.tsv.gz")
        )
        if len(molecule_candidates) != 1:
            raise ValueError(f"{chrom}: expected one molecule sparse-call file")
        for row in read_tsv_gz(molecule_candidates[0]):
            hp = row["hp_family"]
            phase_set = row["phase_set"]
            if hp not in {"1", "2"} or not phase_set:
                continue
            observed = parse_positions(row["positions1"])
            codes = row["call_codes"]
            if len(observed) != len(codes):
                raise ValueError(f"{chrom}: molecule call vector mismatch")
            projected: dict[BlockKey, dict[int, str]] = defaultdict(dict)
            for position, code in zip(observed, codes):
                block_key = route.get((chrom, hp, phase_set, position))
                if block_key is not None:
                    projected[block_key][position] = code if code in {"R", "A"} else "X"
            for block_key, call_map in projected.items():
                vector = "".join(
                    call_map.get(position, "X")
                    for position in blocks[block_key]["positions"]
                )
                if set(vector) != {"X"}:
                    patterns[block_key][vector] += 1

    block_counts = Counter()
    tree_input_examples: dict[str, list[dict[str, object]]] = defaultdict(list)
    pattern_pair_counts = Counter()
    source_pattern_counts: dict[str, Counter[str]] = defaultdict(Counter)
    supported_pattern_counts = Counter()
    tree_input_k = Counter()
    tree_input_memberships = Counter()

    for block_key, block in blocks.items():
        positions = block["positions"]
        block_counts["all"] += 1
        if block["k"] < 2:
            block_counts["k1"] += 1
            continue
        block_counts["k_ge2"] += 1
        supported = {
            pattern: count
            for pattern, count in patterns.get(block_key, {}).items()
            if count >= args.min_read
        }
        strict_pairs = internal_strict_edges.get(block_key, set())
        strict_components = connected_components(positions, strict_pairs)
        component_id = str(block["component_id"])
        pattern_pair_counts["all_within_block_strict_edges"] += len(strict_pairs)
        source_pattern_counts[component_id]["strict_internal_edges"] += len(strict_pairs)
        if len(strict_components) == 1:
            block_counts["strict_induced_connected"] += 1
        else:
            block_counts["strict_induced_disconnected"] += 1
        if not supported:
            block_counts["adapter_pattern_unsupported"] += 1
            block_counts["adapter_unsupported_internal_strict_edges"] += len(strict_pairs)
            pattern_pair_counts[
                "unsupported_block_strict_edges_not_solver_visible"
            ] += len(strict_pairs)
            source_pattern_counts[component_id][
                "strict_internal_edges_not_solver_visible"
            ] += len(strict_pairs)
            continue

        block_counts["tree_input"] += 1
        tree_input_k[block["k"]] += 1
        tree_input_memberships["total"] += block["k"]
        fixed_sets: list[set[int]] = []
        covered: set[int] = set()
        alt_positions: set[int] = set()
        pattern_pairs: set[Pair] = set()
        has_full = False
        max_fixed = 0
        for pattern, count in supported.items():
            fixed = {
                position
                for position, code in zip(positions, pattern)
                if code in {"R", "A"}
            }
            alts = {
                position
                for position, code in zip(positions, pattern)
                if code == "A"
            }
            fixed_sets.append(fixed)
            covered.update(fixed)
            alt_positions.update(alts)
            max_fixed = max(max_fixed, len(fixed))
            has_full |= "X" not in pattern
            supported_pattern_counts["patterns"] += 1
            supported_pattern_counts[
                "fixed_single" if len(fixed) == 1 else "fixed_multi"
            ] += 1
            supported_pattern_counts["molecule_block_weight"] += count
            fixed_ordered = sorted(fixed)
            for left_index, left in enumerate(fixed_ordered):
                for right in fixed_ordered[left_index + 1 :]:
                    pattern_pairs.add((left, right))

        pattern_components = connected_components(positions, pattern_pairs)
        has_multisite = max_fixed >= 2
        coverage_all = covered == set(positions)
        pattern_connected = len(pattern_components) == 1
        mutation_bearing = bool(alt_positions)
        strict_missing = strict_pairs - pattern_pairs
        pattern_extra = pattern_pairs - strict_pairs
        block_counts["has_multisite_supported_pattern"] += has_multisite
        block_counts["only_single_fixed_supported_patterns"] += not has_multisite
        block_counts["supported_patterns_cover_all_sites"] += coverage_all
        block_counts["supported_patterns_miss_some_sites"] += not coverage_all
        block_counts["pattern_graph_connected"] += pattern_connected
        block_counts["pattern_graph_disconnected"] += not pattern_connected
        block_counts["pattern_graph_disconnected_coverage_all"] += (
            not pattern_connected and coverage_all
        )
        block_counts["pattern_graph_disconnected_coverage_incomplete"] += (
            not pattern_connected and not coverage_all
        )
        block_counts["multisite_but_pattern_graph_disconnected"] += (
            has_multisite and not pattern_connected
        )
        block_counts["has_full_supported_pattern"] += has_full
        block_counts["partial_only_supported_patterns"] += not has_full
        block_counts["mutation_bearing"] += mutation_bearing
        block_counts["reference_only"] += not mutation_bearing
        block_counts["alt_universe_size_0"] += len(alt_positions) == 0
        block_counts["alt_universe_size_1"] += len(alt_positions) == 1
        block_counts["alt_universe_size_ge2"] += len(alt_positions) >= 2
        block_counts["strict_connected_but_pattern_disconnected"] += (
            len(strict_components) == 1 and not pattern_connected
        )
        block_counts["pattern_connected_alt_universe_ge2"] += (
            pattern_connected and len(alt_positions) >= 2
        )
        block_counts["pattern_connected_alt_universe_1"] += (
            pattern_connected and len(alt_positions) == 1
        )
        block_counts["pattern_connected_reference_only"] += (
            pattern_connected and len(alt_positions) == 0
        )
        block_counts["pattern_disconnected_alt_universe_ge2"] += (
            not pattern_connected and len(alt_positions) >= 2
        )
        block_counts["only_single_fixed_coverage_all"] += (
            not has_multisite and coverage_all
        )
        block_counts["only_single_fixed_coverage_incomplete"] += (
            not has_multisite and not coverage_all
        )
        block_counts["only_single_fixed_alt_universe_ge2"] += (
            not has_multisite and len(alt_positions) >= 2
        )
        block_counts["only_single_fixed_alt_universe_1"] += (
            not has_multisite and len(alt_positions) == 1
        )
        block_counts["only_single_fixed_reference_only"] += (
            not has_multisite and len(alt_positions) == 0
        )
        block_counts["all_internal_strict_edges_solver_visible"] += not strict_missing
        block_counts["some_internal_strict_edge_not_solver_visible"] += bool(strict_missing)
        block_counts["unexpected_pattern_pair_without_strict_edge"] += bool(pattern_extra)
        tree_input_memberships["covered"] += len(covered)
        tree_input_memberships["not_covered"] += len(set(positions) - covered)

        pattern_pair_counts["tree_input_strict_internal_edges"] += len(strict_pairs)
        pattern_pair_counts["solver_visible_strict_edges"] += len(strict_pairs & pattern_pairs)
        pattern_pair_counts[
            "tree_input_strict_internal_edges_not_solver_visible"
        ] += len(strict_missing)
        pattern_pair_counts["all_within_block_strict_edges_not_solver_visible"] += len(
            strict_missing
        )
        pattern_pair_counts["unexpected_pattern_pairs"] += len(pattern_extra)
        pattern_pair_counts["unique_solver_pattern_pairs"] += len(pattern_pairs)
        source_pattern_counts[component_id]["solver_visible_strict_edges"] += len(
            strict_pairs & pattern_pairs
        )
        source_pattern_counts[component_id]["strict_internal_edges_not_solver_visible"] += len(
            strict_missing
        )

        categories = []
        if not has_multisite:
            categories.append("only_single_fixed")
        if not coverage_all:
            categories.append("coverage_incomplete")
        if not pattern_connected:
            categories.append("pattern_disconnected")
        if strict_missing:
            categories.append("strict_edge_missing")
        for category in categories:
            if len(tree_input_examples[category]) < 8:
                tree_input_examples[category].append(
                    {
                        "chrom": block["chrom"],
                        "block_id": block["block_id"],
                        "component_id": component_id,
                        "k": block["k"],
                        "positions": positions,
                        "supported_patterns": supported,
                        "strict_internal_edges": len(strict_pairs),
                        "strict_edges_not_solver_visible": len(strict_missing),
                        "pattern_components": pattern_components,
                        "covered_sites": len(covered),
                        "alt_universe_sites": len(alt_positions),
                    }
                )

    all_component_ids = {str(block["component_id"]) for block in blocks.values()}
    source_counts = Counter()
    for component_id in all_component_ids:
        edge_info = source_edge_counts.get(component_id, Counter())
        pattern_info = source_pattern_counts.get(component_id, Counter())
        source_counts["source_W"] += 1
        source_counts["source_W_k_gt12"] += len(source_sites[component_id]) > 12
        source_counts["with_cross_block_strict_edge"] += edge_info["cross_block"] > 0
        source_counts["k_gt12_with_cross_block_strict_edge"] += (
            len(source_sites[component_id]) > 12 and edge_info["cross_block"] > 0
        )
        source_counts["k_le12_with_cross_block_strict_edge"] += (
            len(source_sites[component_id]) <= 12 and edge_info["cross_block"] > 0
        )
        source_counts["with_internal_strict_edge_not_solver_visible"] += (
            pattern_info["strict_internal_edges_not_solver_visible"] > 0
        )
        source_counts["with_any_strict_edge_not_solver_visible"] += (
            edge_info["cross_block"] > 0
            or pattern_info["strict_internal_edges_not_solver_visible"] > 0
        )
        source_counts["all_strict_edges_solver_visible"] += (
            edge_info["cross_block"] == 0
            and pattern_info["strict_internal_edges_not_solver_visible"] == 0
        )
        source_counts["with_at_least_one_solver_visible_strict_edge"] += (
            pattern_info["solver_visible_strict_edges"] > 0
        )
        source_counts["with_zero_solver_visible_strict_edges"] += (
            pattern_info["solver_visible_strict_edges"] == 0
        )

    tree_input = block_counts["tree_input"]
    output = {
        "parameters": {
            "min_read": args.min_read,
            "strict_edge_threshold": 3,
            "pattern_projection": "R/A retained; O/D/S/L/X and absent become X",
            "pattern_edge": "two block sites fixed R/A in the same MINREAD-supported exact vector",
        },
        "inputs": {
            "partition_root": str(root),
            "chromosomes": list(AUTOSOMES),
        },
        "receipt_checks": dict(receipt_checks),
        "block_counts": dict(block_counts),
        "tree_input_rates_pct": {
            "of_k_ge2_blocks": pct(tree_input, block_counts["k_ge2"]),
            "has_multisite_supported_pattern": pct(
                block_counts["has_multisite_supported_pattern"], tree_input
            ),
            "only_single_fixed_supported_patterns": pct(
                block_counts["only_single_fixed_supported_patterns"], tree_input
            ),
            "supported_patterns_cover_all_sites": pct(
                block_counts["supported_patterns_cover_all_sites"], tree_input
            ),
            "pattern_graph_connected": pct(
                block_counts["pattern_graph_connected"], tree_input
            ),
            "strict_connected_but_pattern_disconnected": pct(
                block_counts["strict_connected_but_pattern_disconnected"], tree_input
            ),
            "mutation_bearing": pct(block_counts["mutation_bearing"], tree_input),
            "alt_universe_size_ge2": pct(
                block_counts["alt_universe_size_ge2"], tree_input
            ),
            "pattern_connected_alt_universe_ge2": pct(
                block_counts["pattern_connected_alt_universe_ge2"], tree_input
            ),
        },
        "tree_input_k_distribution": dict(sorted(tree_input_k.items())),
        "tree_input_memberships": dict(tree_input_memberships),
        "supported_pattern_counts": dict(supported_pattern_counts),
        "strict_edge_counts": dict(strict_edge_counts),
        "pattern_pair_counts": dict(pattern_pair_counts),
        "strict_edge_visibility_pct": {
            "cross_block_of_all_strict_edges": pct(
                strict_edge_counts["cross_block"], strict_edge_counts["primary_edges"]
            ),
            "solver_visible_of_within_block_strict_edges": pct(
                pattern_pair_counts["solver_visible_strict_edges"],
                pattern_pair_counts["all_within_block_strict_edges"],
            ),
            "not_solver_visible_of_within_block_strict_edges": pct(
                pattern_pair_counts["all_within_block_strict_edges_not_solver_visible"]
                + pattern_pair_counts[
                    "unsupported_block_strict_edges_not_solver_visible"
                ],
                pattern_pair_counts["all_within_block_strict_edges"],
            ),
            "solver_visible_of_all_strict_edges": pct(
                pattern_pair_counts["solver_visible_strict_edges"],
                strict_edge_counts["primary_edges"],
            ),
        },
        "source_W_counts": dict(source_counts),
        "examples": dict(tree_input_examples),
        "cross_block_edge_examples": cross_block_examples,
        "invariants": {
            "all_strict_endpoints_routed": strict_edge_counts["endpoint_route_missing"] == 0,
            "strict_edge_partition_conservation": strict_edge_counts["primary_edges"]
            == strict_edge_counts["within_block"] + strict_edge_counts["cross_block"],
            "pattern_pairs_subset_of_internal_strict_edges": pattern_pair_counts[
                "unexpected_pattern_pairs"
            ]
            == 0,
            "tree_input_requires_supported_pattern": tree_input
            + block_counts["adapter_pattern_unsupported"]
            == block_counts["k_ge2"],
        },
    }
    print(json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
