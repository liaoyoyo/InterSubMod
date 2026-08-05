#!/usr/bin/env python3
"""Independently audit strict HCC1395 child-block tree eligibility.

This script is intentionally read-only.  It rebuilds each bounded block's
threshold-3 endpoint-induced graph and independently reproduces the downstream
adapter's MINREAD=3 partial-pattern projection from authoritative molecule rows.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path


AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--partition-root", required=True, type=Path)
    parser.add_argument("--strict-root", required=True, type=Path)
    parser.add_argument("--molecule-root", required=True, type=Path)
    parser.add_argument("--min-read", type=int, default=3)
    return parser.parse_args()


def read_tsv_gz(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def positions(value: str) -> tuple[int, ...]:
    return tuple(int(token) for token in value.split(",") if token)


def connected_components(
    vertices: tuple[int, ...], edges: list[tuple[int, int, int]]
) -> list[list[int]]:
    adjacency: dict[int, set[int]] = {vertex: set() for vertex in vertices}
    for left, right, _support in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(vertices)
    components: list[list[int]] = []
    while remaining:
        seed = min(remaining)
        stack = [seed]
        component: list[int] = []
        remaining.remove(seed)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(adjacency[current], reverse=True):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return sorted(components, key=lambda item: (item[0], len(item)))


def main() -> int:
    args = parse_args()
    blocks: dict[str, dict[str, object]] = {}
    route_by_chrom: dict[
        str, dict[tuple[str, str, int], str]
    ] = defaultdict(dict)

    for chrom in AUTOSOMES:
        block_path = (
            args.partition_root
            / "chromosomes"
            / chrom
            / "python_partition"
            / "blocks.tsv.gz"
        )
        for row in read_tsv_gz(block_path):
            block_positions = positions(row["positions1"])
            block_id = row["block_id"]
            block = {
                "chrom": chrom,
                "block_id": block_id,
                "unit_id": row["unit_id"],
                "component_id": row["component_id"],
                "linkage_basis": row["linkage_basis"],
                "hp_family": row["hp_family"],
                "phase_set": row["phase_set"],
                "k": int(row["k"]),
                "positions": block_positions,
                "retained_molecule_weight": int(row["retained_molecule_weight"]),
                "retained_pattern_count": int(row["retained_pattern_count"]),
            }
            if block_id in blocks:
                raise ValueError(f"duplicate block_id: {block_id}")
            blocks[block_id] = block
            for position in block_positions:
                route_key = (row["hp_family"], row["phase_set"], position)
                if route_key in route_by_chrom[chrom]:
                    raise ValueError(f"duplicate route: {chrom} {route_key}")
                route_by_chrom[chrom][route_key] = block_id

    edge_index: dict[
        tuple[str, str, str, str], list[tuple[int, int, int]]
    ] = defaultdict(list)
    for chrom in AUTOSOMES:
        candidates = tuple(
            (args.strict_root / "chromosomes" / chrom).glob("*.endpoint_edges.tsv.gz")
        )
        if len(candidates) != 1:
            raise ValueError(f"{chrom}: expected one endpoint edge table")
        for row in read_tsv_gz(candidates[0]):
            if row["passes_primary_threshold"] != "true":
                continue
            key = (
                chrom,
                row["linkage_basis"],
                row["hp_family"],
                row["phase_set"],
            )
            edge_index[key].append(
                (
                    int(row["pos_i1"]),
                    int(row["pos_j1"]),
                    int(row["support_total"]),
                )
            )

    for block in blocks.values():
        vertex_set = set(block["positions"])
        key = (
            block["chrom"],
            block["linkage_basis"],
            block["hp_family"],
            block["phase_set"],
        )
        internal_edges = [
            edge
            for edge in edge_index.get(key, ())
            if edge[0] in vertex_set and edge[1] in vertex_set
        ]
        components = connected_components(block["positions"], internal_edges)
        block["internal_edges"] = internal_edges
        block["induced_components"] = components
        block["induced_connected"] = len(components) == 1

    patterns: dict[str, Counter[str]] = defaultdict(Counter)
    for chrom in AUTOSOMES:
        molecule_path = (
            args.molecule_root
            / "chromosomes"
            / chrom
            / "extraction"
            / f"HCC1395.{chrom}.molecule_sparse_calls.tsv.gz"
        )
        route = route_by_chrom[chrom]
        for row in read_tsv_gz(molecule_path):
            hp = row["hp_family"]
            phase_set = row["phase_set"]
            if hp not in {"1", "2"} or not phase_set:
                continue
            calls = row["call_codes"]
            observed = positions(row["positions1"])
            projected: dict[str, dict[int, str]] = defaultdict(dict)
            for position, code in zip(observed, calls):
                block_id = route.get((hp, phase_set, position))
                if block_id is not None:
                    projected[block_id][position] = code if code in {"R", "A"} else "X"
            for block_id, call_map in projected.items():
                vector = "".join(
                    call_map.get(position, "X")
                    for position in blocks[block_id]["positions"]
                )
                if set(vector) != {"X"}:
                    patterns[block_id][vector] += 1

    for block_id, block in blocks.items():
        supported = {
            pattern: count
            for pattern, count in patterns.get(block_id, {}).items()
            if count >= args.min_read
        }
        supported_multisite = {
            pattern: count
            for pattern, count in supported.items()
            if sum(code in {"R", "A"} for code in pattern) >= 2
        }
        pattern_edges: list[tuple[int, int, int]] = []
        covered_positions: set[int] = set()
        for pattern, count in supported.items():
            fixed_positions = [
                position
                for position, code in zip(block["positions"], pattern)
                if code in {"R", "A"}
            ]
            covered_positions.update(fixed_positions)
            for left_index, left in enumerate(fixed_positions):
                for right in fixed_positions[left_index + 1 :]:
                    pattern_edges.append((left, right, count))
        pattern_components = connected_components(block["positions"], pattern_edges)
        block["supported_patterns"] = supported
        block["adapter_pattern_supported"] = bool(supported)
        block["supported_multisite_patterns"] = supported_multisite
        block["adapter_multisite_pattern_supported"] = bool(supported_multisite)
        block["supported_pattern_covered_all_sites"] = (
            covered_positions == set(block["positions"])
        )
        block["supported_pattern_graph_connected"] = len(pattern_components) == 1

    block_values = list(blocks.values())
    k_ge_2 = [block for block in block_values if block["k"] >= 2]
    zero_weight = [
        block
        for block in k_ge_2
        if block["retained_molecule_weight"] == 0
    ]
    disconnected = [
        block for block in k_ge_2 if not block["induced_connected"]
    ]

    def summarize(values: list[dict[str, object]]) -> dict[str, int]:
        return {
            "total": len(values),
            "adapter_pattern_supported": sum(
                bool(value["adapter_pattern_supported"]) for value in values
            ),
            "adapter_pattern_unsupported": sum(
                not bool(value["adapter_pattern_supported"]) for value in values
            ),
            "adapter_multisite_pattern_supported": sum(
                bool(value["adapter_multisite_pattern_supported"]) for value in values
            ),
            "supported_pattern_covered_all_sites": sum(
                bool(value["supported_pattern_covered_all_sites"]) for value in values
            ),
            "supported_pattern_graph_connected": sum(
                bool(value["supported_pattern_graph_connected"]) for value in values
            ),
        }

    def exception(block: dict[str, object]) -> dict[str, object]:
        repartitioned_components: list[dict[str, object]] = []
        position_to_index = {
            position: index for index, position in enumerate(block["positions"])
        }
        for component in block["induced_components"]:
            projected_counts: Counter[str] = Counter()
            indices = [position_to_index[position] for position in component]
            for pattern, count in patterns.get(block["block_id"], {}).items():
                vector = "".join(pattern[index] for index in indices)
                if set(vector) != {"X"}:
                    projected_counts[vector] += count
            supported_component = {
                pattern: count
                for pattern, count in projected_counts.items()
                if count >= args.min_read
            }
            repartitioned_components.append(
                {
                    "positions": component,
                    "k": len(component),
                    "supported_patterns_after_reprojection": supported_component,
                    "tree_input_eligible_after_reprojection": (
                        len(component) >= 2 and bool(supported_component)
                    ),
                }
            )
        return {
            "chrom": block["chrom"],
            "block_id": block["block_id"],
            "unit_id": block["unit_id"],
            "component_id": block["component_id"],
            "hp_family": block["hp_family"],
            "phase_set": block["phase_set"],
            "k": block["k"],
            "positions": block["positions"],
            "retained_molecule_weight": block["retained_molecule_weight"],
            "retained_pattern_count": block["retained_pattern_count"],
            "internal_edges": block["internal_edges"],
            "induced_components": block["induced_components"],
            "adapter_pattern_supported": block["adapter_pattern_supported"],
            "supported_patterns": block["supported_patterns"],
            "supported_multisite_patterns": block["supported_multisite_patterns"],
            "supported_pattern_covered_all_sites": block[
                "supported_pattern_covered_all_sites"
            ],
            "supported_pattern_graph_connected": block[
                "supported_pattern_graph_connected"
            ],
            "repartitioned_components": repartitioned_components,
        }

    output = {
        "parameters": {
            "min_read": args.min_read,
            "edge_threshold": 3,
            "partial_policy": "R/A retained; O/D/S/L/X and absent calls become X",
        },
        "counts": {
            "blocks_total": len(block_values),
            "blocks_k_eq_1": sum(block["k"] == 1 for block in block_values),
            "blocks_k_ge_2": len(k_ge_2),
            "blocks_k_ge_2_induced_connected": sum(
                bool(block["induced_connected"]) for block in k_ge_2
            ),
            "blocks_k_ge_2_induced_disconnected": len(disconnected),
            "blocks_retained_weight_zero_all": sum(
                block["retained_molecule_weight"] == 0 for block in block_values
            ),
            "blocks_retained_weight_zero_k_ge_2": len(zero_weight),
        },
        "adapter_projection": {
            "all_k_ge_2": summarize(k_ge_2),
            "zero_retained_weight_k_ge_2": summarize(zero_weight),
            "disconnected_k_ge_2": summarize(disconnected),
            "zero_weight_and_disconnected": summarize(
                [
                    block
                    for block in disconnected
                    if block["retained_molecule_weight"] == 0
                ]
            ),
        },
        "disconnected_blocks": [exception(block) for block in disconnected],
        "zero_weight_supported_examples": [
            exception(block)
            for block in zero_weight
            if block["adapter_pattern_supported"]
        ][:10],
        "zero_weight_unsupported_examples": [
            exception(block)
            for block in zero_weight
            if not block["adapter_pattern_supported"]
        ][:10],
    }
    print(json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
