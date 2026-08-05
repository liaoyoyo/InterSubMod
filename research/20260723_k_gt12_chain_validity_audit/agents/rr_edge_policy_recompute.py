#!/usr/bin/env python3
"""Read-only HCC1395 sensitivity audit for RR-only strict graph edges."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Callable, Iterable


Container = tuple[str, str, str]
NodeSet = frozenset[int]


class DisjointSet:
    def __init__(self, nodes: Iterable[int]) -> None:
        self.parent = {node: node for node in nodes}
        self.rank = {node: 0 for node in nodes}

    def find(self, node: int) -> int:
        if self.parent[node] != node:
            self.parent[node] = self.find(self.parent[node])
        return self.parent[node]

    def union(self, left: int, right: int) -> None:
        left = self.find(left)
        right = self.find(right)
        if left == right:
            return
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1


def read_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def make_components(
    nodes_by_container: dict[Container, set[int]],
    edges: list[dict[str, object]],
    predicate: Callable[[dict[str, object]], bool],
) -> dict[Container, tuple[NodeSet, ...]]:
    dsu_by_container = {
        container: DisjointSet(nodes) for container, nodes in nodes_by_container.items()
    }
    for edge in edges:
        if predicate(edge):
            container = edge["container"]
            dsu_by_container[container].union(edge["left"], edge["right"])

    result: dict[Container, tuple[NodeSet, ...]] = {}
    for container, nodes in nodes_by_container.items():
        groups: dict[int, set[int]] = defaultdict(set)
        dsu = dsu_by_container[container]
        for node in nodes:
            groups[dsu.find(node)].add(node)
        result[container] = tuple(
            sorted(
                (frozenset(group) for group in groups.values()),
                key=lambda group: (min(group), len(group), tuple(sorted(group))),
            )
        )
    return result


def summarize(
    components: dict[Container, tuple[NodeSet, ...]],
    edge_count: int,
) -> dict[str, int]:
    all_components = [
        component for values in components.values() for component in values
    ]
    tree_components = [component for component in all_components if len(component) >= 2]
    return {
        "retained_edges": edge_count,
        "components_all": len(all_components),
        "singleton_components": sum(len(component) == 1 for component in all_components),
        "tree_eligible_regions": len(tree_components),
        "tree_linked_memberships": sum(len(component) for component in tree_components),
        "max_k": max(map(len, all_components), default=0),
        "k_gt12_regions": sum(len(component) > 12 for component in tree_components),
    }


def compare_to_primary(
    primary: dict[Container, tuple[NodeSet, ...]],
    alternative: dict[Container, tuple[NodeSet, ...]],
    *,
    source_min_k: int = 2,
) -> dict[str, int]:
    primary_w = [
        (container, component)
        for container, components in primary.items()
        for component in components
        if len(component) >= source_min_k
    ]
    unchanged = 0
    fully_lost = 0
    split_multiple = 0
    partial_one = 0
    affected = 0
    demoted_memberships = 0
    sources_retaining_component_gt12 = 0
    for container, source in primary_w:
        alt_tree = [
            component
            for component in alternative[container]
            if len(component) >= 2 and component <= source
        ]
        retained_nodes = frozenset().union(*alt_tree) if alt_tree else frozenset()
        if any(len(component) > 12 for component in alt_tree):
            sources_retaining_component_gt12 += 1
        if len(alt_tree) == 1 and alt_tree[0] == source:
            unchanged += 1
            continue
        affected += 1
        demoted_memberships += len(source - retained_nodes)
        if not alt_tree:
            fully_lost += 1
        elif len(alt_tree) >= 2:
            split_multiple += 1
        else:
            partial_one += 1
    return {
        "primary_W_total": len(primary_w),
        "primary_W_unchanged": unchanged,
        "primary_W_affected": affected,
        "primary_W_fully_lost": fully_lost,
        "primary_W_split_into_multiple_tree_regions": split_multiple,
        "primary_W_partially_reduced_to_one_tree_region_plus_singletons": partial_one,
        "primary_W_memberships_demoted_to_singleton": demoted_memberships,
        "primary_W_retaining_at_least_one_k_gt12_component": sources_retaining_component_gt12,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve(strict=True)
    chromosome_dirs = sorted(
        (path for path in (root / "chromosomes").glob("chr*") if path.is_dir()),
        key=lambda path: int(path.name.removeprefix("chr")),
    )
    if len(chromosome_dirs) != 22:
        raise RuntimeError(f"expected 22 chromosome directories, found {len(chromosome_dirs)}")

    nodes_by_container: dict[Container, set[int]] = defaultdict(set)
    edges: list[dict[str, object]] = []
    state_mass = {"RR": 0, "RA": 0, "AR": 0, "AA": 0, "total": 0}
    state_conservation_failures = 0
    threshold_flag_failures = 0
    receipt_all_pass_failures = 0
    receipt_output_path_failures = 0
    receipt_output_sha256_failures = 0

    for chrom_dir in chromosome_dirs:
        chrom = chrom_dir.name
        membership_path = chrom_dir / f"HCC1395.{chrom}.site_component_membership.tsv.gz"
        edge_path = chrom_dir / f"HCC1395.{chrom}.endpoint_edges.tsv.gz"
        receipt_path = chrom_dir / "receipt.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("all_pass") is not True:
            receipt_all_pass_failures += 1
        for output_name, output_path in (
            ("edges", edge_path),
            ("membership", membership_path),
        ):
            expected = receipt["outputs"][output_name]
            if Path(expected["path"]).resolve() != output_path.resolve():
                receipt_output_path_failures += 1
            if expected["sha256"] != sha256_path(output_path):
                receipt_output_sha256_failures += 1
        for row in read_tsv_gz(membership_path):
            if int(row["threshold"]) != 3:
                continue
            container = (chrom, row["linkage_basis"], row["phase_set"])
            nodes_by_container[container].add(int(row["site_index"]))
        for row in read_tsv_gz(edge_path):
            support = {
                "RR": int(row["support_RR"]),
                "RA": int(row["support_RA"]),
                "AR": int(row["support_AR"]),
                "AA": int(row["support_AA"]),
            }
            total = int(row["support_total"])
            alt = support["RA"] + support["AR"] + support["AA"]
            state_mass["total"] += total
            for state in ("RR", "RA", "AR", "AA"):
                state_mass[state] += support[state]
            if total != sum(support.values()):
                state_conservation_failures += 1
            if (row["passes_primary_threshold"] == "true") != (total >= 3):
                threshold_flag_failures += 1
            edges.append(
                {
                    "container": (chrom, row["linkage_basis"], row["phase_set"]),
                    "left": int(row["site_i_index"]),
                    "right": int(row["site_j_index"]),
                    "total": total,
                    "alt": alt,
                    **support,
                }
            )

    policies: dict[str, Callable[[dict[str, object]], bool]] = {
        "primary_total_ge3": lambda edge: edge["total"] >= 3,
        "alt_any_with_total_ge3": lambda edge: edge["total"] >= 3 and edge["alt"] >= 1,
        "alt_support_ge3": lambda edge: edge["alt"] >= 3,
    }
    components = {
        name: make_components(nodes_by_container, edges, predicate)
        for name, predicate in policies.items()
    }
    retained_edge_counts = {
        name: sum(predicate(edge) for edge in edges)
        for name, predicate in policies.items()
    }
    primary_qualifying = [edge for edge in edges if edge["total"] >= 3]
    primary_state_mass = {
        state: sum(edge[state] for edge in primary_qualifying)
        for state in ("RR", "RA", "AR", "AA", "total")
    }
    edge_classes = {
        "observed_edges_all_support_levels": len(edges),
        "primary_qualifying_edges": len(primary_qualifying),
        "primary_RR_only_qualifying_edges": sum(
            edge["alt"] == 0 for edge in primary_qualifying
        ),
        "primary_ALT_informative_qualifying_edges": sum(
            edge["alt"] >= 1 for edge in primary_qualifying
        ),
        "primary_edges_ALT_support_ge3": sum(
            edge["alt"] >= 3 for edge in primary_qualifying
        ),
        "primary_qualifying_state_mass": primary_state_mass,
    }
    output = {
        "input_root": str(root),
        "chromosomes": len(chromosome_dirs),
        "containers": len(nodes_by_container),
        "active_node_memberships": sum(map(len, nodes_by_container.values())),
        "edge_state_conservation": {
            "state_mass": state_mass,
            "failure_rows": state_conservation_failures,
            "primary_threshold_flag_failure_rows": threshold_flag_failures,
            "pass": state_conservation_failures == 0 and threshold_flag_failures == 0,
        },
        "input_receipt_checks": {
            "receipts": len(chromosome_dirs),
            "all_pass_failures": receipt_all_pass_failures,
            "output_path_failures": receipt_output_path_failures,
            "output_sha256_failures": receipt_output_sha256_failures,
            "pass": (
                receipt_all_pass_failures == 0
                and receipt_output_path_failures == 0
                and receipt_output_sha256_failures == 0
            ),
        },
        "edge_classes": edge_classes,
        "policies": {
            name: summarize(components[name], retained_edge_counts[name])
            for name in policies
        },
        "impact_vs_primary": {
            name: {
                "all_primary_W": compare_to_primary(
                    components["primary_total_ge3"], components[name]
                ),
                "primary_W_with_k_gt12": compare_to_primary(
                    components["primary_total_ge3"],
                    components[name],
                    source_min_k=13,
                ),
            }
            for name in ("alt_any_with_total_ge3", "alt_support_ge3")
        },
    }
    print(json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
