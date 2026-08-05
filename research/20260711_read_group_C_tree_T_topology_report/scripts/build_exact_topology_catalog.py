#!/usr/bin/env python3
"""Build an exact rooted-unlabeled topology catalog from layered-v2 artifacts.

The script trusts stored trees only when the stored list is complete.  For units
whose display list is truncated, it re-enumerates the full candidate set with
the frozen solver contract and verifies tree, hidden-node, and shape counts.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
SOLVER_DIR = (
    REPO_ROOT
    / "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"
)
sys.path.insert(0, str(SOLVER_DIR))
import tree_enumeration_solver as solver  # noqa: E402


def canon_shape(edges):
    """Return a rooted, directed, unlabeled, sibling-order-free encoding."""
    children = collections.defaultdict(list)
    child_nodes = set()
    all_nodes = set()
    for parent, child in edges:
        children[parent].append(child)
        child_nodes.add(child)
        all_nodes.update((parent, child))
    roots = sorted(all_nodes - child_nodes)

    def visit(node):
        return "(" + "".join(sorted(visit(child) for child in children[node])) + ")"

    return "|".join(sorted(visit(root) for root in roots)) if roots else "()"


def shape_metrics(edges):
    children = collections.defaultdict(list)
    child_nodes = set()
    all_nodes = set()
    for parent, child in edges:
        children[parent].append(child)
        child_nodes.add(child)
        all_nodes.update((parent, child))
    roots = list(all_nodes - child_nodes)
    if not roots:
        return {
            "n_nodes_including_root": 1,
            "n_mutation_state_nodes": 0,
            "n_edges": 0,
            "n_leaves": 1,
            "max_depth": 0,
            "root_degree": 0,
            "n_internal_branch_nodes": 0,
            "max_outdegree": 0,
            "coarse_shape": "root_only",
        }
    stack = [(root, 0) for root in roots]
    depths = []
    degrees = []
    node_depth = {}
    while stack:
        node, depth = stack.pop()
        depths.append(depth)
        node_depth[node] = depth
        degree = len(children[node])
        degrees.append(degree)
        stack.extend((child, depth + 1) for child in children[node])
    root_degree = len(children[roots[0]]) if len(roots) == 1 else None
    internal_branches = sum(
        len(children[node]) >= 2 for node in all_nodes if node not in roots
    )
    max_outdegree = max(degrees)
    n_nonroot = len(all_nodes) - len(roots)
    if n_nonroot == 1:
        coarse = "single"
    elif internal_branches == 0 and root_degree == 1 and max_outdegree <= 1:
        coarse = "linear"
    elif internal_branches == 0 and (root_degree or 0) >= 2:
        coarse = "root_star"
    elif internal_branches > 0 and root_degree == 1:
        coarse = "internal_branch"
    else:
        coarse = "mixed_root_and_internal_branch"
    return {
        "n_nodes_including_root": len(all_nodes),
        "n_mutation_state_nodes": n_nonroot,
        "n_edges": len(edges),
        "n_leaves": sum(degree == 0 for degree in degrees),
        "max_depth": max(depths),
        "root_degree": root_degree,
        "n_internal_branch_nodes": internal_branches,
        "max_outdegree": max_outdegree,
        "coarse_shape": coarse,
    }


def shape_id(signature):
    return "TS-" + hashlib.sha1(signature.encode()).hexdigest()[:10]


def load_raw_groups(sample_dir):
    groups = {}
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        doc = json.loads(path.read_text(encoding="utf-8"))
        for group in doc["groups"]:
            region = f"{group['chrom']}:{group['start']}-{group['end']}"
            if region in groups:
                raise RuntimeError(f"duplicate raw region: {region}")
            groups[region] = group
    return groups


def all_trees_for_unit(unit, raw_group):
    if unit["n_trees"] <= len(unit["trees"]):
        return unit["trees"], False
    family = unit["family"]
    full = (raw_group.get("populations_by_hp") or {}).get(family) or {}
    partial = list(
        ((raw_group.get("subread_groups_by_hp") or {}).get(family) or {}).keys()
    )
    k = len(next(iter(full))) if full else len(partial[0])
    result = solver.enumerate_min_trees(full, partial, k, tree_cap=0)
    checks = {
        "n_trees": result["n_trees"] == unit["n_trees"],
        "n_hidden": result["n_hidden"] == unit["n_hidden"],
        "not_capped": not result["capped"],
        "complete": result["trees_complete"],
        "shape_n": len({canon_shape(tree["edges"]) for tree in result["trees"]})
        == unit["n_distinct_shapes_exact"],
    }
    if not all(checks.values()):
        raise RuntimeError(
            f"rerun mismatch {unit['region']} HP{family}: {checks}"
        )
    return result["trees"], True


def analyze_sample(sample_dir):
    sample = sample_dir.name
    layered_path = next(sample_dir.glob("layered_reconstruction_*.json"))
    doc = json.loads(layered_path.read_text(encoding="utf-8"))
    primary = [unit for unit in doc["detail"] if unit.get("is_primary_lineage")]
    complete = [
        unit for unit in primary if unit.get("analysis_candidate_set_complete")
    ]
    incomplete = [
        unit for unit in primary if not unit.get("analysis_candidate_set_complete")
    ]
    need_rerun = [unit for unit in complete if unit["n_trees"] > len(unit["trees"])]
    raw_groups = load_raw_groups(sample_dir) if need_rerun else {}

    candidate_counts = collections.Counter()
    unit_counts = collections.Counter()
    region_sets = collections.defaultdict(set)
    by_k = collections.defaultdict(collections.Counter)
    by_hp = collections.defaultdict(collections.Counter)
    by_hidden = collections.defaultdict(collections.Counter)
    shape_examples = {}
    total_candidates = 0
    total_incidence = 0
    rerun_units = 0
    unit_rows = []

    for unit in complete:
        trees, reran = all_trees_for_unit(unit, raw_groups.get(unit["region"]))
        rerun_units += int(reran)
        signature_counts = collections.Counter(
            canon_shape(tree["edges"]) for tree in trees
        )
        if len(signature_counts) != unit["n_distinct_shapes_exact"]:
            raise RuntimeError(
                f"shape mismatch {sample} {unit['region']} HP{unit['family']}"
            )
        if sum(signature_counts.values()) != unit["n_trees"]:
            raise RuntimeError(f"tree conservation mismatch {sample} {unit['region']}")
        total_candidates += unit["n_trees"]
        total_incidence += len(signature_counts)
        for signature, count in signature_counts.items():
            stable_id = shape_id(signature)
            candidate_counts[stable_id] += count
            unit_counts[stable_id] += 1
            region_sets[stable_id].add(unit["region"])
            by_k[stable_id][str(unit["n_sSNV"])] += 1
            by_hp[stable_id][unit["family"]] += 1
            by_hidden[stable_id][str(unit["n_hidden"])] += 1
            if stable_id not in shape_examples:
                representative = next(
                    tree
                    for tree in trees
                    if canon_shape(tree["edges"]) == signature
                )
                shape_examples[stable_id] = {
                    "signature": signature,
                    **shape_metrics(representative["edges"]),
                    "example_region": unit["region"],
                    "example_hp": unit["family"],
                    "example_edges": representative["edges"],
                }
        unit_rows.append(
            {
                "region": unit["region"],
                "hp": unit["family"],
                "k": unit["n_sSNV"],
                "n_hidden": unit["n_hidden"],
                "n_trees": unit["n_trees"],
                "shape_candidate_counts": {
                    shape_id(signature): count
                    for signature, count in sorted(signature_counts.items())
                },
            }
        )

    catalog = []
    for stable_id, metadata in shape_examples.items():
        catalog.append(
            {
                "shape_id": stable_id,
                **metadata,
                "candidate_trees": candidate_counts[stable_id],
                "unit_incidence": unit_counts[stable_id],
                "region_incidence": len(region_sets[stable_id]),
                "unit_incidence_by_k": dict(
                    sorted(by_k[stable_id].items(), key=lambda item: int(item[0]))
                ),
                "unit_incidence_by_hp": dict(sorted(by_hp[stable_id].items())),
                "unit_incidence_by_hidden": dict(
                    sorted(
                        by_hidden[stable_id].items(), key=lambda item: int(item[0])
                    )
                ),
            }
        )
    catalog.sort(
        key=lambda row: (-row["candidate_trees"], row["n_nodes_including_root"], row["signature"])
    )
    checks = {
        "candidate_sum_equals_units_ntrees": total_candidates
        == sum(unit["n_trees"] for unit in complete),
        "unit_incidence_equals_exact_shape_sum": total_incidence
        == sum(unit["n_distinct_shapes_exact"] for unit in complete),
        "catalog_candidate_sum": sum(row["candidate_trees"] for row in catalog)
        == total_candidates,
        "catalog_unit_incidence_sum": sum(row["unit_incidence"] for row in catalog)
        == total_incidence,
    }
    return {
        "sample": sample,
        "source": str(layered_path),
        "primary_units": len(primary),
        "complete_units": len(complete),
        "incomplete_units": len(incomplete),
        "rerun_units": rerun_units,
        "candidate_trees": total_candidates,
        "unit_shape_incidence": total_incidence,
        "distinct_shapes": len(catalog),
        "checks": checks,
        "catalog": catalog,
        "unit_rows": unit_rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samples", nargs="*")
    args = parser.parse_args()
    sample_dirs = sorted(
        path for path in (args.run_root / "samples").iterdir() if path.is_dir()
    )
    if args.samples:
        wanted = set(args.samples)
        sample_dirs = [path for path in sample_dirs if path.name in wanted]
    results = []
    for sample_dir in sample_dirs:
        print(f"processing {sample_dir.name}", file=sys.stderr, flush=True)
        results.append(analyze_sample(sample_dir))
    output = {
        "schema_version": "1.0",
        "grain": "primary HP unit; complete non-capped candidate sets",
        "shape_definition": "rooted directed unlabeled unordered AHU-style parentheses",
        "samples": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if not all(all(sample["checks"].values()) for sample in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
