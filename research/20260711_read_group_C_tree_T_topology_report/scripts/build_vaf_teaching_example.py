#!/usr/bin/env python3
"""Build one reproducible HCC1395 multi-tree read-AF ranking example."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def canonical_shape(edges):
    children = defaultdict(list)
    for parent, child in edges:
        children[parent].append(child)

    def visit(node):
        return "(" + "".join(sorted(visit(child) for child in children[node])) + ")"

    return visit("ROOT")


def load_group(sample_dir, target_region):
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        for group in document.get("groups", []):
            region = f"{group['chrom']}:{group['start']}-{group['end']}"
            if region == target_region:
                return group
    raise RuntimeError(f"region not found: {target_region}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--sample", default="HCC1395")
    parser.add_argument("--region", required=True)
    parser.add_argument("--hp", required=True)
    parser.add_argument("--method-script-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.method_script_dir))
    import tree_enumeration_solver as solver
    from read_af_tree_ordering_multisample import (
        ordering_score,
        posterior,
        read_af_from_colcov,
    )

    group = load_group(args.sample_dir, args.region)
    layered = json.loads(
        (args.sample_dir / f"layered_reconstruction_{args.sample}.json").read_text(
            encoding="utf-8"
        )
    )
    unit = next(
        unit
        for unit in layered["detail"]
        if unit["region"] == args.region
        and unit.get("is_primary_lineage")
        and unit["family"] == args.hp
    )
    full = ((group.get("populations_by_hp") or {}).get(args.hp) or {})
    partial = list(((group.get("subread_groups_by_hp") or {}).get(args.hp) or {}).keys())
    positions = group["positions"]
    result = solver.enumerate_min_trees(full, partial, len(positions), tree_cap=0)
    read_af = read_af_from_colcov(
        ((group.get("col_coverage_by_hp") or {}).get(args.hp) or {}),
        positions,
        len(positions),
    )
    if read_af is None:
        raise RuntimeError("target example has incomplete read-AF")

    temperature = 0.05
    threshold = 0.60
    margin = 0.05
    scored = [ordering_score(tree["edges"], read_af) for tree in result["trees"]]
    weights = posterior([score for score, _ in scored], temperature)
    top_index = max(range(len(weights)), key=weights.__getitem__)
    candidates = []
    for index, (tree, (score, comparisons), weight) in enumerate(
        zip(result["trees"], scored, weights), 1
    ):
        candidates.append(
            {
                "candidate": f"T_{index}",
                "edges": tree["edges"],
                "shape_signature": canonical_shape(tree["edges"]),
                "score": score,
                "softmax_weight": weight,
                "ancestry_VAF_deltas": comparisons,
                "selected": index - 1 == top_index,
            }
        )

    colcov = ((group.get("col_coverage_by_hp") or {}).get(args.hp) or {})
    sites = []
    for index, position in enumerate(positions):
        counts = colcov[str(position)]
        sites.append(
            {
                "index": index,
                "position": position,
                "REF_reads": counts[0],
                "ALT_reads": counts[1],
                "read_AF": read_af[index],
            }
        )

    output = {
        "schema_version": "1.0",
        "sample": args.sample,
        "region": args.region,
        "HP": args.hp,
        "CN": unit.get("cn"),
        "k": len(positions),
        "full_read_groups": full,
        "partial_constraints": partial,
        "C": sum("A" in genotype for genotype in full),
        "n_T": result["n_trees"],
        "n_Topo": len({candidate["shape_signature"] for candidate in candidates}),
        "sites": sites,
        "params": {
            "temperature": temperature,
            "top_weight_threshold": threshold,
            "ancestry_violation_margin": margin,
        },
        "candidates": candidates,
        "winner": candidates[top_index]["candidate"],
        "winner_weight": weights[top_index],
        "winner_reaches_threshold": weights[top_index] >= threshold,
        "winner_direction_consistent": all(
            delta >= -margin for delta in scored[top_index][1]
        ),
        "claim_ceiling": "VAF-supported most-likely candidate; exploratory heuristic, not independently confirmed true tree",
        "checks": {
            "candidate_count_matches_layered": result["n_trees"] == unit["n_trees"],
            "candidate_set_complete": bool(result.get("trees_complete")) and not result.get("capped"),
            "unique_top": weights.count(max(weights)) == 1,
        },
    }
    if not all(output["checks"].values()):
        raise RuntimeError(f"example checks failed: {output['checks']}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"VAF EXAMPLE -> {args.output}; T={output['n_T']} Topo={output['n_Topo']} "
        f"winner={output['winner']} weight={output['winner_weight']:.6f}"
    )


if __name__ == "__main__":
    main()
