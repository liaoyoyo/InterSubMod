#!/usr/bin/env python3
"""Exploratory family-specific read-AF ordering over every candidate tree."""

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import tree_enumeration_solver as S


def read_af_from_colcov(colcov, positions, k):
    if not colcov or len(positions) < k:
        return None
    values = {}
    for index in range(k):
        counts = colcov.get(str(positions[index])) or colcov.get(positions[index])
        if not counts:
            return None
        n_ref, n_alt = counts[:2]
        total = n_ref + n_alt
        if not total:
            return None
        values[index] = n_alt / total
    return values


def unlabel(node):
    if node == "ROOT":
        return frozenset()
    value = node[2:] if node.startswith("H_") else node
    return frozenset(index for index, allele in enumerate(value) if allele == "A")


def ordering_score(edges, read_af):
    score = 0.0
    comparisons = []
    for parent, child in edges:
        parent_set, child_set = unlabel(parent), unlabel(child)
        acquired = child_set - parent_set
        if len(acquired) != 1:
            continue
        child_mutation = next(iter(acquired))
        for ancestor_mutation in parent_set:
            delta = read_af[ancestor_mutation] - read_af[child_mutation]
            score += delta
            comparisons.append(delta)
    return score, comparisons


def posterior(scores, temperature):
    maximum = max(scores)
    weights = [math.exp((score - maximum) / temperature) for score in scores]
    total = sum(weights) or 1.0
    return [weight / total for weight in weights]


def cn_class(value):
    return value if value in {"neutral", "gain", "loss", "loh"} else "unavailable"


def load_groups(sample_dir):
    groups = {}
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        doc = json.loads(path.read_text(encoding="utf-8"))
        for group in doc.get("groups", []):
            key = f"{group['chrom']}:{group['start']}-{group['end']}"
            groups[key] = group
    return groups


def analyze_sample(meta, sample_dir, temperatures, thresholds, margins):
    sample = meta["sample"]
    layered_path = sample_dir / f"layered_reconstruction_{sample}.json"
    layered = json.loads(layered_path.read_text(encoding="utf-8"))
    groups = load_groups(sample_dir)
    units = [u for u in layered["detail"]
             if u.get("is_primary_lineage") and not u.get("capped")
             and u.get("L1_base_class") == "ambiguous"]
    prepared = []
    candidate_mismatch = 0
    missing_read_af = 0
    for unit in units:
        group = groups.get(unit["region"])
        if not group:
            missing_read_af += 1
            continue
        family = unit["family"]
        full = (group.get("populations_by_hp", {}) or {}).get(family, {}) or {}
        partial = list(((group.get("subread_groups_by_hp", {}) or {}).get(family, {}) or {}).keys())
        k = len(group.get("positions", []))
        result = S.enumerate_min_trees(full, partial, k, tree_cap=0)
        if result.get("capped") or not result.get("trees_complete"):
            candidate_mismatch += 1
            continue
        if result["n_trees"] != unit["n_trees"]:
            candidate_mismatch += 1
            continue
        values = read_af_from_colcov(
            ((group.get("col_coverage_by_hp", {}) or {}).get(family, {})),
            group.get("positions", []),
            k,
        )
        if values is None:
            missing_read_af += 1
            continue
        scored = [ordering_score(tree["edges"], values) for tree in result["trees"]]
        prepared.append({
            "cn": cn_class(unit.get("cn")),
            "scores": [item[0] for item in scored],
            "comparisons": [item[1] for item in scored],
            "n_trees": result["n_trees"],
        })

    grid = {}
    for temperature in temperatures:
        for threshold in thresholds:
            for margin in margins:
                reached = 0
                strict_neutral_reached = 0
                winner_consistent = 0
                by_cn = Counter()
                for unit in prepared:
                    post = posterior(unit["scores"], temperature)
                    top = max(post)
                    if top < threshold:
                        continue
                    reached += 1
                    by_cn[unit["cn"]] += 1
                    winner_index = post.index(top)
                    if all(delta >= -margin for delta in unit["comparisons"][winner_index]):
                        winner_consistent += 1
                    if unit["cn"] == "neutral":
                        strict_neutral_reached += 1
                key = f"temp={temperature:g}|posterior={threshold:g}|margin={margin:g}"
                grid[key] = {
                    "reached": reached,
                    "reach_fraction": round(reached / len(prepared), 6) if prepared else None,
                    "strict_neutral_reached": strict_neutral_reached,
                    "strict_neutral_reach_fraction": round(strict_neutral_reached / len(prepared), 6) if prepared else None,
                    "winner_order_consistent": winner_consistent,
                    "winner_order_consistent_fraction": round(winner_consistent / reached, 6) if reached else None,
                    "by_cn": dict(by_cn),
                }
    default_key = "temp=0.05|posterior=0.6|margin=0.05"
    return {
        "sample": sample,
        "biological_id": meta["biological_id"],
        "n_ambiguous_primary_units": len(units),
        "n_units_analyzed_all_candidates": len(prepared),
        "n_candidate_mismatch_or_incomplete": candidate_mismatch,
        "n_missing_read_af": missing_read_af,
        "default": grid.get(default_key),
        "sensitivity_grid": grid,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    temperatures = [0.025, 0.05, 0.1]
    thresholds = [0.6, 0.75, 0.9]
    margins = [0.03, 0.05, 0.1]
    samples = [analyze_sample(meta, args.run_root / "samples" / meta["sample"],
                              temperatures, thresholds, margins)
               for meta in manifest["samples"]]
    output = {
        "schema_version": "2.0",
        "method_name": "family_specific_read_af_tree_ordering",
        "epistemic_status": "exploratory heuristic; not independent validation",
        "prohibited_claim": "This is not purity/CN-corrected cancer cell fraction (CCF)",
        "candidate_contract": "all non-capped candidate trees are re-enumerated with tree_cap=0",
        "primary_scope": "mutation-bearing HP1/HP2 units only",
        "default_params": {"temperature": 0.05, "posterior_threshold": 0.6, "violation_margin": 0.05},
        "sensitivity": {"temperature": temperatures, "posterior_threshold": thresholds, "violation_margin": margins},
        "samples": samples,
        "all_candidate_sets_complete": all(
            s["n_candidate_mismatch_or_incomplete"] == 0
            and s["n_missing_read_af"] == 0
            and s["n_units_analyzed_all_candidates"] == s["n_ambiguous_primary_units"]
            for s in samples
        ),
    }
    path = args.output or args.run_root / "read_af_tree_ordering.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"READ-AF ORDERING -> {path}")
    for sample in samples:
        default = sample["default"] or {}
        print(f"  {sample['sample']}: analyzed={sample['n_units_analyzed_all_candidates']}/"
              f"{sample['n_ambiguous_primary_units']} reached={default.get('reached')} "
              f"strict-neutral={default.get('strict_neutral_reached')} mismatch={sample['n_candidate_mismatch_or_incomplete']}")
    raise SystemExit(0 if output["all_candidate_sets_complete"] else 1)


if __name__ == "__main__":
    main()
