#!/usr/bin/env python3
"""Separate exact-tree, topology-shape, and region determinacy denominators."""

import argparse
import json
from pathlib import Path


def ratio(numerator, denominator):
    return {"n": numerator, "denominator": denominator,
            "pct": round(100 * numerator / denominator, 3) if denominator else None}


def sample_summary(sample, layered_path, region_path):
    layered = json.loads(layered_path.read_text(encoding="utf-8"))
    region = json.loads(region_path.read_text(encoding="utf-8"))
    primary = [u for u in layered["detail"] if u.get("is_primary_lineage")]
    noncapped = [u for u in primary if not u.get("capped")]
    complete = [u for u in noncapped if u.get("analysis_candidate_set_complete")]
    exact = sum(u.get("n_trees") == 1 for u in complete)
    shape = sum(u.get("n_distinct_shapes_exact") == 1 for u in complete)
    regions_with_primary = [r for r in region["regions"] if r.get("n_primary_lineages", 0) > 0]
    all_determined = sum(r["region_determinacy"] == "all_determined" for r in regions_with_primary)
    return {
        "sample": sample,
        "primary_lineage_units_all": len(primary),
        "capped_not_exact_eligible": len(primary) - len(noncapped),
        "candidate_complete_noncapped": len(complete),
        "exact_tree_unique": ratio(exact, len(complete)),
        "shape_determined": ratio(shape, len(complete)),
        "region_all_determined": ratio(all_determined, len(regions_with_primary)),
        "regions_without_primary_lineage": sum(r.get("n_primary_lineages", 0) == 0 for r in region["regions"]),
        "shape_minus_exact_pp": round(100 * (shape - exact) / len(complete), 3) if complete else None,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.input_manifest.read_text(encoding="utf-8"))
    samples = []
    for meta in manifest["samples"]:
        sample = meta["sample"]
        root = args.run_root / "samples" / sample
        samples.append(sample_summary(sample, root / f"layered_reconstruction_{sample}.json",
                                      root / f"layered_region_view_{sample}.json"))
    output = {"schema_version": "2.0",
              "denominator_contract": {"exact_tree_unique": "non-capped primary lineage units with complete candidate set",
                                       "shape_determined": "same units; all labels/order variants collapsed to canonical shape",
                                       "region_all_determined": "regions with >=1 primary lineage; every primary lineage exact determined"},
              "samples": samples}
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"DETERMINACY LAYERS -> {args.output}")
    for row in samples:
        print(f"  {row['sample']}: exact={row['exact_tree_unique']['pct']}% "
              f"shape={row['shape_determined']['pct']}% region={row['region_all_determined']['pct']}%")


if __name__ == "__main__":
    main()
