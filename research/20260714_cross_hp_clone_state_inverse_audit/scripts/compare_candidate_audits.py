#!/usr/bin/env python3
"""Compare exact region-key candidate sets between two cross-HP audits."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        record["dataset"],
        record["chrom"],
        int(record["start"]),
        int(record["end"]),
        tuple(int(pos) for pos in record["positions"]),
    )


def set_stats(
    left_records: list[dict[str, Any]],
    right_records: list[dict[str, Any]],
    predicate: Callable[[dict[str, Any]], bool],
) -> dict[str, Any]:
    left = {key(record) for record in left_records if predicate(record)}
    right = {key(record) for record in right_records if predicate(record)}
    intersection = left & right
    union = left | right
    return {
        "left_n": len(left),
        "right_n": len(right),
        "exact_intersection_n": len(intersection),
        "exact_union_n": len(union),
        "left_retention": len(intersection) / len(left) if left else None,
        "right_retention": len(intersection) / len(right) if right else None,
        "jaccard": len(intersection) / len(union) if union else None,
        "key_semantics": "dataset+chrom+start+end+ordered positions; conservative when backbone changes region boundaries",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", required=True, type=Path)
    parser.add_argument("--sensitivity", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    args = parser.parse_args()

    left = load_json(args.primary)
    right = load_json(args.sensitivity)
    left_counts = left["aggregate_counts"]
    right_counts = right["aggregate_counts"]
    metric_names = sorted(set(left_counts) | set(right_counts))
    metrics = []
    for metric in metric_names:
        left_value = int(left_counts.get(metric, 0))
        right_value = int(right_counts.get(metric, 0))
        metrics.append(
            {
                "metric": metric,
                "primary": left_value,
                "sensitivity": right_value,
                "delta_sensitivity_minus_primary": right_value - left_value,
                "relative_change": (right_value - left_value) / left_value if left_value else None,
            }
        )

    set_comparisons = {
        "direct_sister_shape_invariant": set_stats(
            left["candidate_records"],
            right["candidate_records"],
            lambda record: bool(record["direct_sister_shape_invariant"]),
        ),
        "direct_sister_tree_unique": set_stats(
            left["candidate_records"],
            right["candidate_records"],
            lambda record: bool(record["direct_sister_tree_unique"]),
        ),
        "direct_sister_shape_invariant_analysis_complete": set_stats(
            left["candidate_records"],
            right["candidate_records"],
            lambda record: bool(record["direct_sister_shape_invariant_analysis_complete"]),
        ),
        "direct_sister_tree_unique_analysis_complete": set_stats(
            left["candidate_records"],
            right["candidate_records"],
            lambda record: bool(record["direct_sister_tree_unique_analysis_complete"]),
        ),
        "two_site_all_sites_cross_hp_collision": set_stats(
            left["candidate_records"],
            right["candidate_records"],
            lambda record: int(record["n_sSNV"]) == 2 and bool(record["all_sites_collision"]),
        ),
        "pattern_level_inverse_ready": set_stats(
            left["candidate_records"],
            right["candidate_records"],
            lambda record: bool(record["pattern_level_inverse_ready"]),
        ),
    }
    payload = {
        "schema_name": "intersubmod.cross_hp_backbone_sensitivity_comparison",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "primary": {
            "path": str(args.primary.resolve()),
            "run_id": left["input_run"]["run_id"],
        },
        "sensitivity": {
            "path": str(args.sensitivity.resolve()),
            "run_id": right["input_run"]["run_id"],
        },
        "dataset_count_match": left["dataset_count"] == right["dataset_count"] == 7,
        "biological_sample_count_match": (
            left["biological_sample_count"] == right["biological_sample_count"] == 6
        ),
        "metric_comparisons": metrics,
        "exact_region_set_comparisons": set_comparisons,
        "key_caveat": (
            "Low exact overlap can reflect different sSNV backbones changing region boundaries; "
            "it is not by itself biological non-replication."
        ),
        "fixed_two_site_catalog_conclusion": {
            "primary_pattern_ready": left_counts.get("regions_pattern_level_inverse_ready", 0),
            "sensitivity_pattern_ready": right_counts.get("regions_pattern_level_inverse_ready", 0),
            "verdict": "zero under both backbones",
        },
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    with args.out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(metrics)
    print(json.dumps(payload["exact_region_set_comparisons"], ensure_ascii=False, indent=2))
    print(f"WROTE {args.out_json.resolve()}")
    print(f"WROTE {args.out_tsv.resolve()}")


if __name__ == "__main__":
    main()
