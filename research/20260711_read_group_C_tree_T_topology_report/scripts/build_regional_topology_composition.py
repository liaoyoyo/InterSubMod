#!/usr/bin/env python3
"""Compose per-HP topology alternatives into ordered regional HP1/HP2 forests."""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    catalog = json.loads(args.topology.read_text(encoding="utf-8"))
    results = []
    for sample_doc in catalog["samples"]:
        sample = sample_doc["sample"]
        source = json.loads(Path(sample_doc["source"]).read_text(encoding="utf-8"))
        primary_by_region = collections.defaultdict(list)
        for unit in source["detail"]:
            if unit.get("is_primary_lineage"):
                primary_by_region[unit["region"]].append(unit)

        rows = {(row["region"], row["hp"]): row for row in sample_doc["unit_rows"]}
        grain = collections.Counter()
        ordered_shape_candidates = collections.Counter()
        ordered_shape_region_support = collections.Counter()
        n_topology_alternatives = 0
        n_exact_joint_candidates = 0
        sums_by_hp_multiplicity = collections.defaultdict(collections.Counter)
        fully_complete_regions = 0
        incomplete_regions = 0

        for region, units in primary_by_region.items():
            units.sort(key=lambda unit: unit["family"])
            if any(not unit.get("analysis_candidate_set_complete") for unit in units):
                incomplete_regions += 1
                grain[f"hp{len(units)}_incomplete"] += 1
                continue

            fully_complete_regions += 1
            grain[f"hp{len(units)}_complete"] += 1
            unit_rows = [rows[(region, unit["family"])] for unit in units]
            topology_alternatives = 1
            joint_exact_t = 1
            for row in unit_rows:
                topology_alternatives *= len(row["shape_candidate_counts"])
                joint_exact_t *= row["n_trees"]
            n_topology_alternatives += topology_alternatives
            n_exact_joint_candidates += joint_exact_t

            multiplicity = str(len(units))
            sums_by_hp_multiplicity[multiplicity]["regions"] += 1
            sums_by_hp_multiplicity[multiplicity]["topology_alternatives"] += topology_alternatives
            sums_by_hp_multiplicity[multiplicity]["exact_joint_tree_candidates"] += joint_exact_t

            if len(unit_rows) == 1:
                for shape_id, count in unit_rows[0]["shape_candidate_counts"].items():
                    key = f"HP{unit_rows[0]['hp']}={shape_id}"
                    ordered_shape_candidates[key] += count
                    ordered_shape_region_support[key] += 1
            elif len(unit_rows) == 2:
                left, right = unit_rows
                for left_shape, left_count in left["shape_candidate_counts"].items():
                    for right_shape, right_count in right["shape_candidate_counts"].items():
                        key = f"HP{left['hp']}={left_shape}|HP{right['hp']}={right_shape}"
                        ordered_shape_candidates[key] += left_count * right_count
                        ordered_shape_region_support[key] += 1
            else:
                raise RuntimeError(
                    f"unexpected primary HP multiplicity {sample} {region}: {len(units)}"
                )

        ranked = sorted(
            ordered_shape_candidates,
            key=lambda key: (-ordered_shape_candidates[key], key),
        )
        checks = {
            "primary_region_conservation": (
                fully_complete_regions + incomplete_regions == len(primary_by_region)
            ),
            "topology_alternative_conservation": (
                n_topology_alternatives
                == sum(
                    values["topology_alternatives"]
                    for values in sums_by_hp_multiplicity.values()
                )
            ),
            "exact_joint_T_conservation": (
                n_exact_joint_candidates
                == sum(
                    values["exact_joint_tree_candidates"]
                    for values in sums_by_hp_multiplicity.values()
                )
            ),
        }
        if not all(checks.values()):
            raise RuntimeError(f"regional topology checks failed for {sample}: {checks}")

        results.append(
            {
                "sample": sample,
                "primary_regions": len(primary_by_region),
                "fully_complete_regions": fully_complete_regions,
                "incomplete_regions": incomplete_regions,
                "region_grain": dict(grain),
                "sum_region_topology_alternatives": n_topology_alternatives,
                "sum_exact_joint_tree_candidates": n_exact_joint_candidates,
                "sums_by_hp_multiplicity": {
                    key: dict(value) for key, value in sums_by_hp_multiplicity.items()
                },
                "distinct_ordered_regional_shape_signatures": len(ordered_shape_candidates),
                "top_ordered_regional_shapes": [
                    {
                        "signature": key,
                        "exact_joint_candidate_occurrences": ordered_shape_candidates[key],
                        "region_support": ordered_shape_region_support[key],
                    }
                    for key in ranked[:50]
                ],
                "checks": checks,
            }
        )

    output = {
        "schema_version": "1.0",
        "grain": "ordered HP1/HP2 regional forest; complete primary candidate sets only",
        "samples": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"REGIONAL TOPOLOGY -> {args.output}; samples={len(results)}")


if __name__ == "__main__":
    main()
