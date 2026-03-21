#!/usr/bin/env python3
"""Summarize three-way comparison runs from significance_summary.csv outputs."""

import argparse
import csv
import math
import os
from statistics import mean, median


PAIRWISE_FILE = "significance_summary.csv"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        default="output/three_way_comparison",
        help="Three-way comparison root directory",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: same as root)",
    )
    return parser.parse_args()


def load_summary(csv_path):
    with open(csv_path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    keyed = {}
    medians = []
    passed_gating = 0
    for row in rows:
        key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
        try:
            pairwise_median = float(row["PairwiseMedianDist"])
        except (TypeError, ValueError, KeyError):
            pairwise_median = math.nan
        keyed[key] = {
            "pairwise_median": pairwise_median,
            "passed_gating": row.get("PassedGating", "").lower() == "true",
        }
        if not math.isnan(pairwise_median):
            medians.append(pairwise_median)
        if row.get("PassedGating", "").lower() == "true":
            passed_gating += 1

    return rows, keyed, medians, passed_gating


def discover_conditions(root):
    conditions = {}
    for entry in sorted(os.listdir(root)):
        summary_path = os.path.join(root, entry, PAIRWISE_FILE)
        if os.path.isfile(summary_path):
            conditions[entry] = summary_path
    return conditions


def choose_pairwise_conditions(condition_names):
    priorities = {
        "tumor": ["tumor_only_full", "tumor_only"],
        "normal": ["normal_only_full", "normal_only"],
        "mixed": ["mixed_t30n20_full", "mixed_t30n20_full_j7", "mixed_t30n20"],
    }
    chosen = {}
    for label, candidates in priorities.items():
        for candidate in candidates:
            if candidate in condition_names:
                chosen[label] = candidate
                break
    return chosen


def main():
    args = parse_args()
    output_dir = args.output_dir or args.root
    os.makedirs(output_dir, exist_ok=True)

    conditions = discover_conditions(args.root)
    if not conditions:
        raise SystemExit(f"No {PAIRWISE_FILE} found under {args.root}")

    summary_rows = []
    loaded = {}
    for condition, csv_path in conditions.items():
        rows, keyed, medians, passed_gating = load_summary(csv_path)
        loaded[condition] = keyed
        summary_rows.append(
            {
                "condition": condition,
                "regions": len(rows),
                "median_of_median": f"{median(medians):.6f}" if medians else "",
                "mean_of_median": f"{mean(medians):.6f}" if medians else "",
                "frac_lt_0_25": f"{sum(x < 0.25 for x in medians) / len(medians):.6f}" if medians else "",
                "frac_gt_0_50": f"{sum(x > 0.50 for x in medians) / len(medians):.6f}" if medians else "",
                "passed_gating": passed_gating,
            }
        )

    summary_path = os.path.join(output_dir, "three_way_comparison_summary.tsv")
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "condition",
                "regions",
                "median_of_median",
                "mean_of_median",
                "frac_lt_0_25",
                "frac_gt_0_50",
                "passed_gating",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    pairwise_path = os.path.join(output_dir, "three_way_pairwise_delta.tsv")
    chosen = choose_pairwise_conditions(set(conditions))
    if len(chosen) == 3:
        tumor = loaded[chosen["tumor"]]
        normal = loaded[chosen["normal"]]
        mixed = loaded[chosen["mixed"]]
        common_keys = sorted(set(tumor) & set(normal) & set(mixed))

        with open(pairwise_path, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "region_key",
                    "tumor_pairwise_median",
                    "normal_pairwise_median",
                    "mixed_pairwise_median",
                    "delta_nt",
                    "delta_mt",
                    "delta_mn",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            for key in common_keys:
                tumor_m = tumor[key]["pairwise_median"]
                normal_m = normal[key]["pairwise_median"]
                mixed_m = mixed[key]["pairwise_median"]
                if any(math.isnan(v) for v in (tumor_m, normal_m, mixed_m)):
                    continue
                writer.writerow(
                    {
                        "region_key": key,
                        "tumor_pairwise_median": f"{tumor_m:.6f}",
                        "normal_pairwise_median": f"{normal_m:.6f}",
                        "mixed_pairwise_median": f"{mixed_m:.6f}",
                        "delta_nt": f"{normal_m - tumor_m:.6f}",
                        "delta_mt": f"{mixed_m - tumor_m:.6f}",
                        "delta_mn": f"{mixed_m - normal_m:.6f}",
                    }
                )
    else:
        with open(pairwise_path, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "region_key",
                    "tumor_pairwise_median",
                    "normal_pairwise_median",
                    "mixed_pairwise_median",
                    "delta_nt",
                    "delta_mt",
                    "delta_mn",
                ],
                delimiter="\t",
            )
            writer.writeheader()

    print(f"[summarize_three_way_comparison] Wrote {summary_path}")
    print(f"[summarize_three_way_comparison] Wrote {pairwise_path}")


if __name__ == "__main__":
    main()
