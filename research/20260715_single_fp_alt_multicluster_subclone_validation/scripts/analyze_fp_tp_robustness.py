#!/usr/bin/env python3
"""Audit whether pooled FP-versus-TP effects survive dataset reweighting and exclusion."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scipy.stats import binomtest


BIOLOGICAL_GROUPS = {
    "HCC1395": ["HCC1395", "HCC1395_DORADO"],
    "COLO829": ["COLO829"],
    "H1437": ["H1437"],
    "H2009": ["H2009"],
    "HCC1937": ["HCC1937"],
    "HCC1954": ["HCC1954"],
}


def aggregate(stats: list[dict[str, Any]]) -> dict[str, Any]:
    n_pairs = sum(int(row["n_pairs"]) for row in stats)
    fp_n = sum(int(row["fp_n"]) for row in stats)
    tp_n = sum(int(row["tp_n"]) for row in stats)
    fp_only = sum(int(row["fp_only_discordant"]) for row in stats)
    tp_only = sum(int(row["tp_only_discordant"]) for row in stats)
    discordant = fp_only + tp_only
    return {
        "n_pairs": n_pairs,
        "fp_n": fp_n,
        "tp_n": tp_n,
        "fp_fraction": fp_n / n_pairs if n_pairs else None,
        "tp_fraction": tp_n / n_pairs if n_pairs else None,
        "risk_difference_fp_minus_tp": (fp_n - tp_n) / n_pairs if n_pairs else None,
        "fp_only_discordant": fp_only,
        "tp_only_discordant": tp_only,
        "exact_mcnemar_p": (
            binomtest(fp_only, discordant, 0.5, alternative="two-sided").pvalue
            if discordant
            else 1.0
        ),
    }


def dataset_risk_difference(summary: dict[str, Any], sample: str, subset: str, endpoint: str) -> float:
    value = summary["per_sample"][sample][subset][endpoint]["paired_risk_difference_fp_minus_tp"]
    return float(value)


def directional_summary(values: dict[str, float]) -> dict[str, Any]:
    ordered = list(values.values())
    return {
        "n_units": len(ordered),
        "n_positive": sum(value > 0 for value in ordered),
        "n_zero": sum(value == 0 for value in ordered),
        "n_negative": sum(value < 0 for value in ordered),
        "equal_weight_mean_risk_difference": sum(ordered) / len(ordered),
        "minimum_risk_difference": min(ordered),
        "maximum_risk_difference": max(ordered),
        "unit_risk_differences": values,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    summary = json.loads(args.comparison_summary.read_text(encoding="utf-8"))
    datasets = sorted(summary["per_sample"])
    expected = sorted(sample for members in BIOLOGICAL_GROUPS.values() for sample in members)
    if datasets != expected:
        raise RuntimeError(f"Unexpected datasets: observed={datasets}, expected={expected}")

    endpoints = [
        "stable_null_multigroup",
        "residual_unexplained_multigroup",
        "phase_anchored_robust_epigenetic_candidate",
    ]
    subsets = ["all_pairs", "both_evaluable"]
    results: dict[str, Any] = {}
    leave_out_rows: list[dict[str, Any]] = []

    for subset in subsets:
        results[subset] = {}
        for endpoint in endpoints:
            dataset_values = {
                sample: dataset_risk_difference(summary, sample, subset, endpoint)
                for sample in datasets
            }
            biological_values = {
                biological_sample: sum(dataset_values[sample] for sample in members) / len(members)
                for biological_sample, members in BIOLOGICAL_GROUPS.items()
            }
            per_dataset_stats = {
                sample: summary["per_sample"][sample][subset][endpoint] for sample in datasets
            }

            leave_one_dataset_out = {}
            for left_out in datasets:
                value = aggregate([stats for sample, stats in per_dataset_stats.items() if sample != left_out])
                leave_one_dataset_out[left_out] = value
                leave_out_rows.append(
                    {
                        "subset": subset,
                        "endpoint": endpoint,
                        "unit_type": "dataset",
                        "left_out": left_out,
                        **value,
                    }
                )

            leave_one_biological_sample_out = {}
            for left_out, members in BIOLOGICAL_GROUPS.items():
                value = aggregate(
                    [stats for sample, stats in per_dataset_stats.items() if sample not in members]
                )
                leave_one_biological_sample_out[left_out] = value
                leave_out_rows.append(
                    {
                        "subset": subset,
                        "endpoint": endpoint,
                        "unit_type": "biological_sample",
                        "left_out": left_out,
                        **value,
                    }
                )

            results[subset][endpoint] = {
                "dataset_equal_weight": directional_summary(dataset_values),
                "biological_sample_equal_weight": directional_summary(biological_values),
                "leave_one_dataset_out_pooled": leave_one_dataset_out,
                "leave_one_biological_sample_out_pooled": leave_one_biological_sample_out,
            }

    output = {
        "schema_name": "intersubmod.fp_matched_tp_reweighting_leave_one_out_audit",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(args.comparison_summary),
        "biological_grouping": BIOLOGICAL_GROUPS,
        "methods": {
            "equal_weight": (
                "Arithmetic mean of within-dataset paired risk differences; HCC1395 and Dorado are "
                "first averaged before the six-biological-sample mean."
            ),
            "leave_one_out": (
                "Pooled site-weighted counts and exact McNemar test after excluding one dataset or "
                "one biological sample. This is a sensitivity audit, not an independent replication test."
            ),
        },
        "results": results,
        "pass": all(math.isfinite(value) for subset in results.values() for endpoint in subset.values()
                    for value in endpoint["dataset_equal_weight"]["unit_risk_differences"].values()),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "fp_matched_tp_robustness_summary.json"
    summary_path.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    table_path = args.output_dir / "fp_matched_tp_leave_one_out.tsv"
    fields = [
        "subset",
        "endpoint",
        "unit_type",
        "left_out",
        "n_pairs",
        "fp_n",
        "tp_n",
        "fp_fraction",
        "tp_fraction",
        "risk_difference_fp_minus_tp",
        "fp_only_discordant",
        "tp_only_discordant",
        "exact_mcnemar_p",
    ]
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(leave_out_rows)

    print(json.dumps({"summary": str(summary_path), "table": str(table_path), "pass": output["pass"]}, indent=2))


if __name__ == "__main__":
    main()
