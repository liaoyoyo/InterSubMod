#!/usr/bin/env python3
"""Build purity_qc.tsv and purity_rule_eval.tsv from purity_status.tsv + metrics.json."""

import argparse
import csv
import json
import os


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status-tsv", required=True, help="Input purity_status.tsv")
    parser.add_argument("--sample", default="unknown", help="Sample name")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: same directory as status TSV)",
    )
    return parser.parse_args()


def parse_purity_bin(purity_value):
    try:
        purity = float(purity_value)
    except (TypeError, ValueError):
        return "unknown"

    if purity < 40.0:
        return "lt40"
    if purity < 60.0:
        return "40to60"
    return "ge60"


def to_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def load_metrics(output_dir):
    metrics_path = os.path.join(output_dir, "metrics.json")
    if not os.path.exists(metrics_path):
        return None
    with open(metrics_path, "r") as handle:
        return json.load(handle)


def main():
    args = parse_args()
    output_dir = args.output_dir or os.path.dirname(args.status_tsv)

    with open(args.status_tsv, "r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    purity_qc_path = os.path.join(output_dir, "purity_qc.tsv")
    purity_rule_eval_path = os.path.join(output_dir, "purity_rule_eval.tsv")

    os.makedirs(output_dir, exist_ok=True)

    with open(purity_qc_path, "w", newline="") as qc_handle, open(
        purity_rule_eval_path, "w", newline=""
    ) as rule_handle:
        qc_writer = csv.DictWriter(
            qc_handle,
            fieldnames=[
                "sample",
                "run_tag",
                "subsample",
                "purity",
                "purity_bin",
                "source_mm_tags",
                "source_ml_tags",
                "has_mm_ml",
                "step01",
                "step02",
                "step03",
                "analyzable",
                "tp_regions_analyzed",
                "fp_regions_analyzed",
                "baseline_f1",
                "filtered_f1",
                "f1_delta",
                "notes",
                "output_dir",
            ],
            delimiter="\t",
        )
        rule_writer = csv.DictWriter(
            rule_handle,
            fieldnames=[
                "sample",
                "run_tag",
                "subsample",
                "purity",
                "purity_bin",
                "rule_id",
                "trigger_count",
                "trigger_precision_fp",
                "tp_removed",
                "fp_removed",
                "tp_removed_rate",
                "fp_removed_rate",
                "baseline_f1",
                "filtered_f1",
                "f1_delta",
                "output_dir",
            ],
            delimiter="\t",
        )
        qc_writer.writeheader()
        rule_writer.writeheader()

        for row in rows:
            purity_bin = parse_purity_bin(row.get("estimated_purity"))
            source_mm = to_int(row.get("source_mm_tags"), default=0)
            source_ml = to_int(row.get("source_ml_tags"), default=0)
            has_mm_ml = source_mm > 0 and source_ml > 0
            metrics = load_metrics(row.get("output_dir", ""))

            tp_regions = to_int(row.get("tp_regions_analyzed"), default=0)
            fp_regions = to_int(row.get("fp_regions_analyzed"), default=0)
            analyzable = (
                has_mm_ml
                and row.get("step02") == "OK"
                and row.get("step03") == "OK"
                and (tp_regions > 0 or fp_regions > 0)
            )

            baseline_f1 = row.get("baseline_f1", "")
            filtered_f1 = row.get("filtered_f1", "")
            f1_delta = row.get("f1_delta", "")
            if metrics:
                baseline_f1 = metrics.get("baseline", {}).get("f1", baseline_f1)
                filtered_f1 = metrics.get("filtered", {}).get("f1", filtered_f1)
                f1_delta = metrics.get("improvement", {}).get("f1_delta", f1_delta)

            qc_writer.writerow(
                {
                    "sample": args.sample,
                    "run_tag": row.get("run_tag", ""),
                    "subsample": row.get("subsample", ""),
                    "purity": row.get("estimated_purity", ""),
                    "purity_bin": purity_bin,
                    "source_mm_tags": source_mm,
                    "source_ml_tags": source_ml,
                    "has_mm_ml": "true" if has_mm_ml else "false",
                    "step01": row.get("step01", ""),
                    "step02": row.get("step02", ""),
                    "step03": row.get("step03", ""),
                    "analyzable": "true" if analyzable else "false",
                    "tp_regions_analyzed": tp_regions,
                    "fp_regions_analyzed": fp_regions,
                    "baseline_f1": baseline_f1,
                    "filtered_f1": filtered_f1,
                    "f1_delta": f1_delta,
                    "notes": row.get("notes", ""),
                    "output_dir": row.get("output_dir", ""),
                }
            )

            improvement = metrics.get("improvement", {}) if metrics else {}
            rule_writer.writerow(
                {
                    "sample": args.sample,
                    "run_tag": row.get("run_tag", ""),
                    "subsample": row.get("subsample", ""),
                    "purity": row.get("estimated_purity", ""),
                    "purity_bin": purity_bin,
                    "rule_id": metrics.get("rule", {}).get("rule_id", "unknown") if metrics else "unknown",
                    "trigger_count": improvement.get("trigger_count", ""),
                    "trigger_precision_fp": improvement.get("trigger_precision_fp", ""),
                    "tp_removed": improvement.get("tp_removed", ""),
                    "fp_removed": improvement.get("fp_removed", ""),
                    "tp_removed_rate": improvement.get("tp_removed_rate", ""),
                    "fp_removed_rate": improvement.get("fp_removed_rate", ""),
                    "baseline_f1": baseline_f1,
                    "filtered_f1": filtered_f1,
                    "f1_delta": f1_delta,
                    "output_dir": row.get("output_dir", ""),
                }
            )

    print(f"[build_purity_aware_tables] Wrote {purity_qc_path}")
    print(f"[build_purity_aware_tables] Wrote {purity_rule_eval_path}")


if __name__ == "__main__":
    main()
