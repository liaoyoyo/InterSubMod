#!/usr/bin/env python3
"""aggregate_results.py - Aggregate benchmark results across samples and modes.

Scans the benchmark output directory for metrics.json files and generates
a cross-sample summary report.

Usage:
    python3 aggregate_results.py --output-root /path/to/benchmark_pipeline
    python3 aggregate_results.py --output-root /path/to/benchmark_pipeline --output cross_sample_report.csv
"""

import argparse
import csv
import json
import os
import sys


def find_metrics_files(output_root):
    """Recursively find all metrics.json files under the output root."""
    metrics_files = []
    for dirpath, _, filenames in os.walk(output_root):
        for fname in filenames:
            if fname == "metrics.json":
                metrics_files.append(os.path.join(dirpath, fname))
    return sorted(metrics_files)


def load_metrics(filepath):
    """Load a metrics.json file and return its contents with the file path."""
    with open(filepath, "r") as f:
        data = json.load(f)
    data["_source"] = filepath
    # Extract date from directory path if possible
    parent = os.path.dirname(filepath)
    dirname = os.path.basename(parent)
    data["_run_dir"] = dirname
    return data


def generate_csv_report(all_metrics, output_path):
    """Generate a CSV report from all metrics."""
    fieldnames = [
        "sample", "mode", "purity", "run_dir",
        "truth_total",
        "baseline_tp", "baseline_fp", "baseline_precision", "baseline_recall", "baseline_f1",
        "filtered_tp", "filtered_fp", "filtered_precision", "filtered_recall", "filtered_f1",
        "f1_delta", "fp_removed", "tp_removed",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for m in all_metrics:
            row = {
                "sample": m.get("sample", ""),
                "mode": m.get("mode", ""),
                "purity": m.get("purity", ""),
                "run_dir": m.get("_run_dir", ""),
                "truth_total": m.get("truth_total", ""),
            }

            baseline = m.get("baseline", {})
            row["baseline_tp"] = baseline.get("tp", "")
            row["baseline_fp"] = baseline.get("fp", "")
            row["baseline_precision"] = baseline.get("precision", "")
            row["baseline_recall"] = baseline.get("recall", "")
            row["baseline_f1"] = baseline.get("f1", "")

            filtered = m.get("filtered", {})
            row["filtered_tp"] = filtered.get("tp", "")
            row["filtered_fp"] = filtered.get("fp", "")
            row["filtered_precision"] = filtered.get("precision", "")
            row["filtered_recall"] = filtered.get("recall", "")
            row["filtered_f1"] = filtered.get("f1", "")

            improvement = m.get("improvement", {})
            row["f1_delta"] = improvement.get("f1_delta", "")
            row["fp_removed"] = improvement.get("fp_removed", "")
            row["tp_removed"] = improvement.get("tp_removed", "")

            writer.writerow(row)


def print_summary_table(all_metrics):
    """Print a formatted summary table to stdout."""
    print()
    header = f"{'Sample':<12} {'Mode':<10} {'Purity':<8} {'Base F1':>8} {'Filt F1':>8} {'Delta':>8} {'FP Rm':>6} {'TP Rm':>6}"
    print(header)
    print("-" * len(header))

    for m in all_metrics:
        baseline = m.get("baseline", {})
        filtered = m.get("filtered", {})
        improvement = m.get("improvement", {})

        print(
            f"{m.get('sample', '?'):<12} "
            f"{m.get('mode', '?'):<10} "
            f"{m.get('purity', '?'):<8} "
            f"{baseline.get('f1', 0):>8.4f} "
            f"{filtered.get('f1', 0):>8.4f} "
            f"{improvement.get('f1_delta', 0):>+8.4f} "
            f"{improvement.get('fp_removed', 0):>6} "
            f"{improvement.get('tp_removed', 0):>6}"
        )

    print()


def main():
    parser = argparse.ArgumentParser(description="Aggregate benchmark results")
    parser.add_argument("--output-root", required=True, help="Root directory of benchmark outputs")
    parser.add_argument("--output", default=None, help="Output CSV path (default: output_root/cross_sample_report.csv)")
    args = parser.parse_args()

    metrics_files = find_metrics_files(args.output_root)

    if not metrics_files:
        print(f"[Aggregate] No metrics.json found under: {args.output_root}")
        sys.exit(0)

    print(f"[Aggregate] Found {len(metrics_files)} metrics file(s):")
    for f in metrics_files:
        print(f"  - {f}")

    all_metrics = [load_metrics(f) for f in metrics_files]

    # Generate CSV
    output_csv = args.output or os.path.join(args.output_root, "cross_sample_report.csv")
    generate_csv_report(all_metrics, output_csv)
    print(f"\n[Aggregate] CSV report written to: {output_csv}")

    # Print summary
    print_summary_table(all_metrics)


if __name__ == "__main__":
    main()
