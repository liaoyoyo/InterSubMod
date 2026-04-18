#!/usr/bin/env python3
"""compute_metrics.py - Extract metrics from ISM experiment output.

Walks an experiment output directory, finds all metrics.json and
significance_summary.csv files, and produces a unified metrics bundle.

Usage:
    python3 scripts/validation/metrics/compute_metrics.py \
        --experiment-dir /path/to/experiment/output \
        --output metrics_bundle.json
"""

import argparse
import csv
import json
import os
import sys
from collections import Counter


def find_metrics_files(root_dir):
    """Find all metrics.json files recursively."""
    results = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "metrics.json" in filenames:
            results.append(os.path.join(dirpath, "metrics.json"))
    return sorted(results)


def find_significance_summaries(root_dir):
    """Find all significance_summary.csv files recursively."""
    results = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "significance_summary.csv" in filenames:
            results.append(os.path.join(dirpath, "significance_summary.csv"))
    return sorted(results)


def load_metrics_json(filepath):
    """Load and parse a metrics.json file."""
    with open(filepath, "r") as f:
        return json.load(f)


def extract_significance_stats(csv_path):
    """Extract distribution statistics from significance_summary.csv.

    Returns a dict with:
      - verification_class_dist: Counter of VerificationClass values
      - quality_score_stats: mean, median, std of Quality_Score
      - hpfine_ngroups_dist: Counter of HPFineNGroups values
      - dominant_label_dist: Counter of DominantLabel values
      - region_count: total regions in CSV
    """
    verification_classes = Counter()
    quality_scores = []
    hpfine_ngroups = Counter()
    dominant_labels = Counter()
    region_count = 0

    try:
        with open(csv_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                region_count += 1

                vc = row.get("VerificationClass", "Unknown")
                if vc:
                    verification_classes[vc] += 1

                qs = row.get("Quality_Score", "")
                if qs:
                    try:
                        quality_scores.append(float(qs))
                    except (ValueError, TypeError):
                        pass

                hpfn = row.get("HPFineNGroups", "")
                if hpfn:
                    try:
                        hpfine_ngroups[int(float(hpfn))] += 1
                    except (ValueError, TypeError):
                        pass

                dl = row.get("DominantLabel", "")
                if dl:
                    dominant_labels[dl] += 1
    except Exception as e:
        print(f"[WARN] Failed to parse {csv_path}: {e}", file=sys.stderr)
        return None

    # Compute quality score statistics
    qs_stats = {}
    if quality_scores:
        sorted_qs = sorted(quality_scores)
        n = len(sorted_qs)
        qs_mean = sum(sorted_qs) / n
        qs_median = sorted_qs[n // 2] if n % 2 else (sorted_qs[n // 2 - 1] + sorted_qs[n // 2]) / 2
        qs_std = (sum((x - qs_mean) ** 2 for x in sorted_qs) / n) ** 0.5
        qs_stats = {
            "mean": round(qs_mean, 4),
            "median": round(qs_median, 4),
            "std": round(qs_std, 4),
            "count": n,
        }

    return {
        "region_count": region_count,
        "verification_class_dist": dict(verification_classes),
        "quality_score_stats": qs_stats,
        "hpfine_ngroups_dist": {str(k): v for k, v in sorted(hpfine_ngroups.items())},
        "dominant_label_dist": dict(dominant_labels),
    }


def infer_sample_mode(metrics_data, filepath):
    """Infer sample and mode from metrics.json data or path."""
    sample = metrics_data.get("sample", "unknown")
    mode = metrics_data.get("mode", "unknown")
    return sample, mode


def build_metrics_bundle(experiment_dir):
    """Build a complete metrics bundle from experiment output."""
    metrics_files = find_metrics_files(experiment_dir)
    sig_files = find_significance_summaries(experiment_dir)

    bundle = {
        "experiment_dir": experiment_dir,
        "samples": {},
        "summary": {
            "total_metrics_files": len(metrics_files),
            "total_significance_files": len(sig_files),
        },
    }

    # Process metrics.json files
    for mf in metrics_files:
        data = load_metrics_json(mf)
        sample, mode = infer_sample_mode(data, mf)

        if sample not in bundle["samples"]:
            bundle["samples"][sample] = {}

        baseline = data.get("baseline", {})
        filtered = data.get("filtered", {})
        improvement = data.get("improvement", {})

        entry = {
            "source_file": mf,
            "baseline_f1": baseline.get("f1"),
            "baseline_precision": baseline.get("precision"),
            "baseline_recall": baseline.get("recall"),
            "baseline_tp": baseline.get("tp"),
            "baseline_fp": baseline.get("fp"),
            "f1": filtered.get("f1"),
            "precision": filtered.get("precision"),
            "recall": filtered.get("recall"),
            "tp": filtered.get("tp"),
            "fp": filtered.get("fp"),
            "f1_delta": improvement.get("f1_delta"),
            "fp_removed": improvement.get("fp_removed"),
            "tp_removed": improvement.get("tp_removed"),
            "truth_total": data.get("truth_total"),
            "tp_regions_analyzed": data.get("tp_regions_analyzed"),
            "fp_regions_analyzed": data.get("fp_regions_analyzed"),
        }

        bundle["samples"][sample][mode] = entry

    # Process significance_summary.csv files (match to TP/FP labels)
    for sf in sig_files:
        sig_stats = extract_significance_stats(sf)
        if sig_stats is None:
            continue

        # Try to match to a sample/mode by path
        # Path pattern: .../SAMPLE/MODE/RUN_ID/intersubmod_tp/...
        parts = sf.split(os.sep)
        label = None
        for i, p in enumerate(parts):
            if p.startswith("intersubmod_"):
                label = p  # intersubmod_tp or intersubmod_fp
                break

        if label:
            # Attach significance stats to the matching sample entry
            for sample in bundle["samples"]:
                for mode in bundle["samples"][sample]:
                    key = f"significance_{label}"
                    if key not in bundle["samples"][sample][mode]:
                        bundle["samples"][sample][mode][key] = sig_stats
                        break

    return bundle


def main():
    parser = argparse.ArgumentParser(description="Extract metrics from ISM experiment output")
    parser.add_argument("--experiment-dir", required=True,
                        help="Root directory of experiment output")
    parser.add_argument("--output", default=None,
                        help="Output path for metrics bundle JSON (default: stdout)")
    args = parser.parse_args()

    if not os.path.isdir(args.experiment_dir):
        print(f"[ERROR] Experiment directory not found: {args.experiment_dir}", file=sys.stderr)
        sys.exit(1)

    bundle = build_metrics_bundle(args.experiment_dir)

    # Print summary
    sample_count = len(bundle["samples"])
    print(f"[Metrics] Extracted metrics for {sample_count} samples from {args.experiment_dir}",
          file=sys.stderr)
    for sample, modes in sorted(bundle["samples"].items()):
        for mode, data in sorted(modes.items()):
            f1 = data.get("f1")
            f1_str = f"{f1:.4f}" if f1 is not None else "N/A"
            print(f"  {sample}/{mode}: F1={f1_str}", file=sys.stderr)

    # Output
    output_json = json.dumps(bundle, indent=2, ensure_ascii=False)
    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w") as f:
            f.write(output_json)
        print(f"[Metrics] Written to {args.output}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
