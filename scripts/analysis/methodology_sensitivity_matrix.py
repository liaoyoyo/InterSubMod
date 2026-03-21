#!/usr/bin/env python3
"""Aggregate run-level methodology sensitivity summaries."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List

from research_common import load_tsv_rows, read_json, to_float, write_tsv_rows


FIELDS = [
    "sample",
    "platform",
    "region_scope",
    "metric_family",
    "distance_method",
    "read_filter_profile",
    "min_reads",
    "window_bp",
    "methylation_mode",
    "cluster_method",
    "label_mode",
    "regions_tested",
    "strong_regions",
    "subclone_regions",
    "weak_regions",
    "noise_regions",
    "agreement_score",
    "f1",
    "tp",
    "fp",
    "fn",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", action="append", required=True, help="Round sample bundle directory")
    parser.add_argument("--output-tsv", required=True, help="Output methodology_sensitivity.tsv")
    return parser.parse_args()


def count_label_classes(rows: List[Dict[str, str]]) -> Dict[str, int]:
    counts = {"Strong": 0, "Subclone": 0, "Weak": 0, "Noise": 0}
    for row in rows:
        label = row.get("label_class", "")
        if label in counts:
            counts[label] += 1
    return counts


def main() -> None:
    args = parse_args()
    output_rows: List[Dict[str, str]] = []

    for run_dir_raw in args.run_dir:
        run_dir = Path(run_dir_raw).resolve()
        context_path = run_dir / "round_context.json"
        metrics_path = run_dir / "metrics.json"
        agreement_path = run_dir / "label_cluster_agreement.tsv"

        if not context_path.exists() or not metrics_path.exists() or not agreement_path.exists():
            continue

        context = read_json(context_path)
        metrics = read_json(metrics_path)
        agreement_rows = load_tsv_rows(agreement_path)
        class_counts = count_label_classes(agreement_rows)
        applicable = [
            row for row in agreement_rows if row.get("agreement_type") not in {"insufficient_label", ""}
        ]
        consistent = [
            row
            for row in agreement_rows
            if row.get("agreement_type") in {"consistent_strong", "consistent_subclone", "consistent_weak_or_noise"}
        ]
        agreement_score = (len(consistent) / len(applicable)) if applicable else 0.0

        filtered = metrics.get("filtered", {})
        truth_total = int(metrics.get("truth_total", 0))
        tp = int(filtered.get("tp", 0))

        output_rows.append(
            {
                "sample": context.get("sample", run_dir.name),
                "platform": context.get("platform", ""),
                "region_scope": context.get("region_scope", metrics.get("mode", "")),
                "metric_family": context.get("metric_family", "paired_pure_research"),
                "distance_method": context.get("distance_metric", ""),
                "read_filter_profile": context.get("read_filter_profile", "stable_default"),
                "min_reads": str(context.get("min_reads", "")),
                "window_bp": str(context.get("window_bp", "")),
                "methylation_mode": context.get("methylation_mode", "analog_full"),
                "cluster_method": context.get("cluster_method", "UPGMA"),
                "label_mode": context.get("label_mode", "hp+allele"),
                "regions_tested": str(len(agreement_rows)),
                "strong_regions": str(class_counts["Strong"]),
                "subclone_regions": str(class_counts["Subclone"]),
                "weak_regions": str(class_counts["Weak"]),
                "noise_regions": str(class_counts["Noise"]),
                "agreement_score": f"{agreement_score:.6f}",
                "f1": f"{to_float(filtered.get('f1'), 0.0):.4f}",
                "tp": str(tp),
                "fp": str(int(filtered.get("fp", 0))),
                "fn": str(max(truth_total - tp, 0)),
            }
        )

    write_tsv_rows(Path(args.output_tsv), FIELDS, output_rows)
    print(f"[methodology_sensitivity_matrix] Wrote {args.output_tsv}")


if __name__ == "__main__":
    main()
