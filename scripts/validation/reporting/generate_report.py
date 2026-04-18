#!/usr/bin/env python3
"""generate_report.py - Generate experiment report in Markdown.

Usage:
    python3 scripts/validation/reporting/generate_report.py \
        --experiment-id EXP_xxx \
        --hypothesis "Description" \
        --metrics-bundle metrics_bundle.json \
        --comparison comparison.json \
        --output experiment_report.md
"""

import argparse
import json
import os
import sys
from datetime import datetime


def load_json_safe(path):
    """Load JSON file, return None if not found."""
    if not path or not os.path.isfile(path):
        return None
    with open(path, "r") as f:
        return json.load(f)


def generate_report(experiment_id, hypothesis, hypothesis_id, git_commit, git_branch,
                    git_dirty, validation_mode, verdict, verdict_reason,
                    metrics_bundle, comparison):
    """Generate Markdown report string."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        f"<!--",
        f"建立時間: {timestamp}",
        f"目標: 自動化驗證實驗報告",
        f"處理範圍: {experiment_id}",
        f"-->",
        f"",
        f"# Experiment Report: {experiment_id}",
        f"",
        f"| Field | Value |",
        f"|-------|-------|",
        f"| Experiment ID | `{experiment_id}` |",
        f"| Hypothesis | {hypothesis} |",
    ]

    if hypothesis_id:
        lines.append(f"| Hypothesis ID | {hypothesis_id} |")

    lines.extend([
        f"| Git | `{git_branch}@{git_commit}` (dirty={git_dirty}) |",
        f"| Validation Mode | {validation_mode} |",
        f"| Timestamp | {timestamp} |",
        f"| **Verdict** | **{verdict.upper()}** |",
        f"| Reason | {verdict_reason} |",
        f"",
    ])

    # Per-sample results
    if comparison and "per_sample" in comparison:
        lines.extend([
            f"## Per-Sample Results",
            f"",
            f"| Sample | Baseline | Experiment | Delta | Status |",
            f"|--------|----------|------------|-------|--------|",
        ])

        per_sample = comparison["per_sample"]
        for sample in sorted(per_sample.keys()):
            data = per_sample[sample]
            bl = data.get("baseline", 0)
            ex = data.get("experiment", 0)
            delta = data.get("delta", 0)

            if delta > 0.001:
                status = "improved"
            elif delta < -0.005:
                status = "REGRESSION"
            elif delta < -0.001:
                status = "degraded"
            else:
                status = "unchanged"

            lines.append(
                f"| {sample} | {bl:.4f} | {ex:.4f} | {delta:+.4f} | {status} |"
            )

        lines.append("")

    # Cross-sample statistics
    if comparison and "cross_sample_stats" in comparison:
        stats = comparison["cross_sample_stats"]
        bootstrap = stats.get("bootstrap", {})
        t_test = stats.get("t_test", {})

        lines.extend([
            f"## Cross-Sample Statistics",
            f"",
            f"| Statistic | Value |",
            f"|-----------|-------|",
        ])

        if bootstrap.get("mean_delta") is not None:
            lines.append(f"| Mean delta | {bootstrap['mean_delta']:+.6f} |")
        if bootstrap.get("median_delta") is not None:
            lines.append(f"| Median delta | {bootstrap['median_delta']:+.6f} |")
        if bootstrap.get("ci_lower") is not None:
            lines.append(
                f"| Bootstrap 95% CI | [{bootstrap['ci_lower']:+.6f}, {bootstrap.get('ci_upper', 0):+.6f}] |"
            )
        if t_test.get("t_statistic") is not None:
            lines.append(f"| Paired t-statistic | {t_test['t_statistic']:.4f} |")
        if t_test.get("p_value") is not None:
            lines.append(f"| p-value | {t_test['p_value']:.6f} |")

        lines.extend([
            f"| Improved | {stats.get('n_improved', 0)} |",
            f"| Regressed | {stats.get('n_regressed', 0)} |",
            f"| Unchanged | {stats.get('n_unchanged', 0)} |",
            f"",
        ])

    # Metrics summary
    if metrics_bundle and "samples" in metrics_bundle:
        lines.extend([
            f"## Metrics Summary",
            f"",
            f"| Sample | Mode | F1 | Precision | Recall | TP Regions | FP Regions |",
            f"|--------|------|----|-----------|--------|------------|------------|",
        ])

        for sample in sorted(metrics_bundle["samples"].keys()):
            for mode in sorted(metrics_bundle["samples"][sample].keys()):
                data = metrics_bundle["samples"][sample][mode]
                if isinstance(data, dict) and "f1" in data:
                    f1 = data.get("f1")
                    p = data.get("precision")
                    r = data.get("recall")
                    tp_r = data.get("tp_regions_analyzed", "")
                    fp_r = data.get("fp_regions_analyzed", "")
                    f1_str = f"{f1:.4f}" if f1 is not None else "N/A"
                    p_str = f"{p:.4f}" if p is not None else "N/A"
                    r_str = f"{r:.4f}" if r is not None else "N/A"
                    lines.append(
                        f"| {sample} | {mode} | {f1_str} | {p_str} | {r_str} | {tp_r} | {fp_r} |"
                    )

        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate experiment report")
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--hypothesis", required=True)
    parser.add_argument("--hypothesis-id", default="")
    parser.add_argument("--git-commit", default="unknown")
    parser.add_argument("--git-branch", default="unknown")
    parser.add_argument("--git-dirty", default="unknown")
    parser.add_argument("--validation-mode", default="unknown")
    parser.add_argument("--verdict", default="unknown")
    parser.add_argument("--verdict-reason", default="")
    parser.add_argument("--metrics-bundle", default=None)
    parser.add_argument("--comparison", default=None)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    metrics_bundle = load_json_safe(args.metrics_bundle)
    comparison = load_json_safe(args.comparison)

    report = generate_report(
        experiment_id=args.experiment_id,
        hypothesis=args.hypothesis,
        hypothesis_id=args.hypothesis_id,
        git_commit=args.git_commit,
        git_branch=args.git_branch,
        git_dirty=args.git_dirty,
        validation_mode=args.validation_mode,
        verdict=args.verdict,
        verdict_reason=args.verdict_reason,
        metrics_bundle=metrics_bundle,
        comparison=comparison,
    )

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        f.write(report)

    print(f"[Report] Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
