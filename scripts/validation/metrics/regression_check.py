#!/usr/bin/env python3
"""regression_check.py - Per-sample regression detection.

Checks each sample's F1 against the baseline and flags regressions
exceeding the threshold. Can be used standalone or as part of validate.sh.

Usage:
    python3 scripts/validation/metrics/regression_check.py \
        --comparison comparison.json

    python3 scripts/validation/metrics/regression_check.py \
        --baseline baseline_store.json \
        --experiment metrics_bundle.json \
        --threshold 0.005
"""

import argparse
import json
import os
import sys


def check_regressions_from_comparison(comparison, threshold=0.005):
    """Check for regressions in a pre-computed comparison result."""
    per_sample = comparison.get("per_sample", {})
    regressions = []

    for sample, data in sorted(per_sample.items()):
        delta = data.get("delta", 0)
        if delta < -threshold:
            regressions.append({
                "sample": sample,
                "baseline": data.get("baseline"),
                "experiment": data.get("experiment"),
                "delta": delta,
                "severity": "critical" if delta < -2 * threshold else "warning",
            })

    return regressions


def check_regressions_direct(baseline_path, experiment_path, mode="paired_pileup",
                             metric="f1", threshold=0.005):
    """Check regressions by directly comparing baseline store and experiment bundle."""
    with open(baseline_path, "r") as f:
        store = json.load(f)
    with open(experiment_path, "r") as f:
        bundle = json.load(f)

    current_id = store.get("current_baseline_id")
    if not current_id:
        return []

    baseline_entry = store.get("baselines", {}).get(current_id, {})
    baseline_samples = baseline_entry.get("samples", {})
    experiment_samples = bundle.get("samples", {})

    regressions = []
    for sample in sorted(set(baseline_samples.keys()) & set(experiment_samples.keys())):
        bl_val = baseline_samples[sample].get(mode, {}).get(metric)
        ex_val = experiment_samples[sample].get(mode, {}).get(metric)

        if bl_val is None or ex_val is None:
            continue

        delta = float(ex_val) - float(bl_val)
        if delta < -threshold:
            regressions.append({
                "sample": sample,
                "baseline": float(bl_val),
                "experiment": float(ex_val),
                "delta": round(delta, 6),
                "severity": "critical" if delta < -2 * threshold else "warning",
            })

    return regressions


def format_regression_report(regressions, threshold):
    """Format regression report as human-readable text."""
    if not regressions:
        return f"No regressions detected (threshold: {threshold})"

    lines = [
        f"REGRESSION ALERT: {len(regressions)} sample(s) exceeded threshold ({threshold})",
        "",
    ]
    for r in regressions:
        severity_marker = "!!!" if r["severity"] == "critical" else " ! "
        lines.append(
            f"  {severity_marker} {r['sample']:20s}  "
            f"{r['baseline']:.4f} -> {r['experiment']:.4f}  "
            f"(delta: {r['delta']:+.4f})"
        )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Per-sample regression detection")
    parser.add_argument("--comparison", default=None,
                        help="Path to comparison.json (pre-computed)")
    parser.add_argument("--baseline", default=None,
                        help="Path to baseline_store.json (direct mode)")
    parser.add_argument("--experiment", default=None,
                        help="Path to experiment metrics_bundle.json (direct mode)")
    parser.add_argument("--mode", default="paired_pileup", help="Pipeline mode")
    parser.add_argument("--metric", default="f1", help="Metric to check")
    parser.add_argument("--threshold", type=float, default=0.005,
                        help="Regression threshold (default: 0.005)")
    parser.add_argument("--output", default=None, help="Output JSON path")
    args = parser.parse_args()

    if args.comparison:
        with open(args.comparison, "r") as f:
            comparison = json.load(f)
        regressions = check_regressions_from_comparison(comparison, args.threshold)
    elif args.baseline and args.experiment:
        regressions = check_regressions_direct(
            args.baseline, args.experiment,
            mode=args.mode, metric=args.metric, threshold=args.threshold,
        )
    else:
        print("[ERROR] Provide either --comparison or both --baseline and --experiment",
              file=sys.stderr)
        sys.exit(1)

    # Print report
    report = format_regression_report(regressions, args.threshold)
    print(report, file=sys.stderr)

    # Output JSON
    result = {
        "threshold": args.threshold,
        "regression_count": len(regressions),
        "regressions": regressions,
        "passed": len(regressions) == 0,
    }

    output_json = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w") as f:
            f.write(output_json)
    else:
        print(output_json)

    # Exit code: non-zero if regressions found
    sys.exit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
