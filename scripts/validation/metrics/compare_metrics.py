#!/usr/bin/env python3
"""compare_metrics.py - Statistical comparison between experiment and baseline.

Computes per-sample deltas, cross-sample paired t-test, bootstrap CI,
and Go/No-Go verdict.

Usage:
    python3 scripts/validation/metrics/compare_metrics.py \
        --baseline baseline_store.json \
        --experiment metrics_bundle.json \
        --output comparison.json

    python3 scripts/validation/metrics/compare_metrics.py \
        --baseline baseline_store.json \
        --experiment metrics_bundle.json \
        --mode paired_pileup \
        --metric f1
"""

import argparse
import json
import os
import random
import sys


def load_json(filepath):
    with open(filepath, "r") as f:
        return json.load(f)


def get_current_baseline(store):
    """Get the current baseline entry from baseline_store.json."""
    current_id = store.get("current_baseline_id")
    if not current_id:
        return None, None
    entry = store.get("baselines", {}).get(current_id)
    return current_id, entry


def extract_metric_pairs(baseline_entry, experiment_bundle, mode, metric):
    """Extract paired (baseline, experiment) metric values for each sample.

    Returns list of (sample_name, baseline_value, experiment_value) tuples.
    """
    pairs = []
    baseline_samples = baseline_entry.get("samples", {})
    experiment_samples = experiment_bundle.get("samples", {})

    for sample in sorted(set(baseline_samples.keys()) & set(experiment_samples.keys())):
        bl_mode = baseline_samples[sample].get(mode, {})
        ex_mode = experiment_samples[sample].get(mode, {})

        bl_val = bl_mode.get(metric)
        ex_val = ex_mode.get(metric)

        if bl_val is not None and ex_val is not None:
            pairs.append((sample, float(bl_val), float(ex_val)))

    return pairs


def paired_t_test(values_a, values_b):
    """Compute paired t-test statistic and p-value.

    Uses scipy if available, otherwise falls back to manual computation.
    """
    n = len(values_a)
    if n < 2:
        return {"t_statistic": None, "p_value": None, "n": n}

    try:
        from scipy import stats
        t_stat, p_val = stats.ttest_rel(values_b, values_a)
        return {
            "t_statistic": round(float(t_stat), 6),
            "p_value": round(float(p_val), 6),
            "n": n,
        }
    except ImportError:
        # Manual paired t-test
        diffs = [b - a for a, b in zip(values_a, values_b)]
        mean_d = sum(diffs) / n
        var_d = sum((d - mean_d) ** 2 for d in diffs) / (n - 1) if n > 1 else 0
        se_d = (var_d / n) ** 0.5 if var_d > 0 else 0

        t_stat = mean_d / se_d if se_d > 0 else 0

        # Approximate p-value using normal distribution for n >= 7
        import math
        if n >= 7 and se_d > 0:
            # Two-tailed p-value approximation
            z = abs(t_stat)
            p_val = 2 * (1 - 0.5 * (1 + math.erf(z / math.sqrt(2))))
        else:
            p_val = None

        return {
            "t_statistic": round(t_stat, 6),
            "p_value": round(p_val, 6) if p_val is not None else None,
            "n": n,
        }


def bootstrap_ci(values_a, values_b, n_iterations=1000, ci_level=0.95, seed=42):
    """Compute bootstrap confidence interval for mean difference.

    Returns (ci_lower, ci_upper, mean_delta).
    """
    n = len(values_a)
    if n < 2:
        return {"ci_lower": None, "ci_upper": None, "mean_delta": None}

    diffs = [b - a for a, b in zip(values_a, values_b)]
    mean_delta = sum(diffs) / n

    rng = random.Random(seed)
    boot_means = []
    for _ in range(n_iterations):
        boot_sample = [rng.choice(diffs) for _ in range(n)]
        boot_means.append(sum(boot_sample) / n)

    boot_means.sort()
    alpha = 1 - ci_level
    lower_idx = int(alpha / 2 * n_iterations)
    upper_idx = int((1 - alpha / 2) * n_iterations) - 1

    return {
        "ci_lower": round(boot_means[lower_idx], 6),
        "ci_upper": round(boot_means[upper_idx], 6),
        "mean_delta": round(mean_delta, 6),
        "median_delta": round(sorted(diffs)[n // 2], 6),
        "n_iterations": n_iterations,
        "ci_level": ci_level,
    }


def determine_verdict(comparison, f1_improvement_threshold=0.001, f1_regression_threshold=0.005,
                      significance_alpha=0.05, validation_mode="full"):
    """Determine Go/No-Go verdict based on comparison results.

    Verdicts:
        positive:       Statistically significant improvement, no regressions
        positive_pilot: Quick mode positive (needs full confirmation)
        neutral:        No meaningful change
        negative:       Regression detected or significant degradation
    """
    per_sample = comparison.get("per_sample", {})
    cross_sample = comparison.get("cross_sample_stats", {})
    regressions = comparison.get("regressions", [])

    mean_delta = cross_sample.get("bootstrap", {}).get("mean_delta")
    ci_lower = cross_sample.get("bootstrap", {}).get("ci_lower")
    p_value = cross_sample.get("t_test", {}).get("p_value")

    # Regression check
    if regressions:
        return "negative", f"Regression in {len(regressions)} sample(s): {', '.join(regressions)}"

    if mean_delta is None:
        return "neutral", "Insufficient data for comparison"

    # Negative: clear degradation
    if mean_delta < -f1_improvement_threshold:
        return "negative", f"Mean delta F1 = {mean_delta:+.6f} (degradation)"

    # Neutral: no meaningful change
    if abs(mean_delta) <= f1_improvement_threshold:
        return "neutral", f"Mean delta F1 = {mean_delta:+.6f} (within noise)"

    # Positive checks
    if validation_mode == "quick":
        return "positive_pilot", f"Mean delta F1 = {mean_delta:+.6f} (quick mode — needs full confirmation)"

    # Full mode: require bootstrap CI lower > 0
    if ci_lower is not None and ci_lower > 0:
        sig_note = ""
        if p_value is not None and p_value < significance_alpha:
            sig_note = f", p={p_value:.4f}"
        return "positive", f"Mean delta F1 = {mean_delta:+.6f}, CI=[{ci_lower:+.6f}, {cross_sample.get('bootstrap', {}).get('ci_upper', 0):+.6f}]{sig_note}"

    # Improvement but CI includes 0
    return "neutral", f"Mean delta F1 = {mean_delta:+.6f} but CI includes 0 (not significant)"


def compare(baseline_store_path, experiment_bundle_path, mode="paired_pileup",
            metric="f1", n_bootstrap=1000, ci_level=0.95,
            f1_improvement_threshold=0.001, f1_regression_threshold=0.005,
            significance_alpha=0.05, validation_mode="full"):
    """Run full comparison and return structured result."""

    store = load_json(baseline_store_path)
    bundle = load_json(experiment_bundle_path)

    baseline_id, baseline_entry = get_current_baseline(store)
    if baseline_entry is None:
        return {"error": "No current baseline found in store"}

    # Extract paired values
    pairs = extract_metric_pairs(baseline_entry, bundle, mode, metric)

    if not pairs:
        return {
            "error": f"No matching samples found for mode={mode}, metric={metric}",
            "baseline_id": baseline_id,
        }

    # Per-sample deltas
    per_sample = {}
    regressions = []
    for sample, bl_val, ex_val in pairs:
        delta = ex_val - bl_val
        per_sample[sample] = {
            "baseline": round(bl_val, 6),
            "experiment": round(ex_val, 6),
            "delta": round(delta, 6),
        }
        if delta < -f1_regression_threshold:
            regressions.append(sample)

    # Cross-sample statistics
    bl_values = [p[1] for p in pairs]
    ex_values = [p[2] for p in pairs]

    t_test_result = paired_t_test(bl_values, ex_values)
    bootstrap_result = bootstrap_ci(bl_values, ex_values, n_bootstrap, ci_level)

    # Count improvements / regressions
    deltas = [ex - bl for _, bl, ex in pairs]
    n_improved = sum(1 for d in deltas if d > f1_improvement_threshold)
    n_regressed = sum(1 for d in deltas if d < -f1_improvement_threshold)
    n_unchanged = len(deltas) - n_improved - n_regressed

    cross_sample_stats = {
        "t_test": t_test_result,
        "bootstrap": bootstrap_result,
        "n_improved": n_improved,
        "n_regressed": n_regressed,
        "n_unchanged": n_unchanged,
    }

    comparison = {
        "baseline_id": baseline_id,
        "mode": mode,
        "metric": metric,
        "n_samples": len(pairs),
        "per_sample": per_sample,
        "cross_sample_stats": cross_sample_stats,
        "regressions": regressions,
    }

    # Verdict
    verdict, reason = determine_verdict(
        comparison,
        f1_improvement_threshold=f1_improvement_threshold,
        f1_regression_threshold=f1_regression_threshold,
        significance_alpha=significance_alpha,
        validation_mode=validation_mode,
    )
    comparison["verdict"] = verdict
    comparison["verdict_reason"] = reason

    return comparison


def main():
    parser = argparse.ArgumentParser(description="Compare experiment metrics against baseline")
    parser.add_argument("--baseline", required=True, help="Path to baseline_store.json")
    parser.add_argument("--experiment", required=True, help="Path to experiment metrics_bundle.json")
    parser.add_argument("--mode", default="paired_pileup", help="Pipeline mode to compare")
    parser.add_argument("--metric", default="f1", help="Metric to compare (default: f1)")
    parser.add_argument("--validation-mode", default="full", choices=["quick", "full"],
                        help="Validation mode (affects verdict logic)")
    parser.add_argument("--bootstrap-n", type=int, default=1000, help="Bootstrap iterations")
    parser.add_argument("--ci-level", type=float, default=0.95, help="CI confidence level")
    parser.add_argument("--f1-improvement", type=float, default=0.001, help="F1 improvement threshold")
    parser.add_argument("--f1-regression", type=float, default=0.005, help="F1 regression threshold")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance alpha")
    parser.add_argument("--output", default=None, help="Output path (default: stdout)")
    args = parser.parse_args()

    result = compare(
        baseline_store_path=args.baseline,
        experiment_bundle_path=args.experiment,
        mode=args.mode,
        metric=args.metric,
        n_bootstrap=args.bootstrap_n,
        ci_level=args.ci_level,
        f1_improvement_threshold=args.f1_improvement,
        f1_regression_threshold=args.f1_regression,
        significance_alpha=args.alpha,
        validation_mode=args.validation_mode,
    )

    # Print verdict to stderr
    verdict = result.get("verdict", "unknown")
    reason = result.get("verdict_reason", "")
    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  VERDICT: {verdict.upper()}", file=sys.stderr)
    print(f"  Reason:  {reason}", file=sys.stderr)
    if result.get("per_sample"):
        print(f"\n  Per-sample {args.metric} deltas:", file=sys.stderr)
        for sample, data in sorted(result["per_sample"].items()):
            d = data["delta"]
            marker = "!!" if d < -args.f1_regression else ("++" if d > args.f1_improvement else "  ")
            print(f"    {marker} {sample:20s}  {data['baseline']:.4f} -> {data['experiment']:.4f}  ({d:+.4f})",
                  file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)

    # Output JSON
    output_json = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w") as f:
            f.write(output_json)
        print(f"[Compare] Written to {args.output}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
