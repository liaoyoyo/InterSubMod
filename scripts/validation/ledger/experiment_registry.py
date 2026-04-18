#!/usr/bin/env python3
"""experiment_registry.py - Record experiments to the evidence ledger.

Appends a validation_experiment entry to the evidence_ledger.jsonl file.
Extends the existing ledger schema with validation-specific fields.

Usage:
    python3 scripts/validation/ledger/experiment_registry.py \
        --experiment-id EXP_xxx \
        --hypothesis "Description" \
        --verdict positive \
        --metrics-bundle metrics_bundle.json \
        --comparison comparison.json
"""

import argparse
import json
import os
import sys
from datetime import datetime


def load_json_safe(path):
    """Load JSON file, return empty dict if not found."""
    if not path or not os.path.isfile(path):
        return {}
    with open(path, "r") as f:
        return json.load(f)


def build_ledger_entry(experiment_id, hypothesis, hypothesis_id, git_commit, git_branch,
                       validation_mode, verdict, wall_time, metrics_bundle, comparison):
    """Build a ledger entry compatible with existing evidence_ledger.jsonl schema."""
    entry = {
        "cycle_id": experiment_id,
        "type": "validation_experiment",
        "hypothesis_id": hypothesis_id or experiment_id,
        "hypothesis": hypothesis,
        "git_commit": git_commit,
        "git_branch": git_branch,
        "validation_mode": validation_mode,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "wall_time_seconds": int(wall_time) if wall_time else 0,
    }

    # Extract sample list and F1 values from metrics bundle
    if metrics_bundle and "samples" in metrics_bundle:
        datasets = sorted(metrics_bundle["samples"].keys())
        entry["datasets_tested"] = datasets

        result_f1 = {}
        for sample in datasets:
            for mode, data in metrics_bundle["samples"][sample].items():
                if isinstance(data, dict):
                    f1 = data.get("f1")
                    if f1 is not None:
                        result_f1[sample] = f1
        entry["result_f1"] = result_f1

    # Extract comparison results
    if comparison:
        entry["baseline_id"] = comparison.get("baseline_id")

        if "per_sample" in comparison:
            baseline_f1 = {}
            delta_f1 = {}
            for sample, data in comparison["per_sample"].items():
                if "baseline" in data:
                    baseline_f1[sample] = data["baseline"]
                if "delta" in data:
                    delta_f1[sample] = data["delta"]
            entry["baseline_f1"] = baseline_f1
            entry["delta_f1"] = delta_f1

        if "cross_sample_stats" in comparison:
            stats = comparison["cross_sample_stats"]
            entry["cross_sample_stats"] = {
                "mean_delta_f1": stats.get("bootstrap", {}).get("mean_delta"),
                "paired_t_pvalue": stats.get("t_test", {}).get("p_value"),
                "bootstrap_ci": [
                    stats.get("bootstrap", {}).get("ci_lower"),
                    stats.get("bootstrap", {}).get("ci_upper"),
                ],
            }

        entry["regressions"] = comparison.get("regressions", [])

    entry["verdict"] = verdict
    entry["human_decision"] = "pending"
    entry["artifacts_path"] = f"experiments/{experiment_id}/"

    return entry


def append_to_ledger(ledger_path, entry):
    """Append entry to JSONL ledger file."""
    os.makedirs(os.path.dirname(ledger_path) or ".", exist_ok=True)
    with open(ledger_path, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Record experiment to evidence ledger")
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--hypothesis", required=True)
    parser.add_argument("--hypothesis-id", default="")
    parser.add_argument("--git-commit", default="unknown")
    parser.add_argument("--git-branch", default="unknown")
    parser.add_argument("--validation-mode", default="unknown")
    parser.add_argument("--verdict", default="unknown")
    parser.add_argument("--wall-time", default="0")
    parser.add_argument("--metrics-bundle", default=None)
    parser.add_argument("--comparison", default=None)
    parser.add_argument("--ledger", required=True, help="Path to evidence_ledger.jsonl")
    args = parser.parse_args()

    metrics_bundle = load_json_safe(args.metrics_bundle)
    comparison = load_json_safe(args.comparison)

    entry = build_ledger_entry(
        experiment_id=args.experiment_id,
        hypothesis=args.hypothesis,
        hypothesis_id=args.hypothesis_id,
        git_commit=args.git_commit,
        git_branch=args.git_branch,
        validation_mode=args.validation_mode,
        verdict=args.verdict,
        wall_time=args.wall_time,
        metrics_bundle=metrics_bundle,
        comparison=comparison,
    )

    append_to_ledger(args.ledger, entry)
    print(f"[Ledger] Appended experiment {args.experiment_id} to {args.ledger}", file=sys.stderr)
    print(f"[Ledger] Verdict: {args.verdict}", file=sys.stderr)


if __name__ == "__main__":
    main()
