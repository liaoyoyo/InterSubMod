#!/usr/bin/env python3
"""snapshot_baseline.py - Capture current metrics into a versioned baseline.

Scans canonical output directories for metrics.json files and creates
a baseline snapshot that subsequent experiments compare against.

Usage:
    python3 scripts/validation/baseline/snapshot_baseline.py
    python3 scripts/validation/baseline/snapshot_baseline.py --output-root /custom/path
    python3 scripts/validation/baseline/snapshot_baseline.py --description "Post PON-only phasing"
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime


SAMPLES = [
    "HCC1395", "HCC1395_DORADO", "COLO829",
    "H1437", "H2009", "HCC1937", "HCC1954",
]

MODES = ["paired_pileup", "paired_full", "to_pileup", "to_full"]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASELINE_STORE = os.path.join(SCRIPT_DIR, "baseline_store.json")


def find_latest_run_dir(sample_root):
    """Find the most recent run directory (sorted lexicographically)."""
    if not os.path.isdir(sample_root):
        return None
    subdirs = [
        d for d in os.listdir(sample_root)
        if os.path.isdir(os.path.join(sample_root, d))
    ]
    if not subdirs:
        return None
    subdirs.sort()
    return os.path.join(sample_root, subdirs[-1])


def extract_metrics_from_run(run_dir):
    """Extract key metrics from a run's metrics.json.

    Checks both run_dir/metrics.json (complete_matrix runs) and
    run_dir/metrics/metrics.json (older pipeline runs).
    """
    metrics_path = os.path.join(run_dir, "metrics.json")
    if not os.path.isfile(metrics_path):
        metrics_path = os.path.join(run_dir, "metrics", "metrics.json")
    if not os.path.isfile(metrics_path):
        return None

    with open(metrics_path, "r") as f:
        data = json.load(f)

    baseline = data.get("baseline", {})
    filtered = data.get("filtered", {})
    improvement = data.get("improvement", {})

    return {
        "source_run": os.path.basename(run_dir),
        "source_path": metrics_path,
        # Baseline metrics (before ISM filter)
        "baseline_f1": baseline.get("f1"),
        "baseline_precision": baseline.get("precision"),
        "baseline_recall": baseline.get("recall"),
        "baseline_tp": baseline.get("tp"),
        "baseline_fp": baseline.get("fp"),
        # Filtered metrics (after ISM filter)
        "f1": filtered.get("f1"),
        "precision": filtered.get("precision"),
        "recall": filtered.get("recall"),
        "tp": filtered.get("tp"),
        "fp": filtered.get("fp"),
        # Improvement
        "f1_delta": improvement.get("f1_delta"),
        "fp_removed": improvement.get("fp_removed"),
        "tp_removed": improvement.get("tp_removed"),
        # Context
        "truth_total": data.get("truth_total"),
        "tp_regions_analyzed": data.get("tp_regions_analyzed"),
        "fp_regions_analyzed": data.get("fp_regions_analyzed"),
        "purity": data.get("purity"),
    }


def get_git_commit(project_root):
    """Get current git short hash."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=project_root, capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def get_git_branch(project_root):
    """Get current git branch name."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=project_root, capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def load_baseline_store(store_path):
    """Load existing baseline store or create empty one."""
    if os.path.isfile(store_path):
        with open(store_path, "r") as f:
            return json.load(f)
    return {"current_baseline_id": None, "baselines": {}}


def save_baseline_store(store_path, store):
    """Save baseline store to disk."""
    os.makedirs(os.path.dirname(store_path), exist_ok=True)
    with open(store_path, "w") as f:
        json.dump(store, f, indent=2, ensure_ascii=False)


def main():
    parser = argparse.ArgumentParser(description="Snapshot current metrics into versioned baseline")
    parser.add_argument("--output-root", default=None,
                        help="Override canonical output root (default: from config.sh)")
    parser.add_argument("--store", default=BASELINE_STORE,
                        help="Path to baseline_store.json")
    parser.add_argument("--description", default="",
                        help="Description for this baseline snapshot")
    parser.add_argument("--set-current", action="store_true", default=True,
                        help="Set this snapshot as the current baseline (default: True)")
    parser.add_argument("--no-set-current", dest="set_current", action="store_false",
                        help="Do not set this snapshot as the current baseline")
    args = parser.parse_args()

    # Determine project root and output root
    project_root = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
    output_root = args.output_root or os.path.join(
        "/big7_disk/liaoyoyo2001/big7_disk_output", "canonical"
    )

    # Generate baseline ID
    git_commit = get_git_commit(project_root)
    git_branch = get_git_branch(project_root)
    timestamp = datetime.now().isoformat(timespec="seconds")
    baseline_id = f"BL_{datetime.now().strftime('%Y%m%d')}_{git_commit}"

    print(f"[Baseline] Scanning canonical output: {output_root}")
    print(f"[Baseline] Git: {git_branch}@{git_commit}")

    # Scan all samples and modes
    samples_data = {}
    found_count = 0
    missing_count = 0

    for sample in SAMPLES:
        sample_data = {}
        for mode in MODES:
            sample_mode_root = os.path.join(output_root, sample, mode)
            latest_run = find_latest_run_dir(sample_mode_root)

            if latest_run is None:
                missing_count += 1
                continue

            metrics = extract_metrics_from_run(latest_run)
            if metrics is None:
                missing_count += 1
                print(f"  [MISS] {sample}/{mode}: no metrics.json in {os.path.basename(latest_run)}")
                continue

            sample_data[mode] = metrics
            found_count += 1
            f1_val = metrics.get("f1")
            f1_str = f"{f1_val:.4f}" if f1_val is not None else "N/A"
            print(f"  [OK]   {sample}/{mode}: F1={f1_str} (from {metrics['source_run']})")

        if sample_data:
            samples_data[sample] = sample_data

    print(f"\n[Baseline] Found {found_count} metrics, {missing_count} missing")

    if found_count == 0:
        print("[ERROR] No metrics found. Cannot create baseline.", file=sys.stderr)
        sys.exit(1)

    # Build baseline entry
    baseline_entry = {
        "timestamp": timestamp,
        "git_commit": git_commit,
        "git_branch": git_branch,
        "description": args.description or f"Snapshot from {output_root}",
        "metrics_found": found_count,
        "metrics_missing": missing_count,
        "samples": samples_data,
    }

    # Load and update store
    store = load_baseline_store(args.store)
    store["baselines"][baseline_id] = baseline_entry
    if args.set_current:
        store["current_baseline_id"] = baseline_id

    save_baseline_store(args.store, store)

    print(f"\n[Baseline] Saved baseline '{baseline_id}' to {args.store}")
    if args.set_current:
        print(f"[Baseline] Set as current baseline")

    # Print summary of primary metrics (paired_pileup F1)
    print("\n" + "=" * 60)
    print("  Baseline Summary (paired_pileup)")
    print("=" * 60)
    for sample in SAMPLES:
        if sample in samples_data and "paired_pileup" in samples_data[sample]:
            m = samples_data[sample]["paired_pileup"]
            f1 = m.get("f1")
            bf1 = m.get("baseline_f1")
            print(f"  {sample:20s}  baseline_F1={bf1}  filtered_F1={f1}")
        else:
            print(f"  {sample:20s}  (no data)")
    print("=" * 60)


if __name__ == "__main__":
    main()
