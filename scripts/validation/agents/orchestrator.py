#!/usr/bin/env python3
"""orchestrator.py - Multi-experiment parallel orchestrator using git worktrees.

Reads experiment definitions from a YAML file and executes each in an isolated
git worktree. Each experiment runs validate.sh independently; results are
collected and compared in a final summary report.

Usage:
    # Run all experiments defined in YAML
    python3 scripts/validation/agents/orchestrator.py \
        --experiments scripts/validation/agents/experiments.yaml

    # Run specific experiments
    python3 scripts/validation/agents/orchestrator.py \
        --experiments experiments.yaml \
        --only normal_methyl_ref,cpg_stratification

    # Dry run (show what would execute)
    python3 scripts/validation/agents/orchestrator.py \
        --experiments experiments.yaml --dry-run
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime


def load_experiments(yaml_path):
    """Load experiment definitions from YAML file.

    Falls back to simple key-value parsing if PyYAML is not installed.
    """
    try:
        import yaml
        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f)
        return data.get("experiments", [])
    except ImportError:
        return _parse_yaml_fallback(yaml_path)


def _parse_yaml_fallback(yaml_path):
    """Minimal YAML parser for experiment definitions (no PyYAML dependency)."""
    experiments = []
    current = None

    with open(yaml_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            if stripped.startswith("- name:"):
                if current:
                    experiments.append(current)
                current = {"name": stripped.split(":", 1)[1].strip().strip('"').strip("'")}
            elif current and ":" in stripped:
                key, val = stripped.split(":", 1)
                key = key.strip()
                val = val.strip().strip('"').strip("'")
                # Type conversion
                if val.isdigit():
                    val = int(val)
                elif val in ("true", "True"):
                    val = True
                elif val in ("false", "False"):
                    val = False
                current[key] = val

    if current:
        experiments.append(current)

    return experiments


def get_project_root():
    """Get the project root directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, "..", "..", ".."))


def create_worktree(project_root, branch, worktree_dir):
    """Create a git worktree for an experiment branch.

    If the branch already exists, checks it out; otherwise creates from HEAD.
    """
    # Check if branch exists
    result = subprocess.run(
        ["git", "branch", "--list", branch],
        cwd=project_root, capture_output=True, text=True,
    )
    branch_exists = bool(result.stdout.strip())

    if branch_exists:
        cmd = ["git", "worktree", "add", worktree_dir, branch]
    else:
        cmd = ["git", "worktree", "add", "-b", branch, worktree_dir, "HEAD"]

    result = subprocess.run(cmd, cwd=project_root, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to create worktree: {result.stderr}")

    return worktree_dir


def remove_worktree(project_root, worktree_dir):
    """Remove a git worktree."""
    subprocess.run(
        ["git", "worktree", "remove", "--force", worktree_dir],
        cwd=project_root, capture_output=True, text=True,
    )
    # Prune stale worktrees
    subprocess.run(
        ["git", "worktree", "prune"],
        cwd=project_root, capture_output=True, text=True,
    )


def run_validate(worktree_dir, experiment, dry_run=False):
    """Run validate.sh in a worktree and return the result.

    Returns dict with:
        status: success | failed | skipped
        experiment_id: from validate.sh output
        output_dir: experiment output directory
        wall_time: seconds
        log: stdout+stderr
    """
    validate_script = os.path.join(worktree_dir, "scripts", "validation", "validate.sh")

    if not os.path.isfile(validate_script):
        return {
            "status": "failed",
            "error": f"validate.sh not found at {validate_script}",
        }

    cmd = [
        "bash", validate_script,
        "--hypothesis", experiment.get("hypothesis", "No hypothesis"),
        "--mode", experiment.get("mode", "quick"),
    ]

    # Optional parameters
    track = experiment.get("track")
    if track:
        cmd.extend(["--track", track])

    pipeline_mode = experiment.get("pipeline_mode")
    if pipeline_mode:
        cmd.extend(["--pipeline-mode", pipeline_mode])

    hypothesis_id = experiment.get("hypothesis_id")
    if hypothesis_id:
        cmd.extend(["--hypothesis-id", str(hypothesis_id)])

    samples = experiment.get("samples")
    if samples:
        cmd.extend(["--samples", samples])

    threads = experiment.get("threads")
    if threads:
        cmd.extend(["--threads", str(threads)])

    skip_build = experiment.get("skip_build", False)
    if skip_build:
        cmd.append("--skip-build")

    if dry_run:
        cmd.append("--dry-run")

    start_time = time.time()

    try:
        result = subprocess.run(
            cmd, cwd=worktree_dir,
            capture_output=True, text=True,
            timeout=experiment.get("timeout", 7200),  # 2 hour default
        )
        wall_time = time.time() - start_time

        # Parse experiment ID from output
        experiment_id = None
        output_dir = None
        verdict = None
        for line in result.stdout.splitlines():
            if line.startswith("EXPERIMENT_ID="):
                experiment_id = line.split("=", 1)[1]
            elif "VERDICT:" in line:
                verdict = line.split("VERDICT:", 1)[1].strip()

        # Check for experiment output directory
        for line in (result.stdout + result.stderr).splitlines():
            if "Output:" in line:
                parts = line.split("Output:", 1)
                if len(parts) > 1:
                    output_dir = parts[1].strip()

        return {
            "status": "success" if result.returncode == 0 else "failed",
            "return_code": result.returncode,
            "experiment_id": experiment_id,
            "output_dir": output_dir,
            "verdict": verdict,
            "wall_time": round(wall_time, 1),
            "log_tail": result.stderr[-2000:] if result.stderr else "",
        }

    except subprocess.TimeoutExpired:
        return {
            "status": "failed",
            "error": f"Timeout after {experiment.get('timeout', 7200)}s",
            "wall_time": time.time() - start_time,
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": str(e),
            "wall_time": time.time() - start_time,
        }


def run_experiment(project_root, experiment, worktree_base, dry_run=False, keep_worktrees=False):
    """Run a single experiment in an isolated worktree."""
    name = experiment["name"]
    branch = experiment.get("branch", f"experiment/{name}")

    worktree_dir = os.path.join(worktree_base, name)

    print(f"\n{'=' * 60}")
    print(f"  EXPERIMENT: {name}")
    print(f"  Branch: {branch}")
    print(f"  Hypothesis: {experiment.get('hypothesis', 'N/A')}")
    print(f"  Mode: {experiment.get('mode', 'quick')} | Track: {experiment.get('track', 'paired')}")
    print(f"{'=' * 60}")

    result = {"name": name, "branch": branch}

    try:
        # Create worktree
        print(f"  Creating worktree: {worktree_dir}")
        if not dry_run:
            create_worktree(project_root, branch, worktree_dir)
        result["worktree"] = worktree_dir

        # Run validation
        print(f"  Running validate.sh...")
        validate_result = run_validate(worktree_dir, experiment, dry_run=dry_run)
        result.update(validate_result)

        status_icon = "OK" if result.get("status") == "success" else "FAIL"
        print(f"  [{status_icon}] {name}: {result.get('verdict', result.get('status', 'unknown'))}")
        if result.get("wall_time"):
            print(f"  Wall time: {result['wall_time']:.0f}s")

    except Exception as e:
        result["status"] = "failed"
        result["error"] = str(e)
        print(f"  [ERROR] {name}: {e}")

    finally:
        if not keep_worktrees and not dry_run and os.path.isdir(worktree_dir):
            print(f"  Cleaning up worktree...")
            remove_worktree(project_root, worktree_dir)

    return result


def generate_comparison_report(results, output_path):
    """Generate a comparison report across all experiments."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines = [
        f"<!--",
        f"建立時間: {timestamp}",
        f"目標: 多實驗平行比較報告",
        f"-->",
        f"",
        f"# Multi-Experiment Comparison Report",
        f"",
        f"Generated: {timestamp}",
        f"",
        f"## Summary",
        f"",
        f"| Experiment | Branch | Mode | Track | Verdict | Wall Time |",
        f"|------------|--------|------|-------|---------|-----------|",
    ]

    for r in results:
        verdict = r.get("verdict", r.get("status", "N/A"))
        wall = f"{r.get('wall_time', 0):.0f}s" if r.get("wall_time") else "N/A"
        lines.append(
            f"| {r['name']} | `{r.get('branch', 'N/A')}` "
            f"| {r.get('mode', 'N/A')} | {r.get('track', 'paired')} "
            f"| **{verdict}** | {wall} |"
        )

    lines.extend([
        "",
        "## Per-Experiment Details",
        "",
    ])

    for r in results:
        lines.extend([
            f"### {r['name']}",
            f"",
            f"- **Branch**: `{r.get('branch', 'N/A')}`",
            f"- **Status**: {r.get('status', 'N/A')}",
            f"- **Verdict**: {r.get('verdict', 'N/A')}",
            f"- **Experiment ID**: `{r.get('experiment_id', 'N/A')}`",
        ])

        if r.get("error"):
            lines.append(f"- **Error**: {r['error']}")
        if r.get("output_dir"):
            lines.append(f"- **Output**: `{r['output_dir']}`")

        lines.append("")

    report_text = "\n".join(lines)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        f.write(report_text)

    return report_text


def main():
    parser = argparse.ArgumentParser(
        description="Orchestrate parallel experiments using git worktrees",
    )
    parser.add_argument(
        "--experiments", required=True,
        help="Path to experiments YAML definition file",
    )
    parser.add_argument(
        "--only", default=None,
        help="Comma-separated list of experiment names to run (default: all)",
    )
    parser.add_argument(
        "--worktree-base", default=None,
        help="Base directory for worktrees (default: /tmp/ism_worktrees_<timestamp>)",
    )
    parser.add_argument(
        "--report-output", default=None,
        help="Output path for comparison report (default: experiments/<timestamp>_comparison.md)",
    )
    parser.add_argument(
        "--keep-worktrees", action="store_true",
        help="Do not remove worktrees after completion",
    )
    parser.add_argument(
        "--sequential", action="store_true",
        help="Run experiments sequentially (default: parallel via subprocesses)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing",
    )
    args = parser.parse_args()

    project_root = get_project_root()
    experiments = load_experiments(args.experiments)

    if not experiments:
        print("[ERROR] No experiments found in YAML file", file=sys.stderr)
        sys.exit(1)

    # Filter experiments
    if args.only:
        selected = set(args.only.split(","))
        experiments = [e for e in experiments if e["name"] in selected]
        if not experiments:
            print(f"[ERROR] No experiments matched --only filter: {args.only}", file=sys.stderr)
            sys.exit(1)

    # Sort by priority
    experiments.sort(key=lambda e: e.get("priority", 99))

    # Setup worktree base
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    worktree_base = args.worktree_base or f"/tmp/ism_worktrees_{ts}"
    os.makedirs(worktree_base, exist_ok=True)

    print(f"[Orchestrator] Project root: {project_root}")
    print(f"[Orchestrator] Worktree base: {worktree_base}")
    print(f"[Orchestrator] Experiments: {len(experiments)}")
    print(f"[Orchestrator] Mode: {'sequential' if args.sequential else 'sequential (parallel planned)'}")

    # Run experiments
    # Note: true parallel execution (multiprocessing) is planned for Wave 5.
    # Currently runs sequentially to avoid resource contention on shared disk.
    results = []
    for exp in experiments:
        result = run_experiment(
            project_root, exp, worktree_base,
            dry_run=args.dry_run,
            keep_worktrees=args.keep_worktrees,
        )
        # Carry over experiment metadata for report
        result["mode"] = exp.get("mode", "quick")
        result["track"] = exp.get("track", "paired")
        results.append(result)

    # Generate comparison report
    report_output = args.report_output or os.path.join(
        project_root, "docs", "experiments", "outputs",
        f"{ts}_orchestrator_comparison.md",
    )
    report = generate_comparison_report(results, report_output)

    # Final summary
    total = len(results)
    succeeded = sum(1 for r in results if r.get("status") == "success")
    positive = sum(1 for r in results if r.get("verdict", "").startswith("positive"))

    print(f"\n{'=' * 60}")
    print(f"  ORCHESTRATOR COMPLETE")
    print(f"  Total: {total} | Succeeded: {succeeded} | Positive: {positive}")
    print(f"  Report: {report_output}")
    print(f"{'=' * 60}")

    # Dump results JSON
    results_json_path = report_output.replace(".md", ".json")
    with open(results_json_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"  Results JSON: {results_json_path}")

    # Cleanup worktree base if empty
    if not args.keep_worktrees and os.path.isdir(worktree_base):
        try:
            os.rmdir(worktree_base)
        except OSError:
            pass  # Not empty, keep it


if __name__ == "__main__":
    main()
