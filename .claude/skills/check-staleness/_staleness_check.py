#!/usr/bin/env python3
"""
Staleness checker for Resilient Waterfall harness P2 PRECHECK gate.

Reads `state/cycles/{cycle_id}/plan.json`, performs three freshness checks
(binary, dataset, upstream reports), and writes `precheck.json` per
`state/schemas/precheck.schema.json`.

Usage:
    python3 _staleness_check.py <cycle_id>

Exit codes:
    0 = PASS or WARN (cycle may advance, possibly with caveats)
    1 = BLOCKED (cycle must not advance without user override)
    2 = error (plan.json missing, etc.)
"""

import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional


REPO_ROOT = Path(__file__).resolve().parents[3]  # .claude/skills/check-staleness/_staleness_check.py → InterSubMod/
STATE_ROOT = REPO_ROOT / "state"
INVALIDATION_DIR = STATE_ROOT / "invalidation"


# Known dataset probes — extend as new pitfalls are discovered.
# Each probe returns dict {status, schema_violations[]}
DATASET_PROBES: dict[str, dict[str, Any]] = {
    "must_have_columns": ["caller_af"],
    "min_loh_coverage_pct": 50.0,
}


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_plan(cycle_id: str) -> dict:
    plan_path = STATE_ROOT / "cycles" / cycle_id / "plan.json"
    if not plan_path.is_file():
        sys.stderr.write(f"ERROR: plan.json not found at {plan_path}\n")
        sys.exit(2)
    return json.loads(plan_path.read_text())


def git_head_sha() -> Optional[str]:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT, capture_output=True, text=True, check=True
        )
        return out.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def git_commit_distance(from_sha: str, to_sha: str) -> Optional[int]:
    """Number of commits in `from_sha..to_sha`. None on error."""
    try:
        out = subprocess.run(
            ["git", "rev-list", "--count", f"{from_sha}..{to_sha}"],
            cwd=REPO_ROOT, capture_output=True, text=True, check=True
        )
        return int(out.stdout.strip())
    except (subprocess.CalledProcessError, ValueError, FileNotFoundError):
        return None


def git_sha_exists(sha: str) -> bool:
    try:
        subprocess.run(
            ["git", "cat-file", "-e", sha],
            cwd=REPO_ROOT, capture_output=True, check=True
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def check_binary(plan: dict) -> dict:
    expected = plan.get("preconditions", {}).get("binary_version")
    if expected is None:
        return {"status": "skipped", "expected": None, "actual_head": None,
                "stale_distance": None, "note": "Pure-analysis cycle, no binary dependency."}

    head = git_head_sha()
    if head is None:
        return {"status": "missing", "expected": expected, "actual_head": None,
                "stale_distance": None, "note": "Could not query git HEAD."}

    if not git_sha_exists(expected):
        return {"status": "missing", "expected": expected, "actual_head": head,
                "stale_distance": None, "note": "Expected SHA not found in git history."}

    if expected == head:
        return {"status": "fresh", "expected": expected, "actual_head": head,
                "stale_distance": 0, "note": "Binary at HEAD."}

    distance = git_commit_distance(expected, head)
    return {"status": "stale", "expected": expected, "actual_head": head,
            "stale_distance": distance,
            "note": f"HEAD is {distance} commits ahead of plan-stated binary."}


def check_dataset(plan: dict) -> dict:
    dataset_id = plan.get("preconditions", {}).get("dataset_id")
    if dataset_id is None:
        return {"status": "skipped", "expected_id": None,
                "schema_violations": [],
                "note": "Pure-binary cycle, no dataset dependency."}

    # Path A: probe is structural (schema sanity); actual data load happens in P3 pilot.
    # Future Path B: integrate with LlamaIndex dataset_index for richer probes.
    violations = []

    # Stale-mark check via invalidation/stale_marks.jsonl
    stale_marks = read_jsonl(INVALIDATION_DIR / "stale_marks.jsonl")
    for entry in stale_marks:
        if entry.get("type") == "dataset" and entry.get("dataset_id") == dataset_id:
            violations.append(f"Dataset marked stale: {entry.get('reason', 'unspecified')}")

    # Heuristic name-based probes (e.g. "merged" datasets known to lack caller_af)
    if "merged" in dataset_id.lower() and "caller_af" not in dataset_id.lower():
        violations.append(
            "Dataset name suggests merged/synthetic source — "
            "verify 'caller_af' column exists (known pitfall: merged AF trap, see "
            "knowledge/05_data_formats/06_merged_dataset_pitfalls.md)"
        )

    if violations:
        return {"status": "schema_mismatch", "expected_id": dataset_id,
                "schema_violations": violations,
                "note": f"{len(violations)} violation(s) detected."}

    return {"status": "fresh", "expected_id": dataset_id,
            "schema_violations": [],
            "note": "No stale marks or known-pitfall name patterns."}


def check_upstream_reports(plan: dict) -> list[dict]:
    upstream = plan.get("preconditions", {}).get("upstream_reports", [])
    if not upstream:
        return []

    stale_marks = read_jsonl(INVALIDATION_DIR / "stale_marks.jsonl")
    stale_paths = {e["report_path"] for e in stale_marks if e.get("type") == "report"}

    results = []
    for report_path in upstream:
        # Normalize: strip InterSubMod/ prefix if present
        rel_path = report_path.removeprefix("InterSubMod/")
        abs_path = REPO_ROOT / rel_path

        if not abs_path.is_file():
            results.append({"report_path": report_path, "status": "missing",
                            "stale_reason": "file not found"})
            continue

        # Stale mark check
        if report_path in stale_paths or rel_path in stale_paths:
            mark = next((e for e in stale_marks if e.get("report_path") in (report_path, rel_path)), {})
            results.append({"report_path": report_path, "status": "stale_marked",
                            "stale_reason": mark.get("reason", "marked in stale_marks.jsonl")})
            continue

        # Retraction marker check (frontmatter or first ~20 lines)
        try:
            head_text = "".join(abs_path.read_text(errors="replace").splitlines(keepends=True)[:30])
            if re.search(r"\[RETRACTED", head_text) or re.search(r"^retraction:\s*true", head_text, re.M):
                results.append({"report_path": report_path, "status": "retracted",
                                "stale_reason": "retraction marker in report header"})
                continue
        except Exception as e:
            results.append({"report_path": report_path, "status": "skipped",
                            "stale_reason": f"unreadable: {e}"})
            continue

        results.append({"report_path": report_path, "status": "fresh", "stale_reason": None})

    return results


def read_jsonl(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    entries = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            entries.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return entries


def determine_verdict(binary: dict, dataset: dict, reports: list[dict]) -> tuple[str, list[str]]:
    blocking = []

    if binary["status"] in ("stale", "missing"):
        blocking.append(f"binary {binary['status']}: {binary['note']}")
    if dataset["status"] in ("stale", "schema_mismatch", "missing"):
        for v in dataset.get("schema_violations", []):
            blocking.append(f"dataset: {v}")
        if not dataset.get("schema_violations"):
            blocking.append(f"dataset {dataset['status']}: {dataset['note']}")
    for r in reports:
        if r["status"] in ("stale_marked", "retracted", "missing"):
            blocking.append(f"upstream {r['status']}: {r['report_path']} ({r['stale_reason']})")

    if blocking:
        return "BLOCKED", blocking

    # WARN conditions could go here (currently none in Path A)
    return "PASS", []


def write_precheck(cycle_id: str, payload: dict) -> Path:
    out_path = STATE_ROOT / "cycles" / cycle_id / "precheck.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    return out_path


def render_summary(cycle_id: str, payload: dict) -> str:
    lines = [f"[/check-staleness] cycle_id={cycle_id}"]
    b = payload["checks"]["binary"]
    icon = "✓" if b["status"] == "fresh" else ("○" if b["status"] == "skipped" else "✗")
    lines.append(f"  Binary:           {icon} {b['status']} — {b['note']}")
    d = payload["checks"]["dataset"]
    icon = "✓" if d["status"] == "fresh" else ("○" if d["status"] == "skipped" else "✗")
    lines.append(f"  Dataset:          {icon} {d['status']} — {d['note']}")
    rs = payload["checks"]["upstream_reports"]
    if rs:
        n_fresh = sum(1 for r in rs if r["status"] == "fresh")
        bad = [r for r in rs if r["status"] not in ("fresh", "skipped")]
        if bad:
            lines.append(f"  Upstream reports: ✗ {n_fresh}/{len(rs)} fresh; issues:")
            for r in bad:
                lines.append(f"    - {r['report_path']}: {r['status']} ({r['stale_reason']})")
        else:
            lines.append(f"  Upstream reports: ✓ {n_fresh}/{len(rs)} fresh")
    else:
        lines.append("  Upstream reports: ○ none declared")
    lines.append(f"verdict: {payload['verdict']}")
    if payload.get("blocking_reasons"):
        lines.append("blocking_reasons:")
        for r in payload["blocking_reasons"]:
            lines.append(f"  - {r}")
    return "\n".join(lines)


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: _staleness_check.py <cycle_id>\n")
        sys.exit(2)

    cycle_id = sys.argv[1]
    plan = load_plan(cycle_id)

    binary = check_binary(plan)
    dataset = check_dataset(plan)
    reports = check_upstream_reports(plan)

    verdict, blocking = determine_verdict(binary, dataset, reports)

    payload = {
        "schema_version": "1.0",
        "cycle_id": cycle_id,
        "checked_at": utcnow_iso(),
        "checks": {
            "binary": binary,
            "dataset": dataset,
            "upstream_reports": reports,
        },
        "verdict": verdict,
        "blocking_reasons": blocking,
    }

    out_path = write_precheck(cycle_id, payload)
    print(render_summary(cycle_id, payload))
    print(f"\nWritten to: {out_path.relative_to(REPO_ROOT)}")

    sys.exit(1 if verdict == "BLOCKED" else 0)


if __name__ == "__main__":
    main()
