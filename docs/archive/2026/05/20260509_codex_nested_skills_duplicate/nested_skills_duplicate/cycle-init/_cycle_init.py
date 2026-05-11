#!/usr/bin/env python3
"""
P0 REGISTER gate of the Resilient Waterfall harness.

Allocates a cycle_id, creates state/cycles/{cycle_id}/state.json skeleton,
and updates state/active.json index.

Usage:
  python3 _cycle_init.py <slug_or_hypothesis_id> --title "<title>" \\
       [--priority N] [--binary-version SHA] [--dataset-id ID] \\
       [--upstream-report PATH ...]

Schema refs:
  state/schemas/state.schema.json
  state/schemas/active.schema.json

Returns:
  exit 0  — success (cycle_id printed to stdout)
  exit 1  — invalid input
  exit 2  — write failure (e.g. corrupt active.json)
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
from pathlib import Path

# REPO_ROOT = parents[3] because file lives at:
#   InterSubMod/.agents/skills/cycle-init/_cycle_init.py
# parents[0] = cycle-init, [1] = skills, [2] = .agents, [3] = InterSubMod
REPO_ROOT = Path(__file__).resolve().parents[3]
STATE_ROOT = REPO_ROOT / "state"
ACTIVE_JSON = STATE_ROOT / "active.json"
CYCLES_DIR = STATE_ROOT / "cycles"
HYPOTHESIS_QUEUE = REPO_ROOT / "research" / "autoresearch" / "hypothesis_queue.json"
PITFALLS_TABLE = REPO_ROOT / ".agents" / "skills" / "run-evaluator" / "pitfalls_table.json"

CYCLE_ID_REGEX = re.compile(r"^[0-9]{8}-[0-9]{4}-[a-z0-9_-]+$")


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def slugify(text: str, max_len: int = 30) -> str:
    """Lowercase + replace non-alphanumeric with dash + collapse + truncate."""
    s = text.lower()
    s = re.sub(r"[^a-z0-9_-]+", "-", s)
    s = re.sub(r"-+", "-", s).strip("-")
    return s[:max_len].rstrip("-")


def allocate_cycle_id(slug: str) -> str:
    """Generate YYYYMMDD-HHMM-{slug}; auto-suffix on collision."""
    now = dt.datetime.now(dt.timezone.utc)
    base = f"{now.strftime('%Y%m%d-%H%M')}-{slug}"
    if not CYCLE_ID_REGEX.match(base):
        raise ValueError(f"Generated cycle_id violates schema regex: {base}")

    candidate = base
    suffix = 2
    while (CYCLES_DIR / candidate).exists():
        candidate = f"{base}-{suffix}"
        suffix += 1
        if suffix > 99:
            raise RuntimeError("Too many cycle_id collisions in one minute")
    return candidate


def load_active() -> dict:
    if not ACTIVE_JSON.exists():
        # Initialize fresh active.json
        return {
            "schema_version": "1.0",
            "updated_at": utc_now_iso(),
            "max_concurrent": 5,
            "cycles": [],
        }
    try:
        data = json.loads(ACTIVE_JSON.read_text())
    except json.JSONDecodeError as e:
        print(f"ERROR: active.json corrupt — {e}", file=sys.stderr)
        sys.exit(2)
    required = {"schema_version", "updated_at", "cycles"}
    if not required.issubset(data.keys()):
        print(f"ERROR: active.json missing required fields {required - set(data.keys())}", file=sys.stderr)
        sys.exit(2)
    return data


def hypothesis_id_in_queue(hypothesis_id: str) -> bool:
    """Soft check: does hypothesis_id exist in hypothesis_queue.json?"""
    if not HYPOTHESIS_QUEUE.exists():
        return False
    try:
        text = HYPOTHESIS_QUEUE.read_text()
        return hypothesis_id in text
    except Exception:
        return False


def dataset_id_pitfall_warn(dataset_id: str | None) -> list[str]:
    """Warn if dataset_id contains known-pitfall keywords."""
    warnings = []
    if not dataset_id:
        return warnings
    lower = dataset_id.lower()
    risky = ["merged", "pileup_symlink", "phase1_new"]
    hit = [k for k in risky if k in lower]
    if hit:
        warnings.append(f"dataset_id contains pitfall keywords: {hit} — P2 PRECHECK will check carefully")
    return warnings


def build_state(
    cycle_id: str,
    title: str,
    hypothesis_id: str | None,
    priority: int,
    binary_version: str | None,
    dataset_id: str | None,
    upstream_reports: list[str],
    main_axis: dict | None = None,
) -> dict:
    now = utc_now_iso()
    state = {
        "schema_version": "1.0",
        "cycle_id": cycle_id,
        "title": title,
        "hypothesis_id": hypothesis_id,
        "phase": "P0_REGISTER",
        "tier": "pending",
        "verdict": "active",
        "priority": priority,
        "started_at": now,
        "last_advanced_at": now,
        "preconditions": {
            "binary_version": binary_version,
            "dataset_id": dataset_id,
            "upstream_reports": upstream_reports,
        },
        "open_gates": [],
        "artifacts": {
            "plan": None,
            "precheck": None,
            "pilot": None,
            "generalize": None,
            "evaluation": None,
        },
        "history": [
            {
                "timestamp": now,
                "from_phase": None,
                "to_phase": "P0_REGISTER",
                "actor": "user",
                "note": "/cycle-init",
            }
        ],
        "interaction_metrics": {
            "user_interventions": 0,
            "auto_advancements": 0,
            "user_corrections": 0,
            "auto_recovery_attempts": 0,
            "auto_recovery_successes": 0,
            "drift_warnings_emitted": 0,
            "drift_warnings_acknowledged": 0,
        },
    }
    if main_axis is not None:
        state["main_axis"] = main_axis
    return state


def update_active(active: dict, cycle_id: str, title: str, hypothesis_id: str | None) -> tuple[dict, list[str]]:
    """Append cycle to active.cycles[]; return updated active + warnings."""
    warnings = []
    now = utc_now_iso()
    entry = {
        "cycle_id": cycle_id,
        "title": title,
        "hypothesis_id": hypothesis_id,
        "phase": "P0_REGISTER",
        "tier": "pending",
        "verdict": "active",
        "started_at": now,
        "last_advanced_at": now,
    }
    active["cycles"].append(entry)
    active["updated_at"] = now
    if len(active["cycles"]) > active.get("max_concurrent", 5):
        warnings.append(
            f"WARNING: {len(active['cycles'])} active cycles > max_concurrent={active.get('max_concurrent', 5)}; "
            "consider archiving completed cycles via /cycle-state."
        )
    return active, warnings


def write_state_and_active(cycle_id: str, state: dict, active: dict) -> None:
    cycle_dir = CYCLES_DIR / cycle_id
    cycle_dir.mkdir(parents=True, exist_ok=True)
    (cycle_dir / "state.json").write_text(json.dumps(state, indent=2) + "\n")
    ACTIVE_JSON.write_text(json.dumps(active, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Initialize a new P0 REGISTER cycle")
    ap.add_argument("slug_or_hypothesis_id", help="hypothesis_id (e.g. H-2026-05-07-001) or slug (e.g. loh-kde-rescale)")
    ap.add_argument("--title", required=True, help="Cycle title (≤80 chars human-readable)")
    ap.add_argument("--priority", type=int, default=50, help="0-100 priority (default 50)")
    ap.add_argument("--binary-version", default=None, help="C++ binary git SHA (or null for pure-analysis)")
    ap.add_argument("--dataset-id", default=None, help="Dataset identifier")
    ap.add_argument("--upstream-report", action="append", default=[], help="Repeatable; upstream report paths")
    ap.add_argument("--hypothesis-id", default=None, help="Override hypothesis_id (default = positional arg if it looks like H-...)")
    ap.add_argument("--main-axis-anchor", default=None, help="M1: InterSubMod/-prefixed path to anchor doc (e.g. CURRENT_FOCUS.md)")
    ap.add_argument("--main-axis-section", default=None, help="M1: section heading within anchor doc")
    ap.add_argument("--main-axis-goal", default=None, help="M1: one-line goal restatement (required if --main-axis-anchor set)")
    ap.add_argument("--main-axis-required", action="append", default=[], help="M1: repeatable; required keywords")
    ap.add_argument("--main-axis-forbidden", action="append", default=[], help="M1: repeatable; forbidden keywords")
    args = ap.parse_args()

    # Validate priority
    if not 0 <= args.priority <= 100:
        print(f"ERROR: priority must be 0-100 (got {args.priority})", file=sys.stderr)
        return 1

    # Determine hypothesis_id vs slug
    pos = args.slug_or_hypothesis_id
    if args.hypothesis_id:
        hypothesis_id = args.hypothesis_id
        slug_source = pos
    elif re.match(r"^H-\d{4}-\d{2}-\d{2}-\d{3}$", pos):
        hypothesis_id = pos
        slug_source = args.title
    else:
        hypothesis_id = None
        slug_source = pos

    # Build slug
    slug = slugify(slug_source) or slugify(args.title) or f"unnamed-{int(dt.datetime.now().timestamp())}"

    # Allocate cycle_id
    try:
        cycle_id = allocate_cycle_id(slug)
    except (ValueError, RuntimeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    # Soft-check hypothesis_id
    warnings: list[str] = []
    if hypothesis_id and not hypothesis_id_in_queue(hypothesis_id):
        warnings.append(f"WARNING: hypothesis_id={hypothesis_id} not found in hypothesis_queue.json (continuing)")

    # Pitfall warnings on dataset_id
    warnings.extend(dataset_id_pitfall_warn(args.dataset_id))

    # Build main_axis dict (M1 主軸鎖) if user provided
    main_axis = None
    if args.main_axis_anchor:
        if not args.main_axis_goal:
            print("ERROR: --main-axis-goal is required when --main-axis-anchor is set", file=sys.stderr)
            return 1
        main_axis = {
            "anchor_doc": args.main_axis_anchor,
            "anchor_section": args.main_axis_section,
            "one_line_goal": args.main_axis_goal,
            "drift_keywords_required": args.main_axis_required,
            "drift_keywords_forbidden": args.main_axis_forbidden,
        }

    # Build state.json
    state = build_state(
        cycle_id=cycle_id,
        title=args.title,
        hypothesis_id=hypothesis_id,
        priority=args.priority,
        binary_version=args.binary_version,
        dataset_id=args.dataset_id,
        upstream_reports=args.upstream_report,
        main_axis=main_axis,
    )

    # Update active.json
    active = load_active()
    active, more_warnings = update_active(active, cycle_id, args.title, hypothesis_id)
    warnings.extend(more_warnings)

    # Write
    try:
        write_state_and_active(cycle_id, state, active)
    except OSError as e:
        print(f"ERROR: write failed — {e}", file=sys.stderr)
        return 2

    # Output summary
    print(f"[/cycle-init] cycle_id={cycle_id}")
    print(f"  hypothesis_id: {hypothesis_id or '(none)'}")
    print(f'  title: "{args.title}"')
    print(f"  phase: P0_REGISTER")
    print(f"  priority: {args.priority}")
    if main_axis:
        print(f"  main_axis: {main_axis['anchor_doc']} → {main_axis['one_line_goal']}")
        if main_axis["drift_keywords_required"]:
            print(f"    required keywords: {main_axis['drift_keywords_required']}")
    else:
        print(f"  main_axis: (not set; M1 drift detection disabled for this cycle)")
    print(f"  state.json: state/cycles/{cycle_id}/state.json")
    print(f"  active.json updated ({len(active['cycles'])} active cycles total)")
    for w in warnings:
        print(f"  {w}")
    print()
    print("Next: invoke /research-loop or research-loop skill to draft plan.json (P1 PLAN)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
