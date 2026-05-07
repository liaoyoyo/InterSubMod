#!/usr/bin/env python3
"""
Master dashboard for the Resilient Waterfall harness.

Read-only. Loads state/active.json + each state/cycles/{id}/state.json,
sorts by priority + last_advanced, emits dashboard with routing recommendations.

Usage:
  python3 _cycle_state.py [--filter CYCLE_ID] [--phase PHASE] \\
       [--min-priority N] [--format markdown|json|plain] [--include-retro]

Exit codes:
  0  — success (dashboard printed)
  1  — invalid CLI args
  2  — active.json corrupt (cannot proceed)
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
STATE_ROOT = REPO_ROOT / "state"
ACTIVE_JSON = STATE_ROOT / "active.json"
CYCLES_DIR = STATE_ROOT / "cycles"
RETRO_CYCLES_DIR = STATE_ROOT / "retro_cycles"

VALID_FORMATS = {"markdown", "json", "plain"}
PHASE_ENUM = ["P0_REGISTER", "P1_PLAN", "P2_PRECHECK", "P3_PILOT",
              "P4_GENERALIZE", "P5_EVALUATE", "P6_COMMIT"]
STALE_DAYS = 7


def utc_now() -> dt.datetime:
    return dt.datetime.now(dt.timezone.utc)


def parse_iso(s: str) -> dt.datetime | None:
    """Parse ISO-8601 (UTC). Return None on failure."""
    if not s:
        return None
    try:
        if s.endswith("Z"):
            s = s[:-1] + "+00:00"
        return dt.datetime.fromisoformat(s)
    except (ValueError, TypeError):
        return None


def humanize_age(start: dt.datetime | None, now: dt.datetime) -> str:
    if start is None:
        return "?"
    delta = now - start
    if delta.total_seconds() < 0:
        return "future?"
    secs = delta.total_seconds()
    if secs < 3600:
        return f"{int(secs / 60)}m ago"
    if secs < 86400:
        return f"{int(secs / 3600)}h ago"
    return f"{delta.days}d ago"


def load_active() -> dict:
    if not ACTIVE_JSON.exists():
        return {"schema_version": "1.0", "updated_at": "", "max_concurrent": 5, "cycles": []}
    try:
        data = json.loads(ACTIVE_JSON.read_text())
    except json.JSONDecodeError as e:
        print(f"ERROR: active.json corrupt — {e}", file=sys.stderr)
        sys.exit(2)
    return data


def load_state(cycle_id: str, dir_root: Path) -> dict | None:
    state_path = dir_root / cycle_id / "state.json"
    if not state_path.exists():
        return None
    try:
        return json.loads(state_path.read_text())
    except json.JSONDecodeError:
        return {"_corrupt": True, "cycle_id": cycle_id}


def collect_cycles(active: dict, include_retro: bool) -> list[dict]:
    """Build list of dicts merging active.json entry + state.json content."""
    cycles = []
    for entry in active.get("cycles", []):
        cycle_id = entry["cycle_id"]
        state = load_state(cycle_id, CYCLES_DIR)
        if state is None:
            cycles.append({**entry, "_orphan": True})
        else:
            merged = {**entry, **state}
            cycles.append(merged)
    if include_retro and RETRO_CYCLES_DIR.is_dir():
        for retro_dir in sorted(RETRO_CYCLES_DIR.iterdir()):
            if not retro_dir.is_dir():
                continue
            state = load_state(retro_dir.name, RETRO_CYCLES_DIR)
            if state and not state.get("_corrupt"):
                state["_retro"] = True
                cycles.append(state)
    return cycles


def filter_cycles(cycles: list[dict], cycle_id: str | None,
                  phase: str | None, min_priority: int) -> list[dict]:
    out = cycles
    if cycle_id:
        out = [c for c in out if c.get("cycle_id") == cycle_id]
    if phase:
        out = [c for c in out if c.get("phase") == phase]
    if min_priority > 0:
        out = [c for c in out if c.get("priority", 50) >= min_priority]
    return out


def sort_key(cycle: dict) -> tuple:
    priority = cycle.get("priority", 50)
    last_adv = parse_iso(cycle.get("last_advanced_at", "")) or dt.datetime.min.replace(tzinfo=dt.timezone.utc)
    return (-priority, -last_adv.timestamp())


def next_action_for(cycle: dict) -> str:
    """Suggest next skill to invoke based on phase + verdict + artifacts."""
    if cycle.get("_orphan"):
        return "[ORPHAN] remove from active.json or recreate state.json"
    if cycle.get("_corrupt"):
        return "[CORRUPT] manually fix state.json"
    verdict = cycle.get("verdict", "active")
    phase = cycle.get("phase", "P0_REGISTER")
    artifacts = cycle.get("artifacts", {})

    if verdict == "blocked":
        return f"review blocked_reason; consider pivot-direction"
    if verdict == "exit_negative":
        return "/conclude-research (NEGATIVE wrap-up)"
    if verdict == "completed":
        return "archive cycle (move to cycles_archived/)"

    actions = {
        "P0_REGISTER": "/research-loop or research-loop skill (write plan.json)",
        "P1_PLAN": "/check-staleness" if artifacts.get("plan") else "/research-loop (plan.json missing)",
        "P2_PRECHECK": _p2_action(cycle),
        "P3_PILOT": "/multi-sample-consistency" if artifacts.get("pilot") else "/feature-layered-observation",
        "P4_GENERALIZE": "/run-evaluator" if artifacts.get("generalize") else "/multi-sample-consistency",
        "P5_EVALUATE": "/conclude-research" if artifacts.get("evaluation") else "/run-evaluator",
        "P6_COMMIT": "archive cycle (final step)",
    }
    return actions.get(phase, "(unknown phase)")


def _p2_action(cycle: dict) -> str:
    """Refine P2 action based on precheck verdict if available."""
    artifacts = cycle.get("artifacts", {})
    if not artifacts.get("precheck"):
        return "/check-staleness (precheck.json missing)"
    # Try reading precheck.json for verdict
    cycle_id = cycle.get("cycle_id", "")
    precheck_path = CYCLES_DIR / cycle_id / "precheck.json"
    if precheck_path.exists():
        try:
            verdict = json.loads(precheck_path.read_text()).get("verdict", "?")
            if verdict == "PASS":
                return "/test-quick or /feature-layered-observation"
            if verdict == "BLOCKED":
                return "review precheck.json reason; rebuild C++ binary if stale"
            if verdict == "WARN":
                return "review precheck.json; proceed with caution to P3"
        except json.JSONDecodeError:
            pass
    return "/check-staleness (re-run; precheck unclear)"


def open_gates_summary(cycle: dict) -> str:
    gates = cycle.get("open_gates", [])
    if not gates:
        return "–"
    return ", ".join(f"{g.get('gate')}={g.get('status')}" for g in gates)


def health_checks(active: dict, cycles: list[dict], now: dt.datetime) -> list[str]:
    checks = []
    n = len([c for c in cycles if not c.get("_retro")])
    cap = active.get("max_concurrent", 5)
    if n <= cap:
        checks.append(f"✅ active count ({n}) within max_concurrent ({cap})")
    else:
        checks.append(f"⚠ active count ({n}) > max_concurrent ({cap}) — consider archiving")

    orphan_count = sum(1 for c in cycles if c.get("_orphan"))
    if orphan_count == 0:
        checks.append("✅ no orphan references")
    else:
        checks.append(f"⚠ {orphan_count} orphan reference(s) — clean active.json")

    stale = []
    for c in cycles:
        if c.get("_retro") or c.get("_orphan") or c.get("_corrupt"):
            continue
        last = parse_iso(c.get("last_advanced_at", ""))
        if last and (now - last).days > STALE_DAYS:
            stale.append(c["cycle_id"])
    if stale:
        checks.append(f"⚠ {len(stale)} cycle(s) stale >{STALE_DAYS}d: {', '.join(stale)}")
    elif n > 0:
        checks.append(f"✅ no stale cycles (>{STALE_DAYS}d threshold)")

    return checks


def routing_recommendations(cycles: list[dict], now: dt.datetime) -> list[str]:
    recs = []
    for c in cycles:
        if c.get("_retro"):
            continue
        if c.get("_orphan"):
            recs.append(f"⚠ {c.get('cycle_id', '?')}: orphan — remove from active.json or recreate state.json")
            continue
        if c.get("_corrupt"):
            recs.append(f"⚠ {c.get('cycle_id', '?')}: state.json corrupt — manual fix required")
            continue
        verdict = c.get("verdict", "active")
        if verdict == "blocked":
            reason = c.get("blocked_reason") or "(no reason)"
            recs.append(f"⚠ {c['cycle_id']}: blocked — {reason}")
        gates = c.get("open_gates", [])
        for g in gates:
            if g.get("status") in ("open", "blocked"):
                checked = parse_iso(g.get("checked_at", ""))
                age_str = humanize_age(checked, now) if checked else "?"
                recs.append(f"⚠ {c['cycle_id']}: open gate {g.get('gate')}={g.get('status')} ({age_str})")
    if not recs:
        recs.append("✅ no urgent routing issues")
    return recs


def render_markdown(active: dict, cycles: list[dict], now: dt.datetime) -> str:
    lines = []
    lines.append("# /cycle-state — Master Dashboard")
    lines.append(f"Generated at: {now.strftime('%Y-%m-%dT%H:%M:%SZ')}")
    n_active = len([c for c in cycles if not c.get("_retro")])
    cap = active.get("max_concurrent", 5)
    lines.append(f"Active cycles: {n_active} / max_concurrent={cap}")
    n_retro = sum(1 for c in cycles if c.get("_retro"))
    if n_retro:
        lines.append(f"Retro cycles included: {n_retro}")
    lines.append("")

    if n_active == 0 and n_retro == 0:
        lines.append("**No active cycles.** Use `/cycle-init` to start one.")
        return "\n".join(lines)

    lines.append("## Cycles (sorted by priority desc, then last_advanced_at)")
    lines.append("")
    lines.append("| cycle_id | title | phase | tier | pri | last_advanced | open_gates | next_action |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for c in cycles:
        cid = c.get("cycle_id", "?")
        title = (c.get("title") or "").replace("|", "\\|")[:40]
        phase = c.get("phase", "?")
        tier = c.get("tier", "?")
        pri = c.get("priority", "-")
        last = parse_iso(c.get("last_advanced_at", ""))
        last_str = humanize_age(last, now) if last else "?"
        gates = open_gates_summary(c)
        action = next_action_for(c)
        retro_mark = " (retro)" if c.get("_retro") else ""
        lines.append(f"| {cid}{retro_mark} | {title} | {phase} | {tier} | {pri} | {last_str} | {gates} | {action} |")
    lines.append("")

    lines.append("## Routing recommendations")
    for rec in routing_recommendations(cycles, now):
        lines.append(f"- {rec}")
    lines.append("")

    lines.append("## Health checks")
    for check in health_checks(active, cycles, now):
        lines.append(f"- {check}")
    return "\n".join(lines)


def render_plain(active: dict, cycles: list[dict], now: dt.datetime) -> str:
    lines = [f"cycle-state\t{now.strftime('%Y-%m-%dT%H:%M:%SZ')}\tactive={len(cycles)}\tcap={active.get('max_concurrent', 5)}"]
    for c in cycles:
        last = parse_iso(c.get("last_advanced_at", ""))
        last_str = humanize_age(last, now) if last else "?"
        lines.append("\t".join([
            c.get("cycle_id", "?"),
            c.get("phase", "?"),
            str(c.get("tier", "?")),
            str(c.get("priority", "-")),
            last_str,
            next_action_for(c),
        ]))
    return "\n".join(lines)


def render_json(active: dict, cycles: list[dict], now: dt.datetime) -> str:
    out = {
        "generated_at": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "active_count": len([c for c in cycles if not c.get("_retro")]),
        "max_concurrent": active.get("max_concurrent", 5),
        "cycles": [
            {
                "cycle_id": c.get("cycle_id"),
                "title": c.get("title"),
                "phase": c.get("phase"),
                "tier": c.get("tier"),
                "priority": c.get("priority"),
                "last_advanced_at": c.get("last_advanced_at"),
                "open_gates": c.get("open_gates", []),
                "next_action": next_action_for(c),
                "is_retro": bool(c.get("_retro")),
                "is_orphan": bool(c.get("_orphan")),
                "is_corrupt": bool(c.get("_corrupt")),
            }
            for c in cycles
        ],
        "routing_recommendations": routing_recommendations(cycles, now),
        "health_checks": health_checks(active, cycles, now),
    }
    return json.dumps(out, indent=2, ensure_ascii=False)


def main() -> int:
    ap = argparse.ArgumentParser(description="Master dashboard for Resilient Waterfall harness")
    ap.add_argument("--filter", help="Only show this cycle_id")
    ap.add_argument("--phase", choices=PHASE_ENUM, help="Filter by phase")
    ap.add_argument("--min-priority", type=int, default=0, help="Minimum priority to include")
    ap.add_argument("--format", choices=sorted(VALID_FORMATS), default="markdown")
    ap.add_argument("--include-retro", action="store_true", help="Include state/retro_cycles/")
    args = ap.parse_args()

    if not 0 <= args.min_priority <= 100:
        print(f"ERROR: --min-priority must be 0-100 (got {args.min_priority})", file=sys.stderr)
        return 1

    active = load_active()
    cycles = collect_cycles(active, args.include_retro)
    cycles = filter_cycles(cycles, args.filter, args.phase, args.min_priority)
    cycles.sort(key=sort_key)
    now = utc_now()

    if args.format == "markdown":
        print(render_markdown(active, cycles, now))
    elif args.format == "json":
        print(render_json(active, cycles, now))
    else:  # plain
        print(render_plain(active, cycles, now))
    return 0


if __name__ == "__main__":
    sys.exit(main())
