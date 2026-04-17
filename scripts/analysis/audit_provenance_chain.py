#!/usr/bin/env python3
"""Audit the evidence provenance chain for InterSubMod.

Cross-checks:
  1. hypothesis_queue.json ↔ evidence_ledger.jsonl ↔ cycles/*/
  2. tier distribution (uses score_evidence_tier logic)
  3. orphans: cycles without ledger, ledger without cycles
  4. staleness: tier ≤ 2 positive entries older than 365 days
  5. MEMORY.md Concluded consistency with ledger verdicts

Usage:
    python3 scripts/analysis/audit_provenance_chain.py                # full report
    python3 scripts/analysis/audit_provenance_chain.py --weekly       # weekly-report format
    python3 scripts/analysis/audit_provenance_chain.py --json         # machine-readable
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timedelta
from pathlib import Path

PROJECT_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
LEDGER = PROJECT_ROOT / "research/autoresearch/evidence_ledger.jsonl"
QUEUE = PROJECT_ROOT / "research/autoresearch/hypothesis_queue.json"
CYCLES_DIR = PROJECT_ROOT / "research/autoresearch/cycles"
MEMORY = Path("/bip7_disk/liaoyoyo2001/.claude/projects/"
              "-big7-disk-liaoyoyo2001-InterSubMod/memory/MEMORY.md")

FLAG_TIER_MIN = {
    "within_group_ols": 3,
    "af_bin_stratified": 3,
    "permutation_tested": 3,
    "nested_cv": 3,
    "multi_method_agreement": 4,
    "multi_sample_consistent": 4,
    "literature_aligned": 5,
    "independent_reproduction": 5,
}

POSITIVE_VERDICTS = {"positive_pilot", "positive_validated"}
CONCLUSIVE_VERDICTS = POSITIVE_VERDICTS | {"negative", "no_go", "conditional_no_go"}


def compute_cap(flags: list[str]) -> int:
    if not flags:
        return 1
    has_t3 = any(FLAG_TIER_MIN.get(f, 99) <= 3 for f in flags)
    has_t4 = any(FLAG_TIER_MIN.get(f, 99) == 4 for f in flags)
    has_t5 = any(FLAG_TIER_MIN.get(f, 99) == 5 for f in flags)
    if has_t5 and has_t4 and has_t3:
        return 5
    if has_t4 and has_t3:
        return 4
    if has_t3:
        return 3
    return 2


def load_ledger() -> list[dict]:
    entries = []
    if not LEDGER.exists():
        return entries
    with LEDGER.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return entries


def load_queue() -> list[dict]:
    if not QUEUE.exists():
        return []
    with QUEUE.open() as f:
        data = json.load(f)
    return data if isinstance(data, list) else data.get("hypotheses", [])


def list_cycles() -> set[str]:
    if not CYCLES_DIR.exists():
        return set()
    return {d.name for d in CYCLES_DIR.iterdir() if d.is_dir()}


def audit_provenance(entries: list[dict]) -> dict:
    """Check each ledger entry has artifacts_path that resolves."""
    missing = []
    present = []
    for e in entries:
        hid = e.get("hypothesis_id", "?")
        path = e.get("artifacts_path", "")
        if not path:
            missing.append({"hid": hid, "reason": "no artifacts_path"})
            continue
        full = PROJECT_ROOT / path
        if not full.exists():
            missing.append({"hid": hid, "reason": f"path missing: {path}"})
        else:
            present.append({"hid": hid, "path": path})
    return {"present": present, "missing": missing}


def audit_orphan_cycles(entries: list[dict]) -> list[str]:
    """Cycles without corresponding ledger entries."""
    ledger_cycle_ids = {e.get("cycle_id") for e in entries if e.get("cycle_id")}
    cycles_on_disk = list_cycles()
    return sorted(cycles_on_disk - ledger_cycle_ids)


def audit_tier_distribution(entries: list[dict]) -> dict:
    """Tier distribution + over-claimed detection."""
    dist = Counter()
    over_claimed = []
    missing_tier = []
    for e in entries:
        hid = e.get("hypothesis_id", "?")
        tier = e.get("tier")
        flags = e.get("tier_flags", [])
        cap = compute_cap(flags)
        if tier is None:
            missing_tier.append(hid)
            dist["missing"] += 1
        else:
            dist[tier] += 1
            if tier > cap:
                over_claimed.append({"hid": hid, "tier": tier, "cap": cap})
    return {"dist": dict(dist), "over_claimed": over_claimed, "missing_tier": missing_tier}


def audit_staleness(entries: list[dict], days: int = 365) -> list[dict]:
    """Positive entries with tier ≤ 2 older than `days`."""
    cutoff = datetime.now() - timedelta(days=days)
    stale = []
    for e in entries:
        ts = e.get("timestamp", "")
        if not ts:
            continue
        try:
            t = datetime.fromisoformat(ts.replace("Z", "+00:00").split("+")[0])
        except ValueError:
            continue
        verdict = e.get("verdict", "")
        tier = e.get("tier") or 1
        if verdict in POSITIVE_VERDICTS and tier <= 2 and t < cutoff:
            stale.append({
                "hid": e.get("hypothesis_id", "?"),
                "verdict": verdict,
                "tier": tier,
                "age_days": (datetime.now() - t).days,
            })
    return stale


def audit_memory_vs_ledger(entries: list[dict]) -> dict:
    """Cross-check MEMORY.md Concluded section against ledger verdicts."""
    if not MEMORY.exists():
        return {"error": f"MEMORY.md not found at {MEMORY}"}

    text = MEMORY.read_text(encoding="utf-8", errors="ignore")
    concluded_markers = re.findall(r"POSITIVE|NEGATIVE|NO-GO|CONDITIONAL NO-GO", text)
    concluded_count = Counter(concluded_markers)

    ledger_verdict_count = Counter(e.get("verdict", "unknown") for e in entries)
    return {
        "memory_markers": dict(concluded_count),
        "ledger_verdicts": dict(ledger_verdict_count),
    }


def format_full_report(audit: dict) -> str:
    lines = ["# Provenance Chain Audit Report", ""]
    lines.append(f"**Generated**: {datetime.now().isoformat()}")
    lines.append(f"**Ledger**: {LEDGER}")
    lines.append("")

    lines.append("## 1. Provenance Coverage")
    prov = audit["provenance"]
    lines.append(f"- Entries with valid artifacts_path: {len(prov['present'])}")
    lines.append(f"- Entries with missing artifacts: {len(prov['missing'])}")
    if prov["missing"]:
        lines.append("")
        lines.append("**Missing artifacts:**")
        for m in prov["missing"]:
            lines.append(f"  - {m['hid']}: {m['reason']}")
    lines.append("")

    lines.append("## 2. Orphan Cycles (no ledger entry)")
    orphans = audit["orphans"]
    if orphans:
        lines.append(f"Found {len(orphans)} orphan cycles:")
        for c in orphans:
            lines.append(f"  - {c}")
    else:
        lines.append("None (all cycles are tracked in ledger)")
    lines.append("")

    lines.append("## 3. Tier Distribution")
    tier = audit["tier"]
    for level in [1, 2, 3, 4, 5, "missing"]:
        count = tier["dist"].get(level, 0)
        if count > 0:
            lines.append(f"  - Tier {level}: {count}")
    if tier["over_claimed"]:
        lines.append("")
        lines.append("**Over-claimed (tier > confidence_cap):**")
        for o in tier["over_claimed"]:
            lines.append(f"  - {o['hid']}: declared tier={o['tier']} but cap={o['cap']}")
    if tier["missing_tier"]:
        lines.append(f"\n**Missing tier field** ({len(tier['missing_tier'])} entries): "
                     + ", ".join(tier["missing_tier"][:10])
                     + ("..." if len(tier["missing_tier"]) > 10 else ""))
    lines.append("")

    lines.append("## 4. Staleness (tier ≤ 2 positive, > 365 days)")
    stale = audit["staleness"]
    if stale:
        lines.append(f"Found {len(stale)} stale entries needing tier upgrade:")
        for s in stale:
            lines.append(f"  - {s['hid']} (verdict={s['verdict']}, tier={s['tier']}, age={s['age_days']}d)")
    else:
        lines.append("None")
    lines.append("")

    lines.append("## 5. MEMORY.md vs Ledger Cross-Check")
    mem = audit["memory"]
    if "error" in mem:
        lines.append(f"  {mem['error']}")
    else:
        lines.append(f"  MEMORY.md markers: {mem['memory_markers']}")
        lines.append(f"  Ledger verdicts:   {mem['ledger_verdicts']}")
    lines.append("")

    return "\n".join(lines)


def format_weekly(audit: dict) -> str:
    lines = ["## 研究證據盤點（本週）", ""]
    tier = audit["tier"]
    total = sum(tier["dist"].values())
    lines.append(f"- **總條目**: {total}")
    for level in [5, 4, 3, 2, 1]:
        count = tier["dist"].get(level, 0)
        if count > 0:
            pct = 100 * count / total if total else 0
            lines.append(f"- Tier {level}: {count} ({pct:.0f}%)")

    issues = len(tier["over_claimed"]) + len(audit["provenance"]["missing"]) + len(audit["staleness"])
    if issues > 0:
        lines.append("")
        lines.append(f"**需處理 issues: {issues}**")
        if tier["over_claimed"]:
            lines.append(f"  - {len(tier['over_claimed'])} 個 tier over-claimed")
        if audit["provenance"]["missing"]:
            lines.append(f"  - {len(audit['provenance']['missing'])} 個 artifacts 失效")
        if audit["staleness"]:
            lines.append(f"  - {len(audit['staleness'])} 個 stale entries（需升級）")
    else:
        lines.append("\n所有 checks 通過。")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weekly", action="store_true", help="週報格式")
    parser.add_argument("--json", action="store_true", help="JSON 格式")
    parser.add_argument("--stale-days", type=int, default=365)
    args = parser.parse_args()

    entries = load_ledger()
    if not entries:
        print(f"[error] no entries in {LEDGER}", file=sys.stderr)
        return 2

    audit = {
        "provenance": audit_provenance(entries),
        "orphans": audit_orphan_cycles(entries),
        "tier": audit_tier_distribution(entries),
        "staleness": audit_staleness(entries, days=args.stale_days),
        "memory": audit_memory_vs_ledger(entries),
        "meta": {"ledger_entries": len(entries), "cycles_on_disk": len(list_cycles())},
    }

    if args.json:
        print(json.dumps(audit, indent=2, default=str))
    elif args.weekly:
        print(format_weekly(audit))
    else:
        print(format_full_report(audit))

    return 0


if __name__ == "__main__":
    sys.exit(main())
