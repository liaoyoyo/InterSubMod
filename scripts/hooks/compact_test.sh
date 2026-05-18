#!/bin/bash
# compact_test.sh — Regression test: /compact 後 CLAUDE.md §7 保留指令是否仍在 transcript.
#
# Purpose: 防 Context Drift (P2 audit C2 challenge)
# Source: P3 audit Compact-Test
# Trigger: Manual / cron / pre-compact hook (optional)
# Usage:
#   bash InterSubMod/scripts/hooks/compact_test.sh                # latest transcript
#   bash InterSubMod/scripts/hooks/compact_test.sh <path.jsonl>   # specific
# Exit: 0 = all preservation directives still present, 1 = drift detected

set -o pipefail

TRANSCRIPT_DIR="/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod"
TARGET="${1:-latest}"

# CLAUDE.md §7 preservation directives (must remain visible after /compact)
# These are the directives Claude should ALWAYS know post-compact
PRESERVE_KEYS=(
    "Hard Gate"          # CLAUDE.md §1
    "暫停判斷矩陣"        # §1
    "Step→Verify"         # AGENTS.md §6
    "Cynefin"             # /confirmation-protocol (2026-05-18 落地)
    "Productive Failure"  # /scientific-rigor §8.3.1
    "Evidence Tier"       # /scientific-rigor §2
    "Pre-registration"    # §7.1
    "scientific-rigor"    # hub skill
)

export CT_TRANSCRIPT_DIR="$TRANSCRIPT_DIR"
export CT_TARGET="$TARGET"
# Use newline separator to preserve multi-word keys
export CT_PRESERVE_KEYS=$(printf '%s\n' "${PRESERVE_KEYS[@]}")

python3 << 'PYEOF'
import json
import os
import sys
from glob import glob

td = os.environ["CT_TRANSCRIPT_DIR"]
target = os.environ["CT_TARGET"]
keys = [k for k in os.environ["CT_PRESERVE_KEYS"].split("\n") if k.strip()]

if target == "latest":
    files = sorted(glob(f"{td}/*.jsonl"), key=os.path.getmtime, reverse=True)
    tpath = files[0] if files else None
else:
    tpath = target

if not tpath or not os.path.exists(tpath):
    print(f"[compact_test] transcript not found: {tpath}", file=sys.stderr)
    sys.exit(2)

# Find compact markers (system-reminder lines containing "compact" or "壓縮")
compact_indices = []
content_keys_found = {k: [] for k in keys}

with open(tpath, encoding="utf-8") as f:
    for i, line in enumerate(f):
        try:
            d = json.loads(line)
            content = json.dumps(d.get("message", {}).get("content", ""), ensure_ascii=False)
            # Only count REAL /compact user commands (not casual mentions)
            # Pattern: <command-name>/compact</command-name> or "/compact" as standalone slash command
            if "<command-name>/compact</command-name>" in content or "<command-name>compact</command-name>" in content:
                compact_indices.append(i)
            # Track which keys appear at which line
            for k in keys:
                if k in content:
                    content_keys_found[k].append(i)
        except Exception:
            pass

print(f"# compact_test — Context Drift regression check")
print(f"Transcript:    {os.path.basename(tpath)}")
print(f"Compact markers found at lines: {compact_indices if compact_indices else '— (no compaction yet)'}")
print()

# If compaction happened, check if preservation directives appear AFTER the last compact
if compact_indices:
    last_compact = compact_indices[-1]
    print(f"Last compact at line {last_compact}. Checking directives presence AFTER:")
    print()
    print("| Directive | Pre-compact | Post-compact | Status |")
    print("|-----------|-------------|--------------|--------|")
    drift_count = 0
    for k in keys:
        pre = sum(1 for x in content_keys_found[k] if x < last_compact)
        post = sum(1 for x in content_keys_found[k] if x >= last_compact)
        status = "✅" if post > 0 else "⚠ DRIFT"
        if post == 0 and pre > 0:
            drift_count += 1
        print(f"| `{k}` | {pre} | {post} | {status} |")
    print()
    if drift_count:
        print(f"⚠ Context drift detected: {drift_count}/{len(keys)} directives missing post-compact")
        print("→ Recommendation: trigger /research-context-loader Tier 1 to reload")
        sys.exit(1)
    else:
        print(f"✅ All {len(keys)} directives preserved post-compact")
        sys.exit(0)
else:
    print("## No compact event detected — full directive presence:")
    print()
    print("| Directive | Total mentions |")
    print("|-----------|---------------:|")
    for k in keys:
        cnt = len(content_keys_found[k])
        status_icon = "✅" if cnt > 0 else "❌"
        print(f"| `{k}` | {cnt} {status_icon} |")
    missing = [k for k in keys if not content_keys_found[k]]
    if missing:
        print(f"\n⚠ {len(missing)} directives never appeared in this session: {missing}")
        print("→ Check CLAUDE.md §7 cross-reference is correctly loaded at session start")
        sys.exit(1)
    sys.exit(0)
PYEOF
