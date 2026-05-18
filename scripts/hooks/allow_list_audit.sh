#!/bin/bash
# allow_list_audit.sh — Manual audit of .claude/settings.local.json permissions.allow list.
#
# Purpose: 158 entries 累積無 review。偵測 (1) 過度具體 (整段 Bash command) (2) 重複 (3) 寬鬆潛在風險。
# Source: P2 audit M3 — Permission audit (plan §F.4 P3 Living-update)
# Usage:  bash InterSubMod/scripts/hooks/allow_list_audit.sh
# Output: stdout report + optional --json mode

set -uo pipefail

SETTINGS="/big7_disk/liaoyoyo2001/InterSubMod/.claude/settings.local.json"
MODE="${1:-report}"  # report | json

if [ ! -f "$SETTINGS" ]; then
    echo "[allow_list_audit] settings.local.json not found at $SETTINGS" >&2
    exit 1
fi

python3 << PYEOF
import json
import re
from collections import Counter
from pathlib import Path

settings = json.loads(Path("$SETTINGS").read_text())
allow = settings.get("permissions", {}).get("allow", [])
deny = settings.get("permissions", {}).get("deny", [])

# Categorize by prefix
categories = Counter()
exact_dupes = Counter()
overly_specific = []   # full Bash commands with specific args (>100 chars)
broad_wildcards = []   # Bash(*) or similar
file_paths_specific = []  # Read(//abs/path/**) referencing specific dirs

for entry in allow:
    exact_dupes[entry] += 1
    # Extract tool prefix
    m = re.match(r"^([A-Z][a-zA-Z]+)\(", entry)
    if m:
        categories[m.group(1)] += 1
    else:
        categories["[other]"] += 1

    # Overly specific (full command with many args)
    if entry.startswith("Bash(") and len(entry) > 120:
        overly_specific.append(entry)

    # Wildcards
    if entry in ("Bash(*)", "Bash(:*)", "WebFetch", "WebSearch"):
        broad_wildcards.append(entry)

    # Specific Read paths
    if entry.startswith("Read(//") and "**" in entry:
        file_paths_specific.append(entry)

# Find duplicates
dupes = {k: v for k, v in exact_dupes.items() if v > 1}

# Patterns: very specific Bash commands (likely one-shot historical)
one_shot_pattern = re.compile(r"^Bash\([^)]{60,}\)$")
one_shot_count = sum(1 for e in allow if one_shot_pattern.match(e))

# Old timestamped paths likely stale
timestamp_pattern = re.compile(r"(20\d{6}|/tmp/[a-z_]+_test)")
timestamp_entries = [e for e in allow if timestamp_pattern.search(e)]

output = {
    "audit_date": "$(date -Iseconds)",
    "total_allow": len(allow),
    "total_deny": len(deny),
    "categories": dict(categories),
    "exact_duplicates": dict(dupes),
    "overly_specific_count": len(overly_specific),
    "overly_specific_top_5": overly_specific[:5],
    "broad_wildcards": broad_wildcards,
    "one_shot_likely_count": one_shot_count,
    "timestamp_specific_count": len(timestamp_entries),
    "timestamp_specific_top_5": timestamp_entries[:5],
}

mode = "$MODE"

if mode == "json":
    print(json.dumps(output, ensure_ascii=False, indent=2))
else:
    print(f"# allow_list_audit — {output['audit_date']}")
    print()
    print(f"**Total allow entries**: {output['total_allow']}")
    print(f"**Total deny entries**:  {output['total_deny']}")
    print()
    print("## Categories")
    for cat, cnt in sorted(output["categories"].items(), key=lambda x: -x[1]):
        print(f"  {cat:<15} {cnt}")
    print()
    print(f"## Exact duplicates: {len(output['exact_duplicates'])}")
    for k, v in output["exact_duplicates"].items():
        print(f"  ×{v}: {k[:80]}")
    print()
    print(f"## Overly specific (Bash with >120 chars, likely one-shot): {output['overly_specific_count']}")
    for e in output["overly_specific_top_5"]:
        print(f"  - {e[:100]}...")
    print()
    print(f"## Likely one-shot (>60 char Bash commands): {output['one_shot_likely_count']}")
    print(f"  → review & remove if no longer used")
    print()
    print(f"## Timestamp / /tmp_test specific (likely stale): {output['timestamp_specific_count']}")
    for e in output["timestamp_specific_top_5"]:
        print(f"  - {e[:100]}")
    print()
    print(f"## Broad wildcards: {len(output['broad_wildcards'])}")
    for e in output["broad_wildcards"]:
        print(f"  - {e}")
    print()
    print("## Recommendations")
    over_pct = output['one_shot_likely_count'] * 100 / max(1, output['total_allow'])
    print(f"- {over_pct:.0f}% of entries likely one-shot — review quarterly")
    if output['exact_duplicates']:
        print(f"- Remove {len(output['exact_duplicates'])} exact duplicates")
    if output['timestamp_specific_count'] > 5:
        print(f"- {output['timestamp_specific_count']} timestamped/tmp entries should be migrated to wildcards or removed")
PYEOF
