#!/bin/bash
# SessionStart hook — inject CURRENT_FOCUS.md latest section into context.
#
# Source:   InterSubMod/docs/CURRENT_FOCUS.md
# Target:   ~500 tokens (~3000 chars) of the most recent dated section
# Triggers: SessionStart (every new conversation)
# Side effect: warn if CURRENT_FOCUS.md mtime > 7 days (CLAUDE.md §6 contract)
#
# CLAUDE.md §6 「Working State Pointer」 + §4「規劃新增」 — 本 hook 落地

set -uo pipefail

python3 << 'PYEOF'
import json
import os
import re
import time

FOCUS_FILE = "/big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md"
MAX_CHARS = 3000

def emit(ctx: str):
    out = {
        "hookSpecificOutput": {
            "hookEventName": "SessionStart",
            "additionalContext": ctx,
        }
    }
    print(json.dumps(out, ensure_ascii=False))

if not os.path.exists(FOCUS_FILE):
    emit(f"[CURRENT_FOCUS] file not found at {FOCUS_FILE} — please create per CLAUDE.md §6.")
    raise SystemExit(0)

with open(FOCUS_FILE, "r", encoding="utf-8") as f:
    content = f.read()

lines = content.splitlines()
# Match "## YYYY-MM-DD ..." allowing an optional leading marker/emoji token
# (e.g. "## 🔴 2026-05-22 ...") — fixes silent-skip of emoji-prefixed headers
# that previously caused the injector to grab an OLDER plain-date section.
date_header = re.compile(r"^##\s+(?:\S+\s+)?\d{4}-\d{2}-\d{2}")

# Find the first dated section "## YYYY-MM-DD" — that's the latest (newest at top per convention)
start = None
for i, line in enumerate(lines):
    if date_header.match(line):
        start = i
        break

if start is None:
    section = "\n".join(lines[:50])
else:
    # Capture from start until next "## " header or 80 lines
    end = min(start + 80, len(lines))
    for j in range(start + 1, end):
        if lines[j].startswith("## ") and date_header.match(lines[j]):
            end = j
            break
    section = "\n".join(lines[start:end])

# Truncate
truncated = section[:MAX_CHARS]
if len(section) > MAX_CHARS:
    truncated += "\n\n[...truncated, full content at InterSubMod/docs/CURRENT_FOCUS.md]"

# Staleness check
mtime = os.path.getmtime(FOCUS_FILE)
age_days = int((time.time() - mtime) / 86400)
stale_warn = ""
if age_days > 7:
    stale_warn = (
        f"\n\n⚠ [STALE WARNING] CURRENT_FOCUS.md last updated {age_days} days ago "
        f"(CLAUDE.md §6 contract: weekly cadence). Trigger /pivot-direction "
        f"or update before assuming current state."
    )

ctx = (
    "[CURRENT_FOCUS — latest section, auto-injected by SessionStart hook]\n\n"
    + truncated
    + stale_warn
    + "\n\n> Full file: InterSubMod/docs/CURRENT_FOCUS.md | "
    "Skip detailed load if: pure code edits / single doc / simple Q&A. "
    "Trigger /research-context-loader for Tier 2/3 depth."
)

emit(ctx)
PYEOF
