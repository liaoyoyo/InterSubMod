#!/bin/bash
# precompact_autosave.sh — PreCompact hook
#
# Purpose: Before Claude Code compacts the conversation, snapshot the volatile
#          working state (active cycles + CURRENT_FOCUS tail + open git status)
#          so that if compaction drops detail, a durable on-disk record remains.
#          Directly supports CLAUDE.md §7 "/compact 必保留" requirements.
#
# Source:  Community gap G6 — claude-research precompact-autosave.py / postcompact-restore.py.
#          InterSubMod previously had only compact_test.sh (manual verification),
#          no automatic pre-compaction state save.
# Trigger: PreCompact event (fires before both manual /compact and auto-compact).
# Exit:    Always 0 (advisory, never blocks compaction).

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
SNAP_DIR="${ROOT}/state/compact_snapshots"
TS=$(date +%Y%m%d_%H%M%S)
SNAP="${SNAP_DIR}/precompact_${TS}.md"

mkdir -p "$SNAP_DIR" 2>/dev/null || exit 0

# Read hook stdin (may contain trigger=manual|auto, custom_instructions)
INPUT=$(cat 2>/dev/null)
TRIGGER=$(echo "$INPUT" | jq -r '.trigger // "unknown"' 2>/dev/null || echo "unknown")

{
    echo "# PreCompact Snapshot — ${TS} (trigger=${TRIGGER})"
    echo ""
    echo "> Auto-saved by precompact_autosave.sh before context compaction."
    echo "> Restore hint: read this file if post-compact context lost active cycle / focus detail."
    echo ""

    echo "## active.json"
    echo '```json'
    if [ -f "${ROOT}/state/active.json" ]; then
        cat "${ROOT}/state/active.json" 2>/dev/null
    else
        echo "(no state/active.json)"
    fi
    echo '```'
    echo ""

    echo "## CURRENT_FOCUS.md (last 60 lines)"
    echo '```markdown'
    if [ -f "${ROOT}/docs/CURRENT_FOCUS.md" ]; then
        tail -n 60 "${ROOT}/docs/CURRENT_FOCUS.md" 2>/dev/null
    else
        echo "(no CURRENT_FOCUS.md)"
    fi
    echo '```'
    echo ""

    echo "## git status (porcelain)"
    echo '```'
    git -C "$ROOT" status --porcelain 2>/dev/null | head -40 || echo "(git unavailable)"
    echo '```'
    echo ""

    echo "## active cycle state files"
    if [ -d "${ROOT}/state/cycles" ]; then
        find "${ROOT}/state/cycles" -maxdepth 2 -name "state.json" 2>/dev/null | while read -r f; do
            echo "- ${f#$ROOT/}"
        done
    else
        echo "(no state/cycles dir)"
    fi
} > "$SNAP" 2>/dev/null

# Retention: keep last 20 snapshots, prune older
ls -1t "${SNAP_DIR}"/precompact_*.md 2>/dev/null | tail -n +21 | while read -r old; do
    rm -f "$old" 2>/dev/null
done

echo "[PreCompact] state snapshot saved → state/compact_snapshots/precompact_${TS}.md (trigger=${TRIGGER})"
exit 0
