#!/bin/bash
# post_cpp_commit_invalidate.sh
# PostToolUse hook (matcher: Bash) — append to state/invalidation/binary_versions.jsonl
# when a `git commit` modifies any C++ source under src/ or include/.
#
# This hook records *that* a binary version changed; it does NOT scan
# downstream reports. Stale-mark scanning happens lazily inside
# `/check-staleness` skill (Phase 2 PRECHECK), which compares plan-stated
# binary_version against the HEAD recorded here.
#
# Triggers: any Bash invocation containing `git commit`.
# Action: if HEAD commit touches src/**/*.{cpp,hpp,h} or include/**/*.{cpp,hpp,h},
#         append a JSON object to InterSubMod/state/invalidation/binary_versions.jsonl.
#
# Exit: always 0 (notify-only, never block).

set -u

REPO_ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
LEDGER="${REPO_ROOT}/state/invalidation/binary_versions.jsonl"

INPUT=$(cat)
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)

# Bail unless this Bash call was a `git commit`
if ! echo "$COMMAND" | grep -qE '(^|[^a-zA-Z])git[[:space:]]+commit'; then
    exit 0
fi

# Resolve current HEAD; bail if unable
SHA=$(git -C "$REPO_ROOT" rev-parse --short=12 HEAD 2>/dev/null || true)
FULL_SHA=$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || true)
[ -z "$FULL_SHA" ] && exit 0

# Get list of files changed in this commit
CHANGED=$(git -C "$REPO_ROOT" show --name-only --format='' HEAD 2>/dev/null || true)
[ -z "$CHANGED" ] && exit 0

# Extract C++ files under src/ or include/ only (ignore tests/, scripts/)
CPP_FILES=$(echo "$CHANGED" | grep -E '^(src|include)/.*\.(cpp|hpp|h)$' || true)

if [ -z "$CPP_FILES" ]; then
    exit 0
fi

# Resolve metadata
TS=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
SUBJECT=$(git -C "$REPO_ROOT" log -1 --pretty=format:'%s' HEAD 2>/dev/null | head -c 200)
N_FILES=$(echo "$CPP_FILES" | wc -l | tr -d ' ')

# Compose JSON line. Use jq for proper escaping.
ENTRY=$(jq -n -c \
    --arg ts "$TS" \
    --arg sha "$FULL_SHA" \
    --arg short_sha "$SHA" \
    --arg subject "$SUBJECT" \
    --argjson n "$N_FILES" \
    --arg files "$CPP_FILES" \
    '{
        timestamp: $ts,
        commit: $sha,
        short_sha: $short_sha,
        subject: $subject,
        n_cpp_files_changed: $n,
        cpp_files: ($files | split("\n"))
    }')

# Ensure directory exists then append
mkdir -p "$(dirname "$LEDGER")"
echo "$ENTRY" >> "$LEDGER"

# Stderr notify (visible to Claude Code as hook output)
echo "[InvalidationHook] Recorded C++ binary commit ${SHA} (${N_FILES} files) → state/invalidation/binary_versions.jsonl" >&2
echo "  /check-staleness will use this to flag stale upstream reports at P2 PRECHECK." >&2

exit 0
