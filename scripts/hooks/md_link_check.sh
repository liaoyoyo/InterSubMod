#!/usr/bin/env bash
# Markdown Link Check Hook — PostToolUse on Edit|Write
# Reads tool_input from stdin (JSON). When the edited file is .md, runs
# check_md_links.py against it and prints any broken relative refs.
# Exit 0 always (advisory only, never blocks).

set -euo pipefail

INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

[ -z "$FILE_PATH" ] && exit 0

# Only scan Markdown files
case "$FILE_PATH" in
    *.md|*.MD|*.markdown) ;;
    *) exit 0 ;;
esac

# Skip if file no longer exists (e.g., rename in-flight)
[ -f "$FILE_PATH" ] || exit 0

python3 /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/check_md_links.py "$FILE_PATH" 2>/dev/null || true

exit 0
