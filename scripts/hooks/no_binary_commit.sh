#!/usr/bin/env bash
# no_binary_commit.sh — PreToolUse hook for git commit, blocks binary files
#
# Purpose: Enforce binary-files-not-in-git policy after 2026-05-11 history rewrite.
# Industry standard: PPTX, PDF, JPG/PNG, large HTML are deliverables/generated,
# not source. They belong in Drive/LFS/CDN, not git.
#
# Triggers when Claude Code attempts `git commit` (PreToolUse matcher: Bash).
# Reads the staged changes; if any blocked extension is staged, exits 2 (block).
#
# Allowed exceptions (small static assets):
#   - logos/icons/favicons (small, stable, in specific dirs)
#   - .gitkeep placeholders
#
# Override: prefix the commit command with `BYPASS_BINARY_CHECK=1` if intentional.

set -uo pipefail

# Read JSON payload from stdin
INPUT=$(cat)

# Extract command (best-effort; if not git commit, skip)
COMMAND=$(echo "$INPUT" | jq -r '.tool_input.command // ""' 2>/dev/null)

# Only trigger on git commit (not amend, not status)
case "$COMMAND" in
  *"git commit"*) ;;
  *) exit 0 ;;
esac

# Honor explicit bypass
if [ "${BYPASS_BINARY_CHECK:-0}" = "1" ]; then
  exit 0
fi

# Find staged binary files matching blocked patterns
# Blocked: pptx, pdf, jpg, jpeg, png (in figures/assets), >500KB anything
#
# Exception (added 2026-08-06): docs/images/ holds the small, hand-made explanatory
# figures embedded in README.md / README.zh-TW.md and the GitHub Wiki. Those MUST be
# version-controlled or GitHub renders broken images — a raw.githubusercontent URL can
# only resolve against a committed file. They are generated deterministically by
# tools/extract_svg_for_github.py from docs/explain/, are kept small (all well under the
# 500KB warn threshold), and the same carve-out exists in .gitignore. Research output
# figures remain blocked everywhere else, which is what this policy is actually for.
BLOCKED=$(git diff --cached --name-only --diff-filter=AM 2>/dev/null | \
  grep -iE '\.(pptx|pdf|jpg|jpeg|png)$' | \
  grep -v '\.gitkeep$' | \
  grep -v '^docs/images/[^/]*$' || true)

if [ -n "$BLOCKED" ]; then
  echo "ERROR: Blocked binary files in staged changes:" >&2
  echo "$BLOCKED" >&2
  echo "" >&2
  echo "Per InterSubMod policy (post 2026-05-11 history cleanup):" >&2
  echo "  - PPTX/PDF: store in Google Drive, link from .md" >&2
  echo "  - PNG/JPG figures: regenerate via scripts; gitignore figures/" >&2
  echo "" >&2
  echo "Override: BYPASS_BINARY_CHECK=1 git commit -m '...'" >&2
  exit 2
fi

# Also block any large file (>500KB) regardless of type
LARGE=$(git diff --cached --name-only --diff-filter=AM 2>/dev/null | while read -r f; do
  if [ -f "$f" ]; then
    size=$(stat -c%s "$f" 2>/dev/null || echo 0)
    if [ "$size" -gt 524288 ]; then
      printf "%s\t%dKB\n" "$f" "$((size/1024))"
    fi
  fi
done)

if [ -n "$LARGE" ]; then
  echo "WARN: Large files staged (>500KB):" >&2
  echo "$LARGE" >&2
  echo "" >&2
  echo "Consider Git LFS or external storage." >&2
  echo "Proceeding (warning only). Override BYPASS_BINARY_CHECK=1 to silence." >&2
  # Warn only, don't block
fi

exit 0
