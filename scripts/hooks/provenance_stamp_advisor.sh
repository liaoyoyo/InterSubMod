#!/usr/bin/env bash
# provenance_stamp_advisor.sh — PostToolUse(Write|Edit) advisory.
#
# WHY (2026-06-15 audit D6-2): scripts/provenance_stamp.sh existed but was manual-only — even
# the audit doc that DEFINES the stamp rule shipped without it (silent skip → P-17 cross-worktree
# hallucination). This NUDGES (advisory, never blocks) when an audit/status/盤點-class doc is
# written WITHOUT a `build_branch:` stamp, printing the ready-to-paste stamp.
#
# Scope: ONLY fires on audit/status/manifest/盤點/稽核/inventory-class docs (by path/name) — NOT
#        every .md (keeps it quiet + cheap). Exit ALWAYS 0 (advisory, fail-OPEN — never blocks a
#        Write, distinct from the exit-2 gates). Reuses scripts/provenance_stamp.sh for stamp text.
# Source: 2026-06-15 workflow audit D6-2; CLAUDE.md §13-A / known-pitfalls P-17.

INPUT=$(cat 2>/dev/null)
FP=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)
[ -z "$FP" ] && exit 0
case "$FP" in *.md|*.html) ;; *) exit 0 ;; esac

# audit/status-class detection: docs/data_specs/ path OR filename signals 盤點/audit/status/manifest
base=$(basename "$FP")
is_class=0
case "$FP" in */docs/data_specs/*) is_class=1 ;; esac
echo "$base" | grep -qiE 'audit|status|manifest|盤點|稽核|inventory' && is_class=1
[ "$is_class" -eq 0 ] && exit 0

# already stamped? (build_branch: in first 20 lines) → silent
[ -f "$FP" ] && head -20 "$FP" 2>/dev/null | grep -q 'build_branch:' && exit 0

REPO="/big7_disk/liaoyoyo2001/InterSubMod"
{
  echo "[provenance-stamp-advisor] ℹ️ 盤點/audit/status 類文件缺 build_branch: stamp（防 P-17 跨-worktree 幻覺）。建議在 frontmatter 加："
  bash "$REPO/scripts/provenance_stamp.sh" 2>/dev/null | sed 's/^/    /'
} >&2
exit 0
