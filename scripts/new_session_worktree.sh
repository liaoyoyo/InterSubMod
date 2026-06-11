#!/bin/bash
# new_session_worktree.sh <topic> [base-branch] — master-workflow §3-B mechanism.
#
# Create an ISOLATED git worktree for a parallel Claude session so concurrent
# sessions never share the main repo dir's HEAD (the 2026-06-09 cross-session
# commit-pollution fix). Run THIS, then start the parallel `claude` session inside
# the printed worktree dir — NOT in the main repo dir.
#
# Why a launcher (not just `git worktree add`): it also copies the gitignored local
# config (.claude/settings.local.json = your allow-list/permissions) into the new
# worktree so the isolated session is actually runnable, not a permission-prompt shell.
#
# (If your Claude Code version has the native `claude --worktree <topic>` flag, that
#  is even simpler — this script is the version-independent, always-works fallback.)

set -euo pipefail
REPO=/big7_disk/liaoyoyo2001/InterSubMod
TOPIC="${1:?usage: new_session_worktree.sh <topic> [base-branch]   e.g. new_session_worktree.sh asm-rerun}"
BASE="${2:-HEAD}"
SLUG="$(echo "$TOPIC" | tr 'A-Z ' 'a-z-' | tr -cd 'a-z0-9-')"
BR="session/${SLUG}-$(date +%Y%m%d)"
WT="${REPO}/.claude/worktrees/${SLUG}"

if git -C "$REPO" worktree list --porcelain 2>/dev/null | grep -q "worktree ${WT}$"; then
  echo "⚠ worktree already exists: $WT"; exit 0
fi

git -C "$REPO" worktree add -b "$BR" "$WT" "$BASE"

# bring gitignored local config so the isolated session is runnable
if [ -f "$REPO/.claude/settings.local.json" ]; then
  cp "$REPO/.claude/settings.local.json" "$WT/.claude/settings.local.json" 2>/dev/null || true
fi

echo ""
echo "✅ isolated worktree ready:"
echo "   dir   : $WT"
echo "   branch: $BR  (base: $BASE)"
echo ""
echo "→ 在這裡起並行 session（不要在主 repo dir 起）："
echo "     cd $WT && claude"
echo "→ 完成後清理： git -C $REPO worktree remove $WT && git -C $REPO branch -d $BR"
