#!/bin/bash
# git_branch_commit_guard.sh — PreToolUse Bash. Guards `git commit` (master-workflow §3-C/§4③).
#
#   1) exit-2 BLOCK a direct commit to a protected trunk (main / master / develop).
#      This repo ALWAYS works on feature branches and merges to trunk only after
#      user confirmation — a direct trunk commit is essentially always a mistake.
#   2) advisory WARN (never blocks) if other Claude sessions are active in the main
#      repo dir — a commit could land on a shared HEAD / wrong branch (2026-06-09).
#
# fail-OPEN (rigorous): any parse error, missing git, or non-commit command -> exit 0.
# The matcher requires `commit` to be the git SUBCOMMAND (immediately after `git`
# or after `-C <dir>`), so `git log ...commit...` etc. are NOT treated as commits.
# Bias is toward NOT blocking: worst case it misses a block, never over-blocks.

set -uo pipefail

REPO=/big7_disk/liaoyoyo2001/InterSubMod
TDIR=/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod

CMD="$(cat 2>/dev/null | python3 -c 'import sys,json
try:
    print(json.load(sys.stdin).get("tool_input",{}).get("command",""))
except Exception:
    print("")' 2>/dev/null || true)"

# Act ONLY when `commit` is the git subcommand (optionally after `-C <dir>`).
echo "$CMD" | grep -qE 'git([[:space:]]+-C[[:space:]]+[^[:space:]]+)?[[:space:]]+commit([[:space:]]|$)' || exit 0

# Resolve the target dir (worktree commits use `git -C <dir> commit`).
DIR="$(echo "$CMD" | grep -oE '\-C[[:space:]]+[^[:space:]]+' | head -1 | awk '{print $NF}')"
[ -z "${DIR:-}" ] && DIR="$REPO"
[ -d "$DIR" ] || DIR="$REPO"
BR="$(git -C "$DIR" branch --show-current 2>/dev/null || echo '')"

# 1) protected-trunk block
case "${BR:-}" in
  main|master|develop)
    echo "[branch-guard] ⛔ 阻擋：直接 commit 到受保護主線 '$BR'。本 repo 一律用 feature branch，主線只在用戶確認後 merge（master-workflow §3-C/§3-D）。請切到/開 feature branch 再 commit。" 1>&2
    exit 2
    ;;
esac

# 2) advisory concurrent-session warn (never blocks; fail-OPEN on any error)
OTHERS="$(python3 -c "
import os,time,glob
now=time.time()
print(sum(1 for f in glob.glob('$TDIR/*.jsonl') if (now-os.path.getmtime(f))<300))
" 2>/dev/null || echo 0)"
if [ "${OTHERS:-0}" -gt 1 ] 2>/dev/null; then
  echo "[branch-guard] ⚠ 偵測到同 repo 約 $((OTHERS-1)) 個並行 session 近期活躍；commit 前請確認你在自己的 worktree/branch（共用主 dir+HEAD 會讓 commit 落錯 branch，見 master-workflow §3-B）。" 1>&2
fi
exit 0
