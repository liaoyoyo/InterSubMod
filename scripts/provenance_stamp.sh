#!/usr/bin/env bash
# provenance_stamp.sh — 印出可貼進報告 metadata 的 provenance stamp。
# 用途：盤點/audit/status 報告首行標 branch+HEAD+worktree，消滅 5-worktree cross-branch 幻覺（改進 ①，pitfalls P-17）。
# Read-only：不寫任何檔，只 echo。用法：bash scripts/provenance_stamp.sh   （在報告所在 repo/worktree 跑）
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
branch="$(git -C "$repo_root" rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'NOT-A-GIT-REPO')"
commit="$(git -C "$repo_root" rev-parse --short HEAD 2>/dev/null || echo 'NONE')"
wt="$repo_root"

# 並行 session 偵測（主 repo transcript dir >1 活躍 .jsonl 的粗略提示；fail-OPEN）
concurrent_note=""
n_wt="$(git -C "$repo_root" worktree list 2>/dev/null | wc -l | tr -d ' ' || echo 1)"
if [ "${n_wt:-1}" -gt 1 ]; then
  concurrent_note="  # ⚠ ${n_wt} worktree 並行 — 本 stamp 的 current 僅相對 ${branch}"
fi

cat <<EOF
build_branch: ${branch}${concurrent_note}
build_commit: ${commit}
worktree: ${wt}
EOF
