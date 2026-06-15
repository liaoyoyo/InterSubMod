#!/bin/bash
# search_scope_guard.sh — PreToolUse hook (matcher: Bash)
#
# Purpose: This repo is large (build/output/data/_deps/research/figures/assets +
#          3 git worktrees, all .gitignore'd). The Grep/Glob TOOLS (ripgrep) respect
#          .gitignore and auto-skip those. But Bash `grep -r` / `find` / `du` do NOT
#          respect .gitignore — a bare repo-root recursive scan walks gigabytes of
#          data/output/worktree copies = slow + token waste ("卡住的動作").
#
#          This guard BLOCKS the 3 unambiguous repo-root-wide patterns and tells the
#          caller to (1) use the Grep/Glob tool, (2) scope to a subdir, or (3) add
#          -maxdepth. Scoped Bash searches (grep -r X src/, find docs/ ...) pass.
#
# Source:  2026-05-30 user request — judge-scope-first + pre-exclude heavy dirs.
# Exit:    2 = block (stronger restriction, per user); 0 = allow.

INPUT=$(cat 2>/dev/null)
CMD=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)
[ -z "$CMD" ] && exit 0

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
block=0
reason=""

# Allowlist escape: explicit opt-in marker lets a deliberate full scan through
echo "$CMD" | grep -q "ALLOW_FULL_SCAN" && exit 0

# 1) `find .` or `find <repo-root>` recursive without -maxdepth 1/2 and without -prune/-path excludes
if echo "$CMD" | grep -qE '\bfind\s+(\.|'"$(echo "$ROOT" | sed 's/\//\\\//g')"')(\s|$)'; then
    if ! echo "$CMD" | grep -qE '\-maxdepth\s+[0-2]\b' && ! echo "$CMD" | grep -qE '\-(prune|path)\b'; then
        block=1
        reason="find 從 repo root 遞迴（不尊重 .gitignore，會掃 build/output/data/research/worktrees 等 GB 級目錄）"
    fi
fi

# 2) `grep -r/-R` targeting '.' (repo root) — ignores .gitignore
if echo "$CMD" | grep -qE '\bgrep\b[^|]*\s-[A-Za-z]*[rR]' && echo "$CMD" | grep -qE '\s\.(\s|$)'; then
    block=1
    reason="grep -r 對 '.'（repo root，不尊重 .gitignore）"
fi

# 3) `du */` or `du` over repo-root glob — walks everything
if echo "$CMD" | grep -qE '\bdu\b[^|]*\s\*/' || echo "$CMD" | grep -qE '\bdu\b\s+-[a-z]*\s+\.(\s|$)'; then
    block=1
    reason="du 掃全 repo 目錄（research/output 達 GB 級，極慢）"
fi

# 4) walking commands (du / find / ls -R / wc -l) over the HEAVY OUTPUT roots.
#    ISM region dirs (big7_disk_output/.../intersubmod_tp, filtered_snv_tp) hold >10,000
#    subdirs each → these HANG. This is the documented "卡住" #1 cause (2026-06-15 audit
#    D3-1: a du on a complete_matrix dir hung this very session; rules 1-3 above only catch
#    repo-root `.`/`du */`, NOT an absolute big7_disk_output / output/ / data/ path).
HEAVY='big7_disk_output|(^|[[:space:]/])output/|(^|[[:space:]/])data/'
if echo "$CMD" | grep -qE "$HEAVY"; then
    if echo "$CMD" | grep -qE '\bdu\b'; then
        block=1
        reason="du 掃 output/data 重目錄（ISM region 子目錄 >1 萬，會 hang）— 改 ls -la *.bam / ls -d ...*/ | wc -l"
    elif echo "$CMD" | grep -qE '\bfind\b' && ! echo "$CMD" | grep -qE '\-maxdepth\s+[0-9]' && ! echo "$CMD" | grep -qE '\-(prune|path)\b'; then
        block=1
        reason="find 遞迴 output/data 重目錄（無 -maxdepth/-prune，會 hang）"
    elif echo "$CMD" | grep -qE '\bls\b[^|]*\s-[A-Za-z]*R'; then
        block=1
        reason="ls -R output/data 重目錄（遞迴列 >1 萬子目錄，會 hang）— 改單層 ls -d .../*/"
    elif echo "$CMD" | grep -qE '\bwc\b[^|]*\s-[a-z]*l' && echo "$CMD" | grep -qE '\*/'; then
        block=1
        reason="wc -l 對 output glob 展開 region 子目錄（會 hang）— 先 ls -d 一層再數"
    fi
fi

if [ "$block" -eq 1 ]; then
    {
        echo "[search-scope-guard] ⛔ 阻擋大規模搜尋：${reason}"
        echo "  改用以下之一（判斷先 → 精準搜尋）："
        echo "  1. **Grep / Glob 工具**（ripgrep，自動尊重 .gitignore，已跳過 9 重目錄）— 廣搜首選"
        echo "  2. scope 到輕目錄：grep -r <pat> .claude/skills src include docs scripts state tests"
        echo "  3. 淺掃加 -maxdepth 1（或 2）"
        echo "  4. 真要全掃：在指令加 'ALLOW_FULL_SCAN' 註記顯式 opt-in"
    } >&2
    exit 2
fi

exit 0
