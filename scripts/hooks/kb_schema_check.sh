#!/bin/bash
# kb_schema_check.sh - PreToolUse hook (matcher: Bash)
# git commit 前若 staged 檔案含 knowledge/**/*.md，強制跑 3 驗證腳本
# 任一失敗 → exit 2 阻擋

INPUT=$(cat)
CMD=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)

[ -z "$CMD" ] && exit 0

# 只攔截 git commit 命令
if ! echo "$CMD" | grep -qE 'git\s+commit'; then
    exit 0
fi

KB_ROOT=/big7_disk/liaoyoyo2001/InterSubMod/knowledge
PYBIN=/usr/bin/python3

# 若當前無 staged 變動 / 無 knowledge 變動，直接放行
STAGED=$(git -C /big7_disk/liaoyoyo2001/InterSubMod diff --cached --name-only 2>/dev/null)
if [ -z "$STAGED" ]; then
    exit 0
fi
if ! echo "$STAGED" | grep -q '^knowledge/'; then
    exit 0
fi

# 跑 3 驗證腳本
FAIL=0
for script in validate_frontmatter.py check_related_ids_symmetry.py check_canonical_paths.py; do
    OUTPUT=$($PYBIN "$KB_ROOT/scripts/$script" "$KB_ROOT" 2>&1)
    EXIT=$?
    if [ $EXIT -ne 0 ]; then
        echo "[KB SCHEMA GUARD - BLOCKED] scripts/$script FAIL:"
        echo "$OUTPUT" | tail -10
        FAIL=1
    fi
done

if [ $FAIL -eq 1 ]; then
    echo ""
    echo "KB 變動含 schema / related_ids / canonical_paths 違規，禁止 commit。"
    echo "修正後再試。可手動跑："
    echo "  $PYBIN $KB_ROOT/scripts/fix_related_ids_symmetry.py $KB_ROOT"
    exit 2
fi

exit 0
