#!/bin/bash
# pre_commit_compile_check.sh - PreToolUse hook (matcher: Bash)
# git commit 前檢查是否有未編譯的 C++ 修改，有則阻擋

INPUT=$(cat)
CMD=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)

[ -z "$CMD" ] && exit 0

# 只攔截 git commit 命令
if echo "$CMD" | grep -qE 'git\s+commit'; then
    MARKER="/tmp/ism_cpp_pending_compile.txt"
    if [ -f "$MARKER" ] && [ -s "$MARKER" ]; then
        echo "[COMPILE GUARD - BLOCKED] 有未編譯的 C++ 修改，禁止 commit。"
        echo ""
        echo "未編譯檔案:"
        while IFS= read -r f; do
            echo "  - $f"
        done < "$MARKER"
        echo ""
        echo "請先執行: cd build && make -j\$(nproc)"
        echo "編譯成功後再重新 commit。"
        exit 2  # exit 2 = 阻擋工具執行
    fi
fi

exit 0
