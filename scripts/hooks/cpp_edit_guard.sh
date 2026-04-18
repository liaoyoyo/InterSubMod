#!/bin/bash
# cpp_edit_guard.sh - PostToolUse hook (matcher: Edit|Write)
# C++ 檔案修改後建立「待編譯」標記，提醒用戶需編譯

INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

[ -z "$FILE_PATH" ] && exit 0

if echo "$FILE_PATH" | grep -qE '\.(cpp|hpp|h)$'; then
    MARKER="/tmp/ism_cpp_pending_compile.txt"
    echo "$FILE_PATH" >> "$MARKER"
    # 去重
    sort -u "$MARKER" -o "$MARKER" 2>/dev/null
    COUNT=$(wc -l < "$MARKER" 2>/dev/null)
    echo "[Compile Guard] C++ 已修改: $(basename "$FILE_PATH") (待編譯: ${COUNT} 檔)"
    echo "  請在 commit 前執行: cd build && make -j\$(nproc)"
fi
