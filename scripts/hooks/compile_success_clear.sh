#!/bin/bash
# compile_success_clear.sh - PostToolUse hook (matcher: Bash)
# make 成功後清除「待編譯」標記

INPUT=$(cat)
CMD=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)
EXIT_CODE=$(echo "$INPUT" | jq -r '.tool_result.exit_code // "1"' 2>/dev/null)

[ -z "$CMD" ] && exit 0

# 偵測 make 命令且成功
if echo "$CMD" | grep -qE '(^|\s|&&\s*)make(\s|$)'; then
    if [ "$EXIT_CODE" = "0" ]; then
        MARKER="/tmp/ism_cpp_pending_compile.txt"
        if [ -f "$MARKER" ]; then
            rm -f "$MARKER"
            echo "[Compile Guard] 編譯成功，待編譯標記已清除。可以 commit。"
        fi
    fi
fi
