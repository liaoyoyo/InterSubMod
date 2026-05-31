#!/bin/bash
# compile_success_clear.sh - PostToolUse hook (matcher: Bash)
# make 成功後清除「待編譯」標記

INPUT=$(cat)
CMD=$(echo "$INPUT" | jq -r '.tool_input.command // empty' 2>/dev/null)

# 2026-05-31 fix: Claude Code PostToolUse 的結果欄位是 `.tool_response`（與
# external_input_sanitizer.sh 一致），先前誤寫 `.tool_result.exit_code` → jq 永遠
# fallback 回 "1" → marker 從不被清 → 自 2026-05-10 起 stale 累積（已知 lifecycle bug）。
# 改用 .tool_response + 多重 success 訊號 cascade，避免再次 stale-forever。
EXIT_CODE=$(echo "$INPUT" | jq -r '.tool_response.exit_code // .tool_response.exitCode // empty' 2>/dev/null)
IS_ERROR=$(echo "$INPUT" | jq -r '(.tool_response.isError // .tool_response.is_error // empty)|tostring' 2>/dev/null)

[ -z "$CMD" ] && exit 0

# 偵測 make 命令
if echo "$CMD" | grep -qE '(^|\s|&&\s*)make(\s|$)'; then
    # 判斷成功：authoritative 訊號優先，無訊號則 best-effort 清（make 已執行；
    # stale-marker-forever 的危害大於偶爾在失敗 make 後清除 — 下次 C++ edit 會重 arm）。
    SUCCESS=0
    if [ "$EXIT_CODE" = "0" ]; then SUCCESS=1
    elif [ -n "$EXIT_CODE" ]; then SUCCESS=0          # 明確 non-zero exit
    elif [ "$IS_ERROR" = "false" ]; then SUCCESS=1
    elif [ "$IS_ERROR" = "true" ]; then SUCCESS=0
    else SUCCESS=1; fi                                # 無可判訊號 → best-effort

    if [ "$SUCCESS" = "1" ]; then
        MARKER="/tmp/ism_cpp_pending_compile.txt"
        if [ -f "$MARKER" ]; then
            rm -f "$MARKER"
            echo "[Compile Guard] 編譯成功，待編譯標記已清除。可以 commit。"
        fi
    fi
fi
