#!/bin/bash
# kb_sot_guard.sh - PostToolUse hook (matcher: Edit|Write)
# 偵測寫入內容含 ΔF1 SoT 受保護數字（+0.0112 / 0.9762 / -0.0206 / 0.9650）
# 且檔案非 03_pipelines/05_f1_baseline_canonical.md → 印提醒（不阻擋）

INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

[ -z "$FILE_PATH" ] && exit 0

# 僅對 knowledge/ 下 .md 檔觸發
case "$FILE_PATH" in
    */knowledge/*.md) ;;
    *) exit 0 ;;
esac

# SoT 本尊放行
case "$FILE_PATH" in
    */03_pipelines/05_f1_baseline_canonical.md) exit 0 ;;
    */CHANGELOG.md) exit 0 ;;
esac

# 取 new_string 或 content（依 tool 類型）
PAYLOAD=$(echo "$INPUT" | jq -r '.tool_input.new_string // .tool_input.content // empty' 2>/dev/null)
[ -z "$PAYLOAD" ] && exit 0

# 偵測受保護數字（ΔF1 locked values）
PROTECTED='[+]?0\.0112|[-]0\.0206|0\.9762|0\.9650'
if echo "$PAYLOAD" | grep -qE "$PROTECTED"; then
    echo "[KB SoT Guard] 偵測到 ΔF1 SoT 受保護數字寫入 $FILE_PATH"
    echo "  受保護：+0.0112 / 0.9762 / -0.0206 / 0.9650（CL-002/003 locked）"
    echo "  建議：連結 knowledge/03_pipelines/05_f1_baseline_canonical.md 而非複製數字"
    echo "  若為合理引用（速查/對照），可忽略此提醒。"
fi

exit 0
