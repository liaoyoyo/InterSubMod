#!/bin/bash
# memory_recall_logger.sh — Track which skill / memory files are actually referenced.
#
# Purpose: 量化 each skill / memory file 的 monthly recall rate；<5%/month → archive 候選
# Source:  P4 Architect E9 / C1 mitigation (44 skills + 16 plugin = 60 trigger space)
# Trigger: PostToolUse on Read (records which skill/memory file was actually read)
# Output:  docs/postmortems/recall_log_YYYYMM.log (tab-separated, monthly rotation)
# Query:   bash InterSubMod/scripts/hooks/memory_recall_logger.sh --report
#          → produces summary of skill/memory recall counts for the month

set -o pipefail

# Report mode
if [ "${1:-}" = "--report" ]; then
    YYYYMM="${2:-$(date +%Y%m)}"
    LOG="/big7_disk/liaoyoyo2001/InterSubMod/docs/postmortems/recall_log_${YYYYMM}.log"
    if [ ! -f "$LOG" ]; then
        echo "No recall log for $YYYYMM at $LOG"
        exit 0
    fi
    echo "# memory_recall report — $YYYYMM"
    echo "Source: $LOG"
    echo "Total reads tracked: $(wc -l < "$LOG")"
    echo ""
    echo "## Top 20 most-recalled targets"
    awk -F'\t' '{print $2}' "$LOG" | sort | uniq -c | sort -rn | head -20
    echo ""
    echo "## Recall rate per skill (skill_count / total_reads)"
    TOTAL=$(wc -l < "$LOG")
    awk -F'\t' '$2 ~ /^skill=/ {print $2}' "$LOG" | sort | uniq -c | sort -rn | \
        awk -v total="$TOTAL" '{rate=$1*100/total; printf "  %s\t%d reads\t%.1f%%\n", $2, $1, rate}'
    echo ""
    echo "## Skills with 0 reads this month (archive candidates)"
    # All skills
    for sk in /big7_disk/liaoyoyo2001/InterSubMod/.claude/skills/*/; do
        name=$(basename "$sk")
        [ "$name" = "README.md" ] && continue
        cnt=$(grep -cF "skill=${name}" "$LOG" 2>/dev/null || true)
        cnt=${cnt:-0}
        if [ "$cnt" = "0" ]; then
            echo "  - $name (0 reads — archive candidate per C1 mitigation)"
        fi
    done
    exit 0
fi

# Log mode (PostToolUse Read)
INPUT=$(cat 2>/dev/null || echo "{}")
TOOL_NAME=$(echo "$INPUT" | jq -r '.tool_name // ""' 2>/dev/null)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // ""' 2>/dev/null)

[ "$TOOL_NAME" != "Read" ] && exit 0
[ -z "$FILE_PATH" ] && exit 0

# Detect target type
TARGET=""
case "$FILE_PATH" in
    */.claude/skills/*/SKILL.md)
        # Extract skill name
        SKILL_NAME=$(echo "$FILE_PATH" | sed -E 's|.*/\.claude/skills/([^/]+)/.*|\1|')
        TARGET="skill=${SKILL_NAME}"
        ;;
    */memory/feedback_*.md)
        MEM_NAME=$(basename "$FILE_PATH" .md)
        TARGET="memory=${MEM_NAME}"
        ;;
    */memory/project_*.md)
        MEM_NAME=$(basename "$FILE_PATH" .md)
        TARGET="memory=${MEM_NAME}"
        ;;
    */memory/reference_*.md)
        MEM_NAME=$(basename "$FILE_PATH" .md)
        TARGET="memory=${MEM_NAME}"
        ;;
    *)
        exit 0  # Not a tracked target
        ;;
esac

# Append
YYYYMM=$(date +%Y%m)
LOG="/big7_disk/liaoyoyo2001/InterSubMod/docs/postmortems/recall_log_${YYYYMM}.log"
mkdir -p "$(dirname "$LOG")" 2>/dev/null
printf '%s\t%s\n' "$(date -Iseconds)" "$TARGET" >> "$LOG"

exit 0
