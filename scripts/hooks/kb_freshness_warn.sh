#!/bin/bash
# kb_freshness_warn.sh - UserPromptSubmit hook
# 掃描 10_research_status/ 與 09_conclusions/05_hypothesis_queue_snapshot.md
# 若 last_verified > 14 天 → 印提醒（不阻擋）

KB_ROOT=/big7_disk/liaoyoyo2001/InterSubMod/knowledge
TODAY=$(date +%s)
FOURTEEN_DAYS=$((14 * 24 * 3600))

# 關注的時間敏感檔案 pattern
TARGETS=(
    "$KB_ROOT/10_research_status/01_current_focus_snapshot.md"
    "$KB_ROOT/10_research_status/02_active_hypotheses.md"
    "$KB_ROOT/10_research_status/04_blockers_and_risks.md"
    "$KB_ROOT/10_research_status/05_next_milestones.md"
    "$KB_ROOT/09_conclusions/05_hypothesis_queue_snapshot.md"
)

STALE=()
for f in "${TARGETS[@]}"; do
    [ ! -f "$f" ] && continue
    LV=$(grep -m1 '^last_verified:' "$f" | sed 's/last_verified:\s*//; s/\s*$//')
    [ -z "$LV" ] && continue
    LV_EPOCH=$(date -d "$LV" +%s 2>/dev/null)
    [ -z "$LV_EPOCH" ] && continue
    AGE=$((TODAY - LV_EPOCH))
    if [ $AGE -gt $FOURTEEN_DAYS ]; then
        DAYS=$((AGE / 86400))
        STALE+=("  - $(basename "$f"): ${DAYS} 天前 ($LV)")
    fi
done

if [ ${#STALE[@]} -gt 0 ]; then
    echo "[KB Freshness 警示] 以下時間敏感文件超過 14 天未驗證，請考慮更新："
    printf '%s\n' "${STALE[@]}"
    echo "  更新後跑: $KB_ROOT/scripts/refresh_last_verified.py $KB_ROOT --file <relative-path>"
fi

exit 0
