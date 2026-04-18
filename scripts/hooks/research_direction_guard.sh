#!/bin/bash
# research_direction_guard.sh - UserPromptSubmit hook
# 偵測用戶提問是否涉及已關閉(NO-GO/NEGATIVE)的研究方向，發出警告但不阻擋
#
# === 擴展指南 ===
# 當新的研究方向被判定 NO-GO 或 NEGATIVE 時，在下方 WARNINGS 區塊新增條目：
#   1. 複製現有的 if-fi 區塊作為模板
#   2. grep -qiE 的 regex 填入該方向的關鍵字（中英文皆可）
#   3. WARNINGS 第一行：[NO-GO] 或 [NEGATIVE] + 一句話結論 + (日期 + 來源)
#   4. WARNINGS 第二行：關鍵數據（AUC、delta、樣本數）
#   5. WARNINGS 第三行：參考 memory 檔名
#   6. 完成後在 validation-protocol 結論生命週期 Checklist 勾選此步驟
# ================

INPUT=$(cat)
PROMPT=$(echo "$INPUT" | jq -r '.user_prompt // empty' 2>/dev/null)
[ -z "$PROMPT" ] && PROMPT="$INPUT"

WARNINGS=()

# 1. 純甲基化特徵篩選（Beyond-AUC 2026-04-09 確認耗盡）
if echo "$PROMPT" | grep -qiE '(pure methylation|純甲基化|methylation feature|甲基化特徵).*(filter|篩選|區分|discriminat|AUC)'; then
    WARNINGS+=("[NO-GO] 純甲基化特徵空間已確認耗盡 (2026-04-09 Beyond-AUC)")
    WARNINGS+=("  所有 pure methylation AUC <= 0.58，7 方法全面驗證")
    WARNINGS+=("  參考: memory/project_beyond_auc_exhaustion_confirmed.md")
fi

# 2. TO Germline FP 鑑別（G1-G7 全 NO-GO）
if echo "$PROMPT" | grep -qiE '(germline FP|TO.*(post-hoc|特徵探索|feature explor)|TO.*FP.*(filter|篩選))'; then
    WARNINGS+=("[NO-GO] TO Germline FP 鑑別已關閉 (G1-G7 七模組)")
    WARNINGS+=("  48 圖 60+ 特徵全 AUC<0.64，TP loss<=2% 下 FP removal=0%")
    WARNINGS+=("  參考: memory/project_germline_fp_identification_nogo.md")
fi

# 3. O11 甲基化 Heterogeneity
if echo "$PROMPT" | grep -qiE '(epipolymorphism|within-group heterogeneity|甲基化異質性.*區分)'; then
    WARNINGS+=("[NEGATIVE] O11 甲基化 Heterogeneity 假說已否決")
    WARNINGS+=("  epipolymorphism AUC 0.845→0.530 after n_reads correction")
    WARNINGS+=("  參考: memory/project_O11_heterogeneity_negative.md")
fi

# 4. O12 LOH 甲基化場景
if echo "$PROMPT" | grep -qiE '(LOH.*methylation.*scenario|LOH.*甲基化場景|AlleleDelta.*AF)'; then
    WARNINGS+=("[NEGATIVE] O12 LOH 甲基化場景不可區分")
    WARNINGS+=("  AlleleDelta=AF confound，L2 collider bias，全 corrected AUC<0.58")
    WARNINGS+=("  參考: memory/project_O12_loh_methylation_scenarios.md")
fi

# 5. O13 跨區域甲基化 Correlation
if echo "$PROMPT" | grep -qiE '(cross-region|跨區域.*correlation|multi-site methylation|多點甲基化)'; then
    WARNINGS+=("[NEGATIVE] O13 跨區域甲基化 Correlation 已關閉")
    WARNINGS+=("  FP>TP 是 shared read count confound，分層+matching 後差異消失")
    WARNINGS+=("  參考: memory/project_O13_cross_region_negative.md")
fi

# 6. Fine-Pairwise Distance
if echo "$PROMPT" | grep -qiE '(fine.*pairwise|pairwise.*distance.*filter|距離.*篩選)'; then
    WARNINGS+=("[NEGATIVE] Fine-Pairwise Distance 已確認無效")
    WARNINGS+=("  6 pairwise 距離全無效，Paired AUC<0.50（反轉）")
    WARNINGS+=("  參考: memory/project_fine_pairwise_negative.md")
fi

# 7. Option C 雙路測試
if echo "$PROMPT" | grep -qiE '(HP-free combo|option.?C|雙路測試|HP.*free.*AUC)'; then
    WARNINGS+=("[NEGATIVE] Option C 雙路測試已確認無效")
    WARNINGS+=("  HP-free combo AUC=0.564，所有信號來自 HP tags")
    WARNINGS+=("  參考: memory/project_option_c_dual_path_negative.md")
fi

# 8. TO QS 隨機
if echo "$PROMPT" | grep -qiE '(TO.*quality.?score|TO.*QS|TO.*品質分數).*(filter|篩選|有效)'; then
    WARNINGS+=("[NO-GO] TO Quality Score AUC=0.497（隨機水準）")
    WARNINGS+=("  LOH penalty 反向，需根本重設計")
    WARNINGS+=("  參考: memory/project_qs_redesign_evidence.md")
fi

# 僅在有匹配時輸出
if [ ${#WARNINGS[@]} -gt 0 ]; then
    echo "================================================================"
    echo "[Research Direction Guard] 偵測到已關閉的研究方向："
    echo ""
    for w in "${WARNINGS[@]}"; do
        echo "  $w"
    done
    echo ""
    echo "如果你正在探索新的變體方向，請明確說明與已關閉方向的差異。"
    echo "================================================================"
fi

exit 0
