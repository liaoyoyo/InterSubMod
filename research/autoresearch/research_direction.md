<!--
建立時間: 2026-03-22
目標: 記錄 InterSubMod AutoResearch 的當前研究焦點與方向
處理範圍: 研究迴圈的方向文件，供 AI 模型定向用
-->

# 研究方向文件 (Research Direction)

> 最後更新：2026-03-22
> 此文件由 AI 與人類共同維護。AI 在每輪 research-loop 中讀取此文件定向。

---

## 當前研究焦點

**主要 Track**：TO（Tumor-Only）
**理由**：TO 有更多 FP（~300-800），三層架構改善空間比 paired 大，甲基訊號在 TO 下尚有未探索空間。

**主要資料集（pilot）**：HCC1395_5kHz_TO
**次要資料集（medium scale）**：HCC1395_DORADO_TO

**當前基線 F1**：
- HCC1395_5kHz_TO: 0.7127（ClairS-TO baseline）
- HCC1395_DORADO_TO: 0.7226（ClairS-TO baseline）

---

## 當前研究優先項目

1. **W1 設計弱點探索**：HPP 未整合進分類邏輯 → 測試 HPP threshold 作為 FP-B 過濾
2. **W2 設計弱點探索**：AlleleDelta 在 HP 不均時虛高 → 測試 HPP-corrected AlleleDelta
3. **FP 三分類驗證**：確認 HCC1395_TO 中 Type-A/B/C FP 的分佈
4. **TO 甲基訊號穩定性**：Quality_Score 與 PairwiseMedianDist 在 TO 的可靠性

---

## 已確認的研究邊界（不重複探索）

- **甲基作為獨立規則（無 caller gate）**：已驗證全數為負，不再探索
- **DORADO paired 甲基 support**：全數不成立，不再作為主要驗證方向
- **Pool B FN 甲基救援**：已 closeout，不作獨立規則

---

## 轉向歷史

*（此區塊由 pivot-direction skill 自動追加）*

---

## 注意事項

- 若連續 3 輪 `|delta_f1| < 0.001`，自動建議切換策略
- 每輪結束後，FEEDBACK 階段的用戶說法若改變方向，立即更新此文件
- 文獻引導特徵使用 `/inject-hypothesis` 加入佇列時，同步記錄來源於此
