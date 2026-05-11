<!--
建立時間: 2026-03-22
目標: 記錄 InterSubMod AutoResearch 的當前研究焦點與方向
處理範圍: 研究迴圈的方向文件，供 AI 模型定向用
-->

# 研究方向文件 (Research Direction)

> 最後更新：2026-04-26
> 此文件由 AI 與人類共同維護。AI 在每輪 research-loop 中讀取此文件定向。

> 📅 **重大更新 2026-04-26**：研究主軸從 TO rescue/filter 切換到 Thread D LOH-constrained phasing。
> 詳見 [Thread D 主軸報告](../../docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md)
> 與 [Thread B 撤回宣告](../../docs/reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md)。

---

## 研究主軸（更新 2026-04-26）

**Thread D：LOH-constrained phasing signatures**
- 證據鏈：X3 / X5 / X6 / B1 / B3（4 層）
- Wilcoxon p=0.0156（n=6 exact）
- Evidence grade：B
- 主軸報告：`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`

---

## Candidate Queue（status=pending，不自動執行）

| # | 候選方向 | 動機 | 預估時程 | status |
|---|---------|------|---------|--------|
| 1 | Phase 2B HPFineN marker re-verify (master + flag=on) | R-SELFREF 解除 | 2-4 hr | pending |
| 2 | External CN/SV pilot (Wakhan + SAVANA on HCC1395) | Thread D §3 補強 | 6-8 hr | pending |
| 3 | Caller-level benchmark (DeepVariant / Strelka2 vs ClairS-TO) | HCC1954 caller ceiling 驗證 | 1-2 day | pending |
| 4 | TO archive rerun（HCC1937 / HCC1954 / H2009 / H1437 / DORADO） | KDE-corrected n>6 樣本 | 6-8 hr 背景 | pending |
| 5 | NormalBaseline.cpp writer bug 修復（R-DATA-GAP） | /cpp-change 流程 | 2-4 hr | pending |

## 啟動規則
- **所有 queue 項目需用戶手動批准才啟動**
- 本檔案僅作為候選列表，不寫 execution trigger
- 自動執行 hook 已停用（如有）；下次研究循環需明確指令

---

## 舊研究方向（deprecated 2026-04-26）

> ⚠ **deprecated 2026-04-26 — see Thread B retraction**
> 連結：`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md`

### 原 TO main track（已停用）

**主要 Track**：TO（Tumor-Only）
**理由**：TO 有更多 FP（~300-800），三層架構改善空間比 paired 大，甲基訊號在 TO 下尚有未探索空間。

**主要資料集（pilot）**：HCC1395_5kHz_TO
**次要資料集（medium scale）**：HCC1395_DORADO_TO

**當前基線 F1**：
- HCC1395_5kHz_TO: 0.7127（ClairS-TO baseline）
- HCC1395_DORADO_TO: 0.7226（ClairS-TO baseline）

### 原 TO rescue / TO filter / S3-S5 whitelist 優先項目（已停用）

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
