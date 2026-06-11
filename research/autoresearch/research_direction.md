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

---
## 轉向記錄 [2026-05-31] — direction_close (queue hygiene)

**類型**: direction_close
**原因**: H013-H018 全部追溯到 concluded-DEAD 方向（methyl filter / caller-F1-headroom / pure-methylation 判別），與 active.json cycle3 P6 NEGATIVE + guardrail 結論一致。
**動作**: H013-H018 status queued → closed（保留 provenance，priority -40），非刪除。
**新焦點**: 不變 — LOH-constrained phasing ⭐3（論文主軸候選）+ ZAR1L/BRCA2 ASM ⭐3（characterization-only）。
**前一輪結果**: Phase2 Cycle3 caller-F1-headroom filter P6_COMMIT NEGATIVE_filter_direction_failed (2026-05-30)。
**注**: 非新 NO-GO 判定（方向早已 concluded）；純 queue ↔ ledger/state 對齊。commit 標 [manual-queue-edit]。

---
## 轉向記錄 [2026-06-11] — focus_consolidation (論文主軸聚焦)

**類型**: focus_consolidation（既有三軸 reframe，非新假設注入）
**原因**: 用戶決定整理放緩所有任務，聚焦單一論文主軸。
**新焦點**: **「Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing」**
  = 既有 G6 三軸的更聚焦 reframe：somatic haplotagging 重建 LOH-constrained 子克隆結構（脊柱，Grade B+ ⭐3）+ 甲基 characterization/誠實負結果（非重建驅動）。
**新主軸 SoT**: `docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`（另一 session 建立，甲基-phasing-assist 面）+ 本 session 補充面 `docs/concepts/2026/06/20260611_Subclonal_Reconstruction_Paper_Focus_整合篇章_01.md`（ASM-characterization + 四道 NEGATIVE + LOH-phasing 脊柱）
**放緩/凍結**: ASM 顯示/工作站/驗證/SEQC2-CN（本 session 完成→凍結為 §Methods-Neg 證據基座）；G1 ASM survey chr5-22；甲基-filter 維持 DEAD。
**決定性卡關**: HD-1（R-SELFREF circularity 跑 or 降 characterization）— GO/NO-GO 由用戶決定。
**queue 對齊**: H019-H023（ASM funnel/CN/LOH）併入 characterization 軸，非新方向；DEAD filter 早已 closed。
**前一輪結果**: 本 session ASM 工作站/驗證/SEQC2-CN 全落地（commit 8dbb931→5e69a99），支撐 characterization+負結果軸。
**注**: paper-scope 聚焦決定（Tier 4），非 NO-GO 判定；commit 標 [manual-queue-edit] 不涉 queue 數值改動。
