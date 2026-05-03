---
title: 週報 20260502 - Self-phasing audit chain 重評
type: weekly_master_draft
date: 2026-05-02
status: ready_for_handoff
report_type: problem
main_statement: "上游 938f0df 將 V5 納為 default，audit chain 需重評"
audience: advisor
target_duration_min: 25
period: 2026-04-24 ~ 2026-05-02
source_artifacts:
  - InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md
  - InterSubMod/docs/reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md
  - InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
  - InterSubMod/docs/reports/validated/2026/04/20260429_supplement_getVote_design_intent_QA_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_getVote_only_ablation_results_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260501_latest_longphase_to_3paths_audit_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260501_v5_force_path2only_ablation_01.md
material_classification:
  facts: 11
  observations: 2
  inferences: 1
  unconfirmed: 1
priority_buckets:
  ppt: 8
  speaker_note: 5
  appendix: 1
  shelf: 1
suggested_pptx_template: improvement_report
estimated_pptx_slides: 18
professor_qa_count: 7
handoff_choice: pending
---

# §0 Highlights (TL;DR — 教授前 30 秒能抓到的本週關鍵)

⭐⭐⭐ **This Week's Verdict**: **V5 flag 已被上游 938f0df 採納為 default 行為，flag 本身失效但 audit chain 結論仍成立；下週需驗證跨樣本影響並決策 flag 留存。**

### Top Findings（5 條，重要度 ⭐ 標）

1. ⭐⭐⭐ **[F]** 5/01 三條路 audit 證實 commit `938f0df` 把 highPurity 觸發閾值 `>0.95` → `>0.9` + 改寫 purity 為三次多項式，**V5 flag 變完全 no-op**（PASS variant set / HP tag / LOH.bed / F1 全部與新 baseline 等價）
2. ⭐⭐⭐ **[F]** 4/28 V5 audit 通過 haplotag sanity + paired concordance（**caveat**：4/28 audit 自身已標「不等於 caller F1 改善證據」，5/01 後 V5 audit 結論仍成立）
3. ⭐⭐ **[F]** 4/26 主軸切換 Thread D LOH-constrained phasing；Thread B (whitelist filter) 撤回
4. ⭐⭐ **[F]** 4/29 longphase-TO vs V5 技術報告產出（structured-tech-report 首例），4/30 V3F purity=0.6 + getVote-only ablation × 2 完成
5. ⭐ **[O]** HPFineNGroups marker 在 master 5/7 樣本通過，但 HCC1395 反向（chr8 hotspot 待 deep-dive）

### Top Asks（教授必判斷決策點，3 條）

1. ⭐⭐⭐ **[U] V5 flag 是否 deprecate？** 我傾向「deprecate flag 但保留 audit suite 作上游回滾偵測」，需教授確認
2. ⭐⭐⭐ **[U] 938f0df 跨樣本副作用是否需先評估？** 下週 Priority 1 投入是否值得 4 hr？
3. ⭐⭐ **[U] 4/28-29 audit 文件是否需立即補 caveat banner？** 還是等下週 priority 4 統一處理？

### ⭐⭐⭐ Decisive Next Step

> **Priority 1（必做，4 hr）**：跑 938f0df 全 7 樣本影響評估（純 baseline vs V5 ablation）。**若不做，無法判斷上游變更是否有 regression，後續 audit chain 重評卡住。**

---

# §1 本週主線

⭐⭐⭐ **上游 938f0df 將 V5 納為 default，audit chain 需重評**（30 字）

### Sub-thread（混合主線標註）
- **主線 (problem)**: 上游 binary 變更觸發 V5 flag 失效，需重評 audit chain 有效範圍
- **Sub-thread (progress)**: 4/26-30 V5 audit 系列工作（4/28 audit + 4/29 技術報告 + 4/30 V3F ablation × 2）已完成，未受 5/01 發現否定
- **教授視角優先序**: 主線追問（5/01 發現衝擊）> Sub-thread 進度（已完成 deliverable）

# §2 一句話重點

4/28-29 完成 V5 audit + 技術報告（haplotag sanity 通過、不等於 caller F1 改善），5/01 latest longphase-to commit `938f0df` 將 V5 行為納為 default 後，V5 flag 變 no-op，整條 audit 需重評有效範圍與 V5 留存決策。

---

## Layer 0 研究脈絡

### Layer 0.1 宏觀問題定位

```mermaid
graph LR
  A[4/26 Thread D 主軸切換] --> B[4/28 V5 audit]
  B --> C[4/29 longphase TO vs V5 技術報告]
  C --> D[4/30 V3F ablation × 2]
  D --> E[5/01 三條路 audit]
  E --> F{V5 變 no-op}
  F -->|是 [F]| G[audit chain 重評]
  F -.還需驗證跨樣本.-> H[全 7 樣本 ablation [U]]
```

**核心數字**：
- Self-phasing 17.3:1 (TO baseline 原始) → 1.0:1 (V5 修補後) → 1.0:1 (5/01 新 baseline，無需 V5 flag)
- V5 flag 在 commit 938f0df 後 PASS variant set / HP tag / LOH.bed / F1 全部與新 baseline 等價 [F]
- HPFineNGroups marker master 5/7 通過 [O]

### Layer 0.2 背景知識

- **longphase-to-mod**：tumor-only 模式 phasing tool，commit 938f0df 為 2026-04-30 build
- **V5 (Somatic Fallback)**：2026-04 期間設計的 flag，用 second round phasing 處理 high-purity 樣本
- **highPurity 觸發閾值**：commit 938f0df 從 `> 0.95` 改為 `> 0.9`，並改寫 purity 公式為三次多項式
- KB: `06_workflows/phasing-workflow.md`, `05_tools/longphase-to.md`

### Layer 0.3 上週前情提要

> 上週週報（20260423）確認 NG2 LOH-constrained phasing 訊號 + Thread B (whitelist filter) 撤回 + Thread D 主軸切換待落地。本週實際完成 4/26 主軸切換決策、4/28 V5 audit、4/29 技術報告、4/30 V3F ablation，最終於 5/01 發現上游 binary 變更影響整條 audit 結論。

---

## Layer 1 已建立知識參考

- 已關閉假說：Self-phasing causal chain CONFIRMED（V3-Fixed 已被 V5 取代；SEQC2 F1=0.7153）— Memory: `project_self_phasing_causal_chain_confirmed`
- **本週需更新**：上述 Memory 條目需補註 V5 已被上游 default 化

---

## Layer 2 本週調查

### Thread A: V5 audit 系列工作（4/28-30）

#### 問題陳述
V5 Somatic Fallback flag 是否為當前最佳 phasing baseline？是否影響下游 ISM 分析？

#### 假說與可否證條件
- H1: V5 通過 haplotag sanity → 4/28 驗證
- H2: V5 對下游 caller F1 改善有貢獻 → 4/29 技術報告審查
- 推翻條件：若 V5 等同上游 default 行為 → flag 失效

#### 證據卡

**Tier 1（PPT）**：
- §3 [F] 4/28 V5 audit 通過 haplotag sanity + paired concordance
  Source: `InterSubMod/docs/reports/validated/2026/04/20260428_Self_Phasing_Baseline_V5_Audit_01.md`
  Caveat: V5 audit 結論本身已標「不等於 caller F1 改善證據」
- §3 [F] 4/29 technical report：V5 為當下 phasing baseline（**已加 caveat**：5/01 後此說法需更新為「V5 行為已 default 化」）
  Source: `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`
- §3 [F] 4/30 V3F purity=0.6 ablation 完成，path 2 行為獨立於 V5 flag default 化
  Source: `InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md`

**Tier 2（講稿）**：
- §4 [F] 4/29 supplement getVote design intent QA 補充細節
- §4 [F] 4/30 V3F getVote-only ablation
- §4 [I] 4/27 purity06_simulation_plan（pilot 計畫，仍 viable）

### Thread B: 5/01 重大發現 — 上游 938f0df V5 變 no-op

#### 問題陳述
最新 longphase-to commit 938f0df 是否影響 V5 flag 行為？

#### 方法
跑三條路 audit：Pass 1 only / Pass 2 only / Pass 1 + Pass 2，比對 PASS variant set / HP tag / LOH.bed / F1。

#### 證據卡

**Tier 1（PPT，本週 main thesis 核心）**：
- §3 [F] commit 938f0df 把 highPurity 觸發閾值 `> 0.95` → `> 0.9` + 改寫 purity 為三次多項式
  Source: `InterSubMod/docs/experiments/in_progress/2026/04/20260501_latest_longphase_to_3paths_audit_01.md`#§1
- §3 [F] 三條路 audit 結果：V5 flag 在新 binary 變完全 no-op（PASS variant set / HP tag / LOH.bed / F1 全部與新 baseline 等價）
- §3 [F] 4/26 主軸切換 Thread D 決策已落地（`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`）
- §3 [F] 4/26 Thread B 撤回（whitelist filter use case）
- §4 [O] HPFineNGroups marker 重驗：master 5/7 樣本通過 Fisher odds > 1.5（HCC1395 反向 + HCC1954 amplicon 已知）
- §3 [F] 5/01 v5_force_path2only_ablation 補充：V5 path2 等同新 baseline path2

**Tier 1（PPT，待教授判斷的 [U]）**：
- §5 [U] V5 flag 是否仍需保留？（建議 deprecate flag、保留 audit suite 作上游回滾偵測；最終決策待教授）

#### 因果鏈圖

```mermaid
graph TB
  A[4/26 主軸切換 Thread D] --> B[4/28 V5 audit haplotag sanity 通過]
  B --> C[4/29 技術報告產出: V5 為當前 baseline]
  C --> D[4/30 V3F ablation × 2 完成]
  D --> E[5/01 三條路 audit]
  E --> F[發現 commit 938f0df 將 V5 納為 default]
  F --> G[V5 flag 變 no-op]
  G --> H[audit chain 需重評]
  H --> I[V5 flag 留存決策 [U]]
```

#### 結論
- **判決**: V5 audit 結論在 4/28-29 仍成立，但 5/01 後 V5 flag 已被上游採納為 default，flag 本身失效。
- **穩定度**: HIGH（commit hash + 三條路測試 reproducible）
- **影響**: audit chain 文件需更新；V5 flag 留存決策；HPFineNGroups marker 需在新 baseline 重驗
- **重新開啟條件**: 若上游後續回滾 938f0df，audit suite 仍可作為偵測工具

### §7 報告重點優先順序

1. ★ 5/01 重大發現（V5 變 no-op）
2. 4/26 主軸切換 Thread D
3. 4/28 V5 audit 通過 haplotag sanity（連同 caveat）
4. 4/29 技術報告（連同 5/01 後的 caveat）
5. 4/30 V3F ablation × 2
6. V5 flag 留存決策 [U]
7. HPFineNGroups marker 5/7 通過 [O]
8. Thread B 撤回（背景，不深入）

### §8 建議報告順序

→ 上週背景（Thread D 切換）→ 4/26-30 audit 系列工作 → 5/01 重大發現 → audit chain 重評 → V5 留存決策 → 下週 priority

---

## Layer 3 整合更新

### §11 需要補充的資料

- **R1**: 938f0df commit 對 HCC1395 / HCC1395_DORADO / HCC1937 / HCC1954 / H2009 / H1437 / COLO829 全 7 樣本影響評估（純 baseline vs V5 比對）
- **R2**: HPFineNGroups marker 在新 baseline（含 938f0df）下重跑全 7 樣本
- **R3**: 938f0df commit message 中是否有 deprecation note for V5（待 GitHub PR review）

### §12 需要製作的圖表

- 三條路 audit 結果視覺化（Pass 1 / Pass 2 / 1+2 of PASS variant + HP tag + LOH.bed + F1）
- 7 樣本 V5 flag on/off comparison（在 938f0df binary 上）
- highPurity 閾值與 purity 公式變化圖（0.95 → 0.9 + 三次多項式）

### §13 需要補充的定義或解釋

- highPurity 三次多項式 purity 公式：commit 938f0df 中具體公式（從 `PhasingProcess.cpp:142-220` 抽出）
- second round phasing：V5 originally 設計的目標行為

### §14 可用於講稿的例子

- 「2026-04-22 Self_Phasing_complete_report_for_PI 報告中，原 self-phasing 17.3:1 數字源自 v3 audit；V5 修補後 1.0:1，5/01 新 baseline 也是 1.0:1」
- 「audit chain 像保險：上游正確時 flag 看似多餘，但若上游回滾，audit suite 立即偵測」

---

## Layer 4 未來方向

### §16 下一步行動清單（依 W5 紅旗 2 修正）

1. **Priority 1 (必做，4 hr)**：938f0df 全 7 樣本影響評估（純 baseline vs V5 ablation）
2. **Priority 2 (4 hr)**：HPFineNGroups marker 在新 baseline 重跑全 7 樣本，比對 4/28 結果是否變化
3. **Priority 3 (2 hr)**：audit chain 文件更新（4/28 audit + 4/29 技術報告補 caveat banner）
4. **Priority 4 (1 hr)**：撰寫 V5 flag deprecation decision document，列「保留 audit suite 但 deprecate flag」建議供教授決策
5. **Priority 5 (0.5 hr)**：Memory 更新（`project_self_phasing_causal_chain_confirmed.md` 補註 V5 default 化）

### §17 教授可能提問 + 回答準備（7 個）

| # | 教授追問 | 預備回答 |
|:-:|---------|---------|
| 1 | 938f0df 真的等同 V5 行為？三條路證據強嗎？| 5/01 audit 跑三條路（Pass 1 / Pass 2 / 1+2）→ PASS variant / HP tag / LOH.bed / F1 全部等價。Evidence: `InterSubMod/docs/experiments/in_progress/2026/04/20260501_latest_longphase_to_3paths_audit_01.md` |
| 2 | V5 flag 還要留嗎？ | 建議 deprecate flag 但**保留 audit suite** 作上游回滾偵測。技術可行（reduce maintenance）。**[U] 待教授決策** |
| 3 | 4/28-29 audit 結論還有效嗎？| 4/28 haplotag sanity 仍成立。4/29 技術報告主結論不變，但「V5 = 當前最佳」需改述為「V5 行為已 default 化」 |
| 4 | purity 三次多項式對其他樣本（DORADO / COLO829）有副作用嗎？| **[U] 未測**。下週 Priority 1 |
| 5 | HPFineNGroups marker 重驗結果？| master 5/7 通過 **[O]**；HCC1395 反向 chr8 hotspot 待 deep-dive。下週 Priority 2 |
| 6 | V3F ablation 在 V5 變 no-op 後還有意義嗎？| 仍有意義 — V3F 測 path 2 行為，與 V5 flag default 化獨立 |
| 7 | 下週投入優先序？| 938f0df 全 7 樣本評估（必做）→ HPFineNGroups 重驗 → audit chain 文件更新 |

### 風險評估
- **R1**: 若 938f0df 對某樣本（如 HCC1395_DORADO）造成 regression，需考慮上游回滾或本地 patch
- **R2**: V5 flag deprecation 決策可能影響既有 audit suite 引用，需 cross-check 所有 references

---

## 附錄

### §6 不建議放入 PPT

- 內部工具細節（commit hash 對應的 PhasingProcess.cpp 行號）
- 4/26 主軸切換的完整推導過程（教授已知）
- 5/02 weekly-report skill 升級工作（不屬 self-phasing 主軸，已於 C0 移除）

### §15 暫存紀錄

- HCC1395 與 HCC1395_DORADO 的 platform 差異（下週 Priority 1 結果出後可能引用）
- v5_force_path2only_ablation 與 V3F ablation 的方法學對比（下次方法學文件可深入）

### §9 建議 PPT 模板

**improvement_report**（從問題發現 → 重評 → 行動的改進敘事；雖然主線是 problem，但結尾是 actionable plan，故 improvement template 較適合）

### §10 建議投影片架構（18 張）

1. Cover + thesis（main statement）
2. Background: Thread D 主軸切換（amber motivation strip）
3. Pipeline: V5 audit 系列工作（pipeline_flowchart 4 階段）
4. 4/28 V5 audit 結果（含 caveat）
5. 4/29 longphase TO vs V5 技術報告 highlights
6. 4/30 V3F ablation × 2 結果（before-after-split）
7. ★ **5/01 重大發現**: 938f0df commit + 三條路 audit（核心 slide，TLDR 大字）
8. 三條路 audit 結果表（PASS variant / HP tag / LOH.bed / F1 等價）
9. highPurity 閾值 + purity 公式變化（before-after）
10. V5 flag 留存決策（red/green decision tree with [U] marker）
11. audit chain 重評清單（紅色 caveat strip × 3）
12. HPFineNGroups marker 5/7 通過（含 HCC1395 反向 caveat）
13. Thread B 撤回背景（neutral grey caveat）
14. Risk: 上游回滾 + V5 deprecation 影響範圍
15. **Future: 下週 Priority 1-5**（5card_grid）
16. **Take-home**: V5 flag 已成 no-op，audit suite 是 safety net
17. Q&A backup: 7 個追問 evidence 預備
18. Acknowledgments + references
