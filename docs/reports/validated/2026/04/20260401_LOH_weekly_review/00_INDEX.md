<!--
建立時間: 2026-04-01 03:30
目標: LOH 研究週報審閱資料夾導覽索引
處理範圍: 2026-03-25 ~ 2026-04-01 所有 LOH 相關研究
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
-->

# LOH 研究週報審閱資料夾 — 導覽索引

**週期**：2026-03-25 ~ 2026-04-01
**數據基準**：all_region_rows.tsv.gz (748,391 rows × 116 cols, post HP-fix)
**產出**：10 份主題文件 + 3 份審查報告 + 1 份修正紀錄 = 15 份文件，~299 KB
**修正狀態**：v2（已根據用戶審閱回饋修正 HP tag 語義、補充 LOH-like/Region/Tier 定義、新增 Cell Line 限制段落）

---

## 閱讀順序建議

```
00_background.md  ← 先讀（術語與背景）
       ↓
01_hp_integer_tag_fix.md  ← 本週最重要的基礎修正
       ↓
02_loh_evidence_panel_rounds1_4.md  ← LOH 四輪系統研究
       ↓
03_post_hp_fix_loh_enrichment.md  ← HP fix 後 LOH 翻轉確認
       ↓
04_mechanism_investigation.md  ← 核心機制解釋（為什麼翻轉）
       ↓
05_systematic_observation_O1_O10.md  ← 82 張圖的系統性觀察摘要
       ↓
06_methylation_hypothesis_negative.md  ← 兩個假說被否決（重要 negative results）
       ↓
07_qs_mode_aware_change.md  ← 基於觀察的程式碼修改
       ↓
08_literature_survey.md  ← 文獻理論支撐
       ↓
09_conclusions_and_actions.md  ← 總結與下週行動
```

---

## 文件清單

### 主題文件（按邏輯順序）

| # | 檔案 | 頁數 | 內容摘要 |
|---|------|------|---------|
| 0 | [00_background.md](00_background.md) | ~260 行 | ISM 工具說明、7 樣本、Paired/TO 定義、術語表 |
| 1 | [01_hp_integer_tag_fix.md](01_hp_integer_tag_fix.md) | ~354 行 | HP:i:11/21/33 bug 發現、修正、影響量化、文件分級 |
| 2 | [02_loh_evidence_panel_rounds1_4.md](02_loh_evidence_panel_rounds1_4.md) | ~500 行 | 四輪 LOH 研究：基線→分層→HP0 filter→深度驗證 |
| 3 | [03_post_hp_fix_loh_enrichment.md](03_post_hp_fix_loh_enrichment.md) | ~287 行 | HP fix 後 TO LOH enrichment 翻轉（0.805× TP-enriched） |
| 4 | [04_mechanism_investigation.md](04_mechanism_investigation.md) | ~274 行 | 同位點 16-52× 過判、HP imbalance 機制、TP rescue 候選 |
| 5 | [05_systematic_observation_O1_O10.md](05_systematic_observation_O1_O10.md) | ~333 行 | 9 觀察、82 圖表、Top 10 發現、行動建議 |
| 6 | [06_methylation_hypothesis_negative.md](06_methylation_hypothesis_negative.md) | ~376 行 | O11 heterogeneity confound、O12 三場景不可區分、L2 collider bias |
| 7 | [07_qs_mode_aware_change.md](07_qs_mode_aware_change.md) | ~171 行 | TO QS LOH penalty=0、verify_bonus=0 修改 |
| 8 | [08_literature_survey.md](08_literature_survey.md) | ~191 行 | 70+ 篇文獻、三場景理論、mQTL、現有工具比較 |
| 9 | [09_conclusions_and_actions.md](09_conclusions_and_actions.md) | ~148 行 | 10 確認結論、4 否決假說、3 待定、P0/P1/P2 行動 |

### 品質審查報告

| 審查面向 | 檔案 | 發現摘要 |
|---------|------|---------|
| 敘述完整性 | [review_completeness.md](review_completeness.md) | 5/10 五環全滿；07(QS)缺 post-implementation 驗證；08(文獻)缺檢索方法論 |
| 邏輯質疑 | [review_logic.md](review_logic.md) | 20 項質疑：4 高嚴重度、11 中、5 低。核心機制僅有間接證據 |
| 認知門檻 | [review_accessibility.md](review_accessibility.md) | Top 10 補充概念；06(collider bias)認知難度最高；建議 10 張示意圖 |

---

## 本週關鍵結論速覽

### 確認的結論（10 項）

| # | 結論 | 關鍵數據 | 來源 Section |
|---|------|---------|-------------|
| 1 | HP integer tag bug 修正完成 | 88% TO regions → Tier A/A+ | §1 |
| 2 | TO LOH enrichment = TP-enriched | 0.805×, 7/7 樣本一致 | §3 |
| 3 | TO 系統性過判 LOH | 同位點 16-52× | §4 |
| 4 | TO 無單一有效特徵 | 所有 AUC < 0.58 | §5(O1,O4,O5,O8) |
| 5 | LOH penalty 是 TO QS 失效根因 | 移除後 +0.045 AUC | §5(O2,O3) |
| 6 | Paired/TO 必須分離模型 | FP rate 1% vs 31% | §5(O1,O5,O7) |
| 7 | GQ 是 paired 最強特徵 | AUC=0.811 | §5(O4) |
| 8 | 樣本異質性 8.6× | TO FP rate 8.7%-74.6% | §5(O8) |
| 9 | VerificationClass 無效 | Cramér's V=0.023 | §5(O6) |
| 10 | QS Mode-Aware 修改已完成 | commit b9eaba7 | §7 |

### 否決的假說（4 項）

| # | 假說 | 否決原因 | 來源 |
|---|------|---------|------|
| 1 | Within-group heterogeneity 可區分 TP/FP | n_reads confound (AUC 0.845→0.530) | §6(O11) |
| 2 | TO LOH 三場景可透過甲基化區分 | AF confound + L2 collider bias | §6(O12) |
| 3 | HP0 高 = LOH 品質差 | 方向相反 (High HP0 TP% > Low) | §2(R3) |
| 4 | LOH 可作為 binary FP filter | 所有 F1 delta < 0 | §2(R3) |

### 待定問題（3 項）

| # | 問題 | 需要什麼 | 優先級 |
|---|------|---------|--------|
| 1 | FN 位點的 LOH 特性 | FN ISM 數據（目前無） | P2 |
| 2 | HCC1395 chr8 LOH+HPSig 特異性 | ASM+LOH block 深度分析 | P3 |
| 3 | Low AF TP LOH rescue 可行性 | FN-level 數據 | P2 |

---

## 審查 agents 發現的高優先問題

### 邏輯質疑（高嚴重度，4 項）

1. **「TO somatic allele 造成 HP imbalance」僅有間接證據**，缺乏 per-read 級別直接驗證
2. **Epipolymorphism 殘差化方法**——n_reads 可能是 mediator 而非 confound，殘差化可能過度移除信號
3. **HP_Ratio Pearson r=0.001 與 concordance 85%+ 表述矛盾**——**已解決**：散點圖呈 "X" 形交叉，haplotype swap 導致正負相關抵消（Pearson r≈0），但 LOH binary 分類不受 swap 影響（concordance 85.5%）。詳見 §5 O7 修正段落
4. **7 個細胞系的 external validity 風險**——結論推廣到臨床樣本需謹慎

### 認知門檻（Top 5 最需補充概念）

1. **Collider Bias**（OLS 殘差化中的虛假信號）— §6 最難理解的段落
2. **Residualization**（控制混淆因素的回歸方法）
3. **Phasing 與 Haplotype Assignment 原理**
4. **AUC < 0.5 的意義**（方向反轉而非無效）
5. **Simpson's Paradox**

---

## 建議閱讀方式

1. **快速掃描（30 分鐘）**：讀本 INDEX → §9(結論) → §0(背景)
2. **重點審閱（2 小時）**：依序讀 §0→§1→§3→§4→§9，聚焦 HP fix 與機制
3. **完整審閱（半天）**：全 10 份 + 3 份審查報告
4. **驗證模式**：每個 Section 末尾的「待驗證問題」逐項確認
