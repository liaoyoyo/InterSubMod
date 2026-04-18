<!--
建立時間: 2026-04-05 23:30
目標: Step 3 — PERMANOVA 在 LOH 區域的適用性評估與方法修正
處理範圍: PERMANOVA 前提審計、cnLOH 分析、Allele 維度有效性、方法修正提案
關聯檔案:
  - docs/plans/2026/04/20260405_LOH定義與多維標籤交叉分析研究計劃/00_總覽與執行順序.md
  - src/core/SignificanceAnalyzer.cpp
  - src/core/StructureTest.cpp
-->

# Step 3: PERMANOVA 方法評估計畫

## 背景

PERMANOVA 在 ISM 中用於檢驗 methylation distance matrix 上的群組結構。但在 LOH 區域，HP balance 消失（一組 < 3 reads），前提條件可能不成立。

ISM PERMANOVA 相關程式碼：

- `StructureTest.cpp`: min_reads=5, 999 permutations, 完整距離矩陣
- `SignificanceAnalyzer.cpp:172-249`: 分別對 cluster labels 和 HP/Allele labels 做 PERMANOVA
- `LabelTest.hpp`: LabelTestConfig.min_reads_per_group = 3

Master dataset 欄位：

- `LabelHPPermanovaF/P/Valid` — HP 標籤 PERMANOVA
- `LabelAllelePermanovaF/P/Valid` — Allele 標籤 PERMANOVA
- `ClusterPermanovaF/P/Valid` — Cluster PERMANOVA
- `HPMergedSig`, `AlleleSig`, `HP1FamilyN`, `HP2FamilyN`

---

## 子任務 3.1: PERMANOVA 適用性審計 (P0)

**優先級**: P0（最高優先）
**並行性**: 可與 Step 1 並行執行
**目標**: 確認 LOH 區域中 PERMANOVA 的前提是否成立

### 分析項目

1. **`LabelHPPermanovaValid` rate 比較**：LOH vs non-LOH
   - 預期 LOH 內 valid rate 更低，因為一側 HP 組 reads 不足 min_reads_per_group=3
   - 統計量：valid rate (%), 95% CI, Fisher's exact test p-value

2. **`HP1FamilyN`, `HP2FamilyN` 最小值分佈**（在 LOH 區域）
   - 計算 `min(HP1FamilyN, HP2FamilyN)` 分佈
   - 觀察有多少 LOH 位點的少數組 < 3（PERMANOVA 前提失敗閾值）
   - 與 non-LOH 區域比較

3. **`HPMergedSig` rate 比較**：LOH vs non-LOH
   - HPMergedSig 整合了 PERMANOVA + Chi-squared，觀察整體顯著性差異
   - 分 Paired/TO 觀察是否有模式差異

4. **PERMANOVA pseudo-F 分佈比較**
   - `LabelHPPermanovaF` 在 LOH 內 vs non-LOH 的分佈
   - 條件：僅取 `LabelHPPermanovaValid=True` 的位點
   - 觀察 LOH 內 pseudo-F 是否系統性偏高或偏低

5. **分 Paired/TO 模式觀察**
   - 所有上述分析均需分 Paired 與 TO 兩組
   - Paired 模式有真實 phasing，TO 模式有 self-phasing 偏差

6. **分 TP/FP 觀察**
   - 在 LOH 與 non-LOH 各子群中，TP 與 FP 的 PERMANOVA 結果是否有差異
   - 計算 AUC（以 `LabelHPPermanovaP` 區分 TP/FP）

### 產出

- **腳本**: `build_permanova_loh_audit.py`
- **輸出目錄**: `{OUTPUT_ROOT}/20260408_permanova_loh_audit/`
- **圖表**: 8 張
  1. LOH vs non-LOH 的 `LabelHPPermanovaValid` rate（分 Paired/TO，bar chart）
  2. `min(HP1FamilyN, HP2FamilyN)` 分佈（LOH vs non-LOH，histogram）
  3. `HPMergedSig` rate 比較（LOH vs non-LOH，分 Paired/TO，bar chart）
  4. `LabelHPPermanovaF` 分佈（LOH vs non-LOH，violin plot）
  5. `LabelHPPermanovaP` 分佈（LOH vs non-LOH，CDF plot）
  6. TP/FP 在 LOH 內的 pseudo-F 比較（box plot）
  7. 7 樣本逐一的 valid rate 熱圖（sample x LOH_status，heatmap）
  8. 綜合摘要表（validity / sig rate / pseudo-F median，table figure）

### 判定標準

| 指標 | 閾值 | 解讀 |
|------|------|------|
| LOH 內 HP PERMANOVA valid rate | < 50% | HP 維度在 LOH 內大部分不可用 |
| LOH 內 min(HP1N, HP2N) < 3 比例 | > 60% | 群組大小前提普遍不成立 |
| LOH 內 HPMergedSig AUC (TP vs FP) | < 0.55 | HP 顯著性在 LOH 內無區分力 |

---

## 子任務 3.2: cnLOH 分析 (P1)

**優先級**: P1
**依賴**: 子任務 1.2（Coverage_Multiple 定義）+ 子任務 3.1
**目標**: 識別 Copy-neutral LOH (cnLOH) 並分析其特殊性質

### 生物學背景

cnLOH（Copy-neutral LOH）保留兩份同源拷貝，但都來自同一親本染色體。其結果是：

- **覆蓋度正常**（不同於缺失型 LOH 的半量覆蓋）
- **失去 heterozygosity**（germline SNV 變為 homozygous）
- **甲基化模式**：兩份拷貝的甲基化可能一致（同源），減少 epigenetic heterogeneity

### cnLOH 定義

```
cnLOH := HP_Imbalance=True AND Coverage_Multiple ∈ [0.8, 1.2]
```

其中 Coverage_Multiple 由 Step 1.2 定義（區域覆蓋 / 預期覆蓋）。

### 分析項目

1. **cnLOH 比例統計**
   - cnLOH 在 TO 中佔 HP Imbalance 的多少比例
   - 分樣本統計，與文獻參考值比較（實體腫瘤 10-30% LOH 為 cnLOH）

2. **cnLOH 的 TP/FP 分佈**
   - cnLOH vs 非 cnLOH（缺失型 LOH）的 TP/FP 比例比較
   - 觀察 cnLOH 是否有更高的 FP enrichment

3. **cnLOH 內甲基化區分能力**
   - `PairwiseMedianDist` 在 cnLOH 內的分佈（預期更低，因為兩份同源拷貝甲基化一致）
   - `AlleleDelta` 在 cnLOH 內的分佈
   - 與非 cnLOH LOH 區域比較

4. **cnLOH 與 LOH.bed 交集**
   - cnLOH 位點是否落在 LOH.bed 內
   - 兩種定義的一致性（HP Imbalance + normal coverage vs LOH.bed）

### 產出

- **腳本**: `build_cnloh_analysis.py`
- **輸出目錄**: `{OUTPUT_ROOT}/20260408_cnloh_analysis/`
- **圖表**: 6 張
  1. cnLOH 比例（分樣本 + 文獻參考帶，bar chart）
  2. cnLOH vs deletion-LOH 的 TP/FP 比例（stacked bar）
  3. `PairwiseMedianDist` 分佈（cnLOH vs deletion-LOH vs non-LOH，violin）
  4. `AlleleDelta` 分佈（cnLOH vs deletion-LOH vs non-LOH，violin）
  5. cnLOH 與 LOH.bed 交集 Venn 圖
  6. Coverage_Multiple 分佈（cnLOH vs deletion-LOH，histogram + 閾值線）

### 判定標準

| 指標 | 閾值 | 解讀 |
|------|------|------|
| cnLOH 佔 HP Imbalance 比例 | 10-30% | 符合文獻預期範圍 |
| cnLOH 內 PairwiseMedianDist | < non-LOH | 同源拷貝甲基化一致性高 |
| cnLOH 內 AlleleDelta AUC | < 0.55 | ALT/REF 維度在 cnLOH 也失效 |

---

## 子任務 3.3: LOH 內 Allele 維度分析 (P1)

**優先級**: P1
**依賴**: 子任務 3.1
**目標**: 當 HP 維度在 LOH 區域失效時，ALT/REF（Allele）維度是否仍保有區分能力

### 核心邏輯

ALT/REF 標籤來自 SNV 位點的鹼基（reference allele 或 alternative allele），不依賴 phasing 結果。即使 HP balance 消失導致 HP-based PERMANOVA 失效，ALT/REF 仍然可以作為群組維度。

關鍵差異：

| 維度 | 來源 | 受 phasing 影響 | 受 LOH 影響 |
|------|------|----------------|-------------|
| HP (1/2) | LongPhase haplotagging | 是 | 是（一側消失） |
| Allele (ALT/REF) | SNV 位點鹼基 | 否 | 部分（ALT 少但仍存在） |

### 分析項目

1. **LOH 內 AlleleDelta AUC（TP vs FP）**
   - 計算 LOH 區域中 `AlleleDelta` 區分 TP/FP 的 AUC
   - 與 non-LOH 區域的 AUC 比較
   - 分 Paired/TO 觀察

2. **LOH 內 `LabelAllelePermanovaP` 有效性**
   - `LabelAllelePermanovaValid` rate 在 LOH 內是否高於 HP 版本
   - `LabelAllelePermanovaP` < 0.05 的比例（顯著性 rate）
   - Allele PERMANOVA pseudo-F 分佈

3. **Allele 有效性條件判斷**
   - 如果 Allele 維度在 LOH 內仍有效 → LOH 區域應改用 Allele-only 分析
   - 如果 Allele 維度也失效 → LOH 區域需要全新方法

4. **與非 LOH 區域的 AUC 交叉比較**
   - 2x2 矩陣：(LOH/non-LOH) x (HP/Allele) 的 AUC
   - 找出最佳的 (區域類型, 維度) 組合

### 產出

- **腳本**: 整合於 `build_permanova_loh_audit.py` 或獨立 `build_allele_loh_analysis.py`
- **輸出目錄**: `{OUTPUT_ROOT}/20260408_allele_loh_analysis/`
- **圖表**: 6 張
  1. AlleleDelta AUC 比較（LOH vs non-LOH，分 Paired/TO，grouped bar）
  2. Allele vs HP PERMANOVA valid rate（LOH 內，paired bar）
  3. `LabelAllelePermanovaF` 分佈（LOH vs non-LOH，violin）
  4. 2x2 AUC 矩陣熱圖（(LOH/non-LOH) x (HP/Allele)，heatmap）
  5. LOH 內 AlleleDelta 分佈（TP vs FP，histogram + KDE）
  6. AlleleSig rate 比較（LOH vs non-LOH，分 Paired/TO，bar chart）

### 判定標準

| 指標 | 閾值 | 解讀 |
|------|------|------|
| LOH 內 Allele PERMANOVA valid rate | > HP valid rate | Allele 維度在 LOH 內更可靠 |
| LOH 內 AlleleDelta AUC | > 0.58 | Allele 維度在 LOH 內仍有區分力 |
| 2x2 AUC 最佳組合 | LOH+Allele > LOH+HP | 確認 LOH 應切換至 Allele-only |

### 決策樹

```
LOH 內 Allele AUC > 0.58?
├─ 是 → LOH 區域改用 Allele-only 分析（進入 3.4 修正提案）
└─ 否 → LOH 內 Allele 也失效
         ├─ cnLOH 佔多數? → LOH 區域標記為 "low-confidence"，不做 PERMANOVA
         └─ deletion-LOH 佔多數? → 考慮 coverage-based 方法替代
```

---

## 子任務 3.4: 方法修正提案 (P2)

**優先級**: P2
**依賴**: 子任務 3.1 + 3.2 + 3.3（需要前三項結論才能制定修正）
**目標**: 根據 3.1-3.3 的發現，提出具體的 ISM 方法修正方案

### 結構框架

每個修正提案遵循以下結構：

```
現象 → 原因 → 修正方式 → 驗證方法 → 結果解釋
```

### 可能修正方向

#### 方向 A: LOH 區域跳過 HP-based 檢定，僅保留 Allele 維度

- **觸發條件**: 3.3 確認 Allele 維度在 LOH 內有效（AUC > 0.58）
- **現象**: LOH 區域 HP balance 消失，HP PERMANOVA valid rate < 50%
- **原因**: LOH 導致一側 HP 組 reads 不足 min_reads_per_group
- **修正方式**: 在 LOH 區域，SignificanceAnalyzer 跳過 HP-based PERMANOVA，改用 Allele-based 結果
- **程式碼位置**: `SignificanceAnalyzer.cpp:269-304`（Multi-stage HP verification）
- **驗證方法**: HCC1395 TO 數據前後比較，確認 TP loss <= 2%
- **結果解釋**: LOH 區域的 significance 改用 AlleleSig 而非 HPMergedSig

#### 方向 B: Coverage-stratified 分析策略

- **觸發條件**: 3.2 確認 cnLOH 與 deletion-LOH 有不同表現
- **現象**: LOH 區域不是同質的，cnLOH 和 deletion-LOH 行為不同
- **原因**: cnLOH 覆蓋正常但同源，deletion-LOH 覆蓋減半
- **修正方式**: 根據 Coverage_Multiple 區分 LOH 類型，採用不同分析策略
  - cnLOH (Coverage_Multiple 0.8-1.2): 嘗試 Allele-only 或降低 min_reads_per_group
  - deletion-LOH (Coverage_Multiple < 0.8): 標記 low-confidence 或跳過 PERMANOVA
- **程式碼位置**: `RegionProcessor.cpp`（需新增 LOH type 判斷邏輯）
- **驗證方法**: 分 LOH type 觀察修正後的 AUC 變化

#### 方向 C: QS 權重調整

- **觸發條件**: 3.1 確認 LOH 區域的 PERMANOVA 結果不可靠
- **現象**: 目前 QS 計算對 LOH 區域施加 penalty，但 penalty 可能反向（見 memory: QS TO 失效證據）
- **原因**: LOH penalty 基於 HP balance 假設，在 TO self-phasing 下反向
- **修正方式**: 根據 LOH type 差異化 QS 權重
  - HP-Imbalance-only: 降低 HP 權重，提高 Allele 權重
  - LOH.bed overlap: 根據 LOH type 決定是否 penalty
  - 兩者皆有: 最保守策略
- **程式碼位置**: `RegionProcessor.cpp:113-128`（QS 計算邏輯）
- **驗證方法**: QS AUC 前後比較（目標：TO QS AUC > 0.55，目前 0.497）

#### 方向 D: ISM C++ 端新增 Allele-primary mode

- **觸發條件**: 3.3 確認 Allele 維度全面優於 HP 維度（在 LOH 內）
- **現象**: ISM 目前以 HP 為主要維度，Allele 為輔助
- **原因**: 設計假設所有區域都有 HP balance
- **修正方式**: 在 LOH 偵測到時自動切換至 Allele-primary mode
  - `SignificanceAnalyzer` 新增 mode 參數
  - LOH flag 由 RegionProcessor 傳入
  - Allele-primary mode 下 HPMergedSig 不計入最終判定
- **程式碼位置**: `SignificanceAnalyzer.cpp` + `RegionProcessor.cpp`
- **驗證方法**: 端到端測試，確認 LOH 區域自動切換且結果改善

### 修正方向優先序

```
3.1-3.3 結論
    │
    ├─ Allele 有效 (AUC > 0.58) → 方向 A (快速) → 方向 D (完整)
    │
    ├─ cnLOH/deletion 有差異 → 方向 B
    │
    └─ 全面失效 → 方向 C (QS 調整為 fallback)
```

### 產出

- **修正提案文件**: `{OUTPUT_ROOT}/20260408_permanova_method_revision_proposal.md`
- **包含**: 每個方向的完整 現象→原因→修正→驗證→解釋 鏈
- **包含**: 修正前後預期效果估算
- **包含**: 程式碼修改 checklist 與測試計畫

---

## 驗證清單

完成所有子任務後，確認以下關鍵預期：

- [ ] `LabelHPPermanovaValid` 在 LOH 內比 non-LOH 低（方向驗證）
- [ ] cnLOH 比例合理（文獻參考：實體腫瘤 10-30% LOH 為 cnLOH）
- [ ] `AlleleDelta` 在 LOH 內仍可計算（ALT/REF 不受 phasing 影響）
- [ ] 方法修正前後比較：HCC1395 TO 數據，TP loss <= 2%

---

## 時程與依賴關係

```
Step 1.2 (Coverage_Multiple)
    │
    ▼
  3.1 (PERMANOVA 審計) ──────────────────┐
    │                                     │
    ├──► 3.3 (Allele 維度分析) ──────────┤
    │                                     │
    └──► 3.2 (cnLOH 分析，需 1.2) ──────┤
                                          │
                                          ▼
                                    3.4 (方法修正提案)
```

| 子任務 | 預計工時 | 前置依賴 | 可並行項目 |
|--------|---------|---------|-----------|
| 3.1 PERMANOVA 審計 | 3h | 無 | Step 1 全部 |
| 3.2 cnLOH 分析 | 2h | 1.2 + 3.1 | — |
| 3.3 Allele 維度分析 | 2h | 3.1 | 3.2 |
| 3.4 方法修正提案 | 3h | 3.1 + 3.2 + 3.3 | — |
| **合計** | **10h** | | |

---

## 風險與備註

1. **cnLOH 定義敏感度**: Coverage_Multiple 閾值 [0.8, 1.2] 可能需要根據實際分佈微調
2. **Allele 維度在低 AF 位點的可靠性**: AF < 0.1 時 ALT reads 極少，Allele PERMANOVA 也可能失效
3. **方向 D 的工程量**: 新增 Allele-primary mode 涉及多個檔案修改，需要完整的迴歸測試
4. **與 QS 重設計的關係**: 方向 C 與現有 QS TO 失效問題（memory: QS AUC=0.497）高度相關，修正需統一考量
