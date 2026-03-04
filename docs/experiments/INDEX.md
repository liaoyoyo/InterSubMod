<!--
建立時間: 2026-03-05 10:00
目標: 實驗研究歷史主索引，聚合所有已試驗方向，供 AI 快速掌握研究脈絡
處理範圍: InterSubMod 專案 2025-11 至今的所有研究方向
關聯檔案:
  - docs/README.md
  - docs/CURRENT_FOCUS.md
  - docs/experiments/validated/
  - docs/experiments/finalized/
  - docs/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md
-->

# InterSubMod 實驗研究索引

> **AI 使用提示**：本索引為 Layer 2（歷史全貌）。從這裡了解哪些方向已探索、結論為何，再決定下一步。

## 狀態圖例

- ✅ 成功完成，結論明確
- ❌ 失敗或無效，不建議重試
- ⏳ 進行中，尚無定論
- 🔄 值得再探索，有改進空間

---

## Layer 1：研究方向總覽

### 01. 甲基化解析（MM/ML 標籤）

- **期間**: 2025-11
- **狀態**: ✅ 成功
- **關鍵結論**: 支援 ONT BAM 的 `MM`/`ML` 標籤解碼，含正反股 CIGAR 座標校正，CpG 定位正確
- **主要文件**: `docs/solutions/debugging/2025/11/20251130_methylation_parser_debug_01.md`
- **建議後續**: 補齊 MethylationParser 獨立單元測試（目前僅有整合測試）

### 02. CIGAR 座標映射

- **期間**: 2025-11
- **狀態**: ✅ 成功
- **關鍵結論**: 反向股需從 SEQ 末端往回迭代匹配 MM delta encoding；`ref[offset-1]='C' && ref[offset]='G'` 驗證邏輯正確但非顯而易見
- **主要文件**: `docs/solutions/debugging/2025/11/20251130_mm_ml_tags_analysis_01.md`
- **建議後續**: 建立 mock BAM record 測試 forward/reverse 各種 CIGAR 邊界情境

### 03. 距離計算與聚類分析

- **期間**: 2025-11 ~ 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: 實作 NHD / L1 / L2 / CORR / BERNOULLI / JACCARD 六種距離度量；UPGMA 為預設聚類方法
- **主要文件**: `docs/archive/2025/12/20251218_final_completion_report_01.md`
- **建議後續**: 考慮 nearest-neighbor chain 加速聚類（n>500 時 O(n³) 可能成瓶頸）

### 04. OpenMP 平行化

- **期間**: 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: Per-region 平行化，thread-local BamReader/FastaReader 避免鎖競爭；單 Region 平均 < 300ms
- **主要文件**: `docs/architecture/system_overview.md`
- **建議後續**: 相鄰 SNV window 重疊 reads 可考慮共享（可減少 30-50% BAM I/O）

### 05. 統計顯著性：Fisher / χ²

- **期間**: 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: Fisher-Freeman-Halton (RxC) 使用 Patefield 算法，Monte Carlo early stopping 有效節省計算；目前主流程的核心統計方法
- **主要文件**: `docs/archive/2025/12/20251202_phase3_report_01.md`
- **建議後續**: 加入多重檢驗校正（FDR/Benjamini-Hochberg），30K+ region 獨立測試存在 false discovery 風險

### 06. 統計顯著性：PERMANOVA

- **期間**: 2025-12 ~ 2026-01
- **狀態**: ✅ 已實作，⚠️ 預設關閉
- **關鍵結論**: 已實作 pseudo-F = SS_Between/SS_Within；目前主流程 `enable_permanova=false`（效能考量）；啟用時需完整距離矩陣（無 NaN）
- **主要文件**: `docs/experiments/validated/2025/12/20251230_bidirectional_verification_strategy_analysis_01.md`
- **建議後續**: 在 Label-First 驗證框架下（HP/Allele 分組）重新評估 PERMANOVA 的適用時機

### 07. Bernoulli 距離度量

- **期間**: 2025-12
- **狀態**: ✅ 成功
- **關鍵結論**: `w(p)=2×|p-0.5|`（信心權重），高信心位點才對距離有大貢獻；比 NHD 更適合 ONT 中間不確定值
- **主要文件**: `docs/plans/2025/12/20251222_Bernoulli_Implementation_Plan_01.md`
- **建議後續**: 在 purity 系列測試中比較 Bernoulli vs NHD 的效能差異

### 08. TP/FP 特徵富集分析

- **期間**: 2026-01
- **狀態**: ✅ 成功（見 Layer 2 深入摘要）
- **關鍵結論**: LOH_like 富集 TP（OR=7.33），HighCluster 富集 FP（OR=0.02），LowAF 富集 FP（OR=0.05）
- **主要文件**: `docs/reports/validated/2026/03/20260301_甲基化欄位對TPFP與subclone驗證比較_01.md`
- **建議後續**: 對全量 TP/FP（而非 448 筆抽樣）重算 OR/F1，補充統計穩健性

### 09. 雙向驗證（Label-First vs Cluster-First）

- **期間**: 2026-01
- **狀態**: ✅ 框架設計完成（見 Layer 2 深入摘要）
- **關鍵結論**: 三角互證框架：PERMANOVA（Label→Structure）+ Fisher（Cluster→Label）；4 類輸出：Strong / Subclone / Weak / Noise
- **主要文件**: `docs/experiments/validated/2025/12/20251230_bidirectional_verification_strategy_analysis_01.md`
- **建議後續**: 補齊 Bootstrap stability check 實作；LabelTest 多階段 HP 驗證需獨立單元測試

### 10. F1 最佳化

- **期間**: 2026-01
- **狀態**: ✅ 成功（見 Layer 2 深入摘要）
- **關鍵結論**: 最優策略 `AF≥0.3 OR VerificationClass=Subclone`，F1=0.8481（基線 0.8155）；LabelDelta>0.3 可微幅提升至 0.8158
- **主要文件**: `docs/experiments/validated/2026/01/20260107_F1_Optimization_Deep_Analysis_01.md`
- **建議後續**: 在 COLO829 + purity series 做交叉驗證；建立 ML 分類器整合多特徵

### 11. Subsample 混樣甲基化偏差

- **期間**: 2026-02 ~ 2026-03
- **狀態**: ⏳ 進行中（見 Layer 2 深入摘要）
- **關鍵結論**: HCC1395 ONT subsample 無 MM/ML，無法做甲基化分析；已在腳本加入 Methylation Guard；根本問題為 tumor/blood-normal 組織來源差異
- **主要文件**: `docs/ai_sessions/2026/03/20260302_subsample混樣甲基化偏差_現況研究推論與驗證路線圖_01.md`
- **建議後續**: 使用 HCC1395 DORADO（有 MM/ML）做 tumor-only / normal-only / mixed 三路對照

### 12. Purity-Aware 過濾策略

- **期間**: 2026-02 ~ 2026-03
- **狀態**: 🔄 值得再探索（見 Layer 2 深入摘要）
- **關鍵結論**: 低 purity（<40%）時固定甲基化門檻嚴重傷害 Recall；需依 purity 分層的自適應策略
- **主要文件**: `docs/plans/2026/02/20260228_InterSubMod再驗證與再實驗執行計畫_01.md`
- **建議後續**: 低 purity 停用甲基化過濾，中/高 purity 保留；建立 purity-aware feature 融合模型

---

## Layer 2：重點方向深入摘要

> 以下 5 個方向為目前最具研究價值或待解問題，提供完整數據與失敗嘗試紀錄。

---

### 深入：TP/FP 特徵富集分析

**研究問題**: InterSubMod 的甲基化特徵是否能有效區分 Somatic SNV 的真陽性 (TP) 與假陽性 (FP)？

**方法**:
- 資料集：448 筆分層抽樣（TP 288、FP 160），HCC1395 + SEQC2 truth set
- 定義特徵：LOH_like, HighCluster（ClusterCount50kb>20）, LowAF（AF<0.3）, MethylSupport_SubcloneLike
- 統計：Fisher exact test（OR + p-value）+ 規則式 F1 評估

**數據結果**:
- `LOH_like` 富集 TP：OR=7.33, p=3.04e-16，TP rate=74.2%（非條件=28.1%）
- `HighCluster` 富集 FP：OR=0.0185, p=5.75e-26，TP rate=4.9%（非條件=73.6%）
- `LowAF` 富集 FP：OR=0.0458, p=4.44e-40，TP rate=34.1%（非條件=91.9%）
- 低 AF 子群下，`VerificationClass=Subclone`：OR=8.13, p=9.88e-05（顯著 TP 富集）

**最佳 F1 策略**:
| 策略 | F1 |
|---|---|
| 基線（全部保留）| 0.8155 |
| `AF≥0.3` | 0.8238 |
| `AF≥0.3 OR VerificationClass=Subclone` | **0.8481**（最優） |
| `SignificantOnly` | 0.1910（無效）|

**失敗嘗試**: 單獨使用 CramersV 過濾（AUC=0.52）無效；SignificantOnly 策略嚴重損失 Recall

**關鍵限制**: 核心比較集為分層抽樣（448筆），非全量 30K+ 位點；全量結論需再驗證

**建議下一步**:
1. 對全量 TP/FP 重算 LOH_like / HighCluster 的 OR 穩健性
2. 建立 logistic / tree 模型，組合 `VerificationClass + HPMergedDelta + ClusterCount50kb`
3. 在 COLO829 做 held-out 交叉驗證

**相關文件**:
- `docs/reports/validated/2026/03/20260301_甲基化欄位對TPFP與subclone驗證比較_01.md`
- `docs/experiments/validated/2026/01/20260107_F1_Optimization_Deep_Analysis_01.md`
- `docs/experiments/validated/2026/01/20260107_F1_and_Data_Optimization_Report_01.md`

---

### 深入：Purity-Aware 過濾策略

**研究問題**: 如何設計甲基化過濾策略，使其在不同腫瘤純度（purity）下都能保持穩定效能？

**方法**:
- 資料集：HCC1395 DORADO subsample（t7_n29 ~ t50_n00，有 MM/ML）
- 比較固定門檻 vs purity 分層策略的 F1 差異

**數據結果**:
- 低 purity（~19.4%，t7_n29）：固定甲基化門檻顯著傷害 Recall（F1 下降 0.05 級）
- 中高 purity（≥60%）：過濾效果趨近中性或微幅正向
- 結論：固定單一門檻跨 purity 直接套用不可行

**建議策略**:
```
< 40% purity  → 不做甲基化過濾
40-60% purity → 保守門檻或僅標註不刪除
≥ 60% purity  → 啟用甲基化過濾
```

**失敗嘗試**: 固定 `Significant=True` 在低 purity 下等同於強制移除大量 TP

**關鍵限制**: 目前缺乏跨 purity + 跨樣本的完整驗證矩陣；COLO829 與其他樣本尚未納入

**建議下一步**:
1. 停用低 purity 甲基化過濾，重跑 F1 曲線確認恢復
2. 建立統一 `purity_qc.tsv`：每個 subsample 記錄 MM_count, ML_count, analyzable
3. 規劃 Experiment-02（多特徵融合）作為固定門檻的替代方案

**相關文件**:
- `docs/plans/2026/02/20260228_InterSubMod再驗證與再實驗執行計畫_01.md`
- `docs/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md`

---

### 深入：Subsample 混樣甲基化偏差

**研究問題**: 為何 HCC1395 subsample 分析效果不佳？是 tumor/normal 組織來源差異（甲基化背景）造成的系統性偏差？

**方法**:
- 比較 HCC1395 ONT（無 MM/ML）與 HCC1395 DORADO（有 MM/ML）兩批資料
- 分析 tumor（乳腺癌細胞株）vs blood-normal（HCC1395BL，B-lymphoblast）的甲基化背景差異

**數據結果**:
- HCC1395 ONT subsample：無 MM/ML → `tp_regions=0, fp_regions=0`，InterSubMod 完全無有效輸出
- HCC1395 DORADO：有 MM/ML，但低 purity 時固定門檻傷害 Recall
- 外部文獻（Cell Reports Methods 2025）確認：HCC1395 vs HCC1395BL 甲基化 overlap 約五成，組織差異顯著

**已採取措施（2026-03-02）**:
- `scripts/pipeline/run_benchmark.sh`：新增 Methylation Guard（MM/ML 前置檢查）
- `scripts/analysis/run_purity_and_standard_verification.sh`：新增 source BAM MM/ML 檢查
- `scripts/pipeline/config.sh`：新增共用函式 `has_mm_ml_tags()`

**推論（依可信度）**:
1. H1（高）：MM/ML 缺失是第一層阻斷，直接造成無有效統計區域
2. H2（高）：低 purity 下固定門檻將大量 TP 錯誤過濾
3. H3（中高）：tumor/blood-normal 甲基化背景差異放大 composition effect
4. H4（中高）：InterSubMod 主程式未實際把 normal 納入核心矩陣，對混合 BAM 無來源分層

**關鍵限制**: 目前缺乏 tumor-only / normal-only / mixed 三路直接對照數據

**建議下一步**:
1. 使用 DORADO 重建 tumor-only / normal-only / mixed 三組，同一批位點比較 distance 分布
2. 在 InterSubMod 加入來源分層標籤輸出（RG tag 或 BAM 來源標記）
3. 對每個 purity 產出固定 `purity_qc.tsv`

**相關文件**:
- `docs/ai_sessions/2026/03/20260302_subsample混樣甲基化偏差_現況研究推論與驗證路線圖_01.md`
- `docs/ai_sessions/2026/02/20260213_HCC1395_subsample_purity_完整驗證報告_01.md`（Knowledge 路徑）

---

### 深入：雙向驗證（Label-First vs Cluster-First）

**研究問題**: 如何更準確地判斷甲基化距離矩陣中的結構是否有生物學意義？

**方法**:
- **Path A (Label-First)**：使用已知 HP/Allele 標籤分組，PERMANOVA 計算 Pseudo-F 與 p-value，Delta Distance = Between - Within
- **Path B (Cluster-First)**：層次聚類後，Bootstrap stability check（80% 重抽 20-50次），Fisher's Test 計算 Cramér's V

**4 類輸出判定矩陣**:

| 類別 | Label-First | Cluster-First | 解釋 |
|---|---|---|---|
| **Strong** | 顯著（高 R²）| 一致且穩定 | 標籤是甲基化變異主因 |
| **Subclone** | 不顯著 | 顯著且穩定 | 真實結構，但非 HP/Allele 主導 |
| **Weak** | 邊緣顯著 | 不穩定 | 效應存在但不明顯 |
| **Noise** | 不顯著 | 不穩定 | 隨機分群假象 |

**目前實作狀態**:
- PERMANOVA (`StructureTest.cpp`)：已實作，但主流程預設 `enable_permanova=false`
- Fisher (`GlobalTest.cpp`, `LocalTest.cpp`)：已啟用
- Bootstrap stability：已設計，實作待補齊
- LabelTest 多階段：已實作（Merged HP, Fine-grained HP, Allele, Unassigned affinity）

**關鍵限制**: 現有 PERMANOVA 呼叫基於 `cluster_labels`（非 HP/Allele 標籤），不符合 Label-First 定義；需修正

**建議下一步**:
1. 修正 PERMANOVA 呼叫，改用 HP/Allele 標籤作為分組依據
2. 補齊 Bootstrap stability check 的 C++ 實作
3. 建立 LabelTest 多階段 HP 驗證的獨立單元測試

**相關文件**:
- `docs/experiments/validated/2025/12/20251230_bidirectional_verification_strategy_analysis_01.md`
- `docs/architecture/system_overview.md`（S7: SignificanceAnalyzer）

---

### 深入：統計顯著性分析（PERMANOVA + Cramér's V）

**研究問題**: 哪些統計指標能有效度量甲基化分群的生物學顯著性？

**方法**:
- 5 層統計框架：GlobalTest（Fisher FFH）→ LocalTest（One-vs-Rest）→ StructureTest（PERMANOVA）→ LabelTest（HP/Allele）→ Bootstrap
- Cramér's V 作為效應量：`V = sqrt(χ²/n × 1/(min(r,c)-1))`

**關鍵數據**:
- Cramér's V 的 TP/FP 區分能力：AUC=0.52（幾乎無效）
- QUAL（VCF 原生）：AUC=0.9668；AF：AUC=0.9235
- 甲基化特徵整體 AUC 低，但在低 AF（<0.3）子群：NumReads AUC=0.84（有條件價值）
- 甲基化顯著位點中 94.8% 為 TP（高信心保留依據）

**失敗嘗試**: 使用 `Significant=True` 作為主要過濾器（F1=0.19）；使用 CramersV>0.1 過濾（損失過多 TP）

**核心發現**: 甲基化特徵不應作為獨立過濾器，而應作為 VCF 品質不足時的輔助依據（分層策略）

**建議下一步**:
1. 實作多重檢驗校正（FDR / Benjamini-Hochberg）
2. 建立分層過濾策略：QUAL≥0.8 直接接受 → QUAL<0.8 且甲基化顯著 → 保留（救援）→ 否則過濾
3. 在 LabelTest 加入 per-HP PERMANOVA（Label-First 框架）

**相關文件**:
- `docs/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md`（§4.2, §7.1）
- `docs/experiments/validated/2026/01/20260107_F1_and_Data_Optimization_Report_01.md`

---

## 附錄：待驗證方向（尚未正式啟動）

| 方向 | 期望收益 | 難度 | 依據 |
|---|---|---|---|
| 5hmC 雙通道距離矩陣 | 可能更有亞克隆特異性 | 高 | ONT 5kHz 同時提供 C+h |
| 跨 Region 亞克隆一致性 | 真正的亞克隆應跨多 SNV 一致 | 高 | 目前 per-SNV 獨立分析 |
| 機器學習組合特徵分類器 | 整合 15 個特徵的 ensemble | 中 | F1 研究發現特徵互補性 |
| PMD/ChromHMM Gating 啟用 | 降低甲基化背景噪聲 | 中 | 架構已有 `is_pmd` 欄位但未生效 |
| 多樣本交叉驗證（COLO829 等）| 泛化能力評估 | 中 | 目前僅 HCC1395 |
