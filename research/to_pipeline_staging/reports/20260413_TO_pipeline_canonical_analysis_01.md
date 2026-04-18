<!--
建立時間: 2026-04-13 18:30
目標: 使用 canonical 數據完整重建 TO pipeline 多階段 TP/FP 特徵分析
處理範圍: ClairS-TO → LongPhase-TO → ISM 三階段特徵分析（HCC1395 canonical + 7 樣本比較）
關聯檔案:
  - research/to_pipeline_staging/scripts/04_canonical_multi_stage_analysis.py
  - research/to_pipeline_staging/scripts/05_multi_sample_to_benchmark.py
  - research/to_pipeline_staging/scripts/06_canonical_plots.py
  - research/to_pipeline_staging/reports/20260412_TO_pipeline_multi_stage_characterization_01.md (v1,已校正)
-->

# TO Pipeline 多階段 TP/FP 特徵分析（Canonical 校正版）

## 0. 數據來源與校正說明

### 0.1 v1 → v2 校正

| 項目 | v1（已棄用） | v2（本報告） |
|------|-------------|-------------|
| **VCF 來源** | Clair 團隊 run (ONT BAM, 無 MM/ML, Sep 2025, 92MB) | canonical TO pipeline (ONT_5kHz BAM, 有 5mCG+5hmCG MM/ML, Mar 2026, 89MB) |
| **VCF 路徑** | `/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz` | `archive/step01_clairs_to/snv.vcf.gz` |
| **VCF PASS 數** | 47,798 | 48,085 |
| **ISM 來源** | `three_way_comparison/tumor_only_full/` (112K rows) | `archive/step05_intersubmod/` (40,213 rows) |
| **ISM 處理方式** | 全 VCF ISM → 再做 truth 交叉 | 先 truth 交叉 → 分別跑 TP/FP ISM |
| **BED 限制** | 未使用（將 outside-BED 當 FP） | 使用 truth BED（僅分析 BED 內 40,239 variants） |
| **ISM implicit filter FP removal** | ❌ 39.9%（錯誤！含 outside-BED） | ✅ 0.11%（正確，26 variants） |
| **F1 baseline** | ❌ 0.649（錯誤！含 outside-BED FP） | ✅ 0.7127（正確） |

**根因**：v1 將 48,085 PASS 中不在 truth BED 的 7,846 個 variants 全當作 FP，導致 FP 從 11,843 膨脹到 19,689。同時 three_way_comparison 的 ISM 只覆蓋部分 variants，造成「ISM 移除大量 FP」的假象。

### 0.2 Canonical 數據鏈

```
HCC1395 Canonical TO Pipeline (20260307_hcc1395_to_pilot_1)
├── Step 01: ClairS-TO → snv.vcf.gz (88MB, 3,203,048 variants)
│   └── PASS: 48,085 variants
├── Step 02: Benchmark vs SEQC2 truth (in HC BED)
│   ├── PASS in BED: 40,239
│   ├── TP: 28,396 | FP: 11,843 | FN: 11,051
│   └── F1: 0.7127
├── Step 03: LongPhase-TO
│   ├── 同樣 48,085 PASS（FILTER 完全不變）
│   ├── 新增 phasing: PS, GT2, PS2, GT3, PS3
│   ├── LOH.bed: 1,100 regions
│   └── GE.bed: 12,394 regions
├── Step 04: Benchmark LongPhase-TO
│   └── TP=28,396, FP=11,843（與 Step 02 完全相同）
└── Step 05: InterSubMod (分別對 TP/FP 執行)
    ├── TP processed: 28,383 / 28,396 (lost 13, 0.046%)
    ├── FP processed: 11,830 / 11,843 (lost 13, 0.110%)
    └── ISM implicit filter: 僅移除 26 variants
```

---

## 1. Stage 1: ClairS-TO PASS 分析

### 1.1 FILTER 分類統計

| FILTER | Count | 比例 |
|--------|------:|-----:|
| NonSomatic | 3,099,931 | 96.78% |
| LowQual+NonSomatic | 48,385 | 1.51% |
| **PASS** | **48,085** | **1.50%** |
| LowQual/Other | 6,340 | 0.20% |
| LowQual+NoAncestry | 200 | 0.01% |
| LowQual+MultiHap | 107 | 0.00% |
| **Total** | **3,203,048** | |

### 1.2 PASS 在 Truth BED 內的 Benchmark

| 指標 | 值 |
|------|---:|
| PASS in BED | 40,239 |
| TP | 28,396 |
| FP | 11,843 |
| FN | 11,051 |
| Precision | 0.7057 |
| Recall | 0.7199 |
| **F1** | **0.7127** |

**關鍵觀察**：ClairS-TO 在 TO 模式下 precision 僅 0.706。每 10 個 PASS call 中有 ~3 個是 FP（主要是 germline variants）。相比 paired mode F1=0.855，差距 0.14。

### 1.3 VCF 特徵對 TP/FP 的區分力

| 特徵 | AUC | TP mean | FP mean | 解讀 |
|------|----:|--------:|--------:|------|
| AF | 0.334 (inv: 0.666) | 0.499 | 0.666 | FP AF 顯著更高（germline ~0.5-1.0） |
| haplotype_flag (H) | 0.740 | 0.583 | 0.624 | FP 更常帶 H flag |
| SB | 0.521 | 0.589 | 0.589 | 幾乎無區分力 |
| DP | 0.514 | 86.5 | 85.6 | 幾乎無區分力 |
| QUAL | 0.472 (inv: 0.528) | 19.3 | 19.8 | 幾乎無區分力 |
| GQ | 0.513 | 18.8 | 19.3 | 幾乎無區分力 |

**重要發現**：
- **Verdict 全為零**：TO 模式不使用 Verdict 系統（Somatic/Germline/SubclonalSomatic 全部為空）
- **PoN 全為零**：PASS variants 不帶 PoN flags（已被 NonSomatic 過濾掉）
- **AF 是最有效的 VCF 特徵**（AUC=0.666 inverted），因 germline FP 的 AF 集中在 0.5-1.0，而 somatic TP 的 AF 在 0.3-0.6

---

## 2. Stage 2: LongPhase-TO 分析

### 2.1 LongPhase-TO 的作用

LongPhase-TO 在 canonical pipeline 中的角色：
1. **重新 phase** ClairS-TO 輸出 VCF → 添加 PS, GT2, PS2, GT3, PS3 FORMAT fields
2. **輸出 LOH/GE BED** → 標記 LOH 和 genotype equilibrium 區域
3. **Haplotag BAM** → 為 ISM 提供 HP-tagged reads

### 2.2 關鍵結論：LongPhase-TO 不改變 FILTER

```
ClairS-TO PASS: 48,085
LongPhase-TO PASS: 48,085 (完全相同)
```

**所有 phasing-related FILTER tags（NoAncestry, MultiHap）都只出現在 LowQual variants 上**，不影響 PASS variants。因此 LongPhase-TO 對 TP/FP 數量沒有任何影響。

### 2.3 Phasing 統計

| 類別 | Phased | Unphased | Phased Rate |
|------|-------:|---------:|------------:|
| TP (28,396) | 28,194 | 202 | **99.3%** |
| FP (11,843) | 11,649 | 194 | **98.4%** |

TP 和 FP 幾乎都被 phased（差異僅 0.9%）。LongPhase phasing 本身不是有效的 TP/FP 區分器（AUC 名義上很高但實際可用信息極少，因為兩者都 >98% phased）。

### 2.4 LongPhase-TO LOH.bed

| 類別 | In LOH | TP Rate |
|------|-------:|--------:|
| In LOH | 17,576 | 0.737 |
| Not in LOH | 22,663 | 0.682 |

LOH 區域的 TP rate（0.737）略高於非 LOH（0.682），差異 +0.055。這反映了 LOH 區域有更多 somatic variants（allelic imbalance 使得低 VAF somatic mutations 更容易被偵測）。

---

## 3. Stage 3: ISM 分析

### 3.1 ISM Implicit Filter（校正後）

**Canonical 數據確認：ISM implicit filter 效果可忽略。**

| 指標 | 值 |
|------|---:|
| ISM 輸入 TP | 28,396 |
| ISM 處理 TP | 28,383 (lost 13, **0.046%**) |
| ISM 輸入 FP | 11,843 |
| ISM 處理 FP | 11,830 (lost 13, **0.110%**) |
| 總移除 | 26 variants |
| Selectivity | 1.0× (TP 和 FP 等量移除) |

**與 v1 的差異**：v1 聲稱 ISM implicit filter 移除 39.9% FP（7,859 個），實際上是將 outside-BED variants 錯誤計為 FP 造成的假象。Canonical 數據顯示 ISM 幾乎能處理所有 BED 內 PASS variants。

**原因**：在 truth BED 高可信度區域內，幾乎所有 PASS variants 都有 ≥10 reads 和 ≥1 CpG，滿足 ISM 處理條件。

### 3.2 ISM SuggestFilter

| 指標 | 值 |
|------|---:|
| SuggestFilter = true (TP) | 124 |
| SuggestFilter = true (FP) | 95 |
| 總移除 | 219 |
| F1 before | 0.7126 |
| F1 after | **0.7114** |
| ΔF1 | **-0.0012**（負面！）|

**ISM SuggestFilter 反而降低 F1**，因為它移除的 TP (124) 多於 FP (95)。

### 3.3 ISM 特徵 AUC（全部 40,213 ISM-processed variants）

| 特徵 | AUC | 解讀 |
|------|----:|------|
| ism_SuggestFilter | 0.992 | 幾乎全零（binary, 極少觸發） |
| ism_HPMergedSig | 0.973 | HP-based，self-phasing confound |
| ism_CramersV | 0.955 | HP-based Cramér's V |
| ism_Significant | 0.957 | 綜合顯著性 flag |
| ism_HPMergedDelta | 0.917 | HP-based 甲基化差異 |
| ism_HPFineNGroups | 0.826 | HP-based 分群數 |
| ism_UnassignedAffinity | 0.824 | Unassigned read affinity |
| ism_PassedGating | 0.791 | Gating flag |
| ism_AlleleSig | 0.790 | Allele-based 顯著性 |
| ism_HP1FamilyN | 0.692 | HP1 family read count |
| ism_HP2FamilyN | 0.657 | HP2 family read count |
| ism_Quality_Score | 0.641 | 質量分數 |
| ism_HP_Ratio | 0.602 | HP read ratio |
| ism_AlleleDelta | 0.569 | Allele 甲基化差異 |
| ism_HeuristicScore | 0.643 | 啟發式分數 |

**重要觀察**：

1. **HP-based 特徵 AUC 很高（0.8-0.97）**——但這是 self-phasing 造成的：somatic TP 有真實 allelic reads → HP tags 品質好 → 甲基化差異大。FP (germline) 的 HP tags 也存在，但方向和模式不同。

2. **不代表這些特徵可用於過濾**——ISM SuggestFilter 已證明：用這些特徵過濾反而降低 F1。AUC 高不代表在操作閾值下有效。

3. **AlleleDelta（非 HP 依賴）AUC 僅 0.569**——幾乎無用。

### 3.4 ISM VerificationClass 分布

| VerificationClass | TP | FP | TP Rate | 解讀 |
|-------------------|---:|---:|--------:|------|
| Weak | 10,062 | 3,281 | 0.754 | 最多 TP，最佳 TP rate |
| Strong | 6,442 | 2,244 | 0.742 | 第二 |
| Noise | 10,527 | 5,363 | 0.663 | FP 最集中 |
| Subclone | 1,352 | 942 | 0.589 | TP rate 最低 |

Noise 類佔所有 FP 的 45.3%，但 Noise TP rate 仍有 0.663（不能簡單過濾掉 Noise）。

### 3.5 ISM Coverage_Category 分布

| Coverage_Category | TP | FP | TP Rate |
|-------------------|---:|---:|--------:|
| Normal | 10,424 | 4,381 | 0.704 |
| Low | 9,139 | 3,808 | 0.706 |
| Elevated | 3,641 | 1,464 | 0.713 |
| CNV_Loss | 2,529 | 1,142 | 0.689 |
| CNV_Gain | 2,161 | 804 | 0.729 |
| High_Copy | 489 | 231 | 0.679 |

TP rate 在所有 Coverage_Category 之間差異很小（0.679-0.729），說明 ISM Coverage_Category 不是有效的分層器。

---

## 4. CNV/LOH 外部註釋分析

### 4.1 CNV Type 分布

| CNV Type | TP | FP | Total | TP Rate |
|----------|---:|---:|------:|--------:|
| gain | 15,406 | 6,821 | 22,227 | 0.693 |
| loh | 11,866 | 4,394 | 16,260 | 0.730 |
| loss | 146 | 68 | 214 | 0.682 |
| none (neutral) | 978 | 560 | 1,538 | 0.636 |

**重要修正**：v1 報告 neutral TP rate = 16.9%（錯誤），canonical 數據顯示 neutral TP rate = 0.636。差異原因：v1 將 outside-BED variants 計為 FP，而 outside-BED 區域大多是 neutral（沒有 CNV），導致 neutral FP 嚴重膨脹。

### 4.2 AF 在 CNV × TP/FP 的交互作用

| CNV Type | TP AF mean | FP AF mean | Δ(FP-TP) |
|----------|----------:|----------:|---------:|
| gain | 0.512 | 0.700 | +0.188 |
| loh | 0.499 | 0.632 | +0.133 |
| none | 0.401 | 0.558 | +0.157 |

在所有 CNV 類型中，FP 的 AF 一致高於 TP。Gain 區域差異最大（+0.188），因為 gain 複製數增加使得 germline heterozygous alleles 的 AF 偏向 0.67 或 0.75。

---

## 5. 跨樣本比較（7 樣本）

### 5.1 ClairS-TO PASS F1 Summary

| Sample | PASS | TP | FP | FN | Precision | Recall | **F1** | 來源 |
|--------|-----:|---:|---:|---:|----------:|-------:|-------:|------|
| H2009 | 153,854 | 125,707 | 11,990 | 16,384 | 0.913 | 0.885 | **0.899** | orthogonal |
| HCC1395_DORADO | 48,072 | 28,861 | 11,576 | 10,586 | 0.714 | 0.732 | **0.723** | SEQC2 |
| HCC1395 | 48,085 | 28,396 | 11,843 | 11,051 | 0.706 | 0.720 | **0.713** | SEQC2 canonical |
| COLO829 | 50,991 | 33,285 | 17,706 | 9,907 | 0.653 | 0.771 | **0.707** | NYGC |
| H1437 | 65,448 | 45,474 | 13,443 | 35,542 | 0.772 | 0.561 | **0.650** | orthogonal |
| HCC1937 | 29,503 | 12,624 | 12,063 | 4,243 | 0.511 | 0.748 | **0.608** | orthogonal |
| HCC1954 | 74,996 | 17,068 | 50,218 | 5,962 | 0.254 | 0.741 | **0.378** | orthogonal |

**觀察**：
- F1 範圍從 0.378（HCC1954）到 0.899（H2009），差異極大
- **H2009**（高突變負荷 168K truth）表現最好，因 somatic:germline 比例高 → precision 高
- **HCC1954** 表現最差：precision 僅 0.254（每 4 個 PASS 有 3 個 FP），因 HCC1954 是高倍體（近四倍體）→ germline FP 極多
- **HCC1937** precision 0.511 → 接近隨機猜測

### 5.2 AF Distribution（跨樣本一致性）

| Sample | TP AF mean | FP AF mean | Δ(FP-TP) | AF AUC (eff) | 解讀 |
|--------|----------:|----------:|---------:|-------------:|------|
| HCC1937 | 0.420 | 0.719 | **+0.299** | 0.760 | 最大差距——低純度 + 高 germline AF |
| HCC1954 | 0.295 | 0.527 | **+0.233** | 0.802 | TP AF 極低（近四倍體）|
| HCC1395_DORADO | 0.500 | 0.684 | +0.184 | 0.672 | 典型 HCC1395 |
| HCC1395 | 0.499 | 0.666 | +0.167 | 0.666 | Canonical |
| H1437 | 0.553 | 0.616 | +0.062 | 0.555 | 小差距 |
| H2009 | 0.558 | 0.611 | +0.054 | 0.561 | 小差距（高精度樣本） |
| COLO829 | 0.602 | 0.599 | **-0.003** | 0.519 | TP/FP AF 幾乎相同！ |

**跨樣本一致模式**：FP AF > TP AF（6/7 樣本一致），但強度差異大。COLO829 例外（AF 無法區分 TP/FP，可能因 COLO829 高純度 + LOH 廣泛）。

### 5.3 Haplotype Flag（H-flag）一致性

H-flag 是 ClairS-TO VCF 中 INFO 欄位 `H` 的存在與否。

| Sample | H-flag AUC (eff) |
|--------|----------------:|
| HCC1954 | 0.790 |
| COLO829 | 0.776 |
| H2009 | 0.767 |
| H1437 | 0.751 |
| HCC1395_DORADO | 0.743 |
| HCC1937 | 0.742 |
| HCC1395 | 0.740 |

**H-flag AUC 在所有 7 樣本中一致高（0.740-0.790）**。這是 ClairS-TO VCF 中最穩定的 TP/FP 區分器。H flag 表示該 variant 在 haplotype-based calling 中被辨識，FP（germline）更常帶此 flag。

### 5.4 TO vs Paired F1 比較

| Sample | Paired F1 | TO F1 | Δ(TO-Paired) |
|--------|----------:|------:|-------------:|
| HCC1395 | 0.855 | 0.713 | **-0.142** |
| HCC1395_DORADO | 0.866 | 0.723 | **-0.143** |
| COLO829 | 0.887 | 0.707 | **-0.180** |
| H1437 | 0.920 | 0.650 | **-0.270** |
| H2009 | 0.971 | 0.899 | -0.072 |
| HCC1937 | 0.369* | 0.608 | +0.239† |
| HCC1954 | 0.908 | 0.378 | **-0.530** |

*HCC1937 paired_pileup F1=0.369 異常低（已知問題），paired_full F1=0.842
†HCC1937 TO > paired_pileup 是因為 paired_pileup 有已知問題

**中位 ΔF1 = -0.143**：TO 模式比 paired 模式平均損失 0.14 F1。HCC1954 損失最嚴重（-0.53），因為近四倍體 germline 在無 normal 的情況下完全無法過濾。

### 5.5 跨樣本結論

1. **TO F1 受樣本特性影響極大**（0.378-0.899），不如 paired mode 穩定（0.842-0.971 排除 HCC1937 pileup 問題）
2. **高倍體/低純度樣本 TO 表現極差**（HCC1954, HCC1937）
3. **AF 是最一致的 TP/FP 區分器**，但 COLO829 例外
4. **H-flag 是最穩定的特徵**（7/7 樣本 AUC > 0.74）
5. **QUAL, GQ, DP, SB 在所有樣本中均無區分力**（AUC ≈ 0.50）

---

## 6. 綜合結論

### 6.1 Pipeline 各階段貢獻（校正版）

| 階段 | 動作 | TP/FP 變化 | F1 變化 |
|------|------|-----------|---------|
| **ClairS-TO** | FILTER 分類 | 3.2M → 48K PASS | 基線 F1=0.713 |
| **LongPhase-TO** | 添加 phasing info | 無變化 | 無變化 |
| **ISM implicit** | 區域處理門檻 | -13 TP, -13 FP | -0.0001 |
| **ISM SuggestFilter** | 甲基化過濾 | -124 TP, -95 FP | **-0.0012** |

### 6.2 v1 vs v2 核心差異

| 結論 | v1（錯誤） | v2（正確） |
|------|-----------|-----------|
| ISM implicit filter 移除 FP 比例 | 39.9% | **0.11%** |
| ISM implicit filter F1 提升 | +0.064 | **-0.0001** |
| neutral TP rate | 16.9% | **63.6%** |
| ISM SuggestFilter 效果 | 未詳細分析 | **負面 (ΔF1=-0.0012)** |
| 最強特徵 | CNV (in_gain AUC=0.84) | **AF (AUC=0.666 inv) + HP-based ISM** |

### 6.3 v1 錯誤的根因

1. **VCF 來源錯誤**：使用了 Clair 團隊的 ClairS-TO run（ONT BAM, 無 MM/ML, Sep 2025）而非 canonical TO pipeline（ONT_5kHz BAM, 有 5mCG+5hmCG, Mar 2026）
2. **ISM 數據不匹配**：three_way_comparison ISM (112K rows) 對應的是不同 VCF
3. **缺少 BED 限制**：將 outside-BED variants 當作 FP 計算
4. **ISM 覆蓋率錯估**：112K ISM rows 無法 1:1 匹配 47K PASS variants → 誤判為 ISM 移除了大量 FP

### 6.4 對研究方向的影響

- **ISM 在 TO 模式下不具 FP 過濾能力**：implicit filter 可忽略，SuggestFilter 反面
- **HP-based ISM 特徵有高 AUC 但不可操作**：self-phasing confound → 高 AUC ≠ 可用過濾器
- **唯一有效的特徵是 AF**（AUC=0.666）和 CNV 外部註釋（LOH TP rate +0.05）
- **需要 Normal BAM 或外部信息才能根本改善 TO FP**

---

## 7. 圖表索引

### HCC1395 Canonical 分析

| 圖號 | 檔名 | 內容 |
|------|------|------|
| C01 | C01_canonical_pipeline_waterfall.png | Pipeline waterfall + TP/FP per stage |
| C02 | C02_canonical_feature_auc_ranking.png | 全特徵 AUC 排名 |
| C03 | C03_canonical_vcf_features.png | VCF 特徵 TP/FP 分布 (AF, DP, QUAL, GQ, SB, Phasing) |
| C04 | C04_canonical_cnv_stratification.png | CNV/LOH TP rate 分層 |
| C05 | C05_canonical_ism_features.png | ISM 特徵分布 (6 panels) |
| C06 | C06_canonical_ism_filter_detail.png | ISM implicit + SuggestFilter 細節 |

### 多樣本比較

| 圖號 | 檔名 | 內容 |
|------|------|------|
| M01 | M01_multi_sample_f1_comparison.png | 7 樣本 F1/Precision/Recall 比較 |
| M02 | M02_multi_sample_af_distribution.png | 7 樣本 AF 分布 (TP vs FP) |
| M03 | M03_multi_sample_feature_auc_heatmap.png | VCF 特徵 AUC 熱圖 |
| M04 | M04_to_vs_paired_f1.png | TO vs Paired F1 比較 |

---

## 8. 數據文件索引

| 檔名 | 內容 | 行數 |
|------|------|------|
| `hcc1395_canonical_multimodal.csv` | 全特徵 DataFrame | 40,239 |
| `hcc1395_canonical_feature_auc.csv` | 特徵 AUC 排名 | 38 |
| `hcc1395_canonical_stage_metrics.json` | 各階段 F1 指標 | — |
| `multi_sample_to_summary.json` | 多樣本 benchmark 結果 | — |
| `multi_sample_to_features.csv` | 多樣本 VCF 特徵 | — |

## 9. 腳本索引

| 腳本 | 用途 | 版本 |
|------|------|------|
| `04_canonical_multi_stage_analysis.py` | HCC1395 canonical 三階段分析 | v2 (校正) |
| `05_multi_sample_to_benchmark.py` | 7 樣本 ClairS-TO benchmark | new |
| `06_canonical_plots.py` | Canonical 圖表生成 | v2 |
| `01_multi_stage_characterization.py` | ❌ 已棄用（使用 zhenyu112 VCF） | v1 (廢棄) |
| `02_generate_plots.py` | ❌ 已棄用 | v1 (廢棄) |
| `03_updated_waterfall_plot.py` | ❌ 已棄用 | v1 (廢棄) |
