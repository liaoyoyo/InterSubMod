# TP/FP Heatmap 視覺化案例比較報告

**專案**: InterSubMod - 甲基化輔助體細胞突變檢測  
**分析日期**: 2026-01-23  
**目的**: 透過 heatmap 視覺化比較 TP 與 FP 的甲基化模式差異

---

## 報告說明

本報告精選代表性的 **FP 高風險案例** 與 **TP 良好案例**，並展示其：
- **Cluster Heatmap**: 甲基化聚類熱圖
- **Distance Heatmap**: 甲基化距離熱圖
- **關鍵數據**: AlleleDelta, CramersV, VAF, QUAL 等

> [!TIP]
> **觀察重點**: 
> - 左側 annotation 三欄分別為 **HP (haplotype)**、**Strand**、**Allele**
> - 比對 Allele 標記（Ref/Alt）與甲基化模式（紅/藍）是否有對應關係
> - CramersV = 0 表示 Allele 與甲基化之間**無統計關聯**

---

## FP 高風險案例 (Strong + 保留)

這些是過濾後仍保留的高風險假陽性，需要特別檢視。

---

### FP 案例 1: chr6:98149813 (C>T)

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | 0.3867 | ⚠️ 極高 |
| **CramersV** | 0.0000 | ❌ 無統計關聯 |
| **VAF** | 0.3571 | 高於 0.24 門檻 |
| **QUAL** | 0.833 | 較高 |
| **NumReads** | 11 | ⚠️ 低覆蓋 |
| **VerificationClass** | Strong | |
| **DominantLabel** | hp | |
| **PassedGating** | True | |

**為何未被過濾**: VAF (0.36) > 0.24 門檻

#### Cluster Heatmap

![chr6:98149813 Cluster Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr6/chr6_98149813/chr6_98144813_98154813/plots/BERNOULLI/cluster_heatmap.png)

#### Distance Heatmap

![chr6:98149813 Distance Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr6/chr6_98149813/chr6_98144813_98154813/plots/BERNOULLI/distance_heatmap.png)

**觀察重點**:
- Allele (Ref/Alt) 與甲基化模式是否對應？
- 高 AlleleDelta 但 CramersV=0，表示差異不具統計意義
- NumReads 僅 11，統計力可能不足

---

### FP 案例 2: chr12:71376990 (A>C)

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | 0.3051 | ⚠️ 高 |
| **CramersV** | 0.0000 | ❌ 無統計關聯 |
| **VAF** | 0.9481 | ⚡ 極高 (接近 1.0) |
| **QUAL** | 0.9659 | 高 |
| **NumReads** | 61 | 適中 |
| **VerificationClass** | Strong | |
| **DominantLabel** | allele | |

**為何未被過濾**: VAF (0.95) >> 0.24 門檻

#### Cluster Heatmap

![chr12:71376990 Cluster Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr12/chr12_71376990/chr12_71371990_71381990/plots/BERNOULLI/cluster_heatmap.png)

#### Distance Heatmap

![chr12:71376990 Distance Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr12/chr12_71376990/chr12_71371990_71381990/plots/BERNOULLI/distance_heatmap.png)

**觀察重點**:
- VAF 接近 1.0，表示幾乎所有 reads 都是 Alt
- 這可能是 germline 變異或技術偽影
- 高 AD + 高 VAF 的組合需額外研究

---

### FP 案例 3: chr8:82168704 (G>T)

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | 0.2584 | 略高於門檻 |
| **CramersV** | **1.0000** | ✅ 強統計關聯 |
| **VAF** | 0.2529 | 略高於門檻 |
| **QUAL** | 0.9804 | 高 |
| **NumReads** | 72 | 良好 |
| **VerificationClass** | Strong | |
| **Significant** | True | ✅ |

**為何未被過濾**: CramersV (1.0) > 0.05 門檻

#### Cluster Heatmap

![chr8:82168704 Cluster Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr8/chr8_82168704/chr8_82163704_82173704/plots/BERNOULLI/cluster_heatmap.png)

#### Distance Heatmap

![chr8:82168704 Distance Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr8/chr8_82168704/chr8_82163704_82173704/plots/BERNOULLI/distance_heatmap.png)

**觀察重點**:
- CramersV = 1.0 表示 Allele 與甲基化有**完美關聯**
- 這可能是真正的 Allele-Specific Methylation (ASM)
- 但仍是 FP (不在 SEQC2 gold standard)，為何？需深入研究

---

### FP 案例 4: chr1:37177717 (G>A)

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | 0.2454 | 接近門檻 |
| **CramersV** | 0.0000 | ❌ 無統計關聯 |
| **VAF** | 0.1333 | 低 |
| **QUAL** | 0.8122 | 中高 |
| **NumReads** | 38 | 中等 |
| **VerificationClass** | Strong | |
| **DominantLabel** | hp | |

**為何未被過濾**: AlleleDelta (0.245) < 0.25 門檻

#### Cluster Heatmap

![chr1:37177717 Cluster Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr1/chr1_37177717/chr1_37172717_37182717/plots/BERNOULLI/cluster_heatmap.png)

#### Distance Heatmap

![chr1:37177717 Distance Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_fp/filtered_snv_fp/chr1/chr1_37177717/chr1_37172717_37182717/plots/BERNOULLI/distance_heatmap.png)

**觀察重點**:
- AD 剛好在門檻邊緣 (0.245 vs 0.25)
- 若放鬆 AD 門檻至 0.24，此位點會被過濾

---

## TP 良好案例 (Strong + 高品質)

這些是典型的真陽性，特徵符合預期。

---

### TP 案例 1: chr8:115093788 (G>A)

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | -0.0089 | ✅ 接近 0 |
| **CramersV** | 0.0000 | 正常 |
| **VAF** | 0.4167 | ✅ 良好 |
| **QUAL** | 0.9931 | ✅ 極高 |
| **NumReads** | 97 | ✅ 高覆蓋 |
| **VerificationClass** | Strong | |
| **DominantLabel** | hp | |

#### Cluster Heatmap

![chr8:115093788 Cluster Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp/chr8/chr8_115093788/chr8_115088788_115098788/plots/BERNOULLI/cluster_heatmap.png)

#### Distance Heatmap

![chr8:115093788 Distance Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp/chr8/chr8_115093788/chr8_115088788_115098788/plots/BERNOULLI/distance_heatmap.png)

**觀察重點**:
- AlleleDelta ≈ 0，表示 Ref/Alt 甲基化無差異
- 這是典型的**真正體細胞突變**特徵
- 高 QUAL + 高 VAF + 低 AD = 高品質 TP

---

### TP 案例 2: chr3:68729618 (G>A)

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | 0.0000 | ✅ 完美 |
| **CramersV** | 0.0000 | 正常 |
| **VAF** | 0.9740 | ✅ 極高 |
| **QUAL** | 0.9973 | ✅ 極高 |
| **NumReads** | 69 | 良好 |
| **VerificationClass** | Strong | |
| **DominantLabel** | hp | |

#### Cluster Heatmap

![chr3:68729618 Cluster Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp/chr3/chr3_68729618/chr3_68724618_68734618/plots/BERNOULLI/cluster_heatmap.png)

#### Distance Heatmap

![chr3:68729618 Distance Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp/chr3/chr3_68729618/chr3_68724618_68734618/plots/BERNOULLI/distance_heatmap.png)

**觀察重點**:
- AlleleDelta = 0，Ref/Alt 甲基化模式完全一致
- VAF 接近 1.0，可能是 homozygous 突變
- 這是**理想的 TP**特徵

---

### TP 案例 3: chr9:30374811 (G>T)

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | 0.0000 | ✅ 完美 |
| **CramersV** | 0.0000 | 正常 |
| **VAF** | 0.9538 | ✅ 極高 |
| **QUAL** | 0.9983 | ✅ 極高 |
| **NumReads** | 54 | 適中 |
| **VerificationClass** | Strong | |

#### Cluster Heatmap

![chr9:30374811 Cluster Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp/chr9/chr9_30374811/chr9_30369811_30379811/plots/BERNOULLI/cluster_heatmap.png)

#### Distance Heatmap

![chr9:30374811 Distance Heatmap](../../output/bip8_output_archive/20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp/chr9/chr9_30374811/chr9_30369811_30379811/plots/BERNOULLI/distance_heatmap.png)

---

## 被過濾移除的 FP 案例 (成功識別)

用於對照，這些是被過濾條件成功移除的 FP。

---

### 移除 FP 案例: chr12:38697849

| 特徵 | 數值 | 說明 |
|------|------|------|
| **AlleleDelta** | 0.7135 | ⚠️ 極高 |
| **CramersV** | 0.0000 | ❌ 無統計關聯 |
| **VAF** | 0.16 | 低於門檻 |
| **QUAL** | 0.54 | 低 |

**為何被過濾**: AD > 0.25, V < 0.05, VAF < 0.24 ✅

> 參考圖片位於: `output/ml_feature_exploration/2026/01/reports/removed_sites_analysis/`

---

## TP vs FP 特徵對比總結

| 特徵 | 良好 TP 特徵 | 高風險 FP 特徵 |
|------|-------------|---------------|
| **AlleleDelta** | ≈ 0 (無差異) | > 0.25 (高差異) |
| **VAF** | > 0.4 (高) | 0.1-0.3 (低中) |
| **QUAL** | > 0.95 (高) | 0.6-0.9 (中) |
| **CramersV** | 0 或低 | 0 或異常高 (1.0) |
| **NumReads** | > 50 | 10-50 |

---

## 欠缺考慮的問題

透過以上視覺化比較，發現以下需要額外考慮的情況：

### 1. 高 VAF FP

- **案例**: chr12:71376990 (VAF = 0.95)
- **問題**: 現有過濾無法識別高 VAF 的 FP
- **建議**: 需結合其他特徵 (如 QUAL、區域複雜度)

### 2. 高 CramersV FP

- **案例**: chr8:82168704 (V = 1.0)
- **問題**: CramersV > 0.05 會讓位點躲過過濾，但有些 FP 也有高 V
- **建議**: CramersV 高不代表一定是 TP，需要額外驗證

### 3. 邊緣 AlleleDelta

- **案例**: chr1:37177717 (AD = 0.245)
- **問題**: 剛好在 0.25 門檻邊緣
- **建議**: 考慮放鬆 AD 門檻至 0.23-0.24

### 4. 低 NumReads

- **案例**: chr6:98149813 (NumReads = 11)
- **問題**: 低覆蓋度的統計可靠性差
- **建議**: 加入 NumReads < 20 的過濾條件

---

## 建議後續行動

1. **調整 AlleleDelta 門檻**: 從 0.25 降至 0.23
2. **加入 NumReads 條件**: NumReads < 15 視為可疑
3. **研究高 VAF FP**: 這類 FP 可能需要不同的識別策略
4. **研究高 CramersV FP**: chr8:82168704 案例值得深入分析

---

**報告生成時間**: 2026-01-23  
**資料目錄**: `/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260119_all-with-w5000_1/`
