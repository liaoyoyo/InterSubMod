<!--
建立時間: 2026-04-12 18:30
目標: 驗證文獻假說 L1（方向性 ASM）在 ISM 748K 區域數據中的表現
處理範圍: 7 樣本 × TP/FP × 340K 區域 aggregate + 5,360 區域 raw reads
關聯檔案:
  - research/literature_validation/scripts/L1_directional_ASM.py
  - research/literature_validation/scripts/L1_supplementary_AF_stratified.py
  - research/literature_validation/figures/01-04_*.png
-->

# L1 方向性 ASM 文獻驗證報告

## 1. 研究背景與假說

### 文獻預測（epiTRACERx / Gaiti et al. / Basilicata et al.）

多篇高影響力論文指出：

1. **epiTRACERx** (Basilicata, Nature 2024): 體細胞突變位點存在 allele-specific methylation (ASM)，ALT allele 的甲基化模式與 REF 不同，可追蹤腫瘤演化
2. **Gaiti et al.** (Nature Biotechnology 2019): scNMT-seq 在 CLL 中觀察到 somatic driver positions 的 ASM
3. **EVOFLUx** (Moran/Gaiti, Nature Biotechnology 2023): fluctuating CpGs (fCpG) 比穩定 CpGs 更能追蹤亞克隆動態

### 測試假說

- **L1a**: TP somatic SNV 位點的 ALT reads 甲基化程度系統性不同於 REF reads（方向性 ASM）
- **L1b**: AlleleDelta (PERMANOVA) 在 TP 與 FP 之間有顯著差異
- **L1c**: CpG 變異度（fCpG proxy）在 TP 與 FP 之間不同
- **L1d**: |ASM| (絕對 allele-specific methylation 差異) 可區分 TP 與 FP

---

## 2. 數據規模

| 類別 | Part A (aggregate) | Part B (raw reads) |
|------|-------------------|--------------------|
| 數據來源 | significance_summary.csv | reads.tsv + methylation.csv |
| 區域數 | 340,173 | 5,360 (抽樣 500/類別/樣本) |
| 樣本數 | 7 paired | 7 paired |
| 分析方式 | AlleleDelta, AlleleSig | 方向性甲基化差 (ALT_mean - REF_mean) |

---

## 3. 核心發現

### 3.1 L1a: 方向性 ASM — **TP 無方向性，FP 有方向性**

| 指標 | TP (n=3,062) | FP (n=2,298) | 差異 |
|------|-------------|-------------|------|
| signed_delta 均值 | **-0.0003** | **-0.0205** | FP 顯著偏向 REF>ALT |
| signed_delta 中位數 | -0.0007 | -0.0056 | 同上 |
| t-test (H0: delta=0) | p=0.854 | **p=1.16e-14** | TP 不顯著，FP 高度顯著 |
| ALT>REF 比例 | 49.4% | 44.3% | TP 無偏向，FP 偏向 REF |

**結論**：TP somatic SNV 的甲基化完全無方向性偏好（49.4% vs 50.6%，完全隨機）。FP（主要為 germline variants）顯示顯著的 REF>ALT 甲基化偏好。

**文獻對照**：epiTRACERx 的 ASM 觀察是在 **已知 driver mutations** 的 CpG island/promoter 區域，而我們的 TP 包含大量 **passenger mutations**，passenger 不預期有方向性 ASM。

### 3.2 L1b: AlleleDelta (PERMANOVA) — **FP > TP**

| 指標 | TP (n=333,452) | FP (n=6,721) | 統計 |
|------|---------------|-------------|------|
| AlleleDelta 均值 | 0.0200 | **0.0312** | FP 更高 |
| Mann-Whitney p | — | — | **4.71e-04** |
| Cohen's d | — | — | **-0.19** (小效果，FP 方向) |

FP 的 allele-based 甲基化距離 **大於** TP → germline variants 有更強的 ASM，與先前 O1-O13 結論一致。

### 3.3 L1c: CpG Variability (fCpG) — **完全無差異**

| 指標 | TP | FP | p 值 |
|------|----|----|------|
| CpG variance mean | 0.0828 | 0.0825 | **0.767** |
| High-var CpG fraction | 0.510 | 0.530 | — |

**結論**：fCpG 概念（EVOFLUx）在 TP/FP 判別上**完全無效**。TP 和 FP 的 CpG 變異度分佈幾乎完全重疊。

### 3.4 L1d: |ASM| Magnitude — **FP >> TP（反向判別）**

| 指標 | TP | FP | p 值 |
|------|----|----|------|
| |ASM| per CpG | 0.0416 | **0.0699** | **1.43e-24** |
| |AlleleDelta| AUC | — | — | **0.5353** |

**結論**：|ASM| 可以區分 TP 和 FP，但方向**相反** — FP（germline）有更大的 ASM。AUC = 0.535，判別力微弱。

---

## 4. 分層驗證

### 4.1 AF 分層：ASM 差異不隨 AF 改變

| AF Bin | TP_mean | FP_mean | TP-FP | p |
|--------|---------|---------|-------|---|
| [0.0, 0.1) | 0.0143 | 0.0295 | -0.015 | 0.21 |
| [0.1, 0.2) | 0.0219 | 0.0492 | -0.027 | 2.1e-19 |
| [0.2, 0.3) | 0.0206 | 0.0344 | -0.014 | 2.4e-03 |
| [0.3, 0.4) | 0.0216 | 0.0241 | -0.003 | 3.4e-05 |
| [0.4, 0.5) | 0.0230 | 0.0327 | -0.010 | 0.14 |

每個 AF bin 中 FP AlleleDelta ≥ TP，但差異不大且不一致。

### 4.2 甲基化水平分層：FP 的 Mid-methylation 最強

| 甲基化水平 | TP signed_delta | FP signed_delta |
|-----------|----------------|----------------|
| Low (<0.3) | -0.006 | -0.009 |
| **Mid (0.3-0.7)** | +0.005 | **-0.041** |
| High (≥0.7) | -0.000 | -0.011 |

FP 在中度甲基化區域的 ASM 最強（delta=-0.041），而 TP 在所有水平幾乎為零。

### 4.3 VerificationClass 分層

| Class | TP AlleleDelta | FP AlleleDelta | Delta |
|-------|---------------|---------------|-------|
| Noise | 0.0003 | 0.0001 | +0.0002 |
| Weak | 0.0237 | 0.0391 | **-0.0154** |
| **Strong** | **0.0693** | **0.1153** | **-0.0460** |
| Subclone | 0.0007 | 0.0057 | -0.0049 |

Strong class 的差異最大：FP Strong (0.1153) >> TP Strong (0.0693)。

### 4.4 LOH 分層

Non-LOH 和 LOH 子集均顯示 FP ≥ TP 的 AlleleDelta，但差異不顯著（p=0.38, 0.21）。

---

## 5. 核心結論

### L1a 方向性 ASM: **NEGATIVE**
TP somatic SNV 位點不存在系統性方向性 ASM。49.4% ALT>REF 完全是隨機噪聲。epiTRACERx 的觀察可能限於 driver mutations / CpG island promoters。

### L1b AlleleDelta TP vs FP: **INVERTED — FP 更高**
germline ASM >> somatic ASM，與 14 個月研究結論一致。不可用於 TP 鑑別。

### L1c fCpG (EVOFLUx): **NEGATIVE — 完全無差異**
TP 和 FP 的 CpG variability 分佈完全重疊。fCpG 選擇策略無法改善 TP/FP 判別。

### L1d |ASM| Discrimination: **WEAK INVERTED — AUC=0.535**
|ASM| 有微弱判別力但方向相反（FP > TP）。不足以作為獨立過濾器。

---

## 6. 對 ISM 研究方向的影響

### 確認的結論
1. **Pure methylation features 的 AUC 天花板 ≤ 0.58 再次驗證** — 方向性 ASM 也無法突破
2. **germline ASM >> somatic ASM** 是根本性生物學限制，不是 ISM 方法缺陷
3. **fCpG 篩選不會幫助 TP/FP 判別** — EVOFLUx 的 fCpG 概念在此場景無效

### 新發現
1. **FP 的方向性 ASM (REF>ALT)** 值得深入研究 — 如果這是 germline variant 的特徵，理論上可以用反向策略：識別有方向性 ASM 的位點標記為 FP 候選
2. **Mid-methylation FP** 的 ASM 最強（delta=-0.041）— 這些可能對應 imprinted regions 或 mQTL

### 建議後續
- **不建議**繼續投入 L1 方向性 ASM 作為 TP 改善方向
- **不建議**實作 fCpG 篩選（L4）用於 TP/FP 判別
- **可考慮**：FP 的 directional ASM 特徵作為 germline filter 輔助信號（需與 AF confound 解耦）
- **維持**：Phase 2A normal methylation baseline 仍是最有潛力方向（與本結果正交）

---

## 7. 圖表索引

| 圖號 | 內容 | 路徑 |
|------|------|------|
| 01 | AlleleDelta TP/FP violin + per-sample bars | `figures/01_allele_delta_aggregate.png` |
| 02 | 方向性 ASM 直方圖 + CpG variance + |ASM| | `figures/02_directional_asm_detail.png` |
| 03 | ALT vs REF 平均甲基化散點圖 | `figures/03_alt_ref_scatter.png` |
| 04 | AF 分層 + 甲基化水平分層補充分析 | `figures/04_supplementary_af_stratified.png` |
