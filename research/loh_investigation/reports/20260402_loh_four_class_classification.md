<!--
建立時間: 2026-04-02 19:00
目標: 以 SEQC2 CNV Benchmark 為 ground truth，將 LongPhase-TO LOH.bed 分類為 TP/FP/FN/TN 四類，並留存 BED 檔與視覺化
處理範圍: HCC1395 (5kHz) LOH.bed vs SEQC2 LOH benchmark 交集分析
關聯檔案:
  - research/loh_investigation/data/seqc2_cnv_benchmark_v4/ngs_benchmark_cnv_gain_loss_loh.bed
  - research/loh_investigation/data/seqc2_vs_to_validation/
  - research/loh_investigation/reports/20260402_seqc2_cnv_benchmark_analysis.md
  - research/loh_investigation/reports/20260402_seqc2_vs_longphase_to_loh_validation.md
-->

# LOH 四類分類報告 — SEQC2 Benchmark vs LongPhase-TO LOH.bed

> **核心結論**：LongPhase-TO LOH.bed 在體染色體上高度準確（F1=96.2%, Jaccard=0.928）。TO-only 假 LOH 主要來自 chrX（女性細胞株半合子），而非 self-phasing artifact。此發現**修正了先前認為 LOH.bed 嚴重受 self-phasing 汙染的假設**。

---

## 1. 四類分類定義

| 類別 | 定義 | 含義 |
|------|------|------|
| **TP** (True Positive) | SEQC2 = LOH **且** TO LOH.bed = LOH | 兩系統一致判定為 LOH |
| **FP** (False Positive) | SEQC2 ≠ LOH **但** TO LOH.bed = LOH | TO 額外判定的 LOH（可能是 artifact） |
| **FN** (False Negative) | SEQC2 = LOH **但** TO LOH.bed ≠ LOH | TO 漏判的 LOH |
| **TN** (True Negative) | SEQC2 ≠ LOH **且** TO LOH.bed ≠ LOH | 兩系統一致判定為非 LOH |

---

## 2. 全基因組分類統計

### 2.1 bp-Level Confusion Matrix

|  | SEQC2=LOH | SEQC2=non-LOH |
|--|-----------|---------------|
| **TO=LOH** | **1,432 Mb** (TP) | **58 Mb** (FP*) |
| **TO=non-LOH** | **200 Mb** (FN*) | **1,340 Mb** (TN) |

> *注：此處 FP/FN 定義以 SEQC2 為 ground truth。

### 2.2 Region-Level 統計

| 類別 | Regions | Total Mb | 佔基因組 % |
|------|---------|----------|-----------|
| TP | 1,066 | 1,432 | 47.3% |
| FP | 524 | 200 | 6.6% |
| FN | 148 | 58 | 1.9% |
| TN | — | 1,340 | 44.2% |

### 2.3 準確性指標

| 指標 | 含 chrX | 體染色體 only |
|------|---------|-------------|
| **Sensitivity** | 96.1% | 96.1% |
| **Precision** | 87.7% | **96.4%** |
| **F1 Score** | 91.7% | **96.2%** |
| **Jaccard** | 84.7% | **92.8%** |

**體染色體 Precision 從 87.7% → 96.4%**（+8.7%），因為 chrX 上 146 Mb 是 FP 的主要來源。

---

## 3. 視覺化

### Fig01 — 染色體 Ideogram（四類顏色編碼）

![LOH Classification Ideogram](../figures/seqc2_loh_ideogram.png)

全基因組染色體視圖，綠色=TP（兩系統一致 LOH）、紅色=FP（TO-only）、橘色=FN（SEQC2-only）、灰色=TN。

**關鍵觀察**：
- chr8 幾乎全綠（~96% LOH，兩系統高度一致）
- chr3、chr5、chr11 大片段綠色
- chrX 幾乎全紅 — **女性細胞株半合子 X，LongPhase-TO 判為 LOH 但 SEQC2 未標記**
- FN（橘色）零星分佈，多為小片段

### Fig02 — 逐染色體 Stacked Bar

![LOH Classification per Chromosome](../figures/seqc2_loh_per_chr_stacked.png)

每條染色體的 TP/FP/FN/TN 佔比。chr8 LOH 佔 ~96%（幾乎全 TP），chrX LOH 佔 ~93%（幾乎全 FP）。

### Fig03 — Confusion Matrix（bp-Level + Variant-Level）

![LOH Confusion Matrix](../figures/seqc2_loh_confusion_matrix.png)

左圖：bp-level confusion matrix（Mb 單位）。右圖：ISM variant-level confusion matrix。

**ISM Variant-Level**：
- 17,186 variants 落在兩系統一致 LOH 區（TP zone）
- 21,749 variants 落在兩系統一致非 LOH 區（TN zone）
- 369 variants 在 TO LOH 但非 SEQC2 LOH（FP zone）
- 792 variants 在 SEQC2 LOH 但非 TO LOH（FN zone）

### Fig04 — 準確性指標摘要

![LOH Accuracy Summary](../figures/seqc2_loh_accuracy_summary.png)

Sensitivity/Precision/Jaccard 的含 chrX vs 體染色體比較。體染色體 Precision 顯著高於含 chrX 版本。

---

## 4. FP 分析（TO-only LOH）

### 4.1 FP 來源分解

| 來源 | Mb | 佔 FP % | 說明 |
|------|-----|---------|------|
| **chrX** | 146 | **73%** | 女性細胞株半合子 X，預期行為 |
| 體染色體碎片 | 54 | 27% | 邊界碎片化 + 少量 self-phasing |

### 4.2 chrX 說明

HCC1395 是女性三陰性乳癌細胞株。chrX 在腫瘤中一條 X 常失活或丟失，LongPhase-TO 的 self-phasing 在半合子區域判定為 LOH（因為只有一個 haplotype）。這是**預期行為而非 artifact**。SEQC2 benchmark 未將 chrX 半合子標記為 LOH，造成分類不一致。

### 4.3 體染色體 FP

54 Mb 體染色體 FP（524 regions 中排除 chrX 後的碎片），主要是：
- LOH 邊界區域的碎片化差異（region 解析度不同）
- 少量 self-phasing 造成的假 LOH（佔比極小）

---

## 5. FN 分析（SEQC2-only LOH）

| 特性 | 數值 |
|------|------|
| Regions | 148 |
| Total | 58 Mb |
| 平均大小 | 392 Kb |

FN 區域多為小片段，分佈於 chr4（17 Mb）、chr9（12 Mb）等。可能原因：
1. LongPhase-TO self-phasing 在這些區域有足夠 heterozygous SNPs → 未判為 LOH
2. SEQC2 benchmark 含 GAINLOH（CN=3 + LOH），LongPhase-TO 可能因 copy gain 未判為 LOH

---

## 6. ISM Variant 在四類區域的分佈

| LOH Zone | n variants | LOH rate (ISM) | VerificationClass=Noise % |
|----------|-----------|----------------|--------------------------|
| TP (both LOH) | 17,186 | ~54% | ~54% |
| FP (TO-only) | 369 | ~42% | — |
| FN (SEQC2-only) | 792 | ~28% | — |
| TN (neither) | 21,749 | ~22% | ~18% |

**TP zone 的 54% Noise VerificationClass** 確認：即使 LOH.bed 準確判定 LOH 區域，ISM 的表觀遺傳分析在 LOH 區仍然大幅失效（因為只有單一 haplotype → 無法做 HP-based 比較）。

---

## 7. 核心結論修正

### 7.1 先前假設（已修正）

| 假設 | 修正前 | 修正後 |
|------|--------|--------|
| LOH.bed 受 self-phasing 嚴重汙染 | LOH.bed 覆蓋 53.9% → 認為過量 | **體染色體 Precision=96.4%，只有 3.6% 假 LOH** |
| LOH.bed 不可信 | 認為需大幅修正 | **LOH.bed 高度準確，F1=96.2%** |
| TO LOH 過量來自 self-phasing | 全歸因 self-phasing | **73% TO-only 來自 chrX 半合子，非 self-phasing** |

### 7.2 關鍵推論

1. **LOH.bed 可直接用於後續分析** — 體染色體準確率 >96%
2. **ISM 的 LOH 問題不在判定準確性** — 而在 LOH 區域內 ISM epigenetic analysis 的有效性
3. **chrX 應從 LOH 分析中排除** — 半合子 X 不是真正的 somatic LOH
4. **先前 read threshold 報告的「雙問題分解」需重新解讀** — self-phasing scaffold 的影響可能主要集中在 ISM site-level HP_Ratio（而非 LOH.bed region-level）

---

## 8. 輸出檔案

### 8.1 分類 BED 檔

| 檔案 | Regions | Mb | 路徑 |
|------|---------|-----|------|
| TP | 1,066 | 1,432 | `data/seqc2_vs_to_validation/loh_classified_TP.bed` |
| FP | 524 | 200 | `data/seqc2_vs_to_validation/loh_classified_FP.bed` |
| FN | 148 | 58 | `data/seqc2_vs_to_validation/loh_classified_FN.bed` |

### 8.2 圖表檔

| 圖 | 檔案 |
|----|------|
| Fig01 | `figures/seqc2_loh_ideogram.png` |
| Fig02 | `figures/seqc2_loh_per_chr_stacked.png` |
| Fig03 | `figures/seqc2_loh_confusion_matrix.png` |
| Fig04 | `figures/seqc2_loh_accuracy_summary.png` |

### 8.3 數據檔

| 檔案 | 說明 |
|------|------|
| `data/seqc2_vs_to_validation/hcc1395_variant_loh_zone.tsv` | ISM variants with LOH zone 分類 |
| `data/seqc2_vs_to_validation/seqc2_loh_only.bed` | SEQC2 LOH regions (315) |
| `data/seqc2_vs_to_validation/to_5k_autosomal.bed` | TO LOH excluding chrX/Y |

---

## 9. 前置報告

- [SEQC2 CNV Benchmark 論文整理](20260402_seqc2_cnv_benchmark_analysis.md)
- [SEQC2 vs LongPhase-TO 交集驗證](20260402_seqc2_vs_longphase_to_loh_validation.md)
- [LOH Read Threshold 視覺化論證](../../docs/reports/validated/2026/04/20260402_loh_read_threshold_visual_argument_01.md)
- [Self-Phasing 因果鏈報告](../../docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md)
