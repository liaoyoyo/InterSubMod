<!--
建立時間: 2026-04-02 17:30
目標: SEQC2 CNV Benchmark LOH 與 LongPhase-TO LOH.bed 的交集驗證分析
處理範圍: HCC1395 (5kHz + Dorado) LOH.bed vs SEQC2 high-confidence LOH benchmark
關聯檔案:
  - research/loh_investigation/data/seqc2_cnv_benchmark_v4/
  - research/loh_investigation/data/seqc2_vs_to_validation/
  - docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md
-->

# SEQC2 vs LongPhase-TO LOH 交集驗證報告

## 1. 驗證目標

LongPhase-TO 輸出的 `LOH.bed` 覆蓋 HCC1395 基因組 53.9%，遠超一般預期。利用 SEQC2 CNV Benchmark（6 caller + 3 orthogonal technology 建立的 high-confidence LOH regions）作為 **orthogonal ground truth**，驗證 LongPhase-TO LOH.bed 的準確性。

## 2. 數據來源

| 數據集 | 來源 | 方法 | Genome Build |
|--------|------|------|-------------|
| **SEQC2 LOH Benchmark** | Zenodo v4 (DOI: 10.5281/zenodo.14619054) | 6 CNV callers × 21 WGS replicates × 3 orthogonal validations | GRCh38 |
| **LongPhase-TO LOH.bed (5kHz)** | 20260315_hcc1395_to_pilot | TO phased genotype ratio, `--loh` flag | GRCh38 |
| **LongPhase-TO LOH.bed (Dorado)** | 20260315_hcc1395_dorado_to_pilot | 同上，不同 basecaller | GRCh38 |
| **ISM Variants** | all_region_rows.tsv.gz (master dataset) | HCC1395 TO mode, 40,096 regions | GRCh38 |

## 3. 基因組層級交集分析

### 3.1 全基因組（含 chrX）

| 指標 | 數值 |
|------|------|
| SEQC2 LOH | 1,490.4 Mb (315 regions, 49.2% genome) |
| TO LOH (5kHz) | 1,632.2 Mb (1,094 regions, 53.8% genome) |
| **Overlap（兩者一致）** | **1,432.0 Mb (47.2% genome)** |
| TO-only（TO 多判） | 200.2 Mb (6.6% genome) |
| SEQC2-only（TO 漏判） | 58.4 Mb (1.9% genome) |
| **Jaccard Index** | **0.847** |

### 3.2 Autosomal-only（排除 chrX/Y — 核心數字）

chrX 上 LongPhase-TO 判了 146.1 Mb LOH，但 SEQC2 benchmark **不含 chrX**（HCC1395 是女性細胞株，chrX hemizygous 不算 LOH）。排除 chrX 後：

| 指標 | 數值 |
|------|------|
| TO autosomal LOH | 1,485.2 Mb |
| SEQC2 LOH | 1,490.4 Mb |
| **Overlap** | **1,432.0 Mb** |
| **Jaccard Index** | **0.928** |
| **Sensitivity（TO 捕捉 SEQC2）** | **96.1%** |
| **Precision（TO 被 SEQC2 確認）** | **96.4%** |
| **F1** | **96.2%** |
| TO-only（autosomal 假 LOH） | 53.2 Mb (3.6%) |
| SEQC2-missed（TO 漏判） | 58.4 Mb (3.9%) |

### 3.3 核心結論

> **LongPhase-TO 的 autosomal LOH.bed 與 SEQC2 benchmark 高度一致（Jaccard=0.928, F1=96.2%）。**
>
> **TO 多出的 LOH 主要來自 chrX (146 Mb)，不是 self-phasing artifact。**
> Autosomal 上僅多判 53.2 Mb (3.6%)，漏判 58.4 Mb (3.9%)。

### 3.4 Dorado 重複驗證

| 指標 | 5kHz | Dorado |
|------|------|--------|
| TO LOH total | 1,632.2 Mb | 1,633.9 Mb |
| Jaccard vs SEQC2 | 0.847 | 0.845 |
| Sensitivity | 96.1% | 96.0% |
| Precision | 87.7% | 87.6% |

**兩個獨立定序 run 結果幾乎相同** — LOH.bed 再現性極高。

---

## 4. Per-Chromosome 比較

| Chr | SEQC2 (Mb) | TO (Mb) | TO-only (Mb) | SEQC2-miss (Mb) | TO-only% |
|-----|-----------|---------|-------------|-----------------|----------|
| chr1 | 112.9 | 112.0 | 4.1 | 4.9 | 3.6% |
| chr2 | 122.7 | 130.5 | 8.4 | 0.6 | 6.5% |
| chr3 | 157.7 | 153.1 | 2.1 | 6.8 | 1.4% |
| chr4 | 75.8 | 74.0 | 1.4 | 3.2 | 1.9% |
| chr5 | 130.2 | 130.5 | 1.4 | 1.1 | 1.1% |
| chr6 | 76.8 | 80.3 | 6.5 | 3.0 | 8.0% |
| chr7 | 25.9 | 26.4 | 0.5 | 0.0 | 1.9% |
| chr8 | 139.3 | 122.1 | 1.4 | 18.6 | 1.1% |
| chr9 | 45.2 | 45.7 | 1.4 | 0.9 | 3.1% |
| chr10 | 104.7 | 107.0 | 2.9 | 0.6 | 2.7% |
| chr11 | 114.6 | 112.8 | 1.2 | 3.1 | 1.1% |
| chr12 | 102.6 | 102.5 | 0.4 | 0.5 | 0.4% |
| chr13 | 56.1 | 63.8 | 9.4 | 1.7 | 14.7% |
| chr14 | 34.3 | 34.6 | 1.3 | 1.0 | 3.8% |
| chr15 | 33.0 | 33.2 | 1.0 | 0.8 | 3.0% |
| chr17 | 77.4 | 73.1 | 1.0 | 5.2 | 1.3% |
| chr18 | 29.0 | 31.9 | 3.2 | 0.4 | 10.0% |
| chr19 | 32.0 | 26.8 | 0.1 | 5.3 | 0.3% |
| chr22 | 19.2 | 22.2 | 3.6 | 0.7 | 16.3% |
| **chrX** | **0.0** | **146.1** | **146.1** | **0.0** | **100%** |

**觀察**：
- 大部分 autosomal chromosomes TO-only < 5%（邊界微調差異）
- **chr13 (14.7%)** 和 **chr22 (16.3%)** 有較多 TO-only LOH — 可能是 self-phasing 產生的局部假 LOH
- **chr8 SEQC2-miss 最大 (18.6 Mb)** — TO 漏判了 chr8 部分 LOH 區域
- **chrX 100% TO-only** — SEQC2 不含 chrX LOH（女性細胞株 hemizygous X 正常）

---

## 5. ISM Variant 層級驗證

將 HCC1395 TO mode 的 40,096 個 ISM regions 映射到 LOH 區域分類：

### 5.1 五區域分類定義

| Zone | 定義 | 含義 |
|------|------|------|
| **both_LOH** | 落在 SEQC2 LOH ∩ TO LOH.bed | 確認的真 LOH 區域 |
| **TO_only_in_gain_loss** | 落在 TO LOH.bed 但在 SEQC2 gain/loss 內 | TO 誤判為 LOH（實為 gain/loss） |
| **TO_only_novel** | 落在 TO LOH.bed 但 SEQC2 無任何標記 | TO 獨有的 LOH（可能假陽性） |
| **SEQC2_only** | 落在 SEQC2 LOH 但不在 TO LOH.bed | TO 漏判的真 LOH |
| **neither** | 兩者都不是 LOH | 非 LOH 區域 |

### 5.2 TP/FP 分佈

| Zone | n | TP | FP | Precision | FP Rate |
|------|---|----|----|-----------|---------|
| **both_LOH** | 17,186 | 12,710 | 4,476 | **0.740** | 0.260 |
| TO_only_in_gain_loss | 292 | 220 | 72 | 0.753 | 0.247 |
| TO_only_novel | 77 | 61 | 16 | 0.792 | 0.208 |
| SEQC2_only | 792 | 600 | 192 | 0.758 | 0.242 |
| **neither** | 21,749 | 14,904 | 6,845 | **0.685** | 0.315 |

### 5.3 關鍵發現

**1. LOH 區域的 Precision 反而較高（0.740 vs 0.685）**

與直覺相反：確認的 LOH 區域內 TP 比例 (74.0%) **高於**非 LOH 區域 (68.5%)。原因：LOH 區域內 germline variants 已被 PON 高效過濾，剩下的多為 somatic。

**2. TO-only LOH 影響極小**

- TO_only_in_gain_loss: 僅 292 個 variants（佔 0.7%）
- TO_only_novel: 僅 77 個 variants（佔 0.2%）
- 合計 369 個 variants（0.9%），且 Precision 不差（0.753-0.792）

**3. SEQC2-only（TO 漏判）也少**

- 792 個 variants（2.0%），Precision 0.758
- 這些 variant 在 TO 模式下**仍能被分析**，只是不在 LOH.bed 內

**4. chrX 不影響 variant 分析**

ISM 的 40,096 個 variants 中 0 個在 chrX — somatic variant calling pipeline 排除了 chrX。chrX 的 146 Mb TO-only LOH 完全不影響下游分析。

### 5.4 ISM Potential_LOH 與 LOH Zone 的對應

| Zone | ISM Potential_LOH Rate |
|------|----------------------|
| both_LOH | **99.0%** |
| TO_only_in_gain_loss | 97.3% |
| TO_only_novel | 97.4% |
| SEQC2_only | **86.7%** |
| neither | 25.4% |

- **both_LOH 區域 99% 被 ISM 判為 LOH** — ISM HP-Ratio 與 LOH.bed 高度一致
- **SEQC2_only 區域 86.7% 也被 ISM 判為 LOH** — ISM 其實有捕捉到，只是 LOH.bed 邊界沒覆蓋
- **neither 區域仍有 25.4% ISM LOH** — 這部分是 ISM 的 site-level 偏差（局部 HP 不平衡但不在結構性 LOH 區域）

### 5.5 VerificationClass 分佈

| Zone | #1 | #2 | #3 |
|------|-----|-----|-----|
| both_LOH | Noise 54% | Weak 20% | Strong 15% |
| SEQC2_only | Noise 43% | Weak 26% | Strong 23% |
| neither | Weak 45% | Strong 28% | Noise 26% |

**LOH 區域 Noise 佔比極高（54%）**：LOH 破壞了 haplotype 結構，ISM 的 epigenetic analysis 在 LOH 區域失效 → 大量 Noise classification。

非 LOH 區域 Strong 佔比 28% — 只有非 LOH 區域才能產生有意義的 epigenetic subclonal signal。

---

## 6. 綜合結論

### 6.1 LongPhase-TO LOH.bed 準確性評估

| 評估項目 | 結果 | 判定 |
|----------|------|------|
| Autosomal Jaccard | **0.928** | **極高** |
| Autosomal Sensitivity | **96.1%** | 幾乎完全捕捉 SEQC2 LOH |
| Autosomal Precision | **96.4%** | 幾乎無假 LOH |
| Autosomal F1 | **96.2%** | 雙方高度一致 |
| 跨 run 再現性 | 5kHz ≈ Dorado (Jaccard 差 < 0.002) | 極穩定 |

### 6.2 先前結論的修正

| 先前推論 | 修正後 |
|---------|--------|
| "LongPhase-TO LOH.bed 覆蓋 53.9% 異常偏高" | **HCC1395 確實有 49.2% autosomal LOH（SEQC2 確認）。TO 多出的 4.6% 主要是 chrX (hemizygous)** |
| "LOH.bed 大量被 self-phasing 汙染" | **Autosomal LOH.bed 與 orthogonal benchmark F1=96.2%，self-phasing 對 LOH.bed 邊界影響極小 (~3.6%)** |
| "兩套系統都被 self-phasing 汙染" | **LOH.bed 本身高度準確；ISM Potential_LOH 的過量 (+15%) 才是 self-phasing 主要影響** |

### 6.3 重要推論

1. **LOH.bed 是可信的 LOH ground truth** — 對於 HCC1395，autosomal F1=96.2%，可直接使用
2. **ISM Potential_LOH 過度判定的 15.2% 差距**：41.8% ISM LOH vs 26.7% LOH.bed hit → 多出的 15.1% 來自 site-level HP ratio noise，不是結構性 LOH
3. **LOH 區域內 Noise 佔 54%** — ISM 在真 LOH 區域幾乎完全失效，應考慮排除 LOH 區域再做 epigenetic analysis
4. **非 LOH 區域才是 ISM 的有效工作範圍** — Strong 28%, Weak 45%，有意義的 subclonal signal 集中在此
5. **chrX LOH 不影響 ISM 分析** — somatic calling pipeline 不包含 chrX variants

### 6.4 後續可行方向

| 方向 | 可行性 | 說明 |
|------|--------|------|
| 使用 SEQC2 LOH BED 排除 LOH 區域後重算 F1 | 高 | 預期 Precision 改善（移除 LOH 區域的 Noise TP） |
| 擴展到其他 6 個樣本 | 中 | SEQC2 benchmark 僅 HCC1395，其他樣本需其他 CNV truth |
| 用 SEQC2 gain/loss BED 做 copy number aware analysis | 高 | 區分 LOH、gain、loss 對 ISM 的不同影響 |

---

## 7. 數據檔案

```
research/loh_investigation/data/seqc2_vs_to_validation/
├── seqc2_loh_only.bed           # SEQC2 LOH regions (315)
├── seqc2_non_loh.bed            # SEQC2 gain+loss regions (373)
├── to_5k_loh.bed                # LongPhase-TO LOH (HCC1395 5kHz, 1094)
├── to_dorado_loh.bed            # LongPhase-TO LOH (HCC1395 Dorado, 1208)
├── to_5k_autosomal.bed          # TO LOH autosomal only
├── overlap_to5k_seqc2.bed       # Intersection regions
├── to5k_only_loh.bed            # TO-only LOH (subtract)
├── seqc2_only_loh.bed           # SEQC2-only LOH (subtract)
├── to5k_in_seqc2_nonloh.bed     # TO LOH overlapping SEQC2 gain/loss
└── hcc1395_variant_loh_zone.tsv  # ISM variants with LOH zone classification
```
