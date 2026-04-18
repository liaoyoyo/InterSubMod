<!--
建立時間: 2026-04-06
目標: Step 1.4 SEQC2 LOH 驗證 + LOH 象限移除的 F1 影響分析
處理範圍: HCC1395 + HCC1395_DORADO, TO mode, autosomal only
關聯檔案:
  - scripts/analysis/build_seqc2_loh_validation.py
-->

# Step 1.4: SEQC2 LOH 驗證與 F1 影響

## chrX/Y 排除說明

| 項目 | 說明 |
|------|------|
| **排除範圍** | chrX, chrY |
| **原因** | SEQC2 LOH benchmark 不含 chrX/Y。HCC1395 是**女性細胞株**，chrX hemizygous ≠ LOH |
| **影響** | Master dataset HCC1395 的 chrX/Y = **0 筆**（建置時已排除）→ 無實際影響 |
| **COLO829** | 唯一有 chrX/Y 的 sample（1,870 TO rows），本分析僅限 HCC1395 |
| **全樣本分析時** | 必須排除 COLO829 chrX/Y 或標記 autosomal_only=True |

---

## B.1: LOH.bed vs SEQC2 LOH Jaccard

使用 base-pair level interval comparison（非 entry-level），將 LongPhase-TO LOH.bed 與 SEQC2 NGS benchmark 的 320 個 autosomal LOH entries (1,490.5 Mb) 比較。

| Sample | Jaccard | Sensitivity | Precision | F1 | LOH.bed (Mb) | Intersect (Mb) |
|--------|---------|-------------|-----------|-----|--------------|----------------|
| HCC1395 5kHz | **0.9277** | 0.9608 | 0.9642 | 0.9625 | 1,485.2 | 1,432.1 |
| HCC1395 Dorado | **0.9290** | 0.9602 | 0.9662 | 0.9632 | 1,481.2 | 1,431.2 |

**驗證基準**: Target 0.928 ± 0.002 → **PASS** (0.9277 在範圍內)

![SEQC2 Jaccard comparison](figures/w2s_01_seqc2_jaccard_bar.png)

**限制**:
- SEQC2 LOH benchmark 基於 NGS short-read，不是 ONT long-read
- LongPhase-TO LOH.bed 基於 phased genotype ratio，方法完全不同
- 兩者高度一致 (Jaccard=0.928) 驗證了 LOH.bed 的可靠性

---

## B.2: HP_Ratio → SEQC2 LOH 預測

使用 |HP_Ratio - 0.5| 作為 per-variant score 預測是否落在 SEQC2 LOH 區域。

| 方法 | AUC | Best F1 | Precision | Recall |
|------|-----|---------|-----------|--------|
| HP_Ratio score | **0.8979** | 0.858 (thr=0.45) | 0.767 | 0.974 |
| LOH.bed (region-level) | — | **0.966** | 0.977 | 0.956 |

![HP_Ratio ROC for SEQC2 LOH prediction](figures/w2s_03_hp_ratio_roc.png)

![HP_Ratio threshold sweep](figures/w2s_04_hp_ratio_threshold_sweep.png)

HP_Ratio per-variant 對 SEQC2 LOH 有強預測力 (AUC=0.90)，但 LOH.bed region-level 預測更好 (F1=0.97)。這合理 — region-level BED 是 phasing 程式的直接輸出，per-variant HP_Ratio 是間接指標。

---

## B.3: F1 影響分析 (S1-S5)

**基線**: HCC1395 ClairS-TO F1=0.7156 (master dataset: TP=28,495, FP=11,601, FN=11,051)

### 情境模擬結果

| 情境 | 描述 | TP Loss | FP Removal | ΔF1 | 安全 |
|------|------|---------|------------|-----|------|
| S1 | Remove Q1 (Both LOH) | **45.1%** | 38.9% | -0.2135 | FAIL |
| S2 | Remove Q2 (ISM-only) | **15.5%** | 15.5% | -0.0596 | FAIL |
| S3 | Remove Q1+Q2 (all HP Imbalance) | **60.6%** | 54.3% | -0.3151 | FAIL |
| S4 | Keep Q4 only | **61.1%** | 54.8% | -0.3185 | FAIL |
| S5 | Remove Q1 + AlleleSig=False | **25.5%** | 35.1% | -0.0940 | FAIL |

**所有 5 情境均 FAIL TP loss ≤ 2%**。HCC1395_DORADO 結果一致。

![F1 impact per scenario](figures/w2s_05_f1_impact_scenarios.png)

![TP loss vs FP removal tradeoff](figures/w2s_06_tp_fp_tradeoff.png)

### 根本原因

LOH 區域 **TP-enriched** (Q1 FP rate = 0.239，全四象限最低)。LOH 是 tumor suppressor 失活的常見機制 — somatic variants 在 LOH 區域更可能是真實的 (TP)。因此，任何基於 LOH 象限的移除策略都會嚴重損失 TP。

### 輸出檔案

- `f1_impact_scenarios.tsv` — 每個情境的 TP/FP 移除數量與 F1 變化
- `removed_variants_per_scenario.tsv.gz` — 每個被移除 variant 的完整特徵 (166,876 rows)

---

## 操作條件與限制

- Master dataset TP/FP (28,495/11,601) 與 ClairS-TO VCF 基線 (28,396/11,843) 有差異（ISM region 處理合併/排除）
- FN 使用 VCF 基線值 11,051
- F1 計算以 master dataset 內部一致性為準
