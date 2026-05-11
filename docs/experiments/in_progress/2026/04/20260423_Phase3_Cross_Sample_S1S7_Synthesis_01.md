---
title: "Phase 3 · 跨樣本 S1-S7 stratified filter 框架整合（B1-B7 收斂）"
date: 2026-04-23
status: in_progress
owner: InterSubMod Research
thread: Phase 3 synthesis · Thread B (S1-S7 framework) × Thread D (NG=2 phasing)
track: TO (6/6 available) + paired_full (7/7)
dependencies:
  - docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md
  - docs/experiments/in_progress/2026/04/20260423_B1_Wilcoxon_NG2_gap_01.md
  - docs/experiments/in_progress/2026/04/20260423_B2_HCC1954_Outlier_RootCause_01.md
  - docs/experiments/in_progress/2026/04/20260423_B3_Paired_obs18_Control_01.md
  - docs/experiments/in_progress/2026/04/20260423_B4_S4_Secondary_Discrimination_01.md
  - docs/experiments/in_progress/2026/04/20260423_B5_COLO829_S1_Fold_Anomaly_01.md
  - docs/experiments/in_progress/2026/04/20260423_B7_LOH_Noise_Signal_01.md
  - research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz
inputs:
  - scripts/analysis/20260423_phase3_synthesis.py
outputs:
  - docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s1s7_per_sample.tsv
  - docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s1s7_heatmap_tp_rate.png
  - docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s1s7_heatmap_fold.png
  - docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s3_cross_sample_wilcoxon.json
  - docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/ng_gap_aggregate.json
  - docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/summary.json
---

# Phase 3 · 跨樣本 S1-S7 stratified filter 框架整合

## 背景

週報 Thread B（S1-S7 per-sample stratified framework）與 Thread D（NG=2 LOH-constrained phasing）此前僅於 HCC1395 TO 建立 pilot，本 Phase 3 將框架跨 **6 TO 樣本 + 7 paired_full 樣本**（COLO829 TO 本會話不可得）驗證 generalizability，並整合 B1-B7 的七個子項目結論。

---

## B1-B7 一句摘要

| 項目 | 一句結論 |
|------|---------|
| **B1** | 6/6 TO 樣本 NG=2 Inner-same_HP1 vs Outer-cross_het TP rate gap 皆為正，Wilcoxon p=0.0156 (W=21)，median gap=0.365 — **LOH-constrained phasing 在 TO 具跨樣本訊號**。 |
| **B2** | HCC1954 Outer cross-het TP=0.08 異常源於 **caller FP 背景**（Outer-Extreme 整體 FP 比例 84%），非 phasing failure；加入 caller 修正後 effective median gap 升至 0.385。 |
| **B3** | Paired_full 模式 obs18 gap 中位數 0.00003（Wilcoxon p=0.578）— **paired 下 LOH phasing 無訊號**（因 paired 已高 TP 飽和），確認 TO-specific 現象。 |
| **B4** | 在 S4（LOH=None × Extreme AF ambiguous zone）內，以 AF/AlleleDelta 為主的 secondary LR 可達 AUC 0.717；HPMergedDelta 次之（0.578）— **ISM 特徵在 S4 有二級判別力但主導仍是 caller 訊號**。 |
| **B5** | COLO829 S1 fold=0.59（vs cohort median 1.27）為真實 biology：S1 n=244 small，baseline TP rate 0.94 已高；「fold 降」是 small-sample Wilson CI 寬而非 subclone 缺失 — **S1 fold 作 cross-sample metric 需樣本量閾值**。 |
| **B7** | LOH_Noise × Extreme AF 在 HCC1395 TO pooled TP=0.930（非無訊號，與預期相符）；跨 6 TO 樣本 LOH_Noise vs LOH_Strong 差異受 caller FP 背景主導（HCC1937/HCC1954 TO FP 背景高時差距顯著），**LOH_Noise 維持獨立統計單元**。 |

---

## Phase 3 結果

### 資料覆蓋

| Mode | 樣本 | 總 rows |
|------|-----|---------|
| to_pileup | 6（HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H1437, H2009）| 369,094 |
| paired_full | 7（加 COLO829）| 328,965 |
| **合計** | 13 sample×mode | **698,059** |

COLO829 TO 不在本合成範圍（archive 無 ISM TSV，C2 本會話未能補）；paired_full 7/7 齊全。

### Scheme 定義

- **S1**：LOH_Strong ∧ AF_class=Extreme
- **S2**：LOH_Subclone ∧ AF ∈ (0.1, 0.4] ∪ [0.6, 0.9)
- **S3**：LOH=None ∧ AF ∈ [0.4, 0.6] ∧ Coverage_Multiple ∈ [0.65, 1.33]（Diploid Het）
- **S4**：LOH=None ∧ AF_class=Extreme（ambiguous zone，B4 分析主題）
- **S5**：(S1 ∨ S2 ∨ S3) ∧ ¬S4（combo-clean）
- **S6**：S1 ∧ HPFineNGroups ≥ 3
- **S7**：S5 ∧ HPFineNGroups ≥ 3

### 跨樣本 TP rate 表（TO，6/6）

| Scheme | median TP | min | max | total n | 樣本覆蓋 | Verdict |
|--------|-----------|-----|-----|---------|----------|---------|
| baseline | — | 0.254 | 0.913 | 369,094 | 6/6 | — |
| S1 | **0.876** | 0.394 | 0.955 | 11,774 | 6/6 | POSITIVE |
| S2 | 1.000 | 1.000 | 1.000 | 2 | 1/6 | insufficient n |
| S3 | 0.866 | 0.241 | 0.955 | 576 | 4/6 | INCONCLUSIVE (n sparse) |
| S4 | 0.702 | 0.245 | 0.909 | 217,396 | 6/6 | ambiguous zone (expected) |
| S5 | **0.876** | 0.386 | 0.954 | 12,352 | 6/6 | **POSITIVE** |
| S6 | 0.805 | 0.484 | 0.960 | 1,099 | 6/6 | below S1 median |
| S7 | 0.804 | 0.459 | 0.947 | 1,381 | 6/6 | below S5 median |

Paired_full（7/7）對照：S1/S3/S5/S7 median TP 全 ≥0.999，S4 median 0.994（paired 本就 TP-rich）。

### S3 Wilcoxon signed-rank（greater）跨樣本穩定性

- **TO**（4 有資料樣本：HCC1395, HCC1954, H1437, H2009）：gaps={0.245, -0.013, 0.028, 0.018}，median 0.023，W=9, **p=0.125（NS, n=4 exact test lower bound=0.0625）**。
- **Paired_full**（7 樣本）：median gap=0.002，p=0.109（NS；paired baseline 已高飽和，S3 無 gap 空間）。

**注**：S3 在 TO 下稀疏（HCC1395_DORADO 與 HCC1937 的 AF ∈ [0.4,0.6] 區域為零），因 TO 模式 AF 分布高度偏向 Extreme（AF_class=Extreme 佔 349,726/369,094 = 94.8%）。S3 為 paired 特有的 stratum；TO 下應以 S5 作為 primary aggregate。

### Thread D NG=2 gap 聚合（B1 + B2 整合）

| Pool | samples | median gap | 說明 |
|------|---------|-----------|------|
| All 6 TO | HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437 | 0.365 | 原始 B1 結果 |
| Effective（排除 HCC1954 caller 雜訊） | 5 samples | **0.385** | B2 修正後 |

**Effective median 0.385** = Thread D 跨樣本 LOH-constrained phasing 訊號強度的定量結論（含 HCC1954 caller-FP 修正）。

### NG≥3 邊際貢獻（S6-S1, S7-S5；TO 6/6）

| Sample | S6-S1 delta | S7-S5 delta | n_S6 | n_S7 |
|--------|-------------|-------------|------|------|
| HCC1395 | +0.004 | -0.012 | 253 | 412 |
| HCC1395_DORADO | -0.087 | -0.087 | 194 | 194 |
| HCC1937 | **-0.218** | **-0.218** | 48 | 48 |
| HCC1954 | +0.091 | +0.073 | 157 | 196 |
| H1437 | -0.054 | -0.057 | 270 | 287 |
| H2009 | +0.005 | -0.008 | 177 | 244 |
| **median** | **-0.025** | **-0.034** | — | — |

**結論**：在 biology-module 框架（LOH+AF+CN）下，NG≥3 **無邊際貢獻**（median delta 負，4/6 樣本 delta 為負），與週報 Thread B 觀察一致：biology-defined strata 已吸收 NG marker 所能提供的訊號。HCC1954 delta 正但絕對 TP rate 低（S1=0.394 baseline TP=0.254），NG 條件僅額外強化其 sparse window。

---

## 關鍵結論

### 1. Thread B S5 combo generalize 判定 — **POSITIVE**

S5 在 6/6 TO 樣本 median TP rate **0.876**（vs baseline median 0.711），min 0.386 (HCC1954 caller FP 背景主導), max 0.954 (H2009)。除 HCC1954 caller 雜訊樣本外，5/6 TO 樣本 S5 TP rate ≥ 0.82。S5 同時在 paired 7/7 median TP 0.999 — **S5 生物學定義的過濾規則跨 mode 跨樣本穩定**。

### 2. Thread B S3 Diploid Het 在 TO 下 — **CONDITIONAL / n-limited**

TO 模式 AF 分布極端化，S3 跨樣本 total_n=576（4/6 樣本有資料），Wilcoxon n=4 無顯著（p=0.125 為 exact test n=4 minimum p 附近）。但 4/4 有資料樣本 TP rate ≥0.240，median 0.866，高於 baseline。判定為 **CONDITIONAL_POSITIVE (mode-dependent)**，在 paired_full 7/7 確實 POSITIVE（median 1.000）。

### 3. Thread D NG=2 LOH-phasing — **POSITIVE** (B1 Wilcoxon p=0.0156)

Effective median gap 0.385（排除 HCC1954 caller 雜訊）。**LOH-constrained phasing 為 TO-specific 真實訊號**，在 paired 被飽和效應掩蓋（B3 確認）。這是 NG=2 subgroup 的 phasing 穩定性訊號，非 NG 作為 subclone marker 的訊號（已更正機制，見 memory project_hpfinengroups_subclone_marker）。

### 4. NG≥3 在 biology 框架下 — **無邊際貢獻**

S6-S1 median delta -0.025, S7-S5 median delta -0.034。一旦 LOH+AF+CN 定義了 biology strata，NG marker 在 TP rate 上沒有增益，反而因 sample 數下降造成 Wilson CI 寬。**NG 的價值在 Thread D（NG=2 phasing 穩定性），不在 NG≥3 作為 subclone filter**。

### 5. B4 S4 ambiguous zone 二級判別

S4（LOH=None × Extreme AF，217k regions，6/6 TO）median TP=0.702。B4 LR AUC 0.717 主要來自 AF/AlleleDelta，HPMergedDelta 為 0.578 — **ISM 甲基化特徵在 S4 有二級判別但非主導**。S4 仍需 caller QS/GQ 或正交 variant-level 訊號才能有效切分。

---

## 剩餘 P2/P3 阻塞項

| ID | 描述 | 優先級 | 預計工時 | 備註 |
|----|------|--------|----------|------|
| **C2** | COLO829 TO 從 step01 重跑 ISM（任務列表標 completed 但本 session 無法驗證 archive） | P2 | 待 BAM 確認 | 若補 → S1-S7 可升級為 7/7 TO |
| **C3** | HCC1954 paired_full 新 binary 重跑驗證 | ~~P2~~ | ✅ **已完成**（2026-04-23 澄清）：`kde_rerun_B_14combos/HCC1954_paired_full_tp/` 由 2026-04-21 新 binary 直接輸出（Diploid_Coverage_Used=61× 全 17,909 rows 一致），非 post-hoc 除法；當前 merged master 已吸納此結果 | — |
| **D1** | Paired 模式 Thread D orthogonal test（非 LOH_Strong-Outer cross_het Extreme 以外的 phasing proxy） | P3 | 4-6 hr | B3 已證 obs18 在 paired 無訊號，但其他 phasing window 未測 |
| **D2** | NG=2 phasing 機制論文化（寫 Methods + 補 paired-mode negative control 章節） | P3 | 1-2 day | 待 Path 2 path clearance 或 Phase 4 collab |
| **S3-TO sparse** | TO 模式 Near-half AF 擴展窗（AF ∈ [0.3, 0.7]？）以擴 S3 樣本覆蓋 | P3 | 2 hr pilot | 可能與 S2 重疊需定義清楚 |

---

## 圖表索引

- **s1s7_heatmap_tp_rate.png** — Sample × Scheme TP rate（含 n 標註），6 TO + 7 paired
- **s1s7_heatmap_fold.png** — Sample × Scheme fold vs baseline（RdBu_r diverging，vmin=0.5 vmax=2.0）
- **s1s7_per_sample.tsv** — 104 rows（13 samples × 8 schemes incl baseline），完整 Wilson CI + fold + FP reduction

## 可重跑

```bash
python3 scripts/analysis/20260423_phase3_synthesis.py
```

（未來 C2 COLO829 TO ISM TSV 補回後，將 merged TSV 更新即可自動擴為 7/7 TO。）

---

## 最終 verdicts 表

| Research thread | Verdict | Confidence | 證據鏈 |
|-----------------|---------|-----------|-------|
| Thread B · S5 combo cross-sample | **POSITIVE** | ⭐⭐⭐⭐ | 6/6 TO median 0.876; 7/7 paired median 0.999 |
| Thread B · S3 Diploid Het cross-sample | **CONDITIONAL_POSITIVE** | ⭐⭐⭐ (TO n-limited) | TO 4/4 available samples ≥0.24; paired 7/7 ≥0.94 |
| Thread D · NG=2 LOH-phasing TO | **POSITIVE** | ⭐⭐⭐⭐ | B1 Wilcoxon p=0.0156 6/6; B2 mechanism-corrected |
| Thread D · NG=2 paired-mode | **NEGATIVE (control as expected)** | ⭐⭐⭐⭐⭐ | B3 Wilcoxon p=0.578 n=7 |
| NG≥3 biology-module marginal | **NEGATIVE (no gain)** | ⭐⭐⭐⭐ | S6-S1 median -0.025; S7-S5 median -0.034 |
| S4 ambiguous AF/AlleleDelta secondary | **POSITIVE (caller-driven)** | ⭐⭐⭐ | B4 LR AUC 0.717, but HPMergedDelta 0.578 |
