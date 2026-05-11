---
title: "B1 Wilcoxon Signed-Rank — NG=2 Inner same_HP1 vs Outer cross_het TP Gap"
date: 2026-04-23
status: in_progress
verdict: POSITIVE
phase: Thread D follow-up
track: TO
samples:
  - HCC1395
  - HCC1395_DORADO
  - HCC1937
  - HCC1954
  - H2009
  - H1437
hypothesis_id: LOH-constrained-phasing-D
source_data: research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv
artifacts:
  - scripts/analysis/20260423_B1_wilcoxon_ng2_gap.py
  - research/tpfp_loh_af_kde_discrimination/data/obs18_wilcoxon_B1.json
tags:
  - wilcoxon
  - LOH
  - NG2
  - same-hap
  - TO-mode
  - statistical-validation
---

# B1 Wilcoxon Signed-Rank — NG=2 Inner same_HP1 vs Outer cross_het TP Gap

## 目標

對 Thread D 觀察（Obs18：NG=2 在 Inner 93-99% same-hap vs Outer cross-het 飽和）做正式
Wilcoxon signed-rank 統計檢定，確認 6 TO 樣本 Inner-Outer TP gap 不是雜訊。

## 方法

- **資料來源**：`research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv`
  （56 rows，已存在）
- **比對組合**：Inner `same_HP1 (HP1 + HP1-1)` 的 `tp_rate` 減 Outer `cross_het (HP1 + HP2-1)` 的 `tp_rate`
- **統計**：`scipy.stats.wilcoxon(gaps, alternative='greater', zero_method='wilcox', method='exact')`
- **CI**：Bootstrap 1000 resamples（seed=20260423）on median gap，取 2.5/97.5 percentile
- **樣本數**：n=6（HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437）

## 結果

### 每樣本 Gap（tsv 原始值）

| Sample | Inner same_HP1 TP rate (n) | Outer cross_het TP rate (n) | Gap |
|---|---|---|---|
| HCC1395 | 0.9589 (219) | 0.5000 (2) | **+0.4589** |
| HCC1395_DORADO | 0.9385 (7220) | 0.5531 (3305) | **+0.3854** |
| HCC1937 | 0.7586 (6901) | 0.2364 (880) | **+0.5222** |
| HCC1954 | 0.4289 (3402) | 0.0842 (4050) | **+0.3447** |
| H2009 | 0.9318 (27086) | 0.8824 (3785) | +0.0494 |
| H1437 | 0.9196 (4614) | 0.6884 (4439) | +0.2311 |

6/6 樣本 gap > 0（同方向），所有樣本 Inner same_HP1 TP rate > Outer cross_het TP rate。

### Wilcoxon signed-rank（alternative='greater'）

- **W statistic**: **21.0**（6 samples 的最大可能秩和 = 6·7/2 = 21）
- **p-value**: **0.015625**（exact method）
- **顯著性 @ α=0.05**: **TRUE**

W=21 代表所有 6 個樣本 gap 皆為正且排名一致，是 n=6 下 exact test 能給出的最小 p 值。

### Bootstrap 95% CI on median gap

- **Median gap**: 0.3650
- **Mean gap**: 0.3320
- **95% CI (percentile)**: [**0.1403**, **0.4906**]
- **Min / Max**: 0.0494 / 0.5222
- CI 下界 > 0，median gap 顯著正向

## 判定

**POSITIVE**：Inner same_HP1 TP rate 顯著高於 Outer cross_het TP rate（p=0.0156，median gap=+0.365，bootstrap 95% CI [0.14, 0.49]）。

- **與週報宣稱對照**：週報列出的 6 個 gap（+0.46, +0.39, +0.52, +0.35, +0.05, +0.23）與本腳本以 tsv 原值重算的 gap（+0.459, +0.385, +0.522, +0.345, +0.049, +0.231）吻合至小數點第二位。週報數值可信。
- **飽和樣本 H2009**：gap=+0.049，Outer cross_het 基線 0.88 已飽和（該樣本 TP rate 全域高），但仍同方向。
- **最小 n 樣本 HCC1395**：Outer cross_het n=2（極小），但 HCC1395_DORADO（同樣本 DORADO re-basecall）n=3305，+0.385 gap 穩定，可互為支持證據。

## 意涵

- NG=2 LOH-constrained phasing 的 Inner-Outer TP 差異**通過正式非參數檢定**，支持「LOH 區內 `{HP1, HP1-1}` bucket 佔優」的 Thread D 機制。
- HCC1954 基線 TP rate 全面低（~0.1），但 gap 仍 +0.345；顯示此現象不依賴於特定 TP rate 範圍。
- n=6 下 Wilcoxon exact p=0.0156 已是最極端值，進一步提升統計檢力需擴樣本或改用層次模型。

## 後續

- 若要強化，可考慮：
  - 納入 `same_HP2 (HP2 + HP2-1)` Inner vs `cross_het_inv (HP1-1 + HP2)` Outer 做互補驗證
  - 加權 Wilcoxon / mixed-effects model 以納入 per-sample variance
  - permutation test（shuffle Inner/Outer label）做 robustness check

## Artifacts

- 分析腳本：`scripts/analysis/20260423_B1_wilcoxon_ng2_gap.py`
- 結果 JSON：`research/tpfp_loh_af_kde_discrimination/data/obs18_wilcoxon_B1.json`
- 來源資料：`research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv`
