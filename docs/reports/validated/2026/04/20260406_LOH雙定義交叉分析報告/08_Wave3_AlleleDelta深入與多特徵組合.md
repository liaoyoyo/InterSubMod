<!--
建立時間: 2026-04-06
目標: Wave 3 M4 結果 — AlleleDelta 深入分析、多特徵組合、cnLOH 驗證
處理範圍: LOH 區域 Allele 維度深入、組合模型、LOH 子區域分析、confound 檢查
關聯檔案:
  - scripts/analysis/build_allele_deep_dive.py
  - big7_disk_output/synthesis/observation_workspaces/20260406_allele_deep_dive/
-->

# Wave 3: AlleleDelta 深入分析 + 多特徵組合

## M4: 分析結果

### F1: Hopeless vs Separable 子群

在 TO LOH 區域中，以 |AlleleDelta| ≤ 0.01 為閾值：

- **Hopeless 子群** (|AlleleDelta| ≤ 0.01): AlleleDelta 完全無區分力
- **Separable 子群** (|AlleleDelta| > 0.01): 有微弱信號

具體數據見 `f1_hopeless_vs_separable.tsv`。

### F3: 多特徵組合模型

使用 AlleleDelta + CramersV + PairwiseMedianDist 三特徵組合：

| 方法 | AUC | 權重 |
|------|-----|------|
| Best linear combo | 0.5601 | (0.17, 0.83, 0.0) |
| Rank-based | 0.5652 | equal |
| **Voting 2/3** | **0.5774** | equal |

**所有組合 AUC < 0.58 門檻** → 多特徵組合不可行。

### F4: LOH 子區域分析

按 Coverage_Multiple 分類 LOH 子區域：

| LOH 子類型 | 定義 | 最佳特徵 | AUC | N |
|-----------|------|---------|-----|---|
| **cnLOH** | Cov 0.8-1.2 | PairwiseMeanDist | **0.5865** | 65,852 |
| deletion_LOH | Cov < 0.8 | AlleleDelta | 0.5656 | 70,797 |
| gain_LOH | Cov > 1.2 | PairwiseMeanDist | 0.5724 | 38,893 |

cnLOH PairwiseMeanDist 表面上 AUC=0.5865 過 0.58 門檻。

### cnLOH PairwiseMeanDist 深入驗證（追加分析）

| 驗證項目 | 結果 | 判定 |
|----------|------|------|
| Overall AUC | 0.5865 | 過 0.58 |
| Per-sample 5/7 TP-favored | COLO829=0.26, H1437=0.49 FP-favored | **FAIL** |
| Per-sample Mean AUC | **0.4987** (隨機) | **FAIL** |
| NumReads residualize | 0.5864 (不變) | CLEAN |
| Double residualize | 0.5902 | CLEAN |
| AF-bin L3 | 4/4 TP-favored | PASS |

**Per-sample 詳細：**

| Sample | N | TP | FP | AUC | 方向 |
|--------|---|----|----|-----|------|
| COLO829 | 58 | 35 | 23 | 0.260 | FP-favored |
| H1437 | 7,888 | 6,197 | 1,691 | 0.492 | FP-favored |
| H2009 | 31,970 | 29,410 | 2,560 | 0.579 | TP-favored |
| HCC1395 | 6,830 | 5,106 | 1,724 | 0.530 | TP-favored |
| HCC1395_DORADO | 7,452 | 5,523 | 1,929 | 0.543 | TP-favored |
| HCC1937 | 5,047 | 2,750 | 2,297 | 0.526 | TP-favored |
| HCC1954 | 6,607 | 1,864 | 4,743 | 0.561 | TP-favored |

**判定**：整體 AUC=0.5865 是 **Simpson's Paradox** — 被 H2009 (N=31,970, 48.6%) 主導。逐樣本 Mean AUC=0.4987 本質上隨機。cnLOH PairwiseMeanDist **NOT VIABLE**。

### F5: Confound 檢查

| 特徵 | Raw AUC | Residualized (NumReads) | Residualized (Coverage) | 判定 |
|------|---------|-------------------------|------------------------|------|
| AlleleDelta | 0.5564 | 0.5560 | 0.5560 | **CLEAN** (confound-free) |
| **CramersV** | 0.5107 | **0.4644** | **0.4644** | **CONFOUNDED** (drop=0.046) |
| PairwiseMedianDist | 0.5332 | 0.5332 | 0.5332 | **CLEAN** |

**新發現**：CramersV AUC 從 0.511 降至 0.464 (residualized) — 原始 AUC 中 NumReads 貢獻了大部分信號。CramersV 在 LOH 區域的區分力是假的。

---

## 判定總結

| # | 判定 | 置信度 |
|---|------|--------|
| **J13** | **LOH 區域多特徵組合不可行** (最佳 Voting AUC=0.577 < 0.58) | **確定** |
| **J14** | **cnLOH PairwiseMeanDist 0.587 是 Simpson's Paradox** (mean per-sample AUC=0.50) | **確定** |
| **J15** | **CramersV 在 LOH 區域被 NumReads confound** (AUC 0.511→0.464) | **高** |
| **J16** | **AlleleDelta 是 LOH 內唯一真實的 confound-free 信號** (AUC=0.556, 7/7 一致) | **確定** |

### 圖表

![AlleleDelta 分佈: Hopeless vs Separable 子群](figures/w3_m4_01_hopeless_vs_separable.png)

![多特徵組合 AUC top-20 grid search](figures/w3_m4_02_combo_top20.png)

![Per-sample AlleleDelta AUC bar chart](figures/w3_m4_03_per_sample_allele_auc.png)

![LOH 子區域 AUC: cnLOH vs deletion vs gain](figures/w3_m4_04_loh_subregion_auc.png)

![Confound check: Raw vs Residualized AUC](figures/w3_m4_05_confound_check.png)
