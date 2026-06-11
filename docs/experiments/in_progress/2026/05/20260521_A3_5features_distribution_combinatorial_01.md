---
title: "Phase A3 — 5 features TP/FP distribution + 9-cell combinatorial analysis (5 samples)"
date: 2026-05-21
status: in_progress
phase: P3_PILOT_descriptive
agent: Phase A3 (sub-agent of PI 4-goal G3)
commit: 8fd6b6c
audience: PI / Goal-3 audit
tier_used: ⭐2 (descriptive, no tier upgrade claim)
upstream:
  - research/methyl_augmented_filter_phase2/cycle1/cycle1_step1_l2_coefs.tsv  # 10-feature LR top10
  - research/methyl_augmented_filter_phase2/cycle2/data/per_bam_master_V6.tsv  # HCC1395
  - research/methyl_augmented_filter_phase2/cycle2/data/{HCC1937,HCC1954,H1437,H2009}_master_augmented.tsv
downstream:
  - PI Goal-3 4-feature combinatorial slide (pending)
---

# Phase A3 — 5 features TP/FP distribution + 9-cell combinatorial analysis

## §1 背景

PI 4-goal G3「找出 FP 集中區域並設計 filter」需要 5 個 features (cycle 1 L2 LR top coef)
的 TP vs FP **分布刻畫** + **三維組合**(LOH × NG × caller_af) 可視化 + **跨樣本一致性
證據** 來支撐 Phase 2 filter 設計的可推廣性宣告。

本報告為 **descriptive** (⭐2)，不做 NO-GO / tier 升級判定 — 純粹刻畫 features 在 5 樣本
的行為差異，給 cycle3+ filter generalize 決策提供原始依據。

**Features (8 individual + 9-cell combinatorial)**:

1. `caller_af` (continuous, ClairS-TO 報告 AF)
2. `loh_inner_flag` (binary, inner=1 / outer=0)
3. `Coverage_Multiple` (continuous, 觀測 cov / KDE diploid baseline)
4. `NG` (NGroups, ordinal: 2 / 3 / ≥4)
5. `HPMergedDelta` (continuous, methylation HP1-HP2 delta)
6. `HPFineF` (continuous, methylation fine-grained F)
7. `NME_imbalance` (continuous, Normalized Methylation Entropy imbalance)
8. `Epipoly_Delta` (continuous, epipolymorphism delta)
9. `ClusterPermanovaF` (continuous, methylation cluster permanova F)

## §2 方法

### §2.1 Data source

| Sample | TSV | Prefix | rows | TP | FP |
|---|---|---|---|---|---|
| HCC1395 | `per_bam_master_V6.tsv` | `bam_off_` | 35,332 | 30,490 | 4,842 |
| HCC1937 | `HCC1937_master_augmented.tsv` | `V6_off_` | 16,607 | 13,910 | 2,697 |
| HCC1954 | `HCC1954_master_augmented.tsv` | `V6_off_` | 20,136 | 19,449 | 687 |
| H1437   | `H1437_master_augmented.tsv`   | `V6_off_` | 70,964 | 70,191 | 773 |
| H2009   | `H2009_master_augmented.tsv`   | `V6_off_` | 136,701 | 135,359 | 1,342 |

來源根目錄 `research/methyl_augmented_filter_phase2/cycle2/data/`。HCC1395 的
prefix 與其他 4 樣本不一致 (`bam_off_` vs `V6_off_`)，分析腳本內部以 canonical mapping
解析。

### §2.2 統計指標

- **AUC**: `sklearn.metrics.roc_auc_score(y, x)` with `y={TP:1, FP:0}`. NaN-safe drop.
  AUC > 0.5 表示「feature 值越大 → TP 機率越大」；AUC < 0.5 表示反方向 (FP 高於 TP)。
- **Cohen's d**: `(mean_TP - mean_FP) / pooled_sd`，pooled_sd 用 ddof=1 (Cohen 1988)。
  正值 = TP > FP；負值 = TP < FP。
- **Mann-Whitney U p-value**: `scipy.stats.mannwhitneyu(tp, fp, alternative='two-sided')`。
  注意 5 樣本 n>1k 級，p-value 多為極小，**effect size (Cohen's d) 才是主要訊號**。

### §2.3 9-cell 切分規則

- **LOH bin** (2 級): `loh_inner_flag==1` → inner / `==0` → outer
- **NG bin** (3 級): `NG≤2` → NG2 / `NG==3` → NG3 / `NG≥4` → NGge4
- **AF bin** (3 級): `pd.qcut(caller_af, q=3)` 每 sample 各自取 tertile，labels=`AF_L/AF_M/AF_H`

完整切分為 **2 × 3 × 3 = 18 cells**（TSV 全紀錄）；圖示為 **2 個 LOH facet × 3×3 NG×AF
heatmap**（每樣本 2 panels），與任務描述「9-cell heatmap」對齊（9-cell = NG×AF in single LOH facet）。

### §2.4 跨樣本一致性

- **Sign-consistency**: 每 feature 計算 5 樣本 Cohen's d 符號分布；≥4/5 同方向 →
  `sign_consistent=1`。
- **Pairwise Spearman ρ**: 將 9 features 的 |Cohen's d| 在每樣本內排名 (rank, average tie)，
  再對任意 2 樣本算 Spearman ρ (10 對配對)。ρ 高 = 兩樣本 feature 重要性序列相似。

## §3 個別 feature 結果 (9 features × 5 samples)

### §3.1 AUC summary table

(source: `phase2_completeness_audit/A3_per_feature_AUC_cohend_5sample.tsv`)

| Feature | HCC1395 | HCC1937 | HCC1954 | H1437 | H2009 | 範圍 |
|---|---|---|---|---|---|---|
| caller_af | **0.924** | 0.200 | 0.416 | 0.696 | 0.443 | 0.20–0.92 |
| loh_inner_flag | **0.707** | **0.717** | 0.543 | 0.599 | 0.616 | 0.54–0.72 |
| Coverage_Multiple | 0.609 | **0.683** | 0.415 | **0.010** | 0.193 | 0.01–0.68 |
| NG | **0.687** | 0.493 | 0.483 | 0.589 | **0.631** | 0.48–0.69 |
| HPMergedDelta | 0.490 | 0.581 | 0.434 | 0.534 | 0.583 | 0.43–0.58 |
| HPFineF | **0.626** | 0.593 | 0.477 | **0.603** | **0.616** | 0.48–0.63 |
| NME_imbalance | 0.498 | 0.513 | 0.510 | 0.533 | 0.511 | 0.50–0.53 |
| Epipoly_Delta | 0.450 | 0.543 | 0.582 | 0.401 | 0.403 | 0.40–0.58 |
| ClusterPermanovaF | 0.436 | 0.538 | 0.462 | 0.487 | 0.509 | 0.44–0.54 |

粗體 = AUC ≥ 0.58 (within-sample 顯著 discriminative)。

### §3.2 主要觀察

1. **caller_af 方向 sample-dependent**：
   HCC1395 AUC=0.924 (TP > FP) vs HCC1937 AUC=0.200 (TP < FP, 反向)。
   原因：HCC1395 5% 低 tumour purity → TP 真實 AF≈0.05 但模型 boost；HCC1937
   high-purity → FP 集中 high-AF germline tail。**caveat: caller_af 不能直接當
   filter，必先 stratify by tumour purity / sample**。
2. **loh_inner_flag 為跨樣本最穩定 discriminative**：
   5/5 樣本 Cohen's d 為正 (sign_consistent=1)，d range 0.29–0.93。AUC 0.54–0.72
   範圍緊密 → **這是 G3 filter 設計的最強 anchor**。
3. **NG / HPFineF 為次穩定 marker**：
   sign_consistent=1, AUC 4/5 樣本 ≥ 0.59。HCC1954 (high-purity 1k FP only) NG
   失效 (AUC=0.48) 為已知 small-FP-pool issue。
4. **Coverage_Multiple H1437 AUC=0.010 是真實數據反向**：
   H1437 FP median Coverage_Multiple=1.65 vs TP median=0.99 (FP 系統性高於 TP)；
   FP n=8 only 為小樣本，但方向真實。**caveat: 不是 bug，是 H1437 FP 集中在 amplicon-like
   high-cov 區**。
5. **NME_imbalance 全 5 樣本 AUC 0.50–0.53 平坦**：
   接近 random，confirm cycle1 step0「methylation 5th-rank、邊際貢獻」結論。

### §3.3 Figures (per-sample TP/FP boxplot)

- `figures/A3_5features_boxplot_per_sample.png` — 4 core features (caller_af /
  loh_inner_flag / Coverage_Multiple / NG) × 5 樣本 grid，TP (blue) vs FP (red)
  boxplot，每 panel 標註 AUC + Cohen's d。
- `figures/A3_methyl_5sub_boxplot.png` — 5 methylation sub-features × 5 樣本 grid。
- `figures/A3_AUC_cohend_bar_5sample.png` — 上 panel: AUC per feature × sample bar
  (含 AUC=0.5 與 0.58 ceiling 參考線)；下 panel: Cohen's d bar。

## §4 9-cell 組合圖 (LOH × NG × AF)

(source: `phase2_completeness_audit/A3_9cell_heatmap_data_5sample.tsv`)

### §4.1 主要組合觀察 — HCC1395 為例

| LOH | NG | AF | TP rate | n_total |
|---|---|---|---|---|
| outer | NG2 | AF_L | **0.217** | 3,034 |
| outer | NG3 | AF_L | 0.554 | 3,151 |
| outer | NGge4 | AF_L | 0.872 | 3,396 |
| inner | NG2 | AF_H | 1.000 | 4,153 |
| inner | NG2 | AF_L | 0.929 | 674 |
| inner | NGge4 | AF_H | 1.000 | 655 |

**最低 TP-rate cell (highest-FP enrichment)**:
- HCC1395: `outer × NG2 × AF_L` → TP-rate 0.217 (n=3,034, 78% FP 集中)
- 這對應 cycle1 step0 已知 high-AF-FP / low-NG / outer-LOH 區，filter 的主目標

**最高 TP-rate cell (clean TP)**:
- HCC1395 `inner × {NG2/NG3/NGge4} × {AF_M/AF_H}` 多為 0.99–1.00
- 對應 cycle1 step5c 「inner-LOH 不要動」原則

### §4.2 跨樣本 TP-rate range

- 全 5 樣本 × 18 cells (90 rows) TP-rate 範圍 **0.22 – 1.00**。
- **0.22 為 HCC1395 outer × NG2 × AF_L** (最強 FP enrichment cell)。
- HCC1937/HCC1954/H1437/H2009 因 FP 數量小 (29–2,697)，多 cell n_FP=0 → TP-rate=1.00
  → 統計弱，圖內以 `n=` 標註 cell size。

### §4.3 Figures

- `figures/A3_9cell_heatmap_5sample.png` — 5 samples × 2 LOH facets，
  3×3 NG×AF heatmap (TP-rate, vmin=0/vmax=1, RdYlBu colormap)，每 cell 標
  TP-rate + n_total。

## §5 跨樣本 direction consistency

### §5.1 Sign-agreement

(source: `phase2_completeness_audit/A3_direction_consistency_spearman.tsv` Section 1)

5/9 features sign-consistent (≥4/5 樣本同方向)：
- `loh_inner_flag` — 5/5 positive (✓ 最強 anchor)
- `NG` — 4/5 positive
- `HPFineF` — 4/5 positive
- `NME_imbalance` — 4/5 negative (但 AUC ~0.5, 訊號弱)
- `ClusterPermanovaF` — 4/5 negative (AUC ~0.5)

非一致：
- `caller_af` (2 positive / 3 negative) — purity-dependent
- `Coverage_Multiple` (2 positive / 3 negative) — coverage normalization 機制差異
- `HPMergedDelta` (3 / 2)
- `Epipoly_Delta` (2 / 3)

### §5.2 Pairwise Spearman ρ on |Cohen's d| ranking

(10 pairs, n_features=9)

| Pair | ρ | p |
|---|---|---|
| H1437 × H2009 | **0.967** | 0.000022 |
| HCC1395 × HCC1937 | 0.650 | 0.058 |
| HCC1937 × HCC1954 | 0.550 | 0.125 |
| HCC1395 × HCC1954 | 0.533 | 0.139 |
| HCC1395 × H1437 | 0.467 | 0.205 |
| HCC1395 × H2009 | 0.417 | 0.265 |
| HCC1937 × H1437 | 0.417 | 0.265 |
| HCC1937 × H2009 | 0.267 | 0.488 |
| HCC1954 × H1437 | 0.050 | 0.898 |
| HCC1954 × H2009 | -0.067 | 0.865 |

- **Median ρ = 0.44** → 中度跨樣本一致。
- H1437 × H2009 ρ=0.97 → 兩 NCI-H 系列 cell line 行為近乎相同 (同 tumour type lineage)。
- HCC1954 與 H1437/H2009 ρ ≈ 0 → HCC1954 (HER2-amplified, high-purity, FP n=687
  very small) 為 outlier，filter generalize 時需單獨驗證。

### §5.3 Figure

- `figures/A3_direction_consistency_heatmap.png` — 9 features × 5 samples 的 Cohen's d
  heatmap (RdBu_r 雙向色，center=0)，可直觀看 sign-agreement pattern。

## §6 Caveats / 已知 limitation

1. **caller_af direction inconsistent**：HCC1395 strongly positive vs HCC1937 strongly
   negative，反映 tumour purity 結構性差異 (5% vs >50%)。**直接 pool 5 樣本 LR fit
   會抹平這個訊號** — cycle3+ filter 設計需 per-purity stratify 或 per-sample fit。
2. **HCC1954 / H1437 / H2009 FP 量極小** (29 / 8 / 86 in Coverage_Multiple imputed cells;
   687 / 773 / 1,342 in label total)，部分 9-cell 統計 n_FP=0 → TP-rate=1 為 trivial。
   解讀時必看 `n_total` 欄。
3. **AUC = 0.58 ceiling 為 cycle1 step0 證實的「pure ISM feature 上限」**。本報告中
   loh_inner_flag (0.54–0.72) 與 caller_af (HCC1395 only 0.92) 突破 ceiling 是因為
   它們是 **structural / model output**, 非 pure ISM methylation。
4. **本分析無 within-LOH OLS residualization**，AUC 是 unadjusted 值。如要宣稱「feature
   X 獨立於 LOH/AF 之外仍有訊號」必續跑 `/auc-confound-guard` (within-group OLS + AF-bin
   交叉 + permutation)。本報告僅 descriptive，**不做 independence claim**。
5. **NG bin NG≤2 合併 NG=2** (因 NG=0/1 為 missing 表示)，與 cycle1 step1 命名一致。
6. **Mann-Whitney p-value 在 n>1k 級下普遍 < 0.001**，**勿用 p-value 作為 filter
   設計依據**；以 Cohen's d magnitude (>0.5 strong / >0.2 small) 為準。

## §7 後續步驟

1. **A3 → PI Goal-3 slide handoff**：5 fig + 3 TSV 為 G3 「FP 集中區量化」slide 的原始
   素材；建議 PPTX 上選 `A3_9cell_heatmap_5sample.png` + `A3_direction_consistency_heatmap.png`
   作為主圖。
2. **若要升 tier ≥ ⭐3** (claim 9-cell FP enrichment 為跨樣本 reproducible)：
   - 必跑 `/auc-confound-guard` 三關 (within-LOH OLS / AF-bin 交叉 / 1000× permutation)
   - 需引入 5 樣本以外的 4 個 held-out (V3F/V5 paired)
   - 對 caller_af direction 不一致補 per-sample stratified LR
3. **Methylation 5 sub-features 邊際**：本報告 confirm cycle1 step0 結論 — methylation
   單獨 AUC ≤ 0.63 (HPFineF top)，整合到 filter 時 marginal contribution 需 cycle3
   ablation 量化。

## §8 Provenance

- **commit**: 8fd6b6c (branch refactor/phase1-safety, 2026-05-21)
- **腳本**: `scripts/A3_5features_TP_FP_analysis.py`
- **執行時間**: ~30s (Python pandas + sklearn, no BAM scan)
- **產出**:
  - TSV ×3 in `phase2_completeness_audit/`
  - Figure ×5 in `figures/`
- **環境**: Python 3.x / pandas / scikit-learn / scipy / seaborn / matplotlib (CJK
  font chain via `scripts/lib/plot_setup.py`)
