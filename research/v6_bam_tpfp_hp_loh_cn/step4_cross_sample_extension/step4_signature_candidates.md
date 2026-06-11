<!--
build_date: 2026-05-15
agent: Step 4 cross-sample extension — signature candidates synthesis
scope: 5 samples (HCC1395 + H1437 + H2009 + HCC1954 + HCC1937), V6-only
decision_rule: n=5 + direction ≥4/5 + Wilcoxon p ≤ 0.0625 (n=5 exact min)
-->

# Step 4 — Cross-Sample Signature Candidates

## 0. TL;DR

- 跨樣本 n=5 signature candidates: **1 cells**
- 4/5+ direction concordant cells (任意 Wilcoxon): **4**
- HCC1937 sensitivity (n=4 排除 HCC1937, p≤0.125 relaxed) candidates: **1**

## 1. Global TP rates (per sample baseline)

| sample | global_TP_rate |
|---|---:|
| HCC1395 | 0.9832 |
| H1437 | 0.9999 |
| H2009 | 0.9994 |
| HCC1954 | 0.9984 |
| HCC1937 | 0.9845 |

## 2. Signature candidates (n=5, direction ≥4/5, Wilcoxon p ≤ 0.0625)

| cell_id | majority_sign | direction (n_above/n_below) | Wilcoxon p | mean Δ vs global | mean n | HCC1395 | H1437 | H2009 | HCC1954 | HCC1937 | n4-sens flag |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Outer|other|cov_high_gain | above_global | 5/0 | 0.0625 | +0.0069 | 370 | 1.000 (n=33) | 1.000 (n=11) | 1.000 (n=1587) | 1.000 (n=186) | 1.000 (n=34) | True |

## 3. n=5 exact Wilcoxon 門檻說明

- n=5 Wilcoxon signed-rank exact 兩側最小 p = 0.0625（5 個全部同號），故本 step 採 p ≤ 0.0625 作門檻。若 cell n_samples_valid < 5（缺樣本資料），改用 direction concordance 加效應量。
- HCC1937 移除後 n=4 exact min p = 0.125（4 個全部同號），故 sensitivity 用 p≤0.125 relaxed。

## 4. 解釋與後續

- **若 candidates ≥ 1**：升級 candidate cells 為 cross-sample characterization signature，
  並依 plan §H7 conf-guard pass cells 進入 Step 5（H7 needs FDR + permutation pass）。
- **若 candidates = 0**：
  - 仍可看 direction ≥4/5 cells 作 sample-specific characterization 集合
  - HCC1395 + HCC1937 雙 outlier 可能在 same_HP1/cov_normal 與 cross_het 異常
  - 跨樣本 ceiling effect（V6 marker rate 全 ≥0.85）使 TP rate 已接近 1.0 → small absolute delta 難達 Wilcoxon 顯著

## 5. 完整檔案

- `step4_consistency.tsv` — 50 cells × {direction, Wilcoxon, per-sample rates}
- `intermediate/HCC1937_outlier_per_cell.tsv` — HCC1937 vs others deviation per cell
- `intermediate/HCC1937_fp_per_chr.tsv` — per-chr FP rate breakdown
- `intermediate/HCC1937_signature_sensitivity.tsv` — n=4 排除 HCC1937 重算
- `step4_per_sample_grid.tsv` — 5 樣本 × 50 cell 完整 grid
- `step4_HCC1937_outlier_analysis.md` — HCC1937 outlier 詳細分析