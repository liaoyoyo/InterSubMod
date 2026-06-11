<!--
build_date: 2026-05-15
agent: Step 4 cross-sample extension — HCC1937 outlier analysis
scope: HCC1937 (BRCA1 mutant, high-ploidy CNV-driven FP) vs other 4 samples
-->

# Step 4 — HCC1937 Outlier Analysis

## 0. TL;DR

- HCC1937 marker rate 0.8165 邊緣低於 0.85 gate（phase D 唯一未過驗收的樣本）。
- 本檔比較 HCC1937 vs 其他 4 樣本（HCC1395/H1437/H2009/HCC1954）的 per-cell TP rate 差異，並用 n=4 sensitivity 評估訊號是否依賴 HCC1937。

## 1. Per-cell deviation (top 10 cells)

| cell_id | HCC1937_TP_rate | other_mean | deviation | HCC1937_n | n_FP | z_score |
|---|---:|---:|---:|---:|---:|---:|
| Inner|same_HP1|cov_normal | 0.3000 | 0.9979 | -0.6979 | 10 | 7 | -362.42 |
| Inner|same_HP2|cov_normal | 0.6364 | 0.9959 | -0.3595 | 44 | 16 | -61.77 |
| Inner|other|cov_loss | 0.9531 | 0.9944 | -0.0413 | 725 | 34 | -5.43 |
| Outer|other|cov_loss | 0.9565 | 0.9825 | -0.0260 | 69 | 3 | -0.88 |
| Inner|other|cov_gain | 1.0000 | 0.9882 | +0.0118 | 60 | 0 | 0.58 |
| Inner|other|cov_elevated | 0.9949 | 0.9855 | +0.0094 | 979 | 5 | 0.64 |
| Inner|other|cov_normal | 0.9812 | 0.9902 | -0.0090 | 5213 | 98 | -0.63 |
| Outer|other|cov_gain | 0.9974 | 0.9915 | +0.0060 | 783 | 2 | 0.48 |
| Outer|other|cov_elevated | 0.9959 | 0.9984 | -0.0025 | 1696 | 7 | -1.93 |
| Outer|other|cov_normal | 0.9939 | 0.9923 | +0.0016 | 2945 | 18 | 0.13 |

## 2. HCC1937 per-chromosome FP rate (top 10 by FP_rate)

| chr | n | n_FP | FP_rate |
|---|---:|---:|---:|
| chr15 | 304 | 13 | 0.0428 |
| chr17 | 330 | 14 | 0.0424 |
| chr14 | 403 | 13 | 0.0323 |
| chr1 | 997 | 30 | 0.0301 |
| chr22 | 175 | 5 | 0.0286 |
| chr2 | 944 | 25 | 0.0265 |
| chr11 | 612 | 14 | 0.0229 |
| chr5 | 786 | 14 | 0.0178 |
| chr3 | 904 | 15 | 0.0166 |
| chr9 | 577 | 9 | 0.0156 |

## 3. Signature sensitivity (with vs without HCC1937)

- n=5 (含 HCC1937) signature candidates: **1**
- n=4 (排除 HCC1937, Wilcoxon p≤0.125 relaxed) signature candidates: **1**
- → HCC1937 不影響 signature candidates 集合。

## 4. BRCA1 / chr17 / 高 ploidy 染色體解釋

| sample | chr17 n | chr17 n_FP | chr17 FP_rate |
|---|---:|---:|---:|
| H1437 | 902 | 0 | 0.0000 |
| H2009 | 3206 | 3 | 0.0009 |
| HCC1395 | 837 | 36 | 0.0430 |
| HCC1937 | 330 | 14 | 0.0424 |
| HCC1954 | 360 | 3 | 0.0083 |

> BRCA1 在 chr17q21。HCC1937 chr17 FP_rate 與其他樣本對比可揭露 mutation hotspot 對 FP 富集是否有貢獻。

## 5. 結論

- 完整數值見 `intermediate/HCC1937_outlier_per_cell.tsv` + `intermediate/HCC1937_fp_per_chr.tsv`
- HCC1937 signature 影響 sensitivity 結論寫入 `step4_signature_candidates.md`