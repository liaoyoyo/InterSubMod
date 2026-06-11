<!--
build_date: 2026-05-15
agent: Step 4 cross-sample extension — 4 samples V6 ISM
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3 §Step 4)
status: in_progress
-->

# Step 4 — V6 ISM 4-sample cross-sample extension

> Characterization-only. Read-only against `phaseD_v6_5sample/` + master.tsv. 不評估 filter / ΔF1。

## 範圍

- 5 樣本（HCC1395 + H1437 + H2009 + HCC1954 + HCC1937）
- V6 only（H1437/H2009/HCC1954/HCC1937 沒有 V3F/V5 三方對照）
- COLO829 deferred（truth set 0600 權限）

## 檔案

| 檔案 | 內容 |
|---|---|
| `per_sample_master/{sample}_v6_master.tsv` | 每樣本 V6 region-level wide-format TSV |
| `per_sample_grid/{sample}_grid.tsv` | 每樣本 50-cell 3 軸 grid |
| `step4_per_sample_grid.tsv` | 5 樣本 × 50 cells 合併 long-format |
| `step4_consistency.tsv` | 跨樣本 per-cell direction + Wilcoxon |
| `step4_HCC1937_outlier_analysis.md` | HCC1937 outlier 分析 + chr17 BRCA1 |
| `step4_signature_candidates.md` | 跨樣本 signature 候選 cells |
| `figures/{sample}_facets.png` | 每樣本 facet heatmap (LOH × HP × cov) |
| `intermediate/` | per-sample build summary, JSON 摘要 |

## 流程

1. `scripts/build_per_sample_master.py` — 4 樣本 wide TSV
2. `scripts/build_grid_per_sample.py` — 50-cell grid + LR + 簡化 confound guard
3. `scripts/cross_sample_consistency.py` — Wilcoxon signed-rank n=5
4. `scripts/hcc1937_outlier_analysis.py` — HCC1937 vs others
5. `scripts/synthesize_signatures.py` — final candidates synthesis

## Decision rule

- n=5 exact Wilcoxon 兩側 min p = 0.0625（5/5 同號）
- Signature candidate: n_samples=5 + direction ≥ 4/5 + Wilcoxon p ≤ 0.0625
- HCC1937 sensitivity: n=4 (排除 HCC1937) 用 relaxed p ≤ 0.125（n=4 exact min）
