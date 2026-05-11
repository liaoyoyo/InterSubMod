# Feature Layered Observation — Methodology SOP

**Date**: 2026-04-23
**Scope**: G1-G10 functional feature groups × (LOH × AF × CN) cells × TP/FP × 7 samples × 2 modes

## 資料源

- **主資料**：`research/feature_layered_observation/data/merged_with_vcf.tsv.gz` (748,676 rows × 60 cols = 28 ISM + 32 vcf_*)
- **Registry**：`research/feature_layered_observation/data/feature_registry.tsv`（137 features）
- **重要陷阱**：master `AF` 欄位 = `|AlleleDelta|`（非 caller VAF）— 下游切分用 `vcf_AF` 作 caller VAF

## Step 1-6 觀察流程（每特徵 / 每群組）

### Step 1 · 全域分佈
- TP vs FP density (violin + KDE overlay)
- Stats: mean±std, median, Cohen's d, Mann-Whitney U p-value, AUC + Wilson 95% CI
- 輸出：`fig01_global_distribution.png`
- 必要圖例：TP/FP 顏色、樣本 n、AUC bar
- 判定：AUC≥0.58 進 Step 5 confound guard；AUC<0.53 標 NEGATIVE

### Step 2 · LOH × AF × CN 32-cell heatmap
- 切分：`LOH_Subtype (5) × AF_class (3, 用 vcf_AF) × cn_tier_F (5, CovM)` = 75 cells，有效 n≥20 ≈ 32
- Heatmap A：TP rate per cell（色階 0-1）
- Heatmap B：feature mean/median for TP only
- Heatmap C：feature mean/median for FP only
- Heatmap D：Δ(TP-FP) + Cohen's d（signed 色階）
- 輸出：`fig02_{tp_rate,feat_tp,feat_fp,delta}.png`
- 每 cell 標註 n

### Step 3 · 跨樣本一致性
- 7 samples × (LOH×AF×CN) grid heatmap
- 固定樣本順序：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829
- Spearman concordance matrix（7×7）on per-cell TP rate
- 輸出：`fig03_per_sample_consistency.png`

### Step 4 · 分層 AUC
- Global AUC / per-LOH AUC / per-AF AUC / per-CN AUC
- Bar chart with 0.50 (random) + 0.58 (Beyond-AUC ceiling) 虛線
- 輸出：`fig04_stratified_auc.png`

### Step 5 · Confound guard (auc-confound-guard)
- Within-cell OLS 殘差化 on NumReads + vcf_AF + Coverage_Multiple
- 比對 raw AUC vs residualized AUC
- 若 residualized AUC 掉到 ≤0.53 → confound NEGATIVE
- AF-bin 交叉驗證（Extreme/Near-half/Intermediate 各自計 AUC）
- 輸出：`fig05_confound_residualized.png`

### Step 6 · Spatial autocorrelation
- chr + pos 5Mb bin 聚合，per-bin AUC
- 若 AUC 只在 high TP rate 區（baseline > 80%）出現 → artifact warning
- 輸出：`fig06_spatial_autocorrelation.png`

## 圖表規範

- **標題**：明確目標（一句話）+ 樣本 + mode
- **座標軸**：帶單位（counts / rate / log10 / etc）
- **圖例**：TP/FP 固定藍/橘配色
- **Heatmap**：colorbar + threshold line（0.58、0.50）
- **n 標註**：每 cell / bar
- **字型**：Latin + CJK fallback（避免方塊）
- **等比縮放**：不擠壓（memory feedback_pptx_screenshot_rendering_rules）

## 樣本固定順序
`HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829`

## Verdict 規則

| AUC (global) | Confound guard | Cross-sample (>=5/7) | Verdict |
|-------------|----------------|----------------------|---------|
| ≥0.65 | Pass (≥0.55) | Yes | POSITIVE |
| 0.58-0.65 | Pass (≥0.55) | Yes | CONDITIONAL_POSITIVE |
| ≥0.58 | Fail (<0.55) | any | CONFOUND_COLLAPSED |
| 0.50-0.58 | any | any | NEGATIVE |
| any | any | No | SAMPLE_SPECIFIC |

## Feature.md 10 章節標準

1. 特徵定義與來源（含 C++ file:line）
2. 觀察目標
3. 全域分佈（Step 1）
4. LOH×AF×CN 分層（Step 2）
5. 跨樣本一致性（Step 3）
6. 分層 AUC（Step 4）
7. Confound 檢查（Step 5）
8. Spatial autocorrelation（Step 6）
9. 論文與知識庫背景
10. 結論與質疑（3 質疑 + 邏輯鏈 + 後續建議）
