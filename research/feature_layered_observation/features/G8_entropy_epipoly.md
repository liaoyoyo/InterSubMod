---
title: G8 Entropy / Epipolymorphism / PerCpgASM — 特徵觀察
group: G8
features_n: 11
features_available: 11
date: 2026-04-23
status: observation_complete
methodology: research/feature_layered_observation/02_methodology.md
verdict_global: NEGATIVE (all features AUC <=0.56 after PerCpgASM_Valid gate; none exceeds Beyond-AUC 0.58 ceiling)
verdict_residualized: CONDITIONAL_POSITIVE (Entropy_Imbalance residualized AUC 0.638 -- but sign reverses from raw)
caveat: Paired_full n_fp <=2273 per sample (H1437=1, H2009=32) -- per-sample CI extremely wide; cross-sample claims unreliable
---

# G8 · Entropy / Epipolymorphism / PerCpgASM 特徵群組觀察

## 1. 特徵定義與來源（C++ file:line）

| Feature | Source | Dtype | 說明 |
|---------|--------|-------|------|
| NME_HP1 / NME_HP2 | PerCpgAsm.cpp:356-357 compute_nme(binary,2) | cont [0,1] | Normalized Methylation Entropy (window=2 CpGs) per HP |
| Entropy_Imbalance | PerCpgAsm.cpp:360 \|NME_HP1-NME_HP2\| | cont [0,1] | 兩 HP 之 entropy 絕對差 |
| Epipoly_HP1 / HP2 | PerCpgAsm.cpp:367-368 compute_epipolymorphism(binary,4) | cont [0,1] | Scherer 2011 epipolymorphism (window=4 CpGs) |
| Epipoly_Delta | PerCpgAsm.cpp:371 \|epipoly_hp1-epipoly_hp2\| | cont [0,1] | HP 間 pattern diversity 差 |
| PerCpgASM_Valid | PerCpgAsm.cpp per_cpg_asm_valid | binary | 充足覆蓋 + 兩 HP 皆存在 (gate) |
| Fisher_N_Sig / Fisher_Frac_Sig | PerCpgAsm.cpp:336 | ordinal / [0,1] | BH-FDR < 0.05 之 CpG 數 / 比例 |
| Fisher_N_Tested | PerCpgAsm.cpp fisher_n_tested | ordinal | Fisher 2x2 測試的 CpG 位點數 |
| Fisher_MaxNegLogFDR | PerCpgAsm.cpp:337 | cont [0,inf) | max(-log10 BH-FDR) across CpGs |

資料：enriched on master（748,676 rows → 730,738 with G8, 449,034 after PerCpgASM_Valid=True gate, 60.0%）。HCC1954 paired_full 無 G8（舊 pipeline），其餘 13 sample-mode cells 完整。

## 2. 觀察目標

1. 熵 / 多態性 / Fisher-ASM 是否提供 TP vs FP 鑑別力？（registry 標稱 Fisher_Frac_Sig paired AUC=0.726 最高）
2. PerCpgASM_Valid gate 是否為合法 eligibility filter？
3. Entropy_Imbalance 是否 ≈ |AlleleDelta|（entropy 變換的重複特徵）？
4. Epipoly_Delta 與 NME_HP1/HP2 是否獨立？
5. 殘差化後（within-group OLS on vcf_AF+NumReads+CovM）原始訊號是否保留？

## 3. 全域分佈（Step 1，PerCpgASM_Valid=True gated）

| Feature | AUC (gated) | 95% CI | Cohen's d (ungated) | 備註 |
|---------|-------------|--------|---------------------|------|
| Epipoly_HP1 | **0.557** | [0.555, 0.559] | +0.20 | 唯一接近 0.58 ceiling |
| Epipoly_HP2 | 0.547 | [0.545, 0.549] | +0.16 | |
| NME_HP2 | 0.505 | [0.503, 0.507] | +0.03 | near constant (median=1.0) |
| NME_HP1 | 0.500 | [0.498, 0.502] | +0.02 | near random |
| Entropy_Imbalance | 0.496 | [0.494, 0.498] | -0.03 | |
| Fisher_MaxNegLogFDR | 0.493 | [0.491, 0.495] | -0.08 | 方向錯（FP>TP） |
| Fisher_Frac_Sig | **0.479** | [0.477, 0.482] | -0.12 | 與 registry 聲稱 0.726 相反 |
| Fisher_N_Sig | 0.471 | [0.469, 0.473] | -0.10 | |
| Epipoly_Delta | 0.462 | [0.460, 0.464] | -0.13 | TP 比 FP 更 homogeneous |
| Fisher_N_Tested | 0.429 | [0.426, 0.431] | -0.26 | FP 測試的 CpG 位點數顯著更多 |

**關鍵發現**：**所有 G8 特徵 gated AUC ≤ 0.557**；無一通過 0.58 Beyond-AUC ceiling。Fisher 相關特徵方向與 registry 聲稱相反（FP > TP）。→ `fig01_G8_global_auc_bar.png`

## 4. LOH × AF × CN 分層（Step 2）

已對 Fisher_Frac_Sig / Entropy_Imbalance / Epipoly_Delta 繪 32-cell heatmap（ungated）：
- Fisher_Frac_Sig: LOH_None × Extreme AF 之 Δ(TP-FP) 接近 0；所有 LOH subtype 內部 ≤ 0.05 分散
- Entropy_Imbalance: 所有 cell Δ ≤ 0.002（接近尺度下限）
- Epipoly_Delta: 少數 LOH_Subclone cell 出現 Δ > 0.05，但 n < 50
Figures: `fig02_fisher_frac_sig.png`, `fig02_entropy_imb.png`, `fig02_epipoly_delta.png`

## 5. 跨樣本一致性（Step 3，paired_full 6 samples after gate）

**Fisher_Frac_Sig paired_full 7-sample AUC（task 關鍵問題 1 答案）**：

| Sample | n_TP | n_FP | AUC (gated) | CI | AUC (ungated) |
|--------|------|------|-------------|----|---------------|
| COLO829 | 24058 | 1463 | 0.500 | [0.48, 0.51] | 0.504 |
| H1437 | 32369 | **1** | 0.546 | [0.00, 1.09] | 0.651 (n_fp=8) |
| H2009 | 93092 | 32 | **0.637** | [0.55, 0.72] | 0.651 |
| HCC1395 | 14983 | 336 | 0.497 | [0.47, 0.53] | 0.482 |
| HCC1395_DORADO | 12830 | 90 | 0.511 | [0.45, 0.57] | 0.530 |
| HCC1937 | 5445 | 33 | **0.635** | [0.55, 0.72] | 0.631 |
| HCC1954 | — | — | N/A (無 paired G8) | | |

**0/6 samples ≥ 0.65，只 2/6 ≥ 0.58**。paired_full FP 數量極小（1-1463）造成 CI 橫跨整個區間。→ `fig_fisher_frac_sig_paired.png`。跨樣本 heatmap → `fig02_G8_per_sample_heatmap.png`。

## 6. 分層 AUC（Step 4）

`step4_stratified_auc_*.tsv` 顯示：LOH_Subclone 僅在 n 極小時 AUC 偏高；AF_class 分層無系統訊號；CN_tier 不顯著分化。Mode 間 paired_full vs to_pileup 方向不一致（Fisher_Frac_Sig paired median=0.51, TO median=0.49）。

## 7. Confound 檢查（Step 5，within-(sample,mode,LOH) OLS on vcf_AF+NumReads+CovM）

| Feature | AUC raw (gated) | AUC residualized | Δ | 結論 |
|---------|-----------------|------------------|----|------|
| Epipoly_HP1 | 0.557 | **0.500** | -0.057 | CONFOUND_COLLAPSED |
| Epipoly_HP2 | 0.547 | **0.497** | -0.050 | CONFOUND_COLLAPSED |
| Entropy_Imbalance | 0.496 | **0.638** | **+0.142** | 殘差化後訊號反向浮現 ⚠ |
| Fisher_Frac_Sig | 0.479 | 0.525 | +0.046 | 方向翻轉但小 |
| Epipoly_Delta | 0.462 | 0.533 | +0.071 | 方向翻轉 |

**task 關鍵問題 2 答案：Entropy_Imbalance 殘差化後 AUC 從 0.496 → 0.638（保留甚至增強訊號）**，但 (a) raw 方向為 TP<FP、resid 方向反轉；(b) 此現象與 G10 `LabelAllelePermanovaF` 的 AF collider bias 一致（O12 警告），非真獨立訊號。AF-bin 交叉：Extreme=0.497, Intermediate=0.481, Near-half=0.560（強異質 → AF collider 警示）。→ `fig05_G8_confound_residualized.png`

## 8. Spatial autocorrelation（Step 6）

Fisher_Frac_Sig 的 per-5Mb-bin AUC 與 bin TP rate 強相關（`fig06_spatial_Fisher_Frac_Sig.png`），mid-TP-rate（0.4-0.7）bin 中位 AUC ≈ 0.5，僅在 TP rate > 0.9 bin 中見 AUC > 0.6 —— 符合 spatial artifact pattern。

## 9. 論文與知識庫背景

- Scherer et al. 2011（epipolymorphism 原始方法）—— 在 bisulfite short-read tumor vs normal 有區分力，但此處 ONT MM/ML + 亞克隆細分後失效
- O11 heterogeneity NEGATIVE（n_reads confound）已收斂於 2026-03，此次 G8 再確認
- PerCpgASM 作為 gate 的 TP rate 差異：paired_full True vs False 差 ≤ 0.02，TO True vs False 差 -0.05 至 +0.04，證明 gate 非強 TP/FP 篩選器
- Entropy_Imbalance vs AlleleDelta 相關性：Global Pearson=0.017、LOH_Subclone Pearson=0.214（n=241）—— 非強線性相依，但殘差化後異常 AUC 提升暗示非線性 AF collider

## 10. 結論與質疑

**Verdict**: NEGATIVE（所有 G8 特徵 gated global AUC ≤ 0.557，無一達 Beyond-AUC 0.58）。Entropy_Imbalance 殘差化 AUC=0.638 為 CONDITIONAL_POSITIVE，但伴隨 AF collider bias 警示與方向反轉。

**Task 關鍵問題三答**：
1. **Fisher_Frac_Sig paired 7 樣本 ≥0.65？** NO — 0/6（gated）/ 2/6（ungated 但 n_fp<100 不可信）。
2. **Entropy_Imbalance 殘差化保留訊號？** YES 但方向翻轉 (raw 0.496 → resid 0.638)，疑 AF collider。
3. **Epipoly_Delta vs 單側 HP1/HP2 誰強？** Epipoly_HP1（AUC 0.557）> Epipoly_HP2（0.547）> Epipoly_Delta（0.462，反向）。單側 > Delta。

**三質疑**：
- Q1: Registry 宣稱 Fisher_Frac_Sig paired AUC=0.726 的原始來源為何？本次分析 gated 後全量 paired AUC=0.479 — registry 可能使用 positive class label 反轉或特定子集。
- Q2: paired_full FP 數量極端稀少（H1437 n_fp=1, H2009 n_fp=32）使 per-sample AUC 幾乎不可詮釋。是否應在 paired track 只保留 TO cross-validation？
- Q3: Entropy_Imbalance 殘差化 AUC 0.638 若確為 AF collider artifact，應配合 G9 ASM 的 collider 診斷（O12 protocol）再確認，不應單獨引用。

**後續建議**：G8 整體不進入 feature ensemble；若要深入 Entropy_Imbalance 殘差化訊號，需 L3 AF-bin cross + permutation test（auc-confound-guard skill）。PerCpgASM_Valid 可作 eligibility gate 傳遞到 G9 ASM，本身非 TP/FP discriminator。

**Artifacts**:
- `data/G8/G8_auc_table.tsv`, `G8_auc_table_gated.tsv`, `step3_per_sample_auc_gated.tsv`, `step5_confound_gated.tsv`
- Figures in `figures/G8_entropy/`: `fig01_G8_global_auc_bar.png`, `fig02_G8_per_sample_heatmap.png`, `fig05_G8_confound_residualized.png`, `fig_fisher_frac_sig_paired.png` (核心4張 + 完整腳本另產16張)
