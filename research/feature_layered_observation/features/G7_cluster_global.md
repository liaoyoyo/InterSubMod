---
title: G7 · Cluster & Global Methylation Structure · Phase C feature observation
date: 2026-04-23
status: DRAFT · Phase C core output
owner: InterSubMod Research
scope: 11 features × 7 samples × 2 modes × LOH×AF×CN stratification (n=748,676 variants)
data: research/feature_layered_observation/data/G7/ (+ G7_*.tsv sidecars)
figures: research/feature_layered_observation/figures/G7_cluster/ (29 PNGs)
script: research/feature_layered_observation/scripts/observe_G7_cluster.py
related:
  - ../02_methodology.md
  - ../01_feature_inventory.md
  - ../features/G6_hp_fine.md
  - ../../../docs/reports/research_landscape/03_ISM分析價值界定.md
  - ../../../docs/reports/research_landscape/05_證據鏈總覽.md
---

# G7 · Cluster & Global Methylation Structure 特徵群組觀察

> **Thread A — 無監督甲基化分群 + Global p-value / PERMANOVA**。本群 11 特徵來自
> `src/core/GlobalTest.cpp` + `SignificanceAnalyzer.cpp` + `LocalTest.cpp`，描述 region
> 內 read-level cluster 結構的顯著性與強度。本報告在 7 樣本 × 2 模式 pooled
> 748,676 variants（TP 616,831 / FP 131,845）上跑完 Step 1-6，並加入 3 個 G7 特有診斷：
> (a) P-value vs -log10 rank-invariance；(b) Pairwise Mean/Median collinearity；
> (c) PermanovaValid gating 的 F-stat AUC delta。

---

## 1. 特徵定義與來源

| Feature | Source | Type | Computation |
|---|---|---|---|
| GlobalP | `GlobalTest.cpp test_clusters()` | cont [0-1] | Fisher FFH p on cluster_labels × methylation contingency |
| CramersV | `MathUtils.hpp:98` | cont [0-1] | `sqrt(chi2 / (n * min(r-1,c-1)))` from GlobalP contingency |
| PairwiseMeanDist | `DistanceMatrix.cpp` | cont [0-1] | Upper-tri mean pairwise distance across valid reads |
| PairwiseMedianDist | `DistanceMatrix.cpp` | cont [0-1] | Upper-tri median |
| ClusterPermanovaF | `SignificanceAnalyzer.cpp:1609` | cont ≥0 | PERMANOVA pseudo-F `distance_matrix ~ cluster_labels` |
| ClusterPermanovaP | `:1610` | cont [0-1] | Permutation p |
| ClusterPermanovaValid | `:1611` | binary | ≥2 clusters with ≥min_n members |
| ClusterDispersionP | `:1612` | cont [0-1] | PERMDISP ANOVA p on within-cluster variance |
| ClusterDispersionWarn | `:1613` | binary | 1 if dispersion inflates PERMANOVA |
| LocalBestCluster | `LocalTest.cpp 1606` | cat {0..K-1,-1} | Best-matching cluster in local enrichment |
| LocalBestP | `LocalTest.cpp 1605` | cont [0-1] | Best p across local Fisher tests |

**派生特徵**（rank-invariant 但 residualization / 視覺化用）：
- `neg_log10_GlobalP`, `neg_log10_ClusterPermP`, `neg_log10_LocalBestP`

**資料源**：`merged_with_vcf.tsv.gz` 未含 G7 raw 特徵，透過 (sample, mode, RegionID, Chr, Pos)
從對應 run 的 `intersubmod_{tp,fp}/significance_summary.csv` merge 回填；raw vs master
tp_label mismatch = 0，verified n=748,676。

---

## 2. 觀察目標（關鍵問題）

1. **無監督 cluster 結構本身（ClusterPermanovaF）是否能區分 TP/FP？**
2. **-log10(GlobalP) 與原 P-value 哪個 AUC 高？**（rank-invariance 檢查）
3. **PairwiseMean/MedianDist 作為連續特徵，AUC 跨 mode 方向一致嗎？**
4. CramersV 是否已被 HPFineN family 取代？（2×2 sparsity prior 93% zero）
5. 是否所有顯著 AUC 都會在 OLS 殘差化後崩潰？

---

## 3. 全域分佈（Step 1）· `fig01_*.png`

### 3.1 AUC ranking（全 11 特徵 + 派生）

| 特徵 | kind | AUC | 95% CI | Cohen's d | 判定門檻 |
|---|---|---|---|---|---|
| **ClusterPermanovaF** (Valid=1 only) | cont | **0.6776** | [0.6744, 0.6808] | 0.13 | ≥0.58 候選 POSITIVE |
| PairwiseMeanDist | cont | 0.6591 | [0.6576, 0.6606] | 0.55 | ≥0.58 候選 POSITIVE |
| PairwiseMedianDist | cont | 0.6280 | [0.6264, 0.6295] | 0.44 | ≥0.58 候選 POSITIVE |
| neg_log10_LocalBestP | cont | 0.5443 | [0.5426, 0.5459] | 0.32 | NEGATIVE |
| LocalBestCluster | cat | 0.5394 | [0.5377, 0.5411] | — | NEGATIVE |
| CramersV | cont | 0.5317 | [0.5300, 0.5334] | 0.15 | NEGATIVE |
| neg_log10_ClusterPermP | cont | 0.5304 | [0.5287, 0.5321] | 0.14 | NEGATIVE |
| ClusterPermanovaValid | binary | 0.5301 | [0.5284, 0.5317] | — | NEGATIVE |
| neg_log10_GlobalP | cont | 0.5102 | [0.5085, 0.5119] | 0.13 | NEGATIVE |
| ClusterDispersionP | cont | 0.5000 | — | 0.00 | NEGATIVE（100% 值=1.0） |
| ClusterDispersionWarn | binary | 0.5000 | — | — | NEGATIVE（全 0） |
| GlobalP | cont | 0.4898 | — | -0.03 | 反向（FP 稍高 p） |
| ClusterPermanovaP | cont | 0.4696 | — | -0.14 | 反向 |
| LocalBestP | cont | 0.4557 | — | -0.16 | 反向 |

**觀察**：
- P-value 特徵 AUC 近 0.5（含負向，因低 p = TP）；取 `-log10` 後 AUC 數值不變（rank-invariance）
  → **答題 2**：-log10 與 raw P-value AUC 理論上等價（verified `fig01b_pvalue_vs_log.png`）；-log10
  僅視覺/殘差化有用，對 AUC 排序無助。
- **ClusterPermanovaF 是唯一通過 0.58 門檻的 cluster-structural 特徵**，且僅在 Valid=1
  子集有效（178,794 / 748,676 = 23.9% regions），見 §4.3。

### 3.2 Violin / histogram
- `fig01_violin_grid.png`：12-panel 連續特徵 TP/FP 分佈。
- `fig01b_pvalue_vs_log.png`：三對 P vs -log10 AUC 比對（恆等，up to NaN）。

---

## 4. LOH × AF × CN 分層（Step 2-4）· `fig02_*.png`, `fig04_*.png`

### 4.1 ClusterPermanovaF 分層 AUC（Valid=1 only）

| Layer | Group | AUC | n_pos | n_neg |
|---|---|---|---|---|
| global | all | 0.678 | 153,838 | 24,956 |
| **mode** | **paired_full** | **0.786** | 90,702 | 369 |
| **mode** | **to_pileup** | **0.506** | 63,136 | 24,587 |
| AF | Extreme | 0.696 | 121,775 | 19,809 |
| AF | Near-half | **0.261** (reversed) | 5,145 | 167 |
| AF | Intermediate | 0.602 | 26,918 | 4,980 |
| CN | CN_Diploid | 0.716 | 63,982 | 8,867 |
| CN | CN_Near1 | 0.696 | 53,338 | 8,310 |
| CN | CN_Gain | 0.651 | 24,144 | 4,555 |
| CN | CN_Loss | 0.603 | 4,142 | 1,177 |
| CN | CN_HighGain | 0.573 | 8,232 | 2,047 |
| LOH | LOH_Subclone | 0.685 | 7,384 | 977 |
| LOH | LOH_Strong | 0.658 | 35,656 | 3,264 |
| LOH | LOH_Noise | 0.607 | 1,762 | 294 |
| LOH | LOH_Weak | 0.553 | 2,537 | 384 |

**兩大隱患**：
1. **paired_full 高 AUC 主要因為 TP:FP ratio 極度不平衡**（90,702 TP vs 369 FP）→ AUC 膨脹；
   to_pileup 63,136 TP vs 24,587 FP 且 AUC=0.506，幾乎隨機。
2. **Near-half AF 層 AUC=0.261（反向）**，與 Extreme AF 相反 → 潛在 interaction 效應或 caller-dependence。

### 4.2 Pairwise distance 分層（`fig04_PairwiseMean/Median_stratified_auc.png`）
- 全局 AUC 0.659 / 0.628；paired_full 與 to_pileup 方向一致（皆 > 0.5），但幅度差（paired_full ~0.71, to_pileup ~0.58）
  → **答題 3**：Pairwise 作為連續特徵 mode 方向一致，但幅度 paired 強、TO 弱。

### 4.3 ClusterPermanovaValid gating 影響（`fig09_permanova_valid_gate.png`）
- 全體（含 Valid=0）F-stat AUC = 0.539
- Valid=1 gated AUC = 0.678 → **ΔAUC = +0.139**
- Valid rate 不均：per-sample 平均 81-86% 通過，但 ≈15-20% regions Permanova 不可用
- 結論：**F-stat 必須 gating，否則被 low-read-count invalid regions 稀釋**

---

## 5. 跨樣本一致性（Step 3）· `fig03_*.png`

- ClusterPermanovaF per-cell AUC Spearman 全 7 樣本 median ≈ 0.4-0.6（中度一致）
- paired_full 7 樣本方向一致，但 to_pileup 中 H2009 / HCC1395 近隨機
- HCC1954 在 paired_full 有最陡下降（因 Near-half AF 反向 + Permanova 僅 <1% FP）

---

## 6. Confound 檢查（Step 5）· `fig05_*.png`, `G7_confound.tsv`

**Within-row OLS 殘差化 on (NumReads, vcf_AF, Coverage_Multiple)**：

| 特徵 | raw AUC | resid AUC | Δ | 判定 |
|---|---|---|---|---|
| ClusterPermanovaF | 0.678 | **0.521** | -0.157 | **CONFOUND_COLLAPSED** |
| PairwiseMeanDist | 0.659 | **0.553** | -0.106 | borderline（仍 >0.55） |
| PairwiseMedianDist | 0.628 | 0.547 | -0.081 | 近 0.55 邊緣 |
| CramersV | 0.532 | 0.389 | -0.143 | NEGATIVE |
| neg_log10_GlobalP | 0.510 | 0.375 | -0.135 | NEGATIVE |
| neg_log10_LocalBestP | 0.544 | 0.366 | -0.178 | NEGATIVE |

**AF-bin 交叉**（ClusterPermanovaF raw→res）：
- Extreme: 0.696 → 0.520
- Near-half: 0.261 → 0.216（仍反向）
- Intermediate: 0.602 → 0.486

**核心結論**：ClusterPermanovaF 的 raw AUC 幾乎全由 NumReads + vcf_AF + Coverage_Multiple 共變解釋；
**PairwiseMeanDist 是 G7 唯一殘差化後仍有微弱訊號（0.553）的特徵**，但也未達 POSITIVE 門檻（0.58）。

---

## 7. CramersV 診斷 & HPFine 取代證據（`fig07_cramersv_diagnostic.png`）

| 指標 | 值 |
|---|---|
| CramersV == 0 rate (全域) | **88.7%** |
| prior 預期（research_landscape/03） | ~93% |
| TP rate at CramersV=0 | ≈ 樣本 baseline |
| TP rate at CramersV>0 | 輕微 enrichment |
| raw CramersV AUC | 0.532 |
| residualized | 0.389 |

**與 HPFine 家族比對**：
- CramersV AUC 0.532 vs **HPFineNGroups AUC 0.73+（G6 報告）**
- CramersV 受限於 Global 2×2 contingency（cluster × methyl-state）稀疏性；HPFine 4-bucket occupancy 提供
  遠高的訊號密度
- **結論**：CramersV（G7）**已被 HPFineNGroups (G6) / HPFine_NGroups_CF 完全取代**，無獨立保留價值

---

## 8. Pairwise 距離 collinearity（`fig08_pairwise_collinearity.png`）

- Spearman ρ(Mean, Median) = **0.957** → 極嚴重共線性
- Mean vs NumReads: ρ = 0.338（mean 更受 read count 影響）
- Median vs NumReads: ρ = 0.246
- 建議：若做 downstream 模型，**只保留 PairwiseMeanDist**（Median 攜帶訊息幾乎全冗餘）

---

## 9. Spatial autocorrelation（Step 6）· `fig06_*.png`

- ClusterPermanovaF、neg_log10_GlobalP、CramersV 三者的 per-bin AUC 與 bin TP rate 皆呈正相關
  （hi_tp bin AUC > lo_tp bin AUC）→ 空間 auto-correlation confound 存在，`feedback_spatial_autocorrelation_confound`
  警告命中
- 進一步證實 §6 的「大部分訊號來自 chr+pos-correlated covariates」診斷

### 9.1 知識庫引用（Phase D）

查詢詞：`Permanova methylation` (top_score 60.6, partial, high)、`hierarchical clustering methylation` (top_score 93.6, full, high)。兩個主題皆高信心。

| kb_path | kb_title | 與 G7 的關聯 |
|---|---|---|
| `05_tools/intersubmod.md` | InterSubMod | **主要來源**：Cluster heatmap 設計（hierarchical clustering of methylation patterns）；CpG pattern clustering 即本 G7 `ClusterPermanovaF` / `CramersV` / distance-based 特徵的上游 |
| `05_tools/methyl-somatic-analysis.md` | MethylSomaticAnalysis (MSA) | Clustering 定義「methylation pattern 的群組結構」；與 G7 distance matrix 設計相通 |
| `06_workflows/methylation-analysis.md` | Methylation 分析工作流程 | Read × CpG methylation 原始矩陣是 PERMANOVA F 與 pairwise distance 的直接輸入；`features/feature_matrix.tsv` 為下游特徵 |
| `03_file_formats/modcall-vcf.md` | Modcall VCF 規格 | per-site DNA methylation 狀態（GT/MD/UD/DP 欄位）是 Read × CpG matrix 填充依據 |

**Cluster/Global Methyl**：2/2 主題高信心命中。知識庫覆蓋 clustering pipeline 與 methylation 原始資料結構，但對 **PERMANOVA F 統計** 與 **Cramér's V 顯著性** 的具體計算方法（permutation schema、2×2 vs more levels）無專文；建議 Phase F 把 `ClusterPermanovaF` / `ClusterPermanovaValid` gate 邏輯回寫 KB。

---

## 10. 結論與質疑

### 10.1 判定結果（per SOP Verdict 規則）

| Feature | Global AUC | Confound guard | Cross-sample 5/7? | Verdict |
|---|---|---|---|---|
| ClusterPermanovaF | 0.678 | **Fail (0.521)** | 偏 paired-only | **CONFOUND_COLLAPSED** |
| PairwiseMeanDist | 0.659 | 0.553（邊緣） | Yes | **CONDITIONAL_POSITIVE**（僅 paired，需進一步分層） |
| PairwiseMedianDist | 0.628 | 0.547 | Yes | **CONFOUND_COLLAPSED**（與 Mean 幾乎冗餘） |
| CramersV / GlobalP / LocalBestP / Dispersion* | 0.49-0.54 | n/a | — | **NEGATIVE** |

### 10.2 三個質疑

1. **`ClusterPermanovaValid` gate 掉 23.9% region 是否洩漏 label？**
   - 若 Valid=0 主要是 low-n / low-coverage region，而 low-coverage 又與 FP rate 相關，
     則 gating 本身就是 selection bias。需下游 H2 假說：先在 Valid==0 子集測 ClusterPermanovaF
     的 imputed-constant vs random AUC（目前 0.539 raw vs 0.678 gated）
   - 建議：對齊 G6 HPFineNGroups，做 `Valid × TP rate` 2D table 檢驗是否獨立

2. **paired_full AUC 0.786 是否純粹是 TP:FP ratio (245:1) 造成？**
   - Wilson CI 雖窄但 AUC 在極不平衡下對 ranking 敏感度低
   - 建議：per-sample balanced downsample（每 sample 取 min(n_TP, n_FP)）重算 AUC

3. **Near-half AF 層 AUC 反向（0.261），意義為何？**
   - 可能是 caller 在 Near-half 偏向保留 germline-like cluster
   - 或是 Near-half FP 本身有 methylation cluster 信號而 TP 無
   - 建議：per-sample 觀察 Near-half 的 ClusterPermanovaF TP vs FP distribution，確認方向

### 10.3 後續建議

1. **降維**：G7 11 特徵中，只保留 **PairwiseMeanDist** 作為 continuous cluster strength 代表
   （殘差化後仍 0.55）；其餘可歸入 "negative/replaced" 清單
2. **CramersV 排除**：HPFine 家族完全取代（AUC 0.73 vs 0.53）；建議 Phase E 彙整時將
   CramersV/GlobalP/LocalBestP 列入 "beyond-AUC ceiling 0.58 已耗盡" 清單
3. **與 G6 交叉**：檢驗 PairwiseMeanDist 是否在 HPFineNGroups 控制後仍有獨立訊號（nested OLS）
4. **不建議繼續投資**：ClusterPermanovaF 在 OLS 殘差化後 AUC=0.521，與 feedback_spatial_autocorrelation_confound
   + feedback_L2_collider_bias 兩個歷史陷阱共振，無額外研究價值

---

## 11. 論文與知識庫背景（Phase D §9.2 等價節）

### 11.1 內部參照

- `feedback_spatial_autocorrelation_confound` — 本 §9 Step 6 直接確認 chr+pos 聚合特徵在 hi-TP bin AUC 放大
- `feedback_L2_collider_bias` — OLS residualization 後 CramersV/PermanovaF 崩塌，與 L2 collider warning 共振
- `project_beyond_auc_exhaustion_confirmed` — G7 全數特徵 <0.58 再次確認 Beyond-AUC 0.58 ceiling

### 11.2 外部文獻（Phase D）

**PERMANOVA、甲基化 clustering 與 cancer heterogeneity 文獻：**

1. **Anderson, M. J. (2001).** "A new method for non-parametric multivariate analysis of variance." *Austral Ecology* 26(1), 32–46. DOI: **10.1111/j.1442-9993.2001.01070.pp.x** — PERMANOVA 原始論文，建立基於 distance matrix 的 non-parametric multivariate ANOVA 與 permutation-based p-value。與 G7 關聯：`ClusterPermanovaF` 和 `LabelAllelePermanovaF` 使用的就是此 pseudo-F 統計量。Anderson 明確指出 pseudo-F 值對 **group size imbalance 敏感**（Caveat in the original paper）— 與本 G7 發現 paired_full TP:FP=245:1 時 AUC 0.786 高度不穩一致。方向：**挑戰** 我們任何依賴 unbalanced pseudo-F 的結論。

2. **Stefansson, O. A., Moran, S., Gomez, A. et al. (2015).** "Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes." *PLOS ONE* 10(2), e0118453. DOI: **10.1371/journal.pone.0118453** — 在 breast cancer 上示範 hierarchical clustering 能區分 subtype，但強調需要 CpG filtering 和 normalization。與 G7 關聯：ISM 的 `ClusterPermanovaF` 基於 read × CpG matrix 的 cluster separation，理論上應能捕捉 bimodal methylation patterns。方向：**支持** PERMANOVA 作為 characterization 的正當性；**挑戰** 將其當作 TP/FP filter（Stefansson 用於 subtype discovery 而非 variant-level filter）。

3. **Carvalho, D. M., Maia, A. C. S., Melo, S. A. et al. (2019).** "DNA methylation profiles capturing breast cancer heterogeneity." *BMC Genomics* 20, 823. DOI: **10.1186/s12864-019-6142-y** — 報告 breast cancer 的 methylation heterogeneity 需要結合 clustering 與 differential methylation，且 **clustering 單獨不足以定義 subclones**。與 G7 關聯：直接支持我們結論「PairwiseMeanDist 殘差化後 0.55 接近 ceiling，ClusterPermanovaF 在 confound 下崩塌」。方向：**支持** G7 negative verdict — clustering strength 作為 variant filter 沒有獨立訊號的這個結論在 breast cancer methylation 文獻中是共識。

**文獻空白**：bulk ONT read-level methylation distance 矩陣（而非 bisulfite window-level）用於 **per-variant** TP/FP classifier 的直接對照尚無；本 G7 的 characterization-only 結論與 Landan 2012 / Carvalho 2019 的 population-level 用法一致。

---

## 附錄 A · 數據檔案

- `data/G7_global_stats.tsv` — Step 1 AUC / Cohen's d / MWU p (14 rows)
- `data/G7_auc_table.tsv` — Step 4 stratified AUC (global/LOH/AF/CN/mode × 6 flagship)
- `data/G7_cell_delta.tsv` — Step 2 32-cell Δ(TP−FP) medians
- `data/G7_confound.tsv` — Step 5 raw vs residualized + AF-bin cross
- `data/G7_pairwise_collinearity.tsv` — Mean/Median/NumReads Spearman
- `data/G7_summary.tsv` — 四項 headline metrics
- `data/G7/G7_enriched.tsv.gz` — RegionID-level enriched subset (reproducibility)

## 附錄 B · 圖表索引

29 PNGs 於 `figures/G7_cluster/`：
- `fig01_global_auc.png` / `fig01_violin_grid.png` / `fig01b_pvalue_vs_log.png`
- `fig02_{ClusterPermF, nlogGlobalP, CramersV, PairwiseMean, PairwiseMedian, nlogLocalBestP}.png`
- `fig03_{top3}_per_sample.png`
- `fig04_{6 flagship}_stratified_auc.png`
- `fig05_{6 flagship}_confound.png`
- `fig06_{top3}_spatial.png`
- `fig07_cramersv_diagnostic.png`
- `fig08_pairwise_collinearity.png`
- `fig09_permanova_valid_gate.png`
