---
title: G10 Quality Summary & Verification — 特徵觀察
group: G10
features_n: 27
features_available: 24
date: 2026-04-23
status: observation_complete
methodology: research/feature_layered_observation/02_methodology.md
verdict_global: NEGATIVE (synthetic), CONFOUND_COLLAPSED (LabelAllelePermanovaF)
---

# G10 · Quality Summary & Verification 特徵群組觀察

## 1. 特徵定義與來源（含 C++ file:line）

| Feature | Source file:line | Dtype | 說明 |
|---------|------------------|-------|------|
| HeuristicScore | SignificanceAnalyzer.cpp:378 compute_heuristic_score() | cont [0,1] | permanova + cluster + label 合成分數 |
| PassedGating | RegionProcessor.cpp:1608 sig_result.passed_gate | binary | global_p ≤ 0.1 |
| LabelHPPermanovaF/P/Valid | SignificanceAnalyzer.cpp label_hp_permanova_* | cont/binary | PERMANOVA F on distance ~ HP_family |
| LabelHPDispersionP/Warn | SignificanceAnalyzer.cpp label_hp_dispersion_* | cont/binary | PERMDISP ANOVA p on HP clusters |
| LabelAllelePermanovaF/P/Valid | SignificanceAnalyzer.cpp label_allele_permanova_* | cont/binary | PERMANOVA F on distance ~ allele (AF proxy 疑慮) |
| LabelAlleleDispersionP/Warn | SignificanceAnalyzer.cpp label_allele_dispersion_* | cont/binary | PERMDISP ANOVA p on allele |
| UnassignedAffinity/P | LabelTest.cpp:68 test_unassigned_affinity() | cont [-1,1] | HP3/HP0 reads 對 HP1/HP2 中心的偏向 |
| UnassignedDir | LabelTest.cpp unassigned_dir | cat {HP1,HP2,None} | 方向標籤 |
| NHP0 / NHP3 | LabelTest.cpp unassigned_n_hp0/hp3 | ordinal | HP0 / HP3 reads 數量 |
| DominantLabel | SignificanceAnalyzer.cpp:281 dominant_label_dimension | cat | {hp, allele, cluster, none} |
| Stability | SignificanceAnalyzer.cpp cluster_stability | cont [0,1] | cluster bootstrap 穩定度 |
| VerificationClass | SignificanceAnalyzer.cpp:330-339 | cat | {Strong, Subclone, Weak, Noise} |
| Quality_Score | RegionProcessor.cpp:241 compute_quality_score() | cont [0,100] | NumReads + NumCpGs + CovM + significance 合成分數 |
| Quality_Tier | RegionProcessor.cpp:296 determine_quality_tier() | cat | {High, Medium, Low, VeryLow} |
| Significant | RegionProcessor.cpp:1143 is_significant | binary | gate AND p≤0.05 AND CramersV≥0.1 AND NumReads≥20 |
| SuggestFilter | RegionProcessor.cpp:1149 suggest_filter | binary | label_delta > 0.3 (F1 optimization) |
| NHP_Somatic11/21/33 | ReadParser.cpp, RegionProcessor.cpp:1122 n_hp_somatic_* | ordinal | **UNAVAILABLE in canonical runs**（commit a61779c 後新增，canonical 皆 pre-fix） |

**資料源**：7 samples × 2 modes significance_summary.csv（canonical paired_full 7 + TO pilot 6）left-join 於 merged_with_vcf.tsv.gz 上，共 697,078/748,676 rows (93.1%) 帶 G10 欄位；COLO829 TO 無 archive，為 NaN。

## 2. 觀察目標

1. **Quality_Score** 跨 7 樣本/2 模式 AUC 是否全隨機 (~0.50)？
2. **LabelAllelePermanovaF** 殘差化 on vcf_AF 後 AUC 是否崩塌（AF proxy 假說驗證）？
3. **VerificationClass** 4 類別跨 7 樣本 TP rate 是否一致？
4. **NHP_Somatic11/21/33** 是否提供獨立訊號？→ 無法驗證（canonical 缺資料）
5. **SuggestFilter** 當作 FP filter 的 precision/recall 分佈？
6. **Quality_Score** 相對 **NumReads** baseline 增益？

## 3. 全域分佈 (Step 1)

Pooled N_TP = 583,029 / N_FP = 114,049。依 AUC 排序前 10:

| Feature | AUC | 95% CI | Cohen's d | mean_TP | mean_FP |
|---------|-----|--------|-----------|---------|---------|
| NHP0 | 0.764 | [0.763, 0.766] | 0.95 | 34.8 | 5.9 |
| LabelAllelePermanovaF | 0.627 | [0.625, 0.628] | 0.17 | 18.9 | 10.3 |
| LabelAllelePermanovaValid | 0.565 | [0.563, 0.566] | — | 0.93 | 0.80 |
| NHP3 | 0.545 | [0.543, 0.547] | 0.23 | 3.3 | 0.64 |
| Significant | 0.531 | random | — | 0.12 | 0.05 |
| PassedGating | 0.531 | random | — | 0.27 | 0.20 |
| HeuristicScore | 0.509 | random | 0.12 | 3.60 | 2.73 |
| **Quality_Score** | **0.497** | **random** | **0.009** | **92.28** | **92.17** |
| Stability / DispersionWarn / DispersionP | 0.500 | 常數 | 0 | — | — |
| UnassignedAffinity | 0.494 | random | -0.03 | — | — |

**關鍵發現**：
- **只有 NHP0 與 LabelAllelePermanovaF 的 AUC ≥ 0.58**（⚠ 兩者皆涉 confound，見 Step 5）
- Quality_Score AUC=0.497（完全隨機）— 全量重新驗證 2026-04 QS NEGATIVE 結論
- HeuristicScore AUC=0.509 — 合成分數全無訊號
- Stability 永遠為 0（看起來未實作 bootstrap）
- DispersionWarn 全部為 0（常數；no warning triggered）

Figure: `figures/G10_quality/fig01_global_distribution.png`

## 4. LOH×AF×CN 分層 (Step 2)

已對 Quality_Score / HeuristicScore / Stability / LabelAllelePermanovaF 繪 32-cell heatmap：
- Quality_Score 所有 cell 的 |Δ(TP−FP)| < 1%（散佈近零）
- Stability 全為常數 0（Δ 全零）
- HeuristicScore Δ 僅在 LOH_Subclone 某些 cell 出現 ~0.05 正方向
- LabelAllelePermanovaF 在 LOH_Weak/Strong/Noise × Extreme AF 分層見 Δ>3.0（TP>FP）

Figures: `fig02_{Quality_Score,HeuristicScore,Stability,LabelAllelePermanovaF}_heatmap.png`

## 5. 跨樣本一致性 (Step 3)

Per-sample AUC（7 × 2 = 14 cells）：

| Feature | paired_full (7 樣本) | to_pileup (6 樣本) |
|---------|---------------------|---------------------|
| Quality_Score | [0.50–0.74], median 0.59 | [0.49–0.51] **全隨機** |
| LabelAllelePermanovaF | [0.17–0.77] 方向不一致 | [0.52–0.62] 穩定微正 |
| HeuristicScore | [0.47–0.55] 近隨機 | [0.53–0.60] 微正 |
| NHP0 | [0.36–0.78] 極度不一致 | [0.50–0.57] 穩定微正 |

**關鍵發現**：
- Quality_Score 在 **paired_full** 某些樣本看似 >0.58（H2009 0.738, HCC1937 0.680），但 **TO mode 全部 ≤0.51**（證實 2026-04 QS TO 失效 ⭐5 結論）
- LabelAllelePermanovaF 在 paired_full HCC1937/HCC1954 居然 AUC=0.17-0.25（**方向翻轉**，意味 F 值高的 FP 比 TP 多）— AF distribution 跨樣本差異導致
- NHP0 在 paired_full HCC1395 AUC=0.355（低於隨機），與 HCC1937 的 0.784 相差 0.43 — **嚴重 SAMPLE_SPECIFIC**

Figure: `fig03_persample_auc_heatmap.png`

## 6. 分層 AUC (Step 4)

Quality_Score stratified：
- Global AUC = 0.497
- LOH strata: LOH_Weak 0.245, LOH_Strong 0.247, LOH_Subclone 0.332, LOH_Noise 0.422 — **全部 <0.50 翻轉方向**
- Mode: paired_full 0.550 / to_pileup 0.514
- CN: Loss 0.560 / 其他皆 ≈0.49

HeuristicScore stratified：
- LOH_Subclone 0.701（唯一 >0.58 cell，n=8,243） 但 LOH_Weak 僅 0.373
- 跨 layer 不一致

LabelAllelePermanovaF stratified：
- Global 0.627，LOH_Weak 0.697, LOH_Strong 0.685, LOH_Noise 0.668
- Mode: paired_full 0.682 / to_pileup 0.593
- 看似穩定，但需 Step 5 confound guard 驗證

Figure: `fig04_stratified_auc.png`

## 7. Confound 檢查 (Step 5)

殘差化 on `(NumReads, vcf_AF, Coverage_Multiple, AlleleDelta)`：

| Feature | raw AUC | resid AUC | Verdict | Extreme | Near-half | Intermediate |
|---------|---------|-----------|---------|---------|-----------|--------------|
| Quality_Score | 0.497 | 0.517 | **COLLAPSED** | 0.513 | 0.272 | 0.440 |
| HeuristicScore | 0.509 | 0.462 | **COLLAPSED** | 0.510 | 0.431 | 0.425 |
| **LabelAllelePermanovaF** | **0.627** | **0.496** | **COLLAPSED** | 0.612 | 0.614 | 0.690 |
| LabelHPPermanovaF | 0.477 | 0.430 | COLLAPSED | 0.492 | 0.202 | 0.393 |
| UnassignedAffinity | 0.494 | 0.578 | ATTENUATED | 0.494 | 0.516 | 0.503 |

**⚠ 重要結論**：
- **LabelAllelePermanovaF AUC 0.627 → 0.496（殘差化後崩塌）**，確認 **AF proxy 假說成立** — 所有訊號可由 `(NumReads, vcf_AF, CovM, AlleleDelta)` 解釋。與 O12 L2 collider bias 警告一致。
- Quality_Score / HeuristicScore：raw AUC 已近隨機，殘差化後仍 ≤0.52 — 合成分數對 variant discrimination 零貢獻
- UnassignedAffinity 殘差化後 AUC 意外上升到 0.578（可能殘留訊號與 NumReads 符號相反），但 AF-bin 全 ≈0.50 — 非真訊號

Figure: `fig05_confound_residualized.png`

## 8. Spatial autocorrelation (Step 6)

Quality_Score / LabelAllelePermanovaF bin = 5 Mb, 566 bins：
- 兩特徵皆 **無 artifact-suspect flag**（high_TP vs low_TP region ΔAUC < 0.08）
- per-bin AUC 分佈集中於 0.48-0.55（Quality_Score）與 0.55-0.70（LabelAllelePermanovaF）

Figures: `fig06_Quality_Score_spatial.png`, `fig06_LabelAllelePermanovaF_spatial.png`

## 9. VerificationClass 專題 (Step 7)

4 類別在 7 samples × 2 modes 的 TP rate：

| Mode | Noise | Weak | Strong | Subclone |
|------|-------|------|--------|----------|
| paired_full（7 樣本 median） | 0.993 | 0.989 | 0.988 | 1.000 |
| to_pileup（6 樣本 median） | 0.553 | 0.663 | 0.684 | 0.605 |

- **paired_full 全部 ≥0.97**，4 類別無差異（因為 paired_full FP 數量極少，TP rate 天花板）
- **to_pileup TP rate 落差大**（HCC1954 to_pileup Noise=0.213, Strong=0.340）但 **Strong/Weak 之間的差異不穩定**（HCC1395 Strong 0.748 vs HCC1937 Strong 0.604 vs HCC1954 Strong 0.340）— **無單一方向趨勢**
- Subclone 類別在各樣本 n=276-15,168，TP rate 範圍 0.21-1.00（完全不穩定）

**結論**：VerificationClass 分類對 TP 區分無一致訊號，不建議作 filter。

Figure: `fig07_verification_class.png`

## 10. SuggestFilter precision/recall (Step 8)

將 SuggestFilter=1 視為預測 FP：

| Stratum | flag_rate | precision (FP) | recall (FP) | TP loss rate |
|---------|-----------|----------------|-------------|--------------|
| GLOBAL | 2.1% | 0.033 | 0.004 | 2.4% |
| paired_full | 4.1% | 0.007 | 0.025 | 4.1% |
| **to_pileup** | 0.4% | **0.295** | 0.004 | 0.4% |
| HCC1937/to_pileup | 0.06% | **0.643** | 0.001 | 0.04% |
| HCC1954/to_pileup | 0.55% | **0.624** | 0.005 | 0.82% |
| LOH_Strong | 23.3% | 0.010 | 0.031 | 24.97% |

**結論**：
- **Global**：flag 2.1% rows, precision=3.3%（幾乎只 flag TP）— 作為 FP filter 無效
- **TO mode 單樣本（HCC1937/HCC1954）precision 達 0.62-0.64**，但 recall 低於 0.01%（只 flag 極少數 FP）— 極小覆蓋率不實用
- **LOH_Strong 區帶**：flag 23%, 但 TP loss rate=25%（flag 幾乎全是 TP）— 絕不可作 filter

Figure: `fig08_suggestfilter_prec_rec.png`

## 11. Quality_Score vs NumReads baseline (Step 9)

| Sample | Mode | AUC_QS | AUC_NumReads | Δ(QS−NR) |
|--------|------|--------|--------------|----------|
| HCC1395 | paired_full | 0.498 | 0.422 | +0.076 |
| HCC1395_DORADO | paired_full | 0.588 | 0.604 | -0.017 |
| HCC1937 | paired_full | 0.680 | 0.683 | -0.003 |
| HCC1954 | paired_full | 0.637 | 0.415 | **+0.222** |
| H2009 | paired_full | 0.738 | 0.720 | +0.018 |
| H1437 | paired_full | 0.556 | 0.674 | **-0.119** |
| COLO829 | paired_full | 0.531 | 0.511 | +0.019 |
| TO mode（6 樣本） | to_pileup | 0.49-0.51 | 0.48-0.51 | |Δ|<0.02 |

- paired_full 7 樣本 **Δ 中位數 = +0.018**（不顯著），方向不一致（H1437 反向 -0.12, HCC1954 正向 +0.22）
- to_pileup 6 樣本 Δ 全 <0.02（無增益）
- **Quality_Score 整體不提供相對 NumReads baseline 的穩定增益**

Figure: `fig09_qs_vs_numreads.png`

## 12. 論文與知識庫背景

- **Quality_Score 隨機結論**（2026-04 QS redesign evidence project_qs_redesign_evidence.md）：本觀察 4 維度複驗，確認 ⭐5 級結論
- **LabelAllelePermanovaF 疑為 AF proxy**（registry prior_conclusion）：Step 5 殘差化崩塌，確認成立
- **O12 L2 collider bias**（AlleleDelta=AF confound）：LabelAllelePermanovaF 同受此陷阱影響
- **Zone-Aware characterization only**（project_zone_aware_framework）：VerificationClass 印證此結論 — 分類能區分 TP rate 景象但不能作 filter

### 12.1 知識庫引用（Phase D）

查詢詞：`quality score variant caller` (top_score 99.1, full, high confidence)、`ClairS quality score filter` (top_score 102.3, high)。兩個主題皆高信心。

| kb_path | kb_title | 與 G10 的關聯 |
|---|---|---|
| `05_tools/variant-callers.md` | Variant Callers | **主要來源**：ClairS / ClairS-TO / DeepSomatic / ClairS-FSC 等 caller 的 QS 語意與版本差異；說明為何跨 caller 的 `Quality_Score` 不可直接比較 |
| `03_file_formats/vcf-clairs-to.md` | ClairS-TO VCF 詳細規格 | 13 種 FILTER 分類（LowAltBQ / LowSeqEntropy 等）與 INFO/FORMAT 欄位；直接對應本 G10 `SuggestFilter` / `verification_code` 設計 |
| `06_workflows/somatic-variant-calling.md` | Somatic Variant Calling 完整工作流程 | Coverage × MQ × QUAL 的 benchmark 閾值（MQ=20、QUAL cutoff）；定義 QS 特徵的可信區間 |
| `05_tools/longphase-to.md` | LongPhase-TO | Caller-specific polynomial 係數（clairs_to_ssrs / clairs_to_ss / deepsomatic_to）；說明 TO 模式 QS 來自 LongPhase-TO post-processing 而非 raw ClairS-TO |
| `04_databases/seqc2-truth-set.md` | SEQC2 Truth Set | HighConf / MedConf 分類是 benchmark 的 truth label；驅動 `verification_code` 與 `SuggestFilter` 的上游依據 |

**Quality_Verification**：1/1 查詢主題高信心命中。涵蓋 QS 計算、FILTER 定義、benchmark 閾值三層面。
**空缺**：KB 對 **跨 caller QS normalization** 方法（Phred-scale vs logit）與 `Quality_Score` 在 paired vs TO track 的可比性無專文；此為本 G10 Q1 質疑的根因，建議 Phase F 補寫。

### 12.2 外部文獻（Phase D）

**Quality score calibration in deep-learning somatic callers：**

1. **Chen, L., Zheng, Z., Su, J. et al. (2025).** "ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling." *Nature Communications* 16, 64547. DOI: **10.1038/s41467-025-64547-z** — ClairS-TO 使用 ensemble of two opposite-task neural networks（likely-vs-not-likely-somatic），QS 是 two-network consistency score。與 G10 關聯：**直接解釋 Quality_Score AUC=0.497** — ClairS-TO 的內部 QS 已經自我校正 consistency，InterSubMod pipeline 下游合成的 Quality_Score（multi-signal weighted sum）**無法超越 ClairS-TO 本身的 QS**。方向：**支持** 我們 Quality_Score NEGATIVE 結論。

2. **Zheng, Z., Li, S., Su, J. et al. (2022).** "Symphonizing pileup and full-alignment for deep learning-based long-read variant calling." *Nature Computational Science* 2, 797–803. DOI: **10.1038/s43588-022-00387-x** — Clair3 的雙軌（pileup + full-alignment）QS 框架，是 ClairS 的上游 germline caller。與 G10 關聯：Clair3/ClairS 的 QS 已具 **depth-calibrated** 特性（low-coverage 區 QS 自動下修），這解釋了 **Quality_Score vs NumReads Δ 中位數 +0.018** — NumReads 訊號已內建於 QS，兩者高度重疊。方向：**支持** 我們「Quality_Score 不提供相對 NumReads baseline 的穩定增益」結論。

3. **Anderson, M. J. (2001).** "A new method for non-parametric multivariate analysis of variance." *Austral Ecology* 26(1), 32–46. DOI: **10.1111/j.1442-9993.2001.01070.pp.x** — PERMANOVA pseudo-F 原始論文；明確指出 pseudo-F 對 **group size imbalance 和 heterogeneous dispersion** 敏感。與 G10 關聯：`LabelAllelePermanovaF` 基於 allele label (REF/ALT) 的 cluster separation，但 **extreme AF 下幾乎全 ALT read** → F 值主要反映 AF → 殘差化崩塌。方向：**支持** LabelAllelePermanovaF 被 L2 collider bias 解釋為 AF proxy 的結論。

**文獻空白**：本 G10 的 `VerificationClass` 四分類（Noise/Weak/Strong/Subclone）是 InterSubMod 內部 taxonomy，無直接對照文獻；`HeuristicScore` 的加權設計也無直接比對物。

## 13. 結論與 3 項質疑

### Verdict

| 特徵 | Verdict | 理由 |
|------|---------|------|
| Quality_Score | **NEGATIVE** (AUC 0.497, confound stable) | 隨機；4 維度複驗 |
| HeuristicScore | **NEGATIVE** (AUC 0.509) | 近隨機；合成權重未優化 |
| LabelAllelePermanovaF | **CONFOUND_COLLAPSED** | raw 0.627 → resid 0.496；AF proxy |
| LabelHPPermanovaF | NEGATIVE (AUC 0.477) | 方向與 TP 反相 |
| UnassignedAffinity | NEGATIVE | AF-bin 交叉 ≈0.50 |
| NHP0 | **SAMPLE_SPECIFIC** (raw AUC 0.764) | 極度不一致（0.36-0.78）；疑為 NumReads 的相關（paired_full 7 vs TO 6） |
| NHP3 | NEGATIVE | AUC 0.545 隨機 |
| VerificationClass | **CHARACTERIZATION_ONLY** | 跨樣本無一致方向 |
| Stability / DispersionWarn / DispersionP | **NOT_IMPLEMENTED** | 全部常數（0 或 1） |
| PassedGating / Significant / SuggestFilter | NEGATIVE as filter | precision/recall 皆差 |

### 3 項質疑

1. **Stability 是否實作？** G10 全樣本 Stability 值皆為 0，但 registry 記載「bootstrap stability」。應檢查 SignificanceAnalyzer.cpp 是否有 bootstrap iteration，或該欄位是否從未被賦值（default constructor）。
2. **NHP0 的強 AUC (0.764) 是否 NumReads proxy？** paired_full 某些樣本 AUC=0.78（HCC1937），但殘差化未做於 NHP0。建議 step5 extension：殘差化 NHP0 on NumReads 後重評。**（若崩塌，則 NHP0 屬 G1 Coverage confound 子類）**
3. **SuggestFilter 在 TO mode 單樣本 precision 0.62 是否可善用？** 雖 recall <0.01%，但若用於「High-confidence FP flag」只會影響極少 row。需評估對 F1 的 isolated contribution（可能低於量測誤差）。

### 邏輯鏈

- Quality_Score / HeuristicScore 為 **多欄位加權合成分數**，但加權係數是經驗設計（非訓練）→ 期望其 AUC 落於子特徵的中間位置 → 本觀察實測僅 0.50 → **合成權重對 TP 區分未優化**
- LabelAllelePermanovaF 測距離矩陣對 allele label 的分離度，但 allele label = {REF, ALT}，對高 AF variants（Extreme class）幾乎所有 reads 都是 ALT → F 值主要反映 AF → **殘差化後失去訊號符合預期**
- VerificationClass 的 Strong/Weak/Noise/Subclone 分類來自 PERMANOVA + cluster 結果組合，但這些底層訊號（LabelHPPermanovaF 等）全部 collapsed → 分類結果可解釋 FP 聚集景象但對預測無幫助

### 後續建議

1. **P2**：對 NHP0 做 Step 5 confound guard（on NumReads + vcf_AF），釐清是否為 Coverage confound
2. **P3**：review SignificanceAnalyzer.cpp 的 Stability 實作邏輯（確認是否 bug）
3. **P3**：G10 所有合成特徵（Quality_Score, HeuristicScore, VerificationClass）應被標記為 **characterization-only / not predictive**，避免再做 filter 實驗
4. **NHP_Somatic11/21/33**：需要 post-fix 重跑 canonical 才能驗證，暫列於 Phase C 下一批次

---

## 檔案清單

- Script: `scripts/observe_G10_quality.py`
- Report: `features/G10_quality_verification.md`（本文件）
- Figures (13): `figures/G10_quality/fig01-09_*.png`
- Tables: `data/G10/G10_{global_stats,persample_auc,auc_table,confound,spatial,verification_class,suggestfilter_prec_rec,qs_vs_numreads}.tsv`
- Master: `data/G10/master_g10.tsv.gz`（697,078 rows × 80 cols）
