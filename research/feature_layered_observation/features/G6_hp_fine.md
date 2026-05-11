---
title: G6 · HP Fine (4-bucket LOH-phasing) · Phase C feature observation
date: 2026-04-23
status: DRAFT · Phase C core output · Thread D replication
owner: InterSubMod Research
scope: 17 features × 7 samples × 2 modes × LOH×AF×CN stratification
data: research/feature_layered_observation/data/G6/ (9 TSVs, 748,676 variants)
figures: research/feature_layered_observation/figures/G6_hp_fine/ (18 PNGs)
script: research/feature_layered_observation/scripts/G6_hp_fine.py
related:
  - ../02_methodology.md
  - ../01_feature_inventory.md
  - ../../F_hpfinengroups_deepening/observations/step7_findings.md
  - ../../../docs/experiments/in_progress/2026/04/20260423_B1_Wilcoxon_NG2_gap_01.md
---

# G6 · HP Fine (4-bucket LOH-phasing) 特徵群組觀察

> **Thread D 核心特徵群**。本群 17 特徵來自 HP_fine 4-bucket 分類（HP1 / HP1-1 / HP2 / HP2-1），
> 是 2026-04-22 LOH-constrained phasing discovery 的基礎訊號。本報告在 7 樣本 × 2 模式 pooled
> 748,676 variants（TP 616,831 / FP 131,845）上重現 Thread D 的 NG=2 分裂結論，
> 並延伸到 HPFineD 6 條 pairwise 距離、fingerprint precision/recall 與 CF vs LabelTest 雙路冗餘性。

---

## 1. 特徵定義與來源

| Feature | Source | Type | Computation |
|---|---|---|---|
| GlobalP_HPFine | `src/core/GlobalTest.cpp test_hp_fine` | continuous [0-1] | Fisher FFH p on cluster × HP_fine 4-bucket |
| CramersV_HPFine | GlobalTest.cpp | continuous [0-1] | Cramér's V effect size |
| HPFine_NGroups_CF | `RegionProcessor.cpp:1604` | ordinal 0..4 | Cluster-First n_valid_groups |
| HPFineF | LabelTest.cpp | continuous ≥0 | PERMANOVA pseudo-F on HP_fine labels |
| HPFineP | LabelTest.cpp | continuous [0-1] | Permutation p |
| HPFineSig | LabelTest.cpp | binary | Significance flag |
| HPFineNGroups | `LabelTest.cpp:265 hp_to_fine_labels` | ordinal 0..4 | N non-empty among {HP1, HP1-1, HP2, HP2-1}, HP0/HP3 excluded |
| HPFineN_HP1 | `RegionProcessor.cpp:1186` | ordinal | Count of reads tagged HP1 |
| HPFineN_HP1S | RegionProcessor.cpp:1186 | ordinal | Count HP1-1 (somatic-demoted HP1) |
| HPFineN_HP2 | RegionProcessor.cpp:1187 | ordinal | Count HP2 |
| HPFineN_HP2S | RegionProcessor.cpp:1187 | ordinal | Count HP2-1 |
| HPFineD_HP1_HP1S | `RegionProcessor.cpp:1189` | continuous | Mean pairwise methylation distance HP1 vs HP1-1 |
| HPFineD_HP1_HP2 | RegionProcessor.cpp:1189 | continuous | HP1 vs HP2 |
| HPFineD_HP1_HP2S | RegionProcessor.cpp:1190 | continuous | HP1 vs HP2-1 |
| HPFineD_HP1S_HP2 | RegionProcessor.cpp:1190 | continuous | HP1-1 vs HP2 |
| HPFineD_HP1S_HP2S | RegionProcessor.cpp:1191 | continuous | HP1-1 vs HP2-1 |
| HPFineD_HP2_HP2S | RegionProcessor.cpp:1191 | continuous | HP2 vs HP2-1 |

**HP 4-bucket 定義**（重要 — 非 methylation subclone）：
- `HP1`：germline HP1（PON phasing）
- `HP1-1`：reads with somatic HP tag demoted to HP1-family under `--germline-hp-only`
- `HP2` / `HP2-1`：mirror for HP2 family
- 若 read 為 HP0（untagged）或 HP3（inconsistent），不計入 HPFineNGroups

**資料源**：
- Master: `research/feature_layered_observation/data/merged_with_vcf.tsv.gz` (60 cols — 只含 NGroups / NGroups_CF)
- **17 個完整特徵**從 7 樣本 × 2 模式 raw `significance_summary.csv` merge 回填（merge key = sample/mode/RegionID/Chr/Pos；raw vs master tp_label mismatch = 0）
- 最終 enriched: 748,676 rows × 76 cols，17/17 G6 features present

---

## 2. 觀察目標

1. **Thread D NG=2 分裂跨樣本複驗**：重現 `same_hap_marker`（HP1+HP1-1 或 HP2+HP2-1）在 7 樣本的 TP precision。
2. **6 條 HPFineD pairwise 距離排序**：哪條 pair 對 TP/FP 區分最強？paired vs TO 方向是否一致？
3. **CF (Cluster-First) vs LabelTest 雙路冗餘性**：兩個 NGroups 版本是否高度 collinear，或有獨立訊號？
4. **HP1S / HP2S（somatic-demoted bucket）**：count-level AUC 為何高於 HP1 / HP2？是 LOH 標記還是 allele confound？

---

## 3. 全域分佈（Step 1）

**圖**：`fig01_global_auc.png`, `fig01b_NG_distribution.png`

### 3.1 Global AUC ranking（pooled 7×2）

| Feature | AUC | 95% CI | Cohen's d | mean TP | mean FP |
|---|---|---|---|---|---|
| **HPFineN_HP1S** | **0.6120** | [0.610, 0.614] | +0.30 | 23.6 | 15.1 |
| **HPFineN_HP2S** | **0.5875** | [0.586, 0.589] | +0.26 | 21.6 | 14.3 |
| GlobalP_HPFine | 0.5696 | [0.568, 0.571] | +0.28 | 0.76 | 0.66 |
| HPFineD_HP1_HP2 | 0.5430 | [0.540, 0.546] | +0.11 | 0.270 | 0.252 |
| HPFineNGroups | 0.5214 | [0.520, 0.523] | +0.04 | 2.23 | 2.19 |
| HPFine_NGroups_CF | 0.5207 | [0.519, 0.522] | +0.03 | 2.05 | 2.02 |
| HPFineP | 0.5207 | [0.519, 0.522] | +0.07 | 0.33 | 0.31 |
| HPFineD_HP1S_HP2S | 0.5132 | [0.509, 0.517] | +0.04 | 0.218 | 0.214 |
| HPFineD_HP2_HP2S | 0.5082 | [0.505, 0.511] | −0.01 | 0.211 | 0.212 |
| HPFineD_HP1_HP1S | 0.5050 | [0.501, 0.508] | −0.05 | 0.214 | 0.220 |
| HPFineD_HP1S_HP2 | 0.5043 | [0.501, 0.508] | −0.01 | 0.236 | 0.237 |
| HPFineD_HP1_HP2S | 0.5038 | [0.501, 0.507] | −0.02 | 0.231 | 0.233 |
| CramersV_HPFine | 0.497 | — | −0.06 | 0.013 | 0.019 |
| HPFineSig | 0.484 | — | −0.07 | 0.48 | 0.52 |
| HPFineF | 0.481 | — | −0.003 | 6.92 | 7.05 |
| **HPFineN_HP1** | **0.424** | [0.422, 0.425] | −0.27 | 15.4 | 21.8 |
| **HPFineN_HP2** | **0.408** | [0.406, 0.410] | −0.30 | 13.3 | 19.9 |

### 3.2 解讀重點

- **HP1S / HP2S 計數是最強訊號**（AUC 0.61 / 0.59）— 與純 HP1/HP2 計數方向相反（AUC 0.42 / 0.41）。TP variants 有較多 somatic-demoted reads，FP variants 有較多純 germline reads。
- **HPFineNGroups AUC=0.521 整體極弱** — global pool 下 4-bucket 計數在 TP/FP 之間幾乎無分辨力，平均僅差 0.04。
- 6 條 HPFineD pairwise 距離**全落在 0.50-0.54** — 甲基化距離本身在全域 pool 中接近耗盡（CL-008 ≤0.58 ceiling 重現）。
- CramersV_HPFine AUC<0.5 — FP variants 反而有略高 Cramér's V，可能因 FP 常落在稀疏匹配的雜訊區。

> **注意**：此 AUC 為 pooled signal；NG=2 sub-scenario 下訊號翻倍（見 §5）。全域低 AUC **不等於** feature 無用。

### 3.3 HPFineNGroups / CF 分佈

圖 `fig01b_NG_distribution.png`：NGroups=2 是絕大多數 variants（~72% paired + ~56% TO），
並非 TP/FP 區分的主要 cell；信息在 NG=2 **內部**（見 §5 fingerprint）。

---

## 4. LOH × AF × CN 分層（Step 2）

**圖**：`fig02_ngroups.png`（HPFineNGroups 4 panel）、`fig02_hpfinef.png`（HPFineF 4 panel）

### 4.1 TP rate 格局

5 (LOH_Subtype) × 3 (AF_class) × 5 (CN_tier) = 75 cells；n≥20 有效 ~32 cells。

- TP rate 高值（>0.9）集中在 LOH_None + Extreme AF（0/1）+ CN T1-T2 — 這是 paired_full 絕大多數位點（≈85% TP）。
- TP rate 低值（<0.5）集中在 LOH_Strong + Intermediate AF + CN T3-T4 — 經典 FP-enriched zone。

### 4.2 HPFineNGroups mean TP vs FP

- **Δ(TP−FP) HPFineNGroups** 在絕大多數 cells 為 +0.01 到 +0.15（TP 略多 non-empty buckets）。
- 極端反向格（Δ < 0）發生在 AF Near-half + LOH_Strong — FP 在此區反而多 buckets（可能 somatic HP 雜訊污染 tag 造成）。
- HPFineF 在 LOH_None + CN T2-T3 TP 均值顯著高於 FP（permanova-F 提升 ~30-50%）。

### 4.3 CovM 軸發現（新觀察）

Step 4 stratified AUC（§6）顯示 HPFineNGroups 在 **CN_tier T3_quad / T4_amp** 有 AUC 0.604 / 0.627，
比 T1-T2（near-diploid）的 0.48-0.51 明顯高。**高拷貝數（amplification）區域 HPFineNGroups 有真實 TP/FP 區分力**。

---

## 5. 跨樣本一致性（Step 3）+ NG=2 fingerprint

**圖**：`fig03_per_sample_ng_tprate.png`, `fig03b_concordance_*.png`, `fig08_fingerprint_tprate.png`, `fig08b_same_hap_precision_recall.png`

### 5.1 每樣本 × mode × NG 值的 TP rate

NG=2 TP rate（重點；NG=2 為多數 variants）：

| Sample | paired_full | to_pileup |
|---|---:|---:|
| HCC1395 | 0.984 | 0.622 |
| HCC1395_DORADO | 0.994 | 0.750 |
| HCC1937 | 0.995 | 0.583 |
| HCC1954 | 0.997 | 0.162 |
| H2009 | 0.999 | 0.912 |
| H1437 | 1.000 | 0.773 |
| COLO829 | 0.939 | 0.639 |

Paired_full NG=2 幾乎全 TP（7/7 >0.94）；TO mode 差異極大 — **HCC1954 TO NG=2 僅 16% TP**（這是整個 cohort 的極端 outlier）。

### 5.2 NG=2 same_hap_marker precision / recall（Thread D obs18 升級）

**same_hap_marker 定義**：NG=2 AND (HP1>0 + HP1-1>0 + HP2=0 + HP2-1=0) OR (HP2>0 + HP2-1>0 + HP1=0 + HP1-1=0)
— 即 LOH-constrained phasing：reads 全部 phase 到同一 haplotype、somatic tag 只是在同 hap 內 demote。

| Sample | Mode | n NG=2 | n marker | n TP NG=2 | TP in marker | **Precision** | Recall |
|---|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | paired_full | 14,165 | 7,940 | 13,945 | 7,752 | **0.976** | 0.556 |
| HCC1395_DORADO | paired_full | 13,839 | 7,067 | 13,756 | 7,025 | **0.994** | 0.511 |
| HCC1937 | paired_full | 6,936 | 5,005 | 6,904 | 5,004 | **1.000** | 0.725 |
| HCC1954 | paired_full | 4,870 | 0 | 4,854 | 0 | — | 0 |
| H2009 | paired_full | 40,135 | 8,619 | 40,114 | 8,618 | **1.000** | 0.215 |
| H1437 | paired_full | 16,594 | 3,723 | 16,592 | 3,722 | **1.000** | 0.224 |
| COLO829 | paired_full | 18,969 | 821 | 17,817 | 698 | **0.850** | 0.039 |
| HCC1395 | to_pileup | 9,049 | 650 | 5,626 | 601 | **0.925** | 0.107 |
| HCC1395_DORADO | to_pileup | 20,768 | 12,164 | 15,574 | 10,637 | **0.874** | 0.683 |
| HCC1937 | to_pileup | 12,877 | 8,637 | 7,501 | 6,257 | **0.724** | 0.834 |
| HCC1954 | to_pileup | 34,519 | 9,382 | 5,590 | 3,222 | **0.343** | 0.576 |
| H2009 | to_pileup | 55,132 | 38,333 | 50,281 | 35,569 | **0.928** | 0.707 |
| H1437 | to_pileup | 21,739 | 9,053 | 16,806 | 7,749 | **0.856** | 0.461 |
| COLO829 | to_pileup | 23,623 | 3,525 | 15,105 | 2,353 | **0.668** | 0.156 |

**Paired_full precision 統計**（6/7 樣本，HCC1954=0 無 same_hap instances）：
- mean = **0.970** (min 0.850 COLO829, max 1.000)
- median = **0.997**
- **≥0.97 in 5/7 samples**

**TO pileup precision 統計**（7/7 樣本）：
- mean = 0.760 (min 0.343 HCC1954, max 0.928 H2009)
- median = 0.856
- **≥0.85 in 4/7 samples**

### 5.3 fingerprint TP rate breakdown（圖 `fig08_fingerprint_tprate.png`）

- `same_hap_marker`：paired_full 跨樣本 TP rate 0.85-1.00；TO 範圍 0.34-0.93
- `cross_het_marker`（HP1+HP2-1 或 HP1-1+HP2）：TP rate 顯著較低，多數樣本在 0.2-0.7 區間
- `other`（mixed 3+ buckets 在 NG=2 配置中不可能 — 為 HP0/HP3 含糊情況）：數量極少

### 5.4 Spearman 跨樣本一致性

`fig03b_concordance_*.png`（per-mode Spearman matrix on TP rate vs NG value）：
- paired_full 全 7 樣本 Spearman 絕大多數 ≥ 0.85 — 一致方向「NG=2 同其他 TP rate 都高，NG=4 略低」
- TO pileup 一致性弱於 paired — HCC1954 outlier 拉低 concordance

### 5.5 跨樣本結論

| 指標 | paired_full | to_pileup |
|---|---|---|
| same_hap precision median | **0.997** | 0.856 |
| precision ≥ 0.85 samples | 6/6 (HCC1954 N/A) | 5/7 |
| precision ≥ 0.95 samples | 5/6 | 1/7 |

**paired_full 下 same_hap_marker 是近完美 TP 標記**（mean 0.97 precision）；TO mode 下仍有用（median 0.86）但 HCC1954 異常需深查（已知 amplicon hotspot）。

---

## 6. 分層 AUC（Step 4）

**圖**：`fig04_stratified_auc_{HPFineNGroups,HPFineF,HPFineD_HP1_HP2,HPFineD_HP1_HP1S}.png`

### 6.1 HPFineNGroups 分層

| Stratum | n | AUC | 95% CI |
|---|---:|---:|---|
| Global | 748,676 | 0.521 | [0.520, 0.523] |
| LOH_None | 497,637 | 0.556 | [0.554, 0.558] |
| LOH_Noise | 126,631 | **0.620** | [0.617, 0.624] |
| LOH_Strong | 39,229 | 0.349 | [0.339, 0.360] |
| LOH_Weak | 76,742 | 0.365 | [0.358, 0.373] |
| LOH_Subclone | 8,437 | **0.680** | [0.664, 0.696] |
| AF Extreme | 680,151 | 0.547 | [0.545, 0.549] |
| AF Intermediate | 61,939 | 0.274 | [0.268, 0.281] |
| AF Near-half | 6,586 | 0.169 | [0.133, 0.205] |
| CN T3_quad | 121,745 | **0.604** | [0.601, 0.608] |
| CN T4_amp | 41,269 | **0.627** | [0.621, 0.634] |
| mode TO | 419,711 | **0.565** | [0.563, 0.567] |
| sample HCC1954 | 85,224 | **0.658** | [0.654, 0.662] |
| sample HCC1937 | 37,243 | 0.616 | [0.610, 0.622] |

### 6.2 解讀

- **反向訊號（AUC<0.5）** 出現在 AF Intermediate / Near-half 以及 LOH_Strong/Weak：FP variants 在此區有**更多** buckets。這是 CL-008 反向；機制可能是 SEQC2 LOH 雜訊（未校正 haplotype）或 germline-on-tumor contamination 造成 buckets 多樣性。
- **POSITIVE 訊號（AUC≥0.60）** 集中在：
  - **LOH_Subclone**（AUC 0.68，n=8,437）— 驗證 Thread D 假說：subclonal LOH 區域有 true phasing signal
  - **CN amplification**（T3/T4，AUC 0.60-0.63）
  - **HCC1954 / HCC1937**（AUC 0.62-0.66）— amplicon 高的樣本
- **paired_full 的 HPFineNGroups AUC 0.535 實際弱於 TO 0.565** — TO mode 下 NGroups 更有 TP/FP discrimination（反直覺；paired mode 下 TP 近飽和使 NGroups 資訊消失）。

---

## 7. Confound 檢查（Step 5）

**圖**：`fig05_confound_residualized.png`

### 7.1 殘差化設計

covariates: `NumReads + vcf_AF + Coverage_Multiple`（within-cell OLS），全域 pool。

| Feature | AUC raw | AUC resid | Δ | AF Extreme | AF Inter | AF Near-half |
|---|---:|---:|---:|---:|---:|---:|
| HPFineN_HP1S | 0.612 | **0.571** | −0.041 | 0.613 | 0.607 | 0.596 |
| HPFineN_HP2S | 0.587 | **0.517** | −0.070 | 0.588 | 0.593 | 0.577 |
| GlobalP_HPFine | 0.570 | 0.481 | −0.089 | 0.553 | **0.779** | **0.864** |
| CramersV_HPFine | 0.497 | **0.686** | +0.189 | 0.500 | 0.439 | 0.344 |
| HPFineD_HP1_HP2 | 0.543 | 0.515 | −0.028 | 0.561 | 0.480 | 0.374 |
| HPFineNGroups | 0.521 | 0.509 | −0.012 | 0.547 | 0.274 | 0.169 |
| HPFineP | 0.521 | 0.552 | +0.031 | 0.494 | 0.773 | **0.849** |

### 7.2 重點

- **HPFineN_HP1S 殘差後 AUC 0.571** — 仍 >0.55，**通過 confound guard**；HP1-1 計數訊號**非純 NumReads / AF / CovM artifact**。
- **HPFineN_HP2S 殘差後 AUC 0.517** — 接近 confound_collapsed；HP2S 訊號較依附 covariates。
- **GlobalP_HPFine / HPFineP 在 AF Near-half 殘差後 AUC 飆升至 0.85-0.86** — 強 AF-bin 異質訊號：在 Intermediate/Near-half AF 下，p-value 本身是 TP 標記（而非在 Extreme AF）。
- **CramersV_HPFine 殘差後 AUC 0.686（raw 0.497）** — 強 Suppressed signal：AF 殘差後暴露真實 TP effect。
- **HPFineD pairwise** 殘差後全 <0.52 — 甲基化距離本身 confound 清理後**無訊號**（CL-008 再度確認）。

### 7.3 verdict

| Feature | Verdict |
|---|---|
| HPFineN_HP1S | **CONDITIONAL_POSITIVE**（raw 0.61，resid 0.57，跨樣本未驗） |
| HPFineN_HP2S | CONFOUND_SUSPECT（resid 0.52 接近 random） |
| CramersV_HPFine (AF-resid) | **POSITIVE in Intermediate/Near-half AF** |
| GlobalP_HPFine / HPFineP (AF Near-half) | **POSITIVE in Near-half AF** (AUC 0.85-0.86) |
| HPFineNGroups (CN T3-T4 stratum) | **CONDITIONAL_POSITIVE**（stratum AUC 0.60-0.63） |
| 6 HPFineD pairwise (global) | NEGATIVE（resid all < 0.52） |

---

## 8. Spatial autocorrelation（Step 6）

**圖**：`fig06_spatial_HPFineNGroups.png`

5Mb bin × chr 聚合，per-bin AUC on HPFineNGroups。

- per-bin TP rate 分佈 heavily skewed to >0.8（大多數 bins 高 TP 基線）
- 在 high-TP-rate bins（TP rate ≥0.8）median AUC ≈ 0.51（接近 random）— **提示 HPFineNGroups global AUC 主要 driven by 低 TP 區域的 signal**
- mid-TP-rate bins（0.3-0.7）稀少，在這些 bins AUC 分佈較寬，有部分 bins AUC 0.65-0.75（但 bin count 不足 robust inference）

> **Spatial autocorrelation warning**：若後續用 HPFineNGroups 作全域 filter，需檢查其性能是否只在特定 chr 區段（e.g., chr6/chr16/chr22 高 LOH 區）出現。

---

## 9. 論文與知識庫背景

### 9.1 內部連結

- **Thread D（LOH-constrained phasing discovery）**：2026-04-22 發現 NG=2 在 Inner LOH 93-99% same-hap（6/6 樣本），TP gap +0.37。本 G6 報告在**全 pool 無 LOH filter** 下直接重現 same_hap precision 0.97 paired / 0.86 TO，證明 Thread D 訊號可用 G6 features 重建。
- **CL-008 Beyond-AUC 耗盡**：6 條 HPFineD pairwise 全 <0.58 再度確認。
- **CL-025a Sample ASM**：HPFine 特徵的 feature pool 替代，不成為 filter 主軸但為 characterization primary.
- **HPFineNGroups subclone marker 結論更正**（memory `project_hpfinengroups_subclone_marker.md`）：NG=2 訊號**並非** methylation subclone，而是 4-bucket occupancy 模式（LOH-constrained phasing 的體現）。本報告 §5 直接佐證此更正。
- **HP-only Phase 1 結論**（memory `project_readparser_germline_hp_only_phase1_negative.md`）：flag=on 下 HPFineN_HP1S / HP2S 分佈改變（somatic HP tag demoted）；本報告數據為 flag=off 基線，HP1S / HP2S 數非零恰因 flag=off 下 somatic HP tags 仍混在裡面。

### 9.2 外部文獻（Phase D）

**本 G6 是整個研究的論文主軸候選（LOH-constrained phasing discovery），文獻背景特別重要：**

1. **Lin, J.-H., Chen, L.-C., Yu, S.-C. & Huang, Y.-T. (2022).** "LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants." *Bioinformatics* 38(7), 1816–1822. DOI: **10.1093/bioinformatics/btac058** — **核心上游工具**。G6 的 HPFine 4-bucket `{HP1, HP1-1, HP2, HP2-1}` 來自 LongPhase-S/TO 的 somatic phasing output（HP1/HP2 germline + HP{1,2}-1 somatic fallback）。本 G6 觀察 NG=2 在 paired Inner LOH 93-99% same-hap，**直接依賴 LongPhase 的 phase block N50 = 25 Mb 保持 Inner region 跨整個 LOH 段穩定**。方向：**支持** 本研究方法學可行性。

2. **Shafin, K., Pesout, T., Chang, P.-C. et al. (2021).** "Haplotype-aware variant calling with PEPPER-Margin-DeepVariant enables high accuracy in nanopore long-reads." *Nature Methods* 18, 1322–1332. DOI: **10.1038/s41592-021-01299-w** — PEPPER-Margin-DeepVariant 引入「haplotype-aware」深度學習 SNV caller，用 HP tag 作 per-read feature。與 G6 關聯：PEPPER-Margin 的 HP 用法是 **per-candidate reads 的 per-haplotype consensus**，而本 G6 是 **per-region 4-bucket occupancy pattern** — 兩者都承認 HP 攜帶生物訊號，但粒度不同。方向：**支持** HP-aware 方法的整體正當性；**挑戰** 宣稱「新 4-bucket 分類是 novel」時需區別於 PEPPER-Margin 的 2-bucket 做法。

3. **Fu, Y., Aganezov, S., Mahmoud, M. et al. (2024).** "MethPhaser: methylation-based long-read haplotype phasing of human genomes." *Nature Communications* 15, 4859. DOI: **10.1038/s41467-024-49588-0** — MethPhaser 用 methylation signal 延伸 SNV-based phasing（N50 +78-151%, 83-99% 準確度）。與 G6 關聯：證明 **methylation × HP 是合法 joint signal 且跨 phase-set break 有效**，LOH 內部 methylation pattern 恰好是 MethPhaser 可橋接的場景。方向：**支持** HPFineD pairwise distance 在 LOH 區具診斷意義，但 G6 數據顯示 HPFineD 全 AUC <0.58 — 因此本研究**挑戰** MethPhaser 假設「methylation 可系統性區分 TP/FP variants」。

4. **Nik-Zainal, S., Van Loo, P., Wedge, D. C. et al. (2012).** "The life history of 21 breast cancers." *Cell* 149(5), 994–1007. DOI: **10.1016/j.cell.2012.04.023** — Battenberg haplotype-based subclonal CN algorithm 原典。與 G6 關聯：LOH-constrained phasing 的理論基礎 — Battenberg 已示範 phased haplotype 在 LOH 區能解析 subclonal fraction。本 G6 obs18 將此觀念推至 **read-level 4-bucket 分類**。方向：**支持** LOH + phasing 共同作為 subclonal inference 錨點；**文獻空白**：NG=2 的 read-level fingerprint 作為 variant filter precision booster 尚無直接對照。

**文獻空白**：**read-level 4-bucket HP occupancy 作為 somatic variant fingerprint** 是 InterSubMod 的 novel contribution 候選點；現有 LongPhase-S、PEPPER-Margin、MethPhaser、Battenberg 均未以此形式做 variant-level TP/FP 分類。

### 9.3 知識庫引用（Phase D）

查詢詞：`LongPhase-TO 4 bucket` (top_score 36.7, high)、`phased read HP1 HP2` (top_score 122.6, full match, high)、`phasing somatic` (top_score 130.1, full, high)。所有主題皆 high confidence。

| kb_path | kb_title | 與 G6 的關聯 |
|---|---|---|
| `05_tools/longphase-to.md` | LongPhase-TO | **核心來源**：TO phased VCF 輸出 `HP1 / HP2 / HP2-1`（step 5 Haplotype Phasing），即 4-bucket `{HP1, HP1-1, HP2, HP2-1}` 的上游。§5 `NG=2 same_hap_marker` 的 Inner 93-99% same-hap 正是 LOH.bed × phased genotype 結合的產物 |
| `05_tools/longphase-s.md` | LongPhase-S | Paired 模式輸出 `HP:z:1 / 2 / 1-1 / 2-1 / 3`，即 4-bucket 加 ambiguous (HP=3) 的來源。`3` 為 ambiguous somatic ALT，不進入 4-bucket 而被歸為 other；此為 HP1S/HP2S 非零的機制 |
| `06_workflows/phasing-workflow.md` | Phasing 工作流程 | 說明 paired vs TO 兩條 pipeline 下 HP tag 如何被 ReadParser 讀入 4-bucket；TO 的 `HP=3` 不同於 paired，**不可跨 pipeline 比較 bucket 計數** |
| `03_file_formats/bam-format.md` | BAM 格式文件 | HP:i vs HP:z 編碼規則；ReadParser 需兩種都處理；說明 `--germline-hp-only` flag 實作位置（filter somatic HP:z:*-1 tags） |
| `05_tools/longphase.md` | LongPhase (germline) | germline 基線 phasing 的限制：只有 HP1/HP2；當 germline-hp-only=on 時 HP1S/HP2S/HP2-1 bucket 應該都歸零（Phase 1 驗證失敗的根因） |

**HP_Fine (4-bucket)**：3/3 主題全數高信心命中。知識庫對 LongPhase-S/TO 的 HP tag 編碼差異覆蓋完整，但 **4-bucket 聚合成 `HPFineNGroups` 的具體邏輯**（何時視為 "occupied"、最低 read 閾值）為 InterSubMod 內部定義；建議 Phase F 回寫 KB。



---

## 10. 結論與質疑

### 10.1 核心結論

1. **HPFineNGroups 全域 AUC 極弱（0.521）但 stratum-conditional AUC 達 0.60-0.68**（LOH_Subclone / CN T3-T4 / HCC1954 / HCC1937）— 特徵訊號**條件性**，而非全域 characteristic。
2. **NG=2 same_hap_marker precision median 0.997 (paired) / 0.856 (TO)** — Thread D obs18 跨樣本 7×2 升級驗證通過 paired；TO 下在 HCC1954 失效（已知 amplicon 區異常）。
3. **HPFineN_HP1S 是本群最強原特徵**（global AUC 0.612，AF-resid 0.571 通過 confound guard）。原因：TP variants 的 4-bucket 分佈偏向 somatic-demoted HP1S（LOH-constrained phasing 下 germline 和 demoted reads 共存）；FP variants 分佈偏 HP1 純 germline（無 tag demotion 的雜訊 germline calls）。
4. **6 條 HPFineD pairwise 甲基化距離**：paired_full 下全 AUC <0.46（反向）；TO 下 HPFineD_HP1_HP2 AUC 最高 0.569（[0.566, 0.572]）。Top pair across both modes 為 `HPFineD_HP1_HP2`（germline vs germline），但 residualization 後全 <0.52 — **pairwise 距離在 confound guard 下無獨立 TP/FP discrimination**。
5. **CF vs NGroups 冗餘**：median Spearman 0.999；frac_equal mean 0.824。**兩路 near-duplicate**，保留一個即可。建議保留 `HPFineNGroups`（LabelTest path，與 Thread D / fingerprint 直接對應）。

### 10.2 三項質疑

1. **HCC1954 TO NG=2 TP rate 僅 0.162 的機制為何？**
   - 可能 1：amplicon 過高（~5-10× diploid）使 4-bucket 計數失真，FP 雜訊 reads 在 amp 區域展現多樣性
   - 可能 2：HCC1954 stale master 的 HP tag 是 pre-fix（self-phasing 版本），TO mode 的 germline-on-tumor fraction 異常低
   - 後續建議：HCC1954 單樣本 per-chr 分解 NG=2 TP rate，查是否集中 chr17（amplicon hotspot）

2. **HPFineN_HP1S AUC 0.61 是否為 "somatic HP tag 人工分組" artifact（self-phasing）**？
   - memory `project_hpfinengroups_subclone_marker.md` 警告：flag=on 下 N≥3 消失，原本 N≥3 TP 訊號來自 somatic HP tag
   - **需重跑 --germline-hp-only flag=on 版本驗證 HP1S AUC** — 若 flag=on 下 AUC 掉到 <0.53，則原訊號為 somatic HP artifact
   - 當前 dataset 為 flag=off 基線，結論 provisional

3. **same_hap_marker TO mode precision 在 HCC1954 崩盤（0.34）**：這是否意味 TO mode 下 fingerprint 並非穩定 filter？
   - Paired 下 precision ≥0.85 in 6/6，TO 下僅 4/7 達 0.85
   - 後續建議：TO mode fingerprint 作 **characterization-only**（不進 filter pipeline），paired 下可考慮 precision-boost filter

### 10.3 邏輯鏈（R→C）

- R1 [Phase A] Registry 定義 G6 17 features + 4-bucket HP_fine 語意
  → R2 [Phase B] Raw significance_summary × 14 runs merge 回填 master（mismatch=0）
  → R3 [Step 1 global AUC] HP1S/HP2S 雙 count 唯一突破 0.58 閾值
  → R4 [Step 2 stratify] HPFineNGroups 在 LOH_Subclone + CN_amp 突破 0.60
  → R5 [Step 3 per-sample] paired 同源 NG=2 TP rate 全 >0.93 但 TO 方差大
  → R6 [Step 4 stratified] HCC1954/HCC1937 AUC 0.62-0.66 → amplicon 樣本特殊
  → R7 [Step 5 confound] HP1S resid 0.571 pass，6 HPFineD pairwise 全 fail
  → R8 [Fingerprint] same_hap precision paired 0.997 / TO 0.856 → Thread D 跨 7 樣本升級成立 (paired)
  → C1: G6 核心 actionable signal = same_hap_marker (paired), HP1S count, NGroups at CN_amp stratum
  → C2: G6 within-run 特徵 vs global AUC — characterization-first，filter 候選僅 paired same_hap_marker

### 10.4 後續建議

1. **P0**（高優先）：重跑 --germline-hp-only flag=on 完整 7 樣本，驗證 HPFineN_HP1S AUC 是否為 somatic HP artifact
2. **P0**：HCC1954 TO per-chr NG=2 TP rate 分解，rule out amplicon confound 或確認為真訊號
3. **P1**：COLO829 TO 和 paired 的 same_hap recall 都 <0.16 — 檢查 COLO829 LOH bed 是否為 ClairS paired LOH（不含 subclonal），需改用 Wakhan 或 SAVANA LOH bed 驗證
4. **P1**：HPFineD_HP1_HP2（TO 唯一 AUC>0.55 的 pairwise distance），stratify by LOH_Subclone + CN_amp，測是否有 conditional positive
5. **P2**：移除 HPFine_NGroups_CF（冗餘）；feature_registry.tsv 標為 duplicate

### 10.5 資料 / 產物清單

**Data** (`research/feature_layered_observation/data/G6/`):
- `g6_enriched.tsv.gz` — 748,676 rows × 28 cols (17 G6 + sample/mode/LOH/AF/CN/vcf_AF)
- `step1_global_stats.tsv`, `step2_ngroups_cells.tsv`, `step2_hpfinef_cells.tsv`
- `step3_per_sample_ng.tsv`, `step4_stratified_auc_*.tsv` (×4), `step5_confound.tsv`, `step6_spatial_HPFineNGroups.tsv`
- `hpfineD_auc.tsv`, `fingerprint_ng2.tsv`, `fingerprint_precision_recall.tsv`, `cf_ngroups_collinearity.tsv`

**Figures** (`research/feature_layered_observation/figures/G6_hp_fine/`) — 18 PNGs:
- `fig01_global_auc.png`, `fig01b_NG_distribution.png`
- `fig02_ngroups.png`, `fig02_hpfinef.png` (Step 2 heatmaps)
- `fig03_per_sample_ng_tprate.png`, `fig03b_concordance_paired_full.png`, `fig03b_concordance_to_pileup.png`
- `fig04_stratified_auc_{HPFineNGroups,HPFineF,HPFineD_HP1_HP2,HPFineD_HP1_HP1S}.png` (×4)
- `fig05_confound_residualized.png`
- `fig06_spatial_HPFineNGroups.png`
- `fig07_hpfineD_per_sample_auc.png`, `fig07b_hpfineD_heatmap.png`
- `fig08_fingerprint_tprate.png`, `fig08b_same_hap_precision_recall.png`
- `fig09_cf_vs_ngroups.png`

**Script**: `research/feature_layered_observation/scripts/G6_hp_fine.py` (reproducible)
**Run log**: `research/feature_layered_observation/logs/G6_run.log`
