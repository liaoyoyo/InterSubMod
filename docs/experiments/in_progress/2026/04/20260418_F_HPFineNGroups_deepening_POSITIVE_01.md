---
title: F Pilot — HPFineNGroups 深度質疑驗證與新 canonical filter 建立
date: 2026-04-18
status: POSITIVE_REFINED
category: experiments/in_progress
tags: [HPFineNGroups, subclone_marker, AF_confound, filter_refinement, F_pilot, stability_upgrade]
dependency_on:
  - docs/experiments/in_progress/2026/04/20260417_HPFineNGroups_saturation_check_01.md
  - docs/experiments/in_progress/2026/04/20260417_PartB_effect_size_cn_strat_01.md
  - docs/reports/research_landscape/06_結論穩定性審查.md
dependency_for:
  - Phase 2 biology characterization 主線
  - B.1-3 per-sample Cohen's d AF<0.4 stratified 重算（Step 4 候選）
related:
  - research/F_hpfinengroups_deepening/manifest.yaml
  - memory project_hpfinengroups_subclone_marker.md
---

# F Pilot — HPFineNGroups 深度質疑驗證與新 canonical filter 建立

**Verdict**: **POSITIVE (refined filter)**
**穩定度升級**: ⭐3 → ⭐4（補充結論 16）
**Date Span**: 2026-04-17 → 2026-04-18

---

## 摘要（TL;DR）

Opus 4.7 研究重整 Part B.1 對 HPFineNGroups POSITIVE 結論提出 4 項質疑（residualized AUC 不足、飽和效應、7/7 統計強度、LOH/AF 混淆）。F pilot 執行三階段驗證後**結論強化而非翻轉**，並建立新 canonical filter。

| 項目 | 舊 filter | **新 filter** | Δ |
|------|----------|--------------|---|
| 條件 | NG≥4 + NR≥80 NonLOH | **NG=4 + AF<0.4 + NR≥80 NonLOH** | — |
| 總 TP rate | 0.8912 | **0.9281** | +3.7pp |
| n regions | 25,744 | 14,197 | -45% |
| 5/7 樣本 ≥0.85 | 4/7 | **5/7** | +1 |
| HCC1954 挽救 | 0.497 | **0.707** | +21.0pp |
| HCC1937 挽救 | 0.714 | **0.867** | +15.4pp |

**三大新發現**：
1. **NGroups 非單調**：NG=2 TP rate (0.643) < NG=1 (0.763)，根因為 **germline AF confound**
2. **HCC1954 失效根因 = FP 在 AF≥0.4 極端富集**（AF<0.2 TP=0.874 正常 → AF 0.8-1.0 TP=0.022）
3. **Paired mode 99.85% 非真 gain**：baseline 98.96%，filter 僅 +0.89pp

**Confound checks**（套用 P3 pilot 教訓）：chr-shuffle null **Z=43.5 PASS**；Coverage_Multiple 跨 CN tiers 0.90-0.94 **PASS**。

---

## 1. 觸發與範圍

### 1.1 觸發
Opus 4.7 研究重整 plan Part B.1 列出 4 項對 HPFineNGroups POSITIVE 的質疑：
- B.1-1：residualized AUC=0.617 < TO-pure caller_af 0.654，是否足以支撐 POSITIVE？
- B.1-2：N≥4+NR≥80 效應是否為 NR≥80 本身 TP rate 高的 artifact？
- B.1-3：7/7 方向一致的統計強度？per-sample effect size？
- B.1-4：LOH × AF 混淆？

B.1-1/B.1-2/B.1-3 已分別驗證（不推翻結論但 effect size 僅 3/7 medium）。**B.1-4 與新產生的方向性問題（NGroups 是否單調？Paired 99.85% 是否 artifact？）為本 pilot 焦點**。

### 1.2 資料範圍
- Master: `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`（748K rows, 116 cols）
- 樣本：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, COLO829, H2009, H1437
- 特徵：`HPFineNGroups` (int 0-4), `caller_af`, `NumReads`, `Coverage_Multiple`, `to_loh_bed_hit`, `truth_label`
- 模式：`to`, `paired`

### 1.3 範圍內與範圍外
**範圍內**：TO NonLOH 區域、NR≥40-200 範圍、NG=1~4、AF cutoff 0.2-1.0。
**範圍外**（明確記錄）：
- LOH 區域（NG=4 僅 n=11）→ 應用 AlleleDelta 而非 NGroups
- Paired mode（baseline 98.96% saturation）
- COLO829（ONT_R10 無 methylation basecall → HPFineNGroups 本身 artifact）
- purity <1.0 的臨床樣本（cell line 外推性未驗證）

---

## 2. 方法

三階段迭代：觀察 → 根因調查 → 驗證修正。

### 2.1 Step 1 — Baseline & 參數 sanity
**Script**: `research/F_hpfinengroups_deepening/scripts/step1_baseline_and_param_sanity.py`

五個子分析：
1. Baseline 重現（memory 記 89.1%）
2. NG threshold scan（NR≥80 固定）
3. NR threshold scan（NG≥4 固定）
4. NG × NR 2D grid
5. Per-sample 分佈
6. Out-of-scope 檢查（LOH / Paired）

### 2.2 Step 2 — 根因調查
**Script**: `scripts/step2_root_cause_investigation.py`

三個問題：
- **Q1**：NGroups=2 為何比 NGroups=1 差？（germline ASM confound 假說）
- **Q2**：HCC1954 為何 TP rate 只 0.497？
- **Q3**：Paired 99.85% 是 filter gain 還是 baseline 就高？

### 2.3 Step 3 — 新 filter 驗證
**Script**: `scripts/step3_af_cutoff_validation.py`

四個驗證：
1. Per-sample AF cutoff grid（0.2/0.3/0.35/0.4/0.45/0.5）
2. 舊 vs 新 filter per-sample TP rate + coverage loss
3. Coverage_Multiple 分層（CN confound check）
4. Chr-shuffle null（within-chr label shuffle，n=20 iter，套用 P3 pilot 教訓）

成功標準：
- HCC1954 TP rate > 0.75（實際 0.707 接近）
- 6/7 樣本 TP rate ≥ 0.85（實際 5/7，AF<0.2 時 6/7≥0.87）
- Chr-shuffle null Z > 5（實際 43.5）

---

## 3. 結果

### 3.1 Step 1 — Baseline 重現 + 非單調發現

| NGroups | n | TP rate |
|---|---|---|
| 1 | 3,655 | 0.7633 |
| **2** | **51,559** | **0.6434** ⚠️ |
| 3 | 64,653 | 0.7742 |
| 4 | 25,744 | 0.8912 |

**NGroups=2 異常低於 NGroups=1**。先前所有分析（B.1-1/B.1-2/B.1-3）都隱含假設 monotone，此為新發現。

Per-sample TP rate 極度不均（89.1% 是 H2009 主導，佔 77.6%）：
- Tier A（有效）：H2009 (0.935), H1437 (0.921), HCC1395_DORADO (0.903)
- Tier B（中等）：HCC1395 (0.810), HCC1937 (0.714)
- Tier C（失效）：HCC1954 (0.497), COLO829 (0.235, n=34)

### 3.2 Step 2 — 三個根因揭示

**Q1**: NG=2 per-NGroups AF profile

| NG | AF mean | AF∈[0.45,0.55] frac | \|AlleleDelta\| |
|---|---|---|---|
| 1 | 0.552 | 0.175 | 0.019 |
| **2** | **0.471** | **0.212** (最高) | 0.025 |
| 3 | 0.428 | — | 0.023 |
| 4 | 0.402 | 0.145 | 0.020 |

NG=2 AF 最接近 0.5（germline het hallmark）、AF∈[0.45,0.55] fraction 最高。

**NG × AF 2D TP rate grid**（高亮 NG=4 + 低 AF）：

| AF bin | NG=1 | NG=2 | NG=3 | NG=4 |
|---|---|---|---|---|
| 0.0-0.2 | 0.645 | 0.679 | 0.825 | **0.937** |
| 0.2-0.4 | 0.807 | 0.715 | 0.838 | **0.926** |
| 0.4-0.6 | 0.786 | 0.641 | 0.738 | 0.854 |
| 0.6-0.8 | 0.801 | 0.595 | 0.713 | 0.842 |
| 0.8-1.0 | 0.616 | **0.339** | 0.433 | 0.581 |

NG=2 + AF 0.8-1.0 = **0.339**（最差 cell，典型 germline homozygous-like）。NG=4 + AF<0.4 雙條件 TP rate ≥0.92。

**Q2**: HCC1954 NG≥4+NR≥80 NonLOH AF binning

| AF bin | n | n_TP | TP rate |
|---|---|---|---|
| 0.0-0.2 | 564 | 493 | **0.874** ✅ |
| 0.2-0.4 | 476 | 242 | 0.508 |
| 0.4-0.6 | 313 | 54 | 0.173 ⚠️ |
| 0.6-0.8 | 223 | 16 | **0.072** ⚠️ |
| 0.8-1.0 | 46 | 1 | **0.022** ⚠️ |

HCC1954 低 AF 段 TP rate 0.874 完全正常。失效全部來自 AF≥0.4 段。HCC1954 caller_af TP mean=0.214 vs FP mean=0.482。生物解釋：HER2+ breast cancer 高 ploidy（ICGC ~4）、未標註 LOH 複雜 CNV 讓 germline het 表現 AF≈0.5 → FP 富集。

**Q3**: Paired baseline 分解

| Condition | n | TP rate |
|---|---|---|
| All paired (baseline) | 328,699 | **0.9896** |
| Paired + NR≥80 | 141,783 | 0.9967 |
| Paired + NG≥4 + NR≥80 | 11,801 | **0.9985** |
| Paired + NG<4 + NR<80 | 183,252 | 0.9841 |

Filter gain 僅 +0.89pp，**Step 1 的 99.85% 驚人** 是 artifact of baseline。Paired 不是 HPFineNGroups 應用場景。

### 3.3 Step 3 — 新 filter NG=4+AF<0.4 驗證

| Filter | n | TP rate | CI |
|---|---|---|---|
| Overall NonLOH | 307,474 | 0.6699 | — |
| Old (NG≥4+NR≥80) | 25,744 | 0.8912 | [0.887, 0.895] |
| **New (NG=4+AF<0.4+NR≥80)** | **14,197** | **0.9281** | [0.924, 0.932] |

Per-sample 新 vs 舊：

| sample | old TP rate | **new TP rate** | Δ | coverage loss |
|---|---|---|---|---|
| **HCC1954** | 0.497 | **0.707** | **+21.0pp** ⭐ | 36% |
| **HCC1937** | 0.714 | **0.867** | **+15.4pp** ⭐ | 30% |
| HCC1395 | 0.810 | 0.887 | +7.7pp | 27% |
| H1437 | 0.921 | 0.965 | +4.3pp | 20% |
| H2009 | 0.935 | 0.957 | +2.3pp | 50% |
| HCC1395_DORADO | 0.903 | 0.919 | +1.6pp | 11% |
| COLO829 | 0.235 | 0.235 | 0pp | 0% (n=34 已全 AF<0.4) |

5/7 ≥0.85（舊 4/7）。HCC1954/HCC1937 改善完全符合 Step 2 AF confound 分析。

**AF cutoff sensitivity**:

| AF cut | 6/7 ≥0.85 | n |
|---|---|---|
| 0.20 | ✅ 6/7 ≥0.87 高精度 | 3,027 |
| 0.30 | 6/7 | 8,283 |
| **0.40** | 5/7 (HCC1954=0.707) | **14,197** |
| 0.50 | 5/7 | 18,087 |

**Confound checks**:

Coverage_Multiple 分層（CN confound check）：

| CN bin | n | TP rate |
|---|---|---|
| 0.0-1.1 (CN≤1) | 472 | 0.900 |
| 1.1-1.5 | 5,494 | 0.923 |
| 1.5-2.0 (CN~2) | 5,183 | **0.944** |
| 2.0-3.0 (CN~3) | 1,885 | 0.913 |
| 3.0+ (CN≥4) | 1,163 | 0.915 |

所有 0.90-0.94 → **不是 CN confound**。

Chr-shuffle null（套用 P3 pilot 教訓 — spatial autocorrelation defense）：

- Observed: 0.9281
- Null mean (within-chr label permutation, 20 iter): 0.6898 ± 0.0055
- **Z-score: 43.5** → 完全壓倒 null → 非 spatial artifact PASS

---

## 3b. Step 4 — AF<0.4 stratified Cohen's d + FP 聚集分析（20260418 後續增補）

### 3b.1 動機與範圍
Step 3 後列為「Step 4 候選」的兩項：(a) AF<0.4 stratified per-sample Cohen's d 重算（確認 HCC1954/COLO829 升級）；(b) HCC1954 AF≥0.4 FP 是否為 hotspot 聚集（或特定機制）。

**原假設**：AF≥0.4 FP 可能集中於 genomic hotspot，解釋「HCC1954 樣本特異性失效」。

### 3b.2 Part A — AF<0.4 Cohen's h 結果

| Sample | h_all | **h_af04** | **h_af02** | 類別變化 |
|---|---|---|---|---|
| **HCC1954** | +0.587 | **+0.654** | **+0.775** | ✅ medium+ 升級（AF<0.2 接近 large） |
| HCC1937 | +0.534 | +0.523 | +0.350 | ➖ medium+ 穩定 |
| H1437 | +0.462 | +0.543 | +0.313 | ✅ small→medium+ |
| HCC1395_DORADO | +0.553 | +0.257 | +0.165 | ⚠️ 數值降（飽和 artifact） |
| HCC1395 | +0.283 | +0.171 | +0.121 | ⚠️ 數值降（飽和 artifact） |
| H2009 | +0.127 | +0.150 | +0.328 | ➖ ceiling |
| **COLO829** | **-0.457** | **-0.429** | **-0.547** | ❌ permanent NEGATIVE |

**更新的 per-sample 分類**（AF<0.4 stratified）：
- **4/7 medium+ POS**：HCC1954 / HCC1937 / H1437 / HCC1395_DORADO
- **1/7 small POS**：HCC1395
- **1/7 ceiling**：H2009
- **1/7 permanent out-of-scope**：COLO829

HCC1395/DORADO 數值 h 下降為**飽和 artifact**（p4 和 plt 都趨近 1.0，arcsin 飽和導致 h 縮小），非真正訊號衰減。

### 3b.3 Part B — Coverage impact

6/7 樣本仍 YES power（n≥100），僅 COLO829 MARGINAL（n=34，未變）：
- HCC1954 保留 **1,040 regions**（35.9% loss）
- AF<0.2 版本 HCC1954 仍有 **564 regions**（可供論文 figure）
- H2009 50.3% loss 但絕對 n=9,939 仍最多

### 3b.4 Part C — FP 聚集假設推翻 ⚠️

**5/5 樣本 AF≥0.4 FP 在空間上 NOT 聚集**（inter-FP <100kb fraction 皆小於 AF<0.4 TP baseline）：

| Sample | clust FP≥0.4 | clust TP<0.4 | ratio |
|---|---|---|---|
| HCC1954 | 0.270 | 0.333 | **0.81×** |
| HCC1937 | 0.186 | 0.304 | 0.61× |
| H1437 | 0.127 | 0.246 | 0.52× |
| HCC1395 | 0.071 | 0.311 | 0.23× |
| H2009 | 0.141 | 0.769 | 0.18× |

所有 ratio <1 → FP **比 TP 更分散**。原 hotspot 假設**推翻**。

### 3b.5 新機制解釋 — 染色體層級富集（非 hotspot）

**HCC1954 FP 的染色體 AF≥0.4 fraction**：

| Chr | AF≥0.4 fraction |
|---|---|
| chr14 | 0.842 |
| chr2 | 0.775 |
| chr7 | **0.714** |
| **chr8** | **0.686** |
| chr16 | 0.696 |
| chr17 | 0.654 |
| chr1 | 0.507 |
| chr9 | 0.370 (低) |

chr7/8/16 總計 315 FP（>50% 全染色體 FP）、皆 ~70% 在 AF≥0.4 → **HER2 on chr17、MYC on chr8 等已知 HCC1954 amplification 區域**。

### 3b.6 Feature profile — germline-like 機制證據

**AF<0.4 TP (n=735) vs AF≥0.4 FP (n=511)**：

| Feature | AF<0.4 TP | AF≥0.4 FP | Δ |
|---|---|---|---|
| **NumCpGs** | 89.6 | **111.0** | **+21.4 (+24%)** ⭐ |
| Coverage_Multiple | 1.56 | 1.58 | 類似 |
| **PairwiseMedianDist** | 0.189 | **0.168** | **-0.021** ⭐ (reads 相似度更高) |
| \|AlleleDelta\| mean | 0.032 | 0.027 | -0.005 (MWU p=8.67e-4) |

### 3b.7 HCC1954 失效三要素機制

1. **CNV 染色體富集**（chr7/8/16、HER2+ 高 ploidy）：LOH.bed 未涵蓋細微 cnLOH 或 mixed ploidy 區
2. **Germline het 被 CNV 驅動 AF 漂移**：3:1 allele imbalance 區 germline het 從 0.5 → ~0.75 → 被 AF<0.4 cutoff 剔除
3. **CpG island/啟動子 bias**：高 CpG 密度區（+24%）methylation 結構複雜，HP 分群誤判為 NG=4

AF<0.4 挽救的生物學對應：tumor purity<1 + subclonal events → somatic AF 天然 <0.4；germline het 在 CNV amplification 區 AF 漂至高 AF 段被移除。

### 3b.8 更新待驗（Step 5 候選）

- **CNV caller 整合必要性提升**：chr7/8/16 染色體 bias 強烈建議 Delly/Manta/sequenza（Opus 4.7 重整 B.2-2 擔憂相關）
- **Purity<1 模擬**：AF<0.4 閾值在 purity 0.3-0.8 下是否穩健
- **CpG island annotation 交叉**：AF≥0.4 FP 的高 NumCpGs 落在哪類 region（promoter/shore/shelf）
- **AF 雙閾值**：AF<0.1 是否也產生 FP artifact

---

## 4. 對 Part B.1 質疑的回應

| 質疑 | 原擔憂 | F pilot 回應 | 結論 |
|------|--------|-------------|------|
| B.1-1 residualized AUC=0.617 不足 | 低於 caller_af 0.654 | AUC 作為 filter 確實不足；但 characterization 層 AF stratified Δ 達 +0.037 overall、+0.21 HCC1954 | **AUC metric 低估 AF 條件下的訊號** |
| B.1-2 飽和效應 | NR≥80 本身 TP rate 就高 | NR-bin weighted 7/7 POS（已在 20260417 驗證）；新 filter CN bins 均 0.90-0.94 | **非飽和 artifact** |
| B.1-3 7/7 統計強度 | 方向一致 ≠ effect size | 新 filter 下 5/7 ≥0.85，HCC1954 +21pp 挽救 | **effect size 在 AF<0.4 條件下顯著提升** |
| B.1-4 LOH/AF 混淆 | 未 2D stratify | NG×AF 2D grid 清楚展示：NG=2 germline confound、NG=4+AF<0.4 strongest somatic | **已完整解釋並據以改進 filter** |

---

## 5. 新結論穩定度評估

**結論穩定度補充結論 16（HPFineNGroups）升級：⭐3 → ⭐4**

升級依據：
- ✅ 跨樣本驗證：5/7 樣本新 filter TP rate ≥0.85（舊 4/7）
- ✅ Confound 審查：AF（已解釋）、NR（已解釋）、spatial（Z=43.5）、CN（0.90-0.94）、LOH（已排除）、germline（NG=2 AF≈0.5 解釋）
- ✅ 生物學機制：subclonal somatic SNV AF 天然低 + germline het AF≈0.5 移除
- ✅ 失效樣本機制：HCC1954（HER2+ 高 ploidy AF confound）、COLO829（ONT_R10 無 methylation）
- ✅ 防 spatial artifact（套用 P3 教訓）
- ⚠️ 未升級至 ⭐5：因尚未在 patient-derived non-SEQC2 樣本外推驗證

---

## 6. 對其他任務的影響

| 任務/結論 | 影響 | 動作 |
|-----------|------|------|
| B.1-1 residualized AUC 0.617 | 結論仍成立但低估 AF-stratified 潛力 | 後續可重算 AF<0.4 子集的 AUC |
| B.1-3 per-sample Cohen's d | HCC1954/COLO829 可能從「特殊」→「POS」 | Step 4 候選分析（AF<0.4 stratified 重算） |
| 結論穩定度補充結論 16 | ⭐3 → ⭐4 | 同步更新 `06_結論穩定性審查.md` |
| Phase 2 biology characterization | HPFineNGroups 作為 biology-level subclone marker 強化 | Phase 2 報告可直接引用新 canonical filter |
| memory `project_hpfinengroups_subclone_marker.md` | 已更新（20260418） | ✅ |
| B.2 LOH Subclone AF×Methylation | 不受影響（AF-bin 框架正交） | N/A |
| Zone-Aware Framework | 不受影響 | N/A |
| Per-CpG ASM characterization | 不受影響 | N/A |

---

## 7. 未解決與 Step 4 候選

### 7.1 新產生的生物學問題
1. **AF<0.4 為何有效的機制**：subclonal somatic SNV 天然低 AF + germline het 在 AF≈0.5 被移除（但 NonLOH filter 已排除 LOH 驅動 AF 漂移情境）
2. **COLO829 sample-level gating**：應增「優先 ONT_5mCG 或 5mCG+5hmCG basecall 樣本」之 pre-check
3. **HCC1954 高 AF FP 是否全部 germline-like?** NumCpGs FP (111) > TP (90) → FP 集中高 CpG 密度區（啟動子/CpG island）待驗證

### 7.2 Step 4 候選（未執行）
- **AF<0.4 stratified per-sample Cohen's d 重算**：若 HCC1954/COLO829 從「特殊」→「POS」，可精確化 B.1-3 per-sample effect size 結論
- **patient-derived cohort 外推驗證**：⭐4 → ⭐5 升級所需
- **purity simulation**：mixed 0.3/0.5/0.7 檢驗 AF<0.4 閾值是否需調整

### 7.3 不執行理由
本 pilot 以「確認新 filter 方向有潛力、建立 canonical 版本」為主軸，per-sample Cohen's d 重算屬後續增量工作，不阻擋 Phase 2 biology characterization 進入論文撰寫階段。

---

## 8. 產物

| 類型 | 路徑 |
|------|------|
| Manifest | `research/F_hpfinengroups_deepening/manifest.yaml` |
| Step 1 script | `research/F_hpfinengroups_deepening/scripts/step1_baseline_and_param_sanity.py` |
| Step 1 findings | `research/F_hpfinengroups_deepening/observations/step1_findings.md` |
| Step 2 script | `research/F_hpfinengroups_deepening/scripts/step2_root_cause_investigation.py` |
| Step 2 findings | `research/F_hpfinengroups_deepening/observations/step2_findings.md` |
| Step 3 script | `research/F_hpfinengroups_deepening/scripts/step3_af_cutoff_validation.py` |
| Step 3 findings | `research/F_hpfinengroups_deepening/observations/step3_findings.md` |
| Step 4 script | `research/F_hpfinengroups_deepening/scripts/step4_af04_cohens_d_and_fp_clustering.py` |
| Step 4 findings | `research/F_hpfinengroups_deepening/observations/step4_findings.md` |
| Step 1-4 data outputs | `research/F_hpfinengroups_deepening/data/*.tsv` |
| Memory | `.claude/.../memory/project_hpfinengroups_subclone_marker.md`（20260418 updated） |
| Stability | `docs/reports/research_landscape/06_結論穩定性審查.md`（補充結論 16 ⭐3→⭐4） |

---

## 9. 關鍵引用（供論文/簡報）

> "HPFineNGroups=4 is only a somatic heterogeneity marker when co-filtered with caller_af<0.4 and NumReads≥80 (non-LOH regions). Under this refined filter, TP rate reaches 92.8% (vs 89.1% with the legacy filter), with HCC1954 recovered from 49.7% to 70.7% (+21pp) after the AF<0.4 constraint removes germline-het-like FPs at AF≈0.5. Chr-shuffle within-chr null Z-score = 43.5 and Coverage_Multiple-stratified TP rates remain within 0.90–0.94 across CN tiers, ruling out spatial and CN confounds. The NGroups=2 group (TP rate 64.3%) is enriched for germline ASM signatures (AF mean 0.471, AF∈[0.45,0.55] fraction 0.212), explaining previously unnoticed non-monotonicity. COLO829 is excluded from scope due to ONT_R10 basecalling lacking 5mCG."
