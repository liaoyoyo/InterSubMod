---
title: Z3 內部特徵與 AF Germline Band 分離探索
date: 2026-04-18
status: NEGATIVE (main) + CONDITIONAL (HCC1954-specific)
authors: InterSubMod Research
tags: [zone-aware, Z3, HCC1954, germline-escape, LOH, TO-mode, NEGATIVE]
related:
  - docs/reports/research_landscape/08_Zone_Aware.md
  - research/F_hpfinengroups_deepening/
  - docs/experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md
trust_level: in_progress
---

# Z3 內部特徵與 AF Germline Band 分離探索（Step 1/2.5/3）

## 摘要（TL;DR）

**主結論：NEGATIVE**。TO mode Z3（LOH + extreme AF + HPFineNGroups≤1）內部不存在 ≥3 樣本可用的二階區分特徵；AF∈[0.4, 0.6] germline band 的 FP 無法由 CN × HPFineNGroups 聯合分層一致性地分離。

**次結論：HCC1954 單一樣本例外（CONDITIONAL）**。HCC1954 Z3 FP 集中在 **chr5 (43%) / chr8 (29%) / chr17 (13%)**（覆蓋 HER2 與 MYC amplicon 區域），FP 的 NumReads（55 vs TP 37，p=4.7e-9）與 Coverage_Multiple（0.73 vs 0.49，p=4.7e-9）顯著偏高，支持這些 FP 為 **CNV amplicon artifact** 而非 somatic LOH。這與 F pilot `NG=4+AF<0.4+NR≥80` filter 已偵測到的 HCC1954 特殊性一致，沒有產生新的跨樣本 F1 提升路徑。

**對 5 個研究目標的影響**：
- 目標 1-4：**零衝擊**，本次探索為 zone-aware 輔助驗證，不改變 per-CpG ASM / clone structure / two-hit order / normal reference 既有路線
- 目標 5（F1）：**HCC1954 sample-specific 已被 F pilot 覆蓋**，無新增全域路線

---

## 一、研究動機

先前 Zone-Aware Framework 驗證確認 Z1 vs Z3 TP rate 差異真實（Z1=0.965 / Z3=0.608 in TO），但 QS 調整 NEGATIVE。新問題：

- Z3 TO TP rate 跨樣本 σ=0.28（0.050–0.868），代表 **樣本特異性 FP 富集**
- Zone-Aware 層面已結案（characterization only），但「Z3 內是否還有未探索的二階區分信號」未測試
- 用戶另要求驗證 H-Z3d：AF∈[0.4, 0.6] 的 FP 可能是 **germline het 逃逸** 或 **CNV artifact**，可用 CN × NGroups 分離

**四項假說**：
- **H-Z3a**：HCC1954 Z3 低 TP rate 主因是 AF confound，AF bin 內分層應揭露訊號
- **H-Z3b**：HP-free 甲基化特徵（PairwiseMedianDist, HeuristicScore 等）避開 self-phasing artifact，應具獨立區分力
- **H-Z3c**：per-sample 特徵組合可能比 global 閾值有效
- **H-Z3d**：AF∈[0.4, 0.6] FP 在 CN≈2 + NGroups=1 條件下應為 germline het

---

## 二、Z3 樣本間變異基線（摘要）

| 樣本 | Z3 n | Z3 TP rate | 備註 |
|------|------|-----------|------|
| H2009 | 8,772 | 0.868 | 高純度 / 高覆蓋 |
| H1437 | 11,342 | 0.730 | 正常 |
| HCC1395 | 6,661 | 0.551 | CEPH cell line |
| HCC1395_DORADO | 9,186 | 0.535 | 同 HCC1395 另一 basecaller |
| HCC1937 | 5,258 | 0.244 | BRCA1-/-  germline-heavy |
| **HCC1954** | **2,030** | **0.050** | **HER2+ amplicon cell line（最極端）** |
| COLO829 | 9,345 | 0.651 | melanoma |

HCC1954 作為最極端低 TP rate，是本研究主要焦點。

---

## 三、Step 1：Z3 內部特徵 AUC 掃描（NEGATIVE）

**方法**：對 7 樣本 × 12 候選特徵，在 TO Z3 內計算 |AUC|（TP vs FP），判定 gate =「任一特徵在 ≥3 樣本 |AUC|≥0.60」。

**全域結果**：**NO 特徵達到 gate 條件**。

### 單樣本最高 |AUC|（僅作 heatmap 對照）

| 樣本 | Top feature | |AUC| | 次 feature | |AUC| |
|------|------------|-----|-----------|-----|
| H2009 | NumCpGs | 0.607 | PairwiseMedianDist | 0.607 |
| H1437 | NumCpGs | 0.590 | — | — |
| HCC1395 | NumReads / Coverage_Multiple | 0.547 | — | — |
| HCC1395_DORADO | NumReads / Coverage_Multiple | 0.565 | — | — |
| HCC1937 | NumReads / Coverage_Multiple | 0.576 | — | — |
| **HCC1954** | **NumReads / Coverage_Multiple** | **0.673** | NumCpGs | 0.581 |
| COLO829 | NumCpGs | 0.585 | PairwiseMedianDist | 0.580 |

**關鍵觀察**：
- 所有甲基化特徵（AlleleDelta, HPFineF, HPMergedDelta, LabelAllelePermanovaF）幾乎全 ~0.50（在 HPFineNGroups≤1 條件下，HP-dependent 特徵結構性退化）
- HCC1954 單一樣本 NumReads/Coverage_Multiple |AUC|=0.673 為唯一顯著信號
- nhp_ratio（自構 self-phasing proxy）全樣本 ~0.50，證實在 NG≤1 約束下 self-phasing 訊號已被 collapsing

**判定**：NEGATIVE gate → 進入 Step 2.5（不進入 Step 2 AF 分層）。

**產出**：
- `research/z3_internal_feature_exploration/data/z3_feature_auc.tsv`
- `research/z3_internal_feature_exploration/figures/z3_feature_auc_heatmap.png`

---

## 四、Step 2.5：AF Germline Band 分層（NEGATIVE）

**方法**：全樣本 TO mode，`caller_af∈[0.4, 0.6]` 子集，按 CN（Coverage_Multiple <0.75 / 0.75-1.5 / ≥1.5）× HPFineNGroups（≤1 / 2-3 / ≥4）聯合分層，Fisher test 比較 `CN≈2+NG≤1` vs `CN≈2+NG≥4`。

**樣本基線（AF∈[0.4, 0.6] 全體 TP rate）**：

| 樣本 | AF[0.4,0.6] n | TP rate |
|------|--------------|---------|
| H2009 | 53,029 | 0.917 |
| H1437 | 13,655 | 0.762 |
| HCC1395 | 10,595 | 0.725 |
| HCC1395_DORADO | 10,888 | 0.727 |
| HCC1937 | 5,612 | 0.497 |
| **HCC1954** | 17,315 | **0.145** |
| COLO829 | 18,806 | 0.669 |

### Critical cell: CN≈2 + NG≤1 + AF∈[0.4, 0.6]

| 樣本 | n | TP rate | loh_pct | 是否符合 germline pattern (TP<0.3) |
|------|---|---------|---------|-----------------------------------|
| H2009 | 1,786 | 0.945 | 0.25 | ✗ |
| H1437 | 383 | 0.799 | 0.15 | ✗ |
| HCC1395 | 127 | 0.756 | 0.43 | ✗ |
| HCC1395_DORADO | 342 | 0.766 | 0.35 | ✗ |
| HCC1937 | 89 | 0.629 | 0.40 | ✗ |
| **HCC1954** | **254** | **0.146** | **0.03** | **✓** |
| COLO829 | 10 | 0.700 | 0.00 | ✗ |

**Fisher test（CN≈2+NG≤1 vs CN≈2+NG≥4 在 AF[0.4,0.6]）**：
- 所有樣本 delta 介於 -0.23 至 +0.08，無一達顯著（最低 p=0.118，HCC1395_DORADO）
- 預期的「NG=1 低 TP / NG≥4 高 TP（germline vs subclonal）」pattern 未成立

**AlleleDelta 在 germline band 的 AUC**：介於 0.506–0.578，無分離力

**判定**：**NEGATIVE（僅 1/7 樣本符合）**。H-Z3d 假說（CN × NGroups 可一致性分離 germline FP）不成立。AF∈[0.4, 0.6] 的 FP 在多數樣本（6/7）**並非由 germline het 主導**，而是其他異質機制。

**產出**：
- `data/z3_af_germline_band_sample_summary.tsv`
- `data/z3_af_germline_band_cn_ng.tsv`
- `data/z3_af_germline_band_fisher.tsv`
- `data/z3_af_germline_band_alleledelta_auc.tsv`
- `figures/z3_af_germline_band_heatmap.png`

---

## 五、Step 3：HCC1954 vs HCC1395 Z3 機制對比（POSITIVE 特殊觀察）

兩樣本同為 breast cancer cell line，Z3 TP rate 差 11×（HCC1954=0.050 vs HCC1395=0.551）。對比揭示**根本機制差異**。

### 5.1 染色體分佈（關鍵差異）

| HCC1954 top chr | n | % | TP rate | HCC1395 top chr | n | % | TP rate |
|----------------|---|---|---------|-----------------|---|---|---------|
| **chr5** | 875 | **43.1%** | 0.024 | chr3 | 712 | 10.7% | 0.542 |
| **chr8** | 583 | **28.7%** | 0.082 | chr2 | 591 | 8.9% | 0.541 |
| **chr17** | 256 | **12.6%** | 0.031 | chr11 | 585 | 8.8% | 0.564 |
| chr12 | 178 | 8.8% | 0.118 | chr5 | 540 | 8.1% | 0.561 |
| chr19 | 50 | 2.5% | 0.000 | chr12 | 538 | 8.1% | 0.563 |

- **HCC1954**：top-4 染色體佔 93%，高度集中（chr17=HER2 amplicon, chr8=MYC amplicon, chr5 亦 HER2+ 高頻 CNV 區）；每染色體 TP rate ~0.02-0.12（極低）
- **HCC1395**：top-5 染色體佔 ~45%，分佈均勻；每染色體 TP rate ~0.54-0.56（接近 overall 0.551）

**解讀**：HCC1954 Z3 低 TP rate 是 **amplicon-region 地域性 FP 富集** 造成，不是 Z3 機制失效。HCC1395 證明 Z3 在 CNV-clean 樣本中能定位真 somatic LOH。

### 5.2 Coverage / NumReads 差異

| 樣本 | NumReads TP median | NumReads FP median | Mann-Whitney p |
|------|-------------------|-------------------|----------------|
| **HCC1954** | **37** | **55** | **4.7e-9** |
| HCC1395 | 45 | 47 | 2.8e-11 (但效應小) |

| 樣本 | Coverage_Multiple TP median | FP median | p |
|------|----------------------------|-----------|---|
| **HCC1954** | **0.49** | **0.73** | **4.7e-9** |
| HCC1395 | 0.60 | 0.63 | 2.8e-11 (但效應小) |

HCC1954 FP 顯著高覆蓋（代表 CNV amplicon），TP 低覆蓋（可能是正常倍體區的真 somatic）。HCC1395 TP/FP 覆蓋幾乎重疊（兩者同源於一般 LOH）。

### 5.3 PairwiseMedianDist 與 NumCpGs

兩樣本 TP/FP 分佈高度重疊（Mann-Whitney p 雖顯著但 median 差 ~5-8，無實用 AUC），不支持 germline-like reads 或 CpG island bias 假設。

**產出**：
- `data/z3_mechanism_chr_dist.tsv`
- `data/z3_mechanism_feature_dist.tsv`
- `figures/z3_hcc1954_vs_hcc1395_mechanism.png`
- `figures/z3_chr_distribution_contrast.png`

---

## 六、最終判定

| 假說 | 結果 | 信心 |
|------|------|------|
| H-Z3a（AF confound → AF bin 內分層） | **NEGATIVE** — Step 1 caller_af 全樣本 |AUC|<0.55 | 高 |
| H-Z3b（HP-free 特徵避 self-phasing） | **NEGATIVE** — PairwiseMedianDist/HeuristicScore 除 H2009/COLO829 邊緣外全 ~0.50 | 高 |
| H-Z3c（per-sample 組合） | **NEGATIVE** — 僅 HCC1954 一樣本有單特徵 |AUC|>0.60 | 高 |
| H-Z3d（CN×NG 分離 germline） | **NEGATIVE** — 1/7 樣本符合，且該樣本已由 F pilot 覆蓋 | 高 |
| 額外：HCC1954 amplicon-driven FP | **POSITIVE（已知）** — Step 3 確認 chr5/8/17 地域性 CNV FP | 高 |

**核心結論**：Z3 的跨樣本變異來源不是「Z3 內缺少二階可分特徵」，而是「**一個離群樣本（HCC1954）以特定 amplicon 區域 CNV artifact 主導 FP**」。Z3 本身對 6/7 樣本仍是合理的高 TP rate 區域（0.244–0.868）。

---

## 七、對 5 研究目標的影響

| 目標 | 影響 | 說明 |
|------|------|------|
| 1. per-CpG ASM | 零 | Z3 作為 non-ASM 對照的定位不變 |
| 2. clone structure | 零 | Z1 vs Z3 的 NGroups 生物學對比仍有效 |
| 3. two-hit order | 零 | 與本次 zone-aware 路線無交集 |
| 4. normal reference | 零 | 與 Z3 無關 |
| 5. F1 | HCC1954-only CONDITIONAL | 已被 F pilot 的 `NG=4+AF<0.4+NR≥80` canonical filter 覆蓋，無新增路徑 |

**Zone-Aware Framework 總體評估維持不變**：僅作為 **annotation / characterization** 使用，不可作為 QS 調整或 filter design 的基礎。HCC1954 的 amplicon-FP 機制需要在 F1 層由 sample-specific filter（如 high-NR + AF<0.4 黑名單）處理。

---

## 八、Next Steps

- **無立即後續**：本研究為探索性 pilot，NEGATIVE 結論已達成決策目的
- **報告升級**：保留 `in_progress`，若 HCC1954 amplicon 機制納入論文 case study 再升級
- **HCC1954 後續方向**（若用戶決定推進）：可針對 chr5/chr8/chr17 區域做 CN-region blacklist，但需配合 SEQC2 F1 panel 重新評估（預期增益 <+0.005）

---

## 九、產出檔案清單

```
research/z3_internal_feature_exploration/
├── scripts/
│   ├── step1_z3_feature_auc_table.py
│   ├── step2_5_af_germline_band.py
│   └── step3_z3_mechanism_contrast.py
├── data/
│   ├── z3_feature_auc.tsv                        (Step 1)
│   ├── z3_af_germline_band_sample_summary.tsv    (Step 2.5)
│   ├── z3_af_germline_band_cn_ng.tsv             (Step 2.5)
│   ├── z3_af_germline_band_fisher.tsv            (Step 2.5)
│   ├── z3_af_germline_band_alleledelta_auc.tsv   (Step 2.5)
│   ├── z3_mechanism_chr_dist.tsv                 (Step 3)
│   └── z3_mechanism_feature_dist.tsv             (Step 3)
└── figures/
    ├── z3_feature_auc_heatmap.png                (Step 1)
    ├── z3_af_germline_band_heatmap.png           (Step 2.5)
    ├── z3_hcc1954_vs_hcc1395_mechanism.png       (Step 3 violin)
    └── z3_chr_distribution_contrast.png          (Step 3 chr bar)
```

---

## 十、風險與限制

1. HCC1954 Z3 n=2030，分 chr 後 chr5 n=875 仍有足夠 power；Coverage_Multiple AUC=0.673 在 chr-shuffle null 下預期 Z>3，但**未執行** chr-shuffle 驗證（單樣本觀察性結果，不升級結論）
2. HCC1937 Z3 n=5258 + TP rate=0.244 被納入主結論 NEGATIVE 群（6/7 不符合 germline pattern），但其 BRCA1-/- 背景可能有獨立機制，未進一步拆解
3. Step 2.5 Fisher test 在 COLO829 等小 n cell 不可解釋（n=10）；聚集於高 n 樣本（HCC1395/1395_DORADO/H2009）均無顯著

---

**文件狀態**：in_progress
**不升級 validated 原因**：探索性 pilot，結論對既有路線無變動；HCC1954 機制觀察性強但未做 null validation
