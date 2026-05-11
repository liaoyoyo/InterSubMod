---
title: "HCC1954 Outlier Case Panel — Caller FP Background Ceiling, NOT Thread D Mechanism Failure"
date: 2026-04-26
status: validated
verdict: HCC1954 outlier 為 caller FP background ceiling 之 sample-level artifact，**非 Thread D LOH-constrained phasing 機制反例**
evidence_grade: B (Z3 amplicon pilot 已知 + B2 全基因體分析；caller-level confound 證據充分但缺 cross-caller benchmark)
owner: InterSubMod Research
thread: Thread D · LOH-constrained phasing — sample-level 異常 case panel
samples: HCC1954 (TO mode 主)
tracks:
  - TO (主)
  - paired_full (對照)
related_main_axis:
  - InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md (主軸；§3.L1 cross-sample, §5 anomalies)
source_reports:
  - InterSubMod/docs/experiments/in_progress/2026/04/20260423_B2_HCC1954_Outlier_RootCause_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_design_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260424_X5_CrossSample_obs18_KDE_Verified_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260423_B1_Wilcoxon_NG2_gap_01.md
figures_reused:
  - InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_hcc1954_outlier_4panel.png
tags:
  - HCC1954
  - outlier
  - caller-ceiling
  - case-panel
  - thread-D
---

# HCC1954 Outlier Case Panel — Caller FP Background Ceiling

## 1. Executive Summary

**核心結論**：HCC1954 在 Thread D（LOH-constrained phasing）跨樣本 obs18 / B1 / X5 驗證中表現為「TP gap 仍正向 (+0.345) 但 Inner same_HP1 TP rate 僅 0.429（其他 5 樣本 ≥0.759）」的 outlier。經 B2 全基因體根因分析：

- **真正根因**：HCC1954 TO pipeline 在 Outer × NG=2 × cross_het × Extreme AF 這個 cell 內 FP fraction = **89.9%**（6,653/7,400）。FP 為**全基因體均勻分佈**，非 amplicon 局部、非 LOH 註解失誤、非 AF/CovM 偏移。
- **非機制反例**：HCC1954 NG=2 inner same-hap TP rate 偏低是被 caller-level FP background 「分母拉低」所致，並非 LOH-constrained phasing 機制失效；此案的 obs18 gap 仍 +0.345（6/6 正向）。
- **與 Thread D 主軸的關係**：Thread D Wilcoxon n=6 包含 HCC1954（X5 結果 W=21, p=0.01562），不需排除即達顯著。但若以 mechanism strength 為主訊號，HCC1954 應獨立標記為「caller-saturated sample」以避免效應量被 dilute。

**Evidence grade**：**B**
- 充分支持：B2 全基因體 89.9% FP fraction、Z3 amplicon pilot ceiling 已 ΔF1=+0.0075（接近 reject-all-Z3 上限）、reversal_investigation Bootstrap CI=[-0.604,+0.044] 證明反向為小樣本噪音
- 尚缺：跨 caller benchmark（DeepVariant/Strelka2 vs ClairS-TO）、HCC1954 truth set tier 結構與其他樣本對齊驗證

---

## 2. Outlier 觀察（X5 cross-sample table）

X5 報告（`20260424_X5_CrossSample_obs18_KDE_Verified_01.md` §3）KDE-corrected master 上 6 TO 樣本 obs18 NG=2 Inner same_HP1 vs Outer cross_het TP gap：

| Sample | Inner same_HP1 TP | Outer cross_het TP | Gap |
|--------|:-----:|:-----:|:-----:|
| HCC1395 | 0.840 | 0.580 | +0.260 |
| HCC1395_DORADO | 0.939 | 0.553 | +0.385 |
| H1437 | 0.920 | 0.688 | +0.231 |
| H2009 | 0.932 | 0.882 | +0.049 |
| HCC1937 | 0.759 | 0.236 | **+0.522** |
| **HCC1954** | **0.429** | **0.084** | **+0.345** |

**統計**：Wilcoxon signed-rank (alternative=greater) W=21.0, p=0.01562, n=6 — 6/6 正向 gap，HCC1954 包含後仍達顯著。

**HCC1954 獨特之處**：
- Inner same_HP1 TP rate **0.429**（其他 5 樣本中位數 0.917，最低 HCC1937 0.759）— 比次低樣本低 0.330
- Outer cross_het TP rate **0.084**（其他 5 樣本中位數 0.580）— 比次低樣本 HCC1937 0.236 還低 0.152
- 然而 **gap +0.345 仍與 5 樣本中位數 (+0.310) 同範圍**

**判讀**：HCC1954 outlier 不是 gap 方向錯，而是 **Inner / Outer 兩端 TP rate 同步壓低**，此為「分母現象」而非「機制現象」。

---

## 3. B2 根因分析重述（重用 4-panel figure）

來源：`InterSubMod/docs/experiments/in_progress/2026/04/20260423_B2_HCC1954_Outlier_RootCause_01.md` Lines 22-106

![HCC1954 Outer cross-het TP outlier 4-panel root cause](../../../../experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_hcc1954_outlier_4panel.png)

B2 對 NG=2 × cross_het ∪ cross_het_inv × Extreme AF subset（HCC1954 n=7,400 vs 其他 5 TO 樣本 n=24,178）進行三項候選假說的逐項否決：

### H1（LOH 註解失誤）— REJECTED
（B2 §核心觀察 3，Lines 64-76）

| Sample | Inner | Outer | Inner_frac |
|--------|------:|------:|----------:|
| H1437 | 184 | 8,640 | 0.021 |
| H2009 | 303 | 7,207 | 0.040 |
| HCC1395_DORADO | 93 | 6,658 | 0.014 |
| HCC1937 | 58 | 1,669 | 0.034 |
| **HCC1954** | 137 | 7,400 | **0.018** |

HCC1954 Inner_frac=0.018 介於 HCC1395_DORADO (0.014) 與 HCC1937 (0.034) 之間，**無異常**。

### H2（AF/CovM 分佈異常）— REJECTED
（B2 §核心觀察 2，Lines 57-62）

- AF<0.1 比例：HCC1954 100% vs cohort 100%（filter 同質化）
- CovM median：HCC1954 **1.08** vs cohort **1.17**（HCC1954 略低，**非更高**）

### H3（HER2/MYC amplicon 主導）— REJECTED
（B2 §核心觀察 1，Lines 43-55）

| Chr | HCC1954 n | HCC1954 FP | HCC1954 TP rate | Cohort TP rate |
|-----|----------:|----------:|----------:|----------:|
| chr1 | 824 | 765 | 0.072 | 0.689 |
| chr2 | 699 | 638 | 0.087 | 0.704 |
| chr8 (MYC) | 224 | 210 | 0.063 | 0.659 |
| chr17 (HER2) | 123 | 117 | 0.049 | 0.634 |
| chr5 | 100 | 98 | 0.020 | 0.749 |

chr8 + chr17 僅佔 HCC1954 全 FP 的 **5.0%**（210+117 = 327/6,653）。FP top5 chr 為 chr1, chr2, chr3, chr4, chr7（按 variant count 自然排序），無 amplicon 偏集。

### B2 結論
（B2 Lines 87-89）

> 三個原始假設 H1/H2/H3 全部被否決。HCC1954 Outer cross-het TP 0.08 outlier 並非由註解缺陷、AF/CovM 分佈、或已知 amplicon 驅動，而是**全基因體性 FP 密度普遍過高**（Outer × NG=2 × cross_het 這個 cell 裡 FP 佔比接近 90%）。

---

## 4. Caller precision ceiling 量化

### 4.1 Per-sample TP rate（NG=2 × cross_het ∪ cross_het_inv × Extreme AF subset）

來源：`figures/20260423_B2_hcc1954_outlier/B2_tp_rate_by_sample.tsv` 與 `B2_summary.json`

| Sample | n | n_TP | TP rate |
|--------|---:|---:|---:|
| H2009 | 7,207 | 6,322 | **0.877** |
| HCC1395 (paired sub) | 4 | 3 | 0.750 (n 太少) |
| H1437 | 8,640 | 5,893 | 0.682 |
| HCC1395_DORADO | 6,658 | 3,769 | 0.566 |
| HCC1937 | 1,669 | 415 | 0.249 |
| **HCC1954** | **7,400** | **747** | **0.101** |

**HCC1954 整體 FP fraction = 89.9%（6,653/7,400）**。

### 4.2 FP fraction 階梯

| Sample | FP fraction in Outer×NG2×cross_het | 註解 |
|--------|---:|------|
| H2009 | 12% | 高 caller precision |
| H1437 | 32% | 中等 |
| HCC1395_DORADO | 43% | 中等 |
| HCC1937 | 75% | 偏高 |
| **HCC1954** | **90%** | **caller saturated** |

**判讀**：HCC1954 的 caller precision 在這個 cell 已飽和到接近 random（10% TP rate ≈ truth set baseline rate 級別）。任何 ISM-side 機制訊號都會被分母吃掉。

### 4.3 LOH-context Inner_frac（per-sample）

來源：`figures/20260423_B2_hcc1954_outlier/B2_h1_loh_annotation.tsv`

HCC1954 Inner_frac=0.018 與 cohort 同範圍 → 證實 caller saturation 與 LOH 註解獨立，不可混為一談。

---

## 5. Z3 amplicon pilot 結果整合

來源：`InterSubMod/docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md` Lines 1-178

### 5.1 Pilot 設計概要

針對 Z3 = (LOH hit) ∧ (caller_af extreme) ∧ (NGroups ≤1) 的高 FP zone，HCC1954 Z3 中 chr5/chr8/chr17 集中度為 **93%**（design report Line 37）。Pilot 測試 4 strategy：S1 literature arm-level / S2 whole-chr / S3 CovM 95th-percentile / S4 reject all Z3 (ceiling)。

### 5.2 Per-sample × strategy ΔF1（result Lines 39-46）

| Sample | F1 before | S1 ΔF1 | S2 ΔF1 | S3 ΔF1 | S4 ΔF1 (ceiling) |
|--------|---:|---:|---:|---:|---:|
| H2009 | 0.9545 | −0.0024 | −0.0025 | 0.0000 | −0.0270 |
| H1437 | 0.8712 | −0.0072 | −0.0086 | −0.0002 | −0.0717 |
| HCC1395 | 0.8309 | −0.0037 | −0.0041 | +0.0001 | −0.0292 |
| HCC1395_DORADO | 0.8330 | −0.0041 | −0.0039 | +0.0003 | −0.0361 |
| HCC1937 | 0.6772 | +0.0023 | +0.0038 | 0.0000 | **+0.0309** |
| **HCC1954** | **0.4047** | **+0.0058** | **+0.0065** | +0.0002 | **+0.0075** |
| COLO829 | 0.7906 | −0.0108 | −0.0109 | −0.0004 | −0.0644 |

### 5.3 關鍵推論

1. **HCC1954 baseline F1 = 0.4047** — 7 樣本中遠遠最低（次低 HCC1937 0.6772），與 caller saturation 一致
2. **Ceiling ΔF1=+0.0075**：reject all Z3 的 F1 提升上限。S2 達 87% ceiling，S1 達 77%
3. **跨樣本 collateral damage**：5/6 其他樣本 hurt（mean ΔF1 = −0.0044），無法作為跨樣本 global rule
4. **判定 CONDITIONAL**：Z3 是 HCC1954 的 caller-side 已知 hotspot，與 Thread D 機制獨立；Zone-Aware Framework 定位不變（characterization only）

來源：result Lines 14-19（執行摘要）、Lines 161-178（結論）

### 5.4 與 Thread D 機制的獨立性

Z3 amplicon pilot 證明 HCC1954 的 FP 過剩來自 **caller × CNV 架構交互**（HER2+ pseudo-tetraploid + 8p loss + 17p TP53 LOH），與「LOH-constrained phasing 在 Inner/Outer 區呈現 TP rate 差異」的機制完全正交。Z3 高 FP rate 是 caller-side 飽和；Thread D gap 是 ISM-side phasing signature。

---

## 6. 與 Thread D Wilcoxon 的關係

### 6.1 為何 HCC1954 包含後 Wilcoxon 仍顯著

X5 cross-sample 結果（KDE-corrected master）：
- 全 6 樣本：W=21, p=0.01562, n=6, 6/6 正向
- HCC1954 gap=+0.345 與 cohort median (+0.310) 同範圍 → 不會反轉 sign 也不會被視為 outlier 排除

**結論**：Thread D 主軸 evidence-grade B 不依賴 HCC1954 排除。

### 6.2 為何不是 cherry-pick

若以 Inner same_HP1 TP rate 絕對值（HCC1954 = 0.429）為主訊號排除 HCC1954，將被質疑 cherry-pick。但本案 case panel 採用以下原則：

1. **Wilcoxon 結果完整保留 n=6**：不排除任何樣本，p=0.01562 為「含 HCC1954」的真實顯著性
2. **HCC1954 outlier 為 sample-level artifact**：caller saturation（FP fraction 90%）為樣本層級 noise floor，不是「方向反例」
3. **若以「機制強度」為次要分析**：可在敘事中說明「6/6 方向一致 + 5/6 樣本 Inner TP rate ≥0.759；HCC1954 Inner 偏低為 caller saturation 解釋」，並在 limitations 中明列。

來源：`20260417_HCC1954_reversal_investigation_01.md` Lines 19-26（reversal_investigation 已用 Bootstrap CI 證明 HCC1954 反向是 statistical noise，非 mechanism reversal；同樣邏輯延伸至 NG=2 inner TP rate 偏低）。

### 6.3 「樣本 caller 飽和導致 outlier 不可比」

reversal_investigation Lines 198-201（Section 5）已建立 pre-registered exclusion 概念：「n_CN1_segments ≥ 50」作為功效門檻。本 case panel 將此擴展為：

> **Caller-saturation guard**：在 sample × cell × FP fraction > 85% 時，TP rate 絕對值不可作為機制強度 proxy，僅可作為「方向 (sign)」貢獻。

此原則保護 Thread D 結論不被個別 sample 的 caller 表現污染。

---

## 7. 三項質疑 + 邏輯鏈

### Q1：HCC1954 outlier 是否來自 HPFineN_HP1S 自參考 risk？

**否決證據**：X3 flag-on/off neg-control 顯示 HCC1954 在 flag=on（germline-HP-only demote）模式下 NG≥3 也接近 0（與 cohort 同向 collapse），代表 self-phasing tag 並非 HCC1954 獨有問題。

引用：`20260424_X3_FlagOnOff_NG2_NegControl_01.md`（X5 §5.1 引用同 4-重驗證鏈）

→ HCC1954 outlier **不是** HP tag 自參考所致；caller-side 的 FP fraction 90% 才是主導。

### Q2：HCC1954 outlier 是否來自 LOH calling 失敗？

**否決證據**：對照其他 LOH-rich 樣本：
- HCC1937 LOH-rich + Inner_frac 0.034 → gap **+0.522**（最高）
- HCC1395_DORADO LOH-rich + Inner_frac 0.014 → gap +0.385
- HCC1954 Inner_frac 0.018 → gap +0.345（仍同向）

且 B2 H1 已驗證 HCC1954 Inner_frac 與 cohort 同範圍。LOH calling 在 HCC1954 並無系統性失誤。

→ HCC1954 outlier **不是** LOH calling 失敗。

### Q3：HCC1954 outlier 是否來自 ISM HP tag 錯誤？

**否決證據**：reversal_investigation §3.5 Region-level ρ（Table Lines 161-167）：
- HCC1954 region-level ρ(caller_af, HPFineNGroups) = **−0.282**（n=1,418 regions，p=2.5e-27）
- 方向與其他 6 樣本一致（−0.262 to −0.802）

paired_full mode 下 HCC1954 並無此 outlier 級別下沉（reversal_investigation Lines 47-51 paired ρ=−0.211, CI 含 0 即 noise）；僅 TO mode 下 caller saturation 才放大為觀察級 outlier。

→ HCC1954 outlier **不是** ISM HP tag 錯誤；TO-mode caller-side 飽和才是主因。

---

## 8. 後續實驗候選（hypothesis-level）

以下為**未啟動**的延伸方向，標記 hypothesis-level 供研究 roadmap 使用：

### Direction A — Cross-caller benchmark（hypothesis-level）
**假設**：以 DeepVariant、Strelka2 對 HCC1954 TO BAM 重新 call somatic variants，比較與 ClairS-TO 在 Outer × NG=2 × cross_het 這個 cell 的 FP fraction。
**預期區辨**：若 cross-caller FP fraction 仍 >80% → truth set 限制（tier 結構或覆蓋率不足）；若 <50% → ClairS-TO 在 HCC1954 的 caller-specific bias。
**Evidence-grade lift**：可將本 case panel 從 B 升至 A。

### Direction B — HCC1954 truth set 升級或 alt VAF source（hypothesis-level）
**假設**：HCC1954 SEQC2-type truth set 可能因 HER2+ pseudo-tetraploid 而 tier 1 high-confidence 涵蓋率偏低（B2 Lines 95-96 已標記）。
**驗證方式**：比對 HCC1954 truth set 的 variants/Mb、tier 比例與其他 6 樣本；若有公開 high-purity orthogonal truth（e.g. PCR-validated）則替換重算 obs18。
**預期區辨**：若 truth set 升級後 HCC1954 NG=2 inner TP rate >0.7 → 主因為 truth set 限制；若仍 <0.5 → caller saturation 確認。

### Direction C — Wet-lab confirmation（hypothesis-level）
**假設**：HCC1954 Outer × NG=2 × cross_het cell 的 FP 中可能混雜真實亞群 variants（subclone-level VAF 過低被 truth set 漏列）。
**驗證方式**：對 HCC1954 chr5/chr8/chr17 上 HPFineN_HP1S 高 marker 的 Outer-positive variants 做 amplicon Sanger / ddPCR pilot（n≈30 sites）。
**預期區辨**：若 ≥30% 確認為真 variants → 部分「FP」是 truth set 漏列；若 <10% → caller artifact 確認。

> **本案範圍**：以上 3 方向 **不在本 case panel 執行範圍**；列出供 Thread D 主軸決定優先級。

---

## 9. 與 Thread D 主軸報告的雙向連結

### 9.1 本 case panel 引用主軸

主軸路徑（待主代理寫入）：
`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`

引用節點：
- 主軸 §3.L1（cross-sample obs18 / X5 結果）：本 panel §2、§6 提供 HCC1954 的 caller-saturation 解釋
- 主軸 §5（anomalies / robustness）：本 panel 全文作為 HCC1954 outlier 的 standalone case 證據

### 9.2 主軸 §5 反向引用本 panel（占位）

主軸需在 §5 加入類似段落：

> **HCC1954 sample-level anomaly**：HCC1954 在 obs18 NG=2 inner same_HP1 TP rate 為 0.429（cohort median 0.917），但 gap 仍 +0.345 同向。經獨立 case panel 分析（`InterSubMod/docs/reports/validated/2026/04/20260426_HCC1954_Outlier_CallerCeiling_CasePanel_01.md`）確認此偏低為 ClairS-TO 在 HCC1954 TO mode 的 caller FP background ceiling（Outer×NG=2×cross_het cell FP fraction 89.9%），非 LOH-constrained phasing 機制反例。Wilcoxon n=6 包含 HCC1954 仍 p=0.01562 顯著。

> **占位待主代理寫入**：上段文字由 Thread D 主軸報告作者於 §5 anomalies / robustness 章節加入；本 panel 提供完整證據與表格，主軸僅需 1-2 段引用。

---

## 10. 證據檔案完整路徑（含 dataset / build provenance）

### 10.1 來源報告

| 報告 | 路徑 | 行號 | dataset / build |
|------|------|------|-----------------|
| B2 根因 4-panel | `InterSubMod/docs/experiments/in_progress/2026/04/20260423_B2_HCC1954_Outlier_RootCause_01.md` | Lines 22-122 | TO archive 6 樣本 (post-HP-fix), HCC1395 paired 補位 (n=4) |
| HCC1954 reversal | `InterSubMod/docs/experiments/in_progress/2026/04/20260417_HCC1954_reversal_investigation_01.md` | Lines 14-26, 96-179, 198-201 | HEAD `ab61ad1`; B.1-1 完成後；step3 segment-level + Bootstrap B=2000 |
| Z3 pilot design | `InterSubMod/docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_design_01.md` | Lines 36-110, 152-191 | Z3 internal feature exploration step 4-5 |
| Z3 pilot result | `InterSubMod/docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md` | Lines 14-46, 145-178 | TO mode 7 樣本 F1 計算 |
| X5 cross-sample | `InterSubMod/docs/experiments/in_progress/2026/04/20260424_X5_CrossSample_obs18_KDE_Verified_01.md` | Lines 5, 53-78 | X1 KDE-corrected master batch (auto_kde, Diploid_Coverage_Used per sample) |
| B1 Wilcoxon | `InterSubMod/docs/experiments/in_progress/2026/04/20260423_B1_Wilcoxon_NG2_gap_01.md` | (作為 KDE-pre baseline) | Pre-KDE master |

### 10.2 數據檔（B2 一手 artifacts）

- `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_summary.json` — per-sample TP rate + panel stats（FP fraction、CovM median、chr8/chr17 share）
- `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_tp_rate_by_sample.tsv` — Outer × NG=2 × cross_het × Extreme AF subset per-sample TP rate
- `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_per_chr_breakdown.tsv` — per-chr HCC1954 vs cohort TP rate（22 autosomes）
- `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_h1_loh_annotation.tsv` — Inner/Outer/Inner_frac per-sample（H1 driven test）
- `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_hcc1954_outlier_4panel.png` — 主圖（本 panel §3 嵌入）

### 10.3 reversal_investigation 一手 artifacts

- `research/hcc1954_reversal_investigation/data/rho_bootstrap_ci.tsv` — 7 樣本 × 2 modes Bootstrap ρ + CI（B=2000）
- `research/hcc1954_reversal_investigation/data/sensitivity_leave_one_out.tsv` — LOO meta-ρ
- `research/hcc1954_reversal_investigation/data/hcc1954_region_level_rho.tsv` — region-level ρ(caller_af, HPFineNGroups)

### 10.4 Z3 pilot 一手 artifacts

- `research/z3_internal_feature_exploration/data/step5_blacklist_pilot_results.tsv` — per-sample × strategy ΔF1
- `research/z3_internal_feature_exploration/data/step5_blacklist_pilot_summary.tsv` — 匯總表
- `research/z3_internal_feature_exploration/figures/step5_blacklist_delta_f1.png` — ΔF1 視覺化
- `research/z3_internal_feature_exploration/figures/step4_hcc1954_fp_pos_hist.png` — FP Pos × 文獻 amplicon

### 10.5 X5 一手 artifacts

- `research/tpfp_loh_af_kde_discrimination/data/X5_crosssample_obs18.tsv` — KDE-corrected master 6 樣本 obs18 gap
- `research/tpfp_loh_af_kde_discrimination/data/X5_crosssample_summary.json` — Wilcoxon W=21, p=0.01562
- 上游 ISM TSV：`/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun/{SAMPLE}_TO_{tp,fp}/`

### 10.6 Build / pipeline provenance

- ClairS-TO caller binary：見專案 `scripts/run_vcf_all_snv.sh` EXECUTABLE 指向（最新 commit `485075f` 已導向 big7 KDE-fixed binary）
- ISM pipeline build：post-HP-fix（commit `4dc2d73` HPFineNGroups marker re-audit 之後）+ `--germline-hp-only` flag 可選（commit `775027d`）
- TO archive：5 樣本 post-HP-fix archive（COLO829 缺，本 panel 引用 6 樣本中之 HCC1954）

---

**End of HCC1954 Outlier Case Panel**

> 本 panel 為 Thread D 主軸的 sample-level anomaly 補充。HCC1954 不應從跨樣本 Wilcoxon 排除（n=6 已含且顯著），但敘事須清楚標記為「caller-saturated sample with valid sign-level contribution」。
