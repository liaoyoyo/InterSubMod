---
title: "B2 · HCC1954 Outer cross-het TP 0.08 outlier 根因專項分析"
date: 2026-04-23
status: in_progress
owner: InterSubMod Research
thread: Thread D · LOH-constrained phasing deepening
track: TO (paired_full HCC1954 暫未具備 post-hoc rescale ISM，本案以 TO archive 為主)
dependencies:
  - docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md (obs18)
  - research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv
  - research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz
inputs:
  - scripts/analysis/20260423_B2_hcc1954_outlier.py
outputs:
  - docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_hcc1954_outlier_4panel.png
  - docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_per_chr_breakdown.tsv
  - docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_tp_rate_by_sample.tsv
  - docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_h1_loh_annotation.tsv
  - docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_summary.json
---

# B2 · HCC1954 Outer cross-het TP 0.08 outlier 根因分析

## 背景

obs18 顯示：在 Outer × NG=2 × cross_het × Extreme AF cells，跨 6 TO 樣本 TP rate 中位數 ≈ 0.67，但 HCC1954 為 **0.08**（本腳本以 cross_het ∪ cross_het_inv 聚合後為 **0.101**），全樣本最低且與次低的 HCC1937 (0.249) 差距大。

三個候選根因：
- **H1**：Potential_LOH 註解在 HCC1954 的可靠性異常（更多 region 被誤歸為 Outer）
- **H2**：AF 分佈（低 AF 比例特別高 → germline 汙染）、CovM 分佈特殊
- **H3**：HER2 (chr17)、MYC (chr8) 已知 HCC1954 amplicon 主導 FP

## 方法

從 6 TO archive 的 `step05_intersubmod/intersubmod_{tp,fp}/significance_summary.csv` 取 NG=2 rows，依 obs18 協議（HPFineN_HP{1,1S,2,2S} bucket occupancy）分類 combo，濾 Outer × cross_het∪cross_het_inv × Extreme AF。比較 HCC1954 與其他 5 樣本聚合的 AF、CovM、chr 分佈。

subset 規模：HCC1954 = 7,400；其他 5 樣本合計 = 24,178。

HCC1395 paired (`/tmp/ism_hp_fix_phase1/...`) Outer cross-het subset 僅 4 rows（該 sample Outer 本就很少），不足以分析；以 HCC1395_DORADO 代表 HCC1395。

## 結果

### 核心觀察 1：HCC1954 FP 是**全基因體**現象，非 amplicon 局部

HCC1954 subset FP fraction = **89.9%**（6,653 / 7,400）。Top 5 FP 染色體為 chr1, chr2, chr3, chr4, chr7，**依 variant count 自然排序，無 amplicon 偏集**。

| 染色體 | HCC1954 n | HCC1954 FP | HCC1954 TP rate | Cohort TP rate |
|--------|-----------|------------|-----------------|----------------|
| chr1   | 824       | 765        | 0.072           | 0.689          |
| chr2   | 699       | 638        | 0.087           | 0.704          |
| chr8 (MYC)   | 224   | 210    | 0.063           | 0.659          |
| chr17 (HER2) | 123   | 117    | 0.049           | 0.634          |
| chr5   | 100       | 98         | 0.020           | 0.749          |

chr8 + chr17 僅佔 HCC1954 全 FP 的 **5.0%**（210+117 = 327/6653），完全無法解釋 89.9% FP 的主因。**H3 REJECTED**。

### 核心觀察 2：AF 與 CovM 分佈無異常

- AF<0.1 比例：HCC1954 100% vs cohort 100%（因 filter 本身限定 Extreme AF，同質化）
- CovM median：HCC1954 **1.08** vs cohort **1.17**（HCC1954 略低，不是更高），CovM 分佈形狀相似

→ **H2 REJECTED**：AF/CovM 分佈無法解釋 TP rate 差異 5-8 倍。

### 核心觀察 3：Potential_LOH 註解分佈正常

NG=2 × cross_het × Extreme AF 中，Inner 比例：

| Sample | Inner | Outer | Inner_frac |
|--------|-------|-------|-----------|
| H1437           |   184 |  8640 | 0.021 |
| H2009           |   303 |  7207 | 0.040 |
| HCC1395_DORADO  |    93 |  6658 | 0.014 |
| HCC1937         |    58 |  1669 | 0.034 |
| **HCC1954**     |   137 |  7400 | **0.018** |

HCC1954 Inner_frac = 0.018，介於 HCC1395_DORADO (0.014) 與 HCC1937 (0.034) 之間，**無任何異常**。→ **H1 REJECTED**：Potential_LOH 註解不是元凶。

### 核心觀察 4：幾乎每一染色體都低 TP rate

逐染色體對照 HCC1954 vs cohort TP rate（n≥50）：HCC1954 在所有 22 autosome 的 TP rate 都 ≤ 10%，cohort 則普遍 60-75%。差距系統性、全基因體均勻：

- chr5：HCC1954 0.020 vs cohort 0.749（Δ = −0.73）
- chr17：HCC1954 0.049 vs cohort 0.634（Δ = −0.59）
- chr8：HCC1954 0.063 vs cohort 0.659（Δ = −0.60）
- chr1：HCC1954 0.072 vs cohort 0.689（Δ = −0.62）

## 結論

**三個原始假設 H1/H2/H3 全部被否決**。HCC1954 Outer cross-het TP 0.08 outlier 並非由註解缺陷、AF/CovM 分佈、或已知 amplicon 驅動，而是**全基因體性 FP 密度普遍過高**（Outer × NG=2 × cross_het 這個 cell 裡 FP 佔比接近 90%）。

### 排序後的真正根因候選

| 排序 | 候選 | 量化證據 | 備註 |
|:----:|------|--------|------|
| 1 | **Truth set 覆蓋率偏低（SEQC2-type 侷限）** | 幾乎所有 chr 同步降低 → 樣本層級系統偏差 | 需比對 HCC1954 truth set 的 variants/Mb 與其他樣本，或該 truth set 是否為 tier1/high-confidence only |
| 2 | **ClairS-TO caller 在 HCC1954 的 precision 天花板低** | FP frac 89.9%，遠高於 HCC1937 (75%)、H1437 (32%)、H2009 (12%) | caller 在這樣本的背景雜訊高；非 ISM/LOH 問題 |
| 3 | **HCC1954 高度 aneuploid 導致全基因體「Outer region 其實也帶 CN 變化」** | CovM median 無明顯偏移但此 subset 才 clipped to 5 | 需檢 unclipped CovM 全分佈、CN=3/4 比例是否異常 |
| 4 | H3 amplicon 局部 | chr8+chr17 僅 5% FP share | DROPPED |
| 5 | H2 AF/CovM | 分佈無差異 | DROPPED |
| 6 | H1 LOH 註解 | Inner_frac 正常 | DROPPED |

### 對 Thread D / LOH-constrained phasing discovery 的意涵

- **結論更正**：先前 obs18 takeaway 描述 HCC1954 為「LOH-constrained phasing 機制在此樣本失效」可能誤導；真相是 **TP denominator 在整個 HCC1954 TO 都被壓低**，cross_het combination 本身的 phasing 邏輯仍成立，只是觀察窗被 caller-level 的低 precision 遮蔽。
- **論文敘事調整**：LOH-constrained phasing 的主要證據應來自其他 5 TO 樣本（TP gap 顯著），HCC1954 歸為「caller-level 高 FP 背景干擾 ISM signal」的獨立 case study，而非機制反例。
- **下游建議**：(a) B5 COLO829 S1 fold 反轉根因分析時，若也出現同樣 caller-level 高 FP，應與此案合併為「caller quality × ISM signal 交互作用」子課題；(b) 跨樣本結論穩定性評分 (`06_結論穩定性審查`) 應將 HCC1954 TO 獨立標記 confidence ↓1。

## 產出檔

- 圖：`docs/experiments/in_progress/2026/04/figures/20260423_B2_hcc1954_outlier/B2_hcc1954_outlier_4panel.png`
- 逐染色體表：`.../B2_per_chr_breakdown.tsv`
- 樣本 TP rate 表：`.../B2_tp_rate_by_sample.tsv`
- H1 LOH 註解診斷表：`.../B2_h1_loh_annotation.tsv`
- 彙整 JSON：`.../B2_summary.json`
- 腳本：`scripts/analysis/20260423_B2_hcc1954_outlier.py`

## 後續（建議，非本案範圍）

1. 檢核 HCC1954 TO truth set 檔案來源、variants/Mb、tier 定義（對照 SEQC2 HCC1395 的 tier 結構）
2. 拉 HCC1954 TO ClairS-TO paired_full 比對：若 paired 模式 TP rate 也低，則 caller 本身是瓶頸；若 paired 正常，則 TO-only pipeline 的額外校正步驟（germline filter）過度激進
3. 跨樣本 FP 密度歸一化後重跑 obs18：以 baseline FP rate 標準化 TP gap，還原 LOH-constrained phasing 的真實效應量
