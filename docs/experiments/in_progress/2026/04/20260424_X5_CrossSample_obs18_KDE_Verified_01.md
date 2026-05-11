---
title: "X5 | Cross-sample obs18 replication on KDE-corrected master (Thread D verification)"
date: 2026-04-24
status: in_progress
verdict: POSITIVE (obs18 6/6 + Wilcoxon p=0.0156 replicated on KDE-corrected master)
phase: Thread D cross-sample robustness check
track: TO
samples:
  - HCC1395
  - HCC1395_DORADO
  - H1437
  - H2009
  - HCC1937
  - HCC1954
hypothesis_id: LOH-constrained-phasing-cross-sample
source_data: /big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun/
artifacts:
  - scripts/analysis/20260424_X5_crosssample_obs18_S3S5.py
  - scripts/analysis/20260424_X5v2_corrected_S3S5.py
  - research/tpfp_loh_af_kde_discrimination/data/X5_crosssample_obs18.tsv
  - research/tpfp_loh_af_kde_discrimination/data/X5_crosssample_summary.json
  - research/tpfp_loh_af_kde_discrimination/data/X5v2_corrected_S3S5.json
related:
  - docs/experiments/in_progress/2026/04/20260423_B1_Wilcoxon_NG2_gap_01.md
  - docs/experiments/in_progress/2026/04/20260424_X3_FlagOnOff_NG2_NegControl_01.md
  - docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md
  - docs/experiments/in_progress/2026/04/20260423_KDE_Smoke_Test_Cross_Sample_Validation_01.md
tags:
  - cross-sample
  - LOH-constrained-phasing
  - KDE-corrected
  - wilcoxon-replication
---

# X5 — KDE-corrected master 下 Thread D obs18 跨樣本複現

> **主結論**：X1 重跑的 6 TO 樣本 KDE-corrected master 上，obs18 的 Inner same_HP1 vs Outer cross_het TP gap **6/6 正向 + Wilcoxon p=0.0156 完全複製 B1**；但 Thread B S3/S5 filter 定義基於 **caller AF（VCF）而非 AlleleDelta（ISM TSV）**，本次以 AlleleDelta 代理驗證失敗，需 VCF AF merge 重算。Thread D 機制結論通過跨樣本 + KDE-corrected 兩重 robustness check。

---

## 1. 目的

驗證兩件事在 KDE-corrected master（X1 batch）上是否仍成立：
1. **Thread D obs18 跨樣本 gap**（原 B1 驗證過 6/6 正向, p=0.0156, median 0.365）
2. **Thread B S3/S5 filter**（原 v2 report HCC1395 TO S3 n=380 TP=95.5%）

## 2. 資料

- X1 batch: 6 samples × TP+FP = 12 runs 全成功
- 全 auto_kde；Diploid_Coverage_Used per sample：53 (HCC1395/DORADO) / 61 (HCC1954) / 71 (H1437) / 79-81 (H2009) / 89-91 (HCC1937)
- HCC1395 TO=53 與 20260420 pilot 完全一致 ✓

## 3. Thread D obs18 跨樣本複現結果

| Sample | Inner same_HP1 TP | n | Outer cross_het TP | n | Gap |
|--------|:-----:|:---:|:-----:|:---:|:-----:|
| HCC1395 | 0.840 | — | 0.580 | — | +0.260 |
| HCC1395_DORADO | 0.939 | — | 0.553 | — | +0.385 |
| H1437 | 0.920 | — | 0.688 | — | +0.231 |
| H2009 | 0.932 | — | 0.882 | — | +0.049 |
| HCC1937 | 0.759 | — | 0.236 | — | **+0.522** |
| HCC1954 | 0.429 | — | 0.084 | — | +0.345 |

**Wilcoxon signed-rank (alternative=greater)**：W=21.0, p=0.01562, n=6

### 與 B1 (pre-KDE master) 對照

| 統計 | B1 (pre-KDE) | X5 (post-KDE) | 差異 |
|------|:---:|:---:|:---:|
| 6/6 正向 | ✓ | ✓ | 一致 |
| Wilcoxon W | 21 | 21 | 一致 |
| Wilcoxon p | 0.01562 | 0.01562 | 一致 |
| median gap | 0.365 | 0.302 | -0.063（輕微下降但同方向）|

**結論**：Thread D 機制通過 **KDE-corrected master 的 robustness check**。gap median 輕微下降但 Wilcoxon p 完全一致，說明 KDE 校正不影響機制 ordering。

個別 sample 差異主要來自：
- HCC1395 Inner same_HP1 從 B1 的 0.959 降為 0.840 — 新 master 全基因體 data 包含更多 regions（特別是 chr19 以外的 low-coverage regions）
- 其他樣本變動 <±0.05

## 4. Thread B S3/S5 filter 驗證（FAILED — 需 caller AF merge）

### 4.1 發現

使用 X1 master 的 `AlleleDelta` 作為 AF_class 分類（`[0.4, 0.6]` 為 Near-half）：

| Sample | total | Near-half AF 佔比 | S3 n | S3 TP rate |
|--------|:-----:|:-----:|:---:|:---:|
| HCC1395 | 40,096 | **0.2%** | 2 | — |
| HCC1395_DORADO | 40,428 | — | 0 | — |
| H1437 | 58,915 | — | 30 | 0.80 |
| H2009 | 137,695 | 0.2% | 87 | 0.93 |
| HCC1937 | 24,655 | — | 0 | — |
| HCC1954 | 67,286 | 0.2% | 79 | 0.24 |

**跨樣本 S3 達 ≥85% n≥20：1/6**（僅 H2009）

### 4.2 根因

- ISM TSV 的 `AlleleDelta` **不是 caller AF**，而是 ISM 計算的 allele-level methylation delta
- 典型 caller AF 分布 Near-half ≈ 20-40%；此處 AlleleDelta 只 0.2%
- Thread B v2 report（`20260422_LOH_AF_KDE_TPFP_Discrimination_02.md`）的 S3 n=380 是用 VCF caller AF 計算
- 正確驗證需 **merge VCF AF 到 ISM TSV**（以 `RegionID` 或 `Chr:Pos:Ref:Alt` 為 key）

### 4.3 Memory feedback 應用

符合 `feedback_feature_name_vs_definition_rule.md` — `AlleleDelta` 字面上像 AF，但實際計算並非，需查源碼確認。已避開更大誤解。

## 5. 結論

### 5.1 成功驗證
- ✅ **Thread D obs18 跨樣本複現**：6/6 正向、Wilcoxon p=0.01562
- ✅ **KDE-corrected master** 不影響機制結論
- ✅ 配合 B1（obs18 原始）+ B3（paired NC）+ X3（flag=on NC）構成 **4 重驗證**

### 5.2 後續需補
- ❌ Thread B S3/S5 filter 跨樣本：需 VCF caller AF merge 後重算
- 建議 task：寫 `merge_vcf_af_to_ism.py`，將 `filtered_snv_{tp,fp}.vcf.gz` 的 AF INFO 欄位 merge 到 ISM TSV

### 5.3 對論文影響
- Thread D（主軸）穩健性 ✓✓ 全 gated
- Thread B（附軸，落地工具）仍需補 — 優先級可降為 P2

## 6. Artifacts

- 6 TO 樣本 full-genome ISM: `/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun/{SAMPLE}_TO_{tp,fp}/`
- obs18 cross-sample table: `research/tpfp_loh_af_kde_discrimination/data/X5_crosssample_obs18.tsv`
- X5 主結果 JSON: `research/tpfp_loh_af_kde_discrimination/data/X5_crosssample_summary.json`
- X5 v2 S3/S5（caller-AF FAIL）: `research/tpfp_loh_af_kde_discrimination/data/X5v2_corrected_S3S5.json`
- Script: `scripts/analysis/20260424_X5_crosssample_obs18_S3S5.py`
- Script (v2): `scripts/analysis/20260424_X5v2_corrected_S3S5.py`
