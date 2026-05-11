---
title: "X6 | Caller AF merge + Thread B S3/S5 跨樣本驗證 (post-KDE)"
date: 2026-04-24
status: in_progress
verdict: MIXED — AF merge 技術成功；S3/S5 跨樣本不穩定，原 HCC1395 TO S3 TP=95.5% 無法複現
phase: Thread B cross-sample robustness on KDE-corrected master
track: TO
samples: [HCC1395, HCC1395_DORADO, H1437, H2009, HCC1937, HCC1954]
hypothesis_id: Thread-B-S3-cross-sample
source_data:
  ism_master: /big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun/
  vcf_source: archive/202603_early_pilots/*/step02_benchmark_clairs_to/filtered_snv_*.vcf.gz
artifacts:
  - scripts/analysis/20260424_X6_merge_caller_af_S3S5.py
  - research/tpfp_loh_af_kde_discrimination/data/X6_caller_AF_S3S5.json
related:
  - docs/experiments/in_progress/2026/04/20260424_X5_CrossSample_obs18_KDE_Verified_01.md
  - docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md
tags:
  - caller-AF
  - Thread-B
  - S3-Diploid-Het
  - cross-sample-fail
---

# X6 — Caller AF merge + S3/S5 跨樣本驗證

> **結論**：技術上 caller AF merge 成功（0 missing / 369K rows）；但 **Thread B S3 Diploid Het filter 在 KDE-corrected 6 TO 樣本上跨樣本不穩定**（S3 TP≥0.85 僅 1/6；Wilcoxon S3 > baseline W=0, p=1）。原 Thread B v2 report 的 HCC1395 TO S3 n=380 TP=95.5% **無法在 KDE-corrected master 複現**（本次 n=2,200 TP=0.58）。Thread B 附軸需重審或降級為 per-sample characterization。

---

## 1. 方法

1. 對 6 TO 樣本，解析 `filtered_snv_{tp,fp}.vcf.gz` 的 `FORMAT:AF` 欄位
2. 以 `(Chr, Pos, Ref, Alt)` merge 到 X1 KDE-corrected master（12 runs）
3. 用 Thread B v2 report 的 scheme 定義重算：
   - **S1** LOH_Strong + Extreme AF
   - **S2** LOH_Weak ∪ LOH_Subclone
   - **S3** LOH_None + Near-half AF + CN_tier ∈ {T1, T2}
   - **S4** LOH_None + Extreme AF (ambiguous)
   - **S5** (S1 ∨ S2 ∨ S3) \\ S4
   - **S6** S1 + NG≥3
   - **S7** S5 + NG≥3
4. CN tiers: T0 <0.65, T1 [0.65,0.99), T2 [0.99,1.33), T3 [1.33,1.82), T4 ≥1.82

## 2. 結果

### 2.1 Caller AF merge 技術驗收

| Sample | rows | AF missing | AF median | AF distribution |
|--------|-----:|-----------:|----------:|-----------------|
| HCC1395 | 40,096 | **0** | 0.479 | Int 20,012 / NH 10,709 / Ext 9,375 |
| HCC1395_DORADO | 40,428 | **0** | 0.485 | Int 19,597 / NH 10,981 / Ext 9,850 |
| H1437 | 58,915 | **0** | 0.532 | Int 32,728 / NH 13,860 / Ext 12,327 |
| H2009 | 137,695 | **0** | 0.506 | Int 57,924 / NH 53,374 / Ext 26,397 |
| HCC1937 | 24,655 | **0** | 0.503 | Int 11,411 / NH 5,635 / Ext 7,609 |
| HCC1954 | 67,286 | **0** | 0.439 | Int 46,176 / NH 17,589 / Ext 3,521 |

**Merge 100% 成功**；AF 分布符合 caller 輸出預期（Near-half 15-39%）。修正 X5 v2 用 AlleleDelta 導致 Near-half 僅 0.2% 的誤差。

### 2.2 S3 Diploid Het 跨樣本 TP rate

| Sample | S3 n | S3 TP rate | S3 FP reduction | vs baseline |
|--------|-----:|-----------:|---------------:|-------------|
| HCC1395 | 2,200 | 0.583 | 0.921 | ↓ from baseline ~0.71 |
| HCC1395_DORADO | 1,805 | 0.597 | 0.937 | ↓ |
| H1437 | 6,588 | 0.738 | 0.872 | flat |
| **H2009** | **31,436** | **0.903** | 0.746 | ≈baseline 0.93（飽和）|
| HCC1937 | 1,920 | 0.358 | 0.898 | **↓↓** |
| HCC1954 | 10,729 | 0.129 | 0.814 | **↓↓↓** |

**S3 TP≥0.85 & n≥20**：**1/6**（僅 H2009，且 H2009 baseline 本就 0.93 → 飽和）
**Wilcoxon S3 > baseline (one-sided)**：W=0.0, **p=1**（S3 系統性**低於** baseline）

### 2.3 S5 Combo 跨樣本

| Sample | S5 n | S5 TP rate | FP reduction |
|--------|-----:|-----------:|-------------:|
| HCC1395 | 10,652 | 0.726 | 0.749 |
| HCC1395_DORADO | 8,330 | 0.760 | 0.827 |
| H1437 | 11,516 | 0.779 | 0.811 |
| H2009 | 43,726 | 0.915 | 0.691 |
| HCC1937 | 5,414 | 0.636 | 0.836 |
| HCC1954 | 15,208 | 0.182 | 0.752 |

**S5 TP≥0.85 & n≥50**：1/6（同 S3，僅 H2009 飽和）

### 2.4 與原 Thread B v2 report 對比

| 指標 | 原 v2 (HCC1395 TO) | X6 post-KDE |
|------|-------------------:|------------:|
| S3 n | 380 | **2,200（5.8×）**|
| S3 TP rate | **95.5%** | **58.3%（-37pp）** |
| S3 FP reduction | 99.85% | 92.1% |
| S5 n | 886 | 10,652 |
| S5 TP rate | 91.8% | 72.6% |
| S5 FP reduction | 99.37% | 74.9% |

**巨大差異**。可能原因：
- Thread B v2 基於 stale-binary master（75.0 default）→ Coverage_Multiple 系統偏低 → CN_tier 分類不同
- 或 v2 使用的 CN tier cutoff 與 memory 記錄的 0.65/0.99/1.33/1.82 不同
- 或 v2 的 LOH_Subtype 或 AF_class 定義有細節差

**無論如何，Thread B 的 "HCC1395 TO S3 = 95.5% TP rate whitelist" 結論在跨樣本 + KDE-corrected 條件下不成立**。

## 3. 特殊觀察

### 3.1 HCC1954 系統性 S3/S5 崩潰
- S1 TP=0.34, S2=0.31, S3=0.13, S5=0.18
- 對照 Thread D obs18 HCC1954 Inner TP=0.43, Outer=0.08
- 再次確認 `project_loh_subclone_af_methylation_positive.md` 的 "HCC1954 subclonal 分析不可用" 警告

### 3.2 HCC1954 S4 反常高 TP（0.81 n=262）
- 原 Thread B 定義 S4 = "None + Extreme AF = ambiguous / 風險"
- HCC1954 S4 反而 TP 0.81 → 高 purity + CNV 環境 Extreme AF 主要是 CN-gain somatic alleles
- 表示 **S4 定義本身在 CNV-heavy 樣本反向**

### 3.3 H2009 唯一飽和
- 所有 scheme TP rate ≥0.90 → orthogonal-tools truth + FP 少 導致
- H2009 不是 S3 成功證據，是 baseline 飽和 artifact

## 4. 對研究優先序的影響

| 原認定 | X6 證據 | 調整 |
|-------|---------|------|
| Thread B S3 filter 可推廣至 TO 全樣本（論文 whitelist）| 6/6 跨樣本不穩，5/6 樣本 TP rate < baseline | **降級為 per-sample characterization**（HCC1395 TO case study 需重寫）|
| Thread B v2 HCC1395 TO S3 n=380 TP=95.5% | 無法複現（n=2,200 TP=0.58） | **標記為 stale-binary artifact**，需 v3 report 或撤回 |
| Thread D 論文主軸（LOH-constrained phasing）| obs18 6/6 在 X5 完全複現 | **無影響**，主軸仍穩固 |

## 5. 結論

- ✅ Caller AF merge pipeline 技術成功，可複用於其他 TSV × VCF merge 需求
- ❌ Thread B S3 "Diploid Het whitelist" 跨樣本不穩定，不符合論文落地工具標準
- ⚠ Thread B v2 report 的 HCC1395 TO S3 結論需加註 stale-binary + CN_tier-sensitive caveat
- ✓ **Thread D 論文主軸不受影響**（X5 obs18 6/6 + Wilcoxon p=0.0156 完全穩定）

## 6. 後續建議

- 建議撤銷或重寫 Thread B v2 report 的 "S3 Diploid Het whitelist" 宣稱
- Thread B 在論文中降為 "HCC1395 TO case study"，不跨樣本推廣
- 資源集中在 Thread D 論文撰寫與 Wakhan/SAVANA 外部 CN pilot（Thread B 替代方案）

## 7. Artifacts

- Script: `scripts/analysis/20260424_X6_merge_caller_af_S3S5.py`
- JSON: `research/tpfp_loh_af_kde_discrimination/data/X6_caller_AF_S3S5.json`
