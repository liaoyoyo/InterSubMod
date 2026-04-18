<!--
建立時間: 2026-03-28
目標: Phase 1A Round 3 — 驗證 LOH features 是否對 read-level ALT/REF 分類有增益
處理範圍: paired-only multi-bio sample637 shard，4-way feature comparison
關聯檔案:
  - scripts/analysis/run_phase1a_round3_loh_feature_test.py
  - docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md
  - docs/experiments/in_progress/2026/03/20260325_Phase1A_incremental_test與multi_bio_validation_round2_01.md
-->

# Phase 1A Round 3：LOH Feature Test

**日期**：2026-03-28
**腳本**：`scripts/analysis/run_phase1a_round3_loh_feature_test.py`
**資料**：`20260325_phase1_manifest_shard_export_paired_multibio_sample637`（86,521 rows, 620 regions, 7 paired samples）

---

## 1. 背景與假設

### 1.1 Round 2 未解決的問題

Round 2 建立了 `methyl+context` 在 paired-pure multi-bio 上的整體增益（external validation delta F1=+0.0112），但 H2009 出現 delta F1=-0.0039，是所有樣本中唯一明確負向。

### 1.2 LOH 研究結論（Round 1-4）的期待

LOH Round 4 發現：H2009 FP 區域中 **68% 有 Potential_LOH=True**（region-level）。假設如下：

> **若將 LOH features (`Potential_LOH`, `LOH_Subtype`) 加入 Phase 1A，應能幫助模型識別 H2009 的 LOH-related FP。**

### 1.3 4-way 比較設計

| 模型 | Numeric features | Categorical features |
|------|---------|------|
| `context_only` | mapq, PairwiseMedianDist, AlleleDelta, CramersV | hp, mode, VerificationClass, PassedGating |
| `methyl+context` | 上述 + 8 個甲基特徵 | 同上 |
| `loh+context` | mapq, PairwiseMedianDist, AlleleDelta, CramersV | 同上 + **Potential_LOH, LOH_Subtype** |
| `methyl+loh+context` | 上述 + 8 個甲基特徵 | 上述 + Potential_LOH, LOH_Subtype |

---

## 2. 結果

### 2.1 整體比較

![Overall F1 Comparison](../../../../../research/loh_investigation/figures/phase1a_round3/fig01_overall_f1.png)

| Split | model | F1 | delta vs context_only | CI | Supported |
|-------|-------|-----|----------------------|----|-----------|
| external_validation | context_only | 0.8611 | 0.0000 | [0,0] | |
| external_validation | methyl+context | **0.8722** | **+0.0112** | [+0.0044,+0.0188] | ✓ |
| external_validation | loh+context | 0.8606 | **-0.0005** | [-0.0022,+0.0007] | |
| external_validation | methyl+loh+context | **0.8723** | **+0.0112** | [+0.0045,+0.0189] | ✓ |

> **關鍵結論**：`loh+context` 整體比 `context_only` 略差（-0.0005）。`methyl+loh+context` 與 `methyl+context` 幾乎相同（差距僅 0.0001）。LOH features 在整體上無增益。

### 2.2 Per-Dataset 拆解

![Per-Dataset Delta F1](../../../../../research/loh_investigation/figures/phase1a_round3/fig02_per_dataset_delta_f1.png)

| 樣本 | LOH-FP% | context F1 | methyl delta | loh delta | methyl+loh delta |
|-----|---------|-----------|-------------|-----------|-----------------|
| COLO829 | 0% | 0.7532 | **+0.0491** ✓ | 0.0000 | **+0.0491** ✓ |
| H1437 | 38% | 0.7622 | +0.0227 | +0.0018 | +0.0227 |
| **H2009** | **68%** | 0.8056 | -0.0039 | **-0.0021** | -0.0038 |
| HCC1395_DORADO | 0% | 0.9597 | **+0.0055** ✓ | 0.0000 | **+0.0056** ✓ |
| HCC1937 | 38% | 0.9406 | +0.0001 | 0.0000 | +0.0001 |
| HCC1954 | 34% | 0.8928 | **+0.0496** ✓ | 0.0000 | **+0.0497** ✓ |

### 2.3 LOH Region 分佈確認

![LOH Region Distribution](../../../../../research/loh_investigation/figures/phase1a_round3/fig03_loh_region_distribution.png)

H2009 的 LOH-FP 最高（68%），但 LOH feature 仍未能改善該樣本的 read-level 分類結果。

### 2.4 H2009 Bucket Analysis

![Bucket Error Rates](../../../../../research/loh_investigation/figures/phase1a_round3/fig04_bucket_error_rates.png)

H2009 的 FP+LOH=Y bucket 在 `loh+context` 下 error rate 與 `context_only` 幾乎相同，未見改善。

---

## 3. 根本原因分析

### 3.1 為何 LOH Feature 在 Read-Level 無效？

`Potential_LOH` 是 **region-level** 特徵：
- 一個 LOH region 中的**所有 reads**（ALT 和 REF）都得到相同的 `Potential_LOH=True`
- Phase 1A 是 **read-level** 分類任務：判斷每個 read 是 ALT-supporting 還是 REF-supporting
- LOH flag 在 region 內對 ALT/REF reads 是常數，**無法區分 ALT vs REF**

因此，LOH feature 對 read-level 分類器（Phase 1A）本質上無用。

### 3.2 H2009 問題的深層原因

H2009 的 FP 主要是 LOH region 中的 germline SNP（68% FP regions 有 LOH）。這些 FP variants 的讀取特徵與真實 somatic variants 相似：
- 甲基化分佈在 LOH region 可能較均勻（單倍型甲基化）
- 添加甲基特徵可能引入噪訊（與 ALT/REF 讀取的相關性低）
- 結果：methyl 反而造成 H2009 regression

### 3.3 架構釐清

```
Phase 1A（read-level ALT/REF）:
  最佳方案 = methyl+context（paired-pure 整體 +0.0112, CI 全正）
  LOH features = 對 Phase 1A 無用（region-level vs read-level 不匹配）

LOH Evidence Panel（region-level FP risk）:
  與 Phase 1A 互補，作為獨立的 region-level 評分層
  Tier A+(≥50) = 2.02× FP enrichment
  不應與 Phase 1A read-level 模型混用
```

---

## 4. 結論

1. **`loh+context` 對 Phase 1A 無整體增益**（-0.0005，CI 跨零）
2. **`methyl+loh+context` ≈ `methyl+context`**（差距 < 0.0001）
3. **H2009 問題未解決**：LOH features 在 read-level 分類中無法幫助，甲基特徵也無效
4. **架構澄清**：LOH evidence panel 是獨立的 region-level 系統，不適合整合進 read-level Phase 1A

### 4.1 Phase 1A 目前最優結論（保守且正確）

| 任務 | 最優模型 | 整體增益 | CI | 說明 |
|------|---------|---------|-----|------|
| paired-pure multi-bio | `methyl+context` | +0.0112 | [+0.0044,+0.0188] | 整體有統計支持 |
| H2009 specifically | 無有效方案 | ~-0.004 | 跨零 | LOH-dominated FP，需別種方法 |

---

## 5. 下一步

1. **接受 LOH 不適合 Phase 1A** 的負面結論
2. **H2009 診斷**：
   - 研究 H2009 FP+Weak bucket 的甲基化模式（LOH region 內的 methylation 分佈）
   - 考慮 LOH-stratified 策略：對 LOH=Y 和 LOH=N regions 分別建模
   - 考慮 H2009 specific normalization（LOH region 甲基化本質上是 haplotype-specific）
3. **將 paired-pure multi-bio 正式升級為 Phase 1A 主要驗證軸**
   - `methyl+context` 是確認的最佳方案（5/6 樣本正向，全局 CI 全正）
   - H2009 是需要獨立解決的異常案例
4. **LOH evidence panel** 作為獨立 region-level 系統，在下游 evidence aggregation 階段使用
