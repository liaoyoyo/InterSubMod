---
id: exp_20260423_B4_S4_secondary_discrimination
date: 2026-04-23
sample: HCC1395
mode: to_pileup
bucket: S4 (LOH_Subtype=None, AF extreme <0.1)
status: POSITIVE (with caveat)
track: Thread B × Thread D cross
hypothesis: H_B4
---

# B4 | HCC1395 TO S4 Bucket Secondary Discrimination Pilot

## 1. 背景與目的

Thread B 以五維度（LOH_Subtype × AF × Coverage_Multiple × HPFineNGroups × Quality_Score）切分 TO 7 樣本，發現 **S4 bucket = LOH=None + Extreme AF** 在 HCC1395 TO 的 TP rate 與整體 baseline 幾乎相同（71.15% vs 71.07%），**一級切分完全無辨別力**。

B4 假說：Thread D 的 "Inner same-hap fingerprint"（NG=2 時 bucket 組合 `{HP1, HP1-1}` 或 `{HP2, HP2-1}`，單一單倍體來源）可能挑出 S4 內部的 high-precision subset，形成**二級判別**。

## 2. 資料與方法

- **來源**：
  - `/tmp/ism_hp_fix_phase1/{tp_off,fp_off}/significance_summary.csv`（HCC1395 TO post-HP-fix ISM 特徵表）
  - `research/tpfp_loh_af_kde_discrimination/data/master_extended.tsv.gz`（AF 欄位 merge）
- **S4 filter**：`LOH_Subtype=="None"` & `(AF<0.1 | AF>0.9)`
  - n = 30,432（TP=21,652, FP=8,780, TP rate=0.7115）
  - 注意：HCC1395 TO S4 全部落在 `AF<0.1`（0 個 AF>0.9）
- **同單倍體指紋 (SameHapFingerprint)**：
  - = 1 if `(HP1>0 & HP1S>0 & HP2==0 & HP2S==0)` OR `(HP2>0 & HP2S>0 & HP1==0 & HP1S==0)`
- **候選特徵（17）**：HPFineNGroups、SameHapFingerprint、HP_family_imbalance、HP1_family_frac、HP_Ratio、Quality_Score、AlleleDelta、HPFine_NGroups_CF、AF、NTumorReads、NNormalReads、NumCpGs、Fisher_Frac_Sig、Entropy_Imbalance、HPMergedDelta、HPFineF、CramersV_HPFine
- **模型**：LogisticRegression (`class_weight='balanced'`)，StandardScaler，5-fold Stratified CV

### 陷阱偵測：run-level baseline leakage

首輪分析 LR AUC = 0.999，經 drop-one ablation 追蹤到 **Coverage_Multiple + NumReads 組合獨立達 AUC 1.000**。根因：`Diploid_Coverage_Used` 在 `tp_off` run 為 115、`fp_off` run 為 113（run-level 常數），`Coverage_Multiple = NumReads / Diploid_Coverage_Used` 因此可被 LR 反解出 label。**`Coverage_Multiple` 與 `Diploid_Coverage_Used` 已從特徵集中排除。**

## 3. 結果

### 3.1 Single-feature AUC (top 10)

| Rank | Feature | AUC | Direction |
|-----:|:--------|----:|:---------:|
| 1 | AlleleDelta | 0.668 | − |
| 2 | AF | 0.661 | − |
| 3 | HPMergedDelta | 0.578 | + |
| 4 | HP_Ratio | 0.554 | + |
| 5 | HP1_family_frac | 0.554 | + |
| 6 | Fisher_Frac_Sig | 0.549 | + |
| 7 | HP_family_imbalance | 0.548 | + |
| 8 | HPFineNGroups | 0.541 | + |
| 9 | HPFine_NGroups_CF | 0.541 | + |
| 10 | HPFineF | 0.532 | + |

**SameHapFingerprint 本身無法辨別**：S4 內部僅 8/30,432 regions 觸發（AF<0.1 使 HP1 族幾乎無 reads），無 cross-sample signal。

### 3.2 LogisticRegression 5-fold CV

| Config | Mean AUC | Std | Verdict |
|:-------|---------:|----:|:-------:|
| 全 17 特徵 | **0.717** | 0.002 | POSITIVE (≥0.65) |
| AlleleDelta 單一 | 0.668 | 0.005 | — |
| AlleleDelta + AF | 0.668 | 0.005 | — |
| NO AlleleDelta/AF（純 methylation+HP） | 0.619 | 0.006 | — |
| NO AlleleDelta（保留 AF） | 0.693 | 0.003 | — |

### 3.3 Top-K precision (LR 機率排序)

| Percentile | n | TP precision | 提升 vs baseline (0.711) |
|:-----------|--:|-------------:|------------------------:|
| Top-10% | 3,035 | **0.980** | +0.269 |
| Top-20% | 6,071 | 0.962 | +0.251 |
| Top-30% | 9,107 | 0.933 | +0.222 |
| Top-50% | 15,179 | 0.846 | +0.135 |

**Top-10% subset TP precision 98%** — 清楚的 high-precision enrichment。

## 4. 判定與警訊

### 判定：POSITIVE（有條件）

- LR 5-fold CV AUC **0.717 ≥ 0.65 閾值** → 二級判別路徑存在
- Top-10% precision 0.98 具備實務價值

### 警訊與限制

1. **Thread D same-hap fingerprint 在 S4 無效**：S4 全部 AF<0.1（低 VAF，多為 FP-enriched caller 候選），SameHapFingerprint 只觸發 0.03% 的 regions，原假說不成立。
2. **主要訊號來自 AlleleDelta + AF**：這兩者在 S4 內本來就有 AUC≈0.66；LR 僅額外提供 +0.05 AUC。純甲基化 + HP feature subset 的 AUC 僅 0.619，接近 ISM beyond-AUC exhaustion 的歷史結論。
3. **Run-level leakage 發現**：這次診斷找到 `Coverage_Multiple + NumReads` 會 exact-label-leak（原因為 TP/FP 兩個 ISM run 使用不同 Diploid_Coverage_Used 常數）。此陷阱需在未來所有 ISM-feature AUC 分析中注意。

### 與既有結論對照

- 與 MEMORY 中 `project_beyond_auc_exhaustion_confirmed.md` 一致：pure methylation AUC ≤ 0.58（本分析 0.619 在 S4 內稍高，因為 S4 是 extreme AF 子集，和 full dataset 相比 HP-ratio 類特徵略有重新分層）。
- 與 `project_clairs_to_verdict_pilot.md` 的 TO S1 ΔF1=0 一致：即使能挑 high-precision subset，實際 F1 提升需驗證後續 filter scheme 下游量化。

## 5. Artifacts

- **Script**：`scripts/analysis/20260423_B4_S4_discrimination.py`
- **Data**：
  - `research/tpfp_loh_af_kde_discrimination/data/B4_S4_auc_summary.tsv`
  - `research/tpfp_loh_af_kde_discrimination/data/B4_S4_lr_cv.json`
- **Figures** (`docs/experiments/in_progress/2026/04/figures/20260423_B4_S4/`)：
  - `feature_auc_bar.png` — Top-10 single-feature AUC
  - `lr_roc_curve.png` — LR 5-fold CV pooled ROC (AUC 0.717)
  - `samehap_fp_tp_stack.png` — SameHapFingerprint 分佈 (nearly all regions fall in "No")

## 6. 下一步建議

1. **跨樣本擴展**：把同樣 S4 secondary discrimination pipeline 跑 6 個 TO 樣本（COLO829/H1437/H2009/HCC1937/HCC1954/HCC1395_DORADO），確認 LR AUC ≥0.65 是否 generalise。
2. **量化 F1 下游影響**：Top-10% 98% precision → 若以此做 filter cutoff，實際在 ClairS-TO 全流程下 ΔF1 多少？回到 Phase 1/F1 pipeline 跑一次。
3. **殘差化驗證**：把 AlleleDelta / AF residualize 掉後，LR AUC = 0.619 是否 permutation-test significant？走一次 `/auc-confound-guard`。
4. **更新 B3 結論**：Thread D same-hap fingerprint 機制在 AF<0.1 區**不成立**，需更新 LOH-constrained phasing discovery 文件的適用範圍（限 Inner bucket，非 Extreme AF bucket）。

## 7. 證據鏈

- Hypothesis ID: H_B4 (Thread B × Thread D cross)
- Validation Level: L2 (single-sample CV)
- Tier: 2 (single sample, not yet replicated)
- Link to: `project_loh_constrained_phasing_discovery.md` (Thread D 機制在 S4 外成立，S4 內失效)
