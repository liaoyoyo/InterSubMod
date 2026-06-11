<!--
建立時間: 2026-05-15
更新時間: 2026-05-15
目標: V6 BAM TPFP HP-LOH-CN 計畫的方法學細節（confound guard 協議 + LR 公式 + power 規則）
資料來源:
  - /big7_disk/liaoyoyo2001/InterSubMod/research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3)
  - Plan agent 8 點質疑審計（2026-05-15）
  - Memory feedback_L2_collider_bias.md / feedback_pooled_ols_residualization_trap.md / feedback_spatial_autocorrelation_confound.md
狀態: in_progress
audience: Step 2-4 執行 agent + reviewer
-->

# 02 — Methodology Notes（方法學細節）

> 本檔提供 Step 2-4 所需的所有統計細節（grid 計算、LR 公式、confound guard、power 規則、多重比較校正）。Step 2 執行 agent 必須在開跑前 Read 本檔。

## 1. Grid 結構與 cell 定義

### 1.1 主 grid（3 軸 50-cell）

```
cell_id = (loh, hp_bucket, cov_bin)
```

| 軸 | Bins | 計算邏輯 |
|----|------|---------|
| `loh` | inner / outer | region 是否在 LOH.bed 內（由 PhasingGraph.cpp:1817 產出，VAF ≥0.8 → HOM）|
| `hp_bucket` | same_HP1 / same_HP2 / cross_het / cross_het_inv / other | Thread D 4-bucket（LabelTest.cpp:265-302）+ 多 bucket occupied 為 other |
| `cov_bin` | cov_loss / cov_normal / cov_elevated / cov_gain / cov_high_gain | Coverage_Multiple 切分：<0.7 / 0.7-1.3 / 1.3-1.8 / 1.8-2.5 / ≥2.5 |

**bin 邊界 caveat**：bin 邊界是 SEQC2 truth CN 推導；HCC1395 TO 純度 ~80% → CN=3 的 Coverage_Multiple 應落 ~1.4 而非 1.8。Step 2 開跑前 agent B 須驗 `Diploid_Coverage_Used` 欄位中位數，若 ≠ 61（HCC1395 KDE baseline），改用 per-sample quantile (q33/q67/q90) 而非 hardcoded 邊界。

### 1.2 HP bucket 計算規則

對每 region 取 4 個 family count：HP1, HP1S (= HP1-1), HP2, HP2S (= HP2-1)。

| Bucket | 條件 | 物理意義 |
|--------|------|---------|
| `same_HP1` | HP1>0 AND HP1S>0 AND HP2==0 AND HP2S==0 | germline + somatic 同 hap |
| `same_HP2` | HP2>0 AND HP2S>0 AND HP1==0 AND HP1S==0 | germline + somatic 同 hap（另側）|
| `cross_het` | HP1>0 AND HP2S>0 AND HP1S==0 AND HP2==0 | germline het + 對側 somatic（FP-rich）|
| `cross_het_inv` | HP2>0 AND HP1S>0 AND HP2S==0 AND HP1==0 | 同上反向 |
| `other` | 其他 | 多 bucket 同時 occupied |

> 注意：HPFineNGroups (NG) 是這 4 個 bucket 中**非零 bucket 數**（值域 0-4），不在主 grid 軸中（避免與 hp_bucket 共線）。NG 改作 LR covariate。

## 2. Power stratification

| Stratum | Cell n 範圍 | 處理 |
|---------|------------|------|
| **Powered** | n ≥ 50 | 主分析：TP rate + Wilson 95% CI + LR + confound guard |
| **Marginal** | 30 ≤ n < 50 | 輔助：列 rate 但加 「marginal」flag，不參與 H 判定 |
| **Underpowered** | n < 30 | 只列 cell n、TP/FP count，不做 rate inference |

**Wilson 95% CI 公式**（避免 n 小時 Wald 區間越界）：
```
p_hat = k / n
z = 1.96
denom = 1 + z²/n
center = (p_hat + z²/(2n)) / denom
margin = z * sqrt(p_hat * (1-p_hat) / n + z²/(4n²)) / denom
[CI_lo, CI_hi] = [center - margin, center + margin]
```

**Power gate (Gate 1, step1.5)**：
- 算 marginal n (從每軸獨立 frequency 推每 cell 預期 n)
- **如果 powered cells (n≥50) ≥ 15 個**（50 cells 的 30%）→ 進 Step 2 完整 grid
- **否則** → 降級：
  - 選項 A: 合併 cov_bin 鄰近 bin（5→3）→ 2 × 5 × 3 = 30 cells
  - 選項 B: 降為 2 軸（LOH × HP）→ 10 cells
  - 視 marginal 分布 agent B 決定，記在 power_dry_run.md

## 3. Cell-level LR + Covariate

每 powered cell 跑：

```
logit(P(TP)) = β0 + β1 · NG + β2 · caller_af + β3 · NumReads + ε
```

- **β0** = adjusted intercept = log-odds of TP **after** controlling NG/AF/NumReads
- 拿 `expit(β0)` 作為「cell 內固有 TP rate」（characterization 重點）
- **β1, β2, β3** 的 LRT p-value（vs 移除單一 covariate 的 reduced model）顯示每個 covariate 的獨立貢獻

**Plan agent §5 警告落實**：
- ❌ 禁止 report 任何 cell 「AUC」（filter-readiness metric）
- ❌ 禁止用 AUC 0.58 / 0.7 等 thresholds
- ✅ 改報 **deviance decomposition**: D(full) − D(reduced_per_covariate)
- ✅ 改報 **fraction of deviance explained**: 1 − D(full) / D(intercept-only)

## 4. 5 道 Confound Guard

對所有進入 top-list 的 cells（5 個 top-20 list 合併後去重）做：

### Guard 1: within-group OLS residualize on `NumReads`
```
TP_label ~ NumReads (within cell)
residual_TP_rate = mean(TP - predicted_TP)
```
- 若 residualized rate vs raw rate 差異 > 0.05 → flag「NumReads-driven」

### Guard 2: within-group OLS residualize on `caller_af`
- 同 Guard 1 邏輯，控 AF gradient

### Guard 3: L3 AF-bin 內生驗證
- 在 cell 內按 AF bin (0-0.3 / 0.3-0.5 / 0.5-0.7 / 0.7-1.0) 拆分
- 若任一 AF bin 內 TP rate 仍與 cell average 同向（差 < 0.05）→ 通過

### Guard 4: Permutation test
- 隨機 shuffle TP/FP label 1000 次
- 算每次 shuffled cell rate 的 distribution
- cell observed rate 必須 > 99-th percentile 或 < 1-st percentile

### Guard 5: Marginal expectation
- 算 `expected_rate = product of marginal rates`
- `surprise = observed - expected`
- top-surprise cells（multi-axis interaction）標記

**附加 Guard 6: chr-stratified Mantel-Haenszel**（Plan agent §4）
- 對每 cell 做 chr-stratified contingency
- 若 cell 信號完全消失於 chr-stratification → flag「chr-confound」（避免 chr8 hotspot spatial artifact）

**附加 Guard 7: HP1/HP2 symmetry test**（Plan agent §4）
- same_HP1 vs same_HP2 應對稱（V6 修補後 ratio 1.838 接近 1）
- 若 cell 大量集中 HP1 或 HP2 一側 → 二項式 p<0.01 flag「HP-asymmetry」

## 5. 多重比較校正

**FDR 控制（依 Plan agent §7）**：
- 對所有 50 cells 的 fisher_p / LRT_p 做 **BH-FDR q<0.05** 全局校正
- 對 5 個 top-20 list 合併後的 cells（unique ~60-80 個）做 5 道 confound guard = 300-400 tests → 仍套 BH-FDR

**Permutation 改 max-statistic null**（FWER control）：
- 每次 shuffle 取 max(observed cell rate) 形成 max-statistic null
- cell observed rate 必須 > 95-th percentile of max-statistic distribution → 通過

**H7 判定（修正版）**：
- 「Top 20 extreme cells 中 ≥ 5 個通過所有 confound guard **且 BH-FDR q < 0.05**」

## 6. Spatial autocorrelation guard

依 memory `feedback_spatial_autocorrelation_confound.md`：
- 對 top-list cells 額外做 **mid-TP-rate chr+pos window 驗證**
- Window 大小 = 1 Mb，TP rate 40-60% 為 mid（避免 cell 信號被 chr8 hotspot 純空間 artifact 驅動）
- 若 cell 信號在 mid-window 內仍同向 → 通過；否則 flag「spatial-driven」

## 7. Trajectory 分類（Step 1）

對每 region 算 V3F → V5 → V6 trajectory（同一 region 在三個 BAM 版本下的核心指標如 HPFineNGroups 變化）：

| Class | 條件 | 物理意義 |
|-------|------|---------|
| **A 兩段都改善** | V5 > V3F AND V6 > V5 | Layer 1.5 + priority fix 兩段都有效 |
| **B 只 V5 改善** | V5 > V3F AND V6 ≤ V5 | Layer 1.5 加上去有效，但 priority fix 移除 Layer 1.5 反而退步 |
| **C 只 V6 改善** | V5 ≤ V3F AND V6 > V5 | Layer 1.5 沒幫助但 priority fix 後突破 |
| **D 無改善** | V3F ≈ V5 ≈ V6 (Δ < 1e-3) | 兩段都不影響此 region |
| **E 反向** | V6 < V3F | V6 比原 baseline 還差（紅旗，需 zoom-in）|

「改善」定義 per indicator：
- HPFineNGroups: 數值上升 = 改善
- hp_1_1:hp_2_1 ratio: 接近 1 = 改善（abs(ratio - 1) 變小）
- Marker coverage: 數值上升 = 改善

## 8. Output schema（重要欄位）

每個 cell 至少含：

| 欄位 | 型別 | 來源 | 用途 |
|------|------|------|------|
| cell_id | str | (loh, hp_bucket, cov_bin) | unique key |
| n | int | groupby count | power gate |
| n_TP, n_FP | int | groupby sum | rate |
| TP_rate, FP_rate | float | n_TP/n | 觀察值 |
| TP_wilson_lo, TP_wilson_hi | float | 公式 §2 | 信賴區間 |
| TP_enrichment, FP_enrichment | float | rate / global_rate | 富集度 |
| log_odds | float | log(TP / FP) | discriminator |
| fisher_p | float | scipy.stats.fisher_exact | hypothesis test |
| fisher_q_bh | float | BH-FDR | 多重校正 |
| lr_beta0 | float | LR intercept | adjusted log-odds |
| lr_beta_NG, lr_beta_AF, lr_beta_NR | float | LR coefficients | covariate effects |
| lr_dev_explained | float | 1 - D(full)/D(null) | model fit |
| lr_lrt_p_NG, lr_lrt_p_AF, lr_lrt_p_NR | float | LRT vs reduced | covariate 顯著性 |
| residual_TP_rate_NR | float | Guard 1 | NumReads-residualized |
| residual_TP_rate_AF | float | Guard 2 | AF-residualized |
| L3_consistent_flag | bool | Guard 3 | AF-bin 內生穩定 |
| permutation_extreme_flag | bool | Guard 4 | shuffle test |
| surprise | float | Guard 5 | marginal vs observed |
| mh_chr_stratified_p | float | Guard 6 | chr-confound |
| hp_symmetry_binom_p | float | Guard 7 | HP1/HP2 對稱 |
| spatial_robust_flag | bool | §6 | mid-window 驗證 |
| powered_flag | bool | n ≥ 50 | gate |
| marginal_flag | bool | 30 ≤ n < 50 | aux |
| underpowered_flag | bool | n < 30 | exclude |

## 9. 執行順序（Agent 共識）

1. **Agent A (Step 1)**: 整合 phaseC 12 ISM 結果至 master TSV，產 trajectory
2. **Agent B 預檢 (Step 1.5)**: dry-run 算 marginal n + collinearity (Cramér's V/VIF)
3. **Agent B (Step 2)**: 主 grid + LR + 5+2 道 confound guard
4. **Agent C (Step 3)**: 4 FP zone deep dive + LR deviance decomposition
5. **Agent D (背景並行)**: prior art notes
6. **跨樣本擴展 (Step 4)**: n=5 Wilcoxon（依賴 Step 2 完成）
7. **Coordinator (Step 7)**: synthesis.md + H1-H7 判定

## 10. 參考

- **Plan agent 8 點質疑全文**: 主 session 對話記錄（v0.3 plan 修正歷史）
- **Researcher agent prior art**: `02_prior_art_notes.md`（Agent D 並行產出）
- **Memory**: `feedback_L2_collider_bias.md`, `feedback_pooled_ols_residualization_trap.md`, `feedback_spatial_autocorrelation_confound.md`
- **Thread D 主報告**: `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`
- **V6 設計**: `InterSubMod/research/paired_priority_bug_audit/02_V6_proposal_evaluation.md` + `07_V6_validation_findings.md` + `08_phaseD_v6_cross_sample_findings.md`
