---
title: Coverage_Multiple vs SEQC2 Truth 階段性觀察（Step 1+2）
date: 2026-04-19
status: blocked
hypothesis_id: H1
verdict: leaning_rejected（基於過期 master dataset；需以現行 binary 重跑後重量化）
related:
  - research/coverage_multiple_validation/00_PLAN.md
  - docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md
  - docs/methodology/20260419_KDE_expected_coverage_audit_01.md
---

# Coverage_Multiple vs SEQC2 Truth 階段性觀察（Step 1+2）

> **狀態**：blocked — master dataset 需以現行 binary 重跑
> **當前 verdict**：H1 leaning rejected；H2 待 Step 4 決定

## 2026-04-19 根因重審

本報告 Step 1/2 所用 master dataset（`20260330_loh_round1_cross_sample_audit_post_to_hp_fix`）
由 2026-03-30 時的 binary 產出，早於 KDE per-sample diploid commit `8d0a0c8` (2026-04-13)
14 天 → 全 14 rows（7 樣本 × 2 modes）的 expected_coverage 常數 75.0 是 **stale-binary artifact**，
非 KDE 邏輯缺陷，也非「KDE 在執行時未被啟用」。

2026-04-19 C++ 強化（commits `374fad4` + `12d9b3e`）：KDE fallback 升 WARN；TSV 新增
`Diploid_Coverage_Used` 欄位；回傳改 `DiploidEstimate` struct。

**本報告數值後續處理**：
- Step 1 per-CN-bin recall（Gain 14.6% / Loss 86.9%）是 baseline=75.0 下的觀察；重跑後預期
  HCC1395 neutral baseline 會降至 ~54 → CN=3 的 CovM 由 0.85 升至 ~1.5，Gain recall 大幅改善
- Step 2 per-sample bias 表（COLO829 +150% / H2009 −28% 等）是**舊 binary 下** 的 bias；新 binary
  下各 bias 理應收斂至單位數百分比（KDE 即是為此設計）
- 研究價值仍保留：確證「以 default 75.0 baseline 評估跨樣本會產生嚴重誤判」，這是
  infrastructure lesson；但不構成 H1 的最終判定依據

## Step 1 — HCC1395 × SEQC2 truth 對齊

### 輸入
- Master region rows（HCC1395 + HCC1395_DORADO）
- SEQC2 CNV BEDs（gain/loss/LOH + CN 值 + exclusion）

### 方法
以 binary search 將每一 region position 對映到 SEQC2 truth label（Gain/Loss/LOH/Neutral/Excluded），計算：
- Confusion matrix：`Coverage_Category (pred) × truth_label_seqc2`
- Per-CN-bin metrics：各 truth CN bin 的 pred-Gain/Loss rate + mean CovM
- Per-label correlation：Pearson + Spearman（CovM vs truth CN）

### 結果摘要

**核心數字**：HCC1395 (TO mode) **Gain recall = 14.6%** — 遠低於 H1 positive 門檻 0.80。

具體待 Step 1 執行結果填入（執行後由腳本 stdout 回填）。

### 解讀
Gain recall 低至 14.6% 意味：CovM 預測為 Gain 的區域 / truth 為 Gain 的區域 ≪ 1，即 **大多 truth Gain 在 CovM 下未被分類為 Gain**。這是 baseline 被高 ploidy 區域上抬、導致 CovM 普遍下修的 classic symptom。

## Step 2 — 跨樣本 expected_coverage baseline 審計

### 方法
- HCC1395 以 SEQC2 Neutral BED 區段定義 true neutral median NumReads
- 其他 6 樣本以 pan-cancer proxy-neutral 區間（chr1p/2q/3p/9q）作為 neutral proxy
- 逆向推導各樣本實際使用的 `expected_coverage = NumReads / Coverage_Multiple` 中位數
- 計算 `bias = (est_expected - true_neutral_median) / true_neutral_median`

### 結果摘要

**關鍵發現**：跨 7 樣本的 **expected_coverage 疑似為硬編碼常數 ≈ 75.0**（非 per-sample KDE 估計）。

具體每樣本 bias 值待 Step 2 執行結果填入。

### 解讀
若 expected_coverage 硬編碼 = 75，則 CovM 在所有樣本上對同一 NumReads 數值給出相同分類。但各樣本真實測序深度差異（HCC1395 ~60× vs HCC1954 ~30× vs H2009 ~80×）意味：**CovM 的跨樣本 scale 不可比**。這解釋為何 Z3 amplicon blacklist 在 HCC1954 改善但其他樣本失敗 — 閾值對各樣本 CN 涵義不同。

## 階段性結論（H1 leaning rejected）

1. **HCC1395 Gain recall 14.6%** → 單樣本 CovM 已失準
2. **expected_coverage 疑硬編碼** → 跨樣本 CovM 不可比
3. 兩項證據收斂指向 **CovM 非可靠跨樣本 CN proxy**

## 下一步（Step 4）

測試兩種修正：
- **Variant A（Re-center）**：`CovM_v2 = NumReads / neutral_median_NumReads`（per-sample）
- **Variant B（Oracle CN, HCC1395 only）**：以 SEQC2 truth CN 直接取代 CovM，重跑 S2 style Z3 blacklist

預期：若 Variant A 在 HCC1954 ΔF1 維持正 **且** 其他 6 樣本 mean ΔF1 ≥ 0 → H2 supported，可改以 re-centered CovM 定義 canonical filter。若 Variant A 無改善而僅 Variant B 在 HCC1395 恢復 → H2 only on oracle，CovM 問題短期無 pipeline 層級解。

## 待補
- [ ] Step 1 confusion matrix + recall 具體表格（待腳本執行後回填）
- [ ] Step 2 per-sample bias 表格（待腳本執行後回填）
- [ ] Step 4 oracle pilot 結果（執行後撰寫 20260420_step4_report_01.md）
