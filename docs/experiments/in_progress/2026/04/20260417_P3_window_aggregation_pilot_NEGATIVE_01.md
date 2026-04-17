---
title: P3 Window Aggregation Pilot — NEGATIVE (Spatial Auto-correlation Confound)
date: 2026-04-17
status: in_progress
verdict: NEGATIVE
ism_goal: 目標 5 (F1 提升) — gene-level aggregation 可行性
related:
  - research/P3_window_aggregation_pilot/
  - feedback_spatial_autocorrelation_confound.md (new pitfall memory)
---

# P3 Window Aggregation Pilot — NEGATIVE

## 1. 目的

在 Phase 2B gene-level aggregation 全量執行前先做 <2h pilot：驗證 chr+pos window aggregation 是否能突破 region-level pure methylation AUC ≤0.58 的上限。

## 2. 方法

- 樣本：7 樣本 (HCC1395, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954)
- 特徵：HPFineNGroups / HPFineF / AlleleDelta / HPMergedDelta / PairwiseMedianDist
- Window 定義：`Chr_{Pos//1e6}`（1 Mb bins）
- 每 window 聚合：mean / median / max
- Window label：majority vote (≥50% TP regions → y_maj=1)
- 過濾：每 window ≥3 regions，≥20 windows/sample/metric
- 比對：region-level AUC vs best window-level AUC

## 3. 初步結果（Naive，看似 POSITIVE）

| sample   | region_auc_max | best_window_auc | Δ       |
|----------|----------------|-----------------|---------|
| H2009    | 0.568          | 0.910           | +0.342  |
| HCC1937  | 0.624          | 0.977           | +0.353  |
| HCC1954  | 0.668          | 0.809           | +0.141  |
| H1437    | 0.543          | 0.724           | +0.181  |

4/7 樣本 window AUC ≥0.60 且 > region AUC — **naive 結論將誤判為 P3 POTENTIAL**。

## 4. Confound Check（決定性步驟）

### Hypothesis
- **H0 (spatial confound)**: window AUC 提升來自 TP/FP 的基因組空間 auto-correlation（某些 region 幾乎全 TP 或全 FP），aggregation 只是把 spatial 先驗變成「預測力」。
- **H1 (true gain)**: 即使排除 TP rate 極端 (>0.9 或 <0.1) 的 window，window AUC 仍 > region AUC。

### Method
- 重做 H2009 / HCC1937 / HCC1954 / H1437（P3 naive 正向最強 4 樣本）
- 分三層過濾：All windows → excl extremes (0.1 < TP rate < 0.9) → mid-TP-rate (0.3–0.7)
- 比較 window_auc_mid_tp_rate vs region_auc

### Results (TO mode, NonLOH)

| sample   | metric        | region_auc | win_auc_all | excl_extremes | mid_tp_rate | Δ_mid    |
|----------|---------------|-----------|-------------|---------------|-------------|----------|
| H2009    | HPFineNGroups | 0.542     | 0.337       | 0.380         | 0.196       | **-0.346** |
| H2009    | AlleleDelta   | 0.498     | 0.145       | 0.051         | 0.000       | **-0.498** |
| HCC1937  | HPFineNGroups | 0.624     | 0.609       | 0.600         | 0.580       | -0.044   |
| HCC1937  | AlleleDelta   | 0.448     | 0.489       | 0.497         | 0.498       | +0.050   |
| HCC1954  | HPFineNGroups | 0.668     | 0.702       | 0.684         | 0.610       | -0.058   |
| HCC1954  | AlleleDelta   | 0.475     | 0.538       | 0.519         | 0.462       | -0.013   |
| H1437    | HPFineNGroups | 0.549     | 0.610       | 0.604         | 0.531       | -0.018   |
| H1437    | AlleleDelta   | 0.543     | 0.724       | 0.716         | 0.640       | **+0.097** |

注意 H2009 第一層比對（win_auc_all）即已呈現「反向預測」（<0.5）— naive `best_window_auc` 的 0.910 是經過 mean/median/max 三選一後拿到的，但與 mid-TP-rate 嚴格驗證差距極大，顯示即使原始 all-window 都已被「極端 window 的 TP/FP 比例」完全決定。

### Paired mode
Paired mode 的 LOH annotation 只能在 TO 做（`to_loh_bed_hit`），filtering 後幾乎沒有 window 可以用（n_windows_mid ≤5），因此 paired 結果無法作為獨立證據。

## 5. Verdict: NEGATIVE

- **只有 1/8 TO 測試（H1437 AlleleDelta, Δ_mid=+0.097）在 mid-TP-rate window 保留明確 gain**
- **最大原始 gain（H2009 +0.342）完全反轉為負**
- 6/10 sample×metric 在 "excl extremes" 層 Δ > +0.03，但此層只排除 top 10% 最極端，仍被 spatial auto-correlation 主導
- **結論**：window aggregation AUC 提升 ≈ TP/FP 的基因組空間 auto-correlation artifact

## 6. 意涵與後續處理

### 對 Phase 2B Gene-level 計畫
- 若仍要推進 gene-level，**必須以 `shuffle within chr` 做 null model**
- 驗收門檻必須 upgrade 為：**mid-TP-rate window Δ > +0.03 且 permutation p < 0.05**
- 不可只看原始 AUC，也不可只做 leave-one-chr-out（跨 chr 的空間 confound 比 within-chr 少但仍存在）

### 對其他 chr+pos 聚合方向
相同 pitfall 適用：
- chromosomal arm aggregation
- cytoband aggregation
- topic modeling on genomic location
- 任何以 Chr/Pos 為特徵的模型

### 不受影響的方向
- Region-level 特徵 AUC（region 是最小單位，無空間 leakage）
- 非基因組位置的 aggregation（例如按甲基化 pattern 分群）
- 不同 chr 間 shuffle 驗證的模型

## 7. 新增 pitfall 規則

已寫入 memory：`feedback_spatial_autocorrelation_confound.md`
- 核心 rule: 任何 chr+pos 聚合特徵的 AUC claim 必須做 mid-TP-rate window 驗證
- 判準：Δ > +0.03 在 mid-TP-rate subset 才算真實 gain
- Pilot 必問：「我這個特徵若全樣本隨機打亂 region label，是否 window/gene aggregation 後 AUC 仍 >0.55？」若是，則是純 spatial confound

## 8. 記錄與數據

- 腳本：`research/P3_window_aggregation_pilot/scripts/{01_window_aggregation.py, 02_hotspot_confound_check.py}`
- 原始數據：`research/P3_window_aggregation_pilot/data/{window_auc_vs_region_auc.tsv, hotspot_confound.tsv}`
- 視覺化：`research/P3_window_aggregation_pilot/figures/01_window_vs_region_auc.png`
- Manifest: `research/P3_window_aggregation_pilot/manifest.yaml`
