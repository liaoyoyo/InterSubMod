# HPFineNGroups 飽和效應質疑驗證（B.1-2）

**建立**：2026-04-17
**質疑來源**：`docs/plans/opus-4-7-big8-disk-liaoyoyo2001-knowled-cryptic-moore.md` Part B.1-2
**原結論**：HPFineNGroups N≥4+NR≥80 NonLOH TP rate=89.1%（25,514 regions）— 確認為 somatic heterogeneity marker

## 背景

- HPFineNGroups 計算：`src/core/LabelTest.cpp:633` `result.fine_n_groups = unique_groups.size()`
- 上限 = 4（HP1 / HP1-1 / HP2 / HP2-1）
- TSV 欄位：`HPFineNGroups` (col 27)，在 `output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`
- 資料規模：748K rows，7 samples × 2 modes

## 質疑核心

**「N≥4+NR≥80 TP rate 89.1%」是 N≥4 的貢獻還是 NR≥80 的貢獻？**

由於 NGroups 上限=4，高 NR 時大部分 region NGroups 已 ≥4（飽和）。若 NR≥80 的 TP rate 本身就高，89.1% 可能只是 NR≥80 NonLOH 的 baseline 而非 N≥4 的增量。

## 驗證步驟（Step → Verify）

1. **S1**：讀取 all_region_rows.tsv.gz，拆分 TO / Paired mode × NonLOH / LOH。
   → 驗證：總 row 數 ≈ 748K，各 bin n 分佈可列印。

2. **S2**：NR bin × NGroups 雙層交叉 TP rate 表格（TO NonLOH）。
   - NR bins: [0,20), [20,40), [40,60), [60,80), [80,100), [100,150), [150,500)
   - NGroups: 1, 2, 3, 4
   → 驗證：輸出 7×4 TP rate 矩陣 + 樣本數 matrix。

3. **S3**：**核心檢查** — 在 NR≥80 bin 內，NGroups=3 vs 4 的 TP rate 差異。
   - 若 Δ ≥ +10pp → N≥4 有真實增量，結論保留
   - 若 Δ ∈ [5pp, 10pp) → 部分真實，標記為 "weak marker"
   - 若 Δ < +5pp → 飽和效應，89.1% 是 NR≥80 NonLOH baseline，結論降級

4. **S4**：NR 匹配（matched sampling）— 在 NR=[80,100] 中，NGroups=3 vs 4 的 sample 數分佈；計算 Fisher exact p-value 與 Δ TP rate 95% CI。
   → 驗證：p < 0.05 且 CI 不含 0 才算顯著。

5. **S5**：Paired mode 並行檢查（作為獨立複驗）。
   → 驗證：Paired NonLOH 趨勢是否一致。

6. **S6**：每樣本（7/7）NR≥80 內 NGroups=3 vs 4 TP rate 分別計算。
   → 驗證：7/7 是否仍一致方向？effect size 分布？

## 成功標準

- **結論維持（⭐3 不變）**：S3 Δ ≥ +10pp 且 S4 p<0.05；至少 5/7 樣本 S6 方向一致
- **結論降級（⭐3 → ⭐2）**：S3 Δ ∈ [5,10)pp；S6 ≤4/7 方向一致
- **結論翻轉（⭐3 → NEGATIVE）**：S3 Δ < +5pp

## 檔案規範

- scripts/01_saturation_check.py — 主分析
- data/nr_ngroups_crosstab.tsv — S2 輸出
- data/nr_matched_test.tsv — S4 輸出
- data/per_sample_nr80_ngroups.tsv — S6 輸出
- figures/01_nr_ngroups_tp_rate_heatmap.png — S2 視覺化
- figures/02_nr_matched_delta_forest.png — S4/S6 forest plot
- reports/20260417_HPFineNGroups_saturation_check_01.md — 最終報告

## 前後比較 (git)

- 起始 commit：c022918（當前 HEAD）
- 所有程式碼與資料透過 `git add -f research/hpfinengroups_saturation_check/scripts` 加入
- 執行完畢 commit 含結果 TSV，報告獨立 commit
