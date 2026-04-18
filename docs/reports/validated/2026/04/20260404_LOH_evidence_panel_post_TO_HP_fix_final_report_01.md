<!--
建立時間: 2026-04-04
目標: 以 TO HP integer tag fix 後的重跑輸出，重寫 LOH evidence panel 最終結論，並保留可追溯的數據與來源紀錄
處理範圍: LOH Round 1-4 post-fix workspace、TO HP fix 影響邊界、paired/TO feature 定位、binary filter 可行性
複寫對象:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260330_TO_HP_fix重跑驗證與validated文件修正分級_01.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix
-->

# LOH Evidence Panel 最終報告：TO HP fix 後重寫版

> 本報告以 2026-03-30 後的 post-fix workspace 為唯一基準，正式取代舊版 [20260328_LOH_evidence_panel_final_report_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md) 中所有 TO downstream 結論。

## 一、開頭重點結論

1. TO HP integer tag fix 影響的是 `step05_intersubmod` 之後的 HP/QC/LOH-like 衍生欄位，不是 caller benchmark，也不是 LongPhase-TO 上游 `LOH.bed` 的定義。
2. paired 主結論維持成立：LOH 在 paired 全域仍略偏 FP enrichment，且在 `effective_hp_reads >= 50` 時最強。
3. TO 主結論正式改寫為：**LOH 在 TO 不是 FP risk，而是 TP enrichment**。這個方向在全域、Tier A(30-49)、Tier A+(>=50) 都一致。
4. HP0 filter 假說仍然不成立。TO LOH-like 中高 HP0 組的 TP% 仍略高於低 HP0 組，但效果比舊版更弱。
5. `LOH + HPMergedSig` 在 paired 中確實對 FP 有高比例訊號，但主要是 `HCC1395 chr8` 主導的樣本特異現象，不應升級成跨樣本通用規則。
6. 所有 LOH 衍生 filter 的最佳情況仍是 `F1_delta <= 0`。LOH 適合作為 annotation / evidence panel 分數，不適合作為 binary filter。

## 二、重寫原因與方法邊界

根據 Knowledge `05_tools/longphase_to.md`、`06_workflows/phasing_workflow.md` 與 `05_tools/InterSubMod.md`：

1. LongPhase-TO 的正式 LOH 來自 phased VCF 與 `LOH.bed`，不是由 BAM 中單一 HP tag 直接定義。
2. 這次 bug 修正的核心，是 TO 端對整數型 HP tag 的讀取，使 `effective_hp_reads`、`hp_assign_rate`、`Quality_Score`、`Potential_LOH`、`LOH_Subtype` 等下游 summary 改變。
3. 因此需要重寫的是 TO downstream 的統計解讀，而不是重做 caller benchmark 或上游 phasing 定義。

本報告的計算口徑：

1. paired / TO benchmark 基線使用 Round 4 post-fix 的 [step1_baseline_f1.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step1_baseline_f1.tsv)。
2. Round 1 全域 LOH enrichment 使用 [loh_enrichment_by_sample_mode.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/loh_enrichment_by_sample_mode.tsv)。
3. Round 2 Tier A(>=30) 與 HP0 方向使用 [core1_tier_enrichment_global.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix/core1_tier_enrichment_global.tsv) 與 [core2_hp0_by_loh_status.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix/core2_hp0_by_loh_status.tsv)。
4. Round 3 HP0 threshold sweep 與 LOH×HPMergedSig 使用 [core1_hp0_threshold_sweep.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix/core1_hp0_threshold_sweep.tsv)、[core2_joint_loh_methyl_sig.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix/core2_joint_loh_methyl_sig.tsv)。
5. Round 4 filter simulation 與 joint table 使用 [step3_joint_loh_hpsig_global.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step3_joint_loh_hpsig_global.tsv)、[step4_full_benchmark_simulation.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step4_full_benchmark_simulation.tsv)。
6. `A(30-49)` / `A+(>=50)` 的 within-tier enrichment 不直接沿用舊 Step 5 表，而是依 [all_region_rows.tsv.gz](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz) 重新手算：
   `enrichment = (FP_loh / FP_tier_total) / (TP_loh / TP_tier_total)`。

## 三、紀錄完整性確認

| 類別 | 已確認紀錄 | 作用 |
| --- | --- | --- |
| Benchmark 基線 | [step1_baseline_f1.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step1_baseline_f1.tsv) | 提供 TP / FP / FN / precision / recall / F1 |
| LOH 覆蓋率 | [step2_loh_coverage.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step2_loh_coverage.tsv) | 驗證 LOH rows 是否足以代表 benchmark universe |
| 全域 LOH enrichment | [loh_enrichment_by_sample_mode.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/loh_enrichment_by_sample_mode.tsv) | Round 1 paired vs TO 主方向 |
| Tier 與 HP0 | [core1_tier_enrichment_global.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix/core1_tier_enrichment_global.tsv)、[core2_hp0_by_loh_status.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix/core2_hp0_by_loh_status.tsv) | Round 2 support quality 與 HP0 方向 |
| HP0 filter 與 LOH×甲基化 | [core1_hp0_threshold_sweep.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix/core1_hp0_threshold_sweep.tsv)、[core2_joint_loh_methyl_sig.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix/core2_joint_loh_methyl_sig.tsv) | Round 3 假說驗證 |
| Deep validation | [step3_joint_loh_hpsig_global.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step3_joint_loh_hpsig_global.tsv)、[step3_fp_per_sample_decomp.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step3_fp_per_sample_decomp.tsv)、[step3_chr_distribution_fp.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step3_chr_distribution_fp.tsv) | 樣本特異性與 chr8 集中 |
| Filter 可行性 | [step4_full_benchmark_simulation.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step4_full_benchmark_simulation.tsv) | 驗證所有變體 F1 delta 是否仍為負 |

結論：本報告中的每一條主結論都有對應 TSV，且手算段落已明確寫出分母與來源。

## 四、Benchmark 基線與 LOH 覆蓋率

| sample | truth_total | TP | FP | FN | F1 | TP coverage | FP coverage |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| HCC1395 | 39,447 | 29,752 | 544 | 9,695 | 0.8532 | 0.9996 | 1.1526 |
| COLO829 | 43,192 | 35,184 | 2,273 | 8,008 | 0.8725 | 0.9939 | 0.9872 |
| H1437 | 81,016 | 67,460 | 8 | 13,556 | 0.9087 | 1.0001 | 1.0000 |
| H2009 | 142,091 | 132,879 | 85 | 9,212 | 0.9662 | 1.0002 | 1.0118 |
| HCC1937 | 16,867 | 12,392 | 195 | 4,475 | 0.8414 | 1.0000 | 1.0000 |
| HCC1954 | 23,030 | 17,864 | 29 | 5,166 | 0.8731 | 1.0025 | 1.0000 |
| HCC1395_DORADO | 39,447 | 29,877 | 238 | 9,570 | 0.8590 | 1.0002 | 1.0084 |

解讀：

1. 這一輪所有 paired 主結論都仍建立在完整 benchmark universe 上，不是舊 Round 3 那種「未計入 FN 的子集 F1」。
2. LOH rows 對 benchmark TP 的覆蓋率幾乎都是 100%，因此 LOH 統計足以支撐後續 paired/TO 對照。
3. `HCC1395` 的 `FP coverage > 1`，表示 LOH rows 中有額外 FP 落在 benchmark truth scope 之外；這是樣本特異的 coverage 注意事項，不改變主結論，但要避免拿它做跨樣本過度延伸。

## 五、四輪證據鏈重寫

### 5.1 Round 1：paired 與 TO 的全域方向正式分流

由 [loh_enrichment_by_sample_mode.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/loh_enrichment_by_sample_mode.tsv) 加總：

| mode | tp_total | fp_total | tp_loh | fp_loh | TP LOH frac | FP LOH frac | FP/TP enrichment |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| paired | 325,270 | 3,429 | 95,404 | 1,201 | 0.2933 | 0.3502 | 1.1941 |
| to | 291,310 | 128,382 | 129,589 | 45,953 | 0.4448 | 0.3579 | 0.8046 |

解讀：

1. paired：LOH 在 FP 中仍略富集，這是後續 paired evidence panel 的基礎。
2. TO：方向和 paired 相反，LOH 在 TO TP 中更常見，不能再把 TO LOH 當 FP risk。
3. 這裡的 TO 方向不是單一樣本例外，而是 post-fix 後全域加總的正式基準。

### 5.2 Round 2：support quality 與 HP0 的正確定位

由 [core1_tier_enrichment_global.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix/core1_tier_enrichment_global.tsv)：

| mode | tier | TP LOH frac | FP LOH frac | FP/TP enrichment | 解讀 |
| --- | --- | ---: | ---: | ---: | --- |
| paired | A(>=30 合併) | 0.3140 | 0.3671 | 1.1689 | paired 高 support LOH 仍偏 FP enrichment |
| paired | B | 0.3862 | 0.3480 | 0.9011 | 中間 support 不穩定 |
| paired | C | 0.6751 | 0.7000 | 1.0369 | 幾乎無區分力 |
| to | A(>=30 合併) | 0.4291 | 0.3299 | 0.7688 | TO 高 support LOH 為 TP enrichment |

由 [core2_hp0_by_loh_status.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix/core2_hp0_by_loh_status.tsv)：

| mode | truth | LOH-like | n | hp0_ratio_mean |
| --- | --- | --- | ---: | ---: |
| paired | TP | False | 229,866 | 0.0897 |
| paired | TP | True | 95,404 | 0.0411 |
| paired | FP | False | 2,228 | 0.0887 |
| paired | FP | True | 1,201 | 0.0660 |
| to | TP | False | 161,721 | 0.0456 |
| to | TP | True | 129,589 | 0.0957 |
| to | FP | False | 82,429 | 0.0466 |
| to | FP | True | 45,953 | 0.1048 |

解讀：

1. paired LOH-like 的 HP0 比 non-LOH-like 低，符合較乾淨的 one-sided phasing。
2. TO LOH-like 的 HP0 比 non-LOH-like 高，方向仍成立，但這只能說明 phasing completeness 差異，不能直接推出「高 HP0 = 假 LOH」。
3. Round 2 的正確定位應是：paired LOH 有 support quality 訊號；TO LOH 則是 TP enrichment 加上較高 HP0 的混合現象。

### 5.3 Round 3：HP0 filter 仍被否定，但效果幅度縮小

由 [core1_hp0_threshold_sweep.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix/core1_hp0_threshold_sweep.tsv)：

| hp0_thresh | low TP% | high TP% | diff |
| ---: | ---: | ---: | ---: |
| 0.02 | 0.745 | 0.764 | +0.019 |
| 0.05 | 0.748 | 0.762 | +0.014 |
| 0.10 | 0.749 | 0.764 | +0.015 |
| 0.15 | 0.749 | 0.767 | +0.018 |
| 0.30 | 0.750 | 0.766 | +0.015 |

解讀：

1. 舊版「高 HP0 反而 TP% 稍高」的方向仍存在。
2. 但 post-fix 後差距只有 1.4-1.9 個百分點，已不足以支持把 HP0 當作 TO LOH quality filter。
3. 因此 Round 3 的正確結論是：HP0 是 diagnostics，不是 filter。

### 5.4 Round 3 / Round 4：LOH × HPMergedSig 的 paired 訊號成立，但具有樣本特異性

由 [step3_joint_loh_hpsig_global.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step3_joint_loh_hpsig_global.tsv)：

| mode | LOH | HPMergedSig | TP | FP | FP% |
| --- | --- | --- | ---: | ---: | ---: |
| paired | True | True | 1,346 | 80 | 0.0561 |
| paired | True | False | 82,826 | 634 | 0.0076 |
| paired | False | True | 102,234 | 620 | 0.0060 |
| paired | False | False | 81,622 | 611 | 0.0074 |
| to | True | True | 1,343 | 505 | 0.2733 |
| to | True | False | 109,831 | 36,286 | 0.2483 |
| to | False | True | 80,884 | 42,142 | 0.3425 |
| to | False | False | 67,014 | 32,578 | 0.3271 |

由 [step3_fp_per_sample_decomp.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step3_fp_per_sample_decomp.tsv) 與 [step3_chr_distribution_fp.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step3_chr_distribution_fp.tsv)：

1. `paired LOH=True + HPSig=True` 的 80 個 FP 中，`HCC1395 = 70`，占 `87.5%`。
2. 這 80 個 FP 中，`chr8 = 66`，占 `82.5%`。

解讀：

1. paired 的 `LOH + HPMergedSig` 的確是高 FP 比例組合，但主要反映 `HCC1395 chr8` 的特殊結構。
2. TO 中 LOH 反而把 FP% 往下拉，這再次說明 TO LOH 不是 FP risk。
3. 因此 `LOH+HPMergedSig` 可作 paired 的樣本特異 annotation，不應寫成跨樣本通用 feature。

## 六、Round 4 最終驗證：LOH 不適合作 binary filter

由 [step4_full_benchmark_simulation.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix/step4_full_benchmark_simulation.tsv) 取每個樣本最佳 filter：

| sample | best_filter | baseline_f1 | new_f1 | f1_delta | rm_tp_pct | rm_fp_pct |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| COLO829 | TierAplus_LOH_HPSig | 0.8725 | 0.8725 | 0.0000 | 0.0000 | 0.0000 |
| H1437 | TierAplus_LOH_HPSig | 0.9087 | 0.9080 | -0.0006 | 0.0012 | 0.0000 |
| H2009 | TierAplus_LOH_HPSig | 0.9662 | 0.9641 | -0.0021 | 0.0042 | 0.0353 |
| HCC1395 | TierAplus_LOH_HPSig | 0.8532 | 0.8491 | -0.0041 | 0.0099 | 0.1213 |
| HCC1395_DORADO | TierAplus_LOH_HPSig | 0.8590 | 0.8573 | -0.0017 | 0.0034 | 0.0126 |
| HCC1937 | TierAplus_LOH_HPSig | 0.8414 | 0.8406 | -0.0009 | 0.0019 | 0.0051 |
| HCC1954 | TierAplus_LOH_HPSig | 0.8731 | 0.8726 | -0.0005 | 0.0010 | 0.0690 |

解讀：

1. 最佳 filter 仍全部是 `TierAplus_LOH_HPSig`，但除了 `COLO829 = 0.0000` 外，其餘樣本一律是負 delta。
2. 這代表 LOH 特徵最多只能用於 risk scoring / annotation，不足以作 hard filter。
3. 這個結論在 paired 與整體四輪研究中保持一致，且不受這次 TO HP fix 推翻。

## 七、A(30-49) / A+(>=50) 的正確 within-tier feature 定位

舊版 Round 4 的 Step 5 以全域 TP/FP 作分母，不適合直接當 within-tier feature 解讀。這裡依 [all_region_rows.tsv.gz](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz) 手算：

| mode | tier | tier TP total | tier FP total | tier TP LOH | tier FP LOH | TP LOH frac | FP LOH frac | within-tier FP/TP enrichment |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| paired | A(30-49) | 48,088 | 1,065 | 23,854 | 227 | 0.4960 | 0.2131 | 0.4297 |
| paired | A+(>=50) | 219,940 | 880 | 60,318 | 487 | 0.2742 | 0.5534 | 2.0179 |
| to | A(30-49) | 44,177 | 25,255 | 27,111 | 10,942 | 0.6137 | 0.4333 | 0.7060 |
| to | A+(>=50) | 214,895 | 86,256 | 84,063 | 25,849 | 0.3912 | 0.2997 | 0.7661 |

這張表是本輪最重要的 feature 定位依據。

解讀：

1. paired `A(30-49)` 是 TP-protective，不是 FP risk。
2. paired `A+(>=50)` 才是最強、最乾淨的 LOH FP risk feature。
3. TO 無論 `A(30-49)` 或 `A+(>=50)`，方向都仍是 TP enrichment。
4. 所以 evidence panel 若要保留 TO feature，應把它視為「弱 TP support / 非 FP risk」，而不是 paired feature 的鏡像版本。

## 八、Feature 清單與升級判定

| feature | 定義 | paired 定位 | TO 定位 | 是否建議升級 |
| --- | --- | --- | --- | --- |
| `paired_tier_a_loh` | `30 <= effective_hp_reads < 50` 且 `core_loh_like` | TP-protective (`0.43x`) | 不適用 | 不建議當 FP risk，僅可當 annotation |
| `paired_tier_aplus_loh` | `effective_hp_reads >= 50` 且 `core_loh_like` | 強 FP risk (`2.02x`) | 不適用 | 可保留在 paired evidence panel |
| `paired_tier_aplus_loh_hpsig` | `paired_tier_aplus_loh` 且 `HPMergedSig` | 樣本特異，偏 HCC1395 chr8 | 不適用 | 不建議升級成全域規則 |
| `to_tier_a_loh` | TO `30 <= effective_hp_reads < 50` 且 `core_loh_like` | 不適用 | TP enrichment (`0.706x`) | 僅能當弱 TP annotation |
| `to_tier_aplus_loh` | TO `effective_hp_reads >= 50` 且 `core_loh_like` | 不適用 | TP enrichment (`0.766x`) | 不可當 FP risk |
| `hp0_ratio` | `NHP0 / NumReads` | diagnostics | diagnostics | 不建議升級為 filter |

## 九、最終結論

1. **paired**：LOH 仍有價值，但正確用法是 annotation / evidence panel，不是 binary filter。
2. **TO**：TO HP fix 後，LOH 的正式定位已經收斂為 TP enrichment。舊版把 TO A/A+ 當近中性或 FP risk 的結論必須作廢。
3. **方法學**：若後續還要引用 `A(30-49)` / `A+(>=50)`，必須使用 within-tier 分母，不可再直接沿用舊 Step 5 全域分母。
4. **研究紀錄**：本報告已把來源 TSV、手算公式、paired/TO 差異與 superseded 關係寫清楚，可直接作為後續 rewrites 的主依據。

## 十、來源索引

1. 舊 final report：[20260328_LOH_evidence_panel_final_report_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md)
2. TO HP fix 影響評估：[20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md)
3. 重跑驗證與分級：[20260330_TO_HP_fix重跑驗證與validated文件修正分級_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260330_TO_HP_fix重跑驗證與validated文件修正分級_01.md)
4. TO LOH post-fix 主報告：[20260330_TO_LOH_enrichment_post_hp_fix_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md)
5. Round 1 post-fix workspace：[20260330_loh_round1_cross_sample_audit_post_to_hp_fix](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix)
6. Round 2 post-fix workspace：[20260330_loh_round2_support_hp0_analysis_post_to_hp_fix](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix)
7. Round 3 post-fix workspace：[20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix)
8. Round 4 post-fix workspace：[20260330_loh_round4_final_validation_post_to_hp_fix](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix)
