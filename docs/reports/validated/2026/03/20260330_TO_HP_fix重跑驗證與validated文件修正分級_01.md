<!--
建立時間: 2026-03-30 00:00
目標: 驗證 TO HP integer tag fix 後需重跑項目的新輸出，並判定 2026-03-27 之後 validated 文件的修正優先級
處理範圍: LOH Round 1-4、TO PS-block audit、Round 1 v2 figure set、2026-03-27 之後 validated 報告
關聯檔案:
  - docs/reports/validated/2026/03/20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md
  - docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round2_support_hp0_analysis_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round3_methyl_hp0_filter_01.md
  - docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md
  - docs/reports/validated/2026/03/20260329_LOH_round1_cross_sample_audit_v2_01.md
  - docs/reports/validated/2026/03/20260330_LOH_round2_ps_export_and_to_block_audit_01.md
-->

# TO HP fix 重跑驗證與 validated 文件修正分級

## 一、開頭結論

1. 先前判定需要重跑的 LOH downstream workspace 已全部重跑完成，且新輸出已落到新的 post-fix 目錄。
2. 受影響最重的是 **Round 1 / Round 2 / Round 3 / Round 4 中所有 TO LOH 衍生結論**；paired 主結論基本不變。
3. `20260328_LOH_evidence_panel_final_report_01.md` 與 `20260329_LOH_round1_cross_sample_audit_v2_01.md` 不適合只加勘誤，應視為 **需要新版正式報告**。
4. `20260330_LOH_round2_ps_export_and_to_block_audit_01.md` 的 **paired PS 邊界結論與 block concentration 主結論可保留**，但 block 內 `LOH-like_frac` 與對應圖說已改變，應出內容修正版。
5. `20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md` 與 `20260330_TO_LOH_enrichment_post_hp_fix_01.md` 可保留，仍是目前最可靠的 post-fix 基準文件。

## 二、本次實際重跑與輸出確認

### 2.1 已完成重跑的 workspace / asset

| 類型 | 狀態 | 新輸出 |
| --- | --- | --- |
| Round 1 corrected rerun | 完成 | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix` |
| Round 2 corrected rerun | 完成 | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix` |
| Round 3 corrected rerun | 完成 | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix` |
| Round 4 corrected rerun | 完成 | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix` |
| TO PS-block audit corrected rerun | 完成 | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit_post_to_hp_fix` |
| Round 1 v2 figure set post-fix | 完成 | `/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260330_loh_round1_cross_sample_audit_v2_figures_post_to_hp_fix` |

### 2.2 結構檢查摘要

- Round 2 新 workspace 已寫出 10 個 TSV/JSON/MD 主檔與 7 張圖。
- Round 3 新 workspace 已寫出 10 個主檔與 7 張圖。
- Round 4 新 workspace 已寫出 6 個 TSV 主檔與 8 張圖。
- PS-block audit 新 round 已寫出 paired PS audit、read-level PS pilot、TO locus/block summary、decision ledger 與 5 張圖。
- Round 1 v2 的 6 張圖都已重新輸出。

## 三、重跑後真正改變的內容

### 3.1 Round 1 / Round 2：TO LOH enrichment 方向被正式翻轉

paired 幾乎不變，但 TO 全面改變。

Round 1 `loh_enrichment_by_sample_mode.tsv`：

| sample | TO enrichment 舊值 | TO enrichment 新值 |
| --- | ---: | ---: |
| COLO829 | 1.0100 | 0.9558 |
| H1437 | 1.0170 | 0.9193 |
| H2009 | 1.0049 | 0.9232 |
| HCC1395 | 1.0166 | 0.8965 |
| HCC1395_DORADO | 0.9609 | 0.9092 |
| HCC1937 | 1.0164 | 0.8823 |
| HCC1954 | 1.2242 | 0.8521 |

Round 2 `core1_tier_enrichment_by_sample_mode.tsv` 也出現同樣現象。以 Tier A 為例：

| sample | TO Tier A enrichment 舊值 | TO Tier A enrichment 新值 |
| --- | ---: | ---: |
| COLO829 | 1.0696 | 0.9594 |
| H1437 | 1.1825 | 0.9060 |
| H2009 | 1.0445 | 0.9213 |
| HCC1395 | 1.1043 | 0.8858 |
| HCC1395_DORADO | 1.3171 | 0.9174 |
| HCC1937 | 1.0298 | 0.8740 |
| HCC1954 | 1.0505 | 0.8302 |

因此：

- 任何把 TO Tier A / A+ 寫成「可能是 FP risk」的敘述都不能保留。
- 正確口徑應統一成：**TO LOH 在 post-fix 後是 TP enrichment，不是 FP filter。**

### 3.2 Round 2：HP0 的方向沒翻，但數量基礎與 tier 解讀已變

`core2_hp0_by_loh_status.tsv` 顯示：

- TO `FP / LOH-like` 的 `hp0_ratio_mean`：`0.0814 -> 0.1048`
- TO `TP / LOH-like` 的 `hp0_ratio_mean`：`0.0763 -> 0.0957`
- TO `FP / non-LOH-like` 的 `hp0_ratio_mean`：`0.0504 -> 0.0466`
- TO `TP / non-LOH-like` 的 `hp0_ratio_mean`：`0.0552 -> 0.0456`

這表示：

- 「TO LOH-like 比 non-LOH-like 有更高 HP0」這個方向仍成立。
- 但依賴舊 `Potential_LOH` 分組數量、舊 TO tier 行數、舊圖表的段落，仍需改寫。

### 3.3 Round 3：HP0 filter 結論沒翻，但效果變弱，且 TO joint table 數字已換軸

`core1_hp0_threshold_sweep.tsv`：

- `hp0_thresh=0.05`：舊 `high-low diff = +0.030`，新 `+0.014`
- `hp0_thresh=0.15`：舊 `+0.025`，新 `+0.018`

也就是：

- 舊版「High HP0 反而 TP% 稍高」仍成立，
- 但效果幅度縮小，不應再沿用舊版的 2.5-3.1pp 表述。

Round 4 / Round 3 共用的 TO joint table 也改變。例如 `LOH=True, HPSig=False`：

- 舊：`TP=75,788, FP=30,765, FP%=0.2887`
- 新：`TP=109,831, FP=36,286, FP%=0.2483`

因此 Round 3 若繼續保留舊表與舊附錄，讀者會同時看到「方向不變」與「數字基礎已全部換過」的混合版本，不適合作為正式定稿。

### 3.4 Round 4：TO dual-tier enrichment 明顯改變，Evidence panel 的 TO feature 定位要重寫

`step5_dual_tier_enrichment.tsv`：

| mode | tier | 舊 enrichment | 新 enrichment |
| --- | --- | ---: | ---: |
| to | A(30-49) | 0.8302 | 0.9158 |
| to | A+(>=50) | 1.0272 | 0.6977 |
| to | A_combined(>=30) | 0.9190 | 0.7509 |

paired 的三列維持不變，但 TO 完全不同，尤其：

- 舊版 `A+(>=50)` 曾接近中性甚至略偏 FP，
- 新版明確變成 **TP enrichment 0.6977**。

因此 `20260328_LOH_evidence_panel_final_report_01.md` 中任何把 TO `tier_aplus_loh` 寫成 `1.03×` 或近中性的描述，都已失效。

### 3.5 PS-block audit：總體結論大致穩定，但 block 內 LOH-like fraction 已改

以下檔案完全一致或主數字一致：

- `paired_variant_ps_audit.tsv`
- `paired_read_ps_pilot_summary.tsv`
- `to_ps_block_concentration_summary.tsv`
- `to_only_fp_pseudo_state_summary.tsv`
- `decision_ledger.tsv`

但以下檔案實際有變：

- `to_locus_status_table.tsv`
- `to_ps_block_summary.tsv`
- 圖 `fig02_paired_read_ps_fallback_pilot.png`
- 圖 `fig04_to_only_fp_block_scatter.png`

最重要的是 `to_ps_block_summary.tsv` 中 top block 的 `core_loh_like_to_only_fp_fraction`：

| sample | block | 舊值 | 新值 |
| --- | ---: | ---: | ---: |
| HCC1954 | 117485658 | 0.9471 | 0.9952 |
| COLO829 | 69630732 | 0.0498 | 0.9978 |
| H1437 | 42002203 | 0.0857 | 0.9943 |
| HCC1395_DORADO | 71213912 | 0.2143 | 1.0000 |

因此：

- 「TO-only FP broadly PS-linked / multi-locus block-linked」這個主結論可以保留。
- 但報告第 6.4 節那張 top-block 表和相關解讀必須更新，否則會低估 post-fix 下 block 內 LOH-like 的佔比。

## 四、文件修正分級

### 4.1 A 級：必須出新版，不建議只加勘誤

| 文件 | 判定 | 原因 |
| --- | --- | --- |
| `20260327_LOH_round2_support_hp0_analysis_01.md` | 重寫 / 出 v2 | 內文已混入局部修正，但仍大量連到舊 workspace、舊圖、舊 TO tier 表。正式版不應維持「半修正」狀態。 |
| `20260327_LOH_round3_methyl_hp0_filter_01.md` | 重寫 / 出 v2 | 雖然 HP0 filter 大方向未翻，但 TO 計數、joint table、附錄 key numbers 都已改。 |
| `20260328_LOH_evidence_panel_final_report_01.md` | 重寫 / 出 post-fix final | 它整合前述所有 round；TO feature 定位、Round 4 Step5、文內 feature table 已受實質影響。 |
| `20260329_LOH_round1_cross_sample_audit_v2_01.md` | 重寫 / 出 v3 或 post-fix 版 | 引用舊 Round 1 workspace，且 6 張 v2 圖在 post-fix 全部改變。 |

### 4.2 B 級：主結論可保留，但內容與圖表要修

| 文件 | 判定 | 原因 |
| --- | --- | --- |
| `20260330_LOH_round2_ps_export_and_to_block_audit_01.md` | 內容修正 / 建議出 v2 | paired PS 邊界與 TO block-linkage 結論大致不變，但 6.4 top-block 表、少數圖與所有引用舊 Round 1 路徑都應修正。 |
| `20260327_LOH_round1_cross_sample_audit_01.md` | 保留歷史版，但需明確降級 | 目前已有 2026-03-30 勘誤聲明。若只作歷史記錄可保留；若還要當主參考文件，應補上 post-fix 指向。 |

### 4.3 C 級：可直接保留

| 文件 | 判定 | 原因 |
| --- | --- | --- |
| `20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md` | 保留 | 本身就是 post-fix 影響評估。 |
| `20260330_TO_LOH_enrichment_post_hp_fix_01.md` | 保留 | 本身就是 post-fix TO LOH 正式結論。 |

## 五、建議的修正順序

1. 先出 `20260328_LOH_evidence_panel_final_report_01.md` 的 post-fix 版，因為這是對外最像「總結定稿」的文件。
2. 再補 `20260329_LOH_round1_cross_sample_audit_v2_01.md` 的 post-fix 版，因為圖資已全部重跑完成。
3. 接著重寫 `20260327_LOH_round2_support_hp0_analysis_01.md` 與 `20260327_LOH_round3_methyl_hp0_filter_01.md`，把 inline 勘誤整理成乾淨版本。
4. 最後修 `20260330_LOH_round2_ps_export_and_to_block_audit_01.md`，更新 block-level `LOH-like_frac`、圖檔與 Round 1 path。

## 六、目前可直接引用的新基準

1. Round 1 post-fix workspace：
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix`
2. Round 2 post-fix workspace：
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round2_support_hp0_analysis_post_to_hp_fix`
3. Round 3 post-fix workspace：
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round3_methyl_hp0_filter_post_to_hp_fix`
4. Round 4 post-fix workspace：
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round4_final_validation_post_to_hp_fix`
5. PS-block audit post-fix workspace：
   `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit_post_to_hp_fix`
6. Round 1 v2 post-fix assets：
   `/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260330_loh_round1_cross_sample_audit_v2_figures_post_to_hp_fix`

## 七、最終判斷

這次重跑後，可以把文件分成兩群：

- **需要重寫的**：所有直接把 TO LOH downstream 當成正式主結論的 round2 / round3 / round4 closeout 與 round1 v2。
- **只需修正文稿的**：PS-block audit，因為它的主問題不是 HP integer tag，但 block 內 `LOH-like` 組成確實受影響。

所以答案不是「2026-03-27 之後的文件全部報廢」，而是：

1. **總結性文件與 TO LOH 主報告要換新版**。
2. **paired 為主或 PS 邊界為主的文件，可保留核心結論，但要修局部內容**。
