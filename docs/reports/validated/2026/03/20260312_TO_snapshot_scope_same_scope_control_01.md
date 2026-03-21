<!--
建立時間: 2026-03-12 16:25
目標: 正式驗證 HCC1395 5kHz tumor-only 的 snapshot scope same-scope control，確認 candidate-window subset BAM 是否會扭曲既有 read-level diagnostics 指標
處理範圍:
  - HCC1395 5kHz TO full tagged BAM vs candidate-window subset tagged BAM
  - 與 20260311 TO support 特徵 diagnostics 相同的 18 個代表性 region
  - read-level snapshot 指標一致性檢查
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md
-->

# TO snapshot scope same-scope control

## 破題結論

本輪已完成 `HCC1395 5kHz TO` 的 same-scope control，結論可以正式定案：

> **對同一批 18 個代表性 region 而言，從 full tagged BAM 抽成 candidate-window subset BAM，並不會改變本研究目前使用的 snapshot 指標。**

因此：

1. `subset BAM` 本身不是先前 `5kHz TO` 與 `DORADO TO` read-level 圖不能硬比的主因。
2. `candidate-window subset snapshot` 可以安全用於 **同 dataset 內** 的 read-level ranking、TP/FP 現象分類與 support/artifact diagnostics。
3. 目前跨 `5kHz TO` 與 `DORADO TO` 的 read-level 限制，應更精確改寫成：
   - **不可做跨 dataset 的絕對值硬比較**
   - **但 subset 技術本身在 same-dataset / same-region 控制下未觀察到偏移**

---

## 1. 研究問題

上一輪 `TO support 特徵 read-level diagnostics` 的限制是：

- `HCC1395 5kHz TO` 用的是 full tagged BAM
- `HCC1395_DORADO TO` 用的是 candidate-window subset tagged BAM

因此當時只能寫：

- 可以做方向比較
- 不可做跨 dataset 絕對值硬比

這輪要回答的問題是：

### `candidate-window subset BAM` 是否會系統性扭曲既有 snapshot 指標？

若答案是「會」，那麼先前 `DORADO TO` 的 subset snapshot 就需要大幅降權。  
若答案是「不會」，那麼先前限制就應從「subset 可能失真」改成「dataset/platform 本身不同，因此不可做絕對值硬比」。

---

## 2. 實驗設計

### 2.1 控制策略

採用最低成本、但最直接的 same-scope control：

1. 取用 `20260311_TO_support特徵_readlevel_diagnostics` 中，`HCC1395 5kHz TO` 已經做過的 **相同 18 個代表性 region**
2. 從既有 full tagged BAM 抽出這 18 個 region 的 `candidate-window subset BAM`
3. 用與原先相同的 snapshot 流程重新計算指標
4. 將 `full snapshot` 與 `subset snapshot` 逐欄比較

### 2.2 比較指標

本輪固定比較以下 read-level 指標：

- `read_count`
- `target_depth`
- `target_ref_count`
- `target_alt_count`
- `target_alt_fraction`
- `target_other_count`
- `target_mpileup_depth`
- `top_hp_fraction`
- `hp1_collapsed_count`
- `hp2_collapsed_count`
- `other_hp_count`
- `collapsed_hp_balance_delta`
- `na_hp_count`
- `na_hp_fraction`
- `top_hp_tag`

其中：

- `target_alt_fraction`：該位點在 `mpileup` 中的 ALT read 比例
- `collapsed_hp_balance_delta`：collapse 後 HP1 與 HP2 的不平衡程度，越高代表越偏單側 haplotype
- `na_hp_fraction`：沒有 HP tag 的 read 比例

### 2.3 成功條件

若 subset 沒有引入偏移，預期應看到：

1. 同一 region 的上述指標完全一致，或至少沒有實質偏差
2. `full vs subset` 的 direction 判讀不變
3. 可正式將 subset snapshot 定位為：
   - `safe_for_within_dataset_ranking`

---

## 3. 使用資料與輸出

### 3.1 輸入資料

- full tagged BAM：`/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam`
- 代表性 region 來源：`/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260311_to_support_feature_diagnostics/diagnostic_summary.tsv`

### 3.2 本輪輸出

- [InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_snapshot_bias_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_snapshot_bias_summary.md)
- [InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_metric_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_metric_summary.tsv)
- [InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_snapshot_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_snapshot_comparison.tsv)
- [InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/hcc1395_5khz_to_selected_regions.bed](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/hcc1395_5khz_to_selected_regions.bed)
- [InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/hcc1395_5khz_to_candidate_windows_subset.bam](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/hcc1395_5khz_to_candidate_windows_subset.bam)

### 3.3 體積

- subset BAM：約 `95M`
- subset BAI：約 `461K`
- region BED：約 `771B`

這也代表後續若只為 diagnostics 重現相同 region，使用 subset BAM 的成本遠低於 full `260G` tagged BAM。

---

## 4. 主要結果

### 4.1 同一批 18 個 region 全部完全一致

根據 [same_scope_snapshot_bias_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_snapshot_bias_summary.md)：

- `compared_regions = 18`
- `fully_identical_regions = 18`

也就是說，這次 same-scope control 中，**18/18 個 region 的比較結果都是完全一致**。

### 4.2 所有核心指標的最大差值都是 0

根據 [same_scope_metric_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_metric_summary.tsv)：

| metric | rows_compared | max_abs_delta | mean_abs_delta | all_identical |
| --- | ---: | ---: | ---: | --- |
| `read_count` | 18 | 0.0 | 0.0 | true |
| `target_depth` | 18 | 0.0 | 0.0 | true |
| `target_ref_count` | 18 | 0.0 | 0.0 | true |
| `target_alt_count` | 18 | 0.0 | 0.0 | true |
| `target_alt_fraction` | 18 | 0.0 | 0.0 | true |
| `target_other_count` | 18 | 0.0 | 0.0 | true |
| `target_mpileup_depth` | 18 | 0.0 | 0.0 | true |
| `top_hp_fraction` | 18 | 0.0 | 0.0 | true |
| `hp1_collapsed_count` | 18 | 0.0 | 0.0 | true |
| `hp2_collapsed_count` | 18 | 0.0 | 0.0 | true |
| `other_hp_count` | 18 | 0.0 | 0.0 | true |
| `collapsed_hp_balance_delta` | 18 | 0.0 | 0.0 | true |
| `na_hp_count` | 18 | 0.0 | 0.0 | true |
| `na_hp_fraction` | 18 | 0.0 | 0.0 | true |
| `top_hp_tag` | 18 | - | - | true |

這表示在本輪控制下，subset 並沒有改變：

- read 數量
- ALT fraction
- HP balance
- 未分派 HP 比例
- 主導 HP tag

### 4.3 region 細表沒有任何變動列

根據 [same_scope_snapshot_comparison.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260312_snapshot_scope_same_control_hcc1395_5khz/same_scope_snapshot_comparison.tsv)，每列都有：

- `identical_all_metrics = true`
- 各指標 `*_abs_delta = 0.0`
- `top_hp_tag_changed = false`

因此這不是只在 summary 層看起來一致，而是在 region 層逐列完全一致。

---

## 5. 結果的意義

### 5.1 subset 技術本身沒有扭曲 same-region snapshot

這輪控制已足夠回答一個先前未關閉的技術疑慮：

> **若 region 集合完全相同，從 full tagged BAM 抽出 candidate-window subset BAM，不會改變本研究目前使用的 read-level snapshot 指標。**

因此 `subset BAM` 的主要作用是：

- 降低磁碟與 I/O 成本
- 保留同一 region 的觀察能力

而不是改變 read-level 指標。

### 5.2 先前 `5kHz TO full` vs `DORADO TO subset` 的限制需要改寫

先前的限制寫法偏向：

- `full` 與 `subset` scope 不同，因此 read-level 只能做方向比較

現在更精確的寫法應改成：

1. **subset 萃取本身在 same-dataset / same-region 控制下未造成偏移**
2. 真正還不能做跨 dataset 絕對值比較的原因是：
   - dataset / platform 本身不同
   - 候選集合與 coverage ceiling 不同
   - `PairwiseMedianDist` 等 aggregate feature 方向本來就有 dataset-dependence

### 5.3 對 DORADO TO diagnostics 的影響

這個結果會提升 `DORADO TO subset diagnostics` 的可信度，但不代表所有跨 dataset 限制都解除。

目前可以升級的結論是：

- `DORADO TO` 的 subset snapshot **可作為同 dataset 內** 的 support/artifact diagnostics 與 TP/FP 排序依據

目前仍不能升級的結論是：

- `5kHz TO` 與 `DORADO TO` 的 `depth / alt_fraction / hp_assign_rate` 可以直接跨 dataset 比大小

---

## 6. 對主線研究的影響

本輪 same-scope control 完成後，主線敘述應調整為：

1. `TO dual-model availability`：已 closeout，現有 source tree 為 `pileup-route only`
2. `snapshot scope bias`：已完成 `5kHz TO` same-scope control；subset 技術本身未觀察到偏移
3. 真正尚未完全收尾的，只剩：
   - `Pool B FN integration closeout`

也就是說，先前 closing plan 的第 2 項已從：

- `尚未驗證`

收斂成：

- `5kHz TO same-scope control 已完成；cross-dataset absolute comparison 仍受 dataset/platform 差異限制`

---

## 7. 後續建議

### 現在可以做

1. 在 `TO support` 與四象限整合文件中，將 `subset snapshot` 升級為：
   - `safe_for_within_dataset_ranking`
2. 保留 `cross-dataset absolute comparison not allowed` 的限制，但理由改寫為：
   - dataset/platform 差異
   - coverage ceiling
   - feature direction dataset-dependence

### 現在不需要立即做

1. 不需要為了驗證 subset 技術本身，再強制為 `DORADO TO` 生成 `200G+` full tagged BAM
2. 不需要把 `subset vs full` 再當成主要 blocker

### 下一個真正的 closing item

1. `Pool B FN integration closeout`

---

## 8. 關聯文件

- [docs/plans/2026/03/20260312_未完項closing_plan_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md)
- [docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260312_TO雙模型可得性closeout_01.md)
- [docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260311_TO_support特徵_readlevel_diagnostics_01.md)
- [docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md)
