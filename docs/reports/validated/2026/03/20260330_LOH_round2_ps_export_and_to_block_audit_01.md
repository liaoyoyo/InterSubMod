<!--
建立時間: 2026-03-30 02:20
更新時間: 2026-03-30 02:20
狀態: validated
目標: 延續 LOH Round 1 v2 的 5.2 與 5.1，確認 paired 的 PS export 邊界、TO-only FP 是否應升級為 PS-block linkage 問題，並整理下一輪需要的決策
處理範圍:
  - 5.2 paired variant-level PS export feasibility
  - 5.2 paired read-level PS fallback pilot
  - 5.1 TO-only FP within-PS block audit
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260329_LOH_round1_cross_sample_audit_v2_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round2_ps_export_and_to_block_audit.py
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit
-->


<!-- PATH WARNING (2026-04-11 驗證):
本報告中部分絕對路徑已失效。類別: 舊 longphase-s 文件行號引用
詳細路徑清單見 docs/data_specs/20260411_path_inventory.tsv
報告結論仍有效，但引用的外部路徑可能需要更新。
-->

# LOH Round 2：paired PS export feasibility 與 TO-only FP PS-block audit

## 一、開頭重點結論

1. paired Round 1 的問題不是「本來有 somatic variant-level PS，但 split 後被 merge 掉」，而是目前 paired somatic artifact 本身就沒有 variant-level `PS`。7 個樣本合計：
   - `germline_phased_merged.vcf.gz`：`12,194,373 / 12,266,176` records 有非空 `PS`（`99.41%`）
   - `somatic_pass.vcf.gz`：`0 / 568,080`
   - `filtered_snv_tp.vcf.gz`：`0 / 325,507`
   - `filtered_snv_fp.vcf.gz`：`0 / 3,458`

2. paired 的 read-level `PS` fallback 現在就能做，但不能寫成「所有樣本都一樣完整」。
   - pilot 共抽 `138` 個 loci（固定 seed=`20260330`；每個 sample × truth 最多 `10` 個 loci；`H1437 FP` 只有 `8` 個）
   - 整體 median `pileup_ps_fraction`：`TP=0.972`、`FP=0.881`
   - 但 `H1437` 例外非常明顯：`TP median=0.000`，且只有 `4/10` loci 有任何 `PS` read

3. TO-only FP 已經足夠升級成 PS-block 問題，而不是單點問題。
   - `TO-only FP = 126,865`
   - 其中 `120,298` 個有 `PS`（`94.82%`）

4. 這輪沒有支持「TO-only FP 全域上比 TO TP 更集中在少數 block」這個簡化說法。
   - `COLO829`、`H1437`、`HCC1395_DORADO`、`HCC1954` 的 `top5 block share` 都是 `TO-only FP <= TO TP`
   - 但同時也看到大量 `multi-locus blocks` 長期帶有 FP 負荷，所以更準確的說法是：
     - TO-only FP 是 **block-linked**
     - 但不一定是 **比 TP 更極端集中**

5. `GT/GT2/GT3` collapse 後，TO-only FP 的 pseudo-state 不是單一型態：
   - sample-level 的主要非 `Other` pattern 多數仍是 `HP2-1-like`
   - 但 `HCC1395_DORADO` 的 sample-level dominant non-`Other` state 已經是 `LOH hap2-side`
   - 而 `COLO829 / H1437 / HCC1395_DORADO` 的 top FP-load blocks 也常由 `LOH hap2-side` 主導

---

## 二、本輪研究問題、假設與驗收條件

### 2.1 研究問題

1. 5.2：paired 為何在 Round 1 看不到 variant-level `PS`，這是 artifact 遺失、工具輸出邊界，還是可以立即補回的 merge 問題？
2. 5.2：若 paired variant-level `PS` 現在做不到，read-level `PS` summary 是否足夠作為 interim schema？
3. 5.1：`TO-only FP` 是否真的在 `PS block` 內有可追蹤的群聚，值得升級成下一輪主問題？

### 2.2 本輪假設

1. paired upstream 還保有 germline phasing 的 `PS`，但 somatic split artifact 沒有 propagate。
2. paired tagged BAM 仍有 read-level `PS`，因此可以先做 `variant_key -> read-level PS summary`。
3. TO-only FP 多數應該能對到 `PS`，且不少 block 會呈現 multi-locus FP 負荷。

### 2.3 成功條件

1. 能用數據確認 paired `PS` 是在哪一層消失。
2. 能量化 paired read-level `PS` fallback 的可行性，而不是只靠工具文件推論。
3. 能量化 TO-only FP 的 `with_ps` 覆蓋率、block 分布與 pseudo-state 組成。

### 2.4 失敗條件

1. paired upstream 與 split artifact 都沒有 `PS`，導致無法區分是工具邊界還是 pipeline 問題。
2. paired tagged BAM 在大多數 loci 也沒有 `PS`，導致 fallback 不可行。
3. TO-only FP 幾乎沒有 `PS`，或只剩單點事件，無法支撐 block linkage 研究。

---

## 三、資料、腳本與輸出路徑

### 3.1 主要輸入

1. Round 1 主表：[all_region_rows.tsv.gz](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz)
2. paired canonical run bundle（由 Round 1 `source_vcf_file` / `source_tagged_bam` 反推）
3. TO phased fields（直接使用 `all_region_rows.tsv.gz` 中已 ingest 的 `to_phase_ps`, `to_phase_gt`, `to_phase_gt2`, `to_phase_gt3`, `to_phase_ps2`, `to_phase_ps3`）

### 3.2 本輪腳本

1. [build_loh_round2_ps_export_and_to_block_audit.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round2_ps_export_and_to_block_audit.py)

### 3.3 本輪正式輸出

1. round2 output dir：[20260330_loh_round2_ps_export_and_to_block_audit](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit)
2. paired variant-level PS audit：[paired_variant_ps_audit.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit/paired_variant_ps_audit.tsv)
3. paired read-level PS pilot summary：[paired_read_ps_pilot_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit/paired_read_ps_pilot_summary.tsv)
4. TO locus table：[to_locus_status_table.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit/to_locus_status_table.tsv)
5. TO block summary：[to_ps_block_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit/to_ps_block_summary.tsv)
6. TO block concentration summary：[to_ps_block_concentration_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit/to_ps_block_concentration_summary.tsv)
7. pseudo-state summary：[to_only_fp_pseudo_state_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit/to_only_fp_pseudo_state_summary.tsv)
8. decision ledger：[decision_ledger.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit/decision_ledger.tsv)

---

## 四、方法、參數與判讀邏輯

### 4.1 paired 5.2：variant-level `PS` audit

每個 paired sample 固定檢查四個 artifact：

1. `germline_phased_merged.vcf.gz`
2. `somatic_pass.vcf.gz`
3. `filtered_snv_tp.vcf.gz`
4. `filtered_snv_fp.vcf.gz`

每個檔案都重算：

1. `total_records`
2. `ps_header_present`
3. `ps_nonempty_records`
4. `ps_nonempty_fraction`

本輪用來判斷 paired 邊界的工具文件依據：

1. [somatic_haplotag.md#L19](/big8_disk/liaoyoyo2001/longphase-s/docs/somatic_haplotag.md#L19) 明寫 LongPhase-S 把 `PS` 存在 **read/BAM tag**
2. [somatic_haplotag.md#L120](/big8_disk/liaoyoyo2001/longphase-s/docs/somatic_haplotag.md#L120) 明寫 somatic result VCF 是以 input tumor VCF 為基礎，只改 filter/header，其他資訊保持不變
3. [README.md#L130](/big8_disk/liaoyoyo2001/longphase-s/README.md#L130) 也再次說明 `PS` 是 read-level tag，不是 somatic VCF 的額外欄位

### 4.2 paired 5.2：read-level `PS` fallback pilot

pilot 設定：

1. `random_seed = 20260330`
2. 每個 `sample × truth_label` 最多抽 `10` 個 loci
3. 真實抽樣結果：
   - 其餘組合皆為 `10`
   - 只有 `H1437 FP = 8`
4. locus 判讀方式：
   - 只看對應 `tagged BAM`
   - 在變異位點做單點 pileup
   - 排除 `unmapped / secondary / supplementary / duplicate`
5. 每個 locus 輸出：
   - `pileup_total_reads`
   - `pileup_ps_reads`
   - `pileup_hp_reads`
   - `pileup_ps_fraction`
   - `pileup_unique_ps`

### 4.3 TO 5.1：same-locus 與 block 定義

same-locus 狀態由 `sample + variant_key` 重建：

1. `paired_status`
2. `to_status`
3. `locus_overlap_class`
4. `paired_to_concordance`

本輪固定使用：

1. `to_only_fp = paired_to_concordance == "to_only_fp"`
2. `PS block = sample + to_phase_ps`
3. block audit 只對 `has_ps == True` 的 loci 做 block 聚合，但會另外統計 `with_ps_frac`

### 4.4 TO pseudo-state collapse

本輪對 `GT / GT2 / GT3` 做 canonicalization：

1. 把 `0/0`, `0/.`, `./0` 這些 slash 形式先正規化成對應的 bar 形式
2. 再做下列 collapse

| canonical `GT` | `GT2` | `GT3` | pseudo-state |
| --- | --- | --- | --- |
| `0|0` | `1|.` | `./.` | `HP1-1-like` |
| `0|0` | `.|1` | `./.` | `HP2-1-like` |
| `0|0` | `1|1` | `./.` | `HP3-like` |
| `.|0` | 任一 LOH 分支組合 | 任一 LOH 分支組合 | `LOH hap1-side` |
| `0|.` | 任一 LOH 分支組合 | 任一 LOH 分支組合 | `LOH hap2-side` |
| 其他 | 其他 | 其他 | `Other/Unknown` |

---

## 五、5.2 paired：結果與觀察

### 5.2.1 variant-level `PS` propagation 結果

| artifact | samples | total_records | ps_nonempty_records | mean_ps_nonempty_fraction |
| --- | ---: | ---: | ---: | ---: |
| `germline_phased_merged` | 7 | 12,266,176 | 12,194,373 | 0.9941 |
| `somatic_pass` | 7 | 568,080 | 0 | 0.0000 |
| `filtered_snv_tp` | 7 | 325,507 | 0 | 0.0000 |
| `filtered_snv_fp` | 7 | 3,458 | 0 | 0.0000 |

![Fig01 paired variant-level PS audit heatmap](assets/20260330_loh_round2_ps_export_and_to_block_audit/fig01_paired_variant_ps_audit_heatmap.png)

`Fig01` 重點：

1. X 軸是 paired 的四個 artifact 層級
2. Y 軸是 7 個 paired sample
3. 色階是 `ps_nonempty_fraction`
4. 每格同時標 `ps_nonempty_records / total_records` 與百分比

結論：

1. paired upstream germline phasing artifact 明確保有 `PS`
2. paired somatic VCF universe 從 `somatic_pass` 開始就已經沒有 `PS`
3. 所以 Round 1 paired 缺 `variant-level PS` 的主因是 **工具輸出邊界 + pipeline artifact 邊界**
4. 這一點比較接近「需要新 export schema」，不是「只要回頭 merge 一次就會回來」

### 5.2.2 read-level `PS` fallback pilot

| sample | truth | pilot_sites | median_reads | median_ps_fraction | loci_with_any_ps_fraction |
| --- | --- | ---: | ---: | ---: | ---: |
| HCC1395 | FP | 10 | 86.0 | 0.9087 | 1.000 |
| HCC1395 | TP | 10 | 88.0 | 0.9725 | 1.000 |
| HCC1395_DORADO | FP | 10 | 74.5 | 0.6524 | 1.000 |
| HCC1395_DORADO | TP | 10 | 76.0 | 0.9723 | 0.900 |
| COLO829 | FP | 10 | 29.5 | 0.9875 | 0.900 |
| COLO829 | TP | 10 | 34.0 | 0.9563 | 1.000 |
| H1437 | FP | 8 | 72.5 | 0.2595 | 0.625 |
| H1437 | TP | 10 | 70.0 | 0.0000 | 0.400 |
| H2009 | FP | 10 | 95.5 | 0.9105 | 1.000 |
| H2009 | TP | 10 | 112.5 | 0.9592 | 1.000 |
| HCC1937 | FP | 10 | 87.5 | 0.8762 | 1.000 |
| HCC1937 | TP | 10 | 131.5 | 0.9892 | 1.000 |
| HCC1954 | FP | 10 | 62.0 | 0.9176 | 1.000 |
| HCC1954 | TP | 10 | 63.0 | 1.0000 | 1.000 |

![Fig02 paired read-level PS fallback pilot](assets/20260330_loh_round2_ps_export_and_to_block_audit/fig02_paired_read_ps_fallback_pilot.png)

`Fig02` 重點：

1. X 軸是 sample
2. Y 軸是 `pileup_ps_fraction`
3. 紅色是 `FP`，藍色是 `TP`
4. boxplot 顯示 sample-level 分布，散點保留每個 pilot locus

觀察：

1. 多數 paired sample 的 read-level `PS` fallback 明顯可用
2. `H2009`, `HCC1937`, `HCC1954`, `HCC1395` 幾乎都在高 `PS` fraction 區
3. `H1437` 是最重要的例外：
   - `TP median_pileup_ps_fraction = 0.0`
   - `loci_with_any_ps_fraction = 0.4`
4. 因此 paired read-level fallback 應該寫成：
   - **可立即上線做 interim summary**
   - **但必須同時輸出 per-sample completeness metric**
   - **不能假設所有樣本都像 H2009 / HCC1954 一樣完整**

---

## 六、5.1 TO-only FP：結果與觀察

### 6.1 TO-only FP 的 PS 覆蓋率

| sample | total_to_only_fp | with_ps | with_ps_frac | top5_share | ge2_block_share | dominant_non_other_state |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| HCC1395 | 11,512 | 9,425 | 0.8187 | 0.0676 | 0.9811 | `HP2-1-like` |
| HCC1395_DORADO | 11,495 | 11,361 | 0.9883 | 0.0841 | 0.8916 | `LOH hap2-side` |
| COLO829 | 16,193 | 15,898 | 0.9818 | 0.1028 | 0.8236 | `HP2-1-like` |
| H1437 | 13,442 | 13,308 | 0.9900 | 0.1012 | 0.9384 | `HP2-1-like` |
| H2009 | 11,978 | 10,849 | 0.9057 | 0.0548 | 0.9574 | `HP2-1-like` |
| HCC1937 | 12,030 | 10,435 | 0.8674 | 0.0806 | 0.9725 | `HP2-1-like` |
| HCC1954 | 50,215 | 49,022 | 0.9762 | 0.0282 | 0.9694 | `HP2-1-like` |

關鍵解讀：

1. `TO-only FP` 中有 `94.82%` 能直接對到 `PS`
2. 幾乎所有 sample 都顯示大量 multi-locus blocks
3. `ge2_block_share` 在所有 sample 都超過 `82%`
4. 所以這個問題可以明確升級成 block-level linkage 問題

### 6.2 但「比 TP 更集中」不是全域結論

![Fig03 TO-only FP vs TP PS-block concentration](assets/20260330_loh_round2_ps_export_and_to_block_audit/fig03_to_ps_block_concentration_curves.png)

`Fig03` 重點：

1. 每個小圖是一個 focus sample
2. X 軸是依 block 負荷排序後的 `top-k PS blocks`
3. Y 軸是 cumulative share of loci
4. 紅線是 `TO-only FP`，藍線是 `TO TP`

這張圖的主要結論不是「紅線一定壓過藍線」，而是：

1. `COLO829`、`H1437`、`HCC1395_DORADO`：兩條曲線很接近，且多數時候 `TO TP` 稍高
2. `HCC1954`：`TO-only FP` 反而比 `TO TP` 更分散
3. 因此下一輪不能直接寫「TO-only FP 集中在少數 block」
4. 更精確的說法應該是：
   - `TO-only FP` 與 `TP` 一樣都大量落在 multi-locus blocks
   - 但 `TO-only FP` block 內常常同時混有 TP，且 block 本身可能帶有穩定 FP 負荷

### 6.3 block 內 TP/FP mix 的結構

![Fig04 PS blocks carrying TO-only FP](assets/20260330_loh_round2_ps_export_and_to_block_audit/fig04_to_only_fp_block_scatter.png)

`Fig04` 重點：

1. X 軸是同一 block 內的 `TO TP count`
2. Y 軸是同一 block 內的 `TO-only FP count`
3. 點大小是 `block_total_loci`
4. 顏色是該 block 中 `TO-only FP` 的 `LOH-like fraction`

觀察：

1. `COLO829`、`H1437`、`HCC1395_DORADO` 多數大 block 位在 `x > y` 區，也就是 TP 還是比 FP 多
2. 但這些大 block 依然能穩定承載 `150~460` 個 `TO-only FP`
3. `HCC1954` 出現幾個 `FP-dominant` 或 `FP-heavy` block，且 `LOH-like fraction` 很高
4. 所以 block audit 的價值不在於找「純 FP block」，而是找「會持續帶高 FP 負荷的混合 block」

### 6.4 focus sample 的 top block 例子

| sample | phase_block_id | block_total_loci | TP_n | TO-only_FP_n | TO-only_FP_frac | LOH-like_frac | dominant_state |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| HCC1954 | 117485658 | 647 | 231 | 416 | 0.6430 | 0.9471 | `Other/Unknown` |
| COLO829 | 69630732 | 1206 | 707 | 462 | 0.3831 | 0.0498 | `LOH hap2-side` |
| H1437 | 42002203 | 1700 | 1350 | 350 | 0.2059 | 0.0857 | `LOH hap2-side` |
| HCC1395_DORADO | 71213912 | 1232 | 949 | 280 | 0.2273 | 0.2143 | `LOH hap2-side` |

這四個 sample 的意義不同：

1. `HCC1954`：FP 負荷最高，且 block 內 `LOH-like` 也最強，適合做最難場景的 case panel
2. `COLO829 / H1437 / HCC1395_DORADO`：TP 仍然較多，但 `LOH hap2-side` block 可以穩定承載大量 `TO-only FP`

### 6.5 pseudo-state 組成

![Fig05 TO-only FP pseudo-state composition](assets/20260330_loh_round2_ps_export_and_to_block_audit/fig05_to_only_fp_pseudo_state_composition.png)

`Fig05` 重點：

1. X 軸是 sample
2. Y 軸是各 pseudo-state 在 `TO-only FP with PS` 中的 fraction
3. 這張圖已將 `0/0`, `0/.`, `./0` 等 slash 形式 canonicalize 後再做分類

主要觀察：

1. `HCC1395_DORADO`：`LOH hap2-side` 很突出，符合 `Fig04` top blocks 的 dominant state
2. `COLO829 / H1437`：`HP1-1-like + HP2-1-like + LOH hap2-side` 三者共同構成主要成分，其中 `LOH hap2-side` 仍很明顯
3. `H2009 / HCC1937 / HCC1954`：`Other/Unknown` 仍佔大宗，表示目前 `GT/GT2/GT3` collapse 還不能涵蓋所有 TO-only FP 情況
4. 因此 pseudo-state 已經可作為 block annotation，但還不能單獨當成分類器

---

## 七、決策與下一步

### 7.1 這輪可以確認的決策

1. paired 若要做真正的 variant-level same-locus PS analysis，必須新增 **工具或 pipeline 層的 export**。
2. paired 在那之前，應先新增 `variant_key -> read-level PS summary`，至少把下列欄位補進 summary：
   - `pileup_total_reads`
   - `pileup_ps_reads`
   - `pileup_ps_fraction`
   - `pileup_unique_ps`
3. TO 的下一輪應改以 **PS block** 為主單位，不再只看 isolated loci。
4. 報告口徑必須修正成：
   - 正確：`TO-only FP is broadly PS-linked and often carried by multi-locus blocks`
   - 不正確：`TO-only FP is globally more concentrated than TP into a few blocks`

### 7.2 建議的下一輪順序

1. paired：
   - 先把 `read-level PS summary` 補到現有 summary schema
   - 不要等 tool-level somatic PS export 才開始下一輪 paired block diagnostics
2. TO：
   - 優先樣本維持 `HCC1954`、`COLO829`、`H1437`、`HCC1395_DORADO`
   - 對 top FP-load blocks 補 region snapshot / case panel
   - 檢查 block 內 `LOH-like`, `hp_ratio_core`, pseudo-state 是否能組成更可解釋的 evidence panel

### 7.3 目前不建議直接下的結論

1. 不建議寫 paired 只要回補 merge 就能找回 somatic variant-level `PS`
2. 不建議把 read-level `PS` fallback 宣稱成所有樣本都穩定
3. 不建議把 `GT/GT2/GT3` pseudo-state 直接升級成 TO-only FP 分類器
4. 不建議把 TO-only FP 寫成「只集中在極少數 block」

---

## 八、來源索引

1. Round 1 v2 報告：[20260329_LOH_round1_cross_sample_audit_v2_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260329_LOH_round1_cross_sample_audit_v2_01.md)
2. Round 1 workspace：[20260327_loh_round1_cross_sample_audit](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit)
3. 本輪腳本：[build_loh_round2_ps_export_and_to_block_audit.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round2_ps_export_and_to_block_audit.py)
4. 本輪 round output：[20260330_loh_round2_ps_export_and_to_block_audit](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260330_loh_round2_ps_export_and_to_block_audit)
5. LongPhase-S read-level PS / output VCF 邊界：[somatic_haplotag.md#L19](/big8_disk/liaoyoyo2001/longphase-s/docs/somatic_haplotag.md#L19)、[somatic_haplotag.md#L120](/big8_disk/liaoyoyo2001/longphase-s/docs/somatic_haplotag.md#L120)、[README.md#L122](/big8_disk/liaoyoyo2001/longphase-s/README.md#L122)
