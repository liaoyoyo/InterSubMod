<!--
建立時間: 2026-03-29
目標: 依使用者對 Fig01-Fig06 的觀察，補上 Round 1 cross-sample audit 的第 2 版說明、事實修正與下一輪建議
處理範圍:
  - 沿用 2026-03-27 Round 1 的同一批 workspace 與圖表
  - 補強 Fig01-Fig06 的量化解讀
  - 回答 HP_Ratio / effective_hp_reads / PS export 的事實問題
  - 定義下一輪 PS-block 與圖表修正方向
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit_v2_figures.py
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260329_loh_round1_cross_sample_audit_v2_figures
-->

# LOH Round 1 Cross-Sample Audit v2 補充報告

> 日期：2026-03-29
> 性質：v2 補充版；保留 v1，這份只補充新的量化解讀與事實修正
> v1 報告：[20260327_LOH_round1_cross_sample_audit_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md)
> 共用 workspace：[20260327_loh_round1_cross_sample_audit](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit)
> 共用腳本：[build_loh_round1_cross_sample_audit.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit.py)
> v2 圖片資產（相對路徑來源）：`assets/20260329_loh_round1_cross_sample_audit_v2_figures/`

---

## 一、開頭重點結論

1. `7` 個樣本 × `2` 個流程的 Round 1 scope 沒有問題，`14/14` datasets 都正確納入；這份 v2 不是重跑 round，而是把既有圖表觀察轉成可追溯的數據解讀。
2. `Fig01` 你指出的 paired-only `FP LOH-like enrichment` 是真訊號，但它是 **sample-specific paired 訊號**，不是全域規則；最穩的是 `H2009`、`HCC1395_DORADO`、`HCC1937`，`HCC1954` 方向很強但 `FP n` 仍偏小。
3. `Fig02` 現在用的是 `density` 疊圖，不是 count，所以看圖會放大或縮小某些印象。改成 count 軸、`TP/FP` 分圖是必要的；不過即使用現有資料量化，paired 的 `0.5` 中心峰確實較常出現在 TP。
4. `Fig03` 的高支持 peak 與單側偏向值得保留。部分樣本在 `0.25/0.75`、`0.1/0.9` 有高支持 TP 峰；同時部分 TO sample 出現明顯單側 edge bias，這支持後續做 `PS block` 內 clustering，而不是只看全域直方圖。
5. `HP_Ratio = 0.5` 不能解讀成 balanced，這點 v1 是對的；但 v2 進一步確認：`effective_hp_reads = 0` 不只是一種情況，TO 還混有大量 **untracked reads**，所以「balanced」與「no effective HP support」可以區分，但 `no support` 內部還不能完全拆乾淨。
6. TO 的 `PS block` 分析已經具備條件。`to_only_fp` 中 `94.82%` 帶 `PS`，且每個 sample 都存在多個「單一 block 內有大量 FP loci」的情形；paired 端缺 `variant-level PS` 則是 **artifact/export 邊界問題**，不是 upstream 完全沒有 `PS`。

---

## 二、沿用範圍與本版定位

本版沿用 v1 的同一批資料：

1. `748,391` 個 region rows
2. `459,782` 個 same-locus paired-vs-TO union rows
3. `5` 個代表案例

本版不改動以下內容：

1. sample scope
2. mode scope
3. `hp_ratio_core` / `core_loh_like` 的主定義
4. v1 的 case panel 與 fixed figures 路徑

本版新增的是：

1. 依 Fig01-Fig06 的量化補充解讀
2. 對 `3.2`, `3.3`, `3.7`, `3.9` 這幾個事實段落的補充回答
3. 下一輪 `PS block` 與圖表改版建議

知識庫交叉來源：

1. Knowledge Base `05_tools/longphase_s.md`
2. Knowledge Base `05_tools/longphase_to.md`
3. Knowledge Base `03_file_formats/vcf_longphase.md`
4. Knowledge Base `03_file_formats/bam_format.md`
5. Knowledge Base `06_workflows/phasing_workflow.md`

---

## 三、Fig01-Fig06 的 v2 量化解讀

### 3.1 Fig01 — LOH-like Fraction Overview

![Fig01 顯示圖片](assets/20260329_loh_round1_cross_sample_audit_v2_figures/fig01_loh_like_fraction_overview.png)

你的核心觀察是對的：在樣本數夠、而且 paired `FP` 母數不至於只剩個位數的 paired datasets 中，`FP LOH-like` 比 `TP LOH-like` 更高的現象，主要發生在 paired。

| sample | paired FP/TP enrichment | paired FP n | TO FP/TP enrichment | TO FP n | v2 解讀 |
| --- | --- | --- | --- | --- | --- |
| H2009 | `2.685×` | `86` | `1.005×` | `11,989` | paired 訊號明確，TO 幾乎被背景吃掉 |
| HCC1395_DORADO | `1.260×` | `240` | `0.961×` | `11,572` | 方向主要出現在 paired，TO 反而略偏 TP |
| HCC1937 | `1.505×` | `195` | `1.016×` | `12,032` | paired 與 TO 明顯不對稱 |
| HCC1954 | `3.185×` | `29` | `1.224×` | `50,218` | paired 方向很強，但 paired `FP n` 仍偏小，需保守解讀 |

補充兩點：

1. `H1437 paired = 1.795×` 方向也一致，但 `FP=8`，不能與 `H2009/HCC1937` 視為同等級證據。
2. 因此 Fig01 最穩定的結論是：
   - paired 確實存在 sample-specific `FP risk` 訊號
   - TO 則多數被高背景 LOH-like 佔比稀釋

### 3.2 Fig02 — HP Ratio Core Distribution

![Fig02 顯示圖片](assets/20260329_loh_round1_cross_sample_audit_v2_figures/fig02_hp_ratio_core_distribution.png)

v2 版面已改成：

1. 左欄 `TP`、右欄 `FP`
2. 同一樣本 `paired` 在上、`to` 在下
3. 同一樣本四格共用相同 `Y` 軸上限，方便單樣本內直接比較

這張圖目前的主要問題，不是你看錯，而是圖本身真的不利解讀。`Fig02` 在 [build_loh_round1_cross_sample_audit.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit.py) 裡是：

1. `TP/FP` 疊圖
2. `density=True`
3. 沒有直接顯示母數

所以 v2 先把你提的「中心峰」改成數量化觀察。

#### paired：`0.5` 中心峰比較容易留在 TP

以 `hp_ratio_core` 落在 `0.45–0.55`、且 `effective_hp_reads > 0` 計：

| sample | paired TP center | paired FP center | v2 解讀 |
| --- | --- | --- | --- |
| H2009 | `18.2%` | `1.2%` | 差異非常明顯 |
| HCC1937 | `9.2%` | `2.1%` | TP 中間峰遠多於 FP |
| HCC1954 | `14.7%` | `6.9%` | 方向一致，但 paired FP 母數小 |
| HCC1395_HKU_5kHz | `11.2%` | `4.4%` | 中心峰主要留在 TP |
| HCC1395_DORADO | `10.8%` | `8.0%` | 方向一致，但對比沒 H2009 明顯 |

這支持你的原始判讀：paired 中央 `0.5` 峰比較像 TP 背景，而不是 paired FP 的主成分。

#### TO：中間峰不是所有 FP 都高，而是 sample-specific

你說「TO 很多 FP 都有中間高峰」這句需要收斂成更精確的版本：

1. `H2009 TO FP`、`HCC1954 TO FP`、`HCC1395_HKU_5kHz TO FP` 的中心峰確實存在
2. 但 `H1437 TO FP`、`HCC1395_DORADO TO FP`、`COLO829 TO FP` 的中心峰其實很弱

若只看 `effective_hp_reads ≥ 20` 的中心峰：

| sample | TO TP center | TO FP center | v2 解讀 |
| --- | --- | --- | --- |
| H1437 | `2.8%` | `0.9%` | 你說的「中間沒出現 FP」成立 |
| HCC1395_DORADO | `3.9%` | `0.9%` | 同樣成立 |
| COLO829 | `2.7%` | `1.6%` | 中心峰整體都很低 |

所以 v2 會把這段改寫成：

1. TO 中間峰不是全域現象
2. 它主要集中在少數 sample
3. `H1437` 與 `HCC1395_DORADO` 的 TO-FP 確實接近「中間缺口」

#### paired vs TO 的 TP 中心分佈差異

你提到 `COLO829`、`H1437`、`HCC1937` 的 TP 在 paired vs TO 間有差異，這點有明顯支持：

| sample | paired TP center | TO TP center |
| --- | --- | --- |
| COLO829 | `19.8%` | `1.2%` |
| H1437 | `11.0%` | `2.0%` |
| HCC1937 | `9.2%` | `5.5%` |

這表示 mode 間差異不是只發生在 FP，TP 的 `hp_ratio_core` 幾何也會跟著改變。

### 3.3 Fig03 — Effective HP vs HP Ratio Scatter

![Fig03 顯示圖片](assets/20260329_loh_round1_cross_sample_audit_v2_figures/fig03_effective_hp_vs_hp_ratio_scatter.png)

#### 圖面問題

你指出的三個問題都成立：

1. X 軸共用原始 `effective_hp_reads` 範圍，深度高的 sample 會壓扁其他樣本
2. 點大小直接綁 `NumReads`，大點會遮住底下分佈
3. `FP` 顏色太深時會蓋掉 `TP`

#### 高支持 TP 峰不只在 `0.5`

若只看 `effective_hp_reads ≥ 50`，部分 paired TP 在 `0.25/0.75` 的肩峰確實明顯：

| sample | `0.20–0.30` or `0.70–0.80` |
| --- | --- |
| HCC1954 paired TP | `24.2%` |
| H1437 paired TP | `19.5%` |
| H2009 paired TP | `12.4%` |
| HCC1395_DORADO paired TP | `11.2%` |

同時，高支持 edge peak 也存在：

| sample | `<0.1 or >0.9` |
| --- | --- |
| H2009 paired TP | `26.9%` |
| HCC1395_DORADO paired TP | `37.7%` |
| HCC1937 paired TP | `54.1%` |

所以你看到的「TP 會在 `0.75/0.25` 或 `0.9/0.1` 出現高支持高峰」不是錯覺。

#### 對稱性問題值得保留

你提到理論上 `hp_ratio 0.5` 上下應該接近對稱，若明顯偏向單側，可能代表 local clustering。這個判斷值得保留，因為 TO 確實有明顯單側偏向。

以 `effective_hp_reads ≥ 20` 的 edge rows 計，右側（`>0.9`）占比：

| sample | TO TP right-side fraction |
| --- | --- |
| HCC1937 | `83.6%` |
| HCC1395_HKU_5kHz | `73.2%` |
| H2009 | `67.0%` |

這類強單側偏向如果不是 plotting artifact，就比較像：

1. phase block 內局部聚集
2. 或同一側 haplotype 的 region-level偏向

而不是完全隨機的遠距離 HP1/HP2 混合。

### 3.4 Fig04 — VerificationClass × LOH Subtype Structure

![Fig04 顯示圖片](assets/20260329_loh_round1_cross_sample_audit_v2_figures/fig04_verificationclass_lohsubtype_structure.png)

這張圖的顏色確實應該重設，因為現在的維度其實是：

1. 主維度：`VerificationClass`
2. 次維度：是否帶 `LOH_*`

數據上，LOH-like 內部的 class 結構也真的不同：

| mode / truth | Noise | Weak | Strong | Subclone |
| --- | --- | --- | --- | --- |
| paired FP within LOH-like | `71.9%` | `11.1%` | `14.7%` | `2.3%` |
| TO FP within LOH-like | `53.4%` | `27.8%` | `15.7%` | `3.2%` |

這表示：

1. paired 的 LOH-like FP 主要還是 `Noise`
2. TO 的 LOH-like FP 則更常擴散到 `Weak/Strong`

所以你建議的配色方式不只是美觀問題，也更符合資料層級。

### 3.5 Fig05 — LOH vs non-LOH Feature Violin Plots

![Fig05 顯示圖片](assets/20260329_loh_round1_cross_sample_audit_v2_figures/fig05_loh_nonloh_feature_boxplots.png)

改成 violin plot 是合理的，尤其對多峰分佈會更好。

你提到 `hp3_ratio` 看起來沒資料、懷疑有問題。v2 的判斷是：

1. **不是單純畫錯**
2. 而是 **欄位本身就不是 paired / TO 對稱的 feature**

根據 Knowledge Base `05_tools/longphase_to.md`、`06_workflows/phasing_workflow.md`：

1. TO 沒有 LongPhase-S 式的 **read-level** `HP3`
2. 但 TO 的 phased VCF 仍可透過 `GT2/GT3` 表達 `hp3-like` 的 ambiguous somatic phasing
3. 實際上本輪 TO 所有 strata 的 `hp3_nonzero_frac = 0`，表示 `hp3_ratio` 在 **summary/BAM 指標層** 是 structural zero
4. paired 雖然有 `HP3`，但多數 strata 的 `hp3_ratio` 中位數仍是 `0`

所以這張圖若要保留 `hp3_ratio`：

1. 至少要把 paired / TO 分開
2. 或直接在 TO panel 標註 `structural zero（summary/BAM layer）`
3. 不應把 TO 的 `hp3_ratio = 0` 解讀成「TO 完全沒有 hp3-like ambiguity」

### 3.6 Fig06 — Sample Bin Heatmap

![Fig06 顯示圖片](assets/20260329_loh_round1_cross_sample_audit_v2_figures/fig06_sample_bin_heatmap.png)

你說「不知道程度表示的意義」，這點是合理的。因為現在 heatmap 的顏色不是 effect size，而是：

1. 某個 sample/mode 在該 bin 的 **fraction**
2. cell 文字另外補了該 bin 的絕對 `n`
3. 也不是「距離正常值多遠」

v2 已改成較接近你要求的分法：

1. `hp_ratio_core`: `no_eff`, `<0.05`, `0.05-0.10`, `0.10-0.20`, `0.20-0.40`, `0.40-0.60`, `0.60-0.80`, `0.80-0.90`, `0.90-0.95`, `>0.95`
2. `effective_hp_reads`: `0`, `1-9`, `10-29`, `30-49`, `50-74`, `75-99`, `100-149`, `>=150`

因此現在這張圖能同時回答：

1. 邊緣極值與中間區段是否由少數 sample/mode 主導
2. `75x` 左右是否接近該 sample/mode 的常見 read support 區間
3. 某個高比例 cell 到底是不是只由很小的 `n` 撐起來

但它仍然回答不了：

1. 每個 sample 的「正常 read depth」是否應該以 sample-normalized 基準重切
2. 倍體、copy-number 與 read support 的對應是否需要一起納入

所以下一輪若要再升級，應該往 sample-normalized / ploidy-aware bins 走，而不是再回到單純固定分箱。

---

## 四、關於事實敘述的 v2 回答

### 4.1 `HP_Ratio` 本身不能單獨解讀：v2 的完整答案

v1 的核心結論是正確的，但 v2 要把它拆得更細。

#### `effective_hp_reads = 0` 到底是什麼

本輪共有 `69,807` 個 `effective_hp_reads = 0` rows，占全部 `748,391` rows 的 `9.33%`。

拆開後：

| bucket | count | 說明 |
| --- | --- | --- |
| paired `hp0_and_hp3` | `27,198` | paired 中沒有 `HP1/HP2 family`，但有 `HP0` + `HP3` |
| paired `hp3_only` | `2,810` | paired 中只有 `HP3` |
| TO `hp0_only` | `8,875` | TO 中只有顯式 `HP0` |
| TO `neither` | `30,924` | `HP1FamilyN=HP2FamilyN=NHP0=NHP3=0`，但 `NumReads > 0` |

因此答案不是：

1. 全部 unphased
2. 全部 HP3

而是：

1. paired 主要是 `HP0/HP3` 組成
2. TO 則大量混有未被目前 summary 顯式追蹤的 reads

#### `HP1-1` / `HP2-1` 是否已併入主 family

是，已經有明確確認。

根據：

1. Knowledge Base `05_tools/longphase_s.md`
2. Knowledge Base `03_file_formats/bam_format.md`
3. [LabelTest.hpp](/big7_disk/liaoyoyo2001/InterSubMod/include/core/LabelTest.hpp)
4. [LabelTest.cpp](/big7_disk/liaoyoyo2001/InterSubMod/src/core/LabelTest.cpp)
5. [build_loh_round1_cross_sample_audit.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit.py)

本輪主定義是：

1. `HP1 family = HP1 + HP1-1`
2. `HP2 family = HP2 + HP2-1`
3. `HP3` 與 `HP0` 不進 `hp_ratio_core` 分母

#### 能不能區分真正 balanced 與 no effective HP support

能，但只到第一層：

1. `balanced-like`：`effective_hp_reads > 0` 且 `hp_ratio_core` 靠近 `0.5`
2. `no effective HP support`：`effective_hp_reads = 0`

不能完全區分的是 `no support` 的內部細類，尤其在 TO。

v2 的新增判讀是：

1. paired 的 `effective_hp_reads = 0` 還算可解釋，因為 `untracked_reads = 0`
2. TO 的 `effective_hp_reads = 0` 很多其實是 `untracked_reads > 0`

這裡的 `untracked_reads = NumReads - (HP1FamilyN + HP2FamilyN + NHP0 + NHP3)`。

這是 v2 從 summary 欄位 **推論** 出來的 schema 缺口，不是 LongPhase-TO 官方明示欄位。

### 4.2 paired 與 TO 的 LOH-like 整體型態不同：是因為 TO 難 phase 還是還有別的因素

答案是：**不只難 phase，還有工具語意與 export schema 的差異。**

可以拆成三層：

1. **support quality 差異**
   - TO 的 `effective_hp_reads < 10` 本來就更高
   - TO 很多 loci 的 HP support 明顯較弱
2. **label space 差異**
   - paired 有 `HP1-1 / HP2-1 / HP3`
   - TO 的 BAM 只有 `HP:i:1 / HP:i:2`，但 phased VCF 另外提供 `GT2 / GT3 / PS2 / PS3` 來表達 somatic branch 與 LOH resolution
   - 因此 TO 缺的是 LongPhase-S 式的 **read-level** `HP1-1 / HP2-1 / HP3` universe，不是完全沒有 `1-1 / 2-1 / 3-like` 資訊
3. **summary/export schema 差異**
   - TO 存在大量 `untracked_reads`
   - paired 沒有 `variant-level PS` artifact

因此 v2 會把這句話收斂成：

1. `paired` 的確較適合幫忙定義「可解釋 evidence panel」
2. `TO` 較適合幫忙驗證「最難的 FP 場景」
3. 但兩者差異不是單一原因，不能只寫成「TO 難 phase」

### 4.3 same-locus paired-vs-TO 的主要不一致不是 `both_fp`，而是 `TO-only FP`

這個敘述 v1 正確，v2 更強化它。

除了總量上 `to_only_fp = 126,865` 本來就遠高於 `both_fp = 1,517`，v2 再補兩個重要點：

1. `to_only_fp` 中 `94.82%` 具有 `PS`
2. 各 sample 都有大量 multi-locus `PS block`

所以 `TO-only FP` 不只是「TO 比 paired 多一些錯」，而是：

1. mode-specific 的主要 failure universe
2. 已經足以升級成 block-level linkage 問題

### 4.4 為何 paired 仍卡在沒有 variant-level PS export

v2 的答案是：

1. paired 上游不是完全沒有 `PS`
2. 真正卡住的是 **Round 1 使用的 paired somatic split VCF 沒把 `PS` 留下來**

具體證據：

1. Knowledge Base `05_tools/longphase_s.md` 指出 LongPhase-S 需要帶 `PS` 的 phased germline VCF 作輸入。
2. 但 Round 1 paired 軸 ingest 的是：
   - [filtered_snv_tp.vcf.gz](/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix/longphase_s/filtered_snv_tp.vcf.gz)
   - [filtered_snv_fp.vcf.gz](/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix/longphase_s/filtered_snv_fp.vcf.gz)
3. 以 `H2009 paired TP` 抽查前 `1000` 筆，`PS non-empty = 0/1000`
4. 相對地，TO 的 [tumor_phased.vcf](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_phased.vcf) 前 `1000` 筆中，`PS non-empty = 379/1000`

因此 paired 的限制目前比較像：

1. phase 資訊在 upstream phasing universe 存在
2. 但沒有被 propagate 到 Round 1 直接合表使用的 variant-level artifact

這也是為什麼 v2 會把「paired 補 variant-level PS export」列為下一輪必要前置。

---

## 五、下一輪建議

### 5.1 研究面

1. 直接做 `TO-only FP` 的 `within-PS block` clustering audit
2. 優先樣本：`HCC1954`、`COLO829`、`H1437`、`HCC1395_DORADO`
3. 優先回答：
   - 同一 `PS block` 內 FP 是否群聚
   - block 內 TP/FP mix ratio 是否偏向 FP
   - `LOH-like` / `VerificationClass` / `hp_ratio_core` 是否共同偏向某一側

### 5.2 schema / export 面

1. paired 補 `variant-level PS export`
2. 或至少先做 `variant_key -> read-level PS summary`
3. 在 summary 表新增：
   - `tracked_reads`
   - `untracked_reads`
   - `untracked_ratio`

### 5.3 圖表面

1. `Fig03`：若後續資料量再上升，可改成 `hexbin` 或 density contour，避免點雲互相遮蓋
2. `Fig04`：可再拆成 paired / TO 分面，減少 `sample|mode|truth` 全堆在同一張的擁擠感
3. `Fig05`：下一輪可直接拆 paired / TO 兩張，讓 `hp3_ratio` structural-zero 的差異更直觀
4. `Fig06`：下一輪優先改成 sample-normalized 或 ploidy-aware bins，而不是再細切固定區間

---

## 六、最終判斷

這份 v2 補充後，Round 1 可以更清楚地寫成：

1. `HP_ratio`、`LOH-like`、`effective_hp_reads` 都是重要 diagnostics，但還不能單獨決策
2. paired 與 TO 會互相輔助，但輔助方向確實不對稱
3. paired 較適合定義「可解釋 evidence panel」
4. TO 較適合驗證「這組 evidence 是否真能處理最難的 FP 場景」
5. 下一步最值得升級的，不是再做一輪全域 density plot，而是 `TO-only FP` 的 `PS-block linkage audit` 與 paired 的 `PS export`

這個方向與你給的總結是一致的，而且目前有足夠數據支持，適合當後續 round 的正式起點。

---

## 七、來源索引

1. v1 報告：[20260327_LOH_round1_cross_sample_audit_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md)
2. Round 1 workspace：[20260327_loh_round1_cross_sample_audit](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit)
3. 主表：[all_region_rows.tsv.gz](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz)
4. same-locus 表：[same_locus_compare.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/same_locus_compare.tsv)
5. enrichment 表：[loh_enrichment_by_sample_mode.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/loh_enrichment_by_sample_mode.tsv)
6. QC 表：[hp_coverage_qc_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/hp_coverage_qc_summary.tsv)
7. subtype 表：[verificationclass_by_loh_subtype.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/verificationclass_by_loh_subtype.tsv)
8. feature 表：[loh_vs_methyl_feature_summary.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/loh_vs_methyl_feature_summary.tsv)
9. heatmap 表：[sample_bin_heatmap.tsv](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/sample_bin_heatmap.tsv)
10. 主腳本：[build_loh_round1_cross_sample_audit.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit.py)
11. v2 圖片腳本：[build_loh_round1_cross_sample_audit_v2_figures.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/build_loh_round1_cross_sample_audit_v2_figures.py)
12. HP family 對應：[LabelTest.hpp](/big7_disk/liaoyoyo2001/InterSubMod/include/core/LabelTest.hpp)、[LabelTest.cpp](/big7_disk/liaoyoyo2001/InterSubMod/src/core/LabelTest.cpp)
13. LOH evidence closeout：[20260328_LOH_evidence_panel_final_report_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md)
