<!--
建立時間: 2026-03-12 16:20
目標: 正式關閉 Pool B FN caller-side rescue 與甲基訊號整合問題，將結論接回四象限整合、annotation narrative 與主線追蹤入口
處理範圍:
  - HCC1395 5kHz tumor-only Pool B FN（ClairS non-PASS ∩ truth set）
  - 原始 caller-only rescue 規則
  - passified VCF 補跑 InterSubMod 後的 combined rescue 規則
  - Pool B 特徵分布與主線整合口徑
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/
-->

# Pool B FN integration closeout

## 破題結論

`Pool B FN` 這條支線現在可以正式關閉，結論很清楚：

1. `Pool B FN = 840` 再次證明 `ClairS-TO` 在 `HCC1395 5kHz tumor-only` 的非 PASS 區仍存在明確的 caller-side rescue 空間。
2. 在這個母體中，**最有效的救援訊號仍是 caller-first**，尤其是 `VariantCluster` 與 `GQ`。原始 caller-only 最佳規則為 `no_varcluster`（`428 TP / 390 FP`, `F1 +0.003391`）；若要最乾淨的實用規則，則是 `no_varcluster_and_gq15`（`115 TP / 45 FP`, `precision=71.9%`, `F1 +0.001452`）。
3. 補跑甲基後，最佳 combined 規則 `gq15_and_allele_delta_low` 為 `431 TP / 473 FP`, `F1 +0.002702`，只比相近 caller gate `gq15_and_af015` 多 `+0.000386`，但仍**落後**原始 caller-only 最佳 `no_varcluster`。
4. `Pool B` 中 `AlleleDelta` 的意義和 TO kept-set / artifact triage 不同。這裡 `AlleleDelta < 0.15` 比較像「較不 artifact-like」的 somatic-rich 訊號；不能把 `low VAF + high AlleleDelta` 的後段 FP 規則直接全域平移到 `Pool B`。
5. 因此 `Pool B FN` 最終再次支持主線結論：**第一層是 caller-first，甲基只提供弱補強，不適合作為獨立主 rescue engine。**

---

## 1. 研究問題

這輪 closeout 要正式回答的是：

1. `Pool B FN` 是否已足夠成熟，可接回主線報告？
2. `caller-only` 與 `with_methyl` 的規則表，哪一套才是 authoritative evidence？
3. `AlleleDelta`、`PairwiseMedianDist`、`Strong/Subclone` 在 `Pool B` 中到底代表什麼？

---

## 2. Pool B 定義與 baseline

### 2.1 Pool B 定義

`Pool B` 定義為：

- `ClairS-TO final VCF` 中 `FILTER != PASS`
- 排除 `NonSomatic`
- 再與 `truth set` 交集

在 `HCC1395 5kHz tumor-only` baseline 中：

| 指標 | 數值 | 意義 |
| --- | ---: | --- |
| `truth_total` | `39,447` | SEQC2 truth set 總數 |
| `TP` | `28,396` | ClairS-TO 已保留的真陽性 |
| `FP` | `11,843` | ClairS-TO 已保留的偽陽性 |
| `FN` | `11,051` | ClairS-TO 未保留的 truth 位點 |
| `baseline_f1` | `0.712697` | 後續所有 rescue 規則的基準 |
| `Pool B FN` | `840` | 非 PASS 但 truth set 包含的位點 |
| `Pool B non-FN` | `2,642` | 同一非 PASS 母體中，truth set 不含的位點 |

來源：

- [docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md)
- [docs/reports/validated/2026/03/assets/20260312_pool_b_fn_integration_summary_01.tsv](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260312_pool_b_fn_integration_summary_01.tsv)

### 2.2 為何 Pool B 值得獨立分析

`Pool B FN = 840`，遠大於原本用來判斷是否值得繼續分析的最低門檻 `30`。  
這代表：

- `ClairS-TO` 的非 PASS 區域不是只有零星噪音
- caller-side rescue 空間是真實存在的
- 但是否值得用甲基訊號補強，必須靠同母體比較，而不能直接套用 PASS kept-set 的 artifact 規則

---

## 3. 證據層級與方法學差異

這一段是本次 closeout 最重要的技術釐清。

### 3.1 原始 caller-only 規則表

原始 caller-only authoritative 來源是：

- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_caller_rescue_rules.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_caller_rescue_rules.tsv)

這張表直接讀原始 `ClairS-TO final VCF` 的 `FILTER` 與 `GQ/AF`，因此：

- `no_varcluster`
- `lowqual_only`
- `no_varcluster_and_gq15`

這類規則的 `FILTER` 語意是有效且 authoritative 的。

### 3.2 with_methyl / passified VCF

補跑甲基的流程使用了：

- `pool_b_fn_passified.vcf.gz`
- `pool_b_non_fn_passified.vcf.gz`

對應 log 可在：

- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_intersubmod_run.log](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_intersubmod_run.log)

這代表在 `with_methyl/` 目錄下，VCF 的 `FILTER` 已被 passify 成 `PASS`，因此：

- `with_methyl/pool_b_caller_rescue_rules.tsv`

中的 caller 規則，**不能**再與原始 caller-only 規則做橫向比較。

### 3.3 authoritative evidence hierarchy

因此本輪 closeout 固定採以下證據層級：

1. **caller-side FILTER / GQ / AF rescue**
   - 一律使用原始：
   - [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_caller_rescue_rules.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_caller_rescue_rules.tsv)
2. **甲基覆蓋率、joined features、caller+methyl combined 規則**
   - 一律使用：
   - [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_fn_stats.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_fn_stats.tsv)
   - [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_combined_rescue_rules.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_combined_rescue_rules.tsv)
   - [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_joined_features.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_joined_features.tsv)

這個層級必須在後續所有總報告中固定下來，避免再把 `passified VCF` 的 `FILTER` 結果混寫成原始 caller 規則。

---

## 4. 原始 caller-side rescue：主力仍在 VariantCluster 與 GQ

### 4.1 最佳 caller-only 規則

| 規則 | `TP+` | `FP+` | `rescue_precision` | `delta_f1` | 意義 |
| --- | ---: | ---: | ---: | ---: | --- |
| `no_varcluster` | `428` | `390` | `0.5232` | `+0.003391` | 最佳整體 F1 增量 |
| `gq15_and_af015` | `448` | `547` | `0.4503` | `+0.002316` | 單純 caller gate 的高 recall 版本 |
| `gq_ge_15` | `449` | `558` | `0.4459` | `+0.002235` | `GQ` 單獨也成立 |
| `lowqual_only` | `204` | `183` | `0.5271` | `+0.001651` | 單純 LowQual-only 仍有可救空間 |
| `no_varcluster_and_gq15` | `115` | `45` | `0.7188` | `+0.001452` | 最乾淨、最實用的規則 |
| `gq_ge_20` | `187` | `212` | `0.4687` | `+0.001119` | 更嚴格但整體增量較小 |

### 4.2 這些數據代表什麼

1. `VariantCluster` 是 Pool B 中最重要的 caller-side區分訊號。  
   `no_varcluster` 能直接帶來本輪最大的整體 F1 改善，表示很多被 `VariantCluster` 一起刷掉的 truth 位點確實值得被重新檢視。
2. `GQ` 也成立，但在 Pool B 中不是唯一主軸。  
   單用 `gq>=15` 可救 `449 TP`，但也重新引入 `558 FP`；因此 `GQ` 需要配合 caller 其他訊號一起用。
3. 若要追求實用 precision，`no_varcluster_and_gq15` 最合理。  
   它雖不是最佳 `delta_f1`，但 `precision=71.9%` 明顯高於其他規則，代表這條規則更接近可實務採用的乾淨 rescue 池。

---

## 5. 補跑甲基後的 combined rescue：只有弱補強，沒有翻盤

### 5.1 補跑甲基覆蓋率

| 指標 | 數值 |
| --- | ---: |
| `fn_with_methyl_coverage` | `0.9988` (`839/840`) |
| `non_fn_with_methyl_coverage` | `0.9875` (`2609/2642`) |

這表示：

- `with_methyl` 的母體覆蓋率足夠高
- combined rescue 可以正式下結論
- 不能再把「coverage 不夠」當作這條支線未定論的主因

### 5.2 最佳 combined 規則

| 規則 | `TP+` | `FP+` | `rescue_precision` | `delta_f1` | 解讀 |
| --- | ---: | ---: | ---: | ---: | --- |
| `gq15_and_allele_delta_low` | `431` | `473` | `0.4768` | `+0.002702` | 最佳 combined，但仍低於 `no_varcluster` |
| `gq15_and_pairwise_ge_015` | `340` | `448` | `0.4315` | `+0.001471` | 高 pairwise 有訊號，但不強 |
| `gq_ge_20` | `187` | `212` | `0.4687` | `+0.001119` | caller-only 基準線 |
| `gq15_and_not_noise` | `259` | `357` | `0.4205` | `+0.000984` | label 可提供弱 support |
| `gq15_and_strong_or_subclone` | `132` | `224` | `0.3708` | `+0.000128` | 強度不足，不適合升級 |

### 5.3 甲基唯一規則為何失敗

| 規則 | `TP+` | `FP+` | `delta_f1` | 判讀 |
| --- | ---: | ---: | ---: | --- |
| `class_strong_or_subclone` | `259` | `1181` | `-0.006265` | `Strong/Subclone` 在 non-FN 也很多 |
| `pairwise_ge_020` | `450` | `1612` | `-0.006968` | 高 pairwise 不是 Pool B 專屬 somatic 訊號 |

這些結果代表：

1. 甲基特徵在 Pool B 中**不能獨立主導 rescue**。
2. Pool B non-FN 母體含大量 germline / 結構化非 truth 訊號，因此：
   - `Strong/Subclone`
   - 高 `PairwiseMedianDist`
   - `HPMergedSig`
   都會在 non-FN 中出現得更多。

---

## 6. Pool B 的特徵分布：為何 `AlleleDelta` 方向會和 kept-set 不同

### 6.1 分布摘要

| 指標 | `pool_b_fn` (`n=840`) | `pool_b_non_fn` (`n=2642`) | 代表意義 |
| --- | ---: | ---: | --- |
| `median GQ` | `16.0` | `6.0` | FN 偏高 GQ，caller-side rescue 空間明確 |
| `median AF` | `0.375` | `0.24345` | FN 整體支持度更高 |
| `median AlleleDelta` | `0.0097` | `0.0648` | FN 反而更低 AD |
| `median PairwiseMedianDist` | `0.2107` | `0.2222` | pairwise 在 non-FN 略高 |
| `Strong/Subclone` | `259 (30.8%)` | `1181 (44.7%)` | 結構強訊號不代表是 truth |
| `HPMergedSig=true` | `185 (22.0%)` | `868 (32.9%)` | HP 顯著在 non-FN 更多 |
| `AlleleDelta < 0.15` | `814 (96.9%)` | `1934 (73.2%)` | 低 AD 較像 somatic-rich |
| `PairwiseMedianDist >= 0.15` | `640 (76.2%)` | `2228 (84.3%)` | 高 pairwise 不是 Pool B 專屬好訊號 |

### 6.2 為何這不和 TO kept-set artifact triage 矛盾

在 TO kept-set / artifact triage 支線中：

- `low VAF + high AlleleDelta` 較像 artifact / suspect FP

但在 Pool B 中：

- `AlleleDelta < 0.15` 反而更接近 truth-rich rescue 子集

這兩者**不矛盾**，原因是母體不同：

1. **TO kept-set**
   - 母體是 caller 已保留在 callset 內的位點
   - 問題是如何從已保留集合中再辨識 artifact FP
2. **Pool B**
   - 母體是 caller 已刷掉的 non-PASS 位點
   - 問題是如何從 non-PASS 中找回更像 truth 的 TP

因此 `AlleleDelta` 在不同母體中的方向不應全域化。  
這也是本次 closeout 的關鍵方法學結論之一。

---

## 7. 最終整合結論

### 7.1 已正式成立

1. `Pool B FN` 再次證明 `caller-first` 是主力 rescue engine。
2. `VariantCluster` 與 `GQ` 在 Pool B 中是比甲基更有判別力的訊號。
3. 補跑甲基後，最佳 combined 規則只提供弱補強，**沒有超過**原始 caller-only 最佳 `no_varcluster`。
4. 甲基唯一規則在 Pool B 中為負效益，不能升級成主 rescue 規則。
5. `AlleleDelta` 在 Pool B 與 kept-set / artifact triage 的意義不同，不能寫成一條跨母體全域規則。

### 7.2 對主線框架的回寫

Pool B closeout 之後，主線分工應維持：

- 第一層：`caller-first`
- 第二層：`methylation-support`
- 獨立旁路：`artifact triage`
- 補充層：`annotation / QC`

其中 Pool B 特別支持的結論是：

- caller-side rescue 空間真實存在
- 甲基不適合作為獨立主 rescue engine
- 甲基特徵的解讀必須先看母體，再談方向

---

## 8. 已回寫的文件

本次 closeout 已同步回寫到：

1. [docs/CURRENT_FOCUS.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md)
2. [docs/experiments/INDEX.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/INDEX.md)
3. [docs/plans/2026/03/20260312_未完項closing_plan_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260312_未完項closing_plan_01.md)
4. [docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md)

---

## 9. 可複查路徑

### 報告與整合入口

- [docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260308_ClairS邊緣FN探勘與甲基救援_01.md)
- [docs/reports/validated/2026/03/assets/20260312_pool_b_fn_integration_summary_01.tsv](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260312_pool_b_fn_integration_summary_01.tsv)

### 原始 caller-only authoritative outputs

- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_caller_rescue_rules.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_caller_rescue_rules.tsv)

### with_methyl authoritative outputs

- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_fn_stats.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_fn_stats.tsv)
- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_combined_rescue_rules.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_combined_rescue_rules.tsv)
- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_joined_features.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_joined_features.tsv)
- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/with_methyl/pool_b_summary.md)

### 補跑甲基流程 log

- [InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_intersubmod_run.log](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/research_rounds/20260308_clairs_borderline_fn_hcc1395_5khz/pool_b_intersubmod_run.log)
