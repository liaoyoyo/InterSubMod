<!--
建立時間: 2026-07-19T09:08:18+08:00
目標: 固定全 sSNV 與 positional-singleton focal-ALT 甲基多群的證據分層、分母、比例及 claim ceiling
處理範圍: 7 datasets / 6 biological samples / chr1-22 / 469,849 LongPhase-S recalibrated FILTER=PASS biallelic sSNVs
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v3.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/reviews/20260718_positional_singleton數值與claim獨立agent稽核_01.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_summary.json
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_audit_summary.json
-->

# Latent molecular substructure 證據分層與比例

> **PRE-RELEASE / Task Type B。** 本文件中的 M1/M2 與 positional-singleton
> 數值已由凍結結果獨立重算；正式 sSNV 共現、matched-normal、CN/CCF 與最終
> result/report signature 尚在執行。不得把本文件單獨引用為 cellular subclone
> 或 linear ancestry 的確認證據。

## 先講結論

1. 全部 `469,849` 個最新 LongPhase-S recalibrated `FILTER=PASS` autosomal
   biallelic sSNV dataset-sites 中，`459,928`（`97.8885%`）可完成 focal-ALT
   甲基評估。
2. `102,842` 個位點通過 M1 operational stable multigroup screen，占全部
   `21.8883%`、占可評估位點 `22.3605%`。這是候選產生率，不是 subclone 盛行率。
3. 全體 M1 中，strict M2 只有 `1,867` 個可評估，`919` 個 PASS。`919/1,867`
   的 `49.2234%` 是條件比例；對全部 M1 為 `0.8936%`，對全部 sSNV 為
   `0.1956%`。
4. positional-singleton 分支有 `50,432` 個位點，其中 `5,961` 個 M1 FLAGGED，
   strict M2 為 `30 PASS + 18 FAIL`；另有 `5,913 NOT_EVALUABLE` 與
   `44,471 NOT_RUN`。
5. 可支持的最高敘述是「共同 ALT reads 內有 residual read-level epigenetic
   partition，可作 latent molecular substructure / subclone candidate」。
   目前 confirmed cellular subclone=`0`，confirmed linear ancestry=`0`。

## 三種 singleton 口徑不得混用

| 口徑 | 定義 | 現在狀態 |
|---|---|---|
| Positional singleton | 同 dataset/chrom，以相鄰距離 `<=50,000 bp` 建 transitive component 後 `component_size=1` | 已完整重算，`50,432` |
| Local partner unavailable | focal window 內沒有符合 frozen PASS universe 的 partner | 由正式 cooccurrence 結果判定中 |
| Read-sharing degree-zero | 實際沒有任何 read 同時覆蓋 focal 與另一 sSNV | 不能用 positional singleton 直接替代，待正式共現表重算 |

因此，本文件的 `50,432` 僅稱 **positional singleton**。它是「局部難以利用
sSNV 共現」的保守操作型分支，不等於所有可能定義下的 read-sharing
degree-zero。

## 全體 sSNV 證據階梯

| 層級 | 條件 | 分子 / 分母 | 比例 | 可支持敘述 | 證據狀態 |
|---|---|---:|---:|---|---|
| L0 分析母體 | 最新 LongPhase-S PASS autosomal biallelic sSNV | 469,849 / 469,849 | 100% | 正式分析 universe | 已確認 |
| L1 可評估 | ALT reads 與 distance 資訊足夠 | 459,928 / 469,849 | 97.8885% | 可完成 focal-ALT methyl screen | 已確認 |
| L2 M1 | `coarse_ng>=2`、非 unstable、modal ARI>=0.8 | 102,842 / 459,928 | 22.3605% | operational stable multigroup | 已確認 |
| L2b residual screen tier | M1 後未由 screen 中 measured axes 解釋 | 88,004 / 102,842 | 85.5720% | screen-level residual signal | 已確認；非 strict M2 |
| L2c phase-anchored legacy tier | screen 的 phase-anchored robust tier | 17,870 / 102,842 | 17.3762% | 較高門檻 epigenetic candidate | 已確認；非 strict M2 |
| L3 strict M2 evaluable | 八個 measured axes 可形成有效判定 | 1,867 / 102,842 | 1.8154% | 可計算 strict residual gate | 已確認 |
| L3 strict M2 PASS | strict M2 eligible | 919 / 1,867 | 49.2234% | measured-axis residual partition | 已確認；高度選擇子集 |
| L4 genetic association | methyl group 與同 read 其他 sSNV R/A 共分離，global BY | 待正式結果 | 待計算 | local methyl-genetic co-segregation | 執行中 |
| L5 replication/control | multi-seed、REF/normal、CN/CCF、platform consistency | 待正式結果 | 待計算 | 較可信 subclone candidate | 執行中 |
| L6 cellular clone/lineage | 多 marker、CCF/CN、跨樣本或正交 truth 支持 clone 與順序 | 0 confirmed | 0% confirmed | cellular subclone / linear ancestry | 尚不支持 |

額外母數：

- `102,842 / 469,849 = 21.8883%`：M1 對全部 sSNV 的 yield。
- `88,004 / 469,849 = 18.7303%`：legacy residual screen tier 對全部 sSNV。
- `17,870 / 469,849 = 3.8033%`：legacy phase-anchored tier 對全部 sSNV。
- `919 / 102,842 = 0.8936%`：strict M2 PASS 對全部 M1。
- `919 / 469,849 = 0.1956%`：strict M2 PASS 對全部 sSNV。
- `100,974 / 102,842` 因 axis indeterminate 無法通過 strict M2 evaluability；
  另有 1 個 group count >10。這是 M2 選擇性很強的主要原因。

## 全體資料集分層

| Dataset | 全部 | 可評估 | M1 stable | M1 / 可評估 | residual screen tier | phase-anchored legacy tier |
|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | 79,687 | 78,629 | 12,838 | 16.3273% | 10,989 | 1,530 |
| HCC1395_DORADO | 79,739 | 78,637 | 14,789 | 18.8067% | 12,893 | 1,700 |
| COLO829 | 37,788 | 35,484 | 3,579 | 10.0862% | 2,742 | 608 |
| H1437 | 77,080 | 74,961 | 10,187 | 13.5897% | 8,525 | 3,560 |
| H2009 | 154,465 | 152,594 | 54,644 | 35.8101% | 47,593 | 9,153 |
| HCC1937 | 18,690 | 17,886 | 1,938 | 10.8353% | 1,575 | 341 |
| HCC1954 | 22,400 | 21,737 | 4,867 | 22.3904% | 3,687 | 978 |
| **合計** | **469,849** | **459,928** | **102,842** | **22.3605%** | **88,004** | **17,870** |

HCC1395 與 HCC1395_DORADO 是同一 biological sample 的技術資料集，不能算成
兩個獨立生物重現。

## Positional-singleton 分支

### 守恆與比例

```text
50,432 singleton
= 44,471 M1 NOT_FLAGGED
+ 5,913 M1 FLAGGED but M2 NOT_EVALUABLE
+ 18 M2 FAIL
+ 30 M2 PASS
```

| 指標 | 分子 / 分母 | 比例 | 正確解讀 |
|---|---:|---:|---|
| Singleton / 全部 sSNV | 50,432 / 469,849 | 10.7337% | positional 定義下的單點區域 |
| M1-evaluable / singleton | 48,347 / 50,432 | 95.8657% | 技術可評估率 |
| M1 FLAGGED / singleton | 5,961 / 50,432 | 11.8199% | 對全部 singleton 的 operational yield |
| M1 FLAGGED / M1-evaluable | 5,961 / 48,347 | 12.3296% | 對可評估 singleton 的 conditional yield |
| M2-evaluable / singleton M1 | 48 / 5,961 | 0.8052% | strict M2 具有有效判定力的比例 |
| M2 PASS / singleton M1 | 30 / 5,961 | 0.5033% | 對 M1 分母的保守 strict yield |
| M2 PASS / 全 singleton | 30 / 50,432 | 0.0595% | 對全部 singleton 的保守 strict yield |
| M2 PASS / M2-evaluable | 30 / 48 | 62.5% | 僅描述高度選擇的 48 個，不得外推 |

### 各資料集

| Dataset | Singleton | M1-evaluable | M1 FLAGGED | FLAGGED / evaluable | M2 PASS | PASS / singleton |
|---|---:|---:|---:|---:|---:|---:|
| COLO829 | 7,830 | 7,347 | 617 | 8.3980% | 0 | 0% |
| H1437 | 6,696 | 6,298 | 813 | 12.9089% | 2 | 0.0299% |
| H2009 | 2,853 | 2,782 | 602 | 21.6391% | 9 | 0.3155% |
| HCC1395 | 8,279 | 8,074 | 734 | 9.0909% | 2 | 0.0242% |
| HCC1395_DORADO | 8,321 | 8,016 | 962 | 12.0010% | 2 | 0.0240% |
| HCC1937 | 8,469 | 8,069 | 843 | 10.4474% | 14 | 0.1653% |
| HCC1954 | 7,984 | 7,761 | 1,390 | 17.9101% | 1 | 0.0125% |
| **合計** | **50,432** | **48,347** | **5,961** | **12.3296%** | **30** | **0.0595%** |

### Truth label 分層

| Truth | Singleton | M1-evaluable | M1 FLAGGED | FLAGGED / evaluable | M2 P / F / NE / NR |
|---|---:|---:|---:|---:|---:|
| TP | 45,193 | 44,171 | 5,494 | 12.4380% | 30 / 16 / 5,448 / 39,699 |
| FP | 1,084 | 514 | 52 | 10.1167% | 0 / 1 / 51 / 1,032 |
| UNASSESSED | 4,155 | 3,662 | 415 | 11.3326% | 0 / 1 / 414 / 3,740 |

全部 30 個 M2 PASS 均為 benchmark TP，但 FP 的 M2-evaluable 只有 `n=1`，
且 FP 的 M1 可評估率只有 `47.42%`。因此這個結果不能估計 specificity，也不能
解讀為 M2 已能區分 TP/FP。

## 技術重現限制

HCC1395 與 HCC1395_DORADO 的 positional-singleton loci Jaccard=`82.10%`，
M1 multigroup loci Jaccard=`31.57%`；兩邊各 2 個 M2 PASS，但 M2 PASS
intersection=`0`、union=`4`。目前沒有 locus-level strict M2 technical
replication，這是把候選提升為 clone claim 前必須保留的負向證據。

## 生物學解釋邊界

與目前結果相容的模型包括：

1. 共同 ancestral ALT 被多個後代 clone 攜帶，甲基狀態可協助拆分後代。
2. 同一 clone 內存在可逆或狀態依賴的 epigenetic heterogeneity。
3. cis allele-specific methylation、局部 CN/purity、read geometry 或其他未測
   technical/biological factors造成分群。

只觀察單一 focal ALT 時，上述模型可能產生相同資料圖樣。因此 M2 PASS 是
**subclone candidate generator**，不是 clone identifier。提升到 L4/L5 至少需要：

- 同一 read 上其他 sSNV 與 methyl group 的 global-BY association；
- RR/AR/RA/AA 完整 focal ancestry，而非只看 ALT reads；
- multi-seed、REF、matched-normal、HP/PS、strand/read geometry 控制；
- CN/CCF 與 HCC1395/DORADO 或其他正交重現；
- 對 linear 與 branching 模型做可識別性比較。

## 待正式鏈補齊

- `102,842` 位點 raw-read identity、HP/PS 與 typed aux-tag preflight。
- 全量 partner/cooccurrence、999 permutations、global BY 與 topology。
- 真正 local-partner-unavailable 與 read-sharing degree-zero 分母。
- tumor-REF、matched-normal、multi-seed、CN/CCF 與 final candidate tiers。
- 最終 result/report/supplemental signatures、HTML QA 與 evidence ledger。
