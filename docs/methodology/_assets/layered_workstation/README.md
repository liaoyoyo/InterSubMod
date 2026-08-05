<!--
建立時間: 2026-07-07
更新時間: 2026-07-27
目標: 說明 exact-PS×HP layered reconstruction 全基因工作站的證據邊界、資料綁定與重生方式
處理範圍: InterSubMod/docs/methodology/_assets/layered_workstation/
關聯檔案:
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_gene_drug_annotation.py
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/results/validation_receipt.json
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/formal_pair_alt_ref_topology_join.tsv
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/focal_alt_ref_control_join.tsv
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/strict_hp_lane_cpp_topology_join.tsv
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py
  - InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/data/current_v5_read_af_topology/current_v5_read_af_topology.index.json
  - InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json
  - InterSubMod/docs/methodology/20260706_layered_data_model_units_proportions_spec_01.md
-->

# layered_workstation/ — exact-PS chr1–22 分層拓撲工作站

> **Primary authority（2026-07-24）**：exact nonmissing PS × primary HP1/HP2 ×
> threshold-3 strict endpoint read-linkage × bounded block。預設 builder 不再以
> `adjacent_gap<=50000` 傳遞合併。
>
> **可攜版**：7 個 dataset HTML + `index.html`；可離線開啟且零外部依賴。
> driver 會核對 cohort/pipeline/MLHP/topology/census receipt、SHA-256、逐列 join
> 與分母守恆，任何 drift 都 fail closed。

## 新版數據範圍與正確分母

- 98,955 final exact-PS topology groups。
- 85,941 mutation-bearing groups。
- 71,955 `ranked_complete` groups；只在此分母計算 determinacy/signature。
- 10,717 `ABSTAIN_RESOURCE_LIMIT`、3,224 zero denominator、45
  recurrence-required，均保留為獨立狀態。
- 680,527 棵 AF-global-best trees 已作 exact signature census：
  - `UNIQUE_TREE` 39,648（55.1011%）；
  - `TIED_SAME_TOPOLOGY` 23,858（33.1568%）；
  - `TIED_CROSS_TOPOLOGY` 8,449（11.7421%）；
  - 單一 root-preserving unlabeled exact shape 63,506（88.2579%）。
- 7 technical datasets = 6 biological samples；HCC1395 與
  HCC1395_DORADO 是同一 biological sample 的 technical comparison。
- 讀者介面顯示 `HCC1395_HKU`／`HCC1395_NYGC`；authority keys、檔名、
  receipt 與 hash 仍是 `HCC1395`／`HCC1395_DORADO`。
- 以全部 98,955 final groups 為分母，3,554 groups 有 actual sSNV locus
  落入 GENCODE v46 gene body 且該 gene 屬 COSMIC CGC v104；其中 1,252
  groups 在同一 HGNC gene 另有 DGIdb `approved=TRUE AND
  anti_neoplastic=TRUE` association。此層只是 gene-level context。
- 甲基 overlay 是從 407,738 evaluated rows 預先篩出的 7 個 formal G1
  positive pairs，不是 7 個樣本、盛行率或方法準確率。它們映射成 10 個
  exact PS×HP lanes（7 signal：638 direct reads；3 paired RR-only
  background：152 direct reads；合計 790）；
  只有 3/7 focal ALT/REF controls 可檢定，且僅 H2009
  `chr4:2,307,521` 通過 allele-axis gate（V=0.618、permutation p=0.002）。

**Claim ceiling**：read-AF 只排序 recurrence-allowed minimum mathematical
graph models；不是 caller VAF、CCF、calibrated posterior、cellular clone 或真實
ancestry。7/7 technical PASS 不等於 topology-complete。

## 檔案

| 檔 | 內容 | 大小 |
|---|---|--:|
| **`index.html`** | Exact-PS cohort command center：authority、分母漏斗、甲基 evidence overlay、7 dataset composition、癌症基因／藥物 locus、7×7 profile similarity、HCC1395 technical concordance 與收合 provenance。**先開這個。** | 約 90 KB |
| 各 dataset `.html` | 完整 compact final-group index、GRCh38 canvas、多選/搜尋、exact PS/HP boundary、locus/read matrix、甲基四層 evidence、癌症基因／藥物 context、完整 candidate union、selected tree 與 signature census。 | 約 14–182 MB |

## 新版互動與資料模型

- **GRCh38 全基因分布**：chr1–22 共用 `0–248,956,422 bp` 的實體座標尺度
  （chr1 為最大基準）。每條染色體骨架寬度為 `L_chr / L_chr1`，region
  點位 x 座標為 `midpoint / L_chr1`；不再把每條染色體各自拉滿同一寬度。
  全部 98,955 final groups 依 midpoint 落點，因此骨架長度與點位位置使用
  同一尺度。
- **圖例多選**：點 A 選 A，再點 B 為 A∪B，再點已選類別取消；沒有選擇
  時顯示全部。切換 ontology 會清空前一 ontology 的 selection。
- **著色維度**：Determinacy、拓撲形態、HP family、判定狀態、active k、
  癌症基因／藥物 annotation、甲基 evidence。七種 mode 都支援多選與
  第二次取消。
- **快速篩選**：全部、癌症基因 locus、癌症基因 ∩ approved
  抗腫瘤藥；它與目前拓撲圖例選擇採 AND，不覆蓋既有 topology
  selection。篩選後不存在的圖例選擇會自動清除，避免零結果孤兒狀態。
- **甲基快速篩選**：全部、formal G1 signal lanes、focal ALT/REF aligned、
  focal→partner relation resolved；與拓撲圖例、癌症基因／藥物和座標搜尋
  採 AND。沒有 overlay 的 region 是「未在 targeted 7-pair subset」，
  不是甲基陰性。
- **座標搜尋**：支援單點與區間 overlap，例如
  `chr10:87818272-87928739`。搜尋只回傳 final groups；source singleton
  abstain 不會假裝成 topology=0。
- **Region grain**：
  `exact_ps_hp_component_bounded_block`；同一 exact PS、同一 primary HP，
  且 threshold-3 strict endpoint read-linkage 支持。PS 是 boundary，
  不是 ancestry edge。
- **Determinacy**：`UNIQUE_TREE`、`TIED_SAME_TOPOLOGY`、
  `TIED_CROSS_TOPOLOGY` 只以 71,955 ranked units 為分母。
- **形態**：Single-only、Direct-only、Sister-only、Sister+direct 是
  rooted mutation-state graph geometry；不是 clone 類別。
- **Network**：factorized candidate sidecar 保存全部 ranked units 的
  minimum vertex sets 與 per-child legal/best parents。Region detail 的
  灰線是完整 minimum edge union、藍線是所有 AF-global-best trees 的
  edge union、綠線是一棵 deterministic selected exemplar；tied rows
  不把 global-best union 誤畫成一棵 selected tree。
- **甲基四層證據**：formal G1 partner association、focal ALT/REF
  independent joint control、exact PS×HP strict lane、current `all7_v2`
  best-tree pair relation 分開顯示。若 focal→partner 關係在所有
  global-best trees 都成立，候選 network 以橘色 halo 投影相關邊；它不改寫
  AF score、candidate incidence、read support 或 selected tree。H2009
  chr4 的 signal lane 為 `ABSTAIN_RESOURCE_LIMIT`，因此只標 locus
  evidence，不畫成已解 branch。
- **跨樣本比較**：state、determinacy、topology/coarse、active-k 與
  對稱 HP balance 分別正規化後比較；profile similarity 不是同一
  clone/tree 的證明。
- **癌症基因／藥物**：primary key 是 actual sSNV position，不是 region
  span。actual locus 必須落在同一 GENCODE gene body；CGC 與 DGIdb
  以 HGNC join。`NO_HIT_EVALUATED` 表示已評估但未命中，不表示該區沒有
  癌症重要性；drug association 也不是療效、適應症或 clinical
  actionability。
- **Provenance**：所有 `.json`/receipt 路徑只在預設收合的證據區出現，
  不插入主要數字敘事。

## Legacy canonical-v5（僅歷史重現）

舊頁面使用 `adjacent_gap<=50000` 傳遞合併，以及
`W_tree/W_primary/C/Topo` 等不同 grain 與分母。它只能由 `--legacy-v5`
重現，不能與新版 98,955/71,955 分母、88.2579% exact-shape 結果混用。

## 重生
```bash
cd InterSubMod
# 預設：重建 exact-PS authority 的 7 個樣本頁與 cohort index。
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py

# 等價的 direct primary renderer。
python3 \
  research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py \
  --verify-only

python3 \
  research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py

# 只重建 index；7 個既有頁必須已綁定同一 census receipt。
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py --index-only

# 資料與跨面板 regression。
python3 -m unittest \
  research/20260724_exact_ps_cpp_topology_signature_census/tests/test_exact_ps_gene_drug_annotation.py \
  research/20260724_exact_ps_cpp_topology_signature_census/tests/test_exact_ps_candidate_factorization.py \
  research/20260724_exact_ps_cpp_topology_signature_census/tests/test_exact_ps_cohort_similarity.py \
  research/20260724_exact_ps_cpp_topology_signature_census/tests/test_exact_ps_layered_workstation.py

# Chromium 桌機 / 手機真實互動與視覺稽核。
python3 \
  research/20260724_exact_ps_cpp_topology_signature_census/scripts/audit_exact_ps_workstation_playwright.py \
  --output-dir \
  research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v5_all7_10

# 僅需重現舊 canonical-v5 / 50-kb 頁面時才使用。
python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py --legacy-v5
```
- Primary renderer：
  `research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py`。
- Legacy renderer：`build_layered_workstation_v5.py`；只由 `--legacy-v5` 觸發。
- Freshness：sample page 必須符合 `intersubmod-authority =
  20260724-exact-ps-hp-strict-read-linkage`、同一 census receipt SHA 與
  `layered-workstation-exact-ps-v5` UI contract，並宣告
  `intersubmod-genome-scale = shared-grch38-bp-v1`；`--index-only` 不接受舊頁。
- Primary authority 是 2026-07-24 topology cohort receipt/all7 summary、
  per-sample MLHP/topology outputs，以及 signature `receipt.v2.json/summary.json`。
- 改 builder 後必做 syntax、跨面板 regression、Chromium runtime、
  七個 mode 多選/再次取消、雙層 AND filter、座標搜尋、零網路請求、
  390/1440px post-detail overflow 與視覺截圖。
- 本輪 audit receipt：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/artifacts/exact_ps_workstation_audit_v5_all7_10/receipt.json`。

## 與舊 topology_workstation/ 關係
舊資料夾是歷史 pooled／三軸展示；本資料夾是 exact-PS × primary-HP
regional mathematical topology 工作站。兩者可並存供方法史對照，但舊工作站
不可當新版 canonical evidence source。
