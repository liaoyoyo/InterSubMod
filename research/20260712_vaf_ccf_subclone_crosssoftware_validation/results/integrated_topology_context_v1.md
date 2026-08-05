<!--
建立時間: 2026-07-12T07:29:59+08:00
目標: 整合 7 dataset final-shape、HCC1395 pair Read/VAF/shape、gene-drug strata 與 clean-v3 lifecycle 證據
處理範圍: Task B comprehensive validation; GRCh38 chr1-22; historical layered-v2 context
服務目標: G4 / G5
關聯檔案: InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/integrated_topology_context_v1.json
-->

# 拓撲、VAF shape 與 gene/drug 證據快照

用 SCQA：同一 HCC1395 cell line 的兩個技術 dataset 確實再現部分區域與粗結構訊號；但越接近 mutation-labeled exact tree，跨來源一致性越低，而且 clean-v3 尚未完成。因此目前結論是 **PARTIAL technical reproducibility**，不是「每區域唯一真實演化樹」或 clone truth 的證明。

> **PARTIAL / historical layered-v2**：以下拓撲數字來自 `20260710_232501_layered_reconstruction_v2` 工程快照。於 2026-07-12T07:29:59+08:00，clean-v3 `20260712_layered_reconstruction_v3_raw_all_lps_pass_v2` 仍為 `RUNNING`，沒有 terminal `_SUCCESS`，不可取代 layered-v2 成為 canonical scientific result。

## 分母先固定

- `primary_regions = complete_regions + incomplete_regions`。
- `final_single_shape_regions` 是 complete regions 中，原本結構為 `Topo=1`，或 `Topo>1` 但 VAF 第一順位集合能縮成一個 ordered HP shape tuple 的區域。
- 四類 `single only / sister only / direct only / sister+direct` 的百分比，都以 `final_single_shape_regions` 為分母。
- `unresolved = unresolved VAF tie + VAF not evaluable`，其百分比以 `complete_regions` 為分母。
- HCC pair 的 Read/VAF/shape 比較另固定在 **5,720 個 exact-coordinate 且兩側 complete 的 regions**；只有明示為 evaluable-subset 的百分比才改用該層可評估數。

## 7 dataset final-shape 全表

| Dataset | Primary | Complete | Final shape | Final/complete | Single | Sister | Direct | Sister＋direct | Unresolved | Unresolved/complete |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 7,590 | 6,940 | 6,798 | 97.95% | 910 (13.39%) | 358 (5.27%) | 4,976 (73.20%) | 554 (8.15%) | 142 | 2.05% |
| HCC1395_DORADO | 7,268 | 6,750 | 6,082 | 90.10% | 1,337 (21.98%) | 154 (2.53%) | 4,395 (72.26%) | 196 (3.22%) | 668 | 9.90% |
| COLO829 | 7,659 | 6,949 | 5,600 | 80.59% | 1,259 (22.48%) | 18 (0.32%) | 4,296 (76.71%) | 27 (0.48%) | 1,349 | 19.41% |
| H1437 | 8,630 | 6,984 | 6,622 | 94.82% | 877 (13.24%) | 92 (1.39%) | 5,430 (82.00%) | 223 (3.37%) | 362 | 5.18% |
| H2009 | 9,581 | 5,882 | 5,661 | 96.24% | 406 (7.17%) | 119 (2.10%) | 4,598 (81.22%) | 538 (9.50%) | 221 | 3.76% |
| HCC1937 | 2,674 | 2,557 | 2,542 | 99.41% | 344 (13.53%) | 332 (13.06%) | 1,694 (66.64%) | 172 (6.77%) | 15 | 0.59% |
| HCC1954 | 3,975 | 3,823 | 3,734 | 97.67% | 894 (23.94%) | 386 (10.34%) | 2,154 (57.69%) | 300 (8.03%) | 89 | 2.33% |
| **Aggregate** | **47,377** | **39,885** | **37,039** | **92.86%** | **6,027 (16.27%)** | **1,459 (3.94%)** | **27,543 (74.36%)** | **2,010 (5.43%)** | **2,846** | **7.14%** |

括號內四類比例是對各列 `Final shape`，不是對 `Complete`。完整欄位（含 structural Topo=1、VAF-resolved Topo>1、tie 與 not-evaluable）在同名 JSON。

## HCC1395 pair：證據粒度愈細，一致率愈低

| Evidence layer | 分子／分母 | 固定 5,720 比例 | 可評估子集比例 | 能說什麼 |
|---|---:|---:|---:|---|
| 粗五分類一致 | 3,969 / 5,720 | 69.39% | 同左 | 粗 mutation-state 結構在兩技術來源部分再現 |
| Read strict exact + true induced | 1,599 / 5,720 | 27.95% | 1,599 / 4,038 = 39.60% | Read-compatible 投影在嚴格或真正子結構意義下相容 |
| VAF official strict exact + true induced | 1,790 / 5,720 | 31.29% | 1,790 / 3,860 = 46.37% | VAF heuristic 增加被選定的相容結果，但不是獨立真值 |
| VAF normalized sensitivity | 626 / 5,720 | 10.94% | 626 / 1,310 = 47.79% | conditional rate 類似，但 fail-closed 後覆蓋大降 |
| 未排名 exact tree-set、允許 HP swap | 2,014 / 5,720 | 35.21% | 同左 | 完整可行候選集合的技術再現性；不是 selected tree |

粗五分類 agreement 的 95% CI 是 **68.18%–70.57%**，Cohen κ=**0.497**。染色體保留的 5,000 次 permutation null mean=39.51%、q97.5=40.51%，empirical p=0.000200；因此訊號明顯高於這個 null，但 p 值不等於 topology accuracy。

Read 層的候選空間幾乎不衝突：4,038 個可評估 regions 中只有 2 個 disjoint。可是兩側都縮成 singleton projection 的只有 361 個，其中 359 相同（99.45%）；其 coverage 只有 **361/5,720=6.31%**。所以「99.45%」是很窄子集的條件結果，不能當全域準確率。VAF official 層則有 **1,234/3,860=31.97%** disjoint conflicts，顯示它提高決斷性時也增加跨來源衝突。

### VAF-selected tree 與 shape

| Endpoint | Pair n | Coverage of 5,720 | Ordered | HP-swap tolerant | Claim ceiling |
|---|---:|---:|---:|---:|---|
| VAF 唯一 mutation-labeled exact forest | 2,543 | 44.46% | 524 (20.61%) | 949 (37.32%) | 推測樹；不是 biological tree accuracy |
| Structure-first＋VAF 單一 unlabeled shape | 5,168 | 90.35% | 2,317 (44.83%) | 3,667 (70.96%) | 只比較去 mutation labels 的 branching skeleton |
| 兩側原 Topo>1、VAF 各縮至單一 shape | 2,089 | 36.52% | 919 (43.99%) | 1,624 (77.74%) | 最乾淨的 shape-rescue subset；77.74% 仍不是 exact-tree accuracy |

HCC pair 在成功選出 final shape 後，`direct only` 為 **73.20% vs 72.26%**，相差僅 0.94 pp；四類 conditional TVD=**0.0860**。但 unresolved 是 **2.05% vs 9.90%**；把 unresolved 納入 complete-region 五狀態後 TVD 上升到 **0.1455**。因此「主要 direct 結構比例相似」成立，「完整 profile 等價」不成立。

## Gene／drug strata：7/7 均未顯著富集

檢定為 5,000 次雙尾 permutation，按 chromosome＋global region-length decile 分層；annotation 是 context，不是 truth。

| Stratum | n present | Present agreement | Delta vs absent | p(two-sided) |
|---|---:|---:|---:|---:|
| GENCODE gene body | 4,459 | 69.61% | +1.02 pp | 0.4271 |
| COSMIC CGC v104 gene body | 268 | 72.39% | +3.15 pp | 0.5855 |
| DGIdb interaction gene body | 1,104 | 70.20% | +1.01 pp | 0.7600 |
| DGIdb approved-drug gene body | 830 | 70.12% | +0.86 pp | 0.8486 |
| DGIdb approved＋antineoplastic gene body | 454 | 72.25% | +3.11 pp | 0.4295 |
| COSMIC CLP HCC1395 allele, all status | 333 | 69.97% | +0.62 pp | 0.3195 |
| COSMIC CLP confirmed somatic allele | 12 | 75.00% (CI 46.77%–91.11%) | +5.62 pp | 0.7313 |

所有 p>0.05。這些欄位只能提供 face validity、基因脈絡與後續優先排序；不能證明 topology 正確、藥物有效或可用藥。

## Clean-v3 lifecycle 快照

於 2026-07-12T07:29:59+08:00：

- `run_state`: `RUNNING`, sequence=4，state timestamp=`2026-07-11T22:41:01.060128Z`。
- `heartbeat`: seq=49，wall time=`2026-07-11T23:29:04.659136Z`，active samples=`COLO829,H1437`。
- launcher PID 1904861 當時存活；terminal success marker 不存在。
- 因此本快照保留 layered-v2 數值，但明示為 historical engineering context；v3 在 verified terminal success 前不能升級結論。

`run_state.json` snapshot SHA-256=`fb15452400c5c2460933d1783882e108dde9cfe3192730f6192ddc592883da45`；heartbeat 是 mutable live file，該時點 SHA-256=`b6025c743f135fb12b8eb0ad9864b64f89de21225c7ba3033bd3bd23373d0d7b`。

## Claim ceiling

可以說：同一 HCC1395 cell line 的兩技術來源出現高於 null 的區域／粗結構訊號，Read-compatible constraints 與部分 VAF-selected shape 能在可評估子集中再現；這支持工程內部自洽與部分跨技術再現性。

不可以說：每區域已得到準確、唯一的真實演化樹；兩份資料是獨立的 biological clone truth；VAF winner 是 calibrated posterior；gene/drug overlap 驗證了演化或療效；layered-v2 已等同 clean-v3 canonical release。

核心原因是：fine-grained mutation-labeled tree identity 明顯低於 coarse／shape；可評估集合有 selection；Read 與 VAF 重用同一批 reads；目前缺少 single-cell、multi-region、lineage 或 synthetic truth tree。

## Provenance 與 QA

主要來源：

- `InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_four_class_census.json` — `4fd23e16b82543bb3a1a6e42d4a69d1a3f1ccb94d3a422a30b90548d5b4fb2dc`
- `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/topology_pair_analysis.json` — `42a6d99c3085ac5269665b7c61f5121813396ff9d2766d0bb884cc026ba26e1a`
- `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_vaf_pair_summary.json` — `eb9206b23dc6908105d9900db24a5b40e84de2e66cfea769c74607bbe434f1e9`
- `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_annotation_reproducibility.json` — `73127c5c3e703e686b042f0169ea933e75f32b3d96190b5dee28c81172f25f64`
- `InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data/hcc1395_site_topology_summary.json` — `0ac61ef1af3bc7d3f37e446a64711a5c61ba150c0220d095a1f6e43810caa7b1`
- `InterSubMod/research/20260712_hcc1395_pair_site_topology_containment_validation/data/hcc1395_topology_containment_summary.json` — `3659018c4795b9c9a7410faae60941ad3a7481492bbb72c2853af3f68a0f5684`
- `InterSubMod/docs/reports/in_progress/2026/07/20260712_HCC1395兩技術資料粗拓撲與癌症基因藥物一致性驗證_01/artifact.json` — `be2e021b5616a6e262420e218f132dbb9bfa9784bba54ac8506bc684dc6f9fe8`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260712_layered_reconstruction_v3_raw_all_lps_pass_v2/run_state.json` — point-in-time SHA 如上。

上游 QA：final-shape 37/37 PASS；coarse pair 64/64 PASS；site topology 18/18 PASS；site containment 23/23 PASS；VAF tree/shape pair 20/20 PASS。完整 source list、分母與算術 invariant 請以同名 JSON 為準。
