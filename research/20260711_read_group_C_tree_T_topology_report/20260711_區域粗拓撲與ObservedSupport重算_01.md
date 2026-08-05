<!--
建立時間: 2026-07-11T23:58:00+08:00
目標: 重算 complete region 的粗拓撲與 full-read-state endpoint sensitivity
處理範圍: chr1-22; 7 dataset rows; historical layered-v2 engineering snapshot
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_regions.tsv
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_census.json
-->

# 區域粗拓撲與 Full-read-state Endpoint Sensitivity 重算

> **TL;DR**：Complete=39,885 中，Topo=1 的 21,976 個可分為無 within-HP 關係 6,027、sister-only 1,374、direct-only 13,899、sister+direct 676；Topo>1 的 17,909 保留未定。另以 T=1/Topo=1 的 10,832 個區域做 full-read-state endpoint sensitivity：sister-only=1,342、direct-only=2,948、both=80。

> **主張邊界**：tree node 是 mutation state；HP1/HP2 是 ordered forest components；H_ 不是 directly observed clone。`single_only` 是內部欄位名，不得解讀為 single clone。

## Complete-region graph topology

| 類別 | 區域數 | 占 Complete |
|---|---:|---:|
| 無 within-HP 關係 (`single_only` internal key) | 6,027 | 15.11% |
| sister-only | 1,374 | 3.44% |
| direct-only | 13,899 | 34.85% |
| sister + direct | 676 | 1.69% |
| Topo>1 未定 | 17,909 | 44.90% |
| **合計** | **39,885** | **100.00%** |

註：6,027 中只有 837 個是單一 primary HP 的單節點 graph；其餘 5,190 個是 HP1/HP2 各有獨立節點，沒有可比的 within-HP 姐妹或直系關係。該 837 個單節點中，143 個是 full-read observed state，694 個為 H_* hidden/partial-supported state。

## 逐 dataset complete-region 粗拓撲

| Dataset | Complete | 無 HP 內關係 | Sister only | Direct only | Both | Topo>1 未定 |
|---|---:|---:|---:|---:|---:|---:|
| HCC1395 | 6,940 | 910 | 334 | 2,053 | 141 | 3,502 |
| HCC1395_DORADO | 6,750 | 1,337 | 140 | 933 | 34 | 4,306 |
| COLO829 | 6,949 | 1,259 | 17 | 2,592 | 8 | 3,073 |
| H1437 | 6,984 | 877 | 91 | 3,520 | 96 | 2,400 |
| H2009 | 5,882 | 406 | 112 | 2,745 | 265 | 2,354 |
| HCC1937 | 2,557 | 344 | 318 | 1,209 | 70 | 616 |
| HCC1954 | 3,823 | 894 | 362 | 847 | 62 | 1,658 |

## T=1/Topo=1 full-read-state endpoint sensitivity

| 類別 | 區域數 | 占 T=1/Topo=1 |
|---|---:|---:|
| no_observed_within_hp_relation | 6,462 | 59.66% |
| observed_sister_only | 1,342 | 12.39% |
| observed_direct_only | 2,948 | 27.22% |
| observed_sister_and_direct | 80 | 0.74% |
| **合計** | **10,832** | **100.00%** |

## 驗證

- Checks：65/65 PASS。
- 逐 dataset graph 五類均守恆回 Complete。
- Full-read-state endpoint 四類均守恆回 T=1/Topo=1；H_* 可連接路徑，但不當成直接觀測 endpoint。
- 本報告是只讀重算，不修改 layered output。
