---
title: 實作藍圖 — 既有 C++ 資產盤點與 LongLineage 整合路徑
date: 2026-08-06
status: 實作依據（P1 起點）
parent: 2026-08-06-io-contract-spec.md
build_branch: research/subclonal-reconstruction-202606
---

# 實作藍圖

## 0. 關鍵發現：C++ 不需從零寫

使用者原話：「ISM layered solver 理論上是完整可執行的，但是用的不好太舊與沒有清楚整理，
所以要完整用 LongLineage 來 C++ 實作」。

盤點確認：**完整的 C++ 實作已存在**，散落在 `InterSubMod/research/` 下，共 **3,152 行**，
且已跑完 7 樣本 × chr1-22 並凍結。

| 檔案（`InterSubMod/research/...`） | 行數 | 職責 | 對應規範模組 |
|---|---|---|---|
| `20260723_production_exact_ps_strict_read_linkage/cpp/strict_endpoint_graph_verify.cpp` | 325 | **PS×HP 連接分群** | M1 核心 |
| `20260722_exact_ps_k12_hcc1395_pilot/cpp/exact_ps_k12_partition.cpp` | 934 | k≤12 切分 | M1 §1.5 |
| `20260724_exact_ps_cpp_topology_af_all_samples/cpp/exact_ps_topology_af.cpp` | 1,437 | 拓撲 + read-AF ranking | M1 §1.7 |
| `20260724_exact_ps_cpp_topology_signature_census/cpp/exact_topology_signature_census.cpp` | 456 | signature census | 統計層 |
| `20260724_exact_ps_cpp_topology_signature_census/cpp/exact_topology_candidate_factorization.cpp` | — | 候選分解 | 拓撲層 |

⇒ **P1 的真實工作 = 整合，不是重寫。**

## 1. 既有 CLI 介面

### strict_endpoint_graph_verify（M1 核心）
```
--input <molecule_calls>
--edges-output EDGES.tsv
--components-output COMPONENTS.tsv
--receipt-output RECEIPT.json
--threshold <T>
```

### Python 對照與 parity 工具（已存在）
| 腳本 | 用途 |
|---|---|
| `InterSubMod/scripts/build_strict_ps_hp_regions.py`（695 行） | Python 參考實作；`--primary-threshold` 預設 **3**，`--thresholds` 預設 `(1,2,3,5)` |
| `InterSubMod/scripts/compare_strict_graph_python_cpp.py` | **Python↔C++ parity 比對**，schema `intersubmod.strict_endpoint_python_cpp_parity` |
| `InterSubMod/scripts/summarize_strict_graph_parity.py` | parity 彙總 |
| `InterSubMod/scripts/summarize_strict_ps_hp_regions.py` | region 彙總 |
| `InterSubMod/scripts/finalize_strict_partition_run.py` | run 收尾 |
| `InterSubMod/scripts/run_layered_v4_strict.py` | 端到端編排（含 molecule calls） |

⇒ **parity 驗證機制已經內建**，P1 的驗證不需新寫。

## 2. 資料鏈的真實起點

```
tagged BAM ──(extraction，不在 InterSubMod repo 內)──▶ molecule_sparse_calls.tsv.gz
                                                              │
                                          ┌───────────────────┴────────────────┐
                                          ▼                                    ▼
                          build_strict_ps_hp_regions.py         strict_endpoint_graph_verify
                                    (Python 參考)                        (C++ 生產)
                                          └──── compare_strict_graph_python_cpp.py ────┘
```

### ✅ P0 已解（2026-08-06）

extractor 就在 repo 內：
`InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/extract_lossless_read_linkage_collapsing.py`
（**1,248 行 Python**，`run_layered_v4_strict.py:55-59` 的 `DEFAULT_EXTRACTOR`）

CLI：
```
--manifest <含 BAM/VCF 路徑>  --sample <DATASETS>  --chrom <AUTOSOMES>
--output-dir  --mapq-min 20  --baseq-min 20  --bridge-thresholds "1,2,3,5"
--samtools <path>  --samtools-threads N
```
**它透過 samtools 直接讀 BAM** ⇒ BAM → molecule_calls 這一步已存在。

### 完整既有鏈條（全部已存在，已跑過 7 樣本 × chr1-22）

```
tagged BAM + PASS VCF（manifest 綁定）
   │
   ▼  extract_lossless_read_linkage_collapsing.py   [1,248 行 Python + samtools]
molecule_sparse_calls.tsv.gz
   │
   ▼  strict_endpoint_graph_verify.cpp              [325 行 C++]
   │  （Python 對照：build_strict_ps_hp_regions.py 695 行）
edges.tsv + components.tsv + receipt.json
   │
   ▼  exact_ps_k12_partition.cpp                    [934 行 C++]
k≤12 blocks
   │
   ▼  exact_ps_topology_af.cpp                      [1,437 行 C++]
topology + read-AF ranking  →  all7_summary.json（groups_total 11,590）
```

編排器：`InterSubMod/scripts/run_layered_v4_strict.py`
（`PRIMARY_THRESHOLD = 3`、`MAX_BLOCK_SIZE = 12`）

### 🔴 這決定了「用 LongLineage C++ 實作」的真實範圍

現況是 **Python + C++ 混合**：唯一還是 Python 的環節是 **extractor（1,248 行）**，
而它正是直接讀 BAM 的效能瓶頸環節。

⇒ 「完整用 LongLineage 來 C++ 實作」的最大工作量在於 **extractor 的 C++ 化**，
而非 graph/partition/topology（那三支已是 C++）。

LongLineage 恰好已有成熟的 BAM 讀取基礎設施（`src/io/block_reader.cpp`、
`alignment.cpp`、`variant_sites.cpp`），可直接複用 —— 這是把它放進 LongLineage 的最大理由。

## 3. P1 對照基準（已凍結，實測）

來源：`.../20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/summary/all7_summary.json`
（`all_pass: true`，8 項守恆檢查全過）

HCC1395：
| 指標 | 值 |
|---|---|
| `groups_total` | **11,590** |
| `mutation_bearing_units` | 9,624 |
| `mutation_family_complete_units` | 9,357 |
| `mutation_family_abstain_units` | 267 |
| `objective_certified_units` | 11,323 |
| `no_active_alt_units` | 1,966 (17.0%) |
| `ranked_units` | 9,130 |
| `best_tree_unique_units` | 7,047 |
| `best_tree_tied_units` | 2,083 |

`active_k`：`0→1966, 1→3259, 2→4250, 3→1279, 4→405, 5→189, 6→89, 7→49, 8→26, 9→27, 10→13, 11→15, 12→23`

guards：`max_active_bits=12`, `max_family_size=100000`, `max_search_nodes=1000`

## 3b. ✅ 既有 C++ 已實測可編譯（2026-08-06）

```bash
g++ -std=c++17 -O2 -o strict_graph_verify \
  research/20260723_production_exact_ps_strict_read_linkage/cpp/strict_endpoint_graph_verify.cpp
# exit 0，無警告
```
`--help` 輸出的規則與規範 §1.3 完全吻合：
> "Only exact non-missing PS and HP1/HP2 rows enter the graph.
> Only fixed R/A endpoint pairs from a unique molecule add one unit of edge support."

輸入 TSV 欄位：`dataset, chrom, molecule_id, hp_family, phase_set, site_indices, positions1, call_codes`

## 3c. 🔑 HCC1395 strict linkage 完整 canonical 基準

來源：`research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/data/all7_report_data.json`
（`all_pass: true`，**25 項守恆檢查全過**，含 `primary_threshold_is_three` ✓、
`all_154_autosome_rows_present` ✓ = 7 樣本 × 22 染色體）

| 指標 | 值 |
|---|---|
| `candidate_loci_S` | **79,687** |
| `canonical_molecule_rows` | 3,201,802 |
| `all_components` | 39,846 |
| `active_node_memberships` | 62,651 |
| `active_unique_loci` | 36,384（45.66% of S） |
| `W_memberships` | 34,267 |
| **`HP1_W`** | **5,798** |
| **`HP2_W`** | **5,664** |
| HP1_W + HP2_W | **11,462** |
| `W_k_gt12` | 90（0.79%） |
| `W_span_gt_50kb` | 1,064（9.28%） |
| `counterfactual_W_after_50kb_cut` | 11,462（delta = 0） |

### 🔴 決定性一致性證據

`candidate_loci_S` = **79,687** 與 LongLineage parity 報告的 **79,687 sites 完全相同**。
⇒ **兩套系統吃的是同一份 sSNV 輸入**，對接在資料層有堅實基礎。

### W 與 groups_total 的關係

### ✅ 128 的差異已解（2026-08-06 實測）

**topology 階段內部完全自洽**（實測驗證）：
```
active_k 分佈總和         = 11,590 = groups_total          ✓
k=0 (1,966)               = no_active_alt_units            ✓
groups_total − k0 = 9,624 = mutation_bearing_units         ✓
objective_certified 11,323 + objective_abstain 267 = 11,590 ✓
```

**128 的來源 = k>12 的 component 被切分**：
- strict linkage 階段：`HP1_W + HP2_W = 11,462`，其中 `W_k_gt12 = 90`（0.79%）
- topology 階段：`active_k` 分佈**最大值恰為 12**（`12→23`），**無任何 k>12 的 group**
- ⇒ 90 個 k>12 的 component 經 `exact_ps_k12_partition` 切分後變成 218 個 block
  （`11,462 − 90 + 218 = 11,590`，平均 2.42 block/component）

**⚠ 2026-08-06 更正：此等式僅對 HCC1395 成立，非普遍守恆。**

逐樣本實測（W = HP1_W + HP2_W）：

| sample | W | k>12 | groups | groups − W |
|---|---:|---:|---:|---:|
| HCC1395 | 11,462 | 90 | 11,590 | +128 |
| HCC1395_DORADO | 6,828 | 34 | 6,865 | +37 |
| **COLO829** | 13,933 | 30 | 13,919 | **−14** |
| H1437 | 16,326 | 904 | 17,598 | +1,272 |
| H2009 | 24,193 | 4,206 | 36,042 | +11,849 |
| HCC1937 | 5,155 | 36 | 5,195 | +40 |
| HCC1954 | 7,724 | 23 | 7,746 | +22 |
| **TOTAL** | **85,621** | — | **98,955** | **+13,334** |

**COLO829 為負** ⇒ W→groups 不是單純的「切分只會增加」。
必有部分 W 未進入 topology 階段（原因未查明）。

**正確敘述**：全 cohort funnel 為 **85,621 W → 98,955 units**
（與既有記錄的 funnel 描述一致）。個別樣本的增減由 k>12 切分與未進入 topology 的 W 共同決定。

🔴 **P1 待解**：COLO829 的 −14 來源。在釐清前，**不可**把 W→groups 當作守恆等式使用。

### 其他已確認的守恆
`W_pct_of_all_components` = 28.7657% ⟺ 11,462 / 39,846 ✓

## 4. 整合順序（修訂版）

| 階段 | 工作 | 驗證 |
|---|---|---|
| **P0** | 確認 `molecule_sparse_calls` 的產生方式；若不可得，實作 BAM→molecule_calls extractor | 對 HCC1395 chr1 產出 molecule calls，與既有快取比對 |
| **P1** | 移植 `strict_endpoint_graph_verify.cpp` → LongLineage `src/linkage/` | `compare_strict_graph_python_cpp.py` parity 全過；HCC1395 groups ≈ 11,590 |
| **P2** | 移植 `exact_ps_k12_partition.cpp` | 既有 7 檔 parity 測試全過 |
| **P3** | 移植 `exact_ps_topology_af.cpp` + 新增 `read_lineage_assignments` writer | `active_k` 分佈與上表一致 |
| **P4** | 新建 `ll-bam-tag` | `samtools quickcheck` + read 數一致 + IGV 可見 4 tag |
| **P5** | ISM 讀 `lc/lu/lv/ls`；先確認 5 項缺陷實況再修 | `run_tests` 258/258 維持 |
| **P6** | Python spec builder + HTML | 4 種缺件情境各產一份 |

## 5. BAM aux tag 規格（已查證定案）

SAM 官方規格（`hts-specs/SAMtags.tex` 原文）：
> "tags starting with `X`, `Y`, or `Z` and tags containing lowercase letters in either position
> are reserved for local use"

56 個標準 tag：`AM AS BC BQ BZ CB CC CG CM CO CP CQ CR CS CT CY E2 FI FS FZ H0 H1 H2 HI IH LB
MC MD MI ML MM MN MQ NH NM OA OC OP OQ OX PG PQ PT PU Q2 QT QX R2 RG RX SA SM TC TS U2 UQ`
（⚠ `MC` 是標準 tag = Mate CIGAR，故 2026-05 構想提議的 `MC:Z` 撞名，不可用）

實測 HCC1395 tagged BAM（10,000 條 read，chr1:1000000-1200000）出現的 34 個 tag：
`cm zd rl sm rn ms fn ns SA PQ AS sp de PS sd qs tp ch mx s1 pi HP du s2 RG ML st ts MM MN dx sv nn NM`

**定案（三重查證皆無衝突）**：

| tag | SAM 型別 | 內容 |
|---|---|---|
| `lc` | `Z` | component_id（lineage component） |
| `lu` | `Z` | unit_id |
| `lv` | `Z` | vertex_label |
| `ls` | `A` | `U`=UNIQUE / `M`=TIE_CLASS / `P`=PARTIAL / `A`=ABSTAIN |

不變式：`lv` 存在 ⟹ `ls` 必存在；`ls != 'U'` 時 `lv` 僅為代表值。

## 6. 已否決的設計

| 設計 | 否決理由 |
|---|---|
| M1 讀 `site_reads.tsv.bgz` | 全磁碟搜尋零結果，該 artifact 從未落地；取得它需先跑受甲基 gate 影響的 dataset-gate |
| 拆 LongLineage 的 M1/M2 甲基 gate | 實測 gate 是全滅（topology_units=0）；拆它需重定義 partner 規則 + validator 2× 成本 + 撞 P3 parity |
| 在 LongLineage run 內寫 BAM | schema `const false` + 4 處 runtime 硬拒 + typed_aux 摘要導致自我不可重現 |
| 單一全域整樹 | ISM 已判 CLOSED NEGATIVE（等機率、不可枚舉） |
| `MC:Z` 當標籤 tag | `MC` 是 SAM 標準 tag（Mate CIGAR） |
