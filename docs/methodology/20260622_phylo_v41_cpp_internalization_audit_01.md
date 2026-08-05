<!--
建立時間: 2026-06-22
問題類型: C++ 方法（切群標籤內建）
影響 track: TO（tumor-only subclone 切群層）
狀態: pending_decision → 用戶已選 B（執行中）
build_branch: feat/summary-nreadsvalid
-->

# phylo-v4.1 切群標籤 C++ 內建 — 審查文件

## 問題描述

C++ 現有切群輸出靠 `TreeCutter::find_optimal_clusters`（silhouette best_k），`src/core/RegionProcessor.cpp:2548`（及 2196/2343）。silhouette **強制 k≥2、無「無結構」選項、+ PERMANOVA double-dip** → 量化證實 **純噪音 FP 74.7-90.7%**（`phylo_methods_compare.py`），30/30 位點全判顯著=系統性過度宣稱。

Python 已定稿 **phylo-v4.1**（descend-quarantine + per-CpG null + 遞迴分裂 + modal instability + coarse/fine + other + hidden-het + 對齊），全基因組驗證 PASS（structure 23.83%、噪音 FP≈0%）。需內建進 C++ 讓 binary 一次輸出標籤、Python 只讀畫圖。

## 量化影響（全 L1，本 session 驗證）

| 指標 | silhouette(現) | phylo-v4.1(目標) |
|---|---|---|
| 純噪音 FP | **74.7-90.7%** | **≈0%**（null95 閘）|
| 「無結構」verdict | ❌ 強制 ≥2 | ✅ modal |
| 全基因組 structure% | 30/30 全顯著 | TP 23.83%（aligned 19.59% 偏 TP 1.29× / unaligned 4.24% 偏 FP 0.51×）|
| 輸出 | cluster label(flat) | coarse/fine/other + modal_frac + 對齊 |

## 整合點（IDENTIFY）

- **取代點**：`RegionProcessor::perform_clustering_and_significance`（`src/core/RegionProcessor.cpp:2451-2623`）— line 2517 建 UPGMA 樹、2548 silhouette、2623 `work_dist`+`work_meth` 皆在此 scope。within-hp(2167)/tumor-only(2308) 另有 2 個 find_optimal_clusters 呼叫。
- **新模組**：`include/core/PhyloLabeler.hpp` + `src/core/PhyloLabeler.cpp`（CMakeLists 87-96 加一行）。
- **輸出**：`RegionWriter::write_region`（`src/io/RegionWriter.cpp:83`）加 `write_phylo_groups()` → `phylo_groups.tsv` + `phylo_groups_summary.json`。
- **既有可重用**：`HierarchicalClustering`(UPGMA)/`Tree`+`MergeRecord`/`StructureTest::filter_reads_for_complete_matrix`(=peel)+`run_permanova`(置換範本)/`DistanceMatrix`(BERNOULLI)。
- **測試**：`tests/test_phylo_labeler.cpp`（參考 test_hierarchical_clustering.cpp）。

## 修改選項

### 方案 A：不改 C++（Python 後處理讀矩陣，現狀）
維持 binary 只輸出矩陣，Python 讀後做 v4.1。**不達成「Python 只畫圖」目標**，每次重畫需 Python 重算 null（慢、非 binary 真值）。

### 方案 B：新 PhyloLabeler C++ 模組（native v4.1 + 統計等價驗收）⭐ 推薦
v4.1 完整邏輯內建 C++，用 C++ `std::mt19937_64` + **pin compile-time seed**（Hard-Gate 決定論）。取代 silhouette。輸出 phylo_groups.tsv。**驗收=統計等價**（非逐位元）：C++ vs Python 在 30 pilot + 100 WG loci 的 coarse_ng 一致率 >90%，分歧只在 unstable 位點。**分 3 phase**：B1 核心(peel+descend+per-CpG null+遞迴+coarse)；B2(modal K=10 + fine + other + hidden-het + 對齊)；B3(RNG harness + 跨樣本回歸)。

### 方案 C：B + 在 C++ 重實作 PCG64（逐位元等同 Python）
B 加上把 numpy PCG64 permutation 在 C++ 複刻，求 C++ 與 Python 逐位元相同。**成本高、效益邊際**（統計等價已足，bit-identical 非必要）。

## SWOT

**方案 B**（推薦）：

| | Helpful | Harmful |
|---|---|---|
| **Internal** | **S** 重用 70% 既有零件(UPGMA/Tree/peel/DistanceMatrix)；一次 binary run 出全標籤；Python 只讀畫(達成目標) | **W** v4.1 邏輯複雜(~200 行 Python→C++ 遞迴+per-subgroup 重分群 null)；新模組 coupling RegionProcessor |
| **External** | **O** 取代證實薄弱的 silhouette；對外可稱「binary 原生 null-validated 切群」；下游 renderer 大簡化 | **T** C++ RNG≠Python PCG64 → null 門檻不逐位元同(靠統計等價驗收化解)；K=10 modal 全基因組 compute(~5× null) |

**方案 C**：

| | Helpful | Harmful |
|---|---|---|
| **Internal** | **S** 逐位元可重現、Python/C++ 完全對拍 | **W** 複刻 PCG64 高工、易錯、維護負擔 |
| **External** | **O** 最強 reproducibility 宣稱 | **T** numpy 升版改 PCG64 內部即破；過度工程 |

**方案 A**（對照）：S=C++ 不動穩定 / W=不解「Python 重算」痛點、不達目標 / O=資源投他處 / T=silhouette 薄弱持續對外。

> W+T(B)=4 vs S+O(B)=4 → 平手，B 可推進（非降級）。C 的 W+T>S+O → 條件性（不選）。

## 驗收標準

- [ ] B1-B3 各 phase `make -j` 編譯通過（Hard-Gate）
- [ ] `test_phylo_labeler` 單元測試通過（含 peel/descend/null-gate/modal 邊界）
- [ ] **RNG harness**：C++ vs Python v4.1 在 30 pilot + 100 WG loci，coarse_ng 一致率 **>90%**，分歧僅 unstable 位點
- [ ] 全基因組 C++ 跑出 phylo_groups.tsv，structure% 與 Python v4.1（23.83%）一致 ±2pt
- [ ] 跨 ≥3 樣本(HCC1395/COLO829/+1)無編譯/執行回退
- [ ] 既有 PERMANOVA/輸出不變（regression）

## 用戶決策

**選擇**：[x] **B**（native PhyloLabeler，統計等價驗收，分 3 phase）
**日期**：2026-06-22
**理由**：達成「binary 原生切群、Python 只畫圖」目標；取代證實 75-91% FP 的 silhouette；統計等價驗收化解 RNG 跨語言。C（bit-identical）過度工程不選；A 不達目標。
