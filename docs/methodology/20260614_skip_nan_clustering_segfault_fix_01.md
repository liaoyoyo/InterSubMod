---
title: "SKIP NaN 距離矩陣 clustering segfault — 根因與修復（A+C）"
date: 2026-06-14
type: methodology / bugfix
status: validated
branch: fix/clustering-nan-skip-segfault
commit: 517ed90
data_sources:
  - ctest run_tests 221/221（本輪）
  - scripts/regression/regression_check.sh PASS 2624（本輪）
  - 端到端 SKIP /tmp/skip_fixed_e2e/significance_summary.csv（本輪）
---

# SKIP NaN 距離矩陣 clustering segfault — 根因與修復

## L0 結論
`--nan-distance-strategy SKIP` 產生的 NaN 距離矩陣被直接餵進 hierarchical clustering，
clustering 在 NaN 退化下建出**有環 tree**，下游 `Tree::collect_internal_nodes` 沿環無限遞迴
→ **stack overflow → SIGSEGV (exit 139)，崩在第一個 region**。修復 = **A（治本：clustering 前
filter NaN reads，對齊下游 PERMANOVA）+ C（安全網：collect 加 cycle guard）**。三層驗證全綠。

## L1 根因鏈（gdb 重現確認）
gdb 重現（chr1 ±5000 BERNOULLI `--nan-distance-strategy SKIP -j1 OMP_NUM_THREADS=1`）：
崩在 `Progress 1/2624` 後，backtrace `collect_internal_nodes` 自我遞迴 #0–#7+，
**PC 全相同 `0x...00cf`（一直走 `node->left`）= stack overflow**。

| # | 位置 | 行為 |
|---|------|------|
| 1 | `DistanceMatrix.cpp:413,436` | SKIP → 共同 CpG < C_min(3) 的 pair `dist = NaN` |
| 2 | `RegionProcessor.cpp:1497`（修前）| `build_tree(all_dist)` 直接餵含 NaN 完整矩陣 |
| 3 | `HierarchicalClustering build_upgma:100` | `D(i,j) < min_dist` 遇 NaN 恆 false → NaN pair 永不被選；退化建出含環 tree |
| 4 | `RegionProcessor.cpp:1527` | `find_optimal_clusters → cut_by_num_clusters → get_internal_nodes` |
| 5 | `TreeStructure.cpp:73` | `collect_internal_nodes` 沿環無限遞迴 → 💥 SIGSEGV |

**為何只 SKIP 崩**：MAX_DIST 把 invalid pair 設 `1.0`（有限值）→ UPGMA 正常 merge 全 N read →
合法無環 tree（regression 2624 PASS 即證）。SKIP 設 `NaN` → clustering 層**無 NaN 處理能力**。

**防護不一致**：下游 PERMANOVA 本有完整 NaN 防護（`run_permanova` 遇 NaN 即 `valid=false`，
`filter_reads_for_complete_matrix` 專門剔 NaN-heavy read），但 **clustering 跑在這些防護之前就崩**。

## L2 誠實標註：環的精確成因
gdb 證實 `collect_internal_nodes` 無限遞迴（PC 固定走 left）；一個 region reads 中位數 118（< 1000），
正常 caterpillar tree 遞迴 < 1000 層**不可能**爆 8MB stack（~13 萬層才爆）→ **排除深度，必然是環**。
但「哪一次 merge 形成環」的精確序列**未經 Debug-build instrument 最終確認**（L2 強推理）。
**修復對此 robust**：A 治本（NaN 不進 clustering → tree 永遠合法）+ C 安全網（即使建出環也不崩），
且端到端用 gdb 重現的**真實崩潰 input** 驗證通過，不依賴推測環機制。

## 修復方案

### A — clustering 前 filter（治本，`RegionProcessor`）
新增 `extract_complete_submatrix_indices()`：greedy 剔除「非-NaN 鄰居最少」的 read 直到子矩陣無 NaN
（與 `StructureTest::filter_reads_for_complete_matrix` **同演算法** → clustering 與 PERMANOVA 用
**同一個 valid read subset**）。整個 clustering+significance 在 valid-read 空間跑（dist/names/labels/
full_labels/meth_raw 全壓縮到 valid 空間，row 對齊）。
- **MAX_DIST（無 NaN）→ valid == 全部 reads → zero-copy → 與現狀 bit-identical**（regression 安全核心保證）。
- SKIP 語意正確：共同 CpG 不足的邊緣 read 無法可靠分群，剔除（正是 SKIP 設計意圖）。

### C — cycle guard（安全網，`TreeStructure`）
`collect_internal_nodes` / `collect_leaves` 加 `std::set<const TreeNode*> visited`，
節點第二次被訪問即 return。遞迴深度 ≤ 不同節點數 ≤ 2N−1，**永不爆 stack**，且不誤傷正常深 tree。

## 驗證（三層，本輪 fresh）
| 層 | 指令 | 結果 |
|----|------|------|
| 單元（C，TDD）| `run_tests --gtest_filter=TreeCycleGuard.*` | 紅 exit 139 → 綠 **3 PASS**（先寫測試見紅再修） |
| 全套回歸 | `ctest` | **221/221 pass, 0 failed** |
| 結果不變（A 安全）| `regression_check.sh` | **PASS — MAX_DIST chr1 ±5000 BERNOULLI 2624 位點與 golden bit-identical** |
| 端到端（A+C 真實修復）| chr1 SKIP -j16（gdb 同 input）| **exit 139 → exit 0**，2624-region significance（wall 51s，invalid pairs 81% valid ratio） |

## 影響 / 解鎖
- **SKIP 路徑現可用** → 解鎖 memory `project_ism_sampling_layer_review` 的 **U1（MAX_DIST vs SKIP 對
  subclone 分群）** — 最關鍵待測點，先前因此 segfault 無法完成。
- **MAX_DIST 預設路徑零影響**（regression bit-identical）。
- 無 F1 影響（crash fix + 新解鎖路徑，非演算法調參）。

## 後續可選
- U1 對比實驗（MAX_DIST vs SKIP clustering 品質：Strong/Noise/保留率）現可跑。
- 若要 100% 確認環機制 → Debug build + ASan instrument `build_upgma`（非阻塞，修復已 robust）。
