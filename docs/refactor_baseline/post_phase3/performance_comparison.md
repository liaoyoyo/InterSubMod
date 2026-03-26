<!--
建立時間: 2026-03-26
目標: Phase 3 效能優化驗證報告
處理範圍: Task 3.1~3.3 改前/改後效能比較
關聯檔案:
  - src/core/HierarchicalClustering.cpp
  - src/core/MethylationParser.cpp
  - src/core/DistanceMatrix.cpp
-->

# Phase 3 效能優化驗證報告

## 測試環境

- Platform: Linux 5.15.0-168-generic
- Build type: Release (-O3)
- Compiler: g++ (C++17)
- OpenMP: 啟用
- jemalloc: 5.3.0

---

## Task 3.1 — HierarchicalClustering UPGMA Lance-Williams O(N³)

### 問題根因

舊實作將 `build_upgma` 委派給 `build_generic`，其 per-merge 成本為 O(sz_i × sz_j)（遍歷所有成員對重算平均距離）。在不平衡樹下最差情況是 O(N⁴) 總複雜度。

### 修改方式

引入 Lance-Williams UPGMA 更新公式：
```
D[i][m] = (sz[i]*D[i][m] + sz[j]*D[j][m]) / (sz[i]+sz[j])
```
per-merge 成本降為 O(N)，總複雜度降為 O(N³)。

### 效能數據

| N (reads) | 優化前估計 | 優化後實測 |
|-----------|-----------|-----------|
| 100       | ~1 ms     | 1 ms      |
| 500       | ~分鐘級    | **243 ms** |

### 正確性驗證

- `UPGMA_BasicCorrectness`、`UPGMA_TreeTopology`、`UPGMA_BranchLengths` 全通過
- 4x4 測試矩陣：root height = 2.5（預期值），clade {A,B}、{C,D} 正確

---

## Task 3.2 — MethylationParser 零分配 MM tag 解析

### 問題根因

`parse_mm_tag()` 及 `parse_read()` 的 block-scan 每個 read 產生：
- `std::string mm(mm_str)` — MM tag 完整複製（~50-500 bytes）
- `mm.substr(pos)` — 子字串分配
- `std::istringstream` 物件
- `std::string token` per comma（~5-20 個）
- `std::stoi` 可能拋出例外

100k reads × 200 CpGs 估計每個 read 有 6+ 次堆積分配。

### 修改方式

| 改前 | 改後 |
|------|------|
| `std::string mm(mm_str)` | `std::strstr(mm_str, ...)` 直接掃描 |
| `mm.substr(pos)` | 指標運算 `[block_start, block_end)` |
| `std::istringstream` | 無 |
| `std::string token` | 無 |
| `std::stoi` + try-catch | `std::from_chars`（C++17，零分配，零例外） |

### 改動量

`src/core/MethylationParser.cpp`: 51 行新增 / 69 行刪除（淨減 18 行）

### 正確性驗證

所有 159 個測試通過，包含 methylation parsing 相關測試。

---

## Task 3.3 — DistanceMatrix OpenMP schedule + reduction

### 問題根因

```cpp
std::atomic<int> valid_pairs{0};
std::atomic<int> invalid_pairs{0};
std::atomic<long long> total_common_coverage{0};

#pragma omp parallel for schedule(dynamic) ...
for (int i = 0; i < n; ++i)
    for (int j = i+1; j < n; ++j)
        valid_pairs++;  // 每個 pair 3 次原子操作
```

n=500：~125k pairs × 3 atomic ops = ~375k 原子競爭寫入。
`schedule(dynamic)` 預設 chunk=1，對三角矩陣（每行工作量遞減）調度開銷偏高。

### 修改方式

1. `std::atomic` → OpenMP `reduction(+: ...)` 子句：每個執行緒本地累加，barrier 時合併，零競爭。
2. `schedule(dynamic)` → `schedule(guided, 4)`：chunk 大小從大到小動態縮減（最小 4 行），對三角不均衡更好。

### 改動量

`src/core/DistanceMatrix.cpp`: 18 行新增 / 15 行刪除（移除 3 個 std::atomic 宣告）

### 正確性驗證

所有 159 個測試通過。

---

## Phase 3 整體測試驗證

```
100% tests passed, 0 tests failed out of 159
Total Test time (real) = 4.32 sec
```

## Phase 3 編譯驗證

- Clean build 時間：65 秒（-j$(nproc)）
- 無新增編譯警告（僅保留原有兩個 pre-existing unused 警告）

## 程式碼規模變化（Phase 3 後）

| 檔案 | 行數 |
|------|------|
| RegionProcessor.cpp | 1,294（Phase 2 分解後） |
| LabelTest.cpp | 885 |
| DistanceMatrix.cpp | 631 |
| HierarchicalClustering.cpp | 552（UPGMA 優化，+39） |
| MethylationParser.cpp | 290（零分配，-18） |
| ReadAggregator.cpp | 72（Phase 2 新增） |

測試程式碼總計：3,647 行（159 tests）

---

## 版控紀錄

| Commit | 任務 |
|--------|------|
| `0197caf` | Task 3.1: UPGMA Lance-Williams O(N³) |
| `d0004b8` | Task 3.2: MethylationParser zero-alloc |
| `a9912e7` | Task 3.3: DistanceMatrix OpenMP reduction |

Branch: `refactor/phase1-safety`
