<!--
建立時間: 2026-03-26
目標: C++ 程式碼重構的量化基準，所有重構前後的比較基準點
處理範圍: 全量 C++ 核心程式碼（src/core/, include/core/, tests/）
關聯檔案:
  - docs/refactor_baseline/test_results_baseline.txt
  - docs/refactor_baseline/code_metrics_baseline.txt
  - docs/refactor_baseline/memory_baseline.txt
-->

# InterSubMod C++ 重構 Baseline 基準

**建立時間：** 2026-03-26
**Git Branch：** `refactor/baseline-checkpoint`
**基準目的：** 所有重構 Phase 完成後，必須與此 baseline 比對，確認功能等效、無效能退化、無記憶體增加。

---

## 1. 測試套件基準

| 指標 | 數值 |
|------|------|
| 總測試數 | **107** |
| 通過 | **107** |
| 失敗 | **0** |
| 執行時間 | **1.27 秒**（wall clock）|
| Peak 記憶體 | **187 MB**（192,020 KB）|

**測試分佈（17 個 test suite）：**
- DistanceMatrixTest
- FisherExactTest
- MathUtilsTest
- GlobalTestTest
- LocalTestTest
- StructureTestTest
- BootstrapTest
- HierarchicalClusteringTest
- BamReaderTest
- ConfigTest
- SNVLoadingTest
- SignificanceAnalyzerTest
- *(詳見 test_results_baseline.txt)*

---

## 2. 程式碼規模基準

### src/core/*.cpp（共 6,693 行）

| 檔案 | 行數 | 重構優先級 |
|------|------|----------|
| RegionProcessor.cpp | **1,358** | 🔴 Phase 2 核心 |
| LabelTest.cpp | **885** | 🔴 Phase 2 |
| DistanceMatrix.cpp | **628** | 🟠 Phase 2 |
| HierarchicalClustering.cpp | **481** | 🟡 Phase 3 |
| FisherExact.cpp | **415** | 🔴 Phase 1（數值） |
| SignificanceAnalyzer.cpp | **410** | 🟠 Phase 2 |
| MethylationParser.cpp | **308** | 🟡 Phase 3 |
| StructureTest.cpp | **280** | ⏸ 暫緩 |
| LocalTest.cpp | **271** | 🟡 Phase 2 |
| MathUtils.cpp | **239** | 🟢 工具類 |
| GlobalTest.cpp | **239** | 🟢 穩定 |
| ReadParser.cpp | **234** | 🔴 Phase 1（型別）|
| Bootstrap.cpp | **196** | 🟢 穩定 |
| TreeStructure.cpp | **191** | 🟢 穩定 |
| SomaticSnv.cpp | **188** | 🟢 穩定 |
| Config.cpp | **138** | 🔴 Phase 1（洩漏）|
| BamReader.cpp | **117** | 🔴 Phase 1（RAII）|
| MatrixBuilder.cpp | **115** | 🟢 穩定 |

### include/core/*.hpp（共 3,560 行）

最大：Stats.hpp (557)、RegionProcessor.hpp (391)、LabelTest.hpp (245)

### tests/*.cpp（共 2,762 行）

| 測試檔案 | 行數 | 涵蓋模組 |
|---------|------|---------|
| test_hierarchical_clustering.cpp | 536 | HierarchicalClustering |
| test_distance_matrix.cpp | 464 | DistanceMatrix（6 種度量）|
| test_significance_stats.cpp | 400 | MathUtils、FisherExact |
| test_global_local.cpp | 366 | GlobalTest、LocalTest |
| test_structure_bootstrap.cpp | 348 | StructureTest、Bootstrap |
| test_significance_analyzer.cpp | 325 | SignificanceAnalyzer（集成）|
| test_bam_reader.cpp | 147 | BamReader |
| test_config.cpp | 89 | Config |
| test_snv_loading.cpp | 58 | SomaticSnv |
| main_test.cpp | 29 | 測試入口 |
| **LabelTest** | **0** | ⚠️ 無測試！Phase 2 補充 |

---

## 3. 記憶體使用基準

| 場景 | Peak RSS |
|------|---------|
| 測試套件全跑（107 tests） | **187 MB** |
| 單一測試（SignificanceAnalyzer） | ~0.37 MB |

---

## 4. 重構前已知的主要問題

### Phase 1（安全性，立即執行）

| # | 問題 | 位置 | 類型 |
|---|------|------|------|
| 1.1 | `fetch_reads()` 返回原始 `bam1_t*`，例外路徑洩漏 | BamReader.cpp | 記憶體洩漏 |
| 1.2 | `validate()` 中 idx 在 hdr_read 失敗時未釋放 | Config.cpp | 記憶體洩漏 |
| 1.3 | `exp(log_p_k)` 直接累加，極端 table 精度損失 | FisherExact.cpp | 數值穩定性 |
| 1.4 | ~~ss_between = ss_total - ss_within 消去風險~~ | StructureTest.cpp | **⏸ 暫緩** |
| 1.5 | `uint32_t → int32_t` 轉換無溢位保護 | ReadParser.cpp | 型別安全 |

### Phase 2（架構，Phase 1 後執行）

| # | 問題 | 位置 | 影響 |
|---|------|------|------|
| 2.1 | 6 個距離函數 ~200 行重複代碼 | DistanceMatrix.cpp | DRY |
| 2.2 | LabelTest **0% 測試覆蓋** | LabelTest.cpp | 測試盲點 |
| 2.3 | RegionProcessor 1358 行，15+ 職責 | RegionProcessor.cpp | God Class |
| 2.4 | LabelTest 4 個統計測試混在一類 | LabelTest.cpp | 單一職責 |
| 2.5 | API 命名不一致（run/test/compute） | 多處 | 可讀性 |

### Phase 3（效能，Phase 2 後執行）

| # | 問題 | 位置 | 影響 |
|---|------|------|------|
| 3.1 | UPGMA O(N³)，1000+ reads 很慢 | HierarchicalClustering.cpp | 效能 |
| 3.2 | `istringstream` + `stoi` 字符串解析 | MethylationParser.cpp | 效能 |
| 3.3 | OpenMP `schedule(dynamic)` 未優化 | DistanceMatrix.cpp | 效能 |

---

## 5. 重構通過條件（每 Phase 完成後核對）

- [ ] `ctest` 通過率 ≥ 100%（不得低於 107/107）
- [ ] `chr19-verification` 輸出 diff 為空（功能等效）
- [ ] Peak RSS ≤ 200 MB（不得超過 baseline 的 107%）
- [ ] Valgrind 無新增 `definitely lost` 記憶體
- [ ] 各 Task 均有對應 git commit 紀錄

---

## 6. 重構 Git 分支策略

```
test/permanova-enable（目前主線）
 └─ refactor/baseline-checkpoint  ← 此分支（建立基準）
     ├─ refactor/phase1-safety
     ├─ refactor/phase2-architecture
     └─ refactor/phase3-performance
```

每個 Phase 的 PR 合併回 `refactor/baseline-checkpoint`，確認穩定後再考慮合回主線。
