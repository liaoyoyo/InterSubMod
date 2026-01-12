關於現在發現已知錯誤或問題列表，確認問題、分析問題、解決問題，並測試確認真實解決

## 2024-11-25 甲基化位點聚合邏輯確認與效能優化

### 1. 問題確認

**用戶疑慮**：懷疑 `MatrixBuilder` 在聚合甲基化矩陣時，是否僅依據第一條 Read 的 CpG 位點來決定矩陣的欄位（Columns），導致後續 Reads 若包含第一條 Read 沒有的位點時，該位點會被忽略。

**調查結果**：

- **程式碼分析** (`src/core/MatrixBuilder.cpp`)：
  - `finalize()` 函式中，遍歷了 **所有** Reads (`read_methyl_data_`) 的所有甲基化呼叫 (`calls`)。
  - 使用 `std::set<int32_t> unique_cpg_positions` 收集所有出現過的唯一 CpG 位點。
  - 因此，邏輯上是取所有 Reads 位點的 **聯集 (Union)**，而非僅取第一條 Read 的位點。
- **數據驗證** (`output/full_vcf_test/chr19_29283968/chr19_29282968_29284968/methylation.csv`)：
  - 觀察 Read 0（第一列）與 Read 3。
  - Read 0 在位點 `29283558` 為 `NA`（表示該 Read 未覆蓋或無數據）。
  - Read 3 在位點 `29283558` 有數值 `0.0078`。
  - 該位點 `29283558` 存在於 CSV 的標題列中（Column index 3）。
  - **結論**：系統已正確包含並聚合了後續 Reads 獨有的位點，未發生忽略情況。
  - **再次驗證 (2025-11-26)**：
    - 檢查了 SNV `chr19:29,283,968` 的輸出。
    - 確認位點 `29,283,631` 和 `29,283,830` 存在於 `methylation.csv` 的標題中。
    - Read 0 包含這些位點的數據，而 Read 3 雖然沒有這些位點的數據（顯示為 NA），但包含其他位點（如 `29,283,558`），這證明了矩陣正確地取了聯集。
    - **確認無誤**：目前的邏輯是正確的，不需要針對「忽略位點」進行修復。

### 2. 效能與流程優化建議

雖然功能邏輯正確，但在檢視程式碼時發現 `MatrixBuilder` 的實作有優化空間，可提升大量 Reads 時的效能。

**優化項目**：

1. **資料結構優化**：
    - **現狀**：`read_methyl_data_` 使用 `std::map<int, vector...>`。由於 `read_id` 是連續整數（0, 1, 2...），使用 `map` 會造成不必要的紅黑樹操作與記憶體破碎。
    - **修正**：改用 `std::vector<std::vector<std::pair<int32_t, float>>>`。存取時間從 $O(\log N)$ 降為 $O(1)$。

2. **排序與去重優化** (`finalize()` 函式)：
    - **現狀**：使用 `std::set` 來收集唯一位點。`set` 的插入操作成本較高（節點配置 + 樹平衡）。
    - **修正**：改用 `std::vector` 收集所有位點，再使用 `std::sort` + `std::unique`。這對快取（Cache）更友善且通常更快。

### 3. 修正計畫

- 修改 `include/core/MatrixBuilder.hpp` 與 `src/core/MatrixBuilder.cpp`，實施上述效能優化。
