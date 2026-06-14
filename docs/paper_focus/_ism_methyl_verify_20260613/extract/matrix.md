# ISM read×CpG 甲基矩陣 — 源碼稽核 (aspect=matrix)

> 稽核日期: 2026-06-13
> 指派檔案: `src/core/MatrixBuilder.cpp` + `include/core/MatrixBuilder.hpp` + `include/core/MethylationMatrix.hpp`
> 補充佐證檔（CSV 寫出 / 二值化 / 閾值）: `src/io/RegionWriter.cpp` · `src/core/RegionProcessor.cpp` · `include/core/Config.hpp`
> 紀律: 只引實際檔案 + 行號 + 原始碼片段; 不臆測; 疑點標 uncertain。

---

## L0 一眼結論

ISM 內部存在**兩種** read×CpG 矩陣表示，靠 `RegionProcessor::build_methylation_matrix()` 串接：

1. **`MatrixBuilder` 的 `matrix_`（`vector<vector<double>>`）= 連續機率矩陣**：值 ∈ [0.0, 1.0] 是甲基化「機率」（P(modified)），未覆蓋 = `-1.0`。**這是寫到 `methylation.csv` / `methylation_forward.csv` / `methylation_reverse.csv` 的源**。
2. **`MethylationMatrix`（Eigen）= 派生雙表示**：`raw_matrix`（double, NaN=missing）＋ `binary_matrix`（int, 1/0/-1）。**`binary_matrix` 在此處才被二值化產生**，由 `--methyl-high`/`--methyl-low`（預設 0.8 / 0.2）切，三態：1 甲基 / 0 非甲基 / -1（未覆蓋 或 中間模糊）。

**關鍵：二值化不在指派的三個檔案裡發生** —— `MatrixBuilder` 只存連續機率；二值切割發生在 `RegionProcessor.cpp:1417-1422`。

---

## L1 重點邏輯（資料流）

```
ReadParser → MethylCall{ref_pos, probability}
   │  (per read)
   ▼
MatrixBuilder::add_read()           [MatrixBuilder.cpp:17-39]  存 (ref_pos, probability) pair
   │
   ▼
MatrixBuilder::finalize()           [MatrixBuilder.cpp:41-104]
   ├ 收所有 read 的 pos → sort+unique → cpg_positions_ (column 軸)
   ├ matrix_[r].resize(cols, -1.0)  // -1.0 = no coverage 初始
   └ row[col] = prob                // 填連續機率
   │
   ├──► RegionWriter::write_matrix_csv()   methylation.csv（連續值，<0 → "NA"）
   │    RegionWriter::write_strand_matrices()  methylation_forward/reverse.csv
   │
   └──► RegionProcessor::build_methylation_matrix()  [RegionProcessor.cpp:1381-1429]
        └ 產 MethylationMatrix{raw_matrix(NaN), binary_matrix(1/0/-1)}  ← 二值化在此
```

---

## L2 細節 + 原始碼佐證

### 1. 矩陣「值」語意 = 連續機率（P(modified)），非二值

`MatrixBuilder` 內部 `matrix_` 存的是 `float probability` 直接放入 `double`，**無二值化**：

`src/core/MatrixBuilder.cpp:34-36`（add_read 存 probability）:
```cpp
for (const auto& call : methyl_calls) {
    calls.emplace_back(call.ref_pos, call.probability);
}
```

`src/core/MatrixBuilder.cpp:89-96`（finalize 填值，直接放 prob）:
```cpp
float prob = p.second;
auto it = pos_to_col.find(pos);
if (it != pos_to_col.end()) {
    row[it->second] = prob;
}
```

`include/core/MatrixBuilder.hpp:25-28` 明文定義值域:
```
* Value meanings:
* - [0.0, 1.0]: Methylation probability
* - -1.0: This read does not cover this CpG site (using -1 to represent NaN)
```

> ⚠ **uncertain（語意精確度）**：`call.probability` 究竟是 P(5mC) 單獨、還是「any-mod」/「5mC+5hmC 合併」機率，**不在本三檔可確定** —— 取決於上游 `MethylationParser`/`MethylCall` 如何由 MM/ML tag 解析（MEMORY 既有紀錄提及 max-collapse any-mod vs 5mC-only 口徑差異）。本檔僅能斷定：矩陣存的是「一個連續甲基化機率」。本檔範圍內無法分辨 5mC vs 5hmC。

### 2. 二值表示存在但**在下游派生**（非 MatrixBuilder）

`MethylationMatrix` header 宣告兩種表示並存 — `include/core/MethylationMatrix.hpp:22-23`:
```cpp
Eigen::MatrixXd raw_matrix;     ///< Raw methylation probabilities (0.0 - 1.0), NaN for missing
Eigen::MatrixXi binary_matrix;  ///< Binary methylation status: 1 (meth), 0 (unmeth), -1 (missing)
```

實際二值化在 `src/core/RegionProcessor.cpp:1410-1423`:
```cpp
double val = raw_matrix[i][j];
if (val < 0) {  // -1.0 indicates no coverage
    meth_mat.raw_matrix(i, j) = NAN;
    meth_mat.binary_matrix(i, j) = -1;
} else {
    meth_mat.raw_matrix(i, j) = val;
    // Binary threshold (configurable via --methyl-high/--methyl-low)
    if (val >= binary_methyl_high_) {
        meth_mat.binary_matrix(i, j) = 1;
    } else if (val <= binary_methyl_low_) {
        meth_mat.binary_matrix(i, j) = 0;
    } else {
        meth_mat.binary_matrix(i, j) = -1;  // Ambiguous
    }
}
```

**三態二值化**：`val ≥ high → 1`；`val ≤ low → 0`；**中間 (low, high) → -1（Ambiguous，與「未覆蓋」共用 -1 編碼）**。

閾值預設值 — `include/core/Config.hpp:33-34`:
```cpp
double binary_methyl_high = 0.8;  ///< Threshold for methylated (1) call
double binary_methyl_low = 0.2;   ///< Threshold for unmethylated (0) call
```
可由 CLI 覆寫（`--methyl-high` / `--methyl-low`，`include/utils/ArgParser.hpp:60,63`）；Config 驗證 `high > low` 且兩者 ∈ [0,1]（`src/core/Config.cpp:108-113`）。

### 3. 維度怎麼定（reads × CpG sites）

- **Row = read**（讀入順序）。`include/core/MatrixBuilder.hpp:22-23`: "Rows: reads (in order of reading)"。`add_read` 用 `reads_.size()` 當 read_id（`MatrixBuilder.cpp:22`）。
- **Col = unique CpG site，依基因組座標排序**。`include/core/MatrixBuilder.hpp:23`: "Cols: unique CpG sites (sorted by genomic coordinates)"。
  - 收集所有 read 的 `ref_pos` → `std::sort` + `std::unique` 去重（`MatrixBuilder.cpp:62-64`）：
    ```cpp
    std::sort(all_positions.begin(), all_positions.end());
    auto last = std::unique(all_positions.begin(), all_positions.end());
    all_positions.erase(last, all_positions.end());
    ```
  - 即 column 軸 = 「該 region 內，被任一 read 覆蓋過的 CpG 位點聯集」，非固定參考 CpG list。
- 維度分配 `MatrixBuilder.cpp:75-83`: `num_rows = reads_.size()`, `num_cols = cpg_positions_.size()`。
- `num_reads()` / `num_sites()`（`MethylationMatrix.hpp:25-30`）= `read_ids.size()` / `cpg_ids.size()`。

### 4. NA（未覆蓋）如何表示與計入

兩層編碼不同，需注意：

| 層 | 未覆蓋編碼 | 出處 |
|----|-----------|------|
| `MatrixBuilder::matrix_` (連續) | `-1.0`（初始填充值，未被任何 call 覆寫即保持） | `MatrixBuilder.cpp:82` `matrix_[r].resize(num_cols, -1.0);` |
| `MethylationMatrix::raw_matrix` (Eigen) | `NAN` | `RegionProcessor.cpp:1412` |
| `MethylationMatrix::binary_matrix` | `-1` | `RegionProcessor.cpp:1413` |
| `methylation.csv` 輸出 | 字串 `"NA"` | `RegionWriter.cpp:233-234` |

關鍵機制：矩陣**先全部初始化為 -1.0**，只有「該 read 在該 pos 有 MethylCall」才覆寫成機率（`MatrixBuilder.cpp:81-83` + `89-97`）。沒有 call 的格子保持 -1.0 = 未覆蓋。

> ⚠ **編碼歧義（重要）**：`binary_matrix` 的 `-1` **同時表示「未覆蓋」與「中間模糊機率 (low<val<high)」兩種情形**（`RegionProcessor.cpp:1413` vs `1422`），下游若要區分「真未覆蓋」vs「覆蓋但模糊」**需回 `raw_matrix`（NaN vs 非NaN）判斷**。`MethylationMatrix.hpp:23` 註解只寫 "-1 (missing)"，未提及 Ambiguous 也被併入 -1。

「計入」：CSV 連續矩陣以 `-1.0` 判斷分流（`matrix[r][c] < 0.0` → 寫 "NA"，否則寫機率，`RegionWriter.cpp:233-237`）。是否在統計時排除 NA，屬下游 `DistanceMatrix` / `PerCpgAsm` 行為，**不在本三檔範圍**（標 uncertain）。

### 5. methylation.csv 怎麼產生

`src/io/RegionWriter.cpp:217-243` `write_matrix_csv()`：
- 輸出檔名固定 `region_dir + "/methylation.csv"`（`RegionWriter.cpp:219`）。
- **分隔符 = 逗號 `,`**（CSV），表頭 `read_id` + 各 CpG position（`RegionWriter.cpp:222-226`）。
- 每列一條 read，第一欄為 row index `r`（非 read_name）：
  ```cpp
  for (size_t r = 0; r < matrix.size(); r++) {
      ofs << r;
      for (size_t c = 0; c < matrix[r].size(); c++) {
          ofs << ",";
          if (matrix[r][c] < 0.0) { ofs << "NA"; }
          else { ofs << std::fixed << std::setprecision(4) << matrix[r][c]; }
      }
  ```
- **連續機率寫出，固定 4 位小數**（`std::setprecision(4)`，`RegionWriter.cpp:236`）。即 methylation.csv **存的是連續 P(mod)，非二值**。
- 餵入的 `matrix` = `MatrixBuilder::get_matrix()`（連續 `matrix_`）。

> ⚠ uncertain：`.claude/rules/output-structure.md` 描述輸出含 `methylation.csv`（與此一致），但同 rule 提到 `distance_matrix_[METRIC].csv`；本稽核未交叉核對全部輸出檔清單，僅確認 methylation 系列三檔。

### 6. forward/reverse 拆分邏輯

`src/io/RegionWriter.cpp:245-316` `write_strand_matrices()`，**僅在 `output_strand_matrices_` 為 true 時觸發**（`RegionWriter.cpp:114-116`）：

- 依 `reads[i].strand` 分桶（`RegionWriter.cpp:252-258`）：
  ```cpp
  for (size_t i = 0; i < reads.size(); i++) {
      if (reads[i].strand == Strand::FORWARD) { forward_indices.push_back(i); }
      else if (reads[i].strand == Strand::REVERSE) { reverse_indices.push_back(i); }
  }
  ```
  → **strand 不是 FORWARD 也不是 REVERSE 的 read 兩桶都不收**（被丟棄，標 uncertain：是否有 UNKNOWN strand）。
- forward 寫 `methylation_forward.csv`、reverse 寫 `methylation_reverse.csv`（`RegionWriter.cpp:262, 291`）。
- **column 軸不變**：兩個 strand 矩陣共用同一份全域 `cpg_positions`（`RegionWriter.cpp:266-268, 295-297`），只是 **row 子集**。
- **多一欄 `original_read_id`**：表頭 `read_id,original_read_id,<positions...>`（`RegionWriter.cpp:265, 294`），`read_id` 是 strand 內重新編號（`new_id`），`original_read_id` 是全域 row index（`orig_id = forward_indices[new_id]`，`RegionWriter.cpp:273-274, 302-303`）。
- NA/機率寫法與主 CSV 相同（`<0.0 → "NA"`，否則 4 位小數機率，`RegionWriter.cpp:277-281, 306-310`）。
- 計數：metadata 另記 `num_forward` / `num_reverse`（`RegionWriter.cpp:90-97`）。

**結論：forward/reverse 拆分 = 純 row-wise 過濾（按 read strand），值語意/column 軸/NA 編碼完全沿用主矩陣，不改機率本身。**

---

## L3 溯源摘要（檔:行）

| 事實 | 來源 |
|------|------|
| 連續值 ∈ [0,1] = 機率 | MatrixBuilder.hpp:25-28; MatrixBuilder.cpp:34-36,89-96 |
| 未覆蓋 = -1.0 (連續) | MatrixBuilder.cpp:82 |
| Row=read, Col=sorted unique CpG | MatrixBuilder.hpp:22-23; MatrixBuilder.cpp:62-66,75-83 |
| 二值 1/0/-1 派生 + 閾值 | RegionProcessor.cpp:1410-1423 |
| 閾值預設 0.8 / 0.2 | Config.hpp:33-34 |
| binary -1 = missing OR ambiguous（歧義） | RegionProcessor.cpp:1413 vs 1422 |
| methylation.csv 連續+4位小數+NA | RegionWriter.cpp:219-240 |
| forward/reverse = row 子集 + original_read_id 欄 | RegionWriter.cpp:245-316 |

---

## 附註：兩 header 的「Eigen vs vector」差異不是矛盾

`MatrixBuilder.hpp:29` 明寫 "Uses vector<vector<double>> instead of Eigen to avoid dependencies"，而 `MethylationMatrix.hpp:3` 用 Eigen。兩者是**不同類別、不同角色**：`MatrixBuilder` = 採集/組裝層（連續機率 + -1.0 missing），`MethylationMatrix` = 下游分析容器（Eigen raw+binary）。串接點 `RegionProcessor::build_methylation_matrix()`（RegionProcessor.cpp:1381-1429）。`MatrixBuilder.cpp:9-15` 頂部註解描述的舊型別（`map<uint32_t,double>` / `uint32_t` positions）已過時，與現行 `vector<pair<int32_t,float>>`（hpp:111）不符 — 屬陳舊註解，非邏輯 bug。
