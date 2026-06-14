# Verify — claim `matrix_semantics`

> 審查員: 獨立對抗審查（預設懷疑 → 原碼 + 實值雙重佐證）
> 日期: 2026-06-13
> claim 全文:「methylation.csv 值 = 連續 P(5mC) 0-1、per-read×per-CpG、NA=未覆蓋」

---

## Verdict: **PARTIAL**（核心三點 SUPPORTED；唯一不成立子句 = 「P(**5mC**)」的單修飾語意）

| 子 claim | 判定 | 一句理由 |
|----------|------|----------|
| 值 = 連續機率 0–1（非二值） | **SUPPORTED** | 原碼 `ml_data/255.0f` + 實測 240 distinct 值 / 25.2% 中間機率 / 全在 [0,1] |
| per-read × per-CpG | **SUPPORTED** | row=read（讀入序）、col=sorted unique CpG genomic pos；原碼 + metadata 57×203 對齊 |
| NA = 未覆蓋 | **SUPPORTED** | 矩陣初值 `-1.0`，`<0.0 → "NA"`；實測 50.2% NA = 稀疏覆蓋 |
| 「P(**5mC**)」單一 5mC 機率 | **NOT ESTABLISHED（uncertain，偏向 over-specified）** | 原碼 `MatrixBuilder`/`RegionWriter`/`MethylationParser` 寫的是一個泛 `probability`；5mC vs 5hmC 取決於 MM/ML 解析口徑，本路徑無法斷定只含 5mC |

**一句話**：claim 的「連續機率 0-1 / per-read×per-CpG / NA=未覆蓋」三點原碼+實值雙重證實；但把那機率標成「P(**5mC**)」是 over-specification —— 原碼層只能斷定是「一個甲基化機率」，無法保證是純 5mC（可能含/合併 5hmC）。故 PARTIAL 而非全 SUPPORTED。

---

## 證據 A — 原始碼（值=連續機率，非二值）

### A1. 機率來源 = ML tag / 255（8-bit → 連續）
`src/core/MethylationParser.cpp:140`（reverse strand）與 `:180`（forward strand）:
```cpp
float prob = ml_data[ml_offset + delta_idx] / 255.0f;
```
→ BAM `ML` tag 是 0–255 的 8-bit；`/255.0f` 得 [0,1] 連續機率（解釋實測 1/255 granularity）。

`MethylCall` 結構 — `include/core/MethylationParser.hpp:15-16`:
```cpp
int32_t ref_pos;    ///< 1-based genomic position (hg38 coordinates)
float probability;  ///< Methylation probability [0.0, 1.0]
```

### A2. MatrixBuilder 直接存機率，**無二值化**
`src/core/MatrixBuilder.cpp:34-35`（add_read 存 (pos, prob)）:
```cpp
for (const auto& call : methyl_calls) {
    calls.emplace_back(call.ref_pos, call.probability);
}
```
`src/core/MatrixBuilder.cpp:91-95`（finalize 直接填 prob，不過閾值）:
```cpp
float prob = p.second;
auto it = pos_to_col.find(pos);
if (it != pos_to_col.end()) {
    row[it->second] = prob;
}
```
header 值域明文 — `include/core/MatrixBuilder.hpp:25-27`:
```
* Value meanings:
* - [0.0, 1.0]: Methylation probability
* - -1.0: This read does not cover this CpG site (using -1 to represent NaN)
```

### A3. methylation.csv 寫出連續值（4 位小數），未覆蓋→"NA"
`src/io/RegionWriter.cpp:229-239`:
```cpp
for (size_t r = 0; r < matrix.size(); r++) {
    ofs << r;
    for (size_t c = 0; c < matrix[r].size(); c++) {
        ofs << ",";
        if (matrix[r][c] < 0.0) {  // -1.0 indicates no coverage
            ofs << "NA";
        } else {
            ofs << std::fixed << std::setprecision(4) << matrix[r][c];
        }
    }
    ofs << "\n";
}
```
餵入的 `matrix` = `MatrixBuilder::get_matrix()`（連續 `matrix_`，`MatrixBuilder.hpp:59-61`），**非** `MethylationMatrix::binary_matrix`。
→ methylation.csv 存連續 P(mod)，二值化（1/0/-1，閾值 0.8/0.2）發生在下游 `RegionProcessor.cpp:1410-1423`，**不寫進 methylation.csv**。

---

## 證據 B — per-read × per-CpG 維度

- **Row = read（讀入序）**：`include/core/MatrixBuilder.hpp:22` "Rows: reads (in order of reading)"；`add_read` 用 `reads_.size()` 當 read_id（`MatrixBuilder.cpp:22`）。CSV 第一欄 = row index `r`（`RegionWriter.cpp:230`，整數 0..N-1，非 read name）。
- **Col = unique CpG genomic position，座標排序**：`MatrixBuilder.cpp:62-66`：
  ```cpp
  std::sort(all_positions.begin(), all_positions.end());
  auto last = std::unique(all_positions.begin(), all_positions.end());
  all_positions.erase(last, all_positions.end());
  cpg_positions_ = std::move(all_positions);
  ```
  CSV header = `read_id` 後接各 genomic pos（`RegionWriter.cpp:222-225`，**是 position 非 cpg_id**）。
- **CpG context-validated**：`MethylationParser.cpp:138-139`/`176-178` 要求 ref `C`+`G` 才算 CpG。

---

## 證據 C — 實值佐證（flagship methylation.csv，與原碼交叉）

檔案: `research/flagship_chr2_18086020_20260612/ism_out/anchor_18086020/chr2/chr2_18086020/chr2_18071020_18101020/methylation/methylation.csv`

實測（Python 全矩陣掃描，2026-06-13）:
| 量 | 實測值 |
|----|--------|
| 維度 | **57 reads × 203 CpG**（header 204 欄含 read_id；58 行含 header） |
| 值範圍 | min=0.0000, max=1.0000，**out-of-[0,1] = 0** |
| distinct 值數 | **240**；min positive step = **0.003900**（≈ 1/255 = 0.003922）→ 8-bit ML 量化 |
| frac ≤0.05 | 0.254 |
| frac ≥0.95 | 0.494 |
| **frac 嚴格中間 (0.05<v<0.95)** | **0.252** ← 1/4 cell 是中間機率，二值 call 不可能出現 |
| NA | **5807 / 11571 cells (50.2%)**，numeric=5764 |
| header 首 6 欄 | `read_id,18071239,18071625,18071954,18072528,18072619`（genomic pos） |
| row read_id=0 首值 | `0.9882,0.9333,0.0549,0.0235,0.9961,0.8941,...`（連續，含中間值如 0.0549） |

metadata.txt 對齊（同區域目錄）: `Matrix Dimensions: 57 × 203`、`Num Reads: 57 (Forward 32, Reverse 25)` → 與實測完全一致。

→ **240 distinct 值 + 1/255 step + 25.2% 嚴格中間機率** = 連續機率的決定性實證；二值化矩陣只會有 {0,1,NA} 三類值。NA 50.2% = read 只覆蓋自身對齊跨度內 CpG 的正常稀疏。

---

## correction — 唯一需校正子句

claim 寫「連續 **P(5mC)**」。**「連續機率 0-1」「per-read×per-CpG」「NA=未覆蓋」三點為真**，但「P(5mC)」的「5mC」是 over-specification：

- 本路徑原碼（`MethylationParser.cpp:140/180` → `MatrixBuilder` → `RegionWriter`）只搬運一個泛型 `probability` / `float prob = ml_data[...]/255.0f`，**不在這層區分 5mC vs 5hmC**。
- 究竟是 P(5mC) 單修飾、還是 any-mod / 5mC+5hmC 合併（如 max-collapse），取決於上游 MM/ML 多修飾解析口徑 —— 本三檔（MatrixBuilder/RegionWriter/MethylationMatrix）+ MethylationParser 的這兩行**無法斷定**。MEMORY 既有紀錄亦提及「max-collapse any-mod vs 5mC-only 口徑差異」存在。
- **校正後正確版**：「methylation.csv 值 = 連續**甲基化機率** P(modified) ∈ [0,1]、per-read × per-CpG、NA=該 read 在該 CpG 未覆蓋/無 ML call」。若要主張「**5mC**」需另行查 MM/ML 多修飾解析口徑（標 uncertain，不在本 claim 已證範圍內）。

---

## 其他標註 uncertain（不影響本 verdict）

1. **binary_matrix 的 -1 歧義**：下游 `MethylationMatrix::binary_matrix`（Eigen，`MethylationMatrix.hpp:23`）的 `-1` 同時表「未覆蓋」與「中間模糊 (low<val<high)」兩種（`RegionProcessor.cpp:1413` vs `1422`）。**但這是下游派生表示，不寫入 methylation.csv**，與本 claim（針對 methylation.csv 連續值）無關 —— 連續層的 NA 明確只表「未覆蓋」（`-1.0` 初值，唯一來源）。
2. **1-based vs 0-based 命名不一致**：`MethylCall.ref_pos` 註解寫 "1-based"（hpp:15），forward strand `emplace_back(ref_pos_0based + 1, ...)`（`MethylationParser.cpp:181`）對齊 1-based，但 reverse strand `emplace_back(ref_pos_0based, ...)`（:142，註解稱 `-1+1`）。不影響「值語意/維度/NA」結論；CpG 位置精確比對才需查證。
3. 5mC vs 5hmC 之外，本稽核未逐值核對 forward/reverse 拆分檔（多一欄 `original_read_id`，`RegionWriter.cpp:265,294`）。

---

## 溯源摘要（檔:行）

| 事實 | 來源 |
|------|------|
| prob = ML/255（連續來源） | MethylationParser.cpp:140, :180 |
| MethylCall.probability [0,1] | MethylationParser.hpp:15-16 |
| 存機率無二值化 | MatrixBuilder.cpp:34-35, 91-95 |
| 值域 [0,1] + -1.0 missing | MatrixBuilder.hpp:25-27 |
| matrix 初值 -1.0 | MatrixBuilder.cpp:82 |
| CSV 連續+4 位小數+NA | RegionWriter.cpp:229-239 |
| header=read_id+genomic pos | RegionWriter.cpp:222-225 |
| Row=read / Col=sorted unique CpG | MatrixBuilder.hpp:22-23; MatrixBuilder.cpp:62-66 |
| 二值化在下游（不入 csv） | RegionProcessor.cpp:1410-1423 |
| 實測 240 distinct / 25.2% mid / 50.2% NA / 57×203 | flagship methylation.csv（本日 Python 掃描）+ metadata.txt |
