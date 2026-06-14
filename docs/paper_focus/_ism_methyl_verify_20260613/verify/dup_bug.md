# 對抗驗證 — claim「dup_bug」: 「5mC+5hmC 雙列砍半 dup-bug」是否存在於本 ISM C++?

> 驗證日期 2026-06-13 · 獨立對抗審查（預設懷疑）· 已直接讀實際 C++ 源檔複核（非僅信 extract）· 引行號 + 原始碼片段

## verdict: **REFUTED**（dup-bug 不存在於本 ISM C++）

待驗 claim：「5mC+5hmC 雙列砍半 dup-bug」在本 ISM C++（`MethylationParser` / `MatrixBuilder`）有無把同一 CpG 的 5mC+5hmC 兩列平均/砍半。

**結論：本 ISM C++ 不存在此 dup-bug。** ISM C++ 解析路徑天生 **5mC-only（只取 `C+m?`）、每 CpG 單一機率、每 position 單一 column**，從不產生 5hmC 列，故物理上無從「兩列平均/砍半」。dup-bug 屬外部抽取工具 MSA（MethylSomaticAnalysis）的 Python Level-1 處理，**非本 repo 的 C++ binary**。

---

## evidence（直接讀源檔，非信 extract）

### E1 — Parser 只鎖定 5mC `C+m?`，5hmC `C+h?` 永不 emit call

`src/core/MethylationParser.cpp`（`parse_read` 為唯一 MethylCall 生產者）：

- **line 60-63**：唯一 target 是 hard-coded `"C+m?"`（5mC）：
  ```cpp
  const char* target = "C+m?";
  const size_t tlen  = 4;
  bool is_target = (static_cast<size_t>(block_end - block_start) >= tlen) &&
                   std::memcmp(block_start, target, tlen) == 0;
  ```
- **line 65-69**：命中 `C+m?` 才取 deltas 並 `break`；非 target block 不取其機率：
  ```cpp
  if (is_target) {
      deltas      = parse_mm_tag(mm_str, "C+m?");
      found_target = true;
      break;
  }
  ```
- **line 71-82**：非 target block（含 `C+h?`）**只數 delta 個數來推進 `ml_offset`**，藉此把 ML 陣列中 5hmC 那段機率「跳過」，從不轉成 call：
  ```cpp
  // Not our target: count deltas (...) to advance ml_offset.
  ... delta_count++; ml_offset += delta_count;
  ```
- **line 85-87**：找不到 `C+m?` 直接 return 空：`if (!found_target) { return calls; }`
- `parse_mm_tag` 預設 `mod_code = "C+m?"`：`include/core/MethylationParser.hpp:76`。
- 註解明列 ML 順序為 `[C+h? probs..., C+m? probs...]`（`MethylationParser.cpp:41`），故需 offset 跳過 5hmC 段。

**grep 佐證**（`C+h` / `5hmC` 在 parser+matrix 全部出現位置）：
```
include/core/MethylationParser.hpp:33  (註解 MM format)
src/core/MethylationParser.cpp:39      (註解)
src/core/MethylationParser.cpp:41      (註解 — ML 順序)
```
→ `C+h?` 在 ISM C++ **僅出現在註解，無任何 active 解析路徑**；`5hmC` 字串本身在 parser/matrix 程式碼中不出現。

### E2 — 每個 CpG 只 emit 一個機率值（單通道，無雙列來源）

`MethylationParser.cpp`：
- forward（line 180-181）：`float prob = ml_data[ml_offset + delta_idx] / 255.0f;  calls.emplace_back(ref_pos_0based + 1, prob);`
- reverse（line 140-142）：`float prob = ml_data[ml_offset + delta_idx] / 255.0f;  calls.emplace_back(ref_pos_0based, prob);`

每命中一個 ref-context-驗證過的 CpG，只 push **一個** `MethylCall{ref_pos, prob}`。C++ 路徑根本沒有第二個 mod 通道，因此不需要、也不存在「5mC+5hmC 取 max / 取平均 / 砍半合併」步驟。

### E3 — Matrix 以「position 為唯一鍵」去重，同一 CpG 不可能變兩 column

`src/core/MatrixBuilder.cpp`（直接讀源檔複核）：
- **line 34-35**：`add_read` 存 `(call.ref_pos, call.probability)` pair（沿用 parser 單機率）。
- **line 55-64**：收所有 read 的 `p.first`（= ref_pos）→ `std::sort` + `std::unique` 去重：
  ```cpp
  for (const auto& p : calls) { all_positions.push_back(p.first); }
  std::sort(all_positions.begin(), all_positions.end());
  auto last = std::unique(all_positions.begin(), all_positions.end());
  all_positions.erase(last, all_positions.end());
  ```
- **line 69-72**：`pos_to_col[cpg_positions_[i]] = i;` → **1 position ↔ 1 column**。
- **line 89-96**：填值 `row[it->second] = prob;`（同 position 同格直接覆寫，非加總/平均）。

column 軸 = unique CpG position 聯集。同一 CpG 座標只會有一個 column，**結構上不可能因 5mC/5hmC 拆成兩列再砍半**。

### E4 — ISM C++ 中 MethylCall 的唯一生產者就是這個 5mC-only parser

`grep` 全 `src/core src/io include/core`：唯一呼叫 `parse_read` 的點是 `src/core/ReadAggregator.cpp:62`：
```cpp
auto methyl_calls = methyl_parser.parse_read(b, ref_seq, region_start);
```
沒有任何其他 5hmC 解析路徑、沒有別的 `MethylCall` 生產者餵進 `MatrixBuilder`。故「dup-bug 在 ISM C++ 別處潛伏」也被排除。

---

## correction（舊理解校正）

**校正點**：本 claim 名稱「dup_bug」易被誤讀成「ISM C++ 有 5mC+5hmC 雙列砍半 bug」。**正確版：該 dup-bug 不屬本 ISM C++**。

- 舊 memory `project_zar1l_brca2_asm_verification` 提到「buggy −0.054 是 MSA Level1 5mC+5hmC 雙列砍半 artifact」——此 artifact 的歸屬是 **外部 MSA（MethylSomaticAnalysis）Python 抽取工具的 Level-1 處理**，與本 repo C++ binary 無關。
- ISM C++ 的口徑是 **5mC-only**（只解析 `C+m?`，把 `C+h?` 機率用 offset 略過），這與 memory `reference_msa_vs_ism_tool`、`project_code_methodology_audit_2026_06_10`（「C++=5mC-only+germline 軸」）一致。
- 因此 BRCA2 Δβ 等量化若用 ISM C++ 直接算，**本就不會碰到雙列砍半**；只有走 MSA Level-1（5mC+5hmC 雙列）路徑才會。兩條口徑不可混淆。

## uncertain（誠實標註，不影響 verdict）

1. 本輪**未直接讀 MSA（MethylSomaticAnalysis）源碼**佐證其雙列砍半機制細節 — dup-bug 屬 MSA 之認定來自既有 memory，非本輪原碼複核。本驗證只負責「dup-bug 不在本 ISM C++」這半邊，已用原碼定案。
2. 5hmC 被 offset 跳過依賴 MM block 內 `C+h?` 在 `C+m?` 之前的順序假設（`MethylationParser.cpp:41` 註解）；對標準 ONT MM tag 成立。此為 offset 正確性議題，與 dup-bug（雙列砍半）無關。

---

## 溯源（檔:行）

| 事實 | 來源 |
|------|------|
| 唯一 target = `C+m?`（5mC） | `src/core/MethylationParser.cpp:60-63` |
| 命中才取 deltas、否則 break | `src/core/MethylationParser.cpp:65-69` |
| `C+h?` 機率只被 offset 略過、不 emit | `src/core/MethylationParser.cpp:71-82` |
| 找不到 5mC 即 return 空 | `src/core/MethylationParser.cpp:85-87` |
| 每 CpG 單一 prob（fwd/rev） | `src/core/MethylationParser.cpp:140-142, 180-181` |
| `C+h?`/5hmC 僅在註解出現 | grep: `MethylationParser.cpp:39,41` + `hpp:33` |
| position 去重 → 1 pos = 1 col | `src/core/MatrixBuilder.cpp:55-72` |
| 填值為覆寫（非平均/砍半） | `src/core/MatrixBuilder.cpp:89-96` |
| MethylCall 唯一生產者 | `src/core/ReadAggregator.cpp:62` |
