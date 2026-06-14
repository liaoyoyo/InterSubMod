# Verify — claim `5mc_only`：「ISM C++ 只讀 5mC，5hmC 不用/被塌縮」

> **獨立對抗審查** · 稽核日期 2026-06-13 · 預設懷疑、回原始碼複核
> **驗法**：直接 Read `src/core/MethylationParser.cpp`（291 行全文）+ `include/core/MethylationParser.hpp` + repo-scoped grep 驗證 emission 路徑與 caller。未跑 binary、未臆測；不確定處標 **[uncertain]**。

---

## VERDICT：**SUPPORTED**（claim 成立，但需精確措辭）

「ISM C++（`MethylationParser`）只解析 5mC（`C+m?`），5hmC（`C+h?`）**不被使用**」屬實。

⚠ **措辭精確化（針對 claim 內「被塌縮 / max-collapse」字眼）**：5hmC **不是被 max-collapse、不是雙列砍半、不是合併**，而是**機率被 offset 純粹跳過、從不產生任何 call**。若 claim 的「被塌縮」意指 collapse/merge → 該描述**對本 C++ 不準確**（collapse/max 是外部 Python/MSA 工具側行為，非 ISM C++）。對「只讀 5mC、5hmC 不用」這核心命題＝SUPPORTED。

---

## 原碼證據（檔:行 + 片段）

### E1 — 唯一硬編碼目標 = `C+m?`（5mC），block scan 只找它
`src/core/MethylationParser.cpp:59-69`：
```cpp
// Check if this block is "C+m?" (target modification)
const char* target = "C+m?";
const size_t tlen  = 4;
bool is_target = (static_cast<size_t>(block_end - block_start) >= tlen) &&
                 std::memcmp(block_start, target, tlen) == 0;
if (is_target) {
    deltas      = parse_mm_tag(mm_str, "C+m?");
    found_target = true;
    break;
}
```
- 命中 `C+m?` 才取其 deltas 並 `break`；找不到 `C+m?` → `return calls;`（空）（`MethylationParser.cpp:85-87`）。
- `parse_mm_tag` 預設 `mod_code = "C+m?"`（`MethylationParser.hpp:76`）。
- **repo-scoped grep 證實**：全 src/include/tests 中 `parse_mm_tag` 唯一非 default-參數呼叫 = `MethylationParser.cpp:66` 傳的也是 `"C+m?"`（grep 輸出：僅 cpp:66 caller + cpp:201 定義 + hpp:76 宣告）。**無任何路徑以 `C+h?` 為參數呼叫解析。**

### E2 — 5hmC（`C+h?`）只用來推進 `ml_offset`，從不 emit call
`src/core/MethylationParser.cpp:71-82`（非 target block 的處理）：
```cpp
// Not our target: count deltas (commas after the first ',') to advance ml_offset.
const char* comma = ... std::memchr(block_start, ',', ...);
if (comma) {
    int delta_count = 0;
    for (const char* q = comma + 1; q < block_end; ++q) {
        if (*q == ',') ++delta_count;
    }
    delta_count++;            // +1 for last value (no trailing comma)
    ml_offset += delta_count; // 只推進 ML 陣列指標，跳過 5hmC 那段機率
}
```
- 設計註解明說 ML 陣列順序 `[C+h? probs..., C+m? probs...]`（`MethylationParser.cpp:41`），故需 `ml_offset` 跳過前段 5hmC 機率，**讀取時直接從 5mC 段起算**（`prob = ml_data[ml_offset + delta_idx]/255.0f`，`:140`/`:180`）。
- `C+h?`/`5hmC` 在整檔**只出現於註解**（`MethylationParser.cpp:39,41,104`；`MethylationParser.hpp:33`），**無任何 active 解析/運算路徑**。

### E3 — MethylCall 的 emission 點只有 2 處，每處只帶「1 個 5mC 機率」
`grep "emplace_back"` 於本檔僅 2 命中：
- forward：`MethylationParser.cpp:181`：`calls.emplace_back(ref_pos_0based + 1, prob);`
- reverse：`MethylationParser.cpp:142-143`：`calls.emplace_back(ref_pos_0based, prob);`

`MethylCall` 結構只有 `{int32_t ref_pos; float probability;}`（`MethylationParser.hpp:14-22`）—— **不帶 modified-base 種類、不帶 strand、不帶第二個 mod 通道**。每個 CpG 位置只輸出一個 5mC 機率值。

### E4 — 無 max-collapse / 無雙列 / 無 5mC+5hmC 合併運算
- C++ 路徑全程**只有一個 mod 通道**（5mC）。`prob` 直接 = `ml_data[...]/255.0f`（線性，無 0.5 偏移、無分箱、無 clip；`:140`/`:180`，header `:39`），**不存在對兩個修飾取 max / sum / 合併的任何程式碼**。
- caller `src/core/ReadAggregator.cpp:62`：`methyl_parser.parse_read(b, ref_seq, region_start);` —— 上游直接收 sparse 的 `vector<MethylCall>`，無再注入 5hmC。

---

## 對 extract 既有描述的核對結果

extract `parser.md`（F2, L0, L145）與 `percpg.md`（L0#3, C 節）對此 claim 的描述**與原碼逐行吻合**：
- `parser.md:37-49`（F2「只解析 5mC，5hmC 被跳過機率、不輸出」）— **核實正確**。
- `parser.md:145`「不 collapse、不雙列、不取 max」— **核實正確**。
- `percpg.md:64-74`（C 節「dup-bug 不存在於本 C++；ISM 天生 5mC-only 單機率單 column」）— **核實正確**（含 `MatrixBuilder` 以位置為唯一鍵去重，本輪未重讀 MatrixBuilder.cpp，沿用 extract，標 [uncertain-downstream]）。

**結論：extract 對 `5mc_only` 的稽核無需校正，原碼複核一致。**

---

## correction（針對「舊理解」可能的誤區）

若有任何先前敘述把 ISM C++ 寫成「5mC+5hmC **max-collapse**」或「5mC+5hmC **雙列砍半**」——

**正確版**：ISM C++（`MethylationParser`）**只解析 `C+m?`(5mC) 一個通道**；5hmC（`C+h?`）的機率僅被 `ml_offset` 跳過（`MethylationParser.cpp:71-82`），**從不參與任何運算、從不產生 call**。所謂「5mC+5hmC max-collapse 取 any-mod」「雙列砍半 dup-bug」是**外部 Python/MSA（MethylSomaticAnalysis）抽取工具側的 Level-1 處理**（見 memory `reference_msa_vs_ism_tool`、`project_zar1l_brca2_asm_verification`：buggy −0.054 是 MSA 雙列 artifact），**不屬本 C++ binary**。故描述 ISM C++ 行為時應寫「**5mC-only**」，不可寫「max-collapse」或「雙列」。

---

## uncertain（讀不確定 / 邊界）

1. **[uncertain]** 5hmC 被 offset 跳過依賴 MM block 順序 `C+h?` 在前、`C+m?` 在後（`:41` 註解假設）+ 非 target block 含 `,` 才推進 offset（`:75`）。順序若反或空 block，offset 對位邏輯未測。對標準 ONT dorado MM tag 正常無影響。**屬對位魯棒性，不改「只用 5mC」此核心結論。**
2. **[uncertain-downstream]** 下游 `MatrixBuilder` 是否仍維持「1 CpG → 1 column」未在本輪重讀 .cpp（沿用 extract `percpg.md` L71 之描述）。但因 parser 層每 CpG 僅 emit 單一 5mC 機率，即使下游也不可能憑空生出 5hmC 列。
3. `is_cpg_site()` helper（`:283-288`）為 dead-ish（主路徑用 inline CpG 檢查 `:176-178`/`:138-139`），與本 claim 無關。
