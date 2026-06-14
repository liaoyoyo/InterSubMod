# PerCpgAsm 甲基矩陣 ASM 計算 — 源碼稽核

> aspect = `percpg` · 稽核日期 2026-06-13 · 只讀實際檔案、引行號、引原始碼片段佐證

## L0 一眼結論

1. `PerCpgAsm.cpp` 用 **per-CpG Fisher's exact test**（HP1-family vs HP2-family 的 2×2 meth/unmeth 列聯表）+ **BH-FDR** 校正算 per-CpG ASM；另有 NME（entropy）與 Epipolymorphism 兩組 epiallele 異質性度量。
2. 引用註解在 **`src/core/PerCpgAsm.cpp:6`**，現寫 `DAMEfinder (Orjuela 2020)` — **不是** `De Waele 2020`。原本（commit `5b78f1a`, 2026-04-19）誤寫 `De Waele 2020`，於 **commit `891e04b`（2026-06-10）已矯正**為 `Orjuela 2020`。
3. ISM（本 C++ binary）直接讀 BAM `MM/ML` tag（`MethylationParser`），**只解析 `C+m?`（5mC），不產生 5hmC 列** → **5mC+5hmC 雙列砍半 dup-bug 不存在於本 C++**；該 dup-bug 屬外部抽取工具 MSA（MethylSomaticAnalysis）。

---

## L1 重點邏輯

### A. per-CpG ASM 怎麼算（統計檢定 = Fisher）

主函式 `compute_per_cpg_asm`（`PerCpgAsm.cpp:228-375`）三步：

**Step 0 — HP-family 分組**（`PerCpgAsm.cpp:238-258`）
- read 的 `hp_tag` 映射成 family：`"1"|"HP1"|"1-1"` → HP1-family；`"2"|"HP2"|"2-1"` → HP2-family（`PerCpgAsm.cpp:244-250`）。
- 兩組各需 ≥ `min_reads_per_group`（default **5**，header `PerCpgAsm.hpp:66`）否則 return 空 result（`PerCpgAsm.cpp:256-258`）。

**Step 1 — per-CpG Fisher's exact test**（`PerCpgAsm.cpp:262-338`）
- 對每個 CpG column `j`，用 `meth_mat.raw_matrix(idx,j)` 的原始機率，**以 0.5 二值化**計 meth/unmeth：`if (val > 0.5) hp1_meth++ else hp1_unmeth++`（`PerCpgAsm.cpp:282-283`、HP2 同 `289-290`）；`NaN` 跳過（`PerCpgAsm.cpp:281`）。
- 組 2×2 列聯表（rows = HP group、cols = meth/unmeth）：`table.a=hp1_meth; b=hp1_unmeth; c=hp2_meth; d=hp2_unmeth`（`PerCpgAsm.cpp:298-302`）。
- 呼叫 `FisherExact::test_2x2(table)` 取 p-value（`PerCpgAsm.cpp:304-305`）；`FisherConfig.seed=42`（`PerCpgAsm.cpp:267`）。
- Delta = HP1 meth 比例 − HP2 meth 比例（`PerCpgAsm.cpp:308-310`），即 Δβ 形式的 HP-axis 差。
- 僅在 `hp1_total>=1 && hp2_total>=1` 才測（`PerCpgAsm.cpp:295`）。

**Step 1b — BH-FDR 校正 + 顯著計數**（`PerCpgAsm.cpp:317-338` + `bh_fdr_correction` 26-66）
- `bh_fdr_correction`：q_i = p_i · m / rank_i，再由右往左 enforce 單調 + clip 至 1.0（`PerCpgAsm.cpp:46-59`）；保留 NaN（未測 CpG）。
- 顯著閾值 **FDR < 0.05**（`PerCpgAsm.cpp:326`）。
- 輸出 4 個 Fisher 欄位：`fisher_n_tested`（315）/ `fisher_n_sig`（335）/ `fisher_frac_sig = n_sig/n_tested`（336）/ `fisher_max_neg_log_fdr`（337，`-log10(FDR)`，FDR clamp 下限 `1e-300` 防 log(0)，329）。

**Step 2 — NME（Normalized Methylation Entropy）**（`PerCpgAsm.cpp:340-361`, fn `compute_nme` 72-145）
- 每 read 把 binary CpG pattern 串成字串（`1`/`0`/`N`），數 unique pattern → Shannon entropy `H = -Σ p log2 p`（`PerCpgAsm.cpp:108-115`）。
- `H_max = min(log2(valid_reads), valid_cpg_count)`（`PerCpgAsm.cpp:136-138`）；NME = H / H_max（`144`）。
- 對 HP1/HP2 各算一次（`356-357`），`entropy_imbalance = |NME_HP1 − NME_HP2|`（`360`）。min_cpgs=2（`356-357`）。

**Step 3 — Epipolymorphism**（`PerCpgAsm.cpp:363-372`, fn `compute_epipolymorphism` 151-222）
- **4-CpG 滑動窗**（window_size default 4，`367-368`），每窗 epipoly = `1 − Σ p_i²`（pattern 頻率，`PerCpgAsm.cpp:206-213`），對所有窗取平均（`221`）；每窗需 ≥2 valid reads（`204`）。
- HP1/HP2 各算，`epipoly_delta = |epipoly_hp1 − epipoly_hp2|`（`371`）。

> 輸出共 10 欄（`PerCpgAsmResult`, header `PerCpgAsm.hpp:36-54`）：Fisher 4 + NME 3 + Epipoly 3。

### B. 引用註解第幾行寫什麼文獻

`src/core/PerCpgAsm.cpp:5-9`（檔頭 doc block）：
```cpp
 * Literature references:
 * - Per-CpG Fisher: DAMEfinder (Orjuela 2020), pycoMeth (Snajder 2023)   // line 6
 * - NME: CPEL (Jenkinson 2020, Nature Communications)                     // line 7
 * - Epipolymorphism: methclone (Li 2014), Metheor (Lee 2023)             // line 8
```
- **第 6 行**現寫 `DAMEfinder (Orjuela 2020), pycoMeth (Snajder 2023)`。
- git blame：line 6 由 commit `891e04b`（2026-06-10 "修 DAMEfinder 誤引"）重寫；其餘 5/7/8 仍是原始 `5b78f1a`（2026-04-19）。
- `git show 5b78f1a:src/core/PerCpgAsm.cpp` 第 6 行**原文 = `DAMEfinder (De Waele 2020)`** → 已被矯正為 `Orjuela 2020`。
- header `PerCpgAsm.hpp:8` 也把 Fisher 歸給 `DAMEfinder, pycoMeth`（無年份，未含誤引人名）。

> ✅ 結論：本輪實測 `PerCpgAsm.cpp:6` = **Orjuela 2020**（正確），舊 memory 記的 "De Waele 2020 誤引" 已修復（commit 891e04b）。

### C. ISM（本 C++）vs MSA（外部抽取工具）差別 + dup-bug 查核

**一句話差別**：ISM 是本 repo 的 C++ binary，`MethylationParser::parse_read` 直接從 BAM `MM/ML` tag 解析甲基化（`MethylationParser.cpp:15-199`），**只取 `C+m?`(5mC) 一種修飾、每個 CpG 位置一個機率值**；MSA（MethylSomaticAnalysis）是另一個外部抽取工具，會把 5mC 與 5hmC 抽成雙列，存在「雙列砍半 dup-bug」（屬 MSA，非本 C++）。

**dup-bug 在本 C++？— 不存在。 證據鏈：**

1. **Parser 只鎖定 5mC**：`MethylationParser.cpp:60` `const char* target = "C+m?";`；line 66 `deltas = parse_mm_tag(mm_str, "C+m?")`；`parse_mm_tag` default mod_code = `"C+m?"`（header `MethylationParser.hpp:76`）。找不到 `C+m?` 直接 `return calls`（空）（`MethylationParser.cpp:85-87`）。
2. **5hmC 只被「跳過」不產生 call**：`C+h?` block 只用來推進 `ml_offset`（數 delta 數，`MethylationParser.cpp:71-82`），**從不 emit MethylCall**。`grep` 證實 `C+h?`/`5hmC` 全 repo 僅出現在註解（`MethylationParser.cpp:39,41,104` 與 header 33），**無任何 active 解析路徑**。
3. **每 CpG 僅一機率**：forward 路徑 `calls.emplace_back(ref_pos_0based+1, prob)`（`MethylationParser.cpp:181`）、reverse `emplace_back(ref_pos_0based, prob)`（`142`），一個 CpG 位置一個 prob（`prob = ml_data[ml_offset+delta_idx]/255.0`，180/140）。
4. **Matrix 以「位置」為唯一鍵、不會雙列**：`MatrixBuilder::finalize` 收集 `call.ref_pos`（`MatrixBuilder.cpp:35,57`），`std::sort + std::unique` 去重 unique CpG 位置當 column（`MatrixBuilder.cpp:62-66`），`pos_to_col` 1 位置 → 1 column（`69-72`），填值 `row[col]=prob`（`95`）。**同一 CpG 不可能因 5mC/5hmC 變兩 column**。
5. **無 max-collapse 於兩 mod 之內**：C++ 路徑根本沒有第二個 mod 通道，故不需要、也沒有「5mC+5hmC 取 max 合併」步驟（該 max-collapse 是 Python/MSA 側 Level-1 處理；見 memory `reference_msa_vs_ism_tool` / `project_zar1l_brca2_asm_verification`）。

> ✅ 結論：**雙列砍半 dup-bug 不存在於本 C++（MethylationParser + MatrixBuilder）**；ISM C++ 天生是 **5mC-only、每 CpG 單機率單 column**，無從產生 5mC/5hmC 雙列。dup-bug 僅屬外部 MSA 抽取工具。

---

## L2 數字 / 常數 / 閾值（附來源）

| 項目 | 值 | 來源 |
|---|---|---|
| meth 二值化閾值 | `val > 0.5` | `PerCpgAsm.cpp:282,283,289,290` |
| Fisher 顯著閾值 | FDR `< 0.05` | `PerCpgAsm.cpp:326` |
| FisherConfig seed | 42 | `PerCpgAsm.cpp:267` |
| min_reads_per_group | default 5 | `PerCpgAsm.hpp:66` |
| Fisher 測試前提 | `n_reads>=2 && n_cpgs>=1`；`hp1_total>=1 && hp2_total>=1` | `PerCpgAsm.cpp:236,295` |
| NME min_cpgs | 2 | `PerCpgAsm.cpp:356-357` |
| Epipolymorphism window | 4 CpG | `PerCpgAsm.cpp:367-368` |
| Epipoly 每窗最少 reads | 2 | `PerCpgAsm.cpp:204` |
| FDR -log10 下限 clamp | `1e-300` | `PerCpgAsm.cpp:329` |
| ML 機率換算 | `ML[i] / 255.0` | `MethylationParser.cpp:140,180` |
| BH q-value 公式 | `p_i * m / (i+1)`，右→左單調 | `PerCpgAsm.cpp:48,52-58` |
| 輸出欄位數 | 10（Fisher 4 + NME 3 + Epipoly 3） | `PerCpgAsm.hpp:36-54` |

---

## L3 溯源 / 版本

- Fisher 引用矯正：commit `891e04b`（2026-06-10, "修 DAMEfinder 誤引"）— line 6 `De Waele 2020` → `Orjuela 2020`。
- 模組初版：commit `5b78f1a`（2026-04-19, "Phase 2 核心 NormalBaseline + PerCpgAsm + Signed HP Delta"）。
- 檔案：`src/core/PerCpgAsm.cpp`(12268B, mtime Jun 10) / `include/core/PerCpgAsm.hpp` / `src/core/MethylationParser.cpp` / `src/core/MatrixBuilder.cpp` / headers。

## uncertain（讀不確定處）

- `FisherExact::test_2x2` 的 p-value 是雙尾/單尾、是否含 mid-p 修正 — 未讀 `FisherExact.cpp`（僅見 header `include/core/FisherExact.hpp` 被 include，`PerCpgAsm.hpp:23`）；`seed=42` 暗示某種隨機/Monte-Carlo 大表路徑，但未驗證。
- `meth_mat.binary_matrix` 的 0/1/-1 是否同樣以 0.5 為界由上游 build — NME/Epipoly 用 binary_matrix（`PerCpgAsm.cpp:348`），Fisher 用 raw_matrix；binary 化規則在 `MethylationMatrix::build`（宣告 header `MethylationMatrix.hpp:38`），本輪未讀其 .cpp 實作。
- MSA 工具本體（5mC+5hmC 雙列砍半 dup-bug 機制細節）不在本 repo C++；僅憑 memory `reference_msa_vs_ism_tool` 認定屬外部工具，本輪未直接讀 MSA 源碼佐證其 dup-bug。
