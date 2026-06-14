<!--
建立時間: 2026-06-13
任務: 既有方法稽核（甲基段）萃取 — 從 2026-06-10 全碼方法學稽核擷取所有甲基相關 finding
來源稽核: docs/methodology/20260610_code_methodology_detail_audit_01.md (HEAD=bbcc7a9, 2026-06-10)
驗證: 本萃取每條 finding 都回 src/include/tests/scripts 實檔對證行號 + git log 核狀態變化
紀律: 只讀實際檔案、引行號、引原始碼片段；不憑記憶；數字附 檔:行
-->

# ISM 甲基讀取源碼稽核萃取（甲基段）

> 對象稽核 = `InterSubMod/docs/methodology/20260610_code_methodology_detail_audit_01.md`（12-reviewer workflow wf_ea727f9a-ef5，HEAD `bbcc7a9` @ 2026-06-10）。
> 本萃取**逐條回實檔對證**，並標出稽核當時狀態 vs 現行 repo 狀態（HEAD `4ef557b`）的差異。

---

## L0 — 一句話

稽核對「甲基讀取/MM-ML/5mC-5hmC/PerCpgAsm」共給出 6 條核心 finding。**稽核當時**：2 個 🔴 必修（TESTS-1 解析零測試 / PYTOOL-1 反捏造洞）、1 個 ⭐ 統計校準（FISHER-1 over-dispersion）、3 個 ⭐/◽ 口徑 caveat（ASM-2 軸別 / ASM-3 5mC-only / ASM-4 二值化口徑）。**現行 repo（本日實測）**：兩個 🔴 必修**已落地修復**（post-audit commit），De Waele→Orjuela 引用**早已修**（pre-audit），其餘口徑/Fisher 議題**仍在**。

---

## 1. TESTS-1 — 甲基化解析「最易錯的碼」零測試（稽核: 🔴 必修, CONFIRMED HIGH）

**稽核原文（audit L2 T1 表 + L0 + 附錄 B序1）**：
> `MethylationParser.cpp`(290 行) 的 MM/ML 解析、反股 **倒序迭代** CpG 座標（碼內自標 "CRITICAL" 103-157）、5mC/5hmC 分支 **完全無單元測試**；現有測試只覆蓋「filter-to-zero」負路徑。對抗驗證三路嘗試反駁全失敗 → **CONFIRMED HIGH**。
> grep `MethylationParser|MatrixBuilder` 在 `tests/*` = **0 命中**。
> 為何最該先修：此碼是資料入口，錯了會**靜默污染每一個甲基化 call → 每一個 Δβ / 距離 / 分群**，且無號誌。

**實檔對證**：
- `src/core/MethylationParser.cpp` = 290 行（實測 `wc -l` 確認）。
- 反股 CRITICAL 倒序迭代註解確實在碼內 — `MethylationParser.cpp:107-110`：「CRITICAL: MM tags list deltas in 5'->3' order of the ORIGINAL read... for reverse strand, we must iterate BAM SEQ BACKWARDS (Len-1 -> 0)」；對應 loop `MethylationParser.cpp:126` `for (int seq_idx = seq_len - 1; seq_idx >= 0; seq_idx--)`。

**狀態變化（🔴 已修，POST-AUDIT）**：
- 稽核 grep 命中數 0 為**稽核當時（HEAD bbcc7a9）真值**：`git cat-file -e bbcc7a9:tests/test_methylation_parser.cpp` → **ABSENT at bbcc7a9**。
- **現行 repo 已新增** `tests/test_methylation_parser.cpp`（6789 bytes），commit `6593f96` "test(methylation-parser): add MM/ML golden characterization tests (audit TESTS-1)"，日期 2026-06-10 13:51（在稽核 HEAD bbcc7a9 04:50 **之後** ~9hr；`git merge-base --is-ancestor bbcc7a9 6593f96` = post-audit）。
- ⇒ **稽核狀態 = 必修；現行狀態 = 已落地修復**（依稽核附錄 B 建議「補 MM/ML golden test 正/反股+多 CpG+邊界」）。

---

## 2. ASM-3 — C++ 為 5mC-only（稽核: ◽ 口徑 caveat, T4/附錄 A）

**稽核原文（audit L2 T4 ASM-3 + 附錄 A 第5列）**：
> `MethylationParser.cpp:41-90`：C++ **5mC-only**（只解析 `C+m?`，`C+h?` 僅數 delta 推進 ml_offset）→ signed Δβ 非 any-mod max-collapse。
> **memory 說的「5mC+5hmC 雙列砍半 dup-bug」其實在 MSA Python 抽取層，不在 C++ binary**。
> 附錄 A 裁決：「✅ 釐清：C++ 層無此 bug；真議題是 C++ **忽略 5hmC**（口徑 caveat）」。

**實檔對證**：
- `MethylationParser.cpp:60` `const char* target = "C+m?";` — target 鎖死 5mC。
- `MethylationParser.cpp:65-69`：`if (is_target) { deltas = parse_mm_tag(mm_str, "C+m?"); found_target = true; break; }` — 只取 `C+m?` block。
- `MethylationParser.cpp:71-82`：非 target block（含 `C+h?`）僅**數 delta count 推進 `ml_offset`**，不產生 call（`ml_offset += delta_count;` line 81）。
- `MethylationParser.cpp:85-87`：`if (!found_target) return calls;`（找不到 `C+m?` 整條 read 丟）。
- 註解 `MethylationParser.cpp:39`：「MM format: "C+h?,<deltas>;C+m?,<deltas>;"」確認 5hmC block 存在但被跳過。

**狀態（◽ caveat 仍在，非 bug）**：C++ binary 至今仍是 5mC-only；稽核明確將「5mC+5hmC 雙列砍半」歸給 MSA Python 抽取層，C++ 層無該 dup-bug。此為口徑 caveat（與論文要的 any-mod max-collapse 有口徑差），非待修 bug。

---

## 3. ASM-2 / REGION-5 — ASM 軸別：germline/ALLELE 軸（1-1 折進 family 1）（稽核: ⭐ UNCERTAIN→MEDIUM）

**稽核原文（audit L2 T4 ASM-2 表）**：
> `RegionProcessor.cpp:903-906`, `PerCpgAsm.cpp:244-249`：預設把 `{"1","HP1","1-1"}→0`、`{"2","HP2","2-1"}→1` 合併 → signed Δβ 與 per-CpG Fisher 量的是 **HP1-family vs HP2-family = germline/ALLELE 軸**（somatic 子 haplotype `1-1` 被折進 family `1`）。**非 somatic-controlled HP 軸**（HP1 vs HP1-1）。
> 與 memory `feedback_asm_allele_axis_baseline_confound` 一致——碼忠實實作了會被 baseline allelic methylation confound 的軸，但輸出欄名（`tumor_hp_signed_delta`）未明示軸別 → 下游易誤用。

**實檔對證**（`PerCpgAsm.cpp:244-250`）：
```cpp
244  if (hp == "1" || hp == "HP1" || hp == "1-1") {
245      hp_labels[i] = 0;
246      hp1_indices.push_back(i);
247  } else if (hp == "2" || hp == "HP2" || hp == "2-1") {
248      hp_labels[i] = 1;
249      hp2_indices.push_back(i);
250  }
```
→ somatic 子 haplotype `1-1` 與 germline `1`/`HP1` 同併入 label 0；`2-1` 與 `2`/`HP2` 併入 label 1。確認軸別 = HP1-family vs HP2-family（germline/ALLELE），非 somatic-controlled HP1-vs-HP1-1。

**狀態（⭐ caveat/口徑治理問題，仍在）**：稽核裁 UNCERTAIN→MEDIUM；建議 = C++ 輸出欄改名/加註明示軸別（germline-family vs somatic-HP），非單點 bug。現行碼未改。

---

## 4. FISHER-1 — per-CpG Fisher 忽略 read over-dispersion（稽核: ⭐ CONFIRMED→MEDIUM）

**稽核原文（audit L1 T3 + L2 T3 FISHER-1）**：
> `PerCpgAsm.cpp:274-313`：每個 read×CpG binary call 當獨立 Bernoulli 餵 Fisher exact，**無 over-dispersion 處理**（全碼無 beta-binomial/quasi/GEE）。
> 合成示意（ILLUSTRATIVE，非真資料）：ρ=0.3 read 內相依下，NULL 的 Fisher 名目 5% FP 膨脹到 **53–68%（~17–20×）**。真量級需 `methylation.csv` 跑 read-level permutation 對照（needs_rerun）。
> **與專案 memory 自己的 Fisher→beta-binomial+shrinkage 建議一致**。
> 共同主旨：讓「顯著」偏多、p 偏小；對描述性 ASM 存在性率（如 Fisher_Frac_Sig）影響量級可觀，**對 TP/FP 方向、對已 concluded 的 NEGATIVE 結論影響小**。

**實檔對證**（`PerCpgAsm.cpp:274-313`）：
- 逐 read 二值計數，每個 read×CpG cell 當獨立 Bernoulli：`PerCpgAsm.cpp:279-290`（`for (int idx : hp1_indices){ ... if (val>0.5) hp1_meth++; else hp1_unmeth++; }`，HP2 同）。
- 直接餵 2×2 Fisher：`PerCpgAsm.cpp:298-305`（`ContingencyTable2x2 table; table.a=hp1_meth; ...; FisherResult fr = fisher.test_2x2(table);`）。
- 全碼無 beta-binomial/quasi/GEE（稽核宣稱；本萃取信任稽核 grep，未獨立重 grep 全碼 — 標 uncertain）。
- BH-FDR **僅** per-region 跨 CpG：`PerCpgAsm.cpp:317-318`（`if (n_tested > 0){ // BH-FDR correction`）→ 確認 region 內有多重檢定校正，但對「read 內相依」無校正（這是 over-dispersion 議題，非 multiple-testing 議題）。

**狀態（⭐ 仍在，未改）**：稽核建議改 beta-binomial / read-level permutation，或維持 Fisher 但輸出明標 "anti-conservative, independence-assumed"。屬 Hard Gate C++ 改動（需 /methodology-audit→/cpp-change），稽核未自行改。現行碼仍為獨立 Bernoulli Fisher。
> ⚠ 53–68% FP 膨脹數字稽核**自標 ILLUSTRATIVE（合成示意，非真資料）**，真量級 needs_rerun。引用此數字必帶此 caveat。

---

## 5. ASM-4 — per-CpG Fisher 二值化 `raw>0.5` 與 binary_matrix 0.2/0.8 口徑不一（稽核: ◽ 口徑）

**稽核原文（audit L2 T4 ASM-4）**：
> `PerCpgAsm.cpp`：per-CpG Fisher 二值化 `raw>0.5` 硬切，與 binary_matrix 的 0.2/0.8（中間=missing）**口徑不一**。

**實檔對證**：
- `PerCpgAsm.cpp:282` `if (val > 0.5) hp1_meth++; else hp1_unmeth++;`、`PerCpgAsm.cpp:288` 同（HP2）→ 確認 0.5 硬切，無「中間視 missing」帶。
- NaN 在二值化前被跳過：`PerCpgAsm.cpp:281` `if (std::isnan(val)) continue;`（缺值 skip，但 0.2~0.8 的「弱訊號」cell 仍被硬分到 meth/unmeth，不像 binary_matrix 設 0.2/0.8 灰帶）。

**狀態（◽ caveat 仍在）**：口徑不一致，非 bug；稽核歸 T4 口徑治理。

---

## 6. PYTOOL-1 — 反捏造工具 fill_report.py 對 null/NaN 失效（稽核: 🔴 必修, CONFIRMED HIGH）

**稽核原文（audit L1 T2 + L2 T2 PYTOOL-1）**：
> `scripts/fill_report.py:38,54`：`resolve({'auc':None},'auc')=(True,None)` → `found=True` 不進 `missing`；`fmt(None)='None'`、`fmt(nan)='nan'`。
> 🔴 §13-A「template+data 注入缺 key 必 refuse」**對 null 值失效**——數據是 null 時報告靜默出現字串 `None`/`nan` 而非拒絕。
> 修法：`fill_report` 把 `None`/NaN 視同 missing 一併 refuse。

**實檔對證 + 狀態變化（🔴 已修，POST-AUDIT）**：
- 現行 `scripts/fill_report.py:37-44` 已新增 `is_nullish(value)` helper，docstring 明寫「These must be refused like a missing key (anti-fabrication)... Note 0 and False are NOT nullish」。
- 已 wire 進 missing 邏輯：`scripts/fill_report.py:94` `if not found or is_nullish(val): missing.append(key)`。
- commit `eed4300` "fix(fill_report): refuse null/NaN values instead of rendering 'None'/'nan' (audit PYTOOL-1)"，日期 2026-06-10 13:51（post-audit，`git merge-base --is-ancestor bbcc7a9 eed4300` = 在 bbcc7a9 之後）。
- ⇒ **稽核狀態 = 必修；現行狀態 = 已落地修復**。
> 註：fill_report.py 非甲基化解析碼，但屬「餵甲基數字進報告的反捏造防線」，與甲基段間接相關；列此以完整呈現稽核兩個 🔴。

---

## 7. De Waele → Orjuela 引用（稽核: ✅ 已修, 附錄 A — PRE-AUDIT 已修）

**稽核原文（audit 附錄 A 第1列）**：
> memory 宣稱 `PerCpgAsm.cpp:6` 誤引「De Waele 2020」應為 Orjuela。
> 現行碼實況：`= DAMEfinder (Orjuela 2020)`；grep「De Waele」**0 命中**。
> 裁決：✅ **已修**（commit 891e04b），非 open finding。

**實檔對證**：
- `src/core/PerCpgAsm.cpp:6` 現行 = `* - Per-CpG Fisher: DAMEfinder (Orjuela 2020), pycoMeth (Snajder 2023)`（實測）。
- `grep -rn "De Waele" src/ include/` = **0 命中**（exit 1）。
- commit `891e04b` "fix(citation): DAMEfinder 作者誤引（一手 PMC7268773 確認 = Orjuela et al. 2020）— src/core/PerCpgAsm.cpp:6 "De Waele 2020" → "Orjuela 2020""，日期 2026-06-10 01:11，`git merge-base --is-ancestor 891e04b bbcc7a9` = **YES（pre-audit）**。
- ⇒ 此修在稽核**之前**完成；稽核如實將其列為 memory-stale「已修」非 open finding。memory 條目 `project_ism_vs_external_methylation_tools_comparison.md` 標的「🔴 PerCpgAsm.cpp:6 誤引」屬該修之前的舊紀錄。

---

## 附：稽核時間軸（pre-audit vs post-audit）

| commit | 日期 | 性質 | 對甲基段 finding |
|---|---|---|---|
| `891e04b` | 2026-06-10 01:11 | **pre-audit** (ancestor of bbcc7a9) | De Waele→Orjuela 引用修（附錄 A 已記為已修）|
| `bbcc7a9` | 2026-06-10 04:50 | **= 稽核 HEAD** | 稽核快照時點 |
| `6593f96` | 2026-06-10 13:51 | **post-audit** | TESTS-1 修：新增 test_methylation_parser.cpp（MM/ML golden）|
| `eed4300` | 2026-06-10 13:51 | **post-audit** | PYTOOL-1 修：fill_report 對 null/NaN refuse |

> ⇒ 稽核兩個 🔴 必修（TESTS-1 / PYTOOL-1）均在稽核**當日稍晚**落地修復；引用 bug 在稽核**前**已修。仍未動的甲基段議題 = FISHER-1（over-dispersion，需 Hard Gate C++ 改）、ASM-2（軸別欄名）、ASM-3（5mC-only）、ASM-4（二值化口徑）—— 皆口徑治理/統計校準，稽核明示**不翻轉任何已發表 NEGATIVE 結論方向**。

## uncertain / 讀不確定處
- 「全碼無 beta-binomial/quasi/GEE」為稽核宣稱；本萃取未獨立 grep 全 23 個 src/core/*.cpp 重核此 negative claim（只確認 PerCpgAsm.cpp 該段是獨立 Bernoulli Fisher）。標 uncertain。
- ASM-2 引的 `RegionProcessor.cpp:903-906` 行號本萃取未開檔對證（只對證了 PerCpgAsm.cpp:244-250 的同軸別合併邏輯）。RegionProcessor 那側標 uncertain。
- FISHER-1 的 53–68% FP 膨脹為稽核自標 ILLUSTRATIVE 合成示意，非真資料 effect；真量級 needs_rerun。
- `tests/test_methylation_parser.cpp` 內容是否真覆蓋「正股+反股+多CpG+邊界」四情境，本萃取只確認檔案存在+commit message 宣稱「MM/ML golden characterization tests」，未逐 test case 對證覆蓋面。標 uncertain。
