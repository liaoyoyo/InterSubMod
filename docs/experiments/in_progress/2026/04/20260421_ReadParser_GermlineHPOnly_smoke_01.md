---
status: in_progress
date: 2026-04-21
type: smoke_test
phase: 0
plan: /bip7_disk/liaoyoyo2001/.claude/plans/streamed-spinning-wilkinson.md
hypothesis_id: H-HP-ONLY-01
priority: P0
pipeline_track: TO
---

# ReadParser Germline-HP-Only — Phase 0 Smoke Test

## 摘要

**目的**：驗證新加入的 `--germline-hp-only` flag 與 `NHP_Somatic11/21/33` audit 欄位在 chr19 TP 子集上行為正確，確認：
1. Flag off 時 output 保持原行為（pass-through）
2. Flag on 時 somatic HP tags (HP:i:11/21/33 → "1-1"/"2-1"/"3") 被 demote 為 "0"
3. `NHP_Somatic*` 欄位不受 flag 影響（audit 獨立性）
4. 下游 HP-dependent 特徵（NHP3 / HP1FamilyN / HP2FamilyN / HPFineNGroups）正確反映 demotion

**結論**：**所有檢核通過** — Phase 0 驗證 PASS，可進入 Phase 0-8 commit，並視需求進入 Phase 1 全 HCC1395 比對。

---

## 測試設定

| 參數 | 值 |
|------|---|
| BAM | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam`（V3-Fixed haplotag TO BAM）|
| Normal BAM | `/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam` |
| VCF | `clairsto_v3fixed_work/clairsto_tp.vcf.gz` 之 **chr19 子集**（615 sites）|
| Reference | `GRCh38_no_alt_analysis_set.fasta` |
| Window | 5000 bp |
| Threads | 32 |
| Metric | BERNOULLI |
| Binary | `/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod`（Phase 0 local build）|
| 兩次執行差異 | 僅 `--germline-hp-only` flag 是否啟用 |

執行時間：off 20.7s、on ~21s；Memory ~27 MB；615/615 regions 全部成功。

---

## 驗證結果

### 1. Output schema 一致性

- 兩 run 皆 615 rows × 118 columns
- Column 順序一致（新增三欄在 tail，不影響既有欄位）
- `NHP_Somatic11/21/33` 位於 col 115/116/117（`SuggestFilter` 之後）

### 2. Audit 欄位獨立性（核心正確性）

| 欄位 | off sum | on sum | identical |
|------|---------|--------|-----------|
| NHP_Somatic11 | 4,894 | 4,894 | ✅ True |
| NHP_Somatic21 | 8,630 | 8,630 | ✅ True |
| NHP_Somatic33 | 2,498 | 2,498 | ✅ True |

**三欄在兩 run 間 per-row 完全相同**（audit 計數 key 為 `hp_tag_raw`，與 flag 無關）。

### 3. Demotion 數學守恆（critical）

期望：flag on 時 NHP0 增量 = NHP_Somatic11 + NHP_Somatic21 + NHP_Somatic33（+ 可能 NHP3→0 部分）

| 欄位 | off sum | on sum | delta |
|------|---------|--------|-------|
| NHP3 | 2,498 | 0 | **-2,498** |
| NHP0 | 13,985 | 30,007 | **+16,022** |
| HP1FamilyN | 26,462 | 21,568 | **-4,894** |
| HP2FamilyN | 33,134 | 24,504 | **-8,630** |

- NHP3 drop 2,498 = NHP_Somatic33 ✅
- HP1FamilyN drop 4,894 = NHP_Somatic11 ✅
- HP2FamilyN drop 8,630 = NHP_Somatic21 ✅
- NHP0 gain 16,022 = 4,894 + 8,630 + 2,498 ✅ **完全守恆**

### 4. HPFineNGroups 分佈收斂

| NGroups | off count | on count |
|---------|-----------|----------|
| 0 | 2 | 2 |
| 1 | 3 | 10 |
| 2 | 110 | **603** |
| 3 | 394 | 0 |
| 4 | 106 | 0 |

Flag off 時 5 類 NGroups（0-4）；flag on 時收斂至 3 類（0-2）— 符合「11/21/33 三個 somatic labels 移除後，精細分組只剩 germline HP1/HP2 + unphased」的預期。

### 5. HP_Ratio 改變（非核心但有參考性）

| 統計量 | off | on | delta |
|--------|-----|-----|-------|
| median | 0.4615 | 0.4588 | -0.003 |
| mean | 0.4508 | 0.4680 | +0.017 |
| n_nonzero | 612 | 612 | 0 |

**與 plan 預期差異**：plan 期待 TP HP_Ratio median 從 0.836 → 0.55-0.65。本次 smoke 為 **chr19 615-site 子集**（混合 LOH / balanced），與 plan 所引的 6,485 個全基因 balanced 位點組成不同，故 median 基線已是 0.46。本欄位差異**不作為 Phase 0 pass/fail 標準**，留待 Phase 1 全 HCC1395 驗證。

---

## Phase 0 Gate 判定

| 檢核項 | 標準 | 結果 |
|-------|------|------|
| 編譯（make -j）| 無 error，warning 僅 pre-existing | ✅ pass |
| ctest | 215/215 passed（baseline 203 + 新增 12 test_read_parser）| ✅ pass |
| Output schema | flag off 時新欄位在 tail（pandas-safe）| ✅ pass |
| Audit 獨立性 | NHP_Somatic* off == on | ✅ pass |
| Demotion 守恆 | Σ(somatic) = Δ(NHP0) | ✅ exact |
| 下游效應 | NHP3/HP1FamilyN/HP2FamilyN/HPFineNGroups 符合預期 | ✅ pass |

**結論**：Phase 0 全部通過，可進入 commit 階段。

---

## 後續

- **Phase 0-8**（緊接）：commit 本次修改（Config/ArgParser/ReadParser/DataStructs/RegionProcessor/CMakeLists/tests + 本報告）
- **Phase 1**（條件性）：HCC1395 全基因 TO 模式 flag on/off 比對；計算 within_dom_alt_frac / HP_Ratio / HPFineNGroups AUC 差異
- **Phase 2**（條件性）：7 樣本 × 2 模式全量重跑，依 plan R5 與 CovM baseline rerun 合併規劃

---

## 檔案清單（Phase 0 改動）

- `include/core/Config.hpp` — `germline_hp_only` flag
- `include/utils/ArgParser.hpp` — `--germline-hp-only` CLI
- `include/core/ReadParser.hpp` — `ReadFilterConfig::germline_hp_only`
- `include/core/DataStructs.hpp` — `ReadInfo::hp_tag_raw`
- `include/core/RegionProcessor.hpp` — `RegionResult::n_hp_somatic_{11,21,33}`
- `src/core/ReadParser.cpp` — switch 改寫（raw 保留 + conditional demote）
- `src/core/RegionProcessor.cpp` — filter_config 傳遞 + audit 計數 + TSV 擴充
- `tests/test_read_parser.cpp` —（新增）12 test cases
- `CMakeLists.txt` — 加入 test_read_parser.cpp

**驗證資料位置**：`/tmp/ism_hp_fix_smoke/{off,on}/significance_summary.csv`（臨時，commit 後可刪除）
