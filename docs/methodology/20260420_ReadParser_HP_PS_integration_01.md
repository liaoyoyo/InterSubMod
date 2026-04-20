---
title: ReadParser HP + PS tag 整合方法論審查（P0-C）
date: 2026-04-20
status: audit_draft（待選方案 → /cpp-change 實作）
owner: InterSubMod Research
related:
  - docs/reports/audit/decisions/CHECKLIST.md#P0-C
  - include/core/ReadParser.hpp
  - src/core/ReadParser.cpp
  - src/core/LabelTest.cpp
  - include/core/DataStructs.hpp
scope: ReadParser 補齊 PS（phase-set）tag 解析；評估對 29 個 HP-dependent 特徵的影響
decision_required: 選 A / B / C
---

# ReadParser HP + PS Tag 整合方法論審查

## 一、問題現況

### 1.1 實際行為 vs 文件宣稱

**hpp 註解（include/core/ReadParser.hpp:40）**：
```cpp
* - Parsing phasing tags (HP, PS)
```

**實作現況（src/core/ReadParser.cpp:121-142）**：
```cpp
// Extract HP tag (haplotype)
uint8_t* hp_aux = bam_aux_get(b, "HP");
if (hp_aux) { ... }
// NOTE: 無 bam_aux_get(b, "PS") 呼叫
```

**DataStructs.hpp:32**：
```cpp
std::string hp_tag;  ///< Haplotype tag (HP): "1", "2", "1-1", "2-1", "unphase", etc.
// NOTE: 無 ps_tag 欄位
```

### 1.2 PS tag 的生物學意義

- **HP (Haplotype)**：read 屬於哪一條 haplotype（1/2/3）
- **PS (Phase-Set)**：read 所屬的 phasing block ID（integer，通常為 block 起點 position）
- **關鍵關係**：不同 PS block 的 HP=1 彼此**不可比**，因為各 block 的 HP1/HP2 指派是獨立產生的
  - block A 的 HP1 可能對應 block B 的 HP2（方向反轉）
  - 跨 block 使用同一 HP 標籤會混淆真實的 allele lineage

### 1.3 為何當前實作仍可運作

在 ISM 分析中，每個 Region（window_size=5000bp）通常落在**單一 PS block** 內：
- LongPhase PS block 在 5kHz ONT 下通常 >10kb
- 區域內 HP1/HP2 一致性高 → HP-only 處理不致於嚴重錯分

但存在以下**潛在失效場景**：

| 場景 | 機率 | 影響 |
|------|------|------|
| Region 橫跨 PS 邊界（e.g., low-coverage 區或 switch error 斷點） | ~5-10%（視 phase block 長度而定） | HP 混合 → ISM 甲基化異質性假陽性 |
| Region 內存在 phase-switch（LongPhase 重新起算 block） | ~2-5% | 同 HP=1 label 實為跨 block，ISM 跨 region 特徵失真 |
| 跨 region 比較（如 O11 heterogeneity、O13 cross-region） | 所有 cross-region 分析 | HP 無法保證 lineage 一致 |

### 1.4 證據：29 HP-dependent 特徵消費鏈

Grep `hp_tag` 跨 9 個 source file（34 處引用）：

| File | 使用處 | 角色 |
|------|-------|------|
| `src/core/ReadParser.cpp` | 8 | HP 抽取（生產者） |
| `src/core/LabelTest.cpp` | 5 | HP→group 映射（0/1/2/3） |
| `src/core/FisherExact.cpp` | 13 | 依 HP 分組進行 Fisher's test |
| `src/core/RegionProcessor.cpp` | 3 | 主 pipeline（normal read 強制 hp="0"） |
| `src/io/RegionWriter.cpp` | 1 | 輸出 read-level TSV |
| `src/core/PerCpgAsm.cpp` | 1 | ASM 分析 |
| `src/core/LocalTest.cpp` | 1 | 局部統計 |
| `src/core/SignificanceAnalyzer.cpp` | 1 | 顯著性 |
| `src/test/test_phase1_2.cpp` | 1 | 單元測試 |

LabelTest 的典型 HP→group 邏輯（src/core/LabelTest.cpp:227-230）：
```cpp
const std::string& hp = full_labels[i].hp_tag;
// HP1 variants -> group 0 (includes HP1 from any sample: 1, 1-1)
if (hp == "1" || hp == "HP1" || hp == "1-1") { group = 0; }
// ...
```

**關鍵問題**：若 read 來自不同 PS block 的 HP=1，仍被映射到 group 0。

---

## 二、影響量化估計

### 2.1 依存此問題的研究結論

| Audit Card | 結論 | 依賴 HP identity | 若 PS 整合會改變？ |
|-----------|------|-----------------|------------------|
| C08 Read-level LOSO AUC=0.721 | POSITIVE（但 FP removal=0%） | 高（per-read HP-specific features） | **可能**：跨 PS read 混入 HP group 會污染特徵空間 |
| C11 HPFineNGroups N≥4 | POSITIVE（TP 89.12%→92.81%） | 高（group 數＝PS 獨立 HP）| **可能**：若同 region 內無跨 PS，影響微；跨 region aggregation 會改變 |
| C15 Option C 雙路 HP-free combo | NEGATIVE（AUC 0.564） | 反向依賴（HP-free 驗證）| **無影響** |
| C16 HPFineNGroups residualized AUC=0.617 | POSITIVE（>0.58） | 高 | **可能**：bootstrap CI 需重跑 |
| C22 Zone-Aware | CONDITIONAL | 中（zone 定義用 HP） | 輕微 |
| 29 HP-dependent features 全體 | 混合 POSITIVE/NEGATIVE | 全依賴 | 全部需重新驗證 |

**預期影響範圍**：若 PS 整合後多數 region 內仍為單一 PS block（~90-95%），大多數結論應**質化不變**；但 CI 邊界 (AUC±0.01) 可能震盪，需 within-group OLS 重跑確認。

### 2.2 計算成本

| 項目 | 成本 |
|------|------|
| C++ 改動 | ReadParser.cpp +15 行 / DataStructs.hpp +1 field |
| 下游改動（依選項） | 選項 A 全面 +5-20 行/file × 7 files；選項 B 零；選項 C 零 |
| 重編譯 | <30 秒 |
| 7 樣本重跑 | 需搭配 P0-A rerun 窗口 |
| 驗證 within-group OLS | P0-B 範圍內同時完成 |

---

## 三、候選方案

### 選項 A — 複合 identity key（HP:PS composite）

**做法**：ReadParser 產生 `hp_tag = "PS:HP"`（e.g., "1234567:1"），下游所有 string comparison 需更新。

**優點**：
- 語義正確：同 HP 不同 PS 被視為不同 group
- 最嚴格 phasing 保證

**缺點**：
- **破壞性極高**：LabelTest/FisherExact/RegionWriter 全需改 string 比對邏輯
- HPFineNGroups 定義變動（原本 N=4 可能變 N=50+，解釋性崩潰）
- 與歷史結論不可比：C11/C16 重跑後數值會劇烈變動
- 單元測試全部失效

**代價**：~3-5 工作日 + 歷史結論重新解讀

**風險**：🔴 高 — 可能打破現有 POSITIVE 結論，需全量 audit card 重審

### 選項 B — 平行欄位（HP + PS 各自獨立）【推薦】

**做法**：
1. DataStructs.hpp 新增 `std::string ps_tag = "0"`
2. ReadParser.cpp 新增 `bam_aux_get(b, "PS")` 解析（類似 HP 邏輯）
3. ReadInfo.ps_tag 輸出到 read-level TSV（P2-E 源碼分析需求）
4. **下游消費者保持不變**：HP→group 邏輯維持現狀
5. 未來若特定分析需要 PS-aware，可選擇性接入（opt-in）

**優點**：
- **零破壞性**：所有現有結論不變
- 為未來 PS-aware 分析預留欄位（e.g., 跨 region lineage tracking）
- 與 P2-E (Self-phasing algorithmic) 源碼分析契合，提供更完整的 read-level metadata
- 單元測試只需補「ps_tag 被正確解析」一項

**缺點**：
- 不立即解決 HP 跨 PS 混合問題
- PS-aware 效益需後續手動接入下游

**代價**：~0.5 工作日（C++ 改動 + 單元測試）

**風險**：🟢 低 — 純 additive，無破壞性

### 選項 C — Canonicalize（PS 缺失時降為 unphased）

**做法**：
1. ReadParser 讀 PS + HP；**若 PS 不存在或無效，強制 hp_tag = "0"（unphased）**
2. DataStructs 維持不變，無新增欄位
3. 下游無感知

**優點**：
- 最小改動
- 修正「無 PS block 仍被當成 phased」的 edge case

**缺點**：
- LongPhase-s/to 通常 PS 與 HP 一起產生；PS-missing 場景罕見
- 無法解決「跨 PS block 的 HP 混合」核心問題
- 收益很小

**代價**：~2 小時

**風險**：🟡 中 — 可能意外將部分 phased read 降為 unphased，影響 HP-dependent 特徵統計

---

## 四、推薦選項：B（平行欄位）

### 4.1 推薦理由

1. **零下游破壞** — 不影響 22 audit cards 的既有結論，避免回溯風暴
2. **為未來擴展預留** — PS-aware 分析可在 P2-E / Phase 2 A+D 階段按需接入
3. **與 P2-E 契合** — Self-phasing algorithmic 源碼分析需要 read-level PS 資訊才能追蹤 phase-switch 事件
4. **符合 read-level epigenetic context 主軸** — ISM characterization 重心是細緻化 read metadata，而非強化 filter
5. **最低風險 / 高邊際效益** — 一次 C++ 改動解鎖未來多項分析

### 4.2 實作規格

**C++ 改動**：

`include/core/DataStructs.hpp` ReadInfo struct：
```cpp
std::string hp_tag;  ///< Haplotype tag (HP): "1", "2", "1-1", "2-1", "0"=unphased
std::string ps_tag;  ///< Phase-set tag (PS): integer block ID as string, "0"=no PS
```

`src/core/ReadParser.cpp` parse() 函式（HP 後新增）：
```cpp
// Extract PS tag (phase-set block ID)
info.ps_tag = "0";  // Default: no phase set
uint8_t* ps_aux = bam_aux_get(b, "PS");
if (ps_aux) {
    char type = ps_aux[0];
    if (type == 'Z' || type == 'H') {
        info.ps_tag = bam_aux2Z(ps_aux);
    } else if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I') {
        int ps_int = bam_aux2i(ps_aux);
        info.ps_tag = std::to_string(ps_int);
    }
}
```

`src/io/RegionWriter.cpp`：TSV 輸出新增 PS column（位於 HP 欄之後）。

### 4.3 驗證標準（Step → Verify）

1. **編譯** → 驗證：`cd build && make -j$(nproc)` 成功，無 warning
2. **單元測試** → 驗證：`./tests/test_read_parser` 新增 `ParsesPSTagIntegerAndString` case 通過
3. **Integration 測試** → 驗證：對 HCC1395 單 region 輸出 reads.tsv，grep `PS` column 有整數值（非全 "0"）
4. **HP 原行為不變** → 驗證：diff 比對新舊 binary 對同一 region 的 HP label 分布（應 100% 一致）
5. **PS 涵蓋率** → 驗證：7 樣本取樣 reads，PS≠"0" 比例應 >90%（與 HP≠"0" 比例相當）

### 4.4 不做的事（避免 scope creep）

- ❌ 不修改 LabelTest / FisherExact / HPFineNGroups 的 HP→group 邏輯
- ❌ 不產生 composite "PS:HP" 標識符
- ❌ 不重跑 C08/C11/C16 驗證（選項 A 才會打破，選項 B 不會）
- ❌ 不新增 PS-based filter 或 zone 定義（留給未來獨立計畫）

---

## 五、依賴關係

### 前置依賴
- **無強依賴**。可與 P0-A rerun 並行（但建議等 rerun 完成後再 `/cpp-change`，避免 binary 版本混淆）

### 被阻塞項
- **P2-E (Self-phasing algorithmic 源碼分析)** — 需要 ps_tag 在 read-level TSV 以追蹤 phase-switch 事件
- **Phase 2 A+D Sample ASM** — normal read 的 PS 可輔助區分 germline vs somatic phased reads

### 建議執行順序

```
P0-A 7樣本 rerun 完成（current）
  ↓
P0-C /cpp-change 選項 B 實作（0.5 day）
  ↓
P0-B within-group OLS 驗證（可用舊 master；P0-C 若無下游改動不需重跑）
  ↓
（未來）P2-E PS-aware phase-switch 追蹤
```

---

## 六、決策點

```
[ ] A — 複合 identity key（🔴 高破壞性，不推薦）
[x] B — 平行欄位（推薦，零破壞）
[ ] C — Canonicalize PS-missing 為 unphased（🟡 收益小）
```

**預設推薦**：B

**用戶選擇**：____________（填寫後即可執行 /cpp-change）

---

## 關聯文件

- `docs/reports/audit/decisions/01_P0_critical_decisions.md#P0-C`
- `docs/reports/audit/decisions/CHECKLIST.md`
- `docs/methodology/20260419_KDE_expected_coverage_audit_01.md`（前例：P0-A 方法論審查）
- `docs/reports/audit/cards/C08_Read_Level_LOSO.md`
- `docs/reports/audit/cards/C11_HPFineNGroups_Filter.md`
- `docs/reports/audit/cards/C16_HPFineNGroups_Residualized.md`
