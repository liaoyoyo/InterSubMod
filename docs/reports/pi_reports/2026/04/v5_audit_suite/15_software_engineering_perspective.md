---
title: 軟體工程視角 — Baseline vs V5 設計差別與重構分析
date: 2026-04-28
author: liaoyoyo2001
tags: [software-engineering, longphase-to, baseline, v5, refactoring, design-patterns]
status: validated_complete
audience: developer + PI
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/14_user_report_integrated.md
---

# 軟體工程視角：Baseline vs V5 設計差別與重構分析

## §0 一句話結論

> **V5 不是「修 bug」而是「重構」**：套用 SOLID（特別是 SRP 與 OCP）+ 防禦性編程 + Fail-Safe Design，把 baseline 的單體決策（monolithic if-else）重構為**可擴展、可測試、可解釋**的三層架構，同時用 `--pon-only-phasing` 開關保留向後相容性。

### V5 應用的 5 個軟工原則

| # | 原則 | Baseline 違反 | V5 應用 |
|:-:|------|------|------|
| 1 | **Separation of Concerns** | phase 與 tag 階段邏輯耦合（somatic 同時影響 phase anchor 與 tag vote） | flag 只影響 phase 階段，tag 階段獨立修法 |
| 2 | **Single Responsibility** | `getVote()` 同時負責 voting + encoding | 拆成 Layer 1 / 1.5 / 2，各層獨立職責 |
| 3 | **Open/Closed** | 改規則需動 core code | `--pon-only-phasing` flag 為擴展點 |
| 4 | **Defensive Programming** | 無 UNDEFINED guard，null vote 污染 countMap | `countSNPHaplotype` + UNDEFINED guard |
| 5 | **Fail-Safe Design** | 找不到 anchor 就 lost (HP:i:0) | Layer 1.5 fallback + HP:i:33 conservative |

![Figure 15a — 5 SE 維度對比](figures/15_se_perspective/fig15a_se_dimensions.png)

---

## §1 程式設計視角的問題分類

### 1.1 從「軟工原則」回看 baseline

長期維護的程式碼基本上是 layer cake — 越底下的 layer 改動成本越高。Baseline 的問題不只是「邏輯錯」，而是**架構性 design smell**：

```
       ┌─────────────────────────────────┐
       │   ClairS-TO 已給 PASS somatic    │ ← Frozen API
       └─────────────────────────────────┘
                      ↓
       ┌─────────────────────────────────┐
       │  longphase-to phase             │ ← Baseline: monolithic
       │  (PhasingProcess.cpp)           │   addEdge() 不分 origin
       │                                 │   → SOMATIC 也當 anchor
       └─────────────────────────────────┘
                      ↓ phased.vcf
                      ↓
       ┌─────────────────────────────────┐
       │  longphase-to haplotag          │ ← Baseline: getVote()
       │  (HaplotagProcess.cpp)          │   單體 if-else
       │                                 │   無 fallback
       └─────────────────────────────────┘
                      ↓
                 tagged.bam
                  17.3 : 1 ❌
```

### 1.2 三類 design smells

| Smell 類型 | Baseline 症狀 | 軟工術語 |
|----------|-------|------|
| **Mixed responsibility** | `getVote()` 一函式做兩件事 | violates SRP |
| **Tight coupling** | phase 與 tag 都依賴「somatic 可作 anchor」假設 | violates Loose Coupling |
| **No graceful degradation** | LOH 區無 germline → 直接 lost | violates Fail-Safe |

→ V5 不是「補丁式」修 bug，而是針對這三類 smell 做系統性重構。

---

## §2 案例 #1：getVote() 的單一職責原則演進

### 2.1 Baseline 的 monolithic 設計

```cpp
// HaplotagProcess.cpp baseline (concept)
void getVote(countMap, &hpResult) {
    int germHP1 = countMap[HAPLOTYPE1];
    int germHP2 = countMap[HAPLOTYPE2];

    if (germHP1 > 0 || germHP2 > 0) {
        if (germHP1 >= germHP2)
            hpResult = 1;     // germline HP1
        else
            hpResult = 2;     // germline HP2
    }
    else {
        hpResult = 0;         // ❌ lost
    }
}
```

**SRP 違反**：函式同時做：
- (a) Determine germline haplotype (voting)
- (b) Convert to BAM tag value (encoding)
- (c) Handle missing data (error case)

→ 三個職責糾纏在一個 if-else，無法獨立改其中一個。

### 2.2 V5 的職責分層

```cpp
// HaplotagProcess.cpp:512-563 (V5)
void getVote(countMap, &hpResult) {
    // Layer 1: Germline determination
    int germlineResult = 0;
    if (germHP1 > 0 || germHP2 > 0) {
        germlineResult = (germHP1 >= germHP2) ? 1 : 2;
    }
    // Layer 1.5: Somatic fallback (V5 NEW)
    else if (somHP1 > 0 || somHP2 > 0) {
        germlineResult = (somHP1 >= somHP2) ? 1 : 2;
    }
    // else: germlineResult stays 0

    // Layer 2: BAM tag encoding
    if (somaticTotal > 0) {
        if      (germlineResult == 1) hpResult = 11;  // HP:i:11
        else if (germlineResult == 2) hpResult = 21;  // HP:i:21
        else                          hpResult = 33;  // HP:i:33
    } else {
        hpResult = germlineResult;
    }
}
```

**三層獨立可測試**：
- Layer 1：只判 germline → 可寫 unit test 給 (germHP1, germHP2) 各種組合
- Layer 1.5：只判 somatic fallback → 可寫 unit test 給 (somHP1, somHP2)
- Layer 2：只 encoding → 可寫 unit test 給 (germlineResult, somaticTotal)

### 2.3 重構演進歷程（3 版本對比）

![Figure 15b — getVote() 重構演進](figures/15_se_perspective/fig15b_refactor_evolution.png)

| 版本 | 主要改動 | SE 改進 |
|------|---------|--------|
| Baseline | 單體 if-else | (起點) |
| V3-Fixed | Encoding 從 voting 拆出 | SRP step 1：voting / encoding 分離 |
| **V5** | + Layer 1.5 fallback | SRP step 2：完整三層；+ Fail-Safe |

→ V5 是漸進式重構（不是一次大改），每步都通過獨立驗證。

---

## §3 案例 #2：--pon-only-phasing flag 的開放/封閉原則

### 3.1 Open/Closed Principle (OCP)

> **Software entities should be open for extension, but closed for modification.**

V5 想增加「PON-only anchor」行為，但又不能破壞既有 baseline。經典的 OCP 解法是 **feature flag**：

```cpp
// PhasingProcess.cpp:154-157 (V5 NEW)
if (params.ponOnlyPhasing) {  // ← Extension point (flag)
    vGraph->convertNonGermlineToSomatic();
}
// 其他 phasing 邏輯完全不動
```

```cpp
// PhasingGraph.cpp:1139-1145 (V5 NEW)
void VairiantGraph::convertNonGermlineToSomatic() {
    for (auto variantIter = variantPosType->begin();
         variantIter != variantPosType->end(); variantIter++) {
        if (variantIter->second.origin != PON) {
            variantIter->second.origin = SOMATIC;
        }
    }
}
```

### 3.2 Strategy Pattern 應用

variant 的 `origin` field 是一個 strategy key — phasing graph 看 origin 決定該 variant 怎麼處理：

![Figure 15c — Strategy Pattern 應用](figures/15_se_perspective/fig15c_strategy_pattern.png)

| Origin | Baseline strategy | V5 strategy |
|--------|------------------|-------------|
| **PON** | ✓ Phase anchor，✓ NonSomatic filter | ✓ Phase anchor (only!)，✓ NonSomatic filter |
| **SOMATIC** (ClairS-TO PASS) | ❌ **誤用為 anchor** → bias | ✓ GT2/GT3 sub-genotype only |
| **ORIGIN_UNDEFINED** | Unphased | Unphased |

V5 的 `convertNonGermlineToSomatic()` 把 **strategy 切換**做得**乾淨**：
- 不修改 baseline 的 `addEdge()` / `edgeConnectResult()`
- 只在 phasing 主流程前**重新分類** variant origin
- 之後 phasing graph 自動因 origin 不同而採用不同策略

### 3.3 為什麼用 flag 而非分支？

**反例（不好的設計）**：
```cpp
// 如果用 if/else 散落 phasing graph 各處：
if (params.ponOnlyPhasing && origin != PON) skip();
if (params.ponOnlyPhasing && origin != PON) ignore_edge();
// → 多個地方需要同步條件，違反 DRY
```

**正解（V5 採用）**：
```cpp
// 一次性轉換 origin tag，後續 logic 自動 dispatch
convertNonGermlineToSomatic();
// → DRY + 容易測試 + 容易回退
```

→ **Single point of behavior change**，最小侵入式設計。

---

## §4 案例 #3：UNDEFINED guard 的防禦性編程

### 4.1 Baseline 的 silent corruption bug

```cpp
// Baseline countSNPHaplotype (concept, simplified)
void countSNPHaplotype(base, haplotypeBase, countMap, &HP) {
    if (base == haplotypeBase.ref) {
        countMap[haplotypeBase.refHaplotype]++;  // ❌ 如果 refHaplotype 是 UNDEFINED?
        HP = haplotypeBase.refHaplotype;
    }
    else if (base == haplotypeBase.alt) {
        countMap[haplotypeBase.altHaplotype]++;  // 同樣的問題
        HP = haplotypeBase.altHaplotype;
    }
}
```

問題：當 variant 的 `refHaplotype == HAPLOTYPE_UNDEFINED`（例如沒有 phasing 訊息的 site），`countMap[UNDEFINED]++` 會污染 vote 計算。

### 4.2 V5 的 Guard Clauses

```cpp
// HaplotagProcess.cpp:484-510 (V5)
void countSNPHaplotype(base, haplotypeBase, countMap, &HP) {
    if (base == haplotypeBase.ref &&
        haplotypeBase.refHaplotype != HAPLOTYPE_UNDEFINED) {  // ← Guard
        countMap[haplotypeBase.refHaplotype]++;
        HP = haplotypeBase.refHaplotype;
    }
    else if (base == haplotypeBase.alt) {
        if (haplotypeBase.altHaplotype != HAPLOTYPE_UNDEFINED) {  // ← Guard
            countMap[haplotypeBase.altHaplotype]++;
            HP = haplotypeBase.altHaplotype;
        }
    }
}
```

### 4.3 Defensive Programming 原則

```
程式碼層級：「永遠不信任輸入是 well-formed」
```

V5 對每個 enum 比對前都加 guard：
- `HAPLOTYPE_UNDEFINED` (enum value 0) 是合法的「無資料」狀態
- 缺乏 guard 等於 **silent corruption**（vote 數字看起來合理但其實污染）
- 加 guard = **fail-fast**：無資料時不投票，下游 `getVote()` 自然回到 `germlineResult = 0`

---

## §5 案例 #4：Fail-Safe Design — HP:i:33 conservative tagging

### 5.1 Fail-Fast vs Fail-Safe

| 設計選擇 | 何時 fail | 行為 |
|---------|---------|------|
| **Fail-Fast** | 第一個錯誤就 abort | 適合 dev / test |
| **Fail-Safe** | 仍生產輸出但標註不確定 | 適合 production |

V5 對 ambiguous reads（無法確定 HP1 還 HP2 的 somatic-bearing reads）採用 **Fail-Safe**：

```cpp
// V5 Layer 2 (HaplotagProcess.cpp:556-559)
if (somaticTotal > 0) {
    if      (germlineResult == 1) hpResult = 11;  // 清晰 HP1
    else if (germlineResult == 2) hpResult = 21;  // 清晰 HP2
    else                          hpResult = 33;  // ⚠ 模糊 (Fail-Safe)
}
```

### 5.2 為什麼不直接 lost？

| 設計選擇 | 結果 | 問題 |
|---------|------|------|
| Baseline (lost) | 模糊 reads → HP:i:0 | 下游無法區分「真的沒 somatic」與「有 somatic 但模糊」 |
| **V5 (HP:i:33)** | 模糊 reads → HP:i:33 | 下游可選：(a) 排除 HP:i:33 (b) imputation (c) 標記為「需審查」 |

→ V5 把**決策權交給下游**，這是 production-grade 設計。

### 5.3 0.6 purity 場景的 Fail-Safe 價值

從 0.6 purity simulation（見 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md`）:

| Scenario | HP33 比例 |
|----------|:--------:|
| BL @ 0.93 | 1.5% |
| V5 @ 0.93 | 8.2% |
| BL @ 0.6 | 2.0% |
| **V5 @ 0.6** | **12.4%** ← Fail-Safe 顯著啟用 |

V5 在低 purity 下把不確定 reads 推到 HP:i:33（**12.4%**），這是設計意圖的正確 behavior — 不是 bug，是 conservative tagging 的合理表現。

---

## §6 案例 #5：syncPhasingResultOrigins — 狀態一致性設計

### 6.1 為什麼需要這個函式？

V5 的 PON-only flag 把所有非 PON 標為 SOMATIC，但 phasing 結束後，**有些 SOMATIC origin 的 variants 實際是 germline**（之前被誤判）。如果直接輸出，會在 VCF 留下不一致狀態：
- origin = SOMATIC
- 但 GT 是 germline 形式（`0|1`/`1|0`）

### 6.2 V5 的解決方案

```cpp
// PhasingGraph.cpp:1155-1180 (V5 NEW)
void VairiantGraph::syncPhasingResultOrigins(PosPhasingResult &posPhasingResult) {
    // 確保 origin tag 與 GT 格式一致
    for (auto it = posPhasingResult.begin(); it != posPhasingResult.end(); ) {
        auto variantIter = variantPosType->find(it->first);
        if (variantIter != variantPosType->end()) {
            // SOMATIC 但實際是 germline (origin 升級為 ORIGIN_UNDEFINED)
            // → 寫 germline GT 格式
            // OR 移除（避免錯誤標記）
            ...
        }
    }
}
```

### 6.3 軟工原則：State Consistency Invariant

```
不變式 (Invariant)：origin tag 與輸出 GT 格式必須對應
- PON → 0|1 / 1|0 (germline het)
- SOMATIC → 0|0 / .|0 / 0|. (somatic)
- ORIGIN_UNDEFINED → 0/1 (unphased)
```

V5 的 `syncPhasingResultOrigins` 是**主動維護不變式**，避免 silent inconsistency。

---

## §7 軟工流程：從 Bug Detection 到 Fix

### 7.1 Detection Phase（如何發現 17.3:1 bias？）

| 階段 | 方法 | 工具 |
|------|------|------|
| 1 | Aggregate analysis | pysam HP family count |
| 2 | Per-site visualization | IGV 截圖 |
| 3 | Compare with paired ground truth | IGV cross-reference |
| 4 | Statistical test | binom test on HP1 vs HP2 |

### 7.2 Root Cause Analysis (RCA)

採用 **5-Why** 技巧：

```
Q1: 為什麼 HP1:HP2 = 17.3:1？
A1: 因為 somatic reads 都進 HP1 phase block

Q2: 為什麼 somatic reads 都進同一 phase？
A2: 因為 phasing graph 把 somatic site 當 anchor

Q3: 為什麼 phasing graph 用 somatic 當 anchor？
A3: 因為 baseline 沒有 origin 過濾

Q4: 為什麼 baseline 不過濾 origin？
A4: 因為原本設計沒考慮 same-clone read 共現問題

Q5: 為什麼沒考慮共現問題？
A5: 因為 paired mode（germline + tumor）下不會出現 → TO mode 的盲點
```

→ Root cause = **TO mode 設計時延用 paired mode 假設**。

### 7.3 Fix Strategy: Minimal Invasive Change

V5 採用 **少改動但精準** 策略：

| Metric | 數字 |
|--------|:----:|
| 修改的 .cpp 檔案 | 3 個 |
| 新增/修改 lines | +68 / -36 |
| 新增 functions | 2 個（`convertNonGermlineToSomatic`, `syncPhasingResultOrigins`）|
| 新增 flags | 1 個（`--pon-only-phasing`）|
| 對 ClairS-TO calling 的影響 | **0**（calling 完全不變） |

→ 修改範圍最小、影響最可控。

### 7.4 Testing Strategy

V5 的測試覆蓋（驗證見 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md`）：

| 測試類型 | 通過率 |
|---------|:----:|
| 守恆律 A：HP33 + (HP11+HP21) Δ = 0 | 15/15 ✓ |
| 守恆律 B：germline (HP1, HP2) 不變 | 15/15 ✓ |
| Layer 1.5 預期 1：新 HP11/HP21 應在 V3F 是 HP33 | 15/15 ✓ |
| Layer 1.5 預期 2：無 germline → HP33 違規 | 0 violations |
| Untagged → directional：V5 不強行 tag 原 HP=0 | 0 violations |

→ **5 項硬性 sanity check 全 PASS**。

### 7.5 Cross-Purity 驗證（V5 的 robustness）

從 0.6 purity simulation 額外驗證：

| 場景 | V5 行為 | 軟工解讀 |
|------|--------|---------|
| @ 0.93 | 修復 17.3:1 bias | 設計目標達成 |
| @ 0.6 | self-phasing 自然弱 → V5 轉為 conservative tagging | **Graceful degradation** |
| @ 0.6 HP33% 12.4% | 標多 ambiguous | **Fail-Safe 啟用** |

→ V5 在不同情境下「該強強、該保守保守」，這是 production-ready 設計。

---

## §8 邏輯與合理性檢核

### 8.1 V5 邏輯閉環檢核

| 邏輯 | 起點 | 終點 | 是否閉環 |
|------|-----|------|:----:|
| PON-only flag → 限制 anchor → 自然修 bias | PhasingProcess.cpp:154 | tagged.bam HP ratio | ✓ |
| Layer 1.5 fallback → 救起 LOH reads | HaplotagProcess.cpp:537-548 | AMB% 17.5% → 8% | ✓ |
| HP:i:33 fail-safe → 下游可選擇 | HaplotagProcess.cpp:556-559 | downstream pipeline | ✓ |
| syncOrigins → VCF state consistent | PhasingGraph.cpp:1155-1180 | phased.vcf | ✓ |

### 8.2 合理性檢核（潛在風險）

| Risk | 評估 | 緩解 |
|------|-----|------|
| Layer 1.5 fallback 可能 over-correct | 已 sanity check 0 violations | 通過 |
| HP:i:33 比例升高需下游適應 | V5 改動是新增 encoding，向後相容 | 通過 |
| --pon-only-phasing 可關掉回到 baseline | flag 設計原意 | 通過 |
| PON 數量在 LOH 區可能不足 | 0.6 simulation 顯示 PS coverage 87.5% 仍足夠 | 通過 |

→ **無顯著 logic gap 或 unhandled case**。

### 8.3 既有結論的穩定性

| audit suite 結論 | 結論依然成立？ |
|-----------------|:-----:|
| V5 修復 17:1 bias | ✓ (0.93 場景) |
| V5 不傷 ClairS-TO calling | ✓ (Δ F1 = -0.0003) |
| V5 在 clean PS blocks 勝出 | ✓ (+13.3pp paired concordance) |
| V5 Sanity check 全 PASS | ✓ (15/15) |
| V5 跨 purity stable | ✓ (0.6 conservative tagging 啟用) |

---

## §9 SE 視角結論

### 9.1 Baseline vs V5 設計成熟度對比

| 維度 | Baseline | V5 |
|------|:--:|:--:|
| 程式碼可讀性 | 中 | **高**（三層分明） |
| 可測試性 | 低（單體 if-else） | **高**（每層獨立可測） |
| 可擴展性 | 低（要改規則需動 core） | **高**（flag 為擴展點） |
| 防禦性 | 低（無 guard） | **高**（UNDEFINED guard） |
| 故障恢復 | 低（lost reads） | **高**（HP33 fallback） |
| 向後相容 | n/a | **高**（flag 默認 false） |
| 跨情境 robust | 低（只 paired 場景驗過） | **高**（0.93 + 0.6 都驗證） |

### 9.2 V5 不是 patch，是重構

V5 對 baseline 的 4 個 bug 修法不是局部補丁，而是**架構性升級**：

| Bug | 補丁式做法 | V5 重構式做法 |
|-----|-----------|--------------|
| getVote() Layer 1 only | 加 `if (no germline) {}` | **三層職責拆分** + Layer 1.5 fallback |
| HP:i:33 enum 比對 | 改一個比對 | enum 比對 + encoding 邏輯獨立 |
| countSNPHaplotype 污染 | 加 if check | **完整 guard pattern** + 防禦性編程 |
| PhasingProcess 無過濾 | 加 if 散落各處 | **convertNonGermlineToSomatic()** 一次性 strategy 切換 |

### 9.3 給後續開發者的 takeaways

1. **重構 > 補丁**：當 bug 來自 architectural smell，補丁只會累積技術債
2. **Flag 不只是 toggle**：`--pon-only-phasing` 是 OCP 擴展點，讓 baseline 與 V5 並存
3. **Fail-Safe > Fail-Fast for production**：HP:i:33 把決策權留給下游
4. **Invariant 必維護**：`syncPhasingResultOrigins` 看似多餘，實則維持狀態一致
5. **Cross-context 驗證**：0.6 simulation 是 V5 robustness 的驗證 — 不只看設計目標的場景

---

## §10 跨檔索引

### 程式碼引用（關鍵行號）

| 檔案 | 行號 | 函式 | SE 主題 |
|------|------|------|--------|
| `HaplotagProcess.cpp` | 484-510 | `countSNPHaplotype` + UNDEFINED guard | §4 防禦性編程 |
| `HaplotagProcess.cpp` | 512-563 | `getVote()` 三層 | §2 SRP, §5 Fail-Safe |
| `PhasingProcess.cpp` | 154-157 | `--pon-only-phasing` 入口 | §3 OCP |
| `PhasingGraph.cpp` | 1139-1145 | `convertNonGermlineToSomatic()` | §3 Strategy Pattern |
| `PhasingGraph.cpp` | 1155-1180 | `syncPhasingResultOrigins()` | §6 State Invariant |

### 相關 audit suite 文件

| # | 文件 | 與本 SE 報告關係 |
|---|------|---------|
| 01 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` | 詳細 commit-level diff |
| 06 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/06_v5_sanity_bug_check.md` | §7.4 testing strategy 引用 |
| 09 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` | §7.5 cross-purity 驗證 |
| 11 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md` | 5 bugs 清單對應 |
| 13 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/13_phase_vs_tag_algorithm_detail.md` | 演算法細節 |
| 14 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/14_user_report_integrated.md` | PI 整合層 |

### 圖檔索引

| 圖 | 路徑 |
|----|------|
| Figure 15a 5 SE 維度對比 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/15_se_perspective/fig15a_se_dimensions.png` |
| Figure 15b 重構演進 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/15_se_perspective/fig15b_refactor_evolution.png` |
| Figure 15c Strategy Pattern | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/15_se_perspective/fig15c_strategy_pattern.png` |
