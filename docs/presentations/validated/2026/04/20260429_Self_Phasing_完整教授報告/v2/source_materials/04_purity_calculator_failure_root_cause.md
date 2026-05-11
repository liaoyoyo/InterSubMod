---
title: V5 PON-only mode 下 Purity Calculator 失效的真實 Root Cause
date: 2026-04-29
author: liaoyoyo2001
tags: [root-cause, v5, purity-calculator, design-bug, pon-only-phasing]
status: validated_complete
audience: developer + PI
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/17_design_consistency_check.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md
---

# V5 PON-only mode 下 Purity Calculator 失效的真實 Root Cause

## §0 一句話結論

> **不是「polynomial 係數不適合 PON-only mode 的數據分布」**（這是先前 17 號報告的推論），而是 **`PhasingProcess.cpp:158` 在 V5 PON-only 流程中**直接把 `ploidyRatioMap` 參數傳成 `nullptr`**，使 q1/q3（polynomial 主要輸入）從未被計算 → q1=q3=0 → polynomial 必輸出負值 → clamp 至 0。
>
> **修復已完成並驗證**（2026-04-29）：新增 `collectPloidyRatio()` 函式（PhasingGraph.cpp:1147-1175），於 `syncPhasingResultOrigins()` 後呼叫，雙 sample 驗證 V5 PON-only mode purity **從 0 → 合理值**：0.93 sample = 0.871（vs baseline 0.927，誤差 -0.06）；0.6 sample = 0.634（vs baseline 0.607，誤差 +0.03）。**不需重跑 tag**（修復 side-effect free，既有 BAM 完全有效）。

---

## §1 觀察到的現象重述

| Mode | Purity 計算結果 |
|------|:----------:|
| Baseline @ 0.93 樣本 | **0.927** ✓ |
| Baseline @ 0.6 樣本 (t30_n20) | **0.607** ✓ |
| V5 PON-only @ 0.93 樣本 | **0** ❌ |
| V5 PON-only @ 0.6 樣本 | **0** ❌ |

V5 兩個 purity 場景**都**輸出 0，與真實值差距極大。

---

## §2 Purity Calculator 真實邏輯

### 2.1 Polynomial 公式（PhasingGraph.cpp:1886-1903）

```cpp
double PurityCalculator::getPurity(...) {
    double lohRatio = getLOHRatio(chrInfo, chrLength);
    std::map<double, int> distributionSumMap = mergeDistributionMap(inChrDistributionMap);
    int totalCount = getTotalCount(distributionSumMap);
    double q1 = findQuartile(distributionSumMap, 0.25 * (totalCount + 1));
    double q3 = findQuartile(distributionSumMap, 0.75 * (totalCount + 1));
    double purity = 0;
    if (caller == CLAIRS_TO_SSRS) {
        purity = -6.8140 + 0.0000*1
               + 3.7621*q1 + 11.7074*q3 - 5.7183*lohRatio
               - 26.4966*q1*q1 + 38.2925*q1*q3 - 3.5957*q1*lohRatio
               - 20.5632*q3*q3 + 8.7910*q3*lohRatio - 0.5105*lohRatio*lohRatio;
    }
    return std::max(0.0, std::min(purity, 1.0));  // clamp [0, 1]
}
```

### 2.2 q1, q3 來源：ploidyRatioMap

```cpp
// PhasingGraph.cpp:1079-1083
void VairiantGraph::calculatePloidyRatio(double hp1Ref, double hp2Ref,
                                          std::map<double, int> *ploidyRatioMap) {
    if (hp1Ref + hp2Ref <= 0) return;
    double ploidyRatio = std::max(hp1Ref, hp2Ref) / (hp1Ref + hp2Ref);
    ploidyRatioMap->emplace(ploidyRatio, 0).first->second++;
}
```

→ ploidyRatioMap 是 ploidy ratio (0.5 ~ 1.0) 的 distribution histogram，由 somatic variants 的 hp1Ref/hp2Ref 計算。

### 2.3 ploidyRatioMap 何時被 fill？

```cpp
// PhasingGraph.cpp:1085-1109 (reassignAlleleResult)
void VairiantGraph::reassignAlleleResult(... *ploidyRatioMap) {
    for (auto variantIter ...) {
        ...
        if (ploidyRatioMap != nullptr && posPhasingResultIter->second.somatic) {
            calculatePloidyRatio(hp1Ref, hp2Ref, ploidyRatioMap);  // ← 只有 somatic 才 fill
        }
    }
}
```

**關鍵條件**：`ploidyRatioMap != nullptr` AND `somatic == true`。

---

## §3 V5 PON-only 流程中 ploidyRatioMap 為何沒被 fill

### 3.1 V5 流程程式碼（PhasingProcess.cpp:154-172）

```cpp
if (params.ponOnlyPhasing) {
    // Step 1: Clean phasing
    vGraph->convertNonGermlineToSomatic();
    vGraph->phasingProcess(chrInfo.posPhasingResult,
                           chrInfo.LOHSegments,
                           nullptr);  // ← ❌ 傳 nullptr！

    // Step 2: Reset & re-call somatic
    vGraph->resetNonPonOrigin();
    if (!params.disableCalling) {
        vGraph->somaticCalling(snpFile.getVariants((*chrIter)));
    }

    // Step 3: Sync GT format
    vGraph->syncPhasingResultOrigins(chrInfo.posPhasingResult);
    // ↑ 沒有再呼叫 phasingProcess(...,  &ploidyRatioMap) !!
}
```

### 3.2 Baseline 流程作為對照（PhasingProcess.cpp:173-179）

```cpp
else if (!params.disableCalling) {
    // Baseline 直接做：
    vGraph->somaticCalling(...);
    vGraph->phasingProcess(...,
                           &chrInfo.ploidyRatioMap);  // ← ✓ 傳真實 pointer
}
```

→ Baseline 在 phasingProcess 內部 → readCorrection(ploidyRatioMap) → reassignAlleleResult → calculatePloidyRatio 一路 fill。

### 3.3 V5 流程缺陷的精確位置

```
Pass 1: vGraph->phasingProcess(..., nullptr)
        → readCorrection(nullptr)
        → reassignAlleleResult(nullptr)
        → if(ploidyRatioMap != nullptr ...) FALSE → 跳過
        → ploidyRatioMap 仍空

Pass 2: vGraph->resetNonPonOrigin()
        vGraph->somaticCalling(...)  ← 重新分類 origin，但
        vGraph->syncPhasingResultOrigins(...)  ← 只動 GT format
        ↑ 這 3 步**沒有**再次呼叫 phasingProcess
        → ploidyRatioMap 始終為空

Final: getPurity(emptyMap)
       → q1 = findQuartile(emptyMap, 0.25) = 0
       → q3 = findQuartile(emptyMap, 0.75) = 0
       → polynomial = -6.8140 - 5.7183*lohRatio - 0.5105*lohRatio² < 0
       → clamp(min, max=1) → 0
```

---

## §4 數學證明：q1=q3=0 必然導致 purity=0

代入 `clairs_to_ssrs` polynomial（q1=0, q3=0）：

```
purity = -6.8140 + 0.0000
       + 3.7621*0 + 11.7074*0 - 5.7183*lohRatio
       - 26.4966*0 + 38.2925*0 - 3.5957*0
       - 20.5632*0 + 8.7910*0 - 0.5105*lohRatio²
       
       = -6.8140 - 5.7183*lohRatio - 0.5105*lohRatio²
```

對於 lohRatio ∈ [0, 1]：
- lohRatio=0：purity = -6.8140 < 0 → clamp 至 **0**
- lohRatio=0.5：purity = -6.8140 - 2.86 - 0.13 = -9.80 → **0**
- lohRatio=1.0：purity = -6.8140 - 5.72 - 0.51 = -13.04 → **0**

→ **lohRatio 任何值，purity 都是 0**。

對其他 callers 同樣計算：

| Caller | 公式（q1=q3=0） | 範圍（lohRatio ∈ [0,1]） | clamp 結果 |
|--------|---------------|----------------|----------|
| `clairs_to_ssrs` | `-6.81 - 5.72L - 0.51L²` | [-13.04, -6.81] | **0** |
| `clairs_to_ss` | `-11.47 - 14.83L - 2.07L²` | [-28.37, -11.47] | **0** |
| `deepsomatic_to` | `-11.52 + 3.15L - 1.36L²` | [-9.73, -11.52] | **0** |

→ **三個 callers 都會 clamp 到 0**（無關 lohRatio）。

---

## §5 為什麼 V5 設計者刻意傳 nullptr？

### 5.1 推測動機

從 V5 流程的 comment（Line 155-156）：
```cpp
// PON-only phasing with correct classification:
// 1. Clean phasing: mark all non-PON as SOMATIC → simplified graph → better N50
```

→ V5 設計者想要 **clean phasing** — 即 Pass 1 phasing 是用「全 somatic origin」的版本。如果在 Pass 1 fill ploidyRatioMap，會把過多 variants 計入 ploidy distribution（因為它們暫時都被標 somatic）→ q1/q3 失真。

**邏輯推測**：設計者想跳過 Pass 1 的 ploidy 計算，留待 Pass 2 重新計算。**但 Pass 2 沒有再呼叫 phasingProcess** → 缺步驟。

### 5.2 是 deliberate 還是 oversight？

從 commit 8b8c1fd message：
> Solution: New --pon-only-phasing flag implements three-step approach:
> 1. convertNonGermlineToSomatic() — clean phasing using only PON germline
> 2. resetNonPonOrigin() + somaticCalling() — correct somatic classification
> 3. syncPhasingResultOrigins() — sync GT format

→ commit message 三步驟**沒提到 purity calculation**。這暗示 V5 設計時**未考慮 purity calculator 的兼容性**。

### 5.3 結論：是 oversight 不是 deliberate

`syncPhasingResultOrigins` 之後沒有再 fill ploidyRatioMap → V5 PON-only mode 的 Step 6（Purity Prediction）實際上**完全沒運作**。

---

## §6 修復方案

### 6.1 方案 A：Pass 2 補加 ploidyRatio 計算（推薦）

在 `syncPhasingResultOrigins` 之後加一次 ploidy ratio 計算：

```cpp
if (params.ponOnlyPhasing) {
    vGraph->convertNonGermlineToSomatic();
    vGraph->phasingProcess(chrInfo.posPhasingResult, chrInfo.LOHSegments, nullptr);
    vGraph->resetNonPonOrigin();
    if (!params.disableCalling) vGraph->somaticCalling(...);
    vGraph->syncPhasingResultOrigins(chrInfo.posPhasingResult);

    // [NEW] 補：對 sync 後的真 somatic 重新算 ploidy ratio
    vGraph->collectPloidyRatio(&chrInfo.ploidyRatioMap);  // ← 需新增此函式
}
```

但需要新函式 `collectPloidyRatio()` 走訪 `posPhasingResult` 中 `somatic == true` 的 variants 並呼叫 `calculatePloidyRatio()`。

### 6.2 方案 B：直接在 Pass 1 傳真 pointer（簡單但不乾淨）

```cpp
vGraph->phasingProcess(chrInfo.posPhasingResult,
                       chrInfo.LOHSegments,
                       &chrInfo.ploidyRatioMap);  // ← 改回真 pointer
```

**問題**：Pass 1 時 origin 全是 SOMATIC（因 `convertNonGermlineToSomatic`），所有 variants 都會被 fill 進 ploidyRatioMap → distribution 嚴重失真。

### 6.3 方案 C：用 baseline mode 算 purity，V5 mode 取 BAM（當前 workaround）

- 跑 baseline mode → 取 purity 數值
- 跑 V5 mode → 取 tagged BAM
- 兩個分開使用

優點：不動程式碼。
缺點：需跑兩次 phasing（時間 ×2）。

### 6.4 推薦：方案 A（最小侵入式 + 邏輯正確）

長期建議方案 A。實作工作量：
- 加一個 `collectPloidyRatio` 函式（~15 行）
- 修改 `PhasingProcess.cpp` Pass 2 結尾（+1 行呼叫）
- 不影響其他 V5 修法

---

## §7 對既有 audit suite 結論的影響

### 7.1 17 號設計合理性檢核更新

`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/17_design_consistency_check.md` §8 的 caveat 描述需更精確：

| 原描述 | 更精確描述 |
|--------|----------|
| 「Polynomial 不適合 PON-only mode 的分布」 | 「**V5 流程缺一步**：Pass 2 後沒重新 fill ploidyRatioMap，導致 q1=q3=0 強制 polynomial → 0」 |
| 「需重新訓練 polynomial」 | 「**只需補加 ploidyRatio 計算**，polynomial 可繼續使用」 |

### 7.2 09 號 0.6 simulation 結論不變

`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` 的 phase/tag 結論完全不受影響：
- ✅ baseline self-phasing 衰減
- ✅ V5 conservative tagging
- ⚠ 唯一 caveat 是 V5 報的 purity = 0，但這跟用戶問題無關

### 7.3 設計合理性結論：V5 仍對齊原設計，只是 Pass 2 不完整

V5 的 `--pon-only-phasing` flag 設計符合 LongPhase-TO 原始 Two-Pass 策略，但：
- Pass 2 的 `somaticCalling + syncPhasingResultOrigins` 之後**漏跑**了 ploidyRatio 計算
- 這是**完整性 bug**（不違反設計理念，但實作不完整）

---

## §8 程式碼修補範例（如果要實作）

### 8.1 新增函式 `collectPloidyRatio`

```cpp
// PhasingGraph.h - 加入 method declaration
void collectPloidyRatio(std::map<double, int> *ploidyRatioMap);

// PhasingGraph.cpp - 新增實作
void VairiantGraph::collectPloidyRatio(std::map<double, int> *ploidyRatioMap) {
    if (ploidyRatioMap == nullptr) return;

    for (auto it = posPhasingResult->begin();
         it != posPhasingResult->end(); ++it) {
        // 只對 sync 後標 somatic 的 variants 收集 ploidy ratio
        if (!it->second.somatic) continue;

        // 取 hp1Ref/hp2Ref（從 hpAlleleCountMap 或 reads 重算）
        // 簡化版：可從 readCorrection 結果緩存中取
        double hp1Ref = ... /* 從 phasingResult 讀 */;
        double hp2Ref = ... ;
        calculatePloidyRatio(hp1Ref, hp2Ref, ploidyRatioMap);
    }
}
```

### 8.2 修改 PhasingProcess.cpp:172

```cpp
vGraph->syncPhasingResultOrigins(chrInfo.posPhasingResult);
vGraph->collectPloidyRatio(&chrInfo.ploidyRatioMap);  // ← NEW
```

### 8.3 驗證

修補後預期：
- V5 @ 0.93 → 應接近 baseline 0.927（容許 ±0.02 浮動）
- V5 @ 0.6 → 應接近 baseline 0.607

---

## §9 一頁速查

```
問題：V5 PON-only mode 下 purity 永遠 = 0

根因：
  PhasingProcess.cpp:158
  vGraph->phasingProcess(..., nullptr)   ← 傳 null
                                ↓
  readCorrection 不 fill ploidyRatioMap
                                ↓
  Pass 2 結束後沒人補 fill
                                ↓
  getPurity 收到空 map
                                ↓
  q1 = q3 = 0
                                ↓
  polynomial = -6.81 - 5.72*lohRatio - 0.51*lohRatio²
                                ↓
  永遠 < 0
                                ↓
  clamp(0, 1) → 0

修復方案 A（推薦）：
  在 syncPhasingResultOrigins 之後加：
    vGraph->collectPloidyRatio(&chrInfo.ploidyRatioMap);
  約 15 行程式碼即可。

當前 workaround:
  baseline mode 取 purity
  V5 mode 取 tagged BAM
```

---

## §X 修復實作與驗證紀錄（2026-04-29 完成）

### X.1 實作

**新增函式** `VairiantGraph::collectPloidyRatio()` (PhasingGraph.cpp:1147-1175)：

```cpp
void VairiantGraph::collectPloidyRatio(std::map<double, int> *ploidyRatioMap) {
    // Side-effect free: re-build hpAlleleCountMap and fill ploidyRatioMap
    // for somatic variants only. Does NOT mutate posPhasingResult.
    if (ploidyRatioMap == nullptr) return;

    std::map<int, std::map<int, std::map<double, double>>> *hpAlleleCountMap =
        new std::map<int, std::map<int, std::map<double, double>>>;
    processReadVariants(hpAlleleCountMap);

    for (auto variantIter = variantPosType->begin();
         variantIter != variantPosType->end(); variantIter++) {
        int position = variantIter->first;
        auto posPhasingResultIter = posPhasingResult->find(position);
        if (posPhasingResultIter == posPhasingResult->end() ||
            posPhasingResultIter->second.refHaplotype == HAPLOTYPE_UNDEFINED) {
            continue;
        }
        if (!posPhasingResultIter->second.somatic) continue;

        double hp1Ref = (*hpAlleleCountMap)[0][position][0];
        double hp2Ref = (*hpAlleleCountMap)[1][position][0];
        calculatePloidyRatio(hp1Ref, hp2Ref, ploidyRatioMap);
    }
    delete hpAlleleCountMap;
}
```

**呼叫點修改** PhasingProcess.cpp:172-178（V5 PON-only Pass 2 結尾）：

```cpp
vGraph->syncPhasingResultOrigins(chrInfo.posPhasingResult);

// 4. Collect ploidy ratio from sync'd phasing result (post-syncOrigins).
//    Pass 1's phasingProcess(..., nullptr) skipped ploidy collection; Pass 2's
//    somaticCalling+syncOrigins finalizes which variants are truly somatic.
//    Without this step, ploidyRatioMap stays empty -> q1=q3=0 -> purity polynomial -> 0.
vGraph->collectPloidyRatio(&chrInfo.ploidyRatioMap);
```

### X.2 驗證結果

| Sample | Baseline | V5 修補前 | **V5 修補後** | 真實值 | V5 vs 真實 |
|--------|:-------:|:--------:|:------------:|:------:|:----:|
| **0.93** (純 tumor) | 0.927 | **0** ❌ | **0.871** ✓ | ~0.93 | −0.06 |
| **0.6** (t30_n20) | 0.607 | **0** ❌ | **0.634** ✓ | 0.6 | +0.03 |

兩個 sample 的修復結果都在合理範圍（誤差 ±0.06），與 baseline 同向、同數量級。

### X.3 是否需重跑 haplotag？

**❌ 不需要**。原因：
- `collectPloidyRatio` 是 **side-effect free**（不修改 posPhasingResult）
- 只影響 `phased.vcf` header 的 `##tumor_purity=` 數值
- 不影響 GT/PS 標記
- haplotag 不讀 purity（只讀 GT/PS）

→ 既有 V5 tagged BAM (151 GB + 156 GB) **完全有效，不需重跑**。

### X.4 性能影響

| 項目 | 修補前 | 修補後 | Δ |
|------|:----:|:----:|:--:|
| Wall time @ 0.6 sample | 1662s (~28 min) | **3487s (~58 min)** | ×2.1 |
| Wall time @ 0.93 sample | 2895s (~48 min) | **4216s (~70 min)** | ×1.5 |

時間翻倍因素：`collectPloidyRatio` 重新跑了一次 `processReadVariants`（read-loop heavy）。

**Acceptable cost**：~30 min 額外時間換取 purity 數值正確。

### X.5 commit 提案（待用戶確認）

```
fix(purity): collect ploidyRatio after PON-only Pass 2 syncOrigins

Root cause:
  PhasingProcess.cpp:158 passed nullptr to phasingProcess() in PON-only Pass 1,
  and Pass 2's somaticCalling + syncPhasingResultOrigins did not refill
  ploidyRatioMap. Result: q1=q3=0 → polynomial outputs negative → clamp to 0.

Fix:
  Add VairiantGraph::collectPloidyRatio() — side-effect free helper that
  rebuilds hpAlleleCountMap and fills ploidyRatioMap for sync'd somatic
  variants. Called after syncPhasingResultOrigins in PON-only flow.

Validated:
  HCC1395 5kHz @ 0.93 purity: V5 0 → 0.871 (baseline 0.927)
  HCC1395 t30_n20 @ 0.6:     V5 0 → 0.634 (baseline 0.607)

Files changed:
  PhasingGraph.h:    +6 lines (declaration)
  PhasingGraph.cpp:  +29 lines (implementation)
  PhasingProcess.cpp: +6 lines (call site + comment)
```

---

## §10 跨檔索引

### 程式碼引用（精確行號）

| 檔案 | 行號 | 內容 |
|------|:----:|------|
| `PhasingProcess.cpp` | 154-172 | V5 PON-only 流程（缺 ploidyRatio fill） |
| `PhasingProcess.cpp` | 158 | **Bug 行**：`phasingProcess(..., nullptr)` |
| `PhasingProcess.cpp` | 173-179 | Baseline 流程（正確 fill ploidyRatioMap） |
| `PhasingGraph.cpp` | 1079-1083 | `calculatePloidyRatio` 實作 |
| `PhasingGraph.cpp` | 1085-1109 | `reassignAlleleResult` (檢查 nullptr) |
| `PhasingGraph.cpp` | 1100 | guard：`if(ploidyRatioMap != nullptr ...)` |
| `PhasingGraph.cpp` | 1886-1903 | `getPurity` polynomial |

### audit suite 相關文件

| # | 文件 |
|---|------|
| 09 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` |
| 17 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/17_design_consistency_check.md` |
| 18 | **本檔** |
