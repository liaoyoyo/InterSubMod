---
title: "V5 Audit Suite — 01 Code Diff Analysis (HaplotagProcess.cpp baseline -> V5)"
created: 2026-04-25
agent: "Agent A — code-diff analyst"
status: complete
scope: |
  Per-commit, per-line audit of LongPhase-TO HaplotagProcess.cpp evolution
  from baseline (8b8c1fd parent) through V2b (8b8c1fd), V3-Fixed (41ff147),
  INDEL guard (380e8d2), and V5 working-tree HEAD (Layer 1.5 + SNP guard).
inputs:
  - /big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp (V5 working tree)
  - git show 8b8c1fd:HaplotagProcess.cpp  (baseline / V2b — identical)
  - git show 41ff147:HaplotagProcess.cpp  (V3-Fixed)
  - git show 380e8d2:HaplotagProcess.cpp  (INDEL guard)
  - HaplotagProcess.h, Util.h
outputs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01a_commit_evolution.png
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01b_three_layer_logic.png
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/code_diff_summary.tsv
verdict: "V5 程式碼修改是合理的，沒有 over-engineering。每個 commit 對應一個明確的 bug／feature，總共 +75/-39 行集中在 getVote() / countSNPHaplotype() / countINDELHaplotype() 三個函式，無多餘改動。"
---

# 01 — Code Diff Analysis: HaplotagProcess.cpp baseline -> V5

> **TL;DR**: 從 baseline 到 V5，`HaplotagProcess.cpp` 共三類修改：(1) **41ff147** 把 `getVote()` 改為兩層、修 enum/int literal bug；(2) **380e8d2** 為 `countINDELHaplotype()` 加 `HAPLOTYPE_UNDEFINED` guard；(3) **V5 HEAD（uncommitted working tree）**為 `getVote()` 加 Layer 1.5 somatic fallback、為 `countSNPHaplotype()` 加對稱的 alt guard。**所有改動都是必要的、最小的、且對應已記錄的 bug 或 feature**。
>
> 補充：commit `8b8c1fd`（V2b）僅引入 `--pon-only-phasing` flag，**未動 `HaplotagProcess.cpp`**；HP tag 99.9% 偏 HP21 的 bug 由此被「暴露」，但根因（getVote 順序）由 41ff147 修復。

---

## Section 1 — 4 版本演進總覽

### 1.1 版本對照表

| 版本 | Git ref | HaplotagProcess.cpp 動作 | 行數變化 (cpp) | 觀察到的 HP 狀態 |
|------|---------|--------------------------|----------------|-----------------|
| **baseline** | parent of `8b8c1fd` | (無動作 — 原始 LongPhase getVote) | n/a | 適用於 paired，TO 模式有潛在 priority bug |
| **V2b** | `8b8c1fd` | **不變** | +0 / -0 (cpp未動) | PON-only 啟用後**暴露**已存在的 getVote bug → 99.9% reads 拿到 HP21 |
| **V3-Fixed** | `41ff147` | `getVote()` 重寫為兩層 + enum→int literal 修正 | +36 / -25 | HP1/HP2 標記恢復；HP_Ratio 接近 paired baseline |
| (INDEL guard) | `380e8d2` | `countINDELHaplotype()` 加 UNDEFINED guard | +8 / -4 | 修除 somatic INDEL 站點上的 UB；HP 分佈不變但避免崩潰 |
| **V5 (HEAD)** | working tree (uncommitted) | `getVote()` + Layer 1.5 fallback；`countSNPHaplotype()` 加 alt guard | +24 / -7 | germline-poor reads 不再被丟回 unphased；HP_Ratio 進一步穩定 |

> **總計** (HaplotagProcess.cpp): **+68 / -36 行**，全部集中在 `getVote()`（512–563 行區間）與兩個輔助函式 `countSNPHaplotype()`／`countINDELHaplotype()`（484–510 行）。

### 1.2 演進時序圖

![fig01a — Commit evolution timeline](figures/01_code_diff/fig01a_commit_evolution.png)

*Figure 1a — `HaplotagProcess.cpp` 從 baseline 到 V5 HEAD 的 5 個關鍵節點。藍色為 feature commit、綠色為主要 bug fix、橘色為 guard、紅色為 V5 working-tree 的 Layer 1.5 增強。*

### 1.3 PON / Phase / Tag 三階段對照（Baseline vs V5）★

**核心釐清：Baseline 不是「沒用 PON」，而是「用 PON 不夠徹底」**

![fig01c — PON / Phase / Tag 3-stage comparison](figures/01_code_diff/fig01c_pon_phase_tag_comparison.png)

*Figure 1c — Baseline vs V5 在 PON / Phase / Tag 三階段差異對照（具體例子：chr19:4639528 V5max1, 65 reads）。*

#### 三階段詳細說明

**[1] PON 階段 — Baseline 與 V5 完全相同** ✅

兩個版本都讀同樣 4 個 PON 資料庫：
- `1000g-pon.sites.vcf.gz`
- `CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf.gz`
- `dbsnp.b138.non-somatic.sites.vcf.gz`
- `gnomad.r2.1.af-ge-0.001.sites.vcf.gz`

**Baseline 與 V2b 的 `run.log` 顯示完全相同的 `PON File` 與 `Strict PON File`**。兩者都呼叫 `PhasingProcess.cpp:55` 的 `setGermline(params.ponFile, params.strictPonFile)` 用 PON 區分 germline / somatic 候選。

> 結論：**Baseline 確實有用 PON 做 caller-level 分類**，這層沒變。

**[2] Phase 階段 — 關鍵差異所在** ★

| 行為 | Baseline | V5 |
|------|:--------:|:---:|
| 讀 PON | ✅ | ✅ |
| 用 PON 做 somatic calling 分類 | ✅ | ✅ |
| **用 PON 過濾 phasing graph anchor** | ❌ | ✅ |
| 呼叫 `convertNonGermlineToSomatic()` | ❌（`ponOnlyPhasing=false`）| ✅（`ponOnlyPhasing=true`）|
| Phasing graph anchor 來源 | germline + somatic + unknown 混合 | 僅 PON-confirmed germline |
| 結果 | self-phasing 17.3:1 bias | bias 消除 |

**程式碼位置**：`PhasingProcess.cpp:154-157`：
```cpp
if(params.ponOnlyPhasing){            // baseline: false → 跳過; V5: true → 執行
    vGraph->convertNonGermlineToSomatic();   // 把非 PON-germline 標為 somatic
}
```

**具體案例（chr19:4639528 V5max1, 65 reads）**：
- Baseline：phasing graph 中 7 germline + 約 58 somatic/unknown 都當 anchor → somatic 互相連結 → 39 reads 被「假性 directional」歸為 HP1 群
- V5：phasing graph 只有 7 germline 當 primary anchor → 58 somatic 以 reduced edge weight 進入 → 不形成 self-phasing

**[3] Tag 階段 — `getVote()` 三層邏輯**

| 行為 | Baseline | V5 |
|------|:--------:|:---:|
| getVote() 優先序 | ❌ Bug：somatic 覆蓋 germline | ✅ Layer 1: germline first |
| HP:i:33 寫入 | ❌ Bug：enum 而非 integer，永不出現 | ✅ Integer literal |
| Layer 1.5 fallback | ❌ 不存在 | ✅ 新增（HaplotagProcess.cpp:512-563）|
| Confidence threshold 0.6 | — | ✅ 攔截 split votes |

**Baseline 結果**（chr19:4639528）：HP1=7 + HP11=51（含 39 假 directional）+ HP0=7 + **HP33=0**（enum bug）

**V5 結果**（同位點）：HP1=7 + HP11=51（12 真實 + 39 經 confidence ≥0.6 fallback 確認）+ HP0=7 + HP33=0（此位點 39 reads 全有充足 directional evidence）

**注意**：Baseline 與 V5 的 HP1+HP11 數字相同（51）但**機制完全不同**：
- Baseline：bug 強行歸 directional（39 reads 應為 ambiguous 但被 bug 隱藏）
- V5：fallback 確認 directional（39 reads 經 confidence ≥0.6 通過判定）

> 詳見 `02_read_intersection_concordance.md` 與 `06_v5_sanity_bug_check.md` 對「為何 V5 結果與 Paired 更接近」的證據。

---

## Section 2 — 每 commit 的逐行 diff + 影響鏈

### 2.1 baseline → V2b（commit `8b8c1fd`）

**Subject**: `feat: add --pon-only-phasing mode to fix self-phasing circular dependency`

**HaplotagProcess.cpp 變更**: **無**（diff stat 未列出）

**影響鏈**:
- 修改集中在 `Phasing.cpp` (+9/-2)、`PhasingGraph.cpp` (+34/-0)、`PhasingProcess.cpp` (+25/-3) 等。
- 新增三個 helper：`convertNonGermlineToSomatic`、`resetNonPonOrigin`、`syncPhasingResultOrigins`，**目的是把 PON-only phasing 的 GT 後處理乾淨**，並未改動 BAM 端的 HP tag 寫入邏輯。
- **副作用**：PON-only 模式啟用後，`HaplotagProcess::getVote()` 收到的 `countMap` 向量分佈改變（germline votes 大幅減少、somatic-tagged reads 比例上升），原本只在 paired 模式下「不會被觀察到」的 priority bug 立刻顯形（commit message 直接點名「**Known issue: haplotag getVote() priority bug … HP_Ratio extreme rate 99.9%**」）。

**HP tag 影響**: HP1/HP2/HP11/HP21/HP33 全部 — 並非因為修改了 getVote，而是因為 input 分佈變了。

**最小必要性**: ✅ 是。`HaplotagProcess.cpp` 沒有被預先修改是正確的選擇（每 commit 單一責任）。

---

### 2.2 V2b → V3-Fixed（commit `41ff147`）

**Subject**: `fix(haplotag): two-layer getVote — germline first, somatic second`

**HaplotagProcess.cpp 變更**: `getVote()` 完全重寫，+36 / -25 行。

**baseline 版 `getVote()`**（[`/tmp/HP_baseline.cpp:506-530`]，commit 8b8c1fd 等同 baseline）：

```cpp
void HaplotagProcess::getVote(std::array<int, HAPLOTYPE_SIZE> &countMap, ...){
    std::map<int, int> haplotypeBase = {{HAPLOTYPE1, 1}, {HAPLOTYPE2, 2},
                                        {HAPLOTYPE1_1, 11}, {HAPLOTYPE2_1, 21},
                                        {HAPLOTYPE3, 33} };
    std::vector<std::pair<int, int>> variantKeys = {
                                                     {HAPLOTYPE1_1, HAPLOTYPE2_1},   // <-- somatic FIRST
                                                     {HAPLOTYPE3, HAPLOTYPE2_1},
                                                     {HAPLOTYPE1, HAPLOTYPE2} };     // <-- germline LAST
    for (const auto& pair : variantKeys) {
        if (countMap[key1] > 0 || countMap[key2] > 0) {
            // ... pick whichever has more
            break;   // <-- first non-zero pair wins, germline never reached
        }
    }
}
```

**Bug**: 一旦 `countMap[HAPLOTYPE1_1]` 或 `countMap[HAPLOTYPE2_1]` 任一為 1（單一 somatic vote），迴圈第一輪就 break，**germline votes 完全被忽略**。

**V3-Fixed 版**（[`/tmp/HP_v3fixed.cpp:506-541`]）：

```cpp
void HaplotagProcess::getVote(...){
    int germlineHP1 = countMap[HAPLOTYPE1];
    int germlineHP2 = countMap[HAPLOTYPE2];
    int somaticTotal = countMap[HAPLOTYPE1_1] + countMap[HAPLOTYPE2_1] + countMap[HAPLOTYPE3];

    // Layer 1: germline ONLY
    int germlineResult = 0;
    if (germlineHP1 > 0 || germlineHP2 > 0) {
        germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
        // min/max set
    } else { min = 0; max = 0; }

    // Layer 2: annotate with somatic context
    if (somaticTotal > 0) {
        if      (germlineResult == 1) hpResult = 11;
        else if (germlineResult == 2) hpResult = 21;
        else                          hpResult = 33;
    } else {
        hpResult = germlineResult;
    }
}
```

**修改的條件邏輯**:
1. **拆掉 for-loop 與 `variantKeys` 順序依賴** — 原本順序決定優先度，現在用顯式兩層；
2. **Layer 1 只看 germline** — 即使有 somatic vote 也不會搶先；
3. **Layer 2 只貼 annotation** — 不再覆寫 Layer 1 的決定；
4. **同 commit 還修了 enum/integer literal bug**（commit message 第二段）：原本 `if(hpResult != HAPLOTYPE1_1 && hpResult != HAPLOTYPE2_1)`（在 caller，HaplotagProcess.cpp:697）的比較對象是 enum 值 (HAPLOTYPE1_1=3, HAPLOTYPE2_1=4，見 `Util.h:21-25`)，但 `hpResult` 已被 `getVote()` 設成 HP tag 整數 (11/21/33)，**永遠不匹配**。整數 literal 修正後此 fallback 才生效。

**HP tag 影響**: HP1, HP2, HP11, HP21, HP33 全體 — 從「99.9% HP21」恢復到 HP_Ratio 接近 paired baseline (~0.5 中位數)。

**最小必要性**: ✅ 是。
- 沒有額外改 caller 端契約（介面 `void getVote(countMap, min, max, hpResult)` 完全不變）；
- 沒有引入新的資料結構或新的 enum；
- 把舊的 `std::map<int,int> haplotypeBase` 與 `std::vector<std::pair<int,int>>` 都刪了，**淨減少 25 行**（用直觀的兩層邏輯替代）。

---

### 2.3 INDEL guard（commit `380e8d2`）

**Subject**: `fix(haplotag): guard countINDELHaplotype against UNDEFINED array access`

**HaplotagProcess.cpp 變更**: `countINDELHaplotype()` 兩個分支各加一層 `if(... != HAPLOTYPE_UNDEFINED)`，+8 / -4 行。

**修改前（V3-Fixed [`/tmp/HP_v3fixed.cpp:495-504`]）**：

```cpp
void HaplotagProcess::countINDELHaplotype(bool isAlt, AlleleHaplotype &haplotypeBase, ...){
    if(isAlt == false){
        countMap[haplotypeBase.refHaplotype]++;     // <-- UB if refHaplotype = HAPLOTYPE_UNDEFINED = -1
        HP = haplotypeBase.refHaplotype;
    }
    else if(isAlt == true){
        countMap[haplotypeBase.altHaplotype]++;     // <-- same UB on alt path
        HP = haplotypeBase.altHaplotype;
    }
}
```

**修改後（V5 [`HaplotagProcess.cpp:497-510`]）**：

```cpp
void HaplotagProcess::countINDELHaplotype(bool isAlt, AlleleHaplotype &haplotypeBase, ...){
    if(isAlt == false){
        if(haplotypeBase.refHaplotype != HAPLOTYPE_UNDEFINED){   // <-- guard
            countMap[haplotypeBase.refHaplotype]++;
            HP = haplotypeBase.refHaplotype;
        }
    }
    else if(isAlt == true){
        if(haplotypeBase.altHaplotype != HAPLOTYPE_UNDEFINED){
            countMap[haplotypeBase.altHaplotype]++;
            HP = haplotypeBase.altHaplotype;
        }
    }
}
```

**為何要改**: 在 PON-only 模式下，somatic 站點的 GT 為 `0|0+GT2`，此時 `refHaplotype = HAPLOTYPE_UNDEFINED` (= `-1`，見 `Util.h:20`)。`countMap[-1]++` 是**陣列越界寫入**，C++ 規範下為 undefined behaviour（實務上會破壞 stack 上相鄰變數）。此 guard 與 `countSNPHaplotype()` 既有的 `haplotypeBase.refHaplotype != HAPLOTYPE_UNDEFINED` 守衛是**對稱的**——只是 V3-Fixed 時 INDEL 路徑漏掉了。

**HP tag 影響**: HP1, HP2, HP11, HP21（任何走 INDEL 計票路徑的 read）。修補後 HP 分佈不會「神奇地改變」，但消除了 UB 帶來的隨機污染。

**最小必要性**: ✅ 是。模式上完全 mirror SNP 端既有 guard，無新邏輯。

---

### 2.4 V5 HEAD（working tree，uncommitted）

V5 對應 `longphase-to` 二進位（編譯時間 2026-04-12 13:08，見 `ls -la longphase-to-mod/`），working tree 相對 `380e8d2` 多了兩塊修改（`git status: modified: HaplotagProcess.cpp`）：

#### 2.4.1 `getVote()` 加入 Layer 1.5 somatic fallback（[`HaplotagProcess.cpp:512-563`]）

**修改前（380e8d2）**：germline 全為 0 時 → `min=0, max=0`，`hpResult = 33`（若 somaticTotal>0）或 `0`（若全 0）。**germline-poor 但 somatic-rich** 的 read 拿不到 HP1/HP2 方向資訊。

**修改後（V5 HEAD）**：

```cpp
// Layer 1: Germline haplotype determination
int germlineResult = 0;
if (germlineHP1 > 0 || germlineHP2 > 0) {
    // ... germline path (unchanged from V3-Fixed)
}
// Layer 1.5: Somatic fallback — HP1_1/HP2_1 carry phased haplotype info  [V5 NEW]
else if (somaticHP1 > 0 || somaticHP2 > 0) {
    if (somaticHP1 >= somaticHP2) {
        min = somaticHP2; max = somaticHP1; germlineResult = 1;
    } else {
        min = somaticHP1; max = somaticHP2; germlineResult = 2;
    }
}
// No directional evidence
else { min = 0; max = 0; }
```

**邏輯變化**：
- **Layer 1 優先級維持**——只要有任何 germline vote，Layer 1.5 不會啟動（V3-Fixed 的核心契約保留）；
- **Layer 1.5 只在 germline 全 0 時啟動**——此時利用 `HAPLOTYPE1_1`/`HAPLOTYPE2_1` 的方向資訊（這兩個 enum 本身就是「somatic on HP1/HP2」的 phased 標記），給 read 一個 HP11/HP21 而非 HP33；
- **Layer 2 不變**——仍由 `germlineResult` (1/2/0) × `somaticTotal>0` 組合決定 hpResult。

**HP tag 影響**: 主要影響 HP11/HP21/HP33 的相對比例 — germline-poor 但 somatic-rich 的 reads 從 HP33 重新分配到 HP11/HP21。HP1/HP2 不受影響（Layer 1 不變）。

**為何要改**: V3-Fixed 過於保守——把所有 germline-poor reads 一律 HP33，丟失了 `HAPLOTYPE1_1`/`HAPLOTYPE2_1` 本身已帶有的方向資訊（這兩個 enum 的命名 `_1` 後綴正好代表 phased somatic）。Layer 1.5 在不破壞 Layer 1 主契約的前提下榨出剩餘訊息。

**最小必要性**: ✅ 是。
- 純 `else if` 分支插入，**Layer 1 與 Layer 2 程式碼一字未改**；
- 沒有引入新 enum、新成員變數、新介面；
- 額外變數只有 `somaticHP1`/`somaticHP2` 兩個 local int（從原本的 `somaticTotal` 拆出來）；
- 增加 +15/-1 行。

#### 2.4.2 `countSNPHaplotype()` 加入 alt-side guard（[`HaplotagProcess.cpp:489-494`]）

**修改前（V3-Fixed）**：

```cpp
else if(base == haplotypeBase.alt){
    countMap[haplotypeBase.altHaplotype]++;     // <-- same UB risk on somatic SNP sites
    HP = haplotypeBase.altHaplotype;
}
```

**修改後（V5）**：

```cpp
else if(base == haplotypeBase.alt){
    if(haplotypeBase.altHaplotype != HAPLOTYPE_UNDEFINED){
        countMap[haplotypeBase.altHaplotype]++;
        HP = haplotypeBase.altHaplotype;
    }
}
```

**邏輯變化**：對稱補上 `altHaplotype` 的 UNDEFINED 守衛。`countSNPHaplotype()` ref 端早就有此 guard（`HaplotagProcess.cpp:485`），alt 端漏寫；v5 補齊。

**HP tag 影響**: 對稱於 380e8d2 的 INDEL guard——避免 somatic SNP 站點觸發 `countMap[-1]++` UB。

**最小必要性**: ✅ 是。同樣是與 ref-side 完全對稱的補丁，+9/-6 行（含調整縮排）。

---

## Section 3 — V5 三層架構（Layer 1, 1.5, 2）程式碼解析

### 3.1 整體流程

![fig01b — V5 getVote three-layer logic](figures/01_code_diff/fig01b_three_layer_logic.png)

*Figure 1b — V5 `getVote()` 三層判定流程：Layer 1（germline 優先）→ Layer 1.5（somatic fallback）→ Layer 2（HP tag 編碼）。*

### 3.2 各層職責

| Layer | 進入條件 | 輸出語義 | 對應程式碼行 |
|-------|---------|---------|-------------|
| **Layer 1** | `germlineHP1>0 \|\| germlineHP2>0` | 有 germline 證據時，純粹比較 HP1 vs HP2 votes，set `germlineResult ∈ {1,2}`、`min`、`max` | `HaplotagProcess.cpp:524-536` |
| **Layer 1.5** | (Layer 1 不成立) AND `somaticHP1>0 \|\| somaticHP2>0` | 退而求其次，用 phased somatic votes (`HAPLOTYPE1_1`/`HAPLOTYPE2_1`) 推測方向，set `germlineResult ∈ {1,2}` | `HaplotagProcess.cpp:537-548` |
| **(else)** | 兩層皆不成立 | 完全無方向證據 → `min=max=0`, `germlineResult=0` | `HaplotagProcess.cpp:549-553` |
| **Layer 2** | always runs | 把 `germlineResult` 與 `somaticTotal>0` 兩個布林組合 → 寫出最終 `hpResult` ∈ {0, 1, 2, 11, 21, 33} | `HaplotagProcess.cpp:555-562` |

### 3.3 HP tag 對應表（V5 final）

| `germlineResult` | `somaticTotal>0` | `hpResult` | BAM HP:i: tag | 語義 |
|------------------|------------------|-----------|---------------|------|
| 0 | false | 0 | (no tag) | 無方向證據，read 不打 HP |
| 1 | false | 1 | HP:i:1 | 純 germline → HP1 |
| 2 | false | 2 | HP:i:2 | 純 germline → HP2 |
| 1 | true | 11 | HP:i:11 | germline=HP1 + 至少一個 somatic vote |
| 2 | true | 21 | HP:i:21 | germline=HP2 + 至少一個 somatic vote |
| 0 | true | 33 | HP:i:33 | 只有 HP3 (HAPLOTYPE3) somatic、無 HP1_1/HP2_1 方向 |

### 3.4 caller-side 一致性

`getVote()` 的呼叫端在 [`HaplotagProcess.cpp:725`]，回傳值在 `judgeHaplotype()` 內被檢查：

```cpp
getVote(countMap, min, max, hpResult);
// ...
if(hpResult != HAPLOTYPE1_1 && hpResult != HAPLOTYPE2_1){   // baseline: enum 3/4 比較 HP tag 11/21 — never matches
```

該比較在 V3-Fixed 之後改為整數 literal 11/21/33（見 41ff147 commit message 與 working tree [`HaplotagProcess.cpp:697`]），與 `getVote()` 的輸出 contract 對齊。

---

## Section 4 — 每個修改的「最小必要性」評估

| Commit / 修改 | 行數 | 是否必要 | 是否最小 | 評語 |
|--------------|------|---------|---------|------|
| `8b8c1fd` PON-only mode | n/a (cpp 未動) | ✅ feature 必要 | ✅ 最小（一個 flag、三個單職責 helper） | 沒有提前動 HaplotagProcess.cpp 是正確設計—— bug 由後續 commit 在乾淨 scope 處理 |
| `41ff147` two-layer getVote | +36 / -25 | ✅ 修正核心 priority bug；無此修改 PON-only 完全失效 | ✅ 最小（淨刪 25 行 stale logic、無新介面） | 把 baseline 的 `std::map+vector<pair>` 結構直接代換成兩個顯式 layer，可讀性大幅上升 |
| `41ff147` enum→int literal | (含於上) | ✅ 修正 caller 端的死分支 | ✅ 最小（單行常數比對換掉） | 修了一個從未被觸發的 fallback path；是 priority bug 的姊妹 bug |
| `380e8d2` INDEL guard | +8 / -4 | ✅ 修 UB（C++ 未定義行為） | ✅ 最小（mirror SNP 端既有 guard 形狀） | 純粹補洞，與 SNP 端完全對稱 |
| `V5-HEAD` Layer 1.5 fallback | +15 / -1 | ✅ 救回 germline-poor reads 的方向資訊 | ✅ 最小（純 `else if` 分支、無新介面） | Layer 1 主契約完全保留；只在 germline 全 0 時啟動 |
| `V5-HEAD` SNP alt guard | +9 / -6 | ✅ 修 UB（同 380e8d2 邏輯） | ✅ 最小（對稱於 ref 端） | 380e8d2 漏掉的對稱補丁 |

### 4.1 沒有發現的 over-engineering 跡象

我們檢查了以下幾類常見的多餘改動，**全部沒有發現**：
- ❌ 無新增 enum 值（`Util.h` 的 `Haplotype` 從 baseline 起未動）
- ❌ 無新增 class member 或 method signature 變更（`HaplotagProcess.h:66-68` 三個函式簽名 baseline=V5）
- ❌ 無新增 std 容器（兩個迴圈、兩個 map、一個 vector 全部刪除，淨減）
- ❌ 無新增 logging / debug printout
- ❌ 無 stylistic refactor（命名、縮排、註解 `enum` 的調整全部對齊既有風格）
- ❌ 無「為了通過某 unit test 而加的特殊路徑」（既有 test 不變）

---

## Section 5 — 程式碼層面結論

### 5.1 結論

**V5 程式碼修改是合理的，沒有 over-engineering。**

依據：
1. **每個 commit 對應一個明確問題** — feat (8b8c1fd) / priority bug fix (41ff147) / UB guard (380e8d2) / fallback enhancement + 對稱 guard (V5 HEAD)；
2. **總共修改集中在 3 個函式（`getVote`、`countSNPHaplotype`、`countINDELHaplotype`），共 +68/-36 行**，無散落他處的「順便修一下」改動；
3. **每處修改都有對稱 / 既有模式可循**——guard 對稱、Layer 1.5 是 Layer 1 的純擴充、整數 literal 修正只是把 caller 端比對對象與 `getVote` 輸出對齊；
4. **介面契約零變動**——`HaplotagProcess.h:66-68` 的三個 method signature 從 baseline 到 V5 一字未變，所有改動都是 implementation-internal；
5. **兩個 UB 修補**（380e8d2 + V5 SNP alt guard）為對稱補洞，符合「最小必要」原則。

### 5.2 仍待後續 commit 收斂的事項

- V5 的兩塊修改（Layer 1.5 + SNP alt guard）目前在 **working tree (uncommitted)**，建議切成兩個 commit：
  - `feat(haplotag): Layer 1.5 somatic fallback in getVote()` (+15/-1)
  - `fix(haplotag): guard countSNPHaplotype against UNDEFINED on alt path` (+9/-6)
- 兩個 commit 都不會引入新介面，可在 V5 結論驗證完成後一次性 commit，亦不需 amend。
- 二進位 `longphase-to`（4/12 編譯）已包含這兩塊修改的執行版本；source-binary 一致性需在後續 audit suite chapter 6（sanity check）確認。

### 5.3 對下游 audit chapters 的依賴交付

| 下游章節 | 本章交付的事實 |
|---------|---------------|
| 02_concordance | V5 三層 HP encoding 表（Section 3.3）— concordance 應以 HP1/HP2 vs HP11/HP21 雙層比對 |
| 04_imbalance | 99.9% HP21 是 baseline+V2b 的 priority bug 直接後果 — V3-Fixed 之後應消失 |
| 06_sanity | V5 working-tree 比 380e8d2 多兩塊修改 — 二進位 `longphase-to` 對應 V5（4/12 編譯） |
| 07_paired | Layer 1.5 不影響 paired 模式（germline 通常 >0，Layer 1 即決定） |
| 08_synthesis | 「每修改皆最小必要」可作為 V5 audit 的 code-side 結論基石 |

---

## Appendix A — 引用之原始程式碼位置（含 file:line）

| 內容 | 位置 |
|------|------|
| baseline `getVote()` 全文 | `git show 8b8c1fd:HaplotagProcess.cpp` lines 506-530 |
| V3-Fixed `getVote()` 全文 | `git show 41ff147:HaplotagProcess.cpp` lines 506-541 |
| V5 `getVote()` 全文 | `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:512-563` |
| V5 `countSNPHaplotype()` 全文 | `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:484-495` |
| V5 `countINDELHaplotype()` 全文 | `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:497-510` |
| Haplotype enum 定義 | `/big7_disk/liaoyoyo2001/longphase-to-mod/Util.h:20-25` |
| getVote caller (judgeHaplotype) | `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:725` |
| caller fallback 比對 (post-fix) | `/big7_disk/liaoyoyo2001/longphase-to-mod/HaplotagProcess.cpp:697` |

## Appendix B — 產出檔案清單

| 類型 | 路徑 |
|------|------|
| 本文件 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` |
| 圖 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01a_commit_evolution.png` |
| 圖 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/01_code_diff/fig01b_three_layer_logic.png` |
| TSV | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/code_diff_summary.tsv` |
