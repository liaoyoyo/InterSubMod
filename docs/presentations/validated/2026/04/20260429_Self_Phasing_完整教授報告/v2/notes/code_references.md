# 程式碼參考（commit hash + 檔案行號對照）

> 完整 4-commit 演進的精確 code references。所有路徑以 longphase-to-mod fork（**獨立 git repo**，非 InterSubMod）為準。

## 1. Repo 位置與切分

| Repo | 絕對路徑 | 與本工作關係 |
|------|---------|-------------|
| longphase-to-mod fork | `/big7_disk/liaoyoyo2001/longphase-to-mod/` | **4-commit 修補在此**（獨立 git repo）|
| InterSubMod 本 repo | `/big7_disk/liaoyoyo2001/InterSubMod/` | **無 C++ 改動**（src/core/ 不含 HaplotagProcess.cpp）|

## 2. 4-commit 漸進演進

| # | Commit Hash | 名稱 | 解的 bug | 主要檔案修改 |
|---|-------------|------|---------|--------------|
| 1 | **`8b8c1fd`** | V2b PON-only | Phase 層 self-phasing scaffold | Phasing.cpp (+9/-2)、PhasingGraph.cpp (+34/-0)、PhasingProcess.cpp (+25/-3)；HaplotagProcess.cpp **未動** |
| 2 | **`41ff147`** | V3-Fixed | Tag 層 priority bug + enum mismatch | HaplotagProcess.cpp:506-541（getVote 重寫）、:697（caller 端）；+36/-25 |
| 3 | **`380e8d2`** | INDEL guard | UB 補洞（countINDELHaplotype）| HaplotagProcess.cpp:497-510；+8/-4 |
| 4 | **(working tree, 未 commit)** | V5 | V3F 過於保守 | HaplotagProcess.cpp:489-494（countSNPHaplotype alt guard）、:512-563（getVote Layer 1.5 + 三層）；+24/-7 |

## 3. 關鍵檔案與行號

### 3.1 HaplotagProcess.cpp（修補核心）

| 行號 | 函數 | V5 邏輯 | 對應 commit |
|------|------|---------|------------|
| **66-68**（.h）| 三函數 method signature | **介面契約零變動**（baseline → V5 一字未變）| 全部 commits |
| 489-494 | `countSNPHaplotype()` | V5 加 `if (haplotypeBase.altHaplotype != HAPLOTYPE_UNDEFINED)` 對稱 alt guard | V5 working tree |
| 497-510 | `countINDELHaplotype()` | INDEL guard（兩分支補 `HAPLOTYPE_UNDEFINED` guard）| `380e8d2` |
| **512-563** | **`getVote()`** | **V5 三層投票邏輯（核心）** | V3F + V5 working tree |
| 557-559 | hpResult assignment | 直接賦值 integer 11/21/33（V3F 修為 integer literal）| `41ff147` |
| **697** | `judgeHaplotype()` caller 端 | V3F 改為 `if(hpResult != 11)` integer literal 比對；fallback 死分支修復；HP:i:33 開始正確出現 | `41ff147` |
| 725 | `getVote()` 在 `judgeHaplotype()` 內的 caller | - | - |

### 3.2 V5 `getVote()` 三層邏輯（行 512-563）

```
Layer 1 (germline first, 行 ~515-525):
    if (countMap[HAPLOTYPE1] > countMap[HAPLOTYPE2])  germlineResult = 1
    elif (countMap[HAPLOTYPE2] > countMap[HAPLOTYPE1]) germlineResult = 2
    else                                              germlineResult = 0

Layer 1.5 (somatic fallback, V5 新增, 行 ~528-545):
    somaticHP1 = countMap[HAPLOTYPE1_1]
    somaticHP2 = countMap[HAPLOTYPE2_1]
    somaticTotal = somaticHP1 + somaticHP2
    if (germlineResult == 0 and somaticTotal > 0):
        if (somaticHP1 > 0.6 * somaticTotal) → 對齊 HP1 方向 (germlineResult = 1)
        elif (somaticHP2 > 0.6 * somaticTotal) → 對齊 HP2 方向 (germlineResult = 2)

Layer 2 (encode hpResult, 行 ~548-563):
    if (germlineResult == 1 and somaticTotal > 0)  hpResult = 11   // somatic on HP1
    elif (germlineResult == 2 and somaticTotal > 0) hpResult = 21   // somatic on HP2
    elif (germlineResult == 0 and somaticTotal > 0) hpResult = 33   // ambiguous
    elif (germlineResult == 1)                     hpResult = 1    // germline HP1
    elif (germlineResult == 2)                     hpResult = 2    // germline HP2
    else                                            hpResult = 0    // unphased
```

### 3.3 Util.h（enum 定義）

| 行號 | 內容 |
|------|------|
| **20-25**（**注意：不是 21-25**）| `Haplotype` enum 定義：`HAPLOTYPE_UNDEFINED=-1, HAPLOTYPE1=1, HAPLOTYPE2=2, HAPLOTYPE1_1=3, HAPLOTYPE2_1=4, HAPLOTYPE3=5` |

**enum vs HP tag integer mapping**:
- enum HAPLOTYPE1 = 1 ↔ HP:i:1
- enum HAPLOTYPE2 = 2 ↔ HP:i:2
- enum HAPLOTYPE1_1 = 3 ↔ HP:i:11（**型別語意失配 bug 來源**）
- enum HAPLOTYPE2_1 = 4 ↔ HP:i:21
- enum HAPLOTYPE3 = 5 ↔ HP:i:33

### 3.4 PhasingProcess.cpp（V2b PON-only + Purity threshold）

| 行號 | 內容 |
|------|------|
| 55 | PON 設定 |
| 154-157 | `convertNonGermlineToSomatic()` 觸發點（V2b 新增邏輯）|
| **197** | **★ Purity threshold 硬編碼 `purity > 0.95`**（決定是否啟動 Two-Pass）|

**Purity 0.95 邏輯（v3 PPT S11 核心）**：

```cpp
// PhasingProcess.cpp:197 (baseline)
if (purity > 0.95) {
    // 兩條路 / Two-Pass 路徑
    vGraph->convertNonGermlineToSomatic();
    vGraph->phasingProcess(..., nullptr);
}
else {
    // 三條路 / Baseline 標準流程
    vGraph->somaticCalling(...);
    vGraph->phasingProcess(..., &chrInfo.ploidyRatioMap);
}
```

HCC1395 5kHz：真實 purity > 0.95 但 baseline 估 0.927（四捨五入 0.93）→ ≤ 0.95 → **未觸發 Two-Pass** → 走三條路 → 暴露 tag 層 bug。

### 3.5 PhasingGraph.cpp（V5 PON-only Purity 修復，2026-04-29）

| 行號 | 內容 |
|------|------|
| **1147-1175** | **★ `collectPloidyRatio()` 函式（新增，修復 V5 PON-only purity=0 bug）**|
| 1079-1083 | `calculatePloidyRatio()` 既有；只在 ploidyRatioMap≠nullptr 時 fill |
| 1085-1109 | `reassignAlleleResult()` 既有；fill ploidyRatioMap 的入口 |
| 1886-1903 | `getPurity()` polynomial 公式（q1/q3 主要輸入）|

**修復前 V5 PON-only**：`phasingProcess(..., nullptr)` → ploidyRatioMap 永不 fill → q1=q3=0 → polynomial = -6.81 - 5.72L - 0.51L² < 0 → clamp 至 **0**

**修復後**：在 `syncPhasingResultOrigins()` 後呼叫 `collectPloidyRatio()` → ploidyRatioMap 正確填入 → purity 計算正確（0.93 sample → 0.871；0.6 sample → 0.634）

### 3.5 InterSubMod 端（下游消費，無修改）

| 檔案 | 行號 | 內容 |
|------|------|------|
| `InterSubMod/src/core/ReadParser.cpp` | line 120-130 | HP tag parsing：`HP:i:1→"1"`、`11→"1-1"`、`21→"2-1"`、`33→"3"` |
| `InterSubMod/src/core/ReadParser.cpp` | (`--germline-hp-only` flag) | 開啟時 `"1-1"/"2-1"/"3"` demote 為 `"0"`；`hp_tag_raw` 保留供 audit |
| `InterSubMod/src/core/RegionProcessor.cpp` | line 1120, 1269 | ISM SuggestFilter（V5 對 F1 的影響進入這裡）|
| `InterSubMod/include/core/Stats.hpp` | line 323 | `d_hp1_hp1s` 註解："Same haplotype, germline vs somatic" |

## 4. 重現指令

### 4.1 longphase-to-mod 重 build

```bash
cd /big7_disk/liaoyoyo2001/longphase-to-mod/
git log --oneline | head -10   # 確認 4-commit 演進
git diff 41ff147..HEAD HaplotagProcess.cpp  # 看 V5 working tree 修改
make -j$(nproc)
```

### 4.2 InterSubMod 重 build（無 C++ 改動，但需 link V5 binary 重跑 ISM）

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod/
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
./scripts/run_vcf_all_snv.sh --mode all-with-w5000  # 確認 SEQC2 F1 不退步
```

### 4.3 V5 sanity check 重跑

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod/
python3 scripts/analysis/v5_imbalance_improvement.py     # Agent C 產出
python3 scripts/analysis/v5_sanity_paired_check.py       # Agent D 產出
# 期望輸出：4 項 sanity check 15/15 PASS, 0 violation
```

## 5. 短期 P0 行動的程式碼層具體步驟

### F1 — commit V5 working-tree 修改（高優先）

建議切 2 獨立 commits：

```bash
cd /big7_disk/liaoyoyo2001/longphase-to-mod/
# Commit 1
git add HaplotagProcess.cpp  # 只 stage Layer 1.5 部分（行 512-563）
git commit -m "feat(haplotag): Layer 1.5 somatic fallback in getVote()"
# Commit 2
git add HaplotagProcess.cpp  # 只 stage countSNPHaplotype alt guard 部分（行 489-494）
git commit -m "fix(haplotag): guard countSNPHaplotype against UNDEFINED on alt path"
```

### F2 — Confidence threshold 0.6 vote log 驗證

需修改 `getVote()` 加入 vote log 輸出：

```cpp
// 在 Layer 1.5 加入：
if (verbose) {
    std::cout << "[V5_VOTE] pos=" << pos
              << " germline=" << germlineResult
              << " somaticHP1=" << somaticHP1
              << " somaticHP2=" << somaticHP2
              << " hpResult=" << hpResult << "\n";
}
```

或產出 debug TSV：每位點 vote count 與最終 hpResult 對照。

### F3 — 7 樣本 V5 全量重跑

樣本清單（KDE-corrected）：
- HCC1395_5kHz（已驗證）
- HCC1395_DORADO（待）
- HCC1937、HCC1954、H1437、H2009、COLO829（待）

預估 ~10 hr parallel（依 archive TO rerun 規劃）。

## 6. ✅ Fact-check 確認

| 項目 | 確認狀態 |
|------|---------|
| Commit hash V2b `8b8c1fd` | ✅ 與 source 03 §10 §13 一致 |
| Commit hash V3F `41ff147` | ✅ 與 source 03 §10 §13 一致 |
| Commit hash INDEL guard `380e8d2` | ✅ 與 source 03 §10 §13 一致 |
| V5 working tree uncommitted | ✅ 與 source 03 §10 caveat 一致 |
| HaplotagProcess.cpp:512-563（V5 getVote）| ✅ 與 source 03 附錄 A.2 一致 |
| HaplotagProcess.h:66-68（介面契約零變動）| ✅ 與 source 03 §6, §10 一致 |
| Util.h:20-25 enum 範圍 | ✅ **以本檔為準**；storyboard 第一稿誤寫 21-25 已更正 |
| InterSubMod src/core/ 無 phasing/haplotag 檔 | ✅ 與 source 03 §10 確認 |
| **V5 修改量 +68/-36 行** | ✅ 與 source 03 §7.2 一致 |
| **3 函式集中**（getVote / countSNPHaplotype / countINDELHaplotype）| ✅ 一致 |

來源：`v5_audit_suite/01_code_diff_analysis.md`、`source_materials/03_longphase_TO_vs_V5_技術報告.md` §10, §13, 附錄 A.2。
