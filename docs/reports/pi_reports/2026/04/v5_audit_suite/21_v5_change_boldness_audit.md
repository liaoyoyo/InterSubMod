---
title: V5 改動大膽程度審計 — 必須改 vs 大膽改 vs 過度修改
date: 2026-04-29
author: liaoyoyo2001
tags: [audit, v5, baseline, change-classification, conservativeness]
status: validated_complete
audience: developer + PI + reviewers
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/16_baseline_subgenotype_clarification.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/17_design_consistency_check.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/18_purity_calculator_failure_root_cause.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/20_whole_genome_paired_audit.md
---

# V5 改動大膽程度審計：必須改 vs 大膽改 vs 過度修改

## §0 一句話結論

> **V5 的 7 個改動中，4 個是「必須改」（明確 bug + 防禦性編程，無爭議）**：getVote 順序反轉、HP:i:33 enum、UNDEFINED guards、collectPloidyRatio。**3 個是「大膽改」（改變 phasing 行為，需要謹慎評估）**：`--pon-only-phasing` flag、`syncPhasingResultOrigins`、Conservative HP33 tagging。**過度修改的關鍵警訊**：低 purity 下 V5 把 11% reads 推到 HP33，使 V5_06 在 aggregate level 比 BL_06 偏離 PA truth +192%。

---

## §1 V5 改動完整清單

### 1.1 7 個改動

| # | 改動 | Commit / Source | 性質 | 影響範圍 |
|:-:|------|:----:|------|------|
| 1 | `getVote()` Layer 1.5 + 順序反轉 | 41ff147 | bug fix | tag 階段 |
| 2 | HP:i:33 enum 比對修正 | 41ff147 | type bug fix | tag 階段 |
| 3 | `countINDELHaplotype` UNDEFINED guard | 380e8d2 | defensive | tag 階段 |
| 4 | `countSNPHaplotype` UNDEFINED guard | working tree | defensive | tag 階段 |
| 5 | `--pon-only-phasing` flag + `convertNonGermlineToSomatic` | 8b8c1fd | **新功能** | phase 階段 |
| 6 | `resetNonPonOrigin` + `syncPhasingResultOrigins` | 8b8c1fd | **state mutation** | phase 結果 |
| 7 | `collectPloidyRatio` | 本 session（未 commit）| completeness fix | purity calculation |

### 1.2 程式碼行數統計

```
HaplotagProcess.cpp | 36 +++ -9 (V5 working tree)  ← #1, #2, #3, #4
PhasingGraph.cpp    | 29 +++ (本 session 新增)       ← #7
PhasingGraph.h      |  6 +++                        ← #7 declaration
PhasingProcess.cpp  | 27 +++ -7 (跨 commits)         ← #5, #6, #7 call site
總計：68 insertions, 9 deletions（不含 8b8c1fd 引入的程式碼）
```

---

## §2 🟢 必須改的改動（4 項，無爭議）

### 2.1 改動 #1：`getVote()` Layer 1.5 + 順序反轉

**問題**：baseline `getVote()` 用 `for-loop` priority list：
```cpp
variantKeys = {
    {HAPLOTYPE1_1, HAPLOTYPE2_1},  // ← somatic first（錯！）
    {HAPLOTYPE3, HAPLOTYPE2_1},
    {HAPLOTYPE1, HAPLOTYPE2}        // ← germline 反而是備胎
};
```

**為何必須改**：
- ✅ commit 41ff147 message 直接記載 root cause：「99.9% of reads to get HP:i:21 in PON-only mode」
- ✅ phase.md 規範明確指出「GT (Germline and Primary Somatic Phasing)」— germline 是 primary
- ✅ baseline 順序違反規範意圖

**評估**：✅ **明確 bug**，必須改。爭議性 = 0。

---

### 2.2 改動 #2：HP:i:33 enum 比對修正

**問題**：
```cpp
// Baseline (錯)：
if (hpResult != HAPLOTYPE1_1)  // enum value = 3
// 但 hpResult 實際是 HP integer tag (1, 2, 11, 21, 33)
// → 永遠不 match，邏輯失效
```

**為何必須改**：
- ✅ 明確 type confusion bug
- ✅ enum value (3) 與 HP tag int (33) 比對毫無意義

**評估**：✅ **明確 type bug**，必須改。爭議性 = 0。

---

### 2.3 改動 #3 + #4：UNDEFINED guards

**問題**：enum `HAPLOTYPE_UNDEFINED = -1` 是合法狀態（無 phasing 資訊），但 baseline `countSNP/INDELHaplotype` 不檢查：
```cpp
// Baseline (錯)：
if (base == haplotypeBase.ref) {
    countMap[haplotypeBase.refHaplotype]++;  // ← 若 refHaplotype = UNDEFINED?
    // → countMap[-1] 寫入未定義位置 → silent corruption
}
```

V5 修法：
```cpp
if (base == haplotypeBase.ref &&
    haplotypeBase.refHaplotype != HAPLOTYPE_UNDEFINED) {  // ← Guard
    countMap[haplotypeBase.refHaplotype]++;
}
```

**為何必須改**：
- ✅ 防禦性編程標準做法
- ✅ enum 已含 UNDEFINED，不檢查是漏洞
- ✅ Silent corruption 比 explicit fail 更危險

**評估**：✅ **defensive programming 標準**，必須加。爭議性 = 0。

---

### 2.4 改動 #7：`collectPloidyRatio()`

**問題**：V5 PON-only Pass 2（`somaticCalling` + `syncPhasingResultOrigins`）後沒人補 fill `ploidyRatioMap` → `getPurity()` 收到空 map → q1=q3=0 → polynomial 必負 → clamp 至 0。

**為何必須改**：
- ✅ 本 session 已 root cause 確認（18 號報告）
- ✅ 雙 sample 驗證：0.93 從 0 → 0.871，0.6 從 0 → 0.634（修復成功）
- ✅ Side-effect free（不影響 GT/PS 標記，不需重跑 BAM）

**評估**：✅ **完整性修復**，必須改。爭議性 = 0。

---

## §3 ⚠ 大膽改動（3 項，有爭議）

### 3.1 改動 #5：`--pon-only-phasing` flag + `convertNonGermlineToSomatic`

**改動內容**：強制把所有非 PON variants 標為 SOMATIC origin，限制 phasing graph 只用 PON germline 當 anchor。

**為何大膽**：

| 爭議點 | 詳情 |
|--------|------|
| **改變 phasing 算法輸入分布** | 原 LongPhase-TO 接受 ClairS-TO PASS somatic 進 graph；V5 強制排除 |
| **違背 Two-Pass 條件啟用設計** | 原設計只在 purity > 0.95 啟用 `convertNonGermlineToSomatic`；V5 把它變成 user-controlled flag（可在任何 purity 下開啟）|
| **可能改變 phase block N50** | 限制 anchor 來源 → 部分 PS block 變短/斷開；論文 §4.1 報的 N50 10-25Mb 在 V5 PON-only mode 是否仍成立**未驗證** |
| **WG 數據顯示 V5 在 0.6 下偏離 PA** | aggregate HP1:HP2 距離 PA 從 0.124 (BL_06) → 0.362 (V5_06) = **+192%** |
| **purity calculator 副作用** | 直到本 session 才發現需要 `collectPloidyRatio` 補完整 |

**證據鏈（從 audit suite）**：
- 17 號 §3.1：原 Two-Pass 機制只在 purity>0.95 啟用
- 18 號 §3：V5 強制啟用導致 purity calculator 失效（已修）
- 20 號 §2.1：V5_06 距離 PA truth +192%

**評估**：⚠ **保留 flag 但默認關閉（已是這樣）**；建議文件強烈建議「purity ≥ 0.85 才開啟」。

**質疑（待驗證）**：
1. ❓ 是否該保留原 Two-Pass 條件啟用邏輯（purity > 0.95），而不是強制 always-on？
2. ❓ 在低 purity (< 0.7) 場景下 V5 反偏離 PA truth 是 V5 設計缺陷還是 metric bias？
3. ❓ 論文 §4.1 phase quality benchmark 在 PON-only mode 下的 N50 / Phased Ratio 是否仍成立？

---

### 3.2 改動 #6：`resetNonPonOrigin` + `syncPhasingResultOrigins`

**改動內容**：Pass 1 phasing 後，把暫時標 SOMATIC 但實際是 germline 的 variants origin reset，再用 `syncPhasingResultOrigins` 校正 GT format。

**為何大膽**：

| 爭議點 | 詳情 |
|--------|------|
| **State mutation on phasing result** | 在已 phase 過的結果上做 post-hoc origin/GT 修改 |
| **可能丟失 phasing 信息** | `syncPhasingResultOrigins` 內含 `posPhasingResult.erase()` 邏輯（PhasingGraph.cpp:1162-1180），這個移除讓某些 variants 變 unphased |
| **與 Pass 1 phasing 結果不一致** | Pass 1 把 variants A 標 SOMATIC + GT=`0\|0`；Pass 2 reset 為 ORIGIN_UNDEFINED → variants A 可能變成 unphased (`0/1`) |
| **未 quantify removal 數量** | 沒有 logging 記錄「sync 步驟移除了多少 variants」 |

**證據**：
- PhasingGraph.cpp:1155-1180 函式註解直接寫「Variants with HAPLOTYPE_UNDEFINED refHaplotype cannot be converted to germline format, so they are removed」

**評估**：⚠ **必要但需 logging**。建議：
1. 加 logging 記錄 sync 移除的 variants 數量
2. 考慮重設計成「Pass 1 不誤標 origin」而非「Pass 2 sync」

**質疑（待驗證）**：
1. ❓ `syncPhasingResultOrigins` 的 erase 是否丟失重要 phasing 信息？
2. ❓ 是否該在 Pass 1 用 origin filter，避免 Pass 2 sync？
3. ❓ 移除的 variants 在 baseline 中是 phased 還是 unphased？（可能 V5 erase 掉 baseline 已 phase 的 variants）

---

### 3.3 改動 #6.5：Conservative HP33 tagging 在低 purity 的副作用

**現象**：

| Sample | HP33% |
|--------|:-:|
| BL_93 | 0.8% |
| V5_93 | 5.0% |
| BL_06 | 1.7% |
| **V5_06** | **11%** |

V5 在 0.6 下 HP33 比例飆到 **11%**，是 BL_06 的 6.5 倍。

**為何大膽**：

| 爭議點 | 詳情 |
|--------|------|
| **HP33 在 paired metric 被當「不 match」** | PA_93 中 HP33 = 0；任何 V5 標 HP33 的 site 在 site-level concordance 算錯（懲罰 V5）|
| **低 purity 下 over-trigger Layer 1.5** | germline anchor 弱 → Layer 1.5 fallback 啟動更頻繁 → somatic 投票分散 → 多數歸 HP33 |
| **下游無法區分「V5 conservative」與「真錯誤」** | HP33 同時意味「V5 設計保守」與「V5 不確定」，下游處理需要 context |

**評估**：🔴 **可能過度大膽**。

**質疑**：
1. ❓ 是否該加 **purity-aware HP33 threshold**？（高 purity 下嚴格，低 purity 下放寬把 HP33 推回 HP1/HP2）
2. ❓ HP33 的 11% 是否該被視為「V5 在 0.6 下 over-correct」的訊號？

---

## §4 改動爭議性矩陣

### 4.1 必要性 × 爭議性矩陣

| 改動 | 必要性 | 爭議性 | 是否該保留 | 待補強行動 |
|------|:-:|:-:|:-:|------|
| #1 getVote 順序反轉 | 🟢 必須 | 低 | ✅ 保留 | — |
| #2 HP:i:33 enum | 🟢 必須 | 低 | ✅ 保留 | — |
| #3 INDEL guard | 🟢 必須 | 低 | ✅ 保留 | — |
| #4 SNP guard | 🟢 必須 | 低 | ✅ 保留 | — |
| #5 PON-only flag | 🟡 conditional | **中-高** | ⚠ 保留 + 默認 off | 文件加 purity ≥ 0.85 推薦 |
| #6 sync/reset Origins | 🟡 conditional | **中** | ⚠ 加 logging | logging + 後續審視 |
| #6.5 Conservative HP33 | 🟡 conditional | **中** | ⚠ 評估 purity-aware | low-purity HP33 threshold 設計 |
| #7 collectPloidyRatio | 🟢 必須 | 低 | ✅ 保留 | commit |

### 4.2 「大膽程度」與證據對應

| 大膽改動 | 直接證據 | 可信度 |
|---------|------|:----:|
| PON-only flag 改變 phasing | 18 號 §3 (purity calculator 失效)、20 號 §2 (V5_06 偏離 PA) | 高 |
| Sync 可能丟失 phasing 信息 | PhasingGraph.cpp:1162 註解明確記載 erase | 高 |
| HP33 過度 conservative | 20 號 §2.3 (V5_06 HP33 11%) | 高 |

---

## §5 過度修改的關鍵警訊

### 5.1 V5 在 0.6 下的 5 個警訊

1. **Aggregate ratio 偏離 PA**：V5_06 距離 PA = 0.362 vs BL_06 0.124 (+192%)
2. **HP33 飆到 11%**（BL_06 1.7%）
3. **Untagged 從 5% → 7%**（V5 多 2pp untagged）
4. **Site-level concordance match%**：V5_06 49.25% vs BL_06 51.81%（V5 −2.56pp）
5. **purity calculator 失效**（已修，但揭露 PON-only flow 設計不完整）

### 5.2 是否真的「過度」？

| 觀點 | 解讀 |
|------|------|
| **支持「過度」** | V5 在低 purity 下犧牲 paired-truth alignment 換取 conservative tagging，但 baseline 在 0.6 自然就接近 1:1，V5 修法**邊際效用低** |
| **反對「過度」** | V5 conservative tagging 是設計意圖（feature, not bug）；HP33 顯式標出可疑 reads 對下游有價值；用 PA_93 (0.93 truth) 評估 0.6 BAM 是 cross-purity bias |

→ **結論**：在「以 PA_93 為 truth」框架下 V5 在 0.6 過度；但若認可「conservative tagging = feature」，則不過度。**視下游使用方式決定**。

---

## §6 推薦行動

### 6.1 立即處理（高優先）

| 行動 | 說明 |
|------|------|
| **commit collectPloidyRatio 修復** | 本 session 新增（4 files modified）尚未 commit |
| **更新 README** | 注明 `--pon-only-phasing` 推薦 purity ≥ 0.85 使用 |
| **加 logging 給 syncPhasingResultOrigins** | 記錄移除的 variants 數量 |

### 6.2 中期審視（建議實驗）

| 行動 | 說明 |
|------|------|
| **0.6 paired truth 補強實驗** | 對 t30_n20 跑 paired-mode longphase 取真 0.6 paired truth |
| **HP33 conservative-aware concordance** | 計算 V5 排除 HP33 後的 conditional match% |
| **0.6 caller-level F1** | 對 0.6 BAM 做 ClairS-TO calling F1（vs SEQC2 truth）|
| **論文 §4.1 phase quality 重驗證** | Block N50 是否仍 10-25 Mb？Phased Ratio 是否仍 0.55-0.62？ |

### 6.3 長期重構（建議）

| 行動 | 說明 |
|------|------|
| **重新設計 Pass 1 origin tagging** | 避免 Pass 2 需要 sync |
| **purity-aware adaptive logic** | 低 purity 下自動放寬 conservative，避免 over-tagging HP33 |
| **重訓練 polynomial regression** | 加入 PON-only mode 的 q1/q3/lohRatio 分布作為 caller-specific 子分支 |

---

## §7 V5 整體評估（基於本審計）

### 7.1 改動價值評估

| 維度 | 評估 |
|------|:--:|
| **修復明確 bugs** | ✅✅ 4/4 必要改動全部正確 |
| **設計擴展邏輯** | ⚠ PON-only flag 從 conditional 改 always-on 是設計爭議 |
| **State mutation 安全性** | ⚠ Sync 邏輯需 logging 補強 |
| **跨 purity robustness** | 🔴 0.6 下 V5 conservative 副作用顯著（HP33 11%）|
| **文件完整性** | ✅ 經本 audit suite 21 份報告補強至完整 |

### 7.2 V5 整體分類

V5 = **「4 個明確 bug 修復」 + 「3 個有爭議的設計擴展」**：
- 4 個 bug fix 應**永遠保留**
- 3 個爭議改動應**保留但加 caveat 文件**

### 7.3 建議的 V5+ 改進方向

```
V5+:
  ✓ 保留 V5 所有 4 個 bug fix
  ✓ 保留 PON-only flag 但默認 off
  + Add logging for syncPhasingResultOrigins removal count
  + Add purity-aware HP33 threshold (low purity 下放寬)
  + Document recommended use cases:
    - purity ≥ 0.85: --pon-only-phasing 強烈推薦
    - 0.6-0.85: 兩者皆可，依下游需求
    - < 0.5: 慎用，HP33% 可能 > 15%
```

---

## §8 結論

### 8.1 直接答覆

> **「V5 改動是否過於大膽？」**

✅ **4/7 改動是必須改（無爭議）**：
- getVote 順序反轉、HP:i:33 enum、UNDEFINED guards、collectPloidyRatio

⚠ **3/7 改動是大膽改（需 caveat）**：
- PON-only flag、syncPhasingResultOrigins、Conservative HP33 tagging

🔴 **0/7 改動需要回退**：
- 無「絕對應該回退」的改動，但「PON-only flag + HP33」在低 purity 下需謹慎使用

### 8.2 V5 在不同場景的可信度

| Purity | V5 vs BL 推薦 | 理由 |
|--------|:-:|------|
| ≥ 0.85 | ✅ V5 | self-phasing 顯著，V5 修復價值高 |
| 0.6-0.85 | ⚠ 兩者皆可 | self-phasing 中等，V5 conservative 副作用浮現 |
| < 0.5 | ⚠ V5 慎用 | HP33% 可能 > 15%，過度 conservative |

### 8.3 對 audit suite 既有結論的影響

本檔不否定既有結論，但**精確化界定 V5 改善的場景與 caveat**：

| 既有結論 | 與本審計關係 |
|---------|---------|
| V5 修復 17.3:1 bias | ✅ 在 0.93 場景成立，0.6 邊際 |
| V5 不傷 ClairS-TO calling F1 | ✅ 完全成立 |
| V5 sanity check 全 PASS | ✅ 完全成立 |
| Aggregate paired +6.65pp | ⚠ conditional only；WG fixed-denom 持平 |
| V5 conservative tagging | ⚠ 低 purity 下副作用顯著 |

---

## §9 跨檔索引

### 對應 audit suite

| # | 文件 | 與本檔關係 |
|---|------|---------|
| 01 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` | 4 commit diff 細節 |
| 11 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md` | 12 issues 分類 |
| 16 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/16_baseline_subgenotype_clarification.md` | baseline GT2/GT3 真實機制 |
| 17 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/17_design_consistency_check.md` | 對齊論文/README/知識庫 |
| 18 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/18_purity_calculator_failure_root_cause.md` | collectPloidyRatio 修復 |
| 19 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/19_per_site_nuance_audit.md` | 13 位點 nuance 分類 |
| 20 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/20_whole_genome_paired_audit.md` | WG-level paired audit |
| **21** | **本檔** | 改動分類與爭議性矩陣 |

### 程式碼位置（按改動編號）

| # | 檔案 | 行號 |
|:-:|------|:-:|
| 1 | `HaplotagProcess.cpp` | 512-563 (getVote) |
| 2 | `HaplotagProcess.cpp` | 556-559 (HP33 encoding) |
| 3 | `HaplotagProcess.cpp` | 497-510 (INDEL guard) |
| 4 | `HaplotagProcess.cpp` | 484-495 (SNP guard) |
| 5 | `PhasingProcess.cpp` | 154-157 + `PhasingGraph.cpp` 1139-1145 |
| 6 | `PhasingProcess.cpp` | 162-172 + `PhasingGraph.cpp` 1147-1180 |
| 7 | `PhasingGraph.cpp` | 1147-1175 (本 session 新增) |
