---
title: "最新 longphase-to GitHub 程式碼（938f0df）三條路審計"
date: 2026-05-01
status: in_progress
sample: HCC1395_5kHz
purity: 0.93
binary_under_test: longphase-to (commit 938f0df, built 2026-04-30 15:34)
related:
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
  - InterSubMod/docs/reports/validated/2026/04/20260429_supplement_getVote_design_intent_QA_01.md
  - InterSubMod/docs/experiments/in_progress/2026/04/20260430_V3F_ablation_purity06_results_01.md
---

# 最新 longphase-to GitHub 程式碼三條路審計

## §0 一句結論

最新 GitHub commit **938f0df** 把 `highPurity` 觸發閾值從 `> 0.95` 改為 `> 0.9` + 改寫 purity 公式為三次多項式，讓新版 baseline 自動觸發 V5 flag 原本要做的 second round phasing；**V5 flag 在新 binary 變成完全 no-op**（PASS variant set / HP tag / LOH.bed / F1 全部與新 baseline 等價）。

## §1 三條路定義（程式碼出處）

`longphase-to-mod/PhasingProcess.cpp:142-220`：

```
Pass 1（per-chr 並行）：
├── 路 1（無 flag）：somaticCalling → phasingProcess
└── 路 2（--pon-only-phasing）：4 步 PON-only flow
    1. convertNonGermlineToSomatic + phasingProcess(..., nullptr)
    2. resetNonPonOrigin + somaticCalling
    3. syncPhasingResultOrigins
    4. collectPloidyRatio  ← d0bcd8c 修補

Pass 2（per-chr 並行）：
└── 路 3（highPurity = purity > 0.9）：
    convertNonGermlineToSomatic + reset posPhasingResult + phasingProcess(..., nullptr)
```

## §2 四版本對應路徑

| 版本 | binary commit | flag | purity | highPurity 觸發 | 走的路 |
|---|---|---|---|---|---|
| OLD_baseline | pre-8b8c1fd（baseline）| off | 0.927 | NO（< 0.95）| 路 1 |
| OLD_v5_flag | V5 work tree | on | 0（bug）| NO（bug）| 路 2 |
| **NEW_baseline_two_pass** | **938f0df** | off | **0.977** | **YES**（> 0.9）| **路 1 + 路 3** |
| **NEW_v5_flag_three_pass** | **938f0df** | on | **0.984** | **YES**（> 0.9）| **路 2 + 路 3** |

## §3 核心驗證結果

### 3.1 HP tag ratio（per_site_hp_counts）

| 版本 | HP1 | HP2 | HP33 | ratio | vs paired truth (0.92:1) |
|---|---|---|---|---|---|
| OLD_baseline | 96,731 | 72,883 | 2,640 | 1.328 | 偏 HP1 |
| OLD_v5_flag | 65,022 | 88,413 | 14,524 | 0.735 | **反轉成功** |
| NEW_baseline_two_pass | **100,869** | **72,075** | **3,468** | **1.400** | 偏 HP1（更甚）|
| NEW_v5_flag_three_pass | **100,869** | **72,075** | **3,468** | **1.400** | **與新 baseline 完全相同** |

來源：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/threshold_compare.tsv`

### 3.2 PASS variant set（caller 過濾）

| 版本 | total variants | PASS | PASS set 與 OLD_baseline 差異 |
|---|---|---|---|
| OLD_baseline | 3,187,275 | 47,798 | 0 |
| OLD_v5_flag | 3,187,275 | 47,798 | 0 |
| NEW_baseline_two_pass | 3,187,275 | 47,798 | 0 |
| NEW_v5_flag_three_pass | 3,187,275 | 47,798 | 0 |

→ **PASS variant set 在所有 4 版本完全相同**（chr/pos/ref/alt 全等）。longphase-to phase 不修 FILTER 欄位，FILTER 由 caller (ClairS-TO) 決定，phasing 不會 swap PASS↔FAIL。

### 3.3 SEQC2 Caller F1（HCC1395 5kHz, 0.93 purity）

| 版本 | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| OLD_baseline_093 | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| OLD_v5_flag_093 | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| NEW_baseline_two_pass | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |
| NEW_v5_flag_three_pass | 28,509 | 11,606 | 10,938 | 0.7107 | 0.7227 | **0.7166** |

→ **F1 在所有 4 版本完全相同 = 0.7166**（PASS set 全等的必然結果）。

來源：`InterSubMod/docs/experiments/in_progress/2026/04/figures/20260501_latest_binary_f1/latest_binary_f1.tsv`

### 3.4 LOH.bed / GE.bed（除 RGB 顏色欄）

| 比對 | LOH.bed diff | GE.bed diff |
|---|---|---|
| OLD_baseline vs OLD_v5_flag | 0 | 0 |
| OLD_baseline vs NEW_baseline | 0 | 0 |
| NEW_baseline vs NEW_v5_flag | 0 | 0 |

→ **LOH/GE 偵測跨所有 4 版本完全相同**（1094 個 LOH region 全等）。LOH 偵測在 `Clip::detectLOHRegion` 階段完成，發生在 phasing graph 之前，不受 self-phasing / V5 / new pass 影響。

### 3.5 BAM 內容差異（byte-level）

| 比對 | 差異 |
|---|---|
| NEW_baseline.bam vs NEW_v5_flag.bam | 2 bytes（287GB 中的 2 bytes，可能 metadata）|
| NEW_baseline.vcf vs NEW_v5_flag.vcf | 15 bytes（655MB 中的 15 bytes，header 微差）|

→ Byte-level 微差不影響任何 phasing/HP/F1 結果。

## §4 回答用戶四個驗證項

### (i) 新版 baseline F1 vs 舊版 baseline 0.7166 — 是否退步？

**沒有退步，完全相同 = 0.7166**。

原因：longphase-to phase 不改 FILTER 欄位。FILTER 由 caller (ClairS-TO) 決定，phasing 後保留原值。所有 4 版本的 PASS set 完全相同（47,798 個 variants），因此 F1 必然相同。

意涵：**caller-level F1 不能用來評估 phasing/V5 改動的優劣** — 必須用 ISM SuggestFilter 下游消費者的 F1 才能看出差別。

### (ii) 新版 baseline LOH.bed 與 GE.bed 是否仍然帶 self-phasing 偏移？

**LOH.bed / GE.bed 本身沒被 self-phasing 影響**。

原因：LOH region 偵測在 `Clip::detectLOHRegion(snpFile, chrInfo.LOHSegments)` 階段完成（`PhasingProcess.cpp:147`），發生在 phasing graph `addEdge` 之前。LOH 邊界僅依 caller VCF 的 AF distribution clipping 決定，與 phasing 結果無關。

但**self-phasing 偏移仍存在於 HP tag 分配**：
- 新 baseline / 新 v5_flag HP1:HP2 = **1.40:1**（vs paired truth 0.92:1）
- 新版 second round 反而比舊版 baseline (1.33) 更偏 HP1
- 表示 **second round 沒有真正解決 self-phasing**（甚至略微惡化）

### (iii) d0bcd8c 的 ploidyRatio 修補是否反而讓 second round 更偏？

**不直接相關**。`d0bcd8c` 修補的是 **V5 flag 路徑下** Pass 2 的 `ploidyRatioMap` 不會被收集，導致 q1=q3=0 → purity polynomial 算出 0 → highPurity 永遠 false → V5 flag 從不觸發 second round。

修補後：
- V5 flag 路徑能正確計算 purity（如新 v5_flag purity = 0.984）
- highPurity 才能正確觸發
- 但 second round 邏輯本身（路 3）沒被 d0bcd8c 改動

新 baseline 偏移 1.40 來自：
- `938f0df` 的新 purity 公式（HCC1395 0.93 真實 purity → 算出 0.977，估計偏高）
- highPurity threshold 從 0.95 降到 0.9
- → 自動觸發 second round（路 3）
- → second round 用 `nullptr` ploidyRatioMap 進 phasingProcess（line 217）
- → 改變 phasing graph 行為，但仍不能解決 self-phasing

### (iv) 新 binary v5_flag（三條路）vs 新 binary baseline（兩條路）— 哪個好？

**完全等價，V5 flag 在新 binary 是純 no-op**。

完整等價證據（所有指標逐項相同）：

| 指標 | 兩條路 baseline | 三條路 v5_flag | 差異 |
|---|---|---|---|
| HP1 / HP2 / HP33 | 100869 / 72075 / 3468 | 100869 / 72075 / 3468 | 0 |
| HP1:HP2 ratio | 1.400 | 1.400 | 0 |
| Total variants | 3,187,275 | 3,187,275 | 0 |
| PASS variants | 47,798 | 47,798 | 0 |
| PASS variant set (chr/pos/ref/alt) | — | — | 0 |
| LOH.bed lines | 1,094 | 1,094 | 0 |
| LOH.bed content (excl RGB) | — | — | 0 |
| GE.bed content (excl RGB) | — | — | 0 |
| Caller F1 (SEQC2 HC) | 0.7166 | 0.7166 | 0 |
| BAM size | 287,287,277,312 bytes | 287,287,277,310 bytes | 2 bytes |
| VCF size | 655,254,528 bytes | 655,254,543 bytes | 15 bytes |

→ V5 flag 在新 binary 不能產生比 baseline 更好的結果。

## §5 與舊 V5 flag 的真實對比

舊 V5 flag（V5 working tree）vs 新 baseline（938f0df 兩條路）的核心差別：

| 項目 | OLD_v5_flag | NEW_baseline |
|---|---|---|
| HP1:HP2 ratio | **0.735**（成功反轉）| 1.400（仍偏 HP1）|
| HP_33 count | **14,524**（正確標 ambiguous）| 3,468（少 4×）|
| 觸發機制 | flag 強制 | purity 自動 |
| second round（路 3）| 走（in Pass 1 路 2）| 走（in Pass 2 路 3）|

→ **舊 V5 flag 比新 baseline 更接近 paired ground truth (0.92:1)**。新版 baseline 自動 Pass 2 走的「路 3」**不等同**舊 V5 flag 走的「路 2」：

- 路 2（V5 flag 4 步）：`convertNonGermlineToSomatic` → phasingProcess → `resetNonPonOrigin` → **somaticCalling**（重新跑 caller 邏輯）→ syncOrigins → collectPloidyRatio
- 路 3（new highPurity）：`convertNonGermlineToSomatic` → reset PosPhasingResult → phasingProcess（**僅單一 phasing**，無 somaticCalling）

差別在 **somaticCalling 是否重跑** + **ploidyRatio 是否重新收集**。

## §6 重要意涵與後續

### 6.1 新版 baseline 的 self-phasing 沒被解決

新 baseline HP1:HP2 = 1.40:1 比舊 baseline 1.33:1 **更偏 HP1**。雖然新版自動跑 second round，但仍未達舊 V5 flag 的反轉效果（0.735）。

**推測原因**：
- 路 3 不重跑 somaticCalling → somatic 變異仍維持原本錯誤分類
- 只重 phase 不重 call → self-phasing 偏移殘留

### 6.2 V5 flag 與 V5 audit suite 結論需修訂

`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/` 過去結論基於 V5 work tree binary（V5 flag 確實有效），但在最新 GitHub binary（938f0df）下：

| 舊結論 | 新證據下需修訂 |
|---|---|
| 「V5 flag 改善 self-phasing」 | 在新 binary 下 V5 flag 完全 no-op |
| 「V5 不改 caller F1」 | 仍正確（caller F1 永遠 0.7166）|
| 「V5 修 HP_33 enum bug」 | 在新 binary 下 baseline 已自動產生 HP_33（3,468）|

### 6.3 ISM SuggestFilter F1 待驗證

caller F1 不能比較 V5 改動的真實效果，必須跑 ISM SuggestFilter 下游消費者：

| 版本 | ISM F1（待跑）|
|---|---|
| OLD_baseline | 0.7157（已知）|
| OLD_v5_flag | 0.7154（已知）|
| NEW_baseline_two_pass | ❌ 未跑 |
| NEW_v5_flag_three_pass | ❌ 未跑（預期 = 兩條路 baseline 因等價）|

### 6.4 跨 purity 驗證待補

| Purity | OLD_baseline | OLD_v5_flag | NEW_baseline | NEW_v5_flag |
|---|---|---|---|---|
| 0.93 | ✅ 完整 | ✅ 完整 | ✅ 完整 | ✅ 完整 |
| 0.6 | ✅ 完整 | ✅ 完整 | ❌ 未做 | ❌ 未做 |

低 purity (< 0.9) 不會觸發 highPurity → 預期新 baseline (0.6) 行為與舊 baseline 相同。

## §7 落點與索引

| 檔案 | 內容 |
|---|---|
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260501_latest_binary_f1/latest_binary_f1.tsv` | 4 版本 SEQC2 F1 |
| `InterSubMod/docs/experiments/in_progress/2026/04/figures/20260501_latest_binary_f1/path_evidence_summary.tsv` | 路徑與 HP/LOH/F1 一覽 |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/threshold_compare.tsv` | HP1/HP2/HP33 跨版本（既有）|

## §8 變更歷史

| 日期 | 動作 |
|---|---|
| 2026-05-01 03:00 | 建立本報告；驗證新 binary V5 flag 完全 no-op；4 版本 PASS set 全等；F1 全 0.7166 |
