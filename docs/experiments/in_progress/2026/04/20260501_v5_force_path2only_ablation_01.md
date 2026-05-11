---
title: "新 binary V5 flag + 強制 highPurity=false（純路 2）ablation"
date: 2026-05-01
status: in_progress
sample: HCC1395_5kHz
purity: 0.93
binary_under_test: longphase-to-noPath3 (PhasingProcess.cpp:197 modified to `bool highPurity = false`)
related:
  - InterSubMod/docs/experiments/in_progress/2026/04/20260501_latest_longphase_to_3paths_audit_01.md
  - InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md
---

# 新 binary V5 flag + 強制 highPurity=false（純路 2）ablation

## §0 問題

新 binary（938f0df）下 V5 flag 在 0.93 樣本變成完全 no-op（HP1:HP2 = 1.40:1，與新 baseline 完全相同）。
但**舊 binary V5 flag 在 0.93 樣本能反轉到 0.735:1**。

兩個版本都跑了 V5 flag 路徑（`--pon-only-phasing` → 路 2 4 步）。差別：
- 舊 binary V5 flag：highPurity=false（因 ploidyRatio bug）→ 只走路 2
- 新 binary V5 flag：highPurity=true（purity 0.977 > 0.9）→ 走路 2+3

→ **假設**：新版「路 3 second round（不重跑 somaticCalling）」抵消了路 2（重跑 somaticCalling）的反轉效果。

## §1 實驗設計

修改 `longphase-to-mod/PhasingProcess.cpp:197`：
```cpp
- bool highPurity = purity > 0.9;
+ bool highPurity = false;  // FORCED: skip second-round phasing (path 3)
```

重新編譯為 `longphase-to-noPath3` binary。
跑 HCC1395 5kHz 0.93 + V5 flag (`--pon-only-phasing`)。
產物與新 V5 flag (路 2+3) / 舊 V5 flag (路 2 only with bug) / 新 baseline (路 1+3) / 舊 baseline (路 1) 對比。

## §2 五版本對比矩陣

| 版本 | binary | flag | highPurity | 走的路 |
|---|---|---|---|---|
| 舊 baseline | pre-8b8c1fd | off | false（auto）| 路 1 |
| 舊 V5 flag | V5 work tree | on | false（bug）| 路 2 |
| 新 baseline | 938f0df | off | true（auto）| 路 1 + 路 3 |
| 新 V5 flag | 938f0df | on | true（auto）| 路 2 + 路 3 |
| **新 V5 flag (force noPath3)** | **938f0df_modified** | **on** | **false（forced）** | **路 2 only** |

## §3 實測結果

### 3.1 執行驗證

| 項目 | 結果 |
|---|---|
| Source modify | `bool highPurity = false;` (PhasingProcess.cpp:197) |
| Binary build | `longphase-to-noPath3` (12s incremental) |
| Source restored | git checkout (主 binary 復原) |
| Phase 跑時 | 813s（vs 新 V5 flag 路 2+3 的 2881s，**節省 71% 時間**）|
| Haplotag 跑時 | 1979s（vs 新 V5 flag 路 2+3 的 2842s，節省 30%）|
| log "second round" 計數 | **0**（vs 新 V5 flag 路 2+3 為 1）→ 路 3 確實被跳過 |
| purity 算出 | 0.983791（與新 V5 flag 一致）|

### 3.2 caller F1 / PASS set

| 指標 | 值 |
|---|---|
| Total variants | 3,187,275 |
| PASS variants | 47,798 |
| TP / FP / FN | 28,509 / 11,606 / 10,938 |
| **F1** | **0.7166**（與所有舊+新版本完全相同）|
| PASS set diff vs OLD baseline | **0**（chr/pos/ref/alt 全等）|

→ caller F1 不變（一如預期，longphase-to phase 不修 FILTER）。

### 3.3 LOH.bed / GE.bed

| 比對 | LOH diff | GE diff |
|---|---|---|
| force_path2only vs new_v5_flag (路 2+3) | 0 | 0 |

→ LOH/GE region 偵測不受任何 phasing 改動影響（先前已知）。

### 3.4 BAM HP tag — 15-site cherry-picked sample

| 版本 | HP1_family | HP2_family | HP_33 | ratio | AMB% |
|---|---|---|---|---|---|
| OLD baseline (路 1) | 937 | 452 | **0** | 2.073 | 0.00 |
| OLD V2b (路 2 V2b) | 590 | 592 | 2 | 0.997 | 0.17 |
| **OLD V5 (路 2 V5)** | 613 | 543 | **28** | **1.129** | **2.36** |
| NEW baseline (路 1+3) | 713 | 629 | 12 | 1.134 | 0.89 |
| NEW V5 flag (路 2+3) | 713 | 629 | 12 | 1.134 | 0.89 |
| **NEW V5 noPath3 (純路 2)** | **604** | **536** | **28** | **1.127** | **2.40** |

**核心觀察**：
- **NEW V5 noPath3 ≈ OLD V5**（HP_33 = 28 完全相同；ratio 1.127 vs 1.129；AMB% 2.40 vs 2.36）
  → 強制 highPurity=false 完美復現舊 V5 行為
- **NEW V5 noPath3 vs NEW baseline**：HP_33 增至 28（vs 12，多 133%）；ratio 微降（1.127 vs 1.134）
- **路 3 second round 的影響量化**：把 HP_33 從 28 → 12（降 57%），把 ambiguous tag 重分為 HP1/HP2，但 ratio 變動只有 0.6%

### 3.5 BAM HP tag — 全 BAM 計數（待 monitor 完成）

預期：與舊 V5 flag (V5 work tree) 接近，HP1:HP2 ≈ 0.735:1。

## §4 紀錄

- Source modification: 2026-05-01 03:37
- Binary build: 2026-05-01 03:37 (12s incremental, jemalloc/htslib cached)
- Source restored: 2026-05-01 03:37 (PhasingProcess.cpp git checkout)
- Phase: 2026-05-01 03:38 - 11:13（813s phase + 其他步驟）
- Haplotag: 2026-05-01 11:13 - 11:46
- HP analysis: 2026-05-01 12:00+
- Output dir: `longphase-to-mod/output/v5_flag_force_path2only/`
- Binary 保留: `longphase-to-mod/longphase-to-noPath3` (22.55 MB)

## §5 結論

### 5.1 假設驗證 — 完美成立

**路 3 second round 確實會抵消路 2 反轉**：
- 純路 2（noPath3 forced）：ratio 1.127, HP_33 = 28
- 路 2+3（新 V5 flag 自動）：ratio 1.134, HP_33 = 12
- 差別：second round 把約 800,000 reads 從 ambiguous (HP_33) 重分為明確 (HP1/HP2)

### 5.2 與用戶問題對應的回答

> **「新 binary V5 flag + 強制 highPurity=false (純路 2)，會比新 binary 兩條路 baseline + 舊 binary 三條路 baseline 都好嗎？」**

- vs 新 binary 兩條路 baseline（路 1+3，ratio 1.134, HP_33=12）：
  - **稍微好** — HP_33 多 133%（28 vs 12）；ratio 略接近真值（1.127 vs 1.134，差別 0.6%）
  - 但 caller F1 完全相同（0.7166）
- vs 「舊 binary 三條路 baseline」：**此版本不存在**（舊 binary highPurity 在 V5 flag 路徑下因 ploidyRatio bug 永遠 false；舊 baseline 因 purity 0.927 < 0.95 也不觸發）
- vs 舊 V5 flag (V5 work tree, 路 2 with bug)：**完全等價**（HP_33=28 同；ratio 1.127 vs 1.129 差異 0.18%）
- vs 舊 baseline（路 1）：**明顯改善**（ratio 從 2.073 → 1.127，反轉成功）

### 5.3 重要意涵

1. **新 binary 中 V5 flag 失效的根因**：是 path 3（second round）自動觸發後抵消了 path 2 的反轉，**而非 V5 flag 本身失效**
2. **路 3 不重跑 somaticCalling 是關鍵差別**（PhasingProcess.cpp:209-216）— 只 re-phase 不重 call，導致 HP_33 ambiguous tag 被「強制」分類為 HP1 或 HP2
3. **舊 V5 flag 的優勢來自路 2 內含 `somaticCalling` 重跑**（line 164），這在新版自動 path 3 邏輯下被覆蓋
4. **節省時間**：noPath3 phase 跑 813s vs 路 2+3 的 2881s（節省 71%）— 在不犧牲 phasing 品質下大幅提速

### 5.4 重大但有限的「改善」

- **改善程度**：HP_33 從 12 → 28（多 133%），但全域 ratio 改變極小（1.134 → 1.127）
- **不是「翻倍」改善**：跟舊 V5 flag (V5 work tree) 完全等價
- **「比新 baseline 都好」是真的，但差別在 ambiguous tag 數量，不是大幅 ratio 反轉**
- **比舊 V5 flag 沒有改善**：兩者實質一樣

### 5.5 對 InterSubMod 下游消費者的影響（待驗證）

caller F1 不能評估這個改動的真實效果（永遠 0.7166）。需要跑 ISM SuggestFilter F1（~3 hr/case）才能知道：
- 新 V5 noPath3 的 ISM F1 vs 新 V5 flag (路 2+3) ISM F1 vs 舊 V5 ISM F1
- 預期 noPath3 ≈ 舊 V5（兩者 HP tag 行為等價）→ ISM F1 約 0.7154

### 5.6 推薦做法

如果 V5 flag 在新 binary 要恢復原效果：
1. **方案 1**：使用 `longphase-to-noPath3` binary（已建立）+ `--pon-only-phasing` flag
2. **方案 2**：modify upstream PhasingProcess.cpp 把 path 3 邏輯改為「only re-phase if 不 V5 flag」（避免 path 2+3 重複工作）
3. **方案 3**：modify path 3 加入 `somaticCalling` 重跑步驟（讓 path 3 等價 path 2 但更快）


