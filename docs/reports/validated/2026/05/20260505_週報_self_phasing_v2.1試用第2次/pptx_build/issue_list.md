---
title: PPTX 18-slide Wave 1+2 整合 Issue 清單
type: pptx_review_summary
date: 2026-05-06
status: pending_user_decision
agent_count: 6
---

# Wave 1 + Wave 2 整合 Issue 清單

## 6 Agent 結果總表

| Agent | PASS | PARTIAL | FAIL | Top critical |
|-------|:-:|:-:|:-:|---|
| T 字體 | 74 | 10 | 6 | S15/S17 內文 ~11pt FAIL；嵌圖內字級偏小；28pt vs 32pt 規格衝突 |
| C 色彩 | 57 | 23 | 10 | S06 紅綠無雙編碼；S15 紅色誤用 (P1 priority 非 caveat)；S10 紅藍緊鄰 |
| L 佈局 | 64 | 20 | 6 | S01/S06/S10/S16 主視覺 < 55%；S06/S10 caveat 對齊偏 0.18in；S18 focal chip 偏左 |
| B 雙語 | 68 | 16 | 6 | □ 方塊 5 slides；S02/S04/S08/S10 中文 > 60 字 |
| **S 講稿** | — | — | **FAIL** | **Total 26.4 min vs 25 min target，+1.4 min** |
| **D 數據** | — | — | **FAIL × 2** | **★ S07/S08 noPath3 數字 scope 錯誤；★ S04 HP:i:33 方向反轉** |

---

## 🔴 P0 — Critical（事實錯誤，已於 2026-05-07 修正）

### #1 — S04 HP:i:33 metric 混淆 + 數字錯誤（已修）

**原狀**（`build_pptx.py:508` 修前）：
```
HP:i:33 從 14524 → 0（−100%，但忽略 Pass 2）
```
**真因**：把 3paths_audit 的 `HP_33` (threshold_compare scope) 跟 PI 報告 `HP:i:33` (whole genome) 混為一談；數字、scope、方向 3 個都錯。

**已修為**（`build_pptx.py:508`）：
```
HP:i:33 從 239,679 → 110,197（−54%，whole genome read 計數，PI 報告口徑）
```

---

### #2 — S07 OLD baseline HP_33 = 0 是錯的（修前 Agent-D 也誤判，已修）

**原狀**（`build_pptx.py:600` 修前）：
```
("OLD baseline", "路 1", "1.328", "0", "❌"),
```
**真實值**（3paths_audit §3.1 line 53）：OLD baseline HP_33 = **2,640**

**已修為**：`("OLD baseline", "路 1", "1.328", "2,640", "[X]")`

---

### #3 — S07/S08 noPath3 HP_33 = 14,524 scope 錯誤（已修）

**原狀**：S07 第 5 行 + S08 Result 卡都寫 `noPath3 HP_33 = 14,524`
**真因**：14,524 是 OLD_v5_flag (3paths_audit threshold_compare 全 BAM) 的數字；noPath3 真實值是 force_path2only §3.4 的 **15-site cherry-picked HP_33 = 28**

**已修為**：
- S07: `("NEW V5 noPath3", "路 2 only (forced)", "1.127*", "28*", "[OK]*")` + footer「* noPath3 數字來自 §3.4 15-site cherry-picked sample」
- S08 Result 卡: `Result (15-site cherry-picked): HP1:HP2 = 1.127 / HP_33 = 28 / 反轉方向與舊 V5 全 BAM (0.735) 一致 / 全 BAM 等價量化待 monitor`

---

## 🟡 P1 — Important（顯著影響但非事實錯誤）

### #3 — Speaker total 時長超 1.4 min（Agent-S）

**現況**：18 張總字數 10,547 字 ≈ 26.4 min vs target 25 min，**FAIL**。Tier 2 字數 0/18 在 300-360 範圍（16 張超 / 2 張不足）。

**修正建議**（拆遷至 [ORAL-OPTIONAL] 或刪除）：
- S04: ploidy 公式段落 60 字 → Tier 3
- S04: Evidence B 重述 35 字 → Tier 3
- S02: meta-commentary 80 字（敘事框架已 Slide 01 建立）→ 刪除
- S07: 「1.127 vs 0.735 差距解釋」65 字 → Tier 3
- S14: 「Memory 兩條」53 字 → 確認跳過
- S14: Thread B 重述 55 字 → 引用 S02 即可

合計可省 ~350-540 字（≈ 50-80 sec），可使實際時長回到 25 min。

---

### #4 — □ 方塊渲染失敗 5 slides（Agent-B）

**現況**：S07 / S10 / S13 / S15 / S17 wireframe 中出現 `□`、`□□□` 方塊（Unicode 符號 ⭐ / ✓ / ⚠ / ① 等在 PIL CJK fallback 下渲染失敗）。

**注意**：實際 PPTX 在 PowerPoint 開啟時可能不會出現方塊（PPT 自帶字型可能 fallback 到 Segoe UI Symbol）。但 PIL wireframe 確實顯示為方塊 → 視覺驗證受阻 + Memory `feedback_matplotlib_cjk_font_rule` 警告過。

**修正建議**：將下列符號全部改為 ASCII 等價：
- ⭐⭐⭐ → `[P0]` or `[High]` or 純文字
- ✅ ⚠ ❌ → `[OK]` `[!]` `[X]`（已部分採用，但未全替換）
- ① ② ③ → `1.` `2.` `3.`
- 具體位置（grep `⭐\|✅\|⚠\|❌\|①\|②\|③` build_pptx.py）

---

### #5 — Slide 內字級 < 14pt 下限（Agent-T）

| Slide | 位置 | 估算字級 | 建議 |
|:-:|---|:-:|---|
| S15 | P1-P5 priority 卡片內輔助說明（「決定性—若不做主結論失基礎」等）| ~11pt | 改 14pt（`add_text_with_fallback font_size=14`）或標 footnote 9pt 例外 |
| S17 | Q&A 雙欄箭頭說明行 | ~11-12pt | 改 14pt（PI 25 min 內需快速辨讀） |
| S03 | 嵌入 timeline figure 內 commit label / date | ~9-10pt | matplotlib `make_fig_s03_*.py` fontsize 從 9 → 12 |
| S06 | 嵌入 codeflow figure 內 step text | ~11-12pt | matplotlib fontsize → 14 |

---

## 🟢 P2 — Improvement（非阻擋但建議改進）

### #6 — 主視覺比例 < 55% 下限（Agent-L）

| Slide | 主視覺佔比 | 建議 |
|:-:|:-:|---|
| S01 | Main Thesis 框中心 62% canvas（偏低）| 上移 0.15in，元數據文字壓縮為 1 行 |
| S06 | 三欄 codeflow ~42% | 三欄高度從 2.9in 拉伸至 3.4in（+0.5in）|
| S10 | 兩框合計 ~45%（留白 35%）| 兩框各增 0.3in 高，或補 insight chip |
| S16 | Take-home 3 框 ~46% | 三框各 +0.15in，gutter 從 0.35 縮至 0.2 |

### #7 — Caveat strip 對齊偏移（Agent-L）

| Slide | 偏移 | 建議 |
|:-:|:-:|---|
| S06 | +0.20 in | caveat top 座標減 0.10 in |
| S10 | +0.18 in | gap 縮 0.08 in |
| S09 | +0.13 in | 邊界值，可選 |
| S12 | +0.12 in | 邊界值，可選 |
| S13 | +0.11 in | 邊界值，可選 |
| S18 | focal chip 偏左 0.30 in | left 座標 +0.30 in |

### #8 — 中文字數超 60 字（Agent-B）

| Slide | 字數 | 超出 | 建議 |
|:-:|:-:|:-:|---|
| S02 | ~130 | +70 | 本週敘事框架 3 行精簡或拆遷 speaker note |
| S04 | ~140 | +80 | 數字快照 4 行降為 2 行 + 移附錄 |
| S08 | ~120 | +60 | Result+Verdict 6 bullets 精簡至 3 |
| S10 | ~120 | +60 | Caveat Banner 完整引述句改 2 行 assertion |

### #9 — Color 語意/雙編碼（Agent-C）

| Slide | 問題 | 建議 |
|:-:|---|---|
| S06 | 紅綠並列無雙編碼 colorblind FAIL | 路 3 改 amber 或加斜線；節點加數字 1/2/3 |
| S15 | P1 紅色誤用（priority 非 caveat）+ 4 accent | P1 改 blue（method）；priority 用數字徽章 |
| S10 | 4 accent + 紅藍緊鄰 | artifacts 框改 grey 邊框 |
| S02/S09/S12/S16 | 語意偏差 PARTIAL | amber/green/red 對應 token 角色嚴格化 |

### #10 — Evidence path 缺失（Agent-S）

S12 / S15 / S17 缺 .md source path 引用，建議補：
- S12 footer: `InterSubMod/docs/experiments/in_progress/2026/04/20260501_latest_longphase_to_3paths_audit_01.md`
- S15: 各 Priority 對應 source 簡稱
- S17: 必問 Q1-Q3 各補 1 個 .md 引用

---

## 📋 修正範圍 3 選項（待用戶決策）

| 範圍 | 包含 | 估時 | 風險 |
|------|------|------|------|
| **R1（最小，必修）** | P0 #1+#2 兩個事實錯誤 | 10 min | 不修則 PI 報告含錯誤數字傳播 |
| **R2（建議）** | R1 + P1 #3+#4+#5（時長 / 方塊 / 字級）| 30-45 min | 講者時長超 + wireframe 視覺驗證受阻 |
| **R3（完整）** | R1 + R2 + 全 P2 | 1.5-2 hr | 視覺比例/對齊/字數/色彩語意全套 polish |

---

## 個人風格累積觸發點（待你反饋）

依 `InterSubMod/.claude/skills/pptx-build/prompts/feedback_classification.md`，下列 issue 若你決定修正，需強制分類「通用 vs 本次特定」：

1. **emoji/Unicode 符號禁用** → 通用候選（v3 case studies 已多次出現相同問題；強烈建議升 active 規則）
2. **Priority 數字標記用顏色 vs 徽章** → 本次特定（V5 audit 主題特殊）or 通用？
3. **嵌圖 fontsize 下限 14pt（matplotlib 預設過小）** → 通用候選
4. **HP:i:33 方向陳述標準** → 本次特定（self-phasing 主題）

每個 issue 修正前我會強制 AskUserQuestion 分類，累積到 `InterSubMod/.claude/skills/pptx-build/style_library/personal_style_log.md`。
