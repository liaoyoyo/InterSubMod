<!--
建立時間: 2026-04-01 22:00
目標: Phase 4 — 三位審閱者交叉審閱結果彙整
處理範圍: 邏輯/版面/講稿三維度審閱
關聯檔案:
  - 01_full_narrative_report.md
  - 02_ppt_slide_outline.md
  - 03_slide_layout_and_script.md
-->

# 交叉審閱結果彙整

## 評分總覽

| 審閱者 | 角色 | 評分 | 焦點 |
|--------|------|:----:|------|
| 嚴格教授 | 學術邏輯 | 7.2/10 | 邏輯鏈、定義、因果推論 |
| PPT 專家 | 版面設計 | 7.0/10 | 認知負荷、視覺層次、圖片 |
| 講稿教練 | 口述流暢 | 7.5/10 | 自然度、時間分配、過渡句 |

---

## 🔴 Critical Issues（必須修正）

### Issue 1：數字不一致 — LOH penalty 觸發率 vs extreme LOH rate
- **問題**：Slide 8 用 44.5%/35.8%（LOH penalty 觸發率），Slide 10 用 44.6%/35.9%（extreme LOH rate），兩組數字差 0.1pp，描述同一現象卻未區分來源
- **影響**：聽眾會困惑為何兩組數字幾乎一樣但不同
- **建議**：確認是否同一指標 → 統一為一組；若不同 → 明確標註定義差異

### Issue 2：否決假說數量不一致 — 3 vs 4
- **問題**：01_report 寫「4 項」（O11+O12+O13+HP0 filter），但 PPT 封面和結尾寫「3 項」，中間 Slide 7/16 又寫「4 項」
- **影響**：自相矛盾，嚴重損害可信度
- **建議**：統一為 4（需在 PPT 中提及 HP0 filter）或統一為 3（刪除 HP0 filter 計數）

### Issue 3：四步邏輯鏈隱含假設（Slide 10）
- **問題**：Step 1→2 假設「somatic allele 被用作 phasing anchor」，但未交代 LongPhase-TO 實際是否以被評估的 variant 做 anchor
- **影響**：因果鏈最關鍵的一環缺乏確認
- **建議**：查閱 LongPhase-TO 演算法確認機制，或加上「假設」標籤

### Issue 4：缺少 Cell Line vs Clinical Limitation
- **問題**：所有結論基於 7 個高 purity 細胞系，未提及臨床低 purity 樣本的 generalizability
- **影響**：LOH enrichment 翻轉機制高度依賴 tumor purity，低 purity 時偏移幅度減弱
- **建議**：在 Slide 18 加 limitation caveat

---

## 🟡 Moderate Issues（建議修正）

### 定義相關

| # | 問題 | 建議 |
|---|------|------|
| M1 | HP_Ratio 邊界條件未說明（HP1+HP2=0 時？HP:i:33 如何處理？） | Section 2.3 補充精確定義 |
| M2 | LOH-like 閾值 0.1/0.9 缺乏 sensitivity analysis | 至少口述提及閾值穩健性 |
| M3 | L2/L3 術語未定義 | Slide 12 加一行定義 |
| M4 | 觀察數量「九項」vs 實際列出八項（缺 O6/O9） | 修正為「八項」或解釋缺失 |

### 版面相關

| # | 問題 | 建議 |
|---|------|------|
| M5 | Slide 10 嚴重過載（四步邏輯鏈+量化驗證+圖） | 拆為 10A（機制）+ 10B（驗證） |
| M6 | Slide 13 兩個不相關發現硬擠一頁 | 拆分或將 Simpson's 移至 appendix |
| M7 | Slide 3 缺正式 pipeline 流程圖（純 ASCII） | 製作圖示或引用現有 |
| M8 | Slide 12 三假說全為純文字，缺圖片 | 至少插入一張代表性圖 |
| M9 | 多處表格可圖形化（Slide 8 gauge/Slide 11 dumbbell/Slide 16 radar） | 降低表格疲勞 |

### 講稿相關

| # | 問題 | 建議 |
|---|------|------|
| M10 | 4 處缺過渡句（4→5, 6→7, 11→12, 14→15） | 各加一句 bridge sentence |
| M11 | 專有術語口語化不足（PERMANOVA、Shapley、L2 collider bias、ICC） | 口述用功能語言，術語留投影片 |
| M12 | 數字密度過高（每頁 4-6 個精確數字） | 口述最多 2-3 個核心數字 |
| M13 | 高潮 Slide 9/16 語氣力度不足 | 加 [停頓] 標記、更斷言式結論句 |

### 邏輯相關

| # | 問題 | 建議 |
|---|------|------|
| M14 | FN rescue 語氣偏樂觀（enrichment 僅 0.805） | 標註為「假說待驗」 |
| M15 | 「必須分離」隱含「分離後能有效」但尚未證明 | 加 caveat：TO 獨立模型可行性待驗 |
| M16 | O11 否決的 n_reads 校正方法未交代 | 補充具體校正方法 |

---

## 🟢 Minor Issues（可選修正）

| # | 問題 | 位置 |
|---|------|------|
| m1 | Slide 13 講稿「6 個中有 7 個」→ 應為「7 個中有 6 個」 | 03_script Slide 13 |
| m2 | Slide 16 Paired 欄「5/9 一致」含義不清 | 統一為「5/9 反轉」 |
| m3 | 封面摘要過長（6 箭頭）→ 精簡為 15 字標語 | Slide 1 |
| m4 | Enrichment 是 pooled 還是 mean 未標註 | Section 2.4 |
| m5 | Truth set 來源及 bias 未討論 | Slide 18 limitation |
| m6 | Slide 6 表格與柱狀圖訊息重疊 | 二選一 |
| m7 | Slide 15 缺 before/after 視覺化 | 加 gauge chart |

---

## 三位審閱者共識的 Top 5 修正優先順序

| 優先 | 修正項 | 來源 | 預期改善 |
|:----:|--------|------|---------|
| **1** | 統一數字（44.5/44.6、3/4 假說、九/八項觀察） | 教授 🔴 | 消除自相矛盾 |
| **2** | 拆分 Slide 10（機制+驗證）+ 補充 LongPhase-TO 機制 | 教授+PPT | 邏輯嚴謹+認知負荷降低 |
| **3** | 補過渡句 4 處 + 高潮 Slide 停頓標記 | 講稿+PPT | 敘事流暢度質變 |
| **4** | 加 Limitation（cell line/clinical、truth set bias） | 教授 🔴 | 學術完整性 |
| **5** | 專有術語口語化 + 降低數字密度 | 講稿 | 聽眾友善度提升 |

---

## 各審閱者肯定的優點

- **邏輯**：9 觀察 + 3 否決的系統性設計嚴謹；Shapley 分解提供獨立驗證；HP bug→修正→重跑的工作流程展示完整
- **版面**：Slide 2/5/8/9 的對比設計出色；熱圖佔 60% 的決策正確；敘事主線（問題→意外→分析→結論→行動）結構清晰
- **講稿**：總時間 ~21 分鐘合理；口語化程度不錯；故事弧完整有升華（FP Filter→FN Rescue）
