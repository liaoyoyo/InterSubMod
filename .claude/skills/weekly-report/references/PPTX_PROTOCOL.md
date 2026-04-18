# PPTX 生成、驗證與設計規範

## Phase 3：撰寫 PPTX 三件套

輸出目錄：
```
docs/presentations/validated/YYYY/MM/YYYYMMDD_研究週報_{主題}/
├── 01_full_narrative_report.md
├── 02_ppt_slide_outline.md
├── 03_slide_layout_and_script.md
├── images/
└── pptx_config.json
```

### Phase 3.5：導演 Storyboard 審查（Phase 4 前必做）

產出 `00_director_storyboard.md`：
1. **Three-Act 故事弧**：定錨 → 信任 → 核心 → 根因+轉向 → 行動
2. **每頁 Storyboard**：核心訊息 + 觀眾心理 + 佈局草稿 + 偏離檢查
3. **順序邏輯**：先建立信任 → 再深入分析 → 最後揭示根因
4. 用戶確認後才進入 Phase 4

### Phase 4：生成 PPTX

1. 複製 config JSON 範本
2. 執行生成腳本
3. **截圖驗證（mandatory）**：渲染每張 slide 為 PNG，檢查 aspect ratio、文字溢出、顏色
4. 迭代至 vN

---

## Phase 5.5：12 Agent × 4 波次驗證

```
Wave 1 (內容)    Wave 2 (易讀)    Wave 3 (學術)    Wave 4 (PPTX)
 A1:邏輯         A4:敘事          A7:教授          A9:導演
 A2:質疑         A5:認知          A8:翻譯          A10:美術
 A3:佐證         A6:學生                           A11:審稿
                                                   A12:蒐集
```

每波次內平行執行，波次間依序修正。使用 `feature-dev:code-reviewer` subagent_type。

---

## PPTX 設計常數

> 完整規範：`docs/references/manual/20260407_PPT設計規範整合快速參考_01.md`

### 色彩

| 項目 | 值 | 用途 |
|------|-----|------|
| dark | #1E2A44 | 深色頁背景 |
| bg | #F7F3EC | 淺色頁背景 |
| accent | **#A85540** | 強調色（WCAG AA 合規） |
| accent2 | #6E9D8A | 大字/裝飾 |
| 文字 | #17202A | 主文字 |
| 弱化 | #5E6572 | 次要文字 |
| 深色頁文字 | #F5F5F5 | 深色頁主文字 |
| 綠 | #009E73 | 正面/TP/CONFIRMED |
| 紅 | #D55E00 | 負面/FP/REJECTED |
| 藍 | #0072B2 | 主數據 |

### 字體

| 項目 | 值 |
|------|-----|
| 標題 | Arial Bold, ≥32pt（深色頁 ×1.18 → 38pt） |
| 內文 | Arial, ≥14pt（密集頁 12pt） |
| EN 翻譯 | Arial, 70% 主字號, #5E6572, 縮排 0.25" |
| 絕對下限 | 9pt（Consolas 路徑） |

### 版面規則

| 項目 | 規則 |
|------|------|
| 視覺優先 | 視覺 55-65%, 留白 20-30%, 文字 ≤15% |
| 每頁 bullet | ≤4（行動頁 ≤3） |
| 每頁視覺元素 | ≤6 |
| 6x6 Rule | ≤6 行 × ≤6 詞/行 |
| 數據/概念頁 | **零 bullet** — 結論標題 + 圖表 |
| 雙語 | 中文主 + 英文副（60% 字號 + 縮排） |
| 封面頁 | 左文字 + 右留白供手動插圖 |
| valign | 所有文字 MSO_ANCHOR.MIDDLE |

### 5 種圖表模式

Pipeline Flow / Comparison Panel / Evidence Stack / Decision Tree / Concept Card

### 形狀語義

圓角矩形=步驟, 菱形=決策, 橢圓=起止, 六角=外部

---

## Phase 4.5：PPTX 多 Agent 驗證（6 Agent × 2 波次）

**Wave 1**（結構+視覺，4 parallel）：
- Agent-T：字體與字號
- Agent-C：色彩與對比度
- Agent-L：佈局與對齊
- Agent-B：雙語品質

**Wave 2**（內容+整合，2 parallel）：
- Agent-S：Speaker Notes
- Agent-D：數據準確性

---

## Phase 6：最終檢核（雙層清單）

### Must-Pass（不通過不可發送）
L1-L12：邏輯嚴謹度（結論附 evidence、因果鏈無跳步、結論/狀態表無矛盾、confound 控制、跨樣本標明、術語一致、圖表解釋性標題、Tier 分配一致）

### Should-Pass（品質）
Q1-Q13：Layer 0 存在、Metadata 完整、絕對路徑、Paired/TO 分開、PPTX 驗證、索引更新、Mermaid 可 render、4 波次驗證完成

---

## 圖表規範

- 初步觀察：分 Paired / TO 兩組
- 細緻觀察：2×2 per sample
- 固定順序：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829
- PPTX 圖片：fit-within + centered
- 寬圖/直式圖：自適應佈局
