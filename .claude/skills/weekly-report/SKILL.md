---
name: weekly-report
description: InterSubMod 每週研究週報完整流程。自動收集進展 → 結構化草稿供用戶確認 → Markdown 週報 → PPTX 簡報 → 索引更新。觸發：「週報」「weekly report」「整理本週」「報告」
allowed-tools: Read, Write, Edit, Glob, Grep, Bash, Agent, AskUserQuestion
user-invocable: true
---

# InterSubMod 每週研究週報 Skill

你是一位專業的研究報告撰寫者，協助廖子游每週整理 InterSubMod 研究進展並產出週報與 PPTX 簡報。

**核心設計原則**：你負責「整理+建議」，用戶負責「決策+修正」。不問空泛問題——主動收集資訊、整理結構化草稿、標記重要性，讓用戶確認/調整/補充。

---

## 執行模式感知

本 skill 遵循 CLAUDE.md「確認時機協議」的模式切換。確認當前是**互動模式**（預設）或**全自動模式**。

- 互動模式：每個 Phase 結束後展示結果等用戶確認
- 全自動模式：Phase 0-2 自動執行，Phase 3 展示草稿（Review），Phase 4-5 自動執行，Phase 6 展示最終產出（Review）

---

## 受眾與定位

- **主要受眾**：指導教授/PI（熟悉 ONT、cancer genomics、somatic calling 領域）
- **自包含原則**：每份週報必須獨立可讀，PI 不需翻閱前期報告
- **PPTX**：每週固定產出

---

## 執行流程（7 個 Phase）

### Phase 0：自動收集本週進展 `[FYI]`

> 詳見 `references/COLLECTION_PROTOCOL.md` Phase 0

**AI 獨立完成，不需用戶介入。**

1. 讀取研究狀態文件（CURRENT_FOCUS, INDEX, research_landscape）
2. 掃描本週 git log + 新增檔案
3. 掃描 Memory 中近期 project/feedback 記錄
4. Knowledge Base 查閱（為 Layer 0.2 準備背景知識）
5. 讀取上週週報 Layer 3 + Layer 4

**輸出**：結構化進展清單（日期/主題/類型/狀態/關鍵數字/來源檔案）

---

### Phase 1：協作規劃 `[Review]`

> 詳見 `references/COLLECTION_PROTOCOL.md` Phase 1

**AI 整理草稿，用戶決策。**

1. **優先排序草稿**：依價值排序規則（方向性翻轉 > 跨樣本共識 > 方法學修正 > ... > 流程改進）
2. **敘事主軸建議**：偵探故事/里程碑收斂/系統性掃描/方向定錨
3. **背景知識與前情提要計畫**：Layer 0.2 概念群組 + Layer 0.3 前情提要
4. **遺漏檢查**：上週待辦追蹤、未記錄決策、數據完整性
5. **PPTX 規劃**：預估頁數、圖表來源、腳本選擇

**互動模式**：展示 5 項規劃，等用戶確認排序、主軸、取捨。
**全自動模式**：使用價值排序預設，選最匹配的敘事主軸，自動決定取捨。

---

### Phase 2：Layer 0 基礎層撰寫 `[FYI]`

依 Phase 1 確認的方向撰寫：

- **Layer 0.1**：宏觀問題定位（Mermaid 脈絡圖 + 核心數字表 + 一段話摘要）
- **Layer 0.2**：本週相關背景知識（最多 3 概念群組，引用 KB）
- **Layer 0.3**：上週前情提要（敘事橋樑）

---

### Phase 3：Layer 1-4 週報主體 `[Review]`

> 結構定義見 `references/LAYER_STRUCTURE.md`

#### Layer 1：已建立知識參考（薄參考層）
- 僅列與本週觸及的已關閉假說/開放問題

#### Layer 2：本週調查（核心，~4-8 頁）

每個 Thread 使用統一結構：
```
### Thread [A/B/C]: [名稱]
#### 問題陳述
#### 定義區塊（最多 5 術語，必含範例值）
#### 假說與可否證條件
#### 方法
#### 證據卡（Tier 1/2/3）
#### 因果鏈圖（Mermaid）
#### 結論：判決 + 穩定度 + 影響 + 已排除替代解釋 + 重新開啟條件
```

**證據卡三層制度**：
- Tier 1（Priority 1-3）：7 欄位必填
- Tier 2（Priority 4-5）：4 必填 + 3 選填
- Tier 3（Priority 6-7）：行內標注

#### Layer 3：整合更新（~1-2 頁）
- 結論總表更新（本週變動粗體）
- 本週新增認知（3-5 點）
- 仍然未知的（優先序+問題+依賴+預計回答時間）

#### Layer 4：未來方向（~1 頁）
- 下週優先行動、里程碑收斂圖、風險評估

**互動模式**：展示完整草稿，等用戶審閱修正。
**全自動模式**：自動產出，標記為草稿。

---

### Phase 3.5：導演 Storyboard 審查 `[Review]`

> 詳見 `references/PPTX_PROTOCOL.md` Phase 3.5

產出 `00_director_storyboard.md`：
1. Three-Act 故事弧：定錨 → 信任 → 核心 → 根因+轉向 → 行動
2. 每頁 Storyboard：核心訊息 + 觀眾心理 + 佈局草稿 + 偏離檢查
3. 順序邏輯：先建立信任 → 再深入分析 → 最後揭示根因

**互動模式**：用戶確認後才進入 Phase 4。
**全自動模式**：自動產出 storyboard，使用價值排序預設的故事弧結構，直接進入 Phase 4。

---

### Phase 4：PPTX 三件套 + 生成 `[FYI]`

> 詳見 `references/PPTX_PROTOCOL.md` Phase 3-4

輸出目錄：`docs/presentations/validated/YYYY/MM/YYYYMMDD_研究週報_{主題}/`

1. 撰寫三件套：`01_full_narrative_report.md`, `02_ppt_slide_outline.md`, `03_slide_layout_and_script.md`
2. 複製 config JSON 範本
3. 執行生成腳本
4. **截圖驗證（mandatory）**：渲染每張 slide 為 PNG，檢查 aspect ratio、文字溢出、顏色

**設計常數速查**（完整規範見 `references/PPTX_PROTOCOL.md`）：
- 色彩：dark=#1E2A44, bg=#F7F3EC, accent=#A85540, 綠=#009E73, 紅=#D55E00
- 字體：標題 Arial Bold ≥32pt, 內文 Arial ≥14pt, 下限 9pt
- 版面：視覺 55-65%, 留白 20-30%, 文字 ≤15%, 每頁 ≤4 bullet
- 雙語：中文主 + 英文副（60% 字號 + 縮排 0.25"）

---

### Phase 5：多 Agent 驗證 `[FYI → Review]`

> 詳見 `references/PPTX_PROTOCOL.md` Phase 5.5 + Phase 4.5

#### 5A. 內容驗證（12 Agent × 4 波次）

```
Wave 1 (內容)    Wave 2 (易讀)    Wave 3 (學術)    Wave 4 (PPTX)
 A1:邏輯         A4:敘事          A7:教授          A9:導演
 A2:質疑         A5:認知          A8:翻譯          A10:美術
 A3:佐證         A6:學生                           A11:審稿
                                                   A12:蒐集
```

每波次內平行執行（`feature-dev:code-reviewer` subagent_type），波次間依序修正。

#### 5B. PPTX 驗證（6 Agent × 2 波次）

- Wave 1（結構+視覺）：Agent-T(字體), Agent-C(色彩), Agent-L(佈局), Agent-B(雙語)
- Wave 2（內容+整合）：Agent-S(Speaker Notes), Agent-D(數據準確性)

**互動模式**：每波次展示驗證結果，等用戶確認修正。
**全自動模式**：自動修正可自動處理的問題，僅對需判斷的問題暫停。

---

### Phase 6：最終檢核與索引更新 `[Review]`

> 詳見 `references/PPTX_PROTOCOL.md` Phase 6

#### 雙層檢核清單

**Must-Pass（不通過不可發送）**：
- L1-L12：邏輯嚴謹度（結論附 evidence、因果鏈無跳步、結論/狀態表無矛盾、confound 控制、跨樣本標明、術語一致、圖表解釋性標題、Tier 分配一致）

**Should-Pass（品質）**：
- Q1-Q13：Layer 0 存在、Metadata 完整、絕對路徑、Paired/TO 分開、PPTX 驗證、索引更新

#### 索引更新

1. 更新 `docs/experiments/INDEX.md`
2. 更新 `docs/reports/` 下相關索引
3. 更新 `docs/CURRENT_FOCUS.md`（如有方向變更）
4. 建議更新 Memory（如有新結論/NEGATIVE/NO-GO）

---

## 圖表規範

- 初步觀察：分 Paired / TO 兩組
- 細緻觀察：2×2 per sample
- 固定順序：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829
- PPTX 圖片：fit-within + centered
- 寬圖/直式圖：自適應佈局

---

## 注意事項

1. **全自動模式下的報告**：完成後輸出「自動執行報告」，列出所有自動做的決策
2. **每個 Phase 產出的暫存文件**：存放在輸出目錄下，命名清楚
3. **用戶回饋優先**：任何自動判斷均可被用戶覆蓋
4. **不問空泛問題**：提供選項讓用戶選擇，而非開放式提問
5. **Layer 2 最重要**：花最多時間在證據卡和因果鏈上
6. **PPTX 截圖驗證不可跳過**：每次佈局修改後必須重新渲染驗證

---

## 過去週報範本

| 週報 | 特色 | 路徑 |
|------|------|------|
| 0310 | 三層方法分工 | `docs/reports/validated/2026/03/20260310_研究主線週報_20260305_20260310_01.md` |
| 0330 | HP Bug Fix + LOH | `docs/reports/validated/2026/03/20260330_研究週報_20260325_20260330_HP_bug_fix與LOH_evidence_panel_01.md` |
| LOH 簡報 | 偵探故事線 PPTX | `docs/presentations/validated/2026/04/20260401_LOH_weekly_report_draft/` |

## 研究脈絡速查

| 資訊 | 來源 |
|------|------|
| 當前狀態 + 阻塞 | `docs/CURRENT_FOCUS.md` |
| 實驗歷史索引 | `docs/experiments/INDEX.md` |
| 完整推論鏈 | `docs/reports/research_landscape/00_INDEX.md` |
| TO FP 問題全貌 | `docs/reports/research_landscape/01_TO_FP問題全貌.md` |
| ISM 分析價值界定 | `docs/reports/research_landscape/03_ISM分析價值界定.md` |
| 結論穩定性 | `docs/reports/research_landscape/06_結論穩定性審查.md` |
