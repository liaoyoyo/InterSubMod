---
title: weekly-report skill 升級為主流工作流（取代既有版本，串接新 myPPT）
date: 2026-05-01
status: outline-for-review
decision_path: 升級既有 weekly-report，命名不變，內部架構以新 W1-W7 為主
---

# 升級計劃：weekly-report → 主流週報工作流

## 0. 用戶決策（2026-05-01）

- 新架構為主、舊版為參考
- 命名保留 `weekly-report`（避免外部依賴中斷 + 觸發 keyword 不分裂）
- 取代既有版本，成為主流
- 串接新版 myPPT（母稿產出後 handoff）
- 先確認架構，再落地

## 1. 舊版 weekly-report 盤點

### 1.1 舊版檔案
```
InterSubMod/.claude/skills/weekly-report/
├── SKILL.md（236 行，7-Phase 架構）
└── references/
    ├── COLLECTION_PROTOCOL.md
    ├── LAYER_STRUCTURE.md
    └── PPTX_PROTOCOL.md
```

### 1.2 舊版 7-Phase 架構評估

| 舊 Phase | 內容 | 評估 | 命運 |
|----------|------|------|------|
| **Phase 0** 自動收集 | git log + Memory + KB + 上週 Layer 3-4 | ✅ 強，要保留 | **保留並強化**（變新 W1） |
| **Phase 1** 協作規劃 | 優先排序 + 敘事主軸 4 選 + 遺漏檢查 + PPTX 規劃 | ⚠ 主軸 4 選對應新「主線類型」，但 PPTX 規劃應移到 myPPT | **拆分**（前半 → 新 W2，後半 → myPPT） |
| **Phase 2** Layer 0 撰寫 | 宏觀問題 + 背景知識 + 前情提要 | ✅ 結構好 | **保留**（變新 W7 母稿的 Layer 0 段） |
| **Phase 3** Layer 1-4 主體 | 已建立知識 + 本週調查 + 整合更新 + 未來方向 | ✅ Layer 2 Thread 結構（問題+定義+假說+方法+證據卡+因果鏈+結論）非常強 | **保留**（變新 W7 母稿的主體 + 加 4 層分類標籤） |
| **Phase 3.5** Storyboard | Three-Act 故事弧 + 每頁 Storyboard | ✅ 強 | **移到 myPPT**（母稿後啟動） |
| **Phase 4** PPTX 三件套 | narrative + outline + script | ✅ 成熟 | **移到 myPPT** |
| **Phase 5** 多 Agent 驗證 | 12 Agent × 4 波次 + 6 Agent × 2 波次 | ✅ 強 | **拆分**（內容驗證留週報；PPTX 驗證移 myPPT） |
| **Phase 6** 檢核+索引 | Must-Pass L1-L12 + Should-Pass Q1-Q13 | ✅ 完整 | **拆分**（母稿 Must-Pass 留；PPTX 部分移 myPPT） |

### 1.3 舊版缺口（對照您的敘述需求）

| 缺口 | 說明 | 解決方式 |
|------|------|---------|
| ❌ 主線類型識別 | 4 選 1：進展/問題/求協助/探索 | 新增 W2 階段 |
| ❌ 內容 4 層分類 [F]/[O]/[I]/[U] | 與 Tier 1/2/3（優先序）維度不同：F/O/I/U 是「真實性」 | 新增 W3 階段（與 Tier 並用） |
| ❌ 5-7 個教授追問預測 | 舊版無專門段落 | 新增 W6 階段 |
| ❌ 過度宣稱紅旗檢查 | 舊版有「不問空泛問題」原則但無語氣紅旗 | 新增 W5 階段 |
| ❌ 每輪 ≤ 5 問規範 | 舊版有「不問空泛問題」但無數量限制 | 強化 prompts/ 設計 |
| ❌ myPPT 銜接協議 | 舊版直接 Phase 4 產 PPTX，無母稿 ↔ PPT 切割 | 新增 handoff 機制 |

## 2. 新架構：W1-W7（取代舊 Phase 0-3.5）

### 2.1 7 階段對照映射

```
舊 Phase 0      → 新 W1 Raw Data 收集 + 用戶補充     (保留 + 強化)
舊 Phase 1 前半 → 新 W2 主線類型識別                 (升級 4 選 1)
[新增]          → 新 W3 內容 4 層分類 [F]/[O]/[I]/[U]  (新增階段)
舊 Phase 1 後半 → 新 W4 重點排序 + 4 桶分流          (升級評分)
[新增]          → 新 W5 邏輯紅旗檢查                  (新增階段)
[新增]          → 新 W6 教授問答預測                  (新增階段)
舊 Phase 2-3    → 新 W7 母稿產出（Layer 0-4 + 17 段 mapping）  (保留 + 重組)

舊 Phase 3.5    → 移到 myPPT (Storyboard)
舊 Phase 4      → 移到 myPPT (PPTX 三件套)
舊 Phase 5      → 拆分（內容驗證留 W7；PPTX 驗證移 myPPT）
舊 Phase 6      → 拆分（母稿索引留 W7；PPTX 索引移 myPPT）
```

### 2.2 母稿格式：Layer 0-4 + 17 段 hybrid

**保留舊版 Layer 結構為主骨架，17 段作為 Layer 內部標籤**：

```
母稿 .md 結構

# §1 主線（≤ 30 字）+ §2 一句話重點   ← Layer 0.1 宏觀問題定位

## Layer 0.2 背景知識（最多 3 概念群組）  ← 保留
## Layer 0.3 前情提要                    ← 保留

## Layer 1 已建立知識參考（薄）           ← 保留

## Layer 2 本週調查（核心）
   每個 Thread 仍用舊格式：
   ### Thread A: 名稱
   #### 問題陳述
   #### 定義區塊
   #### 假說
   #### 方法
   #### 證據卡（Tier 1/2/3）+ 4 層分類標籤 [F]/[O]/[I]/[U]   ← 新增 4 層標籤
   #### 因果鏈（Mermaid）
   #### 結論

   對應 17 段：§3 [F]、§4 [O][I]、§5 [U]

## §7-§8 重點排序與順序  ← Layer 2 結尾整合，新增段落

## Layer 3 整合更新
   §11 需補資料、§12 需製圖、§13 需補定義、§14 講稿例子

## Layer 4 未來方向
   §16 下一步行動、§17 教授可能提問與回答準備   ← 新增 §17

## §6 不放 PPT 暫存內容 + §15 暫存紀錄  ← 新增段落
```

**4 層分類 vs Tier 並用規則**：

| 維度 | 含義 | 範圍 |
|------|------|------|
| **Tier 1/2/3**（舊有） | 證據卡的優先序 / 重要性 | 決定該證據在 Layer 2 的呈現深度 |
| **[F]/[O]/[I]/[U]**（新增） | 證據的真實性層級 | 決定描述語氣（事實 vs 推論 vs 待確認） |

範例：一筆素材可同時是 Tier 1（最重要）+ [O]（初步觀察，N=3 樣本未驗證至 7）。
- 在 Layer 2 用 Tier 1 完整呈現
- 但描述用 [O] 語氣（「目前觀察到 X，需 7 樣本驗證後確認」），不寫成 [F]

### 2.3 與 myPPT 銜接點

```
W7 母稿 C6 確認 → AskUserQuestion：
  ├─ A. 銜接 myPPT 產 PPTX（自動讀母稿，跳過 main thesis 鎖定）
  ├─ B. 母稿留檔，下次再用 /myPPT --from-draft <path>
  └─ C. 母稿即終點（如 internal status check）

myPPT 從 W7 母稿讀取的欄位：
  - frontmatter: report_type, main_statement, audience
  - Layer 0.1: thesis bar 主視覺
  - Layer 2 Thread × N: 各對應 1-2 張 slide
  - Layer 3: integration slide
  - Layer 4: future tree slide
  - §17 教授追問: backup slides / Q&A slides
```

## 3. 新檔案結構（升級後）

```
InterSubMod/.claude/skills/weekly-report/      # 命名不變
├── SKILL.md                       # 重寫（< 250 行）
├── playbook.md                    # 新增（~800-1000 行，含 W1-W7 細則）
├── README.md                      # 新增（skill 入口）
├── references/
│   ├── COLLECTION_PROTOCOL.md     # 保留（舊版，作 W1 細則）
│   ├── LAYER_STRUCTURE.md         # 升級（保留 Layer 0-4 + 加 17 段 mapping + 加 4 層分類規則）
│   ├── PPTX_PROTOCOL.md           # **DEPRECATED**（內容移到 myPPT）
│   └── HANDOFF_TO_MYPPT.md        # 新增（母稿 schema + 跳過項目 + AskUserQuestion 模板）
├── templates/                     # 新增
│   ├── progress_focus.md          # 進展型主線
│   ├── problem_focus.md           # 問題型主線
│   ├── advisor_consult.md         # 求協助型主線
│   └── new_direction_explore.md   # 探索型主線
├── prompts/                       # 新增
│   ├── raw_data_collect.md        # W1 → C0
│   ├── main_thread_identify.md    # W2 → C1
│   ├── content_classify_4tier.md  # W3 → C2
│   ├── priority_sort_4bucket.md   # W4 → C3
│   ├── logic_check_redflag.md     # W5 → C4
│   ├── professor_qa_predict.md    # W6 → C5
│   ├── master_draft_finalize.md   # W7 → C6
│   └── handoff_to_myppt.md        # C6 後選擇繼續 PPT
└── examples/                      # 新增
    └── master_draft_example.md    # 1 份完整範例母稿
```

**對比舊版：保留 1 + 升級 1 + 棄用 1 + 新增 N**

| 檔案 | 命運 | 說明 |
|------|------|------|
| 舊 SKILL.md | **重寫** | 從 7-Phase → W1-W7 |
| 舊 references/COLLECTION_PROTOCOL.md | **保留** | W1 細則 |
| 舊 references/LAYER_STRUCTURE.md | **升級** | + 17 段 mapping + 4 層分類 |
| 舊 references/PPTX_PROTOCOL.md | **DEPRECATED** | 內容遷移到 myPPT/playbook |
| 新增 14 檔 | **新增** | playbook + README + 4 templates + 8 prompts + 1 example + handoff_protocol |

## 4. 主流流程圖（升級版）

```
[U: 我要做週報] → trigger /weekly-report
    ↓
W1 Raw Data 收集（自動掃描 + 用戶補充）  → C0 確認
    ↓
W2 主線類型識別（4 選 1 + main statement ≤ 30 字）  → C1 確認
    ↓
W3 內容 4 層分類（每筆素材標 [F]/[O]/[I]/[U]）  → C2 確認
    ↓
W4 重點排序 + 4 桶分流（PPT/講稿/備註/暫存）  → C3 確認
    ↓
W5 邏輯紅旗檢查（過度宣稱/流水帳/教授視角缺）  → C4 確認
    ↓
W6 教授問答預測（5-7 個追問 + 預備回答）  → C5 確認
    ↓
W7 17 段母稿產出（Layer 0-4 + 17 段 mapping）  → C6 確認
    ↓
[Output: InterSubMod/docs/reports/validated/YYYY/MM/YYYYMMDD_週報_週主題/
              master_draft.md]
    ↓
[Q: 銜接 myPPT 產 PPTX？]
    ├─ A. 是 → /myPPT --from-draft <path>
    │         ↓
    │         (myPPT 讀母稿，跳過 main thesis；保留 outline/section/slide checkpoint)
    │         ↓
    │         產出 PPTX + speaker script + audit
    │
    ├─ B. 母稿留檔
    └─ C. 母稿即終點
```

## 5. 待您確認的細節（10 個關鍵問題）

請逐項確認：

| # | 問題 | 我的預設 | 您的選擇 |
|:-:|------|--------|----------|
| Q1 | 命名保留 `weekly-report` 還是改 `weekly-master-draft`？ | 保留 weekly-report | _ |
| Q2 | 母稿主骨架用 Layer 0-4 還是 17 段？ | **Layer 0-4 為主，17 段為內部標籤** | _ |
| Q3 | 4 層分類 [F]/[O]/[I]/[U] 與 Tier 1/2/3 並用？還是只留一個？ | **並用**（不同維度：真實性 vs 重要性） | _ |
| Q4 | 舊 Phase 5 多 Agent 驗證（12 + 6 Agent）保留多少？ | 內容驗證 12 Agent 保留；PPTX 驗證 6 Agent 移 myPPT | _ |
| Q5 | 舊 references/PPTX_PROTOCOL.md 移到 myPPT 哪個檔？ | 拆進 myPPT/playbook.md §5/§6/§7 + myPPT/style_library | _ |
| Q6 | 母稿輸出路徑沿用舊版 `docs/reports/validated/YYYY/MM/` 還是新建 `docs/weekly_reports/`？ | **沿用舊版**（避免文件樹分裂） | _ |
| Q7 | C0-C6 七個 checkpoint 是否過多？要合併嗎？ | 7 個剛好對應 7 階段，**不合併** | _ |
| Q8 | Storyboard 階段（舊 Phase 3.5）整體移到 myPPT 還是留週報？ | **移到 myPPT**（簡報視覺敘事屬 myPPT 範疇） | _ |
| Q9 | fast-track（全自動）模式跳過哪些 checkpoint？ | C0/C2/C3 用 AI 預設；C1/C5/C6 必停 | _ |
| Q10 | 升級時舊 SKILL.md 要備份還是直接覆蓋？ | **備份為 SKILL.md.v1.bak**，覆蓋實作 | _ |

## 6. 落地步驟（用戶確認 Q1-Q10 後執行）

```
Step 1：備份舊版（如 Q10=備份）
  cp SKILL.md SKILL.md.v1.bak

Step 2：寫 4 份大綱（已寫 → InterSubMod/.claude/skills/weekly-master-draft/outlines/）
  搬移到 InterSubMod/.claude/skills/weekly-report/outlines/
  並依用戶答案修改

Step 3：依答案落地實際檔案
  - 重寫 SKILL.md
  - 升級 references/LAYER_STRUCTURE.md
  - 移走 references/PPTX_PROTOCOL.md
  - 新增 14 檔 (playbook + README + templates + prompts + example + handoff)

Step 4：myPPT 補 --from-draft 接口
  - myPPT/prompts/outline_confirm.md 加母稿讀取分支
  - myPPT/playbook.md §1 加「若有母稿，跳過 main thesis 鎖定」

Step 5：reviewer 抽查
  - skill-reviewer 檢查升級後 SKILL.md 觸發強度不變
  - 模擬 conversation 測試「整理本週」自動載入

Step 6：刪除中介 outlines/ 子資料夾（成熟後）
```

## 7. 風險

| ID | 風險 | 緩解 |
|----|------|------|
| R1 | 舊版用戶習慣被打斷 | Phase mapping 表清楚標示舊→新對應；保留 Layer 0-4 結構 |
| R2 | W3 4 層分類與 Tier 混淆 | playbook §4 明確區分維度：真實性 vs 重要性 |
| R3 | C0-C6 七個 checkpoint 太繁瑣 | fast-track 模式跳過 C0/C2/C3 |
| R4 | myPPT 還未落地，handoff 規格可能變 | handoff schema 在兩 skill outlines/ 雙向引用，落地時對齊 |
| R5 | 舊 PPTX_PROTOCOL 內容大，搬遷可能出錯 | 搬遷分 2 步：先 myPPT/playbook 引用 → 確認內容齊全後刪舊檔 |
