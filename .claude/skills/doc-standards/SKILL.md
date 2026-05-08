---
name: doc-standards
description: InterSubMod 文檔管理規範。檔案命名格式、元數據模板、圖片規則、多步驟目錄結構、AI session 報告規則。觸發：建立新 .md 檔案或重組文件目錄時。DO NOT USE WHEN：寫 SKILL.md（用 plugin-dev:skill-development）、寫 13 段技術報告（用 structured-tech-report）、寫週報（用 weekly-report）。
allowed-tools: Read, Write, Glob, Grep
user-invocable: true
paths: ["docs/**/*.md", "research/**/*.md", "**/*.md"]
---

# InterSubMod 文檔管理規範

## 文檔目錄結構

```
docs/
├── architecture/        # 專案主軸架構說明
├── concepts/            # 構想紀錄（研究大圖景、發展樹、理論基礎）
├── plans/               # 計劃文件（YYYY/MM 分層）
├── solutions/           # 問題解決報告（任務目標/YYYY/MM 分層）
├── experiments/         # 實驗紀錄（in_progress / validated / finalized）
├── reports/             # 正式說明文件（validated / finalized，索引用 00_INDEX.md）
│   └── {topic}/         # 多檔案說明文件
│       └── figures/     # 說明文件圖表
├── methodology/         # 方法學審查與 closeout 報告
├── standards/           # 命名規範、狀態管理、output 規範等治理文件
├── decisions/           # 重整與治理決策紀錄
├── presentations/       # 簡報 PDF
├── provenance/          # AI 對話 provenance
│   └── ai_sessions/     # AI 對話紀錄（YYYY/MM 分層）
├── data_specs/          # 數據規格
├── references/          # 參考資料（含 manual/、external/）
├── research/            # 長期研究主題工作區
│   └── {study_name}/    # 單一研究主題
│       ├── README.md    # 研究概述與索引
│       ├── figures/     # 圖表（統一命名，gitignore）
│       ├── data/        # 中間數據（gitignore）
│       ├── scripts/     # 分析腳本（git 追蹤）
│       └── reports/     # 研究報告
├── refactor_baseline/   # 重構基準數據
├── archive/             # 歸檔
│   └── deep/            # 深度歸檔（需查詢歷史）
└── trash/               # 暫存待刪除
```

## 檔案命名格式

```
{YYYYMMDD}_{中文說明目標}_{流水號}.md
```

範例：`20260111_文檔庫重整計劃_01.md`

## 多步驟研究專案目錄結構

當研究/任務包含多個 Step 或子計劃時，在對應目錄下建立**專案子資料夾**：

```
docs/plans/YYYY/MM/{YYYYMMDD}_{專案主題}/
  ├── 00_總覽與執行順序.md        # 索引：波次、優先序、依賴關係
  ├── Step1_{子計劃主題}.md
  ├── Step2_{子計劃主題}.md
  └── ...

docs/architecture/{YYYYMMDD}_{專案主題}/
  ├── {YYYYMMDD}_{架構說明}_{流水號}.md
  └── {YYYYMMDD}_{資料追蹤表}_{流水號}.md
```

**規則**：
- 資料夾名稱格式：`{YYYYMMDD}_{中文專案主題}`
- 索引檔用 `00_` 開頭
- Step 計劃用 `StepN_` 前綴（方便排序）
- 同一研究的架構文件、追蹤表集中在同一專案子目錄
- 可按日期或任務主題搜尋定位

## 檔案元數據

每個 Markdown 檔案開頭需包含：

```markdown
<!--
建立時間: YYYY-MM-DD HH:MM
目標: [本檔案的目標或用途]
處理範圍: [涵蓋的工作範圍]
關聯檔案:
  - [相關檔案路徑 1]
  - [相關檔案路徑 2]
-->
```

## 圖片存放與引用規則

- 純圖片子目錄統一命名 `figures/`（不用 `images/`、`plots/`）；混合類型資源目錄可用 `assets/`
- .md 引用圖片用相對路徑：`figures/xxx.png`
- 相對路徑最多 2 層（禁止 `../../../`）；引用 `output/` 或 `research/` 圖片時允許超過 2 層
- 圖片命名：`{NN}_{英文描述}.png`

## AI 對話紀錄撰寫

每次 AI 對話完成重要任務後，撰寫執行報告（可使用 `/report` skill）：

1. **報告位置**：`docs/provenance/ai_sessions/YYYY/MM/`
2. **檔案格式**：`{YYYYMMDD}_{對話主題}_{流水號}.md`
3. **必要內容**：
   - 對話目標
   - 關鍵決策
   - 修改的檔案清單
   - 後續行動建議
