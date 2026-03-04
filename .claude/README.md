# Claude Code 自動化配置說明

本目錄包含 InterSubMod 專案的 Claude Code 自動化配置，包括 SubAgent、Commands、Skills 和 Hooks。

## 目錄結構

```
.claude/
├── settings.json              # 全域設定
├── settings.local.json        # 本地設定（權限、Hooks）
├── CLAUDE.md                  # 專案開發指南
├── README.md                  # 本說明文件
├── agents/                    # 7 個 SubAgent 子代理
│   ├── researcher.md          # 資料收集
│   ├── architect.md           # 任務規劃
│   ├── developer.md           # 程式撰寫
│   ├── optimizer.md           # 程式優化
│   ├── tester.md              # 測試驗證
│   ├── reviewer.md            # 數據審查
│   └── release.md             # 發布管理
├── commands/                  # Slash Commands
│   ├── test-quick.md          # 快速測試
│   ├── test-data.md           # 數據測試
│   ├── test-full.md           # 完整測試
│   ├── git-start.md           # 開始新任務
│   ├── git-commit.md          # 提交變更
│   ├── git-finish.md          # 完成任務
│   └── build.md               # 編譯專案
└── skills/                    # Skills 技能
    └── report/
        └── SKILL.md           # 報告撰寫技能
```

---

## SubAgent 子代理

### 可用的子代理

| Agent | 名稱 | 用途 | 輸出目錄 |
|-------|------|------|----------|
| researcher | 資料收集 | 搜尋學術論文、GitHub、技術文件 | docs/references/ |
| architect | 任務規劃 | 拆分任務、撰寫計劃書 | docs/plans/ |
| developer | 程式撰寫 | 根據計劃書撰寫程式碼 | src/ |
| optimizer | 程式優化 | 審查程式碼品質、效能 | docs/solutions/ |
| tester | 測試驗證 | 執行各級別測試 | docs/experiments/ |
| reviewer | 數據審查 | 分析數據結果、驗證假設 | docs/experiments/ |
| release | 發布管理 | Git 操作、Docker、版本管理 | - |

### 使用方式

在對話中提及子代理名稱，Claude 會自動調用對應的專業角色。

例如：
- "請使用 researcher 搜尋關於甲基化距離計算的論文"
- "請 architect 幫我規劃這個新功能的實作計劃"

---

## Slash Commands

### 可用命令

| 命令 | 說明 |
|------|------|
| `/test-quick` | 快速單點測試 (< 30 秒) |
| `/test-data` | 完整數據測試 (約 1 分鐘) |
| `/test-full` | 完整流程測試 (約 5 分鐘) |
| `/git-start [type] [desc]` | 開始新任務，建立功能分支 |
| `/git-commit [type] [msg]` | 提交變更 |
| `/git-finish` | 完成任務，合併分支 |
| `/build` | 編譯專案 |

### 使用範例

```
/test-quick
/git-start feature add-new-distance
/git-commit feat "add bernoulli distance calculation"
/git-finish
```

---

## Skills 技能

### /report - 報告撰寫

用於撰寫 AI 對話執行報告，記錄重要決策和修改內容。

**使用時機**：
- 完成重要開發任務後
- 做出關鍵技術決策後
- 會話結束前

**輸出位置**：`docs/ai_sessions/{YYYY}/{MM}/`

---

## Hooks 自動化

### 現有 Hooks

| Hook 類型 | 觸發條件 | 動作 |
|-----------|----------|------|
| PreToolUse | 執行危險命令 | 顯示警告 |
| PostToolUse | 編輯 C++ 檔案 | 提醒編譯和測試 |
| PostToolUse | git commit | 提醒確認測試和文檔 |
| SubagentStop | 子代理完成 | 提醒檢查輸出 |
| Stop | 會話結束 | 提醒撰寫執行報告 |

---

## MCP 整合

### PubMed MCP

用於搜尋生物醫學文獻。

**配置環境變數**：
```bash
export NCBI_API_KEY="your_api_key_here"
```

**測試**：
```bash
python3 scripts/mcp/pubmed_server.py --test methylation cancer
```

---

## 開發工作流程

### 典型開發流程

1. **開始任務**
   ```
   /git-start feature add-new-feature
   ```

2. **規劃任務**
   - 使用 architect 子代理規劃實作步驟

3. **撰寫程式碼**
   - 使用 developer 子代理撰寫程式碼
   - 自動觸發編譯提醒

4. **測試驗證**
   ```
   /test-quick    # 快速驗證
   /test-data     # 完整數據測試
   ```

5. **程式碼審查**
   - 使用 optimizer 子代理審查程式碼品質

6. **提交變更**
   ```
   /git-commit feat "add new feature"
   ```

7. **完成任務**
   ```
   /git-finish
   ```

8. **撰寫報告**
   - 使用 /report 技能撰寫執行報告

---

## 注意事項

1. **API Key 安全**：NCBI API Key 應存放在環境變數中，不要直接寫入配置文件
2. **測試順序**：建議按照 quick → data → full 的順序執行測試
3. **報告撰寫**：每次重要對話結束前，使用 /report 技能撰寫執行報告
4. **分支管理**：遵循 Git 分支命名規範 (feature/, fix/, docs/, refactor/)
