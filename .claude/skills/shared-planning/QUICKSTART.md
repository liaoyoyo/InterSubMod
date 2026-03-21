# 快速啟動指南 - Shared Planning

5 分鐘快速上手跨平台任務規劃系統。

## 🎯 目標

讓 Claude Code 和 Antigravity 共享同一套任務規劃文件，實現無縫協作。

## ⚡ 三步驟啟動

### 步驟 1: 驗證安裝

```bash
# 確認檔案存在
ls .claude/skills/shared-planning/

# 應該看到:
# ├── SKILL.md
# ├── antigravity_agent.md
# ├── config.json
# ├── bridge.py (可執行)
# ├── README.md
# └── QUICKSTART.md (本檔案)

# 測試 bridge script
./.claude/skills/shared-planning/bridge.py --help
```

### 步驟 2: 建立第一個任務

#### 方法 A: 使用 Python Bridge (推薦初學者)

```bash
# 建立新任務
./.claude/skills/shared-planning/bridge.py init "my_first_task" \
  --type feature \
  --objective "學習如何使用跨平台規劃系統"

# 檢查建立結果
ls docs/shared_planning/my_first_task/
# 應該看到: task_plan.md, findings.md, progress.md
```

#### 方法 B: 在 Claude Code 中

```
User: "使用 shared-planning skill 建立新任務：學習跨平台協作"

Claude Code 會自動調用 skill 並建立任務
```

### 步驟 3: 更新任務

```bash
# 更新進度
./.claude/skills/shared-planning/bridge.py update "my_first_task" \
  --progress "完成系統安裝和測試"

# 記錄發現
./.claude/skills/shared-planning/bridge.py update "my_first_task" \
  --finding-title "系統運作正常" \
  --finding-context "初次測試" \
  --finding-discovery "所有功能都正常運作"

# 查看任務狀態
./.claude/skills/shared-planning/bridge.py list
```

## 🔄 跨平台協作流程

### 場景：Claude Code → Antigravity 協作

```bash
# === 在 Claude Code 中 ===
User: "建立新任務：實作新功能 X"
Claude: [建立任務並生成計劃]

# 提交到 Git
git add docs/shared_planning/
git commit -m "Planning: Initialize feature X"
git push

# === 切換到 Antigravity ===
# 在 Antigravity 中設定 Custom Agent:
# 1. 開啟 Antigravity
# 2. 建立新 Custom Agent
# 3. 複製 .claude/skills/shared-planning/antigravity_agent.md 內容
# 4. 貼上作為 system prompt
# 5. 命名為 "Planning Companion"

# 使用 Planning Companion
User: "Read the feature X task plan"
Agent: [讀取並理解計劃]

User: "Research implementation approaches"
Agent: [執行研究，更新 findings.md]

# 提交研究結果
git add docs/shared_planning/
git commit -m "Planning: Complete research for feature X"
git push

# === 回到 Claude Code ===
git pull
User: "查看 feature X 的研究結果並開始實作"
Claude: [讀取 findings.md，開始實作]
```

## 📋 常用指令速查

### Python Bridge 指令

```bash
# 初始化
./.claude/skills/shared-planning/bridge.py init <task_name> [--type TYPE] [--objective TEXT]

# 更新進度
./.claude/skills/shared-planning/bridge.py update <task_name> --progress "message"

# 新增發現
./.claude/skills/shared-planning/bridge.py update <task_name> \
  --finding-title "title" \
  --finding-context "context" \
  --finding-discovery "discovery"

# 完成任務項
./.claude/skills/shared-planning/bridge.py update <task_name> \
  --complete-task "task_pattern"

# 列出所有任務
./.claude/skills/shared-planning/bridge.py list

# 任務統計
./.claude/skills/shared-planning/bridge.py stats [task_name]

# 驗證檔案
./.claude/skills/shared-planning/bridge.py validate <task_name>
```

### Git 工作流程

```bash
# 開始工作前
git pull

# 定期提交
git add docs/shared_planning/
git commit -m "Planning: [描述]"

# 切換平台前
git push
```

## 🎓 實戰範例

### 範例 1: 簡單功能開發

```bash
# Day 1 - Claude Code: 規劃
./.claude/skills/shared-planning/bridge.py init "add_dark_mode" \
  --type feature \
  --objective "為應用程式新增深色模式"

# 編輯 task_plan.md 分解任務
# 提交
git add docs/shared_planning/add_dark_mode/
git commit -m "Planning: Initialize dark mode feature"
git push

# Day 2 - Antigravity: 研究
# 使用 Planning Companion agent
# 研究深色模式最佳實踐
# 更新 findings.md
# 提交研究結果

# Day 3 - Claude Code: 實作
git pull
# 根據研究結果實作
# 更新 progress.md
```

### 範例 2: Bug 修復

```bash
# 發現 Bug
./.claude/skills/shared-planning/bridge.py init "fix_memory_leak" \
  --type bugfix \
  --objective "修復甲基化分析中的記憶體洩漏"

# 記錄問題
./.claude/skills/shared-planning/bridge.py update "fix_memory_leak" \
  --finding-title "發現記憶體洩漏位置" \
  --finding-context "執行大型測試時記憶體持續增長" \
  --finding-discovery "RegionProcessor::process_single_region 未釋放暫存資料"

# 完成修復後
./.claude/skills/shared-planning/bridge.py update "fix_memory_leak" \
  --progress "完成記憶體洩漏修復，通過所有測試"

./.claude/skills/shared-planning/bridge.py update "fix_memory_leak" \
  --complete-task "修復記憶體洩漏"
```

### 範例 3: 實驗性研究

```bash
# 啟動實驗
./.claude/skills/shared-planning/bridge.py init "f1_score_optimization" \
  --type experiment \
  --objective "優化甲基化顯著性分析的 F1 分數"

# Antigravity: 並行參數掃描
# 使用 Manager View 同時測試多組參數
# 記錄所有結果到 findings.md

# Claude Code: 根據實驗結果調整實作
# 實作最佳參數配置
```

## 🔧 配置 Antigravity

### 一次性設定步驟

1. **開啟 Antigravity**

2. **建立 Custom Agent**
   - 點擊 "New Agent" 或類似按鈕
   - Agent 名稱: `Planning Companion`
   - 模型選擇: Gemini 3 Pro 或 Claude Sonnet 4.5

3. **匯入 Agent Prompt**
   ```bash
   # 複製檔案內容
   cat .claude/skills/shared-planning/antigravity_agent.md

   # 貼到 Antigravity 的 System Prompt 欄位
   ```

4. **測試 Agent**
   ```
   User: "List all planning tasks"
   Agent should read and list tasks from docs/shared_planning/
   ```

5. **完成！** 現在可以在 Antigravity 中使用 Planning Companion

## ❓ 常見問題

### Q: 為什麼要用這個系統而不是直接使用各平台的原生功能？

A:
- ✅ 統一格式：兩個平台共享相同的規劃文件
- ✅ 版本控制：所有計劃都在 Git 中追蹤
- ✅ 可移植：不依賴特定平台，隨時可切換
- ✅ 人類可讀：純 Markdown 格式，任何人都能閱讀編輯

### Q: Claude Code 和 Antigravity 同時編輯會衝突嗎？

A:
- 使用 Git 工作流程可以避免大部分衝突
- 建議切換平台前先 `git push`
- 如果發生衝突，手動合併（通常很簡單）

### Q: 可以只用其中一個平台嗎？

A:
- 可以！系統完全支援單一平台使用
- Python bridge 在任何環境都能運作
- 只是跨平台協作是設計的主要優勢

### Q: 如何遷移現有的計劃？

A:
```bash
# 建立新任務
./.claude/skills/shared-planning/bridge.py init "existing_task"

# 手動編輯生成的檔案，填入現有計劃內容
# 或使用腳本批次轉換
```

## 📚 下一步

完成快速啟動後，建議閱讀：

1. **[README.md](./README.md)** - 完整功能和範例
2. **[SKILL.md](./SKILL.md)** - Claude Code skill 詳細說明
3. **[antigravity_agent.md](./antigravity_agent.md)** - Antigravity agent 完整指南

## 🎉 你已經準備好了！

現在你可以：
- ✅ 在 Claude Code 中建立任務
- ✅ 在 Antigravity 中繼續工作
- ✅ 透過 Git 同步所有變更
- ✅ 使用 Python bridge 進行批次操作

開始建立你的第一個跨平台任務吧！

```bash
./.claude/skills/shared-planning/bridge.py init "my_awesome_project" \
  --type feature \
  --objective "建立一個令人驚豔的專案"
```

---

**需要幫助？** 查看 [README.md](./README.md) 的故障排除章節
**發現問題？** 在專案中建立 issue 或直接修改文檔
