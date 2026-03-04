# Shared Planning - 跨平台任務規劃系統

**Version**: 1.0.0
**Last Updated**: 2026-01-22

一個可以同時被 **Claude Code** 和 **Google Antigravity** 使用的任務規劃系統。

## 🎯 核心特性

- ✅ **跨平台兼容**: Claude Code 和 Antigravity 共享同一套規劃文件
- ✅ **標準化格式**: 統一的 Markdown 格式，易於版本控制
- ✅ **自動同步**: 兩個平台都能讀寫相同的文件
- ✅ **環境偵測**: 自動識別當前運行平台
- ✅ **Python Bridge**: 獨立的 Python 工具，兩個平台都能使用

## 📁 檔案結構

```
.claude/skills/shared-planning/
├── README.md                   # 本檔案
├── SKILL.md                    # Claude Code skill 定義
├── antigravity_agent.md        # Antigravity agent prompt
├── config.json                 # 配置檔案
├── bridge.py                   # Python bridge script (可執行)
└── templates/                  # 範本目錄 (可選)

docs/shared_planning/{TASK_NAME}/
├── task_plan.md               # 任務計劃
├── findings.md                # 研究發現
└── progress.md                # 進度日誌
```

## 🚀 快速開始

### 方法 1: 在 Claude Code 中使用 (推薦)

```
User: "建立新任務規劃：實作 Bhattacharyya distance"

Claude Code 會自動調用 /shared-planning skill
```

### 方法 2: 使用 Python Bridge (兩個平台通用)

```bash
# 初始化新任務
python .claude/skills/shared-planning/bridge.py init "bhattacharyya_distance" \
  --type feature \
  --objective "實作 Bhattacharyya distance 作為新的甲基化距離度量"

# 更新進度
python .claude/skills/shared-planning/bridge.py update "bhattacharyya_distance" \
  --progress "完成數學公式推導"

# 新增研究發現
python .claude/skills/shared-planning/bridge.py update "bhattacharyya_distance" \
  --finding-title "找到最佳化實作方法" \
  --finding-context "研究現有 C++ 實作" \
  --finding-discovery "Eigen library 提供高效向量運算"

# 標記任務完成
python .claude/skills/shared-planning/bridge.py update "bhattacharyya_distance" \
  --complete-task "文獻調研"

# 列出所有任務
python .claude/skills/shared-planning/bridge.py list

# 查看任務統計
python .claude/skills/shared-planning/bridge.py stats "bhattacharyya_distance"

# 驗證檔案格式
python .claude/skills/shared-planning/bridge.py validate "bhattacharyya_distance"
```

### 方法 3: 在 Antigravity 中使用

#### 選項 A: 使用 Custom Agent

1. 開啟 Antigravity
2. 建立新的 Custom Agent
3. 複製 `antigravity_agent.md` 的內容作為 system prompt
4. 命名為 "Planning Companion"
5. 開始使用：

```
User: "Create a new task for implementing Bhattacharyya distance"

Antigravity Planning Companion will:
- Read the agent instructions
- Create task files in docs/shared_planning/
- Follow the same format as Claude Code
```

#### 選項 B: 直接使用 Python Bridge

在 Antigravity 的終端中：

```bash
python .claude/skills/shared-planning/bridge.py init "task_name"
```

Antigravity 可以透過 terminal access 或 code execution 功能使用 bridge.py。

## 📖 使用情境範例

### 情境 1: Claude Code 開始，Antigravity 接手研究

```bash
# 在 Claude Code 中
User: "建立新任務：優化 PERMANOVA 演算法"
Claude: [調用 /shared-planning skill，建立任務]

# 提交到 Git
git add docs/shared_planning/
git commit -m "Planning: Initialize PERMANOVA optimization task"
git push

# 切換到 Antigravity
# Antigravity 使用 Planning Companion agent
User: "Read the PERMANOVA optimization task plan"
Agent: [讀取 task_plan.md，了解任務]

User: "Research PERMANOVA optimization techniques"
Agent: [使用 Manager View 並行研究]
       [更新 findings.md 記錄發現]
       [更新 progress.md 記錄進度]

# 提交研究成果
git add docs/shared_planning/permanova_optimization/
git commit -m "Planning: Complete PERMANOVA research phase"
git push

# 回到 Claude Code 實作
User: "查看 PERMANOVA 優化的研究結果"
Claude: [讀取 findings.md]
        [根據研究結果實作程式碼]
```

### 情境 2: 兩個平台並行工作不同任務

```bash
# Claude Code: 實作核心功能
User: "實作新的距離度量函數"
Claude: [建立任務: distance_metrics]
        [實作 C++ 程式碼]
        [更新 progress.md]

# 同時，Antigravity: 實驗參數調整
User: "Run parameter sweep for F1 optimization"
Agent: [建立任務: f1_parameter_sweep]
       [使用 Manager View 並行測試多個參數組合]
       [記錄結果到 findings.md]

# 兩個任務互不干擾，各自在 docs/shared_planning/ 下有獨立目錄
```

### 情境 3: 純 Python Bridge 工作流程

```bash
# 適用於需要腳本化或批次處理的場景

# 初始化多個相關任務
for task in "data_preprocessing" "model_training" "evaluation"; do
  python .claude/skills/shared-planning/bridge.py init "$task" --type research
done

# 查看所有任務狀態
python .claude/skills/shared-planning/bridge.py list

# 在腳本中自動更新進度
python .claude/skills/shared-planning/bridge.py update "data_preprocessing" \
  --progress "Completed data cleaning: $(date)"

# 產生統計報告
python .claude/skills/shared-planning/bridge.py stats > task_summary.txt
```

## 🔧 配置說明

### config.json 主要設定

```json
{
  "paths": {
    "base_dir": "docs/shared_planning",     # 任務根目錄
    "task_plan_file": "task_plan.md",       # 計劃檔名
    "findings_file": "findings.md",         # 發現檔名
    "progress_file": "progress.md"          # 進度檔名
  },

  "platforms": {
    "claude_code": {
      "identifier": "Claude Code",          # 平台識別字
      "env_var": "CLAUDE_CODE"              # 環境變數
    },
    "antigravity": {
      "identifier": "Antigravity",
      "env_var": "ANTIGRAVITY_SESSION"
    }
  }
}
```

### 自訂配置

如果需要修改配置，編輯 `.claude/skills/shared-planning/config.json`。

## 📋 檔案格式說明

### task_plan.md - 任務計劃

```markdown
# Task Plan: {TASK_NAME}

**Created**: 2026-01-22 10:00 (via Claude Code)
**Updated**: 2026-01-22 15:30 (via Antigravity)
**Status**: In Progress
**Type**: feature

## Objective
清楚描述任務目標

## Background
背景資訊與前置條件

## Tasks
### Phase 1: 規劃
- [x] 完成的任務
  - Completed: 2026-01-22 (via Claude Code)
- [ ] 待辦任務
  - Assigned to: Antigravity
  - Priority: High

## Dependencies
- Depends on: [其他任務]
- Blocks: [被阻擋的任務]

## Resources
- Documentation: [連結]
- Related Files: `src/file.cpp`
```

### findings.md - 研究發現

```markdown
# Findings: {TASK_NAME}

**Created**: 2026-01-22
**Last Updated**: 2026-01-22 (via Antigravity)

## Key Discoveries

### 2026-01-22: 發現標題

**Context**: 研究背景

**Discovery**: 具體發現

**Evidence**:
- Code: `src/core/Stats.cpp:123`
- Test: 測試結果連結
- Data: 數據分析連結

**Implications**:
- 這如何影響我們的方法？
- 可以做出什麼決策？

**Action Items**:
- [ ] 後續行動 1
- [ ] 後續行動 2
```

### progress.md - 進度日誌

```markdown
# Progress Log: {TASK_NAME}

**Started**: 2026-01-22
**Current Status**: 60% complete | In Progress
**Last Updated**: 2026-01-22 (via Claude Code)

---

## 2026-01-22 Wednesday

**Session Info**:
- Platform: Claude Code
- Session time: 14:00 - 17:00
- Focus area: 實作核心功能

### ✅ Completed Today
- **距離計算函數**: 完成 Bhattacharyya distance 實作
  - Files: `src/core/DistanceMetrics.cpp`, `include/core/DistanceMetrics.hpp`
  - Tests: `tests/test_distance_metrics.cpp`
  - Commit: abc1234

### 🔄 In Progress
- **單元測試**: 撰寫邊界條件測試
  - Progress: 70%
  - Next: 完成零向量和負值測試案例

### 📋 Planned Next
- **整合測試**: 整合到 RegionProcessor
- **效能測試**: Benchmark 與現有度量比較

### 💡 Key Insights
- 使用 Eigen 的 vectorized operations 提升效能 3x
- 參考 findings.md "2026-01-22: Eigen 最佳化"

### ⏱️ Time Allocation
- Implementation: 2h
- Testing: 1h
```

## 🔄 Git 工作流程建議

### 日常工作流程

```bash
# 開始工作前
git pull

# 工作中（定期）
git add docs/shared_planning/
git commit -m "Planning: Update progress for {task_name}"

# 切換平台前
git push

# 重要里程碑
git add docs/shared_planning/
git commit -m "Planning: Complete research phase for {task_name}"
git push
git tag -a "planning/{task_name}/research-complete" -m "Research phase completed"
```

### 衝突處理

如果兩個平台同時編輯導致衝突：

```bash
# 1. 拉取遠端變更
git pull

# 2. 手動解決衝突
# - progress.md: 合併兩邊的進度更新（按時間順序）
# - task_plan.md: 保留最新的任務狀態
# - findings.md: 合併所有發現

# 3. 標記為已解決
git add docs/shared_planning/
git commit -m "Planning: Merge updates from both platforms"
git push
```

## 🛠️ 故障排除

### 問題 1: bridge.py 無法執行

```bash
# 確認執行權限
chmod +x .claude/skills/shared-planning/bridge.py

# 確認 Python 版本
python3 --version  # 需要 Python 3.7+

# 直接使用 python3 執行
python3 .claude/skills/shared-planning/bridge.py --help
```

### 問題 2: 平台偵測錯誤

```bash
# 手動設定環境變數
export CLAUDE_CODE=1
# 或
export ANTIGRAVITY_SESSION=1

# 然後執行 bridge.py
```

### 問題 3: 檔案格式驗證失敗

```bash
# 驗證任務檔案
python .claude/skills/shared-planning/bridge.py validate "task_name"

# 如果檔案損壞，從 Git 恢復
git checkout docs/shared_planning/{task_name}/

# 或重新初始化
rm -rf docs/shared_planning/{task_name}/
python .claude/skills/shared-planning/bridge.py init "task_name"
```

### 問題 4: Claude Code skill 無法調用

```bash
# 確認 skill 已正確安裝
ls .claude/skills/shared-planning/SKILL.md

# 在 Claude Code 中重新載入 skills
/reload-skills  # (如果有這個指令)

# 或直接使用 Python bridge 作為替代
python .claude/skills/shared-planning/bridge.py init "task_name"
```

## 📊 進階功能

### 任務統計報告

```bash
# 查看所有任務統計
python .claude/skills/shared-planning/bridge.py stats

# 輸出為 JSON 方便後續處理
python .claude/skills/shared-planning/bridge.py stats --output json > stats.json

# 查看特定任務詳細統計
python .claude/skills/shared-planning/bridge.py stats "task_name"
```

### 批次操作

```bash
# 列出所有進行中的任務
python .claude/skills/shared-planning/bridge.py list | grep "In Progress"

# 驗證所有任務
for task in $(ls docs/shared_planning/); do
  python .claude/skills/shared-planning/bridge.py validate "$task"
done

# 產生每日摘要
cat <<EOF > daily_summary.sh
#!/bin/bash
echo "Daily Planning Summary - $(date)"
echo "================================"
python .claude/skills/shared-planning/bridge.py list
echo ""
echo "Statistics:"
python .claude/skills/shared-planning/bridge.py stats
EOF
chmod +x daily_summary.sh
./daily_summary.sh
```

### 與其他工具整合

#### Slack 通知

```bash
# 當完成任務時發送 Slack 通知
python .claude/skills/shared-planning/bridge.py update "task_name" \
  --complete-task "實作核心功能" && \
curl -X POST -H 'Content-type: application/json' \
  --data '{"text":"Task completed: 實作核心功能"}' \
  YOUR_SLACK_WEBHOOK_URL
```

#### GitHub Issues 同步

```bash
# 從 GitHub issue 建立任務
gh issue view 123 --json title,body | \
  jq -r '"Issue #123: " + .title' | \
  xargs -I {} python .claude/skills/shared-planning/bridge.py init "{}"
```

## 🎨 最佳實踐

### 1. 命名規範

```bash
# 任務名稱使用底線分隔，小寫
✅ good: "bhattacharyya_distance"
✅ good: "f1_score_optimization"
❌ bad: "Bhattacharyya Distance"
❌ bad: "f1-score-optimization"
```

### 2. 提交訊息格式

```bash
# 使用統一的前綴
git commit -m "Planning: {action} for {task_name}"

# 範例
git commit -m "Planning: Initialize bhattacharyya_distance task"
git commit -m "Planning: Update progress for f1_optimization"
git commit -m "Planning: Complete research phase for permanova_fix"
```

### 3. 任務分解原則

```markdown
# 好的任務分解
- [ ] Phase 1: Research (1-2 days)
  - [ ] Literature review
  - [ ] Code analysis
- [ ] Phase 2: Implementation (3-5 days)
  - [ ] Core algorithm
  - [ ] Unit tests
- [ ] Phase 3: Integration (1-2 days)
  - [ ] Integration tests
  - [ ] Documentation

# 避免過於籠統
❌ - [ ] Complete the task
❌ - [ ] Do research
❌ - [ ] Write code
```

### 4. 發現記錄準則

```markdown
# 好的發現記錄
✅ 包含具體的程式碼位置: `src/core/Stats.cpp:456`
✅ 附上證據連結或數據
✅ 說明影響和後續行動
✅ 使用標籤分類: #performance #bug #research

# 避免模糊記錄
❌ "發現了一個問題"
❌ "這個方法比較快"
❌ 沒有程式碼參考
```

### 5. 平台協作溝通

```markdown
# 在 progress.md 中加入協作備註
### 🔗 Collaboration Notes
- Hand-off to Claude Code: 研究完成，請實作 src/core/NewMetric.cpp
- Coordinating with Antigravity: 參數掃描進行中，預計明天完成
- Files recently updated by Antigravity: findings.md, progress.md
```

## 📚 延伸閱讀

- [SKILL.md](./SKILL.md) - Claude Code skill 完整文檔
- [antigravity_agent.md](./antigravity_agent.md) - Antigravity agent 詳細指南
- [InterSubMod CLAUDE.md](../../../.claude/CLAUDE.md) - 專案整體指南

## 🤝 貢獻

如果你發現任何問題或有改進建議，歡迎：

1. 在專案中建立 issue
2. 直接修改並提交 PR
3. 在 AI 對話中提出建議

## 📝 版本歷程

- **v1.0.0** (2026-01-22)
  - 初始版本發布
  - 支援 Claude Code 和 Antigravity
  - 提供 Python bridge script
  - 完整的文檔和範例

## 📄 授權

本專案遵循 MIT 授權條款。

---

**Last Updated**: 2026-01-22
**Maintained by**: InterSubMod Project Team
**Supported Platforms**: Claude Code, Google Antigravity
