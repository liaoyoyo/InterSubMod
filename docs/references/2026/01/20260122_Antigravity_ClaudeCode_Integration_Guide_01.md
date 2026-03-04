<!--
建立時間: 2026-01-22 17:35
目標: 說明 Antigravity 與 Claude Code 的整合方案
處理範圍: 跨平台 skill 移植與配置
關聯檔案:
  - .claude/skills/shared-planning/SKILL.md
  - .claude/skills/shared-planning/antigravity_agent.md
  - .claude/skills/shared-planning/bridge.py
  - .claude/skills/shared-planning/README.md
  - .claude/skills/shared-planning/QUICKSTART.md
-->

# Antigravity 與 Claude Code 整合指南

## 執行摘要

本文檔說明如何讓 **Google Antigravity** 和 **Claude Code** 共同使用同一套任務規劃系統 (`shared-planning`)，實現跨平台協作。

## 問題背景

使用者希望：
1. 同時使用 Antigravity 和 Claude Code 進行開發
2. 兩個平台能共享計劃文件
3. 功能互補使用
4. 在工作內容與技能切換時無縫銜接

### 挑戰

- Claude Code 的 **skills** 是其特有的系統
- Antigravity 有自己的 **Artifacts** 和 **agent** 系統
- planning-with-files skill 無法直接移植到 Antigravity
- 兩個平台需要共享相同的規劃文件格式

## 解決方案

### 核心設計理念

**不移植 skill 本身，而是建立跨平台的規劃系統**

```
planning-with-files (原始 skill)
        ↓
shared-planning (跨平台版本)
        ↓
    ┌───────┴───────┐
    ↓               ↓
Claude Code     Antigravity
(使用 SKILL.md)  (使用 agent prompt)
    ↓               ↓
    └───→ 共享 ←────┘
  docs/shared_planning/
```

### 實作架構

```
.claude/skills/shared-planning/
├── SKILL.md                    # Claude Code skill 定義
├── antigravity_agent.md        # Antigravity custom agent prompt
├── config.json                 # 通用配置
├── bridge.py                   # Python bridge (兩平台通用)
├── README.md                   # 完整文檔
└── QUICKSTART.md               # 快速入門

docs/shared_planning/{TASK_NAME}/
├── task_plan.md               # 任務計劃 (共享)
├── findings.md                # 研究發現 (共享)
└── progress.md                # 進度日誌 (共享)
```

### 關鍵組件

#### 1. SKILL.md - Claude Code Skill

**用途**: Claude Code 透過 skill 系統調用

**功能**:
- 初始化新任務
- 更新進度和發現
- 維護三個核心文件：task_plan.md, findings.md, progress.md

**調用方式**:
```
User: "建立新任務規劃"
Claude Code: [自動調用 /shared-planning skill]
```

#### 2. antigravity_agent.md - Antigravity Agent Prompt

**用途**: Antigravity custom agent 的 system prompt

**功能**:
- 相同的規劃文件格式
- 相同的工作流程
- 平台識別標記

**設定方式**:
1. 在 Antigravity 中建立 Custom Agent
2. 複製 `antigravity_agent.md` 全部內容
3. 貼上作為 system prompt
4. 命名為 "Planning Companion"

#### 3. bridge.py - Python Bridge Script

**用途**: 獨立的 Python 工具，兩個平台都能使用

**功能**:
```bash
# 初始化任務
./bridge.py init <task_name> --type feature

# 更新進度
./bridge.py update <task_name> --progress "message"

# 新增發現
./bridge.py update <task_name> --finding-title "title"

# 列出所有任務
./bridge.py list

# 任務統計
./bridge.py stats
```

**特點**:
- ✅ 自動偵測平台 (Claude Code / Antigravity)
- ✅ 統一的檔案格式
- ✅ 可以被兩個平台的終端/腳本調用

#### 4. 共享文件格式

**task_plan.md**:
```markdown
# Task Plan: {TASK_NAME}

**Created**: 2026-01-22 10:00 (via Claude Code)
**Updated**: 2026-01-22 15:30 (via Antigravity)  ← 平台標記
**Status**: In Progress

## Objective
...

## Tasks
- [x] Completed task
  - Completed: 2026-01-22 (via Claude Code)  ← 平台標記
- [ ] Pending task
  - Assigned to: Antigravity  ← 分配給哪個平台
```

**findings.md**:
```markdown
# Findings: {TASK_NAME}

**Last Updated**: 2026-01-22 (via Antigravity)  ← 平台標記

### 2026-01-22: Finding Title

**Discovery**: 具體發現
**Evidence**: 證據連結
**Implications**: 影響分析
```

**progress.md**:
```markdown
# Progress Log: {TASK_NAME}

**Last Updated**: 2026-01-22 (via Claude Code)  ← 平台標記

## 2026-01-22 Wednesday

**Session Info**:
- Platform: Antigravity  ← 記錄哪個平台
- Focus area: 研究實作方法

### ✅ Completed Today
- 任務描述
  - Files: src/file.cpp
```

## 使用方式

### 情境 1: Claude Code 建立 → Antigravity 研究

```bash
# === Claude Code ===
User: "建立任務：實作 Bhattacharyya distance"
Claude: [調用 /shared-planning skill]
        [建立 docs/shared_planning/bhattacharyya_distance/]

git add docs/shared_planning/
git commit -m "Planning: Initialize bhattacharyya_distance"
git push

# === Antigravity ===
# 使用 Planning Companion agent
User: "Read the bhattacharyya_distance task plan"
Agent: [讀取 task_plan.md]

User: "Research Bhattacharyya distance implementations"
Agent: [使用 Manager View 並行研究]
       [更新 findings.md 記錄發現]
       [更新 progress.md 記錄進度]

git add docs/shared_planning/bhattacharyya_distance/
git commit -m "Planning: Complete research phase"
git push

# === 回到 Claude Code ===
git pull
User: "根據研究結果實作 Bhattacharyya distance"
Claude: [讀取 findings.md]
        [實作 C++ 程式碼]
```

### 情境 2: 純 Python Bridge 工作流程

```bash
# 可以在任何環境（包括腳本）中使用

# Claude Code 終端
./.claude/skills/shared-planning/bridge.py init "task_name"

# Antigravity 終端
./.claude/skills/shared-planning/bridge.py update "task_name" --progress "..."

# CI/CD 腳本
python .claude/skills/shared-planning/bridge.py stats > report.json
```

### 情境 3: 並行工作不同任務

```bash
# Claude Code: 實作功能 A
User: "實作新的距離度量"
Claude: [建立並實作 distance_metrics 任務]

# 同時，Antigravity: 實驗功能 B
User: "Run parameter optimization"
Agent: [建立並執行 parameter_optimization 任務]

# 兩個任務在不同目錄，互不干擾
docs/shared_planning/
├── distance_metrics/        ← Claude Code
└── parameter_optimization/  ← Antigravity
```

## 平台特性比較與分工建議

| 平台 | 優勢 | 建議用途 |
|------|------|----------|
| **Claude Code** | • C++ 程式碼實作<br>• 測試執行<br>• Git 操作<br>• 快速迭代 | • 核心邏輯實作<br>• 單元測試<br>• Bug 修復<br>• 重構 |
| **Antigravity** | • Manager View 並行<br>• Artifacts 視覺化<br>• Browser automation<br>• Gemini 3 Pro | • 文獻調研<br>• 參數掃描<br>• 實驗分析<br>• 多工並行 |

## 技術細節

### 平台偵測機制

```python
# bridge.py 自動偵測
def _detect_platform(self) -> str:
    if os.environ.get('CLAUDE_CODE'):
        return 'Claude Code'
    elif os.environ.get('ANTIGRAVITY_SESSION'):
        return 'Antigravity'
    return 'Unknown'
```

### 檔案更新標記

每次更新都會自動標記平台：

```python
timestamp = datetime.now().strftime('%Y-%m-%d %H:%M')
platform = self._detect_platform()

update_marker = f"**Updated**: {timestamp} (via {platform})"
```

### 衝突避免策略

1. **Git 工作流程**: 切換平台前先 push
2. **時間戳標記**: 每個更新都有時間戳
3. **平台標識**: 清楚標記誰做了什麼變更
4. **獨立目錄**: 不同任務使用不同目錄

## 檔案清單與用途

| 檔案 | 用途 | 使用者 |
|------|------|--------|
| `SKILL.md` | Claude Code skill 定義 | Claude Code |
| `antigravity_agent.md` | Antigravity agent prompt | Antigravity |
| `config.json` | 通用配置 | 兩者共用 |
| `bridge.py` | Python 工具 | 兩者共用 |
| `README.md` | 完整文檔 | 開發者 |
| `QUICKSTART.md` | 快速入門 | 新使用者 |

## 測試驗證

已完成的測試：

```bash
# ✅ 測試 bridge.py 基本功能
./.claude/skills/shared-planning/bridge.py --help

# ✅ 測試任務初始化
./.claude/skills/shared-planning/bridge.py init "test_task"

# ✅ 測試進度更新
./.claude/skills/shared-planning/bridge.py update "test_task" --progress "..."

# ✅ 測試發現記錄
./.claude/skills/shared-planning/bridge.py update "test_task" --finding-title "..."

# ✅ 測試任務列表
./.claude/skills/shared-planning/bridge.py list

# ✅ 測試統計功能
./.claude/skills/shared-planning/bridge.py stats
```

## 後續擴充建議

### 短期 (1-2 週)

- [ ] 在實際專案中使用並收集回饋
- [ ] 優化驗證函數（修復 Objective 偵測問題）
- [ ] 新增更多範本類型
- [ ] 建立範例專案展示

### 中期 (1-2 月)

- [ ] 整合 GitHub Issues 同步
- [ ] 自動產生任務統計圖表
- [ ] Slack/Discord 通知整合
- [ ] Web UI 查看介面

### 長期 (3-6 月)

- [ ] 支援更多 AI 開發平台 (Cursor, Windsurf, etc.)
- [ ] 任務時間追蹤與估算
- [ ] 團隊協作功能
- [ ] CI/CD 整合

## 常見問題

### Q: 為什麼不直接移植 planning-with-files skill？

A: planning-with-files 是 Claude Code 特有的 skill 系統，無法直接在 Antigravity 中使用。我們建立了一個跨平台的解決方案，兩個平台都能使用相同的文件格式。

### Q: Antigravity 一定要用 custom agent 嗎？

A: 不一定。也可以：
1. 直接使用 Python bridge script
2. 手動編輯 Markdown 文件
3. 透過終端執行指令

Custom agent 只是提供了更好的整合體驗。

### Q: 如果兩個平台同時編輯會怎樣？

A:
1. Git 會偵測到衝突
2. 手動合併（通常很簡單）
3. 建議：切換平台前先 `git push`

### Q: 可以用在其他專案嗎？

A: 可以！整個 `.claude/skills/shared-planning/` 目錄可以複製到任何專案。只需要：
```bash
# 複製 skill 目錄
cp -r .claude/skills/shared-planning /path/to/other/project/.claude/skills/

# 建立 shared_planning 目錄
mkdir -p /path/to/other/project/docs/shared_planning
```

## 結論

我們成功建立了一個跨平台的任務規劃系統，讓 Claude Code 和 Antigravity 能夠：

✅ **共享文件**: 相同的 Markdown 格式
✅ **統一流程**: 一致的工作流程
✅ **平台標記**: 清楚追蹤誰做了什麼
✅ **獨立工具**: Python bridge 兩邊都能用
✅ **版本控制**: Git 友善的純文字格式

這個方案不是簡單的 skill 移植，而是一個**真正的跨平台協作系統**。

## 參考資源

### 官方文檔
- [Google Antigravity 官網](https://antigravity.google/)
- [Claude Code 文檔](https://docs.anthropic.com/claude-code)

### 專案文件
- [Shared Planning README](../../../../.claude/skills/shared-planning/README.md)
- [快速入門指南](../../../../.claude/skills/shared-planning/QUICKSTART.md)
- [Claude Code Skill 定義](../../../../.claude/skills/shared-planning/SKILL.md)
- [Antigravity Agent Prompt](../../../../.claude/skills/shared-planning/antigravity_agent.md)

### 外部資源
- [Build with Google Antigravity](https://developers.googleblog.com/build-with-google-antigravity-our-new-agentic-development-platform/)
- [Antigravity Review](https://dev.to/fabianfrankwerner/an-honest-review-of-google-antigravity-4g6f)

---

**建立日期**: 2026-01-22
**作者**: Claude Sonnet 4.5 (via Claude Code)
**專案**: InterSubMod - 跨平台任務規劃系統整合
**版本**: 1.0.0
