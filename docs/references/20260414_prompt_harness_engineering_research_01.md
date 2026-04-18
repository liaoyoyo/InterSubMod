# Prompt Context Harness Engineering 資料彙整報告

<!--
建立時間: 2026-04-14 
目標: 資料收集與分析 — AI coding agent 的 prompt/context/harness engineering 框架、技術與架構
處理範圍: 2025-2026 年間的框架、論文、開源專案、實務案例
關聯檔案:
  - .claude/CLAUDE.md (InterSubMod 現有配置)
  - .claude/settings.local.json (現有 Hooks 配置)
  - AGENTS.md (現有 agent 指引)
-->

## 1. 搜尋概述

- **關鍵詞**: prompt harness engineering, context engineering, CLAUDE.md, AGENTS.md, cursor rules, windsurf rules, aider conventions, multi-agent orchestration, agent memory framework, progressive disclosure, anti-hallucination rules, meta-prompting, codified context
- **搜尋時間**: 2026-04-14
- **資料來源數**: 40+ 來源（論文、技術部落格、GitHub repos、官方文件、業界指南）
- **重要發現**: 2025 年底至 2026 年初，業界用語從 "prompt engineering" 演進為 **"context engineering"** 再到 **"harness engineering"**，反映了從單次提示到系統級架構設計的範式轉移

---

## 2. 分類概述：六大主題

| # | 類別 | 核心趨勢 | InterSubMod 現有程度 |
|---|------|---------|---------------------|
| 1 | Prompt Engineering 框架 | CLAUDE.md 層級 + 漸進揭露 | **已建置** — CLAUDE.md 約 300 行 |
| 2 | Context Management 架構 | 三層分級載入 + 壓縮策略 | **部分建置** — Knowledge 目錄 + MCP |
| 3 | Agent Orchestration 模式 | 子代理分工 + 驗證迴圈 | **已建置** — researcher/skill agents |
| 4 | Harness Engineering | Hook + 感應器 + 導引系統 | **已建置** — 7 類 hooks |
| 5 | Meta-Prompt 模式 | 反幻覺規則 + 自我修正 + 記憶整合 | **部分建置** — 假設陳述規則 |
| 6 | 實務案例 | 生物資訊 + 科研 + 大型 codebase | **參考價值** |

---

## 3. 類別一：Prompt Engineering 框架

### 3.1 CLAUDE.md 最佳實踐

**來源**: [Claude Code Best Practices](https://code.claude.com/docs/en/best-practices) | [shanraisshan/claude-code-best-practice](https://github.com/shanraisshan/claude-code-best-practice)

**核心概念**: CLAUDE.md 是 Claude Code 的指令核心，直接注入系統提示。其內容會被標記為 "OVERRIDE any default behavior"，具有最高指令優先序。

**指令層級體系** (由高到低):
1. 使用者明確指令 (CLAUDE.md + 直接請求)
2. 自定義系統提示追加
3. 預設系統提示
4. 工具定義

**CLAUDE.md 載入位置** (合併順序):
1. `~/.claude/CLAUDE.md` — 個人全域指令
2. `./CLAUDE.md` — 專案根目錄 (git 共享)
3. 子目錄 CLAUDE.md — 針對特定程式碼區域的範圍指令

**關鍵建議**:
- **控制在 150-200 行以內**以確保可靠遵循 (shanraisshan 建議)
- 讓工具 (ESLint, TypeScript, clang-format) 處理格式規則，CLAUDE.md 專注於不可自動驗證的慣例
- 指向專門文件而非在 CLAUDE.md 重複內容
- 每個修改檔案分別 commit，而非打包

**優勢**: 高優先序、自動載入、跨工具相容 (`ln -s CLAUDE.md agents.md`)
**限制**: 過長會降低遵循品質；靜態載入無法條件分支

**InterSubMod 適用性**: 當前 CLAUDE.md 約 300 行，超出建議上限。可考慮將文檔管理規範、確認時機協議等移入 skills 或 `/docs` 目錄，CLAUDE.md 僅保留核心路由邏輯。

---

### 3.2 跨工具配置標準

**來源**: [Anti-Hallucination Rules Gist](https://gist.github.com/mingrath/7e292d9ca976f63e499db971f21b6bbe) | [DeployHQ Config Guide](https://www.deployhq.com/blog/ai-coding-config-files-guide)

| 工具 | 配置檔 | 格式 | 範圍 |
|------|--------|------|------|
| Claude Code | `.claude/CLAUDE.md` + `AGENTS.md` | Markdown | 專案級 |
| Cursor | `.cursor/rules/*.mdc` (原 `.cursorrules`) | MDC | 可條件觸發 (glob) |
| Windsurf | `.windsurfrules` | Markdown | 專案級 |
| Aider | `.aider.conf.yml` + `CONVENTIONS.md` | YAML + Markdown | 專案/全域 |
| GitHub Copilot | `.github/copilot-instructions.md` | Markdown | 專案級 |

**Cursor 的條件觸發機制 (MDC 格式)**: Cursor 的 `.cursor/rules/*.mdc` 支援以 glob 模式篩選適用檔案（如 `*.ts` 只在 TypeScript 檔案被操作時載入），實現按需載入。這比 Claude Code 的 CLAUDE.md 靜態全量載入更精細。

**AGENTS.md 作為通用標準**: 由 OpenAI Codex CLI 普及後，Linux Foundation 管理，已被 60,000+ 開源專案採用。格式簡單且跨工具相容。

**InterSubMod 適用性**: 現已有 AGENTS.md。若未來需跨工具協作，可維護一份核心 AGENTS.md 作為通用基準。

---

### 3.3 Claude Code 系統提示組裝機制

**來源**: [How Claude Code Builds a System Prompt](https://www.dbreunig.com/2026/04/04/how-claude-code-builds-a-system-prompt.html) | [Piebald-AI/claude-code-system-prompts](https://github.com/Piebald-AI/claude-code-system-prompts)

**架構**: 系統提示由 `getSystemPrompt()` 動態組裝，非靜態字串。包含:
- **Always-included**: 身份設定、系統規則、工具使用、風險評估指引
- **Conditional**: 輸出樣式、語言偏好、記憶系統、MCP 伺服器指令
- **Cache Boundary Marker**: 分離全域可快取與會話特定內容

**內建子代理類型**:
| 子代理 | Token 數 | 用途 |
|--------|---------|------|
| Explore | ~494 | 檔案搜尋專家，程式碼導航 |
| Plan | ~636 | 軟體架構師，設計實作計劃 |
| Task | 可配置 | 通用任務執行 |
| Dream (記憶整合) | 可配置 | 記憶清理與鞏固 |

**InterSubMod 適用性**: 了解內部機制有助於設計更精確的子代理指令。InterSubMod 的 researcher agent 指令約為 Dream agent 規模，合理。

---

## 4. 類別二：Context Management 架構

### 4.1 Codified Context 三層架構

**來源**: [Codified Context (arXiv 2602.20478)](https://arxiv.org/abs/2602.20478) — Vasilopoulos, 2026/02

**核心概念**: 在 108,000 行 C# 分散式系統開發中，建立三層 "codified context" 基礎設施來解決 LLM 跨會話失憶問題。

**三層架構**:

| 層級 | 名稱 | 載入時機 | 大小 | 內容 |
|------|------|---------|------|------|
| Tier 1 | 熱記憶 Constitution | 每次會話自動 | ~660 行 Markdown | 慣例、路由表、架構摘要、已知失敗模式 |
| Tier 2 | 領域專家代理 | 按觸發表呼叫 | 19 個 (115-1233 行) | 子系統知識 + 行為指令 |
| Tier 3 | 冷記憶知識庫 | MCP 按需檢索 | 34 份 (~16,250 行) | 完整子系統規格文件 |

**觸發路由表 (Trigger Table)**:
- 修改網路檔案 → 自動路由至 network-protocol-designer
- 修改座標系統 → 觸發 coordinate-wizard
- 自動路由消除「開發者需記憶呼叫哪個代理」的負擔

**量化結果 (283 場 sessions)**:
- 80%+ 人類提示 < 100 字 (預載入上下文減少了解釋需求)
- 知識/程式碼比 = 24.2%
- 57% 的代理呼叫使用專案特定專家
- 元基礎設施維護僅占 4.3% 互動

**關鍵教訓**:
- 規格過時是首要失敗模式 — 實作變更但文件未更新，導致代理生成衝突程式碼
- 為代理寫規格，非為人類 — 包含檔案路徑、函數名稱、明確的 "do/don't" 指令
- 規格視為承重構件 (load-bearing artifacts)，需 Git 漂移檢測

**InterSubMod 適用性**: **高度相關**。InterSubMod 已有:
- Tier 1 對應: CLAUDE.md (~300 行)
- Tier 2 對應: `.claude/agents/` + skills
- Tier 3 對應: `/big8_disk/liaoyoyo2001/Knowledge/` (MCP 服務)

可改進方向: 加入觸發路由表 (改 `.cpp` 觸發 C++ 專家、改 Python 腳本觸發分析專家)。

---

### 4.2 OpenViking 分級載入系統

**來源**: [volcengine/OpenViking](https://github.com/volcengine/OpenViking)

**核心概念**: ByteDance 開源的 AI Agent 上下文資料庫，以檔案系統範式統一管理 Agent 所需的上下文 (記憶/資源/技能)，實現分級投遞與自演化。

**L0/L1/L2 三層結構**:
| 層級 | 內容 | 載入時機 |
|------|------|---------|
| L0 (Abstract) | 一句話摘要 | 快速檢索/索引 |
| L1 (Overview) | 核心資訊 + 使用情境 | 規劃階段決策 |
| L2 (Details) | 完整原始資料 | 深度閱讀必需時 |

**虛擬檔案系統**: `viking://` 協議，頂層目錄映射到上下文域 (resources/user/agent)。支援目錄遞歸檢索 + 語義搜尋組合。

**InterSubMod 適用性**: Knowledge MCP 已實作類似概念 (search → get_doc)。可增加 L0 摘要層加速初次搜尋判斷。

---

### 4.3 漸進揭露 (Progressive Disclosure)

**來源**: [Stop Bloating Your CLAUDE.md](https://alexop.dev/posts/stop-bloating-your-claude-md-progressive-disclosure-ai-coding-tools/) | [Skill Authoring Best Practices](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/best-practices)

**核心概念**: CLAUDE.md 只放通用上下文 (~50 行)；領域知識、陷阱記錄分散到 `/docs` 和 `.claude/agents/`；Skills 的元資料在啟動時載入但 body 按需展開，實現「安裝多個 Skills 不產生上下文代價」。

**實作模式**:

```
CLAUDE.md (50 行) — 專案概述 + 指令路由
  ├── /docs/nuxt-content-gotchas.md — 領域陷阱
  ├── /docs/testing-strategy.md — 測試策略
  ├── .claude/agents/specialist.md — 領域專家
  └── .claude/skills/SKILL.md — 按需技能
```

**`/learn` Skill 模式**: 當 AI 掙扎於已解決問題時，分析對話提取可重用見解，提議存入 `/docs`，經人核准後保存。形成回饋迴圈，文件庫成為「針對 AI 弱點微調」的知識庫。

**SKILL.md 建議**: body < 500 行；超過則拆分為多檔案漸進揭露。

**關鍵洞見**: "AI forgets. Your documentation doesn't." — 無狀態限制反而成為優勢，前提是系統化捕獲學習成果。

**InterSubMod 適用性**: 已有 skills 目錄。可考慮:
- 將 CLAUDE.md 必讀清單改為 skill 觸發 (會話開始時自動讀取)
- 為反覆犯的錯誤建立 "gotchas" 文件 (如已知的 L2 collider bias 陷阱)

---

### 4.4 Anthropic 官方 Context Engineering 指引

**來源**: [Effective Context Engineering for AI Agents](https://www.anthropic.com/engineering/effective-context-engineering-for-ai-agents)

**核心原則**:
1. **有限注意力預算**: Context 是 "precious, finite resource"，token 數增加會產生 n^2 注意力配對，品質衰減 (context rot)
2. **高信號 Token 選取**: 找到「最小化高信號 token 集合」以最大化期望結果
3. **即時載入 (Just-In-Time)**: 維護輕量識別符 (路徑/URL)，動態載入；非預載所有資料
4. **漸進發現**: 透過元資料信號 (檔案大小/命名/時間戳/目錄層次) 引導自主資訊組裝

**壓縮策略層級**:
| 策略 | 侵入程度 | 適用場景 |
|------|---------|---------|
| 工具結果清除 | 最低 | 安全的輕量壓縮 |
| 對話歷史壓縮 | 中等 | 保留架構決策、未解問題 |
| 結構化筆記 (Agentic Memory) | 中等 | 多小時任務跨上下文重置的一致性 |
| 子代理架構 | 最高 | 大規模探索，壓縮摘要回傳 (1-2k token) |

**InterSubMod 適用性**: 已使用子代理 + Knowledge MCP。可增加:
- 明確的 compact 保留指令 (架構決策、活動 bug、檔案範圍、活動約束)
- 在 40-60% context 使用率時主動 compact

---

## 5. 類別三：Agent Orchestration 模式

### 5.1 多代理編排模式分類

**來源**: [Addy Osmani - Code Agent Orchestra](https://addyosmani.com/blog/code-agent-orchestra/) | [Multi-Agent Orchestration Guide](https://www.codebridge.tech/articles/mastering-multi-agent-orchestration-coordination-is-the-new-scale-frontier)

**六大核心模式**:

| 模式 | 描述 | 適用場景 |
|------|------|---------|
| **Supervisor** | 中心節點分配並收集子代理工作 | 一般任務分解 |
| **Sequential Pipeline** | 前一步輸出為下一步輸入 | 線性依賴流程 |
| **Parallel Fan-Out** | 多代理同時執行獨立任務 | 無狀態共享的獨立工作 |
| **Router** | 根據任務類型路由至專家代理 | 多領域混合查詢 |
| **Hierarchical** | 功能負責人再衍生自己的專家 | 深層任務分解 |
| **Evaluator-Optimizer** | 執行 → 評估 → 回饋 → 修訂迴圈 | 品質要求高的產出 |

**Ralph Loop 模式** (原子化工作):
```
pick task → implement → validate → commit if passing → reset context
```
避免 context overflow，透過外部記憶 (git history, task files) 維持連續性。

**Agent Teams 協調原語**:
- 共享任務清單 + 依賴追蹤
- 點對點代理間訊息傳遞
- 檔案鎖定防衝突

**InterSubMod 適用性**: 現為 Supervisor 模式 (主代理 → researcher/task agents)。未來 Phase 2 多樣本分析可考慮 Parallel Fan-Out (每個樣本一個子代理)。

---

### 5.2 長時運行代理的有效控制裝置

**來源**: [Anthropic - Effective Harnesses for Long-Running Agents](https://www.anthropic.com/engineering/effective-harnesses-for-long-running-agents)

**雙階段代理系統**:
1. **Initializer Agent**: 執行一次，建立基礎設施 (Feature Registry JSON、環境腳本、進度文件)
2. **Coding Agent**: 跨多場 session，漸進推進

**Feature Registry 模式**:
- 建立 JSON 列出所有需求功能 (200+)，初始標記為 "failing"
- 防止代理過早宣告完成
- 提供明確實作目標

**Session Startup 儀式**:
```
1. 確認工作目錄
2. 檢閱進度 + git log
3. 讀取功能需求
4. 執行環境啟動腳本
5. 跑 baseline 功能測試
6. 開始單一功能實作
```

**關鍵規則**: 每次 session 只針對一個 feature。防止 "one-shot" 失敗模式 (嘗試過多工作，context 耗盡於實作中途)。

**InterSubMod 適用性**: 已有類似機制 (必讀清單 = Session Startup)。可增加 Feature Registry 概念用於追蹤 Phase 2 多步驟實作。

---

### 5.3 ACE-FCA 三階段工作流程

**來源**: [Advanced Context Engineering for Coding Agents (HumanLayer)](https://github.com/humanlayer/advanced-context-engineering-for-coding-agents/blob/main/ace-fca.md)

**三階段流程**:

```
Research (探索) → Planning (規劃) → Implementation (實作)
     ↑ 子代理      ↑ 人工審查       ↑ 分階段驗證
```

**Context 最佳化層次** (優先序):
1. **Correctness** — 防止錯誤資訊進入 context
2. **Completeness** — 包含必要資訊
3. **Size** — 最小化雜訊和 token 用量
4. **Trajectory** — 維持問題的邏輯流向

**高槓桿人工審查點**:
- **Research 驗證** — 確認 codebase 理解正確再繼續 (錯誤研究 = 數千行壞程式碼)
- **Plan 審查** — 確保方法與架構對齊再實作 (錯誤計劃 = 數百行壞程式碼)
- **Phase 驗證** — 驗證每階段按預期運作 (錯誤程式碼 = 錯誤程式碼)

**結果驗證**: 非 Rust 專家使用此框架在 300k LOC Rust codebase 貢獻可合併 PR；7 小時內完成 35k LOC 功能新增。

**關鍵洞見**: "Specs are the real code" — 研究文件和計劃是永久構件，實作是次要驗證。

**InterSubMod 適用性**: **高度對齊**。現有的必讀清單 → 實驗設計 → 實作驗證已是類似流程。可將 Step→Verify 格式更明確對齊此三階段模型。

---

### 5.4 驗證瓶頸與品質閘門

**來源**: [Addy Osmani](https://addyosmani.com/blog/code-agent-orchestra/) | [VMAO Framework (arXiv 2603.11445)](https://arxiv.org/html/2603.11445v1)

**核心論述**: "The bottleneck is no longer generation. It's verification."

**品質閘門實作**:
1. **AGENTS.md**: 人工策劃的制度記憶，每次 session 更新發現的模式
2. **Token 預算管理**: 每代理硬限制 + 85% 消耗時自動暫停
3. **Kill Criteria**: 代理卡在同一錯誤 3+ 迭代時重新分配
4. **Plan Approval**: 代理在實作前提案方法，及早捕捉架構問題

**VMAO 框架** (Verified Multi-Agent Orchestration):
```
1. 將查詢分解為 DAG 子問題
2. 分配至領域專家代理
3. 依賴感知排程下並行執行
4. LLM 評估驗證集體完整性
5. 自適應重規劃填補缺口
```

**InterSubMod 適用性**: 已有確認時機協議 (Hard Gate / Gate / Review / FYI)，與 Plan Approval 對齊。可增加 Token 預算管理和 Kill Criteria。

---

## 6. 類別四：Harness Engineering

### 6.1 Harness Engineering 定義與框架

**來源**: [Martin Fowler - Harness Engineering](https://martinfowler.com/articles/harness-engineering.html) | [OpenAI - Harness Engineering](https://openai.com/index/harness-engineering/)

**核心公式**: `Agent = Model + Harness`

Harness 是 AI 代理中除模型以外的一切: 系統提示、程式碼檢索、工具排程、安全執行、上下文管理。

**雙向控制系統**:

| 控制類型 | 名稱 | 時機 | 作用 |
|---------|------|------|------|
| Feedforward | **Guides** (導引) | 代理行動前 | 預期行為並引導，提高首次嘗試成功率 |
| Feedback | **Sensors** (感應器) | 代理行動後 | 觀察結果，啟用自我修正 |

**關鍵洞見**: "單獨使用任一方，要麼得到一個不斷重複錯誤的代理 (feedback-only)，要麼得到一個編碼規則但不驗證效果的代理 (feedforward-only)"。

**執行分類**:
| 類型 | 特性 | 範例 |
|------|------|------|
| **Computational** | 確定性、快速、CPU | 測試、linter、型別檢查 |
| **Inferential** | 語義分析、GPU/NPU、非確定性 | 程式碼審查代理、LLM-as-judge |

---

### 6.2 三大規範範疇

**來源**: [Martin Fowler - Harness Engineering](https://martinfowler.com/articles/harness-engineering.html)

| 範疇 | 導引 (Feedforward) | 感應器 (Feedback) | 覆蓋缺口 |
|------|-------------------|------------------|---------|
| **Maintainability** | 程式碼風格規範、架構規則 | Linter、duplication check、coverage | 誤診、不必要功能、誤解指令 |
| **Architecture Fitness** | 效能需求 skill、可觀察性標準 | 效能測試、日誌慣例驗證 | 跨模組交互效果 |
| **Behaviour** | 功能規格文件 | AI 生成測試套件 | **最弱** — 測試品質依賴 AI 判斷 |

**生命週期整合 ("Keep Quality Left")**:
- **Pre-commit**: 快速 computational controls (linter, 基本 code review)
- **Post-integration pipeline**: 昂貴 sensors (mutation testing, 全面 review)
- **Continuous monitoring**: 漂移檢測 (dead code, 依賴掃描, SLO 退化)

**InterSubMod 適用性**: 已有:
- Feedforward: CLAUDE.md 規範、假設陳述規則、Step→Verify 格式
- Feedback: PreToolUse 阻擋 (未編譯不許 commit)、PostToolUse 標記
- 可增加: 架構適配感應器 (如 benchmark 回歸檢測)

---

### 6.3 Claude Code Hooks 系統

**來源**: [Claude Code Hooks Guide](https://code.claude.com/docs/en/hooks-guide) | [Building Guardrails for Claude Code](https://dev.to/mikelane/building-guardrails-for-ai-coding-assistants-a-pretooluse-hook-system-for-claude-code-ilj) | [Skills and Hooks Starter Kit](https://medium.com/@davidroliver/skills-and-hooks-starter-kit-for-claude-code-c867af2ace32)

**Hook 類型與觸發點**:

| Hook | 觸發時機 | 可做動作 |
|------|---------|---------|
| `PreToolUse` | 工具 (Bash/Write/Edit) 執行前 | 阻擋 (exit 2)、修改輸入、通知 |
| `PostToolUse` | 工具完成後 | 設定標記、觸發後處理、通知 |
| `UserPromptSubmit` | 用戶提交提示時 | 提醒、追加上下文、警告 |
| `SubagentStop` | 子代理完成時 | 檢查結果、清理 |
| `Stop` | 會話結束時 | 報告撰寫提醒 |
| `Notification` | Claude 等待中 | 外部通知 (Telegram, Slack) |

**優先序**: `policySettings (enterprise) > projectSettings > userSettings > localSettings > plugin hooks`；多層 hooks 合併全部執行。

**安全**: 工作空間需 trusted 才執行 hooks；未信任工作空間靜默跳過。

**阻擋模式**: `exit 2` 或回傳 `{"decision": "block"}` 可阻止工具執行。

**InterSubMod 現狀**: 已有 7 類 hooks (見 CLAUDE.md)。這是業界相對完整的配置。

---

### 6.4 Skills 與 Commands 統一

**來源**: [Claude Code Setup After 4 Months](https://okhlopkov.com/claude-code-setup-mcp-hooks-skills-2026/) | [Agent Skills Overview](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/overview)

**2026 現狀**: `.claude/commands/` 仍可用，但建議遷移至 `.claude/skills/`。每個 skill 自動獲得 `/slash-command` 介面。

**Skill 元資料 (YAML frontmatter)**:
```yaml
---
name: skill-name
description: One-line description for discovery
argument_hints:
  - hint1
  - hint2
allowed_tools:
  - Bash
  - Read
---
```

**關鍵設計**: Skill 元資料在啟動時載入系統提示 (僅知道存在)，body 在觸發時才展開。安裝多個 Skills 不產生 context 代價。

**InterSubMod 適用性**: 已有 3 個 skills。可考慮為常見工作流程建立更多 skills (如 benchmark 執行、實驗記錄、文件歸檔)。

---

## 7. 類別五：Meta-Prompt 模式

### 7.1 反幻覺規則

**來源**: [Anti-Hallucination Rules Gist](https://gist.github.com/mingrath/7e292d9ca976f63e499db971f21b6bbe) | [Claude Reduce Hallucinations Guide](https://platform.claude.com/docs/en/test-and-evaluate/strengthen-guardrails/reduce-hallucinations)

**五大核心規則**:

| # | 規則 | 說明 |
|---|------|------|
| 1 | **承認不確定** | 不確定/缺乏資訊時必須說 "I don't know"，而非編造 |
| 2 | **驗證後才聲稱** | 使用工具檢查實際狀態，而非依賴記憶。"Your memory of how code works is often wrong — the file is always right." |
| 3 | **防止假設堆疊** | 不在未驗證基礎上疊加答案。一個不確定假設疊另一個 = 自信但錯誤 |
| 4 | **立即自我修正** | 發現不確定/錯誤時停下承認，而非為流暢完成句子 |
| 5 | **引用所有主張來源** | 標註資訊來自哪個檔案/行號/工具輸出。"If you can't point to where you learned it, you're probably making it up." |

**實作要點**:
- 每條 "do not" 規則對應一個**曾經發生過的幻覺**。規則檔是模型過去錯誤的活記錄
- 模型對明確否定指令反應良好 (如 "do not use getServerSession()")
- 指定技術棧版本號是最有效的單一反幻覺手段

**InterSubMod 適用性**: 已有「假設陳述規則」和「不要憑記憶回答可以查證的事實」。可考慮:
- 增加「已知幻覺記錄」文件，記錄 AI 曾犯的具體錯誤 (如 L2 collider bias 事件)
- 在 CLAUDE.md 中加入版本號 (如 HTSlib 版本、Eigen3 版本)

---

### 7.2 Auto Dream 記憶整合

**來源**: [Claude Code Dreams](https://claudefa.st/blog/guide/mechanics/auto-dream) | [dream-skill (GitHub)](https://github.com/grandamenium/dream-skill)

**核心概念**: 如同人腦 REM 睡眠的記憶整合。Auto Dream 在觸發條件滿足 (24h + 5 sessions) 時自動執行四階段清理:

1. **時間標準化**: 相對日期轉絕對日期 ("Yesterday we decided..." → "On 2026-03-15 we decided...")
2. **矛盾刪除**: 移除被新事實取代的舊記憶 (Express → Fastify)
3. **過時清理**: 刪除指向已移除檔案的除錯筆記
4. **組織整理**: 結構化分類

**四大記憶機制對照**:
| 機制 | 控制方 | 內容 |
|------|--------|------|
| CLAUDE.md | 人工 | 權威規則 |
| Auto Memory | Claude | 工作中學到的知識 |
| Session Memory | Claude | 對話連續性 |
| Auto Dream | Claude (自動) | 保持知識乾淨、最新、有組織 |

**InterSubMod 適用性**: 已有 MEMORY.md 手動管理。可考慮:
- 啟用 Auto Dream 功能 (如已可用)
- 或在 Stop hook 中加入記憶清理提醒
- 定期檢閱 MEMORY.md 中的 "Concluded" 區塊是否需清理

---

### 7.3 Meta-Prompting 框架

**來源**: [Meta Prompting for AI Systems (arXiv 2311.11482)](https://arxiv.org/abs/2311.11482) | [Prompting Guide - Meta Prompting](https://www.promptingguide.ai/techniques/meta-prompting)

**核心概念**: Meta-prompting 專注於任務的**形式結構**而非具體內容，讓 LLM 透過結構性推理提升能力。

**應用於自我感知**:
- 迭代精煉提示 (附回饋)，調整查詢以更好反映真正目標
- 加入偏見規避指令和事實核查要求
- **警告**: 無外部回饋/評估指標時，信任模型自評估可能導致「自信但次優」的提示

**2026 演進**: "Context engineering 已取代 prompt engineering 成為關鍵技術學科"。Meta-prompting 需整合進 harness 架構，而非獨立存在。

**InterSubMod 適用性**: 已有的「假設陳述規則」和「目標驅動驗證格式」本質上是 meta-prompt pattern。可考慮更明確的結構化推理模板。

---

### 7.4 防止重複錯誤的機制

**來源**: 綜合多來源

**已驗證的防錯機制**:

| 機制 | 實作方式 | 效果 |
|------|---------|------|
| **結論索引** | `docs/experiments/INDEX.md` 記錄成功/失敗 | 防止重探已失敗方向 |
| **已關閉研究警告 Hook** | UserPromptSubmit hook 警告 | 自動提醒不要重走老路 |
| **Memory "Concluded" 區塊** | MEMORY.md 分離活躍/結論 | 快速辨識已結案項目 |
| **Anti-hallucination 規則** | "do not" 明確否定清單 | 模型可靠避開指定行為 |
| **Spec drift detection** | Git commit 監控 + 規格更新 | 防止規格過時導致衝突 |

**InterSubMod 現狀**: 前四項**已完整實作**，這在業界屬於領先配置。第五項 (規格漂移檢測) 可通過 PostToolUse hook 自動化。

---

## 8. 類別六：實務案例

### 8.1 Codified Context: 108k LOC C# 分散式系統

**來源**: [arXiv 2602.20478](https://arxiv.org/abs/2602.20478)

| 指標 | 數值 |
|------|------|
| 程式碼規模 | 108,000 行 C# |
| 開發期 | 70 天 (兼職) |
| Session 數 | 283 |
| 人類提示 | 2,801 |
| 專家代理數 | 19 |
| 知識庫文件 | 34 份 |
| 知識/程式碼比 | 24.2% |

**與 InterSubMod 對比**: InterSubMod 的 C++ 核心約 15-20k LOC + Python 腳本，知識庫 (Knowledge/) 內容豐富。知識/程式碼比可能更高，符合科研專案特性。

---

### 8.2 ACE-FCA: 300k LOC Rust Codebase

**來源**: [HumanLayer ACE-FCA](https://github.com/humanlayer/advanced-context-engineering-for-coding-agents)

- 非 Rust 專家使用 ACE-FCA 框架在 300k LOC Rust codebase 貢獻可合併 PR
- 7 小時內完成 35k LOC 功能新增 (cancellation + WASM)
- 關鍵: 頻繁 context reset + 結構化構件 (research → plan → implement)

**與 InterSubMod 對比**: InterSubMod 的 Research → Plan → Implement 流程已隱含此模式。可更明確地將研究文件視為 "primary review target"。

---

### 8.3 生物資訊 AI Agent 系統

**來源**: [PromptBio (bioRxiv)](https://www.biorxiv.org/content/10.1101/2025.07.05.663295v1) | [BioAgents (Nature Scientific Reports)](https://www.nature.com/articles/s41598-025-25919-z) | [From Prompt to Agent Engineering (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12613637/)

**PromptBio**: 多代理生物資訊平台
- DataAgent / OmicsAgent / AnalysisAgent / QAgent 協作
- 逐步人機協作工作流程
- 使用預驗證的領域標準工具

**BioAgents**: 基於小型語言模型 + RAG 的多代理系統
- 在概念基因組學任務上達到人類專家水準
- 支援本地運行和私有資料個人化

**From Prompt to Agent Engineering**: 
- Agent Engineering 定義為「系統化設計、實作和評估 AI 代理以解決生物醫學研究中動態多學科挑戰」
- 代理需具備: 設定目標、分解任務、檢索資訊、與其他代理協作、適應行為

**InterSubMod 適用性**: 這些系統側重於通用生物資訊分析管線。InterSubMod 更接近 Codified Context 案例 — 單一專精工具的深度開發。但 PromptBio 的多代理協作模式可啟發 Phase 2 多樣本分析流程。

---

### 8.4 Claude Code 日常使用 4 個月經驗

**來源**: [Claude Code Setup After 4 Months](https://okhlopkov.com/claude-code-setup-mcp-hooks-skills-2026/)

**實踐者配置**:
- 3 個 MCP Server (Coolify 部署、Telegram 訊息、Codex 雙模型審查)
- Pre-commit hooks 阻擋敏感資料 (.env, .key, .pem)
- Skills: 區塊鏈分析 SQL (~200行)、錢包調查、多語言發布
- 24/7 無人值守運作 + 安全護欄

**生產力影響**: 30 頁區塊鏈報告 + 15 張圖表 = 一個晚上 (vs 一週手動工作)

**與 InterSubMod 對比**: InterSubMod 已有更複雜的 hook 系統和更豐富的 skill 配置，在架構成熟度上已超越此案例。

---

## 9. 衝突觀點與未解決議題

### 9.1 衝突觀點

| 觀點 A | 觀點 B | 可能原因 |
|--------|--------|----------|
| CLAUDE.md 應 < 150 行 (shanraisshan) | Codified Context constitution ~660 行 (arXiv) | 專案複雜度差異；660 行含路由表而非全部是指令 |
| Skills 自動觸發可靠 (官方) | Skills 56% 從未被觸發 (Vercel 評估) | 條件匹配的精確度因場景而異 |
| 全自動模式可行 (多來源) | "Active human engagement" 不可省 (ACE-FCA) | 任務類型差異：常規 vs 創新研究 |
| 記憶框架值得投資 (Mem0/Letta) | 檔案系統 + Git 夠用 (Codified Context) | 專案規模和團隊大小差異 |

### 9.2 未解決議題

1. **行為控制裝置 (Behaviour Harness) 仍最弱** — AI 生成的測試品質尚不足以取代人工驗證
2. **規格漂移檢測** — 自動化程度不足，目前依賴人工紀律
3. **Inferential sensors 的非確定性** — LLM-as-judge 的可靠性仍待提升
4. **Context 壓縮的資訊損失** — 過度壓縮會丟失微妙但關鍵的資訊，無自動品質評估
5. **跨工具一致性** — 不同 AI coding 工具的配置格式碎片化

---

## 10. 資料來源評估

| 來源 | 類型 | 可信度 | 備註 |
|------|------|--------|------|
| [Anthropic Engineering Blog](https://www.anthropic.com/engineering/) | 官方文件 | 高 | 一手來源，但有商業利益 |
| [arXiv 2602.20478 (Codified Context)](https://arxiv.org/abs/2602.20478) | 學術預印本 | 高 | 283 session 量化驗證，方法論嚴謹 |
| [arXiv 2603.05344 (OPENDEV)](https://arxiv.org/abs/2603.05344) | 學術預印本 | 中高 | 工作中論文，持續更新 |
| [Martin Fowler](https://martinfowler.com/articles/harness-engineering.html) | 業界權威部落格 | 高 | 系統性框架分析，含 OpenAI/Stripe 案例 |
| [Addy Osmani](https://addyosmani.com/blog/code-agent-orchestra/) | 業界權威部落格 | 高 | Google 工程師，實戰經驗豐富 |
| [HumanLayer ACE-FCA](https://github.com/humanlayer/advanced-context-engineering-for-coding-agents) | 開源指南 | 中高 | 有量化結果 (35k LOC)，但單一團隊經驗 |
| [shanraisshan/claude-code-best-practice](https://github.com/shanraisshan/claude-code-best-practice) | 社群最佳實踐 | 中 | 實用但缺乏大規模驗證 |
| [Piebald-AI/claude-code-system-prompts](https://github.com/Piebald-AI/claude-code-system-prompts) | 逆向工程 | 中高 | 精確但隨版本變動 |
| [okhlopkov.com](https://okhlopkov.com/claude-code-setup-mcp-hooks-skills-2026/) | 個人經驗分享 | 中 | 單一使用者經驗，但實戰導向 |
| [alexop.dev](https://alexop.dev/posts/stop-bloating-your-claude-md-progressive-disclosure-ai-coding-tools/) | 技術部落格 | 中 | 實用建議，含 Vercel 評估數據 |

---

## 11. InterSubMod 現有配置評估

### 11.1 已具備的強項

| 面向 | InterSubMod 現狀 | 業界對標 |
|------|-----------------|---------|
| Instruction hierarchy | CLAUDE.md + AGENTS.md + agents/ + skills/ | 完整四層 |
| Hook system | 7 類 hooks (阻擋/標記/提醒) | 業界領先 |
| Knowledge tiering | Knowledge MCP (search/get_doc) + docs/ | 對齊 Codified Context Tier 1-3 |
| Research state management | INDEX.md + MEMORY.md (Active/Concluded) | 獨特優勢 |
| Anti-repetition | 已關閉研究警告 hook + Concluded 區塊 | 少見的完整實作 |
| Decision gates | Hard Gate / Gate / Review / FYI 四級 | 超越多數開源實踐 |
| Session startup | 必讀清單 (強制 5 檔案依序讀取) | 對齊 Anthropic 長時運行代理建議 |

### 11.2 可改進的方向

| 改進項目 | 優先度 | 預期效益 | 對應來源 |
|---------|--------|---------|---------|
| CLAUDE.md 瘦身 (300→150 行) | 中 | 提升遵循品質 | shanraisshan, alexop.dev |
| 觸發路由表 (檔案類型→專家代理) | 中 | 減少人工指派代理的負擔 | Codified Context |
| `/learn` Skill (捕獲已知陷阱) | 中 | 系統化累積 AI 弱點知識 | alexop.dev |
| Compact 保留指令 | 低 | 改善長對話品質 | Anthropic, ACE-FCA |
| Feature Registry (Phase 2) | 低 | 追蹤多步驟實作進度 | Anthropic long-running agents |
| Spec drift detection hook | 低 | 防止知識庫過時 | Codified Context |
| 已知幻覺記錄文件 | 低 | 防止重複犯錯 | Anti-hallucination gist |

---

## 12. 建議行動

### 短期 (可立即實施)

1. **建立 "gotchas" 文件** (`docs/references/known_pitfalls.md`) — 記錄 AI 曾犯的具體錯誤 (L2 collider bias, pooled OLS 陷阱等)，CLAUDE.md 中加入「新分析前先讀此檔」指令
2. **在 CLAUDE.md 加入工具版本號** — HTSlib, Eigen3, jemalloc 等精確版本，減少幻覺

### 中期 (Phase 2 開發期間)

3. **CLAUDE.md 瘦身重構** — 將文檔管理規範、確認時機協議移入 skills，CLAUDE.md 保留核心路由 (<150 行)
4. **觸發路由表** — 在 CLAUDE.md 或 Hook 中建立檔案類型→專家代理映射 (`.cpp/.hpp` → C++ specialist, `.py` → analysis specialist)
5. **Feature Registry** — 為 Phase 2 Normal Methylation Reference 建立 JSON 功能追蹤表

### 長期 (持續改進)

6. **`/learn` Skill** — 自動從對話中提取可重用見解，經人核准存入 `/docs`
7. **Spec drift detection** — PostToolUse hook 監測 Knowledge/ 文件最後更新時間，超過 N 天未更新的相關規格提出警告

---

## 附錄 A：關鍵 URL 索引

### 官方文件
- [Claude Code Best Practices](https://code.claude.com/docs/en/best-practices)
- [Claude Code Hooks Guide](https://code.claude.com/docs/en/hooks-guide)
- [Create Custom Subagents](https://code.claude.com/docs/en/sub-agents)
- [Skill Authoring Best Practices](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/best-practices)
- [Compaction](https://platform.claude.com/docs/en/build-with-claude/compaction)
- [Reduce Hallucinations](https://platform.claude.com/docs/en/test-and-evaluate/strengthen-guardrails/reduce-hallucinations)

### 業界框架與分析
- [Anthropic - Effective Context Engineering](https://www.anthropic.com/engineering/effective-context-engineering-for-ai-agents)
- [Anthropic - Effective Harnesses for Long-Running Agents](https://www.anthropic.com/engineering/effective-harnesses-for-long-running-agents)
- [Martin Fowler - Harness Engineering](https://martinfowler.com/articles/harness-engineering.html)
- [Addy Osmani - Code Agent Orchestra](https://addyosmani.com/blog/code-agent-orchestra/)
- [Stack Builders - Beyond AGENTS.md](https://www.stackbuilders.com/insights/beyond-agentsmd-turning-ai-pair-programming-into-workflows/)

### 學術論文
- [Codified Context (arXiv 2602.20478)](https://arxiv.org/abs/2602.20478)
- [OPENDEV (arXiv 2603.05344)](https://arxiv.org/abs/2603.05344)
- [VMAO Framework (arXiv 2603.11445)](https://arxiv.org/html/2603.11445v1)
- [Meta Prompting (arXiv 2311.11482)](https://arxiv.org/abs/2311.11482)
- [Context Engineering: From Prompts to Architecture (arXiv 2603.09619)](https://arxiv.org/pdf/2603.09619)

### 開源專案與工具
- [HumanLayer ACE-FCA](https://github.com/humanlayer/advanced-context-engineering-for-coding-agents)
- [Piebald-AI/claude-code-system-prompts](https://github.com/Piebald-AI/claude-code-system-prompts)
- [shanraisshan/claude-code-best-practice](https://github.com/shanraisshan/claude-code-best-practice)
- [volcengine/OpenViking](https://github.com/volcengine/OpenViking)
- [grandamenium/dream-skill](https://github.com/grandamenium/dream-skill)
- [Meirtz/Awesome-Context-Engineering](https://github.com/Meirtz/Awesome-Context-Engineering)
- [Aider](https://aider.chat/)

### 實務經驗
- [Claude Code Setup After 4 Months](https://okhlopkov.com/claude-code-setup-mcp-hooks-skills-2026/)
- [Stop Bloating Your CLAUDE.md](https://alexop.dev/posts/stop-bloating-your-claude-md-progressive-disclosure-ai-coding-tools/)
- [Anti-Hallucination Rules Gist](https://gist.github.com/mingrath/7e292d9ca976f63e499db971f21b6bbe)
- [How Claude Code Builds a System Prompt](https://www.dbreunig.com/2026/04/04/how-claude-code-builds-a-system-prompt.html)

### 生物資訊 AI Agent
- [PromptBio (bioRxiv)](https://www.biorxiv.org/content/10.1101/2025.07.05.663295v1)
- [BioAgents (Nature Scientific Reports)](https://www.nature.com/articles/s41598-025-25919-z)
- [From Prompt to Agent Engineering (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12613637/)

### Agent Memory 框架
- [mem0.ai - State of AI Agent Memory 2026](https://mem0.ai/blog/state-of-ai-agent-memory-2026)
- [OSS Insight - Agent Memory Race 2026](https://ossinsight.io/blog/agent-memory-race-2026)
- [The New Stack - Memory for AI Agents](https://thenewstack.io/memory-for-ai-agents-a-new-paradigm-of-context-engineering/)
