# Production Harness Engineering Patterns for Long-Running AI Coding Agents (2025-2026)

<!--
建立時間: 2026-04-14 21:00
目標: 收集生產級 AI 代理的 harness engineering 進階模式
處理範圍: guardrails、self-healing、observability、config-as-code、prompt caching、MCP patterns
-->

## 1. Advanced Guardrail Engineering

### 1.1 Three-Tier Handler Architecture (Claude Code Hooks, 2026)

Claude Code hooks 已演化至 **21 lifecycle events + 4 handler types**，形成語義感知的分級防護：

| Handler Type | 機制 | 適用場景 | 延遲成本 |
|---|---|---|---|
| **Command** | Shell script, exit code 控制 | Lint、格式化、路徑黑名單 | < 1s |
| **HTTP** | JSON POST 到遠端端點 | 集中式政策服務、audit logging | 1-5s |
| **Prompt** | 單輪 LLM 語義評估 | 語義安全分類（碰到 auth/payment/DB 則 DENY） | 2-5s |
| **Agent** | 子代理 + 工具存取（Read/Grep/Glob） | 跨檔案一致性驗證、測試覆蓋檢查 | 5-30s |

**關鍵設計原則**：

- **Exit code 2 = 強制阻擋**：Unix 慣例的 exit 1 不會阻擋，必須用 exit 2
- **`if` 條件過濾**：PreToolUse 支援 `"if": "Bash(rm *)"` 語法精確匹配危險命令
- **updatedInput**：hooks 可修改工具輸入參數再放行（不只是 allow/deny）
- **asyncRewake**：背景 hook 完成後喚醒 Claude 繼續（長時間建置/測試）

### 1.2 Graduated Response Model

來自 Galileo 和多個生產系統的共識模式：

```
Level 0: 靜默記錄（FYI）
Level 1: 警告 + 繼續（stderr 訊息）
Level 2: 阻擋 + 要求替代方案（exit 2 + reason）
Level 3: 終止會話（continue: false + stopReason）
Level 4: 人工介入升級（circuit breaker open）
```

### 1.3 AGENTS.md 標準化（Linux Foundation AAIF）

2025 年 8 月 OpenAI 發布，12 月捐贈至 Linux Foundation Agentic AI Foundation。已被 60,000+ 開源專案採用，支援 Cursor、Copilot、Gemini CLI、Devin 等。核心理念：專案配置即 Markdown，目錄樹中最近的檔案優先。

## 2. Self-Healing & Error Recovery

### 2.1 Replit Snapshot Engine — Copy-on-Write 沙箱

Replit 的生產實作是目前最成熟的 AI agent checkpoint/rollback 系統：

- **底層機制**：16 MiB chunk 不可變儲存於 GCS，manifest 指向 chunk 組合
- **複製成本**：複製整個環境 = 複製 manifest = O(1) 常數時間
- **資料庫分支**：PostgreSQL 跑在同一儲存基礎設施上，checkpoint 同時包含 DB 狀態
- **Parallel Sampling**：利用 LLM 非決定性，fork 多個環境讓不同 agent 嘗試同一問題，選最佳方案原子性合併
- **安全實驗**：在隔離環境中放寬 guardrails，允許 agent 執行破壞性測試

### 2.2 SagaLLM — 分散式交易模式套用至多代理系統

VLDB 2025 論文，將 Saga Pattern 應用於 LLM 多代理工作流：

- 每個工作流節點有 **regular agent**（正向執行）+ **compensation agent**（回滾執行）
- 放寬 ACID 為 workflow-wide consistency：模組化 checkpoint + 可補償執行
- LLM 自動生成：狀態追蹤、依賴分析、log schema、recovery orchestration
- 實際效果：一致性、驗證精確度、不確定環境下的自適應協調顯著提升

### 2.3 五層生產安全模式

綜合 Galileo、Temporal、多來源的共識：

1. **Circuit Breaker for Quality**：監控錯誤率/延遲/行為異常，三態切換（Closed/Open/Half-Open）
2. **Validation Gate before Tool Execution**：PreToolUse 驗證層
3. **Idempotent Workflow + Saga Rollback**：每步可重試或回滾
4. **Token & Cycle Budget Guardrail**：硬性上限防止失控
5. **Human Escalation for High-Risk**：明確升級路徑

## 3. Observability for AI Agents

### 3.1 Three Pillars Adapted for Agents

| 傳統支柱 | Agent 對應 | 關鍵差異 |
|---|---|---|
| **Traces** | Decision path reconstruction | 每個 LLM call + tool invocation + retrieval + 中間決策都帶完整 context |
| **Metrics** | Token consumption, latency, cost attribution | Per-agent per-model 成本歸因到 token 級別 |
| **Logs** | Durable audit record | 每個 agent action 記錄 context + policies + data assets |

### 3.2 Token Budget & Cost Tracking

生產級成本監控模式（來自 Vantage、Langfuse、Braintrust）：

- **警報閾值**：「平均 token/request > 上週 1.5 倍」自動觸發
- **5% 請求消耗 50% token** 的 Pareto 分析
- **Decision-tree 視覺化**：標記昂貴的工作流分支
- **Per-feature cost attribution**：哪個功能/agent 花最多

### 3.3 Observability Gap

McKinsey 2025: 88% 組織使用 AI，62% 實驗 AI agents，但僅 < 10% 成功規模化。PwC: 79% 採用 AI agents 但無法追蹤多步驟工作流中的失敗。差距核心 = observability 不足。

### 3.4 工具生態

| 工具 | 定位 | 特色 |
|---|---|---|
| **Langfuse** | 開源 LLM 工程平台 | Span-level cost、multi-step agent workflow |
| **Braintrust** | Evaluation-first | 自動評分、即時監控、production feedback loops |
| **Galileo** | Agent 可靠性 | Free tier、failure recovery 可視化 |

## 4. Configuration as Code

### 4.1 Agent Config 版本化原則

- **Prompts = versioned assets**：與程式碼同等對待，存入 Git
- **Semantic version tags**：連接實驗、code commit、data snapshot、hyperparameter
- **環境隔離**：dev/staging/prod 各有獨立 agent 配置

### 4.2 CI/CD for Evals

AI agent 的測試本質問題：輸出是 probabilistic（Y-ish），不是 deterministic（Y）。

| 評估策略 | 適用場景 | 實作方式 |
|---|---|---|
| **Code-based assert** | 有客觀正確答案 | pytest + exact match |
| **LLM-as-Judge** | 主觀品質評估 | Secondary model scoring |
| **Regression suite** | PR 品質門檻 | CI 自動跑 eval，低於閾值則 block merge |
| **Span-level agent eval** | Tool selection / planning / reasoning | Confident AI + pytest 整合 |

**注意**：LLM-as-Judge 用於客觀任務時不可靠，應保留給真正主觀的評估。

### 4.3 Hook 配置分層

Claude Code 的 hooks 配置有明確的作用域分層：

| 位置 | 作用域 | 可分享 |
|---|---|---|
| `~/.claude/settings.json` | 全域（所有專案） | 否 |
| `.claude/settings.json` | 專案級（所有開發者） | 是（Git） |
| `.claude/settings.local.json` | 專案級（個人） | 否（gitignore） |
| Plugin `hooks/hooks.json` | 啟用時 | 是（plugin 打包） |

## 5. Prompt Caching Strategies (Claude)

### 5.1 Cache-Aware Architecture（Claude Code 核心設計）

Anthropic 工程師 Thariq 的核心洞見：「你必須從 prompt caching 角度出發設計 agent，幾乎每個功能都會碰到它。」

**Cache prefix 建立順序**：`tools → system → messages`

**黃金規則**：靜態內容在前，動態內容在後。

| 層級 | 內容 | 共享範圍 |
|---|---|---|
| System prompt | 所有 Claude Code 使用者共享 | 全球 |
| CLAUDE.md | 同專案使用者共享 | 專案 |
| Tool definitions | 同配置使用者共享 | 配置 |
| Conversation history | 單一會話 | 會話 |

### 5.2 三個反直覺設計決策

1. **Deferred Loading（延遲載入）**：MCP 工具不全部展開 schema，改用輕量 stub（僅名稱 + `defer_loading: true`），需要時透過 ToolSearch 載入。因為新增一個 tool 會改變 prefix → 整個對話 cache 失效。

2. **Plan Mode 不換工具集**：直覺設計是切換 read-only tools，但這會破壞 cache。實際做法是保留全部工具 + 新增 EnterPlanMode/ExitPlanMode 兩個工具 + 以 user message 傳遞模式指令。

3. **Compaction 保留相同 prefix**：壓縮對話時，system prompt + tool definitions + CLAUDE.md 完全不變，只替換舊訊息為摘要。prefix 相同 → 18K token 的 KV cache 直接重用。

### 5.3 實務最佳化

- Cache TTL = 5 分鐘，持續使用自動延長
- 2026-02 起 workspace-level 隔離（非 organization-level）
- 最大效果：長 prompt 成本降低 90%、延遲降低 85%
- 自動 caching 模式：設定單一 `cache_control` 欄位，系統自動管理 breakpoint

## 6. MCP Server Advanced Patterns

### 6.1 Dynamic Tool Registration

MCP 2025 規範支援執行時動態新增/移除工具：

- `McpSyncServer.addTool()` 在連線初始化後動態註冊
- Notification 機制通知 client 工具可用性變更
- 配合 Claude Code 的 deferred loading 避免 cache 失效

### 6.2 Tool Chaining

| 模式 | 說明 | 適用場景 |
|---|---|---|
| **Sequential** | 前一步輸出提供後一步 context | 有依賴的多步驟操作 |
| **Parallel** | 獨立來源併行執行 | 多源資料收集 |
| **Conditional** | 根據結果分支 | 錯誤處理、fallback |

### 6.3 Error Handling in Chains

- **Rich context**：每個錯誤包含操作詳情、輸入、失敗原因、是否為暫時性、建議替代方案
- **Exponential backoff + jitter**：暫時性錯誤
- **Circuit breaker**：持續性失敗快速失敗
- **Fallback to secondary**：主系統失敗時自動切換

### 6.4 State Management

優先 **stateless design**（每次呼叫包含全部 context），必要時才用 session-based（time-based expiration + size limits + 使用者間隔離）。

### 6.5 2026 MCP Roadmap

MCP 規範持續演化中：非同步能力、現代化授權、企業級結構對齊、tool output schema、elicitation 能力。

## 7. 對本專案的適用建議

| 模式 | 目前狀態 | 建議行動 |
|---|---|---|
| Graduated guardrails | 已有 command hooks | 考慮加入 prompt hook 做語義安全檢查 |
| Checkpoint/rollback | 依賴 Git | 長時間實驗前手動 checkpoint 即可 |
| Token/cost tracking | 無 | 低優先；單人專案暫不需要 |
| Config versioning | `.claude/settings.json` 已 Git 追蹤 | 維持現狀 |
| Cache-aware design | CLAUDE.md 已在前 | 注意 MCP tool 數量不要頻繁變動 |
| Deferred loading | 未使用 | 若 MCP server 增加至 > 5 個，應啟用 |
| CI/CD for evals | 無自動化 eval | 可考慮 `run_batch_vcf_analysis.sh` 作為 PostToolUse hook |

---

## Sources

### Guardrail Engineering
- [Claude Code Hooks Reference](https://code.claude.com/docs/en/hooks)
- [Claude Code Hooks: All 12 Events with Examples](https://www.pixelmojo.io/blogs/claude-code-hooks-production-quality-ci-cd-patterns)
- [AGENTS.md Specification](https://agents.md/)
- [Linux Foundation AAIF Announcement](https://www.linuxfoundation.org/press/linux-foundation-announces-the-formation-of-the-agentic-ai-foundation)
- [AI Agent Guardrails Production Guide 2026](https://authoritypartners.com/insights/ai-agent-guardrails-production-guide-for-2026/)
- [Guardrails for Agentic Orchestration (Camunda)](https://camunda.com/blog/2026/01/guardrails-and-best-practices-for-agentic-orchestration/)

### Self-Healing & Error Recovery
- [Replit Snapshot Engine](https://blog.replit.com/inside-replits-snapshot-engine)
- [SagaLLM (VLDB 2025)](https://arxiv.org/abs/2503.11951)
- [Multi-Agent Failure Recovery (Galileo)](https://galileo.ai/blog/multi-agent-ai-system-failure-recovery)
- [Checkpoint/Restore Systems (Eunomia)](https://eunomia.dev/blog/2025/05/11/checkpointrestore-systems-evolution-techniques-and-applications-in-ai-agents/)

### Observability
- [AI Agent Observability Buyer's Guide (Braintrust)](https://www.braintrust.dev/articles/best-ai-observability-tools-2026)
- [AI Cost Observability (Vantage)](https://www.vantage.sh/blog/finops-for-ai-token-costs)
- [Agent Cost Optimization (Galileo)](https://galileo.ai/blog/ai-agent-cost-optimization-observability)

### Configuration as Code
- [CI/CD for Evals (Kinde)](https://www.kinde.com/learn/ai-for-software-engineering/ai-devops/ci-cd-for-evals-running-prompt-and-agent-regression-tests-in-github-actions/)
- [LLM Evaluation Guide (Braintrust)](https://www.braintrust.dev/articles/llm-evaluation-guide)
- [AGENTS.md Explained (Particula)](https://particula.tech/blog/agents-md-ai-coding-agent-configuration)

### Prompt Caching
- [Claude Prompt Caching Docs](https://platform.claude.com/docs/en/build-with-claude/prompt-caching)
- [Lessons from Building Claude Code (Thariq)](https://x.com/trq212/status/2024574133011673516)
- [How Prompt Caching Works in Claude Code](https://www.claudecodecamp.com/p/how-prompt-caching-actually-works-in-claude-code)

### MCP Patterns
- [Advanced MCP Patterns and Tool Chaining](https://dev.to/techstuff/part-4-advanced-mcp-patterns-and-tool-chaining-4ll7)
- [MCP 2026 Roadmap](https://blog.modelcontextprotocol.io/posts/2026-mcp-roadmap/)
- [MCP Tools Specification](https://modelcontextprotocol.io/specification/2025-06-18/server/tools)
