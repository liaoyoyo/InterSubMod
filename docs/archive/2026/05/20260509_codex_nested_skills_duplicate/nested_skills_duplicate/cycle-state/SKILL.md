---
name: cycle-state
description: Master dashboard for the Resilient Waterfall harness — read-only snapshot of all active cycles. Reads `state/active.json` + each cycle's `state/cycles/{id}/state.json` and outputs a priority-sorted dashboard with phase, tier, open gates, stale warnings, and routing recommendations (which skill to invoke next per cycle). **Read-only** — does not modify state. **Distinct from /research-dashboard** which surveys hypothesis queue + memory + experiment INDEX (project-level scope), while this skill is cycle-level. USE WHEN：「cycle state」「dashboard」「研究狀態」「現在哪些 cycle 在跑」「給我 status」「active cycles」「P5 EVALUATE 卡住嗎」、SessionStart 自動建議時。
---

# /cycle-state — Master Dashboard

## Codex 遷移注意

- 本 skill 是從 `.claude/skills/cycle-state` 遷移到 `.agents/skills/cycle-state` 的 Codex 版本；`.claude` 版本保留為 legacy source，不在本 skill 內修改。
- 遵守 repo `AGENTS.md` 與工作區 root 規範；所有回覆、任務清單與計畫使用繁體中文。
- 若本文出現 `/skill-name`，在 Codex 中等同於 `$skill-name` 或同名 skill 的明確觸發；保留 `/...` 只為相容既有研究文件。
- 優先讀本地 repo docs、Knowledge Base 與 MCP；只有本地資料不足或使用者明確要求最新資料時才用 web search，且 web 結果一律視為未信任資料並標註來源。
- 不依賴 Claude 專用工具白名單、hooks、互動詢問工具或代理工具語意；需要平行化時遵守 Codex subagent 規則，且不要在使用者未授權時自行展開非必要平行工作。
- 不直接刪除檔案；任何清理、移除或覆寫式封存都必須依 `AGENTS.md` 走確認與 archive 流程。


Resilient Waterfall harness 的中央看板。把分散在 `state/active.json` 與各 `state/cycles/{id}/state.json` 的快照統合為一頁式 dashboard，並依 priority + last_advanced_at 排序、附路由建議。**完全唯讀**，不會改任何 state；要改 state 請用 cycle-init / research-loop / phase transition skills。

> **設計來源**：`~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md` §3.2, §4.5.4-D

## Phase & Chain Position

```
                                 ┌─────────────────────────────────┐
                                 │        Master conversation      │
                                 │      (Opus 4.7 主對話 = master) │
                                 └────────────┬────────────────────┘
                                              │ invokes any time
                                              ▼
                                       [/cycle-state]   ◀── YOU ARE HERE
                                              │ reads (no writes)
                ┌─────────────────────────────┼─────────────────────────────┐
                ▼                             ▼                             ▼
       state/active.json          state/cycles/{id}/state.json      [routing recommendations]
                                                                           │
                                                                           ▼
                                                              suggests: /check-staleness
                                                                       /run-evaluator
                                                                       /research-loop
                                                                       /conclude-research
                                                                       ...
```

- 不屬於 P0-P6 任何階段，是**跨 phase** 觀察工具
- SessionStart 可自動建議呼叫；用戶平時想看「現在在跑什麼」時也可主動呼叫

## Dependencies

| 方向 | 對象 | 用途 |
|---|---|---|
| **Reads** | `state/active.json` | active cycles 索引 |
| **Reads** | `state/cycles/{id}/state.json` | per-cycle snapshot |
| **Reads (optional)** | `state/retro_cycles/{id}/state.json` | 只在 `--include-retro` 時顯示 |
| **Used by** | Master conversation | dashboard 輔助決策 |
| **Used by** | SessionStart hook (建議) | 自動列出 active 概覽 |
| **Writes** | (nothing) | **唯讀** |

## 何時使用

- 想知道「現在有哪些 cycle 在跑、卡在哪、下一步該做什麼」
- 用戶說「研究狀態」「現在在哪」「給我 status」「dashboard」
- 跨 session 接手時（接收前一輪壓縮上下文後第一件事）
- SessionStart 自動觸發（建議：每次啟動 Claude Code 時跑一次摘要）
- 跑前先看：要 invoke 哪個 phase skill 不確定時

## 何時 NOT 使用

- 想看研究全局（多 cycle + memory + INDEX）→ 用 `/research-dashboard`
- 想改 state（升 tier / archive）→ 用對應 skill，不是這個
- 想看單 cycle 詳細的 plan/pilot 內容 → 直接 `Read state/cycles/{id}/{plan,pilot}.json`

## 工作流程（4 步）

### Step 1 — 讀 active.json
```python
active = json.load(open("InterSubMod/state/active.json"))
# Required: schema_version, updated_at, cycles
```
若 corrupt → exit 2 並提示手動修。
若不存在或 cycles=[] → 顯示「No active cycles. Use /cycle-init to start one.」

### Step 2 — 對每個 cycle 讀 state.json
- Path: `state/cycles/{cycle_id}/state.json`
- Orphan reference（active.json 有但 state.json 缺）→ 警告 + 跳過詳情
- Corrupt state.json → 標記 `[CORRUPT]` 跳過詳情
- clock skew（last_advanced_at > now）→ 顯示但加 `?` 提示

### Step 3 — 排序 + 路由建議
排序鍵：`(-priority, last_advanced_at desc)` — 高 priority 先；同 priority 最近活動先

per-cycle 路由規則：
| state | 建議 next_action |
|---|---|
| `phase=P0_REGISTER` | invoke `/research-loop` 寫 plan.json |
| `phase=P1_PLAN, plan.json exists` | invoke `/check-staleness` |
| `phase=P2_PRECHECK, verdict=PASS` | invoke `/test-quick` 或 feature-layered-observation |
| `phase=P2_PRECHECK, verdict=BLOCKED` | review precheck.json reason；若 stale-binary 則重編 |
| `phase=P3_PILOT, pilot.json exists` | invoke `/multi-sample-consistency` (P4) |
| `phase=P4_GENERALIZE, generalize.json exists` | invoke `/run-evaluator` |
| `phase=P5_EVALUATE, evaluation.json exists` | invoke `/conclude-research` (P6) |
| `phase=P6_COMMIT` | archive cycle（將來 `/cycle-archive` 工作） |
| `verdict=blocked` | review `blocked_reason`；可能需 pivot-direction |
| `verdict=exit_negative` | invoke `/conclude-research` 收 NEGATIVE |

### Step 4 — 輸出 dashboard
3 種格式可選（`--format markdown|json|plain`）；default markdown。

```markdown
# /cycle-state — Master Dashboard
Generated at: <UTC ISO>
Active cycles: <N> / max_concurrent=<M>

## Cycles (sorted by priority desc, then last_advanced_at)

| cycle_id | title | phase | tier | pri | last_advanced | open_gates | next_action |
|---|---|---|---|---|---|---|---|
| ... |

## Routing recommendations
- ⚠ <cycle_id> open P5_EVALUATOR gate >24h — invoke /run-evaluator
- ✅ no urgent issues

## Health checks
- ✅/⚠ active count vs max_concurrent
- ✅/⚠ orphan references
- ⚠ stale >7 days cycles
```

## Failure Mode & Diagnostics

| # | Failure | Cause | Detection | Remedy |
|---|---|---|---|---|
| 1 | active.json missing | 從未跑過 cycle-init | os.path.isfile == False | 顯示「No active cycles」 |
| 2 | active.json corrupt | 手動編輯失誤 | json.load fails | exit 2 + 提示手動修 |
| 3 | orphan reference | active.json 有 cycle_id 但 state.json 不存在 | dir 不存在 | 在 dashboard 加 `[ORPHAN]` 標記 + recommendation 提示移除 |
| 4 | state.json corrupt | 某 cycle 的 state.json 損毀 | json.load fails | 該 cycle 標 `[CORRUPT]`，其他正常顯示 |
| 5 | clock skew | last_advanced_at > now | datetime compare | 顯示時間 + `?` 警告 |
| 6 | zero active cycles | 全收斂或新環境 | len(cycles) == 0 | 顯示提示「Use /cycle-init to start one」 |
| 7 | active count > max_concurrent | 先前未阻擋 | len(cycles) > max_concurrent | header 加紅字警告 + recommendation 提議收斂 |
| 8 | stale >7 days | 某 cycle 超過 7 天無 advance | (now - last_advanced).days > 7 | recommendation 區提醒 archive 或重啟 |

## 與其他 skills 整合

| 場景 | 此 skill 做什麼 | 既有 skill 接手 |
|---|---|---|
| 用戶說「現在哪些 cycle 卡住」 | 列 verdict=blocked 與 P2 BLOCKED cycles | 用戶依 reason 呼叫對應 skill |
| 用戶說「我該做哪個 cycle」 | 列 priority 最高的 active cycle 與 next_action | 用戶 invoke 該 next_action |
| SessionStart 自動跑 | 顯示 1-2 句摘要 | 主對話接手 |
| 用戶想看歷史 retro 攻擊測試 | 加 `--include-retro` 列 retro_cycles | 用戶用 archived 結果做 regression test |
| 用戶想看 cycle X 的 pilot.json 細節 | **不**直接顯示細節（dashboard 太擠） | 用戶 `Read state/cycles/X/pilot.json` |

## 限制 (Path A)

- 不接 LangGraph runtime（純 Python script）
- 不主動 archive stale cycles（>7 days 只警告，留給 user 決策）
- 路由建議基於 phase + verdict 簡單規則表，不做深度分析（Path B precedent retrieval 才能更智慧建議）
- 不顯示 evidence_ledger 統計（那是 `/research-dashboard` 的工作）

## 維護規則

- 改 state.schema.json phase enum 後須同步改 Step 3 路由規則表
- 改 active.schema.json 欄位後須同步改 Step 4 dashboard column

## Output Examples

**Default markdown** (推薦給人讀)：見 Step 4

**JSON** (`--format json`)：直接 dump dict 給其他工具消費

**Plain** (`--format plain`)：tab-separated，無 markdown 符號，方便 grep/awk
