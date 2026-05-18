# Agent Cost Matrix — InterSubMod 16 Agents

> **Purpose**: 業界對齊 Walking Labs L11 observability + OpenAI cost accountability。priori 預期 cost ribbon（actual cost 由 `scripts/hooks/subagent_completion_logger.sh` SubagentStop hook 紀錄）。
>
> **Last updated**: 2026-05-18 (T3 audit 補強)
> **Source of truth (actual cost)**: `InterSubMod/docs/postmortems/subagent_completion_*.log`

## Cost Ribbon 定義

| Ribbon | $ Range | Token Range | 適用情境 |
|--------|---------|-------------|---------|
| 🟢 light | < $0.50 | < 5K out | 短任務、單檔 review、路由建議 |
| 🟡 medium | $0.50 - $2.00 | 5K-20K out | 多檔 audit、benchmark run、報告撰寫 |
| 🔴 heavy | > $2.00 | > 20K out | 全 codebase 掃描、長 trace 研究、cross-sample 分析 |

## 16 Agents Cost Matrix

| Agent | Category | Ribbon | Tokens (in/out typical) | Latency | Files Touched | Isolation |
|-------|----------|--------|------------------------|---------|--------------|-----------|
| **architect** | Planner | 🟡 | 5-15K / 3-8K | 2-5 min | `docs/plans/{YYYY}/{MM}/*.md` | — |
| **developer** | Generator | 🟡 | 10-30K / 5-15K | 5-15 min | `src/core/*.cpp`, `include/core/*.hpp`, `build/` | worktree |
| **optimizer** | Generator | 🟡 | 8-20K / 3-10K | 3-10 min | `src/core/*.cpp`, `docs/solutions/{YYYY}/{MM}/*.md` | worktree |
| **tester** | Verifier | 🟡 | 5-15K / 2-8K | 1-5 min | `docs/experiments/outputs/testing/{YYYY}/{MM}/*.md` | worktree |
| **release** | Generator | 🟢 | 3-8K / 1-5K | 1-5 min | `docs/provenance/ai_sessions/{YYYY}/{MM}/*.md`, `Dockerfile`, `CHANGELOG.md` | — (git ops 需主 branch) |
| **researcher** | Generator | 🟡 | 5-20K / 3-10K | 3-10 min | `docs/references/*.md` | — |
| **research-orchestrator** | Router | 🟢 | 2-5K / 1-3K | < 1 min | read-only (advice 不寫檔) | — |
| **headless-research** | Long-running | 🔴 | 50K-200K / 30K-100K | 1-4 hr | `research/autoresearch/headless/{ts}/*` | — |
| **paper-miner** | Long-running | 🔴 | 30K-100K / 10K-50K | 10-60 min | global `paper-miner-writing-memory.md` | — |
| **parallel-analysis** | Dispatcher | 🟡-🔴 | 10K-50K / 5K-30K | 5-30 min | `research/<topic>/figures/*.png`, `output/synthesis/research_rounds/*` | worktree |
| **parallel-benchmark** | Dispatcher | 🟡-🔴 | 10K-50K / 5K-30K | 10-60 min | benchmark TSV outputs | worktree |
| **evaluator** | Evaluator | 🟢 | 5-15K / 1-5K | 1-3 min | read-only verdict (no Write) | worktree |
| **reviewer** | Evaluator | 🟢 | 5-15K / 3-8K | 1-5 min | `docs/experiments/outputs/analysis/{YYYY}/{MM}/*.md` | worktree |
| **methodology-reviewer** | Evaluator | 🟡 | 8-20K / 3-10K | 2-10 min | JSON verdict (consumed by methodology-audit skill) | worktree |
| **literature-reviewer** | Evaluator | 🟡-🔴 | 20K-80K / 10K-40K | 10-40 min | `literature-review.md`, `references.bib`, Zotero updates | — |

## 業界對齊 Cost Accountability

| 框架 | 對應點 |
|------|------|
| Walking Labs L11 observability | SubagentStop hook 紀錄 actual cost，本表為 priori 預期 |
| OpenAI cost accountability | Per-agent ribbon 公示讓主 agent 在 dispatch 前評估 cost |
| Anthropic 1hr cache TTL | 🔴 agents 應在 1hr session window 內完成（避免 cache miss）|
| Cursor 2.0 worktree default | 11/16 agents 有 isolation（除 release 等 git 操作 + router/long-running 不適合）|

## 觸發決策參考

**Spawn 前評估 5 個維度**:
1. **Cost ribbon**: 是否 budget 容許？預算超 $5 任務應拆 sub-tasks
2. **Latency tolerance**: 用戶當輪等不等得？> 10 min 改用 `run_in_background`
3. **Isolation needed?**: 多 agent 平行寫同檔 → 必 worktree
4. **Cache hit prospect**: 重複任務（>5 次/天）→ cache 維護重要
5. **Fresh-context needs**: PASS/NEEDS_WORK verdict → 用 evaluator 類

## 紀錄與監控

- **Actual cost log**: `InterSubMod/docs/postmortems/subagent_completion_{YYYYMM}.log`
- **Hook script**: `InterSubMod/scripts/hooks/subagent_completion_logger.sh`
- **Cache telemetry**: `InterSubMod/scripts/hooks/cache_telemetry.sh`（96.8% hit baseline）

## 變更 cadence

- 每 30 天 update 一次（依 actual cost log 校正預期 ribbon）
- Cost > 預期 2× 時主動 update + 加入 postmortem
- 加入新 agent → 同時 update 本表 (avoid drift)
