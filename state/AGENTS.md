<!--
建立時間: 2026-05-09 03:09
目標: 提供 Codex 在 state/ 目錄內操作 Resilient Waterfall runtime state 的局部規則
處理範圍: state/active.json、state/cycles、state/retro_cycles、state/invalidation、state/schemas
關聯檔案:
  - ../AGENTS.md
  - ../state/CLAUDE.md
  - ../.agents/skills/README.md
-->

# State Directory - AGENTS.md

> Resilient Waterfall harness 的執行時 state 與 audit 紀錄根目錄。
> 本檔是 `state/CLAUDE.md` 的 Codex 版本；`CLAUDE.md` 保留給 Claude Code legacy，Codex 以本檔與 repo root `AGENTS.md` 為準。

## 結構與職責

```
state/
├── active.json              <- 索引：當前 active cycles（<=5）
├── cycles/                  <- 真實 active cycle 工作目錄
│   └── {cycle_id}/
│       ├── state.json       <- cycle 主狀態（schema 驗證）
│       ├── plan.json
│       ├── precheck.json
│       ├── pilot.json
│       ├── generalize.json
│       ├── evaluation.json
│       └── reflection.log
├── retro_cycles/            <- Drill 1 retro fixture（不在 active.json）
├── cycles_archived/         <- P6 完成後搬移（gitignored 由 cycle）
├── invalidation/            <- stale_marks.jsonl + binary_versions.jsonl
└── schemas/                 <- JSON schema 真理；不可隨意改
```

## 不可手改清單

- `state.json` / `plan.json` 等 cycle artifact：必須由對應 Codex skill 寫入（`cycle-init` / `research-loop` / `check-staleness` / `run-evaluator` / `conclude-research`）。手改會破壞 schema validation 與 history audit 鏈。
- `active.json`：只能由 cycle 初始化或 cycle 終結流程更新；schema 限 maxItems=5。
- `schemas/*.schema.json`：改 schema 等於改現有 cycle 行為。需有 plan v1.x 明確 Decision Log，且不得破壞既有 cycle 解析。

## 可手動操作

- 新建 retro_cycles fixture：可手寫 5 件 artifacts 到 `retro_cycles/{event_id}/`，用於 regression test（例如 Drill 1 retro 模式）。
- `invalidation/stale_marks.jsonl` append-only：可手動加 stale 紀錄；若未來有 hook 或 script，也只能 append。

## 何時碰 / 不該碰

| 情境 | 應做 | 不該 |
|---|---|---|
| 啟動新 cycle | 使用 `$cycle-init` 或同名 Codex skill | 直接 mkdir `cycles/{id}/` |
| 手動推進 phase | 使用對應 phase skill | 直接編輯 `state.json` 的 `phase` |
| cycle 撤回或 retract | 使用 `$conclude-research` 寫 NEGATIVE | 直接刪 cycle 目錄 |
| 想看現狀 | 使用 `$cycle-state` | 自己 cat 多個 `state.json` |

## 與 evidence_ledger 同步規則

`research/autoresearch/evidence_ledger.jsonl` 是歷史 SoT，`state.json` 是現況快照。崩潰恢復時可由 ledger 重建 state。寫入順序採「先寫 ledger，再更新 state」。

## Codex 套用方式

- Claude Code slash commands 在 Codex 端不直接搬移；舊文件中的 `/cycle-init`、`/cycle-state` 等同於 `$cycle-init`、`$cycle-state` 或自然語言觸發同名 skill。
- Claude hooks 不直接搬移；需要 gate 時，優先做成 skill 內 checklist、明確驗證腳本，或在執行前後由 Codex 主流程呼叫 deterministic script。
- Claude agents 不直接搬移；需要平行化時，只在使用者明確要求或任務本身已授權平行工作時，依 Codex subagent 規則拆成獨立輸出目錄的 bounded subtasks。

## 未來 paths 提示

若 Codex 後續支援 repo-local skill path filter，state 相關 skill 可使用：

```yaml
paths: ["state/cycles/**/*.json", "state/active.json"]
```
