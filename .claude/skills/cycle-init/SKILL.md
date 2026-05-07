---
name: cycle-init
description: Phase 0 REGISTER gate of the Resilient Waterfall harness — initialize a new short-validation research cycle. Allocates `cycle_id` (YYYYMMDD-HHMM-{slug}), creates `state/cycles/{cycle_id}/state.json` skeleton + updates `state/active.json` index. Hands off to research-loop for P1 PLAN. **Distinct from /init-research** which scaffolds multi-week project trees in `research/`. USE WHEN：「cycle init」「new cycle」「P0 REGISTER」「啟動新假說 cycle」「state/cycles 建立」、收到 hypothesis_id 後啟動驗證流程時。
allowed-tools: ["Bash", "Read", "Glob", "Grep", "Write"]
user-invocable: true
tags: ["harness", "register", "cycle-init", "p0"]
---

# /cycle-init — Phase 0 REGISTER Gate

短驗證 cycle 的入口 skill。把一個假說候選（已注入 hypothesis_queue.json 或臨時 slug）轉為 harness state machine 中的 active cycle，產出 `state/cycles/{cycle_id}/state.json` 與更新 `state/active.json`。

> **設計來源**：`~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md` §3.7, §4.5.5

## Phase & Chain Position

```
                  hypothesis_queue.json
                       │
                       ▼
[problem-framing-ideation] ──optional──▶ [inject-hypothesis] ──H-id──▶ [/cycle-init] ◀── YOU ARE HERE (P0)
                                                                          │ writes state.json
                                                                          ▼
                                                                  [research-loop] (P1 PLAN)
                                                                          │ writes plan.json
                                                                          ▼
                                                                  [/check-staleness] (P2)
                                                                          │
                                                                          ▼
                                                                  [feature-layered-observation] (P3)
                                                                          │
                                                                          ▼
                                                                  [multi-sample-consistency] (P4)
                                                                          │
                                                                          ▼
                                                                  [/run-evaluator] (P5)
                                                                          │
                                                                          ▼
                                                                  [conclude-research] (P6)
```

- **P0 REGISTER** entry → next is **P1 PLAN**（research-loop）
- 與 `/init-research` 互補不重疊：`/init-research` = `research/{project_name}/` 多週級 project；本 skill = `state/cycles/{cycle_id}/` 單假說短驗證

## Dependencies

| 方向 | 對象 | 用途 |
|---|---|---|
| **Uses** | `state/schemas/state.schema.json` | 寫 state.json 須符合 schema |
| **Uses** | `state/schemas/active.schema.json` | 更新 active.json 須符合 schema |
| **Used by** | `research-loop` (P1) | 從 state.json 讀 cycle_id 與 hypothesis_id |
| **Used by** | `/cycle-state` (dashboard) | 從 active.json + 各 state.json 顯示概覽 |
| **Reads** | `research/autoresearch/hypothesis_queue.json` (optional) | 驗證 hypothesis_id 存在；若缺則警告但繼續 |
| **Writes** | `state/cycles/{cycle_id}/state.json` | 全新 cycle 骨架 |
| **Writes** | `state/active.json` | append cycle entry，updated_at refresh |

## 何時使用

- 用戶決定為某個假說啟動驗證 cycle（例：「驗證 LOH 反 KDE 重建假說」）
- problem-framing-ideation 收斂出具體假說後 → inject-hypothesis 註冊 → /cycle-init 啟動 cycle
- 全自動模式（headless / `auto`）：從 hypothesis_queue.json 取 priority 最高的 pending 假說自動啟動

## 何時 NOT 使用

- 多週級研究**專案** scaffolding（用 `/init-research`，輸出到 `research/`）
- 已有 active cycle 在跑同樣假說（先用 `/cycle-state` 查現況）
- 寫 plan.json（research-loop 的工作；本 skill 只建空骨架）
- 跑 pilot / generalize / evaluator（下游 phase 的工作）

## 工作流程（5 步）

### Step 1 — 收集輸入
| 必填 | 說明 |
|---|---|
| `hypothesis_id` 或 `slug` | 來自 hypothesis_queue.json 或臨時字串（如 `loh-kde-rescale`） |
| `title` | cycle 標題（≤60 字） |

| 選填 | default |
|---|---|
| `priority` | 50 (錨點 80+/60-79/40-59/20-39/<20) |
| `binary_version` | null（pure-analysis cycle） |
| `dataset_id` | null（待 P1 PLAN 填入） |
| `upstream_reports` | `[]` |

### Step 2 — 產生 cycle_id
格式：`YYYYMMDD-HHMM-{slug}`（schema regex `^[0-9]{8}-[0-9]{4}-[a-z0-9_-]+$`）
- date/time = 當前 UTC
- slug 從 title 或 hypothesis_id 推導：lowercase, 替換非 `[a-z0-9_-]` 為 `-`, 去除連續 `-`, 截斷至 ≤30 字
- 衝突偵測：若 `state/cycles/{cycle_id}/` 已存在 → 加 `-2` / `-3` 後綴重試

### Step 3 — 建立 state.json
依 `state.schema.json` 寫入：
```json
{
  "schema_version": "1.0",
  "cycle_id": "<generated>",
  "title": "<input>",
  "hypothesis_id": "<input or null>",
  "phase": "P0_REGISTER",
  "tier": "pending",
  "verdict": "active",
  "priority": <input or 50>,
  "started_at": "<UTC ISO>",
  "last_advanced_at": "<same as started_at>",
  "preconditions": {
    "binary_version": <input or null>,
    "dataset_id": <input or null>,
    "upstream_reports": <input or []>
  },
  "open_gates": [],
  "artifacts": {
    "plan": null, "precheck": null, "pilot": null,
    "generalize": null, "evaluation": null
  },
  "history": [
    {"timestamp": "<UTC>", "from_phase": null, "to_phase": "P0_REGISTER", "actor": "user", "note": "/cycle-init"}
  ]
}
```

### Step 4 — 更新 active.json
- Append cycle entry to `cycles[]`
- Refresh `updated_at`
- **Soft warn** if `len(cycles) > 5`（schema cap maxItems=5）— 提示用戶先收斂或主動 archive，不阻擋 write

### Step 5 — 輸出摘要
```
[/cycle-init] cycle_id=20260507-1430-loh-kde-rescale
  hypothesis_id: H-2026-05-07-001
  title: "LOH 反 KDE 重建假說"
  phase: P0_REGISTER
  priority: 50
  state.json: state/cycles/20260507-1430-loh-kde-rescale/state.json
  active.json updated (3 active cycles total)

Next: invoke /research-loop or research-loop skill to draft plan.json (P1 PLAN)
```

## Failure Mode & Diagnostics

| # | Failure | Cause | Detection | Remedy |
|---|---|---|---|---|
| 1 | cycle_id 衝突 | 同分鐘建 ≥2 個 cycle | dir already exists | 自動加 `-2` 後綴重試 |
| 2 | active.json corrupt | schema 不合（手動編輯導致） | json.load 失敗或缺欄位 | 拒絕寫入；提示用戶手動修正後重跑 |
| 3 | hypothesis_id 不在 queue | 用戶輸入 typo 或臨時假說 | grep hypothesis_queue.json 找不到 | 警告但繼續；state.json 仍寫 hypothesis_id（允許探索性） |
| 4 | active cycles ≥ 5 | 先前 cycles 未收斂 | len(active.cycles) >= 5 | Soft warn 列現有 cycle_id，建議先 archive；不阻擋 |
| 5 | state/cycles/ 不存在 | 首次跑 cycle-init | os.path.isdir | 自動 mkdir -p |
| 6 | binary_version 對不上 git log | 用戶輸入錯 SHA | git cat-file -e 失敗 | 警告，記入 state.json `preconditions.binary_version_unverified=true` 但寫入 |
| 7 | dataset_id 含 known-pitfall keyword | 如 "merged" / "pileup_symlink" | substring match against pitfalls_table.json keywords | 警告，建議改名或預期 P2 PRECHECK 會 catch |
| 8 | slug 為空（title 非 ASCII 全部過濾掉） | 全中文 title 無 ASCII 字元 | re.sub 後字串為空 | 用 hypothesis_id 或 fallback 到 `unnamed-{N}` |

## 與其他 skills 整合

| 場景 | 此 skill 做什麼 | 既有 skill 接手 |
|---|---|---|
| 用戶說「啟動 LOH KDE 重建驗證」 | 建 state.json + 加 active.json | research-loop（P1 PLAN）寫 plan.json |
| 全自動模式取 queue 最高 priority 假說 | 同上但 hypothesis_id 從 queue 拉 | research-loop 或 /research-dashboard 顯示新 cycle |
| 用戶說「先腦力激盪再決定」 | **不**啟動本 skill | problem-framing-ideation 收斂後再回來呼叫 /cycle-init |
| 用戶說「啟動長期專案 X」 | **不**啟動本 skill | /init-research 建 research/{project} |

## 限制 (Path A)

- 不接 LangGraph runtime（純 Python script）
- 不自動驗 binary_version 對應的 git commit（只警告，留給 P2 PRECHECK 嚴格驗）
- active.json 寫入無 file lock（單研究者場景；多人寫時可能 race condition）
- 不寫 evidence_ledger.jsonl（那是 P5/P6 的工作）

## 維護規則

- 改 state.schema.json 後須同步改本 skill 的 Step 3 範本
- 改 active.schema.json 後須同步改 Step 4 邏輯
- 新增 phase enum（如未來加 P0.5_FRAMING）須回頭改本 skill 的 history 初始化
