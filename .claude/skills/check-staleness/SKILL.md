---
name: check-staleness
description: Phase 2 PRECHECK gate of the Resilient Waterfall harness — verify that a cycle's preconditions (C++ binary, dataset, upstream reports) are still fresh before pilot. Reads `state/cycles/{cycle_id}/plan.json`, compares against `state/invalidation/*.jsonl`, and writes `precheck.json`. USE WHEN：「check staleness」「P2 precheck」「驗 precondition」「驗證前置條件」「跑 pilot 前」、cycle 從 P1 → P2 transition 時。涉及 InterSubMod/state/ 內 plan.json 與 invalidation 紀錄。 SKIP WHEN cycle 仍在 P0/P1 PLAN 階段（plan.json 未 ready）、純 build / commit、已通過 P2 進入 P3 PILOT、不涉及 cycle state 變更、純 docs 寫作。
allowed-tools: ["Bash", "Read", "Glob", "Grep"]
user-invocable: true
tags: ["harness", "precheck", "validation", "freshness"]
---

# /check-staleness — Phase 2 PRECHECK Gate

驗證一個 cycle 的前置條件是否仍 fresh，是 Resilient Waterfall harness 的 P2 gate。

> **設計來源**：`~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md` §3.4

## 何時使用

- cycle 進入 P2 PRECHECK 階段（`state.json.phase` 從 `P1_PLAN` → `P2_PRECHECK`）
- 用戶手動跑：`/check-staleness <cycle_id>`
- 主對話 SessionStart 自動掃描所有 active cycles 是否 BLOCKED

## 何時 NOT 使用

- 沒有 plan.json 的 cycle（先做 `/cycle-init` 與 P1 PLAN）
- 純讀檔分析、不依賴 binary / dataset 的查詢任務

## 工作流程（3 個 freshness check）

### Step 1 — 讀取 plan.json
```bash
CYCLE_ID="$1"
PLAN="InterSubMod/state/cycles/${CYCLE_ID}/plan.json"
test -f "$PLAN" || { echo "ERROR: plan.json not found at $PLAN"; exit 2; }
```
從 plan.json 取出：
- `preconditions.binary_version` — git commit SHA
- `preconditions.dataset_id` — dataset 識別字串
- `preconditions.upstream_reports[]` — InterSubMod/-prefixed 路徑

### Step 2 — Binary freshness check
```bash
python3 InterSubMod/.claude/skills/check-staleness/_staleness_check.py "$CYCLE_ID"
```
比對 plan 標的 binary SHA 與目前 HEAD：
- **fresh**: 兩者相同
- **stale**: HEAD 與 plan 標的之間有新 commit（`stale_distance` 為 commit 數）
- **missing**: plan 標的 SHA 在 git 找不到
- **skipped**: plan.binary_version=null（純分析 cycle）

### Step 3 — Dataset schema check
對每個 `dataset_id`，跑 schema sanity check：
- 必須有 `caller_af` 欄位（避免 merged AF 陷阱）
- LOH 覆蓋率 ≥ 50%（避免 HCC1395 phase1_new 殘缺類問題）
- 其他 known pitfalls 由 `_staleness_check.py` 內 `DATASET_PROBES` 表驅動

### Step 4 — Upstream report freshness check
對 `upstream_reports[]` 中每個路徑：
- 讀 frontmatter 的 `binary_version` / `last_verified` 欄位（如有）
- 檢查 `state/invalidation/stale_marks.jsonl` 是否標記此報告為 stale
- 檢查報告是否含 `[RETRACTED` 標記

### Step 5 — 寫出 precheck.json
依 `state/schemas/precheck.schema.json` 格式輸出至 `state/cycles/{cycle_id}/precheck.json`。

verdict 邏輯：
- 任一 check 為 `stale` / `schema_mismatch` / `retracted` → `BLOCKED`
- 全部 `fresh` 或 `skipped` → `PASS`
- 含 `warn`-級別問題（如 LOH 覆蓋 50-58%）→ `WARN`

## 用戶 override 路徑

若 verdict=`BLOCKED` 但用戶確認可繼續（例：知道 stale 影響有限），**不直接編輯 precheck.json**，而是：
1. 用戶提供 override reason
2. skill 在 precheck.json 加 `user_override` 區段（schema 已預留）
3. 將 cycle 從 BLOCKED 推進到 P3，但 evaluation 在 P5 必動（不論 tier）

## 輸出格式（人類可讀摘要）

```
[/check-staleness] cycle_id=20260503-2330-loh-kde-quantify
  Binary:           ✓ fresh (HEAD=8d0a0c8, plan=8d0a0c8)
  Dataset:          ✗ schema_mismatch — caller_af column missing in master_2026-04-01
  Upstream reports: ✓ 3/3 fresh
verdict: BLOCKED
blocking_reasons:
  - dataset master_2026-04-01 lacks caller_af column (known pitfall: merged AF trap)
next steps:
  - Replace dataset_id with master_extended_2026-04-22_KDE-corrected
  - Or override with reason recorded in precheck.user_override
```

## 與既有 skills 整合

| 場景 | 此 skill 做什麼 | 既有 skill 接手 |
|---|---|---|
| 發現 stale binary | 標記 `BLOCKED` | 用戶呼叫 `/cpp-change` 或重跑 benchmark |
| 發現 dataset schema mismatch | 標記 `BLOCKED` + 列具體欄位 | 用戶查 `/known-pitfalls` 確認影響 |
| 發現 upstream retracted | 標記 `BLOCKED` | 用戶用 `/pivot-direction` 重定 cycle |
| 全 fresh | 標記 `PASS` | cycle 推進到 P3，呼叫 `/feature-layered-observation` 或 pilot 流程 |

## 限制 (Path A)

- 不對 LlamaIndex 索引查 precedent — 由 `/run-evaluator` (P5) 負責
- 不主動建議「下一個 fresh 的 dataset」— 只報告 stale，由用戶決定
- 不重新跑 benchmark — 只比對已存在的 artifacts metadata

---

## 與 /scientific-rigor 元方法論的關係

本 skill 為 `/scientific-rigor §7.2 可重現性 checklist` 的**研究實驗層級具體驗證**:
- binary / dataset / upstream 新鮮度檢核 對應 `/scientific-rigor §7.2` 7 項中的「環境記錄 + 數據版本 + 中間產出存檔」
- **場景分流明示**: 研究實驗 reproducibility（seed / data version / commit hash）→ 本 skill；程式碼 build/test 驗證 → `/verification-loop`
- `/scientific-rigor §11 協作圖 step 4 → 5` 之間若 cycle plan.json 已存在，必經本 skill PRECHECK

**級聯觸發**: `/cycle-init` P0 → `/research-loop` P1 → 本 skill P2 → `/feature-layered-observation` P3

---

## Phase Chain Position & Dependencies

- **Phase**: P2 STALENESS（在 /research-loop P1 之後、/feature-layered-observation P3 之前）
- **Upstream**: `/research-loop`（plan.json 已生成）
- **Downstream**: `/feature-layered-observation`（pass）或 `/pivot-direction`（fail）
- **Reads**: `hypothesis_queue.json` + `evidence_ledger.jsonl` + `MEMORY.md` Concluded
- **Writes**: Staleness verdict + reopen-eligibility check

## Failure Mode & Diagnostics

| 症狀 | 可能原因 | 修法 |
|------|---------|------|
| Cycle 啟動但 hypothesis 已 NEGATIVE 30d | 未跑 staleness check | 強制 §8.3.1 reopen 三條件對齊 |
| Reopen 條件模糊 | C1/C2/C3 未具體化 | 要求填 specific 新數據 / 新方法 / 新前置 |

