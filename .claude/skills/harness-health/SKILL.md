---
name: harness-health
description: "Harness 自我稽核固定流程 — 跑 6 顆 L0 健康燈（count drift / Hard-Gate 真實性 / tier-gate wired / state↔doc 同步 / ledger 新鮮度 / queue hygiene）+ 偵測 exit-2 hook 被 '|| exit 0' neuter + stale compile marker，產出 snapshot JSON + standalone HTML 儀表板（localStorage 可 diff 上次）。read-only、無 agent、無 web、純磁碟 grep。USE WHEN：改了 skill/agent/hook/settings 後想確認沒漂移、月度自我稽核、「harness health」「儀表板」「哪些要改進」「漂移檢查」「config drift」、SessionStart 建議跑、commit harness 變動前。SKIP WHEN：研究 cycle 狀態（用 /cycle-state）、專案假說看板（用 /research-dashboard）、跨 cycle 證據鏈審計（用 /provenance-tier-audit）、純 build / commit / docs 寫作、無 harness 變動的研究分析。"
---

# /harness-health — Harness 自我稽核固定流程

每次「更新或任務改動」跑一次（<5 秒），確認「哪些是何、哪些要改進」。**不改任何檔案**（snapshot/HTML 寫入 `state/health_snapshots/` 例外）。

## Phase & Chain Position

- **層級**：元方法論 / harness 自我觀測（與 7-Phase 研究 cycle 正交）。
- **不在研究 waterfall 內**：研究 cycle 狀態看 `/cycle-state`；本 skill 看的是 **harness 本身**（skill/agent/hook/settings/doc 的一致性）。
- **進入點**：(1) 改了 skill/agent/hook/settings 後；(2) 月度 cadence；(3) commit harness 變動前；(4) SessionStart 建議；(5) 用戶問「哪些要改進」。

## 執行

```bash
python3 scripts/harness_health.py
```

產出：
- 終端 6 顆 L0 燈表（GREEN/YELLOW/RED + 每燈可驗證證據）
- `state/health_snapshots/YYYYMMDD-HHMM.json`（snapshot，保留最近 ~10 份）
- `state/health_snapshots/dashboard.html`（standalone，開瀏覽器看；localStorage 自動 diff 上次 snapshot → 「哪些燈變了」一眼可見）

**6 顆 L0 燈**：
1. **COUNT DRIFT** — disk 實際 skill/agent/hook 計數 vs CLAUDE.md §3/§8/§4 宣稱
2. **HARD-GATE TRUTH** — wired hook 中真 `exit 2`（未被 `|| exit 0` neuter）數 vs 宣稱；偵測 neutered bug + 區分刻意 soft（verify_gate）+ stale compile marker
3. **TIER-GATE WIRED** — `pre_tier_upgrade_check` 是否真 wired
4. **STATE↔DOC SYNC** — `active.json` recently_concluded vs CURRENT_FOCUS 待辦描述衝突
5. **LEDGER FRESH** — `evidence_ledger` 最新 entry date vs CURRENT_FOCUS header（ledger 領先 = tier over-claim 風險）
6. **QUEUE HYGIENE** — `hypothesis_queue` queued 數 + 是否引用 concluded-DEAD 方向

**紅旗 → 建議動作**（終端與 HTML 都會列）：
- 燈 2 紅（neutered Hard Gate）→ 修 settings wiring 移除 `|| exit 0`（**先確認 stale marker，否則恢復 gate 會擋 commit**）
- 燈 4 紅 → `/provenance-tier-audit` 或 `/cycle-state` 核對；**勿手改 state.json**
- 燈 5 黃/紅 → 回寫最新 ledger 到 CURRENT_FOCUS
- 燈 6 紅 → `/pivot-direction` 降權 stale queue

## Dependencies

- **Uses**：`scripts/harness_health.py`（唯一執行體；純 stdlib，無外部套件）
- **Reads**（read-only）：`.claude/CLAUDE.md`、`.claude/settings.local.json`、`scripts/hooks/*.sh`、`docs/CURRENT_FOCUS.md`、`research/autoresearch/evidence_ledger.jsonl`、`state/active.json`、`research/autoresearch/hypothesis_queue.json`、`/tmp/ism_cpp_pending_compile.txt`
- **Writes**：僅 `state/health_snapshots/{stamp}.json` + `dashboard.html`
- **Used-by**：SessionStart 建議；commit harness 變動前的人工 checklist
- **與既有 skill 分工**：`/cycle-state`=研究 cycle 級；`/research-dashboard`=專案假說級；`/provenance-tier-audit`=跨 cycle 證據鏈；**本 skill=harness 設定/結構級**（互補不重疊）
- **設計依據**：`docs/references/migration/20260531_harness_audit_dashboard_design_01.md`

## Failure Mode & Diagnostics

| 症狀 | 可能原因 | 診斷 |
|------|----------|------|
| 燈 1 紅但你剛同步過 | claim 正則沒抓到新格式 | 檢查 CLAUDE.md §3「NN 個 SKILL.md」格式；改 `harness_health.py:light_count` 正則 |
| 燈 2 一直紅 | `pre_commit_compile_check`/`kb_schema_check` 仍被 `|| exit 0` masked（已知 2026-05-31 bug，待用戶決定修） | 見 design doc §「neutered Hard Gate」；修前先清 stale marker |
| 燈 4 誤報 | heuristic 比對 cycle id 過寬 | 人工核對 active.json vs CURRENT_FOCUS；必要時調 `light_state` 正則 |
| HTML 不更新 diff | localStorage 被清 / 換瀏覽器 | 第一次跑會建 baseline，第二次起才有 diff |
| script 報錯 | 某來源檔不存在 | 各 `_read` 已 try/except 回空字串，不會 crash；檢查路徑常數 |
