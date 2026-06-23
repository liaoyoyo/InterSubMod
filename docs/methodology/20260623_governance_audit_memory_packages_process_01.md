---
title: 治理稽核紀錄 — 記憶 / 可重用套件 / 流程控管（insights 觸發）
date: 2026-06-23
status: done
type: governance_record
build_branch: research/subclonal-reconstruction-202606
audit_source: workflow wj8dawt6j（4 agent / 3 面）+ 主回合 §13.7 逐項 grep 驗證
data_sources:
  - state/active.json
  - scripts/hooks/script_lint_advisor.sh
  - scripts/tasks/_build_full_tree.py
---

# 治理稽核紀錄（2026-06-23）

/ insights 報告（172 sessions）+ 用戶要求「確認需治理處，特別是記憶與重複使用套件，對流程控管」觸發。
multi-agent 稽核 + 主回合逐項第一手驗證 → 健康度**良好**（多數 insights 建議其實已覆蓋）；落地 5 項真修正，restraint 排除 ~8 項已覆蓋/不值得建。

## A. 記憶治理（memory）

| 項 | 動作 | 驗證 |
|---|---|---|
| MEMORY.md 索引超限 | 33,804 → **24,968 bytes**（-26%，低於 24,986 上限）→ partial-load 解除 | `wc -c`；14 條超長索引 trim 為 hook |
| **分層揭露驗證** | 14 條被 trim 的 detail **逐條 grep 確認已在 topic 檔**（不丟資訊）| 14/14 ✓ |
| 索引完整性 | 所有非-archived topic 檔仍在索引 | grep 全 ✓ |
| review-format 合併 | `feedback_html_for_conclusions_md_for_final_check` 併入 `feedback_review_format_html_ask_md`（merged_from / merged_into，**never 刪**，B 標 status:archived）| — |
| double-dip 交叉連結 | `project_tumor_only_axis_negative...` 加 SUPERSEDED → `project_phylo_subcluster_labeling_doubledip_fix` | — |

**分層揭露機制確認可用**：index = ≤200 byte hook（load-bearing verdict + 🔴 + topic 連結）；detail 在 topic 檔。trim 前先驗 detail 在下層才動，故無損。

## B. 可重用套件（reusable packages）

**唯一新建（restraint 通過）**：`scripts/hooks/script_lint_advisor.sh`
- PostToolUse(Edit|Write) → `.py` py_compile / `.sh` bash -n（+shellcheck if available）；**advisory exit 0 + fail-OPEN**；scope 限 `scripts/`。
- WHY：C++ 有編譯 Hard Gate，但 .py/.sh 零 lint；harness 自家治理 helper（task_graph/number_provenance/harness_health/builder）全未 lint = insights friction #4「自釀 bug」根因。
- 自我測試 4/4 通過（好檔靜默 / 壞 .py 警告+exit0 / scope 外 skip / 壞 .sh 警告）。
- ⚠ wiring 已加 settings.local.json **live**，但 settings 有未提交 churn（§4 既載），**暫不 commit settings**；hook script 本體已 commit。

**restraint 排除（已覆蓋，不重造）**：verified-deck→/html-report-build+/verify-workstation｜confound skill→/auc-confound-guard｜fan-out hook→§8 routing-class｜completion-gate/skill-registry memory→§13.7 家族+harness_health 燈#1｜數字溯源/CJK/anti-overclaim→既有 §13+feedback。

## C. 流程控管（process）

| 項 | 動作 | 驗證 |
|---|---|---|
| **builder-vs-helper 雙寫** | `_build_full_tree.py` 加 guard：graph.json 已存在且無 `--force` → refuse exit 2。宣告 graph.json canonical（由 task_graph.py helper 維護），builder 為一次性 migration | 無 --force refuse ✓ / --force 放行 ✓ |
| active.json stale | G6/G1（06-02 backfill placeholder，hypothesis_id=null）→ 移 `recently_concluded`（verdict=superseded_by_subclonal_axis），cycles 清空；**cycle 目錄未刪** | JSON valid；**T-C1 DRIFT 已清** |

> active.json 手改的正當性：audit 建議「manual + user confirm（非 hook）」+ 用戶「全部改進」=confirm + git-tracked 可回復 + §6 明文這是 placeholder。schema 保持有效（move between arrays，非破壞）。**不加 auto-archiver hook**（會重蹈「假裝多週主軸=live cycle」反模式）。

## 紅線 / 不做（restraint 紀錄）

- ❌ output-token main-turn hook（無可靠 per-turn signal）／fan-out estimator hook（§8 routing-class，一次自我修正事故不到門檻）／number_provenance bare-integer（假陽性風險，metric-shaped 是刻意 tradeoff）。
- ❌ active.json auto-cleanup hook（反模式）。

## 可回復性

全部 repo 變更 git-tracked（active.json / builder / lint hook / 本檔）；memory 變更在 `.claude/projects/` 持久層（review-format B 檔保留為 archived provenance）。settings.local.json live-only 未 commit。
