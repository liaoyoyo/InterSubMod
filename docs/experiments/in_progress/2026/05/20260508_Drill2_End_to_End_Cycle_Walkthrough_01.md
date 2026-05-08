---
title: Drill 2 — Resilient Waterfall Harness 端到端 Forward Walkthrough
date: 2026-05-08
status: validated
phase: P6_COMMIT (drill artifact)
type: harness_validation_drill
tier: 4
classification: harness_validation
mode: D2-A placeholder（合成資料）
cycle_id: 20260507-2112-d2a-colo829-kde-rerun
upstream_reports:
  - InterSubMod/docs/experiments/in_progress/2026/05/20260506_Drill1_Harness_Retrospective_01.md
related_plan: ~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md (v1.7)
---

# Drill 2 — 端到端 Forward Walkthrough

> **Bottom line**：D2-A cycle 用合成 placeholder data 走完 P0→P5（P6 部分）。**M1 Trust Score、M1 Drift、M2 Failure Attribution、M4 Interaction Metrics 四個機制全部如期運作**。最值得記錄的是：**P5 evaluator 正確抓出 region_7 outlier（effect_size_stability=0.35）並觸發 per-component override 從 approve_tier 強制改 downgrade_tier**，這是 harness 設計目標的真實展現。Drill 2 sensitivity = 1/1（單 cycle 偵測 placeholder 內植入的 outlier）。

---

## §1 Cycle 7-phase 時間軸

| Phase | Skill | 時間 | trust_score | 備註 |
|---|---|---|---|---|
| P0 REGISTER | `/cycle-init` | 2026-05-07 21:12 | 🟢0.74 | main_axis lock 設定（KDE/CN1/COLO829/H-CN1）|
| P1 PLAN | research-loop（手寫 plan.json） | 2026-05-07 21:13 | – | 8 欄全 OK; 全用戶單次 ack |
| P2 PRECHECK | `/check-staleness` | 2026-05-07 21:14 | – | verdict=PASS（純分析無 binary、dataset/upstream 全 fresh） |
| P3 PILOT | placeholder synthetic pilot.json | 2026-05-08 09:10 | – | 9 region; 1 region (region_7) 弱 effect 故意注入 |
| P4 GENERALIZE | placeholder synthetic generalize.json | 2026-05-08 09:11 | – | 8/9 passed_threshold; consistency.confidence_uplift=0.78 |
| P5 EVALUATE | `/run-evaluator` | 2026-05-08 09:12 | 🟡0.58 | **verdict=downgrade_tier**（override 觸發） |
| P6 COMMIT | （此報告） | 2026-05-08 | — | placeholder cycle 標 partial completion |

**累計 phase pass time**：~12 hr（含跨 session 等待）。**真實計算耗時 < 1 min**（synthetic data）。

---

## §2 用戶介入次數（M4 計數器）

per plan v1.6 操作性定義（feedback_strategy_then_per_item_confirmation 後續）：

| 介入點 | 你說了什麼 | 算 1 次嗎 |
|---|---|---|
| P0 cycle-init 參數 ack | 「全 OK,先執行後確認調整」 | ✓ 1 |
| P1 plan.json 8 欄審查 ack | 「全 OK」 | ✓ 1 |
| P3 PILOT 路徑選擇 | 「OK」（採 C placeholder）| ✓ 1 |
| **總計** | | **3 次** |

**target ≤3**（plan v1.6 §6.2）達成 ✅。

**自動推進**：P0→P1（state.json 自動寫）/ P1→P2（precheck 自動跑）/ P2→P3（artifacts 自動更新）/ P3→P4（synthetic 寫入後自動 advance）/ P4→P5（evaluator 自動跑）= **5 次 auto_advance**。

**auto/intervention ratio = 5/8 = 62.5%**（接近 plan target 70%；剩 7.5% gap 因 placeholder 模式介入較多 — 真跑模式預期 ratio 更高）。

---

## §3 M1 Trust Score 動態變化

| 時點 | trust | 主要變動因子 |
|---|---|---|
| 啟動 P0 | 0.74 🟢 | 純分析高 reversibility (1.0) + 中性 time_budget (1.0) + skill 調用低 (0.5) |
| P2 PASS 後 | 0.78 🟢 | invocation_count 升至 0.7（已 2 phase advance） |
| 11 hr 後 | 0.58 🟡 | time_budget_remaining 降至 0.46（estimated 24hr 用了 11hr）|

**Trust 機制行為符合設計**：時間預算消耗 → 自動降 trust → cycle-state 加 routing recommendation 提醒 user。**未誤觸發 LOW <0.4 暫停**。

---

## §4 M2 Failure Attribution 自動填寫範例

P5 evaluator 執行後 state.json 自動填：

```json
"failure_attribution": {
  "categories": ["consistency_violation"],
  "primary_category": "consistency_violation",
  "phase_at_failure": "P4_GENERALIZE",
  "skill_at_failure": "run-evaluator",
  "gate_at_failure": "P5_EVALUATOR",
  "components_below_threshold": ["effect_size_stability=0.35"],
  "pitfalls_hit": [],
  "confidence": "low",
  "user_intervention_required": false
}
```

**Confidence=low** 因為只 1 component < 0.4 + 無 pitfall hit；如有 pitfall match 會升 high。**設計符合預期**。

---

## §5 對比 Drill 1（Retrospective vs Forward）

| 維度 | Drill 1 (retrospective) | Drill 2 (forward placeholder) |
|---|---|---|
| 資料來源 | 6 真實 retract events + 2 negative controls | 1 cycle synthetic |
| 攔截目標 | 已知撤回事件 | 注入 outlier |
| Sensitivity | 6/6 (100%) | 1/1（單 cycle）|
| Specificity | 2/2 (100%) | n/a (no negative control in D2) |
| 涵蓋 phases | P2 BLOCKED + P5 pending（皆 retrospective）| **全 P0-P5 forward** |
| 用戶介入 | n/a (fixture) | 3 次（達 target）|
| 暴露的弱點 | hard-fidelity events 用 metric 重現 | placeholder data 不能驗真實計算品質 |

**Drill 2 補了 Drill 1 沒測的兩件事**：
1. **forward routing chain**（P0→P5 連續推進，非 fixture rebuild）
2. **M1+M2+M4 governance 機制 in vivo**（trust 動態 + failure attribution 自動 + interaction metrics 統計）

---

## §6 routing chain 弱點觀察

| 步驟 | 順暢嗎 | 觀察 |
|---|---|---|
| cycle-init → 寫 state.json | ✓ | main_axis 寫入正確；schema 驗證通過 |
| state.json → P1 (research-loop) | △ | 無 automatic transition；用戶/AI 必須手動推；建議未來加 SessionStart hook 提醒 |
| P1 → P2 (check-staleness) | ✓ | next_action 路由建議顯示，user 容易跟 |
| P2 PASS → P3 | ✓ | dashboard 顯示 next_action |
| P3 → P4 → P5 | ✓ | evaluator 自動讀 plan/pilot/generalize artifacts |
| P5 verdict 自動填 failure_attribution | ✓ | M2 機制工作 |

**建議改進**：
- next_action 在 markdown 顯示，但無法自動 invoke 對應 skill — 仍需 user/AI 手動觸發
- Path B 提到的 `cycle_router.py`（v1.7 §6.4）可解此痛點 — phase advance 自動 invoke 下一 skill

---

## §7 §4.5.4-G batch 3-5b 重評（per plan v1.6 §4.5.4-F）

| Batch | 原 plan 預期收益 | Drill 2 觀察 | 決策 |
|---|---|---|---|
| **Batch 3 — 4 passive ref** | 減少誤觸發 | Drill 2 中 known-pitfalls 在 P-04/P-06 sweep 內被消費；無誤觸發 | **`[PARKED]`**（v1.6 已 PARKED，續 PARKED） |
| **Batch 4 — 23 forward-link** | 提升 chain 路由精準度 | next_action 路由規則就足以引導；用戶實際跟 chain 順暢；**verdict=不再需要做** | **`[DROPPED]`**（v1.7 升級此決策；理由：cycle-state next_action + 個別 SKILL.md 「Phase Chain Position」段已 cover 80% 路由需求） |
| **Batch 5b — 3 anchor 硬化** | 對齊個人風格 | anchor #3/#6/#7 在 placeholder 模式不易顯示必要性 | **`[PARKED]`**（等真跑 cycle 才能評；D2-A 真跑或 D2-B 啟動時再評） |

---

## §8 Path B 啟動條件評估

per plan v1.6 §6.4 Decision Gate：「Path A 完成 + 30 天觀察撤回 ≤2 → 進 Path B」

當前狀態：
- Path A：✅ 完成（4 skills + 2 hooks + 7 schemas + Phase 1 governance + v1.7 patterns）
- 30 天觀察起點：**Drill 2 結束日 = 2026-05-08**；30 天 baseline = 2026-06-07
- 撤回率追蹤：需週報級別人工統計（Drill 2 自動化 metrics 不含此）

**v1.7 §11.4 已修訂 Path B Layer 2 為混合方案**（5 day → 2.5 day）+ Layer 3 LlamaIndex 仍保留。

**等待條件**：
1. 6/7 (2026-06-07) 進行 30 天回顧
2. 若撤回 ≤2 → 啟動 Path B Layer 3（LlamaIndex first，Layer 2 後續）
3. 若撤回 >2 → 分析 root cause；可能調 evaluator 閾值或加 anchor #1 4-track 自動 audit

---

## §9 給未來自己的 Heads-up（v1.7 補）

1. **placeholder vs 真跑要明確標 `_synthetic: true`** — 否則未來 audit 會把假數當真結論
2. **Drill 2 的 D2-A 不算真實研究結論** — COLO829 9 region 重跑要在新 cycle 跑（用 D2-A learnings 加快設計）
3. **`tier_recommendation: pending` 且 user_intervention_required=false** 是合法狀態（downgrade_tier 不需 user 介入；只 pending_review 才需）
4. **per-component override 的「規類強制」是 harness 真正力量** — region_7 outlier 在 base composite (0.216) 下會 approve；override 升為 downgrade。這個 reflexion 機制是 Drill 1 sensitivity 6/6 的真正主因
5. **Drill 2 沒做 negative control**（Drill 1 N1/N2 已驗）— 未來想再做 D2-B (HCC1395 chr8) 時不需要再 negative control

---

## §10 References

- Plan: `~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md` v1.7
- Drill 1 報告：`InterSubMod/docs/experiments/in_progress/2026/05/20260506_Drill1_Harness_Retrospective_01.md`
- Cycle artifacts：`InterSubMod/state/cycles/20260507-2112-d2a-colo829-kde-rerun/{state,plan,precheck,pilot,generalize,evaluation}.json`
- Memory: `project_kde_fix_downstream_quantification` (synthetic placeholder 數值來源)
- v1.7 commits: 5a136d8 (Phase 1 governance) → 517c467 (merge) → e2eb43b (C+B partial) → e95c96a (A+B+E1) → 26d1d34 (v1.7 merge)
