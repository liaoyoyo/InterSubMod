---
name: run-evaluator
description: Phase 5 EVALUATE gate of the Resilient Waterfall harness — compute retraction risk score and tier recommendation from a cycle's complete artifacts. Reads `state/cycles/{cycle_id}/{state,plan,precheck,pilot,generalize}.json`, sweeps known-pitfalls, computes 6-component risk + composite + per-component override, writes `evaluation.json`. MANDATORY before tier ⭐4 / ⭐5 upgrade. USE WHEN：「run evaluator」「P5 evaluate」「retraction risk」「tier 升級前」「evaluator」、cycle 從 P4 → P5 transition 時。
allowed-tools: ["Bash", "Read", "Glob", "Grep"]
user-invocable: true
tags: ["harness", "evaluate", "retraction-risk", "tier-gate"]
---

# /run-evaluator — Phase 5 EVALUATE Gate

對一個已完成 P4 GENERALIZE 的 cycle 進行 retraction-risk 評估，是 Resilient Waterfall harness 的 P5 gate。產出的 `evaluation.json` 是 tier ⭐4 / ⭐5 升級的**必要前置**（由 `pre_tier_upgrade_check.sh` hook 強制）。

> **設計來源**：`~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md` §3.5, §4.5.1

## 何時使用

- cycle 進入 P5 EVALUATE 階段（state.json.phase 從 `P4_GENERALIZE` → `P5_EVALUATE`）
- 用戶手動跑：`/run-evaluator <cycle_id>`
- ⭐4 / ⭐5 升級**前**強制呼叫（即使你覺得「結果很穩」也必須跑）
- 全自動模式（headless / `auto`）：不論 tier 都跑

## 何時 NOT 使用

- cycle 還在 P3 PILOT — 等 P4 GENERALIZE 完成再跑
- 純 ⭐1-3 描述性結論且非全自動模式 — user-on-demand
- 已 retracted 的 archived cycle — evaluator 不能逆轉撤回

## 工作流程（5 步）

### Step 1 — 讀 cycle 全部 artifacts
```bash
CYCLE_ID="$1"
CYCLE_DIR="InterSubMod/state/cycles/${CYCLE_ID}"
test -f "${CYCLE_DIR}/state.json" || { echo "ERROR: cycle not found"; exit 2; }
```
必讀：`state.json`, `plan.json`, `pilot.json`, `generalize.json`
建議讀：`precheck.json`（若不存在仍可繼續，但 precondition_freshness 算 0.5 中性）
可選讀：`research/autoresearch/evidence_ledger.jsonl`（為 tier_support_alignment 計分）

### Step 2 — 計算 6 個 risk components

| # | Component | 來源欄位 | 公式 |
|---|---|---|---|
| 1 | multi_sample_consistency | `generalize.consistency.n_samples_passed / n_samples_total` | 0–1 ratio |
| 2 | effect_size_stability | generalize.samples[].metric_value | min(\|v\|) / max(\|v\|), if max>0 else 0 |
| 3 | precondition_freshness | precheck.verdict | PASS=1.0 / WARN=0.7 / BLOCKED=0.3 / missing=0.5 |
| 4 | subgroup_homogeneity | pilot.metric_results 跨 subgroup stddev/mean | clip(1 - stddev/mean, 0, 1) |
| 5 | pitfall_coverage_score | pitfalls_table.json sweep | 1 - n_unchecked_relevant / n_relevant |
| 6 | tier_support_alignment | ledger stability vs tier_used | min(stability, tier) / tier |

7th component (`precedent_similarity`) 為 Path B（LlamaIndex）才填，Path A 設 null 跳過權重重分配。

### Step 3 — Composite formula + 兩階段 verdict

```python
weights = {1: 0.27, 2: 0.22, 3: 0.17, 4: 0.17, 5: 0.17}
risk_base = sum(weights[i] * (1 - components[i]) for i in 1..5)

# Stage 1: base verdict
if risk_base > 0.7:    base_verdict = "pending_review"
elif risk_base > 0.4:  base_verdict = "downgrade_tier"
else:                  base_verdict = "approve_tier"

# Stage 2: 規類強制 override（避免單點崩潰被平均掩蓋）
n_low = count(c < 0.4 for c in components)
n_critical = count(c < 0.2 for c in components)

if n_critical >= 1 or n_low >= 3:
    final_verdict = "pending_review"
elif n_low >= 1:
    final_verdict = max(base_verdict, "downgrade_tier")
else:
    final_verdict = base_verdict
```

### Step 4 — Pitfalls sweep
讀 `.claude/skills/run-evaluator/pitfalls_table.json`（6 條 P-01..P-06）。對每條：
- 比對 `trigger_keywords` 與 plan.hypothesis（substring match, case-insensitive）
- 若命中：套用 `auto_check_rule` 評估 severity
- 收集為 `pitfall_hits[]`，每項 {pitfall_id, severity, evidence}

### Step 5 — 寫 evaluation.json
依 `state/schemas/evaluation.schema.json` 格式寫到 `state/cycles/{cycle_id}/evaluation.json`，含：
- `retraction_risk` (= risk_base)
- `risk_components` (6 個 component 數值；precedent_similarity = null)
- `pitfall_hits[]`
- `tier_recommendation` (依 final_verdict 轉換：approve→keep current, downgrade→tier-1, pending→null)
- `verdict` (final_verdict)
- `required_followups[]`（依 verdict 自動填）

## 用戶 override 路徑

用戶不滿意 evaluator 結論時，**不直接編輯 evaluation.json**，而是：
1. 用戶提供 reviewer name + decision_mode (`override_higher` / `override_lower` / `defer`) + reason
2. skill 在 evaluation.json 寫入 `human_reviewer` 物件（schema 已預留）
3. `pre_tier_upgrade_check.sh` hook 看到 human_reviewer.decision_mode 為 override 時，依其 override_tier 而非 tier_recommendation 判定

## 輸出格式（人類可讀摘要）

```
[/run-evaluator] cycle_id=20260504-1430-loh-kde
  Components:
    [1] multi_sample_consistency:  0.86 ✓
    [2] effect_size_stability:     0.72 ✓
    [3] precondition_freshness:    1.00 ✓
    [4] subgroup_homogeneity:      0.34 ⚠ (low)
    [5] pitfall_coverage_score:    0.83 ✓
  composite risk_base: 0.215
  base_verdict: approve_tier
  override: 1 component < 0.4 → escalate to downgrade_tier
  final_verdict: downgrade_tier
  pitfall_hits: 0
  tier_recommendation: 3 (was plan.tier=4)

Written to: state/cycles/20260504-1430-loh-kde/evaluation.json
```

## 與其他 skills 整合

| 場景 | 此 skill 做什麼 | 既有 skill 接手 |
|---|---|---|
| risk > 0.7 (pending_review) | 寫入 evaluation.json + verdict=pending_review | 用戶呼叫 `/grill-me` 或 P4 重做 |
| 1 component < 0.2 critical | escalate to pending_review | 用戶查 `/known-pitfalls` 判斷是否該 component 真的崩了 |
| pitfall_hits 含 P-04 (pileup symlink) | 寫入 hit + severity=block | 用戶必須回 P2 PRECHECK 重驗 dataset |
| 全部 components > 0.6 | approve_tier | 用戶可呼叫 `/conclude-research` 進 P6 |

## 限制 (Path A)

- 不接 LlamaIndex precedent retrieval（component 7 為 null）
- `tier_support_alignment` 仰賴 evidence_ledger stability 欄位是字串（如「4 (HCC1395 pilot, awaits cross-sample)」）需做 grade extraction
- pitfalls_table.json 7 條規則（P-01..P-07）與 known-pitfalls SKILL.md 同步維護（**新增 pitfall 時需手動同步兩處**）
- 不重跑 statistics — 只消化 cycle 已寫入的 artifacts
- **4 軌證據鏈僅部分自動化**（個人風格 anchor #1 — 詳見 validation-protocol L4 mandatory 區段）：
  - (i) Statistical：人工填 4-track table
  - (ii) Cross-sample：`multi_sample_consistency` 是此軌 surrogate（自動）
  - (iii) Mechanism：人工填
  - (iv) Orthogonal：`pitfall_coverage_score` + P-07 偵測 single-track artifact（自動）
  - **Path B 將新增 component 7 `multi_track_corroboration`** 自動掃 4 軌 artifact reference 覆蓋率
- ⭐4/⭐5 升級的 hook gate (`pre_tier_upgrade_check.sh`) 目前**只檢查 evaluation.json 存在**，不驗 4-track table 完整度；validation-protocol L4 要求 PR 附表，缺軌 reviewer 應拒升

## 維護規則

當 known-pitfalls SKILL.md 有新增 P-XX 條目時，必須**同步更新** pitfalls_table.json：
1. 加 `pitfalls[]` 條目
2. 更新 `last_synced_with_skill` 日期
3. 不可只更新 SKILL.md 不更新 JSON（否則 evaluator 看不到新 pitfall）

維護備忘：`pitfalls_table.json` 是給機器讀；`known-pitfalls/SKILL.md` 是給人讀。內容語意必須一致。

**Path B 候選增強清單**（v1.6 plan §4.5.4-G 5b 候選；未實作）：
- Component 7 `multi_track_corroboration`：掃 cycle 的 plan/pilot/generalize/structured-tech-report 中 4 軌引用，回 0-1 覆蓋率（4/4 → 1.0；3/4 → 0.75 等）
- Component 7 觸發條件：cycle 申請升 ⭐4/⭐5 時必算；⭐3 cycle 可選
- 接入 LlamaIndex precedent retrieval 後，`precedent_similarity` 與 `multi_track_corroboration` 都取代部分人工 reviewer 工作

## DO NOT USE WHEN（v1.7 batch A 採納）

- **cycle 還在 P3 PILOT 或 P4 GENERALIZE** — 等 generalize.json 寫入 P5 phase 後才用；過早跑 evaluator → 多 component 為 None → risk_base 偏低不可信
- **cycle 沒 plan.json 或沒 pilot.json** — 必要 artifact 缺；evaluator 會 silent skip components → 結果無意義
- **想計算 retraction risk 給「未啟動 cycle」評分** — evaluator 不是 a priori 風險評分工具，是 P5 gate
- **想做跨 cycle 全域審計** — 用 `/provenance-tier-audit` 而非本 skill；evaluator 只看單 cycle
- **想要 user-friendly tier upgrade ceremony** — 本 skill 純 mechanical risk score；tier 升級流程交給 conclude-research / weekly-report
- **plan.json 是手寫 mock** — schema 校驗不嚴格，但 risk 計算依賴實際 metric_value 數字；mock 數會誤導

## Quality Checklist — 交付 evaluation.json 前自我檢查（v1.7 batch B）

執行 main() 結束、給用戶結果之前，跑過這 8 條：

- [ ] 5 個 risk_components 都有非 null 值（precedent_similarity 預期 null 為 Path A 正常）
- [ ] composite risk_base 與 verdict 對齊：< 0.4 → approve / 0.4-0.7 → downgrade / > 0.7 → pending_review
- [ ] per-component override 邏輯生效：若任一 component < 0.2 critical → 強制 pending_review；若 ≥3 component < 0.4 → 強制 pending_review
- [ ] tier_recommendation 與 verdict 一致：approve → 不變 / downgrade → tier-1 / pending → null
- [ ] failure_attribution 寫入規則：verdict ≠ approve_tier → 必填；verdict = approve_tier → 不應有此欄位
- [ ] state.json 同步更新（denormalized snapshot mirror）
- [ ] evaluation.json 寫入路徑正確（`state/cycles/{id}/` 或 `state/retro_cycles/{id}/`）
- [ ] pitfall_hits 中每個 entry 含 pitfall + severity + evidence 三欄

## Failure Mode 排查（既有，整理）

| 症狀 | 可能原因 | 排查 |
|---|---|---|
| risk_base 永遠 < 0.05 | 所有 component 全 1.0 → cycle 太完美？或 generalize.consistency 沒讀對 | 印 components dict 看是否有 None / 1.0 |
| pending_review 但 components 都 ≥ 0.4 | per-component override n_low ≥ 3 觸發 | 確認 override 規則照 plan §4.5.1 |
| failure_attribution 永遠 unknown | components 全 null（pilot/generalize 缺）| 先確保 P3-P4 完成 |
| state.json 未鏡寫 | update_state_failure_attribution 例外吞掉 | 加 logging
