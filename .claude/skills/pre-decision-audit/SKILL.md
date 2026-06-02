---
name: pre-decision-audit
description: |
  快速 (≤30min) pre-decision evidence audit — 在進入 cycle / spec / 大改之前做
  §0 Cynefin front-gate (domain classification) + §1-§7 七個 core outputs：
  §1 觀察完整性 checklist L1-L5 ｜ §2 假設可信度評分 5 維 0-100 (含 Pre-Mortem falsifier observable + WRAP Reality-test 三反例)
  §3 Assumption Map 2x2 Bland (importance × known) ｜ §4 Quick pilot guide <30min + checkpoint
  §5 Gap diagnosis (缺什麼+影響+工時+優先級) ｜ §6 Evidence conflict scan (cross-link 既有 NEGATIVE)
  §7 Decision path GO / PROBE / NO-GO + 後續 action.
  與 /implementation-notes (process-time) + /run-evaluator (post) 形成「假說驗證三層樓」之首層；
  對齊 ADR pre-decision separation 原則 (exploration BEFORE decision).
  USE WHEN: cycle 啟動前 (P0), 新 spec 實作前, 跨樣本推廣前 (P3->P4), tier 升級前 (⭐3->⭐4),
  用戶說「實作 <X>」「研究 <Y>」「探索 <Z>」「考慮 <Z>」「我想 <W>」, NO-GO 判定前.
  SKIP WHEN: 純 code edit / typo / build / commit (無 research decision),
  假說已在 active cycle 驗證中 (audit 只在 entry 點), tier ≤ ⭐2 純探索 (過早 audit 浪費),
  retrospective summary (用 /report), 已簽 pre-registration (框架已鎖無須前驗),
  trivial 1-line config 變更.
allowed-tools: Read, Write, Edit, Glob, Grep, Bash
user-invocable: true
paths:
  - "research/**/pre-decision-audit.md"
  - "research/**/audit.md"
  - "state/cycles/**/audit.json"
  - "docs/reports/validated/**/*pre_decision_audit*.md"
---

# Pre-Decision Audit Skill

> **快速 (≤30min) 決策前審計** — 與 `/implementation-notes` (process-time) + `/run-evaluator` (post) 形成「假說驗證三層樓」之首層

## §1 Phase & Chain Position

**Cross-cutting**：P0 / P1 / P3 / P4 / P5 邊界（在 ADR / implementation-notes / cycle-init **之前**呼叫）。

**對齊 ADR pre-decision separation 原則**：exploration BEFORE decision；spike / benchmark / prototype 必在 ADR 寫入前完成。

### 強制觸發點

| 節點 | 觸發強度 | 原因 |
|---|---|---|
| **P0 `/cycle-init`** 註冊前 | 🔴 必觸 | 新 cycle 啟動前驗證假設站不站得住 |
| **P1 `/research-loop`** Step 1 ORIENT 前 | 🔴 必觸 | 定向前審查證據基礎 |
| **P3->P4 跨樣本推廣** 前 | 🔴 Gate | 從 1 樣本擴 N 樣本需先 audit |
| **⭐3 → ⭐4 tier 升級** 前 | 🔴 Gate | tier 升級門檻必跑 |
| **Hard NO-GO 判定** 前 | 🔴 必觸 | 補 §6 conflict scan 證據 |
| **`/inject-hypothesis`** 新 H 註冊前 | 🟡 推薦 | 避免無效假說入隊 |
| **`/cpp-change`** C++ 決策前 | 🟡 推薦 | 大改 C++ 前評估 |
| P2 / P6 | 🟢 可選 | freshness / conclude phase 非關鍵 |

**頻率上限**：1 次 / hour / cycle（避免 nagging）；SessionStart 1 次 / session。

## §2 Dependencies

### Uses（本 skill 引用）

- `/scientific-rigor` §2-§7 — L1-L5 evidence tier + DAG + pre-registration（每項觀察 + 假設標 tier）
- `/confirmation-protocol §Cynefin` — 5 域分類（§0 front-gate 共用 SoT，不重複定義）
- `/research-context-loader` Tier 1/2/3 — 既有報告 context 載入
- `/known-pitfalls` P-01-P-14 — 陷阱交叉檢核（§6 conflict scan）
- `/doc-standards` — filename + tier ribbon 規範

### Used by（誰呼叫本 skill）

- `/cycle-init`（P0）— state.json 寫入前
- `/research-loop`（P1）— Step 1 ORIENT 前
- `/inject-hypothesis` — 新 H append 前快速 audit
- `/cpp-change` — 第 1 步「方案評估」前
- `/implementation-notes` — `init` 時讀取 audit assumptions 填入 §🟡 折衷
- `/run-evaluator`（P5）— tier 升級前 reference §2 credibility score

### Reads

- 既有 reports: `InterSubMod/docs/reports/validated/**/*.md`
- MEMORY index: `MEMORY.md`（過往 NEGATIVE / Concluded section）
- Topic research dir: `InterSubMod/research/<topic>/*.md`
- Active cycle state: `InterSubMod/state/cycles/<cycle_id>/`

### Writes

- Markdown source: `InterSubMod/research/<topic>/pre-decision-audit.md`
- Machine-readable verdict: `InterSubMod/state/cycles/<cycle_id>/audit.json`
- HTML companion (via `standalone_trigger.sh` hook): `pre-decision-audit.standalone.html`
- Finalized 移檔: `InterSubMod/docs/reports/validated/{YYYY}/{MM}/{YYYYMMDD}_<topic>_pre_decision_audit_01.md`

## §3 Three sub-commands

### `init <topic>`

1. 建 dir `InterSubMod/research/<topic>/`（若不存在）
2. 從 `InterSubMod/templates/pre-decision-audit-template.md` 拷貝到 `<dir>/pre-decision-audit.md`
3. 填 frontmatter: `cycle_id` / `topic` / `status=in_progress` / `audit_version=0.1` / `建立時間`
4. 提醒用戶逐 section 填（§0 → §7）
5. 預期 ≤30min 完成

### `score`

1. 讀 `<dir>/pre-decision-audit.md` §2 Credibility Score Table
2. 計算 TOTAL（5 維加總 0-100）
3. 寫入 §7 「TOTAL: __ / 100」
4. 同步寫 `state/cycles/<cycle_id>/audit.json`:
   ```json
   {
     "cycle_id": "...",
     "topic": "...",
     "score": 65,
     "dimensions": {"theory": 20, "observation": 15, "mechanism": 10, "counter_risk": 10, "resource": 10},
     "falsifier_observable": "...",
     "scored_at": "ISO8601"
   }
   ```

### `verdict`

1. 讀 `audit.json` score
2. 套門檻：
   - ≥ 60 → **GO** (status=verdict_GO)
   - 30-60 → **PROBE** (status=verdict_PROBE) + 必執行 §4 pilot
   - < 30 → **NO-GO** (status=verdict_NO-GO)
3. 寫 §7「Next action」：
   - GO → 建議用戶 run `/cycle-init` + `/implementation-notes init <topic>`
   - PROBE → 列 §4 pilot 步驟 + checkpoint
   - NO-GO → 建議寫 MEMORY.md Concluded entry / postmortem（若 mid-cycle）

### Red-team gate（2026-06-02 借鑑 Google co-scientist Reflection agent — 正名 G7 devils-advocate）

**GO verdict 前必跑一次獨立紅隊**（不可與 `score` 同一思路 — generator/evaluator 分離精神）：
1. **最強反方**：「這假說最可能怎麼錯？」列 ≥2 failure mode（confound / circularity / 已 concluded-dead 變體 / 樣本量不足 / collider）。
2. **對照 MEMORY Concluded**：與已 NEGATIVE 方向重疊？重疊 → 需 Productive-Failure reopen 3 條件（C1 新數據 / C2 新方法 / C3 新前置）至少一項，否則 **NO-GO**。
3. **falsifier 檢查**：§7 falsifier_observable 真能否證？不能否證 → 降 **PROBE** 不 GO。

紅隊未過 → verdict 降一級（GO→PROBE / PROBE→NO-GO）。與 4 個 post-數據 read-only verifier **正交**（這是 pre-cycle 假說層紅隊）。

## §4 §0 Front-gate + §1-§7 Seven core output sections

### §0 🟤 Cynefin Domain Gate (front-gate)

決定後續 audit 走 best-practice 還是 probe-first。

- Test：「相同行動是否曾重複產生**可預測**結果？」
  - Yes → **Clear / Complicated** → 可用 best-practice checklist
  - Maybe → **Complex** → **強制 probe-first**（禁套 best-practice）
  - No → **Chaotic** → 先穩定（行動 → 感知 → 反應）再 audit
- 共用 `/confirmation-protocol §Cynefin` 定義（單一 SoT）

### §1 🟢 Observation Completeness Checklist

| Observation | 狀態 (✓/□/✗) | Tier (L1-L5) | 來源 file:line |

- ✓ 已驗 ｜ □ 待補 ｜ ✗ 反例
- 每項標 L1（強實證）→ L5（純假說）evidence tier ribbon
- 對應 `/scientific-rigor §2`

### §2 🔵 Credibility Score (5-dim, 0-100)

5 維 × 0/10/20 → TOTAL 0-100：

1. **理論基礎**：0 speculative / 10 plausible / 20 well-grounded
2. **觀察支撐**：0 none / 10 partial / 20 complete（≥3 樣本一致）
3. **機制清晰度**：0 speculative / 10 plausible / 20 clear DAG
4. **反例風險**：0 high / 10 medium / 20 low（≥3 預期反例皆無）
5. **所需資源**：0 >6h / 10 1-6h / 20 <1h pilot

**Falsifier observable (WRAP)**：「若假設錯，我會看到什麼？」**必填**（避免 unfalsifiable claim）。

### §3 🟡 Assumption Map 2×2 (Bland)

```
                | known            | unknown
HIGH importance | (1) verify quick | (2) MUST validate FIRST ⚠
LOW importance  | (3) document     | (4) defer / ignore
```

- 右上象限（HIGH × UNKNOWN）= **MUST validate before proceed**
- 避免「重要但未知」假設被誤當「已知」推進

### §4 🟠 Quick Pilot Guide (<30min)

3-step minimum：

1. **Step 1**: Read `<file:line>`（context）
2. **Step 2**: Run `<script>` with 1 sample（initial signal）
3. **Step 3**: Check `<metric>` 與 `<threshold>`

**Checkpoint**:
- ≥ threshold → **GO** full verify
- 0.5×threshold ~ threshold → **PROBE** 再 1 樣本
- < 0.5×threshold → **NO-GO** archive

### §5 🔴 Gap Diagnosis

| Missing | Impact | Effort | Priority |

精確告知「為何無法現在決策」 + 缺哪些 + 工時估算 + 排優先級（P0/P1/P2）。

### §6 🟣 Evidence Conflict Scan

grep + cross-link 既有 NEGATIVE / Concluded / dependent reports：

| Prior conclusion | Tier | Relation (support/conflict/dependent) | Source |

**強制**：
- `grep -i "concluded\|negative" InterSubMod/MEMORY.md` → 列所有相關 entries
- `find InterSubMod/docs/reports/validated -name "*NEGATIVE*"` → 列檔案
- 若 conflict 找到 ≥ 1 → §2 反例風險維度自動 −10（如 10→0 或 20→10）

### §7 ⚫ Decision Threshold + Path

- **TOTAL** from §2
- **Verdict**: GO / PROBE / NO-GO
- **Next action**: 對應 `/cycle-init` / pilot / postmortem
- **Decision lock**: Y/N（若 Y，後續不可任意 reopen 除 Productive Failure 3 條件 C1/C2/C3 至少一項，對應 `feedback_productive_failure_reopen_threshold`）

## §5 業界對齊（5 個權威來源）

| 來源 | 用法 |
|---|---|
| [Pre-Mortem (Klein HBR 2007)](https://hbr.org/2007/09/performing-a-project-premortem) | §2 **Falsifier observable** — "假定失敗反推" prospective hindsight |
| [Assumption Mapping (Bland Precoil)](https://www.precoil.com/assumptions-mapping) | §3 2×2 importance × known/unknown — 右上 quadrant MUST validate |
| [WRAP (Heath ModelThinkers)](https://modelthinkers.com/mental-model/wrap-decision-process) | §2 **Reality-test 三反例觀察** — Heath "R" step（"列若假設成立應看到/不該看到"），與 Pre-Mortem falsifier 互補 |
| [Cynefin (Snowden Untools)](https://untools.co/cynefin-framework/) | §0 front-gate domain classification — Complex 域強制 probe-first |
| [ADR pre-decision (AWS)](https://docs.aws.amazon.com/prescriptive-guidance/latest/architectural-decision-records/adr-process.html) | 整 skill 定位 — exploration BEFORE decision separation |

額外 reference：
- [Microsoft Well-Architected ADR (含 confidence level)](https://learn.microsoft.com/en-us/azure/well-architected/architect-role/architecture-decision-record)
- [Building Effective AI Agents — Anthropic](https://www.anthropic.com/research/building-effective-agents)
- [Preregistration — OSF/AsPredicted](https://en.wikipedia.org/wiki/Preregistration_(science))

## §6 與既有 skill 分工

| Skill | 時機 | 形式 | Trigger | Verdict |
|---|---|---|---|---|
| **`/pre-decision-audit`** (本) | 進入 cycle / spec **前** | process-time pre-verdict | 「實作 X」/「研究 Y」/ P0/P1/P3/P4 邊界 | GO / PROBE / NO-GO |
| `/implementation-notes` | spec 實作**中** | process-time live tracking | 4 markers / AI advisory | 4 sections live |
| `/scientific-rigor` | 完整方法論深驗 | rigorous post-evidence | 元方法論 | L1-L5 tier |
| `/run-evaluator` (P5) | cycle 結束**後** | rigorous post-evaluation | tier 升級 | ⭐1-5 + risk components |
| `/report` | session-end | retrospective summary | 「寫報告」 | summary |

**核心切割**：
- vs `/scientific-rigor` — 前者 ≤30min 快速前驗，後者 ≥1h 深驗；前後互補
- vs `/implementation-notes` — 前者 entry 點 verdict，後者 in-progress live；時序串接
- vs `/confirmation-protocol` — Cynefin 共用，但前者額外做 7 outputs；確認矩陣不重複

## §7 嚴謹度繼承（/scientific-rigor §2-§7）

- 每項觀察 + 假設標 L1-L5 evidence tier ribbon
- §2 Credibility Score 對應 §scientific-rigor §3 信度量化
- §6 Conflict Scan 對應 §scientific-rigor §4 DAG / §5 reproducibility
- §7 Decision Lock 對應 `feedback_productive_failure_reopen_threshold` C1/C2/C3

## §8 Failure Mode & Diagnostics

### 1. Score 過度樂觀（每維皆 max）

**症狀**：TOTAL 90-100，但實際無觀察支撐。

**修補**：
- §2 必含 **falsifier observable**（WRAP）— 若填不出 → 觀察維度自動降 0
- §6 必跑 conflict scan — 若有 ≥1 NEGATIVE → tier 降 1 階 → score 重算

### 2. Cynefin 誤分（Complex → Complicated）

**症狀**：Complex 問題誤套 best-practice → 後續 pilot 偏離。

**修補**：
- §0 強制問「**相同行動是否曾重複產生可預測結果**？」
- 若 Maybe → Complex → 禁套 best-practice，必須 safe-to-fail probe

### 3. Conflict scan 漏 NEGATIVE 結論

**症狀**：§6 空表，但 MEMORY.md 「Concluded」section 有相關 NEGATIVE。

**修補**：
- §6 強制執行 `grep -i "concluded\|negative\|no-go" InterSubMod/MEMORY.md`
- 強制執行 `find InterSubMod/docs/reports/validated -name "*NEGATIVE*"`

### 4. Quick pilot 評估失準

**症狀**：§4 pilot 預期 30min 跑完 8h，或 metric threshold arbitrary。

**修補**：
- §4 必含明確 checkpoint + 中止條件
- 若 pilot 跑超時 → 自動標 OVER_BUDGET + 升級為 full cycle

### 5. Audit 過度頻繁（每決策都跑）

**症狀**：1 小時內被觸發 ≥ 3 次，用戶開始 ignore。

**修補**：
- rate limit：每 cycle 1 次 / hour
- audit.json mtime 檢核 — 若 < 1h 前剛跑，return cached verdict
- SessionStart 1 次 / session 上限

## §9 與 /implementation-notes 時序

```
用戶說「實作 <SPEC>」
    ↓
/pre-decision-audit init <topic>   ← Phase 1 PRE
    ↓ (≤30min)
score + verdict
    ↓
verdict_GO ──→ /cycle-init (P0) + /implementation-notes init  ← Phase 2 PROCESS
                    ↓
              spec 實作 + 4 markers append
                    ↓
              /run-evaluator (P5)  ← Phase 3 POST
                    ↓
              /implementation-notes finalize → docs/reports/validated/
verdict_PROBE ──→ §4 Quick Pilot (<30min) → re-run score
verdict_NO-GO ──→ MEMORY.md Concluded + postmortem (if mid-cycle)
```

`/implementation-notes init` 將自動讀取 `state/cycles/<cycle_id>/audit.json` 並把 assumptions 預填到「🟡 折衷考量」section。

## §10 FAQ

**Q: audit 跑 > 30min 怎麼辦？**
A: 自動標 OVER_BUDGET。可能此 topic 不適合 audit（太大 / 太模糊）→ 拆 sub-topic，每個 sub 各跑 audit。

**Q: TOTAL score 邊界 30/60 是 arbitrary？**
A: 是。V3 dogfood + V5 累積後 calibrate。目前依 4 維皆 weak (~25-35) 視為 NO-GO，4 維皆 strong (~60+) 視為 GO，中間 PROBE。

**Q: 我已 inject hypothesis 了，還要 audit 嗎？**
A: 不用。本 skill 只在 entry 點呼叫（cycle 啟動 / 大改前），active cycle 中用 `/implementation-notes` 即可。

**Q: Cynefin 分到 Complex 域怎麼辦？**
A: 禁止套 §1-§7 best-practice checklist。改 safe-to-fail probe — 設計 3-5 個小實驗 parallel 跑，看哪個浮現訊號再 audit。
