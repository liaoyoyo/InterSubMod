<!--
建立時間: {YYYY-MM-DD HH:MM}
目標: <topic> 進入 cycle / spec / 大改前 pre-decision audit
處理範圍: <topic slug>
cycle_id: <cycle_YYYYMMDD-HHMM-slug>
topic: <slug>
status: in_progress | scored | verdict_GO | verdict_PROBE | verdict_NO-GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/.claude/skills/pre-decision-audit/SKILL.md (skill source)
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md §2-§7 (evidence tier)
  - InterSubMod/.claude/skills/confirmation-protocol/SKILL.md §Cynefin (domain SoT)
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md (下游)
-->

# Pre-Decision Audit: <topic>

> **Purpose**: 進入 cycle / spec / 大改前，**≤30min** 完成 7 類審查
> **Verdict**: in_progress → scored → GO / PROBE / NO-GO
> **HTML companion**: `html-report-build` standalone 自動觸發 `<filename>.standalone.html`
> **Cap**: ≤400 行（>500 行 ReadTool 會 truncate）

## Frontmatter

- **Topic**: `<topic slug>`
- **Triggered by**: cycle-init / new-spec / cross-sample / tier-upgrade / NO-GO / ad-hoc
- **AI session**: <session_id 或日期>
- **Last updated**: {YYYY-MM-DD HH:MM}
- **Line count**: <N> / 400
- **Cycle ref**: `InterSubMod/state/cycles/<cycle_id>/`

---

## §0 🟤 Cynefin Domain Gate (front-gate)

> 決定後續 audit 走 best-practice 還是 probe-first。共用 `/confirmation-protocol §Cynefin` 定義。

- [ ] **Domain**: Clear / Complicated / **Complex** / Chaotic / Disorder
- [ ] **Test**: 「相同行動是否曾重複產生**可預測**結果？」
  - **Yes** → Clear / Complicated → 可用 §1-§7 best-practice checklist
  - **Maybe** → **Complex** → **強制 probe-first**（禁套 best-practice，改 safe-to-fail probe）
  - **No** → Chaotic → 先穩定（act → sense → respond）再 audit

**Domain decision**: `<填空>`
**Rationale**: `<為何分到此域>`

---

## §1 🟢 Observation Completeness Checklist

> 與此假設相關的觀察 — ✓ 已驗 / □ 待補 / ✗ 反例
> 每項標 L1-L5 evidence tier（對應 `/scientific-rigor §2`）

| Observation | 狀態 (✓/□/✗) | Evidence tier (L1-L5) | 來源 file:line |
|---|---|---|---|
| (e.g. 7-sample paired AUC ≥ 0.55) | ✓ | L1 ⭐⭐⭐⭐⭐ | `report.md:42` |
| (e.g. cross-sample LOH_AF scatter) | □ | — (not yet) | (not yet) |
| (e.g. within-group residual) | ✗ | L4 ⭐⭐ | `confound_report.md:120` |

---

## §2 🔵 Credibility Score (5-dim, 0-100)

| 維度 | 評分 | 0 (weak) | 10 (moderate) | 20 (strong) |
|---|---|---|---|---|
| **理論基礎** | __ | speculative | plausible | well-grounded |
| **觀察支撐** | __ | none | partial | complete (≥3 樣本一致) |
| **機制清晰度** | __ | speculative | plausible | clear DAG |
| **反例風險** | __ | high | medium | low (≥3 預期反例皆無) |
| **所需資源** | __ | >6h | 1-6h | <1h pilot |
| **TOTAL** | __ / 100 | | | |

**Falsifier observable (WRAP)** ⚠ **必填**：
> 「**若假設錯，我會看到什麼？**」
>
> Answer: `<填空 — 若填不出，§觀察支撐 自動降 0>`

**Reality-test 三個反例觀察**（Heath WRAP "Prepare to be wrong"）：
1. `<若假設成立應該看到 / 不該看到的觀察 1>`
2. `<觀察 2>`
3. `<觀察 3>`

---

## §3 🟡 Assumption Map 2×2 (Bland)

```
                  | known                | unknown
  ----------------|----------------------|--------------------------
  HIGH importance | (1) verify quickly   | (2) MUST validate FIRST ⚠
  ----------------|----------------------|--------------------------
  LOW importance  | (3) document only    | (4) defer / ignore
```

| Assumption | Importance | Known | Quadrant |
|---|---|---|---|
| (e.g. methylation independent of AF) | HIGH | UNKNOWN | (2) ⚠ |
| (e.g. paired pipeline F1 SoT) | HIGH | KNOWN | (1) |
| ... | | | |

**右上 (2) quadrant 假設**（MUST validate first）：
- `<列出>`

---

## §4 🟠 Quick Pilot Guide (<30min)

> 若 verdict 落 PROBE 區（30-60），執行此 pilot 後 re-run score

**Minimum steps**:

1. **Step 1** (5 min): Read `<file:line>` — context
2. **Step 2** (15 min): Run `<script.py>` with sample = `<HCC1395_5kHz or 其他>`
3. **Step 3** (5 min): Check `<metric>` vs `<threshold>`（e.g. AUC > 0.55 / Cohen d > 0.3）

**Checkpoint**:
- `≥ threshold` → **GO** full verify（升 §2 觀察支撐 → 10 or 20）
- `0.5 × threshold ~ threshold` → **PROBE** 再 1 樣本 pilot
- `< 0.5 × threshold` → **NO-GO** archive

**中止條件**:
- Pilot > 30min → 自動標 `OVER_BUDGET`，升級為 full cycle
- Pilot 缺資料 → 列入 §5 Gap Diagnosis

---

## §5 🔴 Gap Diagnosis

> 精確告知「為何無法現在決策」 + 缺什麼 + 優先級

| Missing | Impact (HIGH/MED/LOW) | Effort (h) | Priority (P0/P1/P2) |
|---|---|---|---|
| (e.g. COLO829 BAM symlink) | HIGH | 8h copy | P1 |
| (e.g. methylation-AF independence test) | HIGH | 2h script | P0 |
| ... | | | |

---

## §6 🟣 Evidence Conflict Scan

> 從 `InterSubMod/MEMORY.md` + `InterSubMod/docs/reports/validated/` grep prior conclusions
> ⚠ 強制執行：
> ```bash
> grep -i "concluded\|negative\|no-go" InterSubMod/MEMORY.md
> find InterSubMod/docs/reports/validated -name "*NEGATIVE*"
> ```

| Prior conclusion | Tier | Relation (support / conflict / dependent) | Source |
|---|---|---|---|
| (e.g. L2 collider bias) | concluded | dependent (must avoid) | `MEMORY.md:feedback_L2_collider_bias` |
| (e.g. ASM POSITIVE 32-66%) | ⭐4 | partial support | `snv_methylation_association.md` |
| (e.g. TO Germline FP NO-GO) | concluded | conflict (避開 same approach) | `germline_fp_identification_nogo.md` |

**Conflict count**: `<N>` — 若 ≥ 1 conflict → §2 反例風險 維度自動降級 1 階。

---

## §7 ⚫ Decision Threshold + Path

- **TOTAL score (from §2)**: `__ / 100`
- **Verdict**:
  - ≥ 60 → **GO** full verify
  - 30-60 → **PROBE** small-scale pilot first
  - < 30 → **NO-GO** reject / archive

### Next action

**If GO**:
- Run `/cycle-init <topic>` → 建立 cycle_id
- Run `/implementation-notes init <topic>` → 預填 §🟡 折衷 from this audit assumptions
- 開始 P1 research-loop

**If PROBE**:
- 跑 §4 Quick Pilot（≤30min）
- Pilot success → re-run `score` → 升 GO
- Pilot ambiguous → 二次 pilot 不同樣本
- Pilot failure → NO-GO

**If NO-GO**:
- Append `MEMORY.md` Concluded section（防 re-investigation）
- 若 mid-cycle → run postmortem template `InterSubMod/templates/postmortem.md`
- 否則直接 archive

### Decision lock

- [ ] **Lock decision**: Y / N
- 若 Y → 後續不可任意 reopen，除非 Productive Failure 3 條件 C1 / C2 / C3 至少一項：
  - **C1**: 新數據（new sample / new replicate）
  - **C2**: 新方法（new feature / new model）
  - **C3**: 新前置條件（upstream pipeline 修復）
- 對應 `MEMORY.md:feedback_productive_failure_reopen_threshold`

---

## Provenance Footer

- **Commit hash**: `<git rev-parse HEAD>`
- **Build time**: {YYYY-MM-DD HH:MM:SS}
- **Skill version**: `/pre-decision-audit` v0.1
- **HTML companion**: `<filename>.standalone.html`
- **Line count**: `<N> / 400`（≤400 cap; >400 → 拆 sub-topic）
- **Audit JSON**: `InterSubMod/state/cycles/<cycle_id>/audit.json`
- **Linked artifacts**:
  - Topic dir: `InterSubMod/research/<topic>/`
  - Active cycle: `InterSubMod/state/cycles/<cycle_id>/`
  - Next step (if GO): `InterSubMod/research/<topic>/implementation-notes.md`
  - Postmortem (if NO-GO): `InterSubMod/docs/postmortems/<YYYYMMDD>_<topic>.md`
