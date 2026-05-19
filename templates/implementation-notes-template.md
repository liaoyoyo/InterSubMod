<!--
建立時間: {YYYY-MM-DD HH:MM}
目標: <spec name> 實作過程 living document
處理範圍: <topic>
cycle_id: <cycle_YYYYMMDD-HHMM-slug>
spec_id: <spec slug>
status: in_progress | validated
advisory: on  # set to "off" 關閉 AI 自我偵測 advisory mode
關聯檔案:
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md (skill source)
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md §2-§7 (evidence tier)
  - InterSubMod/.claude/skills/html-report-build/SKILL.md (standalone HTML 生成)
-->

# Implementation Notes: <spec name>

> **Purpose**: 在 spec 實作過程中持續記錄 4 類訊息給用戶隨時翻閱。
> **Trigger**: Hybrid (manual marker `[決策]`/`[偏離]`/`[折衷]`/`[未決]` + AI 自我偵測 advisory ≥2 markers)。
> **HTML companion**: html-report-build standalone 模式自動觸發 `{filename}.standalone.html`。
> **Cap**: ≤400 行（>500 行 ReadTool 會 truncate）。

## Frontmatter

- **Spec source**: `<path/to/spec.md>` 或對話原文（若無正式 spec 文件）
- **AI session**: <session_id 或日期>
- **Last updated**: {YYYY-MM-DD HH:MM}
- **Line count**: <N> / 400

---

## 🔵 設計決定（Design Decisions）

> 規格含糊處的選擇 — AI 在無明確 spec 時做的判斷
> ADR Nygard style: active voice "We will..."

### [YYYY-MM-DD HH:MM] <decision title>
- **Status**: Proposed | Accepted | Superseded by [entry-id] | Closed
- **背景** (Context): <spec 含糊處 / forces at play>
- **決定** (Decision, active voice): "We will use X (not Y)..."
- **理由** (Consequences): <positive + negative + neutral 影響>
- **影響範圍**: <檔案 / 行為 / 後續任務>
- **Revisit if**: <重評條件，如「Cohen's d > 0.5 時重評」>
- **Evidence tier**: L1/L2/L3/L4/L5 ⭐

<!-- 若為 PI / 用戶 lock 不可改項：用 Protected-Decision marker -->
<!-- BEGIN USER-SPECIFIED -->
<!-- **Decision**: <non-negotiable 內容> -->
<!-- **DO NOT change**: <理由 + lock date> -->
<!-- **Rationale**: <PI / user audit reference> -->
<!-- END USER-SPECIFIED -->

---

## 🟠 偏離之處（Deviations）

> 有意偏離規範的部分 + 原因

### [YYYY-MM-DD HH:MM] <deviation title>
- **Status**: Proposed | Accepted | Superseded | Closed
- **規範要求**: <原 spec 要求>
- **實作偏離**: <實際做的>
- **理由**: <為何偏離>
- **風險評估**: <對 spec 完整性影響>
- **回退路徑**: <如何修回原規範>
- **Revisit if**: <重評條件>
- **Evidence tier**: L1/L2/L3/L4/L5 ⭐

---

## 🟡 折衷考量（Trade-offs）

> 替代方案 + 為何選目前做法
> ADR Decision Log 3 欄: Decision + Reason + Alternative

### [YYYY-MM-DD HH:MM] <tradeoff title>
- **Status**: Proposed | Accepted | Superseded | Closed
- **方案 A** (採用，active voice): "We will X..." + Pros
- **方案 B** (備選 — rejected): <做法 + 為何 rejected>
- **方案 C** (備選 — deferred): <做法 + 為何 deferred>
- **採用 A 理由**: <Why A>
- **Tier check**: <Cohen / CI / DAG 對應>
- **Revisit if**: <重評條件>
- **Evidence tier**: L1/L2/L3/L4/L5 ⭐

---

## 🔴 未決問題（Open Questions）

> 待用戶確認或修改的事項

### [YYYY-MM-DD HH:MM] <question title>
- **Status**: open | answered (with answer + timestamp) | closed (no action)
- **Question**: <要問用戶什麼>
- **Context**: <為何問>
- **Options**: <候選答案 A / B / C>
- **Default if no answer**: <若用戶不回的 fallback>
- **Revisit if**: <重評條件>
- **Priority**: critical | major | minor

---

## 📚 Lore（Prior Gotchas / Non-obvious Constraints）

> Cross-session institutional memory — 累積 prior gotchas / edge cases / 不明顯約束
> 對應 TPP Lore section + InterSubMod known-pitfalls

### [YYYY-MM-DD] <gotcha title>
- **Constraint**: <不明顯的約束 / 約束來源>
- **Why it matters**: <對 spec 影響>
- **Evidence**: <file:line 或 prior session 紀錄>

---

## Provenance Footer

- **Commit hash**: `<git rev-parse HEAD>`
- **Build time**: {YYYY-MM-DD HH:MM:SS}
- **Skill version**: `/implementation-notes` v0.1
- **HTML companion**: `<filename>.standalone.html`
- **Line count**: <N> / 400 (≤400 強制 cap; >400 → 拆 spec 或 trim redundancy)
- **Linked artifacts**:
  - Spec: `<path/to/spec.md>`
  - Active cycle: `InterSubMod/state/cycles/<cycle_id>/`
  - Pre-registration: `InterSubMod/research/<topic>/00_INDEX.md`
  - Postmortem (if NEGATIVE): `InterSubMod/docs/postmortems/<YYYYMMDD>_<topic>.md`
