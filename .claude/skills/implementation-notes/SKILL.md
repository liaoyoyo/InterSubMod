---
name: implementation-notes
description: |
  Per-spec living document — 在 spec/feature 實作過程中持續記錄 4 類訊息（設計決定 / 偏離 / 折衷 / 未決問題）到 .md source，自動觸發 html-report-build standalone 模式生成 .html 給用戶隨時翻閱。Hybrid 觸發：用戶輸入 [決策]/[偏離]/[折衷]/[未決] 4 個 marker → AI 立即 append；AI 自我偵測 design decision moment → 主動問用戶是否記錄。與 /report (session-end retrospective) / /structured-tech-report (spec 完成後 13 段) / /weekly-report (週彙整) 互補：本 skill = process-time prospective live tracking。
  USE WHEN: 開始新 spec 實作、用戶說「實作 <SPEC>」、用戶輸入 4 markers ([決策]/[偏離]/[折衷]/[未決])、AI 偵測「我選 A 不選 B 因為 X」/「考量到 Y」/「但這假設 Z」/「不確定 W」等 design decision linguistic markers。
  SKIP WHEN: 純 code edit / build / commit / typo fix（無 design decision moment）、單一 trivial 改動、純 docs 寫作無實作偏離、retrospective session summary（用 /report）、spec 完成後 13 段技術報告（用 structured-tech-report）。
allowed-tools: Read, Write, Edit, Glob, Grep
user-invocable: true
paths: ["research/**/implementation-notes.md", "docs/concepts/**/implementation-notes.md", "docs/reports/validated/**/*implementation_notes*.md"]
---

# Implementation Notes — Per-Spec Living Document Skill

> **Purpose**: 在 spec 實作過程中持續記錄 4 類訊息，讓用戶隨時翻閱 .standalone.html 了解 AI 如何詮釋 spec。
> **不是**: session-end retrospective（用 /report）/ spec 完成後 13 段技術報告（用 /structured-tech-report）/ 週彙整（用 /weekly-report）。

## Phase & Chain Position

**Phase**: cross-cutting — P0-P6 任何階段都可觸發
**Position**: process-time prospective live tracking（介於 /report retrospective 與 /weekly-report 週級之間）

## Dependencies

- **Uses**:
  - `/html-report-build` (standalone 模式) — 生成 .standalone.html companion
  - `/doc-standards` — filename + metadata block + evidence tier ribbon
  - `/scientific-rigor` §2-§7 — 每個 entry 標 L1-L5 evidence tier ribbon
  - `/pre-decision-audit`（**上游 — 假說驗證三層樓首層**）— `init` 時讀取 `state/cycles/<cycle_id>/audit.json` 並把 assumptions 預填到「🟡 折衷考量」section（若 audit verdict=GO）
- **Used by**（預計）: `/cycle-init` (P0) / `/research-loop` (P1) / `/cpp-change` (C++ 修改前) / 任何 design decision moment
- **Reads**: `InterSubMod/templates/implementation-notes-template.md`、`state/cycles/<cycle_id>/audit.json`（若 pre-decision-audit 已跑）
- **Writes**: `InterSubMod/research/<topic>/implementation-notes.md`（進行中）或 `InterSubMod/docs/reports/validated/{YYYY}/{MM}/{YYYYMMDD}_<topic>_implementation_notes_01.md`（完成後）

**假說驗證三層樓位置**：`/pre-decision-audit` (pre, ≤30min) → **`/implementation-notes` (process, live)** → `/run-evaluator` (post, P5 tier)

## Failure Mode & Diagnostics

| 症狀 | 根因 | 修法 |
|------|-----|------|
| marker 衝突（用戶寫 `[決策]` AI 解為 `[折衷]`） | 內容跨 section | append 到主 section + cross-link 其他 |
| AI advisory 過度敏感（每個「但」「考量」都問用戶） | threshold 未啟用 | ≥2 markers 才問；單一片語不問 |
| HTML 未自動生成 | standalone_trigger 路徑不命中 | 確認 `research/**/implementation-notes.md` 在 hook path 內 |
| cross-session lose context | 新 session 不知文件存在 | SessionStart hook 注入 + finalize 時 update CURRENT_FOCUS |
| 文件 >400 行 | append-only 累積 | 拆 spec 或 trim redundancy（移 stale 對話紀錄）|

---

## §1 三個 sub-command

### 1.1 `init <topic>`

從 `InterSubMod/templates/implementation-notes-template.md` 複製到 `InterSubMod/research/<topic>/implementation-notes.md`，填入：
- `cycle_id`（從 active cycle 或新建）
- `spec_id`（topic slug）
- `status: in_progress`
- `建立時間`（current YYYY-MM-DD HH:MM）

### 1.2 `append <section> <text>`

在指定 section 加 timestamped entry：

| Section | Marker | Badge |
|---------|--------|-------|
| Design Decisions | `[決策]` | 🔵 |
| Deviations | `[偏離]` | 🟠 |
| Trade-offs | `[折衷]` | 🟡 |
| Open Questions | `[未決]` | 🔴 |
| Lore | `[gotcha]` | 📚 |

Entry 格式依 ADR Nygard standard（active voice "We will..."）+ 5 個業界補強（Status / Revisit if / Tier / Protected marker / Active voice）。

### 1.3 `finalize`

完成 spec 後：
1. 設 `status: validated`
2. 移檔到 `InterSubMod/docs/reports/validated/{YYYY}/{MM}/{YYYYMMDD}_<topic>_implementation_notes_01.md`
3. trim redundancy（移 stale append-only 對話紀錄）
4. 驗證 ≤400 行
5. commit + update CURRENT_FOCUS.md 標記完成
6. 觸發 html-report-build standalone 生成最終 .standalone.html

---

## §2 Hybrid 觸發規則

### 2.1 Manual marker（用戶輸入）→ AI 立即 append（無需確認）

```
用戶: [決策] 用 KDE bandwidth 0.3 不選 0.5 因為 small sample
AI: 已 append 到 ## 🔵 設計決定 section（無需再問）
```

4 個 markers 互斥但 entry 可 cross-link 多 section（如「決策牽涉折衷」）。

### 2.2 AI 自我偵測 → 主動 prompt（≥2 linguistic markers 才問）

**Linguistic markers**（closed enumeration — 擴充需用戶 opt-in）:
- 設計選擇類：「考量因素」「為何不選」「我選 A 不選 B」「trade-off」
- 假設類：「但這假設」「假定 X」「前提是 Y」
- 不確定類：「不確定」「TBD」「待確認」「mock 替代」
- 暫存類：「目前 hardcode 等以後改」「skip 這部分為了 MVP」「之後可能要重評」

**Threshold（operationalized）**:
- **「同一段」定義 = within one AI assistant response (single turn)** — 一個 AI 回覆 turn 內出現 ≥2 markers 才啟動 advisory
- 單一片語不問（避免噪音）
- ≥2 markers within one turn → AI prompt：「要記下這個 [類型] 嗎？」
- 用戶可關閉 advisory mode（在 implementation-notes.md frontmatter 加 `advisory: off` — 純 documentation marker, AI 讀檔時識別）

### 2.3 Open Questions 的 evidence tier 例外規則

- **🔵 Design Decisions / 🟠 Deviations / 🟡 Trade-offs entry**: 必標 L1-L5 tier ribbon（依 /scientific-rigor §2）
- **🔴 Open Questions entry**: 預設 L5（pre-evidence by definition — questions are unresolved by nature）
- 一旦 Open Question 被 answered → 升 tier 到 L1-L4（依答案 evidence 強度）+ Status 改 `answered`

### 2.3 Cross-Session continuity

- SessionStart hook 注入 CURRENT_FOCUS.md 含 active implementation-notes 路徑
- 新 session 啟動時 AI 主動 Read 該檔，恢復 context
- 若 active implementation-notes mtime > 7 天 → 提示用戶 finalize 或 update

---

## §3 路徑慣例

| 狀態 | 路徑 |
|------|------|
| 進行中 | `InterSubMod/research/<topic>/implementation-notes.md` |
| 概念探索 | `InterSubMod/docs/concepts/<topic>/implementation-notes.md` |
| 完成 (validated) | `InterSubMod/docs/reports/validated/{YYYY}/{MM}/{YYYYMMDD}_<topic>_implementation_notes_01.md` |
| HTML companion | 同路徑 `.standalone.html`（standalone_trigger.sh 自動生成）|

---

## §4 與既有 4 個 skill 分工

| Skill | 時機 | 形式 | Trigger | 互補性 |
|-------|------|------|---------|--------|
| `/report` | session-end | retrospective summary | 「寫報告」 | implementation-notes finalize 後可呼叫 /report 寫 session log |
| `/weekly-report` | 週末 | 週級多主題 | 「週報」 | implementation-notes 多筆彙整入週報 |
| `/structured-tech-report` | spec 完成後 | 13-段技術報告 | 「技術報告」 | 從 implementation-notes finalize 結果寫 13 段 |
| **`/implementation-notes`** (本) | spec **進行中** | **process-time live** | 「實作 <SPEC>」/ 4 markers / AI advisory | live tracking source-of-truth |

**角色區隔**:
- /report: 回顧本 session 做了什麼（retrospective）
- /implementation-notes: 即時記錄 spec 詮釋偏離（prospective live）
- /structured-tech-report: spec 完成後完整技術敘事

---

## §5 嚴謹度繼承 (`/scientific-rigor` §2-§7)

每個 entry 必繼承 `InterSubMod/.claude/skills/scientific-rigor/SKILL.md`：

- **§2 證據分級**: 標 L1-L5 ⭐⭐⭐⭐⭐ ribbon
- **§3 Effect Size**: 數字 metric 必含 Cohen ribbon + CI（如 "+0.0112 F1, marginal < Cohen's small 0.2"）
- **§4 DAG**: 因果 claim 必 reference DAG mermaid 路徑
- **§7 Pre-registration**: 設計決定 entry 必對照事先註冊 H_預測（research/<topic>/00_INDEX.md）
- **§9.2 Postmortem**: 若 spec 走向 NEGATIVE 收尾 → finalize 時 + postmortem.md（互補不重複）

---

## §6 業界 5 個補強落地（2026-05-19 webspec audit）

### 6.1 Status 欄（ADR Nygard）

每個 entry 標：
- `Proposed` (剛記，未 ack)
- `Accepted` (用戶/PI ack)
- `Superseded by [entry-id]` (被新 entry 取代)
- `Closed` (open question 已回答)

### 6.2 Review Triggers（ADR + Augment Code）

每個設計決定 / 折衷 entry 加 `**Revisit if**: <條件>` —
範例：「Cohen's d > 0.5 時重評」「migration 完成後重評 fallback」

### 6.3 Decision 句式 "We will..." active voice

- ✓ 「We will use KDE bandwidth 0.3 not 0.5 because small sample」
- ✗ 「可能採用 0.3」「考慮用 0.3」（hedged 用 Status=Proposed 表達）

### 6.4 Protected-Decision Markers（Augment Code）

非協商項用 HTML comment 標記：

```markdown
<!-- BEGIN USER-SPECIFIED -->
**Decision**: F1 SoT = paired-pure pipeline only
**DO NOT change**: PI 已定，非實驗變項
**Rationale**: 2026-05-XX user audit confirmed
<!-- END USER-SPECIFIED -->
```

AI 看到此標記**不可** override；只能 propose superseding entry 並標 `Status: Proposed (pending user-specified review)`。

### 6.5 Lore Section + ≤400 行 Cap（TPP）

- 4 sections 之外加 `## 📚 Lore`（累積 prior gotchas / 不明顯約束 / edge cases）
- 強制 ≤400 行（ReadTool 在 500 行 truncate）
- 超過 400 行 → 拆 spec 或 trim redundancy
- finalize 必跑 trim：移 stale append-only 對話紀錄

---

## §7 業界對齊

| 框架 | 對應點 |
|------|------|
| [ADR Nygard standard](https://adr.github.io/) | §6.1 Status + §6.3 Active voice + §6.2 Review trigger |
| [Augment Code Living Specs](https://www.augmentcode.com/guides/living-specs-for-ai-agent-development) | §6.4 Protected-Decision Markers + bidirectional spec flow |
| [Cline Memory Bank activeContext.md](https://docs.cline.bot/prompting/cline-memory-bank) | per-task scope + 3-欄 Decision Log |
| [PhotoStructure TPP Pattern](https://photostructure.com/coding/claude-code-tpp/) | §6.5 Lore + ≤400 行 cap + `_active/_done/` folder workflow |
| [Microsoft Azure Well-Architected ADR](https://learn.microsoft.com/en-us/azure/well-architected/architect-role/architecture-decision-record) | Status / Context / Decision / Consequences |
| [Anthropic Sprint Contracts](https://www.anthropic.com/engineering/harness-design-long-running-apps) | sprint-level spec evolution + Memory of work |
| [MkDocs Living Documentation](https://www.mkdocs.org/) | source `.md` → companion `.html` |
| [Diátaxis Explanation framework](https://diataxis.fr/) | implementation-notes ≈ Explanation type |

---

## §8 啟動聲明（觸發時必說）

當用戶輸入 4 markers 或 AI 啟動本 skill 時，**首句必明確告知**：

> 「啟動 /implementation-notes — 記錄到 `InterSubMod/research/<topic>/implementation-notes.md` 的 [section]。HTML companion 將由 standalone_trigger 自動生成。」

避免用戶誤以為 AI 隨口閒聊未記錄。

---

## §9 反 pattern（看到立即修）

| 反 pattern | 為何 | 正確做法 |
|-----------|------|----------|
| 4 markers entry 用 hedged 語言（「可能」「也許」） | 違反 §6.3 active voice | 用 `Status: Proposed` 表達不確定 |
| 同一 spec 多個 implementation-notes 散落 | per-task 原則違反 | 用戶選 per-task 1 個；split spec 才 split notes |
| append-only 不 trim → >400 行 | ReadTool truncate 500 行 | finalize 必 trim |
| Override Protected-Decision marker 內容 | 違反 §6.4 | 只 propose superseding entry，不直接改 |
| AI 每個「但」都 prompt 用戶記錄 | advisory 過度敏感 | ≥2 markers 才問 |
| Implementation-notes 與 postmortem 重複 | 角色重疊 | implementation-notes = in-progress live；postmortem = NEGATIVE 後 SRE blameless |

---

## §10 Quick Reference Card

```
用戶觸發:
  [決策] <text> → AI append to 🔵 Design Decisions
  [偏離] <text> → AI append to 🟠 Deviations
  [折衷] <text> → AI append to 🟡 Trade-offs
  [未決] <text> → AI append to 🔴 Open Questions
  [gotcha] <text> → AI append to 📚 Lore
  /implementation-notes finalize → status=validated + move + commit

AI advisory:
  detect ≥2 linguistic markers → prompt 「要記下這個 [類型] 嗎？」

HTML 生成:
  source .md edit → standalone_trigger.sh → .standalone.html

Status lifecycle:
  Proposed → Accepted → (Superseded by | Closed)

Protected-Decision lifecycle:
  <!-- BEGIN USER-SPECIFIED --> ... <!-- END USER-SPECIFIED -->
  AI 不可 override，只能 propose superseding entry

Cap:
  ≤400 行（>500 行 ReadTool truncate）
```

---

## 參考來源

- ADR Nygard: <https://adr.github.io/>
- Augment Code Living Specs: <https://www.augmentcode.com/guides/living-specs-for-ai-agent-development>
- Cline Memory Bank: <https://docs.cline.bot/prompting/cline-memory-bank>
- PhotoStructure TPP: <https://photostructure.com/coding/claude-code-tpp/>
- Martin Fowler ADR: <https://martinfowler.com/bliki/ArchitectureDecisionRecord.html>
- Microsoft Azure ADR: <https://learn.microsoft.com/en-us/azure/well-architected/architect-role/architecture-decision-record>
- Diátaxis: <https://diataxis.fr/>
- InterSubMod 內部:
  - `InterSubMod/.claude/skills/scientific-rigor/SKILL.md` §2-§7
  - `InterSubMod/.claude/skills/html-report-build/SKILL.md` standalone 模式
  - `InterSubMod/.claude/skills/doc-standards/SKILL.md`
  - `InterSubMod/templates/postmortem.md`
  - `InterSubMod/templates/research_index.md`
