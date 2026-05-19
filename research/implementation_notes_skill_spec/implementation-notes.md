<!--
建立時間: 2026-05-19 16:50
目標: /implementation-notes skill 本身的實作 living document (dogfooding)
處理範圍: implementation-notes skill + template + standalone_trigger 擴展
cycle_id: cycle_20260519-1650-impl_notes_skill
spec_id: implementation_notes_skill_spec
status: in_progress
advisory: on
關聯檔案:
  - InterSubMod/.claude/skills/implementation-notes/SKILL.md (本 skill 自身)
  - InterSubMod/templates/implementation-notes-template.md
  - InterSubMod/scripts/hooks/standalone_trigger.sh
  - /bip7_disk/liaoyoyo2001/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md (本輪 plan)
-->

# Implementation Notes: /implementation-notes Skill (Dogfooding)

> **Purpose**: 用 /implementation-notes 機制紀錄機制本身的實作 — 完美 dogfooding 驗證 skill 可用性。
> **Trigger**: 用戶 plan ack 後啟動實作。
> **HTML companion**: html-report-build standalone 模式自動觸發 `implementation-notes.standalone.html`。
> **Cap**: ≤400 行。

## Frontmatter

- **Spec source**: `/bip7_disk/liaoyoyo2001/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md`（plan v2）
- **AI session**: 2026-05-19 16:00-17:00
- **Last updated**: 2026-05-19 16:50
- **Line count**: ~280 / 400

---

## 🔵 設計決定（Design Decisions）

### [2026-05-19 16:30] Per-task scope（拒絕 Hybrid global+per-task）
- **Status**: Accepted（用戶 AskUserQuestion ack）
- **背景** (Context): 4 個選項 (Per-task / Per-session / Global / Hybrid)；trade-off 是「navigation 集中 vs scope 清晰」
- **決定** (Decision): "We will use Per-task scope — 每個 spec 1 個 implementation-notes.md，存在 research/<topic>/ 或 docs/concepts/<topic>/"
- **理由** (Consequences):
  - ＋scope 清晰，避免 global 文件過長 (>400 lines TPP cap)
  - ＋對應 ADR per-decision unit + Cline activeContext per-task
  - −跨 task 知識散落（無 global index）→ 用戶接受此 trade-off
- **影響範圍**: 路徑慣例 `research/<topic>/` + finalize 移檔到 `docs/reports/validated/`
- **Revisit if**: 跨 task 知識搜尋需求出現 → 評估補 global index（未來輪次）
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐（用戶明示 ack）

### [2026-05-19 16:30] Hybrid 觸發策略（manual marker + AI 自我偵測 advisory）
- **Status**: Accepted（用戶 ack）
- **背景**: 3 選 (Pure manual / Hybrid / Full auto hook)
- **決定**: "We will use Hybrid — manual marker P1 + AI advisory P2 (≥2 linguistic markers threshold)"
- **理由**:
  - ＋manual priority 保證用戶意圖控制
  - ＋AI advisory 避免用戶遺漏記錄
  - ＋threshold 避免 noise
  - −Full auto hook 留待後續（本輪不做）
- **Revisit if**: 用戶 V3 fresh session 試用 friction 高 → 加 advisory hook
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐

### [2026-05-19 16:45] 重用 html-report-build standalone（不新建 HTML middleware）
- **Status**: Accepted
- **背景** (Context): plan 設計時 explore 結果顯示 html-report-build standalone 模式已內建 sticky TOC + collapsible cards + tier badge
- **決定**: "We will use html-report-build standalone mode + standalone_trigger.sh hook — no new HTML middleware"
- **理由**:
  - ＋零新增 HTML 生成 code
  - ＋既有 hook 已 PostToolUse Edit|Write registered
  - ＋對齊 Anthropic LLM-direct HTML pattern
  - −standalone_trigger.sh 需擴展 path 規則覆蓋 `research/*/implementation-notes.md`（修補後完成）
- **影響範圍**: `scripts/hooks/standalone_trigger.sh` 加 4 行 case 規則
- **Revisit if**: standalone 模式無法表達 4 sections 視覺差異 → 考慮新 mode
- **Evidence tier**: L1 ⭐⭐⭐⭐⭐（dry-run 驗證觸發成功）

---

## 🟠 偏離之處（Deviations）

### [2026-05-19 16:40] AGENTS.md 加 §14.1 而非 §11（plan 估錯位置）
- **Status**: Accepted
- **規範要求**: plan 寫「§11 補 implementation-notes living document 工作流提示」
- **實作偏離**: 改加到 §14 任務切割段（§14.1 子段），不在 §11
- **理由**: §11 是「執行檔案 IO 顯示規則」純 IO 規範，不適合放 living document 工作流；§14 是「任務切割與 Agent 啟動」含「紀錄成文件」原則，正是 implementation-notes 對應位置
- **風險評估**: 低 — 內容定位正確，只是 section 編號不同
- **回退路徑**: 若用戶想要 §11 — 可移過去，但語義不合
- **Revisit if**: 用戶 review 後要求調整位置

---

## 🟡 折衷考量（Trade-offs）

### [2026-05-19 16:35] 不做 advisory hook（implementation_notes_advisory.sh）
- **Status**: Accepted（本輪 scope 外，依用戶選 Hybrid 推薦）
- **方案 A** (採用): "We will skip hook in this round — 純 manual + AI 自我偵測 advisory 已足"
  - ＋零 hook overhead
  - ＋AI 自我偵測 threshold (≥2 markers) 已有效
- **方案 B** (rejected): 立即加 hook 自動偵測 PostToolUse marker
  - −false positive 風險
  - −用戶尚未證明 manual+advisory 不足
- **方案 C** (deferred): V3 fresh session 試用後再評估
- **採用 A 理由**: plan Risk & Mitigation 列了「advisory 過度敏感」風險，hook 會放大此風險；先 ship 純 skill + advisory，後續再 hook
- **Tier check**: Cohen 無對應（非數值決策）；DAG: hook → false positive → noise → advisory mode disabled
- **Revisit if**: 用戶實際使用 1 週後反饋「忘記記錄太多次」
- **Evidence tier**: L3 ⭐⭐⭐（hypothetical 預測 hook 不必要）

### [2026-05-19 16:48] 加 5 個業界補強到 template B（plan 後期加入）
- **Status**: Accepted
- **方案 A** (採用): "We will integrate 5 industry best practices (Status / Review trigger / Active voice / Protected-decision marker / Lore + 400 cap) into template directly"
  - ＋一次到位達業界對齊
  - ＋template 立即 enterprise-grade
- **方案 B** (rejected — KISS): minimal 4 sections，業界補強留 SKILL.md 描述
  - −template 與 SKILL.md spec 不一致
  - −用戶必須讀 SKILL.md 才知補強規則
- **方案 C** (deferred): V0.1 minimal → V0.2 加業界補強
  - −第一版就應有 best practices baseline
- **採用 A 理由**: 用戶明示要求「同時查看網路上相關大神對著類 skills 的細節」— 業界對齊是用戶需求
- **Revisit if**: 用戶覺得 template 過 verbose
- **Evidence tier**: L2 ⭐⭐⭐⭐（業界 4 來源 cited）

---

## 🔴 未決問題（Open Questions）

### [2026-05-19 16:50] 是否 dogfood 本實作為 V3 替代？
- **Status**: open
- **Question**: 本 dogfood file 是否可以替代 plan V3 (Fresh session E2E)？
- **Context**: V3 規格要求用戶開新 session 跑 demo_spec；本 dogfood 是當前 session 自我紀錄
- **Options**:
  - A: dogfood 計入 V3 ✓（本實作即驗證；用戶可看 HTML companion 確認）
  - B: dogfood 算 V0 sanity check，V3 仍需新 session
  - C: dogfood 後跑 reviewer agent 視為 V3 替代
- **Default if no answer**: B（保守，未來 fresh session V3 仍跑）
- **Revisit if**: 用戶 commit 後想立即看 HTML companion 確認
- **Priority**: minor（不阻塞 commit）

### [2026-05-19 16:50] reviewer agent V1 何時跑？
- **Status**: open
- **Question**: pr-review-toolkit:code-reviewer 在 commit 前跑還是 commit 後跑？
- **Context**: plan V1 規格說「spawn reviewer agent ≥ 4/5」但未指定時機
- **Options**:
  - A: commit 前跑 — 風險: reviewer findings 需立即修補才 commit
  - B: commit 後跑 — 風險: findings 需 follow-up commit
- **Default if no answer**: B（commit 先確保 baseline，findings 走 commit 0e701cb 模式）
- **Priority**: minor

---

## 📚 Lore（Prior Gotchas / Non-obvious Constraints）

### [2026-05-19] standalone_trigger.sh `set -euo pipefail` 嚴格模式
- **Constraint**: hook 用 `set -euo pipefail` — 任何 unset var 或 pipe failure 會 exit non-zero
- **Why it matters**: 加新 path case 時不可引入 unset var；測試時要 dry-run 確認 exit 0
- **Evidence**: `scripts/hooks/standalone_trigger.sh:23`

### [2026-05-19] ReadTool 500 行 truncate 邊界
- **Constraint**: Claude Code ReadTool 預設 read 2000 行但 hook 提示在 500 行 warning
- **Why it matters**: 400 行 cap 是緩衝；超過 400 → 拆 spec 或 trim
- **Evidence**: [PhotoStructure TPP](https://photostructure.com/coding/claude-code-tpp/) 「The ReadTool quietly truncates files beyond 500 lines」

### [2026-05-19] AGENTS.md plan-prescribed §11 vs 實作 §14.1
- **Constraint**: Plan 可能標錯 section 編號 — 實作時應依語義適配
- **Why it matters**: 盲從 plan section 編號會放錯位置；應 grep AGENTS.md 找最佳語義位置
- **Evidence**: 本實作 deviation entry 已紀錄

---

## Provenance Footer

- **Commit hash**: `<待 commit 後填>`
- **Build time**: 2026-05-19 16:50
- **Skill version**: `/implementation-notes` v0.1
- **HTML companion**: `research/implementation_notes_skill_spec/implementation-notes.standalone.html`（standalone_trigger 將觸發）
- **Line count**: ~280 / 400 ✅
- **Linked artifacts**:
  - Plan: `/bip7_disk/liaoyoyo2001/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md`
  - SKILL: `InterSubMod/.claude/skills/implementation-notes/SKILL.md`
  - Template: `InterSubMod/templates/implementation-notes-template.md`
  - Hook patch: `InterSubMod/scripts/hooks/standalone_trigger.sh` (新加 4 行 case)
  - CLAUDE.md update: `InterSubMod/.claude/CLAUDE.md` §3 (42→43)
  - AGENTS.md update: `InterSubMod/AGENTS.md` §14.1
