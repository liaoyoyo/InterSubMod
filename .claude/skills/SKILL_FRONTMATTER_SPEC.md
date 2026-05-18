<!--
建立時間: 2026-05-18
parent: P4 audit SKILL.md 依賴宣告統一
status: validated spec
purpose: 對齊 harness/harness-skills + Anthropic skill loader spec — InterSubMod skill frontmatter 統一規範
-->

# SKILL.md Frontmatter Spec v1 (2026-05-18)

> **來源**: P4 Architect agent E8 + Researcher agent 結論「對齊 harness/harness-skills 規範」
> **業界對照**: [harness/harness-skills](https://github.com/harness/harness-skills) + Anthropic skill loader spec + [HKUDS/OpenHarness](https://github.com/HKUDS/OpenHarness)
> **適用**: `InterSubMod/.claude/skills/<name>/SKILL.md`

---

## §1 必要欄位（required）

| 欄位 | 型別 | 範例 | 說明 |
|------|------|------|------|
| `name` | string | `scientific-rigor` | kebab-case，對應 skill 目錄名 |
| `description` | string OR YAML block scalar | `"科學嚴謹度元方法論..."` | LLM trigger 用；含 USE WHEN + SKIP WHEN；建議 ≤500 chars；中文 markdown emphasis 用 quoted string 或避免 `**` 開頭（YAML alias bug） |

## §2 推薦欄位（recommended for new skills）

| 欄位 | 型別 | 範例 | 說明 |
|------|------|------|------|
| `user-invocable` | bool | `true` / `false` | 用戶可否直接 `/skill-name` 觸發 |
| `allowed-tools` | string OR array | `Read, Write, Edit` 或 `["Read","Bash"]` | 限制可用工具（Anthropic spec）|
| `paths` / `globs` | array | `["src/**/*.cpp"]` | path-scoped 條件式載入（節省 always-loaded context；對齊 Cursor `.cursor/rules` + Anthropic skill spec）|

## §3 業界對照延伸欄位（optional, future）

| 欄位 | 型別 | 對齊 spec | 範例 |
|------|------|-----------|------|
| `version` | string | InterSubMod 既有 | `0.2.0` |
| `tags` | array | InterSubMod 既有 | `["feature","observation","P3-PILOT"]` |
| `phase` | enum | harness/harness-skills | `"P3-PILOT"` / `"meta"` / `"reporting"` |
| `mcp_required` | array | harness/harness-skills | `["knowledge", "context7"]`（依賴的 MCP server）|
| `dependencies` | object | OpenHarness | `{uses: ["/scientific-rigor"], used_by: ["/cycle-init"]}` |
| `model` | string | Anthropic | `"inherit"` / `"opus-4-7"` / `"sonnet-4-6"` |
| `isolation` | string | Anthropic | `"worktree"`（agent only，skill 不適用） |

## §4 body 必要章節（per feedback_skill_md_must_state_dependencies_and_diagnostics）

每個 SKILL.md body 應含以下 3 段（順序自由）：

1. `## Phase Chain Position` — Upstream / Downstream / Reads / Writes
2. `## 與 /scientific-rigor 元方法論的關係` — 引用 §X-Y hub 章節
3. `## Failure Mode & Diagnostics` — 症狀 × 可能原因 × 修法 表

✅ P1 audit fix 已落地：13 個 D2 缺項 skill 全補完（commit a7c1495）

## §5 Anti-patterns（業界禁忌）

| 反例 | 為何禁 | 修法 |
|------|-------|------|
| `description: **markdown emphasis 開頭**` | YAML 規範下 `**` 視為 alias reference → strict parser fail | 用 `"..."` quoted string 包，或移除開頭 `**` |
| description 不含 USE WHEN | LLM 不知何時 trigger | 加 USE WHEN: 「觸發詞 1」「觸發詞 2」 |
| description 不含 SKIP WHEN | trigger 過寬 | 加 SKIP WHEN: 「不該用於 X」 |
| description >1000 chars | 對齊 OpenAI 「map not 1000-page manual」+ Anthropic 「LLM-friendly 描述要簡」 | 縮到 ≤500，詳細放 body |
| frontmatter 含 `allowed-tools: [...]` array | Anthropic spec 不一致（部分 parser 只接受 string） | 用 string 形式 `"Read, Write, Edit"` |

## §6 audit script

```bash
# 跑 frontmatter audit（含 YAML valid / required fields / anti-pattern 檢核）
python3 InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/scan_script.py
```

輸出: D1-D5 評分 + YAML invalid 清單。

## §7 InterSubMod 現況（2026-05-18）

| 維度 | 現況 |
|------|------|
| YAML valid | 42/42 ✅（P0 已修 5 個）|
| 必要欄位完整 | 42/42 ✅ |
| user-invocable 明示 | 33/42 (78%) ✅ |
| paths-scoped | 3/42（其他多為 always-loaded — cache 94.6% 命中證明可接受）|
| Body §Phase Chain | 14 non-tool skill ✅（P1-D2 補完）|
| Body §scientific-rigor cross-ref | 18 skill 已 cross-ref ✅ |
| Body §Failure Mode | 14 non-tool skill ✅ |

## §8 未來升級路徑

優先級：
1. **P5 未做**: 對 21 工具類 skill 評估補 `paths`（如 image-gen / pptx-build 等明顯 file-type-scoped）
2. **P6 未做**: 引入 `mcp_required` / `phase` / `dependencies` 結構化欄位 — 評估 ROI 後決定
3. **P7 未做**: 引入 SkillTester benchmark（[ai-boost/awesome-harness-engineering](https://github.com/ai-boost/awesome-harness-engineering) 收錄）

## §9 元層說明

本 spec 是 P4 Architect E8 + Researcher 「SKILL.md 依賴宣告統一」結論。**現況已對齊業界 80%**（YAML valid + required fields + body 3 段全綠）；剩 20% 是 frontmatter 細部優化（paths / mcp_required / dependencies）非急性問題。

不強制現有 44 skills 升級到 §3 延伸欄位 — 業界 spec 仍在 evolve；evidence-based 決定何時升級（per `/scientific-rigor §8.3.1 reopen threshold C1/C2/C3`）。
