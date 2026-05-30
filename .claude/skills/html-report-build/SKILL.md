---
name: html-report-build
description: |
  LLM-direct HTML 報告生成（無 Python middleware）— 3 模式：report (Tailwind prose 排版) / slide (16:9 PPT 風) / standalone (PI 終版 + sticky TOC + 折疊 cards + SVG 圖)。Output 為 .md 的 companion 不取代。
  USE WHEN: 「想看 HTML 排版」「給 PI 看 preview」「PPT HTML」「standalone HTML」「html-report-build」、weekly/structured-tech/results-report 結束自動觸發、validated 報告需 PI 終版閱讀體驗。
  standalone PI 終版前必繼承 InterSubMod/.claude/skills/scientific-rigor/SKILL.md §2-§7（claim card 必有 tier badge、§8.4 provenance footer 必含 commit hash）。
  SKIP WHEN: README.md (GitHub 原生 render)、給 LLM 消費的 .md (CLAUDE.md / *.json / MEMORY.md)、純筆記 / CI 文件、in-progress 草稿（用 report companion）。
---

# html-report-build

## Phase & Chain Position

- **Phase 3** of HTML preview + image-gen 3-skill design (supersedes Phase 2 html-preview).
- Stand-alone, callable directly or via myPPT / weekly-report / structured-tech-report / pptx-build orchestration.
- **Upstream**: any .md report or in-context content + design tokens + few-shot examples.
- **Downstream**: PI / collaborators open `.html` in browser; companion to source `.md` (D11: never replaces .md).
- **Replaces**: `html-preview` (deprecated 2026-05-13, Python middleware → LLM-direct).

## Dependencies

**Uses (runtime)**:
- **NO Python**, **NO markdown lib**, **NO jinja2**, **NO build step**.
- Only: Claude LLM reading `.md` + applying `prompts/*.md` + injecting `templates/design_tokens.css`.
- Optional: `base64` CLI for PNG → inline (if image inlining needed; LLM cannot read binary directly).

**Used by**:
- 6 Tier A skills (auto-trigger on Stop hook): `weekly-report`, `structured-tech-report`, `results-report`, `feature-layered-observation`, `methodology-audit`, `conclude-research`.
- `pptx-build` for slide HTML preview (replaces self_phasing-style custom `build_html.py`).
- `myPPT` orchestration (Phase 3).

**Reads**:
- Source `.md` file (full text via Read tool).
- Local PNG referenced in `.md` (resolved relative to .md dir; need base64 for inline).
- `templates/design_tokens.css` (copy into output `<style>` block).
- `prompts/report_prompt.md` or `prompts/slide_prompt.md` depending on mode.
- `examples/*.html` for few-shot grounding.

**Writes**:
- **Report mode (single-file)**: `{md_basename}.html` next to source `.md`.
- **Report mode (topic folder, >=200 lines OR >=5 figures)**: `{md_basename}/` containing `index.html` + `ch{NN}_{slug}.html` + `README.md`.
- **Slide mode**: `preview/index.html` + `preview/slide_XX_{slug}.html` × N + `preview/shared/style.css` (optional).
- **NEVER touches source `.md`** (D11 companion rule).

## Failure Modes & Diagnostics

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Output HTML > 1 MB | Many large PNGs inlined as base64 | Skip inline for large PNG; use external `<img src="../figures/...">` |
| Visual diff > 5% from reference | Prompt missing few-shot example or design tokens | Re-read `examples/slide_pair_example.html`; ensure `:root` CSS variables copied |
| 6-taboo audit fail (gradient / glass / etc.) | LLM hallucinated decorative style | Re-prompt: "Strictly use design_tokens.css palette; no gradients except `.conclusion-arrow`" |
| Speaker note missing | Slide mode skipped `<details class="note-section">` | Re-prompt with `components/speaker_note.html` as required structure |
| 16:9 aspect broken | `aspect-ratio: 16/9` missing on `.slide-canvas` | Verify `<article class="slide-canvas" data-section="...">` exact class names |
| `data-section` color not applied | Section attribute missing or typo | Check S0-S7 + Q&A enum; ensure `<article ... data-section="S2">` exact match |
| Companion replaced .md | Wrong filename collision | Output **MUST** be `{basename}.html`, never `{basename}.md` rewrite |
| Multi-language CJK font fallback | System font missing on PI machine | Use full font stack: Noto Sans CJK TC + Source Han Sans TW + Droid Sans Fallback |
| **Figure broken in browser** ⚠ **2026-05-20 Issue #1** | Relative path `../` off-by-one (LLM 易 miscount) | Before Write: count depth_of(source_dir to repo_root) and verify `../` 數量等於該 depth。grep `<img src="..">` + post-Write `ls path` verify exists |
| **Fabricated metric** ⚠ **2026-05-20 Issue #2** | LLM 內插「合理範圍」數字而非嚴格 source-grep | Before Write 每張 slide: list 所有 numerical values + grep source .md confirm verbatim 或 rounded ≤ 2 sf；多版迭代鎖 source set 不可新增數字 |
| **Slide chars/rows over scenario limit** ⚠ **2026-05-20 Issue #4** | 缺場景對齊量化標準（PI 1-on-1 / Lab meeting / Conf）| Read frontmatter `audience_scenario` 套 memory `reference-pi-scenario-quantitative-standards`；超限自動 flag for compression |
| **Path prefix not InterSubMod/** ⚠ **2026-05-20 Issue #9** | LLM 生成 footer source 引用未自審路徑前綴 | grep `href=` and source citations: any `.md` path must start with `InterSubMod/...` not absolute `/big7_disk/...` or relative `docs/...` |

## Files Manifest (self-check before invocation)

The skill is incomplete without these 23 files. Verify presence via `ls .claude/skills/html-report-build/`:

```
.claude/skills/html-report-build/
├── SKILL.md                                  ← this file
├── evals.json                                ← 6 test cases
├── prompts/
│   ├── report_prompt.md                      ← 8 rules for report mode
│   ├── slide_prompt.md                       ← 8 rules for slide mode
│   └── standalone_prompt.md                  ← 10 rules for standalone PI-end mode (incl. SVG)
├── templates/
│   ├── design_tokens.css                     ← :root CSS variables
│   ├── report_skeleton.html                  ← prose wrapper
│   ├── slide_skeleton.html                   ← 16:9 canvas wrapper
│   ├── slide_index_skeleton.html             ← topbar + tab nav + iframe
│   └── standalone_skeleton.html              ← sticky TOC + main + section cards layout
├── components/
│   ├── callout_critical.html                 ← red-frame warning
│   ├── stat_box.html                         ← big-number visual
│   ├── metric_table.html                     ← comparison table
│   ├── igv_zoom_panel.html                   ← PanZoom IGV
│   ├── speaker_note.html                     ← required per slide
│   ├── conclusion_arrow.html                 ← gradient verdict
│   ├── svg_flow_diagram.html                 ← inline SVG causal chain (standalone)
│   ├── svg_compare_bar.html                  ← inline SVG proportion bar (standalone)
│   └── svg_icon_set.html                     ← inline SVG icons (replaces emoji)
├── examples/
│   ├── README.md                             ← pattern lookup + redirects
│   ├── report_short_example.html             ← report mode few-shot
│   └── slide_pair_example.html               ← slide mode pointer to legacy
└── references/
    └── design_principles.md                  ← 12 rules from Tufte/Reynolds/Duarte/CRAP/AE
```

If any file missing → skill setup is corrupt; abort and report to user.

## Quick Usage

This skill is **prompt-based, not script-based**. To invoke:

1. **Manual**: User types `/html-report-build <md_path>` or asks "convert this .md to HTML".
2. **Auto (Stop hook)**: Tier A skill finishes → hook detects new .md → triggers this skill.
3. **Direct via Skill tool**: `Skill html-report-build` with args `--mode report|slide <path>`.

Workflow:
```
1. Read source .md (full content)
2. Decide mode: report (default) or slide (if .md is presentation outline)
3. Read prompts/{report,slide}_prompt.md for rule set
4. Read templates/design_tokens.css → inline into output <style>
5. Read 1-2 examples/*.html for few-shot grounding
6. Generate output .html using Write tool
7. **Self-audit suite (2026-05-20 Issue #1/#2/#4/#7/#9 升級)**:
   a) 6-taboo grep + structural checks (see Failure Modes table)
   b) **Number-source-grep audit**: list all numerical values in slides; grep source .md to confirm verbatim or rounded ≤ 2 sf. ANY un-sourced number → FAIL, ask user OR use range.
   c) **Path off-by-one audit**: count source_dir → repo_root depth; verify all `<img src="../...">` 數量符合 depth. Post-Write `ls` each path to confirm exists.
   d) **Path prefix audit**: any `.md` path reference in footer/source must start with `InterSubMod/...`.
   e) **Scenario chars/rows audit**: read frontmatter `audience_scenario`; per-slide count chars (excl speaker note) + table rows; flag any exceeding scenario threshold (memory `reference-pi-scenario-quantitative-standards`).
8. Report path to user + audit summary
```

## Multi-Version Output (2026-05-20 Issue #6 升級)

When user requests N versions for comparison (e.g. "produce 3 styles"):

1. **Define N differentiation axes explicitly** before any Write (e.g. minimal / standard / dense / hybrid)
2. Each version picks ONE axis — NOT N=N independent designs
3. **Common source-of-truth lock**: shared metadata block (cycle_id + commit + source .md lines) at top of each HTML; numbers across versions must come from this locked set
4. After all versions: produce 1-page comparison `index.html` with per-version use-case + slide-count + chars-avg
5. **Multi-version fabrication prevention**: when iterating version N+1, **MUST grep version N for all numerical values** and reuse the exact set; new numbers require new source citation

## Evaluator Polish Loop (2026-05-20 Issue #10 升級)

After evaluator audit returns `polish notes`:

- If user selects `ship as-is` → log notes to `known-pitfalls.md` + skill MEMORY.md as future prevention guidance
- If user selects `apply polish` → batch edit + re-run evaluator to confirm no regression
- Do NOT silently drop polish notes — they encode reusable lessons

## Mode Selection (concrete signatures)

| Source signature | Mode | Output |
|------------------|------|--------|
| Plain prose `.md`: no `## Slide` headings, in-progress / draft | **report** | `{basename}.html` (single) or `{basename}/` folder |
| `.md` in `docs/reports/validated/` AND `status: validated` | **standalone** | `{basename}.standalone.html` (PI-end view) |
| `.md` in `docs/experiments/concluded/` AND tier ⭐4-5 | **standalone** | same |
| Source is PI report / errata / verdict document | **standalone** | same |
| Frontmatter `audience: PI / advisor / collaborator` | **standalone** | same |
| User explicitly says「終版」/「給 PI」/「single file」/「standalone」 | **standalone** | same |
| `.md` has `## Slide NN` or `### slide_XX` headings | **slide** | `{basename}/preview/index.html` + slide files |
| `.md` frontmatter has `slide_count:` or `data-section:` keys | **slide** | same as above |
| `.md` has `speaker_note:` / `timing:` repeated per section | **slide** | same as above |
| Filename contains `_outline` / `_slides` / `_deck` / `_presentation` | **slide** (suggestive) | same; ask user if mixed signals |
| Mixed signals (some narrative, some slide-shaped) | **ask user** | Stop and clarify; do not improvise |

**Threshold for topic-folder (report mode)**: `>=200 lines` OR `>=5 figures` → folder; else single-file. See `prompts/report_prompt.md` Rule 4 for concrete path examples + CJK slug derivation rule.

**Mode coexistence**: All three modes can output for the same `.md` if user wants both grep-view (`.md`) and PI-view (`.standalone.html`). Filenames distinguish:
- `{basename}.md` (source, git)
- `{basename}.html` (report mode, prose companion)
- `{basename}/` (topic-folder mode)
- `{basename}.standalone.html` (PI-end terminal view, with sticky TOC + SVG + cards)

## Path-based Routing Table (`USE WHEN` 分流)

This table defines **which mode to invoke for which `.md` path**. Used by `scripts/hooks/standalone_trigger.sh` (PostToolUse Edit|Write) to auto-suggest the right mode.

| `.md` path | Mode | 為什麼 |
|------------|------|--------|
| `InterSubMod/docs/reports/validated/{YYYY}/{MM}/*.md` | **standalone** ★ | Validated reports — PI / advisor 終版閱讀 |
| `InterSubMod/docs/reports/research_landscape/*.md` | **standalone** ★ | 主軸大圖景 — PI 反覆讀，sticky TOC 必要 |
| `InterSubMod/docs/reports/pi_reports/*.md` | **standalone** ★ | 給 PI 終版 — 顧名思義 |
| `InterSubMod/docs/experiments/concluded/{YYYY}/{MM}/*.md` | **standalone** ★ | tier ⭐4-5 收尾，已不再 diff |
| `InterSubMod/docs/concepts/*.md` | **standalone** ★ | 技術說明文章 — SVG flow + sticky TOC 大受益 |
| `InterSubMod/docs/data_specs/*.md` | **standalone** ★ | 規格類技術文章 — 多人來回查 |
| `InterSubMod/docs/references/manual/*.md` | **standalone** ★ | Handbook — 長期參考 |
| `InterSubMod/docs/experiments/in_progress/{YYYY}/{MM}/*.md` | **report** (companion) | 草稿，仍會多次 diff |
| `InterSubMod/docs/plans/*.md` | **report** (companion) | 實驗計畫草稿 |
| `InterSubMod/docs/presentations/*/02_slide_outline.md` | **slide** | 明確的 slide outline |
| `InterSubMod/docs/presentations/*/03_slide_layout_script.md` | **slide** | Layout 腳本 |
| `**/CLAUDE.md` / `**/AGENTS.md` / `**/MEMORY.md` / `**/README.md` | **SKIP** ⛔ | LLM 消費 / GitHub 原生渲染 |
| `state/*.json`, `**/*.jsonl`, `**/*.yaml`, `**/*.csv`, `**/*.tsv` | **SKIP** ⛔ | 結構化資料，非報告 |
| `InterSubMod/docs/presentations/.../self_phasing_synthesis_PI/preview/*` | **SKIP** ⛔ | Legacy custom PPT，凍結保留 |

★ = Auto-suggested by `standalone_trigger.sh` hook on Edit|Write

## Automatic Trigger via Hooks

This skill is auto-suggested by **PostToolUse on Edit|Write** via `InterSubMod/scripts/hooks/standalone_trigger.sh`. The hook:

1. Reads the file path being written/edited.
2. Matches against the Routing Table above.
3. If match → prints `[html-report-build]` reminder lines (3 lines: file modified / suggested mode / output path).
4. If SKIP → exits silently (no noise).
5. If `.standalone.html` already exists AND is newer than `.md` → exits silently (already up-to-date).

**Claude reads the reminder lines** and decides whether to invoke this skill. Hook does NOT auto-execute the skill — it only suggests.

To **disable auto-suggestion** for a specific report (e.g., draft phase): rename source to `*.draft.md` until ready.

To **manually invoke** without the hook: `Skill html-report-build` with `--mode standalone <path>`.

Hook registration: `.claude/settings.local.json` → `hooks.PostToolUse` → entry with `command: "bash .../standalone_trigger.sh"`.

Related hooks already in place:
- `trigger_routing.sh` — general C++/CMake/skill/CLAUDE.md routing reminders
- `kb_sot_guard.sh` — knowledge base source-of-truth guard
- `md_link_check.sh` — markdown cross-link validation
- `evidence_ledger_sync.sh` — evidence ledger entry reminder

These cooperate: when you Edit a validated report, you'll see (a) routing reminder, (b) link check, (c) evidence ledger reminder, (d) **standalone HTML suggestion** — all from PostToolUse hooks, all advisory.

## SVG Diagrams (standalone mode capability)

Standalone mode supports **inline SVG** for small diagrams (Rule 4 in `prompts/standalone_prompt.md`):

| Diagram type | Use SVG? | Component reference |
|--------------|----------|---------------------|
| Flow / causal chain (≤10 nodes) | ✅ YES | `components/svg_flow_diagram.html` |
| Proportion bar / stacked comparison | ✅ YES | `components/svg_compare_bar.html` |
| Icon (replaces emoji-in-heading taboo) | ✅ YES | `components/svg_icon_set.html` |
| Architecture / sparkline / small chart | ✅ YES | hand-write following few-shot |
| Real data scatter (>50 points) | ❌ NO | use matplotlib → PNG → base64 inline |
| IGV / photograph | ❌ NO | PNG only |
| Interactive (zoom / hover tooltip) | ❌ NO | out of scope; use D3 elsewhere or skip |

**SVG hard rules**: `viewBox` + `role="img"` + `<title>` + `<desc>` mandatory. Colors via design tokens (CSS vars). NO animations, NO filters, NO gradients (taboo). Max 5 SVG per report; max viewBox 800×400.

Report mode + slide mode may also use small SVG (icons + arrows), but standalone mode is the primary SVG consumer.

## Companion Rule (D11)

**NEVER** rewrite or delete the source `.md`. Output `.html` lives alongside as companion. Both go into git:
- `.md` = source of truth (human-readable diff, grep-able, INDEX referenced)
- `.html` = visual companion (PI / collaborator browser view, print-ready)

## Design Tokens Reference

All CSS variables defined in `templates/design_tokens.css`:
- **Palette**: `--c-accent: #D97757` / `--c-text: #141413` / `--c-bg: #FAF9F5` / `--c-border: #E3DACC`
- **Spacing**: 8-point scale `--sp-1: 4px` → `--sp-9: 96px`
- **Typography**: system-ui + Noto Sans CJK TC + JetBrains Mono fallbacks
- **Print mode**: `@media print` expands all `<details>`, prints URLs after links, removes shadow

**Inline rule**: Copy `:root { ... }` block from `templates/design_tokens.css` into every output `<style>`. Do NOT link external `.css`.

## 6-Taboo Audit (run on every output)

Reject output that contains:
1. Multiple gradients (>1 occurrence of `linear-gradient` / `radial-gradient`; one allowed for `.conclusion-arrow`)
2. Glass morphism (`backdrop-filter: blur`)
3. Multi-indigo (>2 `#1E3A8A`-family shades; use accent palette instead)
4. Emoji-in-h1/h2/h3 (emoji acceptable inline, not as heading prefix)
5. `text-shadow` (any usage outside `.conclusion-arrow.green`)
6. `box-shadow` with `0 0 *` glow (only use directional shadows `--shadow-sm/md/lg`)

Audit command (manual grep):
```bash
grep -E '(linear-gradient|radial-gradient|backdrop-filter|text-shadow|0 0 [0-9]+px)' output.html
```

## Design Principles (canonical reference)

**Before generating ANY HTML, Read `references/design_principles.md`** — 12 rules synthesized from industry sources:

| # | Rule | Source |
|---|------|--------|
| 1 | 5-second test (first viewer 抓到 takeaway) | NN/g + Berkeley dataviz checklist |
| 2 | 3-second glance (slide is billboard) | Nancy Duarte《slide:ology》 |
| 3 | Assertion-Evidence (title = claim, body = evidence) | Penn State / NSF-backed |
| 4 | Data-ink ratio (remove chartjunk) | Edward Tufte |
| 5 | CRAP (Contrast / Repetition / Alignment / Proximity) | Robin Williams 1993 |
| 6 | Whitespace as design (留白 = signal) | Garr Reynolds Presentation Zen |
| 7 | Hierarchy (重要元素最強視覺權重) | NN/g + Berkeley BPM |
| 8 | 1-2 primary + 1 accent (色彩濫用 = 噪音) | dataviz BP 2026 |
| 9 | Colorblind-safe (pattern + label, not color only) | Nature 2026 + WCAG |
| 10 | WCAG contrast ≥ 4.5:1 | WCAG 2.1 |
| 11 | Vector first (SVG > PNG for diagrams) | Nature / Cell figure guideline |
| 12 | Pre-publish 6-item checklist (must all ✓) | Synthesis |

**The Two Audit Loops** (use both):
1. **Positive** — `references/design_principles.md` Rule 1-12 (how to do it right)
2. **Negative** — `SKILL.md §6-Taboo Audit` (what not to do)

Both must pass before Write. If a rule conflicts (e.g., Rule 8 "1-2 primary" vs slide 9-section coloring) — slide colors are **semantic encoding** (section identity), not decorative; the spirit is preserved.

## 嚴謹度繼承（/scientific-rigor）

HTML 報告（特別 standalone PI 終版）必繼承 `InterSubMod/.claude/skills/scientific-rigor/SKILL.md`：

- **§2 證據分級**: 每個 claim card 必標 ⭐⭐⭐⭐⭐ tier badge（HTML badge style: `<span class="tier-l2">⭐⭐⭐⭐ L2</span>`）
- **§3 Effect Size**: 報告內 metric table 必含 Cohen ribbon + CI 欄
- **§4 DAG**: standalone 模式必嵌 SVG DAG 圖（mermaid CDN 或 inline SVG）
- **§7 Pre-registration**: validated 報告必含「Pre-registration 對照表」欄
- **§8.4 Provenance**: standalone 報告必有 footer 含 commit hash + cycle_id + 生成時間

**最小可用子集**:
- standalone PI 終版: §2 + §3 + §4 + §7 + §8.4 全跑
- report 一般技術報告: §2 + §3
- slide HTML preview: §2 + §3 + slide-level ribbon

**6-Taboo Audit 擴展**: §6-Taboo 加第 7 條「Evidence Tier missing」— 任何 claim card 無 tier badge → 紅旗。

## See Also

- **Spec**: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md` (D11 / D15-D20)
- **SOP**: `InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md` (§3 排除 / §5 模板 / §6 evaluator-optimizer)
- **Plan** (outside InterSubMod repo, design rationale only): `~/.claude/plans/frolicking-tinkering-hopcroft.md` — not required at runtime; reference for design lineage only
- **Deprecated predecessor**: `html-preview` (Python middleware version; **git rm'd 2026-05-30**, see commit history)
- **Legacy reference**: `InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/build_html.py` (PPT custom generator, kept frozen)
- **Companion skills**: `image-gen`, `image-vision-check` (Phase 1, unchanged)
