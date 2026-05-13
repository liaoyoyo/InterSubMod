<!--
build_date: 2026-05-13
audience: PI / lab member / 自己未來
report_class: skill-system-documentation
status: validated
verdict: html-report-build skill 系統完整 shipped — 3 mode (report / slide / standalone) + LLM-direct + 12 design principles + hook auto-trigger + 23 files / 200 KB
last_verified: 2026-05-13
-->

# html-report-build Skill 系統文件 — 2026-05-13

## 0. TL;DR — 3 句話

1. **`html-report-build` skill 取代 `html-preview`**：零 Python 中介，LLM-direct 生成 HTML，遵循 Anthropic Thariq 2026-05 主張的「介面語言從 Markdown 轉向 HTML」。
2. **3 mode 涵蓋 95% 場景**：`report` (in-progress 草稿 companion) / `slide` (PPT deck) / `standalone` (PI 終版單檔 + sticky TOC + inline SVG + 卡片化)。
3. **5 個機制聯動自動化**：path routing → standalone_trigger hook → design principles → skill prompts → output audit。

## 1. 系統組件總圖

| 組件 | 位置 | 功能 |
|------|------|------|
| **新 skill** | `InterSubMod/.claude/skills/html-report-build/` | 3 mode HTML 生成（23 檔案） |
| **舊 skill (deprecated)** | `InterSubMod/.claude/skills/html-preview/` | 標 DEPRECATED；過渡保留 |
| **Hook 自動觸發** | `InterSubMod/scripts/hooks/standalone_trigger.sh` | PostToolUse Edit/Write 偵測路徑提示 |
| **Settings 註冊** | `InterSubMod/.claude/settings.local.json` | PostToolUse 第 10 個 hook entry |
| **Design Principles** | `InterSubMod/.claude/skills/html-report-build/references/design_principles.md` | 12 條業界準則 + 5 秒測試 + 6-item checklist |
| **Memory canonical** | `feedback_design_principles_canonical.md` | 跨 session 記憶（MEMORY.md 索引） |

## 2. 3 個 prompt mode 對照

| Mode | 觸發路徑 | 輸出 | 特色 |
|------|---------|------|------|
| **report** | `docs/experiments/in_progress/` / 一般 .md | `{basename}.html` 或 `{basename}/` folder | Tailwind prose 風格，輕量 |
| **slide** | `docs/presentations/*/02_slide_outline.md` | `preview/index.html` + slide_XX.html × N | 16:9 canvas + 9 section colors + speaker note |
| **standalone** ★ | `docs/reports/validated/` / `pi_reports/` / `concepts/` / `data_specs/` / `references/manual/` | `{basename}.standalone.html` | sticky TOC + section cards + inline SVG + 列印友善 |

## 3. 12 設計準則（業界 source）

| # | Rule | Source |
|---|------|--------|
| 1 | 5 秒測試 | NN/g + Berkeley dataviz checklist |
| 2 | 3 秒 glance / billboard | Nancy Duarte《slide:ology》 |
| 3 | Assertion-Evidence（標題=結論句）| Penn State, NSF-backed |
| 4 | Data-ink ratio | Edward Tufte《Visual Display》 |
| 5 | CRAP (Contrast/Repetition/Alignment/Proximity) | Robin Williams 1993 |
| 6 | Whitespace as design | Garr Reynolds Presentation Zen |
| 7 | Hierarchy（最重要元素最強權重）| NN/g |
| 8 | 1-2 primary + 1 accent | dataviz BP 2026 |
| 9 | Colorblind-safe（pattern + label 雙重編碼）| Nature 2026 + WCAG |
| 10 | WCAG contrast ≥ 4.5:1 | WCAG 2.1 |
| 11 | Vector first（SVG > PNG）| Nature / Cell figure guideline |
| 12 | Pre-publish 6-item checklist | Synthesis |

## 4. SVG 啟用規則

### ✅ 用 inline SVG（≤5 per report）
- 流程圖 / 因果鏈 (≤10 nodes) — `svg_flow_diagram.html`
- 比例條 / 對比 bar — `svg_compare_bar.html`
- Icon (替代 emoji-in-heading 禁忌) — `svg_icon_set.html`

### ❌ 用 PNG（base64 inline 或 external src）
- 真實數據 scatter (>50 points) → matplotlib
- 連續分佈 histogram → matplotlib
- IGV 截圖 / 照片
- 互動 chart → 超出範圍

### SVG 硬規則
- `viewBox` + `role="img"` + `<title>` + `<desc>` 必含
- 顏色用 CSS variables 不寫死 hex
- NO `<animate>` / NO `<filter>` / NO SVG `<linearGradient>`
- viewBox max 800×400

## 5. Auto-trigger 機制

```
1. 用戶 Edit/Write docs/reports/validated/.../foo.md
   ↓
2. PostToolUse hook 觸發 standalone_trigger.sh
   ↓
3. Hook 檢查路徑：
   - SKIP: README.md / CLAUDE.md / state.json / in_progress 草稿
   - SKIP: 若 .standalone.html 已比 .md 新
   - TRIGGER: validated reports / concepts / data_specs / manual
   ↓
4. Hook 印 3 行 [html-report-build] reminder
   ↓
5. Claude 看到 reminder，呼叫 Skill html-report-build
   ↓
6. Skill Read：
   - source .md
   - references/design_principles.md (12 rules)
   - prompts/standalone_prompt.md (10 rules)
   - templates/standalone_skeleton.html
   - components/svg_*.html few-shot
   ↓
7. 生成 {basename}.standalone.html (sticky TOC + cards + SVG)
   ↓
8. 6-Taboo audit + 12 rules self-check → ship
```

## 6. 兩個 pilot 成果

| Pilot | Source | Output | 大小 |
|-------|--------|--------|------|
| Errata standalone | `InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md` (281 行) | `InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.standalone.html` | 44 KB / 724 行 / 4 SVG |
| Synthesis standalone | `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md` (1451 行) | `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.standalone.html` | 59 KB / 931 行 / 4 SVG |
| Slide pilot | `InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/02_slide_outline.md` | `InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/preview_v2/` (4 files) | 38 KB / 808 行 |

## 7. Skill 23 個檔案清單

```
InterSubMod/.claude/skills/html-report-build/
├── SKILL.md                          # 269 lines
├── evals.json                        # 7 test cases
├── prompts/
│   ├── report_prompt.md              # 報告 mode (8+ rules)
│   ├── slide_prompt.md               # PPT mode (8+ rules)
│   └── standalone_prompt.md          # PI 終版 mode (10+ rules, SVG)
├── templates/
│   ├── design_tokens.css             # :root variables
│   ├── report_skeleton.html
│   ├── slide_skeleton.html
│   ├── slide_index_skeleton.html
│   └── standalone_skeleton.html      # sticky TOC + cards layout
├── components/
│   ├── callout_critical.html
│   ├── stat_box.html
│   ├── metric_table.html
│   ├── igv_zoom_panel.html
│   ├── speaker_note.html
│   ├── conclusion_arrow.html
│   ├── svg_flow_diagram.html         # 因果鏈 + 兩層 bug
│   ├── svg_compare_bar.html          # 比例對比
│   └── svg_icon_set.html             # 8 icons (warning/check/info/...)
├── examples/
│   ├── README.md
│   ├── report_short_example.html
│   └── slide_pair_example.html
└── references/
    └── design_principles.md          # 12 rules canonical
```

## 8. 業界 sources

- Anthropic Thariq 2026-05：Claude Code 預設輸出從 MD 切 HTML
- Anthropic Claude Artifacts pattern：self-contained 單檔 HTML
- Edward Tufte《Visual Display of Quantitative Information》
- Garr Reynolds《Presentation Zen》系列
- Nancy Duarte《slide:ology》/《Resonate》
- Robin Williams《Non-Designer's Design Book》(CRAP 1993)
- Michael Alley Assertion-Evidence (Penn State, NSF-backed)
- Nature / Cell / Science figure guidelines (2025-2026)
- WCAG 2.1 contrast
- Nielsen Norman Group dataviz heuristics

## 9. 引用文件

- `InterSubMod/.claude/skills/html-report-build/SKILL.md` — 主 skill 入口
- `InterSubMod/.claude/skills/html-report-build/references/design_principles.md` — 12 規則
- `InterSubMod/.claude/skills/html-preview/SKILL.md` — deprecated skill (referenced for migration)
- `InterSubMod/scripts/hooks/standalone_trigger.sh` — auto-trigger hook
- `InterSubMod/.claude/settings.local.json` — hook registration

## 10. 未來方向

- **Stage 2 (Tier A skill rewrite)**：6 Tier A skill 的 Stop hook 切到 html-report-build
- **Stage 3 (2 週 sunset)**：刪 html-preview Python tools
- **跨樣本 standalone**：對既有 7 個 validated reports 補產 standalone HTML
- **Skill-reviewer agent 整合**：產出後自動跑 plugin-dev:skill-reviewer 評分
