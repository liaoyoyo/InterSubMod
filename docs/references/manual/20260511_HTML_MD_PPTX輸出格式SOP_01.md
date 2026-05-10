---
title: "HTML / Markdown / PPTX 輸出格式 SOP — agent skill 設計與營運手冊"
date: 2026-05-11
status: active
version: 1.0.0
authors: [liaoyoyo2001 (PI), Claude Opus 4.7 (assistant)]
supersedes_partial:
  - InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md (D3, D4, D12 修正)
related:
  - InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md
  - InterSubMod/.claude/skills/image-gen/SKILL.md
  - InterSubMod/.claude/skills/image-vision-check/SKILL.md
  - InterSubMod/.claude/CLAUDE.md
references_external:
  - https://thariqs.github.io/html-effectiveness/
  - https://platform.claude.com/docs/en/agents-and-tools/agent-skills/best-practices
  - https://www.anthropic.com/research/building-effective-agents
  - https://andybrewer.github.io/mvp/
  - https://github.com/anthropics/anthropic-cookbook/tree/main/patterns/agents
tags: [sop, skill_design, html_preview, output_format, agent_workflow]
---

# HTML / Markdown / PPTX 輸出格式 SOP

> **Scope**: 規範 InterSubMod 專案內所有 agent / skill / human 產出文件的輸出格式選擇、技術棧、SKILL.md 撰寫規範、antipattern 防範、evaluator-optimizer 命名。
>
> **Status**: 本 SOP 1.0 取代 2026-05-10 design doc 的 D3 / D4 / D12 部分決策（其餘 D1, D2, D5–D14 仍有效）。新增 D15–D20。

## TL;DR — 決策速查

| 場景 | 用什麼 | 理由 |
|------|--------|------|
| README / repo 入口 | **Markdown** | GitHub 原生渲染、grep 友善 |
| 給 LLM 消費的文件（CLAUDE.md / state.json / hypothesis_queue） | **Markdown / JSON** | 模型解析效率高 |
| Git history-heavy（INDEX.md、MEMORY.md） | **Markdown** | 細粒度 diff |
| 個人 AI session 紀錄 | **Markdown** | 量大、自己看 |
| **PI 報告（status / 週報 / 技術報告 / 實驗結果）** | **MD source + HTML companion** | 視覺化、可分享、可列印 |
| **研究 archive（concluded research）** | **MD source + HTML companion** | 永續可讀、自包含 |
| **教授級正式報告 / 對外簡報** | **PPTX** | 學術/商業慣例 |
| **互動式探索（playground / editor）** | **HTML 單檔** | 無法用其他格式達成 |

---

## 1. 新增決策（D15–D20，補充 2026-05-10 design doc）

| # | 項目 | 決策 |
|---|------|------|
| **D15** | **HTML 技術 stack 改 MVP.css 預設**（取代 Tailwind L2-bake） | 預設用 [MVP.css](https://andybrewer.github.io/mvp/) 8 KB classless framework + 內聯 design tokens；Tailwind 改為 opt-in (`--tailwind`) 對「設計系統文件 / 元件庫」場景使用。理由：thariqs 實際範例（11/14/16）均極簡風格；MVP.css 對齊「語意 HTML 即可」哲學；省 build step；省 17 KB |
| **D16** | **Design tokens 採 Anthropic 主色 + thariqs 範例 spacing** | Palette: `#D97757`（accent）、`#141413`（text）、`#FAF9F5`（bg）、`#E3DACC`（border）；Spacing: 8-point scale (`--sp-1`=4px → `--sp-8`=64px)；Typography: system-ui + 6 級 size scale（48/32/24/16/14/12 px）。**全部以 CSS variables 內聯**，不外部檔案。 |
| **D17** | **主題資料夾條件化** | 報告 ≥200 lines OR ≥5 figures 才建主題資料夾；其餘預設單檔 `{name}.preview.html`（修正 D4 的「永遠主題資料夾」）。理由：thariqs 20 範例中 19/20 為單檔；資料夾增加認知負擔。 |
| **D18** | **Evaluator-Optimizer 模式正式命名** | 我們的 image-gen + image-vision-check 為此 pattern 實例。新增 reference doc 標準化 pattern：generator → evaluator → 反饋 → max 1-2 retry。命名後可重用於其他驗證迴圈（如 PPTX vision check / 報告品質 check）。 |
| **D19** | **frontend-design plugin opt-in 整合** | `/html-preview --polished` flag 觸發 `frontend-design:frontend-design` plugin 生「精緻設計」HTML；預設不啟動（極簡 MVP.css 即可）。 |
| **D20** | **既有 22 skill 全補 evals**（3/skill = 66 evals total）+ description audit | 對齊 Anthropic best-practices「evals-driven development」。Phase 3 期間補完，每 skill 加 `tests/evals.json`（3 case：典型觸發 / SKIP scenario / edge case）。Description 同步 audit 為 directive form。 |

---

## 2. 何時用 HTML — Anthropic 9 category（thariqs）

引用 [thariqs/html-effectiveness](https://thariqs.github.io/html-effectiveness/) + 對應本專案場景：

| 類別 | 對應本專案場景 | Tier |
|------|-------------|------|
| Exploration & Planning | 候選假說比較、3 方案 mockup | A |
| Code Review & PR | gh PR diff（已用 gh CLI） | △ |
| Design & Prototypes | （不太用） | C |
| Diagrams & Illustrations | image-gen 4 類圖（concept/flow/data/icon） | A |
| Reports & Research | weekly-report / structured-tech-report / results-report / feature-layered-observation / methodology-audit / conclude-research | **A++** |
| Decks | 已由 pptx-build 處理 | △ |
| Custom Editors | 假說 reorder / hypothesis_queue 編輯器 | B |
| Matching Style | 設計品味（design tokens） | A |
| Concept Explainers | LOH / phasing / Inter-Subclonal 概念互動 | A |

## 3. 何時保留 Markdown — 6 排除場景（thariqs）

| 場景 | 對應本專案 |
|------|----------|
| README.md | repo 根 + 子目錄 README，**永遠 MD** |
| Slack/Discord 片段 | 我們不用 Slack |
| 給其他 LLM 消費 | `CLAUDE.md` / `AGENTS.md` / `CURRENT_FOCUS.md` / `MEMORY.md` / `state/cycles/*/state.json` |
| Git history-heavy | `docs/experiments/INDEX.md` / `evidence_ledger.jsonl` |
| 個人備忘錄 | `docs/provenance/ai_sessions/*.md`（report skill 輸出） |
| RSS / 電子報 | 不用 |

---

## 4. SKILL.md 撰寫規範（per Anthropic best-practices）

### 4.1 Frontmatter 必需欄位

```yaml
---
name: <gerund-form>           # 新 skill 用 gerund (processing-X / analyzing-Y); 舊 skill 維持原名
description: |                # 必需; max 1024 chars; third-person; directive form
  <verb> <subject> <output>. Use when <trigger phrases>. Skip when <exclusions>.
---
```

### 4.2 Description 寫法 — Directive form（D20 audit 標準）

**❌ 描述式（弱觸發）**：
> "Helps with documents"
> "I can help you process Excel files"
> "用於整理週報"

**✅ Directive form（強觸發）**：
> "Generate descriptive commit messages by analyzing git diffs. Use when the user asks for help writing commit messages or reviewing staged changes."
> "Audit AI-generated images against 6-dimension checklist; output quality.json with verdict pass/partial/fail. Use when image-gen produces output OR PI asks 圖夠不夠好/這張圖能用嗎."

**結構模板**：
```
<Verb 動詞> <受詞> by <方法>. 
USE WHEN: <具體觸發詞 1>, <觸發詞 2>, <觸發詞 3>, <情境 1>, <情境 2>.
SKIP WHEN: <排除場景 1>, <排除場景 2>.
```

### 4.3 Body 規範

| 規則 | 數值 |
|------|------|
| SKILL.md body 上限 | **500 行** |
| Reference 檔案深度 | **1 level deep**（不允許 chained refs） |
| 子檔案 TOC 觸發 | **>100 行** 必加 TOC |
| 命名 | **forward-slash paths**（不允 Windows `\`） |
| Terminology | 一致（不混用同義詞） |
| Time-sensitive info | 移到 `## Old patterns` 折疊區 |

### 4.4 Progressive disclosure 3 種 pattern（依複雜度選）

**Pattern 1**：高層 guide + 平行 references（適合 1 大主題 + 4-5 子主題）
```
SKILL.md → reference/finance.md / sales.md / product.md / marketing.md
```

**Pattern 2**：domain-specific organization（適合多 domain skill）

**Pattern 3**：conditional details（適合「90% 簡單 + 10% 進階」場景）

### 4.5 Workflow 寫法

複雜 skill 必含 checklist 讓 Claude 抄寫追蹤：
```markdown
## Workflow

Copy this checklist:
- [ ] Step 1: <action>
- [ ] Step 2: <action>
- [ ] Step 3: <action>

**Step 1**: <details>
...
```

### 4.6 Feedback loop（evaluator-optimizer 模式）

```markdown
## Validation loop
1. Generate output
2. Run validator (script or Claude self-review)
3. If fail: read suggestions → patch → re-generate
4. Max 1-2 rounds (per Anthropic playbook)
5. Final fail → escalate to PI with checklist
```

---

## 5. HTML Preview 技術規範（D15 + D16 落實）

### 5.1 預設層級：MVP.css 極簡

**完整單檔 HTML 模板**（複製即用，無依賴）：

```html
<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  /* === Design Tokens (D16) === */
  :root {
    --color-accent: #D97757;
    --color-text: #141413;
    --color-bg: #FAF9F5;
    --color-border: #E3DACC;
    --color-muted: #87867F;
    --color-success: #788C5D;
    --color-warning: #C78E3F;
    --color-danger: #B04A4A;
    --sp-1: 4px;  --sp-2: 8px;  --sp-3: 16px;
    --sp-4: 24px; --sp-5: 32px; --sp-6: 48px; --sp-7: 64px;
    --font-sans: system-ui, -apple-system, "Noto Sans CJK TC", sans-serif;
    --font-mono: ui-monospace, "JetBrains Mono", monospace;
  }
  * { box-sizing: border-box; }
  body {
    font-family: var(--font-sans);
    color: var(--color-text);
    background: var(--color-bg);
    max-width: 900px;
    margin: var(--sp-5) auto;
    padding: 0 var(--sp-3);
    line-height: 1.55;
  }
  h1 { font-size: 32px; line-height: 1.2; font-weight: 500;
       border-bottom: 2px solid var(--color-text); padding-bottom: var(--sp-2); }
  h2 { font-size: 24px; line-height: 1.3; font-weight: 500;
       margin-top: var(--sp-5); color: var(--color-text); }
  h3 { font-size: 20px; line-height: 1.4; font-weight: 500; }
  code { font-family: var(--font-mono); font-size: 14px;
         background: var(--color-border); padding: 2px 4px; border-radius: 3px; }
  pre { background: var(--color-border); padding: var(--sp-3); border-radius: 4px;
        overflow-x: auto; }
  pre code { background: none; padding: 0; }
  table { border-collapse: collapse; width: 100%; margin: var(--sp-3) 0; }
  th, td { border: 1px solid var(--color-border); padding: var(--sp-2); text-align: left; }
  th { background: var(--color-border); font-weight: 500; }
  details { border: 1px solid var(--color-border); border-radius: 4px;
            padding: var(--sp-2); margin: var(--sp-2) 0; }
  summary { cursor: pointer; font-weight: 500; }
  details[open] { padding: var(--sp-3); }
  .badge { display: inline-block; padding: 2px 8px; border-radius: 4px;
           font-size: 14px; font-weight: 500; }
  .badge-pass { background: #D1FAE5; color: #065F46; }
  .badge-partial { background: #FEF3C7; color: #92400E; }
  .badge-fail { background: #FEE2E2; color: #991B1B; }
  .figure { margin: var(--sp-4) 0; padding: var(--sp-3);
            border: 1px solid var(--color-border); border-radius: 6px; }
  .figure img { display: block; max-width: 100%; height: auto; margin: 0 auto; }
  .caption { margin-top: var(--sp-2); font-size: 14px; color: var(--color-muted); }
  @media print {
    body { max-width: none; background: white; }
    .figure { break-inside: avoid; }
    details { border: none; }
    details > summary { display: none; }
    details > * { display: block !important; }
  }
</style>
</head>
<body>
  <!-- Content here -->
</body>
</html>
```

特性：
- ~3 KB CSS（vs Tailwind 25 KB）
- 0 build step、0 dependency、零網路呼叫
- 內含 print mode（`@media print` 自動展開所有 details）
- 已遵守 Anthropic 設計禁忌（no gradient / no glass morphism / no multi-indigo / no emoji headers）
- 響應式（max-width 900px + viewport meta）
- 中文友善（Noto Sans CJK TC fallback）

### 5.2 升級觸發條件

| 偵測到 | 升級到 | 加什麼 |
|-------|-------|-------|
| `<!-- interactive: tab|filter|slider -->` 標記 | L3 | Alpine.js 15 KB inline |
| 報告含「設計系統」「元件庫」「prototype」關鍵字 | Tailwind | `--tailwind` flag |
| `--polished` flag (D19) | frontend-design plugin | invoke `frontend-design:frontend-design` skill |

### 5.3 主題資料夾條件（D17 修正 D4）

```python
def should_use_topic_folder(md_path: Path) -> bool:
    content = md_path.read_text()
    line_count = content.count('\n')
    figure_needed_count = content.count('<!-- figure-needed:')
    img_count = content.count('![](')  # markdown images
    total_figures = figure_needed_count + img_count
    return line_count >= 200 or total_figures >= 5
```

| 案例 | 結果 |
|------|------|
| 短報告 (50 行 + 1 圖) | 單檔 `foo.preview.html` 旁邊 |
| 中報告 (180 行 + 3 圖) | 單檔 `foo.preview.html` |
| 大報告 (250 行 + 4 圖) | **主題資料夾** `foo/` |
| feature-layered-observation 10 章節 + 6 圖 | **主題資料夾** + 章節分檔 |

---

## 6. Evaluator-Optimizer 模式規範（D18）

### 6.1 模式定義（per [Anthropic agent playbook](https://www.anthropic.com/research/building-effective-agents)）

```
Generator LLM → Output → Evaluator LLM → Verdict + Suggestions
                              ↓
                          {pass}? → Done
                          {fail}? → Patch input → max 1-2 retry → Final escalate
```

### 6.2 本專案實例（已實作）

| Generator | Evaluator | Optimizer 規則 |
|-----------|-----------|--------------|
| `image-gen` (codex `$imagegen` + cairo) | `image-vision-check` (Claude Read 6 維 checklist) | D8: 自動重生 1 次 → 失敗 escalate PI |
| `pptx-build` slide draft | `pptx-build/prompts/visual_review.md` (Claude Vision 10-checkpoint) | 重 render 1 次 → 失敗 escalate |
| (Phase 2) `html-preview` 渲染 | `html-preview/prompts/design_taste_check.md` (Anthropic 設計禁忌自審) | 修 1 次 → escalate |

### 6.3 新 evaluator-optimizer skill 模板

```markdown
---
name: <generating-X> + <evaluating-X>     # 雙 skill 對
description: |
  <generator>: Generate X. Use when ...
  <evaluator>: Score X against N-dim checklist; output quality.json. Use when X produced.
---

## Workflow
1. <generator> produces output → save to staging
2. <evaluator> reads output + checklist → JSON verdict (pass/partial/fail) + suggestions
3. If pass: commit
4. If partial/fail: max 1 retry with patched input
5. If still fail: escalate to PI with checklist
```

---

## 7. 既有 22 skill audit 計畫（D20 落實）

### 7.1 Audit 三項目

對每個 skill 檢核：

| Audit | Check | 修正動作 |
|-------|-------|---------|
| **A. Description directive** | 是否 directive 動詞開頭？是否含 USE WHEN / SKIP WHEN？是否第三人稱？ | 改寫 description |
| **B. Body 規範** | <500 lines？refs 1-level-deep？>100 行 子檔 TOC？無 Windows path？無 voodoo constants？ | 拆檔 / 加 TOC / 移除 |
| **C. Evals 補件** | tests/evals.json 是否存在？是否 ≥3 cases？ | 寫 3 evals |

### 7.2 22 skill 優先順序

| Tier | Skill | 優先 | 補做時點 |
|------|-------|-----|---------|
| **Tier A 高使用率** | weekly-report / structured-tech-report / results-report / feature-layered-observation / pptx-build / myPPT | P0 | Phase 3 同步 |
| **Tier B 中使用率** | research-loop / cycle-init / check-staleness / run-evaluator / multi-sample-consistency / methodology-audit / conclude-research | P1 | Phase 3 後 |
| **Tier C 低使用率/工具型** | data-audit / observation-analysis / inject-hypothesis / problem-framing-ideation / pivot-direction / review-evidence / research-dashboard / cycle-state / known-pitfalls / report / citation-verification / auc-confound-guard / memory-consolidation / provenance-tier-audit | P2 | 視需要 |

### 7.3 Evals JSON 範本

```json
{
  "skill": "weekly-report",
  "version": "1.0",
  "evals": [
    {
      "id": "eval_01_typical_trigger",
      "query": "整理本週 self-phasing 進度給教授看",
      "expected_behavior": [
        "Skill activates (description 觸發)",
        "Output produces W1-W7 master_draft.md",
        "audience 識別為 PI"
      ]
    },
    {
      "id": "eval_02_skip_scenario",
      "query": "幫我寫 git commit message",
      "expected_behavior": [
        "Skill does NOT activate (應由 git-commit 觸發)"
      ]
    },
    {
      "id": "eval_03_edge_case",
      "query": "整理 4-5 月份 multi-week 進度",
      "expected_behavior": [
        "Skill activates",
        "Asks PI for week range clarification (>7 days non-trivial)",
        "Does not silently span multi-week without ack"
      ]
    }
  ]
}
```

---

## 8. Anti-pattern Lint Checklist（D20 + 通用）

對所有 skill / SOP / 文件，定期 (每月) audit：

```bash
# Anti-pattern lint script (to be created at scripts/lint/skill_lint.sh)

# A. Windows-style paths
grep -rn '\\\\' .claude/skills/ && echo "❌ Windows paths found"

# B. Time-sensitive info (without 'Old patterns' wrap)
grep -rn -E "(after|before) \d{4}" .claude/skills/ | grep -v "Old patterns"

# C. SKILL.md > 500 lines
for f in .claude/skills/*/SKILL.md; do
  lines=$(wc -l < "$f")
  [ "$lines" -gt 500 ] && echo "❌ $f: $lines lines (>500)"
done

# D. Description first-person violation
grep -rn -E "^description:.*( I | I'd | I'll | I am )" .claude/skills/

# E. Voodoo constants in scripts (numbers without comment)
grep -rn -E "(TIMEOUT|DELAY|RETRIES) = [0-9]+\s*$" .claude/skills/ | head -10
```

---

## 9. SOP 對照本專案 22 個現有 skill 的修正建議

（簡略表，完整 audit 報告於 Phase 3 個別 skill 修改時補）

| Skill | description 是否 directive？ | gerund 命名？ | evals 已有？ | 優先補件 |
|-------|--------------------------|------------|-----------|---------|
| `image-gen` | ✓ (Phase 1 已修) | ✗ (動詞) | △ (demo_topic 算非正式) | 加正式 evals.json |
| `image-vision-check` | ✓ (Phase 1 已修) | ✗ | △ | 加正式 evals.json |
| `weekly-report` | △ (改 directive) | ✗ | ✗ | description rewrite + evals |
| `structured-tech-report` | △ | ✗ | ✗ | description rewrite + evals |
| `results-report` | △ | ✗ | ✗ | description rewrite + evals |
| `feature-layered-observation` | △ | ✗ | ✗ | description rewrite + evals |
| `pptx-build` | △ | ✗ | ✗ | description rewrite + evals |
| `myPPT` | △ | ✗ (camelCase!) | ✗ | description rewrite + evals + 不改名 |
| (其他 14 skill) | TBD audit | (大多動詞型) | ✗ | Tier B/C 排程 |

---

## 10. 開發新 skill 的 SOP 檢核（pre-flight）

新 skill 開發前必須 satisfy：

```markdown
## New skill checklist

### Naming
- [ ] Gerund form: 使用 verb+ing 命名（processing-X / analyzing-Y）
- [ ] Lowercase + hyphens only
- [ ] Not 'helper' / 'utils' / 'tools' / 含 'anthropic' 或 'claude'

### Frontmatter
- [ ] description ≤1024 chars, third-person, directive form
- [ ] description 含 USE WHEN + SKIP WHEN 雙列表
- [ ] description 含 ≥5 個具體觸發詞

### Body
- [ ] <500 lines
- [ ] References 1-level deep only
- [ ] >100 lines 子檔含 TOC
- [ ] All paths use forward slashes
- [ ] 一致 terminology
- [ ] No time-sensitive info（或 wrapped in `<details>` Old patterns）

### Workflow（如複雜）
- [ ] Checklist for Claude to copy
- [ ] Validation loop（如有 quality gate）
- [ ] Max retry rules（per evaluator-optimizer pattern）

### Scripts（如有）
- [ ] No Windows paths
- [ ] Solve, don't punt
- [ ] Documented constants (no voodoo numbers)
- [ ] Required packages listed in SKILL.md
- [ ] Pre-flight check 與 dependency declaration

### Evals
- [ ] tests/evals.json with ≥3 cases
- [ ] eval_01: 典型觸發 (must activate)
- [ ] eval_02: SKIP scenario (must NOT activate)
- [ ] eval_03: edge case (correct decision)

### Cross-references
- [ ] SKILL.md 含 "See Also" 指向 spec doc / related skills
- [ ] design tokens 用 CSS variables（如產 HTML）
- [ ] 不重新發明輪子（檢查 frontend-design / cookbook 範例是否已涵蓋）
```

---

## 11. 開發 workflow（per Anthropic best-practices）

### 11.1 Evals-driven development（取代 plan-driven）

```
1. 識別 gaps：跑 Claude 在代表性任務上 → 紀錄 failure
2. 寫 3 evals 涵蓋主要 fail 場景
3. 跑 baseline (without skill) → 紀錄分數
4. 寫最小 SKILL.md 通過 evals
5. 迭代直到 evals 全綠
```

### 11.2 Two-Claude iteration（取代 single-author）

```
Claude A (designer) ←→ Claude B (user)
         ↑                  ↓
         └── observe gaps ──┘

1. Claude A 撰寫初版 SKILL.md
2. Claude B 用 fresh session 跑 evals
3. 觀察 Claude B 行為 → return feedback to Claude A
4. Claude A refine → re-test with B
5. Repeat
```

對應 superpowers:subagent-driven-development workflow，可直接用既有 framework。

---

## 12. 後續落地任務

| Task | Owner | 時點 |
|------|-------|------|
| 修 design doc 加 D15-D20 + 點到本 SOP | 本回合 | 已寫完，下一個 commit |
| Phase 2 html-preview skill 改 MVP.css 預設（D15） | Phase 2 開發者 | Phase 2 plan 撰寫前 |
| Phase 3 開始時，6 Tier A skill 同步補 evals + description audit（D20 P0） | Phase 3 開發者 | Phase 3 plan |
| 寫 anti-pattern lint script `scripts/lint/skill_lint.sh` | TBD | Phase 3 ship 後 |
| Tier B/C 16 skill 補 evals | TBD | Phase 4（暫定） |
| frontend-design opt-in 整合（D19） | Phase 2 期間 | Phase 2 |

---

## Appendix A — 關鍵外部參考清單

1. **[thariqs/html-effectiveness](https://thariqs.github.io/html-effectiveness/)** — 9 category + 20 範例 HTML 檔案。
2. **[Anthropic Skills best-practices](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/best-practices)** — 完整 SKILL.md 規範。
3. **[Anthropic Building Effective Agents](https://www.anthropic.com/research/building-effective-agents)** — 6 agent design pattern 含 evaluator-optimizer。
4. **[Anthropic Cookbook agents/](https://github.com/anthropics/anthropic-cookbook/tree/main/patterns/agents)** — 實作範例。
5. **[MVP.css](https://andybrewer.github.io/mvp/)** — 我們選用的 minimal CSS framework。
6. **[Simple.css](https://simplecss.org/)** — alternative，類似 MVP.css。
7. **[Pico CSS](https://picocss.com/)** — alternative，較完整 minimal framework。
8. **[Simon Willison: Unreasonable Effectiveness of HTML](https://simonwillison.net/2026/May/8/unreasonable-effectiveness-of-html/)** — 業界 May 2026 觀察。

## Appendix B — 22 既有 skill 完整清單（FYI，audit 對象）

```
auc-confound-guard, check-staleness, citation-verification, conclude-research,
confirmation-protocol, cpp-change, cycle-init, cycle-state, data-audit,
doc-standards, fast-learning-coach, feature-layered-observation, grill-me,
init-research, inject-hypothesis, known-pitfalls, memory-consolidation,
methodology-audit, multi-sample-consistency, myPPT, observation-analysis,
pivot-direction, pptx-build, problem-framing-ideation, provenance-tier-audit,
research-dashboard, research-loop, results-analysis, results-report,
review-evidence, run-evaluator, structured-tech-report, validation-protocol,
weekly-report
```

實際 25 skill。Phase 1 新增 image-gen + image-vision-check = 27 total。

## Appendix C — 文件版本與 supersedence

- v1.0 (2026-05-11): 初版，整合研究發現後寫成。
- 取代 design doc 的 D3 / D4 / D12（仍保留 design doc 為歷史紀錄）。
- 新增 D15–D20 為 design doc 第 15-20 號決策。

---

**End of SOP v1.0.**
