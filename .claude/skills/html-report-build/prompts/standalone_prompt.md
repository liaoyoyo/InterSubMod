# standalone_prompt.md — LLM-direct standalone HTML report (PI-end view)

> **Read this prompt fully before generating HTML.** A standalone report is NOT a markdown companion; it's the **terminal interface language** for a one-shot delivery (PI / advisor / collaborator / publication-supplement). The reader gets ONE `.html` file and that's the end of the chain.

---

## Your role

Given a finalized `.md` source (validated, not in-progress), produce **ONE self-contained `.html`** alongside the source. The HTML must give the reader an experience that `.md` rendered in GitHub **cannot** match — sticky navigation, collapsible cards, language-of-color semantics, inline SVG diagrams, print-ready.

**This is the third mode** of html-report-build:
- `report` mode = prose companion (Tailwind-prose style, light restructure)
- `slide` mode = 16:9 PPT deck
- **`standalone` mode (this)** = PI-end terminal HTML, max interactive features

**Output filename**: `{md_basename}.standalone.html` (sibling of source `.md`).

---

## When to use standalone vs report mode

| Signal | Mode |
|--------|------|
| Source in `docs/reports/validated/` AND status = "validated" | **standalone** |
| Source in `docs/experiments/concluded/` AND tier ⭐4-5 | **standalone** |
| Source is PI report / errata / verdict document | **standalone** |
| Source has `audience: PI / advisor / collaborator` frontmatter | **standalone** |
| Source is in-progress / draft / will re-diff | report (companion) |
| Source <100 lines OR no `<h2>` structure | report (single-file) |
| User explicitly says "終版" / "給 PI" / "single file" / "standalone" | **standalone** |

If ambiguous → ask user (don't improvise).

---

## 10 Rules

### Rule 1 — Layout: sticky left TOC + main content (only HTML can do this)
Mandatory page structure:

```html
<div class="layout">
  <aside class="toc">
    <h3>章節導航</h3>
    <ol>
      <li><a href="#tldr">TL;DR</a></li>
      <li><a href="#section-1">Section 1</a></li>
      ...
    </ol>
  </aside>
  <main class="report">
    <article>...</article>
  </main>
</div>
```

CSS:
```css
.layout { max-width: 1320px; margin: 0 auto; display: grid;
          grid-template-columns: 220px 1fr; gap: var(--sp-5); padding: var(--sp-5); }
aside.toc { position: sticky; top: var(--sp-4); align-self: start;
            max-height: calc(100vh - var(--sp-7)); overflow-y: auto; }
@media (max-width: 900px) { .layout { grid-template-columns: 1fr; }
                            aside.toc { position: static; } }
```

TOC must include `id`-anchored links to ALL `<h2>` sections + key sub-anchors (errata IDs, callout IDs).
Mark core/critical entries with class `core` (visual weight + ★ suffix).

### Rule 2 — Article header with metadata panel (not just `<h1>`)
Top of report: structured metadata, not loose paragraphs.

```html
<header class="report-header">
  <div class="pretitle">{report-class · category}</div>
  <h1>{title}</h1>
  <dl class="meta">
    <dt>Build date</dt>      <dd>{ISO date}</dd>
    <dt>Audience</dt>        <dd>{PI / lab / future-self}</dd>
    <dt>Parent report</dt>   <dd><code>{path}</code></dd>
    <dt>Status</dt>          <dd>{validated · last verified ...}</dd>
  </dl>
  <div class="verdict-banner">{one-sentence verdict — uses ONE gradient}</div>
</header>
```

### Rule 3 — Section cards for repeated logical units
When source `.md` has **5+ similar units** (errata E1-E5, findings F1-Fn, samples 1-7), use `<details>` cards:

```html
<details class="errata-card core" id="e4" open>
  <summary>
    <span class="errata-id">E4</span>
    <span class="errata-section">§6.4 / §6.5</span>
    <span class="errata-title">{one-line summary}</span>
    <span class="errata-badge core">CORE</span>
  </summary>
  <div class="errata-body">{full content}</div>
</details>
```

Rules:
- **Core / critical** items: `open` by default, distinct border color
- **Minor / supporting**: collapsed by default
- Each card has unique `id` matching TOC anchor
- Summary fits in one line (≤80 chars)
- Body can contain tables, callouts, SVG, code

### Rule 4 — Inline SVG for small diagrams (NEW capability vs report mode)

**When to use SVG (inline, no external file)**:

| Diagram type | Use SVG? | Why |
|--------------|----------|-----|
| Flow / causal chain (≤10 nodes) | ✅ YES | LLM can write SVG directly; precise, scalable |
| Proportion bar / horizontal stacked | ✅ YES | Smaller than PNG, semantic, accessible |
| Architecture diagram (boxes + arrows, ≤8 boxes) | ✅ YES | No external tool needed |
| Icon (warning / check / info, instead of emoji) | ✅ YES | Avoids emoji-in-heading taboo |
| Sparkline / mini trend (≤20 points) | ✅ YES | Inline data viz |
| Real data plot (continuous, >50 points) | ❌ NO | Use matplotlib → PNG → base64 inline |
| IGV screenshot / photograph | ❌ NO | PNG only |
| Interactive chart (zoom, hover tooltip) | ❌ NO | Out of LLM scope; use D3.js or skip |

**SVG requirements**:
- Use `<svg viewBox="..." role="img" aria-labelledby="title-id desc-id">` with `<title>` + `<desc>` for accessibility
- Colors via design tokens: `<rect fill="var(--c-accent)">` (not hardcoded hex)
- `font-family: var(--ff-body)` for `<text>` elements
- No animations (taboo: flash/glow effects)
- No SVG filters (`<filter>`, `<feGaussianBlur>`) — too heavy + dated aesthetic
- Max viewBox 800×400; if larger, use PNG instead
- Embed at most 5 SVG diagrams per report (keep file <100 KB)

**SVG starter examples** — see `components/svg_*.html`:
- `svg_flow_diagram.html` — boxes + arrows causal chain
- `svg_compare_bar.html` — horizontal proportion bar
- `svg_icon_set.html` — warning / check / info icons

### Rule 5 — Language-of-color row semantics
Tables use class-based row coloring (NOT inline `style`):

| Class | Meaning | Background |
|-------|---------|------------|
| `tr.row-bug` | Failure / negative / blocker | red-tint #FEF2F0 |
| `tr.row-safe` | Pass / verified / positive | green-tint #DCFCE7 |
| `tr.row-warn` | Partial / improving / caveat | yellow-tint #FEF3C7 |

Apply to comparison tables, version diffs, validation results. Reader scans the column of colors and gets the verdict in 2 seconds.

### Rule 6 — Stat grid for focal numbers
Below TL;DR / verdict, place 3-6 key numbers in a grid:

```html
<div class="stat-grid">
  <div class="stat-box">
    <div class="number bug">17.3 : 1</div>
    <div class="label">HP1:HP2 全基因組偏移 (baseline)</div>
  </div>
  ...
</div>
```

Variants: `.number` (default accent), `.number.bug` (red), `.number.success` (green), `.number.warning` (yellow).

This is the **5-second skim layer** — PI reads only these 4 numbers + verdict banner if short on time.

### Rule 7 — Print mode is first-class (PI may print for advisor)
Mandatory `@media print` rules:

```css
@media print {
  body { background: white; }
  .layout { display: block; max-width: 100%; padding: 0; }
  aside.toc { display: none; }                              /* sticky TOC hidden */
  main.report { box-shadow: none; border: none; padding: 0; }
  .errata-card { page-break-inside: avoid; border-color: #999; }
  details[open] > summary, details > summary { list-style: none; }
  details:not([open]) > * { display: block !important; }   /* force expand */
  details:not([open]) > summary { background: white !important; }
  a[href^="http"]::after { content: " (" attr(href) ")"; font-size: 0.85em; color: #666; }
  a[href^="../"]::after { content: " [" attr(href) "]"; font-size: 0.8em; color: #999; }
  pre, table, blockquote, .errata-card, svg { page-break-inside: avoid; }
  h1, h2, h3 { page-break-after: avoid; }
}
```

Print test: when user prints, no content should be cut, no JS interactivity needed, all `<details>` forced open, URLs printed after links.

### Rule 8 — Inline everything, ONE file, ZERO external dep
- All CSS inline in single `<style>` block (no `<link rel="stylesheet">`)
- All SVG inline (no `<img src="*.svg">`)
- Real data PNG images base64-inline ONLY IF total HTML stays <500 KB
- If file would exceed 500 KB: convert to topic folder OR use external `<img src="../figures/...">`
- No `<script src="...cdn...">` — standalone reports do NOT need JS interactivity (TOC is `<a href="#">` anchors, details are HTML-native)
- Reader can save the file to disk, email it, open offline — must work fully.

### Rule 9 — 6-Taboo audit (canonical list in SKILL.md)
Standalone mode allows **TWO** `linear-gradient` exceptions (one more than report mode):
1. `.verdict-banner` (header verdict, ONE gradient — green=POSITIVE / blue=NEUTRAL / yellow=CAUTION)
2. `.conclusion-arrow` (final section verdict, optional, if separate from header banner)

All other taboos strict: no `backdrop-filter`, no `text-shadow`, no glow `box-shadow`, no multi-indigo, no emoji-in-headings.

Audit command:
```bash
grep -cE '(linear-gradient|radial-gradient|backdrop-filter|text-shadow|0 0 [0-9]+px)' output.standalone.html
# Expected: ≤2 (gradients for verdict-banner + optional conclusion-arrow)
```

### Rule 10 — Companion rule (D11) — NEVER replace .md
Source `.md` is **always preserved**. Output filename suffix `.standalone.html` distinguishes from:
- `{basename}.html` (report mode, prose companion)
- `{basename}/` (topic folder mode)
- `{basename}.standalone.html` (this mode)

All three modes can coexist for the same `.md` if user wants both PI-view and grep-view.

INDEX.md continues to reference the `.md`. The standalone HTML is a delivery artifact, not the SoT.

---

## Standard skeleton structure

Reference `templates/standalone_skeleton.html` for the canonical layout. Key macros to fill:

| Placeholder | Source |
|-------------|--------|
| `{{REPORT_TITLE}}` | `.md` `<h1>` or frontmatter `title:` |
| `{{REPORT_CLASS}}` | frontmatter `report_class:` (e.g. "Errata Companion · Self-Phasing") |
| `{{BUILD_DATE}}` | frontmatter `build_date:` |
| `{{AUDIENCE}}` | frontmatter `audience:` |
| `{{PARENT_REPORT}}` | frontmatter `parent_report:` |
| `{{STATUS}}` | frontmatter `status: validated · last_verified: ...` |
| `{{VERDICT_LINE}}` | frontmatter `verdict:` (one sentence, fits in banner) |
| `{{TOC_ENTRIES}}` | All `<h2>` headings + key card IDs |
| `{{TLDR_BLOCK}}` | First section content + stat-grid |
| `{{CARD_GRID}}` | Section cards if source has repeated units |
| `{{TIMELINE_TABLE}}` | If source has "修訂歷程" / "Timeline" section |
| `{{REFERENCES_LIST}}` | If source has references section |

---

## Design Principles (Read first)

Before generating, **Read `references/design_principles.md`** — 12 rules from Tufte / Garr Reynolds / Nancy Duarte / Robin Williams CRAP / Assertion-Evidence / Nature figure guidelines / WCAG. Self-audit against the **6-item pre-publish checklist** before Write.

Key rules to apply for standalone mode:
- **Rule 1 (5-second test)**: `.stat-grid` + `.verdict-banner` must surface takeaway in 5 sec
- **Rule 3 (Assertion-Evidence)**: each `<h2>` heading should be a complete claim, not a topic label
- **Rule 5 (CRAP)**: Contrast (font size hierarchy) + Repetition (consistent badge/card classes) + Alignment (grid layout) + Proximity (card summary+body tight, cards separated)
- **Rule 7 (Hierarchy)**: verdict-banner > stat-grid > core cards (open) > minor cards (collapsed) > footer refs
- **Rule 9 (Colorblind-safe)**: row-bug/safe/warn classes must ALSO have text verdict ✓/✗/⚠ — never color-only
- **Rule 11 (Vector first)**: SVG > PNG for diagrams ≤10 nodes

## Workflow (your steps)

1. **Read** source `.md` fully (use Read tool).
2. **Read** `references/design_principles.md` (Rule 1-12 internalized).
3. **Parse frontmatter** (lines between `<!--` and `-->` or `---` and `---`).
3. **Decide layout**:
   - 5+ similar units → use section cards
   - 3-6 key numbers → use stat-grid
   - Has flow/causal chain → consider inline SVG (≤5 diagrams)
4. **Read** `templates/standalone_skeleton.html` for canonical structure.
5. **Read** `components/svg_*.html` if planning SVG diagrams (few-shot).
6. **Generate** HTML: fill skeleton, restructure body per Rules 1-7.
7. **Self-audit** (Rule 9 taboo grep + structure checks).
8. **Write** output to `{md_basename}.standalone.html` (sibling of source).
9. **Report** path + file size to user.

---

## What NOT to do

- ❌ Do not call `python3` / `pandoc` / `markdown` lib — LLM-direct only.
- ❌ Do not use external CDN (no `<script src="https://...">`).
- ❌ Do not use D3.js / Chart.js / any JS chart lib — write inline SVG or skip.
- ❌ Do not split into multiple files — standalone = ONE file.
- ❌ Do not put real data scatter plots inline (>50 points) — use matplotlib → PNG.
- ❌ Do not rewrite the source `.md`.
- ❌ Do not add tracking / analytics / external font CDN.
- ❌ Do not animate SVG (no `<animate>`, no `transform: scale` keyframes).
- ❌ Do not use SVG `<filter>` for visual effects (drop-shadow / blur) — dated aesthetic.
- ❌ Do not exceed 500 KB total file size (compress images, defer to external if needed).
