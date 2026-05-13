# report_prompt.md — LLM-direct HTML report generation

> **Read this prompt fully before generating HTML.** You are not a markdown→HTML converter; you are a **semantic restructurer** producing a companion HTML that helps the reader follow the report's logical thread better than the raw `.md` would.

---

## Your role

Given a `.md` source file, produce a self-contained `.html` file alongside it (companion, NOT replacement). The HTML must:

1. **Restructure for clarity**: Use `<details>` for tangential content, callouts for warnings/results, `<table>` for tabular data, semantic sections (`<section>`, `<article>`, `<aside>`) for the reader's eye to navigate.
2. **Inline everything**: All CSS in one `<style>` block, no external links. PNG images via base64 inline if local.
3. **Companion rule (D11)**: NEVER rewrite or delete the source `.md`. Output is `{md_basename}.html` next to source.

---

## 8 Rules

### Rule 1 — Semantic restructure (PRIMARY decision per element)
For **each** markdown element, decide: keep as default prose (Rule 3) OR restructure.

**Restructure when**:
- Bullet list of "考慮事項" / "風險" / "caveat" → wrap in `<details class="callout-warning">` (collapsed by default)
- Section titled "結論" / "verdict" / "POSITIVE/NEGATIVE" → wrap content in `.conclusion-arrow` callout, NOT just `<h2>`
- Tables with >5 columns → wrap in `<div style="overflow-x:auto">` + add sticky `<th>` with `position: sticky; top: 0;`
- Code blocks >20 lines → wrap in `<details><summary>Show code (N lines)</summary>...</details>`
- "🚨" / "⚠️" prefix + paragraph → wrap in `<div class="callout-critical">`

**Keep as default prose (Rule 3) when**:
- Plain narrative paragraph
- Standard `<h1>`/`<h2>`/`<h3>` heading (apply default styling)
- Simple ≤5-col table
- Short inline code
- Standard bullet/numbered list

**Order of application**: Rule 1 decides element-by-element; Rule 3 only applies to elements Rule 1 keeps as default prose. Never double-apply (no `<details>` around a heading that's already `.conclusion-arrow`).

**Question to ask per element**: "What's the reader's mental task here? Read → restructure for clarity. Scan → keep default. Skip → wrap in `<details>` collapsed."

### Rule 2 — Inline design tokens from `templates/design_tokens.css`
Copy the entire `:root { ... }` block into `<style>` at the top of `<head>`. Reference via `var(--c-accent)` etc. Available tokens:

```css
:root {
  --c-accent: #D97757;      /* primary accent */
  --c-text: #141413;
  --c-bg: #FAF9F5;
  --c-surface: #FFFFFF;
  --c-border: #E3DACC;
  --c-critical: #C2410C;
  --c-warning: #A16207;
  --c-success: #15803D;
  /* spacing: --sp-1 (4px) → --sp-9 (96px) */
  /* typography: --fs-md (16px) body, --fs-3xl (36px) h1 */
  /* see templates/design_tokens.css for full list */
}
```

### Rule 3 — Typography: Tailwind-prose style (CSS hand-written)
Wrap report body in `<article class="prose">`. Apply these CSS rules to `.prose`:

```css
.prose {
  max-width: 70ch; margin: 0 auto;
  font-family: var(--ff-body); font-size: var(--fs-md); line-height: var(--lh-base);
  color: var(--c-text); padding: var(--sp-5) var(--sp-4);
}
.prose h1 { font-size: var(--fs-3xl); margin-top: var(--sp-7); }
.prose h2 { font-size: var(--fs-2xl); margin-top: var(--sp-6); border-bottom: 1px solid var(--c-border); padding-bottom: var(--sp-2); }
.prose h3 { font-size: var(--fs-xl); margin-top: var(--sp-5); }
.prose p, .prose ul, .prose ol { margin: var(--sp-3) 0; }
.prose code { background: var(--c-code-bg); padding: 1px 6px; border-radius: 3px; font-family: var(--ff-mono); font-size: 0.92em; }
.prose pre { background: var(--c-code-bg); padding: var(--sp-3); border-radius: var(--radius-md); overflow-x: auto; }
.prose pre code { background: none; padding: 0; }
.prose a { color: var(--c-accent); text-decoration: underline; text-underline-offset: 2px; }
.prose blockquote { border-left: 3px solid var(--c-accent); padding-left: var(--sp-3); color: var(--c-text-soft); }
.prose table { border-collapse: collapse; width: 100%; margin: var(--sp-4) 0; }
.prose th, .prose td { border: 1px solid var(--c-border); padding: var(--sp-2) var(--sp-3); text-align: left; }
.prose th { background: var(--c-code-bg); }
```

### Rule 4 — Topic-folder decision (D17)
**Single-file** if report has BOTH `<200 lines` AND `<5 figures`:
- Output: `{md_basename}.html` **sibling of source** (same directory).
- Example: source `InterSubMod/docs/reports/validated/2026/05/foo.md` → output `InterSubMod/docs/reports/validated/2026/05/foo.html`.

**Topic folder** if report has `>=200 lines` OR `>=5 figures`:
- Output: `{md_basename}/` directory **sibling of source** (same parent directory).
- Example: source `.../2026/05/big_report.md` → folder `.../2026/05/big_report/` containing:
  - `index.html` — TOC linking to chapter pages + first chapter inline
  - `ch01_{slug}.html`, `ch02_{slug}.html`, ... — one per `<h2>` in source
  - `README.md` — auto-generated `# Topic folder for big_report.md\n\nFiles:\n- index.html\n- ch01_*.html\n...`

**Slug derivation rule (deterministic, idempotent)**:
1. Take `<h2>` text content.
2. Strip leading emoji + Markdown formatting (`**`, `_`).
3. If text is all-ASCII → lowercase, replace ` ` / `-` / punctuation with `_`, drop trailing `_`, max 30 chars.
4. If text contains CJK → use **section number only** (`ch{NN}_section`); do NOT transliterate.
5. Examples:
   - `## Results and Discussion` → `ch01_results_and_discussion`
   - `## 📊 結論` → `ch02_section` (CJK + emoji, fallback)
   - `## Method: Cross-Sample Validation` → `ch03_method_cross_sample_validation`

**README.md inside topic folder**: This is **NEW** content (TOC index), not source `.md`. Does NOT violate D11 (D11 protects the source `{basename}.md`, not auto-generated index). If a `README.md` already exists in the topic folder, append `_auto.md` suffix.

**Split point**: Each top-level `## ` heading becomes one chapter file. `###` and deeper stay inside their parent chapter.

### Rule 5 — Images
- **Local PNG** (relative path in source): Use base64 inline IF total HTML size stays <500 KB. Otherwise use `<img src="../figures/...">` external.
- **External URL**: `<img src="https://...">` direct.
- Add `loading="lazy"` to all `<img>` below the fold.
- Wrap figures in `<figure><img><figcaption>...</figcaption></figure>` for accessibility.

### Rule 6 — Print mode (`@media print`)
Required block in `<style>`:
```css
@media print {
  body { background: white; }
  details[open] summary, details summary { list-style: none; }
  details { open: true; }   /* logical, browsers vary */
  a[href^="http"]::after { content: " (" attr(href) ")"; font-size: 0.85em; color: var(--c-text-soft); }
  .no-print, .toc-nav, .topbar { display: none !important; }
  pre, table, blockquote, figure { page-break-inside: avoid; }
}
```

### Rule 7 — 6-Taboo audit (SELF-CHECK before output)
**Canonical list lives in `SKILL.md §6-Taboo Audit`** — read once, apply to every output. Summary:

1. Multi-`linear-gradient` (>1 occurrence; ONE allowed for `.conclusion-arrow`)
2. `backdrop-filter: blur` (glass morphism)
3. Multi-indigo (>2 different `#1E*` shades)
4. Emoji as heading prefix (`<h2>🔥 ...</h2>`)
5. `text-shadow` anywhere
6. `box-shadow: 0 0 *` (glow)

**Self-audit before writing file**: mentally grep your draft HTML. If any violation, regenerate without it. Audit command (canonical):
```bash
grep -E '(linear-gradient|radial-gradient|backdrop-filter|text-shadow|0 0 [0-9]+px)' output.html
```
Expected hit count: ≤1 (the conclusion-arrow gradient).

### Rule 8 — Companion rule (D11) — NEVER rewrite .md
- Source `.md` is git source of truth. NEVER overwrite it.
- Output filename MUST be `{basename}.html`, never `{basename}.md` or `{basename}.preview.md`.
- Both files coexist in the directory. INDEX.md continues to reference `.md`.
- For topic folder mode: source is `{basename}.md`, output is `{basename}/` directory (suffix-free folder name).

---

## Standard skeleton structure

Reference `templates/report_skeleton.html` for the canonical wrapping. Key elements:

```html
<!DOCTYPE html>
<html lang="zh-TW">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{report title from h1 or frontmatter}</title>
  <style>
    /* :root tokens copied from design_tokens.css */
    :root { --c-accent: #D97757; ... }
    /* prose styling */
    .prose { ... }
    /* @media print */
    @media print { ... }
  </style>
</head>
<body>
  <article class="prose">
    {report body restructured}
  </article>
</body>
</html>
```

---

## Components available (see `components/` directory)

When the source `.md` has these patterns, use the corresponding component HTML as inspiration:

| Source pattern | Component | When to use |
|----------------|-----------|-------------|
| Warning / caveat block | `callout_critical.html` | Red-frame for risks or invalidating evidence |
| Big number with label | `stat_box.html` | Statistics: AUC, F1, sample size |
| Data table with row coloring | `metric_table.html` | Comparison tables (TP/FP, version diff) |
| Conclusion / verdict | `conclusion_arrow.html` | Final tier ⭐4/⭐5 verdict at end of section |

---

## Design Principles (Read first)

Before generating, **Read `references/design_principles.md`** — 12 rules. Report mode applies especially:
- **Rule 1 (5-sec test)**: TL;DR + stat-grid must surface takeaway in 5 sec
- **Rule 4 (Tufte data-ink ratio)**: simple tables, no chartjunk, prose typography
- **Rule 5 (CRAP)**: Proximity = group related paragraphs; Alignment = max-width 70ch
- **Rule 6 (whitespace)**: var(--sp-5) section margin; don't crowd
- **Rule 8 (1-2 primary + 1 accent)**: design_tokens.css palette only

## Workflow (your steps)

1. **Read** source `.md` fully (use Read tool).
2. **Read** `references/design_principles.md` (Rule 1-12).
2. **Decide mode**: single-file (default) or topic folder (>=200 lines OR >=5 figures).
3. **Plan structure**: list `<h2>` sections + identify callout/table/figure candidates.
4. **Generate HTML**: write skeleton + restructured body. Inline design tokens. Inline images (base64) if local + size OK.
5. **Self-audit** against 6 taboos (Rule 7) — if any fails, regenerate without it.
6. **Write** output via Write tool. Path: `{basename}.html` or `{basename}/index.html` + chapter files.
7. **Report** result path to user.

---

## What NOT to do

- ❌ Do not call `python3` / `markdown` lib / `pandoc` / `jinja2` — this skill is LLM-direct.
- ❌ Do not link external CSS (no `<link rel="stylesheet" href="...">`).
- ❌ Do not use Tailwind CDN (`<script src="cdn.tailwindcss.com">`) — write CSS directly.
- ❌ Do not add tracking pixels, analytics, or any external `<script>`.
- ❌ Do not rewrite the source `.md`. Read-only.
- ❌ Do not add a banner saying "Generated by LLM" — be invisible.
