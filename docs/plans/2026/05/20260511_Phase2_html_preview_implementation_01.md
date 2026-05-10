---
title: "Phase 2 Implementation Plan — html-preview skill (MVP.css + conditional topic folder)"
date: 2026-05-11
status: ready_to_execute
version: 0.1.0
authors: [liaoyoyo2001 (PI), Claude Opus 4.7 (assistant)]
spec:
  - InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md (D1-D20)
  - InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md (SOP v1.0)
phase: 2 of 3
estimated_effort: 3-4 工作天 (15 tasks × 20-40 min/task)
tags: [implementation_plan, phase2, html_preview, mvp_css, evaluator_optimizer]
---

# Phase 2 Implementation Plan: html-preview skill (MVP.css + conditional topic folder)

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `html-preview` skill that converts markdown source `.md` → companion HTML viewer using **MVP.css minimal framework + Anthropic design tokens**, with **conditional topic folder** (only ≥200 lines OR ≥5 figures), evaluator-optimizer self-audit against design taboos, and 3 evals per Anthropic best-practices.

**Architecture:**
- Skill at `InterSubMod/.claude/skills/html-preview/`
- Core converter: `tools/dispatch.py` → reads source `.md`, decides single-file vs topic-folder, applies MVP.css template
- Conversion engine: Python `markdown` library (pure-python, no pandoc) + `jinja2` templates + `beautifulsoup4` post-processing
- Three rendering levels: **L1** (raw HTML, fallback), **L2** (MVP.css + design tokens, **default**), **L3** (+ Alpine.js for interactive markers)
- Evaluator-Optimizer pattern (D18): generated HTML → `design_taste_check.py` Claude self-audit → max 1 retry → escalate
- All design tokens **inline** as CSS variables (D16); zero CDN, zero network, zero build-step (D15)

**Tech Stack:** Python 3.10+ + `markdown` + `jinja2` + `beautifulsoup4` + `pyyaml` (already installed) + `Pillow` (already from Phase 1) · MVP.css 1.15 inline · Alpine.js 3.x inline (only for L3)

**Reference:**
- Design doc decisions D15-D20 — `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md`
- SOP §4-§6 (SKILL.md規範 + HTML 模板 + Evaluator-Optimizer 模式) — `InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md`
- Phase 1 plan (architectural patterns, TDD discipline) — `InterSubMod/docs/plans/2026/05/20260510_Phase1_image_gen_vision_check_implementation_01.md`

---

## Naming pragmatic exception (per N5 user decision)

User decided: "新 skill 強制 gerund". Per Anthropic best-practices, this skill ideally would be `previewing-html`. **However**, the design doc (3954 lines) and SOP refer to it as `html-preview` throughout. Renaming would invalidate dozens of cross-references.

**Resolution**: Keep skill name `html-preview` (noun-phrase form, still ACCEPTABLE per Anthropic's "Acceptable alternatives" list). Apply gerund rule strictly to Phase 3+ skills. Document this exception in skill body.

---

## Pre-Flight Check (run once before Task 1)

Phase 2 has dep installs the previous phase didn't have. Run these and ensure all `[OK]`:

```bash
# 1. Python markdown library
python3 -c "import markdown; print(f'[OK] markdown {markdown.__version__}')" 2>&1 \
  || pip install --user 'markdown>=3.5'

# 2. Jinja2 templating
python3 -c "import jinja2; print(f'[OK] jinja2 {jinja2.__version__}')" 2>&1 \
  || pip install --user 'jinja2>=3.0'

# 3. beautifulsoup4 (HTML manipulation)
python3 -c "import bs4; print(f'[OK] bs4 {bs4.__version__}')" 2>&1 \
  || pip install --user 'beautifulsoup4>=4.12'

# 4. From Phase 1: pyyaml + Pillow (should already be installed)
python3 -c "import yaml, PIL; print(f'[OK] phase1 deps')" 2>&1
```

Expected: 4 `[OK]` lines after install (if missing). Running the install commands automatically is acceptable — these are pure-python pip installs with no system dependencies.

---

## File Structure

```
InterSubMod/
├── .claude/skills/html-preview/                      # NEW skill
│   ├── SKILL.md                                       # ≤500 lines, directive description
│   ├── playbook.md
│   ├── templates/
│   │   ├── design_tokens.css                          # CSS variables (D16) — Anthropic palette
│   │   ├── base_l1.html.jinja                         # raw HTML (fallback, no framework)
│   │   ├── base_l2.html.jinja                         # MVP.css inline (DEFAULT, D15)
│   │   ├── base_l3.html.jinja                         # + Alpine.js (interactive)
│   │   └── topic_index.html.jinja                     # main index for topic folder
│   ├── components/
│   │   ├── readme_template.md.jinja                   # auto-gen folder README
│   │   └── design_taboos_audit.md                     # Anthropic禁忌 ref
│   ├── tools/
│   │   ├── preflight.sh                                # check pip deps
│   │   ├── md_to_html.py                              # markdown → HTML body
│   │   ├── topic_folder_decider.py                    # D17: should-split logic
│   │   ├── topic_folder_builder.py                    # build dir + README
│   │   ├── interactivity_detect.py                    # L2 vs L3 (D3 modified by D15)
│   │   ├── inline_assets.py                           # PNG → base64
│   │   ├── design_taste_check.py                      # eval-optimizer audit (D18)
│   │   └── dispatch.py                                # main entry
│   ├── prompts/
│   │   ├── design_taste_check_prompt.md               # for Claude self-audit
│   │   └── frontend_design_polish_prompt.md           # opt-in (D19)
│   ├── tests/
│   │   ├── test_md_to_html.py
│   │   ├── test_topic_folder_decider.py
│   │   ├── test_topic_folder_builder.py
│   │   ├── test_interactivity_detect.py
│   │   ├── test_inline_assets.py
│   │   ├── test_dispatch.py
│   │   └── fixtures/
│   │       ├── short_report.md                        # 50 lines, 1 fig (single-file case)
│   │       ├── medium_report.md                       # 180 lines, 3 figs
│   │       ├── large_report.md                        # 250 lines, 4+ figs (topic-folder case)
│   │       ├── interactive_marker.md                  # contains <!-- interactive: tab --> (L3 trigger)
│   │       └── sample_fig.png                         # for inline test
│   └── evals.json                                     # 3 evals per D20
│
└── examples/phase2_demo/                              # NEW demo
    └── (real PI report converted at Task 15)
```

Total new files: ~25 (templates + tools + tests + skill docs + demo).

---

### Task 1: Project scaffolding + dep install

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/{templates,components,tools,prompts,tests/fixtures}/.gitkeep`
- Create: `InterSubMod/examples/phase2_demo/.gitkeep`

- [ ] **Step 1: Create directories**

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
mkdir -p .claude/skills/html-preview/{templates,components,tools,prompts,tests/fixtures}
mkdir -p examples/phase2_demo
touch .claude/skills/html-preview/{templates,components,tools,prompts,tests/fixtures}/.gitkeep
touch examples/phase2_demo/.gitkeep
```

- [ ] **Step 2: Install Python deps**

```bash
pip install --user 'markdown>=3.5' 'jinja2>=3.0' 'beautifulsoup4>=4.12'
```

Verify each prints OK:
```bash
python3 -c "import markdown, jinja2, bs4; print('[OK] markdown', markdown.__version__); print('[OK] jinja2', jinja2.__version__); print('[OK] bs4', bs4.__version__)"
```

- [ ] **Step 3: Verify directories**

```bash
find .claude/skills/html-preview examples/phase2_demo -type d | sort
```
Expected: 7 directories.

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/html-preview/ examples/phase2_demo/
git commit -m "feat(html-preview): scaffolding for Phase 2

Phase 2 of HTML preview + image-gen 3-skill design.
Spec: InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md (D15-D20)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 2: preflight.sh dependency check

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/preflight.sh`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_preflight.sh`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_preflight.sh`:

```bash
#!/usr/bin/env bash
# Test that preflight.sh exits 0 when all deps present.
set -euo pipefail
SCRIPT="$(dirname "$0")/../tools/preflight.sh"

[ -x "$SCRIPT" ] || { echo "FAIL: $SCRIPT not executable"; exit 1; }

if "$SCRIPT" >/dev/null 2>&1; then
    echo "PASS: preflight.sh returns 0 in healthy env"
else
    echo "FAIL: preflight.sh returned non-zero"
    "$SCRIPT" 2>&1 | head -10
    exit 1
fi

echo "all preflight tests passed"
```

```bash
chmod +x .claude/skills/html-preview/tests/test_preflight.sh
```

- [ ] **Step 2: Run test (expect FAIL)**

```bash
bash .claude/skills/html-preview/tests/test_preflight.sh
```
Expected: `FAIL: ... preflight.sh not executable`.

- [ ] **Step 3: Write preflight.sh**

Create `InterSubMod/.claude/skills/html-preview/tools/preflight.sh`:

```bash
#!/usr/bin/env bash
# preflight.sh — verify html-preview skill dependencies.
set -euo pipefail
errors=0

# 1. Python 3.10+
PY_VER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
if python3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)"; then
    echo "[OK] python3 ${PY_VER}"
else
    echo "[FAIL] python3 ${PY_VER} (need >= 3.10)" >&2
    errors=$((errors+1))
fi

# 2. markdown library
if python3 -c "import markdown" 2>/dev/null; then
    echo "[OK] markdown $(python3 -c 'import markdown; print(markdown.__version__)')"
else
    echo "[FAIL] markdown missing. Install: pip install --user 'markdown>=3.5'" >&2
    errors=$((errors+1))
fi

# 3. jinja2 library
if python3 -c "import jinja2" 2>/dev/null; then
    echo "[OK] jinja2 $(python3 -c 'import jinja2; print(jinja2.__version__)')"
else
    echo "[FAIL] jinja2 missing. Install: pip install --user 'jinja2>=3.0'" >&2
    errors=$((errors+1))
fi

# 4. beautifulsoup4 library
if python3 -c "import bs4" 2>/dev/null; then
    echo "[OK] bs4 $(python3 -c 'import bs4; print(bs4.__version__)')"
else
    echo "[FAIL] bs4 missing. Install: pip install --user 'beautifulsoup4>=4.12'" >&2
    errors=$((errors+1))
fi

# 5. pyyaml (Phase 1 dep)
if python3 -c "import yaml" 2>/dev/null; then
    echo "[OK] pyyaml"
else
    echo "[FAIL] pyyaml missing. Install: pip install --user pyyaml" >&2
    errors=$((errors+1))
fi

# 6. Pillow (Phase 1 dep)
if python3 -c "from PIL import Image" 2>/dev/null; then
    echo "[OK] Pillow"
else
    echo "[FAIL] Pillow missing. Install: pip install --user Pillow" >&2
    errors=$((errors+1))
fi

if [ "$errors" -gt 0 ]; then
    echo "" >&2
    echo "$errors dependency check(s) failed." >&2
    exit 1
fi
```

```bash
chmod +x .claude/skills/html-preview/tools/preflight.sh
```

- [ ] **Step 4: Run test (expect PASS)**

```bash
bash .claude/skills/html-preview/tests/test_preflight.sh
```
Expected: `PASS: preflight.sh returns 0 in healthy env`.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/html-preview/tools/preflight.sh .claude/skills/html-preview/tests/test_preflight.sh
git commit -m "feat(html-preview): preflight.sh — 6-dep check (python3.10+ markdown jinja2 bs4 pyyaml Pillow)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: design_tokens.css — Anthropic palette + 8-point spacing (D16)

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/templates/design_tokens.css`

This file is 100% standalone (no test needed beyond visual review at Task 15).

- [ ] **Step 1: Write design_tokens.css**

Create `InterSubMod/.claude/skills/html-preview/templates/design_tokens.css`:

```css
/* design_tokens.css — Anthropic palette + 8-point spacing (D16, per SOP §5.1)
 * To be inlined inside <style> blocks by base_lN templates. */

:root {
  /* === Color palette (from thariqs design-system 05) === */
  --color-accent: #D97757;
  --color-text: #141413;
  --color-bg: #FAF9F5;
  --color-border: #E3DACC;
  --color-muted: #87867F;
  --color-success: #788C5D;
  --color-warning: #C78E3F;
  --color-danger: #B04A4A;
  --color-link: #5C7CA3;

  /* === Spacing (8-point) === */
  --sp-1: 4px;
  --sp-2: 8px;
  --sp-3: 16px;
  --sp-4: 24px;
  --sp-5: 32px;
  --sp-6: 48px;
  --sp-7: 64px;
  --sp-8: 96px;

  /* === Typography === */
  --font-sans: system-ui, -apple-system, "Noto Sans CJK TC", sans-serif;
  --font-mono: ui-monospace, "JetBrains Mono", "Cascadia Mono", monospace;
  --fs-display: 48px;
  --fs-h1: 32px;
  --fs-h2: 24px;
  --fs-h3: 20px;
  --fs-body: 16px;
  --fs-small: 14px;
  --fs-caption: 12px;

  /* === Layout === */
  --max-content: 900px;
  --radius-sm: 3px;
  --radius-md: 6px;
}
```

- [ ] **Step 2: Verify CSS syntax**

```bash
python3 -c "
import re
with open('.claude/skills/html-preview/templates/design_tokens.css') as f:
    content = f.read()
# Basic syntactic validity: balanced braces
opens = content.count('{')
closes = content.count('}')
assert opens == closes, f'Brace mismatch: {opens} open vs {closes} close'
# All --color-* and --sp-* and --fs-* and --font-* defined
required = ['--color-accent', '--color-text', '--color-bg', '--sp-1', '--sp-3', '--fs-h1', '--font-sans']
for r in required:
    assert r in content, f'Missing token: {r}'
print('[OK] design_tokens.css valid')
"
```
Expected: `[OK] design_tokens.css valid`.

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/html-preview/templates/design_tokens.css
git commit -m "feat(html-preview): design_tokens.css — D16 Anthropic palette + 8-point spacing

Tokens: 9 colors (#D97757 accent / #141413 text / #FAF9F5 bg / 6 semantic),
8-point spacing (4-96px), system-ui + JetBrains Mono fonts, 6 size scale.
All inlined as CSS variables (no external file dependency at runtime).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: base_l2.html.jinja — MVP.css template (DEFAULT per D15)

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/templates/base_l2.html.jinja`

- [ ] **Step 1: Write base_l2.html.jinja**

Create `InterSubMod/.claude/skills/html-preview/templates/base_l2.html.jinja`:

```jinja
<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ title | e }}</title>
<style>
{{ design_tokens_css }}

/* === MVP.css base (semantic HTML, D15) === */
* { box-sizing: border-box; }
body {
  font-family: var(--font-sans);
  color: var(--color-text);
  background: var(--color-bg);
  max-width: var(--max-content);
  margin: var(--sp-5) auto;
  padding: 0 var(--sp-3);
  line-height: 1.55;
}
h1 {
  font-size: var(--fs-h1); line-height: 1.2; font-weight: 500;
  border-bottom: 2px solid var(--color-text);
  padding-bottom: var(--sp-2);
  margin-top: var(--sp-6);
}
h2 {
  font-size: var(--fs-h2); line-height: 1.3; font-weight: 500;
  margin-top: var(--sp-5); color: var(--color-text);
}
h3 { font-size: var(--fs-h3); line-height: 1.4; font-weight: 500; margin-top: var(--sp-4); }
h4, h5, h6 { font-weight: 500; margin-top: var(--sp-3); }
p { margin: var(--sp-3) 0; }
a { color: var(--color-link); text-decoration: underline; }
a:hover { color: var(--color-accent); }
code {
  font-family: var(--font-mono); font-size: var(--fs-small);
  background: var(--color-border); padding: 2px 4px; border-radius: var(--radius-sm);
}
pre {
  background: var(--color-border); padding: var(--sp-3);
  border-radius: var(--radius-md); overflow-x: auto;
  font-size: var(--fs-small); line-height: 1.4;
}
pre code { background: none; padding: 0; font-size: inherit; }
table { border-collapse: collapse; width: 100%; margin: var(--sp-3) 0; }
th, td { border: 1px solid var(--color-border); padding: var(--sp-2); text-align: left; }
th { background: var(--color-border); font-weight: 500; }
blockquote {
  border-left: 4px solid var(--color-accent);
  margin: var(--sp-3) 0; padding-left: var(--sp-3);
  color: var(--color-muted);
}
ul, ol { padding-left: var(--sp-4); }
li { margin: var(--sp-1) 0; }
img { max-width: 100%; height: auto; }

/* === Collapsible (Markdown <details> survival) === */
details {
  border: 1px solid var(--color-border); border-radius: var(--radius-md);
  padding: var(--sp-2); margin: var(--sp-2) 0;
}
summary { cursor: pointer; font-weight: 500; padding: var(--sp-1); }
details[open] { padding: var(--sp-3); }

/* === Vision check / status badges === */
.badge { display: inline-block; padding: 2px 8px; border-radius: var(--radius-sm);
         font-size: var(--fs-small); font-weight: 500; }
.badge-pass { background: #D1FAE5; color: #065F46; }
.badge-partial { background: #FEF3C7; color: #92400E; }
.badge-fail { background: #FEE2E2; color: #991B1B; }
.badge-info { background: var(--color-border); color: var(--color-text); }

/* === Figure (image + caption) === */
.figure {
  margin: var(--sp-4) 0; padding: var(--sp-3);
  border: 1px solid var(--color-border); border-radius: var(--radius-md);
}
.figure img { display: block; margin: 0 auto; }
.caption { margin-top: var(--sp-2); font-size: var(--fs-small); color: var(--color-muted); }

/* === Print mode (D15 print-friendly) === */
@media print {
  body { max-width: none; background: white; margin: 1cm; }
  .figure { break-inside: avoid; border: 1px solid #ccc; }
  details { border: none; padding: 0; }
  details > summary { display: none; }
  details > * { display: block !important; }
  a { color: var(--color-text); text-decoration: none; }
  a::after { content: " (" attr(href) ")"; font-size: var(--fs-caption); color: var(--color-muted); }
}
</style>
</head>
<body>
<header>
  <p style="font-size: var(--fs-small); color: var(--color-muted);">
    Source: <a href="{{ source_md_relpath | e }}"><code>{{ source_md_relpath | e }}</code></a>
    · Generated by <code>html-preview</code> v{{ skill_version }}
    · {{ generated_at }}
  </p>
</header>

<main>
{{ body_html }}
</main>

<footer style="margin-top: var(--sp-7); padding-top: var(--sp-3); border-top: 1px solid var(--color-border); font-size: var(--fs-small); color: var(--color-muted);">
  Auto-generated companion HTML. Edit the source <code>.md</code> and re-run <code>/html-preview {{ source_md_relpath | e }}</code>.
</footer>
</body>
</html>
```

- [ ] **Step 2: Verify Jinja2 template renders**

```bash
python3 -c "
from jinja2 import Template
from pathlib import Path
tmpl_text = Path('.claude/skills/html-preview/templates/base_l2.html.jinja').read_text()
tokens = Path('.claude/skills/html-preview/templates/design_tokens.css').read_text()
t = Template(tmpl_text)
out = t.render(
    title='Test Report',
    design_tokens_css=tokens,
    source_md_relpath='docs/test.md',
    skill_version='2.0.0',
    generated_at='2026-05-11T10:00:00Z',
    body_html='<h1>Hello</h1><p>world</p>',
)
assert '<title>Test Report</title>' in out
assert '#D97757' in out, 'design tokens not inlined'
assert '<h1>Hello</h1>' in out, 'body not rendered'
assert '@media print' in out, 'print mode missing'
print(f'[OK] base_l2 template renders, output {len(out)} chars')
"
```
Expected: `[OK] base_l2 template renders, output ~5500-6000 chars`.

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/html-preview/templates/base_l2.html.jinja
git commit -m "feat(html-preview): base_l2.html.jinja — MVP.css default template (D15)

~3 KB inline CSS using design tokens, 0 external dep, 0 build step.
Includes: semantic typography, table styles, <details> collapse, badge
helpers (pass/partial/fail/info), figure with caption, print mode
(@media print expands all details, adds URL after links).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 5: md_to_html.py — markdown → HTML body (Python markdown lib)

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/md_to_html.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_md_to_html.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/fixtures/short_report.md`

- [ ] **Step 1: Write fixture short_report.md**

Create `InterSubMod/.claude/skills/html-preview/tests/fixtures/short_report.md`:

```markdown
---
title: "Short Test Report"
date: 2026-05-11
---

# Short Test Report

This is a short report for unit testing the markdown converter.

## Background

Some background text with **bold** and *italic* and `code`.

## Results

| Metric | Value |
|--------|-------|
| AUC | 0.78 |
| F1 | 0.65 |

## Code

```python
def hello():
    print("hello")
```

## Conclusion

A [link to docs](https://example.com).
```

- [ ] **Step 2: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_md_to_html.py`:

```python
"""Test md_to_html: markdown → HTML body."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from md_to_html import convert, parse_frontmatter  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_frontmatter():
    md = (FIXTURES / "short_report.md").read_text()
    fm, body = parse_frontmatter(md)
    assert fm.get("title") == "Short Test Report", f"Bad fm: {fm}"
    assert fm.get("date") == "2026-05-11", f"Bad date: {fm.get('date')}"
    assert "# Short Test Report" in body, "Body missing H1"


def test_convert_basic():
    md = (FIXTURES / "short_report.md").read_text()
    html = convert(md)
    assert "<h1>Short Test Report</h1>" in html
    assert "<table>" in html, "Table not rendered"
    assert "<code>" in html, "Inline code not rendered"
    assert "<pre>" in html, "Code block not rendered"
    assert "<strong>bold</strong>" in html, "Bold not rendered"
    assert '<a href="https://example.com">link to docs</a>' in html


def test_convert_strips_frontmatter():
    md = (FIXTURES / "short_report.md").read_text()
    html = convert(md)
    # Frontmatter (date, title in YAML form) should NOT appear in body
    assert "---" not in html, "YAML separator leaked into HTML"


def test_convert_handles_tables():
    md = "| A | B |\n|---|---|\n| 1 | 2 |\n"
    html = convert(md)
    assert "<table>" in html
    assert "<th>A</th>" in html
    assert "<td>1</td>" in html


def test_convert_handles_details():
    # GitHub-flavored details survive
    md = "<details><summary>Title</summary>\n\nbody\n\n</details>"
    html = convert(md)
    assert "<details>" in html
    assert "<summary>Title</summary>" in html


def main():
    tests = [
        test_parse_frontmatter, test_convert_basic, test_convert_strips_frontmatter,
        test_convert_handles_tables, test_convert_handles_details,
    ]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 3: Run test (expect FAIL — ImportError)**

```bash
python3 .claude/skills/html-preview/tests/test_md_to_html.py
```
Expected: `ImportError: cannot import name 'convert'`.

- [ ] **Step 4: Write md_to_html.py**

Create `InterSubMod/.claude/skills/html-preview/tools/md_to_html.py`:

```python
"""md_to_html — convert markdown source to HTML body fragment.

Uses the Python `markdown` library with extensions for tables, fenced code,
and pass-through HTML (so `<details>` survives).
Strips YAML frontmatter (returned separately for template variables).
"""
from __future__ import annotations

import re
from typing import Tuple

import markdown
import yaml

# Markdown extensions: tables, fenced code, no-implicit-paragraphs around HTML
_EXTENSIONS = ["tables", "fenced_code", "md_in_html"]
_FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?)\n---\s*\n", re.DOTALL)


def parse_frontmatter(md_text: str) -> Tuple[dict, str]:
    """Return (frontmatter_dict, body_without_frontmatter)."""
    m = _FRONTMATTER_RE.match(md_text)
    if not m:
        return {}, md_text
    fm_yaml = m.group(1)
    body = md_text[m.end():]
    try:
        fm = yaml.safe_load(fm_yaml) or {}
    except yaml.YAMLError:
        fm = {}
    if not isinstance(fm, dict):
        fm = {}
    # Normalize date to ISO string if datetime
    if "date" in fm and not isinstance(fm["date"], str):
        fm["date"] = str(fm["date"])
    return fm, body


def convert(md_text: str) -> str:
    """Convert markdown source to HTML body. Strips YAML frontmatter."""
    _, body = parse_frontmatter(md_text)
    return markdown.markdown(body, extensions=_EXTENSIONS, output_format="html")


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: md_to_html.py <input.md> [<output.html>]")
        return 2
    from pathlib import Path
    src = Path(argv[1])
    md_text = src.read_text()
    fm, _ = parse_frontmatter(md_text)
    html = convert(md_text)
    if len(argv) >= 3:
        Path(argv[2]).write_text(html)
        print(f"[md_to_html] {len(html)} chars → {argv[2]} (title={fm.get('title')!r})")
    else:
        print(html)
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 5: Run test (expect PASS)**

```bash
python3 .claude/skills/html-preview/tests/test_md_to_html.py
```
Expected: 5 PASS lines.

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/html-preview/tools/md_to_html.py \
        .claude/skills/html-preview/tests/test_md_to_html.py \
        .claude/skills/html-preview/tests/fixtures/short_report.md
git commit -m "feat(html-preview): md_to_html.py — Python markdown lib converter

Extensions: tables, fenced_code, md_in_html (preserves <details>).
parse_frontmatter strips YAML; convert returns body HTML only.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 6: topic_folder_decider.py — D17 conditional logic

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/topic_folder_decider.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_topic_folder_decider.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/fixtures/medium_report.md`
- Create: `InterSubMod/.claude/skills/html-preview/tests/fixtures/large_report.md`

- [ ] **Step 1: Write fixtures**

Create `InterSubMod/.claude/skills/html-preview/tests/fixtures/medium_report.md` (180 lines, 3 figs):

```markdown
---
title: Medium Test Report
---

# Medium Test Report

Stub body to reach 180 lines and contain 3 figures.

![](figures/fig1.png)
![](figures/fig2.png)
<!-- figure-needed: concept_diagram, slug=fig3 -->
```
(and 175 more lines of `text\n` to reach 180)

Use this script to generate fixture deterministically (no manual line-counting needed):
```bash
{ printf -- "---\ntitle: Medium Test Report\n---\n\n# Medium Test Report\n\nStub body.\n\n![](figures/fig1.png)\n![](figures/fig2.png)\n<!-- figure-needed: concept_diagram, slug=fig3 -->\n\n"; for i in $(seq 1 170); do printf "Line %d filler text.\n" "$i"; done; } > .claude/skills/html-preview/tests/fixtures/medium_report.md
```

Create `InterSubMod/.claude/skills/html-preview/tests/fixtures/large_report.md` (250 lines, 5 figures):

```bash
{ printf -- "---\ntitle: Large Test Report\n---\n\n# Large Report\n\n## Section 1\n![](figures/a.png)\n## Section 2\n![](figures/b.png)\n## Section 3\n![](figures/c.png)\n## Section 4\n![](figures/d.png)\n## Section 5\n<!-- figure-needed: icon, slug=e -->\n\n"; for i in $(seq 1 235); do printf "Line %d filler text.\n" "$i"; done; } > .claude/skills/html-preview/tests/fixtures/large_report.md
```

- [ ] **Step 2: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_topic_folder_decider.py`:

```python
"""Test topic_folder_decider applies D17 logic."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from topic_folder_decider import should_use_topic_folder  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_short_report_single_file():
    short = FIXTURES / "short_report.md"
    assert should_use_topic_folder(short) is False, "Short report (~30 lines, 0 figs) should be single-file"


def test_medium_report_single_file():
    medium = FIXTURES / "medium_report.md"
    # 180 lines, 3 figures — should be SINGLE FILE per D17 (need ≥200 lines OR ≥5 figs)
    assert should_use_topic_folder(medium) is False, \
        "Medium (180 lines, 3 figs) should be single-file"


def test_large_report_topic_folder():
    large = FIXTURES / "large_report.md"
    # 250 lines, 5 figures — exceeds both thresholds
    assert should_use_topic_folder(large) is True, \
        "Large (250 lines, 5 figs) should be topic-folder"


def test_lines_threshold_only():
    # Hypothetical: 220 lines, 0 figures
    text = "---\ntitle: x\n---\n\n# x\n\n" + "\n".join(f"line {i}" for i in range(220))
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(text)
        f.flush()
        path = Path(f.name)
    try:
        assert should_use_topic_folder(path) is True, ">200 lines should trigger"
    finally:
        path.unlink()


def test_figures_threshold_only():
    # Hypothetical: 30 lines, 6 figures
    text = "---\ntitle: x\n---\n\n# x\n\n" + "\n".join(f"![](f{i}.png)" for i in range(6))
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(text)
        f.flush()
        path = Path(f.name)
    try:
        assert should_use_topic_folder(path) is True, ">=5 figures should trigger"
    finally:
        path.unlink()


def main():
    tests = [test_short_report_single_file, test_medium_report_single_file,
             test_large_report_topic_folder, test_lines_threshold_only, test_figures_threshold_only]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 3: Run test (expect FAIL — ImportError)**

```bash
python3 .claude/skills/html-preview/tests/test_topic_folder_decider.py
```
Expected: `ImportError`.

- [ ] **Step 4: Write topic_folder_decider.py**

Create `InterSubMod/.claude/skills/html-preview/tools/topic_folder_decider.py`:

```python
"""topic_folder_decider — D17 logic: should_use_topic_folder?

Rule: report >= 200 lines OR >= 5 figures → topic folder.
Otherwise: single companion .preview.html file.

Figure count = markdown ![]() images + <!-- figure-needed: ... --> markers.
"""
from __future__ import annotations

from pathlib import Path

LINE_THRESHOLD = 200
FIGURE_THRESHOLD = 5


def count_figures(md_text: str) -> int:
    """Count both markdown !() images and <!-- figure-needed: ... --> markers."""
    img_count = md_text.count("![](") + md_text.count("![ ")  # ![alt](src) or ![ alt
    # More robust: regex
    import re
    img_count = len(re.findall(r"!\[[^\]]*\]\([^)]+\)", md_text))
    figure_needed_count = md_text.count("<!-- figure-needed:")
    return img_count + figure_needed_count


def should_use_topic_folder(md_path: Path) -> bool:
    """Return True if the report meets D17 conditions for topic folder."""
    text = Path(md_path).read_text()
    line_count = text.count("\n") + (0 if text.endswith("\n") else 1)
    fig_count = count_figures(text)
    return line_count >= LINE_THRESHOLD or fig_count >= FIGURE_THRESHOLD


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: topic_folder_decider.py <md_file>")
        return 2
    p = Path(argv[1])
    if not p.exists():
        print(f"Not found: {p}", file=sys.stderr)
        return 2
    text = p.read_text()
    lc = text.count("\n") + (0 if text.endswith("\n") else 1)
    fc = count_figures(text)
    decision = should_use_topic_folder(p)
    print(f"lines={lc} figures={fc} → topic_folder={decision}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 5: Run test (expect 5 PASS)**

```bash
python3 .claude/skills/html-preview/tests/test_topic_folder_decider.py
```
Expected: 5 PASS lines.

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/html-preview/tools/topic_folder_decider.py \
        .claude/skills/html-preview/tests/test_topic_folder_decider.py \
        .claude/skills/html-preview/tests/fixtures/medium_report.md \
        .claude/skills/html-preview/tests/fixtures/large_report.md
git commit -m "feat(html-preview): topic_folder_decider.py — D17 conditional split

Rule: lines >= 200 OR figures >= 5 → topic folder; else single file.
Figure count includes both ![]() and <!-- figure-needed: --> markers.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 7: interactivity_detect.py — L2 vs L3 routing

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/interactivity_detect.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_interactivity_detect.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/fixtures/interactive_marker.md`

- [ ] **Step 1: Write fixture**

```bash
cat > .claude/skills/html-preview/tests/fixtures/interactive_marker.md <<'EOF'
---
title: Interactive Demo
---

# Interactive Demo

<!-- interactive: tab -->

This report uses interactive tabs to compare options.
EOF
```

- [ ] **Step 2: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_interactivity_detect.py`:

```python
"""Test interactivity_detect: returns L2 (default) or L3 (interactive)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from interactivity_detect import detect_level  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_short_report_is_l2():
    p = FIXTURES / "short_report.md"
    assert detect_level(p) == "l2", "Default should be L2"


def test_interactive_marker_is_l3():
    p = FIXTURES / "interactive_marker.md"
    assert detect_level(p) == "l3", "<!-- interactive: tab --> should trigger L3"


def test_explicit_l1_override():
    # Tests that --force-l1 override works
    import tempfile
    text = "<!-- html-preview: force-l1 -->\n# x"
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(text)
        f.flush()
        p = Path(f.name)
    try:
        assert detect_level(p) == "l1"
    finally:
        p.unlink()


def main():
    tests = [test_short_report_is_l2, test_interactive_marker_is_l3, test_explicit_l1_override]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 3: Run test (expect FAIL)**

```bash
python3 .claude/skills/html-preview/tests/test_interactivity_detect.py
```
Expected: `ImportError`.

- [ ] **Step 4: Write interactivity_detect.py**

Create `InterSubMod/.claude/skills/html-preview/tools/interactivity_detect.py`:

```python
"""interactivity_detect — choose L1/L2/L3 rendering level.

L1: raw HTML (only when explicitly forced via <!-- html-preview: force-l1 -->)
L2: MVP.css default (D15)
L3: + Alpine.js when <!-- interactive: tab|filter|slider --> is present
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

Level = Literal["l1", "l2", "l3"]

_INTERACTIVE_MARKER = re.compile(r"<!--\s*interactive\s*:", re.IGNORECASE)
_FORCE_L1_MARKER = re.compile(r"<!--\s*html-preview\s*:\s*force-l1\s*-->", re.IGNORECASE)


def detect_level(md_path: Path) -> Level:
    text = Path(md_path).read_text()
    if _FORCE_L1_MARKER.search(text):
        return "l1"
    if _INTERACTIVE_MARKER.search(text):
        return "l3"
    return "l2"


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: interactivity_detect.py <md_file>")
        return 2
    print(detect_level(Path(argv[1])))
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 5: Run test (expect 3 PASS)**

```bash
python3 .claude/skills/html-preview/tests/test_interactivity_detect.py
```

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/html-preview/tools/interactivity_detect.py \
        .claude/skills/html-preview/tests/test_interactivity_detect.py \
        .claude/skills/html-preview/tests/fixtures/interactive_marker.md
git commit -m "feat(html-preview): interactivity_detect.py — L1/L2/L3 routing

Default L2 (MVP.css). L3 when <!-- interactive: ... --> marker found.
L1 only via explicit override marker. Aligned with D15 (L2 default) + D3 (L3 trigger).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 8: inline_assets.py — PNG → Base64 (reuses Phase 1 PIL)

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/inline_assets.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_inline_assets.py`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_inline_assets.py`:

```python
"""Test inline_assets: replace <img src="..."> with data: URIs."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from inline_assets import inline_images  # noqa
from PIL import Image


def make_test_png(path: Path):
    Image.new("RGB", (16, 16), (50, 100, 200)).save(path, "PNG")


def test_inlines_existing_relative_png():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        png = d / "fig.png"
        make_test_png(png)
        html = '<img src="fig.png" alt="x">'
        out = inline_images(html, base_dir=d)
        assert "data:image/png;base64," in out, f"Not inlined: {out}"
        assert 'src="fig.png"' not in out, "Original src not replaced"
        assert 'alt="x"' in out, "Other attributes preserved"


def test_keeps_external_url_unchanged():
    with tempfile.TemporaryDirectory() as tmp:
        html = '<img src="https://example.com/img.png">'
        out = inline_images(html, base_dir=Path(tmp))
        assert out == html, "External URL should be untouched"


def test_keeps_data_uri_unchanged():
    with tempfile.TemporaryDirectory() as tmp:
        html = '<img src="data:image/png;base64,abc">'
        out = inline_images(html, base_dir=Path(tmp))
        assert out == html, "data: URI should be untouched"


def test_missing_file_keeps_src_with_warning_comment():
    with tempfile.TemporaryDirectory() as tmp:
        html = '<img src="missing.png">'
        out = inline_images(html, base_dir=Path(tmp))
        assert 'src="missing.png"' in out, "Should keep src on miss"
        assert "<!--" in out and "missing" in out.lower(), "Should warn via comment"


def main():
    tests = [test_inlines_existing_relative_png, test_keeps_external_url_unchanged,
             test_keeps_data_uri_unchanged, test_missing_file_keeps_src_with_warning_comment]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Run test (expect FAIL)**

```bash
python3 .claude/skills/html-preview/tests/test_inline_assets.py
```
Expected: `ImportError`.

- [ ] **Step 3: Write inline_assets.py**

Create `InterSubMod/.claude/skills/html-preview/tools/inline_assets.py`:

```python
"""inline_assets — replace <img src="local.png"> with data: URIs.

Skips external URLs (http/https) and existing data: URIs. Inlines relative
paths resolved from base_dir. Missing files: keeps src + adds HTML comment.
"""
from __future__ import annotations

import base64
import re
from pathlib import Path

_IMG_RE = re.compile(r'<img\s+[^>]*?src="([^"]+)"[^>]*?>', re.IGNORECASE)


def _inline_one(match: re.Match, base_dir: Path) -> str:
    full_tag = match.group(0)
    src = match.group(1)
    if src.startswith(("http://", "https://", "data:")):
        return full_tag
    src_path = (base_dir / src).resolve()
    if not src_path.exists():
        return f"<!-- inline_assets: file not found at {src} relative to {base_dir} -->" + full_tag
    encoded = base64.b64encode(src_path.read_bytes()).decode("ascii")
    suffix = src_path.suffix.lower().lstrip(".")
    mime = {"png": "png", "jpg": "jpeg", "jpeg": "jpeg", "gif": "gif", "svg": "svg+xml"}.get(suffix, suffix)
    data_uri = f"data:image/{mime};base64,{encoded}"
    return full_tag.replace(f'src="{src}"', f'src="{data_uri}"', 1)


def inline_images(html: str, *, base_dir: Path) -> str:
    return _IMG_RE.sub(lambda m: _inline_one(m, base_dir), html)


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 3:
        print("Usage: inline_assets.py <input.html> <base_dir>")
        return 2
    html = Path(argv[1]).read_text()
    base = Path(argv[2])
    out = inline_images(html, base_dir=base)
    print(out)
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 4: Run test (expect 4 PASS)**

```bash
python3 .claude/skills/html-preview/tests/test_inline_assets.py
```

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/html-preview/tools/inline_assets.py .claude/skills/html-preview/tests/test_inline_assets.py
git commit -m "feat(html-preview): inline_assets.py — PNG/JPG/SVG → base64 data URI

Skips external URLs and existing data URIs. Missing files: warn comment.
SVG → image/svg+xml mime; others → image/<ext>.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 9: topic_folder_builder.py — D17 folder + README auto-gen

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/topic_folder_builder.py`
- Create: `InterSubMod/.claude/skills/html-preview/components/readme_template.md.jinja`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_topic_folder_builder.py`

- [ ] **Step 1: Write README template**

Create `InterSubMod/.claude/skills/html-preview/components/readme_template.md.jinja`:

```jinja
# {{ basename }}

> Auto-generated by `html-preview` skill v{{ skill_version }} at {{ generated_at }}

## Source

Markdown source: [`../{{ basename }}.md`](../{{ basename }}.md)

## Contents

- Main entry: [`index.html`](./index.html)
{% for ch in chapters -%}
- Chapter {{ ch.num }}: [`{{ ch.filename }}`]({{ ch.filename }}) — {{ ch.title }}
{% endfor %}

## Figures

{% if figures -%}
| # | File | Status |
|---|------|--------|
{% for fig in figures -%}
| {{ loop.index }} | `figures/{{ fig.filename }}` | {{ fig.status }} |
{% endfor %}
{% else -%}
No figures present in this report.
{%- endif %}

## Rebuild

```bash
/html-preview ../{{ basename }}.md --rebuild
```

To re-generate any AI figures:
```bash
/image-gen ./prompts/ ./figures/
```
```

- [ ] **Step 2: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_topic_folder_builder.py`:

```python
"""Test topic_folder_builder creates folder + README."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from topic_folder_builder import build_folder  # noqa


def test_creates_folder_with_readme():
    with tempfile.TemporaryDirectory() as tmp:
        md = Path(tmp) / "test_report.md"
        md.write_text("---\ntitle: Test\n---\n\n# Test\n")
        folder = build_folder(md, chapters=[
            {"num": "01", "filename": "ch01_intro.html", "title": "Intro"},
            {"num": "02", "filename": "ch02_methods.html", "title": "Methods"},
        ], figures=[])
        assert folder == md.parent / "test_report"
        assert folder.exists() and folder.is_dir()
        readme = folder / "README.md"
        assert readme.exists()
        readme_text = readme.read_text()
        assert "test_report" in readme_text
        assert "ch01_intro.html" in readme_text
        assert "Methods" in readme_text


def test_creates_subdirs_prompts_figures():
    with tempfile.TemporaryDirectory() as tmp:
        md = Path(tmp) / "x.md"
        md.write_text("# x")
        folder = build_folder(md, chapters=[], figures=[])
        assert (folder / "prompts").is_dir()
        assert (folder / "figures").is_dir()


def test_includes_figures_in_readme():
    with tempfile.TemporaryDirectory() as tmp:
        md = Path(tmp) / "x.md"
        md.write_text("# x")
        folder = build_folder(md, chapters=[], figures=[
            {"filename": "fig1.png", "status": "✓ generated"},
            {"filename": "fig2.png", "status": "pending"},
        ])
        readme_text = (folder / "README.md").read_text()
        assert "fig1.png" in readme_text
        assert "pending" in readme_text


def main():
    tests = [test_creates_folder_with_readme, test_creates_subdirs_prompts_figures, test_includes_figures_in_readme]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 3: Run test (expect FAIL)**

```bash
python3 .claude/skills/html-preview/tests/test_topic_folder_builder.py
```

- [ ] **Step 4: Write topic_folder_builder.py**

Create `InterSubMod/.claude/skills/html-preview/tools/topic_folder_builder.py`:

```python
"""topic_folder_builder — build {basename}/ topic folder structure (D17).

Creates: foo/ + foo/prompts/ + foo/figures/ + foo/README.md
README is rendered from components/readme_template.md.jinja.
"""
from __future__ import annotations

import datetime
from pathlib import Path
from typing import Iterable

from jinja2 import Template

_SKILL_DIR = Path(__file__).resolve().parent.parent
_README_TEMPLATE_PATH = _SKILL_DIR / "components" / "readme_template.md.jinja"
SKILL_VERSION = "2.0.0"


def build_folder(md_path: Path, *, chapters: Iterable[dict], figures: Iterable[dict]) -> Path:
    """Build topic folder next to md_path. Returns folder Path.

    chapters: list of {"num": "01", "filename": "ch01_x.html", "title": "..."}
    figures: list of {"filename": "fig1.png", "status": "..."}
    """
    basename = md_path.stem
    folder = md_path.parent / basename
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "prompts").mkdir(exist_ok=True)
    (folder / "figures").mkdir(exist_ok=True)
    (folder / "figures" / ".gitkeep").touch(exist_ok=True)
    (folder / "prompts" / ".gitkeep").touch(exist_ok=True)

    tmpl_text = _README_TEMPLATE_PATH.read_text()
    rendered = Template(tmpl_text).render(
        basename=basename,
        skill_version=SKILL_VERSION,
        generated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        chapters=list(chapters),
        figures=list(figures),
    )
    (folder / "README.md").write_text(rendered)
    return folder


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: topic_folder_builder.py <md_file>")
        return 2
    md = Path(argv[1])
    folder = build_folder(md, chapters=[], figures=[])
    print(f"[topic_folder_builder] built {folder}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 5: Run test (expect 3 PASS)**

```bash
python3 .claude/skills/html-preview/tests/test_topic_folder_builder.py
```

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/html-preview/tools/topic_folder_builder.py \
        .claude/skills/html-preview/components/readme_template.md.jinja \
        .claude/skills/html-preview/tests/test_topic_folder_builder.py
git commit -m "feat(html-preview): topic_folder_builder.py + README jinja template (D17)

Creates {basename}/ next to source .md, with prompts/ + figures/ subdirs
and auto-gen README.md (lists source link, chapters, figures, rebuild cmds).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 10: design_taste_check.py — Anthropic taboo audit (D14, D18 evaluator)

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/design_taste_check.py`
- Create: `InterSubMod/.claude/skills/html-preview/prompts/design_taste_check_prompt.md`
- Create: `InterSubMod/.claude/skills/html-preview/components/design_taboos_audit.md`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_design_taste_check.py`

- [ ] **Step 1: Write design_taboos_audit.md reference**

Create `InterSubMod/.claude/skills/html-preview/components/design_taboos_audit.md`:

```markdown
# Anthropic / thariqs Design Taboos (audit reference)

Per [thariqs/html-effectiveness](https://thariqs.github.io/html-effectiveness/) FAQ:

| # | Taboo | Detection regex / signal |
|---|-------|--------------------------|
| 1 | gradient overuse (>=2 gradients) | `linear-gradient` or `radial-gradient` count >= 2 |
| 2 | glass morphism | `backdrop-filter:` or `backdrop-blur` |
| 3 | multi-indigo stacking (multiple indigo shades) | `#4F46E5` `#6366F1` `#818CF8` (>=2 of) or `--color-indigo` chains |
| 4 | emoji-decorated headers | `<h[1-6][^>]*>[^<]*[\u{1F300}-\u{1FAFF}]` |
| 5 | drop shadows on text | `text-shadow:` (any non-`none`) |
| 6 | glow effects | `filter: drop-shadow` outside icons |

Also flagged: garish neon colors (#FF00FF, #00FFFF), 3D effects (`transform: rotateY`), comic sans / decorative fonts.
```

- [ ] **Step 2: Write design_taste_check_prompt.md**

Create `InterSubMod/.claude/skills/html-preview/prompts/design_taste_check_prompt.md`:

```markdown
# Design Taste Check Prompt (used by Claude when auditing rendered HTML)

You are auditing a rendered HTML document for Anthropic / thariqs design taboos.
Read the inline `<style>` block AND the body markup. Score each taboo present (1) or absent (0).

## Output format (JSON ONLY)

```json
{
  "taboos_detected": {
    "gradient_overuse": <true|false>,
    "glass_morphism": <true|false>,
    "multi_indigo_stacking": <true|false>,
    "emoji_decorated_headers": <true|false>,
    "drop_shadow_text": <true|false>,
    "glow_effects": <true|false>
  },
  "n_taboos": <integer>,
  "verdict": <"clean"|"minor_issues"|"redesign">,
  "suggestions": [
    "<concrete edit, e.g., 'Remove gradient on .hero — use --color-bg flat fill'>"
  ]
}
```

## Decision rules

- `n_taboos == 0` → `verdict: clean`
- `n_taboos == 1` → `verdict: minor_issues`
- `n_taboos >= 2` → `verdict: redesign`

For each detected taboo, write a 1-sentence suggestion.
```

- [ ] **Step 3: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_design_taste_check.py`:

```python
"""Test design_taste_check static detection (regex-based, no Claude call)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from design_taste_check import detect_taboos_static  # noqa


def test_clean_html_passes():
    html = """
    <style>body { color: #141413; background: #FAF9F5; }</style>
    <h1>Title</h1>
    """
    result = detect_taboos_static(html)
    assert result["n_taboos"] == 0
    assert result["verdict"] == "clean"


def test_gradient_overuse_flagged():
    html = """
    <style>
    .a { background: linear-gradient(red, blue); }
    .b { background: linear-gradient(green, yellow); }
    </style>
    """
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["gradient_overuse"] is True
    assert result["n_taboos"] >= 1


def test_glass_morphism_flagged():
    html = '<style>.x { backdrop-filter: blur(10px); }</style>'
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["glass_morphism"] is True


def test_emoji_header_flagged():
    html = '<h1>🚀 Welcome 🎉</h1>'
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["emoji_decorated_headers"] is True


def test_drop_shadow_text_flagged():
    html = '<style>h1 { text-shadow: 2px 2px 4px black; }</style>'
    result = detect_taboos_static(html)
    assert result["taboos_detected"]["drop_shadow_text"] is True


def test_minor_issues_verdict_one_taboo():
    html = '<style>.x { backdrop-filter: blur(10px); }</style>'
    result = detect_taboos_static(html)
    assert result["verdict"] == "minor_issues"


def test_redesign_verdict_two_taboos():
    html = """
    <style>
    .x { backdrop-filter: blur(10px); }
    h1 { text-shadow: 2px 2px 4px black; }
    </style>
    """
    result = detect_taboos_static(html)
    assert result["verdict"] == "redesign"
    assert result["n_taboos"] >= 2


def main():
    tests = [test_clean_html_passes, test_gradient_overuse_flagged, test_glass_morphism_flagged,
             test_emoji_header_flagged, test_drop_shadow_text_flagged,
             test_minor_issues_verdict_one_taboo, test_redesign_verdict_two_taboos]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test (expect FAIL)**

```bash
python3 .claude/skills/html-preview/tests/test_design_taste_check.py
```

- [ ] **Step 5: Write design_taste_check.py**

Create `InterSubMod/.claude/skills/html-preview/tools/design_taste_check.py`:

```python
"""design_taste_check — regex-based static detection of Anthropic design taboos.

Per design doc D14 + D18 (evaluator-optimizer pattern). For deeper audit,
use Claude vision with `prompts/design_taste_check_prompt.md` separately.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import TypedDict


class TabooResult(TypedDict):
    taboos_detected: dict
    n_taboos: int
    verdict: str
    suggestions: list


_GRADIENT_RE = re.compile(r"(linear|radial)-gradient", re.IGNORECASE)
_GLASS_RE = re.compile(r"backdrop-filter\s*:|backdrop-blur", re.IGNORECASE)
_INDIGO_HEXES = ("#4F46E5", "#6366F1", "#818CF8", "#A5B4FC", "#C7D2FE")
_EMOJI_HEADER_RE = re.compile(
    r"<h[1-6][^>]*>[^<]*[\U0001F300-\U0001FAFF\U00002600-\U000027BF]",
    re.IGNORECASE,
)
_TEXT_SHADOW_RE = re.compile(r"text-shadow\s*:\s*[^;n]", re.IGNORECASE)
_GLOW_RE = re.compile(r"filter\s*:\s*drop-shadow", re.IGNORECASE)


def detect_taboos_static(html: str) -> TabooResult:
    detected: dict = {}
    suggestions: list = []

    n_gradients = len(_GRADIENT_RE.findall(html))
    detected["gradient_overuse"] = n_gradients >= 2
    if detected["gradient_overuse"]:
        suggestions.append(f"Remove or consolidate gradients (found {n_gradients}); use flat var(--color-bg) fills.")

    detected["glass_morphism"] = bool(_GLASS_RE.search(html))
    if detected["glass_morphism"]:
        suggestions.append("Remove backdrop-filter / glass morphism; use solid backgrounds.")

    indigo_hits = sum(1 for hex_ in _INDIGO_HEXES if hex_.lower() in html.lower())
    detected["multi_indigo_stacking"] = indigo_hits >= 2
    if detected["multi_indigo_stacking"]:
        suggestions.append(f"Stop stacking indigo shades (found {indigo_hits}); pick one accent color from design tokens.")

    detected["emoji_decorated_headers"] = bool(_EMOJI_HEADER_RE.search(html))
    if detected["emoji_decorated_headers"]:
        suggestions.append("Remove emoji from <h1>-<h6> headers; reserve emojis for body content if any.")

    detected["drop_shadow_text"] = bool(_TEXT_SHADOW_RE.search(html))
    if detected["drop_shadow_text"]:
        suggestions.append("Remove text-shadow on text; use weight or color for emphasis instead.")

    detected["glow_effects"] = bool(_GLOW_RE.search(html))
    if detected["glow_effects"]:
        suggestions.append("Remove drop-shadow filters from non-icon elements.")

    n_taboos = sum(1 for v in detected.values() if v)
    if n_taboos == 0:
        verdict = "clean"
    elif n_taboos == 1:
        verdict = "minor_issues"
    else:
        verdict = "redesign"

    return {
        "taboos_detected": detected,
        "n_taboos": n_taboos,
        "verdict": verdict,
        "suggestions": suggestions,
    }


def main(argv: list[str]) -> int:
    import json
    import sys
    if len(argv) < 2:
        print("Usage: design_taste_check.py <html_file>")
        return 2
    html = Path(argv[1]).read_text()
    result = detect_taboos_static(html)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["verdict"] == "clean" else 1


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 6: Run test (expect 7 PASS)**

```bash
python3 .claude/skills/html-preview/tests/test_design_taste_check.py
```

- [ ] **Step 7: Commit**

```bash
git add .claude/skills/html-preview/tools/design_taste_check.py \
        .claude/skills/html-preview/prompts/design_taste_check_prompt.md \
        .claude/skills/html-preview/components/design_taboos_audit.md \
        .claude/skills/html-preview/tests/test_design_taste_check.py
git commit -m "feat(html-preview): design_taste_check.py — D14+D18 Anthropic taboo audit

Static regex detection of 6 taboos: gradient overuse, glass morphism,
multi-indigo stacking, emoji headers, text-shadow, glow effects.
Verdict: clean (0) / minor_issues (1) / redesign (>=2).
Companion Claude prompt for deeper visual audit.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 11: dispatch.py — main entry, end-to-end orchestration

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/tools/dispatch.py`
- Create: `InterSubMod/.claude/skills/html-preview/tests/test_dispatch.py`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/html-preview/tests/test_dispatch.py`:

```python
"""Test dispatch end-to-end: md → companion HTML (single-file or topic-folder)."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from dispatch import build_preview  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_short_report_produces_single_file():
    with tempfile.TemporaryDirectory() as tmp:
        # Copy short_report.md into tmp
        src = FIXTURES / "short_report.md"
        md = Path(tmp) / "short_report.md"
        md.write_text(src.read_text())

        result = build_preview(md)
        assert result["mode"] == "single_file", f"Got mode: {result['mode']}"
        out = result["output_path"]
        assert out.suffix == ".html"
        assert out.name == "short_report.preview.html"
        assert out.exists()
        html = out.read_text()
        assert "<title>Short Test Report</title>" in html
        assert "#D97757" in html, "design tokens missing"
        assert "<table>" in html, "table from md missing"


def test_large_report_produces_topic_folder():
    with tempfile.TemporaryDirectory() as tmp:
        src = FIXTURES / "large_report.md"
        md = Path(tmp) / "large_report.md"
        md.write_text(src.read_text())

        result = build_preview(md)
        assert result["mode"] == "topic_folder", f"Got mode: {result['mode']}"
        folder = result["output_path"]
        assert folder.is_dir()
        assert (folder / "index.html").exists()
        assert (folder / "README.md").exists()
        assert (folder / "figures").is_dir()
        assert (folder / "prompts").is_dir()


def test_taste_check_clean_passes():
    with tempfile.TemporaryDirectory() as tmp:
        src = FIXTURES / "short_report.md"
        md = Path(tmp) / "short_report.md"
        md.write_text(src.read_text())
        result = build_preview(md)
        assert result["taste_check"]["verdict"] == "clean", \
            f"Default template should be clean, got {result['taste_check']}"


def main():
    tests = [test_short_report_produces_single_file, test_large_report_produces_topic_folder,
             test_taste_check_clean_passes]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Run test (expect FAIL)**

```bash
python3 .claude/skills/html-preview/tests/test_dispatch.py
```

- [ ] **Step 3: Write dispatch.py**

Create `InterSubMod/.claude/skills/html-preview/tools/dispatch.py`:

```python
"""dispatch — main entry for html-preview skill.

Workflow:
  1. Detect topic-folder vs single-file (D17)
  2. Detect L2 vs L3 interactivity level (D15 default L2)
  3. Convert .md body → HTML (md_to_html)
  4. Render base_lN.html.jinja with design tokens inline
  5. Inline images (PNG → base64)
  6. Static taste audit (D14/D18)
  7. Write output (single .preview.html OR topic folder index.html)
"""
from __future__ import annotations

import datetime
import sys
from pathlib import Path

from jinja2 import Template

_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS))

import design_taste_check  # noqa: E402
import inline_assets  # noqa: E402
import interactivity_detect  # noqa: E402
import md_to_html  # noqa: E402
import topic_folder_builder  # noqa: E402
import topic_folder_decider  # noqa: E402

_SKILL_DIR = _TOOLS.parent
_TEMPLATES_DIR = _SKILL_DIR / "templates"
_DESIGN_TOKENS_PATH = _TEMPLATES_DIR / "design_tokens.css"
SKILL_VERSION = "2.0.0"


def _render(md_path: Path, level: str) -> str:
    """Render md to full HTML string at given level."""
    md_text = md_path.read_text()
    fm, _ = md_to_html.parse_frontmatter(md_text)
    body_html = md_to_html.convert(md_text)

    # Inline any local <img src="...">
    body_html = inline_assets.inline_images(body_html, base_dir=md_path.parent)

    template_path = _TEMPLATES_DIR / f"base_{level}.html.jinja"
    if not template_path.exists():
        # Fallback to L2 if requested level doesn't exist yet
        template_path = _TEMPLATES_DIR / "base_l2.html.jinja"
    tmpl_text = template_path.read_text()
    tokens_css = _DESIGN_TOKENS_PATH.read_text()

    rendered = Template(tmpl_text).render(
        title=fm.get("title", md_path.stem),
        design_tokens_css=tokens_css,
        source_md_relpath=md_path.name,
        skill_version=SKILL_VERSION,
        generated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        body_html=body_html,
    )
    return rendered


def build_preview(md_path: Path) -> dict:
    """Build companion HTML for md_path. Returns:
        {
          "mode": "single_file" | "topic_folder",
          "output_path": Path (file for single_file, folder for topic_folder),
          "level": "l1"|"l2"|"l3",
          "taste_check": {...}
        }
    """
    md_path = Path(md_path).resolve()
    use_folder = topic_folder_decider.should_use_topic_folder(md_path)
    level = interactivity_detect.detect_level(md_path)
    html = _render(md_path, level)
    taste = design_taste_check.detect_taboos_static(html)

    if use_folder:
        folder = topic_folder_builder.build_folder(
            md_path, chapters=[], figures=[],
        )
        index_path = folder / "index.html"
        index_path.write_text(html)
        return {
            "mode": "topic_folder",
            "output_path": folder,
            "level": level,
            "taste_check": taste,
        }
    else:
        out = md_path.parent / (md_path.stem + ".preview.html")
        out.write_text(html)
        return {
            "mode": "single_file",
            "output_path": out,
            "level": level,
            "taste_check": taste,
        }


def main(argv: list[str]) -> int:
    import argparse
    parser = argparse.ArgumentParser(prog="html-preview dispatch")
    parser.add_argument("md_path", type=Path)
    parser.add_argument("--rebuild", action="store_true",
                        help="overwrite existing companion (default: skip if exists)")
    parser.add_argument("--polished", action="store_true",
                        help="(D19, opt-in) use frontend-design plugin (Phase 2.5 not implemented yet)")
    args = parser.parse_args(argv[1:])

    if not args.md_path.exists():
        print(f"Source not found: {args.md_path}", file=sys.stderr)
        return 2
    if args.polished:
        print("[html-preview] --polished requires frontend-design plugin invocation; not yet wired in Phase 2 base. "
              "Falling back to MVP.css default.", file=sys.stderr)

    result = build_preview(args.md_path)
    print(f"[html-preview] mode={result['mode']} level={result['level']} → {result['output_path']}")
    print(f"[html-preview] taste_check: verdict={result['taste_check']['verdict']} "
          f"n_taboos={result['taste_check']['n_taboos']}")
    if result["taste_check"]["suggestions"]:
        for s in result["taste_check"]["suggestions"]:
            print(f"  - {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
```

- [ ] **Step 4: Run test (expect 3 PASS)**

```bash
python3 .claude/skills/html-preview/tests/test_dispatch.py
```

- [ ] **Step 5: Smoke test — convert short_report fixture**

```bash
mkdir -p /tmp/htmlprev_smoke
cp .claude/skills/html-preview/tests/fixtures/short_report.md /tmp/htmlprev_smoke/
python3 .claude/skills/html-preview/tools/dispatch.py /tmp/htmlprev_smoke/short_report.md
ls -la /tmp/htmlprev_smoke/
file /tmp/htmlprev_smoke/short_report.preview.html
```
Expected: produces `short_report.preview.html` ~5-10 KB; `file` reports `HTML document`.

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/html-preview/tools/dispatch.py .claude/skills/html-preview/tests/test_dispatch.py
git commit -m "feat(html-preview): dispatch.py — main entry, conditional split + taste check

Decides single-file vs topic-folder per D17, picks L2 default per D15,
inlines design tokens and PNG assets, runs static taste audit (D14/D18).
--polished flag stub for D19 frontend-design integration (Phase 2.5).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 12: SKILL.md + playbook (per Anthropic best-practices, SOP §4)

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/SKILL.md`
- Create: `InterSubMod/.claude/skills/html-preview/playbook.md`

- [ ] **Step 1: Write SKILL.md per directive form (SOP §4.2)**

Create `InterSubMod/.claude/skills/html-preview/SKILL.md`:

```markdown
---
name: html-preview
description: |
  Generate companion HTML viewer next to any markdown report by converting source to MVP.css-styled single-file HTML or, for large reports (>=200 lines OR >=5 figures), a topic folder with index.html + chapter files + README. Inlines design tokens (Anthropic palette + 8-point spacing), inlines PNG assets as base64, audits result against 6 design taboos (gradient overuse / glass morphism / multi-indigo / emoji headers / text-shadow / glow), and provides print-friendly @media print rules. Companion (not replacement) for the source .md.
  USE WHEN: 「想看看排版」「給 PI 看 preview」「快速確認感覺」「html preview」「報告 HTML 預覽」「companion HTML」「.md 預覽」「給人看的版本」、寫完報告後 PI 提到「想先看看」、weekly-report / structured-tech-report / results-report / feature-layered-observation / methodology-audit / conclude-research 結束時自動觸發。
  SKIP WHEN: README.md（GitHub 原生渲染即可）、Slack 片段、給其他 LLM 消費的 .md（如 CLAUDE.md / state.json / hypothesis_queue.json）、純個人筆記、CI 自動化文件、JSON / YAML / CSV state 檔。
---

# html-preview

## Phase & Chain Position

- Phase 2 of HTML preview + image-gen 3-skill design (Phase 1 = image-gen + image-vision-check shipped).
- Stand-alone, callable directly or via myPPT / weekly-report / structured-tech-report.
- Upstream: any .md report (especially Tier A 6 skills per SOP §7.2).
- Downstream: PI opens output in browser; Phase 3 接入會把這 skill 嵌進 6 個既有 skill 自動觸發。

## Dependencies

**Uses (runtime)**:
- `python3` >= 3.10 with `markdown >= 3.5`, `jinja2 >= 3.0`, `beautifulsoup4 >= 4.12`, `pyyaml`, `Pillow`
- No pandoc, no Node.js, no external CDN. All deps pure-python pip installable.

**Used by**:
- `myPPT` orchestration (Phase 3)
- 6 Tier A skills via Phase 3 接入 (opt-in flag)

**Reads**:
- Source `.md` file
- Any local `<img src="...">` referenced in the .md (resolved relative to .md dir)

**Writes**:
- Single-file mode: `{basename}.preview.html` next to source .md
- Topic-folder mode: `{basename}/` directory with index.html + README.md + figures/ + prompts/

## Failure Modes & Diagnostics

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `[FAIL] markdown missing` | pip dep not installed | `pip install --user 'markdown>=3.5'` |
| Output HTML > 1 MB | many large PNG inlines | reduce image size before run, or resize via `image-gen --postprocess-size` |
| `taste_check verdict=redesign` | template has taboos (gradient/glass/etc) | normally only happens if user added inline `<style>` to .md; review |
| Topic folder built when expecting single file | report >=200 lines OR >=5 figures | verify with `tools/topic_folder_decider.py <md>` |
| Frontmatter not parsed | non-YAML or missing `---` delimiters | ensure .md starts with `---\n<key>: <val>\n---\n\n` |

## Quick Usage

```bash
# Default: build companion HTML next to .md
python3 InterSubMod/.claude/skills/html-preview/tools/dispatch.py \
  InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md

# Force rebuild over existing
python3 InterSubMod/.claude/skills/html-preview/tools/dispatch.py <md> --rebuild

# Use frontend-design plugin (D19, Phase 2.5; currently logs warning + falls back)
python3 InterSubMod/.claude/skills/html-preview/tools/dispatch.py <md> --polished

# Just check what level / mode would be chosen
python3 InterSubMod/.claude/skills/html-preview/tools/topic_folder_decider.py <md>
python3 InterSubMod/.claude/skills/html-preview/tools/interactivity_detect.py <md>
```

## Pre-flight Required

```bash
bash InterSubMod/.claude/skills/html-preview/tools/preflight.sh
```

## Output Layout

**Single-file mode (D17 default for short reports)**:
```
docs/reports/validated/2026/05/
├── 20260511_short_report.md          # source (in git)
└── 20260511_short_report.preview.html # companion (in git)
```

**Topic-folder mode (>=200 lines OR >=5 figures)**:
```
docs/reports/validated/2026/05/
├── 20260511_large_report.md         # source
└── 20260511_large_report/            # topic folder
    ├── README.md                     # auto-gen, lists contents
    ├── index.html                    # main entry
    ├── prompts/                      # in git (figure prompts)
    └── figures/                      # .gitignore (PNG output)
```

## Design Tokens (D16)

All tokens inline as CSS variables (no external file at runtime):
- Palette: `#D97757` accent / `#141413` text / `#FAF9F5` bg / `#E3DACC` border
- Spacing: 8-point scale (--sp-1=4px → --sp-8=96px)
- Typography: system-ui + Noto Sans CJK TC + JetBrains Mono fallbacks
- Print mode: @media print expands all `<details>`, prints URLs after links

## Cost Model

| Operation | Cost |
|-----------|------|
| md → HTML conversion | 0 (local Python) |
| Static taste check | 0 (regex only) |
| Optional Claude vision audit (manual invoke) | ~$0.005 / page |

## See Also

- Spec: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md` (D1-D20)
- SOP: `InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md` (§5 templates, §6 evaluator-optimizer)
- Companion skills: `image-gen`, `image-vision-check` (Phase 1)
```

- [ ] **Step 2: Verify SKILL.md ≤500 lines (SOP §4.3)**

```bash
wc -l .claude/skills/html-preview/SKILL.md
```
Expected: well under 500.

- [ ] **Step 3: Write playbook.md**

Create `InterSubMod/.claude/skills/html-preview/playbook.md`:

```markdown
# html-preview Playbook

## Step-by-step workflow

### Step 1: Preflight
Run `tools/preflight.sh`. Halt on any `[FAIL]`.

### Step 2: Check if source exists
The argument must be an existing `.md` file. If folder, error.

### Step 3: Decide mode (D17)
Call `topic_folder_decider.should_use_topic_folder(md_path)`:
- True → topic-folder mode
- False → single-file mode (default)

### Step 4: Decide level (D15 modified by D3)
Call `interactivity_detect.detect_level(md_path)`:
- L2 → MVP.css (default)
- L3 → + Alpine.js (when `<!-- interactive: tab -->` present)
- L1 → raw HTML (only via explicit `<!-- html-preview: force-l1 -->`)

### Step 5: Convert + render
- `md_to_html.convert(text)` → body HTML
- `inline_assets.inline_images(html, base_dir)` → embed PNGs as base64
- `Template(base_lN).render(title, design_tokens_css, body_html, ...)` → full HTML

### Step 6: Audit (D14 + D18 evaluator-optimizer)
- `design_taste_check.detect_taboos_static(html)` → verdict
- If `clean` → write output
- If `minor_issues` → write output but warn PI in stdout
- If `redesign` → write output, escalate full checklist to PI

### Step 7: Write output
- single_file → write `{basename}.preview.html`
- topic_folder → call `topic_folder_builder.build_folder()` + write `index.html` + `README.md`

### Step 8: Suggest follow-up
- If figures still pending: suggest `/image-gen ./prompts/ ./figures/`
- If taste_check verdict != clean: list suggestions

## Anthropic Design Taboos (D14 enforcement)

Static audit checks for:
1. gradient overuse (>=2 gradient-* uses)
2. glass morphism (backdrop-filter)
3. multi-indigo stacking (>=2 indigo hex codes)
4. emoji-decorated headers
5. text-shadow on text
6. drop-shadow filters (non-icon)

Default L2 template guarantees `clean` verdict (no taboos in design tokens).
Verdict only fails if user/upstream injected custom inline styles.

## Known Edge Cases

- **Empty .md**: produces empty companion, no error
- **Malformed YAML frontmatter**: parse_frontmatter returns `{}`, body still renders
- **External `<img src="https://...">`**: NOT inlined; passed through as-is
- **SVG inline already in .md**: passed through unchanged
- **Nested folders for source .md**: companion always created next to source (same dir)
- **`<details>` blocks in .md**: preserved via `md_in_html` extension; styled by template
```

- [ ] **Step 4: Commit**

```bash
git add .claude/skills/html-preview/SKILL.md .claude/skills/html-preview/playbook.md
git commit -m "feat(html-preview): SKILL.md + playbook per Anthropic best-practices (SOP §4)

SKILL.md: directive description with USE WHEN / SKIP WHEN, third-person,
phase chain position, dependencies, failure modes, quick usage. <500 lines.
playbook.md: 8-step workflow, taboo audit explanation, edge cases.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 13: evals.json — 3 evaluation cases (per D20)

**Files:**
- Create: `InterSubMod/.claude/skills/html-preview/evals.json`

- [ ] **Step 1: Write evals.json**

Create `InterSubMod/.claude/skills/html-preview/evals.json`:

```json
{
  "skill": "html-preview",
  "version": "2.0.0",
  "schema_version": "1.0",
  "created_at": "2026-05-11",
  "evals": [
    {
      "id": "eval_01_typical_trigger",
      "category": "activation",
      "query": "幫我把這份報告轉成 HTML 給 PI 看 preview",
      "context": {
        "files_in_scope": ["InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md"]
      },
      "expected_behavior": [
        "Skill activates (USE WHEN matched: 'preview', '給 PI 看')",
        "dispatch.py invoked with the .md path",
        "Output is companion HTML (single-file or topic-folder per D17)",
        "Result printed includes mode + level + taste_check verdict"
      ],
      "must_not": [
        "Replace the source .md (companion only, D11)",
        "Use Tailwind by default (D15 says MVP.css)"
      ]
    },
    {
      "id": "eval_02_skip_scenario",
      "category": "skip",
      "query": "幫我整理 CLAUDE.md",
      "context": {
        "files_in_scope": ["InterSubMod/.claude/CLAUDE.md"]
      },
      "expected_behavior": [
        "Skill does NOT activate (CLAUDE.md is for LLM consumption per SKIP WHEN)",
        "Either no skill triggers, or doc-standards / claude-md-management:claude-md-improver activates instead"
      ],
      "must_not": [
        "Generate a CLAUDE.md.preview.html",
        "Auto-create a CLAUDE.md/ topic folder"
      ]
    },
    {
      "id": "eval_03_edge_case_large_report_topic_folder",
      "category": "edge_case",
      "query": "convert this 250-line report to HTML",
      "context": {
        "files_in_scope": ["any .md with >=200 lines OR >=5 figures"]
      },
      "expected_behavior": [
        "topic_folder_decider returns True",
        "build_folder creates {basename}/ next to source",
        "README.md auto-generated with source link + chapter list + figure table",
        "index.html created inside folder",
        "figures/ and prompts/ subdirs created with .gitkeep"
      ],
      "must_not": [
        "Create single .preview.html for large report (D17 violation)",
        "Forget to create the README.md (D12 explicit requirement)"
      ]
    }
  ],
  "notes": [
    "Evals are per Anthropic best-practices (SOP §10).",
    "Run manually via Claude Code with the skill loaded; observe behavior matches expected_behavior.",
    "Future: add scripts/run_evals.py for automated harness once eval framework standardized."
  ]
}
```

- [ ] **Step 2: Validate JSON**

```bash
python3 -c "import json; d=json.load(open('.claude/skills/html-preview/evals.json')); print(f'[OK] evals.json valid: {len(d[\"evals\"])} evals')"
```
Expected: `[OK] evals.json valid: 3 evals`.

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/html-preview/evals.json
git commit -m "feat(html-preview): evals.json — 3 cases per D20 (typical / skip / edge)

eval_01: typical trigger ('preview', '給 PI 看') → activate + companion HTML
eval_02: skip CLAUDE.md (LLM-consumption per SKIP WHEN)
eval_03: edge case — 250-line report routes to topic-folder mode

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 14: Run all tests + lint check before demo

- [ ] **Step 1: Run all tests**

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
echo "=== preflight ===" && bash .claude/skills/html-preview/tests/test_preflight.sh
echo "=== md_to_html ===" && python3 .claude/skills/html-preview/tests/test_md_to_html.py
echo "=== topic_folder_decider ===" && python3 .claude/skills/html-preview/tests/test_topic_folder_decider.py
echo "=== interactivity_detect ===" && python3 .claude/skills/html-preview/tests/test_interactivity_detect.py
echo "=== inline_assets ===" && python3 .claude/skills/html-preview/tests/test_inline_assets.py
echo "=== topic_folder_builder ===" && python3 .claude/skills/html-preview/tests/test_topic_folder_builder.py
echo "=== design_taste_check ===" && python3 .claude/skills/html-preview/tests/test_design_taste_check.py
echo "=== dispatch ===" && python3 .claude/skills/html-preview/tests/test_dispatch.py
```
Expected: All PASS lines, no FAIL or ERROR.

- [ ] **Step 2: Lint check (SOP §8 anti-pattern checklist)**

```bash
# A. Windows-style paths
echo "--- Windows paths ---"
grep -rn '\\\\' .claude/skills/html-preview/ && echo "❌ Windows paths found" || echo "✓ no Windows paths"

# B. Time-sensitive info (skip "Old patterns" wrap)
echo "--- Time-sensitive ---"
grep -rn -E "(after|before) [0-9]{4}" .claude/skills/html-preview/ | grep -v "Old patterns" | head -3 || echo "✓ no time-sensitive"

# C. SKILL.md > 500 lines
echo "--- SKILL.md size ---"
wc -l .claude/skills/html-preview/SKILL.md

# D. Description first-person violation
echo "--- first-person ---"
grep -E "description:.*( I | I'd | I'll | I am )" .claude/skills/html-preview/SKILL.md && echo "❌ first-person" || echo "✓ third-person"

# E. Voodoo constants in scripts
echo "--- voodoo constants ---"
grep -rn -E "(TIMEOUT|DELAY|RETRIES) = [0-9]+" .claude/skills/html-preview/tools/ | grep -v "#" | head -3 || echo "✓ no undocumented magic"
```
Expected: No `❌` lines.

- [ ] **Step 3: Commit (no changes if all green)**

If any lint findings produced fixes, commit:
```bash
git add -A .claude/skills/html-preview/
git commit -m "chore(html-preview): lint pass per SOP §8 anti-pattern checklist

(Optional task — only if any fixes were needed.)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```
If no fixes needed, skip this commit step.

---

### Task 15: End-to-end demo on real PI report + ship

**Goal**: Convert a real existing report to HTML preview, visually confirm output, commit demo.

- [ ] **Step 1: Pick demo source**

Use `InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md` (existing, ~150-300 lines).

```bash
wc -l docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md
# Note: the line count determines if it's single-file or topic-folder mode
```

- [ ] **Step 2: Run dispatch**

```bash
python3 .claude/skills/html-preview/tools/dispatch.py \
  docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.md
```

Expected output prints mode (`single_file` or `topic_folder`), level (`l2`), taste_check verdict (`clean`), and output path.

- [ ] **Step 3: Verify output exists and is valid HTML**

```bash
# If single_file mode:
ls -la docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.preview.html 2>/dev/null

# If topic_folder mode:
ls -la docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01/

# Either way: validate HTML syntax (head/tail check, basic structural)
python3 -c "
from pathlib import Path
candidates = [
    Path('docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.preview.html'),
    Path('docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01/index.html'),
]
for p in candidates:
    if p.exists():
        h = p.read_text()
        assert h.startswith('<!DOCTYPE html>'), 'Missing doctype'
        assert '#D97757' in h, 'Design tokens missing'
        assert '<title>' in h, 'No title'
        print(f'[OK] {p}: {len(h)} chars')
        break
else:
    print('[FAIL] no output found')
"
```
Expected: `[OK] <path>: ~10000-200000 chars` (depends on figure inlines).

- [ ] **Step 4: PI visual confirmation (manual)**

Open the output HTML in a browser:
```bash
# On local machine after scp:
# scp <user>@<server>:.../{filename}.preview.html /tmp/ && open /tmp/{filename}.preview.html

# Or print a Python sanity check of section count
python3 -c "
from pathlib import Path
import re
candidates = list(Path('docs/reports/validated/2026/05').glob('20260509_PI_Report_4_29_Errata_01*'))
for p in candidates:
    if p.suffix == '.html':
        h = p.read_text()
        h2_count = len(re.findall(r'<h2[ >]', h))
        print(f'{p.name}: {h2_count} sections')
"
```

PI checklist:
- [ ] Title bar shows correct .md title
- [ ] All sections present and styled (not raw markdown leaking)
- [ ] Tables render correctly with borders
- [ ] Code blocks have monospace font and grey background
- [ ] Links are blue (`#5C7CA3`); hover changes to accent (`#D97757`)
- [ ] Print preview (Ctrl+P) shows expanded `<details>` blocks
- [ ] No browser console errors

- [ ] **Step 5: Commit demo**

```bash
git add docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01*
# Note: this glob picks up either {basename}.preview.html OR {basename}/ folder
git commit -m "demo(phase2): html-preview applied to 20260509_PI_Report_4_29_Errata_01

Phase 2 closed-loop validation. Mode: <single_file|topic_folder> (per D17).
Level: l2 (MVP.css default per D15). Taste check: clean (D14/D18).

Source .md unchanged (companion-only per D11).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

- [ ] **Step 6: Update CURRENT_FOCUS.md**

Append (or insert) to `docs/CURRENT_FOCUS.md`:

```markdown

## 2026-05-11 — html-preview skill shipped (Phase 2)

✅ `InterSubMod/.claude/skills/html-preview/` shipped. End-to-end demo at
`InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01.preview.html`
(or `.../20260509_PI_Report_4_29_Errata_01/index.html` if topic-folder mode).

Spec: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md` (D1-D20)
SOP: `InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md`
Phase 2 plan: `InterSubMod/docs/plans/2026/05/20260511_Phase2_html_preview_implementation_01.md`

Ready to plan Phase 3 (Tier A 6 skill 接入 html-preview as companion).
```

```bash
git add docs/CURRENT_FOCUS.md
git commit -m "docs(focus): mark Phase 2 (html-preview) complete

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

- [ ] **Step 7: Add memory file**

Use Write tool to create `/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/project_phase2_html_preview_shipped.md`:

```markdown
---
name: Phase 2 html-preview skill shipped
description: 2026-05-11 Phase 2 完成。html-preview skill 用 MVP.css + design tokens 把 .md 轉 companion HTML; 大報告 (>=200 lines OR >=5 figs) 自動建主題資料夾; 內建設計禁忌 audit。
type: project
---

**2026-05-11 Phase 2 ship**

`InterSubMod/.claude/skills/html-preview/` 上線：
- 預設 MVP.css 內聯 (~3 KB CSS, D15 取代 Tailwind)
- 設計 tokens: Anthropic 主色 #D97757 / #141413 / #FAF9F5 + 8-point spacing (D16)
- 條件化主題資料夾 (>=200 lines OR >=5 figures) (D17)
- 6 維設計禁忌靜態 audit (D14/D18 evaluator-optimizer 模式)
- 3 evals 內建 (D20)
- frontend-design plugin opt-in (--polished, D19) — 暫 stub fallback

**Why**: PI 痛點 — .md 在 terminal 看不出排版、要分享難、列印不友善。
**How to apply**: 寫完 .md 後 `python3 .claude/skills/html-preview/tools/dispatch.py <md>`，產 companion HTML。
**Phase 3 待做**: 6 個 Tier A skill 接 html-preview 為自動 companion (opt-in flag); 補 evals + description audit per D20。

**Spec**: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md`
**SOP**: `InterSubMod/docs/references/manual/20260511_HTML_MD_PPTX輸出格式SOP_01.md`
**Plan**: `InterSubMod/docs/plans/2026/05/20260511_Phase2_html_preview_implementation_01.md`
```

Add pointer line to `MEMORY.md` under `## Active Research` or new `## Skills` section:
```
- [Phase 2 html-preview shipped](project_phase2_html_preview_shipped.md) — 2026-05-11 MVP.css + 條件化主題資料夾 + 6 維 taboo audit
```

---

## Self-Review

**1. Spec coverage** (D1-D20):
- ✓ D1, D2 covered by Phase 1 (image-gen + image-vision-check); html-preview is the third skill
- ✓ D3 (interactivity) updated by D15: L2 default + L3 trigger via marker — Task 7
- ✓ D4 (HTML lifecycle companion) updated by D17 conditional — Task 6, 9, 11
- ✓ D5, D6, D8, D14 image-gen related (Phase 1)
- ✓ D7 .gitignore for figures inherited from Phase 1
- ✓ D9 first-batch scope (Phase 2 = html-preview alone)
- ✓ D10, D11 (companion not replacement) — built into dispatch.py never touching source .md
- ✓ D12 folder naming/files — Tasks 9, 11 + readme template
- ✓ D13 codex OAuth (image-gen, Phase 1, not relevant here)
- ✓ D15 MVP.css default — Task 3 (tokens) + Task 4 (template) + Task 11 (dispatch picks L2)
- ✓ D16 design tokens — Task 3 (full CSS variables)
- ✓ D17 conditional topic folder — Task 6 (decider) + Task 9 (builder) + Task 11 (dispatch routing)
- ✓ D18 Evaluator-Optimizer — Task 10 (taste check) + Task 11 (dispatch invokes audit)
- ✓ D19 frontend-design opt-in — Task 11 dispatch has `--polished` flag stub (full integration deferred to Phase 2.5)
- ✓ D20 evals — Task 13 (3 evals.json)
- ✓ SOP §4 SKILL.md規範 — Task 12 (directive description, ≤500 lines, 3rd person)
- ✓ SOP §10 pre-flight checklist — implicit in plan structure

Missing (acceptable for Phase 2 scope):
- L3 Alpine.js template (`base_l3.html.jinja`) — not yet created; only L2 in Task 4
- frontend-design `--polished` full integration (Phase 2.5)
- SVG inline (currently relies on `md_in_html` extension; works for inline `<svg>` in source)

**2. Placeholder scan**: No TBD/TODO/handwave found. All test code, all impl code, all commands shown.

**3. Type / signature consistency**:
- `convert(md_text)` Task 5 → consumed by `dispatch._render` Task 11 ✓
- `parse_frontmatter(md_text)` Task 5 → consumed by `dispatch._render` Task 11 ✓
- `should_use_topic_folder(md_path)` Task 6 → consumed by `dispatch.build_preview` Task 11 ✓
- `detect_level(md_path) -> Literal["l1","l2","l3"]` Task 7 → consumed by `dispatch._render` Task 11 ✓
- `inline_images(html, *, base_dir)` Task 8 → consumed by `dispatch._render` Task 11 ✓
- `build_folder(md_path, *, chapters, figures) -> Path` Task 9 → consumed by `dispatch.build_preview` Task 11 ✓
- `detect_taboos_static(html) -> TabooResult` Task 10 → consumed by `dispatch.build_preview` Task 11 ✓
- `build_preview(md_path) -> dict` Task 11 → consumed by `dispatch.main` + tests ✓

All keyword-only args used consistently; all return types match consumers.

**4. Scope**: Phase 2 only. Phase 2.5 items (L3 template + frontend-design integration) explicitly deferred. Phase 3 (Tier A 接入) is out of scope.

---

## Execution Handoff

Plan complete. Two execution options:

**1. Subagent-Driven (recommended)** — Dispatch fresh subagent per task, two-stage review (spec + quality), continuous execution. Same pattern as Phase 1.

**2. Inline Execution** — Execute tasks in this session using `superpowers:executing-plans`, batch with checkpoints.

Which approach?
