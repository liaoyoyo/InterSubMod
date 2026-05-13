---
name: html-preview
description: |
  **[DEPRECATED 2026-05-13]** Superseded by `html-report-build` (LLM-direct, no Python middleware). New invocations should use html-report-build. This skill kept for 1 release as fallback only — see migration note below.
  Legacy purpose: Generate companion HTML viewer next to any markdown report by converting source to MVP.css-styled single-file HTML or, for large reports (>=200 lines OR >=5 figures), a topic folder with index.html + chapter files + README. Uses Python pipeline (markdown lib + jinja2 + bs4 + Pillow). Companion (not replacement) for the source .md.
  USE WHEN: **DO NOT** — use `html-report-build` instead. Only invoke explicitly with `--legacy` flag if html-report-build LLM-direct output unexpectedly fails and Python fallback needed.
  SKIP WHEN: All normal use cases — route to html-report-build.
---

# html-preview ⚠ DEPRECATED

> **[2026-05-13] SUPERSEDED BY `html-report-build`**
>
> This skill is **deprecated** in favor of `html-report-build`, which removes the entire Python middleware chain (markdown lib / jinja2 / bs4 / Pillow) and uses LLM-direct HTML generation following the Claude Artifacts pattern (Anthropic's HTML-first default since 2026-05).
>
> **Why deprecated**:
> 1. Python middleware is no longer needed — LLM can directly read `.md` and produce semantically-restructured HTML, not just transform markdown syntax.
> 2. Actual adoption count was **0** (`.preview.html` files never produced; only 1 manual self_phasing topic folder).
> 3. User feedback (2026-05-13): "直接產生適合釐清邏輯思路的 HTML 比較好" — LLM-direct beats Python transformation.
>
> **Migration**:
> - Old: `python3 .claude/skills/html-preview/tools/dispatch.py <md>`
> - New: Invoke `html-report-build` skill (manual, auto-trigger, or via `Skill` tool)
>
> **Sunset timeline**:
> - **Stage 1 (now)**: Banner added; skill remains callable but discouraged.
> - **Stage 2 (Tier A skill rewrite)**: 6 Tier A skills (weekly-report / structured-tech-report / results-report / feature-layered-observation / methodology-audit / conclude-research) switch Stop hook to `html-report-build`.
> - **Stage 3 (2 weeks, if no regression)**: Delete `tools/*.py`; keep this SKILL.md as stub.
>
> **References**:
> - New skill: `InterSubMod/.claude/skills/html-report-build/SKILL.md`
> - Plan: `/bip7_disk/liaoyoyo2001/.claude/plans/frolicking-tinkering-hopcroft.md`

---

# html-preview (legacy reference content below)

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
