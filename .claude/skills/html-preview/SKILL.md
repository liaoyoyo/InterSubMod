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
