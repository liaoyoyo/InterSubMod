---
name: image-gen
description: |
  生成示意圖 / 流程圖 / 概念示意 / icon 圖片。雙軌：concept_diagram & data_mockup → AI 軌 (codex CLI $imagegen via ChatGPT OAuth, 免 USD)；flow_diagram & icon → 程式軌 (pycairo, 零成本、不錯字、deterministic)。輸出 PNG 到指定 figures/ 目錄。
  USE WHEN：「需要示意圖」「畫個流程圖」「補個 icon」「LOH 概念圖怎麼示意」「報告補圖」「生成圖片」「生圖」「畫示意」「畫個概念圖」「補張圖」「需要圖」「插圖」、PI 寫報告時提到圖片需求、weekly-report / structured-tech-report / results-report 結束時偵測 figure-needed 標記。
  SKIP WHEN：matplotlib/seaborn 真實數據圖（用 InterSubMod/scripts/lib/plot_setup.py）、UI 截圖（用瀏覽器開發工具）、純文字流程（用 mermaid markdown）、PPTX 內已 review 過的圖。
---

# image-gen

## Phase & Chain Position

- Stand-alone skill, callable directly or via myPPT routing.
- Upstream: `weekly-report` / `structured-tech-report` / `results-report` etc. emit `<!-- figure-needed: type, slug=foo -->` markers in `.md`.
- Downstream: `image-vision-check` reads produced PNGs and audits.

## Dependencies

**Uses (runtime)**:
- `codex` CLI v0.125+ (must be `Logged in using ChatGPT` via `codex login`).
- `python3` ≥ 3.10 with `pycairo`, `Pillow`, `pyyaml`.

**Used by**:
- `image-vision-check` (consumes PNG output)
- `html-preview` (Phase 2; consumes PNG via PIL base64 encoding)
- `myPPT` (orchestration)

**Reads**:
- `<prompts_dir>/*.yaml` (each conforms to `templates/{type}.yaml` schema)

**Writes**:
- `<figures_dir>/<basename>.png` for each `<basename>.yaml`

## Failure Modes & Diagnostics

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `[FAIL] codex login` | OAuth session expired | run `codex login` |
| `[FAIL] pycairo missing` | system pycairo not installed | `pip install --user pycairo` or `apt install python3-cairo` |
| AI track produces blank/empty PNG | ChatGPT quota exhausted | wait until next month OR enable deprecated API key path; cairo track unaffected |
| cairo track raises `ValueError: backend mismatch` | YAML `backend: ai` but routed to cairo | check YAML `type` field (only `flow_diagram` & `icon` are cairo) |
| Output PNG not found in target dir | codex wrote to `~/.codex/generated_images/` | `codex_output_collector.py` should rescue; if not, check `since_epoch` filter |

## Quick Usage

```bash
# Dry-run (preview only, no codex call, no cairo render)
python3 InterSubMod/.claude/skills/image-gen/tools/dispatch.py \
  <prompts_dir> <figures_dir> --dry-run

# Full run
python3 InterSubMod/.claude/skills/image-gen/tools/dispatch.py \
  <prompts_dir> <figures_dir>

# Re-render even if target exists
python3 InterSubMod/.claude/skills/image-gen/tools/dispatch.py \
  <prompts_dir> <figures_dir> --force

# Resize all outputs to standard size after gen
python3 InterSubMod/.claude/skills/image-gen/tools/dispatch.py \
  <prompts_dir> <figures_dir> --postprocess-size 1024x1024
```

## Pre-flight Required

Before any run:
```bash
bash InterSubMod/.claude/skills/image-gen/tools/preflight.sh
```

## Cost Model

| Backend | Per-image cost | Quota source |
|---------|----------------|--------------|
| AI (codex `$imagegen`) | ~30k agent tokens (low quality), 3-5× consumption rate | ChatGPT subscription monthly |
| cairo (program render) | 0 | local CPU only |

## See Also

- Spec: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md`
- Playbook: `InterSubMod/.claude/skills/image-gen/playbook.md`
- Companion skill: `InterSubMod/.claude/skills/image-vision-check/SKILL.md`
