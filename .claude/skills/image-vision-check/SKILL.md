---
name: image-vision-check
description: |
  Audit AI/cairo 產出的圖片品質，依 6 維度 checklist 評分（prompt fidelity / 可讀性 / 重點突出 / 設計禁忌避免 / 列印友善 / 跨圖風格一致），輸出 quality.json + verdict (pass / partial / fail) + suggestions。閥值 4/6 通過。fail 時可觸發 image-gen 自動重生 1 次。
  USE WHEN：「圖夠不夠好」「圖達標了嗎」「檢核圖片」「這張圖能用嗎」「vision check」「品檢」「review 圖」、image-gen 產出後自動觸發。
  SKIP WHEN：真實數據圖 (用人眼確認 + statistical test)、PPTX 內已 review 過的圖、純色 placeholder。
---

# image-vision-check

## Phase & Chain Position

- Stand-alone skill, callable directly or via myPPT routing.
- Upstream: `image-gen` (consumes its PNG output).
- Downstream: `html-preview` (Phase 2; uses verdict to decide which images to inline).

## Dependencies

**Uses (runtime)**:
- Claude `Read` tool (built-in vision) — reads PNG and scores 6 dimensions.
- `python3` ≥ 3.10 with `pyyaml`, `jsonschema` (optional, for strict schema check).

**Used by**:
- `image-gen` (auto-retry loop on fail)
- `html-preview` (Phase 2)
- `myPPT` (orchestration)

**Reads**:
- `<figures_dir>/*.png`
- (paired) `<prompts_dir>/*.yaml` for context (subject, constraints) — optional but improves Claude's judgment.

**Writes**:
- `<figures_dir>/quality.json` (one entry per image)

## Failure Modes & Diagnostics

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Claude `Read` fails on PNG | corrupt PNG | re-run image-gen on that prompt |
| All scores 6/6 across all images | over-permissive Claude | review checklist criteria; tighten language |
| All scores ≤2/6 | under-permissive Claude OR genuinely bad images | inspect a few PNGs visually |
| `quality.json` schema invalid | runner code drift | run `tools/vision_check_runner.py --validate-schema-only` |

## Quick Usage

```bash
# Check single image
python3 InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py \
  <figures_dir>/<basename>.png

# Batch check whole dir
python3 InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py \
  <figures_dir>

# With paired prompts for context
python3 InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py \
  <figures_dir> --prompts <prompts_dir>
```

## Verdict Rules (D6 + D8)

| Score | Verdict | Action |
|-------|---------|--------|
| 6/6 | `pass` | OK, proceed |
| 4-5/6 | `partial` | Auto-retry 1× via image-gen; if still partial, surface to PI with suggestions |
| ≤3/6 | `fail` | Auto-retry 1× via image-gen; if still fail, halt + surface to PI |

## See Also

- Spec: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md`
- Companion: `InterSubMod/.claude/skills/image-gen/SKILL.md`
- Quality JSON schema: `InterSubMod/.claude/skills/image-vision-check/schemas/quality.schema.json`
