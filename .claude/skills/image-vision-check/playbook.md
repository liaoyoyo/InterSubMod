# image-vision-check Playbook

## Step-by-step workflow

### Step 1: Inputs
Caller provides:
- `figures_dir` (required) — directory of PNGs to check.
- `prompts_dir` (optional) — paired YAMLs for context. If absent, runner falls back to `figures_dir/../prompts/`.

### Step 2: Per-image audit
For each PNG in `figures_dir`:
1. Find paired YAML at `<prompts_dir>/<basename>.yaml`. If missing, log warning and use generic checklist.
2. Determine `type` from YAML.
3. Load `checklists/<type>_check.md`.
4. Use Claude `Read` tool on PNG. Send to Claude with `prompts/vision_review_master.md` + checklist content + YAML context.
5. Claude returns JSON. Validate against `schemas/quality.schema.json`.
6. Append entry to `<figures_dir>/quality.json`.

### Step 3: Verdict aggregation
After all images checked, compute:
- `n_pass` / `n_partial` / `n_fail`
- Overall verdict: `all_pass` if all pass; `mixed` otherwise.

### Step 4: Auto-retry coordination (D8)
If any verdict ∈ {`partial`, `fail`}:
- For each fail/partial image, locate paired YAML.
- Apply suggestions to YAML (string-level edit; mark `version: 2` in retry history).
- Re-run image-gen for those YAMLs only.
- Re-check.
- If still fail → set verdict `needs_human` and write to `quality.json`. Do not retry again (D8: max 1 retry).

### Step 5: Hand off
- If `all_pass`, suggest `/html-report-build` or report DONE.
- If `needs_human`, surface checklist + suggestions to PI for manual prompt tuning.

## quality.json structure (cumulative)

```json
{
  "checked_at": "2026-05-10T14:30:00Z",
  "figures_dir": "InterSubMod/examples/phase1_demo/demo_topic/figures",
  "summary": { "n_pass": 3, "n_partial": 1, "n_fail": 0, "overall": "mixed" },
  "entries": [
    { /* per-image entry, see schema */ }
  ]
}
```

## Cost note

Each Claude `Read` call ~ 500 tokens prompt + ~ 200 tokens output ≈ ~$0.005/image at Opus 4.7 pricing.
