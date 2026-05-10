# image-gen Playbook

## Step-by-step workflow

### Step 1: Preflight
Run `tools/preflight.sh`. Halt on any `[FAIL]`. Provide install instruction to user.

### Step 2: Plan
Call `dispatch.plan_run(prompts_dir, output_dir)` to enumerate prompts and classify by backend.
Show cost preview to user via `dispatch.cost_preview(plan)`.
If user types `n` or hesitates, exit.

### Step 3: Dispatch each prompt
For each prompt in plan:
- If target PNG exists and `--force` not set, skip.
- Else `route_one(prompt_path, target_path)`.

### Step 4: Post-process (optional)
If `--postprocess-size <WxH>` given, run PIL resize for all outputs.

### Step 5: Hand off to image-vision-check
After all prompts complete, suggest:
> `[image-gen] All N images generated. Run vision check?`
> `→ /image-vision-check <figures_dir>`

## Prompt YAML Authoring Tips

- One YAML per image. Filename `{slug}.yaml` → output `{slug}.png`.
- Pick `backend` based on `type`:
  - `concept_diagram` → `ai`
  - `data_mockup` → `ai`
  - `flow_diagram` → `cairo`
  - `icon` → `cairo`
- For AI prompts, write `subject` as a single declarative sentence ("LOH region triggers haplotype phasing degradation"). Avoid imperative ("draw a chromosome"). Prefer scientific accuracy in subject text.
- For cairo prompts, define `nodes` and `edges` (flow) or `shape` + `glyph` (icon) precisely; cairo will not infer.

## Anthropic Design Taboos (D14 enforcement)

In every YAML's `constraints` field, include at minimum:
- `no gradient overuse`
- `no glass morphism`
- `no multi-indigo stacking`
- `colorblind-friendly`

These propagate into the AI prompt text via `extract_prompt_text()`.

## Known Edge Cases

- **codex CLI prompts for approval**: `--ask-for-approval never` should suppress; if not, run interactively first to OAuth refresh.
- **Cairo glyph clipped**: increase `padding_px` in icon YAML, re-run.
- **AI generates wrong text**: AI tracks have inherent randomness for text; if labels critical, switch to `flow_diagram` (cairo).
