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
