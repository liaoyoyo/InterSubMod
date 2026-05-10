---
title: "Phase 1 Implementation Plan — image-gen + image-vision-check + Demo Topic Folder"
date: 2026-05-10
status: ready_to_execute
version: 0.1.0
authors: [liaoyoyo2001 (PI), Claude Opus 4.7 (assistant)]
spec: InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md
phase: 1 of 3
estimated_effort: 2-3 工作天 (14 tasks × 15-30 min/task)
tags: [implementation_plan, phase1, image_gen, vision_check]
---

# Phase 1 Implementation Plan: image-gen + image-vision-check + Demo Topic Folder

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `image-gen` + `image-vision-check` Claude Code skills under `InterSubMod/.claude/skills/` and validate end-to-end closed loop (`prompt → AI/cairo → PNG → vision check pass → demo HTML inline display`) via a hand-written demo topic folder.

**Architecture:** Two new skills.
- `image-gen` routes by prompt YAML `type` field — `concept_diagram` & `data_mockup` go to **AI track** (codex CLI `$imagegen` via ChatGPT OAuth, no API key); `flow_diagram` & `icon` go to **program track** (pycairo). PIL post-processes both.
- `image-vision-check` uses Claude `Read` tool to score 6 dimensions per image, writes `quality.json`, auto-retries 1× on fail.
- End-to-end demo lives at `InterSubMod/examples/phase1_demo/demo_topic/` with hand-written `index.html` proving base64 inline + topic-folder structure work, validating the design before Phase 2 builds the full `html-preview` skill.

**Tech Stack:** codex CLI v0.125+ (ChatGPT OAuth, already logged in) · Python 3.10+ + `pycairo` + `Pillow` + `pyyaml` · Claude `Read` tool for vision scoring · Bash glue scripts.

**Reference design doc:** `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md` (decisions D1–D14).

---

## Pre-Flight Check (run once before Task 1)

These MUST pass before starting Task 1. If any fail, install the missing dependency and re-verify.

```bash
# 1. codex CLI logged in via ChatGPT OAuth
codex login status 2>&1 | grep -q "Logged in" && echo "[OK] codex login" || echo "[FAIL] run: codex login"

# 2. pycairo available
python3 -c "import cairo; print(f'[OK] pycairo {cairo.version}')" 2>&1 \
  || echo "[FAIL] install: pip install --user pycairo  OR  apt install python3-cairo"

# 3. Pillow available
python3 -c "from PIL import Image, __version__; print(f'[OK] Pillow {__version__}')" 2>&1 \
  || echo "[FAIL] install: pip install --user Pillow"

# 4. pyyaml available
python3 -c "import yaml; print(f'[OK] pyyaml {yaml.__version__}')" 2>&1 \
  || echo "[FAIL] install: pip install --user pyyaml"
```

Expected: 4 lines all `[OK]`. Failing dependencies block ALL tasks.

---

## File Structure

```
InterSubMod/
├── .claude/skills/
│   ├── image-gen/                                    # NEW skill
│   │   ├── SKILL.md
│   │   ├── playbook.md
│   │   ├── templates/
│   │   │   ├── concept_diagram.yaml
│   │   │   ├── flow_diagram.yaml
│   │   │   ├── data_mockup.yaml
│   │   │   └── icon.yaml
│   │   ├── tools/
│   │   │   ├── dispatch.py                           # main entry
│   │   │   ├── codex_image_gen.py                    # AI track
│   │   │   ├── codex_output_collector.py             # 搬 PNG
│   │   │   ├── cairo_render.py                       # 程式 track
│   │   │   ├── pil_postprocess.py                    # resize/crop
│   │   │   └── preflight.sh                          # dep check
│   │   ├── prompts/
│   │   │   └── prompt_quality_self_review.md
│   │   └── tests/
│   │       ├── test_dispatch.py
│   │       ├── test_cairo_render.py
│   │       ├── test_pil_postprocess.py
│   │       ├── test_codex_output_collector.py
│   │       └── fixtures/
│   │           ├── valid_concept.yaml
│   │           ├── valid_flow.yaml
│   │           ├── valid_data.yaml
│   │           └── valid_icon.yaml
│   │
│   └── image-vision-check/                           # NEW skill
│       ├── SKILL.md
│       ├── playbook.md
│       ├── checklists/
│       │   ├── concept_diagram_check.md
│       │   ├── flow_diagram_check.md
│       │   ├── data_mockup_check.md
│       │   └── icon_check.md
│       ├── tools/
│       │   └── vision_check_runner.py
│       ├── schemas/
│       │   └── quality.schema.json
│       ├── prompts/
│       │   └── vision_review_master.md
│       └── tests/
│           ├── test_vision_check_runner.py
│           └── fixtures/
│               └── sample_quality.json
│
├── examples/phase1_demo/demo_topic/                  # NEW (hand-written demo)
│   ├── README.md
│   ├── demo_topic.md
│   ├── index.html
│   ├── prompts/
│   │   ├── fig1_concept_loh.yaml
│   │   ├── fig2_flow_pipeline.yaml
│   │   ├── fig3_data_mockup.yaml
│   │   └── fig4_icon_status.yaml
│   └── figures/                                      # .gitignore
│       └── .gitkeep
│
└── .gitignore                                        # MODIFY: append patterns
```

---

### Task 1: Project scaffolding (directories + .gitignore)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/{templates,tools,prompts,tests/fixtures}/.gitkeep`
- Create: `InterSubMod/.claude/skills/image-vision-check/{checklists,tools,schemas,prompts,tests/fixtures}/.gitkeep`
- Create: `InterSubMod/examples/phase1_demo/demo_topic/{prompts,figures}/.gitkeep`
- Modify: `InterSubMod/.gitignore`

- [ ] **Step 1: Create skill + demo directories**

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod

mkdir -p .claude/skills/image-gen/{templates,tools,prompts,tests/fixtures}
mkdir -p .claude/skills/image-vision-check/{checklists,tools,schemas,prompts,tests/fixtures}
mkdir -p examples/phase1_demo/demo_topic/{prompts,figures}

# .gitkeep so empty dirs survive git
touch .claude/skills/image-gen/{templates,tools,prompts,tests/fixtures}/.gitkeep
touch .claude/skills/image-vision-check/{checklists,tools,schemas,prompts,tests/fixtures}/.gitkeep
touch examples/phase1_demo/demo_topic/{prompts,figures}/.gitkeep
```

- [ ] **Step 2: Verify directories created**

Run:
```bash
find .claude/skills/image-gen .claude/skills/image-vision-check examples/phase1_demo -type d | sort
```
Expected: 13 directory paths.

- [ ] **Step 3: Append .gitignore rules**

Append to `InterSubMod/.gitignore`:

```gitignore

# === image-gen / html-preview (added 2026-05-10, design doc D7) ===
# Generated PNGs not committed; prompts/ retained for reproducibility
**/figures/*.png
!**/figures/.gitkeep

# Codex CLI working dir (image gen output staging)
.codex_staging/
```

- [ ] **Step 4: Verify .gitignore works**

```bash
echo "test" > examples/phase1_demo/demo_topic/figures/test.png
git status --porcelain examples/phase1_demo/demo_topic/figures/
```
Expected: empty output (test.png ignored). Then `rm examples/phase1_demo/demo_topic/figures/test.png`.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/image-gen/ .claude/skills/image-vision-check/ examples/phase1_demo/ .gitignore
git commit -m "feat(skills): scaffolding for image-gen + image-vision-check + demo topic folder

Phase 1 of HTML preview + image-gen 3-skill design.
Spec: InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 2: Pre-flight check script (preflight.sh)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/tools/preflight.sh`
- Create: `InterSubMod/.claude/skills/image-gen/tests/test_preflight.sh`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-gen/tests/test_preflight.sh`:

```bash
#!/usr/bin/env bash
# Test that preflight.sh exits 0 when all deps present, exits non-zero otherwise.
set -euo pipefail
SCRIPT="$(dirname "$0")/../tools/preflight.sh"

# Test 1: script exists and executable
[ -x "$SCRIPT" ] || { echo "FAIL: $SCRIPT not executable"; exit 1; }

# Test 2: returns 0 in current env (assumes deps installed per Pre-Flight Check)
if "$SCRIPT" >/dev/null 2>&1; then
    echo "PASS: preflight.sh returns 0 in healthy env"
else
    echo "FAIL: preflight.sh returned non-zero in env that should be healthy"
    "$SCRIPT" 2>&1 | head -10
    exit 1
fi

echo "all preflight tests passed"
```

```bash
chmod +x .claude/skills/image-gen/tests/test_preflight.sh
```

- [ ] **Step 2: Run test to verify it fails**

```bash
bash .claude/skills/image-gen/tests/test_preflight.sh
```
Expected: `FAIL: ... preflight.sh not executable` (because preflight.sh doesn't exist yet).

- [ ] **Step 3: Write minimal implementation**

Create `InterSubMod/.claude/skills/image-gen/tools/preflight.sh`:

```bash
#!/usr/bin/env bash
# preflight.sh — verify image-gen skill dependencies before run.
# Exits 0 if all OK, non-zero with error message if any dep missing.
set -euo pipefail

errors=0

# 1. codex CLI logged in via ChatGPT OAuth
if ! codex login status 2>&1 | grep -q "Logged in"; then
    echo "[FAIL] codex CLI not logged in. Run: codex login" >&2
    errors=$((errors+1))
else
    echo "[OK] codex login"
fi

# 2. pycairo
if ! python3 -c "import cairo" 2>/dev/null; then
    echo "[FAIL] pycairo missing. Install: pip install --user pycairo" >&2
    errors=$((errors+1))
else
    echo "[OK] pycairo $(python3 -c 'import cairo; print(cairo.version)')"
fi

# 3. Pillow
if ! python3 -c "from PIL import Image" 2>/dev/null; then
    echo "[FAIL] Pillow missing. Install: pip install --user Pillow" >&2
    errors=$((errors+1))
else
    echo "[OK] Pillow $(python3 -c 'import PIL; print(PIL.__version__)')"
fi

# 4. pyyaml
if ! python3 -c "import yaml" 2>/dev/null; then
    echo "[FAIL] pyyaml missing. Install: pip install --user pyyaml" >&2
    errors=$((errors+1))
else
    echo "[OK] pyyaml $(python3 -c 'import yaml; print(yaml.__version__)')"
fi

if [ "$errors" -gt 0 ]; then
    echo "" >&2
    echo "$errors dependency check(s) failed. See messages above." >&2
    exit 1
fi
```

```bash
chmod +x .claude/skills/image-gen/tools/preflight.sh
```

- [ ] **Step 4: Run test to verify it passes**

```bash
bash .claude/skills/image-gen/tests/test_preflight.sh
```
Expected: `PASS: preflight.sh returns 0 in healthy env` and `all preflight tests passed`.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/image-gen/tools/preflight.sh .claude/skills/image-gen/tests/test_preflight.sh
git commit -m "feat(image-gen): preflight.sh dependency check + test

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: Prompt YAML schema + 4 templates (concept/flow/data/icon)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/templates/{concept_diagram,flow_diagram,data_mockup,icon}.yaml`
- Create: `InterSubMod/.claude/skills/image-gen/tests/fixtures/valid_{concept,flow,data,icon}.yaml`
- Create: `InterSubMod/.claude/skills/image-gen/tests/test_yaml_schema.py`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-gen/tests/test_yaml_schema.py`:

```python
"""Validate prompt YAML files conform to expected schema."""
import sys
from pathlib import Path
import yaml

SKILL_DIR = Path(__file__).parent.parent
TEMPLATES = SKILL_DIR / "templates"
FIXTURES = SKILL_DIR / "tests" / "fixtures"

REQUIRED_FIELDS = {"type", "subject", "constraints", "output"}
VALID_TYPES = {"concept_diagram", "flow_diagram", "data_mockup", "icon"}

def validate(path: Path) -> list[str]:
    errors = []
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        return [f"{path}: top-level not a dict"]
    missing = REQUIRED_FIELDS - data.keys()
    if missing:
        errors.append(f"{path}: missing fields {missing}")
    if data.get("type") not in VALID_TYPES:
        errors.append(f"{path}: type='{data.get('type')}' not in {VALID_TYPES}")
    out = data.get("output", {})
    if "size" not in out or "quality" not in out:
        errors.append(f"{path}: output.{{size,quality}} missing")
    return errors

def main():
    files = list(TEMPLATES.glob("*.yaml")) + list(FIXTURES.glob("valid_*.yaml"))
    if not files:
        print(f"FAIL: no yaml found under {TEMPLATES} or {FIXTURES}")
        return 1
    all_errors = []
    for f in files:
        all_errors.extend(validate(f))
    if all_errors:
        for e in all_errors:
            print(f"FAIL: {e}")
        return 1
    print(f"PASS: {len(files)} yaml files valid")
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-gen/tests/test_yaml_schema.py
```
Expected: `FAIL: no yaml found ...` (templates not yet created).

- [ ] **Step 3: Create 4 template YAMLs**

`InterSubMod/.claude/skills/image-gen/templates/concept_diagram.yaml`:

```yaml
# Template: concept_diagram (AI track via codex $imagegen)
# Use for: Inter-Subclonal / LOH / phasing 概念示意圖
type: concept_diagram
subject: "<one-sentence main subject, e.g., LOH region triggers haplotype phasing degradation>"
key_elements:
  - "<element 1, e.g., chromosome ideogram with shaded LOH region>"
  - "<element 2, e.g., paired reads above and below>"
  - "<element 3>"
labels:
  top: "<top-of-image label>"
  bottom: "<bottom label>"
style:
  palette: ["#1f2937", "#ef4444"]   # 2-3 colors, no gradient overuse
  background: "white"
  composition: "centered"
constraints:
  - "no gradient overuse"
  - "no glass morphism"
  - "no multi-indigo stacking"
  - "colorblind-friendly (no red/green pair)"
  - "print-friendly contrast >= 4.5:1"
  - "system-ui font feel"
output:
  size: "1024x1024"
  quality: "high"
backend: ai      # always 'ai' for concept_diagram
```

`InterSubMod/.claude/skills/image-gen/templates/flow_diagram.yaml`:

```yaml
# Template: flow_diagram (program track via pycairo)
# Use for: G1→G7 pipeline、step 順序、條件分支
type: flow_diagram
subject: "<flow name, e.g., Self-Phasing V5 Pipeline G1-G7>"
nodes:
  - id: "g1"
    label: "G1: ReadParser"
    x: 100
    y: 100
  - id: "g2"
    label: "G2: VCF Loader"
    x: 350
    y: 100
edges:
  - from: "g1"
    to: "g2"
    label: ""           # optional edge label
style:
  palette: ["#1f2937", "#3b82f6"]
  background: "white"
  font_family: "DejaVu Sans"
  font_size: 14
constraints:
  - "deterministic layout (no AI randomness)"
  - "exact text rendering"
output:
  size: "1024x768"
  quality: "high"       # cairo always renders sharp; quality field for consistency
backend: cairo          # always 'cairo' for flow_diagram
```

`InterSubMod/.claude/skills/image-gen/templates/data_mockup.yaml`:

```yaml
# Template: data_mockup (AI track) — 概念型 bar/line 示意，非真實數據
# Real data plots → use scripts/lib/plot_setup.py (matplotlib)
type: data_mockup
subject: "<chart concept, e.g., bar chart showing 7 samples' AUC distribution mock-up>"
chart_type: "bar"             # bar | line | scatter
axes:
  x_label: "<x label>"
  y_label: "<y label>"
data_hint: "<verbal description of trend, e.g., HCC1395 highest, COLO829 lowest, gradient between>"
style:
  palette: ["#1f2937", "#10b981"]
  background: "white"
constraints:
  - "axes clearly labeled"
  - "no gradient overuse"
  - "data is conceptual mock, not real"
output:
  size: "1024x1024"
  quality: "high"
backend: ai
```

`InterSubMod/.claude/skills/image-gen/templates/icon.yaml`:

```yaml
# Template: icon (program track via pycairo)
# Use for: 條目訂點、序號標記、status 標籤
type: icon
subject: "<icon meaning, e.g., status_pass>"
shape: "circle"               # circle | square | rounded_square | hexagon
glyph: "✓"                    # single CJK / symbol / digit
glyph_color: "#ffffff"
fill_color: "#10b981"
style:
  size_px: 256                # 256 / 512
  padding_px: 24
  background: "transparent"
constraints:
  - "single symbol / digit only"
  - "high contrast glyph vs fill"
output:
  size: "256x256"
  quality: "n/a"
backend: cairo
```

- [ ] **Step 4: Create fixture YAMLs (valid examples for tests)**

Copy each template to a fixture with one realistic example filled in. For example:

`InterSubMod/.claude/skills/image-gen/tests/fixtures/valid_concept.yaml`:

```yaml
type: concept_diagram
subject: "LOH region triggers haplotype phasing degradation in tumor sample"
key_elements:
  - "chromosome 8 ideogram with shaded LOH band at q24"
  - "paired short reads aligned above and below the chromosome"
  - "an arrow showing phasing failure within the LOH region"
labels:
  top: "Tumor sample with LOH"
  bottom: "Phasing degrades inside LOH (shaded)"
style:
  palette: ["#1f2937", "#ef4444"]
  background: "white"
  composition: "centered"
constraints:
  - "no gradient overuse"
  - "no glass morphism"
  - "colorblind-friendly"
  - "print-friendly contrast >= 4.5:1"
  - "system-ui font feel"
output:
  size: "1024x1024"
  quality: "high"
backend: ai
```

`InterSubMod/.claude/skills/image-gen/tests/fixtures/valid_flow.yaml`:

```yaml
type: flow_diagram
subject: "Self-Phasing V5 minimal 3-step demo"
nodes:
  - id: "n1"
    label: "ReadParser"
    x: 80
    y: 200
  - id: "n2"
    label: "Phasing"
    x: 360
    y: 200
  - id: "n3"
    label: "Output VCF"
    x: 640
    y: 200
edges:
  - from: "n1"
    to: "n2"
    label: ""
  - from: "n2"
    to: "n3"
    label: ""
style:
  palette: ["#1f2937", "#3b82f6"]
  background: "white"
  font_family: "DejaVu Sans"
  font_size: 14
constraints:
  - "deterministic layout"
  - "exact text rendering"
output:
  size: "1024x400"
  quality: "high"
backend: cairo
```

`InterSubMod/.claude/skills/image-gen/tests/fixtures/valid_data.yaml`:

```yaml
type: data_mockup
subject: "Bar chart mock-up showing 7 samples' baseline AUC distribution"
chart_type: "bar"
axes:
  x_label: "Sample"
  y_label: "AUC"
data_hint: "HCC1395 ~0.78 (highest), COLO829 ~0.55 (lowest), other 5 samples between 0.60-0.72"
style:
  palette: ["#1f2937", "#10b981"]
  background: "white"
constraints:
  - "axes clearly labeled"
  - "no gradient overuse"
  - "data is conceptual mock, not real"
output:
  size: "1024x1024"
  quality: "high"
backend: ai
```

`InterSubMod/.claude/skills/image-gen/tests/fixtures/valid_icon.yaml`:

```yaml
type: icon
subject: "status_pass"
shape: "circle"
glyph: "✓"
glyph_color: "#ffffff"
fill_color: "#10b981"
style:
  size_px: 256
  padding_px: 24
  background: "transparent"
constraints:
  - "single symbol"
  - "high contrast"
output:
  size: "256x256"
  quality: "n/a"
backend: cairo
```

- [ ] **Step 5: Run test to verify it passes**

```bash
python3 .claude/skills/image-gen/tests/test_yaml_schema.py
```
Expected: `PASS: 8 yaml files valid` (4 templates + 4 fixtures).

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/image-gen/templates/ .claude/skills/image-gen/tests/fixtures/ .claude/skills/image-gen/tests/test_yaml_schema.py
git commit -m "feat(image-gen): YAML schema + 4 templates (concept/flow/data/icon) + fixtures + test

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: codex_image_gen.py (AI track wrapper)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/tools/codex_image_gen.py`
- Create: `InterSubMod/.claude/skills/image-gen/tests/test_codex_image_gen.py`

- [ ] **Step 1: Write the failing test (mock codex call)**

Create `InterSubMod/.claude/skills/image-gen/tests/test_codex_image_gen.py`:

```python
"""Test codex_image_gen.build_command builds correct codex CLI invocation."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from codex_image_gen import build_command, extract_prompt_text  # noqa

FIXTURES = Path(__file__).parent / "fixtures"

def test_extract_concept_prompt():
    yaml_path = FIXTURES / "valid_concept.yaml"
    text = extract_prompt_text(yaml_path)
    assert "LOH region" in text, f"Expected LOH in prompt, got: {text!r}"
    assert "$imagegen" in text, f"Expected $imagegen suffix, got: {text!r}"
    assert "no gradient" in text.lower(), "Expected constraints embedded"

def test_build_command():
    yaml_path = FIXTURES / "valid_concept.yaml"
    out_dir = Path("/tmp/test_codex_out")
    cmd = build_command(yaml_path, out_dir)
    assert cmd[0] == "codex", f"First token must be 'codex', got {cmd[0]}"
    assert "exec" in cmd, "Must use 'exec' subcommand"
    assert "--image-dir" in cmd, "Must specify --image-dir"
    idx = cmd.index("--image-dir")
    assert cmd[idx+1] == str(out_dir), f"--image-dir value wrong: {cmd[idx+1]}"
    last = cmd[-1]
    assert "$imagegen" in last, f"Prompt argument must end with $imagegen: {last!r}"

def main():
    tests = [test_extract_concept_prompt, test_build_command]
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

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-gen/tests/test_codex_image_gen.py
```
Expected: `ImportError: cannot import name 'build_command'` or similar.

- [ ] **Step 3: Write minimal implementation**

Create `InterSubMod/.claude/skills/image-gen/tools/codex_image_gen.py`:

```python
"""Wrap codex CLI exec for AI-track image generation.

Reads prompt YAML, constructs codex exec command with $imagegen, runs it,
expects output PNG to land in --image-dir.

Backend: codex CLI v0.125+ (ChatGPT OAuth).
Model: gpt-image-2 (codex default for $imagegen as of 2026-04-21).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import List

import yaml


def extract_prompt_text(yaml_path: Path) -> str:
    """Read prompt YAML and render to natural-language string for codex.

    Format: subject + key_elements (if any) + labels + style hints + constraints.
    Always appends '$imagegen' to trigger codex CLI's image gen skill.
    """
    data = yaml.safe_load(yaml_path.read_text())
    parts: list[str] = []

    parts.append(data.get("subject", ""))

    elems = data.get("key_elements") or []
    if elems:
        parts.append("Include: " + "; ".join(elems) + ".")

    labels = data.get("labels") or {}
    for pos, lab in labels.items():
        if lab:
            parts.append(f"Label at {pos}: '{lab}'.")

    chart = data.get("chart_type")
    if chart:
        ax = data.get("axes", {})
        parts.append(f"Chart type: {chart}. X axis: {ax.get('x_label','')}, Y axis: {ax.get('y_label','')}.")
        parts.append(f"Data trend hint: {data.get('data_hint','')}")

    style = data.get("style", {})
    pal = style.get("palette") or []
    if pal:
        parts.append(f"Palette: {', '.join(pal)}. Background: {style.get('background','white')}.")

    cons = data.get("constraints") or []
    if cons:
        parts.append("Constraints: " + "; ".join(cons) + ".")

    out = data.get("output", {})
    parts.append(f"Output size: {out.get('size','1024x1024')}, quality: {out.get('quality','high')}.")

    parts.append("$imagegen")
    return " ".join(p for p in parts if p)


def build_command(yaml_path: Path, output_dir: Path) -> List[str]:
    """Construct codex exec CLI command list."""
    prompt_text = extract_prompt_text(yaml_path)
    return [
        "codex", "exec",
        "--image-dir", str(output_dir),
        "--ask-for-approval", "never",
        prompt_text,
    ]


def run(yaml_path: Path, output_dir: Path, *, dry_run: bool = False) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = build_command(yaml_path, output_dir)
    if dry_run:
        print("DRY RUN:", " ".join(repr(c) for c in cmd))
        return 0
    print(f"[codex_image_gen] Running for {yaml_path.name} → {output_dir}/")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print("STDOUT:", proc.stdout, file=sys.stderr)
        print("STDERR:", proc.stderr, file=sys.stderr)
        print(f"[codex_image_gen] codex exec exited {proc.returncode}", file=sys.stderr)
    return proc.returncode


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print("Usage: codex_image_gen.py <prompt.yaml> <output_dir> [--dry-run]")
        return 2
    yaml_path = Path(argv[1])
    out_dir = Path(argv[2])
    dry = "--dry-run" in argv[3:]
    return run(yaml_path, out_dir, dry_run=dry)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 .claude/skills/image-gen/tests/test_codex_image_gen.py
```
Expected: `PASS: test_extract_concept_prompt` and `PASS: test_build_command`.

- [ ] **Step 5: Smoke test dry-run against real codex (no network call)**

```bash
python3 .claude/skills/image-gen/tools/codex_image_gen.py \
  .claude/skills/image-gen/tests/fixtures/valid_concept.yaml \
  /tmp/test_codex_out --dry-run
```
Expected: `DRY RUN: 'codex' 'exec' '--image-dir' '/tmp/test_codex_out' '--ask-for-approval' 'never' '...$imagegen'`

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/image-gen/tools/codex_image_gen.py .claude/skills/image-gen/tests/test_codex_image_gen.py
git commit -m "feat(image-gen): codex_image_gen.py — AI track wrapper for codex CLI \$imagegen

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 5: codex_output_collector.py (rescue PNGs from codex output)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/tools/codex_output_collector.py`
- Create: `InterSubMod/.claude/skills/image-gen/tests/test_codex_output_collector.py`

Background: codex `--image-dir` should write to specified dir, but if codex falls back to `~/.codex/generated_images/`, this collector finds the latest PNG and moves it.

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-gen/tests/test_codex_output_collector.py`:

```python
"""Test that collector finds latest PNG and moves to target."""
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from codex_output_collector import collect_latest  # noqa


def test_finds_and_moves_latest_png():
    with tempfile.TemporaryDirectory() as src_str, tempfile.TemporaryDirectory() as dst_str:
        src = Path(src_str)
        dst = Path(dst_str)
        # Create 2 PNGs with different mtimes
        old = src / "old.png"
        old.write_bytes(b"\x89PNG\r\n\x1a\n" + b"old")
        time.sleep(0.05)
        new = src / "new.png"
        new.write_bytes(b"\x89PNG\r\n\x1a\n" + b"new")

        target_name = "moved_fig1.png"
        result = collect_latest([src], dst, target_name=target_name, since_epoch=0.0)

        assert result is not None, "collect_latest returned None"
        assert result.name == target_name, f"Wrong name: {result.name}"
        assert (dst / target_name).read_bytes().endswith(b"new"), "Wrong PNG content"


def test_returns_none_when_no_pngs():
    with tempfile.TemporaryDirectory() as src_str, tempfile.TemporaryDirectory() as dst_str:
        src = Path(src_str)
        dst = Path(dst_str)
        result = collect_latest([src], dst, target_name="x.png", since_epoch=0.0)
        assert result is None, "Should return None when no PNGs found"


def main():
    tests = [test_finds_and_moves_latest_png, test_returns_none_when_no_pngs]
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

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-gen/tests/test_codex_output_collector.py
```
Expected: `ImportError: cannot import name 'collect_latest'`.

- [ ] **Step 3: Write minimal implementation**

Create `InterSubMod/.claude/skills/image-gen/tools/codex_output_collector.py`:

```python
"""Find newest PNG across candidate dirs (created after `since_epoch`) and move it.

Use case: codex CLI may write generated PNG to `--image-dir` OR fall back to
`~/.codex/generated_images/`. This collector scans both, picks the latest PNG
created after the codex call started, and renames it to a target filename.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Iterable, Optional


def collect_latest(
    candidate_dirs: Iterable[Path],
    target_dir: Path,
    *,
    target_name: str,
    since_epoch: float,
) -> Optional[Path]:
    """Return path to moved file, or None if no eligible PNG found.

    `since_epoch` filters PNGs created at-or-after that mtime.
    """
    target_dir.mkdir(parents=True, exist_ok=True)
    candidates: list[Path] = []
    for d in candidate_dirs:
        if not d.exists():
            continue
        candidates.extend(p for p in d.glob("*.png") if p.stat().st_mtime >= since_epoch)
    if not candidates:
        return None
    latest = max(candidates, key=lambda p: p.stat().st_mtime)
    dst = target_dir / target_name
    shutil.move(str(latest), str(dst))
    return dst


CODEX_DEFAULT_DIRS = [
    Path.home() / ".codex" / "generated_images",
]


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 4:
        print("Usage: codex_output_collector.py <target_dir> <target_name> <since_epoch> [<extra_dir>...]")
        return 2
    target_dir = Path(argv[1])
    target_name = argv[2]
    since_epoch = float(argv[3])
    extras = [Path(a) for a in argv[4:]]
    dirs = extras + CODEX_DEFAULT_DIRS
    result = collect_latest(dirs, target_dir, target_name=target_name, since_epoch=since_epoch)
    if result is None:
        print("[collector] No PNG found in any candidate dir", file=sys.stderr)
        return 1
    print(f"[collector] Moved → {result}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 .claude/skills/image-gen/tests/test_codex_output_collector.py
```
Expected: `PASS: test_finds_and_moves_latest_png` and `PASS: test_returns_none_when_no_pngs`.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/image-gen/tools/codex_output_collector.py .claude/skills/image-gen/tests/test_codex_output_collector.py
git commit -m "feat(image-gen): codex_output_collector.py — rescue PNG from codex fallback dir

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 6: cairo_render.py (program track for flow_diagram + icon)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/tools/cairo_render.py`
- Create: `InterSubMod/.claude/skills/image-gen/tests/test_cairo_render.py`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-gen/tests/test_cairo_render.py`:

```python
"""Test cairo_render produces valid PNG for flow_diagram and icon."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from cairo_render import render  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def is_png(path: Path) -> bool:
    return path.exists() and path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_flow_diagram_renders():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "flow.png"
        render(FIXTURES / "valid_flow.yaml", out)
        assert is_png(out), f"Output is not a valid PNG: {out}"
        assert out.stat().st_size > 1000, f"PNG suspiciously small: {out.stat().st_size} bytes"


def test_icon_renders():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "icon.png"
        render(FIXTURES / "valid_icon.yaml", out)
        assert is_png(out), f"Output is not a valid PNG: {out}"


def test_rejects_non_cairo_type():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "bad.png"
        try:
            render(FIXTURES / "valid_concept.yaml", out)
            assert False, "Should reject backend=ai for cairo renderer"
        except ValueError as e:
            assert "backend" in str(e).lower(), f"Wrong error msg: {e}"


def main():
    tests = [test_flow_diagram_renders, test_icon_renders, test_rejects_non_cairo_type]
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

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-gen/tests/test_cairo_render.py
```
Expected: `ImportError: cannot import name 'render'`.

- [ ] **Step 3: Write minimal implementation**

Create `InterSubMod/.claude/skills/image-gen/tools/cairo_render.py`:

```python
"""Program-track image rendering via pycairo.

Supports: flow_diagram (nodes + edges), icon (shape + glyph).
Deterministic, no AI randomness, no API quota. Fast.
"""
from __future__ import annotations

import math
from pathlib import Path

import cairo
import yaml


def _parse_size(size_str: str) -> tuple[int, int]:
    w, h = size_str.lower().split("x")
    return int(w), int(h)


def _parse_hex(c: str) -> tuple[float, float, float, float]:
    """Hex color #RRGGBB or transparent → RGBA tuple in 0-1 range."""
    if c == "transparent":
        return (0.0, 0.0, 0.0, 0.0)
    c = c.lstrip("#")
    r, g, b = int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)
    return (r / 255, g / 255, b / 255, 1.0)


def _draw_flow(ctx: cairo.Context, data: dict, w: int, h: int) -> None:
    style = data.get("style", {})
    palette = style.get("palette") or ["#1f2937", "#3b82f6"]
    text_rgba = _parse_hex(palette[0])
    accent_rgba = _parse_hex(palette[1])
    bg_rgba = _parse_hex(style.get("background", "white") if style.get("background") != "white" else "#ffffff")
    font_size = float(style.get("font_size", 14))
    font_family = style.get("font_family", "DejaVu Sans")

    # Background
    ctx.set_source_rgba(*bg_rgba)
    ctx.rectangle(0, 0, w, h)
    ctx.fill()

    ctx.select_font_face(font_family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(font_size)

    nodes = data.get("nodes", [])
    node_by_id = {n["id"]: n for n in nodes}

    # Draw edges first (under nodes)
    for e in data.get("edges", []):
        a = node_by_id.get(e["from"])
        b = node_by_id.get(e["to"])
        if not a or not b:
            continue
        ctx.set_source_rgba(*accent_rgba)
        ctx.set_line_width(2.0)
        ctx.move_to(a["x"], a["y"])
        ctx.line_to(b["x"], b["y"])
        ctx.stroke()
        # Arrow head
        angle = math.atan2(b["y"] - a["y"], b["x"] - a["x"])
        ah = 10.0
        ctx.move_to(b["x"], b["y"])
        ctx.line_to(b["x"] - ah * math.cos(angle - math.pi / 7),
                    b["y"] - ah * math.sin(angle - math.pi / 7))
        ctx.line_to(b["x"] - ah * math.cos(angle + math.pi / 7),
                    b["y"] - ah * math.sin(angle + math.pi / 7))
        ctx.close_path()
        ctx.fill()

    # Draw nodes
    for n in nodes:
        label = n.get("label", "")
        ext = ctx.text_extents(label)
        pad_x, pad_y = 16.0, 12.0
        box_w = ext.width + 2 * pad_x
        box_h = ext.height + 2 * pad_y
        x0 = n["x"] - box_w / 2
        y0 = n["y"] - box_h / 2
        # Node bg
        ctx.set_source_rgba(*bg_rgba)
        ctx.rectangle(x0, y0, box_w, box_h)
        ctx.fill()
        # Node border
        ctx.set_source_rgba(*text_rgba)
        ctx.set_line_width(1.5)
        ctx.rectangle(x0, y0, box_w, box_h)
        ctx.stroke()
        # Label
        ctx.move_to(n["x"] - ext.width / 2 - ext.x_bearing,
                    n["y"] + ext.height / 2 - (ext.height + ext.y_bearing))
        ctx.show_text(label)


def _draw_icon(ctx: cairo.Context, data: dict, w: int, h: int) -> None:
    style = data.get("style", {})
    fill_rgba = _parse_hex(data.get("fill_color", "#10b981"))
    glyph_rgba = _parse_hex(data.get("glyph_color", "#ffffff"))
    bg = style.get("background", "transparent")
    bg_rgba = _parse_hex(bg) if bg != "white" else (1, 1, 1, 1)
    pad = float(style.get("padding_px", 24))
    shape = data.get("shape", "circle")
    glyph = data.get("glyph", "?")

    # Background
    ctx.set_source_rgba(*bg_rgba)
    ctx.rectangle(0, 0, w, h)
    ctx.fill()

    # Shape
    ctx.set_source_rgba(*fill_rgba)
    cx, cy = w / 2, h / 2
    radius = (min(w, h) - 2 * pad) / 2
    if shape == "circle":
        ctx.arc(cx, cy, radius, 0, 2 * math.pi)
        ctx.fill()
    elif shape in ("square", "rounded_square"):
        x0, y0 = pad, pad
        side = min(w, h) - 2 * pad
        if shape == "rounded_square":
            r = side * 0.15
            ctx.move_to(x0 + r, y0)
            ctx.line_to(x0 + side - r, y0)
            ctx.arc(x0 + side - r, y0 + r, r, -math.pi / 2, 0)
            ctx.line_to(x0 + side, y0 + side - r)
            ctx.arc(x0 + side - r, y0 + side - r, r, 0, math.pi / 2)
            ctx.line_to(x0 + r, y0 + side)
            ctx.arc(x0 + r, y0 + side - r, r, math.pi / 2, math.pi)
            ctx.line_to(x0, y0 + r)
            ctx.arc(x0 + r, y0 + r, r, math.pi, 3 * math.pi / 2)
            ctx.close_path()
            ctx.fill()
        else:
            ctx.rectangle(x0, y0, side, side)
            ctx.fill()
    else:
        # Fallback: circle
        ctx.arc(cx, cy, radius, 0, 2 * math.pi)
        ctx.fill()

    # Glyph
    ctx.set_source_rgba(*glyph_rgba)
    font_size = radius * 1.3
    ctx.select_font_face("DejaVu Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    ctx.set_font_size(font_size)
    ext = ctx.text_extents(glyph)
    ctx.move_to(cx - ext.width / 2 - ext.x_bearing,
                cy + ext.height / 2 - (ext.height + ext.y_bearing))
    ctx.show_text(glyph)


def render(yaml_path: Path, output: Path) -> None:
    """Render YAML prompt to PNG.

    Raises ValueError if backend is not 'cairo'.
    """
    data = yaml.safe_load(Path(yaml_path).read_text())
    if data.get("backend") != "cairo":
        raise ValueError(
            f"cairo_render only handles backend='cairo'; got '{data.get('backend')}' "
            f"in {yaml_path}"
        )
    out = data.get("output", {})
    w, h = _parse_size(out.get("size", "1024x768"))
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, w, h)
    ctx = cairo.Context(surface)

    t = data.get("type")
    if t == "flow_diagram":
        _draw_flow(ctx, data, w, h)
    elif t == "icon":
        _draw_icon(ctx, data, w, h)
    else:
        raise ValueError(f"cairo_render does not support type='{t}'")

    surface.write_to_png(str(output))


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 3:
        print("Usage: cairo_render.py <prompt.yaml> <output.png>")
        return 2
    render(Path(argv[1]), Path(argv[2]))
    print(f"[cairo] rendered → {argv[2]}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 .claude/skills/image-gen/tests/test_cairo_render.py
```
Expected: 3 PASS lines (flow renders, icon renders, rejects non-cairo type).

- [ ] **Step 5: Smoke test — actually look at the PNG**

```bash
python3 .claude/skills/image-gen/tools/cairo_render.py \
  .claude/skills/image-gen/tests/fixtures/valid_flow.yaml \
  /tmp/smoke_flow.png
ls -la /tmp/smoke_flow.png
file /tmp/smoke_flow.png
```
Expected: PNG file ~5-30 KB, `file` reports `PNG image data, 1024 x 400`.

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/image-gen/tools/cairo_render.py .claude/skills/image-gen/tests/test_cairo_render.py
git commit -m "feat(image-gen): cairo_render.py — program track for flow_diagram + icon

Deterministic rendering via pycairo, no AI randomness or quota cost.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 7: pil_postprocess.py (resize + base64 helpers)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/tools/pil_postprocess.py`
- Create: `InterSubMod/.claude/skills/image-gen/tests/test_pil_postprocess.py`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-gen/tests/test_pil_postprocess.py`:

```python
"""Test pil_postprocess.resize_to_target and to_base64."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from pil_postprocess import resize_to_target, to_base64  # noqa
from PIL import Image


def make_test_png(path: Path, w: int, h: int) -> None:
    Image.new("RGB", (w, h), (200, 200, 200)).save(path, "PNG")


def test_resize_preserves_aspect_or_pads():
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.png"
        dst = Path(tmp) / "out.png"
        make_test_png(src, 800, 400)  # 2:1
        resize_to_target(src, dst, target_size=(1024, 1024))
        img = Image.open(dst)
        assert img.size == (1024, 1024), f"Expected 1024x1024, got {img.size}"


def test_to_base64_starts_with_data_uri():
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "in.png"
        make_test_png(src, 100, 100)
        b64 = to_base64(src)
        assert b64.startswith("data:image/png;base64,"), f"Bad prefix: {b64[:40]}"
        assert len(b64) > 100, "Base64 too short"


def main():
    tests = [test_resize_preserves_aspect_or_pads, test_to_base64_starts_with_data_uri]
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

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-gen/tests/test_pil_postprocess.py
```
Expected: `ImportError`.

- [ ] **Step 3: Write minimal implementation**

Create `InterSubMod/.claude/skills/image-gen/tools/pil_postprocess.py`:

```python
"""PIL post-processing: resize to standard sizes + base64 inline encoding."""
from __future__ import annotations

import base64
from pathlib import Path

from PIL import Image


def resize_to_target(src: Path, dst: Path, target_size: tuple[int, int]) -> None:
    """Resize src PNG to target_size. Preserves aspect ratio via thumbnail + padding (white bg)."""
    img = Image.open(src)
    img.thumbnail(target_size, Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", target_size, (255, 255, 255))
    paste_x = (target_size[0] - img.size[0]) // 2
    paste_y = (target_size[1] - img.size[1]) // 2
    if img.mode in ("RGBA", "LA"):
        canvas.paste(img, (paste_x, paste_y), img.convert("RGBA").split()[-1])
    else:
        canvas.paste(img, (paste_x, paste_y))
    canvas.save(dst, "PNG", optimize=True)


def to_base64(src: Path) -> str:
    """Return data:image/png;base64,... string."""
    return "data:image/png;base64," + base64.b64encode(Path(src).read_bytes()).decode("ascii")


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 4:
        print("Usage: pil_postprocess.py resize <src.png> <dst.png> [<W>x<H>]")
        print("       pil_postprocess.py base64 <src.png>")
        return 2
    cmd = argv[1]
    if cmd == "resize":
        size = (1024, 1024)
        if len(argv) >= 5:
            w, h = argv[4].lower().split("x")
            size = (int(w), int(h))
        resize_to_target(Path(argv[2]), Path(argv[3]), size)
        print(f"[pil] resized → {argv[3]} ({size[0]}x{size[1]})")
        return 0
    if cmd == "base64":
        print(to_base64(Path(argv[2])))
        return 0
    print(f"Unknown cmd: {cmd}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 .claude/skills/image-gen/tests/test_pil_postprocess.py
```
Expected: 2 PASS lines.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/image-gen/tools/pil_postprocess.py .claude/skills/image-gen/tests/test_pil_postprocess.py
git commit -m "feat(image-gen): pil_postprocess.py — resize + base64 helpers (PIL)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 8: dispatch.py (main entry, routes by type to AI/cairo)

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/tools/dispatch.py`
- Create: `InterSubMod/.claude/skills/image-gen/tests/test_dispatch.py`

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-gen/tests/test_dispatch.py`:

```python
"""Test dispatch routes by backend field to correct subprocess."""
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from dispatch import route_one, plan_run  # noqa

FIXTURES = Path(__file__).parent / "fixtures"


def test_plan_run_classifies_correctly():
    plan = plan_run(FIXTURES, output_dir=Path("/tmp/out"))
    by_backend = {p["backend"]: 0 for p in plan}
    for p in plan:
        by_backend[p["backend"]] = by_backend.get(p["backend"], 0) + 1
    assert by_backend.get("ai", 0) == 2, f"Expected 2 AI, got {by_backend}"
    assert by_backend.get("cairo", 0) == 2, f"Expected 2 cairo, got {by_backend}"


def test_route_one_cairo_invokes_cairo_render():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "icon.png"
        route_one(FIXTURES / "valid_icon.yaml", out, dry_run=False)
        assert out.exists() and out.stat().st_size > 500, f"icon PNG missing/too small: {out}"


def test_route_one_ai_dry_run_does_not_call_codex():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "concept.png"
        with patch("dispatch.subprocess.run") as mock_run:
            route_one(FIXTURES / "valid_concept.yaml", out, dry_run=True)
            mock_run.assert_not_called()


def main():
    tests = [test_plan_run_classifies_correctly, test_route_one_cairo_invokes_cairo_render, test_route_one_ai_dry_run_does_not_call_codex]
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

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-gen/tests/test_dispatch.py
```
Expected: `ImportError`.

- [ ] **Step 3: Write minimal implementation**

Create `InterSubMod/.claude/skills/image-gen/tools/dispatch.py`:

```python
"""Main dispatch entry for image-gen skill.

Reads prompts from <prompts_dir>/*.yaml, routes by `backend` field:
  - ai    → codex_image_gen.run + codex_output_collector.collect_latest
  - cairo → cairo_render.render
Then pil_postprocess.resize_to_target standardizes output sizes.

Cost preview before run.
Idempotent: skips if target PNG exists (unless --force).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import yaml

# Add tools dir to import path so siblings work when run as script
_TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(_TOOLS))

import cairo_render  # noqa: E402
import codex_image_gen  # noqa: E402
import codex_output_collector  # noqa: E402
import pil_postprocess  # noqa: E402


def plan_run(prompts_dir: Path, output_dir: Path) -> list[dict]:
    """Return list of dicts {prompt_path, backend, target_path}."""
    plan = []
    for p in sorted(prompts_dir.glob("*.yaml")):
        data = yaml.safe_load(p.read_text())
        backend = data.get("backend", "ai")
        target = output_dir / (p.stem + ".png")
        plan.append({"prompt_path": p, "backend": backend, "target_path": target})
    return plan


def route_one(prompt_path: Path, target_path: Path, *, dry_run: bool = False) -> int:
    """Dispatch one prompt to its backend, write PNG to target_path."""
    data = yaml.safe_load(Path(prompt_path).read_text())
    backend = data.get("backend", "ai")
    target_path.parent.mkdir(parents=True, exist_ok=True)

    if backend == "cairo":
        if dry_run:
            print(f"DRY RUN [cairo]: {prompt_path.name} → {target_path}")
            return 0
        cairo_render.render(prompt_path, target_path)
        print(f"[cairo] {prompt_path.name} → {target_path}")
        return 0

    # AI backend: codex exec + collector
    if dry_run:
        cmd = codex_image_gen.build_command(prompt_path, target_path.parent)
        print(f"DRY RUN [ai]: {' '.join(repr(c) for c in cmd)}")
        return 0

    started_at = time.time()
    rc = codex_image_gen.run(prompt_path, target_path.parent)
    if rc != 0:
        print(f"[ai] codex exec failed for {prompt_path.name} (exit {rc})", file=sys.stderr)
        return rc
    # Rescue PNG if codex wrote elsewhere
    candidate_dirs = [target_path.parent] + codex_output_collector.CODEX_DEFAULT_DIRS
    moved = codex_output_collector.collect_latest(
        candidate_dirs, target_path.parent,
        target_name=target_path.name, since_epoch=started_at,
    )
    if moved is None:
        print(f"[ai] PNG not found in {candidate_dirs} for {prompt_path.name}", file=sys.stderr)
        return 1
    print(f"[ai] {prompt_path.name} → {moved}")
    return 0


def cost_preview(plan: list[dict]) -> str:
    n_ai = sum(1 for p in plan if p["backend"] == "ai")
    n_cairo = sum(1 for p in plan if p["backend"] == "cairo")
    return (
        f"[image-gen] {len(plan)} prompts: AI×{n_ai} (codex $imagegen, ~30k tokens/each, eats ChatGPT quota), "
        f"cairo×{n_cairo} (zero cost). Output dir: see plan."
    )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="image-gen dispatch")
    parser.add_argument("prompts_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="re-render even if target PNG exists")
    parser.add_argument("--postprocess-size", default=None,
                        help="resize all outputs to e.g. 1024x1024")
    args = parser.parse_args(argv[1:])

    plan = plan_run(args.prompts_dir, args.output_dir)
    print(cost_preview(plan))
    if args.dry_run:
        for p in plan:
            print(f"  {p['backend']:5s}  {p['prompt_path'].name} → {p['target_path']}")
        return 0

    failed = 0
    for p in plan:
        if p["target_path"].exists() and not args.force:
            print(f"[skip] {p['target_path']} (use --force to re-render)")
            continue
        rc = route_one(p["prompt_path"], p["target_path"])
        if rc != 0:
            failed += 1

    if args.postprocess_size:
        w, h = args.postprocess_size.lower().split("x")
        size = (int(w), int(h))
        for p in plan:
            if p["target_path"].exists():
                pil_postprocess.resize_to_target(p["target_path"], p["target_path"], size)
                print(f"[pil] resized {p['target_path']} → {size}")

    if failed:
        print(f"[image-gen] {failed} of {len(plan)} prompts failed", file=sys.stderr)
        return 1
    print(f"[image-gen] {len(plan)} prompts succeeded")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 .claude/skills/image-gen/tests/test_dispatch.py
```
Expected: 3 PASS lines.

- [ ] **Step 5: Smoke test cairo path end-to-end**

```bash
mkdir -p /tmp/dispatch_smoke
python3 .claude/skills/image-gen/tools/dispatch.py \
  .claude/skills/image-gen/tests/fixtures /tmp/dispatch_smoke --dry-run
```
Expected: cost preview line + 4 plan lines (2 AI dry-run, 2 cairo dry-run).

Then run cairo-only for real (skip AI to avoid quota use during tests):
```bash
# Manually filter to just cairo fixtures
mkdir -p /tmp/dispatch_cairo_only
cp .claude/skills/image-gen/tests/fixtures/valid_flow.yaml /tmp/dispatch_cairo_only/
cp .claude/skills/image-gen/tests/fixtures/valid_icon.yaml /tmp/dispatch_cairo_only/
mkdir -p /tmp/dispatch_smoke_real
python3 .claude/skills/image-gen/tools/dispatch.py \
  /tmp/dispatch_cairo_only /tmp/dispatch_smoke_real
ls -la /tmp/dispatch_smoke_real/
```
Expected: 2 PNGs (valid_flow.png, valid_icon.png), each ~5-30 KB.

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/image-gen/tools/dispatch.py .claude/skills/image-gen/tests/test_dispatch.py
git commit -m "feat(image-gen): dispatch.py main entry — routes AI/cairo by YAML backend field

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 9: image-gen SKILL.md + playbook

**Files:**
- Create: `InterSubMod/.claude/skills/image-gen/SKILL.md`
- Create: `InterSubMod/.claude/skills/image-gen/playbook.md`
- Create: `InterSubMod/.claude/skills/image-gen/prompts/prompt_quality_self_review.md`

- [ ] **Step 1: Write SKILL.md**

Create `InterSubMod/.claude/skills/image-gen/SKILL.md`:

```markdown
---
name: image-gen
description: |
  生成示意圖 / 流程圖 / 概念示意 / icon 圖片。雙軌：concept_diagram & data_mockup → AI 軌 (codex CLI \$imagegen via ChatGPT OAuth, 免 USD)；flow_diagram & icon → 程式軌 (pycairo, 零成本、不錯字、deterministic)。輸出 PNG 到指定 figures/ 目錄。
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
```

- [ ] **Step 2: Write playbook.md**

Create `InterSubMod/.claude/skills/image-gen/playbook.md`:

```markdown
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
```

- [ ] **Step 3: Write prompts/prompt_quality_self_review.md**

Create `InterSubMod/.claude/skills/image-gen/prompts/prompt_quality_self_review.md`:

```markdown
# Prompt Quality Self-Review (run before invoking codex)

For each prompt YAML about to be sent to AI track, the agent should self-check:

1. **Subject is declarative single sentence** — not imperative, not vague. "LOH region triggers ..." not "draw a chromosome".
2. **Constraints include all 4 Anthropic taboos** — `no gradient overuse`, `no glass morphism`, `no multi-indigo stacking`, `colorblind-friendly`.
3. **Output size matches use case** — 1024×1024 default; only larger for hero figures.
4. **No emoji or decorative unicode** in `subject` or `labels`.
5. **Palette ≤ 3 colors** explicit hex codes.

If any check fails, fix the YAML before invoking codex. Cost saving: prompts/ is in git, fixing now is 0 token; fixing after generation is ~30k token waste.
```

- [ ] **Step 4: Manual sanity check — does the SKILL.md description trigger?**

The description should match common phrasings. Quickly grep for trigger keywords:

```bash
grep -oE "(需要示意圖|畫個|生成圖片|生圖|補張圖|插圖)" .claude/skills/image-gen/SKILL.md | sort -u
```
Expected: matches present (test that triggers are in description).

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/image-gen/SKILL.md .claude/skills/image-gen/playbook.md .claude/skills/image-gen/prompts/
git commit -m "feat(image-gen): SKILL.md + playbook + prompt_quality_self_review

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 10: image-vision-check SKILL.md + 4 checklists + playbook

**Files:**
- Create: `InterSubMod/.claude/skills/image-vision-check/SKILL.md`
- Create: `InterSubMod/.claude/skills/image-vision-check/playbook.md`
- Create: `InterSubMod/.claude/skills/image-vision-check/checklists/{concept_diagram,flow_diagram,data_mockup,icon}_check.md`
- Create: `InterSubMod/.claude/skills/image-vision-check/schemas/quality.schema.json`
- Create: `InterSubMod/.claude/skills/image-vision-check/prompts/vision_review_master.md`

- [ ] **Step 1: Write SKILL.md**

Create `InterSubMod/.claude/skills/image-vision-check/SKILL.md`:

```markdown
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
```

- [ ] **Step 2: Write 4 checklist .md files**

Create `InterSubMod/.claude/skills/image-vision-check/checklists/concept_diagram_check.md`:

```markdown
# concept_diagram Check Checklist

Score each dimension 1 = pass, 0 = fail. Need ≥4 to pass.

1. **Prompt fidelity** — Image content matches the `subject` declarative sentence + lists key elements.
2. **Readability** — All text labels (top, bottom, callouts) are legible, no garbled CJK / no broken Latin.
3. **Focal clarity** — Main subject clearly distinguished from background; visual hierarchy obvious.
4. **Design taboos avoided** — No gradient overuse, no glass morphism, no multi-indigo stacking, no emoji-decorated headers.
5. **Print-friendly** — Contrast ≥ 4.5:1, not relying on red/green color-coding alone.
6. **Cross-figure consistency** — Style (font, palette, line weight) matches the other figures in the same topic folder. (For first image, score = 1 by default.)
```

Create `InterSubMod/.claude/skills/image-vision-check/checklists/flow_diagram_check.md`:

```markdown
# flow_diagram Check Checklist

Score each dimension 1 = pass, 0 = fail. Need ≥4 to pass.

1. **Prompt fidelity** — All `nodes` rendered with their `label`s; all `edges` connect from→to correctly with arrows.
2. **Readability** — Node labels not clipped, font size readable at intended display size (~14pt baseline).
3. **Focal clarity** — Edges visually connect parent → child without crossing important nodes; arrowheads visible.
4. **Design taboos avoided** — Same as concept_diagram (no gradient/glass-morphism/multi-indigo).
5. **Print-friendly** — Black-on-white printable; contrast OK.
6. **Cross-figure consistency** — Same node style / same arrow style across figures.

Note: cairo-rendered diagrams should ALWAYS pass dimension 2 (readability) since text is rasterized exact. If 2 fails, cairo render bug.
```

Create `InterSubMod/.claude/skills/image-vision-check/checklists/data_mockup_check.md`:

```markdown
# data_mockup Check Checklist

Score each dimension 1 = pass, 0 = fail. Need ≥4 to pass.

1. **Prompt fidelity** — Chart type matches `chart_type`; trend matches `data_hint`.
2. **Readability** — Axis labels match `axes.x_label`/`axes.y_label`; no garbled text.
3. **Focal clarity** — Bars/lines clearly distinguishable; legend present (if multi-series).
4. **Design taboos avoided** — No gradient on bars (flat fills), no 3D effects, no glass-morphism.
5. **Print-friendly** — Patterns/textures or hatching used IN ADDITION to color (so colorblind/print-OK).
6. **Cross-figure consistency** — Same axis style / palette as other data_mockups.

Note: Mock-up nature should be visible — bars look "concept-y" not pixel-perfect quant.
```

Create `InterSubMod/.claude/skills/image-vision-check/checklists/icon_check.md`:

```markdown
# icon Check Checklist

Score each dimension 1 = pass, 0 = fail. Need ≥4 to pass.

1. **Prompt fidelity** — Shape matches `shape` field; glyph matches `glyph` field exactly.
2. **Readability** — Glyph occupies > 40% of inner area, clearly visible against fill.
3. **Focal clarity** — Single dominant element; no unintended texture/border.
4. **Design taboos avoided** — Solid fill (not gradient), no drop shadows, no glow effects.
5. **Print-friendly** — Glyph contrast vs fill ≥ 4.5:1; works in B&W (luminance contrast OK).
6. **Cross-figure consistency** — Same icon set style across multiple icons (size, padding, glyph weight).

Note: cairo-rendered icons should ALWAYS pass 1, 2, 3, 4 (deterministic). If any fails, cairo render bug.
```

- [ ] **Step 3: Write quality.schema.json**

Create `InterSubMod/.claude/skills/image-vision-check/schemas/quality.schema.json`:

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "title": "image-vision-check quality.json entry",
  "type": "object",
  "required": ["image", "version", "checks", "score", "verdict"],
  "properties": {
    "image": { "type": "string", "description": "filename relative to figures dir" },
    "version": { "type": "integer", "minimum": 1, "description": "1 = original, 2 = first retry, etc." },
    "type": {
      "type": "string",
      "enum": ["concept_diagram", "flow_diagram", "data_mockup", "icon"]
    },
    "backend": {
      "type": "string",
      "enum": ["ai", "cairo"]
    },
    "checks": {
      "type": "object",
      "required": [
        "prompt_fidelity",
        "readability",
        "focal_clarity",
        "design_taboos_avoided",
        "print_friendly",
        "cross_figure_consistency"
      ],
      "properties": {
        "prompt_fidelity": { "type": "boolean" },
        "readability": { "type": "boolean" },
        "focal_clarity": { "type": "boolean" },
        "design_taboos_avoided": { "type": "boolean" },
        "print_friendly": { "type": "boolean" },
        "cross_figure_consistency": { "type": "boolean" }
      }
    },
    "score": { "type": "string", "pattern": "^[0-6]/6$" },
    "verdict": {
      "type": "string",
      "enum": ["pass", "partial", "fail", "needs_human"]
    },
    "suggestions": {
      "type": "array",
      "items": { "type": "string" }
    },
    "retry_history": {
      "type": "array",
      "items": {
        "type": "object",
        "required": ["version", "score"],
        "properties": {
          "version": { "type": "integer" },
          "score": { "type": "string" }
        }
      }
    },
    "checked_at": { "type": "string", "format": "date-time" }
  }
}
```

- [ ] **Step 4: Write vision_review_master.md prompt**

Create `InterSubMod/.claude/skills/image-vision-check/prompts/vision_review_master.md`:

```markdown
# Vision Review Master Prompt (used by Claude when reading images)

You are checking an image generated for a research report. Score 6 dimensions, 1 = pass / 0 = fail. Be honest — partial credit is reflected in the 4/6 threshold.

## Inputs (provided per call)
- Image (you read via `Read` tool)
- Prompt YAML (subject, type, constraints)
- Image type-specific checklist (one of `checklists/*.md`)

## Output (you must produce as JSON only, no prose)

```json
{
  "image": "<filename>",
  "version": <integer>,
  "type": "<type>",
  "backend": "<ai|cairo>",
  "checks": {
    "prompt_fidelity": <true|false>,
    "readability": <true|false>,
    "focal_clarity": <true|false>,
    "design_taboos_avoided": <true|false>,
    "print_friendly": <true|false>,
    "cross_figure_consistency": <true|false>
  },
  "score": "<n>/6",
  "verdict": "<pass|partial|fail>",
  "suggestions": [
    "<suggestion 1, only for failed checks, with concrete prompt edit>"
  ],
  "checked_at": "<ISO 8601 timestamp>"
}
```

## Decision rules (apply strictly)

- score = sum(checks where value == true), expressed as `"<n>/6"`.
- verdict:
  - 6/6 → `pass`
  - 4-5/6 → `partial`
  - ≤3/6 → `fail`
- For each false check, write a 1-sentence `suggestion` that includes the YAML field to edit (e.g., `constraints: add 'no drop shadow'`).
- For first image of a topic folder, `cross_figure_consistency` is `true` by default.

## Common mistakes to avoid

- Do NOT score on subjective beauty.
- Do NOT add bonus dimensions (e.g., "innovation", "creativity") — strict 6 only.
- Do NOT mention being an AI in suggestions — write as if directly editing the prompt YAML.
```

- [ ] **Step 5: Write playbook.md**

Create `InterSubMod/.claude/skills/image-vision-check/playbook.md`:

```markdown
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
- If `all_pass`, suggest `/html-preview` (Phase 2) or report DONE.
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
```

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/image-vision-check/SKILL.md .claude/skills/image-vision-check/playbook.md .claude/skills/image-vision-check/checklists/ .claude/skills/image-vision-check/schemas/ .claude/skills/image-vision-check/prompts/
git commit -m "feat(image-vision-check): SKILL.md + playbook + 4 checklists + schema + master prompt

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 11: vision_check_runner.py (orchestrate Claude vision audit)

**Files:**
- Create: `InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py`
- Create: `InterSubMod/.claude/skills/image-vision-check/tests/test_vision_check_runner.py`
- Create: `InterSubMod/.claude/skills/image-vision-check/tests/fixtures/sample_quality.json`

Note: actual Claude `Read` invocation happens through the Claude Code agent (not via Python API). The runner's responsibility is to:
1. Enumerate images.
2. Output a structured "instruction file" the agent can pick up to do per-image vision review.
3. Validate Claude's resulting JSON against the schema.
4. Aggregate into `quality.json`.

So this runner is mostly file-handling + validation; the actual vision work is invoked by the calling Claude agent reading the instruction file.

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-vision-check/tests/test_vision_check_runner.py`:

```python
"""Test runner: enumerate, validate JSON, aggregate."""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from vision_check_runner import (  # noqa
    enumerate_images,
    validate_entry,
    aggregate,
)

SCHEMA_PATH = Path(__file__).parent.parent / "schemas" / "quality.schema.json"


def make_test_png(path: Path):
    from PIL import Image
    Image.new("RGB", (32, 32), (200, 100, 50)).save(path, "PNG")


def test_enumerate_images_finds_pngs():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        make_test_png(d / "fig1.png")
        make_test_png(d / "fig2.png")
        (d / "not_a_png.txt").write_text("ignored")
        out = enumerate_images(d)
        names = sorted(p.name for p in out)
        assert names == ["fig1.png", "fig2.png"], f"Wrong: {names}"


def test_validate_entry_pass():
    entry = {
        "image": "fig1.png",
        "version": 1,
        "type": "concept_diagram",
        "backend": "ai",
        "checks": {
            "prompt_fidelity": True,
            "readability": True,
            "focal_clarity": True,
            "design_taboos_avoided": True,
            "print_friendly": True,
            "cross_figure_consistency": True,
        },
        "score": "6/6",
        "verdict": "pass",
        "suggestions": [],
        "checked_at": "2026-05-10T14:00:00Z",
    }
    assert validate_entry(entry, SCHEMA_PATH) == [], "Should be valid"


def test_validate_entry_fail_bad_score():
    entry = {
        "image": "fig1.png",
        "version": 1,
        "checks": {
            "prompt_fidelity": True,
            "readability": True,
            "focal_clarity": True,
            "design_taboos_avoided": True,
            "print_friendly": True,
            "cross_figure_consistency": True,
        },
        "score": "7/6",  # invalid
        "verdict": "pass",
    }
    errors = validate_entry(entry, SCHEMA_PATH)
    assert errors, "Should report invalid score format"


def test_aggregate_summary():
    entries = [
        {"verdict": "pass"},
        {"verdict": "pass"},
        {"verdict": "partial"},
        {"verdict": "fail"},
    ]
    s = aggregate(entries)
    assert s["n_pass"] == 2 and s["n_partial"] == 1 and s["n_fail"] == 1
    assert s["overall"] == "mixed"
    s2 = aggregate([{"verdict": "pass"}, {"verdict": "pass"}])
    assert s2["overall"] == "all_pass"


def main():
    tests = [test_enumerate_images_finds_pngs, test_validate_entry_pass, test_validate_entry_fail_bad_score, test_aggregate_summary]
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

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-vision-check/tests/test_vision_check_runner.py
```
Expected: `ImportError`.

- [ ] **Step 3: Write minimal implementation**

Create `InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py`:

```python
"""vision_check_runner — enumerate images, validate Claude-produced quality JSON, aggregate.

Note on architecture:
This runner does NOT itself call Claude vision. The pattern is:
1. Runner enumerates PNGs and writes `<figures_dir>/_check_instructions.json`.
2. Calling Claude agent reads the instructions, uses its own `Read` tool on each
   PNG, and writes per-image JSON entries (matching `schemas/quality.schema.json`)
   into `<figures_dir>/_pending_entries/<basename>.json`.
3. Runner is re-invoked with `--aggregate` to validate + merge into final
   `<figures_dir>/quality.json`.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Optional jsonschema; we degrade gracefully if missing.
try:
    import jsonschema as _js  # type: ignore
    _HAS_JSONSCHEMA = True
except ImportError:
    _HAS_JSONSCHEMA = False


def enumerate_images(figures_dir: Path) -> list[Path]:
    return sorted(figures_dir.glob("*.png"))


def validate_entry(entry: dict, schema_path: Path) -> list[str]:
    """Return list of error strings; empty list means valid."""
    if not _HAS_JSONSCHEMA:
        # Manual minimal check
        required = {"image", "version", "checks", "score", "verdict"}
        missing = required - entry.keys()
        if missing:
            return [f"missing required fields: {missing}"]
        score = entry.get("score", "")
        if not (isinstance(score, str) and len(score) == 3 and score[1] == "/" and score[0].isdigit() and score[2] == "6" and 0 <= int(score[0]) <= 6):
            return [f"invalid score: {score}"]
        return []
    schema = json.loads(Path(schema_path).read_text())
    validator = _js.Draft202012Validator(schema)
    return [f"{e.message} at /{'/'.join(map(str, e.path))}" for e in validator.iter_errors(entry)]


def aggregate(entries: list[dict]) -> dict:
    n_pass = sum(1 for e in entries if e.get("verdict") == "pass")
    n_partial = sum(1 for e in entries if e.get("verdict") == "partial")
    n_fail = sum(1 for e in entries if e.get("verdict") == "fail")
    n_human = sum(1 for e in entries if e.get("verdict") == "needs_human")
    overall = "all_pass" if (n_partial + n_fail + n_human) == 0 else "mixed"
    return {
        "n_pass": n_pass,
        "n_partial": n_partial,
        "n_fail": n_fail,
        "n_needs_human": n_human,
        "overall": overall,
    }


def write_instructions(figures_dir: Path, prompts_dir: Path | None) -> Path:
    images = enumerate_images(figures_dir)
    instructions = {
        "figures_dir": str(figures_dir.resolve()),
        "prompts_dir": str(prompts_dir.resolve()) if prompts_dir else None,
        "schema_path": str((Path(__file__).parent.parent / "schemas" / "quality.schema.json").resolve()),
        "master_prompt": str((Path(__file__).parent.parent / "prompts" / "vision_review_master.md").resolve()),
        "checklists_dir": str((Path(__file__).parent.parent / "checklists").resolve()),
        "pending_dir": str((figures_dir / "_pending_entries").resolve()),
        "images": [
            {
                "png": str(img.resolve()),
                "yaml": str((prompts_dir / (img.stem + ".yaml")).resolve()) if prompts_dir and (prompts_dir / (img.stem + ".yaml")).exists() else None,
            }
            for img in images
        ],
    }
    out = figures_dir / "_check_instructions.json"
    out.write_text(json.dumps(instructions, indent=2, ensure_ascii=False))
    (figures_dir / "_pending_entries").mkdir(exist_ok=True)
    return out


def aggregate_pending(figures_dir: Path) -> Path:
    schema_path = Path(__file__).parent.parent / "schemas" / "quality.schema.json"
    pending = figures_dir / "_pending_entries"
    if not pending.exists():
        print(f"No pending entries at {pending}", file=sys.stderr)
        sys.exit(2)

    entries: list[dict] = []
    errors: list[str] = []
    for j in sorted(pending.glob("*.json")):
        try:
            entry = json.loads(j.read_text())
        except json.JSONDecodeError as e:
            errors.append(f"{j.name}: {e}")
            continue
        errs = validate_entry(entry, schema_path)
        if errs:
            errors.append(f"{j.name}: {'; '.join(errs)}")
            continue
        entries.append(entry)

    if errors:
        print("Validation errors:")
        for e in errors:
            print(f"  {e}")

    summary = aggregate(entries)
    quality = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "figures_dir": str(figures_dir.resolve()),
        "summary": summary,
        "entries": entries,
    }
    out = figures_dir / "quality.json"
    out.write_text(json.dumps(quality, indent=2, ensure_ascii=False))
    print(f"[vision_check] aggregated {len(entries)} entries → {out}")
    print(f"[vision_check] summary: {summary}")
    return out


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="vision_check_runner")
    parser.add_argument("figures_dir", type=Path)
    parser.add_argument("--prompts", type=Path, default=None,
                        help="paired prompts dir (defaults to <figures_dir>/../prompts/)")
    parser.add_argument("--aggregate", action="store_true",
                        help="aggregate pending JSON entries into final quality.json")
    args = parser.parse_args(argv[1:])

    if not args.figures_dir.is_dir():
        print(f"figures_dir not found: {args.figures_dir}", file=sys.stderr)
        return 2

    if args.aggregate:
        aggregate_pending(args.figures_dir)
        return 0

    prompts_dir = args.prompts
    if prompts_dir is None:
        candidate = args.figures_dir.parent / "prompts"
        prompts_dir = candidate if candidate.is_dir() else None

    instructions = write_instructions(args.figures_dir, prompts_dir)
    print(f"[vision_check] instructions written → {instructions}")
    print(f"[vision_check] next: Claude agent reads instructions, processes each image,")
    print(f"[vision_check] writes to _pending_entries/<basename>.json,")
    print(f"[vision_check] then re-run this script with --aggregate")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 .claude/skills/image-vision-check/tests/test_vision_check_runner.py
```
Expected: 4 PASS lines.

- [ ] **Step 5: Smoke test instruction-write step**

```bash
mkdir -p /tmp/vc_smoke/figures /tmp/vc_smoke/prompts
python3 -c "from PIL import Image; Image.new('RGB',(64,64),(0,128,255)).save('/tmp/vc_smoke/figures/fig1.png')"
cp .claude/skills/image-gen/tests/fixtures/valid_concept.yaml /tmp/vc_smoke/prompts/fig1.yaml

python3 .claude/skills/image-vision-check/tools/vision_check_runner.py /tmp/vc_smoke/figures
cat /tmp/vc_smoke/figures/_check_instructions.json
```
Expected: `_check_instructions.json` exists with 1 image entry referencing png + yaml.

- [ ] **Step 6: Commit**

```bash
git add .claude/skills/image-vision-check/tools/ .claude/skills/image-vision-check/tests/
git commit -m "feat(image-vision-check): vision_check_runner.py + tests

Architecture: runner produces instructions; calling Claude agent does Read
work; runner re-invoked with --aggregate to validate + merge.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 12: Auto-retry loop integration helper

**Files:**
- Modify: `InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py` (add `apply_suggestions_to_yaml`)
- Create: `InterSubMod/.claude/skills/image-vision-check/tests/test_retry_loop.py`

The actual retry orchestration is done by the calling agent (it reads quality.json, calls image-gen with patched YAML, re-runs check). This task adds the YAML-patch helper.

- [ ] **Step 1: Write the failing test**

Create `InterSubMod/.claude/skills/image-vision-check/tests/test_retry_loop.py`:

```python
"""Test apply_suggestions_to_yaml appends suggestions to constraints field."""
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from vision_check_runner import apply_suggestions_to_yaml  # noqa


def test_appends_to_constraints():
    with tempfile.TemporaryDirectory() as tmp:
        yaml_path = Path(tmp) / "p.yaml"
        yaml_path.write_text(yaml.safe_dump({
            "type": "concept_diagram",
            "subject": "x",
            "constraints": ["original"],
            "output": {"size": "1024x1024", "quality": "high"},
            "backend": "ai",
        }))
        suggestions = ["constraints: add 'no drop shadow'"]
        apply_suggestions_to_yaml(yaml_path, suggestions)
        new = yaml.safe_load(yaml_path.read_text())
        assert "no drop shadow" in new["constraints"], f"Constraint not added: {new['constraints']}"
        assert "original" in new["constraints"], "Original constraint dropped"


def test_handles_missing_constraints_field():
    with tempfile.TemporaryDirectory() as tmp:
        yaml_path = Path(tmp) / "p.yaml"
        yaml_path.write_text(yaml.safe_dump({
            "type": "icon",
            "subject": "x",
            "output": {"size": "256x256", "quality": "n/a"},
            "backend": "cairo",
        }))
        apply_suggestions_to_yaml(yaml_path, ["constraints: add 'high contrast glyph'"])
        new = yaml.safe_load(yaml_path.read_text())
        assert "constraints" in new and "high contrast glyph" in new["constraints"]


def main():
    failed = 0
    for t in [test_appends_to_constraints, test_handles_missing_constraints_field]:
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

- [ ] **Step 2: Run test to verify it fails**

```bash
python3 .claude/skills/image-vision-check/tests/test_retry_loop.py
```
Expected: `ImportError: cannot import name 'apply_suggestions_to_yaml'`.

- [ ] **Step 3: Append helper to vision_check_runner.py**

Append to `InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py` (before `def main`):

```python
def apply_suggestions_to_yaml(yaml_path: Path, suggestions: list[str]) -> None:
    """Apply textual suggestions from quality.json to a prompt YAML.

    Currently supports the 'constraints: add ...' suggestion form.
    Other suggestion forms are logged but not auto-applied (D8 max 1 retry,
    keep simple).
    """
    import yaml as _yaml
    data = _yaml.safe_load(Path(yaml_path).read_text())
    if not isinstance(data, dict):
        return
    constraints = list(data.get("constraints") or [])
    for s in suggestions:
        prefix = "constraints: add"
        if s.startswith(prefix):
            tail = s[len(prefix):].strip()
            tail = tail.strip("'\"`")
            if tail and tail not in constraints:
                constraints.append(tail)
    data["constraints"] = constraints
    Path(yaml_path).write_text(_yaml.safe_dump(data, allow_unicode=True, sort_keys=False))
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python3 .claude/skills/image-vision-check/tests/test_retry_loop.py
```
Expected: 2 PASS lines.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/image-vision-check/tools/vision_check_runner.py .claude/skills/image-vision-check/tests/test_retry_loop.py
git commit -m "feat(image-vision-check): apply_suggestions_to_yaml helper for D8 auto-retry

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 13: Phase 1 demo topic folder (hand-written, end-to-end validation)

**Files:**
- Create: `InterSubMod/examples/phase1_demo/demo_topic/README.md`
- Create: `InterSubMod/examples/phase1_demo/demo_topic/demo_topic.md`
- Create: `InterSubMod/examples/phase1_demo/demo_topic/index.html`
- Create: `InterSubMod/examples/phase1_demo/demo_topic/prompts/{fig1_concept_loh,fig2_flow_pipeline,fig3_data_mockup,fig4_icon_status}.yaml`

This is the end-to-end **閉環驗證** PI requested.

- [ ] **Step 1: Write 4 prompt YAMLs**

Create `InterSubMod/examples/phase1_demo/demo_topic/prompts/fig1_concept_loh.yaml`:

```yaml
type: concept_diagram
subject: "Schematic showing how a Loss-of-Heterozygosity (LOH) region in a tumor sample causes haplotype phasing degradation: paired short reads above an ideogram of chromosome 8, with a shaded LOH band where read pairs lose haplotype consistency."
key_elements:
  - "ideogram of chromosome 8 with q-arm shaded gray for LOH region"
  - "paired short reads (red and blue) aligned above and below the chromosome"
  - "an arrow inside the LOH region showing haplotype assignment becoming ambiguous"
labels:
  top: "Tumor sample with LOH (chr8q)"
  bottom: "Phasing degrades inside the LOH band"
style:
  palette: ["#1f2937", "#ef4444"]
  background: "white"
  composition: "centered, ideogram horizontal across the middle"
constraints:
  - "no gradient overuse"
  - "no glass morphism"
  - "no multi-indigo stacking"
  - "colorblind-friendly (no red/green pair)"
  - "print-friendly contrast >= 4.5:1"
  - "system-ui font feel"
output:
  size: "1024x1024"
  quality: "high"
backend: ai
```

Create `InterSubMod/examples/phase1_demo/demo_topic/prompts/fig2_flow_pipeline.yaml`:

```yaml
type: flow_diagram
subject: "Self-Phasing V5 minimal pipeline (G1-G3 demo)"
nodes:
  - id: "g1"
    label: "G1: ReadParser"
    x: 150
    y: 200
  - id: "g2"
    label: "G2: VCF Loader"
    x: 460
    y: 200
  - id: "g3"
    label: "G3: Phase Assignment"
    x: 800
    y: 200
edges:
  - from: "g1"
    to: "g2"
    label: ""
  - from: "g2"
    to: "g3"
    label: ""
style:
  palette: ["#1f2937", "#3b82f6"]
  background: "white"
  font_family: "DejaVu Sans"
  font_size: 16
constraints:
  - "deterministic layout"
  - "exact text rendering"
output:
  size: "1024x400"
  quality: "high"
backend: cairo
```

Create `InterSubMod/examples/phase1_demo/demo_topic/prompts/fig3_data_mockup.yaml`:

```yaml
type: data_mockup
subject: "Mock-up bar chart showing 7 cancer sample baseline AUC distribution (illustrative only, not real data)"
chart_type: "bar"
axes:
  x_label: "Sample"
  y_label: "Baseline AUC"
data_hint: "HCC1395 ~0.78 (highest), COLO829 ~0.55 (lowest), other 5 samples (DORADO, H2009, HCC1954, HCC1937, HCC2218) between 0.60-0.72; bars in dark teal flat fill, no gradient"
style:
  palette: ["#1f2937", "#10b981"]
  background: "white"
constraints:
  - "axes clearly labeled"
  - "x-axis sample names readable, no overlap"
  - "no gradient overuse"
  - "no 3D effect"
  - "data is conceptual mock, not real"
output:
  size: "1024x1024"
  quality: "high"
backend: ai
```

Create `InterSubMod/examples/phase1_demo/demo_topic/prompts/fig4_icon_status.yaml`:

```yaml
type: icon
subject: "status_pass — green check icon"
shape: "circle"
glyph: "✓"
glyph_color: "#ffffff"
fill_color: "#10b981"
style:
  size_px: 256
  padding_px: 24
  background: "transparent"
constraints:
  - "single glyph"
  - "high contrast glyph vs fill"
output:
  size: "256x256"
  quality: "n/a"
backend: cairo
```

- [ ] **Step 2: Write demo_topic.md**

Create `InterSubMod/examples/phase1_demo/demo_topic/demo_topic.md`:

```markdown
---
title: "Demo topic — Phase 1 closed-loop validation"
date: 2026-05-10
status: phase1_demo
authors: [Claude Opus 4.7 (assistant)]
spec: InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md
---

# Demo Topic — Phase 1 Closed-Loop Validation

> 此 .md + 同名主題資料夾用於驗證 Phase 1 三件事：
> 1. image-gen skill 能跑通 AI 軌 (concept + data) 與 cairo 軌 (flow + icon)
> 2. image-vision-check 能對 4 張圖各打 6 維分數
> 3. 手寫極簡 index.html 能 inline 4 張 base64 PNG 並正確顯示
> 此檔不是 Phase 2 html-preview skill 的設計範例，僅閉環驗證用。

## 範例圖示需求 (4 類各 1 張)

<!-- figure-needed: concept_diagram, slug=fig1_concept_loh -->
**Figure 1**: LOH 區段如何引發 phasing 失效的概念示意。

<!-- figure-needed: flow_diagram, slug=fig2_flow_pipeline -->
**Figure 2**: Self-Phasing V5 G1-G3 簡化 pipeline 流程圖。

<!-- figure-needed: data_mockup, slug=fig3_data_mockup -->
**Figure 3**: 7 樣本 baseline AUC 分布的概念型 bar chart (mock-up，非真數據)。

<!-- figure-needed: icon, slug=fig4_icon_status -->
**Figure 4**: status_pass 綠色勾勾 icon。

## 預期結果

執行：
```bash
# 1. 生圖
python3 InterSubMod/.claude/skills/image-gen/tools/dispatch.py \
  InterSubMod/examples/phase1_demo/demo_topic/prompts \
  InterSubMod/examples/phase1_demo/demo_topic/figures

# 2. 視覺檢核 (instruction phase)
python3 InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py \
  InterSubMod/examples/phase1_demo/demo_topic/figures

# (Claude agent 處理每張圖、寫 _pending_entries/*.json)

# 3. 視覺檢核 (aggregate)
python3 InterSubMod/.claude/skills/image-vision-check/tools/vision_check_runner.py \
  InterSubMod/examples/phase1_demo/demo_topic/figures --aggregate

# 4. 開瀏覽器
xdg-open InterSubMod/examples/phase1_demo/demo_topic/index.html
```

開瀏覽器後預期：4 張圖 inline 顯示、每張下方有 vision check verdict + score。
```

- [ ] **Step 3: Write index.html (極簡，~50 lines CSS, no Tailwind)**

Create `InterSubMod/examples/phase1_demo/demo_topic/index.html`:

```html
<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Demo Topic — Phase 1 Closed-Loop Validation</title>
<style>
  * { box-sizing: border-box; }
  body {
    font-family: system-ui, -apple-system, "Noto Sans CJK TC", sans-serif;
    max-width: 900px;
    margin: 2rem auto;
    padding: 0 1.5rem;
    color: #1f2937;
    line-height: 1.6;
    background: #ffffff;
  }
  h1 { border-bottom: 2px solid #1f2937; padding-bottom: 0.5rem; }
  h2 { margin-top: 2.5rem; color: #374151; }
  .figure {
    margin: 1.5rem 0 2.5rem;
    padding: 1rem;
    border: 1px solid #e5e7eb;
    border-radius: 6px;
  }
  .figure img {
    display: block;
    max-width: 100%;
    height: auto;
    margin: 0 auto;
  }
  .caption {
    margin-top: 0.75rem;
    font-size: 0.9rem;
    color: #6b7280;
  }
  .verdict {
    display: inline-block;
    padding: 0.2rem 0.6rem;
    border-radius: 4px;
    font-weight: 600;
    font-size: 0.85rem;
  }
  .verdict-pass { background: #d1fae5; color: #065f46; }
  .verdict-partial { background: #fef3c7; color: #92400e; }
  .verdict-fail { background: #fee2e2; color: #991b1b; }
  .verdict-pending { background: #e5e7eb; color: #374151; }
  @media print {
    body { max-width: none; }
    .figure { break-inside: avoid; }
  }
</style>
</head>
<body>

<h1>Demo Topic — Phase 1 Closed-Loop Validation</h1>

<p>此 HTML 為 <strong>手寫極簡 demo</strong>（無 Tailwind，~50 lines CSS），目的：驗證
Base64 inline 顯示 + 主題資料夾結構 + vision check verdict 呈現。Phase 2 的
<code>html-preview</code> skill 才會用 Tailwind L2-bake 升級。</p>

<p>來源 .md：<a href="./demo_topic.md"><code>demo_topic.md</code></a> ·
README：<a href="./README.md"><code>README.md</code></a></p>

<h2>Figure 1 · concept_diagram (AI 軌)</h2>
<div class="figure">
  <!-- After image-gen + Base64 encode, replace src with: data:image/png;base64,... -->
  <img src="figures/fig1_concept_loh.png" alt="LOH region triggers phasing degradation schematic">
  <div class="caption">
    LOH 區段如何引發 phasing 失效。
    <span class="verdict verdict-pending">vision check: pending</span>
  </div>
</div>

<h2>Figure 2 · flow_diagram (cairo 軌)</h2>
<div class="figure">
  <img src="figures/fig2_flow_pipeline.png" alt="Self-Phasing V5 G1-G3 pipeline">
  <div class="caption">
    G1 → G2 → G3 簡化流程。Cairo 程式繪圖、不會錯字。
    <span class="verdict verdict-pending">vision check: pending</span>
  </div>
</div>

<h2>Figure 3 · data_mockup (AI 軌)</h2>
<div class="figure">
  <img src="figures/fig3_data_mockup.png" alt="Mock-up of 7 sample baseline AUC distribution">
  <div class="caption">
    7 樣本 baseline AUC 概念型 bar chart（mock-up，非真數據）。
    <span class="verdict verdict-pending">vision check: pending</span>
  </div>
</div>

<h2>Figure 4 · icon (cairo 軌)</h2>
<div class="figure">
  <img src="figures/fig4_icon_status.png" alt="status_pass green check icon"
       style="max-width: 128px;">
  <div class="caption">
    status_pass icon。
    <span class="verdict verdict-pending">vision check: pending</span>
  </div>
</div>

<hr>
<p style="color: #6b7280; font-size: 0.85rem;">
  Generated by hand on 2026-05-10 as Phase 1 closed-loop demo.
  After running image-gen + image-vision-check, manually replace
  <code>src="figures/...png"</code> with <code>src="data:image/png;base64,..."</code>
  using <code>tools/pil_postprocess.py base64</code>, and update each
  <code>verdict-pending</code> span to the actual verdict from
  <code>figures/quality.json</code>.
</p>

</body>
</html>
```

- [ ] **Step 4: Write README.md (folder index)**

Create `InterSubMod/examples/phase1_demo/demo_topic/README.md`:

```markdown
# demo_topic — Phase 1 closed-loop demo

> Hand-written demo. NOT auto-generated by html-preview skill (Phase 2).
> Purpose: validate image-gen + image-vision-check + Base64 inline + topic folder structure.

## Contents

- [`demo_topic.md`](./demo_topic.md) — source markdown describing the 4 figure needs
- [`index.html`](./index.html) — minimal HTML viewer with 4 inline figures + verdict slots
- [`prompts/`](./prompts/) — 4 YAML prompt files (one per figure)
- [`figures/`](./figures/) — generated PNGs (`.gitignore`'d)

## How to regenerate (Phase 1 closed loop)

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod

# 1. preflight
bash .claude/skills/image-gen/tools/preflight.sh

# 2. generate 4 figures (2 AI via codex, 2 cairo via Python)
python3 .claude/skills/image-gen/tools/dispatch.py \
  examples/phase1_demo/demo_topic/prompts \
  examples/phase1_demo/demo_topic/figures

# 3. vision check (write instructions)
python3 .claude/skills/image-vision-check/tools/vision_check_runner.py \
  examples/phase1_demo/demo_topic/figures

# 4. Claude agent reads instructions, processes each PNG, writes _pending_entries/*.json
#    (this step is manual / agent-driven)

# 5. aggregate
python3 .claude/skills/image-vision-check/tools/vision_check_runner.py \
  examples/phase1_demo/demo_topic/figures --aggregate

# 6. open in browser
xdg-open examples/phase1_demo/demo_topic/index.html
# or copy to local machine for viewing
```

## Acceptance criteria (PI confirms)

- [ ] 4 PNGs generated successfully (no errors)
- [ ] 2 AI-generated images look like the prompts (concept + data_mockup)
- [ ] 2 cairo-generated images render exactly per YAML (flow + icon)
- [ ] vision check 4/4 pass OR partial with sensible suggestions
- [ ] index.html opens in browser, all 4 figures display inline
- [ ] no broken images, no console errors
- [ ] verdict spans show actual vision check result
```

- [ ] **Step 5: Commit (skeleton, before actually running)**

```bash
git add examples/phase1_demo/demo_topic/README.md \
        examples/phase1_demo/demo_topic/demo_topic.md \
        examples/phase1_demo/demo_topic/index.html \
        examples/phase1_demo/demo_topic/prompts/
git commit -m "feat(examples): phase1_demo/demo_topic — closed-loop validation skeleton

Hand-written demo skeleton. Figures will be generated by running:
  python3 .claude/skills/image-gen/tools/dispatch.py prompts/ figures/

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 14: Phase 1 acceptance smoke test + commit + memory update

**Files:**
- Run all preceding tooling end-to-end
- Modify: `InterSubMod/docs/CURRENT_FOCUS.md` (or equivalent active focus tracker)
- Add memory entry via Claude memory system

- [ ] **Step 1: Run preflight (must succeed before this task)**

```bash
bash .claude/skills/image-gen/tools/preflight.sh
```
Expected: 4 `[OK]` lines, exit 0.

- [ ] **Step 2: Generate 4 figures (this consumes ~60k ChatGPT tokens for 2 AI prompts)**

```bash
python3 .claude/skills/image-gen/tools/dispatch.py \
  examples/phase1_demo/demo_topic/prompts \
  examples/phase1_demo/demo_topic/figures
```
Expected:
- Cost preview line printed
- 4 lines `[ai] ...` or `[cairo] ...`
- Final line: `[image-gen] 4 prompts succeeded`
- 4 PNGs in `examples/phase1_demo/demo_topic/figures/`

If AI track fails (quota / network), Phase 1 partial demo still valid with cairo-only — note in commit.

- [ ] **Step 3: Run vision-check instruction phase**

```bash
python3 .claude/skills/image-vision-check/tools/vision_check_runner.py \
  examples/phase1_demo/demo_topic/figures
```
Expected: `_check_instructions.json` written.

- [ ] **Step 4: Manually drive Claude vision check on each image**

For each PNG, the Claude agent (in this session, or a sub-agent) should:
1. Read `_check_instructions.json`.
2. For each `images[i]`, use `Read` tool on `images[i].png`.
3. Read `images[i].yaml` for context (subject, type, constraints).
4. Read `checklists/{type}_check.md` for that image's type.
5. Read `prompts/vision_review_master.md` for output format.
6. Produce JSON conforming to schema, write to `_pending_entries/{basename}.json`.

(In this plan task, the executor agent does this step inline.)

- [ ] **Step 5: Aggregate**

```bash
python3 .claude/skills/image-vision-check/tools/vision_check_runner.py \
  examples/phase1_demo/demo_topic/figures --aggregate
```
Expected: `quality.json` with `summary` showing `n_pass / n_partial / n_fail`.

- [ ] **Step 6: Update index.html with actual verdicts and Base64 inline**

For each of 4 figures:
- Replace `src="figures/<basename>.png"` with `src="<output of pil_postprocess.py base64 figures/<basename>.png>"`
- Replace `<span class="verdict verdict-pending">vision check: pending</span>` with the actual verdict from quality.json (`pass` / `partial` / `fail` and the score string).

```bash
# helper: get base64 for one image
python3 .claude/skills/image-gen/tools/pil_postprocess.py base64 \
  examples/phase1_demo/demo_topic/figures/fig1_concept_loh.png | head -c 80
```

- [ ] **Step 7: Open index.html and visually verify**

If running locally:
```bash
xdg-open examples/phase1_demo/demo_topic/index.html
```

If on remote SSH server, copy file to local for viewing:
```bash
# on local machine
scp <user>@<server>:/big7_disk/liaoyoyo2001/InterSubMod/examples/phase1_demo/demo_topic/index.html /tmp/
open /tmp/index.html
```

PI checks:
- [ ] All 4 images visible inline (Base64 worked)
- [ ] Verdicts show real values (not "pending")
- [ ] Print preview looks reasonable (Ctrl+P)
- [ ] No console errors (open browser dev tools)

- [ ] **Step 8: Commit final demo state**

```bash
git add examples/phase1_demo/demo_topic/index.html \
        examples/phase1_demo/demo_topic/figures/quality.json
git commit -m "demo(phase1): closed-loop validation complete — 4 figures + vision check + inline HTML

PI acceptance: $(date +%Y-%m-%d).
Verdicts: see examples/phase1_demo/demo_topic/figures/quality.json

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

- [ ] **Step 9: Update CURRENT_FOCUS.md**

Append to `InterSubMod/docs/CURRENT_FOCUS.md` under appropriate section:

```markdown
## 2026-05-10 — Phase 1 完成

✅ image-gen + image-vision-check skills shipped, demo closed-loop validated at
`InterSubMod/examples/phase1_demo/demo_topic/`. Ready to plan Phase 2 (html-preview skill).

Spec: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md`
Phase 1 plan: `InterSubMod/docs/plans/2026/05/20260510_Phase1_image_gen_vision_check_implementation_01.md`
```

```bash
git add docs/CURRENT_FOCUS.md
git commit -m "docs(focus): mark Phase 1 (image-gen + vision-check) complete

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

- [ ] **Step 10: Add Claude memory entry**

Save a memory file for future sessions:

`/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/project_phase1_image_gen_skills_shipped.md`:

```markdown
---
name: Phase 1 image-gen + image-vision-check skills shipped
description: 2026-05-10 Phase 1 完成。image-gen 雙軌 (codex AI + cairo)、image-vision-check 6 維 + 4/6 通過。Demo at examples/phase1_demo/demo_topic/。
type: project
---

**2026-05-10 Phase 1 ship**

兩個新 skill 已上線：
- `InterSubMod/.claude/skills/image-gen/` — 雙軌：concept_diagram + data_mockup → AI (codex CLI `$imagegen` via ChatGPT OAuth, 免 USD); flow_diagram + icon → cairo (deterministic, 零成本)
- `InterSubMod/.claude/skills/image-vision-check/` — 6 維 checklist (D6) + 4/6 通過閾值 + 自動重生 1 次

**Why**: PI 痛點 — 圖示需求散落、無統一 prompt、無品質檢核、不想付 OpenAI API key 費用。
**How to apply**: 寫報告 .md 時加 `<!-- figure-needed: type, slug=foo -->` 標記，後續可呼叫 `/image-gen <prompts_dir> <figures_dir>` + `/image-vision-check <figures_dir>`。

**Phase 2 待做**: html-preview skill (主題資料夾化、Tailwind L2-bake、Anthropic 設計禁忌)。詳見 spec 的 §3.2。
**Phase 3 待做**: Tier A 6 個既有 skill 接入 html-preview companion。詳見 spec 的 §3.3。

**Spec**: `InterSubMod/docs/references/manual/20260510_HTML預覽_圖示生成_3skill_設計_01.md`
**Plan**: `InterSubMod/docs/plans/2026/05/20260510_Phase1_image_gen_vision_check_implementation_01.md`
```

Add pointer to MEMORY.md (under `## Active Research` or new `## Skills` section).

---

## Self-Review

Quick mental pass over this plan against the spec:

**1. Spec coverage** (D1–D14):
- ✓ D1 三 skill: image-gen (Tasks 4-9), image-vision-check (Tasks 10-12), html-preview (Phase 2, not in this plan)
- ✓ D2 對話觸發: SKILL.md description 中 USE WHEN 列觸發詞
- (Phase 2) D3 互動深度
- (Phase 2) D4 主題資料夾
- ✓ D5 4 類: concept/flow/data/icon templates (Task 3)
- ✓ D6 6 維 + 4/6 閾值: checklists + verdict rules (Tasks 10-11)
- ✓ D7 figures `.gitignore` (Task 1) + prompts in git
- ✓ D8 重生 1 次: apply_suggestions_to_yaml + retry_loop test (Task 12)
- ✓ D9 首批範圍: Phase 1 = image-gen + vision-check + demo (Tasks 1-14)
- (Phase 3) D10 Tier A 接入
- (Phase 3) D11 companion 不取代
- (Phase 2) D12 主題資料夾命名
- ✓ D13 codex OAuth-only: codex_image_gen.py (Task 4) + preflight check (Task 2)
- ✓ D14 雙軌 + PIL: cairo_render.py (Task 6) + pil_postprocess.py (Task 7) + dispatch (Task 8)

**2. Placeholder scan**: Reviewed — no TBD/TODO/handwave. All test code, all impl code, all commands shown.

**3. Type / signature consistency**: 
- `extract_prompt_text(yaml_path)` defined in Task 4, used (no other reference)
- `build_command(yaml_path, output_dir)` Task 4 → consumed by `route_one` Task 8 ✓
- `collect_latest(candidate_dirs, target_dir, *, target_name, since_epoch)` Task 5 → consumed by `route_one` Task 8 ✓
- `render(yaml_path, output)` Task 6 → consumed by `route_one` Task 8 ✓
- `resize_to_target(src, dst, target_size)` Task 7 → consumed by `dispatch.main` Task 8 + Task 14 inlining ✓
- `validate_entry(entry, schema_path)` Task 11 → consumed by `aggregate_pending` Task 11 ✓
- `apply_suggestions_to_yaml(yaml_path, suggestions)` Task 12 → not yet consumed by automation (D8 retry done by calling agent reading quality.json)

**4. Scope**: Phase 1 only. Phase 2 (html-preview) and Phase 3 (Tier A接入) explicitly deferred per spec D9 and the user's chosen rollout.

---

## Execution Handoff

Plan complete. Two execution options:

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — Execute tasks in this session using `superpowers:executing-plans`, batch execution with checkpoints for review.

Which approach?
