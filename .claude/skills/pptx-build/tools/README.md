# myPPT — Tools Index

This page indexes the Python helpers used by the myPPT skill. The actual
implementation lives in `InterSubMod/tools/ppt_toolkit/` (so it can be
reused outside the skill, e.g. by future presentation scripts in
`docs/presentations/...`).

For full API documentation see the package README:
[`InterSubMod/tools/ppt_toolkit/README.md`](../../../../tools/ppt_toolkit/README.md)

---

## Module Map

| Module | Purpose | Key API |
|--------|---------|---------|
| `text_helpers` | Per-character Latin/CJK font fallback, image fitter, speaker-note bound checker | `add_text_with_fallback`, `fit_image_within`, `set_speaker_note` |
| `style_library_loader` | YAML `objects/*.yaml` and `layouts/*.yaml` deserialization | `load_object`, `load_layout`, `load_palette`, `ObjectSpec`, `LayoutSpec` |
| `assertion_title` | Thesis-style title bar with generic-label rejection (PLOS Rule 1 + Alley AE) | `assertion_title`, `is_assertion_sentence` |
| `tier_aware_speaker_note` | Tier 2 / Tier 3 split + speaking-time estimator (CN 400 字/min, EN 150 wpm) | `set_tier_aware_note`, `build_tier_note`, `estimate_speaking_seconds` |
| `claude_vision_review` | 10-checkpoint markdown template for Claude Vision review (no API call) | `discover_slide_pngs`, `build_review_template`, `build_summary_table` |

---

## Quick Start

```python
import sys
from pathlib import Path

# Add tools/ to import path. In production, install ppt_toolkit as a package.
sys.path.insert(0, str(Path("/big7_disk/liaoyoyo2001/InterSubMod/tools").resolve()))

from pptx import Presentation
from pptx.util import Inches

from ppt_toolkit import (
    add_text_with_fallback, fit_image_within,
    assertion_title, set_tier_aware_note,
    load_object, load_layout,
    discover_slide_pngs, build_review_template,
)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)

slide = prs.slides.add_slide(prs.slide_layouts[6])

# 1. Thesis-style title (rejects generic labels)
assertion_title(slide, "V5 三層投票把 17.3:1 降回 1:1",
                subtitle="clean PS blocks, 全基因組")

# 2. Render a layout from YAML
layout = load_layout("before_after_split")
layout.populate(slide, {
    "title.text": "Baseline vs V5",
    "left.heading": "Baseline (somatic-first)",
    "left.image": "figures/baseline.png",
    "right.heading": "V5 (germline-first + Layer 1.5)",
    "right.image": "figures/v5.png",
    "delta_strip.content": "Δ AMB 17.5% → 8.0%",
})

# 3. Tier-split speaker note with timing budget
estimate = set_tier_aware_note(slide,
    tier2_text="V5 三層投票把 17.3:1 降回 1:1，AMB% −9.5pp。",
    tier3_text="守恆律 A 在 15 sites 全 PASS；cherry-picked 需 F8 隨機抽樣補強。",
    target_seconds=75,
)
print(estimate)
```

---

## Versioning

| Version | Date | Notes |
|---------|------|-------|
| 0.1.0 | 2026-04-30 | Initial extraction from v2 build_pptx.py + new myPPT helpers |

---

## Source Provenance

- `add_text_with_fallback`, `fit_image_within`, `set_speaker_note`
  — extracted from `docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/scripts/build_pptx.py:109-271`.
- `assertion_title`, `set_tier_aware_note`, `style_library_loader`,
  `claude_vision_review`
  — new modules per myPPT plan §13, §17, §22.

The v2 file is **frozen** (the source-of-truth deck for the 2026-04-29
Self-Phasing presentation). `ppt_toolkit` is a *copy-style* extraction —
v2 keeps its inline helpers; new decks (v3 onward) should `import ppt_toolkit`.
