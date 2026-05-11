# ppt_toolkit

Reusable building blocks for the [`myPPT`](../../.claude/skills/myPPT/) skill.

`ppt_toolkit` is a small Python package that consolidates the helpers proven in
v2/v3 of the Self-Phasing 教授報告 build (`build_pptx.py`) and adds the new
helpers required by the myPPT playbook (`assertion_title`, `tier_aware_speaker_note`,
`claude_vision_review`, YAML-driven `style_library_loader`).

The toolkit is **framework-agnostic at the deck level** — every helper operates
on python-pptx `Slide` / `TextFrame` objects. It can be imported by any future
presentation script (weekly reports, conference talks, ad-hoc review decks).

---

## Install / Import

The package is intended to be imported directly from `InterSubMod/tools/`:

```python
import sys
from pathlib import Path

sys.path.insert(0, str(Path("/big7_disk/liaoyoyo2001/InterSubMod/tools").resolve()))

from ppt_toolkit import (
    add_text_with_fallback, fit_image_within, set_speaker_note,
    load_object, load_layout,
    assertion_title, set_tier_aware_note,
    build_review_template, discover_slide_pngs,
)
```

External dependency: `pptx`, `Pillow`, `lxml`, `pyyaml`.

---

## Module Map

| Module | Public symbols | Source |
|--------|----------------|--------|
| `text_helpers` | `add_text_with_fallback`, `fit_image_within`, `set_speaker_note`, `PALETTE`, `DEFAULT_LATIN_FONT`, `DEFAULT_CJK_FONT` | Extracted from `v2/scripts/build_pptx.py:109-271` |
| `style_library_loader` | `load_object`, `load_layout`, `load_palette`, `ObjectSpec`, `LayoutSpec` | New (myPPT plan §13, §17) |
| `assertion_title` | `assertion_title`, `is_assertion_sentence` | New (myPPT plan §5) |
| `tier_aware_speaker_note` | `set_tier_aware_note`, `build_tier_note`, `estimate_speaking_seconds`, `SpeakerNoteEstimate` | New (myPPT plan §4 + §6) |
| `claude_vision_review` | `discover_slide_pngs`, `build_review_template`, `build_summary_table`, `CHECKPOINTS` | New (myPPT plan §7 + §28) |

---

## Quick Reference

### 1. `add_text_with_fallback(text_frame, text, ...)`

Per-character Latin / CJK font fallback. Injects `<a:latin>` and `<a:ea>`
elements so PowerPoint paints ASCII glyphs in `latin_font` (default `Arial`)
and CJK glyphs in `cjk_font` (default `Droid Sans Fallback`).

```python
from pptx import Presentation
from pptx.util import Inches
from ppt_toolkit import add_text_with_fallback

prs = Presentation()
slide = prs.slides.add_slide(prs.slide_layouts[6])
box = slide.shapes.add_textbox(Inches(1), Inches(1), Inches(8), Inches(1)).text_frame
add_text_with_fallback(box, "V5 把 17.3:1 降回 1:1", font_size=28, bold=True)
```

### 2. `fit_image_within(slide, path, x, y, max_w, max_h)`

Equal-ratio scaling — never forces both width AND height (avoids the squashed-
IGV failure mode). Renders a light-grey fallback rectangle with the filename
when the source PNG is missing.

```python
from ppt_toolkit import fit_image_within
fit_image_within(slide, "figures/v5max1.png",
                 Inches(1), Inches(1.5), Inches(11), Inches(5))
```

### 3. `set_speaker_note(slide, text, min_chars=350, max_chars=None)`

Mandatory bound check on speaker notes. v3 enhancement: also enforces an
**upper** bound — long notes (> 360 CJK chars ≈ 90 sec) overshoot the
PLOS 1 min/slide rule.

```python
from ppt_toolkit import set_speaker_note
set_speaker_note(slide, "...", min_chars=500, max_chars=900)
```

### 4. `assertion_title(slide, thesis_sentence, subtitle=None)`

Renders a title bar **and** rejects generic labels ("Results", "Discussion",
"結果"...) unless `allow_label=True`. Heuristic check: requires a verb / claim
marker or a number / comparison operator in the title.

```python
from ppt_toolkit import assertion_title
assertion_title(slide, "Phasing 沒問題，是 Tag 有問題")          # OK
assertion_title(slide, "結果")                                     # raises
assertion_title(slide, "Acknowledgments", allow_label=True)        # bypassed
```

### 5. `set_tier_aware_note(slide, tier2, tier3, target_seconds=75)`

Auto-splits Tier 2 (must-say) and Tier 3 (oral-optional, tagged
`[ORAL-OPTIONAL]`) and estimates speaking time. Returns a
`SpeakerNoteEstimate` for build reports.

```python
from ppt_toolkit import set_tier_aware_note

est = set_tier_aware_note(slide,
    tier2_text="V5 三層投票把 17.3:1 降回 1:1，AMB% −9.5pp。",
    tier3_text="守恆律 A 在 15/15 sites PASS；cherry-picked 樣本可信但需 F8 隨機抽樣補強。",
    target_seconds=75,
)
print(est)  # [OK] tier2=24c (~4s) tier3=58c (~9s) total=13s vs target=75s
```

### 6. `load_object` / `load_layout` (style_library)

```python
from ppt_toolkit import load_object, load_layout

caveat = load_object("caveat_red_strip")
caveat.render(slide, Inches(1), Inches(5.5),
              {"caveat_text": "Source mis-reported 38%, real 16.5%"})

layout = load_layout("before_after_split")
layout.populate(slide, {
    "title.text": "Baseline tag voting vs V5 三層投票",
    "left.heading": "Baseline (somatic-first)",
    "left.image": "figures/baseline_getvote.png",
    "right.heading": "V5 (germline-first + Layer 1.5)",
    "right.image": "figures/v5_getvote.png",
    "delta_strip.content": "Δ AMB 17.5% → 8.0%",
})
```

YAML lookup order:
1. `$MYPPT_STYLE_LIBRARY` env var (testing override)
2. `<repo>/.claude/skills/myPPT/style_library/{objects,layouts}/`

### 7. `build_review_template` (Claude Vision)

Generates a 10-checkpoint markdown template per slide PNG. Does **not** call
any Claude API — the parent conversation reads each PNG with the Read tool
and fills in PASS / FAIL / PARTIAL verdicts.

```python
from pathlib import Path
from ppt_toolkit import discover_slide_pngs, build_review_template

pngs = discover_slide_pngs("/path/to/v3/output/screenshots/")
md = build_review_template(pngs, deck_name="v3 Self-Phasing")
Path("audits/v3_vision_review_template.md").write_text(md, encoding="utf-8")
```

---

## Coding standards

- **Style**: Black-compatible (88 cols) + type hints + docstrings.
- **Comments**: English only.
- **Tests**: not bundled here; the parent skill exercises the toolkit via
  v3 `build_pptx.py` and `screenshot_all.py`.

## Relation to v2 build_pptx.py

The v2 file is **frozen** (it is the source-of-truth deck for the 2026-04-29
Self-Phasing presentation). `ppt_toolkit` is a *copy-style* extraction — v2
keeps its inline helpers, but new decks (v3 onward) should import from
`ppt_toolkit` instead of duplicating logic.
