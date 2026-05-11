"""ppt_toolkit — Reusable building blocks for the myPPT skill.

Re-exports the public API so callers can write::

    from ppt_toolkit import (
        add_text_with_fallback, fit_image_within, set_speaker_note,
        load_object, load_layout, assertion_title,
        set_tier_aware_note, build_review_template,
    )

The toolkit is intentionally framework-agnostic at the deck level — every
helper operates on python-pptx slide / text-frame objects so it can be
reused by future presentations (weekly reports, conference talks, etc.).

Module map:
    - text_helpers              — per-character font fallback, image fitter,
                                  speaker-note bound checker.
    - style_library_loader      — YAML object/layout deserialization.
    - assertion_title           — thesis-style title bar guard.
    - tier_aware_speaker_note   — Tier 2/3 split + speaking-time estimator.
    - claude_vision_review      — 10-checkpoint markdown template generator.
"""

from ppt_toolkit.text_helpers import (
    DEFAULT_CJK_FONT,
    DEFAULT_CJK_FONT_PATH,
    DEFAULT_LATIN_FONT,
    PALETTE,
    add_text_with_fallback,
    fit_image_within,
    set_speaker_note,
)
from ppt_toolkit.style_library_loader import (
    LayoutSpec,
    ObjectSpec,
    clear_cache,
    load_layout,
    load_object,
    load_palette,
)
from ppt_toolkit.assertion_title import (
    assertion_title,
    is_assertion_sentence,
)
from ppt_toolkit.tier_aware_speaker_note import (
    SpeakerNoteEstimate,
    build_tier_note,
    estimate_speaking_seconds,
    set_tier_aware_note,
)
from ppt_toolkit.claude_vision_review import (
    CHECKPOINTS,
    CheckpointResult,
    SlidePng,
    SlideReview,
    build_review_template,
    build_summary_table,
    discover_slide_pngs,
    make_default_checkpoints,
)

__version__ = "0.1.0"

__all__ = [
    # text_helpers
    "DEFAULT_CJK_FONT",
    "DEFAULT_CJK_FONT_PATH",
    "DEFAULT_LATIN_FONT",
    "PALETTE",
    "add_text_with_fallback",
    "fit_image_within",
    "set_speaker_note",
    # style_library_loader
    "LayoutSpec",
    "ObjectSpec",
    "clear_cache",
    "load_layout",
    "load_object",
    "load_palette",
    # assertion_title
    "assertion_title",
    "is_assertion_sentence",
    # tier_aware_speaker_note
    "SpeakerNoteEstimate",
    "build_tier_note",
    "estimate_speaking_seconds",
    "set_tier_aware_note",
    # claude_vision_review
    "CHECKPOINTS",
    "CheckpointResult",
    "SlidePng",
    "SlideReview",
    "build_review_template",
    "build_summary_table",
    "discover_slide_pngs",
    "make_default_checkpoints",
]
