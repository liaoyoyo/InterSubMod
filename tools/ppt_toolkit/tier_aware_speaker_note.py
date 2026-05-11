"""tier_aware_speaker_note.py — Auto-tier-split speaker notes (myPPT §4 + §6).

Every slide carries a speaker note split into:
    Tier 2 (must-say)        — core narration the speaker cannot drop.
    Tier 3 (oral-optional)   — depth detail tagged ``[ORAL-OPTIONAL]``,
                               speaker decides at runtime whether to deliver.

This helper:
    1. Concatenates Tier 2 + ``[ORAL-OPTIONAL]`` block + Tier 3.
    2. Estimates speaking time using language-aware char-per-min rates
       (Mandarin 400 chars/min, English 150 wpm ≈ 750 chars/min).
    3. Surfaces a warning when Tier 2 alone exceeds ``target_seconds``
       (PLOS Rule: 1 min/slide).

Usage:

    >>> from ppt_toolkit.tier_aware_speaker_note import set_tier_aware_note
    >>> set_tier_aware_note(slide,
    ...     tier2_text="V5 三層投票把 17.3:1 降回 1:1，AMB% −9.5pp。",
    ...     tier3_text="守恆律 A 在 15/15 sites PASS，cherry-picked 樣本...",
    ...     target_seconds=75,
    ... )
    SpeakerNoteEstimate(tier2_seconds=12.0, tier3_seconds=45.0,
                       total_seconds=57.0, over_target=False, ...)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from ppt_toolkit.text_helpers import set_speaker_note

# Language-aware speaking rate constants (chars per second).
# Mandarin: ~400 字/min ≈ 6.67 chars/sec (conservative; rapid speakers reach 8).
# English:  ~150 wpm ≈ 750 chars/min ≈ 12.5 chars/sec (1 word ≈ 5 chars).
CN_CHARS_PER_SEC = 6.67
EN_CHARS_PER_SEC = 12.5

ORAL_OPTIONAL_DELIMITER = "\n\n[ORAL-OPTIONAL]\n"


@dataclass
class SpeakerNoteEstimate:
    """Result of tier-aware speaker note construction."""

    tier2_chars: int
    tier3_chars: int
    tier2_seconds: float
    tier3_seconds: float
    total_seconds: float
    over_target: bool
    target_seconds: float
    warnings: List[str] = field(default_factory=list)

    @property
    def total_chars(self) -> int:
        return self.tier2_chars + self.tier3_chars

    def __str__(self) -> str:
        flag = "OVER" if self.over_target else "OK"
        return (
            f"[{flag}] tier2={self.tier2_chars}c (~{self.tier2_seconds:.0f}s) "
            f"tier3={self.tier3_chars}c (~{self.tier3_seconds:.0f}s) "
            f"total={self.total_seconds:.0f}s vs target={self.target_seconds:.0f}s"
        )


def _is_cjk(ch: str) -> bool:
    cp = ord(ch)
    # Subset of full CJK_RANGES sufficient for chars/sec estimation.
    return (
        0x3000 <= cp <= 0x303F
        or 0x3040 <= cp <= 0x30FF
        or 0x3400 <= cp <= 0x4DBF
        or 0x4E00 <= cp <= 0x9FFF
        or 0xF900 <= cp <= 0xFAFF
        or 0xFF00 <= cp <= 0xFFEF
    )


def estimate_speaking_seconds(text: str) -> float:
    """Estimate seconds to speak *text* using mixed CN/EN rates.

    CJK characters use ``CN_CHARS_PER_SEC``; everything else (including
    punctuation, digits, ASCII letters) uses ``EN_CHARS_PER_SEC``. Whitespace
    is excluded from the count to avoid double-charging spaces between words.
    """
    if not text:
        return 0.0
    cn_chars = 0
    en_chars = 0
    for ch in text:
        if ch.isspace():
            continue
        if _is_cjk(ch):
            cn_chars += 1
        else:
            en_chars += 1
    cn_sec = cn_chars / CN_CHARS_PER_SEC if cn_chars else 0.0
    en_sec = en_chars / EN_CHARS_PER_SEC if en_chars else 0.0
    return cn_sec + en_sec


def build_tier_note(
    tier2_text: str,
    tier3_text: Optional[str] = None,
    target_seconds: float = 75.0,
    min_chars: int = 350,
) -> tuple[str, SpeakerNoteEstimate]:
    """Build the combined note text + estimate without writing it.

    Useful for dry-run reporting (e.g. "scripts/screenshot_all.py" can call
    this to print a per-slide budget table before final commit).

    Returns
    -------
    (combined_text, estimate)
    """
    tier2_text = (tier2_text or "").strip()
    tier3_text = (tier3_text or "").strip()

    if not tier2_text:
        raise ValueError("tier2_text (must-say content) is required and was empty.")

    if tier3_text:
        combined = tier2_text + ORAL_OPTIONAL_DELIMITER + tier3_text
    else:
        combined = tier2_text

    tier2_seconds = estimate_speaking_seconds(tier2_text)
    tier3_seconds = estimate_speaking_seconds(tier3_text)
    total_seconds = tier2_seconds + tier3_seconds
    warnings: List[str] = []

    if tier2_seconds > target_seconds:
        warnings.append(
            f"Tier 2 alone exceeds target ({tier2_seconds:.0f}s > {target_seconds:.0f}s); "
            f"split must-say content across two slides."
        )
    if total_seconds > target_seconds * 1.5:
        warnings.append(
            f"Total speaker note ({total_seconds:.0f}s) is >1.5x target; "
            f"demote some Tier 3 content to backup slide / Q&A prep."
        )
    if len(combined) < min_chars:
        warnings.append(
            f"Combined note has {len(combined)} chars, below min_chars={min_chars}; "
            f"add Tier 3 oral-optional detail to satisfy storyboard requirement."
        )

    estimate = SpeakerNoteEstimate(
        tier2_chars=len(tier2_text),
        tier3_chars=len(tier3_text),
        tier2_seconds=tier2_seconds,
        tier3_seconds=tier3_seconds,
        total_seconds=total_seconds,
        over_target=tier2_seconds > target_seconds,
        target_seconds=target_seconds,
        warnings=warnings,
    )
    return combined, estimate


def set_tier_aware_note(
    slide,
    tier2_text: str,
    tier3_text: Optional[str] = None,
    target_seconds: float = 75.0,
    min_chars: int = 350,
    max_chars: Optional[int] = None,
    raise_on_over_target: bool = False,
) -> SpeakerNoteEstimate:
    """Write a tier-split speaker note onto *slide*.

    Parameters
    ----------
    slide : pptx.slide.Slide
    tier2_text : str
        Must-say narrative.
    tier3_text : str, optional
        Oral-optional depth content; appended after ``[ORAL-OPTIONAL]`` marker.
    target_seconds : float
        Per-slide speaking budget (PLOS rule recommends 60s; v3 storyboard
        uses 75-90s for high-density slides).
    min_chars, max_chars : passed through to :func:`set_speaker_note`.
    raise_on_over_target : bool
        When True, raises ValueError if Tier 2 alone exceeds target. Default
        False — emits warnings via the SpeakerNoteEstimate.warnings list,
        letting the caller decide.

    Returns
    -------
    SpeakerNoteEstimate
        Diagnostic estimate (also returned for inclusion in build reports).
    """
    combined, estimate = build_tier_note(
        tier2_text=tier2_text,
        tier3_text=tier3_text,
        target_seconds=target_seconds,
        min_chars=min_chars,
    )
    if raise_on_over_target and estimate.over_target:
        raise ValueError(
            f"Slide note Tier 2 over target: {estimate}. Set raise_on_over_target=False "
            f"to allow with warnings, or split must-say content."
        )
    set_speaker_note(slide, combined, min_chars=min_chars, max_chars=max_chars)
    return estimate
