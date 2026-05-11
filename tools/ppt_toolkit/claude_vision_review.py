"""claude_vision_review.py — Module-level helper for Claude Vision review.

This module does **not** call any Claude API directly. Per the myPPT plan
(§7, §28), the visual review is performed by the parent conversation that
already has Claude Vision read access. This helper only:

    1. Discovers slide PNGs produced by ``scripts/screenshot_all.py``.
    2. Emits a 10-checkpoint markdown template per slide that the parent
       conversation can fill in by reading the PNGs with the Read tool.
    3. Produces a per-deck summary table for inclusion in audit reports.

The 10 checkpoints follow myPPT playbook §5:
    1. Heading is an assertion sentence (not generic label).
    2. One idea per slide.
    3. ≤ 6 elements.
    4. Text density (EN ≤30 words / CN ≤60 chars, excluding title).
    5. Distracted takeaway: main point clear without speaker audio.
    6. Visual contrast (colour, font, size).
    7. Figure-to-text ratio ≥ 60% visual.
    8. Citation / data source visible.
    9. Speaker note ≤ ~75s estimated.
   10. Disaster fallback (PDF readable, no animation dependency).

Usage:

    >>> from ppt_toolkit.claude_vision_review import (
    ...     discover_slide_pngs, build_review_template, build_summary_table,
    ... )
    >>> pngs = discover_slide_pngs("/path/to/v3/output/screenshots/")
    >>> md = build_review_template(pngs, deck_name="v3 Self-Phasing")
    >>> Path("audits/v3_vision_review_template.md").write_text(md)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

CHECKPOINTS = [
    ("CP1", "Heading is an assertion sentence (not generic label)"),
    ("CP2", "One idea per slide"),
    ("CP3", "≤ 6 visual elements"),
    ("CP4", "Text density: EN ≤ 30 words / CN ≤ 60 chars (excl. title)"),
    ("CP5", "Distracted takeaway: main point clear without audio"),
    ("CP6", "Visual contrast (colour / font / size hierarchy)"),
    ("CP7", "Figure-to-text ratio ≥ 60% visual"),
    ("CP8", "Citation / data source visible"),
    ("CP9", "Speaker note ≤ ~75s (estimated chars)"),
    ("CP10", "Disaster fallback: PDF readable, no animation dependency"),
]

_SLIDE_PNG_RE = re.compile(r"slide[_-]?(\d+)", re.IGNORECASE)


@dataclass
class SlidePng:
    """Discovered slide screenshot."""

    slide_index: int
    path: Path
    title_hint: Optional[str] = None  # filled by caller if storyboard mapping known


@dataclass
class CheckpointResult:
    """Single checkpoint result for a slide (filled by parent conversation)."""

    cp_id: str
    description: str
    verdict: str = "PENDING"  # PASS / FAIL / PARTIAL / PENDING
    notes: str = ""


@dataclass
class SlideReview:
    """Per-slide vision review result."""

    slide_index: int
    title_hint: Optional[str] = None
    checkpoints: List[CheckpointResult] = field(default_factory=list)
    overall_verdict: str = "PENDING"

    def pass_count(self) -> int:
        return sum(1 for cp in self.checkpoints if cp.verdict == "PASS")

    def fail_count(self) -> int:
        return sum(1 for cp in self.checkpoints if cp.verdict == "FAIL")

    def partial_count(self) -> int:
        return sum(1 for cp in self.checkpoints if cp.verdict == "PARTIAL")


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------
def discover_slide_pngs(directory, pattern: str = "*.png") -> List[SlidePng]:
    """Walk *directory* and return SlidePng entries sorted by slide_index.

    Matches filenames containing ``slide<N>`` / ``slide_<N>`` / ``slide-<N>``
    case-insensitively. Files without a slide number are skipped.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"Screenshot directory not found: {directory}")

    found: List[SlidePng] = []
    for png in sorted(directory.glob(pattern)):
        m = _SLIDE_PNG_RE.search(png.stem)
        if not m:
            continue
        idx = int(m.group(1))
        found.append(SlidePng(slide_index=idx, path=png))
    found.sort(key=lambda s: s.slide_index)
    return found


# ---------------------------------------------------------------------------
# Template generation
# ---------------------------------------------------------------------------
def build_review_template(
    slides: List[SlidePng],
    deck_name: str = "Deck",
    title_map: Optional[Dict[int, str]] = None,
) -> str:
    """Return a markdown template the parent conversation should fill in.

    Each slide section lists the 10 checkpoints with a PENDING verdict and
    an empty notes field. The parent conversation reads each PNG using the
    Read tool, evaluates the checkpoints, and replaces the PENDING markers.

    Parameters
    ----------
    slides : list of SlidePng
        Output of :func:`discover_slide_pngs`.
    deck_name : str
        Human-readable deck identifier for the markdown title.
    title_map : dict[int, str], optional
        Map of slide_index -> storyboard title hint.

    Returns
    -------
    str
        Markdown body ready to write to ``audits/<deck>_vision_review.md``.
    """
    title_map = title_map or {}
    lines: List[str] = []
    lines.append(f"# {deck_name} — Claude Vision 10-Checkpoint Review")
    lines.append("")
    lines.append(
        "> Replace each `PENDING` verdict with `PASS` / `FAIL` / `PARTIAL` after "
        "reading the slide PNG with the `Read` tool. Add 1-2 sentences in the "
        "`notes` column when verdict is `FAIL` or `PARTIAL`."
    )
    lines.append("")

    for slide in slides:
        title = title_map.get(slide.slide_index) or slide.title_hint or ""
        header = f"## Slide {slide.slide_index}"
        if title:
            header += f" — {title}"
        lines.append(header)
        lines.append("")
        lines.append(f"PNG: `{slide.path}`")
        lines.append("")
        lines.append("| CP | Check | Verdict | Notes |")
        lines.append("|----|-------|---------|-------|")
        for cp_id, desc in CHECKPOINTS:
            lines.append(f"| {cp_id} | {desc} | PENDING |  |")
        lines.append("")
        lines.append("**Overall verdict (PASS / PARTIAL / FAIL)**: PENDING")
        lines.append("")

    return "\n".join(lines)


def build_summary_table(reviews: List[SlideReview]) -> str:
    """Return a markdown summary table aggregating per-slide verdicts.

    Format::

        | Slide | Title | PASS | FAIL | PARTIAL | Overall |
    """
    lines: List[str] = []
    lines.append("| Slide | Title | PASS | FAIL | PARTIAL | Overall |")
    lines.append("|-------|-------|------|------|---------|---------|")
    for r in sorted(reviews, key=lambda x: x.slide_index):
        title = r.title_hint or ""
        lines.append(
            f"| S{r.slide_index} | {title} | "
            f"{r.pass_count()} | {r.fail_count()} | "
            f"{r.partial_count()} | {r.overall_verdict} |"
        )
    return "\n".join(lines)


def make_default_checkpoints() -> List[CheckpointResult]:
    """Return a fresh list of CheckpointResult entries (all PENDING)."""
    return [CheckpointResult(cp_id=cp, description=desc) for cp, desc in CHECKPOINTS]
