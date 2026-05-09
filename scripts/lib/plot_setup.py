"""
plot_setup.py — Centralized matplotlib rcParams for CJK font fallback.

Source: prior art E-1 (Matplotlib >=3.6 per-glyph fallback).
Plan ref: ~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md §12.2 T1-4.
Lessons: feedback_matplotlib_cjk_font_rule.md (S20/S23 中文亂碼 incidents).

Usage:
    from scripts.lib.plot_setup import setup_plot_style
    setup_plot_style()
    # then proceed with normal matplotlib code; CJK characters render correctly.

Design notes:
    - Never raise on missing fonts (would break unrelated plot scripts).
    - rcParams font.sans-serif fallback chain: DejaVu Sans -> Noto Sans CJK TC ->
      Source Han Sans TW -> Droid Sans Fallback -> sans-serif.
    - matplotlib 3.6+ honors per-glyph fallback automatically (issue #20740 fixed).
    - Avoid LaTeX math + CJK in same figure (issue #29173 still open).
"""

from __future__ import annotations

import logging
from typing import Iterable

logger = logging.getLogger(__name__)

# Ordered preference: Latin first, then CJK fallbacks. Ordering matters: matplotlib
# walks the chain and picks the first font that contains each glyph.
DEFAULT_FONT_CHAIN = (
    "DejaVu Sans",
    "Noto Sans CJK TC",
    "Noto Sans CJK SC",
    "Source Han Sans TW",
    "Source Han Sans CN",
    "WenQuanYi Zen Hei",
    "Droid Sans Fallback",
    "sans-serif",
)


def _available_fonts() -> set[str]:
    """Return set of font family names matplotlib can find."""
    try:
        from matplotlib import font_manager
    except ImportError:
        logger.warning("matplotlib not installed; plot_setup is a no-op")
        return set()
    return {f.name for f in font_manager.fontManager.ttflist}


def setup_plot_style(
    font_chain: Iterable[str] = DEFAULT_FONT_CHAIN,
    *,
    base_size: int = 11,
    dpi: int = 120,
    enable_unicode_minus: bool = False,
) -> dict:
    """
    Apply central rcParams for InterSubMod plots.

    Returns dict of {applied_fonts, missing_fonts} for diagnostics.
    Never raises — missing fonts only emit a warning so unrelated plot scripts
    are not broken by font availability.
    """
    try:
        import matplotlib
    except ImportError:
        logger.warning("matplotlib unavailable; setup_plot_style is a no-op")
        return {"applied_fonts": [], "missing_fonts": list(font_chain)}

    available = _available_fonts()
    requested = list(font_chain)
    applied = [f for f in requested if f in available or f == "sans-serif"]
    missing = [f for f in requested if f not in available and f != "sans-serif"]

    if not applied:
        # Extreme fallback: keep matplotlib defaults; warn loudly.
        logger.warning(
            "No requested fonts available (%s). Falling back to matplotlib default; "
            "CJK characters may render as boxes.",
            requested,
        )
        applied = ["sans-serif"]

    if missing:
        logger.info(
            "plot_setup: requested fonts not found: %s. Using available chain: %s",
            missing,
            applied,
        )

    # Surface a louder warning when all CJK fonts are missing — figures with
    # Chinese text will render glyph boxes otherwise. This is an env gap, not
    # a code defect; the suggested install commands resolve it on Debian/Ubuntu.
    cjk_keywords = ("CJK", "Han", "WenQuanYi", "Droid Sans Fallback")
    cjk_applied = [f for f in applied if any(k in f for k in cjk_keywords)]
    if not cjk_applied:
        logger.warning(
            "plot_setup: NO CJK font available in current environment. "
            "Chinese characters will render as boxes. "
            "Install one of: `sudo apt-get install fonts-noto-cjk` "
            "or `sudo apt-get install fonts-wqy-zenhei` "
            "or via conda: `conda install -c conda-forge font-ttf-noto-sans-cjk-jp` (and -tc/-sc variants)."
        )

    matplotlib.rcParams["font.family"] = ["sans-serif"]
    matplotlib.rcParams["font.sans-serif"] = applied
    matplotlib.rcParams["font.size"] = base_size
    matplotlib.rcParams["figure.dpi"] = dpi
    matplotlib.rcParams["savefig.dpi"] = dpi
    matplotlib.rcParams["axes.unicode_minus"] = enable_unicode_minus

    return {"applied_fonts": applied, "missing_fonts": missing}


__all__ = ["setup_plot_style", "DEFAULT_FONT_CHAIN"]
