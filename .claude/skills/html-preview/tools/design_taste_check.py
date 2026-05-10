"""design_taste_check — regex-based static detection of Anthropic design taboos.

Per design doc D14 + D18 (evaluator-optimizer pattern). For deeper audit,
use Claude vision with `prompts/design_taste_check_prompt.md` separately.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import TypedDict


class TabooResult(TypedDict):
    taboos_detected: dict
    n_taboos: int
    verdict: str
    suggestions: list


_GRADIENT_RE = re.compile(r"(linear|radial)-gradient", re.IGNORECASE)
_GLASS_RE = re.compile(r"backdrop-filter\s*:|backdrop-blur", re.IGNORECASE)
_INDIGO_HEXES = ("#4F46E5", "#6366F1", "#818CF8", "#A5B4FC", "#C7D2FE")
_EMOJI_HEADER_RE = re.compile(
    r"<h[1-6][^>]*>[^<]*[\U0001F300-\U0001FAFF\U00002600-\U000027BF]",
    re.IGNORECASE,
)
_TEXT_SHADOW_RE = re.compile(r"text-shadow\s*:\s*[^;n]", re.IGNORECASE)
_GLOW_RE = re.compile(r"filter\s*:\s*drop-shadow", re.IGNORECASE)


def detect_taboos_static(html: str) -> TabooResult:
    detected: dict = {}
    suggestions: list = []

    n_gradients = len(_GRADIENT_RE.findall(html))
    detected["gradient_overuse"] = n_gradients >= 2
    if detected["gradient_overuse"]:
        suggestions.append(f"Remove or consolidate gradients (found {n_gradients}); use flat var(--color-bg) fills.")

    detected["glass_morphism"] = bool(_GLASS_RE.search(html))
    if detected["glass_morphism"]:
        suggestions.append("Remove backdrop-filter / glass morphism; use solid backgrounds.")

    indigo_hits = sum(1 for hex_ in _INDIGO_HEXES if hex_.lower() in html.lower())
    detected["multi_indigo_stacking"] = indigo_hits >= 2
    if detected["multi_indigo_stacking"]:
        suggestions.append(f"Stop stacking indigo shades (found {indigo_hits}); pick one accent color from design tokens.")

    detected["emoji_decorated_headers"] = bool(_EMOJI_HEADER_RE.search(html))
    if detected["emoji_decorated_headers"]:
        suggestions.append("Remove emoji from <h1>-<h6> headers; reserve emojis for body content if any.")

    detected["drop_shadow_text"] = bool(_TEXT_SHADOW_RE.search(html))
    if detected["drop_shadow_text"]:
        suggestions.append("Remove text-shadow on text; use weight or color for emphasis instead.")

    detected["glow_effects"] = bool(_GLOW_RE.search(html))
    if detected["glow_effects"]:
        suggestions.append("Remove drop-shadow filters from non-icon elements.")

    n_taboos = sum(1 for v in detected.values() if v)
    if n_taboos == 0:
        verdict = "clean"
    elif n_taboos == 1:
        verdict = "minor_issues"
    else:
        verdict = "redesign"

    return {
        "taboos_detected": detected,
        "n_taboos": n_taboos,
        "verdict": verdict,
        "suggestions": suggestions,
    }


def main(argv: list[str]) -> int:
    import json
    import sys
    if len(argv) < 2:
        print("Usage: design_taste_check.py <html_file>")
        return 2
    html = Path(argv[1]).read_text()
    result = detect_taboos_static(html)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["verdict"] == "clean" else 1


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
