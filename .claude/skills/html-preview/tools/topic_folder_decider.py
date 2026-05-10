"""topic_folder_decider — D17 logic: should_use_topic_folder?

Rule: report >= 200 lines OR >= 5 figures → topic folder.
Otherwise: single companion .preview.html file.

Figure count = markdown ![]() images + <!-- figure-needed: ... --> markers.
"""
from __future__ import annotations

from pathlib import Path

LINE_THRESHOLD = 200
FIGURE_THRESHOLD = 5


def count_figures(md_text: str) -> int:
    """Count both markdown !() images and <!-- figure-needed: ... --> markers."""
    import re
    img_count = len(re.findall(r"!\[[^\]]*\]\([^)]+\)", md_text))
    figure_needed_count = md_text.count("<!-- figure-needed:")
    return img_count + figure_needed_count


def should_use_topic_folder(md_path: Path) -> bool:
    """Return True if the report meets D17 conditions for topic folder."""
    text = Path(md_path).read_text()
    line_count = text.count("\n") + (0 if text.endswith("\n") else 1)
    fig_count = count_figures(text)
    return line_count >= LINE_THRESHOLD or fig_count >= FIGURE_THRESHOLD


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: topic_folder_decider.py <md_file>")
        return 2
    p = Path(argv[1])
    if not p.exists():
        print(f"Not found: {p}", file=sys.stderr)
        return 2
    text = p.read_text()
    lc = text.count("\n") + (0 if text.endswith("\n") else 1)
    fc = count_figures(text)
    decision = should_use_topic_folder(p)
    print(f"lines={lc} figures={fc} → topic_folder={decision}")
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
