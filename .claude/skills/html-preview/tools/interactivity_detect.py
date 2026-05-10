"""interactivity_detect — choose L1/L2/L3 rendering level.

L1: raw HTML (only when explicitly forced via <!-- html-preview: force-l1 -->)
L2: MVP.css default (D15)
L3: + Alpine.js when <!-- interactive: tab|filter|slider --> is present
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

Level = Literal["l1", "l2", "l3"]

_INTERACTIVE_MARKER = re.compile(r"<!--\s*interactive\s*:", re.IGNORECASE)
_FORCE_L1_MARKER = re.compile(r"<!--\s*html-preview\s*:\s*force-l1\s*-->", re.IGNORECASE)


def detect_level(md_path: Path) -> Level:
    text = Path(md_path).read_text()
    if _FORCE_L1_MARKER.search(text):
        return "l1"
    if _INTERACTIVE_MARKER.search(text):
        return "l3"
    return "l2"


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: interactivity_detect.py <md_file>")
        return 2
    print(detect_level(Path(argv[1])))
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
