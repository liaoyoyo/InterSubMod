"""md_to_html -- convert markdown source to HTML body fragment.

Uses the Python `markdown` library with extensions for tables, fenced code,
and pass-through HTML (so `<details>` survives).
Strips YAML frontmatter (returned separately for template variables).
"""
from __future__ import annotations

import re
from typing import Tuple

import markdown
import yaml

# Markdown extensions: tables, fenced code, no-implicit-paragraphs around HTML
_EXTENSIONS = ["tables", "fenced_code", "md_in_html"]
_FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?)\n---\s*\n", re.DOTALL)


def parse_frontmatter(md_text: str) -> Tuple[dict, str]:
    """Return (frontmatter_dict, body_without_frontmatter)."""
    m = _FRONTMATTER_RE.match(md_text)
    if not m:
        return {}, md_text
    fm_yaml = m.group(1)
    body = md_text[m.end():]
    try:
        fm = yaml.safe_load(fm_yaml) or {}
    except yaml.YAMLError:
        fm = {}
    if not isinstance(fm, dict):
        fm = {}
    # Normalize date to ISO string if datetime
    if "date" in fm and not isinstance(fm["date"], str):
        fm["date"] = str(fm["date"])
    return fm, body


def convert(md_text: str) -> str:
    """Convert markdown source to HTML body. Strips YAML frontmatter."""
    _, body = parse_frontmatter(md_text)
    return markdown.markdown(body, extensions=_EXTENSIONS, output_format="html")


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 2:
        print("Usage: md_to_html.py <input.md> [<output.html>]")
        return 2
    from pathlib import Path
    src = Path(argv[1])
    md_text = src.read_text()
    fm, _ = parse_frontmatter(md_text)
    html = convert(md_text)
    if len(argv) >= 3:
        Path(argv[2]).write_text(html)
        print(f"[md_to_html] {len(html)} chars -> {argv[2]} (title={fm.get('title')!r})")
    else:
        print(html)
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
