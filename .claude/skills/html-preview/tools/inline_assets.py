"""inline_assets — replace <img src="local.png"> with data: URIs.

Skips external URLs (http/https) and existing data: URIs. Inlines relative
paths resolved from base_dir. Missing files: keeps src + adds HTML comment.
"""
from __future__ import annotations

import base64
import re
from pathlib import Path

_IMG_RE = re.compile(r'<img\s+[^>]*?src="([^"]+)"[^>]*?>', re.IGNORECASE)


def _inline_one(match: re.Match, base_dir: Path) -> str:
    full_tag = match.group(0)
    src = match.group(1)
    if src.startswith(("http://", "https://", "data:")):
        return full_tag
    src_path = (base_dir / src).resolve()
    if not src_path.exists():
        return f"<!-- inline_assets: file not found at {src} relative to {base_dir} -->" + full_tag
    encoded = base64.b64encode(src_path.read_bytes()).decode("ascii")
    suffix = src_path.suffix.lower().lstrip(".")
    mime = {"png": "png", "jpg": "jpeg", "jpeg": "jpeg", "gif": "gif", "svg": "svg+xml"}.get(suffix, suffix)
    data_uri = f"data:image/{mime};base64,{encoded}"
    return full_tag.replace(f'src="{src}"', f'src="{data_uri}"', 1)


def inline_images(html: str, *, base_dir: Path) -> str:
    return _IMG_RE.sub(lambda m: _inline_one(m, base_dir), html)


def main(argv: list[str]) -> int:
    import sys
    if len(argv) < 3:
        print("Usage: inline_assets.py <input.html> <base_dir>")
        return 2
    html = Path(argv[1]).read_text()
    base = Path(argv[2])
    out = inline_images(html, base_dir=base)
    print(out)
    return 0


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv))
