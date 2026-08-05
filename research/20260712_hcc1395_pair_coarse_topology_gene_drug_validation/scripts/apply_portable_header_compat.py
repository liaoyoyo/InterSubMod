#!/usr/bin/env python3
"""Prevent 100vw report headers from overflowing beside classic scrollbars."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


STYLE_ID = "hcc1395-portable-header-scrollbar-compat"
STYLE = f"""<style id=\"{STYLE_ID}\">
/* Keep the portable reader's full-width sticky bars inside the layout viewport
   when the browser reserves physical width for a vertical scrollbar. */
.portable-page-header,
.analytics-top-bar {{
  width: 100% !important;
  max-width: 100% !important;
  margin-right: 0 !important;
  margin-left: 0 !important;
  inset-inline: 0 auto !important;
}}
</style>"""


def digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", type=Path)
    args = parser.parse_args()

    html = args.html.read_text(encoding="utf-8")
    before = digest(html)
    html = html.replace('<html lang="en"', '<html lang="zh-TW"', 1)
    if f'id="{STYLE_ID}"' not in html:
        marker = "</head>"
        if marker not in html:
            raise RuntimeError("portable HTML has no </head> marker")
        html = html.replace(marker, f"{STYLE}\n{marker}", 1)
    after = digest(html)
    if after != before:
        args.html.write_text(html, encoding="utf-8")
    print(f"input_output={args.html}")
    print(f"sha256_before={before}")
    print(f"sha256_after={after}")
    print(f"style_id={STYLE_ID}")
    print("document_language=zh-TW")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
