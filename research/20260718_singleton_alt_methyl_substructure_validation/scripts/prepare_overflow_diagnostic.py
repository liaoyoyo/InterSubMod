#!/usr/bin/env python3
"""Prepare a self-reporting portable HTML file for classic-scrollbar diagnosis."""

from __future__ import annotations

import argparse
from pathlib import Path


UPSTREAM = (
    "width:100vw;height:48px;min-height:48px;"
    "margin-right:calc(50% - 50vw);margin-left:calc(50% - 50vw);"
)
PATCHED = "width:100%;height:48px;min-height:48px;margin-right:0;margin-left:0;"
RUNTIME_PATCH = (
    '<style data-hcc1395-scrollbar-compat="true">'
    ".analytics-top-bar{width:100%!important;margin-right:0!important;margin-left:0!important}"
    "</style>"
)

PROBE = r"""
<script>
setTimeout(() => {
  const clientWidth = document.documentElement.clientWidth;
  const visible = (el) => {
    const style = getComputedStyle(el);
    const rect = el.getBoundingClientRect();
    return style.display !== "none" && style.visibility !== "hidden" &&
      Number(style.opacity || 1) !== 0 && rect.width > 0 && rect.height > 0;
  };
  const rows = Array.from(document.querySelectorAll("body *"))
    .filter(visible)
    .map((el) => {
      const rect = el.getBoundingClientRect();
      const style = getComputedStyle(el);
      return {
        tag: el.tagName,
        id: el.id || "",
        className: typeof el.className === "string" ? el.className.slice(0, 180) : "",
        left: Math.round(rect.left * 10) / 10,
        right: Math.round(rect.right * 10) / 10,
        width: Math.round(rect.width * 10) / 10,
        clientWidth: el.clientWidth,
        scrollWidth: el.scrollWidth,
        overflowX: style.overflowX,
        position: style.position,
        text: (el.innerText || el.textContent || "").replace(/\s+/g, " ").slice(0, 100),
      };
    })
    .filter((row) => row.right > clientWidth + 1 || row.left < -1)
    .sort((a, b) => (b.right - clientWidth) - (a.right - clientWidth));
  const result = {
    innerWidth,
    clientWidth,
    documentScrollWidth: document.documentElement.scrollWidth,
    bodyClientWidth: document.body.clientWidth,
    bodyScrollWidth: document.body.scrollWidth,
    offenderCount: rows.length,
    offenders: rows.slice(0, 100),
  };
  const output = document.createElement("pre");
  output.id = "overflow-diagnostics-json";
  output.textContent = JSON.stringify(result);
  document.body.prepend(output);
}, 10000);
</script>
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    html = args.input.read_text(encoding="utf-8")
    count = html.count(UPSTREAM)
    if count not in (0, 1):
        raise ValueError(f"Expected zero or one upstream header rules, found {count}")
    html = html.replace(UPSTREAM, PATCHED)
    if "</head>" not in html:
        raise ValueError("Input HTML has no closing head tag")
    html = html.replace("</head>", f"{RUNTIME_PATCH}</head>", 1)
    if "</body>" not in html:
        raise ValueError("Input HTML has no closing body tag")
    html = html.replace("</body>", f"{PROBE}</body>", 1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html, encoding="utf-8")
    print(f"wrote={args.output.resolve()} patched_header_rules={count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
