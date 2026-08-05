#!/usr/bin/env python3
"""Report visible DOM elements that extend outside the desktop viewport."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

HEADLESS_SHELL = (
    "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/"
    "chromium_headless_shell-1223/chrome-headless-shell-linux64/chrome-headless-shell"
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", type=Path)
    parser.add_argument("--width", type=int, default=1440)
    args = parser.parse_args()
    uri = args.html.resolve().as_uri()
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True, executable_path=HEADLESS_SHELL)
        page = browser.new_page(viewport={"width": args.width, "height": 1100}, device_scale_factor=1)
        page.goto(uri, wait_until="load", timeout=60000)
        time.sleep(8)
        probe = """() => {
              const clientWidth = document.documentElement.clientWidth;
              const visible = (el) => {
                const style = getComputedStyle(el);
                const rect = el.getBoundingClientRect();
                return style.display !== 'none' && style.visibility !== 'hidden' &&
                  Number(style.opacity || 1) !== 0 && rect.width > 0 && rect.height > 0;
              };
              const rows = Array.from(document.querySelectorAll('body *'))
                .filter(visible)
                .map((el) => {
                  const r = el.getBoundingClientRect();
                  return {
                    tag: el.tagName,
                    id: el.id || '',
                    className: typeof el.className === 'string' ? el.className.slice(0, 180) : '',
                    left: Math.round(r.left * 10) / 10,
                    right: Math.round(r.right * 10) / 10,
                    width: Math.round(r.width * 10) / 10,
                    scrollWidth: el.scrollWidth,
                    clientWidth: el.clientWidth,
                    overflowX: getComputedStyle(el).overflowX,
                    text: (el.innerText || el.textContent || '').replace(/\\s+/g, ' ').slice(0, 180),
                  };
                })
                .filter((row) => row.right > clientWidth + 1 || row.left < -1)
                .sort((a, b) => (b.right - clientWidth) - (a.right - clientWidth));
              return {
                clientWidth,
                documentScrollWidth: document.documentElement.scrollWidth,
                bodyScrollWidth: document.body.scrollWidth,
                offenderCount: rows.length,
                offenders: rows.slice(0, 50),
              };
            }"""
        results = []
        for index, frame in enumerate(page.frames):
            try:
                result = frame.evaluate(probe)
                result["frame_index"] = index
                result["url"] = frame.url
                results.append(result)
            except Exception as error:  # pragma: no cover - diagnostic fallback
                results.append({"frame_index": index, "url": frame.url, "error": str(error)})
        browser.close()
    print(json.dumps(results, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
