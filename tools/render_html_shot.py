#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""render_html_shot.py — headless HTML -> full-page PNG screenshot (playwright chromium).

Implements the verify-workstation W7 "playwright 截圖檢查無溢出/錯位" step as a reusable
tool, and closes the artifact-QA loop (AI-Scientist-v2 VLM-feedback pattern, human-gated):
  render -> screenshot -> Claude Reads the PNG (the "VLM") -> detect tofu/overflow/broken-img
  / unverified numbers -> fix -> re-render until pass.

Vision review itself is done by Claude reading the PNG with the Read tool (same pattern as
tools/ppt_toolkit/claude_vision_review.py). This tool only produces the PNG(s).

Usage:
  python3 tools/render_html_shot.py <input.html> -o <out.png> [--width 1400] [--full]
  # relative <img src="figs/..."> resolve against the HTML file location (file:// base).
"""
import argparse
import os
import sys
from pathlib import Path


def render(html_path, out_png, width=1400, full=True, wait_ms=1200):
    from playwright.sync_api import sync_playwright
    url = Path(html_path).resolve().as_uri()  # file:// so relative figs/ resolve
    with sync_playwright() as p:
        b = p.chromium.launch(args=["--no-sandbox"])
        pg = b.new_page(viewport={"width": width, "height": 1000},
                        device_scale_factor=1)
        errors = []
        pg.on("pageerror", lambda e: errors.append(str(e)))
        pg.goto(url, wait_until="networkidle", timeout=30000)
        pg.wait_for_timeout(wait_ms)
        # broken-image detection (naturalWidth==0 on loaded <img>)
        broken = pg.evaluate(
            "Array.from(document.images).filter(i=>i.complete && i.naturalWidth===0)"
            ".map(i=>i.getAttribute('src'))")
        pg.screenshot(path=out_png, full_page=full)
        b.close()
    return {"errors": errors, "broken_images": broken}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("html")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--width", type=int, default=1400)
    ap.add_argument("--full", action="store_true", default=True)
    ap.add_argument("--viewport-only", dest="full", action="store_false")
    a = ap.parse_args()
    if not os.path.exists(a.html):
        print(f"ERROR: no such html {a.html}", file=sys.stderr)
        sys.exit(1)
    r = render(a.html, a.out, a.width, a.full)
    print(f"[ok] wrote {a.out}")
    if r["errors"]:
        print(f"⚠ pageerror x{len(r['errors'])}: {r['errors'][:3]}")
    if r["broken_images"]:
        print(f"🔴 broken images x{len(r['broken_images'])}: {r['broken_images'][:5]}")
    if not r["errors"] and not r["broken_images"]:
        print("✅ no pageerror / no broken images (still Read the PNG for tofu/overflow)")


if __name__ == "__main__":
    main()
