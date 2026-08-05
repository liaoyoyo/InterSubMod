#!/usr/bin/env python3
"""Locate rendered elements that cause horizontal overflow in a portable report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--browser", required=True, type=Path)
    parser.add_argument("--screenshot", type=Path)
    parser.add_argument("--width", type=int, default=1440)
    parser.add_argument("--height", type=int, default=1000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    html = args.html.resolve(strict=True)
    browser_path = args.browser.resolve(strict=True)

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=str(browser_path),
            headless=True,
        )
        page = browser.new_page(viewport={"width": args.width, "height": args.height})
        page.goto(html.as_uri(), wait_until="networkidle")
        page.wait_for_timeout(5_000)
        diagnostics = page.evaluate(
            """() => {
                const viewportWidth = document.documentElement.clientWidth;
                const rows = [];
                for (const element of document.querySelectorAll("body *")) {
                    const rect = element.getBoundingClientRect();
                    const style = getComputedStyle(element);
                    const excessRight = rect.right - viewportWidth;
                    const excessLeft = -rect.left;
                    if (excessRight <= 1 && excessLeft <= 1) continue;
                    rows.push({
                        tag: element.tagName.toLowerCase(),
                        id: element.id || null,
                        className: typeof element.className === "string"
                            ? element.className
                            : null,
                        blockId: element.closest("[data-portable-block-id]")
                            ?.getAttribute("data-portable-block-id") || null,
                        rect: {
                            left: Math.round(rect.left * 100) / 100,
                            right: Math.round(rect.right * 100) / 100,
                            width: Math.round(rect.width * 100) / 100,
                        },
                        scrollWidth: element.scrollWidth,
                        clientWidth: element.clientWidth,
                        overflowX: style.overflowX,
                        whiteSpace: style.whiteSpace,
                        excessRight: Math.round(excessRight * 100) / 100,
                        excessLeft: Math.round(excessLeft * 100) / 100,
                        text: (element.textContent || "").trim().slice(0, 180),
                    });
                }
                rows.sort((a, b) =>
                    Math.max(b.excessRight, b.excessLeft) -
                    Math.max(a.excessRight, a.excessLeft)
                );
                return {
                    clientWidth: viewportWidth,
                    documentScrollWidth: document.documentElement.scrollWidth,
                    bodyScrollWidth: document.body.scrollWidth,
                    offenders: rows.slice(0, 100),
                };
            }"""
        )
        if args.screenshot is not None:
            args.screenshot.parent.mkdir(parents=True, exist_ok=True)
            page.screenshot(path=str(args.screenshot), full_page=False)
        browser.close()

    print(json.dumps(diagnostics, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
