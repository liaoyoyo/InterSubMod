#!/usr/bin/env python3
"""Diagnose page-level overflow in a portable Data Analytics report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--width", required=True, type=int)
    parser.add_argument("--height", default=900, type=int)
    parser.add_argument("--executable-path", type=Path)
    parser.add_argument("--show-scrollbars", action="store_true")
    parser.add_argument("--screenshot", type=Path)
    parser.add_argument("--wait-ms", default=5_000, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    html_path = args.html.resolve()
    if not html_path.is_file():
        raise SystemExit(f"Missing HTML: {html_path}")

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True,
            executable_path=str(args.executable_path.resolve())
            if args.executable_path
            else None,
            ignore_default_args=["--hide-scrollbars"] if args.show_scrollbars else None,
        )
        page = browser.new_page(viewport={"width": args.width, "height": args.height})
        console_errors: list[str] = []
        page.on(
            "console",
            lambda message: console_errors.append(message.text)
            if message.type == "error"
            else None,
        )
        page.goto(html_path.as_uri(), wait_until="networkidle")
        page.wait_for_timeout(args.wait_ms)
        if args.screenshot:
            args.screenshot.parent.mkdir(parents=True, exist_ok=True)
            page.screenshot(path=str(args.screenshot.resolve()), full_page=True)
            args.screenshot.resolve().chmod(0o444)

        diagnostics = page.evaluate(
            """
            () => {
              const viewportWidth = document.documentElement.clientWidth;
              const describe = (element) => {
                const rect = element.getBoundingClientRect();
                const style = getComputedStyle(element);
                const ancestors = [];
                let current = element;
                for (let depth = 0; current && depth < 5; depth += 1) {
                  const id = current.id ? `#${current.id}` : '';
                  const classes = current.classList?.length
                    ? `.${Array.from(current.classList).slice(0, 3).join('.')}`
                    : '';
                  ancestors.push(`${current.tagName.toLowerCase()}${id}${classes}`);
                  current = current.parentElement;
                }
                return {
                  element: ancestors[0],
                  ancestors,
                  left: Math.round(rect.left * 10) / 10,
                  right: Math.round(rect.right * 10) / 10,
                  width: Math.round(rect.width * 10) / 10,
                  clientWidth: element.clientWidth,
                  scrollWidth: element.scrollWidth,
                  overflowX: style.overflowX,
                  whiteSpace: style.whiteSpace,
                  text: (element.innerText || element.textContent || '').trim().slice(0, 180),
                };
              };
              const all = Array.from(document.querySelectorAll('body *'));
              const outside = all
                .filter((element) => {
                  const rect = element.getBoundingClientRect();
                  return rect.width > 0 && (rect.right > viewportWidth + 1 || rect.left < -1);
                })
                .map(describe)
                .sort((a, b) => b.right - a.right)
                .slice(0, 50);
              const internalOverflow = all
                .filter((element) => element.scrollWidth > element.clientWidth + 1)
                .map(describe)
                .sort((a, b) => (b.scrollWidth - b.clientWidth) - (a.scrollWidth - a.clientWidth))
                .slice(0, 50);
              const layoutRows = Array.from(
                document.querySelectorAll('.analytics-layout-row')
              ).slice(0, 50).map((row) => {
                const rowRect = row.getBoundingClientRect();
                const descendants = Array.from(row.querySelectorAll('*'))
                  .map((element) => {
                    const rect = element.getBoundingClientRect();
                    const style = getComputedStyle(element);
                    return {
                      ...describe(element),
                      boxSizing: style.boxSizing,
                      marginLeft: style.marginLeft,
                      marginRight: style.marginRight,
                      paddingLeft: style.paddingLeft,
                      paddingRight: style.paddingRight,
                    };
                  })
                  .filter((entry) => entry.right > rowRect.right + 1)
                  .sort((a, b) => b.right - a.right)
                  .slice(0, 10);
                return {
                  blockId: row.querySelector('[data-layout-block-id]')
                    ?.getAttribute('data-layout-block-id') || null,
                  row: describe(row),
                  descendants,
                };
              }).filter((entry) => entry.descendants.length > 0);
              const specialElements = Array.from(document.querySelectorAll(
                'header, #data-analytics-portable-reader, #data-analytics-portable-fallback, ' +
                '#data-analytics-portable-reader-root, [data-portable-fallback]'
              )).map((element) => {
                const style = getComputedStyle(element);
                return {
                  ...describe(element),
                  ariaHidden: element.getAttribute('aria-hidden'),
                  hidden: element.hidden,
                  display: style.display,
                  position: style.position,
                  visibility: style.visibility,
                  boxSizing: style.boxSizing,
                  marginLeft: style.marginLeft,
                  marginRight: style.marginRight,
                  paddingLeft: style.paddingLeft,
                  paddingRight: style.paddingRight,
                };
              });
              const edgeElements = all
                .filter((element) => {
                  const rect = element.getBoundingClientRect();
                  return rect.width > 0 &&
                    rect.right > document.documentElement.clientWidth + 1 &&
                    rect.right <= document.documentElement.scrollWidth + 1;
                })
                .map(describe)
                .sort((a, b) => b.right - a.right)
                .slice(0, 50);
              const blockRows = Array.from(
                document.querySelectorAll('.analytics-layout-row')
              ).filter((element) => {
                const rect = element.getBoundingClientRect();
                const style = getComputedStyle(element);
                return rect.width > 0 && rect.height > 0 && style.visibility !== 'hidden';
              });
              const overlapPairs = [];
              for (let leftIndex = 0; leftIndex < blockRows.length; leftIndex += 1) {
                const left = blockRows[leftIndex];
                const leftRect = left.getBoundingClientRect();
                for (let rightIndex = leftIndex + 1; rightIndex < blockRows.length; rightIndex += 1) {
                  const right = blockRows[rightIndex];
                  const rightRect = right.getBoundingClientRect();
                  const overlapWidth = Math.min(leftRect.right, rightRect.right) -
                    Math.max(leftRect.left, rightRect.left);
                  const overlapHeight = Math.min(leftRect.bottom, rightRect.bottom) -
                    Math.max(leftRect.top, rightRect.top);
                  if (overlapWidth > 1 && overlapHeight > 1) {
                    overlapPairs.push({
                      leftBlockId: left.querySelector('[data-layout-block-id]')
                        ?.getAttribute('data-layout-block-id') || null,
                      rightBlockId: right.querySelector('[data-layout-block-id]')
                        ?.getAttribute('data-layout-block-id') || null,
                      overlapWidth: Math.round(overlapWidth * 10) / 10,
                      overlapHeight: Math.round(overlapHeight * 10) / 10,
                    });
                  }
                }
              }
              return {
                viewportWidth,
                documentClientWidth: document.documentElement.clientWidth,
                documentScrollWidth: document.documentElement.scrollWidth,
                bodyClientWidth: document.body.clientWidth,
                bodyScrollWidth: document.body.scrollWidth,
                outside,
                internalOverflow,
                layoutRows,
                specialElements,
                edgeElements,
                overlapCount: overlapPairs.length,
                overlapPairs: overlapPairs.slice(0, 50),
              };
            }
            """
        )
        diagnostics["consoleErrors"] = console_errors
        diagnostics["pass"] = (
            diagnostics["documentScrollWidth"] <= diagnostics["documentClientWidth"] + 1
            and diagnostics["bodyScrollWidth"] <= diagnostics["bodyClientWidth"] + 1
            and diagnostics["overlapCount"] == 0
            and not console_errors
        )
        print(json.dumps(diagnostics, ensure_ascii=False, indent=2))
        browser.close()


if __name__ == "__main__":
    main()
