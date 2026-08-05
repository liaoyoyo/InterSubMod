#!/usr/bin/env python3
"""Audit a rendered portable report for horizontal overflow."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from playwright.sync_api import sync_playwright


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--width", type=int, default=1440)
    parser.add_argument("--height", type=int, default=1000)
    parser.add_argument("--delay-ms", type=int, default=32)
    args = parser.parse_args()

    html = args.html.resolve()
    if not html.exists():
        raise FileNotFoundError(html)

    console_errors: list[str] = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": args.width, "height": args.height})
        page.emulate_media(reduced_motion="reduce", color_scheme="light")
        page.on(
            "console",
            lambda message: console_errors.append(message.text)
            if message.type == "error"
            else None,
        )
        page.goto(html.as_uri(), wait_until="networkidle")
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            if page.evaluate(
                "document.documentElement.dataset.dataAnalyticsPortableReader"
            ) == "ready":
                break
            page.wait_for_timeout(50)
        else:
            raise RuntimeError("portable reader did not reach ready state")
        page.wait_for_timeout(args.delay_ms)
        result = page.evaluate(
            """
            () => {
              const vw = document.documentElement.clientWidth;
              const visible = (el) => {
                const style = getComputedStyle(el);
                const rect = el.getBoundingClientRect();
                return style.display !== 'none' && style.visibility !== 'hidden' && rect.width > 0 && rect.height > 0;
              };
              const label = (el) => {
                const id = el.id ? `#${el.id}` : '';
                const classes = [...el.classList].slice(0, 4).map(c => `.${c}`).join('');
                const text = (el.textContent || '').replace(/\\s+/g, ' ').trim().slice(0, 120);
                return `${el.tagName.toLowerCase()}${id}${classes} :: ${text}`;
              };
              const offenders = [...document.querySelectorAll('body *')]
                .filter(visible)
                .map(el => {
                  const rect = el.getBoundingClientRect();
                  return {
                    label: label(el),
                    left: Math.round(rect.left * 10) / 10,
                    right: Math.round(rect.right * 10) / 10,
                    width: Math.round(rect.width * 10) / 10,
                    clientWidth: el.clientWidth,
                    scrollWidth: el.scrollWidth,
                    overflowX: getComputedStyle(el).overflowX,
                  };
                })
                .filter(row => row.right > vw + 1 || row.left < -1 || row.scrollWidth > row.clientWidth + 1)
                .sort((a, b) => Math.max(b.right - vw, b.scrollWidth - b.clientWidth) - Math.max(a.right - vw, a.scrollWidth - a.clientWidth))
                .slice(0, 40);
              const edgeOffenders = [...document.querySelectorAll('body *')]
                .filter(visible)
                .map(el => {
                  const rect = el.getBoundingClientRect();
                  return {
                    label: label(el),
                    left: Math.round(rect.left * 10) / 10,
                    right: Math.round(rect.right * 10) / 10,
                    width: Math.round(rect.width * 10) / 10,
                    overflowX: getComputedStyle(el).overflowX,
                  };
                })
                .filter(row => row.right > vw + 1 && row.right < vw + 50)
                .sort((a, b) => b.right - a.right)
                .slice(0, 80);
              return {
                viewportWidth: vw,
                documentClientWidth: document.documentElement.clientWidth,
                documentScrollWidth: document.documentElement.scrollWidth,
                bodyClientWidth: document.body.clientWidth,
                bodyScrollWidth: document.body.scrollWidth,
                state: document.documentElement.dataset.reportState || document.body.dataset.reportState || null,
                offenders,
                edgeOffenders,
              };
            }
            """
        )
        browser.close()

    print(json.dumps({"html": str(html), "layout": result, "console_errors": console_errors}, indent=2))


if __name__ == "__main__":
    main()
