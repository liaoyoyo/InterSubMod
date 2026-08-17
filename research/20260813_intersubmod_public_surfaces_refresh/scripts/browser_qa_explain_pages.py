#!/usr/bin/env python3
"""Chromium desktop/mobile/no-JS/print smoke for the 17 explain pages."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

from playwright.sync_api import Browser, BrowserContext, Page, sync_playwright


BASE = "http://127.0.0.1:8765/docs/explain"
LOCAL_ORIGIN = "http://127.0.0.1:8765/"
PAGES = [
    "index.standalone.html",
    "01_background-glossary.standalone.html",
    "02_ism-core.standalone.html",
    "03_methylation-read-filter.standalone.html",
    "04_subclone-reconstruction-chr2-18M.standalone.html",
    "05_subclone-correction-audit-chr2-18M.standalone.html",
    "06_ism-subclone-pipeline-concept.standalone.html",
    "07_subclone-judgment-workstation-chr2-18M.standalone.html",
    "08_subclone-logic-chain-chr2-18M.standalone.html",
    "09_three-stats-division-of-labor.standalone.html",
    "10_ism-cpp-vs-chr2-subclone-capability.standalone.html",
    "11_system-map-overview.standalone.html",
    "12_intersubmod-io.standalone.html",
    "13_longlineage-io.standalone.html",
    "14_upstream-data.standalone.html",
    "15_python-html-layer.standalone.html",
    "16_how-to-run.standalone.html",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default=BASE)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def inspect_page(page: Page, name: str, profile: str) -> dict[str, Any]:
    console_errors: list[str] = []
    page_errors: list[str] = []
    failed_requests: list[str] = []
    page.on(
        "console",
        lambda message: console_errors.append(message.text) if message.type == "error" else None,
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.on(
        "requestfailed",
        lambda request: failed_requests.append(f"{request.url}: {request.failure}"),
    )
    response = page.goto(f"{BASE}/{name}", wait_until="networkidle", timeout=30_000)
    status = response.status if response else None
    title = page.title().strip()
    h1_count = page.locator("h1").count()
    body_text_length = len(page.locator("body").inner_text().strip())
    overflow = page.evaluate(
        "Math.max(document.documentElement.scrollWidth, document.body.scrollWidth) "
        "> document.documentElement.clientWidth + 1"
    )
    overflow_elements = page.evaluate(
        """() => {
          const viewport = document.documentElement.clientWidth;
          return Array.from(document.querySelectorAll('body *'))
            .map((element) => {
              const rect = element.getBoundingClientRect();
              const style = getComputedStyle(element);
              return {
                tag: element.tagName.toLowerCase(),
                id: element.id || '',
                class_name: String(element.className || '').slice(0, 100),
                left: Math.round(rect.left),
                right: Math.round(rect.right),
                width: Math.round(rect.width),
                scroll_width: element.scrollWidth,
                client_width: element.clientWidth,
                overflow_x: style.overflowX,
                text: (element.innerText || '').replace(/\\s+/g, ' ').slice(0, 120),
              };
            })
            .filter((item) => item.right > viewport + 1 || item.left < -1)
            .slice(0, 20);
        }"""
    ) if overflow else []
    visible_svgs = page.locator("svg:visible").count()
    return {
        "page": name,
        "profile": profile,
        "status": status,
        "title": title,
        "h1_count": h1_count,
        "body_text_length": body_text_length,
        "visible_svg": visible_svgs,
        "horizontal_overflow": bool(overflow),
        "overflow_elements": overflow_elements,
        "console_errors": console_errors,
        "page_errors": page_errors,
        "failed_requests": failed_requests,
    }


def context_rows(browser: Browser, profile: str, **context_options: Any) -> list[dict[str, Any]]:
    context: BrowserContext = browser.new_context(**context_options)
    rows: list[dict[str, Any]] = []
    try:
        for name in PAGES:
            page = context.new_page()
            try:
                if profile == "print":
                    page.emulate_media(media="print")
                rows.append(inspect_page(page, name, profile))
            finally:
                page.close()
    finally:
        context.close()
    return rows


def row_errors(row: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    if row["status"] != 200:
        errors.append(f"HTTP status {row['status']}")
    if not row["title"]:
        errors.append("empty title")
    if row["h1_count"] < 1:
        errors.append("missing h1")
    if row["body_text_length"] < 200:
        errors.append(f"body text too short: {row['body_text_length']}")
    if row["horizontal_overflow"]:
        errors.append("horizontal overflow")
    errors.extend(f"console: {item}" for item in row["console_errors"])
    errors.extend(f"pageerror: {item}" for item in row["page_errors"])
    # Do not fail on intentionally external links/images in this smoke; failed local requests do fail.
    errors.extend(
        f"requestfailed: {item}"
        for item in row["failed_requests"]
        if item.startswith(LOCAL_ORIGIN)
    )
    return errors


def main() -> int:
    args = parse_args()
    global BASE, LOCAL_ORIGIN
    BASE = args.base.rstrip("/")
    parsed = urlsplit(BASE)
    LOCAL_ORIGIN = f"{parsed.scheme}://{parsed.netloc}/"
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            rows = []
            rows.extend(
                context_rows(
                    browser,
                    "desktop",
                    viewport={"width": 1440, "height": 1000},
                    device_scale_factor=1,
                )
            )
            rows.extend(
                context_rows(
                    browser,
                    "mobile",
                    viewport={"width": 390, "height": 844},
                    device_scale_factor=1,
                    is_mobile=True,
                )
            )
            rows.extend(
                context_rows(
                    browser,
                    "no_js",
                    viewport={"width": 1280, "height": 900},
                    java_script_enabled=False,
                )
            )
            rows.extend(
                context_rows(
                    browser,
                    "print",
                    viewport={"width": 1280, "height": 900},
                    device_scale_factor=1,
                )
            )
        finally:
            browser.close()

    failures = []
    for row in rows:
        errors = row_errors(row)
        if errors:
            failures.append({"page": row["page"], "profile": row["profile"], "errors": errors})
    receipt = {
        "pages": len(PAGES),
        "profiles": ["desktop", "mobile", "no_js", "print"],
        "page_profile_checks": len(rows),
        "failures": failures,
        "verdict": "PASS" if not failures else "FAIL",
        "rows": rows,
    }
    rendered = json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
