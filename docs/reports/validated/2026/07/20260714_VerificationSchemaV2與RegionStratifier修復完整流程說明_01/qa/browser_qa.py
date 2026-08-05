#!/usr/bin/env python3
"""Browser QA for the standalone Verification schema repair report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


VIEWPORTS = {
    "desktop": {"width": 1440, "height": 1000},
    "mobile": {"width": 390, "height": 844},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    html_path = args.html.resolve(strict=True)
    output_dir = args.output_dir.resolve(strict=True)
    results = {}
    console_errors = []
    page_errors = []
    external_requests = []

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for name, viewport in VIEWPORTS.items():
            page = browser.new_page(viewport=viewport)
            page.on(
                "console",
                lambda message: console_errors.append(message.text)
                if message.type == "error"
                else None,
            )
            page.on("pageerror", lambda error: page_errors.append(str(error)))
            page.on(
                "request",
                lambda request: external_requests.append(request.url)
                if not request.url.startswith(("file:", "data:"))
                else None,
            )
            page.goto(html_path.as_uri(), wait_until="networkidle")

            title = page.title()
            h1_count = page.locator("h1").count()
            nav_links = page.locator("aside nav a")
            nav_count = nav_links.count()
            missing_targets = []
            for index in range(nav_count):
                href = nav_links.nth(index).get_attribute("href") or ""
                if not href.startswith("#") or page.locator(href).count() != 1:
                    missing_targets.append(href)

            dimensions = page.evaluate(
                """() => ({
                    innerWidth: window.innerWidth,
                    bodyWidth: document.body.scrollWidth,
                    htmlWidth: document.documentElement.scrollWidth,
                    svgWidth: document.querySelector('.flow svg').getBoundingClientRect().width,
                    svgHeight: document.querySelector('.flow svg').getBoundingClientRect().height,
                    tables: document.querySelectorAll('table').length,
                    wrappedTables: document.querySelectorAll('.table-wrap > table').length,
                    overflowElements: [...document.querySelectorAll('body *')]
                      .map(element => {
                        const rect = element.getBoundingClientRect();
                        return {
                          tag: element.tagName,
                          className: String(element.className || ''),
                          id: element.id || '',
                          left: Math.round(rect.left),
                          right: Math.round(rect.right),
                          width: Math.round(rect.width)
                        };
                      })
                      .filter(row => row.left < -1 || row.right > window.innerWidth + 1)
                      .slice(0, 30)
                })"""
            )
            body_overflow = max(dimensions["bodyWidth"], dimensions["htmlWidth"]) > (
                dimensions["innerWidth"] + 1
            )
            placeholder_count = page.locator(
                "text=/P2_VALIDATOR|VALIDATED_PENDING|TODO|TBD/"
            ).count()
            status_text = page.locator("#report-status").inner_text()

            screenshot = output_dir / (name + ".png")
            page.screenshot(path=str(screenshot), full_page=True)
            results[name] = {
                "viewport": viewport,
                "title": title,
                "h1_count": h1_count,
                "nav_count": nav_count,
                "missing_nav_targets": missing_targets,
                "body_overflow": body_overflow,
                "svg_width": round(dimensions["svgWidth"], 2),
                "svg_height": round(dimensions["svgHeight"], 2),
                "tables": dimensions["tables"],
                "wrapped_tables": dimensions["wrappedTables"],
                "overflow_elements": dimensions["overflowElements"],
                "placeholder_count": placeholder_count,
                "report_status": status_text,
                "screenshot": screenshot.name,
            }
            page.close()
        browser.close()

    checks = {
        "console_errors_empty": not console_errors,
        "page_errors_empty": not page_errors,
        "external_requests_empty": not external_requests,
        "desktop_no_body_overflow": not results["desktop"]["body_overflow"],
        "mobile_no_body_overflow": not results["mobile"]["body_overflow"],
        "all_nav_targets_exist": all(not row["missing_nav_targets"] for row in results.values()),
        "single_h1": all(row["h1_count"] == 1 for row in results.values()),
        "svg_visible": all(row["svg_width"] > 0 and row["svg_height"] > 0 for row in results.values()),
        "all_tables_wrapped": all(row["tables"] == row["wrapped_tables"] for row in results.values()),
        "no_placeholders": all(row["placeholder_count"] == 0 for row in results.values()),
        "final_status": all(row["report_status"] == "VALIDATED" for row in results.values()),
    }
    payload = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "html": str(html_path),
        "viewports": results,
        "checks": checks,
        "console_errors": console_errors,
        "page_errors": page_errors,
        "external_requests": external_requests,
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
