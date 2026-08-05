#!/usr/bin/env python3
"""Browser-level QA for the standalone optimized backend HTML report."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import time
from typing import Any

from playwright.sync_api import sync_playwright


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--html", type=pathlib.Path, required=True)
    parser.add_argument("--screenshot", type=pathlib.Path, required=True)
    parser.add_argument("--receipt", type=pathlib.Path, required=True)
    args = parser.parse_args()

    html_path = args.html.resolve(strict=True)
    args.screenshot.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, Any]] = []
    console_errors: list[str] = []
    page_errors: list[str] = []

    def check(name: str, condition: bool, observed: Any) -> None:
        checks.append(
            {
                "name": name,
                "status": "PASS" if condition else "FAIL",
                "observed": observed,
            }
        )

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        browser_version = browser.version
        page = browser.new_page(
            viewport={"width": 1440, "height": 1200},
            device_scale_factor=1,
            color_scheme="light",
        )
        page.on(
            "console",
            lambda message: console_errors.append(message.text)
            if message.type == "error"
            else None,
        )
        page.on("pageerror", lambda error: page_errors.append(str(error)))
        page.goto(html_path.as_uri(), wait_until="networkidle")

        check(
            "document_title",
            page.title()
            == "Hypercube Optimized Backend｜完整流程、正確性與效能稽核",
            page.title(),
        )
        check("single_h1", page.locator("h1").count() == 1, page.locator("h1").count())
        check(
            "ten_report_sections",
            page.locator("main > section").count() == 10,
            page.locator("main > section").count(),
        )
        check(
            "scope_ribbon_is_non_production",
            "NON-PRODUCTION" in page.locator(".scope-ribbon").inner_text(),
            page.locator(".scope-ribbon").inner_text(),
        )
        check(
            "bounded_verdict_visible",
            page.get_by_text("PASS FOR BOUNDED DUAL PILOT", exact=True).count() == 1,
            page.get_by_text("PASS FOR BOUNDED DUAL PILOT", exact=True).count(),
        )
        evidence = json.loads(
            page.locator("#report-evidence").text_content() or "{}"
        )
        check(
            "embedded_oracle_count",
            evidence.get("oracle_cases") == 27355,
            evidence.get("oracle_cases"),
        )
        check(
            "embedded_objective_certificate",
            evidence.get("objective_status") == "PROVEN_OPTIMAL",
            evidence.get("objective_status"),
        )
        check(
            "embedded_family_certificate",
            evidence.get("family_status") == "COMPLETE",
            evidence.get("family_status"),
        )

        link_results = page.evaluate(
            """() => Array.from(document.querySelectorAll('nav.toc a')).map(a => ({
                href: a.getAttribute('href'),
                targetExists: Boolean(document.querySelector(a.getAttribute('href')))
            }))"""
        )
        check(
            "toc_targets",
            len(link_results) == 10
            and all(result["targetExists"] for result in link_results),
            link_results,
        )

        initial_theme = page.locator("html").get_attribute("data-theme")
        page.locator("#themeToggle").click()
        toggled_theme = page.locator("html").get_attribute("data-theme")
        check(
            "theme_toggle",
            initial_theme == "paper" and toggled_theme == "night",
            {"initial": initial_theme, "after_click": toggled_theme},
        )
        page.locator("#themeToggle").click()

        page.locator("#collapseAll").click()
        collapsed_count = page.locator("details[open]").count()
        page.locator("#expandAll").click()
        expanded_count = page.locator("details[open]").count()
        total_details = page.locator("details").count()
        check(
            "details_controls",
            collapsed_count == 0 and expanded_count == total_details,
            {
                "collapsed_open_count": collapsed_count,
                "expanded_open_count": expanded_count,
                "total_details": total_details,
            },
        )
        page.evaluate(
            """() => {
                document.documentElement.style.scrollBehavior = "auto";
                document.body.style.scrollBehavior = "auto";
                if (document.activeElement) document.activeElement.blur();
                window.scrollTo(0, 0);
            }"""
        )
        page.keyboard.press("Home")
        page.wait_for_timeout(100)
        check("desktop_screenshot_at_top", page.evaluate("window.scrollY") == 0, page.evaluate("window.scrollY"))
        page.screenshot(path=str(args.screenshot), full_page=False)

        desktop_overflow = page.evaluate(
            "() => ({body: document.body.scrollWidth, viewport: window.innerWidth})"
        )
        check(
            "desktop_no_body_overflow",
            desktop_overflow["body"] <= desktop_overflow["viewport"] + 1,
            desktop_overflow,
        )
        check(
            "substantial_rendered_content",
            len(page.locator("main").inner_text()) >= 5000,
            len(page.locator("main").inner_text()),
        )

        page.set_viewport_size({"width": 390, "height": 844})
        page.reload(wait_until="networkidle")
        mobile_overflow = page.evaluate(
            "() => ({body: document.body.scrollWidth, viewport: window.innerWidth})"
        )
        check(
            "mobile_no_body_overflow",
            mobile_overflow["body"] <= mobile_overflow["viewport"] + 1,
            mobile_overflow,
        )
        check(
            "mobile_scope_ribbon_visible",
            page.locator(".scope-ribbon").is_visible(),
            page.locator(".scope-ribbon").bounding_box(),
        )
        browser.close()

    check("console_errors", not console_errors, console_errors)
    check("page_errors", not page_errors, page_errors)
    receipt = {
        "schema": "intersubmod.optimized_backend_html_browser_qa.v1",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "html_path": str(html_path),
        "html_sha256": sha256_file(html_path),
        "browser": {"engine": "chromium", "version": browser_version},
        "screenshot": {
            "path": str(args.screenshot.resolve()),
            "sha256": sha256_file(args.screenshot),
            "size_bytes": args.screenshot.stat().st_size,
        },
        "checks": checks,
        "summary": {
            "total": len(checks),
            "passed": sum(row["status"] == "PASS" for row in checks),
            "failed": sum(row["status"] == "FAIL" for row in checks),
            "all_pass": all(row["status"] == "PASS" for row in checks),
        },
    }
    args.receipt.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt["summary"], ensure_ascii=False, sort_keys=True))
    return 0 if receipt["summary"]["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
