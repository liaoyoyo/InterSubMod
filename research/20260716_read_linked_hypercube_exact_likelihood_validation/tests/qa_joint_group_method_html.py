#!/usr/bin/env python3
"""Browser QA for the joint-group method teaching HTML."""

from __future__ import annotations

import json
from pathlib import Path

from playwright.sync_api import Browser, Page, sync_playwright


ROOT = Path(__file__).resolve().parents[1]
HTML = ROOT / "20260717_聯合GroupConstraints到Likelihood排序方法教學_01.html"
OUT = ROOT / "results" / "20260717_joint_group_method_html_qa"


def assert_no_document_overflow(page: Page, label: str) -> dict[str, int]:
    metrics = page.evaluate(
        """() => ({
          documentScrollWidth: document.documentElement.scrollWidth,
          documentClientWidth: document.documentElement.clientWidth,
          bodyScrollWidth: document.body.scrollWidth,
          bodyClientWidth: document.body.clientWidth
        })"""
    )
    if metrics["documentScrollWidth"] > metrics["documentClientWidth"] + 1:
        raise AssertionError(f"{label}: document horizontal overflow: {metrics}")
    return metrics


def open_checked_page(browser: Browser, viewport: dict[str, int]) -> tuple[Page, list[str], list[str]]:
    page = browser.new_page(viewport=viewport)
    console_errors: list[str] = []
    page_errors: list[str] = []
    page.on(
        "console",
        lambda message: console_errors.append(message.text) if message.type == "error" else None,
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.goto(HTML.as_uri(), wait_until="load")
    return page, console_errors, page_errors


def main() -> None:
    if not HTML.is_file():
        raise FileNotFoundError(HTML)
    OUT.mkdir(parents=True, exist_ok=True)

    receipt: dict[str, object] = {
        "input_html": str(HTML),
        "scope": "Task Type F — method teaching demo; not full-data validation evidence",
        "checks": {},
        "artifacts": {},
    }

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)

        desktop, desktop_console, desktop_page_errors = open_checked_page(
            browser, {"width": 1440, "height": 1000}
        )
        if desktop.title() != "聯合 Group Constraints → Minimum-extra 候選 → Likelihood 排序｜方法教學":
            raise AssertionError(f"Unexpected title: {desktop.title()}")
        if desktop.locator("html").get_attribute("lang") != "zh-Hant":
            raise AssertionError("Document language must be zh-Hant")
        if desktop.locator("h1").count() != 1:
            raise AssertionError("Expected exactly one h1")
        if desktop.locator(".status-strip").count() != 1:
            raise AssertionError("Missing visible scope/status strip")
        status_text = desktop.locator(".status-strip").inner_text()
        for required in ("audited", "Demo", "NO-GO", "Not full-data"):
            if required.lower() not in status_text.lower():
                raise AssertionError(f"Missing status token: {required}")
        if desktop.locator(".evidence-card").count() != 4:
            raise AssertionError("Expected four evidence-state cards")
        if desktop.locator("svg[role='img'] title").count() < 1:
            raise AssertionError("SVG must have a title")
        if desktop.locator("svg[role='img'] desc").count() < 1:
            raise AssertionError("SVG must have a description")
        if desktop.locator("table caption").count() != desktop.locator("table").count():
            raise AssertionError("Every data table must have a caption")
        if desktop.locator("thead th[scope='col']").count() != desktop.locator("thead th").count():
            raise AssertionError("Every table header must declare column scope")

        desktop_metrics = assert_no_document_overflow(desktop, "desktop-1440")
        details_count = desktop.locator("details").count()
        if details_count < 5:
            raise AssertionError(f"Expected collapsible technical notes, found {details_count}")
        desktop.locator("#toggle-details").click()
        open_count = desktop.locator("details[open]").count()
        if open_count != details_count:
            raise AssertionError(f"Expand-all failed: {open_count}/{details_count}")
        desktop.locator("#toggle-details").click()
        if desktop.locator("details[open]").count() != 0:
            raise AssertionError("Collapse-all failed")

        desktop.locator(".topnav a[href='#step2']").click()
        desktop.wait_for_timeout(200)
        if desktop.locator("#step2-title").count() != 1:
            raise AssertionError("Step 2 anchor target is missing")
        external_resources = desktop.evaluate(
            """() => performance.getEntriesByType('resource')
              .map((entry) => entry.name)
              .filter((name) => /^https?:/i.test(name))"""
        )
        if external_resources:
            raise AssertionError(f"Unexpected external resources: {external_resources}")

        desktop.evaluate("window.scrollTo(0, 0)")
        desktop.wait_for_timeout(100)
        desktop_png = OUT / "desktop_1440.png"
        desktop.screenshot(path=str(desktop_png), full_page=True)

        desktop.emulate_media(media="print")
        print_detail_visibility = desktop.evaluate(
            """() => Array.from(document.querySelectorAll('details')).map((item) => {
              const body = item.querySelector('.detail-body');
              const summary = item.querySelector('summary');
              return {
                bodyHeight: body ? body.getBoundingClientRect().height : 0,
                summaryHeight: summary ? summary.getBoundingClientRect().height : 0
              };
            })"""
        )
        if any(item["bodyHeight"] <= 0 for item in print_detail_visibility):
            raise AssertionError(f"Print media did not reveal detail bodies: {print_detail_visibility}")
        print_pdf = OUT / "method_teaching_A4.pdf"
        desktop.pdf(
            path=str(print_pdf),
            format="A4",
            print_background=True,
            margin={"top": "12mm", "right": "11mm", "bottom": "14mm", "left": "11mm"},
        )
        desktop.close()

        viewport_results: dict[str, object] = {}
        for width, height in ((1024, 900), (390, 844), (320, 760)):
            page, console_errors, page_errors = open_checked_page(
                browser, {"width": width, "height": height}
            )
            metrics = assert_no_document_overflow(page, f"viewport-{width}")
            if not page.locator(".status-strip").is_visible():
                raise AssertionError(f"Status strip not visible at width {width}")
            if width <= 700 and page.locator(".evidence-card").count() != 4:
                raise AssertionError(f"Evidence cards missing at mobile width {width}")
            if width == 390:
                page.screenshot(path=str(OUT / "mobile_390.png"), full_page=True)
            viewport_results[str(width)] = {
                "metrics": metrics,
                "console_errors": console_errors,
                "page_errors": page_errors,
            }
            page.close()

        browser.close()

    if desktop_console or desktop_page_errors:
        raise AssertionError(
            f"Desktop runtime errors: console={desktop_console}, page={desktop_page_errors}"
        )
    for width, result in viewport_results.items():
        if result["console_errors"] or result["page_errors"]:
            raise AssertionError(f"Viewport {width} runtime errors: {result}")

    receipt["checks"] = {
        "desktop_1440_no_document_overflow": desktop_metrics,
        "responsive_no_document_overflow": viewport_results,
        "single_h1": True,
        "visible_demo_and_no_go_status": True,
        "four_evidence_state_cards": True,
        "svg_title_and_description": True,
        "table_caption_and_header_scopes": True,
        "external_resources": "0",
        "expand_collapse_all": f"PASS ({details_count} details)",
        "print_details_visible": True,
        "console_and_page_errors": "0",
    }
    receipt["artifacts"] = {
        "desktop_screenshot": str(OUT / "desktop_1440.png"),
        "mobile_screenshot": str(OUT / "mobile_390.png"),
        "print_pdf": str(OUT / "method_teaching_A4.pdf"),
    }
    receipt_path = OUT / "qa_receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
