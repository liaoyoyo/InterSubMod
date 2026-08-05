#!/usr/bin/env python3
"""Render and validate the standalone research report at desktop and mobile sizes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def validate_viewport(page, html_uri: str, name: str, width: int, height: int, output_dir: Path) -> dict:
    console_errors: list[str] = []
    page_errors: list[str] = []
    page.on("console", lambda message: console_errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.set_viewport_size({"width": width, "height": height})
    page.goto(html_uri, wait_until="load")
    page.wait_for_function("Array.from(document.images).every(img => img.complete)")

    image_state = page.locator("img").evaluate_all(
        "imgs => imgs.map(img => ({src: img.getAttribute('src'), width: img.naturalWidth, "
        "height: img.naturalHeight, complete: img.complete}))"
    )
    missing_images = [item for item in image_state if not item["complete"] or item["width"] <= 0 or item["height"] <= 0]
    layout = page.evaluate(
        """() => ({
          viewportWidth: window.innerWidth,
          documentWidth: document.documentElement.scrollWidth,
          bodyWidth: document.body.scrollWidth,
          documentHeight: document.documentElement.scrollHeight,
          overflow: Array.from(document.querySelectorAll('h1,h2,h3,.metric,.claim,.case,.badge,summary'))
            .filter(el => el.scrollWidth > el.clientWidth + 1)
            .map(el => ({tag: el.tagName, cls: el.className, text: el.textContent.trim().slice(0, 100),
                         scrollWidth: el.scrollWidth, clientWidth: el.clientWidth}))
        })"""
    )

    page.locator('a[href="#strict"]:visible').first.click()
    page.wait_for_function("location.hash === '#strict'")
    strict_target_ok = page.locator("#strict").count() == 1
    first_details = page.locator("details").first
    first_details.locator("summary").click()
    details_open = first_details.evaluate("el => el.open")

    page.goto(html_uri, wait_until="load")
    screenshot = output_dir / f"{name}_full.png"
    page.screenshot(path=str(screenshot), full_page=True)

    horizontal_overflow = max(layout["documentWidth"], layout["bodyWidth"]) > layout["viewportWidth"] + 1
    passed = not console_errors and not page_errors and not missing_images and not horizontal_overflow
    passed = passed and not layout["overflow"] and strict_target_ok and details_open
    return {
        "viewport": {"name": name, "width": width, "height": height},
        "layout": layout,
        "n_images": len(image_state),
        "missing_images": missing_images,
        "console_errors": console_errors,
        "page_errors": page_errors,
        "strict_target_ok": strict_target_ok,
        "details_open": details_open,
        "screenshot": str(screenshot),
        "pass": passed,
    }


def main() -> int:
    args = parse_args()
    html_path = args.html.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if not html_path.exists():
        raise FileNotFoundError(html_path)

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            results = [
                validate_viewport(browser.new_page(), html_path.as_uri(), "desktop_1440x1000", 1440, 1000, output_dir),
                validate_viewport(browser.new_page(), html_path.as_uri(), "mobile_390x844", 390, 844, output_dir),
            ]
        finally:
            browser.close()

    summary = {
        "schema_name": "intersubmod.standalone_html_browser_qa",
        "schema_version": "1.0.0",
        "html": str(html_path),
        "viewports": results,
        "pass": all(result["pass"] for result in results),
    }
    summary_path = output_dir / "browser_qa_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
