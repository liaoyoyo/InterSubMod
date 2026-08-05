#!/usr/bin/env python3
"""Browser QA for the standalone cross-HP observation report."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import sync_playwright


def inspect_page(page, url: str, viewport: dict[str, int], screenshot: Path) -> dict:
    console_errors: list[str] = []
    page_errors: list[str] = []
    page.on("console", lambda msg: console_errors.append(msg.text) if msg.type == "error" else None)
    page.on("pageerror", lambda exc: page_errors.append(str(exc)))
    page.set_viewport_size(viewport)
    response = page.goto(url, wait_until="networkidle")
    page.screenshot(path=str(screenshot), full_page=True)
    details = page.locator("details").first
    details.locator("summary").click()
    checks = {
        "http_status_200": response is not None and response.status == 200,
        "title": page.title() == "跨 HP clone-state 反推｜可行性與觀察紀錄",
        "h1_visible": page.locator("h1").is_visible(),
        "section_count_11": page.locator("main section").count() == 11,
        "details_opens": details.get_attribute("open") is not None,
        "key_51815_present": page.get_by_text("51,815", exact=True).count() >= 1,
        "key_22779_present": page.get_by_text("22,779", exact=True).count() >= 1,
        "flow_svg_visible": page.locator(".flow-svg").is_visible(),
        "candidate_cards_3": page.locator(".candidate-card").count() == 3,
        "no_console_errors": not console_errors,
        "no_page_errors": not page_errors,
        "no_external_subresources": page.evaluate(
            """() => performance.getEntriesByType('resource')
              .filter(x => !x.name.startsWith(location.origin)).length === 0"""
        ),
        "no_body_horizontal_overflow": page.evaluate(
            "() => document.documentElement.scrollWidth <= window.innerWidth + 1"
        ),
        "all_nav_targets_exist": page.evaluate(
            """() => [...document.querySelectorAll('nav a[href^="#"]')]
              .every(a => document.querySelector(a.getAttribute('href')))"""
        ),
    }
    return {
        "viewport": viewport,
        "checks": checks,
        "all_pass": all(checks.values()),
        "console_errors": console_errors,
        "page_errors": page_errors,
        "screenshot": str(screenshot.resolve()),
        "body_dimensions": page.evaluate(
            "() => ({scrollWidth:document.documentElement.scrollWidth,"
            "scrollHeight:document.documentElement.scrollHeight,"
            "innerWidth:window.innerWidth,innerHeight:window.innerHeight})"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page()
        desktop = inspect_page(
            page,
            args.url,
            {"width": 1440, "height": 1000},
            args.out_dir / "desktop_full.png",
        )
        mobile = inspect_page(
            page,
            args.url,
            {"width": 390, "height": 844},
            args.out_dir / "mobile_full.png",
        )
        page.emulate_media(media="print")
        print_checks = {
            "nav_hidden": page.locator("nav").evaluate(
                "(node) => getComputedStyle(node).display === 'none'"
            ),
            "hero_not_full_viewport": page.locator(".hero").evaluate(
                "(node) => node.getBoundingClientRect().height < 800"
            ),
        }
        browser.close()
    payload = {
        "schema_name": "intersubmod.cross_hp_observation_html_qa",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "url": args.url,
        "desktop": desktop,
        "mobile": mobile,
        "print_checks": print_checks,
        "all_pass": desktop["all_pass"] and mobile["all_pass"] and all(print_checks.values()),
    }
    receipt = args.out_dir / "browser_qa_receipt.json"
    receipt.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if not payload["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
