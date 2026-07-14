#!/usr/bin/env python3
"""Capture the archived reference and current canonical workstation before redesign."""

import json
from pathlib import Path

from playwright.sync_api import sync_playwright


ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
SHOTS = OUT / "screenshots"
OLD = ROOT / "docs/methodology/_assets/topology_workstation/HCC1395.html"
NEW = ROOT / "docs/methodology/_assets/layered_workstation/HCC1395.html"


def capture_surface(browser, label, path, viewport, targets):
    context = browser.new_context(viewport=viewport, device_scale_factor=1)
    page = context.new_page()
    console_errors = []
    page_errors = []
    page.on("console", lambda message: console_errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.goto(path.as_uri(), wait_until="load", timeout=120_000)
    page.wait_for_load_state("networkidle", timeout=120_000)
    page.wait_for_timeout(300)
    records = []
    for name, selector in targets:
        locator = page.locator(selector).first
        locator.wait_for(state="visible", timeout=30_000)
        locator.scroll_into_view_if_needed()
        page.wait_for_timeout(150)
        target = SHOTS / f"{label}_{name}.png"
        locator.screenshot(path=str(target), animations="disabled")
        records.append(
            {
                "name": name,
                "selector": selector,
                "path": str(target.relative_to(ROOT)),
                "box": locator.bounding_box(),
            }
        )
    overflow = page.evaluate(
        "Math.max(document.body.scrollWidth, document.documentElement.scrollWidth) - window.innerWidth"
    )
    context.close()
    return {
        "label": label,
        "source": str(path.relative_to(ROOT)),
        "viewport": viewport,
        "captures": records,
        "body_overflow_px": overflow,
        "console_errors": console_errors,
        "page_errors": page_errors,
    }


def main():
    SHOTS.mkdir(parents=True, exist_ok=True)
    runs = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        runs.append(
            capture_surface(
                browser,
                "old_desktop",
                OLD,
                {"width": 1440, "height": 1000},
                [
                    ("overview", ".archive-banner"),
                    ("statistics", ".stats"),
                    ("three_facets", ".facets"),
                ],
            )
        )
        runs.append(
            capture_surface(
                browser,
                "new_desktop",
                NEW,
                {"width": 1440, "height": 1000},
                [
                    ("overview", ".hero"),
                    ("genome_and_concepts", "#genome-overview"),
                    ("filters", "#region-browser .filters"),
                ],
            )
        )
        runs.append(
            capture_surface(
                browser,
                "old_mobile",
                OLD,
                {"width": 390, "height": 844},
                [("three_facets", ".facets")],
            )
        )
        runs.append(
            capture_surface(
                browser,
                "new_mobile",
                NEW,
                {"width": 390, "height": 844},
                [
                    ("genome_and_concepts", "#genome-overview"),
                    ("filters", "#region-browser .filters"),
                ],
            )
        )
        version = browser.version
        browser.close()

    metrics = {
        "chromium_version": version,
        "old_reference_is_archived": True,
        "runs": runs,
        "fatal_errors": sum(len(run["page_errors"]) for run in runs),
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
