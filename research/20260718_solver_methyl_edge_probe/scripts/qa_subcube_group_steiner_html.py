#!/usr/bin/env python3
"""Browser QA for the standalone Subcube Group-Steiner method explainer."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


VIEWPORTS = (
    (1440, 1000),
    (1024, 768),
    (390, 844),
    (320, 568),
)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} HTML_PATH")
    html_path = Path(sys.argv[1]).resolve()
    if not html_path.is_file():
        raise SystemExit(f"missing HTML: {html_path}")

    results: dict[str, object] = {
        "html": str(html_path),
        "viewports": [],
        "console_errors": [],
        "page_errors": [],
        "no_js": {},
        "print": {},
    }

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for width, height in VIEWPORTS:
            page = browser.new_page(viewport={"width": width, "height": height})
            page.on(
                "console",
                lambda message: (
                    results["console_errors"].append(message.text)
                    if message.type == "error"
                    else None
                ),
            )
            page.on("pageerror", lambda error: results["page_errors"].append(str(error)))
            page.goto(html_path.as_uri(), wait_until="networkidle")
            metrics = page.evaluate(
                """() => ({
                    viewportWidth: window.innerWidth,
                    documentWidth: document.documentElement.scrollWidth,
                    bodyWidth: document.body.scrollWidth,
                    title: document.title,
                    h1: document.querySelectorAll('h1').length,
                    h2: document.querySelectorAll('h2').length,
                    svg: document.querySelectorAll('svg').length,
                    details: document.querySelectorAll('details').length,
                    tables: document.querySelectorAll('table').length,
                    externalAssets: document.querySelectorAll(
                      'script[src], link[rel="stylesheet"], img[src], iframe[src]'
                    ).length,
                    mainTextLength: (document.querySelector('main')?.innerText || '').length,
                    horizontalOverflow:
                      document.documentElement.scrollWidth > window.innerWidth + 1
                })"""
            )
            metrics["screenshot"] = f"/tmp/subcube_group_steiner_{width}.png"
            page.screenshot(path=metrics["screenshot"], full_page=False)
            results["viewports"].append(metrics)
            page.close()

        no_js_context = browser.new_context(
            java_script_enabled=False,
            viewport={"width": 390, "height": 844},
        )
        no_js_page = no_js_context.new_page()
        no_js_page.goto(html_path.as_uri(), wait_until="load")
        results["no_js"] = no_js_page.evaluate(
            """() => ({
                title: document.title,
                mainTextLength: (document.querySelector('main')?.innerText || '').length,
                svg: document.querySelectorAll('svg').length,
                horizontalOverflow:
                  document.documentElement.scrollWidth > window.innerWidth + 1
            })"""
        )
        no_js_context.close()

        print_page = browser.new_page(viewport={"width": 1240, "height": 900})
        print_page.goto(html_path.as_uri(), wait_until="networkidle")
        print_page.emulate_media(media="print")
        results["print"] = print_page.evaluate(
            """() => {
                const detailBodies = [...document.querySelectorAll('details > .details-body')];
                return {
                  visibleDetailBodies: detailBodies.filter(
                    node => getComputedStyle(node).display !== 'none'
                  ).length,
                  totalDetailBodies: detailBodies.length,
                  bodyBackground: getComputedStyle(document.body).backgroundColor,
                  navDisplay: getComputedStyle(document.querySelector('nav')).display
                };
            }"""
        )
        print_page.pdf(
            path="/tmp/subcube_group_steiner_print.pdf",
            format="A4",
            print_background=True,
        )
        results["print"]["pdf"] = "/tmp/subcube_group_steiner_print.pdf"
        print_page.close()
        browser.close()

    passed = (
        not results["console_errors"]
        and not results["page_errors"]
        and all(not item["horizontalOverflow"] for item in results["viewports"])
        and all(item["externalAssets"] == 0 for item in results["viewports"])
        and results["no_js"]["mainTextLength"] > 5_000
        and not results["no_js"]["horizontalOverflow"]
        and results["print"]["visibleDetailBodies"] == results["print"]["totalDetailBodies"]
        and results["print"]["navDisplay"] == "none"
    )
    results["pass"] = passed
    print(json.dumps(results, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
