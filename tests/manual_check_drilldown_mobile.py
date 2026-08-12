#!/usr/bin/env python3
"""Browser smoke check for drilldown mobile/sticky/overflow contracts.

This is intentionally not auto-collected by pytest because Playwright is an
optional development dependency.  Run it against a freshly rendered bundle:

    python3 tests/manual_check_drilldown_mobile.py \
        --url http://127.0.0.1:8765/index.html --screenshot /tmp/drilldown-mobile.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--url", required=True)
    p.add_argument("--screenshot", type=Path)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    browser_errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 390, "height": 844})
        page.on("console", lambda msg: browser_errors.append(
            f"console:{msg.type}:{msg.text}") if msg.type == "error" else None)
        page.on("pageerror", lambda err: browser_errors.append(f"pageerror:{err}"))
        page.goto(args.url)
        page.wait_for_load_state("networkidle")

        # Check the real v1 status and descriptive cooccurrence behavior.
        page.locator("button[data-view='selfcheck']").click()
        selfcheck_heading = page.locator("#view-selfcheck h2").inner_text()
        if "未通過" not in selfcheck_heading:
            raise AssertionError(f"expected blocking v1 selfcheck heading: {selfcheck_heading!r}")
        page.locator("button[data-view='cooccur']").click()
        filter_cells = page.locator("#view-cooccur td[data-x]").count()
        if filter_cells:
            raise AssertionError(f"descriptive matrices exposed {filter_cells} fake filter controls")
        if args.screenshot:
            args.screenshot.parent.mkdir(parents=True, exist_ok=True)
            page.screenshot(path=str(args.screenshot), full_page=True)

        # Exercise sticky positioning after scroll, when the mobile authority bar wraps.
        page.evaluate("window.scrollTo(0, 700)")
        page.wait_for_timeout(100)
        sticky = page.evaluate("""() => {
            const a = document.querySelector('.authority').getBoundingClientRect();
            const l = document.querySelector('.levelbar').getBoundingClientRect();
            return {authorityBottom: a.bottom, levelbarTop: l.top};
        }""")
        if sticky["levelbarTop"] + 1 < sticky["authorityBottom"]:
            raise AssertionError(f"sticky bars overlap: {sticky}")

        # Synthetic panel uses the production CSS.  In fit mode neither the figure nor
        # body may overflow; at 8x only .panel-scroll may own horizontal overflow.
        page.evaluate("""() => {
            window.scrollTo(0, 0);
            const host = document.createElement('section');
            host.className = 'section';
            host.id = 'overflow-contract-fixture';
            host.innerHTML = `<div class="panel-scroll" data-mode="fit">
                <div class="panel-fig"><div style="width:100%;height:120px"></div></div>
            </div>`;
            document.querySelector('.shell').prepend(host);
        }""")
        fit = page.evaluate("""() => ({
            viewport: innerWidth,
            body: document.documentElement.scrollWidth,
            panelClient: document.querySelector('#overflow-contract-fixture .panel-scroll').clientWidth,
            panelScroll: document.querySelector('#overflow-contract-fixture .panel-scroll').scrollWidth
        })""")
        if fit["body"] > fit["viewport"]:
            raise AssertionError(f"fit mode leaks horizontal body overflow: {fit}")

        zoom = page.evaluate("""() => {
            const panel = document.querySelector('#overflow-contract-fixture .panel-scroll');
            const fig = panel.querySelector('.panel-fig');
            panel.dataset.mode = 'zoom';
            fig.style.width = '2400px';
            return {
                viewport: innerWidth,
                body: document.documentElement.scrollWidth,
                panelClient: panel.clientWidth,
                panelScroll: panel.scrollWidth
            };
        }""")
        if zoom["body"] > zoom["viewport"]:
            raise AssertionError(f"zoom leaks horizontal body overflow: {zoom}")
        if zoom["panelScroll"] <= zoom["panelClient"]:
            raise AssertionError(f"zoom did not create local panel overflow: {zoom}")

        browser.close()

    result = {
        "status": "PASS",
        "selfcheck_heading": selfcheck_heading,
        "descriptive_filter_cells": filter_cells,
        "sticky": sticky,
        "fit": fit,
        "zoom": zoom,
        "browser_errors": browser_errors,
        "screenshot": str(args.screenshot) if args.screenshot else None,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
