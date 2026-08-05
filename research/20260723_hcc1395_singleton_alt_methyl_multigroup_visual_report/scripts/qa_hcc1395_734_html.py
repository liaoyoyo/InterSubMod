#!/usr/bin/env python3
"""Independent Playwright QA for the self-contained HCC1395 HTML report."""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import sync_playwright


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
WORK = REPO / "research/20260723_hcc1395_singleton_alt_methyl_multigroup_visual_report"
ARTIFACT = WORK / "artifact.json"
HTML = WORK / "20260723_HCC1395_734穩定甲基多群統計與視覺報告_01.standalone.html"
RESULTS = WORK / "results"
RECEIPT = RESULTS / "20260723_HCC1395_734位點HTML_Playwright_QA_01.json"
TOP_SCREENSHOT = RESULTS / "03_HCC1395_734_HTML_top_preview.png"
CASE_SCREENSHOT = RESULTS / "04_HCC1395_chr14_UPGMA_preview.png"
CHROMIUM = Path("/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome")


def expected_counts(artifact: dict) -> dict[str, int]:
    blocks = artifact["manifest"]["blocks"]
    metric_cards = sum(
        len(block.get("cardIds", [])) for block in blocks if block["type"] == "metric-strip"
    )
    return {
        "blocks": sum(block["type"] != "metric-strip" for block in blocks) + metric_cards,
        "charts": sum(block["type"] == "chart" for block in blocks),
        "html": sum(block["type"] == "html" for block in blocks),
        "metrics": metric_cards,
        "tables": sum(block["type"] == "table" for block in blocks),
    }


def qa_viewport(browser, artifact: dict, viewport: dict, capture: bool) -> dict:
    console_errors: list[str] = []
    page_errors: list[str] = []
    external_requests: list[str] = []
    context = browser.new_context(
        viewport={"width": viewport["width"], "height": viewport["height"]},
        color_scheme="light",
        reduced_motion="reduce",
    )
    page = context.new_page()
    page.set_default_timeout(60_000)
    page.on("console", lambda message: console_errors.append(message.text) if message.type == "error" else None)
    page.on("pageerror", lambda error: page_errors.append(str(error)))

    def capture_request(request) -> None:
        url = request.url
        if url.startswith(("http://", "https://", "ws://", "wss://")):
            external_requests.append(url)

    page.on("request", capture_request)
    page.goto(HTML.as_uri(), wait_until="load", timeout=120_000)
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        state = page.evaluate("document.documentElement.dataset.dataAnalyticsPortableReader || ''")
        if state == "ready":
            break
        if state in {"failed", "missing-runtime", "unsupported"}:
            raise AssertionError(f"Portable reader entered terminal state: {state}")
        page.wait_for_timeout(100)
    else:
        raise TimeoutError("Portable reader did not reach ready within 120 seconds")
    page.wait_for_timeout(500)

    reader_state = page.evaluate("document.documentElement.dataset.dataAnalyticsPortableReader")
    title = page.locator("#data-analytics-portable-reader h1").first.text_content().strip()
    expected_title = artifact["manifest"]["title"]
    if title != expected_title:
        raise AssertionError(f"Title mismatch at {viewport['name']}: {title!r} != {expected_title!r}")

    actual_counts = {
        "blocks": page.locator("#data-analytics-portable-reader [data-analytics-layout-item]").count(),
        "charts": page.locator("#data-analytics-portable-reader .chart-frame").count(),
        "html": page.locator("#data-analytics-portable-reader iframe.report-html-frame").count(),
        "metrics": page.locator("#data-analytics-portable-reader .report-metric-card").count(),
        "tables": page.locator("#data-analytics-portable-reader .table-card").count(),
    }
    expected = expected_counts(artifact)
    if actual_counts != expected:
        raise AssertionError(f"Rendered counts differ at {viewport['name']}: {actual_counts} != {expected}")

    geometry = page.evaluate(
        """
        () => {
          const selectors = ['[data-analytics-layout-item]', '.chart-frame', '.table-card', '.report-metric-card'];
          const bad = [];
          for (const selector of selectors) {
            for (const node of document.querySelectorAll('#data-analytics-portable-reader ' + selector)) {
              const rect = node.getBoundingClientRect();
              if (!(rect.width > 0 && rect.height > 0)) bad.push({selector, width: rect.width, height: rect.height});
            }
          }
          return {
            bad,
            clientWidth: document.documentElement.clientWidth,
            scrollWidth: document.documentElement.scrollWidth,
          };
        }
        """
    )
    if geometry["bad"]:
        raise AssertionError(f"Zero-size content at {viewport['name']}: {geometry['bad'][:5]}")
    if geometry["scrollWidth"] > geometry["clientWidth"] + 1:
        raise AssertionError(f"Horizontal page overflow at {viewport['name']}: {geometry}")

    chart_svgs = page.locator("#data-analytics-portable-reader .chart-frame svg.recharts-surface").count()
    chart_render_details = page.locator("#data-analytics-portable-reader .chart-frame").evaluate_all(
        """
        frames => frames.map(frame => {
          const plot = frame.querySelector('.chart-plot');
          const rect = plot?.getBoundingClientRect();
          return {
            childElements: plot?.childElementCount ?? 0,
            height: rect?.height ?? 0,
            svgs: plot?.querySelectorAll('svg').length ?? 0,
            width: rect?.width ?? 0,
          };
        })
        """
    )
    if any(
        detail["childElements"] < 1 or detail["width"] <= 0 or detail["height"] <= 0
        for detail in chart_render_details
    ):
        raise AssertionError(f"A chart plot did not render at {viewport['name']}: {chart_render_details}")

    iframe_locator = page.locator("iframe.report-html-frame")
    image_checks = []
    for index in range(iframe_locator.count()):
        frame = iframe_locator.nth(index).content_frame
        frame.locator("img").wait_for(state="visible", timeout=60_000)
        status = frame.locator("img").evaluate(
            "img => ({complete: img.complete, naturalWidth: img.naturalWidth, naturalHeight: img.naturalHeight, alt: img.alt})"
        )
        if not status["complete"] or status["naturalWidth"] < 1 or status["naturalHeight"] < 1:
            raise AssertionError(f"Embedded figure {index} did not load: {status}")
        image_checks.append(status)

    source_interaction = "not_run"
    if viewport["name"] == "desktop":
        source_button = page.locator(
            '#data-analytics-portable-reader button[data-artifact-action="open-options"]'
            '[data-artifact-has-source="true"]'
        ).first
        source_button.scroll_into_view_if_needed()
        source_button.evaluate("button => button.click()")
        source_action = page.locator('[role="menu"] [data-artifact-action="view-source"]').first
        source_action.wait_for(state="visible")
        source_action.evaluate("button => button.click()")
        dialog = page.locator('[data-artifact-dialog="source"]').first
        dialog.wait_for(state="visible")
        overview = dialog.get_by_role("tab", name="Overview", exact=True)
        if overview.count() != 1:
            raise AssertionError("Source dialog lacks its Overview tab")
        page.keyboard.press("Escape")
        source_interaction = "passed"

    if capture:
        page.evaluate("window.scrollTo(0, 0)")
        page.wait_for_timeout(250)
        page.screenshot(path=str(TOP_SCREENSHOT), full_page=False)
        iframe_locator.first.evaluate("node => node.scrollIntoView({block: 'start'})")
        page.wait_for_timeout(500)
        page.screenshot(path=str(CASE_SCREENSHOT), full_page=False)

    if external_requests:
        raise AssertionError(f"External requests at {viewport['name']}: {external_requests}")
    if console_errors or page_errors:
        raise AssertionError(
            f"Browser errors at {viewport['name']}: console={console_errors}, page={page_errors}"
        )

    result = {
        "name": viewport["name"],
        "viewport": {"width": viewport["width"], "height": viewport["height"]},
        "reader_state": reader_state,
        "title": title,
        "counts": actual_counts,
        "chart_svgs": chart_svgs,
        "chart_render_details": chart_render_details,
        "embedded_images": image_checks,
        "horizontal_overflow": geometry,
        "source_interaction": source_interaction,
        "external_requests": external_requests,
        "console_errors": console_errors,
        "page_errors": page_errors,
    }
    context.close()
    return result


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    if not HTML.exists() or not CHROMIUM.exists():
        raise FileNotFoundError("HTML or Chromium is missing")
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True,
            executable_path=str(CHROMIUM),
            args=["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"],
        )
        try:
            viewports = [
                qa_viewport(browser, artifact, {"name": "desktop", "width": 1440, "height": 1000}, True),
                qa_viewport(browser, artifact, {"name": "mobile", "width": 390, "height": 844}, False),
            ]
        finally:
            browser.close()

    receipt = {
        "ok": True,
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "html": str(HTML),
        "artifact": str(ARTIFACT),
        "browser": str(CHROMIUM),
        "viewports": viewports,
        "screenshots": [str(TOP_SCREENSHOT), str(CASE_SCREENSHOT)],
        "note": (
            "Canonical portable builder passed validation/package/structural verification. Its dump-dom chart extraction "
            "timed out in this host; this independent Playwright run verifies the actual dynamic reader, charts, tables, "
            "sandboxed figures, source dialog, errors, external requests, and responsive geometry."
        ),
    }
    RECEIPT.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
