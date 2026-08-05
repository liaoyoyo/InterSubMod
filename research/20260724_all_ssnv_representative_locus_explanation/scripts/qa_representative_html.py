#!/usr/bin/env python3
"""Playwright QA for the representative-locus explanatory HTML."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from playwright.sync_api import sync_playwright


REPO = Path("/big7_disk/liaoyoyo2001/InterSubMod")
WORK = REPO / "research/20260724_all_ssnv_representative_locus_explanation"
HTML = WORK / "20260724_單位點sSNV甲基多群代表位點定義與停止判定_01.standalone.html"
RESULTS = WORK / "results"
RECEIPT = RESULTS / "20260724_HTML_Playwright_QA_02.json"
TOP_SCREENSHOT = RESULTS / "04_report_top_desktop_v2.png"
CHR14_SCREENSHOT = RESULTS / "05_chr14_figures_desktop_v2.png"
MOBILE_SCREENSHOT = RESULTS / "06_report_top_mobile_v2.png"
CHROMIUM = Path("/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def qa_viewport(browser, name: str, width: int, height: int, capture: bool) -> dict:
    console_errors: list[str] = []
    page_errors: list[str] = []
    external_requests: list[str] = []
    context = browser.new_context(
        viewport={"width": width, "height": height},
        color_scheme="light",
        reduced_motion="reduce",
    )
    page = context.new_page()
    page.on(
        "console",
        lambda message: console_errors.append(message.text) if message.type == "error" else None,
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))

    def capture_request(request) -> None:
        if request.url.startswith(("http://", "https://", "ws://", "wss://")):
            external_requests.append(request.url)

    page.on("request", capture_request)
    page.set_default_timeout(120_000)
    page.goto(HTML.as_uri(), wait_until="load")
    image_locator = page.locator("main img")
    if image_locator.count() != 6:
        raise AssertionError(f"Expected 6 images, found {image_locator.count()}")
    for index in range(image_locator.count()):
        image = image_locator.nth(index)
        image.scroll_into_view_if_needed()
        image.evaluate("img => { img.loading = 'eager'; }")
        page.wait_for_function(
            "img => img.complete && img.naturalWidth > 0 && img.naturalHeight > 0",
            arg=image.element_handle(),
        )

    title = page.locator("h1").inner_text().strip()
    expected_title = "單位點 sSNV 甲基多群：代表位點定義、圖解與停止判定"
    if title != expected_title:
        raise AssertionError(f"Title mismatch: {title!r}")

    counts = {
        "h2_sections": page.locator("main h2").count(),
        "toc_links": page.locator(".toc a").count(),
        "images": page.locator("main img").count(),
        "tables": page.locator("main table").count(),
        "stat_boxes": page.locator(".stat-box").count(),
        "locus_cards": page.locator(".locus-card").count(),
        "inline_svgs": page.locator("main svg").count(),
    }
    expected_counts = {
        "h2_sections": 10,
        "toc_links": 10,
        "images": 6,
        "tables": 3,
        "stat_boxes": 8,
        "locus_cards": 6,
        "inline_svgs": 1,
    }
    if counts != expected_counts:
        raise AssertionError(f"Structural counts differ: {counts} != {expected_counts}")

    geometry = page.evaluate(
        """
        () => {
          const root = document.documentElement;
          const bad = [];
          for (const node of document.querySelectorAll(
            'main, article, .stat-box, .locus-card, figure, .figure-frame, .evidence-step, svg'
          )) {
            const rect = node.getBoundingClientRect();
            if (!(rect.width > 0 && rect.height > 0)) {
              bad.push({tag: node.tagName, cls: node.className, width: rect.width, height: rect.height});
            }
          }
          return {
            bad,
            clientWidth: root.clientWidth,
            scrollWidth: root.scrollWidth,
          };
        }
        """
    )
    if geometry["bad"]:
        raise AssertionError(f"Zero-size content: {geometry['bad']}")
    if geometry["scrollWidth"] > geometry["clientWidth"] + 1:
        raise AssertionError(f"Page-level horizontal overflow: {geometry}")

    images = page.locator("main img").evaluate_all(
        """
        imgs => imgs.map(img => ({
          alt: img.alt,
          complete: img.complete,
          naturalWidth: img.naturalWidth,
          naturalHeight: img.naturalHeight,
          renderedWidth: img.getBoundingClientRect().width,
          renderedHeight: img.getBoundingClientRect().height,
        }))
        """
    )
    for item in images:
        if (
            not item["alt"]
            or not item["complete"]
            or item["naturalWidth"] < 1
            or item["naturalHeight"] < 1
            or item["renderedWidth"] < 1
            or item["renderedHeight"] < 1
        ):
            raise AssertionError(f"Broken or inaccessible image: {item}")

    flow = page.locator("svg.flow").evaluate(
        """
        svg => {
          const rect = svg.getBoundingClientRect();
          return {
            width: rect.width,
            height: rect.height,
            title: svg.querySelector('title')?.textContent || '',
            desc: svg.querySelector('desc')?.textContent || '',
          };
        }
        """
    )
    if flow["width"] < 1 or flow["height"] < 1 or not flow["title"] or not flow["desc"]:
        raise AssertionError(f"Accessible flow diagram did not render: {flow}")

    page.locator('a[href="#chr14"]').click()
    page.wait_for_timeout(150)
    if page.evaluate("location.hash") != "#chr14":
        raise AssertionError("TOC anchor navigation failed")

    if capture and name == "desktop":
        page.evaluate("window.scrollTo(0, 0)")
        page.wait_for_timeout(150)
        page.screenshot(path=str(TOP_SCREENSHOT), full_page=False)
        page.locator("#chr14").scroll_into_view_if_needed()
        page.wait_for_timeout(150)
        page.screenshot(path=str(CHR14_SCREENSHOT), full_page=False)
    if capture and name == "mobile":
        page.evaluate("window.scrollTo(0, 0)")
        page.wait_for_timeout(150)
        page.screenshot(path=str(MOBILE_SCREENSHOT), full_page=False)

    if external_requests:
        raise AssertionError(f"Unexpected external requests: {external_requests}")
    if console_errors or page_errors:
        raise AssertionError(
            f"Browser errors: console={console_errors}, page={page_errors}"
        )

    result = {
        "name": name,
        "viewport": {"width": width, "height": height},
        "title": title,
        "counts": counts,
        "geometry": geometry,
        "images": images,
        "flow": flow,
        "anchor_navigation": "passed",
        "external_requests": external_requests,
        "console_errors": console_errors,
        "page_errors": page_errors,
    }
    context.close()
    return result


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not HTML.is_file() or not CHROMIUM.is_file():
        raise FileNotFoundError("HTML or Chromium executable is missing")

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            headless=True,
            executable_path=str(CHROMIUM),
            args=["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"],
        )
        try:
            viewports = [
                qa_viewport(browser, "desktop", 1440, 1000, True),
                qa_viewport(browser, "mobile", 390, 844, True),
            ]
        finally:
            browser.close()

    receipt = {
        "ok": True,
        "created_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "html": str(HTML),
        "html_sha256": sha256(HTML),
        "browser": str(CHROMIUM),
        "viewports": viewports,
        "screenshots": [
            str(TOP_SCREENSHOT),
            str(CHR14_SCREENSHOT),
            str(MOBILE_SCREENSHOT),
        ],
        "claim_boundary": (
            "Representative visuals support stable focal-ALT read-level methylation "
            "substructure candidates; they do not establish cellular clones or lineage."
        ),
    }
    RECEIPT.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
