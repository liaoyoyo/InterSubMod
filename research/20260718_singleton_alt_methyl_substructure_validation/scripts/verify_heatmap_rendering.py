#!/usr/bin/env python3
"""Verify that all six portable HTML heatmaps render as colored matrices."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from playwright.sync_api import sync_playwright


HEADLESS_SHELL = (
    "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/"
    "chromium_headless_shell-1223/chrome-headless-shell-linux64/chrome-headless-shell"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_viewport(page: Any, html_path: Path, width: int, screenshot_dir: Path) -> dict[str, Any]:
    page.set_viewport_size({"width": width, "height": 1000 if width > 600 else 844})
    page.goto(html_path.resolve().as_uri(), wait_until="load", timeout=60_000)
    page.wait_for_function(
        "() => document.querySelectorAll('iframe.report-html-frame').length === 6",
        timeout=30_000,
    )
    page.wait_for_timeout(1_000)

    parent_geometry = page.evaluate(
        """() => ({
          clientWidth: document.documentElement.clientWidth,
          scrollWidth: document.documentElement.scrollWidth,
        })"""
    )
    frames = page.locator("iframe.report-html-frame")
    frame_results = []
    for index in range(frames.count()):
        locator = frames.nth(index)
        handle = locator.element_handle()
        frame = handle.content_frame()
        frame.wait_for_selector("[data-true-heatmap]", timeout=10_000)
        result = frame.evaluate(
            """() => {
              const root = document.querySelector("[data-true-heatmap]");
              const cells = Array.from(document.querySelectorAll("[data-heatmap-cell]"));
              const colors = [...new Set(cells.map((cell) => getComputedStyle(cell).backgroundColor))];
              const values = cells.map((cell) => cell.getAttribute("data-value"));
              const legend = document.querySelector(".hm-ramp");
              const table = document.querySelector(".hm-grid");
              const rect = table?.getBoundingClientRect();
              return {
                id: root?.getAttribute("data-true-heatmap"),
                title: root?.querySelector("h2")?.textContent?.trim(),
                cellCount: cells.length,
                finiteValueCount: values.filter((value) => value && value !== "nan").length,
                uniqueColorCount: colors.length,
                colors,
                groupALabels: document.querySelectorAll(".hm-group-a,.hm-row-a").length,
                groupBLabels: document.querySelectorAll(".hm-group-b,.hm-row-b").length,
                legendBackground: legend ? getComputedStyle(legend).backgroundImage : "",
                tableWidth: rect?.width || 0,
                tableHeight: rect?.height || 0,
                clientWidth: document.documentElement.clientWidth,
                scrollWidth: document.documentElement.scrollWidth,
                clientHeight: document.documentElement.clientHeight,
                scrollHeight: document.documentElement.scrollHeight,
              };
            }"""
        )
        result["pass"] = bool(
            result["id"]
            and result["cellCount"] == 64
            and result["finiteValueCount"] > 0
            and result["uniqueColorCount"] >= 3
            and result["groupALabels"] > 0
            and result["groupBLabels"] > 0
            and "gradient" in result["legendBackground"]
            and result["tableWidth"] > 0
            and result["tableHeight"] > 0
            and result["scrollWidth"] <= result["clientWidth"] + 1
            and result["scrollHeight"] <= result["clientHeight"] + 2
        )
        frame_results.append(result)
        if width == 1440:
            screenshot_dir.mkdir(parents=True, exist_ok=True)
            locator.screenshot(path=str(screenshot_dir / f"{index + 1:02d}_{result['id']}.png"))

    return {
        "viewport_width": width,
        "parent_geometry": parent_geometry,
        "parent_horizontal_overflow": parent_geometry["scrollWidth"] > parent_geometry["clientWidth"] + 1,
        "heatmaps": frame_results,
        "pass": (
            not (parent_geometry["scrollWidth"] > parent_geometry["clientWidth"] + 1)
            and len(frame_results) == 6
            and all(result["pass"] for result in frame_results)
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--screenshot-dir", type=Path, required=True)
    args = parser.parse_args()

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True, executable_path=HEADLESS_SHELL)
        page = browser.new_page(viewport={"width": 1440, "height": 1000}, device_scale_factor=1)
        results = [
            inspect_viewport(page, args.html, 1440, args.screenshot_dir),
            inspect_viewport(page, args.html, 390, args.screenshot_dir),
        ]
        browser.close()

    receipt = {
        "schema_name": "intersubmod.heatmap_rendering_qa",
        "schema_version": "1.0.0",
        "created_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "html": str(args.html.resolve()),
        "html_sha256": sha256_file(args.html),
        "expected_heatmaps": 6,
        "expected_cells_per_heatmap": 64,
        "viewports": results,
        "pass": all(result["pass"] for result in results),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "pass": receipt["pass"],
                "html_sha256": receipt["html_sha256"],
                "viewports": [
                    {
                        "width": result["viewport_width"],
                        "heatmaps": len(result["heatmaps"]),
                        "all_colored": all(item["uniqueColorCount"] >= 3 for item in result["heatmaps"]),
                        "pass": result["pass"],
                    }
                    for result in results
                ],
                "receipt": str(args.output.resolve()),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0 if receipt["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
