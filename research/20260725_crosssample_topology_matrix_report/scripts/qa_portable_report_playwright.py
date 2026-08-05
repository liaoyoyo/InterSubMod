#!/usr/bin/env python3
"""Fallback Playwright QA for the self-contained portable HTML report."""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from playwright.sync_api import Page, sync_playwright


VIEWPORTS = (
    {"name": "desktop", "width": 1440, "height": 1000},
    {"name": "mobile", "width": 390, "height": 844},
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chromium-executable", required=True, type=Path)
    return parser.parse_args()


def wait_for_reader(page: Page) -> str:
    page.wait_for_load_state("networkidle")
    page.wait_for_function(
        """() => {
          const state = document.documentElement.dataset.dataAnalyticsPortableReader || "";
          return ["ready", "failed", "missing-runtime", "unsupported"].includes(state);
        }""",
        timeout=60_000,
    )
    return page.evaluate(
        "() => document.documentElement.dataset.dataAnalyticsPortableReader || ''"
    )


def inspect_page(
    page: Page, viewport: dict[str, Any], screenshot: Path, html_uri: str
) -> dict[str, Any]:
    page_errors: list[str] = []
    console_errors: list[str] = []
    external_requests: list[str] = []
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.on(
        "console",
        lambda message: (
            console_errors.append(message.text) if message.type == "error" else None
        ),
    )
    page.on(
        "request",
        lambda request: (
            external_requests.append(request.url)
            if re.match(r"^https?://", request.url)
            else None
        ),
    )
    page.goto(html_uri, wait_until="domcontentloaded")
    reader_state = wait_for_reader(page)
    page.screenshot(path=str(screenshot), full_page=True)
    layout = page.evaluate(
        """() => {
          const root = document.documentElement;
          const body = document.body;
          const charts = document.querySelectorAll(
            '#data-analytics-portable-reader section[data-artifact-kind="chart"]'
          );
          const chartSurfaces = document.querySelectorAll(
            '#data-analytics-portable-reader .chart-plot svg.recharts-surface'
          );
          const chartPlots = [
            ...document.querySelectorAll(
              '#data-analytics-portable-reader section[data-artifact-kind="chart"] .chart-plot'
            )
          ];
          const tables = document.querySelectorAll(
            '#data-analytics-portable-reader section[data-artifact-kind="table"]'
          );
          const blocks = document.querySelectorAll(
            '#data-analytics-portable-reader [data-artifact-block-id]'
          );
          const visibleDialogs = [...document.querySelectorAll('[role="dialog"]')]
            .filter((node) => {
              const style = getComputedStyle(node);
              const rect = node.getBoundingClientRect();
              return style.display !== 'none' && style.visibility !== 'hidden' &&
                rect.width > 0 && rect.height > 0;
            }).length;
          const sourceIds = document.querySelectorAll('[data-source-id]');
          const sourceHosts = document.querySelectorAll('[data-portable-source-host]');
          const sourceTooltips = document.querySelectorAll(
            '[data-portable-source-host][aria-describedby] + [role="tooltip"], ' +
            '[data-portable-source-host][aria-describedby] [role="tooltip"]'
          );
          const overflowRows = [...document.querySelectorAll('body *')]
            .map((node) => {
              const rect = node.getBoundingClientRect();
              return {
                tag: node.tagName.toLowerCase(),
                classes: String(node.className || '').slice(0, 180),
                artifactId:
                  node.getAttribute('data-artifact-id') ||
                  node.getAttribute('data-artifact-block-id') ||
                  node.getAttribute('data-chart-id') ||
                  node.getAttribute('data-table-id') ||
                  '',
                ariaLabel: node.getAttribute('aria-label') || '',
                clientWidth: node.clientWidth,
                scrollWidth: node.scrollWidth,
                left: Math.round(rect.left),
                right: Math.round(rect.right),
                width: Math.round(rect.width),
                insideScrollableTable: Boolean(node.closest('.table-wrap')),
              };
            })
            .filter((row) =>
              row.right > window.innerWidth + 1 ||
              row.left < -1 ||
              row.scrollWidth > row.clientWidth + 1
            )
            .sort((a, b) => Math.max(b.right - window.innerWidth, b.scrollWidth - b.clientWidth) -
              Math.max(a.right - window.innerWidth, a.scrollWidth - a.clientWidth))
          const overflowingElements = overflowRows.slice(0, 30);
          const outerOverflowingElements = overflowRows
            .filter((row) => !row.insideScrollableTable)
            .slice(0, 30);
          return {
            bodyScrollWidth: body.scrollWidth,
            chartCards: charts.length,
            chartPlots: chartPlots.length,
            chartPlotDimensions: chartPlots.map((node) => {
              const rect = node.getBoundingClientRect();
              return {
                width: Math.round(rect.width),
                height: Math.round(rect.height),
                childCount: node.children.length,
              };
            }),
            chartSurfaces: chartSurfaces.length,
            documentScrollWidth: root.scrollWidth,
            fallbackHidden:
              document.querySelector('#data-analytics-portable-fallback')?.hidden ?? null,
            readerAriaHidden:
              document.querySelector('#data-analytics-portable-reader')?.getAttribute('aria-hidden'),
            tableCards: tables.length,
            renderedBlocks: blocks.length,
            sourceIds: sourceIds.length,
            sourceHosts: sourceHosts.length,
            sourceTooltips: sourceTooltips.length,
            visibleDialogs,
            overflowingElements,
            outerOverflowingElements,
          };
        }"""
    )
    buttons = page.locator("button").all_inner_texts()
    return {
        "viewport": viewport,
        "reader_state": reader_state,
        "layout": layout,
        "button_count": len(buttons),
        "button_text_sample": buttons[:20],
        "external_requests": sorted(set(external_requests)),
        "page_errors": page_errors,
        "console_errors": console_errors,
        "screenshot": str(screenshot),
    }


def main() -> int:
    args = parse_args()
    if not args.html.is_file():
        raise FileNotFoundError(args.html)
    if not args.chromium_executable.is_file():
        raise FileNotFoundError(args.chromium_executable)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    results: list[dict[str, Any]] = []
    html_uri = args.html.resolve().as_uri()
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=str(args.chromium_executable),
            headless=True,
            args=["--no-sandbox", "--disable-gpu"],
        )
        try:
            for viewport in VIEWPORTS:
                context = browser.new_context(
                    viewport={"width": viewport["width"], "height": viewport["height"]},
                    color_scheme="light",
                    reduced_motion="reduce",
                )
                try:
                    page = context.new_page()
                    screenshot = args.output_dir / f"{viewport['name']}.png"
                    results.append(inspect_page(page, viewport, screenshot, html_uri))
                finally:
                    context.close()
        finally:
            browser.close()

    checks = {
        "all_reader_states_ready": all(row["reader_state"] == "ready" for row in results),
        "chart_card_count_5": all(row["layout"]["chartCards"] == 5 for row in results),
        "chart_plot_count_5": all(row["layout"]["chartPlots"] == 5 for row in results),
        "all_chart_plots_nonempty": all(
            all(
                plot["width"] > 0 and plot["height"] > 0 and plot["childCount"] > 0
                for plot in row["layout"]["chartPlotDimensions"]
            )
            for row in results
        ),
        "table_card_count_4": all(row["layout"]["tableCards"] == 4 for row in results),
        "rendered_block_count_24": all(row["layout"]["renderedBlocks"] == 24 for row in results),
        "no_document_overflow": all(
            row["layout"]["documentScrollWidth"] <= row["viewport"]["width"] + 1
            for row in results
        ),
        "no_external_requests": all(not row["external_requests"] for row in results),
        "no_page_errors": all(not row["page_errors"] for row in results),
        "no_console_errors": all(not row["console_errors"] for row in results),
        "source_traceability_embedded": all(
            row["layout"]["sourceIds"] >= 9
            and row["layout"]["sourceHosts"] >= 5
            and row["layout"]["sourceTooltips"] >= 5
            for row in results
        ),
    }
    receipt = {
        "schema_name": "intersubmod.portable_report_playwright_qa",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "html": str(args.html.resolve()),
        "checks": checks,
        "all_pass": all(checks.values()),
        "viewports": results,
        "note": (
            "Fallback QA after the official portable static-SVG Chromium dump timed out; "
            "the official structural verifier passed independently."
        ),
    }
    receipt_path = args.output_dir / "playwright_qa_receipt.json"
    receipt_path.write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
