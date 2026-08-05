#!/usr/bin/env python3
"""Run deterministic desktop/mobile visual QA for the portable region report."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Mapping, Sequence


VIEWPORTS = (
    ("desktop", 1440, 1100),
    ("mobile", 390, 844),
)
EXPECTED_EVIDENCE_FILES = frozenset(
    {
        "desktop_overview.png",
        "mobile_overview.png",
        "desktop_method_flow.png",
        "desktop_distance_qc.png",
        "desktop_distance_cut_impact.png",
        "nojs_desktop_overview.png",
        "nojs_mobile_overview.png",
    }
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: JSON root must be an object")
    return value


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"expected a regular file: {path}")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def expected_visual_ids(
    artifact: Mapping[str, Any],
) -> tuple[list[str], list[str]]:
    """Return exact visible chart/table IDs after validating the manifest graph."""
    manifest = artifact.get("manifest")
    if not isinstance(manifest, Mapping):
        raise ValueError("artifact manifest is absent")
    blocks = manifest.get("blocks")
    charts = manifest.get("charts")
    tables = manifest.get("tables")
    if not all(isinstance(value, list) for value in (blocks, charts, tables)):
        raise ValueError("artifact blocks/charts/tables must be arrays")

    def require_unique_ids(
        values: Sequence[Any],
        *,
        key: str,
        label: str,
    ) -> list[str]:
        ids = [
            value.get(key) if isinstance(value, Mapping) else None
            for value in values
        ]
        if any(not isinstance(value, str) or not value for value in ids):
            raise ValueError(f"every {label} must declare a nonempty {key}")
        if len(set(ids)) != len(ids):
            raise ValueError(f"{label} {key} values must be unique")
        return sorted(ids)

    chart_ids = require_unique_ids(charts, key="id", label="manifest chart")
    table_ids = require_unique_ids(tables, key="id", label="manifest table")
    block_chart_ids = require_unique_ids(
        [
            block
            for block in blocks
            if isinstance(block, Mapping) and block.get("type") == "chart"
        ],
        key="chartId",
        label="visible chart block",
    )
    block_table_ids = require_unique_ids(
        [
            block
            for block in blocks
            if isinstance(block, Mapping) and block.get("type") == "table"
        ],
        key="tableId",
        label="visible table block",
    )
    if block_chart_ids != chart_ids:
        raise ValueError(
            "visible chart block IDs do not exactly match manifest chart IDs"
        )
    if block_table_ids != table_ids:
        raise ValueError(
            "visible table block IDs do not exactly match manifest table IDs"
        )
    return chart_ids, table_ids


def make_output_dir(path: Path) -> Path:
    if path.exists():
        if not path.is_dir() or next(path.iterdir(), None) is not None:
            raise ValueError(f"output directory must be new or empty: {path}")
    else:
        path.mkdir(parents=True, exist_ok=False)
    return path.resolve()


def write_json(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--artifact", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    from playwright.sync_api import sync_playwright

    args = parse_args(argv)
    html_path = args.html.resolve(strict=True)
    artifact_path = args.artifact.resolve(strict=True)
    output_dir = make_output_dir(args.output_dir)
    artifact = read_json(artifact_path)
    manifest = artifact["manifest"]
    expected_title = manifest["title"]
    expected_chart_ids, expected_table_ids = expected_visual_ids(artifact)
    expected_charts = len(expected_chart_ids)
    expected_tables = len(expected_table_ids)
    expected_svg_variants = expected_charts * 2
    results: list[dict[str, Any]] = []

    with sync_playwright() as playwright:
        browser_executable = Path(playwright.chromium.executable_path).resolve(
            strict=True
        )
        browser_identity = file_identity(browser_executable)
        browser = playwright.chromium.launch(
            headless=True,
            executable_path=str(browser_executable),
        )
        browser_version = browser.version
        for label, width, height in VIEWPORTS:
            context = browser.new_context(
                viewport={"width": width, "height": height},
                color_scheme="light",
                reduced_motion="reduce",
                device_scale_factor=1,
            )
            page = context.new_page()
            console_errors: list[str] = []
            page_errors: list[str] = []
            page.on(
                "console",
                lambda message: console_errors.append(message.text)
                if message.type == "error"
                else None,
            )
            page.on("pageerror", lambda error: page_errors.append(str(error)))
            page.goto(html_path.as_uri(), wait_until="networkidle")
            page.locator(
                'html[data-data-analytics-portable-reader="ready"]'
            ).wait_for(state="attached", timeout=30_000)
            page.wait_for_timeout(300)

            metrics = page.evaluate(
                """() => {
                  const root = document.documentElement;
                  const visible = (element) => {
                    const style = getComputedStyle(element);
                    const rect = element.getBoundingClientRect();
                    return style.display !== "none" && style.visibility !== "hidden" &&
                      rect.width > 0 && rect.height > 0;
                  };
                  const reader = document.querySelector(
                    "#data-analytics-portable-reader-root"
                  );
                  const chartPanels = [...reader.querySelectorAll(
                    '.chart-panel[data-artifact-id]'
                  )].filter(visible);
                  const tablePanels = [...reader.querySelectorAll(
                    '.table-card[data-artifact-id]'
                  )].filter(visible);
                  const staticChartPanels = [...document.querySelectorAll(
                    '#data-analytics-portable-fallback [data-artifact-kind="chart"]'
                  )];
                  const staticTablePanels = [...document.querySelectorAll(
                    '#data-analytics-portable-fallback [data-artifact-kind="table"]'
                  )];
                  const svgVariants = [...document.querySelectorAll(
                    '.portable-static-chart-svg'
                  )];
                  const brokenCharts = staticChartPanels.filter((panel) =>
                    panel.querySelectorAll('.portable-static-chart-svg').length !== 2
                  ).map((panel) => panel.getAttribute("data-artifact-id"));
                  const methodBlocks = [...reader.querySelectorAll(
                    '.analytics-layout-item[data-layout-block-id="method_flow"]'
                  )].filter(visible);
                  const visibleHeading = [...reader.querySelectorAll("h1")].find(visible);
                  const ids = (elements) => elements
                    .map((element) => element.getAttribute("data-artifact-id"))
                    .sort();
                  return {
                    clientWidth: root.clientWidth,
                    scrollWidth: root.scrollWidth,
                    horizontalOverflowPx: root.scrollWidth - root.clientWidth,
                    heading: visibleHeading?.textContent?.trim(),
                    liveReaderVisible: visible(reader),
                    fallbackVisible: visible(document.querySelector(
                      "#data-analytics-portable-fallback"
                    )),
                    chartPanels: chartPanels.length,
                    liveChartIds: ids(chartPanels),
                    staticChartPanels: staticChartPanels.length,
                    staticChartIds: ids(staticChartPanels),
                    tablePanels: tablePanels.length,
                    liveTableIds: ids(tablePanels),
                    staticTablePanels: staticTablePanels.length,
                    staticTableIds: ids(staticTablePanels),
                    svgVariants: svgVariants.length,
                    brokenCharts,
                    methodBlock: methodBlocks.length === 1,
                    distanceChart: Boolean(reader.querySelector(
                      '.chart-panel[data-artifact-id="distance_w_qc_chart"]'
                    )),
                    distanceCutImpactChart: Boolean(reader.querySelector(
                      '.chart-panel[data-artifact-id="distance_cut_impact_chart"]'
                    )),
                  };
                }"""
            )
            # The report contains a 154-row table and may exceed Chromium's
            # maximum single-image height. Capture the deterministic viewport
            # plus targeted evidence cards instead of an unreliable full-page
            # bitmap.
            screenshot_path = output_dir / f"{label}_overview.png"
            page.screenshot(path=str(screenshot_path), full_page=False)
            if label == "desktop":
                method_path = output_dir / "desktop_method_flow.png"
                distance_path = output_dir / "desktop_distance_qc.png"
                distance_cut_path = output_dir / "desktop_distance_cut_impact.png"
                page.locator(
                    '.analytics-layout-item[data-layout-block-id="method_flow"]'
                ).screenshot(path=str(method_path))
                page.locator(
                    '.chart-panel[data-artifact-id="distance_w_qc_chart"]'
                ).screenshot(path=str(distance_path))
                page.locator(
                    '.chart-panel[data-artifact-id="distance_cut_impact_chart"]'
                ).screenshot(path=str(distance_cut_path))

            checks = {
                "title_matches": metrics["heading"] == expected_title,
                "no_horizontal_overflow": metrics["horizontalOverflowPx"] <= 1,
                "live_reader_visible": metrics["liveReaderVisible"],
                "fallback_hidden_after_hydration": not metrics["fallbackVisible"],
                "chart_count_matches": metrics["chartPanels"] == expected_charts,
                "live_chart_ids_match": (
                    metrics["liveChartIds"] == expected_chart_ids
                ),
                "table_count_matches": metrics["tablePanels"] == expected_tables,
                "live_table_ids_match": (
                    metrics["liveTableIds"] == expected_table_ids
                ),
                "static_chart_count_matches": (
                    metrics["staticChartPanels"] == expected_charts
                ),
                "static_chart_ids_match": (
                    metrics["staticChartIds"] == expected_chart_ids
                ),
                "static_table_count_matches": (
                    metrics["staticTablePanels"] == expected_tables
                ),
                "static_table_ids_match": (
                    metrics["staticTableIds"] == expected_table_ids
                ),
                "svg_variant_count_matches": (
                    metrics["svgVariants"] == expected_svg_variants
                ),
                "every_chart_has_light_and_dark_svg": not metrics["brokenCharts"],
                "method_flow_present": metrics["methodBlock"],
                "distance_chart_present": metrics["distanceChart"],
                "distance_cut_impact_chart_present": metrics[
                    "distanceCutImpactChart"
                ],
                "console_errors_zero": not console_errors,
                "page_errors_zero": not page_errors,
            }
            results.append(
                {
                    "mode": "javascript_enabled",
                    "viewport": {
                        "label": label,
                        "width": width,
                        "height": height,
                    },
                    "all_pass": all(checks.values()),
                    "checks": checks,
                    "metrics": metrics,
                    "console_errors": console_errors,
                    "page_errors": page_errors,
                    "screenshot": str(screenshot_path),
                }
            )
            context.close()

        for label, width, height in VIEWPORTS:
            context = browser.new_context(
                viewport={"width": width, "height": height},
                color_scheme="light",
                reduced_motion="reduce",
                device_scale_factor=1,
                java_script_enabled=False,
            )
            page = context.new_page()
            console_errors: list[str] = []
            page_errors: list[str] = []
            page.on(
                "console",
                lambda message: console_errors.append(message.text)
                if message.type == "error"
                else None,
            )
            page.on("pageerror", lambda error: page_errors.append(str(error)))
            page.goto(html_path.as_uri(), wait_until="load")
            page.locator("#data-analytics-portable-fallback").wait_for(
                state="visible",
                timeout=30_000,
            )

            metrics = page.evaluate(
                """() => {
                  const root = document.documentElement;
                  const visible = (element) => {
                    if (!element) return false;
                    const style = getComputedStyle(element);
                    const rect = element.getBoundingClientRect();
                    return style.display !== "none" && style.visibility !== "hidden" &&
                      rect.width > 0 && rect.height > 0;
                  };
                  const fallback = document.querySelector(
                    "#data-analytics-portable-fallback"
                  );
                  const reader = document.querySelector(
                    "#data-analytics-portable-reader-root"
                  );
                  const chartPanels = [...fallback.querySelectorAll(
                    '[data-artifact-kind="chart"][data-artifact-id]'
                  )].filter(visible);
                  const tablePanels = [...fallback.querySelectorAll(
                    '[data-artifact-kind="table"][data-artifact-id]'
                  )].filter(visible);
                  const ids = (elements) => elements
                    .map((element) => element.getAttribute("data-artifact-id"))
                    .sort();
                  const brokenCharts = chartPanels.filter((panel) =>
                    panel.querySelectorAll('.portable-static-chart-svg').length !== 2
                  ).map((panel) => panel.getAttribute("data-artifact-id"));
                  const visibleHeading = [...fallback.querySelectorAll("h1")].find(visible);
                  return {
                    clientWidth: root.clientWidth,
                    scrollWidth: root.scrollWidth,
                    horizontalOverflowPx: root.scrollWidth - root.clientWidth,
                    heading: visibleHeading?.textContent?.trim(),
                    fallbackVisible: visible(fallback),
                    liveReaderAbsent: reader === null,
                    chartPanels: chartPanels.length,
                    chartIds: ids(chartPanels),
                    tablePanels: tablePanels.length,
                    tableIds: ids(tablePanels),
                    svgVariants: fallback.querySelectorAll(
                      '.portable-static-chart-svg'
                    ).length,
                    brokenCharts,
                  };
                }"""
            )
            screenshot_path = output_dir / f"nojs_{label}_overview.png"
            page.screenshot(path=str(screenshot_path), full_page=False)
            checks = {
                "title_matches": metrics["heading"] == expected_title,
                "no_horizontal_overflow": metrics["horizontalOverflowPx"] <= 1,
                "fallback_visible": metrics["fallbackVisible"],
                "live_reader_absent": metrics["liveReaderAbsent"],
                "chart_count_matches": metrics["chartPanels"] == expected_charts,
                "chart_ids_match": metrics["chartIds"] == expected_chart_ids,
                "table_count_matches": metrics["tablePanels"] == expected_tables,
                "table_ids_match": metrics["tableIds"] == expected_table_ids,
                "svg_variant_count_matches": (
                    metrics["svgVariants"] == expected_svg_variants
                ),
                "every_chart_has_light_and_dark_svg": not metrics["brokenCharts"],
                "console_errors_zero": not console_errors,
                "page_errors_zero": not page_errors,
            }
            results.append(
                {
                    "mode": "javascript_disabled",
                    "viewport": {
                        "label": label,
                        "width": width,
                        "height": height,
                    },
                    "all_pass": all(checks.values()),
                    "checks": checks,
                    "metrics": metrics,
                    "console_errors": console_errors,
                    "page_errors": page_errors,
                    "screenshot": str(screenshot_path),
                }
            )
            context.close()
        browser.close()

    evidence_files = {
        path.name: file_identity(path)
        for path in sorted(output_dir.glob("*.png"))
    }
    evidence_files_exact = set(evidence_files) == EXPECTED_EVIDENCE_FILES
    receipt = {
        "schema_name": "intersubmod.strict_region_html_visual_qa",
        "schema_version": "2.0.0",
        "created_at_utc": utc_now(),
        "all_pass": (
            all(result["all_pass"] for result in results)
            and evidence_files_exact
        ),
        "inputs": {
            "html": file_identity(html_path),
            "artifact": file_identity(artifact_path),
        },
        "browser": {
            "engine": "chromium",
            "version": browser_version,
            "executable": browser_identity,
        },
        "expected": {
            "title": expected_title,
            "visible_charts": expected_charts,
            "visible_tables": expected_tables,
            "static_svg_variants": expected_svg_variants,
            "chart_ids": expected_chart_ids,
            "table_ids": expected_table_ids,
            "modes": ["javascript_enabled", "javascript_disabled"],
        },
        "evidence": {
            "all_expected_files_present": evidence_files_exact,
            "files": evidence_files,
        },
        "viewports": results,
    }
    write_json(output_dir / "visual_qa_receipt.json", receipt)
    print(
        json.dumps(
            {
                "all_pass": receipt["all_pass"],
                "viewports": [
                    {
                        "mode": result["mode"],
                        "label": result["viewport"]["label"],
                        "horizontal_overflow_px": result["metrics"][
                            "horizontalOverflowPx"
                        ],
                        "charts": result["metrics"]["chartPanels"],
                        "tables": result["metrics"]["tablePanels"],
                        "svg_variants": result["metrics"].get("svgVariants"),
                    }
                    for result in results
                ],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
